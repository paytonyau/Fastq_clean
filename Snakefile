# ==============================================================================
# Fastq_clean (modernised) — Snakemake workflow
#
# Reimplements the three conceptual stages of the original 2014 Fastq_clean
# Perl pipeline (illumina_clean.pl / sRNA_clean.pl) with current, actively
# maintained tools:
#
#   1. adapter / quality trim   : illumina_clean.pl pattern-matching  -> fastp
#   2. rRNA / virus filtering   : bin/bwa (bundled, uncompiled binary) -> bwa mem (apt/bioconda) + samtools
#   3. stats / clean report     : tools/qcStatistics.pl + Fq_statistics.R -> fastp JSON + samtools flagstat,
#                                  aggregated by workflow/scripts/summarize_stats.py
#
# Usage:
#   snakemake -c4 --configfile config/config.yaml
#   snakemake -c4 --configfile config/config.yaml -n     # dry run
# ==============================================================================

import pandas as pd

configfile: "config/config.yaml"

OUTDIR = config["outdir"]
samples_df = pd.read_csv(config["samples"], sep="\t", dtype=str).set_index("sample_id", drop=False)
SAMPLES = list(samples_df.index)
CONTAM_ENABLED = config["contaminant_filter"]["enabled"]

def is_paired(sample):
    r2 = samples_df.loc[sample, "fastq_2"]
    return isinstance(r2, str) and r2.strip() != "" and r2.lower() != "nan"

def fq1(sample):
    return samples_df.loc[sample, "fastq_1"]

def fq2(sample):
    return samples_df.loc[sample, "fastq_2"]

rule all:
    input:
        f"{OUTDIR}/{config['stats']['summary_tsv']}",
        f"{OUTDIR}/{config['stats']['summary_html']}"

# ------------------------------------------------------------------ #
# Stage 1: adapter / quality trimming (fastp)
#   One rule, branches PE/SE internally -> avoids ambiguous-rule clashes.
# ------------------------------------------------------------------ #

rule fastp_trim:
    input:
        r1 = lambda wc: fq1(wc.sample),
        r2 = lambda wc: fq2(wc.sample) if is_paired(wc.sample) else []
    output:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        json = f"{OUTDIR}/trimmed/{{sample}}.fastp.json",
        html = f"{OUTDIR}/trimmed/{{sample}}.fastp.html"
    log:
        f"{OUTDIR}/logs/{{sample}}.fastp.log"
    params:
        q = config["trimming"]["quality_phred"],
        n = config["trimming"]["n_base_limit"],
        l = config["trimming"]["length_required"],
        detect_adapter = "--detect_adapter_for_pe" if config["trimming"]["detect_adapter_for_pe"] else "",
        extra = config["trimming"]["extra_args"],
        paired = lambda wc: is_paired(wc.sample)
    threads: 4
    run:
        common = (
            f"--qualified_quality_phred {params.q} "
            f"--n_base_limit {params.n} "
            f"--length_required {params.l} "
            f"{params.extra} --thread {threads} "
            f"--json {output.json} --html {output.html}"
        )
        if params.paired:
            shell(
                "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
                "{params.detect_adapter} " + common + " > {log} 2>&1"
            )
        else:
            # single-end: still declare an (empty) R2 output so downstream rules
            # that expect two files per sample slot always find one
            shell(
                "fastp -i {input.r1} -o {output.r1} " + common + " > {log} 2>&1 "
                "&& : > {output.r2}"
            )

# ------------------------------------------------------------------ #
# Stage 2: contaminant (rRNA / virus) filtering via bwa mem + samtools
#   Original pipeline aligned trimmed reads against rRNA/virus references
#   and discarded any read that mapped. Same logic here: align, keep only
#   reads where the mate pair (or single read) is UNMAPPED (-f 12/-f 4).
# ------------------------------------------------------------------ #

rule bwa_index_contaminant:
    input:
        config["contaminant_filter"]["reference_fasta"]
    output:
        config["contaminant_filter"]["reference_fasta"] + ".bwt"
    log:
        f"{OUTDIR}/logs/bwa_index.log"
    shell:
        "bwa index {input} > {log} 2>&1"

rule filter_contaminants:
    input:
        r1 = f"{OUTDIR}/trimmed/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{OUTDIR}/trimmed/{{sample}}_R2.trimmed.fastq.gz",
        ref = config["contaminant_filter"]["reference_fasta"],
        idx = config["contaminant_filter"]["reference_fasta"] + ".bwt"
    output:
        r1 = f"{OUTDIR}/clean/{{sample}}_R1.clean.fastq.gz",
        r2 = f"{OUTDIR}/clean/{{sample}}_R2.clean.fastq.gz",
        bam = f"{OUTDIR}/contaminant_screen/{{sample}}.bam",
        flagstat = f"{OUTDIR}/contaminant_screen/{{sample}}.flagstat.txt"
    log:
        f"{OUTDIR}/logs/{{sample}}.bwa.log"
    params:
        threads = config["contaminant_filter"]["bwa_threads"],
        paired = lambda wc: is_paired(wc.sample)
    threads: 4
    run:
        if params.paired:
            shell(
                "bwa mem -t {params.threads} {input.ref} {input.r1} {input.r2} 2> {log} "
                "| samtools sort -@ {params.threads} -o {output.bam} - "
                "&& samtools index {output.bam} "
                "&& samtools flagstat {output.bam} > {output.flagstat} "
                "&& samtools fastq -f 12 -F 256 "
                "-1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {output.bam} "
                ">> {log} 2>&1"
            )
        else:
            shell(
                "bwa mem -t {params.threads} {input.ref} {input.r1} 2> {log} "
                "| samtools sort -@ {params.threads} -o {output.bam} - "
                "&& samtools index {output.bam} "
                "&& samtools flagstat {output.bam} > {output.flagstat} "
                "&& samtools fastq -f 4 -n {output.bam} 2>> {log} | gzip > {output.r1} "
                "&& : > {output.r2}"
            )

# ------------------------------------------------------------------ #
# Stage 3: per-sample stats -> aggregated summary
# ------------------------------------------------------------------ #

def all_fastp_json():
    return [f"{OUTDIR}/trimmed/{s}.fastp.json" for s in SAMPLES]

def all_flagstat():
    if not CONTAM_ENABLED:
        return []
    return [f"{OUTDIR}/contaminant_screen/{s}.flagstat.txt" for s in SAMPLES]

def all_clean_r1():
    return [f"{OUTDIR}/clean/{s}_R1.clean.fastq.gz" if CONTAM_ENABLED
            else f"{OUTDIR}/trimmed/{s}_R1.trimmed.fastq.gz" for s in SAMPLES]

rule summarize:
    input:
        fastp_json = all_fastp_json(),
        flagstat = all_flagstat(),
        clean = all_clean_r1()
    output:
        tsv = f"{OUTDIR}/{config['stats']['summary_tsv']}",
        html = f"{OUTDIR}/{config['stats']['summary_html']}"
    params:
        samples = ",".join(SAMPLES)
    shell:
        """
        python workflow/scripts/summarize_stats.py \
            --samples {params.samples} \
            --outdir {OUTDIR} \
            --tsv {output.tsv} \
            --html {output.html}
        """
