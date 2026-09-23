# Fastq_clean (modernised)

A Snakemake reimplementation of the 2014 `Fastq_clean` Perl pipeline
(Zhang et al. 2014, BIBM), using current, actively-maintained tools in place
of the bundled/vendored ones from the original repo. **Tested end-to-end**
against synthetic data (see "Testing" below) — every stage's output count
matches the known composition of the test data exactly.

## Stage mapping (old → new)

| Stage | Original (2014) | Modernised |
|---|---|---|
| Adapter / quality trim | `illumina_clean.pl` pattern-matching + `tools/fastq_clipper.pl` / `bin/trim_ends.pl` | [`fastp`](https://github.com/OpenGene/fastp) |
| QC report | Bundled FastQC v0.11.2 (2014) | current `fastqc` (or skip — fastp's own HTML/JSON report already covers this) |
| rRNA / virus filtering | `bin/bwa` (uncompiled binary checked into the repo) aligned to `rRNA_reference` / `virus_reference` | `bwa mem` (current, via apt/conda) + `samtools` — same align-and-discard logic |
| SRA download | `tools/fastq-dump` (vendored binary) | current `sra-tools` (`prefetch` / `fasterq-dump` via conda) — not wired into this workflow, add if you pull directly from SRA |
| Stats / clean report | `tools/qcStatistics.pl`, `tools/extract_clean_report.pl`, `tools/Fq_statistics.R` | `workflow/scripts/summarize_stats.py` (parses fastp JSON + samtools flagstat into one TSV/HTML) |
| Orchestration | shell scripts + `conf/*.config` files, no resumability | Snakemake DAG — automatic parallelism, resumable, dry-run support |

## Parameter mapping

The original `illumina_clean.pl` defaults (tuned for 100 bp Illumina reads)
are carried over as fastp equivalents in `config/config.yaml`:

| Original | Value | fastp equivalent |
|---|---|---|
| `quality_Cutoff` | 20 | `--qualified_quality_phred 20` |
| `n_Cutoff` | 2 | `--n_base_limit 2` |
| `read_Length` (min kept) | 25 | `--length_required 25` |

The original's `bwa aln`-era parameters (`max_dist`, `max_open`,
`max_extension`, `len_seed`, `dist_seed`) don't map 1:1 onto `bwa mem`
(a different, newer alignment algorithm with its own well-tuned defaults for
short-read contamination screening) — `bwa mem`'s defaults are used as-is
rather than forcing the old seed/mismatch settings through.

## Requirements

```bash
# Ubuntu/Debian (what this was tested with)
sudo apt install fastp bwa samtools
pip install snakemake pandas

# or, reproducibly, via conda/mamba:
mamba env create -f envs/environment.yml
conda activate fastq_clean_modern
```

## Usage

1. Edit `config/samples.tsv` — one row per sample (`sample_id`, `fastq_1`,
   `fastq_2`; leave `fastq_2` blank for single-end data).
2. Point `contaminant_filter.reference_fasta` in `config/config.yaml` at your
   own rRNA + virus reference FASTA (concatenate multiple references into one
   file — that's what the original pipeline's `rRNA_reference` /
   `virus_reference` did too).
3. Dry-run to check the plan, then run for real:

```bash
snakemake -c4 --configfile config/config.yaml -n   # dry run
snakemake -c4 --configfile config/config.yaml       # actual run
```

Outputs land in `results/`:
- `results/trimmed/` — fastp-trimmed reads + per-sample JSON/HTML reports
- `results/clean/` — final reads with contaminants removed
- `results/contaminant_screen/` — BAMs + flagstat from the contaminant screen
- `results/summary.tsv` / `results/summary.html` — aggregated stats across all samples

## Testing

`test_data/make_test_data.py` generates two synthetic paired-end samples with
a known composition (random "clean" reads, reads drawn verbatim from a
synthetic contaminant reference, and reads with adapter read-through +
degraded quality tails). Regenerate and re-run with:

```bash
python3 test_data/make_test_data.py
snakemake -c4 --configfile config/config.yaml --forceall
```

On the bundled test data, every stage's output count matched the injected
composition exactly (e.g. sample1: 550 injected pairs → 1054/2=527 pairs
survive trimming → 427 remain after removing exactly the 100 injected
contaminant pairs). That's a synthetic sanity check, not a substitute for
running against your own real reference and reads.

## Notes / things to adapt for your own data

- The synthetic contaminant reference here is random sequence, not a real
  rRNA/virus database — swap in SILVA/RefSeq rRNA + a virus reference for
  actual use.
- Poly-G trimming (`--trim_poly_g`) is off by default; turn it on via
  `trimming.extra_args` in the config if you're on NovaSeq/NextSeq
  (2-colour chemistry can produce spurious high-quality G-runs at read ends).
- `sRNA_clean.pl`'s small-RNA-specific logic (short, single adapter,
  narrow length window) isn't reproduced here — for that, use fastp with
  `--length_required` set low (e.g. 15–18) and no polyG/overlap trimming, or
  use `cutadapt` directly, which is more common for sRNA-seq specifically.
