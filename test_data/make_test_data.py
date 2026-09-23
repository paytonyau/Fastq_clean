#!/usr/bin/env python3
"""Generate small synthetic test data for the modernised Fastq_clean workflow:
   - a "contaminant" reference (stand-in for rRNA/virus references)
   - 2 paired-end samples, each a mix of:
       * clean reads (random, not from the contaminant reference)
       * contaminant reads (drawn from the reference -> should be filtered out)
       * low-quality / adapter-contaminated reads (should be trimmed/dropped)
"""
import gzip
import random

random.seed(42)
BASES = "ACGT"
ADAPTER = "AGATCGGAAGAGC"  # standard Illumina TruSeq adapter


def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def rand_qual(n, mean_q=35):
    # Phred+33 quality string, mostly high quality
    return "".join(chr(33 + max(2, min(40, int(random.gauss(mean_q, 4))))) for _ in range(n))


def revcomp(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def write_fastq(path, records):
    with gzip.open(path, "wt") as fh:
        for name, seq, qual in records:
            fh.write(f"@{name}\n{seq}\n+\n{qual}\n")


def make_sample(sample_id, n_clean=400, n_contam=100, n_lowqual=50, read_len=100):
    contam_ref = CONTAM_REF
    r1_records, r2_records = [], []

    # clean genomic-like reads (random sequence, high quality)
    for i in range(n_clean):
        frag = rand_seq(read_len + 50)
        r1 = frag[:read_len]
        r2 = revcomp(frag)[:read_len]
        name = f"{sample_id}:clean:{i}"
        r1_records.append((name, r1, rand_qual(read_len, 36)))
        r2_records.append((name, r2, rand_qual(read_len, 36)))

    # contaminant reads: sampled straight from the reference -> bwa should map these
    for i in range(n_contam):
        start = random.randint(0, len(contam_ref) - read_len - 50)
        frag = contam_ref[start:start + read_len + 50]
        r1 = frag[:read_len]
        r2 = revcomp(frag)[:read_len]
        name = f"{sample_id}:contam:{i}"
        r1_records.append((name, r1, rand_qual(read_len, 36)))
        r2_records.append((name, r2, rand_qual(read_len, 36)))

    # low-quality / adapter-through reads: short insert so adapter reads through
    # into the 3' end, plus poor quality tail -> fastp should trim/drop these
    for i in range(n_lowqual):
        insert_len = random.randint(20, 40)
        frag = rand_seq(insert_len)
        r1 = (frag + ADAPTER + rand_seq(read_len))[:read_len]
        r2 = (revcomp(frag) + ADAPTER + rand_seq(read_len))[:read_len]
        # degrade quality after the insert to simulate real adapter-readthrough reads
        q1 = rand_qual(insert_len, 36) + rand_qual(read_len - insert_len, 8)
        q2 = rand_qual(insert_len, 36) + rand_qual(read_len - insert_len, 8)
        name = f"{sample_id}:lowqual:{i}"
        r1_records.append((name, r1, q1))
        r2_records.append((name, r2, q2))

    random.shuffle(r1_records)
    # keep r2 order matching r1 by name
    order = {name: idx for idx, (name, _, _) in enumerate(r1_records)}
    r2_records.sort(key=lambda rec: order[rec[0]])

    write_fastq(f"test_data/{sample_id}_R1.fastq.gz", r1_records)
    write_fastq(f"test_data/{sample_id}_R2.fastq.gz", r2_records)
    print(f"{sample_id}: {len(r1_records)} read pairs "
          f"({n_clean} clean, {n_contam} contaminant, {n_lowqual} low-qual/adapter)")


CONTAM_REF = rand_seq(3000)  # stand-in for a concatenated rRNA+virus reference

if __name__ == "__main__":
    with open("test_data/contaminant_reference.fasta", "w") as fh:
        fh.write(">synthetic_rRNA_virus_reference\n")
        for i in range(0, len(CONTAM_REF), 70):
            fh.write(CONTAM_REF[i:i + 70] + "\n")

    make_sample("sample1", n_clean=400, n_contam=100, n_lowqual=50)
    make_sample("sample2", n_clean=350, n_contam=150, n_lowqual=60)
