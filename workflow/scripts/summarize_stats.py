#!/usr/bin/env python3
"""
Aggregate per-sample fastp JSON + samtools flagstat output into one summary
table. Replaces the original pipeline's tools/qcStatistics.pl,
tools/extract_clean_report.pl and tools/Fq_statistics.R.
"""
import argparse
import gzip
import json
import os
import sys


def read_fastp_json(path):
    with open(path) as fh:
        d = json.load(fh)
    before = d.get("summary", {}).get("before_filtering", {})
    after = d.get("summary", {}).get("after_filtering", {})
    filtering = d.get("filtering_result", {})
    return {
        "reads_before_trim": before.get("total_reads", 0),
        "bases_before_trim": before.get("total_bases", 0),
        "q20_rate_before": before.get("q20_rate", 0),
        "q30_rate_before": before.get("q30_rate", 0),
        "reads_after_trim": after.get("total_reads", 0),
        "bases_after_trim": after.get("total_bases", 0),
        "q20_rate_after": after.get("q20_rate", 0),
        "q30_rate_after": after.get("q30_rate", 0),
        "reads_passed_filter": filtering.get("passed_filter_reads", 0),
        "reads_low_quality": filtering.get("low_quality_reads", 0),
        "reads_too_many_N": filtering.get("too_many_N_reads", 0),
        "reads_too_short": filtering.get("too_short_reads", 0),
        "adapter_trimmed_reads": d.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
    }


def read_flagstat(path):
    total = mapped = 0
    with open(path) as fh:
        for line in fh:
            if " in total" in line:
                total = int(line.split()[0])
            elif " mapped (" in line and "primary mapped" not in line and "mate mapped" not in line:
                mapped = int(line.split()[0])
    contaminant_reads = mapped
    pct_contaminant = (100.0 * mapped / total) if total else 0.0
    return {
        "screened_reads_total": total,
        "contaminant_reads_removed": contaminant_reads,
        "pct_contaminant": round(pct_contaminant, 3),
    }


def count_clean_reads(fastq_gz_path):
    if not os.path.exists(fastq_gz_path) or os.path.getsize(fastq_gz_path) == 0:
        return 0
    n = 0
    opener = gzip.open if fastq_gz_path.endswith(".gz") else open
    with opener(fastq_gz_path, "rt") as fh:
        for i, _ in enumerate(fh):
            pass
    try:
        n = (i + 1) // 4
    except UnboundLocalError:
        n = 0
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, help="comma-separated sample IDs")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--html", required=True)
    args = ap.parse_args()

    samples = args.samples.split(",")
    rows = []
    for s in samples:
        row = {"sample": s}
        fastp_json = os.path.join(args.outdir, "trimmed", f"{s}.fastp.json")
        if os.path.exists(fastp_json):
            row.update(read_fastp_json(fastp_json))

        flagstat = os.path.join(args.outdir, "contaminant_screen", f"{s}.flagstat.txt")
        clean_path = os.path.join(args.outdir, "clean", f"{s}_R1.clean.fastq.gz")
        if os.path.exists(flagstat):
            row.update(read_flagstat(flagstat))
            row["final_clean_reads"] = count_clean_reads(clean_path)
        else:
            trimmed_path = os.path.join(args.outdir, "trimmed", f"{s}_R1.trimmed.fastq.gz")
            row["final_clean_reads"] = count_clean_reads(trimmed_path)

        rows.append(row)

    fieldnames = []
    for r in rows:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    with open(args.tsv, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in fieldnames) + "\n")

    html = ["<html><head><title>Fastq_clean summary</title>",
            "<style>table{border-collapse:collapse;font-family:sans-serif;font-size:13px}",
            "td,th{border:1px solid #ccc;padding:4px 8px;text-align:right}",
            "th{background:#f0f0f0}td:first-child,th:first-child{text-align:left}</style>",
            "</head><body><h2>Fastq_clean (modernised) — run summary</h2><table><tr>"]
    html.append("".join(f"<th>{f}</th>" for f in fieldnames))
    html.append("</tr>")
    for r in rows:
        html.append("<tr>" + "".join(f"<td>{r.get(k, '')}</td>" for k in fieldnames) + "</tr>")
    html.append("</table></body></html>")

    with open(args.html, "w") as fh:
        fh.write("\n".join(html))

    print(f"Wrote {args.tsv} and {args.html} for {len(rows)} sample(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
