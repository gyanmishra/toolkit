#!/usr/bin/env python3

"""Aggregate per-sample Qualimap RNA-seq reports into a TSV."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_report(report_path: Path) -> dict[str, str]:
    metrics: dict[str, str] = {}
    with report_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if "=" not in line:
                continue
            key, value = line.strip().split("=", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate a combined TSV from Qualimap RNA-seq reports.")
    parser.add_argument("-d", "--directory", required=True, help="Qualimap parent directory")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    qualimap_dir = Path(args.directory).resolve()

    if not qualimap_dir.is_dir():
        raise FileNotFoundError(f"Qualimap directory not found: {qualimap_dir}")

    report_dict: dict[str, dict[str, str]] = {}
    for child in sorted(qualimap_dir.iterdir()):
        if not child.is_dir():
            continue
        report_file = child / "rnaseq_qc_results.txt"
        if not report_file.exists():
            continue
        report_dict[child.name] = parse_report(report_file)

    if not report_dict:
        raise ValueError(f"No rnaseq_qc_results.txt files found under: {qualimap_dir}")

    df = pd.DataFrame.from_dict(report_dict, orient="index")
    output_path = qualimap_dir / "Qualimap_report.tsv"
    df.to_csv(output_path, sep="\t")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
