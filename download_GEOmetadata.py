#!/usr/bin/env python3

"""Download GEO sample metadata and optional FASTQ files for selected GSM IDs.

Examples
--------
python download_GEOmetadata.py -GEO GSE144825 -GSM All -e you@example.com -key <API_KEY>
python download_GEOmetadata.py -GEO GSE144825 -GSM GSM4297142,GSM4297154 -d ./downloads -e you@example.com -key <API_KEY>
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import GEOparse
import pandas as pd
from Bio import Entrez

logging.getLogger("GEOparse").setLevel(logging.ERROR)


def fetch_geo_samples(geo_accession: str, gsm_ids: str) -> list[tuple[str, str]]:
    """Fetch sample names and GSM IDs for a GEO accession."""
    handle = Entrez.esearch(db="gds", term=geo_accession)
    record = Entrez.read(handle)
    handle.close()

    if not record.get("IdList"):
        raise ValueError(f"No GEO series found for accession: {geo_accession}")

    geo_id = record["IdList"][0]
    handle = Entrez.esummary(db="gds", id=geo_id)
    summary = Entrez.read(handle)
    handle.close()

    samples = summary[0].get("Samples", [])
    if "All" not in gsm_ids:
        target_ids = set(g.strip() for g in gsm_ids.split(",") if g.strip())
        samples = [item for item in samples if item.get("Accession") in target_ids]

    return [(sample["Title"], sample["Accession"]) for sample in samples]


def map_gsm_to_sra(gsm_id: str) -> str | None:
    """Map a GSM ID to corresponding SRA run accession."""
    handle = Entrez.esearch(db="sra", term=gsm_id)
    record = Entrez.read(handle)
    handle.close()

    sra_ids = record.get("IdList", [])
    if not sra_ids:
        return None

    handle = Entrez.efetch(db="sra", id=sra_ids[0], rettype="xml")
    xml_data = handle.read()
    handle.close()

    root = ET.fromstring(xml_data)
    for run in root.iter("RUN"):
        return run.attrib.get("accession")
    return None


def map_gsm_to_org(gsm_id: str, gse: GEOparse.GEO.GSE) -> str:
    """Return semicolon-joined organism field from GEO sample metadata."""
    sample = gse.gsms.get(gsm_id)
    if sample is None:
        return ""
    organisms = sample.metadata.get("organism_ch1", [])
    return ";".join(organisms)


def download_sra(sra_id: str, directory: Path) -> None:
    """Download FASTQ files for an SRA run if not already present."""
    paired_gz = (directory / f"{sra_id}_1.fastq.gz", directory / f"{sra_id}_2.fastq.gz")
    single_gz = directory / f"{sra_id}.fastq.gz"
    if (paired_gz[0].exists() and paired_gz[1].exists()) or single_gz.exists():
        return
    print(f"Downloading {sra_id} to {directory} ...")
    subprocess.run(["fasterq-dump", sra_id, "--outdir", str(directory)], check=True)


def maybe_gzip(path: Path) -> Path:
    """Gzip file if plain file exists; return gz path."""
    gz_path = path.with_suffix(path.suffix + ".gz")
    if gz_path.exists():
        return gz_path
    if path.exists():
        subprocess.run(["gzip", str(path)], check=True)
    return gz_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download SRA metadata from GEO accession and optionally download FASTQ files."
    )
    parser.add_argument("-GEO", "--GEOaccession", type=str, required=True, help="GEO accession (e.g. GSE144825)")
    parser.add_argument(
        "-GSM",
        "--GSMid",
        type=str,
        required=True,
        help='GSM IDs comma-separated, or "All"',
    )
    parser.add_argument("-d", "--outDir", type=str, help="Output directory to download FASTQ files")
    parser.add_argument("-e", "--EntrezEmail", type=str, required=True, help="Entrez email")
    parser.add_argument("-key", "--EntrezKey", type=str, required=True, help="Entrez API key")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    Entrez.email = args.EntrezEmail
    Entrez.api_key = args.EntrezKey

    out_dir = Path(args.outDir).resolve() if args.outDir else None
    if out_dir and not out_dir.is_dir():
        raise FileNotFoundError(f"Output directory does not exist: {out_dir}")

    geo_samples = fetch_geo_samples(args.GEOaccession, args.GSMid)
    if not geo_samples:
        raise ValueError(f"No matching GSM samples found for {args.GEOaccession} using {args.GSMid}")

    gse = GEOparse.get_GEO(geo=args.GEOaccession, destdir=".")

    data: list[dict[str, str]] = []
    for sample_name, gsm_id in geo_samples:
        sra_id = map_gsm_to_sra(gsm_id)
        organism = map_gsm_to_org(gsm_id, gse)

        row = {
            "SampleName": sample_name,
            "GSM_ID": gsm_id,
            "SRA_ID": sra_id or "",
            "Organism": organism,
            "Read_1": "",
            "Read_2": "",
        }

        if out_dir and sra_id:
            download_sra(sra_id, out_dir)
            fastq_1 = maybe_gzip(out_dir / f"{sra_id}_1.fastq")
            fastq_2 = maybe_gzip(out_dir / f"{sra_id}_2.fastq")
            single = maybe_gzip(out_dir / f"{sra_id}.fastq")

            if fastq_1.exists() and fastq_2.exists():
                row["Read_1"] = fastq_1.name
                row["Read_2"] = fastq_2.name
            elif single.exists():
                row["Read_1"] = single.name

        data.append(row)

    pd.DataFrame(data).to_csv("data.csv", index=False)


if __name__ == "__main__":
    main()
