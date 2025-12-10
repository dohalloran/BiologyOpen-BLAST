#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Filter DrugBank compounds to those whose targets (UniProt_ID) "
            "have BLAST hits against Ancylostoma."
        )
    )
    p.add_argument(
        "--blast",
        default="MajorityVote_targets_vs_ancylostoma.tsv",
        help="BLAST tabular (outfmt 6) results file (default: %(default)s)",
    )
    p.add_argument(
        "--fasta",
        default="MajorityVote_targets.fasta",
        help="FASTA file used as BLAST queries (default: %(default)s)",
    )
    p.add_argument(
        "--targets_csv",
        default="MajorityVote_DrugBank_to_UniProt.csv",
        help=(
            "CSV with DrugBank targets. Must contain columns "
            "'UniProt_ID' and 'DrugBank_ID' (default: %(default)s)"
        ),
    )
    p.add_argument(
        "--evalue_cutoff",
        type=float,
        default=1e-4,
        help="E-value cutoff to consider a query as having a hit (default: %(default)g)",
    )
    p.add_argument(
        "--out_drugs",
        default="MajorityVote_filtered_drugs.csv",
        help="Output CSV with one DrugBank_ID per row (default: %(default)s)",
    )
    p.add_argument(
        "--out_mapping",
        default="MajorityVote_filtered_drug_to_uniprot.csv",
        help=(
            "Output CSV mapping DrugBank_ID -> UniProt_IDs with BLAST hits "
            "(default: %(default)s)"
        ),
    )
    return p.parse_args()


# ------------------- BLAST / FASTA helpers ------------------- #

def load_passing_queries(blast_path, evalue_cutoff):
    """
    Read BLAST outfmt 6 and return:
      dict: qseqid -> best_evalue (only if best_evalue <= cutoff).
    Assumes evalue is the 11th column (index 10).
    """
    passing = {}
    with open(blast_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            qseqid = parts[0]
            try:
                evalue = float(parts[10])
            except ValueError:
                continue

            if evalue <= evalue_cutoff:
                if qseqid not in passing or evalue < passing[qseqid]:
                    passing[qseqid] = evalue
    return passing


def parse_fasta_ids(fasta_path):
    """
    Parse FASTA and return:
      id_to_header: qseqid -> full header (without '>')
    qseqid is taken as the header up to the first whitespace.
    """
    id_to_header = {}
    current_id = None
    current_header = None

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:].strip()
                qid = header.split()[0]
                current_id = qid
                current_header = header
                id_to_header[current_id] = current_header
            else:
                # sequence itself is not needed for this filtering step
                continue

    return id_to_header


def extract_uniprot_from_header(header):
    """
    Extract UniProt ID from a FASTA header string.

    Examples:
      'drugbank_target|P45059 Peptidoglycan synthase FtsI (DB00303)'
        -> 'P45059'
      'P19113|something'
        -> 'P19113'
      'Q9UI32'
        -> 'Q9UI32'
    """
    token = header.split()[0]
    if "|" in token:
        return token.split("|", 1)[1]
    return token


# ------------------- DrugBank / UniProt mapping ------------------- #

def build_uniprot_to_drugbank(targets_csv_path):
    """
    Build:
      uni_to_drugs: UniProt_ID -> set of DrugBank_IDs
      all_drugs:    set of all DrugBank_IDs in the file

    Requires columns 'UniProt_ID' and 'DrugBank_ID'. If 'DrugBank_ID'
    contains multiple IDs separated by ';', they are all captured.
    """
    df = pd.read_csv(targets_csv_path)

    if "UniProt_ID" not in df.columns or "DrugBank_ID" not in df.columns:
        raise ValueError(
            f"{targets_csv_path} must contain 'UniProt_ID' and 'DrugBank_ID' columns."
        )

    uni_to_drugs = defaultdict(set)
    all_drugs = set()

    for _, row in df.iterrows():
        uni = str(row["UniProt_ID"]).strip()
        if not uni or uni.lower() == "nan":
            continue

        raw = str(row["DrugBank_ID"])
        if not raw or raw.lower() == "nan":
            continue

        for dbid in raw.split(";"):
            dbid = dbid.strip()
            if not dbid:
                continue
            uni_to_drugs[uni].add(dbid)
            all_drugs.add(dbid)

    return uni_to_drugs, all_drugs


# ------------------- main ------------------- #

def main():
    args = parse_args()

    print(f"Reading DrugBank targets from: {args.targets_csv}")
    uni_to_drugs, all_drugs = build_uniprot_to_drugbank(args.targets_csv)
    print(f"  Unique UniProt_IDs in targets CSV:   {len(uni_to_drugs)}")
    print(f"  Unique DrugBank_IDs in targets CSV:  {len(all_drugs)}")

    print(f"\nReading BLAST results from: {args.blast}")
    passing_queries = load_passing_queries(args.blast, args.evalue_cutoff)
    print(f"  Queries with hits at e <= {args.evalue_cutoff:g}: {len(passing_queries)}")

    print(f"\nReading FASTA headers from: {args.fasta}")
    id_to_header = parse_fasta_ids(args.fasta)
    print(f"  Total query sequences in FASTA: {len(id_to_header)}")

    # Map passing qseqid -> UniProt_ID
    passing_unis = set()
    missing_in_fasta = 0
    for qid in passing_queries:
        header = id_to_header.get(qid)
        if header is None:
            missing_in_fasta += 1
            continue
        uni = extract_uniprot_from_header(header)
        passing_unis.add(uni)

    print(f"\n  Unique UniProt_IDs with BLAST hits:  {len(passing_unis)}")
    if missing_in_fasta:
        print(f"  WARNING: {missing_in_fasta} passing qseqid(s) not found in FASTA headers.")

    # Collect DrugBank_IDs for those UniProt targets
    filtered_drugs = set()
    drug_to_unis = defaultdict(set)

    for uni in passing_unis:
        for dbid in uni_to_drugs.get(uni, []):
            filtered_drugs.add(dbid)
            drug_to_unis[dbid].add(uni)

    print(f"\nUnique DrugBank_IDs after filtering by BLAST: {len(filtered_drugs)}")
    print(f"Reduction from original {len(all_drugs)} -> {len(filtered_drugs)}")

    # 1) Write one-drug-per-row CSV
    print(f"\nWriting filtered drug list to: {args.out_drugs}")
    with open(args.out_drugs, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["DrugBank_ID"])
        for dbid in sorted(filtered_drugs):
            w.writerow([dbid])

    # 2) Write mapping DrugBank_ID -> UniProt_IDs_with_hits
    print(f"Writing DrugBank_ID -> UniProt_IDs mapping to: {args.out_mapping}")
    with open(args.out_mapping, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(["DrugBank_ID", "UniProt_IDs_with_hits"])
        for dbid in sorted(filtered_drugs):
            unis = sorted(drug_to_unis[dbid])
            w.writerow([dbid, ";".join(unis)])

    print("\nDone.")


if __name__ == "__main__":
    main()
