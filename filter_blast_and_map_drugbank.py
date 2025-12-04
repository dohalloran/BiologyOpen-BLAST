#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Filter BLAST hits of MajorityVote_targets.fasta by e-value and "
            "map passing sequences back to DrugBank IDs using DrugBank_Targets.csv"
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
        default="DrugBank_Targets.csv",
        help="DrugBank targets CSV with a 'UniProt ID' and 'Drug IDs' column (default: %(default)s)",
    )
    p.add_argument(
        "--evalue_cutoff",
        type=float,
        default=1e-4,
        help="E-value cutoff to consider a query as having a hit (default: %(default)g)",
    )
    p.add_argument(
        "--out_fasta",
        default="MajorityVote_targets_filtered.fasta",
        help="Output FASTA of query sequences that passed the e-value cutoff (default: %(default)s)",
    )
    p.add_argument(
        "--out_csv",
        default="MajorityVote_targets_filtered_mapping.csv",
        help="Output CSV mapping UniProt ID -> DrugBank ID for passing sequences (default: %(default)s)",
    )
    return p.parse_args()


def load_passing_queries(blast_path, evalue_cutoff):
    """
    Read BLAST outfmt 6 and return:
      - dict: qseqid -> best_evalue (only if best_evalue <= cutoff)
    Assumes evalue is the 11th column (index 10) in standard outfmt 6.
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
      - id_to_header: qseqid -> full header (without leading '>')
      - id_to_seq:    qseqid -> sequence (no line breaks)

    Here qseqid is taken as the header up to the first whitespace.
    """
    id_to_header = {}
    id_to_seq = {}
    current_id = None
    current_header = None
    seq_chunks = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # flush previous
                if current_id is not None:
                    id_to_header[current_id] = current_header
                    id_to_seq[current_id] = "".join(seq_chunks)
                header = line[1:].strip()
                qid = header.split()[0]
                current_id = qid
                current_header = header
                seq_chunks = []
            else:
                if current_id is not None:
                    seq_chunks.append(line.strip())

    if current_id is not None:
        id_to_header[current_id] = current_header
        id_to_seq[current_id] = "".join(seq_chunks)

    return id_to_header, id_to_seq


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
    token = header.split()[0]  # first token
    if "|" in token:
        return token.split("|", 1)[1]
    return token


def build_uniprot_to_drugbank(targets_csv_path):
    """
    Build UniProt ID -> set of DrugBank IDs from DrugBank_Targets.csv.
    Assumes columns 'UniProt ID' and 'Drug IDs' exist.
    'Drug IDs' can contain multiple IDs separated by ';'.
    """
    df = pd.read_csv(targets_csv_path)

    possible_uni_cols = ["UniProt ID", "Uniprot ID", "UniProtID", "UniProt"]
    uni_col = None
    for c in possible_uni_cols:
        if c in df.columns:
            uni_col = c
            break
    if uni_col is None:
        raise ValueError(
            f"Could not find a UniProt ID column in {targets_csv_path}. "
            f"Looked for: {possible_uni_cols}"
        )

    possible_drug_cols = ["Drug IDs", "Drug IDs ", "Drugs", "DrugIDs"]
    drug_col = None
    for c in possible_drug_cols:
        if c in df.columns:
            drug_col = c
            break
    if drug_col is None:
        raise ValueError(
            f"Could not find a 'Drug IDs' column in {targets_csv_path}. "
            f"Looked for: {possible_drug_cols}"
        )

    mapping = defaultdict(set)
    for _, row in df.iterrows():
        uni = str(row[uni_col]).strip()
        if not uni or uni.lower() == "nan":
            continue
        drugs_raw = str(row[drug_col])
        if not drugs_raw or drugs_raw.lower() == "nan":
            continue
        for dbid in drugs_raw.split(";"):
            dbid = dbid.strip()
            if dbid:
                mapping[uni].add(dbid)

    return mapping


def main():
    args = parse_args()

    print(f"Reading BLAST results from: {args.blast}")
    passing_queries = load_passing_queries(args.blast, args.evalue_cutoff)
    print(f"  Queries with hits at e <= {args.evalue_cutoff:g}: {len(passing_queries)}")

    print(f"Reading FASTA queries from: {args.fasta}")
    id_to_header, id_to_seq = parse_fasta_ids(args.fasta)
    print(f"  Total sequences in FASTA: {len(id_to_header)}")

    print(f"Reading DrugBank targets from: {args.targets_csv}")
    uni_to_drugbank = build_uniprot_to_drugbank(args.targets_csv)
    print(f"  Unique UniProt IDs in DrugBank targets: {len(uni_to_drugbank)}")

    print(f"Writing filtered FASTA to: {args.out_fasta}")
    with open(args.out_fasta, "w") as fasta_out, open(args.out_csv, "w", newline="") as csv_out:
        csv_writer = csv.writer(csv_out)
        csv_writer.writerow(
            [
                "QueryID",
                "BestEvalue",
                "FASTA_Header",
                "UniProt_ID",
                "DrugBank_IDs",
            ]
        )

        kept = 0
        for qid, best_e in sorted(passing_queries.items(), key=lambda x: x[0]):
            if qid not in id_to_header:
                continue

            header = id_to_header[qid]
            seq = id_to_seq[qid]
            uni = extract_uniprot_from_header(header)
            db_ids = sorted(uni_to_drugbank.get(uni, []))

            fasta_out.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fasta_out.write(seq[i : i + 60] + "\n")

            csv_writer.writerow(
                [qid, f"{best_e:.3g}", header, uni, ";".join(db_ids)]
            )
            kept += 1

    print("Done.")
    print(f"  Filtered sequences written: {kept}")
    print(f"  Mapping CSV written to: {args.out_csv}")


if __name__ == "__main__":
    main()
