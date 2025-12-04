#!/usr/bin/env python3
"""
Extract UniProt targets and protein sequences for DrugBank compounds
predicted active by MajorityVote_Class == 1.

Inputs (default filenames):
  - DrugBank_approved_predictions_all_models.csv
  - DrugBank_Targets.csv
  - DrugBank_protein.fasta

Outputs:
  - MajorityVote_DrugBank_to_UniProt.csv
  - MajorityVote_targets.fasta
"""

import argparse
import sys
from typing import Dict, List, Set, Tuple

import pandas as pd


def load_active_drugbank_ids(predictions_path: str) -> Tuple[str, List[str]]:
    """
    Load predictions file and return:
      (drugbank_id_column_name, list_of_active_drugbank_ids)
    Active = MajorityVote_Class == 1
    """
    print(f"Loading predictions from: {predictions_path}")
    df = pd.read_csv(predictions_path, sep=None, engine="python")

    if "MajorityVote_Class" not in df.columns:
        raise KeyError(
            "Column 'MajorityVote_Class' not found in predictions file. "
            f"Available columns: {list(df.columns)}"
        )

    # Try to find the DrugBank ID column (e.g. 'DrugBank ID')
    dbid_candidates = [c for c in df.columns if "DrugBank" in c and "ID" in c]
    if not dbid_candidates:
        raise KeyError(
            "Could not find a DrugBank ID column in predictions file. "
            f"Available columns: {list(df.columns)}"
        )
    drugbank_col = dbid_candidates[0]

    active_ids = (
        df.loc[df["MajorityVote_Class"] == 1, drugbank_col]
        .dropna()
        .astype(str)
        .str.strip()
        .unique()
        .tolist()
    )

    print(f"  Found {len(active_ids)} DrugBank IDs with MajorityVote_Class == 1")
    return drugbank_col, active_ids


def map_drugbank_to_uniprot(
    targets_path: str, active_db_ids: List[str]
) -> List[Tuple[str, str, str]]:
    """
    From DrugBank_Targets.csv, map active DrugBank IDs to UniProt IDs.

    Returns list of tuples: (DrugBank_ID, UniProt_ID, Target_Name)
    """
    print(f"\nLoading targets from: {targets_path}")
    df = pd.read_csv(targets_path, sep=None, engine="python")

    # Column names in the sample:
    # 'UniProt ID' and 'Drug IDs' and 'Name'
    required_cols = ["UniProt ID", "Drug IDs"]
    for col in required_cols:
        if col not in df.columns:
            raise KeyError(
                f"Column '{col}' not found in targets file. "
                f"Available columns: {list(df.columns)}"
            )

    active_set: Set[str] = set(active_db_ids)
    mappings_set: Set[Tuple[str, str, str]] = set()

    for _, row in df.iterrows():
        drug_ids_field = str(row["Drug IDs"])
        if not drug_ids_field or drug_ids_field.lower() == "nan":
            continue

        # Split multi-valued field like "DB11300; DB11311; DB11571"
        raw_ids = [x.strip() for x in drug_ids_field.split(";") if x.strip()]
        uniprot = str(row["UniProt ID"]).strip()
        target_name = str(row.get("Name", "")).strip()

        for dbid in raw_ids:
            if dbid in active_set and uniprot:
                mappings_set.add((dbid, uniprot, target_name))

    mappings = sorted(mappings_set)
    print(f"  Found {len(mappings)} DrugBank–UniProt mappings for active compounds")
    return mappings


def write_mapping_csv(
    mappings: List[Tuple[str, str, str]], out_csv: str
) -> None:
    """
    Write DrugBank–UniProt mapping to CSV.
    """
    if not mappings:
        print("WARNING: no mappings to write to CSV.")
        return

    df_out = pd.DataFrame(
        mappings, columns=["DrugBank_ID", "UniProt_ID", "Target_Name"]
    )
    df_out.to_csv(out_csv, index=False)
    print(f"  Mapping table written to: {out_csv}")


def parse_fasta_and_select(
    fasta_path: str, wanted_uniprots: Set[str], out_fasta: str
) -> None:
    """
    Read DrugBank_protein.fasta and write only the sequences whose
    UniProt ID is in wanted_uniprots to out_fasta.

    FASTA headers look like:
      >drugbank_target|P45059 Peptidoglycan synthase FtsI (DB00303)

    We treat the token after 'drugbank_target|' and before the first
    space as the UniProt ID (e.g. 'P45059').
    """
    print(f"\nScanning FASTA file: {fasta_path}")
    print(f"  Number of UniProt IDs of interest: {len(wanted_uniprots)}")

    n_written = 0
    with open(fasta_path, "r") as fin, open(out_fasta, "w") as fout:
        write_seq = False

        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()

                # Default: don't write this record
                write_seq = False

                # Extract UniProt ID from header
                uniprot_id = None
                if "|" in header:
                    after_pipe = header.split("|", 1)[1]
                    # UniProt ID is up to first whitespace
                    uniprot_id = after_pipe.split()[0].strip()

                if uniprot_id and uniprot_id in wanted_uniprots:
                    write_seq = True
                    n_written += 1

                if write_seq:
                    fout.write(line)
            else:
                if write_seq:
                    fout.write(line)

    print(f"  Wrote {n_written} matching FASTA entries to: {out_fasta}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract UniProt targets and protein sequences for DrugBank compounds "
            "predicted active by MajorityVote_Class == 1."
        )
    )
    parser.add_argument(
        "--predictions",
        default="DrugBank_approved_predictions_all_models.csv",
        help="Predictions file with MajorityVote_Class column",
    )
    parser.add_argument(
        "--targets",
        default="DrugBank_Targets.csv",
        help="DrugBank targets metadata file (with UniProt ID and Drug IDs columns)",
    )
    parser.add_argument(
        "--fasta",
        default="DrugBank_protein.fasta",
        help="FASTA file with DrugBank protein sequences",
    )
    parser.add_argument(
        "--out-mapping",
        default="MajorityVote_DrugBank_to_UniProt.csv",
        help="Output CSV mapping DrugBank_ID to UniProt_ID and target name",
    )
    parser.add_argument(
        "--out-fasta",
        default="MajorityVote_targets.fasta",
        help="Output FASTA with sequences for targets of MajorityVote active compounds",
    )

    args = parser.parse_args()

    try:
        _, active_db_ids = load_active_drugbank_ids(args.predictions)
        mappings = map_drugbank_to_uniprot(args.targets, active_db_ids)
        write_mapping_csv(mappings, args.out_mapping)

        wanted_uniprots = {u for _, u, _ in mappings}
        parse_fasta_and_select(args.fasta, wanted_uniprots, args.out_fasta)

        print("\nDone.")
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
