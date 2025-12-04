# README – From Majority-Vote DrugBank Predictions to Filtered BLAST Targets

This folder contains everything needed to go from a table of DrugBank majority-vote
activity predictions to a final set of DrugBank target proteins that:

1. Were predicted **active** by a majority of machine-learning models, and  
2. Have at least one **BLASTP hit with e-value ≤ 1e-4** against the
   *Ancylostoma caninum* proteome.

All steps described here start from the point where the script
`extract_majorityvote_targets.py` is available in this directory.

---

## 1. Directory layout

Inside the `BLAST` directory you should see:

```text
BLAST/
├── DrugBank_Targets.csv
├── DrugBank_approved_predictions_all_models.csv
├── DrugBank_protein.fasta
├── MajorityVote_DrugBank_to_UniProt.csv
├── MajorityVote_targets.fasta
├── MajorityVote_targets_filtered.fasta
├── MajorityVote_targets_filtered_mapping.csv
├── MajorityVote_targets_vs_ancylostoma.tsv
├── db/
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.pdb
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.phr
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.pin
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.pot
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.psq
│   ├── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.ptf
│   └── ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa.pto
├── extract_majorityvote_targets.py
├── filter_blast_and_map_drugbank.py
└── run_blast_and_filter.sh
````

**Key inputs**

* `DrugBank_approved_predictions_all_models.csv`
  Machine-learning predictions for DrugBank compounds. Includes a `DrugBank ID`
  column and a `MajorityVote_Class` column (1 = predicted active, 0 = inactive).

* `DrugBank_Targets.csv`
  DrugBank target annotation table. Must contain at least:

  * `UniProt ID` – target UniProt accession
  * `Drug IDs` – one or more DrugBank IDs (e.g. `DB00114; DB01017`).

* `DrugBank_protein.fasta`
  FASTA containing sequences for DrugBank targets. FASTA headers include the
  UniProt accession, for example:

  ```text
  >drugbank_target|P19113 Histidine decarboxylase (DB00114)
  ```

* `db/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa*`
  Pre-formatted BLASTP database for the *A. caninum* proteome
  (created previously with `makeblastdb`).

**Scripts**

* `extract_majorityvote_targets.py`
  Extracts targets of majority-vote active drugs and builds the query FASTA
  for BLAST.

* `run_blast_and_filter.sh`
  Convenience shell script that runs BLASTP and then calls
  `filter_blast_and_map_drugbank.py`.

* `filter_blast_and_map_drugbank.py`
  Filters BLAST results at a chosen e-value threshold and maps the passing
  sequences back to DrugBank IDs.

Generated intermediate and final files are described in the sections below.

---

## 2. Software requirements

* **Python** 3.7+
* **Python packages**

  * `pandas`
* **NCBI BLAST+** (here: `ncbi-blast-2.10.0+`), already unpacked in this folder.
  The script `run_blast_and_filter.sh` points directly to
  `./ncbi-blast-2.10.0+/bin/blastp`.

If `python` on your system points to Python 2, use `python3` instead in the
commands below.

---

## 3. Step 1 – Extract majority-vote DrugBank targets

This step uses `extract_majorityvote_targets.py` to:

1. Find all DrugBank compounds predicted active by the **majority** of ML models.
2. Map those DrugBank IDs to their UniProt targets using `DrugBank_Targets.csv`.
3. Pull out the corresponding protein sequences from `DrugBank_protein.fasta`.
4. Write out:

   * a mapping table (`MajorityVote_DrugBank_to_UniProt.csv`), and
   * a FASTA file of the unique target sequences
     (`MajorityVote_targets.fasta`).

### 3.1 Inputs

* `DrugBank_approved_predictions_all_models.csv`
* `DrugBank_Targets.csv`
* `DrugBank_protein.fasta`

### 3.2 Command

From inside the `BLAST` directory:

```bash
cd /path/to/BLAST
python extract_majorityvote_targets.py
```

### 3.3 What the script does (logic overview)

1. **Select majority-vote active drugs**

   * Reads `DrugBank_approved_predictions_all_models.csv`.
   * Keeps only rows where `MajorityVote_Class == 1`.
   * Collects the set of “active” DrugBank compound IDs.

2. **Map DrugBank IDs → UniProt IDs**

   * Reads `DrugBank_Targets.csv`.
   * For each row, splits the `Drug IDs` field on `;` to get individual
     DrugBank IDs.
   * If any of these IDs belong to the majority-vote active set, the
     corresponding `UniProt ID` is recorded.
   * The script writes this mapping to
     **`MajorityVote_DrugBank_to_UniProt.csv`**, which typically contains
     at least:

     * `DrugBank ID`
     * `UniProt ID`
     * (optionally additional columns from `DrugBank_Targets.csv`).

3. **Extract protein sequences**

   * Opens `DrugBank_protein.fasta`.
   * For each FASTA header, extracts the UniProt accession
     (e.g. in `drugbank_target|P19113 ...` it picks `P19113`).
   * If the UniProt ID matches one of the majority-vote targets from step 2,
     the sequence is written to **`MajorityVote_targets.fasta`**.
   * Only unique (UniProt) sequences are kept to avoid duplication.

### 3.4 Outputs

* **`MajorityVote_DrugBank_to_UniProt.csv`**
  Table mapping majority-vote active DrugBank compounds to their UniProt
  target IDs.

* **`MajorityVote_targets.fasta`**
  FASTA file containing the target protein sequences for all majority-vote
  active compounds. This file is used as the query set for BLASTP in the
  next step.

---

## 4. Step 2 – Run BLASTP against the *A. caninum* proteome

In this stage we compare each majority-vote target protein against the
*A. caninum* proteome using BLASTP and store the results in a tabular file.

This is automated by the shell script `run_blast_and_filter.sh`.

### 4.1 BLAST database

The BLAST database is already present in `db/`:

* Base name: `db/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa`

BLASTP uses this base name plus the associated index files (`.pin`, `.psq`,
etc.) that have been created previously with `makeblastdb`.

### 4.2 Command

From inside the `BLAST` directory:

```bash
cd /path/to/BLAST

# Make sure the script is executable (only needed once)
chmod +x run_blast_and_filter.sh

# Run BLASTP and subsequent filtering/mapping
./run_blast_and_filter.sh
```

The script assumes the BLAST+ binaries are in
`./ncbi-blast-2.10.0+/bin/`. If your BLAST+ installation is elsewhere,
edit the `BLASTP=` line near the top of `run_blast_and_filter.sh`.

### 4.3 What `run_blast_and_filter.sh` does

1. **BLASTP search**

   It runs:

   ```bash
   ./ncbi-blast-2.10.0+/bin/blastp \
     -query MajorityVote_targets.fasta \
     -db db/ancylostoma_caninum.PRJNA72585.WBPS19.protein.fa \
     -evalue 1e-2 \
     -max_target_seqs 10 \
     -num_threads 4 \
     -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
     -out MajorityVote_targets_vs_ancylostoma.tsv
   ```

   * **Query**: `MajorityVote_targets.fasta`
   * **Database**: *A. caninum* proteome (`db/ancylostoma_caninum...`)
   * **E-value threshold (for reporting)**: `1e-2`
   * **Output format**: tabular (BLAST outfmt 6)

   This produces **`MajorityVote_targets_vs_ancylostoma.tsv`** with one
   line per hit, containing (in order):
   `qseqid, sseqid, pident, length, mismatch, gapopen,
    qstart, qend, sstart, send, evalue, bitscore`.

2. **Call the Python filter & mapping script**

   After BLASTP finishes, the script calls:

   ```bash
   python3 filter_blast_and_map_drugbank.py \
     --blast MajorityVote_targets_vs_ancylostoma.tsv \
     --fasta MajorityVote_targets.fasta \
     --targets_csv DrugBank_Targets.csv \
     --evalue_cutoff 1e-4 \
     --out_fasta MajorityVote_targets_filtered.fasta \
     --out_csv MajorityVote_targets_filtered_mapping.csv
   ```

   This filters BLAST hits at **e ≤ 1e-4** and maps passing queries back
   to DrugBank IDs (details in the next section).

### 4.4 Output

* **`MajorityVote_targets_vs_ancylostoma.tsv`**
  Raw BLASTP results in tabular format (all hits with e ≤ 1e-2 as per
  the BLAST run parameters).

This file is kept as an intermediate, and the final filtered set is
produced by `filter_blast_and_map_drugbank.py`.

---

## 5. Step 3 – Filter BLAST hits and map back to DrugBank IDs

The script `filter_blast_and_map_drugbank.py` takes the BLAST results,
applies a stringent e-value filter, and links each passing query sequence
back to its UniProt ID and DrugBank compound(s).

### 5.1 Inputs

* `MajorityVote_targets_vs_ancylostoma.tsv` – BLAST outfmt 6 results.
* `MajorityVote_targets.fasta` – query FASTA used for BLASTP.
* `DrugBank_Targets.csv` – to map UniProt IDs back to DrugBank IDs.

### 5.2 Command (run inside `BLAST`)

Normally you do not need to call this directly because
`run_blast_and_filter.sh` already runs it. If you want to run it by
hand or change parameters:

```bash
python filter_blast_and_map_drugbank.py \
  --blast MajorityVote_targets_vs_ancylostoma.tsv \
  --fasta MajorityVote_targets.fasta \
  --targets_csv DrugBank_Targets.csv \
  --evalue_cutoff 1e-4 \
  --out_fasta MajorityVote_targets_filtered.fasta \
  --out_csv MajorityVote_targets_filtered_mapping.csv
```

To use a different e-value threshold, modify the `--evalue_cutoff` value.

### 5.3 Logic overview

1. **Identify passing queries (per-query best e-value)**

   * Reads `MajorityVote_targets_vs_ancylostoma.tsv` line by line.
   * For each query (`qseqid`), it tracks the **best (lowest) e-value**
     observed across all hits.
   * A query is considered to **pass** if its best e-value is
     `≤ evalue_cutoff` (default `1e-4`).

2. **Map query IDs → FASTA header, sequence, UniProt ID**

   * Parses `MajorityVote_targets.fasta`.
   * Uses the **first token of each FASTA header** (without the leading `>`)
     as the query ID (`qseqid`) to match the BLAST file.
   * Stores:

     * full header (without `>`)
     * concatenated sequence
   * Extracts the UniProt accession from the header. Typical patterns:

     * `drugbank_target|P45059 ...` → `P45059`
     * `P19113|something` → `P19113`
     * `Q9UI32` → `Q9UI32`

3. **Build UniProt ID → DrugBank ID mapping**

   * Reads `DrugBank_Targets.csv` using pandas.
   * Detects the UniProt column (preferring `UniProt ID` but also
     tolerating similar variants such as `Uniprot ID`).
   * Detects the DrugBank column (e.g. `Drug IDs`).
   * Splits `Drug IDs` on `;` to obtain individual DrugBank IDs and
     builds a mapping `UniProt ID → {DrugBank IDs}`.

4. **Write filtered FASTA and mapping table**

   For each query that passes the e-value cutoff:

   * Writes the full sequence to **`MajorityVote_targets_filtered.fasta`**
     using the original FASTA header.
   * Writes one row to **`MajorityVote_targets_filtered_mapping.csv`**
     with the columns:

     * `QueryID` – FASTA query ID (first header token)
     * `BestEvalue` – best BLAST e-value (formatted)
     * `FASTA_Header` – full header from `MajorityVote_targets.fasta`
     * `UniProt_ID` – UniProt accession parsed from header
     * `DrugBank_IDs` – semicolon-separated DrugBank IDs associated with
       that UniProt ID (as found in `DrugBank_Targets.csv`)

### 5.4 Outputs

* **`MajorityVote_targets_filtered.fasta`**
  FASTA file containing only those majority-vote target sequences that
  have at least one BLASTP hit at **e ≤ 1e-4** against the
  *A. caninum* proteome.

* **`MajorityVote_targets_filtered_mapping.csv`**
  A concise table directly linking each filtered target to its UniProt ID,
  DrugBank ID(s), and best BLAST e-value.

These two files represent the final filtered set of targets and their
drug associations.

---

## 6. Quick “one-liner” workflow summary

From the `BLAST` directory:

```bash
# 1) Extract majority-vote target sequences
python extract_majorityvote_targets.py

# 2) Run BLASTP against A. caninum and filter + map results
chmod +x run_blast_and_filter.sh   # only needed once
./run_blast_and_filter.sh
```

After these steps you will have:

* `MajorityVote_targets_filtered.fasta` – protein sequences of interest,
* `MajorityVote_targets_filtered_mapping.csv` – mapping of those sequences
  to UniProt IDs, DrugBank IDs, and BLAST statistics.

These are ready for downstream analyses (e.g. pathway enrichment, manual
inspection, prioritization of targets, etc.).

