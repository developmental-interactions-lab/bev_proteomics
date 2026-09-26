# Gut peptide coverage across the human gut microbiome

Reproducible workflow for `gut_peptide_coverage_uhgg_v2_1_1.py`

## What this script does

This script estimates how broadly a short bacterial peptide sequence is represented across abundant human-gut bacterial communities.

It combines:

1. **Human gut bacterial relative-abundance profiles** from the Unified Human Gut Virome (UHGV) read-mapping resource.
2. **Genome metadata** from the Unified Human Gastrointestinal Genome collection (UHGG v2.0.2).
3. **Protein sequences and cluster membership** from the Unified Human Gastrointestinal Protein catalogue (UHGP v2.0.2).

The analysis is intentionally **orthology-independent**. A candidate peptide is considered potentially relevant wherever the same or a sufficiently similar sequence occurs in any protein encoded by a bacterial genome. This is appropriate when the biological question is the potential taxonomic breadth of a linear antibody epitope rather than conservation of the peptide's original parent protein.

The headline quantity is **breadth-weighted community coverage**. For each gut metagenome and candidate peptide,

```text
coverage(sample) = sum over species [
    species relative abundance in that sample
    ×
    fraction of high-quality UHGG genomes of that species carrying the sequence
]
```

For example, if a species represents 20% of a sample and the qualifying peptide sequence occurs in 75% of its high-quality UHGG reference genomes, that species contributes 15% to the estimated community coverage.

The script reports the median of this quantity across human-gut bulk metagenomes.

---

## Sequence-similarity tiers

The script evaluates six peptide-matching tiers:

| Tier | Definition |
|---|---|
| `exact` | Exact peptide match |
| `cons_le1` | ≤1 substitution, all substitutions conservative |
| `cons_le2` | ≤2 substitutions, all substitutions conservative |
| `cons_le3` | ≤3 substitutions, all substitutions conservative |
| `unrestricted_le1` | ≤1 substitution of any type |
| `unrestricted_le2` | ≤2 substitutions of any type |

By default, a substitution is considered **conservative when its Grantham distance is ≤50**.

No insertions or deletions are allowed. Candidate-sized protein windows containing non-standard or ambiguous amino-acid symbols are skipped.

The main summary emphasizes:

- exact matches;
- ≤1 conservative substitution;
- ≤2 conservative substitutions.

The remaining tiers are sensitivity analyses.

---

# 1. Requirements

## Software

Recommended:

- Linux or another Unix-like environment
- Python ≥3.9
- `pandas`
- `numpy`
- `openpyxl`
- `pyahocorasick`

Install into a virtual environment, for example:

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install pandas numpy openpyxl pyahocorasick
```

`pyahocorasick` is strongly recommended and is effectively required for searching the full UHGP-100 catalogue efficiently.

The script itself uses only Python packages plus standard Unix utilities. GNU `grep` is used opportunistically to accelerate filtering of the large UHGP membership table; the script has a Python fallback.

## Hardware

The abundance-processing and summary stages are modest. The **UHGP-100 protein scan is the expensive step** and is best run on an HPC system or a workstation with substantial memory and storage.

A reasonable HPC starting point is:

```text
16 CPUs
64 GB RAM
up to 24 h wall time
```

Actual requirements depend on filesystem speed, CPU performance, and the number of candidate peptides.

The UHGP-100 archive and extracted files are large, so allow substantial local/scratch storage.

---

# 2. Download the external datasets

Create a working directory:

```bash
mkdir -p gut_peptide_coverage/{uhgv,uhgg,results}
cd gut_peptide_coverage
```

The commands below use `wget`. `curl -O` can be used instead.

## 2.1 UHGV abundance and metadata files

The script uses three files from the UHGV resource:

```text
UHGV/read_mapping/relative_abundance.tsv
UHGV/read_mapping/sample_metadata.tsv
UHGV/host_predictions/host_genomes_info.tsv
```

Download them:

```bash
cd uhgv

wget https://portal.nersc.gov/UHGV/read_mapping/relative_abundance.tsv
wget https://portal.nersc.gov/UHGV/read_mapping/sample_metadata.tsv
wget https://portal.nersc.gov/UHGV/host_predictions/host_genomes_info.tsv

cd ..
```

The UHGV download page is:

```text
https://uhgv.jgi.doe.gov/downloads
```

Why these files are needed:

- `relative_abundance.tsv` contains read-mapping-derived relative abundances of viruses and prokaryotes across bulk metagenomes and viromes.
- `sample_metadata.tsv` provides sample and BioProject information.
- `host_genomes_info.tsv` provides updated GTDB r207 taxonomy for prokaryotic reference genomes, including UHGG genomes.

The analysis uses **bulk metagenomes** and the **prokaryotic component** of the UHGV read-mapping data. It does not analyze viral sequences.

## 2.2 UHGG v2.0.2 genome metadata

Download the authoritative UHGG v2.0.2 genome metadata from MGnify:

```bash
cd uhgg

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv

cd ..
```

This file provides the UHGG genome denominator, species-cluster membership, completeness, contamination, and identifiers needed to link UHGP proteins to genomes.

## 2.3 UHGP-100 v2.0.2 protein catalogue

Download the 100%-identity clustered UHGP protein catalogue:

```bash
cd uhgg

wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-100.tar.gz
```

Inspect the archive if desired:

```bash
tar -tzf uhgp-100.tar.gz | head
```

Extract it:

```bash
tar -xzf uhgp-100.tar.gz
```

The analysis requires:

```text
uhgp-100.faa
uhgp-100.tsv
```

If the archive extracts into a subdirectory, use the actual paths to those two files in the command line.

Return to the project root:

```bash
cd ..
```

## 2.4 Optional: record file checksums

For reproducibility, it is useful to record checksums for the downloaded inputs:

```bash
sha256sum \
  uhgv/relative_abundance.tsv \
  uhgv/sample_metadata.tsv \
  uhgv/host_genomes_info.tsv \
  uhgg/genomes-all_metadata.tsv \
  uhgg/uhgp-100.faa \
  uhgg/uhgp-100.tsv \
  > input_sha256.txt
```

Adjust paths if the UHGP archive extracted into a subdirectory.

---

# 3. Prepare the candidate peptide FASTA

Candidate peptides are supplied as ordinary FASTA.

Example:

```fasta
>candidate_1
PEERREVIEWISPAIN
>candidate_2
PLEASECITEMYPAPER
>candidate_3
SIGNIFICANCEISFAKE
```

Requirements:

- FASTA identifiers must be unique.
- Only the 20 standard amino-acid letters are accepted.
- Default candidate length range: **12–25 aa**.

The default limits can be changed with:

```text
--min-peptide-length
--max-peptide-length
```

Save the file, for example, as:

```text
candidate_peptides.fasta
```

---

# 4. Recommended directory layout

A convenient layout is:

```text
gut_peptide_coverage/
├── gut_peptide_coverage_uhgg_v2_1_1.py
├── candidate_peptides.fasta
├── uhgv/
│   ├── relative_abundance.tsv
│   ├── sample_metadata.tsv
│   └── host_genomes_info.tsv
├── uhgg/
│   ├── genomes-all_metadata.tsv
│   ├── uhgp-100.faa
│   └── uhgp-100.tsv
└── results/
```

The exact layout is not required; all input paths are supplied explicitly.

---

# 5. Run the complete analysis

Activate the Python environment:

```bash
source venv/bin/activate
```

Then run:

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py all \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --relative-abundance uhgv/relative_abundance.tsv \
  --sample-metadata uhgv/sample_metadata.tsv \
  --host-genomes-info uhgv/host_genomes_info.tsv \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv \
  --uhgp-faa uhgg/uhgp-100.faa \
  --uhgp-membership uhgg/uhgp-100.tsv \
  --threads 16
```

This uses the publication-oriented defaults described below.

---

# 6. Important default parameters

The main defaults in v2.1.1 are:

```text
Candidate length                     12–25 aa
Species abundance threshold          ≥1% bacterial relative abundance
Abundance prevalence threshold       ≥1% of bulk samples
Minimum independent BioProjects      2
High-quality genome completeness     ≥90%
High-quality genome contamination    ≤5%
UHGG cluster→GTDB mapping purity      ≥90%
Conservative substitution            Grantham distance ≤50
Maximum substitutions searched       3
Secondary broad-carrier threshold    ≥75% of HQ genomes
Minimum HQ genomes for broad QC      1
```

The **75% broad-carrier threshold is not used to calculate the headline coverage metric**. It is retained as a secondary all-or-none QC/sensitivity view.

The principal coverage metric uses the actual high-quality-genome breadth continuously.

---

# 7. What happens internally

The `all` command runs three stages:

```text
prepare → search → summarize
```

They can also be run separately.

## 7.1 `prepare`

The preparation stage:

1. reads the candidate peptide FASTA;
2. reads UHGV prokaryotic abundance profiles;
3. keeps bulk-metagenome samples;
4. aggregates prokaryotic abundance to GTDB r207 species;
5. renormalizes within the evaluated species-labelled bacterial component;
6. identifies ecologically relevant species;
7. maps UHGG v2.0.2 species clusters to updated GTDB r207 species;
8. selects authoritative UHGG genomes belonging to the retained species;
9. marks high-quality genomes using completeness and contamination.

By default, a species is retained if it reaches **≥1% bacterial relative abundance in ≥1% of included bulk samples and does so in at least two BioProjects**.

UHGG species-cluster assignments are linked to updated GTDB r207 species using overlap with the UHGG genomes represented in `host_genomes_info.tsv`. A cluster mapping is accepted only when mapping purity is ≥90% by default.

## 7.2 `search`

The search stage scans the UHGP-100 representative protein FASTA for candidate-sized sequence windows.

Matches are classified by:

- total substitutions;
- number of conservative substitutions;
- number of non-conservative substitutions;
- Grantham distance of each substitution.

No gaps or indels are allowed.

Matching UHGP-100 representative proteins are expanded to their member genomes using `uhgp-100.tsv`.

For each candidate × species combination, the script calculates the fraction of:

- all mapped UHGG genomes carrying the qualifying sequence;
- high-quality UHGG genomes carrying the qualifying sequence.

This is the **reference-genome breadth**.

## 7.3 `summarize`

For each sample, peptide, and similarity tier, the script calculates:

```text
breadth-weighted coverage
=
Σ [species abundance × high-quality reference-genome breadth]
```

Species with no evaluable high-quality-genome breadth contribute zero.

The principal output is the **median breadth-weighted coverage across bulk gut metagenomes**.

A secondary QC calculation also treats species as all-or-none carriers when the qualifying sequence is found in at least 75% of their high-quality reference genomes.

---

# 8. Running the stages separately

This can be useful for debugging or for avoiding unnecessary repetition.

## Preparation only

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py prepare \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --relative-abundance uhgv/relative_abundance.tsv \
  --sample-metadata uhgv/sample_metadata.tsv \
  --host-genomes-info uhgv/host_genomes_info.tsv \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv
```

## Protein search and genome breadth

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py search \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv \
  --uhgp-faa uhgg/uhgp-100.faa \
  --uhgp-membership uhgg/uhgp-100.tsv \
  --threads 16
```

## Summarization only

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py summarize \
  --peptides candidate_peptides.fasta \
  --outdir results
```

The later stages expect the earlier-stage files to exist in the same output directory.

---

# 9. Resuming a run and `--reuse`

The expensive part is scanning UHGP-100 and mapping hit representatives back to selected UHGG genomes.

The `search` and `all` commands support:

```text
--reuse
```

Use it only when all of the following are unchanged:

- candidate peptide FASTA;
- UHGP-100 FASTA;
- UHGP membership table;
- species/genome selection;
- sequence-matching parameters relevant to the search.

Example:

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py search \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv \
  --uhgp-faa uhgg/uhgp-100.faa \
  --uhgp-membership uhgg/uhgp-100.tsv \
  --threads 16 \
  --reuse
```

Do **not** reuse old UHGP hit files after changing the candidate peptide set. The existing hit table cannot contain hits for newly added candidates.

---

# 10. Main output files

## `concise_summary.tsv`

The principal result table.

One row per candidate with:

```text
median_coverage_exact
median_coverage_le1_conservative
median_coverage_le2_conservative
```

These are median breadth-weighted community coverages across the analyzed bulk gut metagenomes.

For example:

```text
median_coverage_exact = 0.62
```

means that the estimated median fraction of the evaluated bacterial community carrying the exact peptide is 62%, using high-quality UHGG reference-genome breadth as a proxy for within-species prevalence.

It is **not** a direct measurement that 62% of bacterial cells carry the peptide.

## `extended_summary.tsv`

Adds:

```text
≤3 conservative substitutions
≤1 unrestricted substitution
≤2 unrestricted substitutions
```

These are sensitivity views and can be useful when considering possible broader antibody recognition.

## `UHGG_peptide_coverage_summary.xlsx`

Publication-oriented workbook containing the headline summary plus the supporting evidence and QC tables required for interpretation.

## `UHGG_peptide_coverage.xlsx`

Full audit workbook. In v2.1.1 its contents are intentionally very similar to the summary workbook so either can be archived for reproducibility.

---

# 11. Supporting and audit outputs

Important additional files include:

### Species and genome selection

```text
selected_species.tsv
selected_species_abundance_matrix.tsv.gz
selected_uhgg_genomes.tsv.gz
```

These document which species and genomes entered the analysis.

### Raw peptide/protein hits

```text
uhgp_candidate_hits.tsv.gz
```

Candidate-sized matching windows in UHGP-100 representatives, including substitution counts and Grantham-distance information.

### UHGP representative → UHGG genome mapping

```text
uhgp_hit_rep_to_selected_genome.tsv.gz
```

Maps matching UHGP representatives to selected UHGG genomes through the UHGP-100 membership table.

### Genome-level peptide presence

```text
peptide_genome_tier_presence.tsv.gz
```

Records whether each selected genome contains at least one qualifying match for each sequence-similarity tier.

### Species-level reference-genome breadth

```text
peptide_species_breadth.tsv.gz
```

For every candidate × selected species, this contains the number and fraction of all and high-quality UHGG genomes carrying each qualifying variant tier.

### Per-sample coverage

```text
peptide_sample_coverage.tsv.gz
```

Contains both:

```text
breadth_weighted_coverage
broad_hq_thresholded_coverage
```

for every candidate, sample, and sequence tier.

The first is the headline continuous metric. The second is the secondary ≥75%-breadth all-or-none QC metric.

### Broad-carrier QC

```text
broad_hq_threshold_qc.tsv
```

Summarizes the older binary carrier interpretation retained as a sensitivity check.

### Human-readable hit examples

```text
uhgp_hit_examples.tsv
```

Provides representative sequence matches for manual sanity checking.

### Diagnostic files

The workflow also writes JSON diagnostic/QC files documenting parameters, counts, software versions, mapping performance, and known limitations.

These should be retained with the results when the analysis is archived.

---

# 12. How to interpret the Excel workbook

The most important sheets are:

## `Concise summary`

Start here.

The three main columns report median breadth-weighted gut-community coverage for:

1. exact peptide matches;
2. ≤1 conservative substitution;
3. ≤2 conservative substitutions.

## `Candidate evaluation`

Provides a compact comparison between candidates and shows how much apparent coverage is gained when conservative or unrestricted substitutions are permitted.

A large jump from exact to ≤1 or ≤2 substitution coverage means that much of the predicted breadth depends on peptide variants rather than the exact immunizing sequence.

## `Species breadth`

Use this to identify which bacterial species drive a candidate's coverage.

For each species, inspect:

- mean abundance;
- number of high-quality UHGG genomes;
- `breadth_hq_*`.

A breadth of 1.0 means all evaluated high-quality UHGG reference genomes of that species are associated with a qualifying peptide sequence.

## `Per-sample coverage`

Shows the underlying sample-by-sample result rather than only the median.

This is useful for detecting candidates that have high average coverage but perform poorly in a subset of microbiomes.

## `≥75% breadth QC`

This is a secondary sensitivity analysis.

It deliberately converts within-species breadth into an all-or-none call: a species contributes its full abundance only when at least 75% of its high-quality reference genomes carry the qualifying variant.

Do not confuse this with the headline breadth-weighted metric.

## `UHGP hit examples`

Use this as a sequence-level sanity check.

Inspect:

- query sequence;
- matched protein window;
- substitutions;
- Grantham distances;
- total number of substitutions.

## `Run & QC`

Contains the analysis thresholds, software versions, diagnostics, and methodological caveats. This sheet should be archived with any published analysis.

---

# 13. Important limitations

## UHGG reference-genome breadth is not measured strain prevalence

UHGG genomes are not a random population sample.

A result such as:

```text
breadth_hq_exact = 0.80
```

means that the exact sequence was associated with 80% of the available high-quality UHGG reference genomes used for that species.

It does **not** establish that 80% of naturally occurring strains, bacterial cells, or hosts carry that peptide.

The breadth-weighted community coverage should therefore be interpreted as an **ecologically weighted reference-catalogue estimate**.

## UHGP-100 cluster expansion can slightly overestimate genome carriage

UHGP-100 is clustered at 100% amino-acid identity with an 80% target-coverage criterion. The distributed membership table does not contain the full amino-acid sequence of every cluster member.

A shorter cluster member can therefore, in principle, be credited with a terminal peptide that is present in the representative sequence but outside the shorter member.

Exact verification would require the individual protein sequences for every UHGG genome, which this workflow does not reconstruct.

The script records this explicitly in the QC output.

## Similar sequence does not guarantee antibody recognition

The conservative and unrestricted substitution tiers are sequence-similarity sensitivity analyses.

They are **not probabilities of antibody cross-reactivity**.

Recognition of a substituted peptide depends on:

- the position of the substitution;
- amino-acid chemistry;
- epitope conformation and accessibility;
- the actual polyclonal antibody repertoire.

For antibody-antigen design, exact matches should therefore be interpreted most directly, while conservative-substitution tiers provide plausible broader sequence-space estimates.

## The analysis is orthology-independent

A qualifying peptide may occur in a protein unrelated to the protein from which the candidate was originally selected.

This is intentional when estimating potential antibody recognition across bacteria.

If the research question instead concerns evolutionary conservation of the original protein or functional orthologues, this workflow is not sufficient by itself.

## The analyzed sample set is not explicitly a "healthy adult" cohort

The workflow uses eligible UHGV **human-gut bulk metagenomes** and does not impose a universal health-status or adult-only filter.

Do not describe the resulting coverage as specifically representing healthy adults unless such a filter is added separately and documented.

---

# 14. Changing the analysis thresholds

Run:

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py all --help
```

or the help for an individual stage:

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py prepare --help
python gut_peptide_coverage_uhgg_v2_1_1.py search --help
python gut_peptide_coverage_uhgg_v2_1_1.py summarize --help
```

Important options include:

```text
--abundant-relative
--abundant-prevalence
--min-abundant-bioprojects
--min-completeness
--max-contamination
--cluster-mapping-min-purity
--conservative-grantham-max
--broad-breadth
--min-hq-genomes-for-broad
--min-peptide-length
--max-peptide-length
```

When reporting modified results, record all non-default parameters.

---

# 15. Example SLURM job

The following is a reasonable starting template for an HPC cluster. Scheduler account and partition names are site-specific.

```bash
#!/bin/bash -l
#SBATCH --job-name=gut_peptide_cov
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=gut_peptide_cov_%j.out
#SBATCH --error=gut_peptide_cov_%j.err

set -euo pipefail

source venv/bin/activate

python gut_peptide_coverage_uhgg_v2_1_1.py all \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --relative-abundance uhgv/relative_abundance.tsv \
  --sample-metadata uhgv/sample_metadata.tsv \
  --host-genomes-info uhgv/host_genomes_info.tsv \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv \
  --uhgp-faa uhgg/uhgp-100.faa \
  --uhgp-membership uhgg/uhgp-100.tsv \
  --threads "${SLURM_CPUS_PER_TASK}"
```

For CSC Roihu, add the appropriate project account and partition directives and load the site's Python environment before activating the virtual environment.

One configuration used during development was based on:

```bash
module load python-data
python3 -m venv --system-site-packages venv
source venv/bin/activate
python -m pip install pyahocorasick openpyxl
```

Cluster modules change over time, so consult the current CSC documentation if these commands no longer match the environment.

---

# 16. Reproducibility checklist

For a publication or archived analysis, retain:

- `gut_peptide_coverage_uhgg_v2_1_1.py`
- the exact candidate FASTA
- downloaded input filenames and preferably SHA-256 checksums
- the script version
- all command-line parameters
- `run_diagnostics.json`
- preparation/search QC JSON files
- `concise_summary.tsv`
- `extended_summary.tsv`
- `peptide_species_breadth.tsv.gz`
- `peptide_sample_coverage.tsv.gz`
- the Excel workbook
- the SLURM/job log if run on HPC

Do not rerun with `--reuse` after changing the peptide FASTA or sequence-search parameters.

---

# 17. Data sources

## UHGV

Unified Human Gut Virome resource and downloads:

- https://uhgv.jgi.doe.gov/
- https://uhgv.jgi.doe.gov/downloads
- https://portal.nersc.gov/UHGV/

Files used here:

- `read_mapping/relative_abundance.tsv`
- `read_mapping/sample_metadata.tsv`
- `host_predictions/host_genomes_info.tsv`

## UHGG / UHGP

MGnify Unified Human Gastrointestinal Genome catalogue v2.0.2:

- https://www.ebi.ac.uk/metagenomics/api/v1/genome-catalogues/human-gut-v2-0-2
- https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/

Files used here:

- `genomes-all_metadata.tsv`
- `protein_catalogue/uhgp-100.tar.gz`, providing `uhgp-100.faa` and `uhgp-100.tsv`

Original UHGG publication:

Almeida A. et al. *A unified catalog of 204,938 reference genomes from the human gut microbiome.* Nature Biotechnology (2021). https://doi.org/10.1038/s41587-020-0603-3

---

# 18. Citation and versioning

If this workflow is used in a publication, cite the UHGV resource used for the abundance profiles, the UHGG/UHGP resource, and the version of this script.

Recommended repository practice:

```text
script: gut_peptide_coverage_uhgg_v2_1_1.py
version: 2.1.1
UHGG/UHGP: v2.0.2
GTDB taxonomy bridge: r207 via UHGV host_genomes_info.tsv
```

If changes are made to the script, increment the version and document changes in a `CHANGELOG.md`.

---

# 19. Minimal quick-start

After installing the Python packages and downloading the required data:

```bash
python gut_peptide_coverage_uhgg_v2_1_1.py all \
  --peptides candidate_peptides.fasta \
  --outdir results \
  --relative-abundance uhgv/relative_abundance.tsv \
  --sample-metadata uhgv/sample_metadata.tsv \
  --host-genomes-info uhgv/host_genomes_info.tsv \
  --uhgg-metadata uhgg/genomes-all_metadata.tsv \
  --uhgp-faa uhgg/uhgp-100.faa \
  --uhgp-membership uhgg/uhgp-100.tsv \
  --threads 16
```

Then begin interpretation with:

```text
results/concise_summary.tsv
results/UHGG_peptide_coverage_summary.xlsx
```

The most important number is the **median breadth-weighted community coverage for the exact peptide**, followed by the ≤1 and ≤2 conservative-substitution sensitivity tiers.
