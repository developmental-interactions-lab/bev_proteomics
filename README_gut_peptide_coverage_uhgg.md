# Gut peptide coverage in human gut bacterial communities

Version **3.1.2**, 27 September 2026. Script: `gut_peptide_coverage_uhgg.py`.

This workflow searches UHGP-100 representative proteins for candidate peptides,
infers genome carriage through cluster membership, and combines within-species
reference-genome breadth with UHGV bacterial abundance profiles. It is deliberately
orthology-independent: a qualifying sequence in any protein counts. It estimates
**encoded sequence carriage**, not protein expression, exposure or antibody binding.

## Quick start

Install Python >=3.9 and the required packages in your existing environment:

```bash
python -m pip install pandas numpy openpyxl pyahocorasick
```

Put the pipeline, SLURM script, your `candidate_peptides.fasta`, and
`controls.fasta` in the project directory. The supplied SLURM job expects:

```text
/scratch/project_2009813/evprots/check_conservation/
  gut_peptide_coverage_uhgg.py
  run_gut_peptide_coverage_uhgg.slurm
  candidate_peptides.fasta
  controls.fasta                 # required by supplied SLURM job
  venv/
  uhgv/relative_abundance.tsv
  uhgv/sample_metadata.tsv
  uhgv/host_genomes_info.tsv
  uhgp/genomes-all_metadata.tsv
  uhgp/uhgp-100.faa
  uhgp/uhgp-100.tsv
```

Keep the actual extracted protein paths if the archive creates a subdirectory;
edit the two paths in the SLURM script accordingly. Verify the account, partition
and environment module on your cluster. The supplied account is inferred from your
scratch path. The existing `python-data` module and virtual environment setup are
retained; they were not tested on CSC in this revision.

```bash
sbatch run_gut_peptide_coverage_uhgg.slurm
```

The job writes to **`results_uhgg_v3`** and requires the supplied blinded controls FASTA.
It creates `input_peptides.fasta` containing the actual combined input. Duplicate
FASTA IDs are rejected. The controls have neutral IDs and a shuffled order; the pipeline is not given
their positive/negative identities. It processes them identically to candidates.
The separate `controls_key.tsv` is for interpretation only and is never read by the
pipeline or the SLURM script.

Start with:

- `UHGG_peptide_coverage_summary.xlsx`: one workbook, with interpretive notes.
- `concise_summary.tsv`: equal-rate extrapolated estimates alongside reference-supported
  shares, evaluability and the extreme missing-reference scenarios, for exact,
  <=1 conservative and <=2 conservative matches.
- `evaluability_summary.tsv`: how much abundance and reference information was evaluated.

No trusted migration of legacy v2 caches is performed. The first v3 run performs
its own search and membership inventory. Adding control peptides requires a new
protein search in any case. Later compatible v3 stages can be reused safely.
When updating an existing v3.0.x result to v3.1.2, run `all --reuse` in the same
output directory. Preparation and species breadth are rebuilt using the new
genome-level taxonomy mapping; the validated, compatible UHGP protein scan and
peptide-independent membership inventory are reused. Save any earlier result
tables you want to compare before the new summaries overwrite them.
If v3.1.0 or v3.1.1 has already completed in that directory, v3.1.2 changes
only the summary: run `summarize` with the same combined FASTA and output directory.
The supplied SLURM job's `all --reuse` also works; compatible earlier stages
are reused.

## Panel selection and denominator

The input profiles are restricted to bulk metagenomes and prokaryotic records,
joined to GTDB r207 taxonomy, then restricted to species-labelled Bacteria.
Abundances are renormalized within that retained bacterial component per sample.
Samples with no positive retained bacterial abundance are excluded and counted.
The published UHGV sample metadata repeat some `sample_name` values. The pipeline
uses only their BioProject assignment: repeated rows with the same
`bioproject_id` are collapsed and counted in `preparation_diagnostics.json`.
Conflicting BioProject assignments for one name stop preparation with examples;
unrelated metadata fields do not affect the abundance analysis.

Thus **100% means species-labelled bacterial abundance represented in the input
reference profiles**. It excludes archaea, references without a species label and
organisms absent from the profiling reference framework. The workflow cannot
recover the last fraction from these tables. It reports the retained labelled
bacterial fraction of the *input prokaryotic abundance* separately; that diagnostic
is not the fraction of the entire original biological community recovered.

Species are ranked by mean relative abundance across all included samples, with
zeros included where absent. The smallest shared panel reaching
`--target-community-coverage 0.95` is retained. Ties are ordered by species name.
The target is met **on average across samples**, not necessarily in every sample.
Mean, median, 5th percentile, minimum and maximum panel coverage are reported.
There is no additional abundance, prevalence or two-BioProject selection gate.

Selection precedes reference-quality filtering: the panel is not silently
renormalized or expanded to hide quality losses. An additional preparation
diagnostic reports abundance belonging to any taxonomy-mappable species with HQ
genomes, including species outside the panel. This is a pre-membership ceiling on
HQ species availability, not a guarantee that all those genomes are in UHGP.

## Reference quality and mapping

UHGG metadata provide the genome denominator. Each UHGG genome receives its own
GTDB r207 species label from the matching `data_source=UHGG` record in the UHGV
host taxonomy. Versionless MGYG matching is attempted only when unambiguous.
For a genome with no direct label, a species label is inferred from its UHGG
cluster only when the modal species among overlapping genomes has at least 90%
purity. A mixed cluster never overrides a directly identified genome, even when
its modal fraction is below 90%. This matters when updated taxonomy splits an
older UHGG species cluster into multiple species. `--cluster-mapping-min-purity`
controls **only this fallback**. The selected-genome table records the mapping
source; species QC counts direct and cluster-inferred genomes separately. Cluster
overlap and modal-purity columns describe cluster composition, not the certainty
of each direct genome label; the 90% purity criterion does not gate direct
genome assignments. Cross-release taxonomy remains an approximation. The v3.1
change alters species membership and HQ breadth denominators, so v3.0.x and
v3.1.x peptide coverage values must not be compared as if only reporting changed.

Operational HQ criteria are completeness >=90% and contamination <=5%. These are
not a claim that all selected genomes satisfy the full MIMAG HQ definition.

Before peptide-specific mapping, the script inventories genome aliases across the
**complete UHGP membership table**. This inventory is peptide-independent and is
cached using the membership/metadata input identities. It identifies genomes with
at least one represented protein. It does not prove that every gene or candidate
locus was recovered. MGYG and legacy GUT_GENOME aliases are supported; unresolved
recognized aliases and lines without any recognizable alias are reported.

A genome is evaluable for the main breadth calculation when it passes HQ filtering
and occurs in that inventory. A species needs at least one such genome. Species
supported by only 1-4 evaluable HQ genomes are flagged; their abundance is reported.
No shrinkage estimator or uncertainty distribution is imposed.

## Coverage, evaluability and the equal-rate estimate

For species s, let:

- `a_js`: its bacterial relative abundance in sample j, before panel restriction;
- `N_s`: all HQ reference genomes assigned to that selected species;
- `n_s`: HQ genomes represented in the membership inventory;
- `k_s`: represented HQ genomes linked to a qualifying peptide hit.

The reported `breadth_hq_*` is **k_s/n_s**, missing when n_s=0. The
reference-supported community share is:

```text
C_j = sum over evaluable selected species [a_js * k_s/n_s]
```

This extrapolates observed reference breadth to the corresponding species'
abundance. `C_j` is the **reference-supported share** of the species-labelled
bacterial abundance, with no assigned support for unevaluable species. It is
not a guaranteed biological lower bound because genome carriage itself is
inferred from protein-cluster membership. For example, 20% species abundance
and 75% breadth contribute 15 percentage points.

Evaluability is reported separately:

```text
E_j = sum abundance of selected species with n_s > 0
R_j = sum over selected species [a_js * n_s/N_s]
```

Use a zero ratio when N_s=0. `E_j` describes species-level evaluability. In
practice, the **unevaluable fraction** `1-E_j` contains species outside the shared
panel and panel species without a mapped UHGG genome, without an HQ genome, or
without a represented HQ genome in UHGP membership. It is a reference-availability
gap, not evidence that the microbiota sample is poor or that the peptide is
absent. Species with only 1-4 evaluable HQ genomes remain evaluable but are
flagged separately as sparse evidence. `R_j` weights evaluability by the share
of all HQ reference genomes represented in UHGP; it is a **reference-based
diagnostic**, not measured missing biological abundance.

The main whole-community **equal-rate extrapolation** is calculated for each
sample and peptide tier as:

```text
T_j = C_j / E_j, if E_j > 0; otherwise undefined
```

It assumes that the **abundance-weighted match rate** in the unevaluable fraction
equals the rate in the evaluable fraction of that same sample. Equivalently,
`T_j = C_j + (1-E_j)*(C_j/E_j)`. Individual missing species need not have equal
rates. This is a transparent working assumption, not an observed carriage
measurement or a testable property of the unevaluable species. Species without
HQ references are not a random sample of taxa. The estimate can be too high
or too low, depending on the peptide and the taxa missed. Report the median
`T_j` **together with** the
median `C_j`, median `E_j`, median `1-E_j`, and the number of samples with
`E_j=0`, and show the extreme missing-reference scenarios alongside the point
estimate. Each median is computed from per-sample quantities; medians cannot be
multiplied or divided to reproduce another median. The main denominator remains
the species-labelled bacterial abundance, not total gut prokaryotes.

## Missing-reference sensitivity

The missing-reference sensitivity endpoints are:

```text
L_j = sum over selected species [a_js * k_s/N_s]
U_j = L_j + (1 - R_j)
```

These **extreme sensitivity scenarios** assign unresolved HQ genomes and
unevaluable/outside-panel species either zero or complete carriage. L_j, C_j
and U_j are calculated for each sample before medians are taken. Never add
separately calculated medians to reconstruct an endpoint. The per-sample file
also provides `unevaluated_species_upper = C_j + 1-E_j`, which changes only
completely unevaluated-species assumptions and holds the within-species breadth
estimate fixed.
For each sample, the columns satisfy
`missing_reference_lower <= breadth_weighted_coverage <= unevaluated_species_upper <= missing_reference_upper`.
When E_j>0, the equal-rate extrapolation also satisfies
`missing_reference_lower <= breadth_weighted_coverage <= equal_rate_extrapolated_coverage <= missing_reference_upper`.
The concise summary orders these four columns from lower to upper for each tier;
median endpoints are computed separately from the sample-level endpoints.

**These are assumption ranges, not confidence intervals or guaranteed biological
bounds.** They do not cover all sampling bias, missing loci within represented
genomes, taxonomy errors, abundance-profile errors or cluster-expansion errors.
When all HQ genomes are represented, L_j=C_j and the extreme range reduces to the
simple missing-species range.

Disjoint sample-level abundance categories are provided for: outside-panel species;
selected species with no mapped genomes; selected species with mapped genomes but
no HQ genomes; HQ species with no membership representation; and evaluable species.
Within the last category, partial reference loss is reported as
`partly_missing_hq_reference_weight`. Do not add that reference weight to a physical
community-abundance estimate without retaining the reference-based interpretation.

## Sequence tiers

| Tier | Definition |
|---|---|
| `exact` | No amino-acid substitutions |
| `cons_le1` | <=1 substitution, each with Grantham distance <=50 |
| `cons_le2` | <=2 substitutions, each with Grantham distance <=50 |
| `cons_le3` | <=3 substitutions, each with Grantham distance <=50 |
| `unrestricted_le1` | <=1 substitution, any chemistry |
| `unrestricted_le2` | <=2 substitutions, any chemistry |

Exact matches are included in all broader tiers. Search windows allow no indels.
The seed search uses four disjoint query segments to retrieve all windows with at
most three substitutions; full candidate-sized windows are then checked. Search
hits with three non-conservative substitutions can appear in the raw audit table
but contribute to none of the six reported tiers.

Candidate peptides default to 12-25 aa and must contain only the 20 standard amino
acids, with unique identifiers. Protein terminal `*` symbols are trimmed; internal
stops and ambiguous residues remain barriers, and windows containing them are
excluded. Candidate FASTAs reject stop symbols rather than silently altering them.

Grantham <=50 is an operational physicochemical cutoff, not an antibody-binding
model. Exact sequence matching applies to the **representative protein**;
per-genome peptide carriage is inferred through membership. Orthology is not required.

A secondary binary QC view counts a species in full only if breadth reaches
`--broad-breadth 0.75`, with at least `--min-hq-genomes-for-broad 1` represented HQ
genomes. This threshold does not affect the continuous main score. Actual settings
are retained in provenance and used in workbook explanations.

## Running stages and changing settings

A complete command matching the SLURM defaults is:

```bash
python gut_peptide_coverage_uhgg.py all \
  --peptides candidate_peptides.fasta \
  --outdir results_uhgg_v3 \
  --relative-abundance uhgv/relative_abundance.tsv \
  --sample-metadata uhgv/sample_metadata.tsv \
  --host-genomes-info uhgv/host_genomes_info.tsv \
  --uhgg-metadata uhgp/genomes-all_metadata.tsv \
  --uhgp-faa uhgp/uhgp-100.faa \
  --uhgp-membership uhgp/uhgp-100.tsv \
  --target-community-coverage 0.95 \
  --threads 16 --reuse
```

`prepare`, `search` and `summarize` can also run separately. Use each command's
`--help` for required arguments. Later stages validate the supplied candidate
FASTA against the prepared table; changing `--peptides` cannot silently select a
new query set without preparation.

| Change | Required work |
|---|---|
| Only re-export summaries or change Excel row budget | `summarize` |
| Upgrade a completed v3.1.0/v3.1.1 run to v3.1.2 reporting | `summarize` using the same combined FASTA and output directory; no new UHGP scan |
| Target abundance, input abundances, HQ criteria or taxonomy bridge | `all --reuse`; compatible protein search is retained, affected mapping/breadth are rebuilt |
| Add/remove/change peptides, including controls | `all --reuse`; incompatible scan and dependent stages are rebuilt |
| Grantham cutoff | `search --reuse` with the new cutoff, then `summarize`; scan is rebuilt |
| Broad-carrier threshold or minimum HQ count for that QC | `search --reuse`, then `summarize`; scan and compatible mapping are retained |
| UHGP FASTA, membership or UHGG metadata | `all --reuse`; source identities invalidate dependent stages |

On `summarize`, optional `--broad-breadth` and `--conservative-grantham-max` are only
consistency assertions. A mismatch with saved search settings stops the run.
For a SLURM run, use `results_uhgg_v3/input_peptides.fasta` when summarizing the
combined candidate/control analysis, not the candidate-only input FASTA.

```bash
python gut_peptide_coverage_uhgg.py summarize \
  --peptides results_uhgg_v3/input_peptides.fasta \
  --outdir results_uhgg_v3
```

`--reuse` validates separate completed protein-search and hit-mapping caches.
Incompatible stages are rebuilt, not silently reused. The peptide-independent
membership inventory is reused automatically when its inputs are unchanged.
Completion manifests and atomic final output writes prevent ordinary incomplete
runs from being treated as completed caches. Do not edit cached TSV/JSON files.
Use a separate output directory for simultaneous jobs: the workflow does not lock
output directories against concurrent writers.

Small source files have full SHA-256 checksums in manifests. Huge FASTA/membership
files use absolute path, size, nanosecond modification time and SHA-256 of sampled
blocks. Completed outputs are guarded by size/mtime. These guards are designed for
ordinary workflow changes, **not adversarial integrity verification**. Record full
external-data checksums once for publication; moving giant files invalidates reuse.

## Outputs

| File | Purpose |
|---|---|
| `concise_summary.tsv` | Three equal-rate extrapolated medians, companion reference-supported shares, evaluability and ordered extreme endpoints |
| `extended_summary.tsv` | Equal-rate extrapolations for all six tiers, reference-supported scores, sample 5th percentiles, extreme endpoints and equal-BioProject sensitivity of the supported score |
| `UHGG_peptide_coverage_summary.xlsx` | Single workbook with README, results and combined QC |
| `evaluability_summary.tsv` | Mean/median/p05/min/max across samples for abundance and reference diagnostics |
| `sample_evaluability.tsv.gz` | Per-sample panel, quality, mapping and partial-reference diagnostics |
| `species_evaluability.tsv` | Species abundance, original/evaluable genome counts and mapping evidence |
| `peptide_species_breadth.tsv.gz` | Hit counts and breadth for each candidate/species/tier |
| `peptide_sample_coverage.tsv.gz` | C, equal-rate extrapolation (undefined when E=0), L, U and secondary metrics per sample/candidate/tier |
| `bioproject_coverage.tsv` | Per-project medians and sample counts |
| `broad_hq_threshold_qc.tsv` | Secondary binary species-carrier interpretation |
| `uhgp_candidate_hits.tsv.gz` | Matching representative windows and substitution details |
| `uhgp_hit_rep_to_selected_genome.tsv.gz` | Hit representative-to-selected-genome mapping |
| `peptide_genome_tier_presence.tsv.gz` | Genome-level presence; header-only if no hits |
| `uhgp_hit_examples.tsv` | Limited hit examples, not the complete evidence table |
| `uhgp_genome_inventory.tsv.gz` | Peptide-independent membership representation of all metadata genomes |
| `selected_genome_membership_qc.tsv.gz` | Selected genomes with quality/representation flags |
| `all_species_abundance.tsv`, `selected_species.tsv` | Panel selection audit |
| `selected_species_abundance_matrix.tsv.gz` | Original-denominator abundance of selected species |
| `included_sample_metadata.tsv` | Sample metadata and project assignments retained for summaries |
| `*_qc.json`, `*_diagnostics.json`, `*_manifest.json` | Counts, settings, source/cache identities and limitations |

The former duplicate `UHGG_peptide_coverage.xlsx` is no longer written. A sheet
exceeding `--excel-max-rows` (default 200,000 data rows) is replaced by an explicit
pointer to the complete TSV; it is not silently truncated. TSVs remain complete.
`excel_export_qc.json` lists omitted sheets. Excel writing streams rows to reduce
memory use. The genome-hit merge can still be memory-intensive for ubiquitous
short peptides, although duplicate windows are collapsed before expansion.

BioProject summaries are descriptive sensitivity checks. Missing BioProject IDs
are excluded from that view and counted. Equal-project weighting does not establish
representativeness or remove repeated-subject dependence. No health-status or
adult-only cohort restriction is applied.

## Resources, data and validation

The supplied SLURM request retains 16 CPUs, 64 GB RAM and 24 hours. These are starting
allocations, not measured guarantees. The first v3 run adds a complete membership
inventory pass. Protein scan cost depends strongly on peptide length/composition;
12-aa queries and ubiquitous controls can increase seed hits and mapping size.
Inspect job logs and scheduler resource reports before reducing or increasing
allocations. Temporary hit files and final compressed tables coexist during writing.

Download the UHGV inputs from:

- https://portal.nersc.gov/UHGV/read_mapping/relative_abundance.tsv
- https://portal.nersc.gov/UHGV/read_mapping/sample_metadata.tsv
- https://portal.nersc.gov/UHGV/host_predictions/host_genomes_info.tsv

Use UHGG/UHGP v2.0.2 from:

- https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/genomes-all_metadata.tsv
- https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.2/protein_catalogue/uhgp-100.tar.gz

Extract the archive and use the actual paths to `uhgp-100.faa` and `uhgp-100.tsv`.
All hit representatives must occur in membership; incomplete representative mapping
stops the workflow rather than producing unqualified coverage results.
If representatives are missing, `membership_mapping_qc.json` records their count
and up to 20 example IDs before the workflow stops. Check the
`hit_representatives_seen_fraction` in a prior v2 QC file if available: a value
below 1 predicts this stop, but does not establish why a representative is absent.
The self-pair behavior of the specific distributed membership file has not been
verified here. A compatible completed protein scan can be reused after diagnosis.
Catalogue membership remains an inference of sequence carriage: 100%-identity clustering
with partial-length coverage does not establish complete per-member peptide identity.
The magnitude and net direction of resulting error are not quantified here.

Run the supplied synthetic regression suite beside the script:

```bash
python -m unittest -v test_gut_peptide_coverage_uhgg.py
```

Nineteen tests passed during development, including 500 brute-force search comparisons,
known end-to-end coverage/evaluability values, stop-symbol barriers, no-hit outputs,
cache invalidation and panel expansion, summary-setting rejection, single/multiple
process agreement, numeric sample-name preservation, missing-representative QC,
and Excel omission guards. The development environment used the
Python seed-search fallback because pyahocorasick was unavailable; the suite also
tests the Aho-Corasick path when that dependency is installed. The production-scale
UHGP run, CSC environment and real-catalogue mapping completeness have not been
validated in this environment. See `VALIDATION.md`.

Positive controls are biological sanity checks, not universal pass/fail standards
or correction factors. Expected-negative controls are also included in the same neutral-ID FASTA.
See `controls_notes.md` for source sequences and
the scope of conservation verification. A low control score warrants investigation
of sequence variants and mapping; it does not by itself diagnose a namespace bug.

## References and publication record

- Almeida A et al. A unified catalog of 204,938 reference genomes from the human
  gut microbiome. *Nature Biotechnology* (2021).
  https://doi.org/10.1038/s41587-020-0603-3
- UHGG/UHGP v2.0.2 release, MGnify.
  https://www.ebi.ac.uk/metagenomics/api/v1/genome-catalogues/human-gut-v2-0-2
- Camargo AP et al. A genomic atlas of the human gut virome elucidates genetic
  factors shaping host interactions. *bioRxiv* (2025).
  https://doi.org/10.1101/2025.11.01.686033
- Nayfach S, Camargo A. Unified Human Gut Virome (UHGV), version 1.0 [data set].
  Zenodo (2025). https://doi.org/10.5281/zenodo.17402089
- Grantham R. Amino acid difference formula to help explain protein evolution.
  *Science* 185:862-864 (1974). https://doi.org/10.1126/science.185.4154.862

Repository intended for archiving:
https://github.com/developmental-interactions-lab/bev_proteomics/

This deliverable was not pushed to that repository. Record its actual commit/tag
when published. Retain input release/checksums, combined peptide FASTA, v3 script,
manifests, QC JSON, TSVs and job logs. A licence should be chosen by the repository
owner; this revision does not assign one. See `CHANGELOG.md` for changes from v2.
