#!/usr/bin/env python3
"""
gut_peptide_coverage_uhgg_v2_1_0.py

Publication-oriented analysis of candidate peptide sequence coverage across abundant
human-gut bacterial species using UHGV abundance profiles and the UHGG/UHGP v2.0.2
reference catalogue.

The analysis is deliberately orthology-independent: a peptide variant can occur in
any protein encoded by a bacterial genome. This is appropriate for estimating the
potential taxonomic breadth of an antibody epitope rather than conservation of the
parent protein.

Written by ChatGPT 5.6 Sol, evaluated by Mikael Niku, 2026-09-26

HEADLINE DEFINITION
-------------------
For each species, reference-genome breadth is the fraction of its high-quality UHGG
reference genomes associated with a qualifying peptide variant (default HQ: >=90%
completeness and <=5% contamination). Community coverage is calculated continuously
for each UHGV bulk metagenome as the sum, across evaluated species, of species
relative abundance multiplied by that species' high-quality reference-genome breadth.
Thus a species at 20% abundance with 75% reference-genome breadth contributes 15%
to the breadth-weighted community coverage. A >=75% breadth carrier classification
is retained only as a secondary QC/sensitivity view.

Peptide variants are grouped into:
  exact                  : 0 substitutions
  cons_le1               : <=1 substitution, all conservative
  cons_le2               : <=2 substitutions, all conservative
  cons_le3               : <=3 substitutions, all conservative
  unrestricted_le1       : <=1 substitution of any physicochemical type
  unrestricted_le2       : <=2 substitutions of any physicochemical type

A substitution is called conservative when its Grantham distance is <=50 by
default. No gaps or indels are allowed.

IMPORTANT LIMITATION
--------------------
UHGP-100 was clustered at 100% amino-acid identity with an 80% target-coverage
criterion. Its distributed membership file does not contain the full sequence of
every cluster member. Consequently, expanding a matching representative protein
to all genomes in its UHGP-100 cluster can, in principle, credit a truncated member
that does not contain a terminal epitope present in the representative. Exact
per-genome verification would require the individual genome protein sequences (or
reconstruction from the UHGG GFF files) and is not performed here. The script
therefore uses only high-quality genomes for the headline breadth estimate and records
this limitation explicitly in the QC output. The continuous breadth-weighted coverage
reduces the all-or-none effect of a hard within-species cutoff, but reference-genome
breadth must still not be interpreted as measured population strain prevalence.

OUTPUT
------
  concise_summary.tsv
      Median breadth-weighted community coverage for exact, <=1 conservative, and
      <=2 conservative variants.
  extended_summary.tsv
      Adds <=3 conservative and <=1/<=2 unrestricted substitution views.
  UHGG_peptide_coverage_summary.xlsx
      Standalone Excel workbook containing only the concise and extended summaries.
  UHGG_peptide_coverage.xlsx
      Full workbook containing the two summaries plus species-level breadth, selected
      species, per-sample coverage, >=75% breadth QC, representative hit examples,
      and run/QC information.
  Additional TSV/TSV.GZ and JSON audit files.

DEPENDENCIES
------------
Python >=3.9, pandas, numpy, openpyxl; pyahocorasick is strongly recommended and
required for large UHGP-100 searches.

Version 2.1.0
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import multiprocessing as mp
import os
import platform
import re
import shutil
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

try:
    import ahocorasick
except ImportError:
    ahocorasick = None

VERSION = "2.1.1"
AA20 = set("ACDEFGHIKLMNPQRSTVWY")
MGYG_RE = re.compile(r"(MGYG\d{9}(?:\.\d+)?)")
GUT_RE = re.compile(r"(GUT_GENOME\d{6,9})")

# Grantham distances (Grantham R. Science. 1974;185:862-864).
# Symmetric matrix; diagonal values are added programmatically as 0.
_GRANTHAM_ROWS = {
    "A": {"R":112,"N":111,"D":126,"C":195,"Q":91,"E":107,"G":60,"H":86,"I":94,"L":96,"K":106,"M":84,"F":113,"P":27,"S":99,"T":58,"W":148,"Y":112,"V":64},
    "R": {"N":86,"D":96,"C":180,"Q":43,"E":54,"G":125,"H":29,"I":97,"L":102,"K":26,"M":91,"F":97,"P":103,"S":110,"T":71,"W":101,"Y":77,"V":96},
    "N": {"D":23,"C":139,"Q":46,"E":42,"G":80,"H":68,"I":149,"L":153,"K":94,"M":142,"F":158,"P":91,"S":46,"T":65,"W":174,"Y":143,"V":133},
    "D": {"C":154,"Q":61,"E":45,"G":94,"H":81,"I":168,"L":172,"K":101,"M":160,"F":177,"P":108,"S":65,"T":85,"W":181,"Y":160,"V":152},
    "C": {"Q":154,"E":170,"G":159,"H":174,"I":198,"L":198,"K":202,"M":196,"F":205,"P":169,"S":112,"T":149,"W":215,"Y":194,"V":192},
    "Q": {"E":29,"G":87,"H":24,"I":109,"L":113,"K":53,"M":101,"F":116,"P":76,"S":68,"T":42,"W":130,"Y":99,"V":96},
    "E": {"G":98,"H":40,"I":134,"L":138,"K":56,"M":126,"F":140,"P":93,"S":80,"T":65,"W":152,"Y":122,"V":121},
    "G": {"H":98,"I":135,"L":138,"K":127,"M":127,"F":153,"P":42,"S":56,"T":59,"W":184,"Y":147,"V":109},
    "H": {"I":94,"L":99,"K":32,"M":87,"F":100,"P":77,"S":89,"T":47,"W":115,"Y":83,"V":84},
    "I": {"L":5,"K":102,"M":10,"F":21,"P":95,"S":142,"T":89,"W":61,"Y":33,"V":29},
    "L": {"K":107,"M":15,"F":22,"P":98,"S":145,"T":92,"W":61,"Y":36,"V":32},
    "K": {"M":95,"F":102,"P":103,"S":121,"T":78,"W":110,"Y":85,"V":97},
    "M": {"F":28,"P":87,"S":135,"T":81,"W":67,"Y":36,"V":21},
    "F": {"P":114,"S":155,"T":103,"W":40,"Y":22,"V":50},
    "P": {"S":74,"T":38,"W":147,"Y":110,"V":68},
    "S": {"T":58,"W":177,"Y":144,"V":124},
    "T": {"W":128,"Y":92,"V":69},
    "W": {"Y":37,"V":88},
    "Y": {"V":55},
}
GRANTHAM: Dict[Tuple[str,str], int] = {}
for a in AA20:
    GRANTHAM[(a,a)] = 0
for a, row in _GRANTHAM_ROWS.items():
    for b, d in row.items():
        GRANTHAM[(a,b)] = d
        GRANTHAM[(b,a)] = d

TIERS = [
    ("exact", "Exact match"),
    ("cons_le1", "≤1 conservative substitution"),
    ("cons_le2", "≤2 conservative substitutions"),
    ("cons_le3", "≤3 conservative substitutions"),
    ("unrestricted_le1", "≤1 substitution (non-conservative allowed)"),
    ("unrestricted_le2", "≤2 substitutions (non-conservative allowed)"),
]


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def die(msg: str) -> None:
    raise SystemExit(f"ERROR: {msg}")


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def clean_col(x: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(x).strip().lower()).strip()


def autodetect_col(columns: Sequence[object], candidates: Sequence[str], regexes: Sequence[str] = ()) -> Optional[str]:
    cols = [str(c) for c in columns]
    norm = {clean_col(c): c for c in cols}
    for c in candidates:
        if clean_col(c) in norm:
            return norm[clean_col(c)]
    for rgx in regexes:
        p = re.compile(rgx, re.I)
        hits = [c for c in cols if p.search(c)]
        if len(hits) == 1:
            return hits[0]
    return None


def open_text(path: str | Path, mode: str = "rt"):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, mode, encoding=None if "b" in mode else "utf-8", errors=None if "b" in mode else "replace")
    return open(p, mode, encoding=None if "b" in mode else "utf-8", errors=None if "b" in mode else "replace")


def parse_fasta(path: str | Path) -> Iterator[Tuple[str, str]]:
    name = None
    seq: List[str] = []
    with open_text(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper().replace("*", "")
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(re.sub(r"\s+", "", line))
        if name is not None:
            yield name, "".join(seq).upper().replace("*", "")


def read_candidate_peptides(path: str | Path, min_len: int, max_len: int) -> pd.DataFrame:
    rows = []
    seen = set()
    for pid, seq in parse_fasta(path):
        if pid in seen:
            die(f"Duplicate FASTA ID: {pid}")
        seen.add(pid)
        bad = sorted(set(seq) - AA20)
        if bad:
            die(f"Candidate {pid} contains non-standard residues: {','.join(bad)}")
        if not min_len <= len(seq) <= max_len:
            die(f"Candidate {pid} has length {len(seq)}; allowed range is {min_len}-{max_len} aa.")
        rows.append({"peptide_id": pid, "sequence": seq, "length": len(seq)})
    if not rows:
        die("No peptide sequences found in candidate FASTA.")
    return pd.DataFrame(rows)


def read_host_genomes(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    use = ["genome_id", "data_source", "completeness", "contamination", "lineage_gtdb_r207_updated"]
    if p.suffix.lower() == ".zip":
        with zipfile.ZipFile(p) as z:
            names = [n for n in z.namelist() if n.endswith(".tsv") and not n.startswith("__MACOSX")]
            if len(names) != 1:
                die(f"Expected one TSV in {p}; found {names}")
            with z.open(names[0]) as fh:
                df = pd.read_csv(fh, sep="\t", usecols=lambda c: c in use, low_memory=False)
    else:
        df = pd.read_csv(p, sep="\t", usecols=lambda c: c in use, low_memory=False)
    needed = {"genome_id", "data_source", "lineage_gtdb_r207_updated"}
    if not needed.issubset(df.columns):
        die(f"host_genomes_info is missing columns: {sorted(needed - set(df.columns))}")
    return df


def parse_gtdb_lineage(lineage: object) -> Dict[str, str]:
    ranks = {r: "" for r in "dpcofgs"}
    for part in str(lineage).split(";"):
        m = re.match(r"^([dpcofgs])__(.*)$", part.strip())
        if m:
            ranks[m.group(1)] = m.group(2).strip()
    return ranks


def add_taxonomy(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    tx = out["lineage_gtdb_r207_updated"].map(parse_gtdb_lineage)
    names = dict(zip("dpcofgs", ["domain","phylum","class","order","family","genus","species"]))
    for r, name in names.items():
        out[name] = [x[r] for x in tx]
    return out


def base_mgyg(x: object) -> str:
    return re.sub(r"\.\d+$", "", str(x).strip())


def read_uhgg_metadata(path: str | Path) -> Tuple[pd.DataFrame, Dict[str, str]]:
    header = pd.read_csv(path, sep="\t", nrows=0)
    cols = header.columns
    genome_col = autodetect_col(cols, ["Genome", "genome", "MGYG accession"], [r"^genome$"])
    species_rep_col = autodetect_col(cols, ["Species_rep", "Species rep", "species representative"], [r"species.*rep"])
    complete_col = autodetect_col(cols, ["Completeness", "completeness"], [r"^completeness$"])
    contam_col = autodetect_col(cols, ["Contamination", "contamination"], [r"^contamination$"])
    original_col = autodetect_col(
        cols,
        ["Original_name", "Original name", "Original_ID", "Original ID", "Original_genome", "Original genome"],
        [r"original.*(name|id|genome|accession)"],
    )
    if genome_col is None or species_rep_col is None or complete_col is None or contam_col is None:
        die(
            "Could not identify required columns in genomes-all_metadata.tsv. "
            f"Detected Genome={genome_col}, Species_rep={species_rep_col}, Completeness={complete_col}, Contamination={contam_col}."
        )
    usecols = [genome_col, species_rep_col, complete_col, contam_col] + ([original_col] if original_col else [])
    m = pd.read_csv(path, sep="\t", usecols=usecols, dtype=str, low_memory=False)
    out = pd.DataFrame({
        "genome_id": m[genome_col].astype(str).str.strip(),
        "species_rep": m[species_rep_col].astype(str).str.strip(),
        "completeness": pd.to_numeric(m[complete_col], errors="coerce"),
        "contamination": pd.to_numeric(m[contam_col], errors="coerce"),
    })
    if original_col:
        out["original_genome_id"] = m[original_col].astype(str).str.strip()
    else:
        out["original_genome_id"] = ""
    out = out.drop_duplicates("genome_id", keep="first")
    detected = {
        "genome": genome_col,
        "species_rep": species_rep_col,
        "completeness": complete_col,
        "contamination": contam_col,
        "original_genome_id": original_col or "",
    }
    return out, detected


def map_uhgg_clusters_to_r207(uhgg: pd.DataFrame, host_bac: pd.DataFrame, min_purity: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    h = host_bac[host_bac["data_source"].astype(str).eq("UHGG")][["genome_id", "species"]].drop_duplicates("genome_id")
    h["genome_base"] = h["genome_id"].map(base_mgyg)
    u = uhgg[["genome_id", "species_rep"]].copy()
    u["genome_base"] = u["genome_id"].map(base_mgyg)

    exact = u.merge(h[["genome_id","species"]], on="genome_id", how="left")
    missing = exact["species"].isna()
    if missing.any():
        # Version-insensitive fallback is safe only when the host base ID is unique.
        hbase = h.groupby("genome_base").filter(lambda x: x["species"].nunique() == 1).drop_duplicates("genome_base")
        fallback = exact.loc[missing, ["genome_id","species_rep","genome_base"]].merge(
            hbase[["genome_base","species"]], on="genome_base", how="left"
        )
        exact.loc[missing, "species"] = fallback.set_index("genome_id").reindex(exact.loc[missing,"genome_id"])["species"].values

    x = exact[exact["species"].fillna("").ne("")].copy()
    counts = x.groupby(["species_rep","species"]).size().rename("n_overlap").reset_index()
    totals = counts.groupby("species_rep")["n_overlap"].sum().rename("n_cluster_overlap")
    counts = counts.merge(totals, on="species_rep", how="left")
    counts["mapping_purity"] = counts["n_overlap"] / counts["n_cluster_overlap"]
    best = counts.sort_values(["species_rep","n_overlap"], ascending=[True,False]).drop_duplicates("species_rep")
    best["mapping_accepted"] = best["mapping_purity"].ge(min_purity)

    out = uhgg.merge(best[["species_rep","species","mapping_purity","n_cluster_overlap","mapping_accepted"]], on="species_rep", how="left")
    out = out.rename(columns={"species":"species_key"})
    # The left merge introduces missing values for UHGG clusters with no r207 overlap,
    # which can coerce this column to object dtype. Normalize it explicitly to bool
    # so Boolean negation cannot become integer bitwise negation (-1/-2).
    out["mapping_accepted"] = out["mapping_accepted"].eq(True).astype(bool)
    return out, best


def build_genome_alias_map(uhgg: pd.DataFrame) -> Dict[str, str]:
    aliases: Dict[str, str] = {}
    ambiguous: Set[str] = set()
    def add(alias: str, current: str):
        alias = str(alias).strip()
        if not alias or alias.lower() in {"nan","none"}:
            return
        if alias in aliases and aliases[alias] != current:
            ambiguous.add(alias)
        else:
            aliases[alias] = current
    for r in uhgg.itertuples(index=False):
        add(r.genome_id, r.genome_id)
        add(base_mgyg(r.genome_id), r.genome_id)
        orig = getattr(r, "original_genome_id", "")
        add(orig, r.genome_id)
        # Original accessions may appear as a prefix in protein IDs.
        if orig:
            add(base_mgyg(orig), r.genome_id)
    for a in ambiguous:
        aliases.pop(a, None)
    return aliases


def prepare_species_and_genomes(args) -> None:
    outdir = ensure_dir(args.outdir)
    peptides = read_candidate_peptides(args.peptides, args.min_peptide_length, args.max_peptide_length)
    peptides.to_csv(outdir / "candidate_peptides.tsv", sep="\t", index=False)

    log("Reading UHGV host-genome taxonomy (updated GTDB r207)...")
    host = add_taxonomy(read_host_genomes(args.host_genomes_info))
    host_bac = host[host["domain"].eq("Bacteria") & host["species"].ne("")].copy()
    host_map = host_bac[["genome_id","species"]].drop_duplicates("genome_id")

    log("Reading UHGV bulk-metagenome relative abundances...")
    ra = pd.read_csv(args.relative_abundance, sep="\t", low_memory=False)
    need = {"sample_name","sample_type","genome_type","genome","relative_abundance"}
    if not need.issubset(ra.columns):
        die(f"relative_abundance.tsv missing columns: {sorted(need-set(ra.columns))}")
    ra = ra[ra["sample_type"].astype(str).eq("bulk") & ra["genome_type"].astype(str).eq("prok")].copy()
    ra["relative_abundance"] = pd.to_numeric(ra["relative_abundance"], errors="coerce").fillna(0.0)
    ra = ra.merge(host_map, left_on="genome", right_on="genome_id", how="inner")
    # The host map above contains species-labelled Bacteria only. Archaea and
    # prokaryotic reference genomes without a species label are therefore excluded.
    ra = ra.groupby(["sample_name","species"], as_index=False)["relative_abundance"].sum()
    denom = ra.groupby("sample_name")["relative_abundance"].transform("sum")
    ra = ra[denom.gt(0)].copy()
    ra["abundance"] = ra["relative_abundance"] / denom[denom.gt(0)]

    sm = pd.read_csv(args.sample_metadata, sep="\t", low_memory=False)
    if "sample_name" not in sm.columns or "bioproject_id" not in sm.columns:
        die("sample_metadata.tsv must contain sample_name and bioproject_id.")
    sm = sm.drop_duplicates("sample_name")
    samples = pd.DataFrame({"sample_name": sorted(ra["sample_name"].unique())}).merge(sm, on="sample_name", how="left")
    sample_ids = set(samples["sample_name"].astype(str))
    ra = ra[ra["sample_name"].astype(str).isin(sample_ids)].copy()
    nsamples = ra["sample_name"].nunique()
    if nsamples == 0:
        die("No eligible bulk bacterial metagenomes remain.")
    sample_to_bp = samples.set_index("sample_name")["bioproject_id"].astype(str).to_dict()
    ra["bioproject_id"] = ra["sample_name"].map(sample_to_bp).fillna("NA")

    g = ra.groupby("species")
    stats = g.agg(
        abundance_sum=("abundance","sum"),
        prevalence_n=("sample_name","nunique"),
        max_abundance=("abundance","max"),
        n_bioprojects_detected=("bioproject_id","nunique"),
    )
    stats["mean_abundance"] = stats["abundance_sum"] / nsamples
    stats["prevalence_any"] = stats["prevalence_n"] / nsamples
    high = ra[ra["abundance"].ge(args.abundant_relative)]
    hs = high.groupby("species").agg(
        n_samples_ge_abundant=("sample_name","nunique"),
        n_bioprojects_ge_abundant=("bioproject_id","nunique"),
    )
    stats = stats.join(hs, how="left").fillna({"n_samples_ge_abundant":0,"n_bioprojects_ge_abundant":0})
    stats["prevalence_ge_abundant"] = stats["n_samples_ge_abundant"] / nsamples
    stats["selected"] = stats["prevalence_ge_abundant"].ge(args.abundant_prevalence) & stats["n_bioprojects_ge_abundant"].ge(args.min_abundant_bioprojects)
    selected = stats[stats["selected"]].copy().sort_values("mean_abundance", ascending=False)
    if selected.empty:
        die("Species-selection thresholds selected zero species.")

    log("Reading authoritative UHGG genomes-all_metadata.tsv denominator...")
    uhgg_raw, detected_cols = read_uhgg_metadata(args.uhgg_metadata)
    uhgg_mapped, cluster_map = map_uhgg_clusters_to_r207(uhgg_raw, host_bac, args.cluster_mapping_min_purity)
    bad_clusters = cluster_map[~cluster_map["mapping_accepted"]]
    if len(bad_clusters):
        log(f"WARNING: {len(bad_clusters):,} UHGG species clusters had r207 mapping purity below {args.cluster_mapping_min_purity:.0%}; they are excluded.")
    uhgg_mapped.loc[~uhgg_mapped["mapping_accepted"], "species_key"] = np.nan

    selected_genomes = uhgg_mapped[uhgg_mapped["species_key"].isin(selected.index)].copy()
    selected_genomes["high_quality"] = selected_genomes["completeness"].ge(args.min_completeness) & selected_genomes["contamination"].le(args.max_contamination)

    gmstats = selected_genomes.groupby("species_key").agg(
        n_uhgg_genomes_all=("genome_id","nunique"),
        n_uhgg_genomes_hq=("high_quality","sum"),
        mean_completeness=("completeness","mean"),
        species_cluster=("species_rep","first"),
        cluster_mapping_purity=("mapping_purity","first"),
    )
    selected = selected.join(gmstats, how="left")
    selected[["n_uhgg_genomes_all","n_uhgg_genomes_hq"]] = selected[["n_uhgg_genomes_all","n_uhgg_genomes_hq"]].fillna(0).astype(int)
    selected["mapped_to_authoritative_uhgg"] = selected["n_uhgg_genomes_all"].gt(0)
    selected["hq_evaluable"] = selected["n_uhgg_genomes_hq"].gt(0)
    selected = selected.reset_index().rename(columns={"species":"species_key"})

    # Abundance matrix is intentionally restricted to the pre-specified common-species panel.
    panel_species = selected.loc[selected["mapped_to_authoritative_uhgg"], "species_key"].tolist()
    mat = ra[ra["species"].isin(panel_species)].pivot_table(index="sample_name", columns="species", values="abundance", aggfunc="sum", fill_value=0.0)
    mat = mat.reindex(sorted(ra["sample_name"].unique()), fill_value=0.0)
    panel_cov = mat.sum(axis=1)

    selected.to_csv(outdir / "selected_species.tsv", sep="\t", index=False)
    selected_genomes.to_csv(outdir / "selected_uhgg_genomes.tsv.gz", sep="\t", index=False, compression="gzip")
    mat.to_csv(outdir / "selected_species_abundance_matrix.tsv.gz", sep="\t", compression="gzip")
    cluster_map.to_csv(outdir / "uhgg_cluster_to_gtdb_r207_mapping.tsv.gz", sep="\t", index=False, compression="gzip")

    diag = {
        "script_version": VERSION,
        "candidate_peptides_n": int(len(peptides)),
        "candidate_length_min": int(peptides["length"].min()),
        "candidate_length_max": int(peptides["length"].max()),
        "bulk_metagenome_samples_n": int(nsamples),
        "bioprojects_n": int(samples["bioproject_id"].astype(str).nunique()),
        "species_selection": {
            "abundant_relative_threshold": args.abundant_relative,
            "minimum_sample_prevalence_at_threshold": args.abundant_prevalence,
            "minimum_bioprojects_at_threshold": args.min_abundant_bioprojects,
            "selected_species_n": int(len(selected)),
            "selected_species_mapped_to_authoritative_uhgg_n": int(selected["mapped_to_authoritative_uhgg"].sum()),
        },
        "panel_sample_weighted_mean_community_coverage": float(panel_cov.mean()),
        "panel_median_community_coverage": float(panel_cov.median()),
        "uhgg_total_genomes_in_metadata": int(uhgg_raw["genome_id"].nunique()),
        "selected_uhgg_genomes_all_n": int(selected_genomes["genome_id"].nunique()),
        "selected_uhgg_genomes_hq_n": int(selected_genomes.loc[selected_genomes["high_quality"], "genome_id"].nunique()),
        "high_quality_definition": {"min_completeness": args.min_completeness, "max_contamination": args.max_contamination},
        "uhgg_metadata_detected_columns": detected_cols,
        "taxonomy_mapping": "UHGG species clusters mapped to UHGV lineage_gtdb_r207_updated by overlapping UHGG genomes; mapping propagated to all genomes in each accepted UHGG species cluster.",
        "cluster_mapping_min_purity": args.cluster_mapping_min_purity,
        "archaea_excluded": True,
        "notes": [
            "UHGV abundances are renormalized over species-labelled bacterial host genomes after archaea and unlabelled prokaryotic genomes are excluded.",
            "The two-BioProject rule requires the species to reach the abundance threshold in at least two BioProjects.",
            "Community coverage cannot exceed the abundance represented by the selected common-species panel; the observed panel mean and median ceilings are reported above.",
        ],
    }
    (outdir / "preparation_diagnostics.json").write_text(json.dumps(diag, indent=2), encoding="utf-8")
    log(f"Prepared {len(selected):,} selected species from {nsamples:,} bulk metagenomes; panel mean coverage={panel_cov.mean():.1%}.")


def grantham_distance(a: str, b: str) -> int:
    """Return Grantham distance for two standard amino acids.

    Ambiguous/non-standard residues must be filtered before this function is
    called. Raise a normal exception rather than SystemExit so that an
    unexpected residue in a multiprocessing worker is propagated to the
    parent process instead of potentially terminating a worker silently.
    """
    try:
        return GRANTHAM[(a,b)]
    except KeyError as exc:
        raise ValueError(f"No Grantham distance for amino-acid pair {a}/{b}") from exc


def classify_variant(query: str, target: str, conservative_cutoff: int) -> Dict[str, object]:
    subs = []
    for i, (a,b) in enumerate(zip(query,target), start=1):
        if a != b:
            d = grantham_distance(a,b)
            subs.append((i,a,b,d,d <= conservative_cutoff))
    n = len(subs)
    ncons = sum(x[4] for x in subs)
    nnon = n - ncons
    return {
        "n_substitutions": n,
        "n_conservative": ncons,
        "n_nonconservative": nnon,
        "substitutions": ";".join(f"{a}{i}{b}:{d}" for i,a,b,d,_ in subs),
        "max_grantham": max([x[3] for x in subs], default=0),
        "grantham_sum": sum(x[3] for x in subs),
        "exact": n == 0,
        "cons_le1": n <= 1 and nnon == 0,
        "cons_le2": n <= 2 and nnon == 0,
        "cons_le3": n <= 3 and nnon == 0,
        "unrestricted_le1": n <= 1,
        "unrestricted_le2": n <= 2,
    }


def partition_query(seq: str, n_parts: int) -> List[Tuple[int,str]]:
    n_parts = min(max(1,n_parts), len(seq))
    base, rem = divmod(len(seq), n_parts)
    out=[]; pos=0
    for i in range(n_parts):
        n = base + (1 if i < rem else 0)
        out.append((pos, seq[pos:pos+n])); pos += n
    return out


def _seed_payloads(peptides: pd.DataFrame, max_substitutions: int):
    payloads: Dict[str, List[Tuple[str,int,int,str]]] = defaultdict(list)
    for r in peptides.itertuples(index=False):
        for off, seed in partition_query(r.sequence, max_substitutions + 1):
            payloads[seed].append((r.peptide_id, off, r.length, r.sequence))
    return payloads


def build_seed_automaton(peptides: pd.DataFrame, max_substitutions: int):
    payloads = _seed_payloads(peptides, max_substitutions)
    if ahocorasick is None:
        return None, payloads
    A = ahocorasick.Automaton()
    for seed, pl in payloads.items():
        A.add_word(seed, (seed,pl))
    A.make_automaton()
    return A, payloads


def search_one_sequence(A, payloads, prot_id: str, prot_seq: str, max_substitutions: int, conservative_cutoff: int) -> Tuple[List[Tuple], int]:
    proposed=set(); rows=[]; n_ambiguous_windows=0
    if A is not None:
        iterator = ((end-len(seed)+1, seed, pl) for end,(seed,pl) in A.iter(prot_seq))
    else:
        def fallback():
            for seed, pl in payloads.items():
                start=0
                while True:
                    i=prot_seq.find(seed,start)
                    if i < 0: break
                    yield i,seed,pl
                    start=i+1
        iterator=fallback()
    for seed_start, seed, pl in iterator:
        for pid,qoff,qlen,qseq in pl:
            start=seed_start-qoff
            key=(pid,start)
            if key in proposed: continue
            proposed.add(key)
            if start < 0 or start+qlen > len(prot_seq): continue
            win=prot_seq[start:start+qlen]
            # Grantham distances are defined only for the 20 standard amino
            # acids. A window containing X/B/Z/J/U/O/* (or any other
            # non-standard symbol) has unknown substitution chemistry and is
            # therefore excluded rather than classified arbitrarily.
            if any(aa not in AA20 for aa in win):
                n_ambiguous_windows += 1
                continue
            mm=sum(a!=b for a,b in zip(qseq,win))
            if mm > max_substitutions: continue
            c=classify_variant(qseq,win,conservative_cutoff)
            rows.append((
                pid,qseq,prot_id,start+1,win,c["n_substitutions"],c["n_conservative"],c["n_nonconservative"],
                c["substitutions"],c["max_grantham"],c["grantham_sum"],
                int(c["exact"]),int(c["cons_le1"]),int(c["cons_le2"]),int(c["cons_le3"]),
                int(c["unrestricted_le1"]),int(c["unrestricted_le2"]),
            ))
    return rows, n_ambiguous_windows


_WORKER_A=None; _WORKER_PAYLOADS=None; _WORKER_MAX=None; _WORKER_CUTOFF=None

def _worker_init(peptide_records, max_substitutions, conservative_cutoff):
    global _WORKER_A,_WORKER_PAYLOADS,_WORKER_MAX,_WORKER_CUTOFF
    pdf=pd.DataFrame(peptide_records)
    _WORKER_A,_WORKER_PAYLOADS=build_seed_automaton(pdf,max_substitutions)
    _WORKER_MAX=max_substitutions; _WORKER_CUTOFF=conservative_cutoff


def _worker_batch(batch):
    rows=[]; n_ambiguous_windows=0
    try:
        for pid,seq in batch:
            hit_rows, n_ambig = search_one_sequence(
                _WORKER_A,_WORKER_PAYLOADS,pid,seq,_WORKER_MAX,_WORKER_CUTOFF
            )
            rows.extend(hit_rows)
            n_ambiguous_windows += n_ambig
    except BaseException as exc:
        # Convert worker-level SystemExit and other unexpected failures into a
        # normal exception that multiprocessing propagates to the parent.
        raise RuntimeError(
            f"UHGP worker failed while processing batch near protein {pid!r}: {exc}"
        ) from exc
    return len(batch),rows,n_ambiguous_windows


def fasta_batches(path: str | Path, batch_size: int):
    b=[]
    for x in parse_fasta(path):
        b.append(x)
        if len(b)>=batch_size:
            yield b; b=[]
    if b: yield b


def scan_uhgp(peptides: pd.DataFrame, fasta: str | Path, out_path: Path, threads: int, batch_size: int, progress_every: int, conservative_cutoff: int) -> Set[str]:
    max_substitutions=3
    if ahocorasick is None and Path(fasta).stat().st_size > 100_000_000:
        die("pyahocorasick is required for a full UHGP run. Install with: python -m pip install pyahocorasick")
    header=["peptide_id","query_sequence","protein_rep_id","protein_start_1based","matched_sequence","n_substitutions","n_conservative","n_nonconservative","substitutions","max_grantham","grantham_sum","exact","cons_le1","cons_le2","cons_le3","unrestricted_le1","unrestricted_le2"]
    matched=set(); nseq=0; nhit=0; n_ambiguous_windows=0; nextlog=progress_every
    with open(out_path,"wt",encoding="utf-8",newline="") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(header)
        if threads<=1:
            A,payloads=build_seed_automaton(peptides,max_substitutions)
            for batch in fasta_batches(fasta,batch_size):
                for pid,seq in batch:
                    rows, n_ambig = search_one_sequence(A,payloads,pid,seq,max_substitutions,conservative_cutoff)
                    n_ambiguous_windows += n_ambig
                    for r in rows: w.writerow(r); matched.add(r[2])
                    nhit += len(rows)
                nseq += len(batch)
                if nseq>=nextlog:
                    log(f"UHGP search: {nseq:,} representative proteins; {nhit:,} qualifying windows; {n_ambiguous_windows:,} ambiguous windows skipped")
                    nextlog += progress_every
        else:
            recs=peptides[["peptide_id","sequence","length"]].to_dict("records")
            ctx=mp.get_context("fork" if sys.platform.startswith("linux") else "spawn")
            with ctx.Pool(threads,initializer=_worker_init,initargs=(recs,max_substitutions,conservative_cutoff)) as pool:
                for nb,rows,n_ambig in pool.imap_unordered(_worker_batch,fasta_batches(fasta,batch_size),chunksize=1):
                    nseq += nb
                    n_ambiguous_windows += n_ambig
                    for r in rows: w.writerow(r); matched.add(r[2])
                    nhit += len(rows)
                    if nseq>=nextlog:
                        log(f"UHGP search: {nseq:,} representative proteins; {nhit:,} qualifying windows; {n_ambiguous_windows:,} ambiguous windows skipped")
                        nextlog += progress_every
    log(
        f"UHGP scan complete: {nseq:,} representatives; {nhit:,} qualifying windows "
        f"in {len(matched):,} representatives; {n_ambiguous_windows:,} ambiguous windows skipped."
    )
    search_qc = {
        "script_version": VERSION,
        "representative_proteins_scanned": int(nseq),
        "qualifying_candidate_sized_windows": int(nhit),
        "hit_representatives_n": int(len(matched)),
        "ambiguous_candidate_windows_skipped": int(n_ambiguous_windows),
        "ambiguous_window_rule": (
            "Candidate-sized windows containing any residue outside the 20 standard "
            "amino acids ACDEFGHIKLMNPQRSTVWY are excluded because Grantham "
            "similarity cannot be assigned unambiguously."
        ),
        "conservative_grantham_max": int(conservative_cutoff),
        "max_substitutions_searched": int(max_substitutions),
    }
    (out_path.parent / "uhgp_search_qc.json").write_text(
        json.dumps(search_qc, indent=2), encoding="utf-8"
    )
    return matched


def extract_genome_alias(text: object) -> Optional[str]:
    s=str(text)
    m=MGYG_RE.search(s)
    if m: return m.group(1)
    m=GUT_RE.search(s)
    if m: return m.group(1)
    return None


def grep_membership_lines(membership: str | Path, hit_reps: Set[str], tmpdir: Path) -> Optional[Path]:
    if not hit_reps: return None
    grep=shutil.which("grep")
    if grep is None or str(membership).endswith(".gz"): return None
    pats=tmpdir/"hit_rep_ids.txt"; pats.write_text("\n".join(sorted(hit_reps))+"\n",encoding="utf-8")
    out=tmpdir/"uhgp_membership_hit_lines.tsv"
    env=os.environ.copy(); env["LC_ALL"]="C"
    log("Prefiltering UHGP membership table with GNU grep -F...")
    with open(out,"wb") as fh:
        p=subprocess.run([grep,"-F","-f",str(pats),str(membership)],stdout=fh,stderr=subprocess.PIPE,env=env)
    if p.returncode not in (0,1):
        log("WARNING: grep prefilter failed; using Python scan.")
        out.unlink(missing_ok=True); return None
    return out


def parse_membership_for_hits(path: str | Path, hit_reps: Set[str], selected_genomes: Set[str], alias_map: Dict[str,str]) -> Tuple[pd.DataFrame, Set[str], Set[str]]:
    rows=set(); seen_reps=set(); unmapped_aliases=set()
    with open_text(path,"rt") as fh:
        for line in fh:
            if not line.strip(): continue
            fields=line.rstrip("\r\n").split("\t")
            tokens=[]
            for f in fields:
                tokens.extend([x for x in re.split(r"[,; ]+",f.strip()) if x])
            reps=[t for t in tokens if t in hit_reps]
            if not reps:
                reps=[f.split()[0] for f in fields if f.split() and f.split()[0] in hit_reps]
            if not reps: continue
            seen_reps.update(reps)
            gids=set()
            for t in tokens + reps:
                alias=extract_genome_alias(t)
                if not alias: continue
                current=alias_map.get(alias) or alias_map.get(base_mgyg(alias))
                if current is None:
                    unmapped_aliases.add(alias); continue
                if current in selected_genomes:
                    gids.add(current)
            for rep in set(reps):
                # fresh copy per representative: avoids cross-representative mutation.
                rep_gids=set(gids)
                alias=extract_genome_alias(rep)
                if alias:
                    current=alias_map.get(alias) or alias_map.get(base_mgyg(alias))
                    if current in selected_genomes:
                        rep_gids.add(current)
                for g in rep_gids:
                    rows.add((rep,g))
    return pd.DataFrame(sorted(rows),columns=["protein_rep_id","genome_id"]), seen_reps, unmapped_aliases


def map_hits(args, hit_reps: Set[str], selected_genomes_df: pd.DataFrame, uhgg_raw: pd.DataFrame) -> pd.DataFrame:
    outdir=Path(args.outdir)
    selected=set(selected_genomes_df["genome_id"].astype(str))
    alias_map=build_genome_alias_map(uhgg_raw)
    if not hit_reps:
        out=pd.DataFrame(columns=["protein_rep_id","genome_id"])
        out.to_csv(outdir/"uhgp_hit_rep_to_selected_genome.tsv.gz",sep="\t",index=False,compression="gzip")
        return out
    tmp=ensure_dir(outdir/"_tmp_membership")
    pre=grep_membership_lines(args.uhgp_membership,hit_reps,tmp)
    source=pre if pre is not None else args.uhgp_membership
    repmap,seen,unmapped=parse_membership_for_hits(source,hit_reps,selected,alias_map)
    if len(seen)==0:
        die("UHGP search produced hits, but none of the hit representative IDs were found in the membership table. Check that uhgp-100.faa and uhgp-100.tsv are from the same release.")
    seen_fraction=len(seen)/len(hit_reps)
    if seen_fraction < 0.5:
        log(f"WARNING: only {seen_fraction:.1%} of hit UHGP representatives were found in membership lines; check file compatibility.")
    if len(repmap)==0 and len(unmapped)>0:
        die(
            "Hit representatives were found in the UHGP membership table, but genome identifiers could not be resolved to authoritative UHGG accessions. "
            "This is consistent with an identifier-namespace mismatch (for example legacy GUT_GENOME IDs without an original-to-MGYG crosswalk in genomes-all_metadata.tsv). "
            "Check the detected original-accession column and ensure the UHGP and UHGG files are from the same release."
        )
    if len(repmap)==0:
        log("WARNING: hit representatives were found in the membership table, but none mapped to genomes in the selected species panel.")
    repmap.to_csv(outdir/"uhgp_hit_rep_to_selected_genome.tsv.gz",sep="\t",index=False,compression="gzip")
    qc={
        "hit_representatives_n":len(hit_reps),
        "hit_representatives_seen_in_membership_n":len(seen),
        "hit_representatives_seen_fraction":seen_fraction,
        "selected_genomes_mapped_n":int(repmap["genome_id"].nunique()) if len(repmap) else 0,
        "unmapped_genome_alias_examples":sorted(unmapped)[:50],
        "unmapped_genome_alias_n":len(unmapped),
    }
    (outdir/"membership_mapping_qc.json").write_text(json.dumps(qc,indent=2),encoding="utf-8")
    return repmap


def run_search_and_breadth(args) -> None:
    outdir=ensure_dir(args.outdir)
    peptides=pd.read_csv(outdir/"candidate_peptides.tsv",sep="\t")
    selected=pd.read_csv(outdir/"selected_species.tsv",sep="\t")
    genomes=pd.read_csv(outdir/"selected_uhgg_genomes.tsv.gz",sep="\t")
    uhgg_raw,_=read_uhgg_metadata(args.uhgg_metadata)

    hits_gz=outdir/"uhgp_candidate_hits.tsv.gz"
    if args.reuse and hits_gz.exists():
        hits=pd.read_csv(hits_gz,sep="\t")
        hit_reps=set(hits["protein_rep_id"].astype(str))
        log("Reusing existing UHGP candidate hit table.")
    else:
        tmp_hits=outdir/"uhgp_candidate_hits.tsv"
        threads=args.threads if args.threads>0 else int(os.environ.get("SLURM_CPUS_PER_TASK","1"))
        hit_reps=scan_uhgp(peptides,args.uhgp_faa,tmp_hits,threads,args.batch_sequences,args.progress_every,args.conservative_grantham_max)
        hits=pd.read_csv(tmp_hits,sep="\t")
        hits.to_csv(hits_gz,sep="\t",index=False,compression="gzip")
        tmp_hits.unlink(missing_ok=True)

    map_gz=outdir/"uhgp_hit_rep_to_selected_genome.tsv.gz"
    if args.reuse and map_gz.exists():
        repmap=pd.read_csv(map_gz,sep="\t")
        log("Reusing existing representative-to-genome mapping.")
    else:
        repmap=map_hits(args,hit_reps,genomes,uhgg_raw)

    if len(hits) and len(repmap):
        x=hits.merge(repmap,on="protein_rep_id",how="inner").merge(
            genomes[["genome_id","species_key","high_quality"]],on="genome_id",how="inner"
        )
    else:
        x=pd.DataFrame(columns=list(hits.columns)+["genome_id","species_key","high_quality"])

    # Genome-level presence per tier: any qualifying hit in any protein is sufficient.
    tier_cols=[t[0] for t in TIERS]
    genome_rows=[]
    if len(x):
        grp=x.groupby(["peptide_id","genome_id","species_key","high_quality"],sort=False)
        for key,g in grp:
            rec={"peptide_id":key[0],"genome_id":key[1],"species_key":key[2],"high_quality":bool(key[3])}
            for tier in tier_cols:
                rec[tier]=bool(pd.to_numeric(g[tier],errors="coerce").fillna(0).astype(bool).any())
            genome_rows.append(rec)
    genome_presence=pd.DataFrame(genome_rows)
    genome_presence.to_csv(outdir/"peptide_genome_tier_presence.tsv.gz",sep="\t",index=False,compression="gzip")

    denom=genomes.groupby("species_key").agg(
        n_uhgg_genomes_all=("genome_id","nunique"),
        n_uhgg_genomes_hq=("high_quality","sum"),
        mean_completeness=("completeness","mean"),
        species_rep=("species_rep","first"),
        cluster_mapping_purity=("mapping_purity","first"),
    ).reset_index()
    selcols=["species_key","mean_abundance","prevalence_ge_abundant","n_bioprojects_ge_abundant"]
    denom=selected[selcols].merge(denom,on="species_key",how="left")
    denom[["n_uhgg_genomes_all","n_uhgg_genomes_hq"]]=denom[["n_uhgg_genomes_all","n_uhgg_genomes_hq"]].fillna(0).astype(int)

    rows=[]
    for p in peptides.itertuples(index=False):
        gp=genome_presence[genome_presence["peptide_id"].eq(p.peptide_id)] if len(genome_presence) else genome_presence
        for d in denom.itertuples(index=False):
            gs=gp[gp["species_key"].eq(d.species_key)] if len(gp) else gp
            rec={
                "peptide_id":p.peptide_id,"sequence":p.sequence,"length":p.length,"species_key":d.species_key,
                "mean_abundance":d.mean_abundance,"prevalence_ge_abundant":d.prevalence_ge_abundant,
                "n_bioprojects_ge_abundant":d.n_bioprojects_ge_abundant,
                "n_uhgg_genomes_all":int(d.n_uhgg_genomes_all),"n_uhgg_genomes_hq":int(d.n_uhgg_genomes_hq),
                "mean_completeness":d.mean_completeness,"species_rep":d.species_rep,"cluster_mapping_purity":d.cluster_mapping_purity,
            }
            for tier in tier_cols:
                if len(gs):
                    ids=set(gs.loc[gs[tier].eq(True),"genome_id"])
                    hqids=set(gs.loc[gs["high_quality"].eq(True) & gs[tier].eq(True),"genome_id"])
                else:
                    ids=set(); hqids=set()
                rec[f"hit_genomes_all_{tier}"]=len(ids)
                rec[f"breadth_all_{tier}"]=len(ids)/d.n_uhgg_genomes_all if d.n_uhgg_genomes_all else np.nan
                rec[f"hit_genomes_hq_{tier}"]=len(hqids)
                rec[f"breadth_hq_{tier}"]=len(hqids)/d.n_uhgg_genomes_hq if d.n_uhgg_genomes_hq else np.nan
                rec[f"carrier_broad_hq_{tier}"]=bool(d.n_uhgg_genomes_hq >= args.min_hq_genomes_for_broad and rec[f"breadth_hq_{tier}"] >= args.broad_breadth)
            rows.append(rec)
    breadth=pd.DataFrame(rows)
    breadth.to_csv(outdir/"peptide_species_breadth.tsv.gz",sep="\t",index=False,compression="gzip")

    # Compact examples for sanity checking; deduplicate by peptide/protein representative.
    if len(hits):
        examples=hits.sort_values(["peptide_id","n_substitutions","max_grantham","protein_rep_id"]).groupby("peptide_id",group_keys=False).head(args.sanity_hits_per_peptide)
    else:
        examples=hits.copy()
    examples.to_csv(outdir/"uhgp_hit_examples.tsv",sep="\t",index=False)

    qc={
        "script_version":VERSION,
        "search_max_total_substitutions":3,
        "conservative_grantham_max":args.conservative_grantham_max,
        "broad_hq_fraction":args.broad_breadth,
        "min_hq_genomes_for_broad":args.min_hq_genomes_for_broad,
        "reference_genome_breadth_is_not_strain_prevalence":True,
        "uhgp_cluster_expansion_limitation":(
            "UHGP-100 cluster membership can include shorter members because clustering used 100% identity with an 80% target-coverage threshold. "
            "The distributed membership table does not include every member sequence, so per-genome epitope presence cannot be re-verified without individual genome protein sequences. "
            "Headline carrier calls therefore use high-quality genomes and a broad within-species threshold, but a residual upward bias is possible for epitopes located outside truncated member sequences."
        ),
    }
    (outdir/"search_breadth_qc.json").write_text(json.dumps(qc,indent=2),encoding="utf-8")
    log("UHGP search and species breadth calculation complete.")


def summarize_coverage(args) -> None:
    outdir=ensure_dir(args.outdir)
    peptides=pd.read_csv(outdir/"candidate_peptides.tsv",sep="\t")
    breadth=pd.read_csv(outdir/"peptide_species_breadth.tsv.gz",sep="\t")
    mat=pd.read_csv(outdir/"selected_species_abundance_matrix.tsv.gz",sep="\t",index_col=0)
    selected=pd.read_csv(outdir/"selected_species.tsv",sep="\t")
    prep=json.loads((outdir/"preparation_diagnostics.json").read_text())

    # Headline metric: continuous high-quality reference-genome breadth weighting.
    # For species s in sample j, contribution = abundance(j,s) * breadth_hq(s).
    # Species without an evaluable HQ breadth contribute zero. This avoids the
    # discontinuity of a hard 75% carrier cutoff while retaining the exact breadth
    # estimate in the detailed outputs.
    sample_rows=[]; ext_rows=[]; threshold_rows=[]
    for p in peptides.itertuples(index=False):
        b=breadth[breadth["peptide_id"].eq(p.peptide_id)].copy()
        er={"peptide_id":p.peptide_id,"sequence":p.sequence,"length":int(p.length)}
        qr={"peptide_id":p.peptide_id,"sequence":p.sequence,"length":int(p.length)}
        for tier,label in TIERS:
            breadth_col=f"breadth_hq_{tier}"
            if breadth_col not in b.columns:
                die(f"Missing expected species-breadth column: {breadth_col}")
            weights=(b.set_index("species_key")[breadth_col]
                       .pipe(pd.to_numeric,errors="coerce")
                       .fillna(0.0).clip(lower=0.0,upper=1.0))
            cols=[c for c in mat.columns if c in weights.index]
            if cols:
                w=weights.reindex(cols).fillna(0.0)
                weighted_cov=mat[cols].mul(w,axis=1).sum(axis=1)
            else:
                weighted_cov=pd.Series(0.0,index=mat.index,dtype=float)

            er[f"median_coverage_{tier}"]=float(weighted_cov.median())
            er[f"mean_coverage_{tier}"]=float(weighted_cov.mean())
            er[f"abundance_weighted_mean_species_breadth_{tier}"]=(
                float((selected.set_index("species_key")["mean_abundance"].reindex(weights.index).fillna(0.0) * weights).sum())
            )

            # Secondary/QC view retained for continuity with the previous analysis:
            # all-or-none species carriage at the configured broad-breadth threshold.
            carrier=set(b.loc[b[f"carrier_broad_hq_{tier}"].eq(True),"species_key"])
            broad_cols=[c for c in mat.columns if c in carrier]
            broad_cov=mat[broad_cols].sum(axis=1) if broad_cols else pd.Series(0.0,index=mat.index,dtype=float)
            qr[f"median_thresholded_coverage_{tier}"]=float(broad_cov.median())
            qr[f"n_broad_carrier_species_{tier}"]=len(carrier)

            for sample,val in weighted_cov.items():
                sample_rows.append({
                    "peptide_id":p.peptide_id,"sample_name":sample,"tier":tier,
                    "breadth_weighted_coverage":float(val),
                    "broad_hq_thresholded_coverage":float(broad_cov.loc[sample]),
                })
        ext_rows.append(er)
        threshold_rows.append(qr)

    ext=pd.DataFrame(ext_rows)
    concise=ext[["peptide_id","sequence","median_coverage_exact","median_coverage_cons_le1","median_coverage_cons_le2"]].copy()
    concise.columns=["peptide_id","sequence","median_coverage_exact","median_coverage_le1_conservative","median_coverage_le2_conservative"]
    concise.to_csv(outdir/"concise_summary.tsv",sep="\t",index=False)

    ext_keep=ext[[
        "peptide_id","sequence","length",
        "median_coverage_exact","median_coverage_cons_le1","median_coverage_cons_le2","median_coverage_cons_le3",
        "median_coverage_unrestricted_le1","median_coverage_unrestricted_le2",
    ]].copy()
    ext_keep=ext_keep.rename(columns={
        "median_coverage_cons_le1":"median_coverage_le1_conservative",
        "median_coverage_cons_le2":"median_coverage_le2_conservative",
        "median_coverage_cons_le3":"median_coverage_le3_conservative",
        "median_coverage_unrestricted_le1":"median_coverage_le1_nonconservative_allowed",
        "median_coverage_unrestricted_le2":"median_coverage_le2_nonconservative_allowed",
    })
    ext_keep.to_csv(outdir/"extended_summary.tsv",sep="\t",index=False)

    threshold_qc=pd.DataFrame(threshold_rows)
    threshold_qc.to_csv(outdir/"broad_hq_threshold_qc.tsv",sep="\t",index=False)
    sample_df=pd.DataFrame(sample_rows)
    sample_df.to_csv(outdir/"peptide_sample_coverage.tsv.gz",sep="\t",index=False,compression="gzip")

    # Run information, including parameters requested for publication/QC.
    run_info={
        **prep,
        "script_version":VERSION,
        "headline_coverage_definition":(
            "For each sample and sequence tier, sum species relative abundance multiplied by the fraction of "
            "high-quality UHGG reference genomes of that species associated with a qualifying variant."
        ),
        "headline_formula":"coverage_j = sum_s abundance_j,s * breadth_hq_s",
        "headline_metric_interpretation":(
            "Approximate expected fraction of evaluated bacterial abundance carrying a qualifying variant, if UHGG high-quality "
            "reference-genome breadth is used as a proxy for within-species variant prevalence; it is not a direct cell count."
        ),
        "secondary_broad_species_carrier_definition":f"qualifying variant present in >= {args.broad_breadth:.0%} of high-quality UHGG reference genomes",
        "conservative_substitution_definition":f"Grantham distance <= {args.conservative_grantham_max}",
        "reported_tiers":[label for _,label in TIERS],
        "headline_summary_tiers":["Exact match","≤1 conservative substitution","≤2 conservative substitutions"],
        "reference_genome_breadth_is_not_population_strain_prevalence":True,
        "python_version":platform.python_version(),
        "pandas_version":pd.__version__,
        "numpy_version":np.__version__,
        "pyahocorasick_available":ahocorasick is not None,
        "code_version":VERSION,
    }
    try:
        import openpyxl
        run_info["openpyxl_version"]=openpyxl.__version__
    except Exception:
        run_info["openpyxl_version"]="not available"
    (outdir/"run_diagnostics.json").write_text(json.dumps(run_info,indent=2),encoding="utf-8")

    write_excel(outdir,concise,ext_keep,breadth,selected,sample_df,threshold_qc,run_info)
    log("Summary stage complete.")


def write_excel(outdir: Path, concise: pd.DataFrame, extended: pd.DataFrame, breadth: pd.DataFrame, selected: pd.DataFrame, sample_df: pd.DataFrame, threshold_qc: pd.DataFrame, run_info: Dict[str,object]) -> None:
    try:
        import openpyxl
        from openpyxl.styles import Font, PatternFill, Alignment
        from openpyxl.utils import get_column_letter
    except Exception:
        log("WARNING: openpyxl not available; Excel workbooks not written. TSV outputs are complete.")
        return

    examples_path=outdir/"uhgp_hit_examples.tsv"
    if examples_path.exists() and examples_path.stat().st_size > 0:
        try:
            examples=pd.read_csv(examples_path,sep="\t")
        except pd.errors.EmptyDataError:
            examples=pd.DataFrame()
    else:
        examples=pd.DataFrame()

    qc_rows=[]
    def flatten(prefix,obj):
        if isinstance(obj,dict):
            for k,v in obj.items():
                flatten(f"{prefix}.{k}" if prefix else k,v)
        elif isinstance(obj,list):
            qc_rows.append((prefix," | ".join(map(str,obj))))
        else:
            qc_rows.append((prefix,obj))
    flatten("",run_info)
    qc=pd.DataFrame(qc_rows,columns=["parameter","value"])

    readme=pd.DataFrame([
        ["Concise summary","Headline result: median breadth-weighted gut-community coverage for exact, ≤1 conservative, and ≤2 conservative variants."],
        ["Candidate evaluation","Compact interpretation table combining headline medians, gains from conservative substitutions, unrestricted sensitivity views, and the older ≥75%-breadth QC values."],
        ["Extended summary","Adds ≤3 conservative and unrestricted ≤1/≤2 substitution views. 'Non-conservative allowed' means the total substitution count is limited but substitution chemistry is unrestricted."],
        ["Breadth-weighted coverage","For every species and sample, contribution = species relative abundance × fraction of HQ UHGG reference genomes with a qualifying variant. A 20% abundant species with 75% breadth contributes 15%."],
        ["Species breadth","Per candidate and species, fraction of authoritative UHGG genomes and high-quality UHGG genomes linked to a qualifying UHGP-100 protein cluster."],
        ["Selected species","Ecologically selected bacterial species, abundance statistics, and authoritative UHGG genome counts."],
        ["Per-sample coverage","Underlying breadth-weighted community coverage for every candidate/sample/tier, plus the older ≥75%-thresholded QC value."],
        ["≥75% breadth QC",f"Secondary all-or-none sensitivity view: a species counts fully only when the qualifying variant occurs in ≥{run_info['secondary_broad_species_carrier_definition'].split('>=')[-1].split('of')[0].strip()} of its HQ genomes."],
        ["UHGP hit examples","Representative sequence windows and substitution chemistry for sanity checking."],
        ["Run & QC","Input counts, thresholds, software versions, panel abundance ceiling, and methodological caveats."],
        ["Conservative substitution",f"Defined as {run_info['conservative_substitution_definition']}."],
        ["Important interpretation","Breadth-weighted coverage is an estimate based on reference-genome breadth, not a direct measurement of the fraction of bacterial cells carrying the peptide. UHGG genomes are not a random population sample."],
        ["Important caveat","UHGP-100 cluster expansion may over-credit truncated cluster members because full member sequences are not distributed in uhgp-100.tsv; this can slightly inflate breadth estimates."],
        ["Antibody interpretation","The similarity tiers are alternative sequence-similarity views, not direct probabilities of antibody binding. Antibody recognition remains position- and context-dependent."],
    ],columns=["item","explanation"])

    concise_excel = concise.rename(columns={
        "peptide_id":"Candidate", "sequence":"Sequence",
        "median_coverage_exact":"Median breadth-weighted coverage – exact",
        "median_coverage_le1_conservative":"Median breadth-weighted coverage – ≤1 conservative",
        "median_coverage_le2_conservative":"Median breadth-weighted coverage – ≤2 conservative",
    })
    extended_excel = extended.rename(columns={
        "peptide_id":"Candidate", "sequence":"Sequence", "length":"Length (aa)",
        "median_coverage_exact":"Median breadth-weighted coverage – exact",
        "median_coverage_le1_conservative":"Median breadth-weighted coverage – ≤1 conservative",
        "median_coverage_le2_conservative":"Median breadth-weighted coverage – ≤2 conservative",
        "median_coverage_le3_conservative":"Median breadth-weighted coverage – ≤3 conservative",
        "median_coverage_le1_nonconservative_allowed":"Median breadth-weighted coverage – ≤1 substitution (non-conservative allowed)",
        "median_coverage_le2_nonconservative_allowed":"Median breadth-weighted coverage – ≤2 substitutions (non-conservative allowed)",
    })
    threshold_excel=threshold_qc.rename(columns={
        "peptide_id":"Candidate","sequence":"Sequence","length":"Length (aa)",
        **{f"median_thresholded_coverage_{tier}":f"Median ≥75%-breadth coverage – {label}" for tier,label in TIERS},
        **{f"n_broad_carrier_species_{tier}":f"Broad carrier species – {label}" for tier,label in TIERS},
    })

    # Candidate-level evaluation: readable comparison without losing the detailed audit sheets.
    candidate_eval=extended.copy()
    candidate_eval["gain_exact_to_le1_conservative"]=(
        candidate_eval["median_coverage_le1_conservative"]-candidate_eval["median_coverage_exact"]
    )
    candidate_eval["gain_exact_to_le2_conservative"]=(
        candidate_eval["median_coverage_le2_conservative"]-candidate_eval["median_coverage_exact"]
    )
    candidate_eval["gain_le1_cons_to_unrestricted"]=(
        candidate_eval["median_coverage_le1_nonconservative_allowed"]-candidate_eval["median_coverage_le1_conservative"]
    )
    candidate_eval["gain_le2_cons_to_unrestricted"]=(
        candidate_eval["median_coverage_le2_nonconservative_allowed"]-candidate_eval["median_coverage_le2_conservative"]
    )
    tq=threshold_qc[[
        "peptide_id",
        "median_thresholded_coverage_exact",
        "median_thresholded_coverage_cons_le1",
        "median_thresholded_coverage_cons_le2",
    ]].copy()
    candidate_eval=candidate_eval.merge(tq,on="peptide_id",how="left")
    candidate_eval=candidate_eval.rename(columns={
        "peptide_id":"Candidate","sequence":"Sequence","length":"Length (aa)",
        "median_coverage_exact":"Median BW coverage – exact",
        "median_coverage_le1_conservative":"Median BW coverage – ≤1 conservative",
        "median_coverage_le2_conservative":"Median BW coverage – ≤2 conservative",
        "median_coverage_le3_conservative":"Median BW coverage – ≤3 conservative",
        "median_coverage_le1_nonconservative_allowed":"Median BW coverage – ≤1 unrestricted",
        "median_coverage_le2_nonconservative_allowed":"Median BW coverage – ≤2 unrestricted",
        "gain_exact_to_le1_conservative":"Gain: exact → ≤1 conservative",
        "gain_exact_to_le2_conservative":"Gain: exact → ≤2 conservative",
        "gain_le1_cons_to_unrestricted":"Extra gain from non-conservative allowance – ≤1",
        "gain_le2_cons_to_unrestricted":"Extra gain from non-conservative allowance – ≤2",
        "median_thresholded_coverage_exact":"Median ≥75%-breadth QC – exact",
        "median_thresholded_coverage_cons_le1":"Median ≥75%-breadth QC – ≤1 conservative",
        "median_thresholded_coverage_cons_le2":"Median ≥75%-breadth QC – ≤2 conservative",
    })

    def style_workbook(wb):
        header_fill=PatternFill("solid",fgColor="1F4E78")
        header_font=Font(color="FFFFFF",bold=True)
        light_fill=PatternFill("solid",fgColor="D9EAF7")
        for ws in wb.worksheets:
            ws.freeze_panes="A2"
            for cell in ws[1]:
                cell.fill=header_fill
                cell.font=header_font
                cell.alignment=Alignment(wrap_text=True,vertical="top")
            if ws.max_row >= 1 and ws.max_column >= 1:
                ws.auto_filter.ref=ws.dimensions
            for col_cells in ws.columns:
                cells=list(col_cells)
                vals=[str(c.value) if c.value is not None else "" for c in cells[:200]]
                width=min(max(10,max((len(v) for v in vals),default=10)+2),48)
                ws.column_dimensions[get_column_letter(cells[0].column)].width=width
            for row in ws.iter_rows(min_row=2):
                for c in row:
                    if isinstance(c.value,(int,float)):
                        hdr=str(ws.cell(row=1,column=c.column).value or "").lower()
                        if "coverage" in hdr or "gain" in hdr:
                            c.number_format="0.0%"
            if ws.title in {"README","Run & QC"}:
                for row in ws.iter_rows(min_row=2):
                    for c in row:
                        c.alignment=Alignment(wrap_text=True,vertical="top")
        if "Concise summary" in wb.sheetnames:
            ws=wb["Concise summary"]
            for cell in ws[1]:
                cell.fill=PatternFill("solid",fgColor="17365D")
        if "Candidate evaluation" in wb.sheetnames:
            ws=wb["Candidate evaluation"]
            for cell in ws[1]:
                cell.fill=PatternFill("solid",fgColor="244062")
            # Visually separate the secondary ≥75%-breadth QC columns.
            for c in range(14,17):
                if c <= ws.max_column:
                    ws.cell(row=1,column=c).fill=light_fill
                    ws.cell(row=1,column=c).font=Font(color="000000",bold=True)

    def write_rich_workbook(path: Path) -> None:
        with pd.ExcelWriter(path,engine="openpyxl") as writer:
            readme.to_excel(writer,sheet_name="README",index=False)
            concise_excel.to_excel(writer,sheet_name="Concise summary",index=False)
            candidate_eval.to_excel(writer,sheet_name="Candidate evaluation",index=False)
            extended_excel.to_excel(writer,sheet_name="Extended summary",index=False)
            breadth.to_excel(writer,sheet_name="Species breadth",index=False)
            selected.to_excel(writer,sheet_name="Selected species",index=False)
            sample_df.to_excel(writer,sheet_name="Per-sample coverage",index=False)
            threshold_excel.to_excel(writer,sheet_name="≥75% breadth QC",index=False)
            examples.to_excel(writer,sheet_name="UHGP hit examples",index=False)
            qc.to_excel(writer,sheet_name="Run & QC",index=False)
            style_workbook(writer.book)

    # The standalone summary is deliberately rich: headline tables first, followed by
    # the same supporting evidence/QC sheets needed to interpret them.
    summary_xlsx=outdir/"UHGG_peptide_coverage_summary.xlsx"
    write_rich_workbook(summary_xlsx)
    log(f"Wrote {summary_xlsx}")

    # Full audit workbook retained for backwards compatibility; contents are identical
    # so either workbook is sufficient for interpretation and archiving.
    xlsx=outdir/"UHGG_peptide_coverage.xlsx"
    write_rich_workbook(xlsx)
    log(f"Wrote {xlsx}")

def run_all(args) -> None:
    prepare_species_and_genomes(args)
    run_search_and_breadth(args)
    summarize_coverage(args)


def common_args(p):
    p.add_argument("--peptides",required=True,help="Candidate peptide FASTA.")
    p.add_argument("--outdir",default="results_uhgg",help="Output directory.")
    p.add_argument("--min-peptide-length",type=int,default=12)
    p.add_argument("--max-peptide-length",type=int,default=25)


def preparation_args(p):
    p.add_argument("--relative-abundance",required=True,help="UHGV read_mapping/relative_abundance.tsv")
    p.add_argument("--sample-metadata",required=True,help="UHGV read_mapping/sample_metadata.tsv")
    p.add_argument("--host-genomes-info",required=True,help="UHGV host_predictions/host_genomes_info.tsv or .zip")
    p.add_argument("--uhgg-metadata",required=True,help="UHGG v2.0.2 genomes-all_metadata.tsv; authoritative genome denominator")
    p.add_argument("--abundant-relative",type=float,default=0.01,help="Species abundance threshold within species-labelled bacterial references.")
    p.add_argument("--abundant-prevalence",type=float,default=0.01,help="Minimum fraction of bulk samples reaching --abundant-relative.")
    p.add_argument("--min-abundant-bioprojects",type=int,default=2,help="Minimum BioProjects in which abundance threshold is reached.")
    p.add_argument("--min-completeness",type=float,default=90.0)
    p.add_argument("--max-contamination",type=float,default=5.0)
    p.add_argument("--cluster-mapping-min-purity",type=float,default=0.90,help="Minimum overlap purity for assigning a UHGG species cluster to one updated GTDB-r207 species.")


def search_args(p):
    p.add_argument("--uhgp-faa",required=True,help="UHGP-100 representative protein FASTA (uhgp-100.faa).")
    p.add_argument("--uhgp-membership",required=True,help="UHGP-100 membership table (uhgp-100.tsv).")
    p.add_argument("--uhgg-metadata",required=True,help="UHGG v2.0.2 genomes-all_metadata.tsv.")
    p.add_argument("--conservative-grantham-max",type=int,default=50,help="Maximum Grantham distance classified as conservative.")
    p.add_argument("--broad-breadth",type=float,default=0.75,help="Fraction of high-quality genomes required for a species to count as a carrier.")
    p.add_argument("--min-hq-genomes-for-broad",type=int,default=1,help="Minimum number of HQ genomes required to evaluate broad carriage. Default 1 preserves singleton species; raise for sensitivity analysis if desired.")
    p.add_argument("--threads",type=int,default=0,help="0 uses SLURM_CPUS_PER_TASK or 1.")
    p.add_argument("--batch-sequences",type=int,default=10000)
    p.add_argument("--progress-every",type=int,default=1_000_000)
    p.add_argument("--sanity-hits-per-peptide",type=int,default=50)
    p.add_argument("--reuse",action="store_true",help="Reuse UHGP hit and membership mapping files. Safe only when candidate FASTA and relevant search/mapping parameters are unchanged.")


def summary_args(p):
    p.add_argument("--conservative-grantham-max",type=int,default=50)
    p.add_argument("--broad-breadth",type=float,default=0.75)


def make_parser():
    ap=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description=__doc__)
    sp=ap.add_subparsers(dest="command",required=True)
    p=sp.add_parser("prepare",formatter_class=argparse.ArgumentDefaultsHelpFormatter,help="Build abundance-defined species panel and authoritative UHGG genome denominator.")
    common_args(p); preparation_args(p); p.set_defaults(func=prepare_species_and_genomes)
    p=sp.add_parser("search",formatter_class=argparse.ArgumentDefaultsHelpFormatter,help="Search UHGP-100 and calculate high-quality reference-genome breadth.")
    common_args(p); search_args(p); p.set_defaults(func=run_search_and_breadth)
    p=sp.add_parser("summarize",formatter_class=argparse.ArgumentDefaultsHelpFormatter,help="Calculate breadth-weighted community coverage and write summaries/Excel.")
    common_args(p); summary_args(p); p.set_defaults(func=summarize_coverage)
    p=sp.add_parser("all",formatter_class=argparse.ArgumentDefaultsHelpFormatter,help="Run prepare + search + summarize.")
    common_args(p); preparation_args(p)
    # search_args without duplicated uhgg-metadata
    p.add_argument("--uhgp-faa",required=True)
    p.add_argument("--uhgp-membership",required=True)
    p.add_argument("--conservative-grantham-max",type=int,default=50)
    p.add_argument("--broad-breadth",type=float,default=0.75)
    p.add_argument("--min-hq-genomes-for-broad",type=int,default=1)
    p.add_argument("--threads",type=int,default=0)
    p.add_argument("--batch-sequences",type=int,default=10000)
    p.add_argument("--progress-every",type=int,default=1_000_000)
    p.add_argument("--sanity-hits-per-peptide",type=int,default=50)
    p.add_argument("--reuse",action="store_true")
    p.set_defaults(func=run_all)
    return ap


def main(argv=None):
    ap=make_parser(); args=ap.parse_args(argv)
    log(f"gut_peptide_coverage_uhgg_v2.py version {VERSION}")
    args.func(args)


if __name__=="__main__":
    main()
