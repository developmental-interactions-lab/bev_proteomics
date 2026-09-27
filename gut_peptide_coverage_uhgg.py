#!/usr/bin/env python3
"""
Gut peptide coverage, version 3.1.0.

Select a shared species panel covering a target cumulative mean bacterial abundance
(default 95%), search UHGP-100 representatives for exact/conservative/unrestricted
peptide variants, and combine evaluable reference-genome breadth with UHGV abundances.
Run --help or see README_gut_peptide_coverage_uhgg.md for interpretation and provenance.
This estimates encoded sequence carriage, not protein expression or antibody binding.
All FASTA records are processed identically; no control-type metadata is consumed.
Original workflow developed using ChatGPT 5.6 Sol; revised 2026-09-27.
"""

from __future__ import annotations

import argparse
import hashlib
import time
import csv
import gzip
import json
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

VERSION = "3.1.2"
CACHE_SCHEMA = 3
TAXONOMY_MAPPING_ALGORITHM = "direct_uhgv_genome_taxonomy_then_pure_cluster_v1"
MAX_SUBSTITUTIONS = 3
EXCEL_DATA_ROW_LIMIT = 200000
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


def parse_fasta(path: str | Path, protein: bool = True) -> Iterator[Tuple[str, str]]:
    name = None
    seq: List[str] = []
    with open_text(path, "rt") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper().rstrip("*") if protein else "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(re.sub(r"\s+", "", line))
        if name is not None:
            yield name, "".join(seq).upper().rstrip("*") if protein else "".join(seq).upper()


def read_candidate_peptides(path: str | Path, min_len: int, max_len: int) -> pd.DataFrame:
    rows = []
    seen = set()
    for pid, seq in parse_fasta(path, protein=False):
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
    if out["genome_id"].duplicated().any():
        die("Duplicate genome IDs in UHGG metadata; resolve them explicitly.")
    detected = {
        "genome": genome_col,
        "species_rep": species_rep_col,
        "completeness": complete_col,
        "contamination": contam_col,
        "original_genome_id": original_col or "",
    }
    return out, detected


def map_uhgg_clusters_to_r207(uhgg: pd.DataFrame, host_bac: pd.DataFrame, min_purity: float) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Prefer the genome's own r207 label; infer only otherwise unlabeled genomes."""
    h = host_bac[host_bac["data_source"].astype(str).eq("UHGG")][["genome_id", "species"]].drop_duplicates("genome_id")
    h["genome_base"] = h["genome_id"].map(base_mgyg)
    u = uhgg[["genome_id", "species_rep"]].copy()
    u["genome_base"] = u["genome_id"].map(base_mgyg)

    exact = u.merge(h[["genome_id","species"]], on="genome_id", how="left")
    exact["mapping_source"] = np.where(exact["species"].notna(), "genome_exact", "unmapped")
    missing = exact["species"].isna()
    if missing.any():
        # Version-insensitive fallback is safe only when the host base ID is unique.
        hbase = h.groupby("genome_base").filter(lambda x: x["species"].nunique() == 1).drop_duplicates("genome_base")
        fallback = exact.loc[missing, ["genome_id","species_rep","genome_base"]].merge(
            hbase[["genome_base","species"]], on="genome_base", how="left"
        )
        exact.loc[missing, "species"] = fallback.set_index("genome_id").reindex(exact.loc[missing,"genome_id"])["species"].values
        exact.loc[missing & exact["species"].notna(), "mapping_source"] = "genome_base"

    x = exact[exact["species"].fillna("").ne("")].copy()
    counts = x.groupby(["species_rep","species"]).size().rename("n_overlap").reset_index()
    totals = counts.groupby("species_rep")["n_overlap"].sum().rename("n_cluster_overlap")
    counts = counts.merge(totals, on="species_rep", how="left")
    counts["mapping_purity"] = counts["n_overlap"] / counts["n_cluster_overlap"]
    best = counts.sort_values(["species_rep","n_overlap"], ascending=[True,False]).drop_duplicates("species_rep")
    best["mapping_accepted"] = best["mapping_purity"].ge(min_purity)

    direct = exact[["genome_id", "species", "mapping_source"]].rename(columns={"species":"species_key"})
    modal = best[["species_rep","species","mapping_purity","n_cluster_overlap","mapping_accepted"]].rename(
        columns={"species":"cluster_modal_species", "mapping_accepted":"cluster_fallback_accepted"})
    out = uhgg.merge(direct, on="genome_id", how="left", validate="one_to_one").merge(
        modal, on="species_rep", how="left", validate="many_to_one")
    fallback_mask = out["species_key"].isna() & out["cluster_fallback_accepted"].eq(True)
    out.loc[fallback_mask, "species_key"] = out.loc[fallback_mask, "cluster_modal_species"]
    out.loc[fallback_mask, "mapping_source"] = "cluster_fallback"
    out["mapping_source"] = out["mapping_source"].fillna("unmapped")
    out["mapping_accepted"] = out["species_key"].notna()
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


def digest_json(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str).encode()).hexdigest()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def source_signature(path, large=False):
    p = Path(path).resolve()
    if not p.is_file():
        die(f"Input file does not exist: {p}")
    stat = p.stat()
    result = {"path": str(p), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if large:
        # Avoid repeated full reads of very large catalogues. This is an identity
        # guard, not a cryptographic checksum of the complete catalogue.
        h = hashlib.sha256()
        with p.open("rb") as f:
            for offset in sorted({0, max(0, stat.st_size // 2 - 524288), max(0, stat.st_size - 1048576)}):
                f.seek(offset); h.update(f.read(1048576))
        result["sampled_sha256"] = h.hexdigest()
    else:
        result["sha256"] = sha256_file(p)
    return result


def artifact_signature(path):
    p = Path(path); st = p.stat()
    return {"size": st.st_size, "mtime_ns": st.st_mtime_ns}


def write_json(path, value):
    path = Path(path); tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")
    tmp.replace(path)


def write_table(df, path, index=False):
    path = Path(path); tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=index, compression="gzip" if str(path).endswith(".gz") else None)
    tmp.replace(path)


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save_cache(path, inputs, outputs):
    write_json(path, {"cache_schema": CACHE_SCHEMA, "script_version": VERSION,
                      "inputs": inputs, "outputs": {Path(f).name: artifact_signature(f) for f in outputs}})


def cache_matches(path, inputs):
    path = Path(path)
    if not path.exists():
        return False
    try:
        manifest = read_json(path)
        return (manifest["cache_schema"] == CACHE_SCHEMA and manifest["inputs"] == inputs
                and all((path.parent / name).is_file() and artifact_signature(path.parent / name) == sig
                        for name, sig in manifest["outputs"].items()))
    except (KeyError, ValueError, OSError):
        return False


def frame_digest(df):
    return hashlib.sha256(df.to_csv(index=False).encode()).hexdigest()


def validate_prepared(args):
    outdir = Path(args.outdir)
    mf = outdir / "prepare_manifest.json"
    if not mf.exists():
        die("No v3 preparation manifest. Run prepare (or all) in a new output directory.")
    manifest = read_json(mf)
    if not cache_matches(mf, manifest["inputs"]):
        die("Prepared files are incomplete or changed; rerun prepare.")
    if manifest["inputs"].get("parameters", {}).get("taxonomy_mapping_algorithm") != TAXONOMY_MAPPING_ALGORITHM:
        die("Prepared taxonomy predates direct genome mapping; rerun prepare or all. The compatible protein scan can be reused.")
    supplied = read_candidate_peptides(args.peptides, args.min_peptide_length, args.max_peptide_length)
    saved = pd.read_csv(outdir / "candidate_peptides.tsv", sep="\t", keep_default_na=False, dtype={"peptide_id": str})
    if supplied.to_dict("records") != saved.to_dict("records"):
        die("Candidate FASTA differs from prepared candidates. Rerun prepare or all; --reuse safely invalidates incompatible caches.")
    return saved


def coverage_stats(values):
    values = pd.Series(values, dtype=float)
    return {"mean": float(values.mean()), "median": float(values.median()),
            "p05": float(values.quantile(.05)), "min": float(values.min()), "max": float(values.max())}


def select_cumulative_panel(stats, target):
    ordered = stats.sort_values(["mean_abundance", "species_key"], ascending=[False, True]).copy()
    cumulative = ordered["mean_abundance"].cumsum().to_numpy()
    n = min(len(ordered), int(np.searchsorted(cumulative, target, side="left")) + 1)
    ordered["selected"] = np.arange(len(ordered)) < n
    ordered["cumulative_mean_abundance"] = cumulative
    return ordered


def unique_sample_projects(metadata):
    """Collapse repeated UHGV sample names only if their BioProject agrees."""
    sm = metadata[["sample_name", "bioproject_id"]].copy()
    if sm["sample_name"].isna().any() or sm["sample_name"].astype(str).str.strip().eq("").any():
        die("Sample metadata contain missing or blank sample_name values.")
    sm["sample_name"] = sm["sample_name"].astype(str)
    sm["bioproject_id"] = sm["bioproject_id"].fillna("").astype(str).str.strip()
    conflicts = sm.groupby("sample_name", sort=False)["bioproject_id"].nunique().gt(1)
    if conflicts.any():
        names = conflicts.index[conflicts].tolist()
        die(f"Conflicting bioproject_id values for {len(names)} sample_name(s) in sample metadata; "
            f"examples: {', '.join(names[:5])}. Resolve the source assignments explicitly.")
    duplicate_rows = int(sm["sample_name"].duplicated().sum())
    if duplicate_rows:
        log(f"Collapsed {duplicate_rows} repeated sample-metadata rows with matching sample_name and bioproject_id.")
    return sm.drop_duplicates(subset="sample_name"), duplicate_rows


def prepare_species_and_genomes(args) -> None:
    outdir = ensure_dir(args.outdir)
    # Invalidate completion marker first, so a failed preparation cannot be reused.
    (outdir / "prepare_manifest.json").unlink(missing_ok=True)
    peptides = read_candidate_peptides(args.peptides, args.min_peptide_length, args.max_peptide_length)
    log("Reading bacterial taxonomy and bulk-metagenome abundances...")
    host = add_taxonomy(read_host_genomes(args.host_genomes_info))
    if host["genome_id"].duplicated().any():
        die("Duplicate genome_id values in host taxonomy; resolve them explicitly.")
    host_bac = host[host["domain"].eq("Bacteria") & host["species"].ne("")].copy()
    ra = pd.read_csv(args.relative_abundance, sep="\t", low_memory=False,
                     dtype={"sample_name": str})
    need = {"sample_name", "sample_type", "genome_type", "genome", "relative_abundance"}
    if not need.issubset(ra):
        die(f"Relative abundance file lacks {sorted(need - set(ra))}")
    ra = ra[ra["sample_type"].eq("bulk") & ra["genome_type"].eq("prok")].copy()
    ra["relative_abundance"] = pd.to_numeric(ra["relative_abundance"], errors="coerce")
    if (~np.isfinite(ra["relative_abundance"]) | ra["relative_abundance"].lt(0)).any():
        die("Bulk prokaryotic abundances contain missing, non-numeric, negative or infinite values.")
    ra["sample_name"] = ra["sample_name"].astype(str)
    all_prok = ra.groupby("sample_name")["relative_abundance"].sum()
    ra = ra.merge(host_bac[["genome_id", "species"]], left_on="genome", right_on="genome_id", how="inner", validate="many_to_one")
    ra = ra.groupby(["sample_name", "species"], as_index=False)["relative_abundance"].sum()
    labelled = ra.groupby("sample_name")["relative_abundance"].sum()
    totals = ra.groupby("sample_name")["relative_abundance"].transform("sum")
    ra = ra[totals.gt(0)].copy()
    ra["abundance"] = ra["relative_abundance"] / totals[totals.gt(0)]
    if ra.empty:
        die("No positive species-labelled bacterial abundances remain.")
    sample_ids = sorted(ra["sample_name"].unique()); nsamples = len(sample_ids)
    sm = pd.read_csv(args.sample_metadata, sep="\t", dtype={"sample_name": str, "bioproject_id": str})
    if not {"sample_name", "bioproject_id"}.issubset(sm):
        die("sample_metadata.tsv requires sample_name and bioproject_id.")
    sm, duplicate_metadata_rows = unique_sample_projects(sm)
    samples = pd.DataFrame({"sample_name": sample_ids}).merge(sm, on="sample_name", how="left", validate="one_to_one")
    samples["bioproject_id"] = samples["bioproject_id"].fillna("").astype(str).str.strip()
    samples["labelled_bacterial_fraction_of_input_prokaryotic_abundance"] = samples["sample_name"].map(labelled / all_prok.replace(0, np.nan))
    ra = ra.merge(samples[["sample_name", "bioproject_id"]], on="sample_name", validate="many_to_one")
    stats = ra.groupby("species").agg(abundance_sum=("abundance", "sum"), max_abundance=("abundance", "max")).reset_index().rename(columns={"species": "species_key"})
    stats["mean_abundance"] = stats["abundance_sum"] / nsamples
    detected = ra[ra["abundance"].gt(0)].groupby("species")["sample_name"].nunique()
    projects = ra[ra["abundance"].gt(0) & ra["bioproject_id"].ne("")].groupby("species")["bioproject_id"].nunique()
    stats["prevalence_any"] = stats["species_key"].map(detected).fillna(0) / nsamples
    stats["n_bioprojects_detected"] = stats["species_key"].map(projects).fillna(0).astype(int)
    stats = select_cumulative_panel(stats, args.target_community_coverage)
    selected = stats[stats["selected"]].copy()
    uhgg_raw, detected_cols = read_uhgg_metadata(args.uhgg_metadata)
    mapped, cluster_map = map_uhgg_clusters_to_r207(uhgg_raw, host_bac, args.cluster_mapping_min_purity)
    genomes = mapped[mapped["species_key"].isin(selected["species_key"])].copy()
    genomes["high_quality"] = genomes["completeness"].ge(args.min_completeness) & genomes["contamination"].le(args.max_contamination)
    counts = genomes.groupby("species_key").agg(n_uhgg_genomes_all=("genome_id", "nunique"), n_uhgg_genomes_hq=("high_quality", "sum"),
        cluster_mapping_purity_min=("mapping_purity", "min"), cluster_overlap_min=("n_cluster_overlap", "min"),
        n_uhgg_species_clusters=("species_rep", "nunique")).reset_index()
    source_counts = genomes.groupby(["species_key", "mapping_source"]).size().unstack(fill_value=0)
    for source in ("genome_exact", "genome_base", "cluster_fallback"):
        counts[f"n_mapped_{source}"] = counts["species_key"].map(source_counts[source] if source in source_counts else {}).fillna(0).astype(int)
    selected = selected.merge(counts, on="species_key", how="left", validate="one_to_one")
    for c in ["n_uhgg_genomes_all", "n_uhgg_genomes_hq", "n_uhgg_species_clusters",
              "n_mapped_genome_exact", "n_mapped_genome_base", "n_mapped_cluster_fallback"]:
        selected[c] = selected[c].fillna(0).astype(int)
    mat = ra[ra["species"].isin(selected["species_key"])].pivot_table(index="sample_name", columns="species", values="abundance", aggfunc="sum", fill_value=0.0)
    mat = mat.reindex(index=sample_ids, columns=selected["species_key"].tolist(), fill_value=0.0)
    sample_qc = pd.DataFrame(index=mat.index)
    sample_qc["panel_abundance"] = mat.sum(axis=1)
    sample_qc["outside_panel_abundance"] = 1 - sample_qc["panel_abundance"]
    no_map = selected.loc[selected["n_uhgg_genomes_all"].eq(0), "species_key"]
    no_hq = selected.loc[selected["n_uhgg_genomes_all"].gt(0) & selected["n_uhgg_genomes_hq"].eq(0), "species_key"]
    sample_qc["panel_no_mapped_genomes_abundance"] = mat[no_map].sum(axis=1)
    sample_qc["panel_no_hq_genomes_abundance"] = mat[no_hq].sum(axis=1)
    sample_qc["hq_species_abundance_before_membership_qc"] = mat[selected.loc[selected["n_uhgg_genomes_hq"].gt(0), "species_key"]].sum(axis=1)
    hq_all = mapped["completeness"].ge(args.min_completeness) & mapped["contamination"].le(args.max_contamination) & mapped["mapping_accepted"]
    available_species = set(mapped.loc[hq_all, "species_key"])
    attainable = ra[ra["species"].isin(available_species)].groupby("sample_name")["abundance"].sum().reindex(sample_ids, fill_value=0)
    diag = {"script_version": VERSION, "bulk_metagenome_samples_n": nsamples, "selected_species_n": len(selected),
        "target_community_coverage": args.target_community_coverage,
        "selection_rule": "Smallest shared panel reaching target cumulative mean abundance; zeros included; ties ordered by species key.",
        "panel_coverage": coverage_stats(sample_qc["panel_abundance"]),
        "all_species_hq_mappable_abundance_before_membership_qc": coverage_stats(attainable),
        "samples_excluded_no_labelled_bacterial_abundance": int(len(all_prok) - nsamples),
        "samples_missing_bioproject": int(samples["bioproject_id"].eq("").sum()),
        "sample_metadata_duplicate_rows_collapsed": duplicate_metadata_rows,
        "high_quality_definition": {"min_completeness": args.min_completeness, "max_contamination": args.max_contamination},
        "cluster_mapping_min_purity": args.cluster_mapping_min_purity,
        "taxonomy_mapping_algorithm": TAXONOMY_MAPPING_ALGORITHM,
        "mapped_genomes_by_source": mapped["mapping_source"].value_counts().to_dict(),
        "uhgg_metadata_detected_columns": detected_cols,
        "denominator": "Species-labelled bacterial abundance assigned to host references, renormalized within each sample. Archaea, unlabelled references and unmapped community are outside this denominator.",
        "hq_definition_note": "Operational completeness/contamination filter, not the full MIMAG high-quality definition."}
    outputs = []
    for df, name, index in [(peptides,"candidate_peptides.tsv",False),(stats,"all_species_abundance.tsv",False),
        (selected,"selected_species.tsv",False),(genomes,"selected_uhgg_genomes.tsv.gz",False),
        (mat,"selected_species_abundance_matrix.tsv.gz",True),(samples,"included_sample_metadata.tsv",False),
        (sample_qc,"preparation_sample_qc.tsv.gz",True),(cluster_map,"uhgg_cluster_to_gtdb_r207_mapping.tsv.gz",False)]:
        write_table(df, outdir / name, index); outputs.append(outdir / name)
    write_json(outdir / "preparation_diagnostics.json", diag); outputs.append(outdir / "preparation_diagnostics.json")
    inputs = {"sources": {k: source_signature(getattr(args,k)) for k in ["peptides","relative_abundance","sample_metadata","host_genomes_info","uhgg_metadata"]},
              "parameters": {**{k: getattr(args,k) for k in ["target_community_coverage","min_completeness","max_contamination","cluster_mapping_min_purity"]},
                             "taxonomy_mapping_algorithm": TAXONOMY_MAPPING_ALGORITHM}}
    save_cache(outdir / "prepare_manifest.json", inputs, outputs)
    log(f"Selected {len(selected):,} species; mean panel abundance {sample_qc.panel_abundance.mean():.2%}; 5th percentile {sample_qc.panel_abundance.quantile(.05):.2%}.")


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
    max_substitutions=MAX_SUBSTITUTIONS
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
    if p.returncode == 1:
        log("No literal hit-representative IDs found by grep; checking the full membership table for diagnosis.")
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


def inventory_membership(args, uhgg_raw):
    outdir = Path(args.outdir)
    inputs = {"membership": source_signature(args.uhgp_membership, large=True),
              "uhgg_metadata": source_signature(args.uhgg_metadata)}
    manifest = outdir / "membership_inventory_manifest.json"
    table = outdir / "uhgp_genome_inventory.tsv.gz"
    qc_path = outdir / "membership_inventory_qc.json"
    if cache_matches(manifest, inputs):
        log("Reusing validated, peptide-independent membership inventory.")
        return pd.read_csv(table, sep="\t"), read_json(qc_path), inputs
    manifest.unlink(missing_ok=True)
    aliases = build_genome_alias_map(uhgg_raw)
    pattern = re.compile(r"MGYG\d{9}(?:\.\d+)?|GUT_GENOME\d{6,9}")
    seen = set(); represented = set(); unknown = set(); no_alias_examples = []
    nlines = n_no_alias = 0
    log("Inventorying genome identifiers across the complete membership table (cached for later runs)...")
    with open_text(args.uhgp_membership) as fh:
        for line in fh:
            if not line.strip():
                continue
            nlines += 1; found = False
            for match in pattern.finditer(line):
                found = True; alias = match.group(0)
                if alias in seen:
                    continue
                seen.add(alias)
                current = aliases.get(alias) or aliases.get(base_mgyg(alias))
                if current is None:
                    unknown.add(alias)
                else:
                    represented.add(current)
            if not found:
                n_no_alias += 1
                if len(no_alias_examples) < 10:
                    no_alias_examples.append(line.strip()[:200])
            if nlines % 10000000 == 0:
                log(f"Membership inventory: {nlines:,} lines; {len(represented):,} genomes resolved.")
    if not represented:
        die("No UHGG genomes resolved anywhere in membership table; check namespace and metadata release.")
    inventory = pd.DataFrame({"genome_id": sorted(uhgg_raw["genome_id"].astype(str).unique())})
    inventory["in_uhgp_membership"] = inventory["genome_id"].isin(represented)
    qc = {"nonempty_membership_lines": nlines, "resolved_genomes_n": len(represented),
          "unresolved_recognized_alias_n": len(unknown), "unresolved_alias_examples": sorted(unknown)[:50],
          "lines_without_recognized_alias_n": n_no_alias, "lines_without_recognized_alias_examples": no_alias_examples,
          "interpretation": "At least one protein identifier represents the genome; this does not prove complete gene recovery. Only MGYG/GUT_GENOME aliases are recognized. Lines lacking recognized identifiers are reported; mixed-format lines may also contain unrecognized identifiers."}
    write_table(inventory, table); write_json(qc_path, qc); save_cache(manifest, inputs, [table, qc_path])
    return inventory, qc, inputs


def map_hits(args, hit_reps, genomes, uhgg_raw):
    outdir = Path(args.outdir)
    aliases = build_genome_alias_map(uhgg_raw)
    selected = set(genomes["genome_id"].astype(str))
    if hit_reps:
        tmp = ensure_dir(outdir / "_tmp_membership")
        pre = grep_membership_lines(args.uhgp_membership, hit_reps, tmp)
        mapping, seen, unresolved = parse_membership_for_hits(pre if pre is not None else args.uhgp_membership, hit_reps, selected, aliases)
        fraction = len(seen) / len(hit_reps)
        if pre is not None:
            pre.unlink(missing_ok=True)
        (tmp / "hit_rep_ids.txt").unlink(missing_ok=True)
    else:
        mapping = pd.DataFrame(columns=["protein_rep_id", "genome_id"])
        seen = set(); unresolved = set(); fraction = None
    qc = {"hit_representatives_n": len(hit_reps), "hit_representatives_seen_n": len(seen),
          "hit_representatives_seen_fraction": fraction,
          "missing_hit_representatives_n": len(hit_reps - seen),
          "missing_hit_representative_examples": sorted(hit_reps - seen)[:20],
          "selected_genomes_with_any_search_hit_n": int(mapping["genome_id"].nunique()),
          "unresolved_recognized_alias_n_on_hit_lines": len(unresolved), "unresolved_alias_examples": sorted(unresolved)[:50]}
    write_json(outdir / "membership_mapping_qc.json", qc)
    if hit_reps - seen:
        die(f"Only {len(seen):,}/{len(hit_reps):,} search-hit representatives occur in membership. "
            f"Examples: {', '.join(qc['missing_hit_representative_examples'])}. "
            "See membership_mapping_qc.json; resolve the mismatch before reporting coverage.")
    write_table(mapping, outdir / "uhgp_hit_rep_to_selected_genome.tsv.gz")
    return mapping


def run_search_and_breadth(args):
    outdir = ensure_dir(args.outdir); peptides = validate_prepared(args)
    selected = pd.read_csv(outdir / "selected_species.tsv", sep="\t")
    genomes = pd.read_csv(outdir / "selected_uhgg_genomes.tsv.gz", sep="\t")
    uhgg_raw, _ = read_uhgg_metadata(args.uhgg_metadata)
    prep = read_json(outdir / "prepare_manifest.json")
    if source_signature(args.uhgg_metadata) != prep["inputs"]["sources"]["uhgg_metadata"]:
        die("UHGG metadata differs from preparation input; rerun prepare.")
    (outdir / "search_manifest.json").unlink(missing_ok=True)
    inventory, inventory_qc, inventory_inputs = inventory_membership(args, uhgg_raw)
    genomes = genomes.merge(inventory, on="genome_id", how="left", validate="one_to_one")
    if genomes["in_uhgp_membership"].isna().any():
        die("Inventory does not cover all selected genome IDs.")
    genomes["evaluable_hq"] = genomes["high_quality"] & genomes["in_uhgp_membership"]
    hits_path = outdir / "uhgp_candidate_hits.tsv.gz"
    scan_inputs = {"peptides": peptides.to_dict("records"), "protein_fasta": source_signature(args.uhgp_faa, large=True),
                   "max_substitutions": MAX_SUBSTITUTIONS, "conservative_grantham_max": args.conservative_grantham_max,
                   "algorithm": CACHE_SCHEMA}
    scan_manifest = outdir / "protein_search_manifest.json"
    if args.reuse and cache_matches(scan_manifest, scan_inputs):
        log("Reusing compatible complete protein search.")
        hits = pd.read_csv(hits_path, sep="\t", dtype={"peptide_id": str, "protein_rep_id": str})
        hit_reps = set(hits["protein_rep_id"])
    else:
        scan_manifest.unlink(missing_ok=True)
        tmp = outdir / "uhgp_candidate_hits.partial.tsv"
        threads = args.threads or int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
        hit_reps = scan_uhgp(peptides, args.uhgp_faa, tmp, threads, args.batch_sequences, args.progress_every, args.conservative_grantham_max)
        hits = pd.read_csv(tmp, sep="\t", dtype={"peptide_id": str, "protein_rep_id": str})
        write_table(hits, hits_path); tmp.unlink(missing_ok=True)
        save_cache(scan_manifest, scan_inputs, [hits_path, outdir / "uhgp_search_qc.json"])
    map_inputs = {"inventory_inputs": inventory_inputs, "selected_genomes": digest_json(sorted(genomes["genome_id"].tolist())),
                  "hit_representatives": digest_json(sorted(hit_reps)), "algorithm": CACHE_SCHEMA}
    map_manifest = outdir / "hit_mapping_manifest.json"; map_path = outdir / "uhgp_hit_rep_to_selected_genome.tsv.gz"
    if args.reuse and cache_matches(map_manifest, map_inputs):
        log("Reusing compatible representative-to-selected-genome mapping.")
        mapping = pd.read_csv(map_path, sep="\t")
    else:
        map_manifest.unlink(missing_ok=True)
        mapping = map_hits(args, hit_reps, genomes, uhgg_raw)
        save_cache(map_manifest, map_inputs, [map_path, outdir / "membership_mapping_qc.json"])
    tiers = [t for t, _ in TIERS]
    # Collapse redundant windows BEFORE expanding representative proteins to genomes.
    protein_presence = hits.groupby(["peptide_id", "protein_rep_id"], sort=False)[tiers].max().reset_index()
    protein_presence = protein_presence[protein_presence[tiers].any(axis=1)]
    x = protein_presence.merge(mapping, on="protein_rep_id", how="inner").merge(
        genomes[["genome_id", "species_key", "high_quality", "evaluable_hq"]], on="genome_id", how="inner", validate="many_to_one")
    groupcols = ["peptide_id", "genome_id", "species_key", "high_quality", "evaluable_hq"]
    if len(x):
        gp = x.groupby(groupcols, sort=False)[tiers].max().reset_index()
    else:
        gp = pd.DataFrame(columns=groupcols + tiers)
    write_table(gp, outdir / "peptide_genome_tier_presence.tsv.gz")
    repcounts = genomes.groupby("species_key").agg(n_genomes_in_uhgp=("in_uhgp_membership", "sum"), n_hq_genomes_in_uhgp=("evaluable_hq", "sum")).reset_index()
    species_qc = selected.merge(repcounts, on="species_key", how="left", validate="one_to_one")
    for c in ["n_genomes_in_uhgp", "n_hq_genomes_in_uhgp"]:
        species_qc[c] = species_qc[c].fillna(0).astype(int)
    species_qc["hq_reference_representation_fraction"] = species_qc["n_hq_genomes_in_uhgp"] / species_qc["n_uhgg_genomes_hq"].replace(0, np.nan)
    species_qc["hq_evaluable"] = species_qc["n_hq_genomes_in_uhgp"].gt(0)
    species_qc["low_hq_reference_n"] = species_qc["n_hq_genomes_in_uhgp"].between(1, 4)
    count_all = gp.groupby(["peptide_id", "species_key"])[tiers].sum() if len(gp) else pd.DataFrame(columns=tiers)
    count_hq = gp[gp["evaluable_hq"].eq(True)].groupby(["peptide_id", "species_key"])[tiers].sum() if len(gp) else pd.DataFrame(columns=tiers)
    breadth = peptides.merge(species_qc, how="cross")
    idx = pd.MultiIndex.from_frame(breadth[["peptide_id", "species_key"]])
    for tier in tiers:
        ka = count_all[tier].reindex(idx, fill_value=0).to_numpy(dtype=int) if len(count_all) else np.zeros(len(breadth), dtype=int)
        kh = count_hq[tier].reindex(idx, fill_value=0).to_numpy(dtype=int) if len(count_hq) else np.zeros(len(breadth), dtype=int)
        breadth[f"hit_genomes_all_{tier}"] = ka; breadth[f"hit_genomes_hq_{tier}"] = kh
        breadth[f"breadth_all_{tier}"] = ka / breadth["n_genomes_in_uhgp"].replace(0, np.nan)
        # Missing catalogue genomes are excluded from the evaluable denominator,
        # and separately enter the missing-reference sensitivity interval.
        breadth[f"breadth_hq_{tier}"] = kh / breadth["n_hq_genomes_in_uhgp"].replace(0, np.nan)
        breadth[f"breadth_hq_all_references_{tier}"] = kh / breadth["n_uhgg_genomes_hq"].replace(0, np.nan)
        breadth[f"carrier_broad_hq_{tier}"] = breadth["n_hq_genomes_in_uhgp"].ge(args.min_hq_genomes_for_broad) & breadth[f"breadth_hq_{tier}"].ge(args.broad_breadth)
    examples = hits.sort_values(["peptide_id", "n_substitutions", "max_grantham", "protein_rep_id"]).groupby("peptide_id", group_keys=False).head(args.sanity_hits_per_peptide)
    outputs = [outdir / name for name in ["selected_genome_membership_qc.tsv.gz", "species_evaluability.tsv", "peptide_species_breadth.tsv.gz", "uhgp_hit_examples.tsv", "peptide_genome_tier_presence.tsv.gz", "search_breadth_qc.json"]]
    for df, path in zip([genomes, species_qc, breadth, examples], outputs[:4]):
        write_table(df, path)
    qc = {"script_version": VERSION, "conservative_grantham_max": args.conservative_grantham_max,
          "max_substitutions_searched": MAX_SUBSTITUTIONS, "broad_breadth": args.broad_breadth,
          "min_hq_genomes_for_broad": args.min_hq_genomes_for_broad,
          "headline_breadth_denominator": "HQ genomes independently shown to have at least one protein in UHGP membership.",
          "cluster_expansion_limitation": "Genome carriage is inferred from matching representatives. Partial-length clustering can misassign peptide presence; direction and magnitude of net error are not established. Genome presence in membership does not establish locus recovery.",
          "selected_hq_genomes": int(genomes["high_quality"].sum()), "selected_evaluable_hq_genomes": int(genomes["evaluable_hq"].sum())}
    write_json(outputs[-1], qc)
    save_cache(outdir / "search_manifest.json", {"preparation_inputs": prep["inputs"], "scan_inputs": scan_inputs, "mapping_inputs": map_inputs,
               "broad_breadth": args.broad_breadth, "min_hq_genomes_for_broad": args.min_hq_genomes_for_broad}, outputs)
    log("Search, membership inventory and species breadth complete.")


def flatten_dict(obj, prefix=""):
    rows = []
    for key, value in obj.items():
        label = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            rows.extend(flatten_dict(value, label))
        else:
            rows.append((label, json.dumps(value, default=str) if isinstance(value, list) else value))
    return rows


def summarize_coverage(args):
    outdir = Path(args.outdir); peptides = validate_prepared(args)
    search_manifest = outdir / "search_manifest.json"
    if not search_manifest.exists():
        die("No completed v3 search. Run search or all first.")
    manifest = read_json(search_manifest)
    if not cache_matches(search_manifest, manifest["inputs"]):
        die("Search/breadth outputs changed or are incomplete; rerun search.")
    if manifest["inputs"]["preparation_inputs"] != read_json(outdir / "prepare_manifest.json")["inputs"]:
        die("Preparation changed since search; rerun search with --reuse to rebuild affected stages.")
    search_qc = read_json(outdir / "search_breadth_qc.json")
    if search_qc["selected_evaluable_hq_genomes"] == 0:
        log("WARNING: no HQ genomes are evaluable; coverage has no observed support and the assumption interval is [0,1].")
    for name in ["broad_breadth", "conservative_grantham_max"]:
        supplied = getattr(args, name, None)
        if supplied is not None and supplied != search_qc[name]:
            die(f"--{name.replace('_','-')} differs from search provenance; rerun search to change it.")
    breadth = pd.read_csv(outdir / "peptide_species_breadth.tsv.gz", sep="\t", dtype={"peptide_id": str})
    species = pd.read_csv(outdir / "species_evaluability.tsv", sep="\t").set_index("species_key")
    mat = pd.read_csv(outdir / "selected_species_abundance_matrix.tsv.gz", sep="\t",
                      index_col="sample_name", dtype={"sample_name": str})
    meta = pd.read_csv(outdir / "included_sample_metadata.tsv", sep="\t", dtype={"sample_name": str, "bioproject_id": str}).set_index("sample_name")
    sample_qc = pd.read_csv(outdir / "preparation_sample_qc.tsv.gz", sep="\t",
                            index_col="sample_name", dtype={"sample_name": str})
    if not mat.index.is_unique or not meta.index.is_unique or not sample_qc.index.is_unique:
        die("Duplicate sample names in prepared matrix or sample metadata; rerun prepare.")
    if not mat.index.equals(meta.index) or set(mat.index) != set(sample_qc.index):
        die("Sample identifiers differ across prepared abundance, metadata and QC files; rerun prepare.")
    sample_qc = sample_qc.reindex(mat.index)
    species = species.reindex(mat.columns)
    if species["n_hq_genomes_in_uhgp"].isna().any():
        die("Species QC does not cover the selected abundance matrix.")
    eval_mask = species["n_hq_genomes_in_uhgp"].gt(0)
    r = species["hq_reference_representation_fraction"].fillna(0)
    sample_qc["panel_hq_species_absent_from_uhgp_abundance"] = mat.loc[:, species["n_uhgg_genomes_hq"].gt(0) & ~eval_mask].sum(axis=1)
    sample_qc["evaluable_species_abundance"] = mat.loc[:, eval_mask].sum(axis=1)
    sample_qc["unevaluable_species_abundance"] = 1 - sample_qc["evaluable_species_abundance"]
    sample_qc["represented_hq_reference_weight"] = mat.mul(r, axis=1).sum(axis=1)
    sample_qc["unresolved_reference_weight"] = 1 - sample_qc["represented_hq_reference_weight"]
    sample_qc["partly_missing_hq_reference_weight"] = sample_qc["evaluable_species_abundance"] - sample_qc["represented_hq_reference_weight"]
    sample_qc["low_n_hq_species_abundance"] = mat.loc[:, species["n_hq_genomes_in_uhgp"].between(1,4)].sum(axis=1)
    sample_qc["bioproject_id"] = meta["bioproject_id"].reindex(mat.index).fillna("")
    sample_qc["labelled_bacterial_fraction_of_input_prokaryotic_abundance"] = meta["labelled_bacterial_fraction_of_input_prokaryotic_abundance"].reindex(mat.index)
    for col in sample_qc.select_dtypes(include="number"):
        if col.endswith(("abundance", "weight")):
            if ((sample_qc[col] < -1e-8) | (sample_qc[col] > 1+1e-8)).any():
                die(f"Invalid fraction in sample QC: {col}")
            sample_qc[col] = sample_qc[col].clip(0,1)
    # Undefined, rather than zero, when no species can be evaluated in a sample.
    evaluable_denominator = sample_qc["evaluable_species_abundance"].where(
        sample_qc["evaluable_species_abundance"].gt(0))
    have_evaluable_samples = evaluable_denominator.notna().any()
    ext_rows = []; sample_parts = []; project_parts = []; threshold_rows = []
    for p in peptides.itertuples(index=False):
        b = breadth[breadth["peptide_id"].eq(p.peptide_id)].set_index("species_key").reindex(mat.columns)
        if b["sequence"].isna().any():
            die(f"Incomplete species breadth grid for {p.peptide_id}.")
        er = {"peptide_id": p.peptide_id, "sequence": p.sequence, "length": p.length}
        qr = dict(er)
        for tier, _ in TIERS:
            w = b[f"breadth_hq_{tier}"].fillna(0)
            if ((w < 0) | (w > 1)).any():
                die("Breadth outside [0,1].")
            score = mat.mul(w, axis=1).sum(axis=1)
            equal_rate = score.div(evaluable_denominator).clip(0,1)
            lower = mat.mul(w*r, axis=1).sum(axis=1)
            upper = (lower + sample_qc["unresolved_reference_weight"]).clip(0,1)
            species_upper = (score + sample_qc["unevaluable_species_abundance"]).clip(0,1)
            broad = mat.loc[:, b[f"carrier_broad_hq_{tier}"].eq(True)].sum(axis=1)
            er[f"median_coverage_{tier}"] = float(score.median())
            er[f"median_equal_rate_extrapolated_coverage_{tier}"] = float(equal_rate.median()) if have_evaluable_samples else np.nan
            er[f"p05_coverage_{tier}"] = float(score.quantile(.05))
            er[f"median_missing_reference_lower_{tier}"] = float(lower.median())
            er[f"median_missing_reference_upper_{tier}"] = float(upper.median())
            qr[f"median_thresholded_coverage_{tier}"] = float(broad.median())
            qr[f"n_broad_carrier_species_{tier}"] = int(b[f"carrier_broad_hq_{tier}"].eq(True).sum())
            part = pd.DataFrame({"peptide_id": p.peptide_id, "sample_name": mat.index, "tier": tier,
                "breadth_weighted_coverage": score.to_numpy(), "equal_rate_extrapolated_coverage": equal_rate.to_numpy(),
                "missing_reference_lower": lower.to_numpy(),
                "missing_reference_upper": upper.to_numpy(), "unevaluated_species_upper": species_upper.to_numpy(),
                "broad_hq_thresholded_coverage": broad.to_numpy(), "bioproject_id": sample_qc["bioproject_id"].to_numpy()})
            sample_parts.append(part)
            known = part[part["bioproject_id"].ne("")]
            if len(known):
                projects = known.groupby("bioproject_id").agg(n_samples=("sample_name","size"), median_coverage=("breadth_weighted_coverage","median"),
                    median_missing_reference_lower=("missing_reference_lower","median"), median_missing_reference_upper=("missing_reference_upper","median")).reset_index()
                projects.insert(0, "tier", tier); projects.insert(0, "peptide_id", p.peptide_id)
                project_parts.append(projects)
                er[f"median_of_bioproject_medians_{tier}"] = float(projects["median_coverage"].median())
            else:
                er[f"median_of_bioproject_medians_{tier}"] = np.nan
        for field in ["panel_abundance", "evaluable_species_abundance", "unevaluable_species_abundance",
                      "represented_hq_reference_weight", "unresolved_reference_weight"]:
            er[f"median_{field}"] = float(sample_qc[field].median())
        ext_rows.append(er); threshold_rows.append(qr)
    extended = pd.DataFrame(ext_rows)
    concise_cols = ["peptide_id", "sequence"]
    for tier in ["exact", "cons_le1", "cons_le2"]:
        concise_cols += [f"median_missing_reference_lower_{tier}", f"median_coverage_{tier}",
                         f"median_equal_rate_extrapolated_coverage_{tier}", f"median_missing_reference_upper_{tier}"]
    concise_cols += ["median_panel_abundance", "median_evaluable_species_abundance",
                     "median_unevaluable_species_abundance", "median_represented_hq_reference_weight"]
    concise = extended[concise_cols].rename(columns={f"{prefix}_{tier}": f"{prefix}_{name}"
        for prefix in ["median_coverage", "median_equal_rate_extrapolated_coverage"]
        for tier,name in [("cons_le1","le1_conservative"),("cons_le2","le2_conservative")]})
    sample_df = pd.concat(sample_parts, ignore_index=True)
    projects = pd.concat(project_parts, ignore_index=True) if project_parts else pd.DataFrame(columns=["peptide_id","tier","bioproject_id","n_samples","median_coverage","median_missing_reference_lower","median_missing_reference_upper"])
    threshold = pd.DataFrame(threshold_rows)
    qc_summary = pd.DataFrame([{"metric": col, **coverage_stats(sample_qc[col])} for col in sample_qc.select_dtypes(include="number")])
    outputs = [(concise,"concise_summary.tsv",False),(extended,"extended_summary.tsv",False),(sample_df,"peptide_sample_coverage.tsv.gz",False),
               (sample_qc,"sample_evaluability.tsv.gz",True),(qc_summary,"evaluability_summary.tsv",False),(projects,"bioproject_coverage.tsv",False),
               (threshold,"broad_hq_threshold_qc.tsv",False)]
    for df,name,index in outputs:
        write_table(df, outdir/name, index)
    run_info = {"script_version": VERSION, "prepared": read_json(outdir / "preparation_diagnostics.json"), "search": search_qc,
        "protein_search_qc": read_json(outdir / "uhgp_search_qc.json"), "membership_inventory": read_json(outdir / "membership_inventory_qc.json"),
        "hit_mapping": read_json(outdir / "membership_mapping_qc.json"),
        "formula": "C_j=sum_s a_js*(k_s/n_evaluableHQ_s); unevaluable species contribute no observed support.",
        "equal_rate_extrapolation": "For E_j=sum_s a_js over evaluable species, T_j=C_j/E_j if E_j>0, otherwise undefined. T_j estimates the total species-labelled bacterial fraction only under the assumption that the unevaluable fraction has the same abundance-weighted peptide match rate as the evaluable fraction. Per-sample ratios precede medians; this assumption is not tested by the missing species.",
        "samples_zero_evaluable_n": int(sample_qc["evaluable_species_abundance"].eq(0).sum()),
        "missing_reference_sensitivity": "L_j=sum_s a_js*k_s/n_allHQ_s; U_j=L_j+1-sum_s a_js*n_evaluableHQ_s/n_allHQ_s. Zero-reference species use ratio 0. These endpoints vary only missing-species/reference assumptions; not confidence intervals or biological bounds.",
        "partial_reference_note": "C extrapolates breadth observed in represented HQ genomes to their species' abundance. The sensitivity range instead allows missing genomes to be all negative or all positive. Catalogue genomes are not population-random samples.",
        "bioproject_note": "Pooled medians are primary descriptive summaries. Equal-project median-of-medians is a sensitivity view; missing project IDs are excluded from that view. No independence or population-representativeness inference is made.",
        "provenance_note": "Small inputs have full SHA-256. Giant FASTA/membership identity uses path, size, nanosecond mtime and sampled SHA-256, not a full-file cryptographic checksum. Completed-cache artifacts have size/mtime guards. Keep the external catalogue release/checksums with the run.",
        "python_version": platform.python_version(), "pandas_version": pd.__version__, "numpy_version": np.__version__,
        "pyahocorasick_available": ahocorasick is not None,
        "excel_data_row_limit": args.excel_max_rows}
    write_json(outdir / "run_diagnostics.json", run_info)
    tables = {"Concise summary": concise, "Extended summary": extended, "Evaluability": qc_summary,
              "Species breadth": breadth, "Selected species QC": species.rename_axis("species_key").reset_index(), "Sample evaluability": sample_qc.reset_index(),
              "Per-sample coverage": sample_df, "BioProject coverage": projects, "Broad carrier QC": threshold,
              "UHGP hit examples": pd.read_csv(outdir / "uhgp_hit_examples.tsv", sep="\t"),
              "Run & QC": pd.DataFrame(flatten_dict(run_info), columns=["parameter","value"])}
    write_excel(outdir, tables, run_info, args.excel_max_rows)
    log("Summary complete. Start with concise_summary.tsv and evaluability_summary.tsv.")


def write_excel(outdir, tables, run_info, max_rows):
    try:
        import openpyxl
        from openpyxl.styles import Font, PatternFill, Alignment
    except ImportError:
        log("WARNING: openpyxl unavailable; all TSV/JSON results are complete.")
        return
    readme = pd.DataFrame([
        ("Estimate and support", "Equal-rate extrapolated coverage = reference-supported community share / evaluable species abundance, formed per sample before the median. It estimates coverage of the species-labelled bacterial fraction if unevaluable species have the same abundance-weighted match rate. Present with supported share, evaluability and the extreme sensitivity range. Undefined when no species are evaluable."),
        ("Supported and unknown", "Median coverage_* is the reference-supported share of the full species-labelled bacterial denominator; median_evaluable_species_abundance and median_unevaluable_species_abundance show how much was assessed and remained unknown. These separate medians must not be multiplied or subtracted to reconstruct a sample."),
        ("Target", "Panel selected BEFORE quality checks to reach 95% mean abundance by default; not a per-sample guarantee. See Evaluability for the actual target/result."),
        ("Missing-reference range", "Extreme scenarios: unresolved species and HQ references are all negative or all positive. Endpoint medians are formed after the per-sample calculations; they are not confidence intervals or guaranteed biological bounds."),
        ("Reference weight", "Abundance times the fraction of HQ catalogue genomes represented in UHGP; this is a reference-based diagnostic, not measured community missingness."),
        ("Broad carrier QC", f"Secondary binary view using breadth >= {run_info['search']['broad_breadth']:.0%} and at least {run_info['search']['min_hq_genomes_for_broad']} represented HQ genomes."),
        ("Similarity", f"Conservative = Grantham <= {run_info['search']['conservative_grantham_max']}; sequence categories are not antibody-binding probabilities."),
        ("Carriage", "Exact refers to the representative protein window; genome carriage is inferred through cluster membership and does not prove locus recovery, protein expression or accessibility."),
        ("Large tables", "Tables exceeding the configured Excel row budget are omitted in full, with an explicit TSV pointer. TSV files always contain all rows."),
        ("Sources", "See README and the provenance manifests for releases, source identity, reuse rules and limitations.")], columns=["item","explanation"])
    sources = {"Species breadth":"peptide_species_breadth.tsv.gz", "Per-sample coverage":"peptide_sample_coverage.tsv.gz",
               "Sample evaluability":"sample_evaluability.tsv.gz", "BioProject coverage":"bioproject_coverage.tsv",
               "UHGP hit examples":"uhgp_hit_examples.tsv", "Selected species QC":"species_evaluability.tsv",
               "Concise summary":"concise_summary.tsv", "Extended summary":"extended_summary.tsv", "Evaluability":"evaluability_summary.tsv",
               "Broad carrier QC":"broad_hq_threshold_qc.tsv", "Run & QC":"run_diagnostics.json"}
    path = outdir / "UHGG_peptide_coverage_summary.xlsx"
    tmp = outdir / "UHGG_peptide_coverage_summary.partial.xlsx"
    # Write-only mode limits workbook memory. Styling is confined to header and
    # numeric cells while streaming; no second full worksheet traversal.
    wb = openpyxl.Workbook(write_only=True)
    from openpyxl.cell import WriteOnlyCell
    def sheet(name, df):
        ws = wb.create_sheet(name); ws.freeze_panes = "A2"
        ws.auto_filter.ref = f"A1:{openpyxl.utils.get_column_letter(max(1,len(df.columns)))}{len(df)+1}"
        for i, col in enumerate(df.columns, 1):
            ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = min(48,max(16,len(str(col))+2))
        header = []
        for col in df.columns:
            c = WriteOnlyCell(ws, str(col)); c.font = Font(color="FFFFFF",bold=True); c.fill = PatternFill("solid",fgColor="1F4E78")
            c.alignment = Alignment(wrap_text=True); header.append(c)
        ws.append(header)
        for row in df.itertuples(index=False, name=None):
            values=[]
            for col, value in zip(df.columns,row):
                if pd.isna(value): value=None
                elif isinstance(value,np.generic): value=value.item()
                if isinstance(value,str): value=value[:32760]
                c=WriteOnlyCell(ws,value)
                if isinstance(value,str): c.data_type="s"
                if isinstance(value,float) and any(k in str(col).lower() for k in ["coverage","abundance","breadth","fraction","weight","missing_reference"]): c.number_format="0.0%"
                values.append(c)
            ws.append(values)
    sheet("README", readme)
    omitted=[]
    for name,df in tables.items():
        if len(df) > max_rows or len(df) > 1048575:
            source = sources.get(name, "see TSV outputs")
            omitted.append({"sheet":name,"rows":len(df),"complete_file":source})
            sheet(name,pd.DataFrame([{"status":"Complete table omitted from Excel due to row budget", "rows":len(df), "complete_file":source}]))
        else:
            sheet(name,df)
    wb.save(tmp); tmp.replace(path)
    write_json(outdir / "excel_export_qc.json", {"workbook":path.name,"max_data_rows_per_sheet":max_rows,"omitted_tables":omitted})
    log(f"Wrote {path}; {len(omitted)} large table(s) retained in TSV only.")


def run_all(args):
    prepare_species_and_genomes(args); run_search_and_breadth(args); summarize_coverage(args)


def common_args(p):
    p.add_argument("--peptides", required=True, help="Candidate FASTA; validated against prepared candidates in later stages.")
    p.add_argument("--outdir", default="results_uhgg_v3")
    p.add_argument("--min-peptide-length", type=int, default=12)
    p.add_argument("--max-peptide-length", type=int, default=25)


def preparation_args(p):
    for name in ["relative-abundance","sample-metadata","host-genomes-info","uhgg-metadata"]:
        p.add_argument("--"+name, required=True)
    p.add_argument("--target-community-coverage", type=float, default=.95, help="Target cumulative mean bacterial abundance before QC; not a per-sample target.")
    p.add_argument("--min-completeness", type=float, default=90)
    p.add_argument("--max-contamination", type=float, default=5)
    p.add_argument("--cluster-mapping-min-purity", type=float, default=.90)


def search_args(p, include_metadata=True):
    if include_metadata: p.add_argument("--uhgg-metadata", required=True)
    p.add_argument("--uhgp-faa", required=True); p.add_argument("--uhgp-membership", required=True)
    p.add_argument("--conservative-grantham-max", type=int, default=50)
    p.add_argument("--broad-breadth", type=float, default=.75)
    p.add_argument("--min-hq-genomes-for-broad", type=int, default=1)
    p.add_argument("--threads", type=int, default=0, help="0 uses SLURM_CPUS_PER_TASK, otherwise 1.")
    p.add_argument("--batch-sequences", type=int, default=5000)
    p.add_argument("--progress-every", type=int, default=250000)
    p.add_argument("--sanity-hits-per-peptide", type=int, default=50)
    p.add_argument("--reuse", action="store_true", help="Reuse only compatible completed stage caches; rebuild changed stages automatically. Legacy v2 caches are not trusted.")


def summary_args(p, provenance_options=True):
    if provenance_options:
        p.add_argument("--conservative-grantham-max", type=int, default=None, help="Optional consistency assertion; does not change results.")
        p.add_argument("--broad-breadth", type=float, default=None, help="Optional consistency assertion; does not change results.")
    p.add_argument("--excel-max-rows", type=int, default=EXCEL_DATA_ROW_LIMIT, help="Maximum data rows per Excel sheet; larger full tables remain in TSV.")


def make_parser():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub=ap.add_subparsers(dest="command",required=True)
    p=sub.add_parser("prepare",formatter_class=argparse.ArgumentDefaultsHelpFormatter); common_args(p); preparation_args(p); p.set_defaults(func=prepare_species_and_genomes)
    p=sub.add_parser("search",formatter_class=argparse.ArgumentDefaultsHelpFormatter); common_args(p); search_args(p); p.set_defaults(func=run_search_and_breadth)
    p=sub.add_parser("summarize",formatter_class=argparse.ArgumentDefaultsHelpFormatter); common_args(p); summary_args(p); p.set_defaults(func=summarize_coverage)
    p=sub.add_parser("all",formatter_class=argparse.ArgumentDefaultsHelpFormatter); common_args(p); preparation_args(p); search_args(p,False); summary_args(p,False); p.set_defaults(func=run_all)
    return ap


def validate_args(args):
    if args.min_peptide_length < 12 or args.max_peptide_length < args.min_peptide_length:
        die("Peptide length bounds must satisfy max >= min >= 12; shorter three-substitution searches are impractical on UHGP.")
    for name in ["target_community_coverage","cluster_mapping_min_purity"]:
        if hasattr(args,name) and not 0 < getattr(args,name) <= 1: die(f"{name} must be in (0,1].")
    if getattr(args,"broad_breadth",None) is not None and not 0 < args.broad_breadth <= 1: die("broad_breadth must be in (0,1].")
    if getattr(args,"conservative_grantham_max",None) is not None and not 0 <= args.conservative_grantham_max <= 215: die("Grantham cutoff must be in [0,215].")
    for name in ["min_completeness","max_contamination"]:
        if hasattr(args,name) and not 0 <= getattr(args,name) <= 100: die(f"{name} must be in [0,100].")
    for name in ["batch_sequences","progress_every","sanity_hits_per_peptide","min_hq_genomes_for_broad","excel_max_rows"]:
        if hasattr(args,name) and getattr(args,name) < 1: die(f"{name} must be positive.")
    if hasattr(args,"threads") and args.threads < 0: die("threads cannot be negative.")
    if getattr(args,"excel_max_rows",0) > 1048575: die("Excel allows at most 1,048,575 data rows plus the header.")


def main(argv=None):
    args=make_parser().parse_args(argv); validate_args(args)
    log(f"gut_peptide_coverage_uhgg.py version {VERSION}")
    started=time.monotonic(); args.func(args)
    log(f"Completed {args.command} in {(time.monotonic()-started)/60:.1f} minutes.")


if __name__ == "__main__":
    main()
