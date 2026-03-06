#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-AlphaFold Analysis Pipeline

For each AlphaFold model, computes per-model metrics (pLDDT, RMSD, SASA),
then aggregates by scaffold design to produce a final per-scaffold summary.

Since each scaffold generates N sequences via ProteinMPNN → N AF2 models,
the per-scaffold summary averages metrics across all N models and identifies
the single best model (lowest mean epitope RMSD + highest pLDDT).

Model name format expected (ProteinMPNN + ColabFold):
    scaffoldID__ep1_ep2__T0.2_sample3_...
    scaffoldID__ep1_ep2_relaxed_rank_001_...

Outputs (in output_dir/):
  af2_metrics_summary_TIMESTAMP.csv   → one row per AF2 model
  af2_per_residue_TIMESTAMP.csv       → one row per residue
  af2_epitope_rmsd_TIMESTAMP.csv      → one row per (model, epitope)
  scaffold_summary_TIMESTAMP.csv      → one row per scaffold design ← KEY OUTPUT

Usage:
    python process_af2_outputs.py \\
        --af2_dir      AlphaFold_results/ \\
        --epitopes_dir Epitopes/ \\
        --output_dir   AF2_analysis/
"""

import os
import re
import glob
import csv
import warnings
import argparse
from datetime import datetime
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser, PPBuilder, Superimposer
from Bio.SeqUtils import seq1

warnings.simplefilter("ignore")

import pyrosetta
from pyrosetta import pose_from_pdb
pyrosetta.init("-mute all")


# ══════════════════════════════════════════════════════════════════
# EPITOPE DICTIONARY
# ══════════════════════════════════════════════════════════════════

def extract_sequence_from_pdb(pdb_path):
    parser    = PDBParser(QUIET=True)
    name      = os.path.splitext(os.path.basename(pdb_path))[0]
    structure = parser.get_structure(name, pdb_path)
    model     = structure[0]
    ppb       = PPBuilder()
    sequences = [str(pp.get_sequence()) for pp in ppb.build_peptides(model)]
    if not sequences:
        for chain in model:
            residues = [r for r in chain.get_residues() if r.get_id()[0] == " "]
            seq = seq1("".join(r.resname for r in residues))
            if seq:
                sequences.append(seq)
            break
    return "".join(sequences) if sequences else None


def build_epitope_dictionary(epitopes_dir):
    pdb_files = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files in: {epitopes_dir}")
    epitope_dict = {}
    print(f"\nBuilding epitope dictionary from {len(pdb_files)} PDB files...")
    for pdb_path in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        seq  = extract_sequence_from_pdb(pdb_path)
        if seq:
            epitope_dict[name] = seq
            print(f"  [{name}]  {seq}  ({len(seq)} aa)")
        else:
            print(f"  WARNING: Could not read sequence from {pdb_path}")
    print(f"  → {len(epitope_dict)} epitopes loaded.\n")
    return epitope_dict


def find_epitope_spans(full_sequence, epitope_dict):
    spans = {}
    for name, seq in epitope_dict.items():
        idx = full_sequence.find(seq)
        if idx != -1:
            spans[name] = (idx, idx + len(seq))
    return spans


# ══════════════════════════════════════════════════════════════════
# SCAFFOLD ID PARSING
# ══════════════════════════════════════════════════════════════════

def parse_scaffold_id(model_name):
    """
    Strips ProteinMPNN/ColabFold suffixes to recover the scaffold design ID.

    Examples:
      relaxed_88__frag1_frag3__T0.2_sample1_score0.18  → relaxed_88__frag1_frag3
      relaxed_88__frag1_frag3_relaxed_rank_001_alpha... → relaxed_88__frag1_frag3
      relaxed_88__frag1_frag3                           → relaxed_88__frag1_frag3
    """
    # ColabFold pattern: _relaxed_rank_ or _unrelaxed_rank_
    m = re.match(r"^(.+?)(?:_relaxed_rank_|_unrelaxed_rank_)", model_name)
    if m:
        return m.group(1)
    # ProteinMPNN pattern: __T0.x or __sample
    m = re.match(r"^(.+?)(?:__T[\d\.]+|__sample\d)", model_name)
    if m:
        return m.group(1)
    return model_name


# ══════════════════════════════════════════════════════════════════
# 1. pLDDT
# ══════════════════════════════════════════════════════════════════

def load_af2_model(pdb_path):
    parser = PDBParser(QUIET=True)
    name   = os.path.splitext(os.path.basename(pdb_path))[0]
    struct = parser.get_structure(name, pdb_path)
    return struct, struct[0]


def get_residue_plddt(biopdb_model):
    residue_data = []
    for chain in biopdb_model:
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue
            if "CA" not in residue:
                continue
            ca      = residue["CA"]
            resname = residue.resname
            aa      = seq1(resname) if resname != "UNK" else "X"
            residue_data.append({
                "chain":       chain.get_id(),
                "res_id":      residue.get_id()[1],
                "resname":     resname,
                "aa":          aa,
                "plddt":       ca.get_bfactor(),
                "residue_obj": residue,
            })
    return residue_data


def get_full_sequence(residue_data):
    return "".join(r["aa"] for r in residue_data)


def compute_plddt_metrics(residue_data, epitope_spans):
    all_plddt = [r["plddt"] for r in residue_data]
    global_mean   = sum(all_plddt) / len(all_plddt) if all_plddt else 0.0
    global_median = float(np.median(all_plddt))      if all_plddt else 0.0

    epitope_mask = [False] * len(residue_data)
    for name, (start, end) in epitope_spans.items():
        for i in range(start, min(end, len(epitope_mask))):
            epitope_mask[i] = True

    ep_plddt   = [r["plddt"] for i, r in enumerate(residue_data) if epitope_mask[i]]
    scaf_plddt = [r["plddt"] for i, r in enumerate(residue_data) if not epitope_mask[i]]

    per_epitope = {}
    for name, (start, end) in epitope_spans.items():
        vals = [residue_data[i]["plddt"] for i in range(start, end)
                if i < len(residue_data)]
        per_epitope[name] = round(sum(vals) / len(vals), 3) if vals else None

    return {
        "global_mean_plddt":   round(global_mean,   3),
        "global_median_plddt": round(global_median, 3),
        "mean_epitope_plddt":  round(sum(ep_plddt)   / len(ep_plddt),   3) if ep_plddt   else None,
        "mean_scaffold_plddt": round(sum(scaf_plddt) / len(scaf_plddt), 3) if scaf_plddt else None,
        "n_epitope_residues":  len(ep_plddt),
        "n_scaffold_residues": len(scaf_plddt),
        "per_epitope_plddt":   per_epitope,
        "epitope_mask":        epitope_mask,
    }


# ══════════════════════════════════════════════════════════════════
# 2. RMSD
# ══════════════════════════════════════════════════════════════════

def compute_epitope_rmsd(af2_residue_data, epitope_spans, epitopes_dir):
    parser  = PDBParser(QUIET=True)
    results = []

    for epitope_name, (start, end) in epitope_spans.items():
        epitope_pdb = os.path.join(epitopes_dir, f"{epitope_name}.pdb")
        if not os.path.exists(epitope_pdb):
            results.append({"epitope_name": epitope_name, "n_ca_atoms": 0,
                             "rmsd": None, "status": "ORIGINAL_PDB_NOT_FOUND"})
            continue
        try:
            af2_res  = af2_residue_data[start:end]
            orig_str = parser.get_structure(epitope_name, epitope_pdb)
            orig_res = [r for chain in orig_str[0] for r in chain
                        if r.get_id()[0] == " " and "CA" in r]

            n = min(len(af2_res), len(orig_res))
            if n == 0:
                raise ValueError("No CA atoms to align")

            af2_ca  = [r["residue_obj"]["CA"] for r in af2_res[:n]
                       if "CA" in r["residue_obj"]]
            orig_ca = [r["CA"] for r in orig_res[:n]]
            n       = min(len(af2_ca), len(orig_ca))

            if n < 3:
                raise ValueError(f"Too few CA atoms ({n})")

            sup = Superimposer()
            sup.set_atoms(orig_ca[:n], af2_ca[:n])
            results.append({"epitope_name": epitope_name, "n_ca_atoms": n,
                             "rmsd": round(sup.rms, 4), "status": "OK"})
        except Exception as e:
            results.append({"epitope_name": epitope_name, "n_ca_atoms": 0,
                             "rmsd": None, "status": f"ERROR: {e}"})
    return results


# ══════════════════════════════════════════════════════════════════
# 3. SASA
# ══════════════════════════════════════════════════════════════════

def compute_sasa_pyrosetta(pdb_path, epitope_spans, full_sequence):
    pose       = pose_from_pdb(pdb_path)
    calc       = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
    total_sasa = calc.calculate(pose)
    rsd_sasa   = calc.get_residue_sasa()

    epitope_mask     = [False] * len(full_sequence)
    epitope_name_map = {}
    for name, (start, end) in epitope_spans.items():
        for i in range(start, min(end, len(full_sequence))):
            epitope_mask[i]     = True
            epitope_name_map[i] = name

    pdb_info    = pose.pdb_info()
    per_residue = []
    for i in range(1, pose.total_residue() + 1):
        idx    = i - 1
        region = "EPITOPE" if (idx < len(epitope_mask) and epitope_mask[idx]) \
                 else "SCAFFOLD"
        per_residue.append({
            "rosetta_index": i,
            "pdb_resnum":    pdb_info.number(i) if pdb_info else i,
            "chain":         pdb_info.chain(i)  if pdb_info else "A",
            "resname":       pose.residue(i).name3(),
            "aa_1letter":    pose.residue(i).name1(),
            "region":        region,
            "epitope_name":  epitope_name_map.get(idx, ""),
            "sasa":          round(rsd_sasa[i], 3),
        })
    return round(total_sasa, 3), per_residue


# ══════════════════════════════════════════════════════════════════
# 4. CSV HELPERS
# ══════════════════════════════════════════════════════════════════

SUMMARY_COLS = [
    "model_name", "scaffold_id", "n_residues",
    "global_mean_plddt", "global_median_plddt",
    "mean_epitope_plddt", "mean_scaffold_plddt",
    "n_epitope_residues", "n_scaffold_residues",
    "epitopes_found", "n_epitopes_found",
    "mean_epitope_rmsd", "per_epitope_rmsd",
    "total_sasa", "mean_epitope_sasa", "mean_scaffold_sasa",
]

PER_RESIDUE_COLS = [
    "model_name", "rosetta_index", "pdb_resnum", "chain",
    "resname", "aa_1letter", "region", "epitope_name", "plddt", "sasa",
]

RMSD_COLS = ["model_name", "epitope_name", "n_ca_atoms", "rmsd", "status"]

SCAFFOLD_SUMMARY_COLS = [
    # Identity
    "scaffold_id",
    "n_models",
    "epitopes_found",
    "n_epitopes",
    # pLDDT across models
    "mean_global_plddt",
    "std_global_plddt",
    "best_global_plddt",
    "mean_epitope_plddt",
    "mean_scaffold_plddt",
    # RMSD across models
    "mean_epitope_rmsd",
    "std_epitope_rmsd",
    "best_epitope_rmsd",
    "per_epitope_mean_rmsd",   # frag1=0.82;frag3=1.14
    # SASA across models
    "mean_epitope_sasa",
    "std_epitope_sasa",
    "mean_total_sasa",
    # Best single model (lowest RMSD + highest pLDDT)
    "best_model_name",
    "best_model_plddt",
    "best_model_rmsd",
    "best_model_epitope_sasa",
]


def init_csv(path, cols):
    with open(path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=cols).writeheader()


def append_csv(path, rows, cols):
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})


def write_csv(path, rows, cols):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})


# ══════════════════════════════════════════════════════════════════
# 5. PER-SCAFFOLD AGGREGATION
# ══════════════════════════════════════════════════════════════════

def aggregate_by_scaffold(model_rows):
    """
    Groups per-model summary rows by scaffold_id and computes:
      - mean / std / best for pLDDT, RMSD, SASA across all N models
      - per-epitope mean RMSD
      - best single model identified by composite score
        (pLDDT_norm + RMSD_norm_inverted + SASA_norm)

    Parameters
    ----------
    model_rows : list of dicts, one per AF2 model, as stored in SUMMARY_COLS

    Returns a list of per-scaffold dicts sorted by mean_epitope_rmsd asc.
    """
    def safe_float(v):
        try:
            return float(v) if v not in (None, "", "None") else None
        except (ValueError, TypeError):
            return None

    def mean(vals):
        v = [x for x in vals if x is not None]
        return sum(v) / len(v) if v else None

    def std(vals):
        v = [x for x in vals if x is not None]
        return float(np.std(v)) if len(v) > 1 else 0.0

    def minmax_norm(vals, invert=False):
        """Normalise list to [0,1]. invert=True means lower raw = higher score."""
        clean = [v for v in vals if v is not None]
        if not clean:
            return [None] * len(vals)
        lo, hi = min(clean), max(clean)
        if hi == lo:
            return [1.0 if v is not None else None for v in vals]
        result = []
        for v in vals:
            if v is None:
                result.append(None)
            elif invert:
                result.append((hi - v) / (hi - lo))
            else:
                result.append((v - lo) / (hi - lo))
        return result

    # Group by scaffold
    by_scaffold = defaultdict(list)
    for row in model_rows:
        sid = parse_scaffold_id(row["model_name"])
        by_scaffold[sid].append(row)

    scaffolds = []

    for scaffold_id, models in by_scaffold.items():

        # ── pLDDT ────────────────────────────────────────────────
        g_plddts  = [safe_float(m["global_mean_plddt"])   for m in models]
        ep_plddts = [safe_float(m["mean_epitope_plddt"])  for m in models]
        sc_plddts = [safe_float(m["mean_scaffold_plddt"]) for m in models]

        # ── RMSD ─────────────────────────────────────────────────
        # Each model row stores mean_epitope_rmsd and per_epitope_rmsd
        model_rmsds = [safe_float(m.get("mean_epitope_rmsd")) for m in models]

        # Accumulate per-epitope RMSD across models
        ep_rmsd_accum = defaultdict(list)
        for m in models:
            raw = m.get("per_epitope_rmsd", "")
            for part in raw.split(";"):
                if "=" in part:
                    ep, val = part.split("=", 1)
                    v = safe_float(val)
                    if v is not None:
                        ep_rmsd_accum[ep.strip()].append(v)

        per_ep_mean_str = ";".join(
            f"{ep}={round(sum(vals)/len(vals), 3)}"
            for ep, vals in sorted(ep_rmsd_accum.items())
        )

        # ── SASA ─────────────────────────────────────────────────
        ep_sasas  = [safe_float(m.get("mean_epitope_sasa")) for m in models]
        tot_sasas = [safe_float(m.get("total_sasa"))        for m in models]

        # ── Best model (composite score) ──────────────────────────
        # Normalise each metric across this scaffold's models then sum:
        #   pLDDT (higher = better), RMSD (lower = better), SASA (higher = better)
        rmsd_fallback = [v if v is not None else 999.0 for v in model_rmsds]
        sasa_fallback = [v if v is not None else 0.0   for v in ep_sasas]
        plddt_fallback= [v if v is not None else 0.0   for v in g_plddts]

        plddt_n = minmax_norm(plddt_fallback)
        rmsd_n  = minmax_norm(rmsd_fallback, invert=True)
        sasa_n  = minmax_norm(sasa_fallback)

        composites = [
            (p or 0) + (r or 0) + (s or 0)
            for p, r, s in zip(plddt_n, rmsd_n, sasa_n)
        ]
        best_idx        = composites.index(max(composites))
        best_model      = models[best_idx]
        best_model_name = best_model["model_name"]

        # ── epitopes found ────────────────────────────────────────
        ep_found = models[0].get("epitopes_found", "")
        n_ep     = len([e for e in ep_found.split(";") if e])

        def r3(v): return round(v, 3) if v is not None else None
        def r4(v): return round(v, 4) if v is not None else None

        scaffolds.append({
            "scaffold_id":           scaffold_id,
            "n_models":              len(models),
            "epitopes_found":        ep_found,
            "n_epitopes":            n_ep,
            "mean_global_plddt":     r3(mean(g_plddts)),
            "std_global_plddt":      r3(std(g_plddts)),
            "best_global_plddt":     r3(max(v for v in g_plddts if v is not None)) if any(g_plddts) else None,
            "mean_epitope_plddt":    r3(mean(ep_plddts)),
            "mean_scaffold_plddt":   r3(mean(sc_plddts)),
            "mean_epitope_rmsd":     r4(mean(model_rmsds)),
            "std_epitope_rmsd":      r4(std(model_rmsds)),
            "best_epitope_rmsd":     r4(min(v for v in model_rmsds if v is not None)) if any(v is not None for v in model_rmsds) else None,
            "per_epitope_mean_rmsd": per_ep_mean_str,
            "mean_epitope_sasa":     r3(mean(ep_sasas)),
            "std_epitope_sasa":      r3(std(ep_sasas)),
            "mean_total_sasa":       r3(mean(tot_sasas)),
            "best_model_name":       best_model_name,
            "best_model_plddt":      r3(safe_float(best_model.get("global_mean_plddt"))),
            "best_model_rmsd":       r4(safe_float(best_model.get("mean_epitope_rmsd"))),
            "best_model_epitope_sasa": r3(safe_float(best_model.get("mean_epitope_sasa"))),
        })

    # Sort: lowest mean RMSD first; ties by highest pLDDT
    scaffolds.sort(key=lambda x: (
        x["mean_epitope_rmsd"] if x["mean_epitope_rmsd"] is not None else 999,
        -(x["mean_global_plddt"] or 0),
    ))
    return scaffolds


# ══════════════════════════════════════════════════════════════════
# 6. MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════

def run_analysis(af2_dir, epitopes_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    summary_csv          = os.path.join(output_dir, f"af2_metrics_summary_{timestamp}.csv")
    per_res_csv          = os.path.join(output_dir, f"af2_per_residue_{timestamp}.csv")
    rmsd_csv             = os.path.join(output_dir, f"af2_epitope_rmsd_{timestamp}.csv")
    scaffold_summary_csv = os.path.join(output_dir, f"scaffold_summary_{timestamp}.csv")

    init_csv(summary_csv, SUMMARY_COLS)
    init_csv(per_res_csv, PER_RESIDUE_COLS)
    init_csv(rmsd_csv,    RMSD_COLS)

    epitope_dict = build_epitope_dictionary(epitopes_dir)

    pdb_files = sorted(glob.glob(os.path.join(af2_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {af2_dir}")

    print(f"Analysing {len(pdb_files)} AlphaFold model(s)...\n")

    all_summary_rows = []

    for i, pdb_path in enumerate(pdb_files):
        model_name  = os.path.splitext(os.path.basename(pdb_path))[0]
        scaffold_id = parse_scaffold_id(model_name)
        print(f"[{i+1}/{len(pdb_files)}] {model_name}")
        print(f"  Scaffold ID     : {scaffold_id}")

        try:
            # ── pLDDT ────────────────────────────────────────────
            _, biopdb_model = load_af2_model(pdb_path)
            residue_data    = get_residue_plddt(biopdb_model)
            full_seq        = get_full_sequence(residue_data)
            print(f"  Sequence length : {len(full_seq)}")

            epitope_spans = find_epitope_spans(full_seq, epitope_dict)
            found         = list(epitope_spans.keys())
            print(f"  Epitopes found  : {found if found else 'none'}")

            plddt_m = compute_plddt_metrics(residue_data, epitope_spans)
            print(f"  Global pLDDT    : {plddt_m['global_mean_plddt']:.1f}")
            if plddt_m["mean_epitope_plddt"]:
                print(f"  Epitope pLDDT   : {plddt_m['mean_epitope_plddt']:.1f}"
                      f"  |  Scaffold: {plddt_m['mean_scaffold_plddt']:.1f}")

            # ── RMSD ─────────────────────────────────────────────
            rmsd_results  = compute_epitope_rmsd(residue_data, epitope_spans, epitopes_dir)
            ok_rmsds      = [r["rmsd"] for r in rmsd_results if r["rmsd"] is not None]
            mean_ep_rmsd  = round(sum(ok_rmsds) / len(ok_rmsds), 4) if ok_rmsds else None
            per_ep_rmsd_str = ";".join(
                f"{r['epitope_name']}={r['rmsd']}"
                for r in rmsd_results if r["rmsd"] is not None
            )
            for r in rmsd_results:
                tag = f"RMSD={r['rmsd']:.3f} Å" if r["rmsd"] is not None else r["status"]
                print(f"    [{r['epitope_name']}] {tag}")

            # ── SASA ─────────────────────────────────────────────
            total_sasa, sasa_per_res = compute_sasa_pyrosetta(
                pdb_path, epitope_spans, full_seq)
            ep_sasa   = [r["sasa"] for r in sasa_per_res if r["region"] == "EPITOPE"]
            scaf_sasa = [r["sasa"] for r in sasa_per_res if r["region"] == "SCAFFOLD"]
            mean_ep_sasa   = round(sum(ep_sasa)   / len(ep_sasa),   3) if ep_sasa   else None
            mean_scaf_sasa = round(sum(scaf_sasa) / len(scaf_sasa), 3) if scaf_sasa else None
            print(f"  Total SASA      : {total_sasa:.1f} Å²"
                  + (f"  |  Epitope: {mean_ep_sasa:.1f}" if mean_ep_sasa else ""))

            # ── Write per-model CSVs ──────────────────────────────
            summary_row = {
                "model_name":           model_name,
                "scaffold_id":          scaffold_id,
                "n_residues":           len(residue_data),
                "global_mean_plddt":    plddt_m["global_mean_plddt"],
                "global_median_plddt":  plddt_m["global_median_plddt"],
                "mean_epitope_plddt":   plddt_m["mean_epitope_plddt"],
                "mean_scaffold_plddt":  plddt_m["mean_scaffold_plddt"],
                "n_epitope_residues":   plddt_m["n_epitope_residues"],
                "n_scaffold_residues":  plddt_m["n_scaffold_residues"],
                "epitopes_found":       ";".join(found),
                "n_epitopes_found":     len(found),
                "mean_epitope_rmsd":    mean_ep_rmsd,
                "per_epitope_rmsd":     per_ep_rmsd_str,
                "total_sasa":           total_sasa,
                "mean_epitope_sasa":    mean_ep_sasa,
                "mean_scaffold_sasa":   mean_scaf_sasa,
            }
            append_csv(summary_csv, [summary_row], SUMMARY_COLS)
            all_summary_rows.append(summary_row)

            per_res_rows = []
            for rd, sd in zip(residue_data, sasa_per_res):
                per_res_rows.append({
                    "model_name":    model_name,
                    "rosetta_index": sd["rosetta_index"],
                    "pdb_resnum":    sd["pdb_resnum"],
                    "chain":         sd["chain"],
                    "resname":       sd["resname"],
                    "aa_1letter":    sd["aa_1letter"],
                    "region":        sd["region"],
                    "epitope_name":  sd["epitope_name"],
                    "plddt":         rd["plddt"],
                    "sasa":          sd["sasa"],
                })
            append_csv(per_res_csv, per_res_rows, PER_RESIDUE_COLS)

            rmsd_rows = [{"model_name": model_name, **r} for r in rmsd_results]
            append_csv(rmsd_csv, rmsd_rows, RMSD_COLS)

            print(f"  ✓ Done\n")

        except Exception as e:
            import traceback
            print(f"  ✗ ERROR: {e}")
            traceback.print_exc()
            print()
            continue

    # ── Per-scaffold aggregation ──────────────────────────────────
    print(f"{'='*55}")
    print(f"Aggregating {len(all_summary_rows)} models into per-scaffold summary...")

    scaffold_rows = aggregate_by_scaffold(all_summary_rows)
    write_csv(scaffold_summary_csv, scaffold_rows, SCAFFOLD_SUMMARY_COLS)

    print(f"\n  Top 5 scaffold designs:")
    print(f"  {'Scaffold':<35} {'N':>3}  {'RMSD':>7}  {'pLDDT':>7}  {'SASA':>7}")
    print(f"  {'-'*65}")
    for row in scaffold_rows[:5]:
        rmsd  = f"{row['mean_epitope_rmsd']:.3f}" if row['mean_epitope_rmsd'] else "N/A"
        plddt = f"{row['mean_global_plddt']:.1f}"  if row['mean_global_plddt'] else "N/A"
        sasa  = f"{row['mean_epitope_sasa']:.1f}"  if row['mean_epitope_sasa'] else "N/A"
        print(f"  {row['scaffold_id']:<35} {row['n_models']:>3}  "
              f"{rmsd:>7} Å  {plddt:>6}  {sasa:>6} Å²")

    print(f"\n{'='*55}")
    print(f"Analysis complete.")
    print(f"  Per-model summary  : {summary_csv}")
    print(f"  Per-residue data   : {per_res_csv}")
    print(f"  Epitope RMSD       : {rmsd_csv}")
    print(f"  Scaffold summary   : {scaffold_summary_csv}  ← main output")
    print(f"{'='*55}")


# ══════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("--af2_dir",      required=True,
                        help="Folder with AlphaFold PDB models")
    parser.add_argument("--epitopes_dir", required=True,
                        help="Folder with original epitope PDB files")
    parser.add_argument("--output_dir",   default="AF2_analysis/",
                        help="Output folder (default: AF2_analysis/)")
    args = parser.parse_args()
    run_analysis(args.af2_dir, args.epitopes_dir, args.output_dir)
