#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-AlphaFold Analysis Pipeline

For each AlphaFold-modelled structure, this script:
  1. Loads the model and extracts pLDDT per residue (stored in B-factor column)
  2. Locates grafted epitopes in the sequence using the epitope dictionary
  3. Computes mean pLDDT for epitope regions vs scaffold regions
  4. Aligns each epitope region against the original epitope PDB → RMSD
  5. Calculates per-residue SASA using PyRosetta
  6. Saves everything to structured CSVs for downstream analysis

Outputs (in output_dir/):
  - af2_metrics_summary_TIMESTAMP.csv   → one row per model (global metrics)
  - af2_per_residue_TIMESTAMP.csv       → one row per residue (pLDDT + SASA + region)
  - af2_epitope_rmsd_TIMESTAMP.csv      → one row per (model, epitope) with RMSD

Usage:
    python post_alphafold_analysis.py \\
        --af2_dir      AlphaFold_results/ \\
        --epitopes_dir Epitopes/ \\
        --output_dir   AF2_analysis/ 
"""

import os
import glob
import csv
import math
import warnings
import argparse
from datetime import datetime
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser, PPBuilder, Superimposer
from Bio.SeqUtils import seq1

warnings.simplefilter("ignore")

# PyRosetta (for SASA)
import pyrosetta
from pyrosetta import pose_from_pdb
pyrosetta.init("-mute all")


# ═══════════════════════════════════════════════════════════
# EPITOPE DICTIONARY  (reused from earlier pipeline stages)
# ═══════════════════════════════════════════════════════════

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
    pdb_files    = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {epitopes_dir}")
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
    """
    Returns a dict: { epitope_name: (start_0based, end_0based_exclusive) }
    Only includes epitopes found in the sequence.
    """
    spans = {}
    for name, seq in epitope_dict.items():
        idx = full_sequence.find(seq)
        if idx != -1:
            spans[name] = (idx, idx + len(seq))
    return spans


# ═══════════════════════════════════════════════════════════
# 1. LOAD AF2 MODEL & EXTRACT pLDDT
# ═══════════════════════════════════════════════════════════

def load_af2_model(pdb_path):
    """
    Loads the AlphaFold PDB. Returns (structure, model, chain_sequence, residues).
    AlphaFold stores pLDDT in the B-factor column of each atom.
    We use CA atoms to get one value per residue.
    """
    parser    = PDBParser(QUIET=True)
    name      = os.path.splitext(os.path.basename(pdb_path))[0]
    structure = parser.get_structure(name, pdb_path)
    biopdb_model = structure[0]
    return structure, biopdb_model


def get_residue_plddt(biopdb_model):
    """
    Extracts pLDDT per residue from B-factor of CA atoms.
    Returns a list of dicts with residue info and pLDDT value,
    in the order they appear in the model (all chains concatenated).
    """
    residue_data = []
    for chain in biopdb_model:
        for residue in chain:
            if residue.get_id()[0] != " ":
                continue  # skip HETATM
            if "CA" not in residue:
                continue
            ca     = residue["CA"]
            plddt  = ca.get_bfactor()
            res_id = residue.get_id()[1]
            resname= residue.resname
            aa     = seq1(resname) if resname != "UNK" else "X"
            residue_data.append({
                "chain":   chain.get_id(),
                "res_id":  res_id,
                "resname": resname,
                "aa":      aa,
                "plddt":   plddt,
                "residue_obj": residue,
            })
    return residue_data


def get_full_sequence(residue_data):
    """Returns the full 1-letter sequence from residue_data list."""
    return "".join(r["aa"] for r in residue_data)


# ═══════════════════════════════════════════════════════════
# 2. pLDDT METRICS
# ═══════════════════════════════════════════════════════════

def compute_plddt_metrics(residue_data, epitope_spans):
    """
    Computes pLDDT statistics globally, per region (EPITOPE vs SCAFFOLD),
    and per individual epitope.

    Returns a dict of metrics.
    """
    all_plddt     = [r["plddt"] for r in residue_data]
    global_mean   = sum(all_plddt) / len(all_plddt) if all_plddt else 0.0
    global_median = float(np.median(all_plddt)) if all_plddt else 0.0

    # Tag each residue with its region
    epitope_mask = [False] * len(residue_data)
    for name, (start, end) in epitope_spans.items():
        for i in range(start, end):
            if i < len(epitope_mask):
                epitope_mask[i] = True

    ep_plddt   = [r["plddt"] for i, r in enumerate(residue_data) if epitope_mask[i]]
    scaf_plddt = [r["plddt"] for i, r in enumerate(residue_data) if not epitope_mask[i]]

    mean_ep_plddt   = sum(ep_plddt)   / len(ep_plddt)   if ep_plddt   else None
    mean_scaf_plddt = sum(scaf_plddt) / len(scaf_plddt) if scaf_plddt else None

    # Per-epitope pLDDT
    per_epitope = {}
    for name, (start, end) in epitope_spans.items():
        vals = [residue_data[i]["plddt"] for i in range(start, end)
                if i < len(residue_data)]
        per_epitope[name] = sum(vals) / len(vals) if vals else None

    return {
        "global_mean_plddt":    round(global_mean,   3),
        "global_median_plddt":  round(global_median, 3),
        "mean_epitope_plddt":   round(mean_ep_plddt,   3) if mean_ep_plddt   else None,
        "mean_scaffold_plddt":  round(mean_scaf_plddt, 3) if mean_scaf_plddt else None,
        "n_epitope_residues":   len(ep_plddt),
        "n_scaffold_residues":  len(scaf_plddt),
        "per_epitope_plddt":    per_epitope,
        "epitope_mask":         epitope_mask,
    }


# ═══════════════════════════════════════════════════════════
# 3. RMSD — EPITOPE vs ORIGINAL
# ═══════════════════════════════════════════════════════════

def get_backbone_atoms(residue_obj, atom_names=("CA", "N", "C", "O")):
    """Returns a list of atoms present in the residue from atom_names."""
    atoms = []
    for name in atom_names:
        if name in residue_obj:
            atoms.append(residue_obj[name])
    return atoms


def compute_epitope_rmsd(af2_residue_data, epitope_spans, epitopes_dir):
    """
    For each epitope found in the AF2 model:
      1. Extracts the corresponding residues from the AF2 model
      2. Loads the original epitope PDB
      3. Superimposes using Bio.PDB Superimposer (CA atoms)
      4. Returns the RMSD

    Returns a list of dicts: [{epitope_name, n_atoms, rmsd, status}]
    """
    parser = PDBParser(QUIET=True)
    results = []

    for epitope_name, (start, end) in epitope_spans.items():
        epitope_pdb = os.path.join(epitopes_dir, f"{epitope_name}.pdb")
        if not os.path.exists(epitope_pdb):
            results.append({
                "epitope_name": epitope_name,
                "n_ca_atoms":   0,
                "rmsd":         None,
                "status":       "ORIGINAL_PDB_NOT_FOUND",
            })
            continue

        try:
            # AF2 residues in the epitope span
            af2_residues = af2_residue_data[start:end]

            # Original epitope structure
            orig_structure = parser.get_structure(epitope_name, epitope_pdb)
            orig_model     = orig_structure[0]
            orig_residues  = [
                r for chain in orig_model
                for r in chain
                if r.get_id()[0] == " " and "CA" in r
            ]

            n_af2 = len(af2_residues)
            n_ori = len(orig_residues)

            # Align on the shorter length (handle minor size mismatches)
            n_align = min(n_af2, n_ori)
            if n_align == 0:
                raise ValueError("No CA atoms to align")

            af2_ca  = [r["residue_obj"]["CA"] for r in af2_residues[:n_align]
                       if "CA" in r["residue_obj"]]
            orig_ca = [r["CA"] for r in orig_residues[:n_align]]

            n_common = min(len(af2_ca), len(orig_ca))
            if n_common < 3:
                raise ValueError(f"Too few CA atoms for superimposition ({n_common})")

            sup = Superimposer()
            sup.set_atoms(orig_ca[:n_common], af2_ca[:n_common])
            rmsd = round(sup.rms, 4)

            results.append({
                "epitope_name": epitope_name,
                "n_ca_atoms":   n_common,
                "rmsd":         rmsd,
                "status":       "OK",
            })

        except Exception as e:
            results.append({
                "epitope_name": epitope_name,
                "n_ca_atoms":   0,
                "rmsd":         None,
                "status":       f"ERROR: {e}",
            })

    return results


# ═══════════════════════════════════════════════════════════
# 4. SASA (PyRosetta)
# ═══════════════════════════════════════════════════════════

def compute_sasa_pyrosetta(pdb_path, epitope_spans, full_sequence):
    """
    Loads the PDB into PyRosetta and computes per-residue SASA.
    Maps SASA values back to regions (EPITOPE / SCAFFOLD).

    Returns (total_sasa, per_residue_list) where per_residue_list
    has one dict per residue with sasa and region.
    """
    pose       = pose_from_pdb(pdb_path)
    calc       = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
    total_sasa = calc.calculate(pose)
    rsd_sasa   = calc.get_residue_sasa()

    # Build epitope mask over the full sequence
    epitope_mask = [False] * len(full_sequence)
    epitope_name_map = {}
    for name, (start, end) in epitope_spans.items():
        for i in range(start, min(end, len(full_sequence))):
            epitope_mask[i] = True
            epitope_name_map[i] = name

    pdb_info = pose.pdb_info()
    per_residue = []

    for i in range(1, pose.total_residue() + 1):
        seq_idx = i - 1  # 0-based
        region  = "EPITOPE" if (seq_idx < len(epitope_mask) and epitope_mask[seq_idx]) \
                  else "SCAFFOLD"
        ep_name = epitope_name_map.get(seq_idx, "")

        per_residue.append({
            "rosetta_index": i,
            "pdb_resnum":   pdb_info.number(i) if pdb_info else i,
            "chain":        pdb_info.chain(i)  if pdb_info else "A",
            "resname":      pose.residue(i).name3(),
            "aa_1letter":   pose.residue(i).name1(),
            "region":       region,
            "epitope_name": ep_name,
            "sasa":         round(rsd_sasa[i], 3),
        })

    return round(total_sasa, 3), per_residue


# ═══════════════════════════════════════════════════════════
# 5. CSV WRITERS
# ═══════════════════════════════════════════════════════════

SUMMARY_COLS = [
    "model_name",
    "n_residues",
    "global_mean_plddt",
    "global_median_plddt",
    "mean_epitope_plddt",
    "mean_scaffold_plddt",
    "n_epitope_residues",
    "n_scaffold_residues",
    "epitopes_found",
    "total_sasa",
    "mean_epitope_sasa",
    "mean_scaffold_sasa",
]

PER_RESIDUE_COLS = [
    "model_name",
    "rosetta_index",
    "pdb_resnum",
    "chain",
    "resname",
    "aa_1letter",
    "region",
    "epitope_name",
    "plddt",
    "sasa",
]

RMSD_COLS = [
    "model_name",
    "epitope_name",
    "n_ca_atoms",
    "rmsd",
    "status",
]


def init_csv(path, cols):
    with open(path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=cols).writeheader()


def append_csv(path, rows, cols):
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in cols})


# ═══════════════════════════════════════════════════════════
# 6. MAIN PIPELINE
# ═══════════════════════════════════════════════════════════

def run_analysis(af2_dir, epitopes_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    summary_csv    = os.path.join(output_dir, f"af2_metrics_summary_{timestamp}.csv")
    per_res_csv    = os.path.join(output_dir, f"af2_per_residue_{timestamp}.csv")
    rmsd_csv       = os.path.join(output_dir, f"af2_epitope_rmsd_{timestamp}.csv")

    init_csv(summary_csv, SUMMARY_COLS)
    init_csv(per_res_csv, PER_RESIDUE_COLS)
    init_csv(rmsd_csv,    RMSD_COLS)

    epitope_dict = build_epitope_dictionary(epitopes_dir)

    pdb_files = sorted(glob.glob(os.path.join(af2_dir, "*.pdb")))
    if not pdb_files:
        raise FileNotFoundError(f"No PDB files found in: {af2_dir}")

    print(f"Analysing {len(pdb_files)} AlphaFold model(s)...\n")

    for i, pdb_path in enumerate(pdb_files):
        model_name = os.path.splitext(os.path.basename(pdb_path))[0]
        print(f"[{i+1}/{len(pdb_files)}] {model_name}")

        try:
            # ── Step 1: load model & pLDDT ──────────────────────────
            structure, biopdb_model = load_af2_model(pdb_path)
            residue_data = get_residue_plddt(biopdb_model)
            full_seq     = get_full_sequence(residue_data)
            print(f"  Sequence length : {len(full_seq)}")

            # ── Step 2: locate epitopes ──────────────────────────────
            epitope_spans = find_epitope_spans(full_seq, epitope_dict)
            found = list(epitope_spans.keys())
            print(f"  Epitopes found  : {found if found else 'none'}")

            # ── Step 3: pLDDT metrics ────────────────────────────────
            plddt_metrics = compute_plddt_metrics(residue_data, epitope_spans)
            print(f"  Global pLDDT    : {plddt_metrics['global_mean_plddt']:.1f}")
            if plddt_metrics["mean_epitope_plddt"] is not None:
                print(f"  Epitope pLDDT   : {plddt_metrics['mean_epitope_plddt']:.1f}")
                print(f"  Scaffold pLDDT  : {plddt_metrics['mean_scaffold_plddt']:.1f}")

            # Per-epitope pLDDT print
            for ep_name, ep_plddt in plddt_metrics["per_epitope_plddt"].items():
                if ep_plddt is not None:
                    print(f"    [{ep_name}] mean pLDDT: {ep_plddt:.1f}")

            # ── Step 4: RMSD vs original epitopes ───────────────────
            rmsd_results = compute_epitope_rmsd(
                residue_data, epitope_spans, epitopes_dir
            )
            for r in rmsd_results:
                status = f"RMSD={r['rmsd']:.3f} Å" if r["rmsd"] is not None \
                         else r["status"]
                print(f"    [{r['epitope_name']}] {status}")

            # ── Step 5: SASA ─────────────────────────────────────────
            total_sasa, sasa_per_res = compute_sasa_pyrosetta(
                pdb_path, epitope_spans, full_seq
            )
            ep_sasa   = [r["sasa"] for r in sasa_per_res if r["region"] == "EPITOPE"]
            scaf_sasa = [r["sasa"] for r in sasa_per_res if r["region"] == "SCAFFOLD"]
            mean_ep_sasa   = sum(ep_sasa)   / len(ep_sasa)   if ep_sasa   else None
            mean_scaf_sasa = sum(scaf_sasa) / len(scaf_sasa) if scaf_sasa else None
            print(f"  Total SASA      : {total_sasa:.1f} Å²")
            if mean_ep_sasa:
                print(f"  Mean epitope SASA : {mean_ep_sasa:.1f} Å²")

            # ── Step 6: Write CSVs ───────────────────────────────────

            # Summary row
            append_csv(summary_csv, [{
                "model_name":           model_name,
                "n_residues":           len(residue_data),
                "global_mean_plddt":    plddt_metrics["global_mean_plddt"],
                "global_median_plddt":  plddt_metrics["global_median_plddt"],
                "mean_epitope_plddt":   plddt_metrics["mean_epitope_plddt"],
                "mean_scaffold_plddt":  plddt_metrics["mean_scaffold_plddt"],
                "n_epitope_residues":   plddt_metrics["n_epitope_residues"],
                "n_scaffold_residues":  plddt_metrics["n_scaffold_residues"],
                "epitopes_found":       ";".join(found),
                "total_sasa":           total_sasa,
                "mean_epitope_sasa":    round(mean_ep_sasa,   3) if mean_ep_sasa   else "",
                "mean_scaffold_sasa":   round(mean_scaf_sasa, 3) if mean_scaf_sasa else "",
            }], SUMMARY_COLS)

            # Per-residue rows (merge pLDDT + SASA)
            per_res_rows = []
            for j, (rd, sd) in enumerate(zip(residue_data, sasa_per_res)):
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

            # RMSD rows
            rmsd_rows = [{"model_name": model_name, **r} for r in rmsd_results]
            append_csv(rmsd_csv, rmsd_rows, RMSD_COLS)

            print(f"  ✓ Done\n")

        except Exception as e:
            print(f"  ✗ ERROR processing {model_name}: {e}\n")
            continue

    print(f"{'='*55}")
    print(f"Analysis complete.")
    print(f"  Summary CSV    : {summary_csv}")
    print(f"  Per-residue CSV: {per_res_csv}")
    print(f"  RMSD CSV       : {rmsd_csv}")
    print(f"{'='*55}")


# ═══════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════

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