#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GRAFFITI — Unified Pipeline Runner
===================================
Grafting Routine for Automated Fragment Fitting,
Insertion, and Targeted Immunogenic epitopes

Pipeline stages:
  Step 1  — Parallel motif grafting (all scaffold × epitope pairs)
  Step 2  — Overlap resolution, SASA scoring, top scaffold selection
  Step 3  — ProteinMPNN sequence design on selected scaffolds
  [Step 4]  AlphaFold2 — run externally, results analysed in Step 5
  Step 5  — Post-AF2 analysis: pLDDT, RMSD, SASA per region

Usage:
  Run full pipeline from scratch:
    python graffiti_run.py --steps 1 2 3

  Resume from Step 2 (grafting already done):
    python graffiti_run.py --steps 2 3

  Run only post-AF2 analysis:
    python graffiti_run.py --steps 5

  Override any config value on the command line:
    python graffiti_run.py --steps 1 2 3 --cpus 16 --top_percent 10

  Dry run (print config and exit):
    python graffiti_run.py --dry_run
"""

import os
import sys
import glob
import argparse
import textwrap
import subprocess
from datetime import datetime

# ══════════════════════════════════════════════════════════════════════
# ██  MASTER CONFIG  — edit everything here, nowhere else
# ══════════════════════════════════════════════════════════════════════

CFG = {
    # ── Directories ──────────────────────────────────────────────────
    "scaffolds_dir":    "few_dummies/",   # input scaffold PDBs
    "epitopes_dir":     "Epitopes/",          # input epitope PDBs
    "grafts_dir":       "Grafts_individual/", # step 1 output
    "optimized_dir":    "Grafts_optimized/",  # step 2 output
    "mpnn_dir":         "MPNN_out/",          # step 3 output
    "af2_dir":          "AlphaFold_results/", # step 4 input (run externally)
    "af2_analysis_dir": "AF2_analysis/",      # step 5 output

    # ── Step 1: Grafting ─────────────────────────────────────────────
    "cpus":             8,       # parallel workers (None = all available)

    # ── Step 2: Overlap resolution + filtering ───────────────────────
    "overlap_buffer":   0,       # residue margin between graft ranges
    "top_percent":      20,      # % of top scaffolds to carry forward
    "min_grafts":       None,    # None = auto (max observed); int = hard floor
    "n_workers_step2":  8,       # workers for parallel PDB generation

    # ── Step 3: ProteinMPNN ──────────────────────────────────────────
    "mpnn_chain":       "A",     # chain to design
    "mpnn_num_seq":     10,      # sequences per scaffold
    "mpnn_temp":        "0.2",   # sampling temperature
    "mpnn_model":       "v_48_020",
    "mpnn_script":      "",      # path to protein_mpnn_run.py (auto-detected if empty)

    # ── Step 5: Post-AF2 analysis ────────────────────────────────────
    "rmsd_threshold":   2.0,     # Å — epitope RMSD cutoff for quality flag
    "plddt_threshold":  70.0,    # global pLDDT minimum for quality flag
}

# ══════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════

def banner(title, width=60):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print()
    print("=" * width)
    print(f"  {title}")
    print(f"  {ts}")
    print("=" * width)


def section(msg):
    print(f"\n── {msg}")


def check_dir(path, label):
    if not os.path.isdir(path):
        print(f"  ERROR: {label} directory not found: {path}")
        sys.exit(1)
    pdbs = glob.glob(os.path.join(path, "*.pdb"))
    print(f"  {label}: {path}  ({len(pdbs)} PDB files)")
    return pdbs


def find_latest_csv(directory, pattern):
    matches = sorted(glob.glob(os.path.join(directory, pattern)))
    if not matches:
        return None
    return matches[-1]


def run(cmd, label):
    """Run a shell command, streaming output, exit on failure."""
    print(f"\n  $ {cmd}\n")
    ret = subprocess.call(cmd, shell=True)
    if ret != 0:
        print(f"\n  ERROR: {label} exited with code {ret}")
        sys.exit(ret)


# ══════════════════════════════════════════════════════════════════════
# PIPELINE STEPS
# ══════════════════════════════════════════════════════════════════════

def step1_grafting(cfg):
    banner("STEP 1 — Parallel Motif Grafting")

    check_dir(cfg["scaffolds_dir"], "Scaffolds")
    check_dir(cfg["epitopes_dir"],  "Epitopes")

    cpus = cfg["cpus"] or ""
    cpus_arg = f"--cpus {cpus}" if cpus else ""

    run(
        f"python step1_paralel.py "
        f"--scaffolds {cfg['scaffolds_dir']} "
        f"--epitopes  {cfg['epitopes_dir']} "
        f"--output    {cfg['grafts_dir']} "
        f"{cpus_arg}",
        "Step 1 — Grafting"
    )

    csv = find_latest_csv(cfg["grafts_dir"], "graft_map_*.csv")
    if csv:
        print(f"\n  Output CSV: {csv}")
    else:
        print("\n  WARNING: no graft_map CSV found after step 1")


def step2_overlap(cfg):
    banner("STEP 2 — Overlap Resolution + Filtering")

    csv = find_latest_csv(cfg["grafts_dir"], "graft_map_*.csv")
    if not csv:
        print(f"  ERROR: no graft_map CSV found in {cfg['grafts_dir']}")
        print("  Run Step 1 first.")
        sys.exit(1)
    print(f"  Using CSV: {csv}")

    min_grafts_arg = f"--min_grafts {cfg['min_grafts']}" if cfg["min_grafts"] else ""
    run(
        f"python step2.py "
        f"--phase1_csv     {csv} "
        f"--scaffolds_dir  {cfg['scaffolds_dir']} "
        f"--epitopes_dir   {cfg['epitopes_dir']} "
        f"--output_dir     {cfg['optimized_dir']} "
        f"--overlap_buffer {cfg['overlap_buffer']} "
        f"--top_percent    {cfg['top_percent']} "
        f"--n_workers      {cfg['n_workers_step2']} "
        f"{min_grafts_arg}",
        "Step 2 — Overlap + Filtering"
    )

    top_dir = os.path.join(cfg["optimized_dir"], "top_pdbs/")
    pdbs = glob.glob(os.path.join(cfg["optimized_dir"], "**/*.pdb"), recursive=True)
    print(f"\n  Optimized PDBs generated: {len(pdbs)}")


def step3_mpnn(cfg):
    banner("STEP 3 — ProteinMPNN Sequence Design")

    # Find the top scaffold PDBs from step 2
    top_dirs = sorted(glob.glob(os.path.join(cfg["optimized_dir"], "top*pct*/")))
    if top_dirs:
        input_dir = top_dirs[-1]   # most recently created top% folder
    else:
        # fallback: use all final PDBs
        input_dir = os.path.join(cfg["optimized_dir"], "final_pdbs/")

    if not os.path.isdir(input_dir):
        print(f"  ERROR: could not find input PDBs for MPNN in {cfg['optimized_dir']}")
        print("  Run Step 2 first.")
        sys.exit(1)

    pdbs = glob.glob(os.path.join(input_dir, "*.pdb"))
    print(f"  Input PDBs: {input_dir}  ({len(pdbs)} files)")

    mpnn_script_arg = (f"--mpnn_script {cfg['mpnn_script']}"
                       if cfg["mpnn_script"] else "")

    run(
        f"python proteinmpnn_design.py "
        f"{input_dir} "
        f"{cfg['mpnn_dir']} "
        f"--epitopes_dir {cfg['epitopes_dir']} "
        f"--chain        {cfg['mpnn_chain']} "
        f"--num_seq      {cfg['mpnn_num_seq']} "
        f"--temp         {cfg['mpnn_temp']} "
        f"--model        {cfg['mpnn_model']} "
        f"{mpnn_script_arg}",
        "Step 3 — ProteinMPNN"
    )

    fas = glob.glob(os.path.join(cfg["mpnn_dir"], "**/*.fa"), recursive=True)
    print(f"\n  .fa files generated: {len(fas)}")
    print(f"\n  ┌─ NEXT STEP ────────────────────────────────────────")
    print(f"  │  Run AlphaFold2 on the sequences in: {cfg['mpnn_dir']}")
    print(f"  │  Then save the PDB models to:        {cfg['af2_dir']}")
    print(f"  │  Then run:  python graffiti_run.py --steps 5")
    print(f"  └────────────────────────────────────────────────────")


def step5_af2_analysis(cfg):
    banner("STEP 5 — Post-AlphaFold2 Analysis")

    check_dir(cfg["af2_dir"],      "AlphaFold results")
    check_dir(cfg["epitopes_dir"], "Epitopes")

    run(
        f"python post_alphafold_analysis.py "
        f"--af2_dir      {cfg['af2_dir']} "
        f"--epitopes_dir {cfg['epitopes_dir']} "
        f"--output_dir   {cfg['af2_analysis_dir']}",
        "Step 5 — Post-AF2 analysis"
    )

    csvs = glob.glob(os.path.join(cfg["af2_analysis_dir"], "*.csv"))
    print(f"\n  Analysis files: {len(csvs)}")


# ══════════════════════════════════════════════════════════════════════
# CLI + MAIN
# ══════════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            GRAFFITI — Unified Pipeline Runner

            Steps:
              1  Parallel motif grafting (scaffold × epitope)
              2  Overlap resolution, SASA scoring, top selection
              3  ProteinMPNN sequence design
              [4  AlphaFold2 — run externally]
              5  Post-AF2 analysis (pLDDT, RMSD, SASA)
        """),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--steps", nargs="+", type=int,
        default=[1, 2, 3],
        metavar="N",
        help="Steps to run, e.g. --steps 1 2 3  or  --steps 2  (default: 1 2 3)"
    )

    # override any CFG key from CLI
    parser.add_argument("--scaffolds_dir",    type=str)
    parser.add_argument("--epitopes_dir",     type=str)
    parser.add_argument("--grafts_dir",       type=str)
    parser.add_argument("--optimized_dir",    type=str)
    parser.add_argument("--mpnn_dir",         type=str)
    parser.add_argument("--af2_dir",          type=str)
    parser.add_argument("--af2_analysis_dir", type=str)
    parser.add_argument("--cpus",             type=int)
    parser.add_argument("--overlap_buffer",   type=int)
    parser.add_argument("--top_percent",      type=int)
    parser.add_argument("--min_grafts",       type=int)
    parser.add_argument("--n_workers_step2",  type=int)
    parser.add_argument("--mpnn_num_seq",     type=int)
    parser.add_argument("--mpnn_temp",        type=str)
    parser.add_argument("--mpnn_chain",       type=str)
    parser.add_argument("--mpnn_script",      type=str)
    parser.add_argument("--dry_run", action="store_true",
                        help="Print config and exit without running anything")

    return parser.parse_args()


def main():
    args = parse_args()

    # Apply CLI overrides into CFG
    for key in CFG:
        val = getattr(args, key, None)
        if val is not None:
            CFG[key] = val

    steps = sorted(set(args.steps))

    # ── print config ─────────────────────────────────────────────
    banner("GRAFFITI — Pipeline Configuration")
    print(f"  Steps to run : {steps}")
    print()
    print(f"  Directories")
    print(f"    Scaffolds        : {CFG['scaffolds_dir']}")
    print(f"    Epitopes         : {CFG['epitopes_dir']}")
    print(f"    Grafts (step 1)  : {CFG['grafts_dir']}")
    print(f"    Optimized (step2): {CFG['optimized_dir']}")
    print(f"    MPNN out (step3) : {CFG['mpnn_dir']}")
    print(f"    AF2 in  (step4)  : {CFG['af2_dir']}")
    print(f"    AF2 out (step5)  : {CFG['af2_analysis_dir']}")
    print()
    print(f"  Step 1 — Grafting")
    print(f"    CPUs             : {CFG['cpus'] or 'all available'}")
    print()
    print(f"  Step 2 — Overlap + Filter")
    print(f"    Overlap buffer   : {CFG['overlap_buffer']} residues")
    print(f"    Top percent      : {CFG['top_percent']}%")
    print(f"    Min grafts       : {CFG['min_grafts'] or 'auto (max observed)'}")
    print(f"    Workers          : {CFG['n_workers_step2'] or 'all available'}")
    print()
    print(f"  Step 3 — ProteinMPNN")
    print(f"    Chain            : {CFG['mpnn_chain']}")
    print(f"    Sequences/target : {CFG['mpnn_num_seq']}")
    print(f"    Temperature      : {CFG['mpnn_temp']}")
    print(f"    Model            : {CFG['mpnn_model']}")
    print()
    print(f"  Step 5 — Post-AF2")
    print(f"    pLDDT threshold  : {CFG['plddt_threshold']}")
    print(f"    RMSD threshold   : {CFG['rmsd_threshold']} Å")

    if args.dry_run:
        print("\n  Dry run — exiting.")
        return

    t_start = datetime.now()

    # ── run steps ─────────────────────────────────────────────────
    step_fns = {
        1: step1_grafting,
        2: step2_overlap,
        3: step3_mpnn,
        5: step5_af2_analysis,
    }

    for s in steps:
        if s == 4:
            banner("STEP 4 — AlphaFold2  [external]")
            print("  AlphaFold2 must be run externally.")
            print(f"  Input sequences : {CFG['mpnn_dir']}")
            print(f"  Save PDB models : {CFG['af2_dir']}")
            continue
        if s not in step_fns:
            print(f"  WARNING: unknown step {s}, skipping")
            continue
        step_fns[s](CFG)

    elapsed = datetime.now() - t_start
    banner(f"Pipeline finished  ({str(elapsed).split('.')[0]} elapsed)")


if __name__ == "__main__":
    main()
