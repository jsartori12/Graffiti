#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GRAFFITI - Phase 1: Parallel Motif Grafting

Tests every (scaffold, epitope) pair independently using multiprocessing.
Each pair runs in its own spawned process with its own PyRosetta instance.

Requires phase1_worker.py in the same directory.

Usage:
    python phase1_grafting.py                         # uses all CPUs
    python phase1_grafting.py --cpus 8                # limit to 8
    python phase1_grafting.py --cpus 1                # sequential (debug)
    python phase1_grafting.py --scaffolds my_pdbs/ --epitopes my_eps/
"""

import os
import glob
import csv
import argparse
import multiprocessing
from datetime import datetime

from phase1_worker import worker


# =====================================================================
# CONFIG
# =====================================================================
SCAFFOLD_DIR = "relaxed_results/"
EPITOPES_DIR = "Epitopes/"
OUTPUT_DIR   = "Grafts_individual/"

CSV_COLUMNS = [
    "scaffold_id", "motif_name", "status",
    "graft_method",       # MOTIFGRAFT | LOOP_INSERT | FAILED_NO_LOOP
    "epitope_ss_class",   # HELIX | SHEET | LOOP | MIXED
    "scaffold_range_start", "scaffold_range_end",
    "pdb_chain", "pdb_start_resnum", "pdb_end_resnum",
    "motif_size", "sequence", "connection_resnums",
    "scaffold_total_residues", "error_message",
]


# =====================================================================
# REAL-TIME PROGRESS CALLBACK
# Called by apply_async as each result arrives.
# =====================================================================
def make_progress_callback(total):
    counter = {"n": 0}
    def callback(result):
        counter["n"] += 1
        n      = counter["n"]
        sid    = result.get("scaffold_id", "?")
        motif  = result.get("motif_name",  "?")
        status = result.get("status",      "?")
        icon   = "OK" if status == "SUCCESS" else "XX"
        pct    = n / total * 100
        print(f"  [{n:>5}/{total}  {pct:5.1f}%]  {icon}  {sid} + {motif}")
    return callback


# =====================================================================
# MAIN
# =====================================================================
def main():
    n_available = multiprocessing.cpu_count()

    parser = argparse.ArgumentParser(
        description="GRAFFITI Phase 1 - Parallel Motif Grafting"
    )
    parser.add_argument(
        "--cpus", type=int, default=n_available,
        help=f"Number of CPUs (default: all = {n_available})"
    )
    parser.add_argument(
        "--scaffolds", type=str, default=SCAFFOLD_DIR,
        help=f"Scaffold PDB directory (default: {SCAFFOLD_DIR})"
    )
    parser.add_argument(
        "--epitopes", type=str, default=EPITOPES_DIR,
        help=f"Epitope PDB directory (default: {EPITOPES_DIR})"
    )
    parser.add_argument(
        "--output", type=str, default=OUTPUT_DIR,
        help=f"Output directory (default: {OUTPUT_DIR})"
    )
    args = parser.parse_args()

    n_workers    = max(1, min(args.cpus, n_available))
    scaffold_dir = args.scaffolds
    epitopes_dir = args.epitopes
    output_dir   = args.output

    os.makedirs(output_dir, exist_ok=True)

    pdb_files   = sorted(glob.glob(os.path.join(scaffold_dir, "*.pdb")))
    motifs_list = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))

    if not pdb_files:
        print(f"Error: no PDB files found in {scaffold_dir}")
        return
    if not motifs_list:
        print(f"Error: no PDB files found in {epitopes_dir}")
        return

    tasks = [
        (pdb, motif, output_dir)
        for pdb   in pdb_files
        for motif in motifs_list
    ]
    total = len(tasks)

    print("=" * 55)
    print("GRAFFITI - Phase 1: Parallel Motif Grafting")
    print("=" * 55)
    print(f"  Scaffolds  : {len(pdb_files)}")
    print(f"  Epitopes   : {len(motifs_list)}")
    print(f"  Total pairs: {total}")
    print(f"  CPUs       : {n_workers} / {n_available} available")
    print(f"  Output     : {output_dir}")
    print("=" * 55)
    print()

    # spawn: each worker gets a clean Python interpreter
    # Required for PyRosetta (not fork-safe due to C++ global state)
    # apply_async + callback gives real-time progress vs silent pool.map
    ctx      = multiprocessing.get_context("spawn")
    callback = make_progress_callback(total)

    with ctx.Pool(processes=n_workers) as pool:
        async_results = [
            pool.apply_async(worker, (task,), callback=callback)
            for task in tasks
        ]
        pool.close()
        pool.join()

    results = [ar.get() for ar in async_results]

    # save CSV
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path  = os.path.join(output_dir, f"graft_map_{timestamp}.csv")

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in results:
            writer.writerow({col: row.get(col, "") for col in CSV_COLUMNS})

    success = sum(1 for r in results if r.get("status") == "SUCCESS")
    failed  = total - success

    print()
    print("=" * 55)
    print("Phase 1 complete.")
    print(f"  Total pairs  : {total}")
    print(f"  Success      : {success}  ({success/total*100:.1f}%)")
    print(f"  Failed       : {failed}")
    print(f"  CSV report   : {csv_path}")
    print("=" * 55)

    failures = [r for r in results if r.get("status") != "SUCCESS"]
    if failures:
        n_show = min(5, len(failures))
        print(f"\nFirst {n_show} failure(s):")
        for r in failures[:n_show]:
            print(f"  {r['scaffold_id']} + {r['motif_name']}")
            print(f"    {r.get('error_message', '')}")


if __name__ == "__main__":
    main()
