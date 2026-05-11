#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase1_worker.py — SS-aware worker for GRAFFITI Phase 1.

MUST be in the same directory as step1_paralel.py and ss_utils.py.

For each (scaffold, epitope) pair:
  - If the epitope has defined secondary structure (HELIX / SHEET / MIXED):
      → Standard MotifGraft (PyRosetta RosettaScripts)
  - If the epitope is predominantly a loop (>= LOOP_FRACTION coil):
      → Find compatible loop regions in the scaffold and mutate their
        sequence to match the epitope (no backbone remodelling)

The graft_method column in the CSV records which path was taken:
  MOTIFGRAFT      — structured epitope, MotifGraft succeeded
  LOOP_INSERT     — loop epitope, sequence insertion succeeded
  FAILED_NO_LOOP  — loop epitope but no compatible loop found in scaffold
  FAILED          — structured epitope, MotifGraft raised an exception
  INIT_FAILED     — PyRosetta could not initialise
"""

import os


def worker(args):
    pdb_path, motif_path, output_dir = args

    scaffold_id = os.path.splitext(os.path.basename(pdb_path))[0]
    motif_name  = os.path.basename(motif_path)

    # ── PyRosetta init ────────────────────────────────────────────
    try:
        import pyrosetta
        from pyrosetta import pose_from_pdb
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
        pyrosetta.init("-ex1 -ex2 -mute all")
    except Exception as e:
        return {
            "scaffold_id":   scaffold_id,
            "motif_name":    motif_name,
            "status":        "INIT_FAILED",
            "graft_method":  "INIT_FAILED",
            "error_message": str(e),
        }

    # ── local helpers ─────────────────────────────────────────────
    def xml_motif_graft(context_path, frag_path):
        return f"""
        <ROSETTASCRIPTS>
          <MOVERS>
            <MotifGraft name="motif_grafting"
              context_structure="{context_path}"
              motif_structure="{frag_path}"
              RMSD_tolerance="3"
              NC_points_RMSD_tolerance="10.0"
              clash_score_cutoff="5"
              clash_test_residue="GLY"
              combinatory_fragment_size_delta="0:0"
              full_motif_bb_alignment="1"
              graft_only_hotspots_by_replacement="0"
              revert_graft_to_native_sequence="0"/>
          </MOVERS>
          <PROTOCOLS><Add mover="motif_grafting"/></PROTOCOLS>
        </ROSETTASCRIPTS>"""

    def get_labeled(pose, labels):
        out      = set()
        pdb_info = pose.pdb_info()
        if pdb_info is None:
            return out
        for i in range(1, pose.total_residue() + 1):
            if any(lbl in labels for lbl in pdb_info.get_reslabels(i)):
                out.add(i)
        return out

    def extract_motif_ranges(pose):
        pdb_info  = pose.pdb_info()
        motif_idx = sorted(get_labeled(pose, {"MOTIF"}))
        conn_idx  = sorted(get_labeled(pose, {"CONNECTION"}))
        if not motif_idx:
            return None
        pdb_nums = [pdb_info.number(i) for i in motif_idx]
        chains   = [pdb_info.chain(i)  for i in motif_idx]
        sequence = "".join(pose.residue(i).name1() for i in motif_idx)
        return {
            "pdb_resnums":          pdb_nums,
            "chains":               chains,
            "start_resnum":         pdb_nums[0],
            "end_resnum":           pdb_nums[-1],
            "start_chain":          chains[0],
            "motif_size":           len(motif_idx),
            "connection_resnums":   [pdb_info.number(i) for i in conn_idx],
            "sequence":             sequence,
            "scaffold_range_start": motif_idx[0],
            "scaffold_range_end":   motif_idx[-1],
        }

    def remove_context(pose):
        ctx   = get_labeled(pose, {"CONTEXT"})
        clean = pose.clone()
        for i in range(clean.total_residue(), 0, -1):
            if i in ctx:
                clean.conformation().delete_residue_slow(i)
        return clean

    # ── load both poses ───────────────────────────────────────────
    try:
        scaffold_pose = pose_from_pdb(pdb_path)
        epitope_pose  = pose_from_pdb(motif_path)
        scaffold_size = scaffold_pose.total_residue()
        epitope_len   = epitope_pose.total_residue()
    except Exception as e:
        return {
            "scaffold_id":   scaffold_id,
            "motif_name":    motif_name,
            "status":        "FAILED",
            "graft_method":  "FAILED",
            "error_message": f"PDB load failed: {e}",
        }

    # ── classify epitope SS ───────────────────────────────────────
    try:
        from ss_utils import classify_epitope_ss, find_compatible_loops, \
                             graft_loop_with_modeling
        ep_ss_class = classify_epitope_ss(epitope_pose)
    except Exception as e:
        # if ss_utils fails, fall back to MotifGraft
        ep_ss_class = "HELIX"

    # ── BRANCH ───────────────────────────────────────────────────
    # ── PATH A: structured epitope → MotifGraft ──────────────────
    if ep_ss_class in ("HELIX", "SHEET", "MIXED"):
        try:
            objs  = XmlObjects.create_from_string(
                        xml_motif_graft(pdb_path, motif_path))
            mover = objs.get_mover("motif_grafting")
            mover.apply(scaffold_pose)

            ranges = extract_motif_ranges(scaffold_pose)
            if ranges is None:
                raise RuntimeError("No MOTIF residues labeled after graft")

            clean_pose = remove_context(scaffold_pose)
            pdb_out    = os.path.join(output_dir, f"{scaffold_id}__{motif_name}")
            clean_pose.dump_pdb(pdb_out)

            return {
                "scaffold_id":             scaffold_id,
                "motif_name":              motif_name,
                "status":                  "SUCCESS",
                "graft_method":            "MOTIFGRAFT",
                "epitope_ss_class":        ep_ss_class,
                "scaffold_range_start":    ranges["scaffold_range_start"],
                "scaffold_range_end":      ranges["scaffold_range_end"],
                "pdb_chain":               ranges["start_chain"],
                "pdb_start_resnum":        ranges["start_resnum"],
                "pdb_end_resnum":          ranges["end_resnum"],
                "motif_size":              ranges["motif_size"],
                "sequence":                ranges["sequence"],
                "connection_resnums":      ";".join(
                                               str(r) for r in
                                               ranges["connection_resnums"]),
                "scaffold_total_residues": scaffold_size,
                "error_message":           "",
            }

        except Exception as e:
            return {
                "scaffold_id":      scaffold_id,
                "motif_name":       motif_name,
                "status":           "FAILED",
                "graft_method":     "MOTIFGRAFT",
                "epitope_ss_class": ep_ss_class,
                "error_message":    str(e).strip().splitlines()[0],
            }

    # ── PATH B: loop epitope → structural loop modeling ───────────
    else:  # ep_ss_class == "LOOP"
        import sys, traceback
        try:
            compatible_loops = find_compatible_loops(
                scaffold_pose, epitope_len)

            if not compatible_loops:
                print(f"  [phase1] FAILED_NO_LOOP: {scaffold_id} + {motif_name} "
                      f"ep_len={epitope_len}", file=sys.stderr, flush=True)
                return {
                    "scaffold_id":      scaffold_id,
                    "motif_name":       motif_name,
                    "status":           "FAILED",
                    "graft_method":     "FAILED_NO_LOOP",
                    "epitope_ss_class": ep_ss_class,
                    "error_message":    (
                        f"No compatible loop (ep_len={epitope_len}, "
                        f"scaffold={scaffold_id})"
                    ),
                }

            loop_start, loop_end = compatible_loops[0]
            print(f"  [phase1] LOOP_MODEL: {scaffold_id} + {motif_name} "
                  f"ep_len={epitope_len} -> loop [{loop_start}-{loop_end}] "
                  f"len={loop_end - loop_start + 1}",
                  file=sys.stderr, flush=True)

            remodelled_pose, meta = graft_loop_with_modeling(
                scaffold_pose, epitope_pose, loop_start, loop_end)

            pdb_out = os.path.join(output_dir, f"{scaffold_id}__{motif_name}")
            remodelled_pose.dump_pdb(pdb_out)

            ep_seq = "".join(
                epitope_pose.residue(i).name1()
                for i in range(1, epitope_pose.total_residue() + 1)
            )

            return {
                "scaffold_id":             scaffold_id,
                "motif_name":              motif_name,
                "status":                  "SUCCESS",
                "graft_method":            "LOOP_MODEL",
                "epitope_ss_class":        ep_ss_class,
                "scaffold_range_start":    meta["scaffold_range_start"],
                "scaffold_range_end":      meta["scaffold_range_end"],
                "pdb_chain":               "A",
                "pdb_start_resnum":        meta["insert_start"],
                "pdb_end_resnum":          meta["insert_end"],
                "motif_size":              epitope_len,
                "sequence":                ep_seq,
                "connection_resnums":      "",
                "scaffold_total_residues": scaffold_size,
                "error_message":           "",
            }

        except Exception as e:
            print(f"  [phase1] LOOP_MODEL exception: {scaffold_id} + {motif_name}: {e}",
                  file=sys.stderr, flush=True)
            traceback.print_exc(file=sys.stderr)
            return {
                "scaffold_id":      scaffold_id,
                "motif_name":       motif_name,
                "status":           "FAILED",
                "graft_method":     "LOOP_MODEL",
                "epitope_ss_class": ep_ss_class,
                "error_message":    str(e).strip().splitlines()[0],
            }
