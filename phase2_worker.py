#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase2_worker.py — SS-aware cumulative grafting worker for GRAFFITI Phase 2.

MUST be in the same directory as step2.py and ss_utils.py.

For each selected motif on a scaffold:
  - HELIX / SHEET / MIXED → MotifGraft (backbone remodelling)
  - LOOP                  → sequence insertion into a compatible scaffold loop

The graft_method per motif is recorded in the result.
"""

import os


def worker(args):
    """
    Processes a full scaffold: cumulative grafting of all selected motifs.
    Returns a serialisable dict (no Rosetta objects).
    """
    (scaffold_id, selected_motifs,
     scaffold_dir, epitopes_dir, final_pdb_dir) = args

    # ── PyRosetta init ────────────────────────────────────────────
    try:
        import pyrosetta
        from pyrosetta import pose_from_pdb
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
        import os as _os
        _devnull      = _os.open(_os.devnull, _os.O_WRONLY)
        _saved_stdout = _os.dup(1)
        _saved_stderr = _os.dup(2)
        _os.dup2(_devnull, 1)
        _os.dup2(_devnull, 2)
        try:
            pyrosetta.init("-mute all", silent=True)
        finally:
            _os.dup2(_saved_stdout, 1)
            _os.dup2(_saved_stderr, 2)
            _os.close(_devnull)
            _os.close(_saved_stdout)
            _os.close(_saved_stderr)
    except Exception as e:
        return {
            "scaffold_id": scaffold_id, "status": "INIT_FAILED",
            "error": str(e), "grafts_ok": [], "graft_methods": {},
            "sasa_records": [], "avg_ep_sasa": None, "out_pdb": None,
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
        out, pdb_info = set(), pose.pdb_info()
        if pdb_info is None:
            return out
        for i in range(1, pose.total_residue() + 1):
            if any(lbl in labels for lbl in pdb_info.get_reslabels(i)):
                out.add(i)
        return out

    def remove_context(pose):
        ctx   = get_labeled(pose, {"CONTEXT"})
        clean = pose.clone()
        for i in range(clean.total_residue(), 0, -1):
            if i in ctx:
                clean.conformation().delete_residue_slow(i)
        return clean

    def calc_sasa(pose):
        c = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
        t = c.calculate(pose)
        return t, c.get_residue_sasa()

    def make_sasa_records(pose, scaffold_id, motifs_applied,
                          rsd_sasa, total_sasa):
        pdb_info = pose.pdb_info()
        motif_i  = get_labeled(pose, {"MOTIF"})
        conn_i   = get_labeled(pose, {"CONNECTION"})
        hot_i    = get_labeled(pose, {"HOTSPOT"})
        rows = []
        for i in range(1, pose.total_residue() + 1):
            if   i in motif_i: region = "MOTIF"
            elif i in conn_i:  region = "CONNECTION"
            elif i in hot_i:   region = "HOTSPOT"
            else:              region = "SCAFFOLD"
            rows.append({
                "scaffold_id":    scaffold_id,
                "motifs_applied": motifs_applied,
                "total_sasa":     f"{total_sasa:.3f}",
                "rosetta_index":  i,
                "pdb_resnum":     pdb_info.number(i) if pdb_info else i,
                "chain":          pdb_info.chain(i)  if pdb_info else "A",
                "resname":        pose.residue(i).name3(),
                "aa_1letter":     pose.residue(i).name1(),
                "region":         region,
                "sasa":           f"{rsd_sasa[i]:.3f}",
            })
        return rows

    # ── load scaffold ─────────────────────────────────────────────
    pdb_path = os.path.join(scaffold_dir, f"{scaffold_id}.pdb")
    if not os.path.exists(pdb_path):
        return {
            "scaffold_id": scaffold_id, "status": "PDB_NOT_FOUND",
            "error": pdb_path, "grafts_ok": [], "graft_methods": {},
            "sasa_records": [], "avg_ep_sasa": None, "out_pdb": None,
        }

    # import ss_utils once per worker
    try:
        from ss_utils import (classify_epitope_ss, find_compatible_loops,
                              graft_loop_with_modeling)
        ss_utils_ok = True
    except Exception:
        ss_utils_ok = False   # fallback: treat everything as structured

    sorted_motifs = sorted(selected_motifs,
                           key=lambda r: r["scaffold_range_start"])
    pose          = pose_from_pdb(pdb_path)
    grafts_ok     = []
    graft_methods = {}          # motif_name → method used
    used_positions = set()      # Rosetta indices already occupied

    for rec in sorted_motifs:
        motif_name = rec["motif_name"]
        motif_path = os.path.join(epitopes_dir, motif_name)
        if not os.path.exists(motif_path):
            continue

        # ── classify this epitope's SS ────────────────────────────
        ep_ss_class = "HELIX"   # default if ss_utils unavailable
        if ss_utils_ok:
            try:
                epitope_pose = pose_from_pdb(motif_path)
                ep_ss_class  = classify_epitope_ss(epitope_pose)
                ep_len       = epitope_pose.total_residue()
            except Exception:
                ep_ss_class = "HELIX"
                ep_len      = rec.get("motif_size", 8)
        else:
            ep_len = rec.get("motif_size", 8)

        # ── PATH A: structured → MotifGraft ──────────────────────
        if ep_ss_class in ("HELIX", "SHEET", "MIXED") or not ss_utils_ok:
            try:
                objs  = XmlObjects.create_from_string(
                            xml_motif_graft(pdb_path, motif_path))
                mover = objs.get_mover("motif_grafting")
                mover.apply(pose)
                pose  = remove_context(pose)

                # record positions used by this graft
                used_positions.update(
                    range(rec["scaffold_range_start"],
                          rec["scaffold_range_end"] + 1))

                grafts_ok.append(motif_name)
                graft_methods[motif_name] = f"MOTIFGRAFT({ep_ss_class})"

            except RuntimeError:
                graft_methods[motif_name] = f"MOTIFGRAFT_FAILED({ep_ss_class})"
                continue

        # ── PATH B: loop → structural loop modeling ─────────────
        else:   # ep_ss_class == "LOOP"
            try:
                compatible = find_compatible_loops(
                    pose, ep_len, already_used=used_positions)

                if not compatible:
                    graft_methods[motif_name] = "FAILED_NO_LOOP"
                    import sys
                    print(f"    [phase2] FAILED_NO_LOOP: {motif_name} "
                          f"ep_len={ep_len} on {scaffold_id}",
                          file=sys.stderr, flush=True)
                    continue

                loop_start, loop_end = compatible[0]
                import sys
                print(f"    [phase2] LOOP_INSERT: {motif_name} "
                      f"ep_len={ep_len} -> loop [{loop_start}-{loop_end}] "
                      f"len={loop_end - loop_start + 1}",
                      file=sys.stderr, flush=True)

                # epitope_pose already loaded above for classify_epitope_ss
                # reuse it; if not available, reload
                if "epitope_pose" not in dir():
                    epitope_pose = pose_from_pdb(motif_path)

                new_pose, meta = graft_loop_with_modeling(
                    pose, epitope_pose, loop_start, loop_end)

                # IMPORTANT: replace working pose with the remodelled one
                pose = new_pose

                # mark inserted positions as used (indices in new pose)
                used_positions.update(
                    range(meta["insert_start"], meta["insert_end"] + 1))

                grafts_ok.append(motif_name)
                graft_methods[motif_name] = "LOOP_MODEL"

            except Exception as e:
                import traceback, sys
                graft_methods[motif_name] = f"LOOP_MODEL_FAILED"
                print(f"    [phase2] LOOP_MODEL exception for {motif_name}: {e}",
                      file=sys.stderr, flush=True)
                traceback.print_exc(file=sys.stderr)
                continue

    if not grafts_ok:
        return {
            "scaffold_id":  scaffold_id, "status": "NO_GRAFTS",
            "error":        "", "grafts_ok": [], "graft_methods": graft_methods,
            "sasa_records": [], "avg_ep_sasa": None, "out_pdb": None,
        }

    # ── save PDB ──────────────────────────────────────────────────
    motifs_tag     = "_".join(os.path.splitext(m)[0] for m in grafts_ok)
    motifs_applied = ";".join(grafts_ok)
    out_pdb        = os.path.join(
                         final_pdb_dir, f"{scaffold_id}__{motifs_tag}.pdb")
    pose.dump_pdb(out_pdb)

    # ── SASA ─────────────────────────────────────────────────────
    sasa_records = []
    avg_ep_sasa  = None
    try:
        total_sasa, rsd_sasa = calc_sasa(pose)
        sasa_records = make_sasa_records(
            pose, scaffold_id, motifs_applied, rsd_sasa, total_sasa)
        ep_vals     = [float(r["sasa"]) for r in sasa_records
                       if r["region"] == "MOTIF"]
        avg_ep_sasa = round(sum(ep_vals) / len(ep_vals), 3) if ep_vals else 0.0
    except Exception:
        pass

    # Build human-readable graft method summary
    methods_summary = ";".join(
        f"{m}={graft_methods.get(m, '?')}" for m in grafts_ok)

    return {
        "scaffold_id":    scaffold_id,
        "status":         "SUCCESS",
        "grafts_ok":      grafts_ok,
        "graft_methods":  graft_methods,
        "methods_summary": methods_summary,
        "out_pdb":        out_pdb,
        "avg_ep_sasa":    avg_ep_sasa,
        "sasa_records":   sasa_records,
        "error":          "",
    }
