#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 16:52:17 2026

@author: joao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ss_utils.py — Secondary structure utilities for GRAFFITI.

Shared by phase1_worker.py and phase2_worker.py.
Must be in the same directory as both workers.

Updated: 
- Epitope LOOP classification now triggers sequence mutation if 
  epitope length <= scaffold loop length.
- Documentation in English.
"""

from pyrosetta import Pose
from pyrosetta.rosetta.protocols.grafting import delete_region, insert_pose_into_pose
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta import *

# ── Thresholds ────────────────────────────────────────────────────────
LOOP_FRACTION   = 0.6    # ≥60% coil → LOOP
MAX_SIZE_DELTA  = 5     # scaffold loop length within ±20 of epitope
FLANK_RESIDUES  = 3      # residues on each side included in CCD/relax
CCD_CYCLES       = 200    # CCD iterations
RELAX_ROUNDS    = 0      # FastRelax rounds (0 = minimal)

# ══════════════════════════════════════════════════════════════════════
# SECONDARY STRUCTURE CLASSIFICATION
# ══════════════════════════════════════════════════════════════════════

def get_dssp(pose):
    """
    Runs DSSP on pose, inserts SS into pose, returns the full SS string.
    """
    from pyrosetta.rosetta.core.scoring.dssp import Dssp
    dssp_obj = Dssp(pose)
    dssp_obj.insert_ss_into_pose(pose)
    return pose.secstruct()


def classify_epitope_ss(epitope_pose):
    """
    Classifies the dominant secondary structure of an epitope pose.
    Returns: "HELIX" | "SHEET" | "LOOP" | "MIXED"
    """
    ss = get_dssp(epitope_pose)
    n  = len(ss)
    if n == 0:
        return "LOOP"

    n_H = ss.count("H")
    n_E = ss.count("E")
    n_L = ss.count("L")

    if n_L / n >= LOOP_FRACTION:
        return "LOOP"
    if n_H > n_E:
        return "HELIX"
    if n_E > n_H:
        return "SHEET"
    return "MIXED"


# ══════════════════════════════════════════════════════════════════════
# HELPER UTILITIES
# ══════════════════════════════════════════════════════════════════════

def mutate_repack(starting_pose, posi, amino, scorefxn):
    """
    Introduce a point mutation at a given Pose position and repack neighbours.

    The function clones the input pose so that the original is never modified.
    The mutation is applied via a design-restricted TaskFactory that:

    1. Allows rotamer sampling only within the neighbourhood of the target
       residue (atoms within the shell defined by
       ``NeighborhoodResidueSelector``).
    2. Restricts the target residue to the single canonical amino acid
       specified by *amino*.
    3. Restricts all other neighbourhood residues to repacking only
       (no sequence change).
    4. Freezes residues outside the neighbourhood entirely.

    Disulfide bonds are preserved (``NoRepackDisulfides``).

    Parameters
    ----------
    starting_pose : pyrosetta.Pose
        Template pose. Not modified.
    posi : int
        Rosetta Pose index (1-based) of the residue to mutate.
    amino : str
        One-letter code of the target amino acid (e.g. ``'A'`` for alanine).
    scorefxn : pyrosetta.ScoreFunction
        Score function used for rotamer packing.

    Returns
    -------
    pyrosetta.Pose
        A new pose carrying the requested mutation with repacked neighbours.
    """
    pose = starting_pose.clone()

    # --- Residue selectors ---------------------------------------------------
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)

    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)

    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)

    # --- Task factory --------------------------------------------------------
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Freeze residues outside the neighbourhood
    prevent_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        prevent_rlt, nbr_selector, True))

    # Repack-only for non-target residues inside the neighbourhood
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),
        not_design))

    # Restrict the target position to the requested amino acid
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        aa_to_design, mut_posi))

    # --- Packing -------------------------------------------------------------
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)

    return pose
        
def model_sequence(pose, mutations, scorefxn, relax=False):
    """
    Apply a set of point mutations to a pose and optionally relax the result.

    Iterates over the *mutations* dictionary and calls :func:`mutate_repack`
    sequentially for each position. If *relax* is ``True``, a Cartesian
    FastRelax pass (via :func:`pack_relax`) is applied to the fully mutated
    pose before returning, allowing the backbone and side chains to accommodate
    all changes simultaneously.

    Parameters
    ----------
    pose : pyrosetta.Pose
        Template pose. Not modified — a clone is used internally by
        :func:`mutate_repack` at each step.
    mutations : dict[int, str]
        Dictionary mapping Pose index (int, 1-based) to the desired one-letter
        amino acid code. Typically produced by :func:`Compare_sequences`.
        Pass an empty dict to skip mutations and only relax (if *relax=True*).
    scorefxn : pyrosetta.ScoreFunction
        Score function used for rotamer packing and relaxation.
    relax : bool, optional
        If ``True`` (default), run :func:`pack_relax` on the final mutated pose.
        Set to ``False`` to skip relaxation when speed is more important than
        structural quality (e.g. during preliminary screening).

    Returns
    -------
    pyrosetta.Pose
        New pose carrying all requested mutations, optionally relaxed.
    """
    new_pose = pose.clone()
    for index, target_aa in mutations.items():
        new_pose = mutate_repack(
            starting_pose=new_pose,
            posi=index,
            amino=target_aa,
            scorefxn=scorefxn,
        )

    return new_pose
def _mutate_to_epitope_seq(pose, insert_start, insert_end, epitope_seq):
    """
    Mutates residues [insert_start, insert_end] in pose to epitope_seq
    using the pack_mutate and model_sequence logic for better structural integrity.
    
    Args:
        pose (pyrosetta.Pose): The scaffold pose to modify.
        insert_start (int): Start index (1-based PyRosetta indexing).
        insert_end (int): End index (1-based PyRosetta indexing).
        epitope_seq (str): The target amino acid sequence.
    """
    # Ensure the sequence length matches the target range
    range_length = (insert_end - insert_start) + 1
    if len(epitope_seq) != range_length:
        print(f"Warning: Sequence length ({len(epitope_seq)}) does not match range ({range_length})")
    
    # Map the epitope sequence to specific residue indices
    # This prepares the data format expected by your model_sequence module
    mutations_to_apply = []
    for offset, aa in enumerate(epitope_seq):
        res_idx = insert_start + offset
        
        # Safety check for pose bounds
        if res_idx > pose.total_residue():
            break
            
        mutations_to_apply.append((res_idx, aa))
    
    # Use your model_sequence module to apply mutations with repacking
    try:
        # model_sequence handles the iterative mutate_repack calls
        # and preserves structural constraints like disulfides
        model_sequence(pose, mutations_to_apply)
    except Exception as e:
        print(f"Error during epitope modeling at residues {insert_start}-{insert_end}: {e}")

def _build_movemap(pose, loop_start, loop_end, flank=FLANK_RESIDUES):
    """
    Returns a MoveMap restricted to the loop region ± flank residues.
    """
    from pyrosetta.rosetta.core.kinematics import MoveMap
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    n_res = pose.total_residue()
    lo = max(1, loop_start - flank)
    hi = min(n_res, loop_end + flank)
    for i in range(lo, hi + 1):
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    return mm


# ══════════════════════════════════════════════════════════════════════
# REPLACEMENT LOGIC
# ══════════════════════════════════════════════════════════════════════

def _replace_region_with_epitope(scaffold_pose, epitope_pose,
                                 loop_start, loop_end, ep_ss_class):
    """
    Determines if it should mutate the sequence or perform a full graft.
    
    If ep_ss_class == "LOOP" and ep_len <= scaffold_loop_len:
        - Mutates the sequence of the scaffold loop.
        - Deletes remaining scaffold residues if epitope is smaller.
    Else:
        - Performs full deletion and insertion (grafting).
    """
    loop_len = loop_end - loop_start + 1
    ep_len   = epitope_pose.total_residue()
    ep_seq   = "".join(epitope_pose.residue(i).name1() for i in range(1, ep_len + 1))

    pose = scaffold_pose.clone()

    # MODIFICATION: Check if we just mutate sequence for loops
    if ep_ss_class == "LOOP" and ep_len <= loop_len:
        # Case 1: Epitope fits inside the existing loop scaffold
        # Mutate the segment starting from loop_start
        actual_end = loop_start + ep_len - 1
        _mutate_to_epitope_seq(pose, loop_start, actual_end, ep_seq)
        
        # If the epitope is strictly smaller, delete the leftover loop residues
        if ep_len < loop_len:
            delete_region(pose, actual_end + 1, loop_end)
            
        return pose, loop_start, actual_end

    else:
        # Case 2: Epitope is a structured region or longer than the loop
        # Full grafting: delete old region, insert new pose
        delete_region(pose, loop_start, loop_end)
        insert_pose_into_pose(pose, epitope_pose, loop_start - 1)
        actual_end = loop_start + ep_len - 1
        return pose, loop_start, actual_end


# ══════════════════════════════════════════════════════════════════════
# MAIN GRAFTING PIPELINE
# ══════════════════════════════════════════════════════════════════════

def graft_loop_with_modeling(scaffold_pose, epitope_pose,
                              loop_start, loop_end,
                              flank=FLANK_RESIDUES,
                              ccd_cycles=CCD_CYCLES,
                              relax_rounds=RELAX_ROUNDS):
    """
    Full grafting pipeline with selective modeling.

    1. Classifies epitope Secondary Structure.
    2. Replaces sequence (if LOOP/shorter) or pose (otherwise).
    3. Runs CCD to close junctions.
    4. Runs FastRelax for clash resolution.
    """
    # Identify Epitope properties
    ep_ss_class = classify_epitope_ss(epitope_pose)
    ep_len  = epitope_pose.total_residue()
    ep_seq  = "".join(epitope_pose.residue(i).name1() for i in range(1, ep_len + 1))

    # Step 1: Replace/Mutate region
    new_pose, insert_start, insert_end = _replace_region_with_epitope(
        scaffold_pose, epitope_pose, loop_start, loop_end, ep_ss_class
    )

    # Step 2: CCD loop closure
    # If it was a sequence mutation and the backbone didn't change, 
    # CCD typically exits quickly. Necessary if trimming/insertion occurred.
    from pyrosetta.rosetta.protocols.loops import Loop, Loops
    from pyrosetta.rosetta.protocols.loops.loop_closure.ccd import CCDLoopClosureMover
    
    try:
        cut = (insert_start + insert_end) // 2
        loop_obj = Loop(insert_start, insert_end, cut)
        mm = _build_movemap(new_pose, insert_start, insert_end, flank)
        ccd = CCDLoopClosureMover(loop_obj, mm)
        ccd.max_cycles(ccd_cycles)
        ccd.apply(new_pose)
    except Exception:
        pass

    # Step 3: FastRelax
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta import get_fa_scorefxn

    try:
        sfxn = get_fa_scorefxn()
        mm_relax = _build_movemap(new_pose, insert_start, insert_end, flank)
        fr = FastRelax(sfxn, relax_rounds)
        fr.set_movemap(mm_relax)
        fr.apply(new_pose)
    except Exception:
        pass

    meta = {
        "loop_start":            loop_start,
        "loop_end":              loop_end,
        "insert_start":          insert_start,
        "insert_end":            insert_end,
        "epitope_seq":           ep_seq,
        "ep_ss_class":           ep_ss_class,
        "scaffold_range_start":  insert_start,
        "scaffold_range_end":    insert_end,
    }

    return new_pose, meta
