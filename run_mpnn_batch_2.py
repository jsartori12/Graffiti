#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ProteinMPNN batch design with automatic epitope locking.

Instead of hardcoding CDR sequences, this script reads the epitope PDB files
directly and builds the sequence dictionary automatically. During design, all
grafted epitope regions are fixed (locked), and only the scaffold regions are
allowed to mutate.

Usage:
    python proteinmpnn_design.py input_dir output_dir --epitopes_dir Epitopes/

    # With all options:
    python proteinmpnn_design.py Grafts_optimized/final_pdbs/ MPNN_out/ \\
        --epitopes_dir Epitopes/ \\
        --chain A \\
        --num_seq 10 \\
        --temp 0.2 \\
        --model v_48_020 \\
        --mpnn_script /path/to/protein_mpnn_run.py
"""

import os
import argparse
import glob
import subprocess
import sys
import json
import warnings
from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1

warnings.simplefilter("ignore")


# ─────────────────────────────────────────────
# BUILD EPITOPE DICTIONARY FROM PDB FILES
# ─────────────────────────────────────────────
def extract_sequence_from_pdb(pdb_path):
    """
    Reads a PDB file and returns the full amino acid sequence
    of the first chain found (1-letter codes).
    Uses PPBuilder to handle multi-residue chains correctly.
    """
    parser  = PDBParser(QUIET=True)
    name    = os.path.splitext(os.path.basename(pdb_path))[0]
    structure = parser.get_structure(name, pdb_path)
    model   = structure[0]

    ppb = PPBuilder()
    sequences = []
    for pp in ppb.build_peptides(model):
        sequences.append(str(pp.get_sequence()))

    if not sequences:
        # Fallback: manual extraction if PPBuilder finds nothing
        for chain in model:
            residues = [r for r in chain.get_residues() if r.get_id()[0] == " "]
            seq = seq1("".join(r.resname for r in residues))
            if seq:
                sequences.append(seq)
            break  # only first chain

    return "".join(sequences) if sequences else None


def build_epitope_dictionary(epitopes_dir):
    """
    Reads all PDB files in epitopes_dir and builds a dictionary:
        { epitope_name : sequence }

    This replaces the hardcoded CDRs_DICTIONARY.
    The epitope name is the filename without extension.

    Returns the dictionary and prints a summary.
    """
    pdb_files = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))

    if not pdb_files:
        raise FileNotFoundError(
            f"No PDB files found in epitopes directory: {epitopes_dir}"
        )

    epitope_dict = {}
    print(f"\nBuilding epitope dictionary from {len(pdb_files)} PDB files...")

    for pdb_path in pdb_files:
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        seq  = extract_sequence_from_pdb(pdb_path)

        if seq:
            epitope_dict[name] = seq
            print(f"  [{name}]  {seq}  ({len(seq)} aa)")
        else:
            print(f"  WARNING: Could not extract sequence from {pdb_path} — skipping.")

    print(f"  → {len(epitope_dict)} epitopes loaded.\n")
    return epitope_dict


# ─────────────────────────────────────────────
# SEQUENCE & INDEX UTILITIES
# ─────────────────────────────────────────────
def get_chain_sequence(chain):
    """Extracts the 1-letter amino acid sequence from a Bio.PDB chain."""
    residues = [r for r in chain.get_residues() if r.get_id()[0] == " "]
    return seq1("".join(r.resname for r in residues))


def find_epitope_indices(scaffold_sequence, epitope_dict):
    """
    Searches for each epitope sequence as a substring of the scaffold sequence.
    Returns a sorted list of 1-based residue indices that belong to any epitope.

    If an epitope is not found, a warning is printed but execution continues —
    it may mean that epitope was not grafted into this particular scaffold.
    """
    all_fixed = []

    for epitope_name, epitope_seq in epitope_dict.items():
        start = scaffold_sequence.find(epitope_seq)

        if start != -1:
            end     = start + len(epitope_seq)
            indices = [i + 1 for i in range(start, end)]  # convert to 1-based
            all_fixed.extend(indices)
            print(f"      -> Locked [{epitope_name}] at positions {start+1}–{end} "
                  f"({len(indices)} residues)")
        else:
            # Not a hard error: this epitope may not be present in this scaffold
            print(f"      -> [{epitope_name}] not found in chain — skipping.")

    fixed = sorted(set(all_fixed))
    return fixed


# ─────────────────────────────────────────────
# JSONL GENERATION
# ─────────────────────────────────────────────
def create_fixed_positions_jsonl(pdb_path, output_dir, target_chain_id, epitope_dict):
    """
    Generates a JSONL file defining which residues are fixed for ProteinMPNN:

    - Target chain (e.g., A): fix ONLY the grafted epitope regions.
      Everything else on this chain is free to be redesigned.
    - Other chains (e.g., context/antigen): fix ENTIRELY to protect them.

    Returns the path to the generated JSONL, or None on failure.
    """
    pdb_filename  = os.path.basename(pdb_path)
    pdb_name_base = os.path.splitext(pdb_filename)[0]

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_name_base, pdb_path)
        model     = structure[0]
    except Exception as e:
        print(f"   -> ERROR reading {pdb_path}: {e}")
        return None

    chain_ids = [c.get_id() for c in model]
    if target_chain_id not in chain_ids:
        print(f"   -> ERROR: chain '{target_chain_id}' not found. "
              f"Available: {chain_ids}")
        return None

    chain_definitions = {}

    for chain in model:
        chain_id = chain.get_id()
        seq      = get_chain_sequence(chain)

        if not seq:
            continue

        if chain_id == target_chain_id:
            # Design chain: lock only the epitope residues
            fixed_indices = find_epitope_indices(seq, epitope_dict)

            if not fixed_indices:
                print(f"      Chain {chain_id} (Design): No epitopes found — "
                      f"all {len(seq)} residues are free to mutate.")
            else:
                free = len(seq) - len(fixed_indices)
                print(f"      Chain {chain_id} (Design): "
                      f"{len(fixed_indices)} residues locked (epitopes), "
                      f"{free} free to mutate.")

            chain_definitions[chain_id] = fixed_indices

        else:
            # Context chain: lock entirely
            fixed_indices = list(range(1, len(seq) + 1))
            chain_definitions[chain_id] = fixed_indices
            print(f"      Chain {chain_id} (Context): fully locked ({len(seq)} residues).")

    final_output = {pdb_name_base: chain_definitions}

    jsonl_dir  = os.path.join(output_dir, "fixed_pos_jsonl_temp")
    os.makedirs(jsonl_dir, exist_ok=True)
    jsonl_path = os.path.join(jsonl_dir, f"{pdb_name_base}.jsonl")

    try:
        with open(jsonl_path, "w") as f:
            f.write(json.dumps(final_output) + "\n")
        return jsonl_path
    except Exception as e:
        print(f"   -> ERROR writing JSONL: {e}")
        return None


# ─────────────────────────────────────────────
# MAIN BATCH RUNNER
# ─────────────────────────────────────────────
def run_proteinmpnn_batch(
    input_dir, output_dir, epitopes_dir,
    mpnn_script_path, num_sequences,
    model_name, sampling_temp, chain_id
):
    os.makedirs(output_dir, exist_ok=True)

    # Build epitope dictionary once — reused for all scaffolds
    epitope_dict = build_epitope_dictionary(epitopes_dir)

    pdb_files = sorted(glob.glob(os.path.join(input_dir, "*.pdb")))
    if not pdb_files:
        print(f"No PDB files found in: {input_dir}")
        return

    print(f"Processing {len(pdb_files)} PDB files...")
    print(f"Target chain for design : {chain_id}")
    print(f"Epitopes to lock        : {list(epitope_dict.keys())}\n")

    success = 0
    failed  = 0

    for i, pdb_path in enumerate(pdb_files):
        pdb_filename = os.path.basename(pdb_path)
        print(f"\n[{i+1}/{len(pdb_files)}] {pdb_filename}")

        # 1. Generate JSONL with fixed positions
        jsonl_path = create_fixed_positions_jsonl(
            pdb_path, output_dir, chain_id, epitope_dict
        )

        if not jsonl_path:
            print("   -> Skipping: JSONL generation failed.")
            failed += 1
            continue

        # 2. Build ProteinMPNN command
        command = [
            sys.executable,
            mpnn_script_path,
            f"--pdb_path={pdb_path}",
            f"--out_folder={output_dir}",
            f"--num_seq_per_target={num_sequences}",
            f"--model_name={model_name}",
            f"--sampling_temp={sampling_temp}",
            "--seed=0",
            "--save_score=1",
            "--batch_size=1",
            f"--fixed_positions_jsonl={jsonl_path}",
        ]

        # 3. Run
        try:
            subprocess.run(command, check=True, capture_output=False)
            print(f"   -> Done.")
            success += 1
        except subprocess.CalledProcessError:
            print(f"   -> ERROR: ProteinMPNN failed for {pdb_filename}.")
            failed += 1
        except FileNotFoundError:
            print(f"   -> CRITICAL: Script not found at {mpnn_script_path}")
            break

    print(f"\n{'='*50}")
    print(f"Batch complete.  Success: {success}  |  Failed: {failed}")
    print(f"{'='*50}")


# ─────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_dir",  help="Folder with grafted PDB files to design")
    parser.add_argument("output_dir", help="Output folder for ProteinMPNN results")

    parser.add_argument(
        "--epitopes_dir", required=True,
        help="Folder containing the original epitope PDB files.\n"
             "Sequences are extracted automatically and used to lock those\n"
             "regions during design."
    )
    parser.add_argument(
        "--mpnn_script",
        default="/home/joao/Downloads/RFdiffusion/ProteinMPNN/protein_mpnn_run.py",
        help="Path to protein_mpnn_run.py"
    )
    parser.add_argument("--num_seq", type=int,   default=10,      help="Sequences per target")
    parser.add_argument("--model",   type=str,   default="v_48_020", help="ProteinMPNN model name")
    parser.add_argument("--temp",    type=str,   default="0.2",   help="Sampling temperature")
    parser.add_argument("--chain",   type=str,   default="A",     help="Chain ID to design (default: A)")

    args = parser.parse_args()

    run_proteinmpnn_batch(
        input_dir        = args.input_dir,
        output_dir       = args.output_dir,
        epitopes_dir     = args.epitopes_dir,
        mpnn_script_path = args.mpnn_script,
        num_sequences    = args.num_seq,
        model_name       = args.model,
        sampling_temp    = args.temp,
        chain_id         = args.chain,
    )
    
    
