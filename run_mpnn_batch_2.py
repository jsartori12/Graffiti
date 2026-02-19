import os
import argparse
import glob
import subprocess
import sys
import json
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import warnings

# Suppress BioPython warnings
warnings.simplefilter('ignore')

# Dictionary of CDRs
CDRs_DICTIONARY = {
    "H1": "KLDII",
    "H2": "IKFTIKFQEFSPNLWGLEFQKNKDYYII",
    "H3": "KILMKV"
}

# --- Helper Functions ---

def get_chain_sequence_and_indices(chain):
    """
    Extracts sequence from a Bio.PDB chain.
    """
    # Filter only amino acid residues (ignore waters/hetero)
    residues = [res for res in chain.get_residues() if res.get_id()[0] == ' ']
    sequence = seq1(''.join(res.resname for res in residues))
    return sequence

def get_formatted_fixed_indices(scaffold_sequence, cdr_dict):
    """
    Finds CDR indices in the sequence and returns 1-based list.
    """
    all_fixed_indices = []
    
    for cdr_name, cdr_sequence in cdr_dict.items():
        start_index = scaffold_sequence.find(cdr_sequence)
        
        if start_index != -1:
            # start_index is 0-based.
            end_position = start_index + len(cdr_sequence)
            # Convert to 1-based for ProteinMPNN (range excludes end, so +1 is implicit in logic)
            # We want indices: start+1, start+2... end
            indices_list = [i + 1 for i in range(start_index, end_position)]
            all_fixed_indices.extend(indices_list)
            print(f"      -> Found {cdr_name} at pos {start_index+1}-{end_position}")
        else:
            print(f"      -> WARNING: CDR '{cdr_name}' not found in sequence.")

    return sorted(list(set(all_fixed_indices)))

def create_fixed_positions_jsonl(pdb_path, output_dir, target_chain_id):
    """
    Generates JSONL defining rules:
    - Target Chain (e.g., B): Fix ONLY the CDRs.
    - Other Chains (e.g., A): Fix EVERYTHING (protect the antigen).
    """
    pdb_filename = os.path.basename(pdb_path)
    pdb_name_base = os.path.splitext(pdb_filename)[0]
    
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_name_base, pdb_path)
        model = structure[0]
    except Exception as e:
        print(f"   -> ERROR reading structure {pdb_path}: {e}")
        return None

    chain_definitions = {}
    
    # Check if target chain exists in structure
    chain_ids_in_structure = [c.get_id() for c in model]
    if target_chain_id not in chain_ids_in_structure:
        print(f"   -> ERROR: Target chain '{target_chain_id}' not found in PDB. Available chains: {chain_ids_in_structure}")
        return None

    for chain in model:
        chain_id = chain.get_id()
        seq = get_chain_sequence_and_indices(chain)
        
        if not seq:
            continue

        if chain_id == target_chain_id:
            # === TARGET CHAIN (Binder) ===
            # Fix only CDRs
            fixed_indices = get_formatted_fixed_indices(seq, CDRs_DICTIONARY)
            chain_definitions[chain_id] = fixed_indices
            print(f"      Chain {chain_id} (Target): Fixed {len(fixed_indices)} residues (CDRs).")
        else:
            # === CONTEXT CHAINS (Antigen) ===
            # Fix everything
            fixed_indices = list(range(1, len(seq) + 1))
            chain_definitions[chain_id] = fixed_indices
            print(f"      Chain {chain_id} (Context): Fixed completely ({len(seq)} residues).")

    final_output = {
        pdb_name_base: chain_definitions
    }

    jsonl_temp_dir = os.path.join(output_dir, "fixed_pos_jsonl_temp")
    os.makedirs(jsonl_temp_dir, exist_ok=True)
    jsonl_path = os.path.join(jsonl_temp_dir, f"{pdb_name_base}.jsonl")
    
    try:
        with open(jsonl_path, 'w') as f:
            f.write(json.dumps(final_output) + '\n')
        return jsonl_path
    except Exception as e:
        print(f"   -> ERROR saving JSONL: {e}")
        return None

# --- Main Execution ---
def run_proteinmpnn_batch(input_dir, output_dir, mpnn_script_path, num_sequences, model_name, sampling_temp, chain_id):
    
    os.makedirs(output_dir, exist_ok=True)
    pdb_files = glob.glob(os.path.join(input_dir, "*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files found in: {input_dir}")
        return

    print(f"Starting processing of {len(pdb_files)} PDB files...")
    print(f"Target Chain for Design: {chain_id}")

    for i, pdb_path in enumerate(pdb_files):
        pdb_filename = os.path.basename(pdb_path)
        print(f"\n[{i+1}/{len(pdb_files)}] Processing {pdb_filename}...")
        
        # 1. Generate JSONL
        fixed_pos_jsonl_path = create_fixed_positions_jsonl(pdb_path, output_dir, chain_id)
        
        if not fixed_pos_jsonl_path:
            print("   -> Skipping: JSONL generation failed.")
            continue

        # 2. Build Command
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
            "--batch_size=1" 
        ]
        
        if fixed_pos_jsonl_path:
            command.append(f"--fixed_positions_jsonl={fixed_pos_jsonl_path}")
        
        # 3. Execute
        try:
            # We use capture_output=False so you can see ProteinMPNN progress bars in real time
            subprocess.run(command, check=True, capture_output=False) 
            print(f"   -> Success!")
        except subprocess.CalledProcessError as e:
            print(f"   -> ERROR running ProteinMPNN.")
        except FileNotFoundError:
            print(f"   -> CRITICAL ERROR: Script not found at {mpnn_script_path}")
            break

    print("\nBatch processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run ProteinMPNN fixing CDRs on a specific chain.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("input_dir", type=str, help="Input PDB folder")
    parser.add_argument("output_dir", type=str, help="Output folder")
    parser.add_argument("--mpnn_script", type=str, default="/home/joao/Downloads/RFdiffusion/ProteinMPNN/protein_mpnn_run.py", help="Path to protein_mpnn_run.py")
    parser.add_argument("--num_seq", type=int, default=10, help="Sequences per target")
    parser.add_argument("--model", type=str, default="v_48_020", help="Model name")
    parser.add_argument("--temp", type=str, default="0.2", help="Sampling temperature")
    
    # NEW ARGUMENT: Chain ID (Defaults to B based on your PDB)
    parser.add_argument("--chain", type=str, default="B", help="Target Chain ID to design (default: B)")

    args = parser.parse_args()

    run_proteinmpnn_batch(
        args.input_dir, 
        args.output_dir, 
        args.mpnn_script, 
        args.num_seq, 
        args.model, 
        args.temp,
        args.chain # Pass the chain argument
    )