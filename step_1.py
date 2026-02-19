#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 14:55:46 2026

@author: joao
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASE 1 — Grafting independente de cada motif em cada scaffold.

Para cada par (scaffold, motif):
  - Tenta o graft isoladamente
  - Se bem-sucedido, salva os ranges de resíduos ocupados (start, end)
  - Registra tudo em um CSV para análise posterior na Fase 2

Saída:
  - graft_map_TIMESTAMP.csv   → tabela com ranges por par (scaffold, motif)
  - PDBs grafted individuais  → Grafts_individual/scaffold_X__motif_Y.pdb
"""

from pyrosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.select.residue_selector import ResiduePDBInfoHasLabelSelector
import glob
import os
import csv
from datetime import datetime

init()


# ─────────────────────────────────────────────
# XML
# ─────────────────────────────────────────────
def xmlobject_MotifGrafting(context_path, frag_path):
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
        <PROTOCOLS>
            <Add mover="motif_grafting"/>
        </PROTOCOLS>
    </ROSETTASCRIPTS>
    """


# ─────────────────────────────────────────────
# UTILITÁRIOS DE POSE
# ─────────────────────────────────────────────
def get_residues_by_label(pose, target_labels):
    result = set()
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        return result
    for i in range(1, pose.total_residue() + 1):
        for label in pdb_info.get_reslabels(i):
            if label in target_labels:
                result.add(i)
                break
    return result


def extract_motif_ranges(pose):
    """
    Extrai os ranges de resíduos ocupados pelo motif grafted.

    Retorna um dict com:
      - rosetta_indices : set de índices internos (1-based)
      - pdb_resnums     : lista de números PDB dos resíduos MOTIF
      - chains          : lista de chains dos resíduos MOTIF
      - start_resnum    : primeiro resnum PDB do motif
      - end_resnum      : último resnum PDB do motif
      - start_chain     : chain do primeiro resíduo
      - end_chain       : chain do último resíduo
      - motif_size      : número de resíduos MOTIF
      - connection_resnums : resnums dos resíduos CONNECTION (junções)
      - sequence        : sequência de aminoácidos inserida
      - scaffold_range_start : índice Rosetta do primeiro resíduo MOTIF no scaffold
      - scaffold_range_end   : índice Rosetta do último resíduo MOTIF no scaffold
    """
    pdb_info = pose.pdb_info()
    motif_idx   = sorted(get_residues_by_label(pose, {"MOTIF"}))
    conn_idx    = sorted(get_residues_by_label(pose, {"CONNECTION"}))

    if not motif_idx:
        return None

    pdb_resnums = [pdb_info.number(i) for i in motif_idx]
    chains      = [pdb_info.chain(i)  for i in motif_idx]
    sequence    = "".join(pose.residue(i).name1() for i in motif_idx)

    conn_resnums = [pdb_info.number(i) for i in conn_idx]

    return {
        "rosetta_indices":      motif_idx,
        "pdb_resnums":          pdb_resnums,
        "chains":               chains,
        "start_resnum":         pdb_resnums[0],
        "end_resnum":           pdb_resnums[-1],
        "start_chain":          chains[0],
        "end_chain":            chains[-1],
        "motif_size":           len(motif_idx),
        "connection_resnums":   conn_resnums,
        "sequence":             sequence,
        "scaffold_range_start": motif_idx[0],
        "scaffold_range_end":   motif_idx[-1],
    }


def remove_context(pose):
    """Remove resíduos CONTEXT da pose."""
    context_idx = get_residues_by_label(pose, {"CONTEXT"})
    clean = pose.clone()
    for i in range(clean.total_residue(), 0, -1):
        if i in context_idx:
            clean.conformation().delete_residue_slow(i)
    return clean


# ─────────────────────────────────────────────
# CSV
# ─────────────────────────────────────────────
CSV_COLUMNS = [
    "scaffold_id",
    "motif_name",
    "status",
    # Ranges no scaffold (índices Rosetta — comparáveis entre motifs do mesmo scaffold)
    "scaffold_range_start",
    "scaffold_range_end",
    # Ranges no PDB (numeração do arquivo de saída)
    "pdb_chain",
    "pdb_start_resnum",
    "pdb_end_resnum",
    # Detalhes
    "motif_size",
    "sequence",
    "connection_resnums",
    "scaffold_total_residues",
    "error_message",
]


def init_csv(path):
    with open(path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=CSV_COLUMNS).writeheader()


def append_csv(path, row):
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writerow({col: row.get(col, "") for col in CSV_COLUMNS})


def make_success_record(scaffold_id, motif_name, ranges, scaffold_total):
    conn_str = ";".join(str(r) for r in ranges["connection_resnums"])
    return {
        "scaffold_id":            scaffold_id,
        "motif_name":             motif_name,
        "status":                 "SUCCESS",
        "scaffold_range_start":   ranges["scaffold_range_start"],
        "scaffold_range_end":     ranges["scaffold_range_end"],
        "pdb_chain":              ranges["start_chain"],
        "pdb_start_resnum":       ranges["start_resnum"],
        "pdb_end_resnum":         ranges["end_resnum"],
        "motif_size":             ranges["motif_size"],
        "sequence":               ranges["sequence"],
        "connection_resnums":     conn_str,
        "scaffold_total_residues": scaffold_total,
        "error_message":          "",
    }


def make_failure_record(scaffold_id, motif_name, error_msg, scaffold_total=""):
    return {
        "scaffold_id":            scaffold_id,
        "motif_name":             motif_name,
        "status":                 "FAILED",
        "scaffold_range_start":   "",
        "scaffold_range_end":     "",
        "pdb_chain":              "",
        "pdb_start_resnum":       "",
        "pdb_end_resnum":         "",
        "motif_size":             "",
        "sequence":               "",
        "connection_resnums":     "",
        "scaffold_total_residues": scaffold_total,
        "error_message":          str(error_msg).strip().splitlines()[0],
    }


# ─────────────────────────────────────────────
# CONFIGURAÇÕES
# ─────────────────────────────────────────────
scaffold_dir  = "relaxed_results/"
epitopes_dir  = "Epitopes/"
output_dir    = "Grafts_individual/"
os.makedirs(output_dir, exist_ok=True)

pdb_files   = sorted(glob.glob(os.path.join(scaffold_dir, "*.pdb")))
motifs_list = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
csv_path  = os.path.join(output_dir, f"graft_map_{timestamp}.csv")
init_csv(csv_path)

total_pairs    = len(pdb_files) * len(motifs_list)
pair_count     = 0
success_count  = 0

print(f"Fase 1 — Grafting independente")
print(f"Scaffolds : {len(pdb_files)}")
print(f"Motifs    : {len(motifs_list)}")
print(f"Total pares: {total_pairs}")
print(f"CSV: {csv_path}\n")


# ─────────────────────────────────────────────
# LOOP PRINCIPAL
# ─────────────────────────────────────────────
for pdb_path in pdb_files:
    scaffold_id   = os.path.splitext(os.path.basename(pdb_path))[0]
    scaffold_size = pose_from_pdb(pdb_path).total_residue()

    print(f"\n{'='*55}")
    print(f"Scaffold: {scaffold_id}  ({scaffold_size} resíduos)")
    print(f"{'='*55}")

    for motif_path in motifs_list:
        motif_name = os.path.basename(motif_path)
        pair_count += 1
        progress   = f"[{pair_count}/{total_pairs}]"

        try:
            print(f"  {progress} → {motif_name}", end=" ... ", flush=True)

            # Carrega scaffold fresco a cada tentativa (grafts independentes)
            pose = pose_from_pdb(pdb_path)

            objs  = XmlObjects.create_from_string(
                        xmlobject_MotifGrafting(pdb_path, motif_path))
            mover = objs.get_mover("motif_grafting")
            mover.apply(pose)

            # Extrair ranges ANTES de remover contexto
            ranges = extract_motif_ranges(pose)

            if ranges is None:
                # Graft aplicado mas nenhum resíduo MOTIF encontrado (caso raro)
                raise RuntimeError("Graft applied but no MOTIF residues labeled")

            # Salvar PDB individual (sem contexto)
            clean_pose  = remove_context(pose)
            pdb_out     = os.path.join(output_dir, f"{scaffold_id}__{motif_name}")
            clean_pose.dump_pdb(pdb_out)

            # Registrar no CSV
            record = make_success_record(scaffold_id, motif_name, ranges, scaffold_size)
            append_csv(csv_path, record)
            success_count += 1

            print(f"OK  | range [{ranges['scaffold_range_start']}–{ranges['scaffold_range_end']}]"
                  f"  seq: {ranges['sequence']}")

        except RuntimeError as e:
            msg = str(e)
            label = "sem região compatível" if "not suitable" in msg else msg.splitlines()[0]
            print(f"FAIL | {label}")
            append_csv(csv_path, make_failure_record(scaffold_id, motif_name, msg, scaffold_size))

        except Exception as e:
            print(f"ERROR | {e}")
            append_csv(csv_path, make_failure_record(scaffold_id, motif_name, str(e), scaffold_size))


# ─────────────────────────────────────────────
# RESUMO
# ─────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"Fase 1 concluída.")
print(f"  Pares testados  : {total_pairs}")
print(f"  Grafts aceitos  : {success_count}")
print(f"  Grafts rejeitados: {total_pairs - success_count}")
print(f"  Taxa de sucesso : {success_count / total_pairs * 100:.1f}%")
print(f"  Mapa salvo em   : {csv_path}")
print(f"{'='*55}")