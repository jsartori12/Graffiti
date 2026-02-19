#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASE 2 — Análise de sobreposição e seleção otimizada.

Lê o CSV gerado pela Fase 1 e para cada scaffold:
  1. Identifica quais pares de motifs têm ranges sobrepostos
  2. Resolve o conflito removendo o motif que causa mais sobreposições
     (estratégia greedy: maximiza o número de motifs compatíveis)
  3. Gera um relatório com:
     - Quais motifs foram selecionados por scaffold
     - Quais foram removidos e por qual conflito
     - Ranking de scaffolds por número de motifs compatíveis

Saídas:
  - overlap_report_TIMESTAMP.csv     → decisão por par (scaffold, motif)
  - scaffold_ranking_TIMESTAMP.csv   → ranking de scaffolds
  - overlap_charts_TIMESTAMP.png     → visualizações
"""

import os
import csv
import glob
from datetime import datetime
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np


# ─────────────────────────────────────────────
# CARREGAR CSV DA FASE 1
# ─────────────────────────────────────────────
def load_graft_map(csv_path):
    """
    Lê o CSV da Fase 1 e retorna apenas os registros SUCCESS
    como lista de dicts com campos convertidos para int onde necessário.
    """
    records = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            if row["status"] != "SUCCESS":
                continue
            row["scaffold_range_start"] = int(row["scaffold_range_start"])
            row["scaffold_range_end"]   = int(row["scaffold_range_end"])
            row["motif_size"]           = int(row["motif_size"])
            records.append(row)
    return records


# ─────────────────────────────────────────────
# DETECÇÃO DE SOBREPOSIÇÃO
# ─────────────────────────────────────────────
def ranges_overlap(start1, end1, start2, end2, buffer=0):
    """
    Verifica se dois ranges se sobrepõem.
    buffer: margem extra de segurança (resíduos) ao redor de cada range.
    """
    return not (end1 + buffer < start2 or end2 + buffer < start1)


def find_overlaps(motif_records, buffer=0):
    """
    Dado um conjunto de registros de motifs para um scaffold,
    retorna uma lista de pares (motif_A, motif_B) que se sobrepõem.
    """
    overlaps = []
    for i in range(len(motif_records)):
        for j in range(i + 1, len(motif_records)):
            a = motif_records[i]
            b = motif_records[j]
            if ranges_overlap(
                a["scaffold_range_start"], a["scaffold_range_end"],
                b["scaffold_range_start"], b["scaffold_range_end"],
                buffer=buffer
            ):
                overlaps.append((a["motif_name"], b["motif_name"]))
    return overlaps


# ─────────────────────────────────────────────
# SELEÇÃO GREEDY (maximiza motifs compatíveis)
# ─────────────────────────────────────────────
def greedy_select(motif_records, buffer=0):
    """
    Algoritmo greedy para selecionar o máximo de motifs sem sobreposição:

    1. Ordena os motifs pelo tamanho do range (menor primeiro)
       → favorece motifs menores que "ocupam menos espaço"
    2. Itera: aceita o próximo motif se não sobrepõe com nenhum já aceito
    3. Registra os rejeitados e com qual motif conflitaram

    Retorna:
      selected : lista de records aceitos
      rejected : lista de dicts {record, conflict_with}
    """
    # Ordenar por tamanho do range (menor = menos invasivo)
    sorted_records = sorted(
        motif_records,
        key=lambda r: r["scaffold_range_end"] - r["scaffold_range_start"]
    )

    selected = []
    rejected = []

    for candidate in sorted_records:
        conflict = None
        for accepted in selected:
            if ranges_overlap(
                candidate["scaffold_range_start"], candidate["scaffold_range_end"],
                accepted["scaffold_range_start"],  accepted["scaffold_range_end"],
                buffer=buffer
            ):
                conflict = accepted["motif_name"]
                break

        if conflict is None:
            selected.append(candidate)
        else:
            rejected.append({
                "record":        candidate,
                "conflict_with": conflict,
            })

    return selected, rejected


# ─────────────────────────────────────────────
# ANÁLISE PRINCIPAL
# ─────────────────────────────────────────────
def analyze(records, buffer=0):
    """
    Para cada scaffold, roda a seleção greedy e coleta estatísticas.

    Retorna:
      results : dict scaffold_id → {selected, rejected, all_motifs, overlaps}
    """
    # Agrupar por scaffold
    by_scaffold = defaultdict(list)
    for r in records:
        by_scaffold[r["scaffold_id"]].append(r)

    results = {}
    for scaffold_id, motif_records in by_scaffold.items():
        overlaps          = find_overlaps(motif_records, buffer=buffer)
        selected, rejected = greedy_select(motif_records, buffer=buffer)

        results[scaffold_id] = {
            "all_motifs": motif_records,
            "selected":   selected,
            "rejected":   rejected,
            "overlaps":   overlaps,
            "n_total":    len(motif_records),
            "n_selected": len(selected),
            "n_rejected": len(rejected),
            "n_overlapping_pairs": len(overlaps),
        }

    return results


# ─────────────────────────────────────────────
# SALVAR CSVs
# ─────────────────────────────────────────────
def save_decision_csv(results, path):
    """
    CSV detalhado: uma linha por (scaffold, motif) com a decisão.
    """
    cols = [
        "scaffold_id", "motif_name", "decision",
        "scaffold_range_start", "scaffold_range_end",
        "motif_size", "sequence",
        "conflict_with", "pdb_chain",
        "pdb_start_resnum", "pdb_end_resnum",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()

        for scaffold_id, data in sorted(results.items()):
            for r in data["selected"]:
                writer.writerow({
                    "scaffold_id":          scaffold_id,
                    "motif_name":           r["motif_name"],
                    "decision":             "KEEP",
                    "scaffold_range_start": r["scaffold_range_start"],
                    "scaffold_range_end":   r["scaffold_range_end"],
                    "motif_size":           r["motif_size"],
                    "sequence":             r["sequence"],
                    "conflict_with":        "",
                    "pdb_chain":            r["pdb_chain"],
                    "pdb_start_resnum":     r["pdb_start_resnum"],
                    "pdb_end_resnum":       r["pdb_end_resnum"],
                })
            for item in data["rejected"]:
                r = item["record"]
                writer.writerow({
                    "scaffold_id":          scaffold_id,
                    "motif_name":           r["motif_name"],
                    "decision":             "REMOVE",
                    "scaffold_range_start": r["scaffold_range_start"],
                    "scaffold_range_end":   r["scaffold_range_end"],
                    "motif_size":           r["motif_size"],
                    "sequence":             r["sequence"],
                    "conflict_with":        item["conflict_with"],
                    "pdb_chain":            r["pdb_chain"],
                    "pdb_start_resnum":     r["pdb_start_resnum"],
                    "pdb_end_resnum":       r["pdb_end_resnum"],
                })
    print(f"  📄 Decisões salvas: {path}")


def save_ranking_csv(results, path):
    """
    CSV de ranking: um scaffold por linha, ordenado por motifs selecionados (desc).
    """
    cols = [
        "rank", "scaffold_id",
        "n_motifs_grafted", "n_selected", "n_removed",
        "n_overlapping_pairs", "selected_motifs", "removed_motifs",
    ]
    rows = []
    for scaffold_id, data in results.items():
        rows.append({
            "scaffold_id":         scaffold_id,
            "n_motifs_grafted":    data["n_total"],
            "n_selected":          data["n_selected"],
            "n_removed":           data["n_rejected"],
            "n_overlapping_pairs": data["n_overlapping_pairs"],
            "selected_motifs":     ";".join(r["motif_name"] for r in data["selected"]),
            "removed_motifs":      ";".join(i["record"]["motif_name"] for i in data["rejected"]),
        })

    rows.sort(key=lambda x: (-x["n_selected"], -x["n_motifs_grafted"]))

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for rank, row in enumerate(rows, start=1):
            row["rank"] = rank
            writer.writerow(row)

    print(f"  📄 Ranking salvo: {path}")
    return rows


# ─────────────────────────────────────────────
# GRÁFICOS
# ─────────────────────────────────────────────
def save_charts(results, ranking_rows, timestamp, output_dir):
    # Limitar a 50 scaffolds por gráfico para não estourar o tamanho da imagem
    MAX_PER_CHART = 50

    # # ── Gráfico 1: Ranking de scaffolds (barras empilhadas KEEP + REMOVE) ──
    # chunks = [ranking_rows[i:i+MAX_PER_CHART]
    #           for i in range(0, len(ranking_rows), MAX_PER_CHART)]

    # for chunk_idx, chunk in enumerate(chunks):
    #     sc_names = [r["scaffold_id"] for r in chunk]
    #     sc_sel   = [r["n_selected"]  for r in chunk]
    #     sc_rem   = [r["n_removed"]   for r in chunk]
    #     n        = len(sc_names)

    #     # Horizontal se muitos scaffolds, vertical se poucos
    #     if n > 20:
    #         fig1, ax1 = plt.subplots(figsize=(8, max(6, n * 0.35)))
    #         y = np.arange(n)
    #         ax1.barh(y, sc_sel, label="Mantidos (KEEP)", color="#2ecc71", edgecolor="white")
    #         ax1.barh(y, sc_rem, left=sc_sel, label="Removidos (REMOVE)",
    #                  color="#e74c3c", edgecolor="white")
    #         for i, (sel, rem) in enumerate(zip(sc_sel, sc_rem)):
    #             total = sel + rem
    #             if total > 0:
    #                 ax1.text(total + 0.05, i, str(total), va="center", fontsize=7, color="gray")
    #             if sel > 0:
    #                 ax1.text(sel / 2, i, str(sel), ha="center", va="center",
    #                          fontsize=7, fontweight="bold", color="white")
    #         ax1.set_yticks(y)
    #         ax1.set_yticklabels(sc_names, fontsize=7)
    #         ax1.invert_yaxis()
    #         ax1.set_xlabel("Número de motifs", fontsize=11)
    #         ax1.grid(axis="x", alpha=0.3)
    #     else:
    #         fig1, ax1 = plt.subplots(figsize=(max(8, n * 0.9), 6))
    #         x = np.arange(n)
    #         ax1.bar(x, sc_sel, label="Mantidos (KEEP)", color="#2ecc71", edgecolor="white")
    #         ax1.bar(x, sc_rem, bottom=sc_sel, label="Removidos (REMOVE)",
    #                 color="#e74c3c", edgecolor="white")
    #         for i, (sel, rem) in enumerate(zip(sc_sel, sc_rem)):
    #             total = sel + rem
    #             if total > 0:
    #                 ax1.text(i, total + 0.1, str(total),
    #                          ha="center", va="bottom", fontsize=8, color="gray")
    #             if sel > 0:
    #                 ax1.text(i, sel / 2, str(sel), ha="center", va="center",
    #                          fontsize=8, fontweight="bold", color="white")
    #         ax1.set_xticks(x)
    #         ax1.set_xticklabels(sc_names, rotation=40, ha="right", fontsize=8)
    #         ax1.set_ylabel("Número de motifs", fontsize=11)
    #         ax1.grid(axis="y", alpha=0.3)

    #     suffix = f"_part{chunk_idx+1}" if len(chunks) > 1 else ""
    #     title  = "Ranking de Scaffolds — Motifs Selecionados vs Removidos"
    #     if len(chunks) > 1:
    #         title += f" ({chunk_idx+1}/{len(chunks)})"
    #     ax1.set_title(title, fontsize=12, fontweight="bold", pad=12)
    #     ax1.legend(fontsize=10)

    #     chart1 = os.path.join(output_dir, f"ranking_chart{suffix}_{timestamp}.png")
    #     plt.tight_layout()
    #     plt.savefig(chart1, dpi=150, bbox_inches="tight")
    #     plt.close()
    #     print(f"  📊 Gráfico ranking salvo: {chart1}")

    # ── Gráfico 2: Distribuição — quantos scaffolds ficaram com N motifs ──
    all_n_selected = [r["n_selected"] for r in ranking_rows]
    max_sel        = max(all_n_selected) if all_n_selected else 0
    dist           = {i: all_n_selected.count(i) for i in range(max_sel + 1)}

    fig2, ax2 = plt.subplots(figsize=(max(7, (max_sel + 1) * 1.1), 5))

    colors = []
    for k in dist:
        if k == 0:              colors.append("#bdc3c7")
        elif k == max_sel:      colors.append("#2ecc71")
        elif k >= max_sel * .5: colors.append("#27ae60")
        elif k >= max_sel * .25:colors.append("#e67e22")
        else:                   colors.append("#e74c3c")

    bars = ax2.bar(list(dist.keys()), list(dist.values()),
                   color=colors, edgecolor="white", linewidth=0.8, width=0.6)

    for bar, val in zip(bars, dist.values()):
        if val > 0:
            ax2.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 0.1, str(val),
                     ha="center", va="bottom", fontsize=11, fontweight="bold")

    ax2.set_xticks(list(dist.keys()))
    ax2.set_xlabel("Motifs selecionados (sem sobreposição)", fontsize=11)
    ax2.set_ylabel("Número de scaffolds", fontsize=11)
    ax2.set_title("Distribuição de Scaffolds por Motifs Compatíveis",
                  fontsize=13, fontweight="bold", pad=12)
    ax2.set_ylim(0, max(dist.values()) + 2)
    ax2.grid(axis="y", alpha=0.3)

    chart2 = os.path.join(output_dir, f"distribution_chart_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(chart2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📊 Gráfico distribuição salvo: {chart2}")

    # # ── Gráfico 3: Mapa de sobreposições por scaffold (heatmap motif × motif) ──
    # # Pega o scaffold com mais sobreposições como exemplo representativo
    # most_conflicts = max(results.items(), key=lambda x: x[1]["n_overlapping_pairs"])
    # sc_id, sc_data = most_conflicts

    # all_m = [r["motif_name"] for r in sc_data["all_motifs"]]
    # n_m   = len(all_m)

    # if n_m > 1:
    #     matrix = np.zeros((n_m, n_m), dtype=int)
    #     overlap_pairs = {frozenset(p) for p in sc_data["overlaps"]}
    #     for i, mi in enumerate(all_m):
    #         for j, mj in enumerate(all_m):
    #             if i != j and frozenset([mi, mj]) in overlap_pairs:
    #                 matrix[i, j] = 1

    #     short_names = [m.replace(".pdb", "").replace("epitope_", "ep") for m in all_m]

    #     fig3, ax3 = plt.subplots(figsize=(max(6, n_m * 0.8), max(5, n_m * 0.8)))
    #     cmap = mcolors.ListedColormap(["#ecf0f1", "#e74c3c"])
    #     ax3.imshow(matrix, cmap=cmap, vmin=0, vmax=1)

    #     for i in range(n_m):
    #         for j in range(n_m):
    #             if i == j:
    #                 ax3.text(j, i, "–", ha="center", va="center", fontsize=9, color="gray")
    #             elif matrix[i, j] == 1:
    #                 ax3.text(j, i, "✗", ha="center", va="center",
    #                          fontsize=10, fontweight="bold", color="white")

    #     ax3.set_xticks(range(n_m))
    #     ax3.set_yticks(range(n_m))
    #     ax3.set_xticklabels(short_names, rotation=40, ha="right", fontsize=8)
    #     ax3.set_yticklabels(short_names, fontsize=8)
    #     ax3.set_title(f"Mapa de Sobreposições — {sc_id}\n(scaffold com mais conflitos)",
    #                   fontsize=11, fontweight="bold", pad=10)

    #     patch_ov  = mpatches.Patch(color="#e74c3c", label="Sobreposição")
    #     patch_ok  = mpatches.Patch(color="#ecf0f1", label="Compatível")
    #     ax3.legend(handles=[patch_ov, patch_ok], fontsize=9, loc="upper right",
    #                bbox_to_anchor=(1.25, 1.0))

    #     chart3 = os.path.join(output_dir, f"overlap_heatmap_{timestamp}.png")
    #     plt.tight_layout()
    #     plt.savefig(chart3, dpi=150, bbox_inches="tight")
    #     plt.close()
    #     print(f"  📊 Heatmap sobreposição salvo: {chart3}")


# ─────────────────────────────────────────────
# CONFIGURAÇÕES
# ─────────────────────────────────────────────

# Aponta para o CSV gerado pela Fase 1
# Altere para o caminho correto ou use glob para pegar o mais recente
phase1_csvs = sorted(glob.glob("Grafts_individual/graft_map_*.csv"))
if not phase1_csvs:
    raise FileNotFoundError("Nenhum CSV da Fase 1 encontrado em Grafts_individual/")

phase1_csv = phase1_csvs[-1]  # usa o mais recente
print(f"Usando CSV da Fase 1: {phase1_csv}")

# Buffer de segurança: ignora sobreposições menores que N resíduos
# 0 = qualquer sobreposição conta; 2-3 = permite pequena adjacência
OVERLAP_BUFFER = 0

output_dir = "Grafts_optimized/"
os.makedirs(output_dir, exist_ok=True)
timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")


# ─────────────────────────────────────────────
# EXECUÇÃO
# ─────────────────────────────────────────────
print(f"\nFase 2 — Análise de sobreposição (buffer={OVERLAP_BUFFER} resíduos)")

records = load_graft_map(phase1_csv)
print(f"Registros SUCCESS carregados: {len(records)}")

results = analyze(records, buffer=OVERLAP_BUFFER)

decision_csv = os.path.join(output_dir, f"overlap_decisions_{timestamp}.csv")
ranking_csv  = os.path.join(output_dir, f"scaffold_ranking_{timestamp}.csv")

save_decision_csv(results, decision_csv)
ranking_rows = save_ranking_csv(results, ranking_csv)
save_charts(results, ranking_rows, timestamp, output_dir)


# ─────────────────────────────────────────────
# GERAÇÃO DOS PDBs FINAIS
# ─────────────────────────────────────────────
from pyrosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

init()

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

def remove_context(pose):
    context_idx = get_residues_by_label(pose, {"CONTEXT"})
    clean = pose.clone()
    for i in range(clean.total_residue(), 0, -1):
        if i in context_idx:
            clean.conformation().delete_residue_slow(i)
    return clean


# ─────────────────────────────────────────────
# SASA
# ─────────────────────────────────────────────
def calculate_sasa(pose):
    """
    Calcula o SASA total e por resíduo da pose.
    Retorna (total_sasa, rsd_sasa_vector).
    """
    calc       = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
    total_sasa = calc.calculate(pose)
    rsd_sasa   = calc.get_residue_sasa()
    return total_sasa, rsd_sasa


def build_sasa_records(pose, scaffold_id, motifs_applied, rsd_sasa, total_sasa):
    """
    Constrói uma lista de dicts — um por resíduo — com:
      - scaffold_id, motifs_applied
      - resíduo: índice, resname, resnum PDB, chain
      - sasa_value (Å²)
      - region: MOTIF / SCAFFOLD / CONNECTION / HOTSPOT
    """
    pdb_info = pose.pdb_info()
    motif_idx = get_residues_by_label(pose, {"MOTIF"})
    conn_idx  = get_residues_by_label(pose, {"CONNECTION"})
    hot_idx   = get_residues_by_label(pose, {"HOTSPOT"})

    records = []
    for i in range(1, pose.total_residue() + 1):
        if i in motif_idx:
            region = "MOTIF"
        elif i in conn_idx:
            region = "CONNECTION"
        elif i in hot_idx:
            region = "HOTSPOT"
        else:
            region = "SCAFFOLD"

        pdb_num = pdb_info.number(i) if pdb_info else i
        chain   = pdb_info.chain(i)  if pdb_info else "A"

        records.append({
            "scaffold_id":    scaffold_id,
            "motifs_applied": motifs_applied,
            "total_sasa":     f"{total_sasa:.3f}",
            "rosetta_index":  i,
            "pdb_resnum":     pdb_num,
            "chain":          chain,
            "resname":        pose.residue(i).name3(),
            "aa_1letter":     pose.residue(i).name1(),
            "region":         region,
            "sasa":           f"{rsd_sasa[i]:.3f}",
        })
    return records


SASA_COLUMNS = [
    "scaffold_id", "motifs_applied", "total_sasa",
    "rosetta_index", "pdb_resnum", "chain",
    "resname", "aa_1letter", "region", "sasa",
]

def init_sasa_csv(path):
    with open(path, "w", newline="") as f:
        csv.DictWriter(f, fieldnames=SASA_COLUMNS).writeheader()

def append_sasa_csv(path, records):
    with open(path, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SASA_COLUMNS)
        writer.writerows(records)


# ─────────────────────────────────────────────
# GERAÇÃO DOS PDBs FINAIS + SASA
# ─────────────────────────────────────────────

import pyrosetta

# Diretórios
scaffold_dir  = "relaxed_results/"
epitopes_dir  = "Epitopes/"
final_pdb_dir = os.path.join(output_dir, "final_pdbs/")
os.makedirs(final_pdb_dir, exist_ok=True)

sasa_csv = os.path.join(output_dir, f"sasa_per_residue_{timestamp}.csv")
init_sasa_csv(sasa_csv)
print(f"  SASA CSV iniciado: {sasa_csv}")

print(f"\n{'='*55}")
print(f"Gerando PDBs finais otimizados + SASA...")
print(f"{'='*55}")

generated = 0
skipped   = 0

for scaffold_id, data in results.items():
    selected_motifs = data["selected"]

    if not selected_motifs:
        print(f"  ⚠  {scaffold_id} — nenhum motif selecionado, pulando")
        skipped += 1
        continue

    pdb_path = os.path.join(scaffold_dir, f"{scaffold_id}.pdb")
    if not os.path.exists(pdb_path):
        print(f"  ✗  {scaffold_id} — PDB não encontrado: {pdb_path}")
        skipped += 1
        continue

    print(f"\n  Scaffold: {scaffold_id}")
    print(f"  Motifs selecionados: {[r['motif_name'] for r in selected_motifs]}")

    selected_sorted = sorted(selected_motifs, key=lambda r: r["scaffold_range_start"])

    scaffold_pose = pose_from_pdb(pdb_path)
    grafts_ok     = []

    for record in selected_sorted:
        motif_name = record["motif_name"]
        motif_path = os.path.join(epitopes_dir, motif_name)

        if not os.path.exists(motif_path):
            print(f"    ✗ Motif não encontrado: {motif_path} — pulando")
            continue

        try:
            objs  = XmlObjects.create_from_string(
                        xmlobject_MotifGrafting(pdb_path, motif_path))
            mover = objs.get_mover("motif_grafting")
            mover.apply(scaffold_pose)
            scaffold_pose = remove_context(scaffold_pose)
            grafts_ok.append(motif_name)
            print(f"    ✓ {motif_name} inserido")

        except RuntimeError as e:
            print(f"    ✗ {motif_name} falhou: {str(e).splitlines()[0]}")
            continue

    if grafts_ok:
        motifs_tag    = "_".join(os.path.splitext(m)[0] for m in grafts_ok)
        motifs_applied = ";".join(grafts_ok)
        out_pdb       = os.path.join(final_pdb_dir, f"{scaffold_id}__{motifs_tag}.pdb")
        scaffold_pose.dump_pdb(out_pdb)
        print(f"  ★  Salvo: {out_pdb}  ({len(grafts_ok)} motifs)")

        # ── Calcular SASA na pose final ──
        try:
            total_sasa, rsd_sasa = calculate_sasa(scaffold_pose)
            sasa_records = build_sasa_records(
                scaffold_pose, scaffold_id, motifs_applied, rsd_sasa, total_sasa
            )
            append_sasa_csv(sasa_csv, sasa_records)

            # Resumo SASA dos epitopos enxertados
            motif_sasa = [r for r in sasa_records if r["region"] == "MOTIF"]
            avg_epitope_sasa = (
                sum(float(r["sasa"]) for r in motif_sasa) / len(motif_sasa)
                if motif_sasa else 0
            )
            print(f"     SASA total: {total_sasa:.1f} Å²  |  "
                  f"SASA médio epitopos: {avg_epitope_sasa:.1f} Å²  |  "
                  f"resíduos epitopo: {len(motif_sasa)}")

        except Exception as e:
            print(f"     ⚠ SASA falhou: {e}")

        generated += 1
    else:
        print(f"  ⚠  {scaffold_id} — nenhum motif aplicado com sucesso")
        skipped += 1


# ─────────────────────────────────────────────
# RESUMO NO TERMINAL
# ─────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"Fase 2 concluída.")
print(f"\nTop 5 scaffolds (por motifs compatíveis):")
for row in ranking_rows[:5]:
    print(f"  #{row['rank']:2d}  {row['scaffold_id']:<25} "
          f"{row['n_selected']} motifs  |  "
          f"{row['n_overlapping_pairs']} conflitos")

total_selected = sum(r["n_selected"] for r in ranking_rows)
total_removed  = sum(r["n_removed"]  for r in ranking_rows)
print(f"\nTotal motifs mantidos  : {total_selected}")
print(f"Total motifs removidos : {total_removed}")
print(f"\nPDBs finais gerados    : {generated}")
print(f"Scaffolds pulados      : {skipped}")
print(f"Diretório final        : {final_pdb_dir}")
print(f"SASA por resíduo       : {sasa_csv}")
print(f"{'='*55}")


