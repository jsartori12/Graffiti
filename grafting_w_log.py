#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Motif Grafting com log CSV pivotado e gráfico de análise.
Colunas: Scaffold | frag1 | frag2 | ... | fragN | total_inserted | success_rate
"""

from pyrosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.select.residue_selector import ResiduePDBInfoHasLabelSelector
from pyrosetta.rosetta.protocols.grafting import return_region
import glob
import os
import csv
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

init()


# ─────────────────────────────────────────────
# XML
# ─────────────────────────────────────────────
def xmlobject_MotifGrafting(context_path, frag_path):
    xml_definition = f"""
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
    return xml_definition


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


def remove_labels_from_pose(pose, labels_to_remove):
    indices = get_residues_by_label(pose, labels_to_remove)
    new_pose = pose.clone()
    for i in range(new_pose.total_residue(), 0, -1):
        if i in indices:
            new_pose.conformation().delete_residue_slow(i)
    return new_pose


def keep_only_scaffold_and_motifs(pose):
    return remove_labels_from_pose(pose, {"CONTEXT"})


# ─────────────────────────────────────────────
# LOG PIVOTADO
# ─────────────────────────────────────────────
class GraftLogger:
    """
    Mantém uma tabela pivotada em memória:
      linhas  = scaffolds
      colunas = fragmentos
      valor   = "YES" / "NO"

    Ao final, salva CSV + gráfico.
    """

    def __init__(self, motifs_list, output_dir):
        self.output_dir  = output_dir
        # Nomes curtos dos fragmentos (sem path, sem extensão)
        self.frag_names  = [os.path.splitext(os.path.basename(m))[0] for m in motifs_list]
        # Dicionário principal: scaffold_name → {frag_name: "YES"/"NO"}
        self.table       = {}

    def init_scaffold(self, scaffold_name):
        """Registra um novo scaffold com todos os fragmentos como 'NO'."""
        self.table[scaffold_name] = {frag: "NO" for frag in self.frag_names}

    def mark_success(self, scaffold_name, motif_path):
        frag_name = os.path.splitext(os.path.basename(motif_path))[0]
        if scaffold_name in self.table and frag_name in self.table[scaffold_name]:
            self.table[scaffold_name][frag_name] = "YES"

    def save_csv(self, timestamp):
        """
        Salva o CSV pivotado a cada chamada (sobrescreve).
        Scaffold | frag1 | ... | fragN | total_inserted | pct_inserted
        """
        csv_path = os.path.join(self.output_dir, f"graft_log_{timestamp}.csv")
        fieldnames = ["scaffold"] + self.frag_names + ["total_inserted", "pct_inserted"]

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for scaffold_name, frags in self.table.items():
                n_yes   = sum(1 for v in frags.values() if v == "YES")
                n_total = len(frags)
                pct     = f"{(n_yes / n_total * 100):.1f}%" if n_total > 0 else "0.0%"
                row = {"scaffold": scaffold_name}
                row.update(frags)
                row["total_inserted"] = n_yes
                row["pct_inserted"]   = pct
                writer.writerow(row)

        print(f"  💾 CSV atualizado: {csv_path}")
        return csv_path

    def save_charts(self, timestamp):
        """
        Gera dois gráficos separados:
        1. Distribuição absoluta — quantos fragmentos foram aceitos no total (por fragmento)
        2. Percentual — % de scaffolds que aceitaram cada fragmento
        """
        if not self.table:
            return

        scaffolds   = list(self.table.keys())
        frags       = self.frag_names
        n_scaffolds = len(scaffolds)
        n_frags     = len(frags)

        # Percentual por fragmento (usado no gráfico 2)
        frag_pct = []
        for frag in frags:
            n_yes = sum(1 for scaf in scaffolds if self.table[scaf].get(frag) == "YES")
            frag_pct.append((n_yes / n_scaffolds * 100) if n_scaffolds > 0 else 0)

        bar_colors_pct = ["#2ecc71" if p >= 50 else "#e67e22" if p >= 25 else "#e74c3c"
                          for p in frag_pct]

        # ────────────────────────────────────────
        # GRÁFICO 1 — Distribuição: quantos scaffolds saíram com N fragmentos
        # ────────────────────────────────────────

        # Contar quantos fragmentos cada scaffold aceitou
        frags_per_scaffold = [
            sum(1 for v in self.table[scaf].values() if v == "YES")
            for scaf in scaffolds
        ]

        # Distribuição: para cada valor de 0..n_frags, quantos scaffolds tiveram aquele total
        dist = {i: 0 for i in range(n_frags + 1)}
        for count in frags_per_scaffold:
            dist[count] += 1

        x_vals   = list(dist.keys())
        y_vals   = list(dist.values())

        # Cor por quantidade: 0 = cinza, resto gradiente verde
        dist_colors = []
        for x in x_vals:
            if x == 0:
                dist_colors.append("#bdc3c7")
            elif x / n_frags >= 0.5:
                dist_colors.append("#2ecc71")
            elif x / n_frags >= 0.25:
                dist_colors.append("#e67e22")
            else:
                dist_colors.append("#e74c3c")

        fig1, ax1 = plt.subplots(figsize=(max(8, (n_frags + 1) * 1.1), 5))

        bars1 = ax1.bar(x_vals, y_vals, color=dist_colors, edgecolor="white",
                        linewidth=0.8, width=0.6)

        for bar, val in zip(bars1, y_vals):
            if val > 0:
                ax1.text(bar.get_x() + bar.get_width() / 2,
                         bar.get_height() + 0.15,
                         str(val),
                         ha="center", va="bottom", fontsize=11, fontweight="bold")

        ax1.set_xticks(x_vals)
        ax1.set_xticklabels([str(x) for x in x_vals], fontsize=10)
        ax1.set_xlabel("Número de fragmentos inseridos", fontsize=11)
        ax1.set_ylabel("Número de scaffolds", fontsize=11)
        ax1.set_title("Distribuição de Scaffolds por Número de Fragmentos Inseridos",
                      fontsize=13, fontweight="bold", pad=12)
        ax1.set_ylim(0, max(y_vals) + 2)
        ax1.grid(axis="y", alpha=0.3)

        leg_patches = [
            mpatches.Patch(color="#bdc3c7", label="0 fragmentos"),
            mpatches.Patch(color="#e74c3c", label="< 25% dos frags"),
            mpatches.Patch(color="#e67e22", label="25–49% dos frags"),
            mpatches.Patch(color="#2ecc71", label="≥ 50% dos frags"),
        ]
        ax1.legend(handles=leg_patches, fontsize=9, loc="upper right")

        chart1_path = os.path.join(self.output_dir, f"graft_chart_dist_{timestamp}.png")
        plt.tight_layout()
        plt.savefig(chart1_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  📊 Gráfico 1 (distribuição) salvo: {chart1_path}")

        # ────────────────────────────────────────
        # GRÁFICO 2 — Percentual por fragmento
        # ────────────────────────────────────────
        fig2, ax2 = plt.subplots(figsize=(max(8, n_frags * 1.1), 5))

        bars2 = ax2.bar(frags, frag_pct, color=bar_colors_pct, edgecolor="white", linewidth=0.8)

        for bar, pct in zip(bars2, frag_pct):
            ax2.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 1.0,
                     f"{pct:.1f}%",
                     ha="center", va="bottom", fontsize=10, fontweight="bold")

        ax2.set_xticks(range(n_frags))
        ax2.set_xticklabels(frags, rotation=35, ha="right", fontsize=9)
        ax2.set_ylabel("% de scaffolds que aceitaram", fontsize=11)
        ax2.set_xlabel("Fragmento (Epitopo)", fontsize=11)
        ax2.set_title("Taxa de Aceitação (%) por Fragmento",
                      fontsize=13, fontweight="bold", pad=12)
        ax2.set_ylim(0, 115)
        ax2.axhline(y=50, color="gray", linestyle="--",
                    linewidth=0.8, alpha=0.5, label="50% de referência")
        ax2.grid(axis="y", alpha=0.3)

        pct_patches = [
            mpatches.Patch(color="#2ecc71", label="≥ 50%"),
            mpatches.Patch(color="#e67e22", label="25–49%"),
            mpatches.Patch(color="#e74c3c", label="< 25%"),
        ]
        ax2.legend(handles=pct_patches, fontsize=9, loc="upper right")

        chart2_path = os.path.join(self.output_dir, f"graft_chart_pct_{timestamp}.png")
        plt.tight_layout()
        plt.savefig(chart2_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  📊 Gráfico 2 (percentual) salvo: {chart2_path}")


# ─────────────────────────────────────────────
# CONFIGURAÇÕES
# ─────────────────────────────────────────────
scaffold_dir = "scaffolds_to_test/"
epitopes_dir = "Epitopes/"
output_dir   = "Grafts_out/"
os.makedirs(output_dir, exist_ok=True)

pdb_files   = sorted(glob.glob(os.path.join(scaffold_dir, "*.pdb")))
motifs_list = sorted(glob.glob(os.path.join(epitopes_dir, "*.pdb")))

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
logger    = GraftLogger(motifs_list, output_dir)


# ─────────────────────────────────────────────
# LOOP PRINCIPAL
# ─────────────────────────────────────────────
for pdb_path in pdb_files:
    base_name   = os.path.splitext(os.path.basename(pdb_path))[0]
    output_path = os.path.join(output_dir, f"grafted_{base_name}.pdb")

    print(f"\n{'='*55}")
    print(f"Scaffold: {base_name}")
    print(f"{'='*55}")

    logger.init_scaffold(base_name)

    scaffold_pose   = pose_from_pdb(pdb_path)
    grafts_inserted = []

    for motif_path in motifs_list:
        motif_name = os.path.basename(motif_path)

        try:
            print(f"\n  → Tentando graft: {motif_name}")
            print(f"     Resíduos na pose atual: {scaffold_pose.total_residue()}")

            objs  = XmlObjects.create_from_string(
                        xmlobject_MotifGrafting(pdb_path, motif_path))
            mover = objs.get_mover("motif_grafting")
            mover.apply(scaffold_pose)

            scaffold_pose = keep_only_scaffold_and_motifs(scaffold_pose)
            grafts_inserted.append(motif_name)
            logger.mark_success(base_name, motif_path)

            print(f"  ✓ Graft aceito! Total inseridos: {len(grafts_inserted)}")

        except RuntimeError as e:
            error_msg = str(e)
            if "not suitable scaffold grafts" in error_msg:
                print(f"  ✗ Sem regiões compatíveis para {motif_name}")
            else:
                print(f"  ✗ Erro Rosetta: {error_msg.splitlines()[0]}")
            continue

        except Exception as e:
            print(f"  ✗ Erro inesperado: {e}")
            continue

    if grafts_inserted:
        scaffold_pose.dump_pdb(output_path)
        print(f"\n  ★ Salvo: {output_path}")
        print(f"     Motifs inseridos ({len(grafts_inserted)}): {grafts_inserted}")
    else:
        print(f"\n  ⚠ Nenhum graft bem-sucedido para: {base_name}")

    # ── Salvar CSV após cada scaffold (evita perda de dados) ──
    logger.save_csv(timestamp)


# ─────────────────────────────────────────────
# GRÁFICOS FINAIS
# ─────────────────────────────────────────────
print(f"\n{'='*55}")
print("Gerando gráficos finais...")
logger.save_charts(timestamp)
print("Processamento concluído.")