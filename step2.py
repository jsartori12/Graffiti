#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASE 2 — Análise de sobreposição + geração paralela de PDBs + SASA.

Cada scaffold é processado em um núcleo isolado via multiprocessing.Pool.
PyRosetta é inicializado dentro de cada worker para evitar conflitos de fork.

Configuração em CAPS na seção CONFIG abaixo.
"""

import os
import csv
import glob
import math
import shutil
import random
import multiprocessing
from datetime import datetime
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ══════════════════════════════════════════════════════════════════
# CONFIG  ← ajuste aqui
# ══════════════════════════════════════════════════════════════════
SCAFFOLD_DIR   = "few_dummies/"
EPITOPES_DIR   = "Epitopes/"
OUTPUT_DIR     = "Grafts_optimized/"
OVERLAP_BUFFER = 0      # 0 = qualquer sobreposição é conflito
TOP_PERCENT    = 5     # salva top X% por mean epitope SASA
N_WORKERS      = 128   # None = usa todos os CPUs disponíveis


# ══════════════════════════════════════════════════════════════════
# CARREGAR CSV DA FASE 1
# ══════════════════════════════════════════════════════════════════
def load_graft_map(csv_path):
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


# ══════════════════════════════════════════════════════════════════
# OVERLAP + GREEDY
# ══════════════════════════════════════════════════════════════════
def ranges_overlap(s1, e1, s2, e2, buffer=0):
    return not (e1 + buffer < s2 or e2 + buffer < s1)


def find_overlaps(motif_records, buffer=0):
    overlaps = []
    for i in range(len(motif_records)):
        for j in range(i + 1, len(motif_records)):
            a, b = motif_records[i], motif_records[j]
            if ranges_overlap(a["scaffold_range_start"], a["scaffold_range_end"],
                              b["scaffold_range_start"], b["scaffold_range_end"],
                              buffer=buffer):
                overlaps.append((a["motif_name"], b["motif_name"]))
    return overlaps


def greedy_one_pass(records, buffer=0):
    selected, rejected = [], []
    for candidate in records:
        conflict = next(
            (a["motif_name"] for a in selected
             if ranges_overlap(candidate["scaffold_range_start"],
                               candidate["scaffold_range_end"],
                               a["scaffold_range_start"],
                               a["scaffold_range_end"], buffer=buffer)),
            None
        )
        if conflict is None:
            selected.append(candidate)
        else:
            rejected.append({"record": candidate, "conflict_with": conflict})
    return selected, rejected


def greedy_select(motif_records, buffer=0, n_random_trials=200, seed=42):
    rng = random.Random(seed)
    best_sel, best_rej, best_score = [], [], -1

    for reverse in (False, True):
        ordered = sorted(motif_records,
                         key=lambda r: r["scaffold_range_end"] - r["scaffold_range_start"],
                         reverse=reverse)
        sel, rej = greedy_one_pass(ordered, buffer)
        if len(sel) > best_score:
            best_sel, best_rej, best_score = sel, rej, len(sel)

    shuffled = list(motif_records)
    for _ in range(n_random_trials):
        rng.shuffle(shuffled)
        sel, rej = greedy_one_pass(shuffled, buffer)
        if len(sel) > best_score:
            best_sel, best_rej, best_score = sel, rej, len(sel)

    return best_sel, best_rej


def analyze(records, buffer=0):
    by_scaffold = defaultdict(list)
    for r in records:
        by_scaffold[r["scaffold_id"]].append(r)

    results = {}
    for scaffold_id, motif_records in by_scaffold.items():
        overlaps           = find_overlaps(motif_records, buffer=buffer)
        selected, rejected = greedy_select(motif_records, buffer=buffer)
        results[scaffold_id] = {
            "all_motifs":          motif_records,
            "selected":            selected,
            "rejected":            rejected,
            "overlaps":            overlaps,
            "n_total":             len(motif_records),
            "n_selected":          len(selected),
            "n_rejected":          len(rejected),
            "n_overlapping_pairs": len(overlaps),
        }
    return results


# ══════════════════════════════════════════════════════════════════
# WORKER — roda em processo filho isolado
# Recebe tudo via args (picklable). Não usa globals do processo pai.
# ══════════════════════════════════════════════════════════════════
def _worker(args):
    """
    Processa um scaffold completo em um processo isolado:
      1. Inicializa PyRosetta silenciosamente
      2. Graft acumulativo dos motifs selecionados (menor range primeiro)
      3. Salva o PDB final
      4. Calcula SASA por resíduo
      5. Retorna dict com resultados serializáveis (sem objetos Rosetta)
    """
    (scaffold_id, selected_motifs,
     scaffold_dir, epitopes_dir, final_pdb_dir) = args

    # PyRosetta DEVE ser inicializado dentro do processo filho.
    # Com 'spawn', cada worker começa do zero — sem herdar estado do pai.
    try:
        import pyrosetta
        from pyrosetta import pose_from_pdb
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
        pyrosetta.init("-mute all")
    except Exception as e:
        return {
            "scaffold_id": scaffold_id, "status": "INIT_FAILED",
            "error": str(e), "grafts_ok": [], "sasa_records": [],
            "avg_ep_sasa": None, "out_pdb": None,
        }

    # ── helpers locais (evita problemas de pickle) ──────────────
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

    def make_sasa_records(pose, scaffold_id, motifs_applied, rsd_sasa, total_sasa):
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

    # ── graft acumulativo ────────────────────────────────────────
    pdb_path = os.path.join(scaffold_dir, f"{scaffold_id}.pdb")
    if not os.path.exists(pdb_path):
        return {
            "scaffold_id": scaffold_id, "status": "PDB_NOT_FOUND",
            "error": pdb_path, "grafts_ok": [], "sasa_records": [],
            "avg_ep_sasa": None, "out_pdb": None,
        }

    sorted_motifs = sorted(selected_motifs, key=lambda r: r["scaffold_range_start"])
    pose          = pose_from_pdb(pdb_path)
    grafts_ok     = []

    for rec in sorted_motifs:
        motif_path = os.path.join(epitopes_dir, rec["motif_name"])
        if not os.path.exists(motif_path):
            continue
        try:
            objs  = XmlObjects.create_from_string(
                        xml_motif_graft(pdb_path, motif_path))
            mover = objs.get_mover("motif_grafting")
            mover.apply(pose)
            pose = remove_context(pose)
            grafts_ok.append(rec["motif_name"])
        except RuntimeError:
            continue   # motif failed — continue with the next one

    if not grafts_ok:
        return {
            "scaffold_id": scaffold_id, "status": "NO_GRAFTS",
            "error": "", "grafts_ok": [], "sasa_records": [],
            "avg_ep_sasa": None, "out_pdb": None,
        }

    # ── salvar PDB ───────────────────────────────────────────────
    motifs_tag     = "_".join(os.path.splitext(m)[0] for m in grafts_ok)
    motifs_applied = ";".join(grafts_ok)
    out_pdb        = os.path.join(final_pdb_dir, f"{scaffold_id}__{motifs_tag}.pdb")
    pose.dump_pdb(out_pdb)

    # ── SASA ─────────────────────────────────────────────────────
    sasa_records = []
    avg_ep_sasa  = None
    try:
        total_sasa, rsd_sasa = calc_sasa(pose)
        sasa_records = make_sasa_records(
            pose, scaffold_id, motifs_applied, rsd_sasa, total_sasa
        )
        ep_vals     = [float(r["sasa"]) for r in sasa_records
                       if r["region"] == "MOTIF"]
        avg_ep_sasa = round(sum(ep_vals) / len(ep_vals), 3) if ep_vals else 0.0
    except Exception:
        pass

    return {
        "scaffold_id":  scaffold_id,
        "status":       "SUCCESS",
        "grafts_ok":    grafts_ok,
        "out_pdb":      out_pdb,
        "avg_ep_sasa":  avg_ep_sasa,
        "sasa_records": sasa_records,
        "error":        "",
    }


# ══════════════════════════════════════════════════════════════════
# CSV HELPERS
# ══════════════════════════════════════════════════════════════════
SASA_COLS = [
    "scaffold_id", "motifs_applied", "total_sasa",
    "rosetta_index", "pdb_resnum", "chain",
    "resname", "aa_1letter", "region", "sasa",
]


def write_csv(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})


def save_decision_csv(results, path):
    cols = ["scaffold_id", "motif_name", "decision",
            "scaffold_range_start", "scaffold_range_end",
            "motif_size", "sequence", "conflict_with",
            "pdb_chain", "pdb_start_resnum", "pdb_end_resnum"]
    rows = []
    for sid, data in sorted(results.items()):
        for r in data["selected"]:
            rows.append({**r, "scaffold_id": sid, "decision": "KEEP",
                         "conflict_with": ""})
        for item in data["rejected"]:
            rows.append({**item["record"], "scaffold_id": sid,
                         "decision": "REMOVE",
                         "conflict_with": item["conflict_with"]})
    write_csv(path, cols, rows)
    print(f"  📄 Decisões: {path}")


def save_ranking_csv(results, path):
    cols = ["rank", "scaffold_id", "n_motifs_grafted", "n_selected",
            "n_removed", "n_overlapping_pairs", "selected_motifs", "removed_motifs"]
    rows = []
    for sid, data in results.items():
        rows.append({
            "scaffold_id":         sid,
            "n_motifs_grafted":    data["n_total"],
            "n_selected":          data["n_selected"],
            "n_removed":           data["n_rejected"],
            "n_overlapping_pairs": data["n_overlapping_pairs"],
            "selected_motifs": ";".join(r["motif_name"] for r in data["selected"]),
            "removed_motifs":  ";".join(i["record"]["motif_name"]
                                        for i in data["rejected"]),
        })
    rows.sort(key=lambda x: (-x["n_selected"], -x["n_motifs_grafted"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
    write_csv(path, cols, rows)
    print(f"  📄 Ranking: {path}")
    return rows


# ══════════════════════════════════════════════════════════════════
# TOP X% BY MEAN EPITOPE SASA
# ══════════════════════════════════════════════════════════════════
def save_top_scaffolds(all_sasa_records, results, final_pdb_dir, output_dir,
                       top_percent, timestamp, min_grafts=None):
    """
    Generates two output folders:

    1. top{N}pct_by_epitope_sasa/
       Top X% of ALL scaffolds ranked by mean epitope SASA.

    2. max_grafts/
       ALL scaffolds that achieved the maximum number of successfully
       grafted epitopes (or >= min_grafts if specified).
       Sorted by mean epitope SASA within the group.
    """
    from collections import Counter

    # ── aggregate SASA + graft count per scaffold ─────────────────
    scaffold_data = defaultdict(lambda: {
        "motif_sasa": [], "motifs_applied": "", "total_sasa": 0.0
    })
    for row in all_sasa_records:
        sid = row["scaffold_id"]
        scaffold_data[sid]["motifs_applied"] = row["motifs_applied"]
        scaffold_data[sid]["total_sasa"]     = float(row["total_sasa"])
        if row["region"] == "MOTIF":
            scaffold_data[sid]["motif_sasa"].append(float(row["sasa"]))

    ranked = []
    for sid, data in scaffold_data.items():
        vals      = data["motif_sasa"]
        mean_ep   = round(sum(vals) / len(vals), 3) if vals else 0.0
        n_grafted = results.get(sid, {}).get("n_selected", 0)
        tag       = "_".join(os.path.splitext(m)[0]
                             for m in data["motifs_applied"].split(";") if m)
        ranked.append({
            "scaffold_id":       sid,
            "motifs_applied":    data["motifs_applied"],
            "n_grafts":          n_grafted,
            "n_motif_residues":  len(vals),
            "mean_epitope_sasa": mean_ep,
            "total_sasa":        round(data["total_sasa"], 3),
            "pdb_path":          os.path.join(
                                     final_pdb_dir,
                                     f"{sid}__{tag}.pdb" if tag else f"{sid}.pdb"),
        })

    # ── graft count distribution ──────────────────────────────────
    all_counts = [e["n_grafts"] for e in ranked]
    max_grafts = max(all_counts) if all_counts else 0
    threshold  = min_grafts if min_grafts is not None else max_grafts

    print(f"\n  Graft count distribution:")
    for count, n in sorted(Counter(all_counts).items(), reverse=True):
        bar  = "#" * min(n, 40)
        flag = "  <-- threshold" if count == threshold else ""
        print(f"    {count:>2} grafts : {n:>4} scaffolds  {bar}{flag}")

    # ── Ranking CSV (all scaffolds) ───────────────────────────────
    ranked_by_sasa = sorted(ranked, key=lambda x: x["mean_epitope_sasa"], reverse=True)
    n_top          = max(1, math.ceil(len(ranked_by_sasa) * top_percent / 100))

    ranking_csv = os.path.join(output_dir, f"sasa_ranking_{timestamp}.csv")
    cols = ["rank", "scaffold_id", "motifs_applied", "n_grafts",
            "n_motif_residues", "mean_epitope_sasa", "total_sasa",
            "in_top_pct", "in_max_grafts", "pdb_found"]
    rows = []
    max_grafts_ids = {e["scaffold_id"] for e in ranked if e["n_grafts"] >= threshold}
    top_pct_ids    = {e["scaffold_id"] for e in ranked_by_sasa[:n_top]}
    for rank, entry in enumerate(ranked_by_sasa, 1):
        rows.append({
            "rank":              rank,
            "scaffold_id":       entry["scaffold_id"],
            "motifs_applied":    entry["motifs_applied"],
            "n_grafts":          entry["n_grafts"],
            "n_motif_residues":  entry["n_motif_residues"],
            "mean_epitope_sasa": entry["mean_epitope_sasa"],
            "total_sasa":        entry["total_sasa"],
            "in_top_pct":        "YES" if entry["scaffold_id"] in top_pct_ids    else "NO",
            "in_max_grafts":     "YES" if entry["scaffold_id"] in max_grafts_ids else "NO",
            "pdb_found":         "YES" if os.path.exists(entry["pdb_path"])      else "NO",
        })
    write_csv(ranking_csv, cols, rows)
    print(f"  📄 SASA ranking: {ranking_csv}")

    # ── Folder 1: top X% by mean epitope SASA ────────────────────
    top_sasa_dir = os.path.join(output_dir, f"top{top_percent}pct_by_epitope_sasa/")
    os.makedirs(top_sasa_dir, exist_ok=True)
    copied = 0
    print(f"\n  📁 top{top_percent}pct_by_epitope_sasa/  ({n_top} scaffolds)")
    for rank, entry in enumerate(ranked_by_sasa[:n_top], 1):
        src = entry["pdb_path"]
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(top_sasa_dir, os.path.basename(src)))
            print(f"    #{rank:3d}  SASA={entry['mean_epitope_sasa']:.1f} Å²"
                  f"  grafts={entry['n_grafts']}"
                  f"  {os.path.basename(src)}")
            copied += 1
        else:
            print(f"    #{rank:3d}  ⚠ not found: {src}")
    print(f"  ✓ {copied}/{n_top} copied")

    # ── Folder 2: all max-grafts scaffolds ────────────────────────
    max_grafts_entries = sorted(
        [e for e in ranked if e["n_grafts"] >= threshold],
        key=lambda x: x["mean_epitope_sasa"], reverse=True
    )
    label = f"max{threshold}grafts"
    top_grafts_dir = os.path.join(output_dir, f"{label}/")
    os.makedirs(top_grafts_dir, exist_ok=True)
    copied2 = 0
    print(f"\n  📁 {label}/  ({len(max_grafts_entries)} scaffolds with >= {threshold} grafts)")
    for rank, entry in enumerate(max_grafts_entries, 1):
        src = entry["pdb_path"]
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(top_grafts_dir, os.path.basename(src)))
            print(f"    #{rank:3d}  grafts={entry['n_grafts']}"
                  f"  SASA={entry['mean_epitope_sasa']:.1f} Å²"
                  f"  {os.path.basename(src)}")
            copied2 += 1
        else:
            print(f"    #{rank:3d}  ⚠ not found: {src}")
    print(f"  ✓ {copied2}/{len(max_grafts_entries)} copied")

    return ranking_csv, top_sasa_dir, top_grafts_dir


# ══════════════════════════════════════════════════════════════════
# GRÁFICO
# ══════════════════════════════════════════════════════════════════
def save_chart(ranking_rows, timestamp, output_dir):
    all_n = [r["n_selected"] for r in ranking_rows]
    max_n = max(all_n) if all_n else 0
    dist  = {i: all_n.count(i) for i in range(max_n + 1)}

    colors = []
    for k in dist:
        if   k == 0:            colors.append("#bdc3c7")
        elif k == max_n:        colors.append("#2ecc71")
        elif k >= max_n * 0.5:  colors.append("#27ae60")
        elif k >= max_n * 0.25: colors.append("#e67e22")
        else:                   colors.append("#e74c3c")

    fig, ax = plt.subplots(figsize=(max(7, (max_n + 1) * 1.1), 5))
    bars = ax.bar(list(dist.keys()), list(dist.values()),
                  color=colors, edgecolor="white", linewidth=0.8, width=0.6)
    for bar, val in zip(bars, dist.values()):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.1, str(val),
                    ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_xticks(list(dist.keys()))
    ax.set_xlabel("Motifs selecionados (sem sobreposição)", fontsize=11)
    ax.set_ylabel("Número de scaffolds", fontsize=11)
    ax.set_title("Distribuição de Scaffolds por Motifs Compatíveis",
                 fontsize=13, fontweight="bold", pad=12)
    ax.set_ylim(0, max(dist.values()) + 2)
    ax.grid(axis="y", alpha=0.3)
    path = os.path.join(output_dir, f"distribution_chart_{timestamp}.png")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  📊 Gráfico: {path}")


# ══════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════
def main():
    import argparse
    n_avail = multiprocessing.cpu_count()

    parser = argparse.ArgumentParser(description="GRAFFITI Step 2 — Overlap + Filtering")
    parser.add_argument("--phase1_csv",    type=str, default=None,
                        help="Path to graft_map CSV from step 1. "
                             "If omitted, uses the latest in --grafts_dir.")
    parser.add_argument("--grafts_dir",    type=str, default="Grafts_individual/",
                        help="Directory containing step 1 graft_map CSVs (default: Grafts_individual/)")
    parser.add_argument("--scaffolds_dir", type=str, default=SCAFFOLD_DIR)
    parser.add_argument("--epitopes_dir",  type=str, default=EPITOPES_DIR)
    parser.add_argument("--output_dir",    type=str, default=OUTPUT_DIR)
    parser.add_argument("--overlap_buffer",type=int, default=OVERLAP_BUFFER)
    parser.add_argument("--top_percent",   type=int, default=TOP_PERCENT)
    parser.add_argument("--n_workers",     type=int, default=N_WORKERS or n_avail)
    parser.add_argument("--min_grafts",    type=int, default=None,
                        help="Minimum grafts to pass filter (default: auto = max observed)")
    args = parser.parse_args()

    # resolve phase1 CSV
    if args.phase1_csv:
        phase1_csv = args.phase1_csv
    else:
        matches = sorted(glob.glob(os.path.join(args.grafts_dir, "graft_map_*.csv")))
        if not matches:
            raise FileNotFoundError(f"No graft_map CSV found in {args.grafts_dir}")
        phase1_csv = matches[-1]

    scaffold_dir   = args.scaffolds_dir
    epitopes_dir   = args.epitopes_dir
    output_dir     = args.output_dir
    overlap_buffer = args.overlap_buffer
    top_percent    = args.top_percent
    n_workers      = args.n_workers
    min_grafts     = args.min_grafts

    os.makedirs(output_dir, exist_ok=True)
    final_pdb_dir = os.path.join(output_dir, "final_pdbs/")
    os.makedirs(final_pdb_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"{'='*55}")
    print(f"GRAFFITI — Step 2: Overlap + Filtering")
    print(f"  CSV Fase 1 : {phase1_csv}")
    print(f"  Scaffolds  : {scaffold_dir}")
    print(f"  Epitopes   : {epitopes_dir}")
    print(f"  Output     : {output_dir}")
    print(f"  Workers    : {n_workers}")
    print(f"  Buffer     : {overlap_buffer} residues")
    print(f"  Top %      : {top_percent}%")
    print(f"  Min grafts : {min_grafts or 'auto (max observed)'}")
    print(f"{'='*55}\n")

    # ── 2a: overlap analysis (single-process, fast) ──────────────
    print("Analysing overlaps...")
    records      = load_graft_map(phase1_csv)
    print(f"  SUCCESS records: {len(records)}")
    results      = analyze(records, buffer=overlap_buffer)

    decision_csv = os.path.join(output_dir, f"overlap_decisions_{timestamp}.csv")
    ranking_csv  = os.path.join(output_dir, f"scaffold_ranking_{timestamp}.csv")
    save_decision_csv(results, decision_csv)
    ranking_rows = save_ranking_csv(results, ranking_csv)
    save_chart(ranking_rows, timestamp, output_dir)

    # ── 2b: build task list ──────────────────────────────────────
    tasks = [
        (sid, data["selected"], scaffold_dir, epitopes_dir, final_pdb_dir)
        for sid, data in results.items()
        if data["selected"]
    ]
    print(f"\n{'='*55}")
    print(f"Generating {len(tasks)} PDB(s) in parallel ({n_workers} workers)...")
    print(f"{'='*55}")

    # ── 2c: parallel pool ────────────────────────────────────────
    ctx = multiprocessing.get_context("spawn")
    with ctx.Pool(processes=n_workers) as pool:
        worker_results = pool.map(_worker, tasks)

    # ── collect results ──────────────────────────────────────────
    all_sasa_records = []
    generated = failed = 0

    for res in worker_results:
        sid = res["scaffold_id"]
        if res["status"] == "SUCCESS":
            ep_str  = (f"  | SASA: {res['avg_ep_sasa']:.1f} Å²"
                       if res["avg_ep_sasa"] is not None else "")
            methods = res.get("methods_summary", "")
            meth_str = f"  | {methods}" if methods else ""
            print(f"  ✓ {sid:<35} {len(res['grafts_ok'])} motif(s){ep_str}{meth_str}")
            all_sasa_records.extend(res["sasa_records"])
            generated += 1
        else:
            print(f"  ✗ {sid:<35} [{res['status']}] {res['error']}")
            failed += 1

    # ── SASA CSV ─────────────────────────────────────────────────
    sasa_csv_path = os.path.join(output_dir, f"sasa_per_residue_{timestamp}.csv")
    write_csv(sasa_csv_path, SASA_COLS, all_sasa_records)
    print(f"\n  📄 SASA per residue: {sasa_csv_path}")

    # ── top X% selection ─────────────────────────────────────────
    print(f"\n{'='*55}")
    print(f"Selecting top {top_percent}% | graft filter: {min_grafts or 'max observed'}...")
    ranking_sasa_csv, top_sasa_dir, top_grafts_dir = save_top_scaffolds(
        all_sasa_records, results, final_pdb_dir,
        output_dir, top_percent, timestamp,
        min_grafts=min_grafts,
    )

    # ── summary ──────────────────────────────────────────────────
    total_sel = sum(r["n_selected"] for r in ranking_rows)
    total_rem = sum(r["n_removed"]  for r in ranking_rows)
    print(f"\n{'='*55}")
    print("Step 2 complete.")
    print(f"\nTop 5 scaffolds (by compatible motifs):")
    for row in ranking_rows[:5]:
        print(f"  #{row['rank']:2d}  {row['scaffold_id']:<28} "
              f"{row['n_selected']} motifs | {row['n_overlapping_pairs']} conflicts")
    print(f"\nMotifs kept      : {total_sel}")
    print(f"Motifs removed   : {total_rem}")
    print(f"PDBs generated   : {generated}")
    print(f"PDBs failed      : {failed}")
    print(f"Final PDB dir    : {final_pdb_dir}")
    print(f"SASA per residue : {sasa_csv_path}")
    print(f"SASA ranking     : {ranking_sasa_csv}")
    print(f"Top {top_percent}% SASA dir : {top_sasa_dir}")
    print(f"Max grafts dir   : {top_grafts_dir}")
    print(f"{'='*55}")


if __name__ == "__main__":
    main()
