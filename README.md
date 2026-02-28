# Graffiti (GRafting Algorithm For Finding, Incorporating and Tailoring Inserts)
<div align="center">
  <img src="grafiti-removebg-preview.png" width="400" />
</div>
![Graffiti](grafiti-removebg-preview.png)

A PyRosetta-based pipeline for multi-motif grafting onto protein scaffolds, with overlap optimization, SASA analysis, and structured logging.

---

## Overview

This pipeline automates the process of grafting multiple epitope fragments (motifs) onto protein scaffolds using Rosetta's [MotifGraft](https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/MotifGraftMover) algorithm. It is structured in two independent phases:

- **Phase 1** — Tests every `(scaffold, motif)` pair independently and maps where each graft lands in the scaffold sequence.
- **Phase 2** — Reads the Phase 1 map, detects overlapping grafts, selects the optimal non-overlapping combination per scaffold using a greedy algorithm, generates the final multi-motif PDB structures, and computes per-residue SASA for each final pose.

The goal is to find scaffolds that can accommodate the maximum number of epitopes simultaneously without structural conflicts, while providing quantitative metrics (graft position, SASA exposure) to support downstream design decisions.

---

## Directory Structure

```
project/
├── scaffolds_to_test/          # Input scaffold PDB files
├── Epitopes/                   # Input motif/epitope PDB files
├── Grafts_individual/          # Phase 1 output: individual graft PDBs + map CSV
├── Grafts_optimized/
│   ├── final_pdbs/             # Phase 2 output: multi-motif grafted structures
│   ├── overlap_decisions_*.csv # KEEP/REMOVE decision per (scaffold, motif)
│   ├── scaffold_ranking_*.csv  # Scaffolds ranked by compatible motif count
│   ├── sasa_per_residue_*.csv  # Per-residue SASA for all final poses
│   ├── ranking_chart_*.png     # Bar chart: motifs kept vs removed per scaffold
│   ├── distribution_chart_*.png# Histogram: scaffolds by number of compatible motifs
│   └── overlap_heatmap_*.png   # Conflict matrix for the most contested scaffold
├── phase1_independent_grafting.py
└── phase2_overlap_analysis.py
```

---

## Requirements

- Python 3.8+
- [PyRosetta](https://www.pyrosetta.org/) (licensed separately)
- `matplotlib`
- `numpy`

Install Python dependencies:
```bash
pip install matplotlib numpy
```

> PyRosetta requires an academic or commercial license. See [pyrosetta.org](https://www.pyrosetta.org/downloads) for installation instructions.

---

## Usage

### Phase 1 — Independent Grafting

Tests every `(scaffold, motif)` pair in isolation. Each scaffold is reloaded fresh for every attempt, so results are fully independent.

```bash
python phase1_independent_grafting.py
```

**Configure** at the top of the script:
```python
scaffold_dir = "scaffolds_to_test/"
epitopes_dir = "Epitopes/"
output_dir   = "Grafts_individual/"
```

**Output:** `Grafts_individual/graft_map_TIMESTAMP.csv`

| Column | Description |
|---|---|
| `scaffold_id` | Scaffold name (without `.pdb`) |
| `motif_name` | Motif filename |
| `status` | `SUCCESS` or `FAILED` |
| `scaffold_range_start` | First Rosetta index of the grafted region |
| `scaffold_range_end` | Last Rosetta index of the grafted region |
| `pdb_chain` | Chain of the inserted motif |
| `pdb_start_resnum` / `pdb_end_resnum` | PDB residue numbers of the motif |
| `motif_size` | Number of residues in the motif |
| `sequence` | Amino acid sequence inserted |
| `connection_resnums` | Junction residues (N-/C-terminal connections) |

---

### Phase 2 — Overlap Analysis, Optimization & SASA

Reads the Phase 1 CSV, resolves conflicts, generates final PDBs, and computes SASA.

```bash
python phase2_overlap_analysis.py
```

**Configure** at the top of the execution section:
```python
OVERLAP_BUFFER = 0   # Safety margin in residues (0 = any overlap counts)
scaffold_dir   = "scaffolds_to_test/"
epitopes_dir   = "Epitopes/"
output_dir     = "Grafts_optimized/"
```

Increasing `OVERLAP_BUFFER` allows small gaps between motifs to be tolerated — useful when scaffold flexibility may accommodate closely spaced grafts.

---

## How the Overlap Optimization Works

For each scaffold, Phase 2 runs a **greedy selection** algorithm:

1. All successfully grafted motifs for that scaffold are sorted by range size (smallest first).
2. Motifs are accepted one by one. A candidate is rejected if its `scaffold_range` overlaps with any already-accepted motif.
3. The result is the largest set of non-overlapping motifs for that scaffold.

This approach prioritizes smaller motifs (less invasive to scaffold structure) and maximizes the total count of compatible grafts.

---

## Output Files

### CSVs

| File | Contents |
|---|---|
| `graft_map_*.csv` | Phase 1 raw results — all pairs, success/failure, residue ranges |
| `overlap_decisions_*.csv` | Phase 2 per-pair decisions (`KEEP` / `REMOVE`) with conflict attribution |
| `scaffold_ranking_*.csv` | Scaffolds ranked by number of compatible motifs |
| `sasa_per_residue_*.csv` | Per-residue SASA (Å²) for every final grafted pose |

### SASA CSV Schema

| Column | Description |
|---|---|
| `scaffold_id` | Scaffold name |
| `motifs_applied` | Semicolon-separated list of applied motifs |
| `total_sasa` | Total SASA of the full pose (Å²) |
| `rosetta_index` | Internal Rosetta residue index (1-based) |
| `pdb_resnum` | PDB residue number |
| `chain` | Chain identifier |
| `resname` | 3-letter amino acid name |
| `aa_1letter` | 1-letter amino acid code |
| `region` | `MOTIF`, `SCAFFOLD`, `CONNECTION`, or `HOTSPOT` |
| `sasa` | Per-residue SASA (Å²) |

Filtering `region == MOTIF` isolates the SASA values of the grafted epitopes, which can be used as a downstream metric to assess epitope surface exposure — a key parameter for immunogenicity predictions.

### Plots

| File | Description |
|---|---|
| `ranking_chart_*.png` | Stacked bar chart — motifs kept vs removed per scaffold (split into parts if >50 scaffolds) |
| `distribution_chart_*.png` | Histogram of scaffolds by number of compatible motifs |
| `overlap_heatmap_*.png` | Conflict matrix (motif × motif) for the scaffold with the most overlaps |

---

## MotifGraft Parameters

Both scripts use the following MotifGraft settings, which can be adjusted directly in `xmlobject_MotifGrafting()`:

| Parameter | Default | Description |
|---|---|---|
| `RMSD_tolerance` | `3.0` Å | Maximum backbone RMSD for fragment alignment |
| `NC_points_RMSD_tolerance` | `10.0` Å | Maximum RMSD at N-/C-terminal junction points |
| `clash_score_cutoff` | `5` | Maximum tolerated atomic clashes |
| `clash_test_residue` | `GLY` | Residue used for clash testing (GLY = smallest) |
| `full_motif_bb_alignment` | `1` | Full backbone alignment (requires exact fragment size match) |

---

## Design Decisions & Notes

**Why test pairs independently in Phase 1?**
Grafting is sensitive to scaffold geometry. Testing each motif in isolation gives a clean, unbiased picture of where each epitope can land before any accumulation effects interfere.

**Why greedy and not exhaustive search?**
With many motifs and scaffolds, exhaustive search of all non-overlapping subsets is NP-hard. The greedy approach (smallest-first) is fast, interpretable, and performs well in practice for typical epitope counts (< 20).

**Why is SASA computed on the final pose only?**
SASA is context-dependent — the presence of other grafted motifs can bury or expose neighboring residues. Computing it on the final multi-motif structure gives the most realistic estimate of epitope accessibility.

**Why keep the labels from MotifGraft?**
Rosetta automatically labels residues as `MOTIF`, `SCAFFOLD`, `CONNECTION`, `CONTEXT`, and `HOTSPOT` after grafting. These labels persist in the pose and are used throughout the pipeline to identify regions without re-parsing PDB files.

---

## Reference

Silva, D., Correia, B.E., and Procko, E. (2016) *Motif-driven Design of Protein-Protein Interactions.* Methods Mol. Biol. 1414:285–304.

---

## License

This project is released under the MIT License. Note that PyRosetta itself is subject to its own separate licensing terms.
