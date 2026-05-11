# GRAFFITI
<div align="center">
  <img src="grafiti-removebg-preview.png" width="400" />
</div>

**G**rafting **R**outine for **A**utomated **F**ragment **F**itting, **I**nsertion, and **T**argeted **I**mmunogenic epitopes

A PyRosetta pipeline for transplanting immunogenic epitopes onto protein scaffolds, with secondary-structure-aware grafting, SASA-based ranking, and ProteinMPNN sequence design.

---

## Overview

GRAFFITI systematically tests every scaffold × epitope combination, resolves positional conflicts, selects top-performing scaffolds, and prepares sequences for structure prediction.

```
Scaffolds (PDB) ─┐
                  ├─► Step 1: Graft ─► Step 2: Filter ─► Step 3: MPNN ─► [AF2] ─► Step 5: Analyse
Epitopes  (PDB) ─┘
```

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `step1_paralel.py` | Parallel grafting of all scaffold × epitope pairs |
| 2 | `step2.py` | Overlap resolution, SASA scoring, top scaffold selection |
| 3 | `graffiti_run.py` | ProteinMPNN sequence design on selected scaffolds |
| 4 | *(external)* | AlphaFold2 structure prediction |
| 5 | `graffiti_run.py` | Post-AF2 analysis: pLDDT, RMSD, SASA per region |

---

## Secondary Structure-Aware Grafting

Epitopes are classified by DSSP into **HELIX**, **SHEET**, **LOOP**, or **MIXED**. The classification drives two completely different grafting strategies:

```
                    ┌─────────────────────────────┐
                    │   Scaffold + epitope PDB     │
                    │   (one job per worker)       │
                    └──────────────┬──────────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │      DSSP classification     │
                    │  Count H (helix), E (sheet), │
                    │  L (coil) residues via DSSP  │
                    └──────────────┬──────────────┘
                                   │
                      ┌────────────▼────────────┐
                      │    ≥ 60% coil residues? │
                      │   (LOOP_FRACTION = 0.6) │
                      └──────┬──────────┬───────┘
                          YES │          │ NO
                              │          │
               ┌──────────────▼─┐   ┌───▼──────────────────┐
               │  LOOP epitope  │   │  HELIX / SHEET / MIXED│
               └──────┬─────────┘   └───────────┬───────────┘
                      │                          │
      ┌───────────────▼──────────┐    ┌──────────▼──────────┐
      │ Find compatible scaffold │    │      MotifGraft      │
      │ loop (ep_len ≤ loop_len  │    │  Full backbone graft │
      │     ± MAX_SIZE_DELTA)    │    │  via RosettaScripts  │
      └──────────┬───────────────┘    └──────────┬──────────┘
                 │                               │
        ┌────────▼────────┐             ┌────────▼────────┐
        │  Loop found?    │             │  graft_method:  │
        └──┬──────────┬───┘             │   MOTIFGRAFT    │
        YES│          │NO               └─────────────────┘
           │          │
┌──────────▼───────┐  └──► FAILED_NO_LOOP
│ Mutate sequence  │
│ → CCD closure    │
│ → FastRelax      │
└──────────┬───────┘
           │
   ┌───────▼────────┐
   │ graft_method:  │
   │   LOOP_MODEL   │
   └────────────────┘
```

**Path A — structured epitope (HELIX / SHEET / MIXED):** PyRosetta's `MotifGraft` mover searches the scaffold for backbone segments that geometrically match the epitope's N- and C-terminal anchors, then grafts the full epitope backbone. The scaffold's 3D geometry changes.

**Path B — loop epitope (LOOP):** A full backbone graft is overkill for flexible loops. GRAFFITI finds an existing scaffold loop of compatible length, mutates its sequence residue-by-residue to match the epitope, then runs CCD to close chain breaks and FastRelax to resolve clashes. Backbone movement is minimal.

This logic lives in `ss_utils.py` and is shared by both worker modules.

---

## Requirements

- Python ≥ 3.8
- [PyRosetta](https://www.pyrosetta.org/) (licensed, install separately)
- NumPy
- Matplotlib
- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) (for Step 3)

Install Python dependencies:

```bash
pip install numpy matplotlib
```

PyRosetta requires a separate license and installation — see [pyrosetta.org](https://www.pyrosetta.org/downloads).

---

## Directory Structure

```
project/
├── graffiti_run.py        # Unified pipeline runner
├── step1_paralel.py       # Phase 1 orchestrator
├── step2.py               # Phase 2 orchestrator
├── phase1_worker.py       # Per-pair grafting worker
├── phase2_worker.py       # Per-scaffold cumulative grafting worker
├── ss_utils.py            # DSSP classification + loop modeling utilities
│
├── few_dummies/           # Input scaffold PDBs
├── Epitopes/              # Input epitope PDBs
├── Grafts_individual/     # Step 1 output (PDBs + graft_map CSV)
├── Grafts_optimized/      # Step 2 output (final PDBs, SASA CSVs, rankings)
├── MPNN_out/              # Step 3 output (designed sequences)
├── AlphaFold_results/     # Step 4 input (run externally)
└── AF2_analysis/          # Step 5 output (pLDDT, RMSD, SASA analysis)
```

---

## Usage

### Run full pipeline

```bash
python graffiti_run.py --steps 1 2 3
```

### Resume from Step 2 (grafting already done)

```bash
python graffiti_run.py --steps 2 3
```

### Run post-AF2 analysis only

```bash
python graffiti_run.py --steps 5
```

### Override config on the command line

```bash
python graffiti_run.py --steps 1 2 3 --cpus 16 --top_percent 10
```

### Dry run (print config and exit)

```bash
python graffiti_run.py --dry_run
```

---

## Configuration

All defaults live in the `CFG` dict at the top of `graffiti_run.py`. Every key can be overridden via CLI flag.

| Key | Default | Description |
|-----|---------|-------------|
| `scaffolds_dir` | `few_dummies/` | Input scaffold PDB directory |
| `epitopes_dir` | `Epitopes/` | Input epitope PDB directory |
| `grafts_dir` | `Grafts_individual/` | Step 1 output directory |
| `optimized_dir` | `Grafts_optimized/` | Step 2 output directory |
| `cpus` | 8 | Parallel workers for Step 1 |
| `overlap_buffer` | 0 | Residue margin between graft ranges |
| `top_percent` | 20 | % of top scaffolds carried to Step 3 |
| `min_grafts` | `None` (auto) | Minimum grafts to pass filter |
| `n_workers_step2` | 8 | Workers for Step 2 PDB generation |
| `mpnn_num_seq` | 10 | Sequences per scaffold (ProteinMPNN) |
| `mpnn_temp` | `0.2` | ProteinMPNN sampling temperature |
| `rmsd_threshold` | 2.0 Å | Epitope RMSD cutoff for quality flag |
| `plddt_threshold` | 70.0 | Minimum pLDDT for quality flag |

---

## Step 1 in isolation

```bash
python step1_paralel.py \
  --scaffolds my_scaffolds/ \
  --epitopes  my_epitopes/ \
  --output    Grafts_individual/ \
  --cpus      8
```

Outputs a timestamped `graft_map_YYYYMMDD_HHMMSS.csv` with one row per pair, recording graft method, scaffold range, sequence, and SASA.

---

## Step 2 in isolation

```bash
python step2.py \
  --phase1_csv    Grafts_individual/graft_map_*.csv \
  --scaffolds_dir few_dummies/ \
  --epitopes_dir  Epitopes/ \
  --output_dir    Grafts_optimized/ \
  --top_percent   20
```

Outputs:
- `final_pdbs/` — all resolved scaffold PDBs
- `top{N}pct_sasa/` — top % by mean epitope SASA
- `max{N}grafts/` — scaffolds with maximum compatible grafts
- `sasa_per_residue_*.csv` — per-residue SASA table
- `scaffold_ranking_*.csv` — ranked scaffold summary
- `distribution_chart_*.png` — motif count distribution plot

---

## Key Parameters in `ss_utils.py`

| Constant | Default | Description |
|----------|---------|-------------|
| `LOOP_FRACTION` | 0.6 | Coil fraction threshold to classify as LOOP |
| `MAX_SIZE_DELTA` | 5 | Max length difference (scaffold loop vs epitope) |
| `FLANK_RESIDUES` | 3 | Flanking residues included in CCD/relax |
| `CCD_CYCLES` | 200 | CCD loop closure iterations |
| `RELAX_ROUNDS` | 0 | FastRelax rounds (0 = minimal) |

---

## Multiprocessing Notes

PyRosetta is not fork-safe due to C++ global state. Both Step 1 and Step 2 use `multiprocessing.get_context("spawn")` so each worker initialises its own clean PyRosetta instance. This is handled automatically — no user action required.

---

## Output CSV Columns (Step 1)

| Column | Description |
|--------|-------------|
| `scaffold_id` | Scaffold PDB stem |
| `motif_name` | Epitope filename |
| `status` | `SUCCESS` / `FAILED` / `FAILED_NO_LOOP` |
| `graft_method` | `MOTIFGRAFT` / `LOOP_MODEL` / `FAILED_NO_LOOP` |
| `epitope_ss_class` | `HELIX` / `SHEET` / `LOOP` / `MIXED` |
| `scaffold_range_start` | Rosetta start index of graft |
| `scaffold_range_end` | Rosetta end index of graft |
| `motif_size` | Epitope length (residues) |
| `sequence` | Grafted sequence (one-letter) |
| `scaffold_total_residues` | Total scaffold residues |
| `error_message` | Failure reason (empty on success) |

---

## License

This project uses PyRosetta, which requires a separate academic or commercial license from the [Rosetta Commons](https://www.pyrosetta.org/downloads#academic). All other code in this repository is released under the MIT License.
