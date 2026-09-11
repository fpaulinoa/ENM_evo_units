# Ecological Niche Modelling of *Caiman latirostris* and its Evolutionary Units

Complete multi-algorithm ecological niche modelling (ENM) pipeline (MaxEnt,
Random Forest, GAM) for the broad-snouted caiman (*Caiman latirostris*) and
the five evolutionary units described by genomic data in Ornellas et al.
(2026): *Caiman_flumi*, *Caiman_gac*, *Caiman_gpi*, *Caiman_nocaa*,
*Caiman_norma*, *Caiman_paran*, and *Caiman_saof* (GAC and GPI are the two
aggregated clades: flumi+norma and nocaa+paran+saof).

This repository contains the full pipeline (8 units). For the version
restricted to the whole species (Article 1), see the companion repository.

## Requirements

Same as Repository 1 (see corresponding README), plus the `rvest` package
for reading MaxEnt HTML reports in `MNE012`.

## Expected data layout

Identical to Repository 1, except `output/Models Maxent/`, `output/Models
Random Forest/`, and `output/models_gam/` each contain one subfolder per
unit (`Caiman_flumi/`, `Caiman_gac/`, ..., `Caiman_saof/`,
`Caiman_latirostris/`).

## Execution order

Scripts `MNE00` through `MNE012` process all 8 units automatically (internal
loop), with no manual per-species editing required. Execution order matches
Repository 1 (see the table in that README), plus the evolutionary-unit
scripts below:

| Script | Description |
|--------|-------------|
| `MNE09b_rf_unidades.R` | Random Forest for the 7 evolutionary units/aggregates (excludes `Caiman_latirostris`) |
| `MNE09c_gam_unidades.R` | GAM for the 7 evolutionary units/aggregates |

These last two scripts support Article 2 (niche comparison across
evolutionary units) and are not required to reproduce the Article 1 results.

## Note on the niche-overlap script (Schoener's D between units)

A script to compute pairwise niche overlap between the 5 evolutionary units
(analogous to the original `ENM013b_Schoener_Hillinger_TODAS.R` design) is
planned but not yet finalised in this repository version. It should use
variable selection by name (not by position) and the already-cleaned
occurrence points in `output/pts_thin/`.

## Evaluation criteria

Identical to Repository 1: `wAUC ≥ 0.5`, `TSS ≥ 0`, `CBI ≥ 0.4`.

## Citation

If you use this code, please cite the associated articles (Article 1:
publication details to be added upon acceptance; Article 2: in preparation)
and this repository. For the evolutionary-unit genomic structure, also cite
Ornellas et al. (2026), *Molecular Phylogenetics and Evolution*.
