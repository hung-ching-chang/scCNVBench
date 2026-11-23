# Changes to scDesign3

This repository adapts the original [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3) R scripts. The table compares the author-provided originals with the versions in this repository. **Line numbers refer to the modified files before attribution headers were added** (and may shift with later edits).

| Script | Modified lines | Changes relative to original |
|---|---|---|
| **`scDesign3_mod.R`** | **1–16, 21–47, 49–138, 140–199, 205–309, 312–334** | **Main executable script adapted from the scDesign3 workflow (`scdesign3.R`).** Sources the supporting scripts; clusters reference cells; builds subclone-specific CNV factors from DNA dosage and Gamma-distributed signal strengths; runs CNV-adjusted scDesign3 simulation; saves synthetic count and cell-label files; and plots a heatmap. This is a study-specific driver script, not a direct line-by-line replacement for the upstream `scdesign3.R`. |
| `scdesign3.R` | 113–114, 202–211 | Added optional `cnv_factor_mat`; multiplies the fitted mean matrix elementwise by CNV factors before simulation (defaults to a matrix of ones). Added dimension diagnostics. |
| `fit_copula.R` | 242, 326 | Replaced group-level `lapply()` with `mclapply()` and set `mc.cores = 5`. |
| `simu_new.R` | 91–181, 200, 241–247, 275–278, 292, 329–331 | Added `sample_copula_group()` helper (not currently called), printed mean-matrix dimensions, switched group sampling to `mclapply()` with five workers, and made minor formatting/comment changes. |
| `extract_para.R` | 101; original line 276 deleted | Removed `as.numeric()` from `data$gene` assignment; deleted one commented-out line. |
| `plot_reduceddim.R` | 121–122, 126–127, 135 | Replaced programmatic ggplot aesthetic mappings with `aes_string()` for PCA/UMAP plots; minor formatting change. |
| `data.R` | Original lines 6–17 deleted | Removed documentation entries for `example_count` and `pseudotime`; retained `example_sce`. |
| `construct_data.R` | — | No change. |
| `fit_marginal.R` | — | No change. |
| `gamlss_fix.R` | — | No change. |
| `perform_lrt.R` | — | No change. |
| `sparse_cov.R` | — | No change. |

## Running the simulation

Run **`scDesign3_mod.R`** as the main script; the supporting `.R` files are loaded automatically with `source()`. Before running, set `scDesign3_path`, the reference count-matrix and CNV input paths, and the output directory in the script. Supply the command-line Gamma parameters `--beta_amp`, `--alpha_amp`, `--beta_del`, and `--alpha_del` as appropriate (the beta parameters have no usable defaults). Example:

```bash
Rscript scDesign3_mod.R --alpha_amp 1 --beta_amp <amp_scale> --alpha_del 1 --beta_del <del_scale>
```

The script writes simulated counts (`*.scDesign3Simulated.txt`), simulated cell labels (`*.scDesign3Simulated.CellTypeLabel.txt`), and a count heatmap (`*.count_heatmap.png`). The example input paths and dataset-specific assumptions must be adapted to the user’s data.

## Attribution

These scripts are derived from the [original scDesign3 project](https://github.com/SONGDONGYUAN1994/scDesign3). Cite the original scDesign3 publication and repository when using this adaptation, and retain the original copyright and license notices when redistributing the code.
