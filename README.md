# Foot morphology statistics

R code for pediatric foot morphology reference curves, age gradients, classification, adjusted associations, and optional international comparisons. No research data are included.

The pipeline includes the Gamma log-link AV mean model used for Figure 2. Fixed-foot-length age gradients are expressed in mm³/year, with participant-cluster robust covariance and simultaneous bands across ages and both sexes. AV stabilisation, cumulative change, and stable-window distributions use this model. Inference uses 5,000 joint Gaussian coefficient perturbations, conditional on the reference foot-length path.

The log-AV GAM remains the conditional reference model for centiles, cross-validation, and classifications. Its original gradients are retained as `av_log_scale_*` outputs. Gamma bilateral means and right-minus-left contrasts are saved in `results/av_gamma/`.

## Run

From the repository directory:

```sh
Rscript INSTALL_PACKAGES.R
Rscript CHECK_CODE.R
Rscript RUN_ALL.R
```

Dependencies: `mgcv`, `geepack`, `gamlss`, `gamlss.dist`, and `sandwich`. Tested versions are listed in `dependencies.csv`. Code checks use synthetic observations.

Before analysis, provide `data/observations.csv` and edit `00_config.R`. Set `stage = "check"` for input validation, `"quick"` for a reduced computation, or `"full"` for the configured analysis. Outputs are written under `results/`.

The original numerical defaults are retained, including the expected sample of 4,929 children and 9,858 feet. For another cohort, update those counts or disable `require_expected_n`. Keep generated data, model objects, and outputs local.

## Input

`data/observations.example.csv` contains the required header only. Each child must have exactly two rows, one per side, with consistent demographic values.

| Column | Definition or coding |
| --- | --- |
| `participant_id` | Child identifier shared by both rows |
| `gender` | 1 = boys; 2 = girls |
| `age` | Years, from 3 inclusive to 19 exclusive |
| `height`, `weight`, `BMI` | cm, kg, kg/m² |
| `side` | 1 = Left; 2 = Right |
| `FL` | Foot length, mm |
| `AI` | Arch index, ratio |
| `hallux` | Angle, degrees |
| `AS`, `S` | Areas, mm² |
| `AV` | Volume, mm³ |
| `AS/S` | Arch area ratio (AAR); internally named `AS_ratio` for compatibility |

All 14 columns are required. Missing anthropometry is allowed; demographic values must agree between the two rows. Identifier, gender, age, and side cannot be missing. Use blank cells or `NA` for missing measurements. BMI and AAR are recomputed and checked against supplied values.