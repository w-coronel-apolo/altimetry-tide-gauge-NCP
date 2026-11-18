# ALT–TG vertical land motion

MATLAB code and supplementary tables for the 131-tide-gauge analysis, January
1993–December 2024.

## Data and complete reproducibility package

The complete package will be archived on Zenodo as a ZIP containing the code,
all five input files, bundled mapping files, and saved outputs.

**Zenodo DOI/link: pending publication.**

For a complete reproduction, download and extract that ZIP and run `run_all`
inside the extracted folder. To use this GitHub copy instead, copy the entire
`data` folder from the matching Zenodo archive into this repository's root folder
before running the code. The required files are `sea_level.mat`,
`altimetry_monthly.h5`, `forcing.mat`, `river_monthly.mat`, and `ancillary.mat`;
their contents and sources are described in [DATA_README.md](DATA_README.md).

The four CSV outputs and bundled M_Map files are included here. The large analysis
inputs and generated MAT outputs are distributed through Zenodo.

## Run

After adding the input data, set this folder in MATLAB path and type in the command window:
`run_all`

The four (a01-a04) scripts rebuild the river indices, prepare sea level, fit the
models, and calculate the results and tables. Figures 1–9 and S1–S6 then open in
MATLAB.

For analysis-only runs or selected figures, type the following in command window:

```matlab
run_all(false)             % Analysis and tables only
main_figures               % All nine main figures from existing outputs
supplementary_figures      % All six supplementary figures
main_figures([7 8])        % Selected figures
```


Tested with MATLAB R2025b and Statistics and Machine Learning Toolbox.
[M_Map](https://www.eoas.ubc.ca/~rich/map.html) 1.4n and the low-, intermediate-
and high-resolution GSHHG 2.3.7 coastlines are included in `external/m_map`.
A numerical run takes approximately six minutes on the development computer.

## Files

| File | What it does |
|---|---|
| `a01_river_indices.m` | Build river indices from the screened monthly discharge records |
| `a02_prepare_sea_level.m` | Remove seasonal cycles, correlate ALT/TG, select offshore cells and form composites |
| `a03_fit_models.m` | Local wind PCA, stepwise regression, correction, component fits, and OLS/HAC trends |
| `a04_results_and_tables.m` | GNSS/GIA comparisons, window and forcing sensitivities, and the supplementary tables |
| `main_figures.m` | Create main figures, each in a labeled case |
| `supplementary_figures.m` | Supplementary figures, each in a labeled case |
| `functions/` | Shared OLS/HAC trend calculation |
| `external/m_map/` | Bundled M_Map code and the three coastline resolutions used by the figures |
| `data/` | Copy the five input files from the Zenodo archive; see `DATA_README.md` |
| `outputs/` | Four CSV files; running the analysis also creates the MAT results |

Use `run_all` to run the steps in order. Rerunning the first step starts a new
`outputs/results.mat`; complete the remaining steps before opening figures.
Existing files under `outputs/` are replaced when their step runs.

The principal outputs are listed below. The CSVs are included in this GitHub
copy; `run_all` creates the MAT files after the input data have been added.

- `results.mat`: river indices (`R`), sea-level observations/composites (`O`),
  regression and trend results (`A`), models and PCA bases (`Models`), component fits
  (`C`), and scientific comparisons/sensitivities (`D`).
- `grids.mat`: deseasonalized offshore ALT (`ALTds`) and ALT/TG correlations (`R_grid`).
- `river_indices.csv`: the 31 monthly river-index series, in m³/s.
- `Table_S1_station_rates.csv`, `Table_S2_model_summary.csv`, and
  `Table_S3_regression_coefficients.csv`: supplementary data, retaining the
  manuscript's selected columns.

## Variables

The four numbered files are ordinary scripts: open one and run it after the
preceding steps. They load their own inputs and leave variables in the workspace.
`run_all` runs them in its own workspace.

| Variable | Meaning |
|---|---|
| `M`, `L`, `N`, `t` | TG matrix (mm), coordinates [lon lat], station names, decimal years |
| `MRD`, `LRD`, `NRD`, `tRD` | River-index matrix (m3/s), coordinates [lon lat], river names, decimal years |
| `ALT` | Deseasonalized SLA+DAC composite for each TG, aligned with `M` and `t` (mm) |
| `Mds`, `Mdtds` | TG series after seasonal removal, or detrending then seasonal removal |
| `Models` | Regression fits and wind PCA bases; `M` is reserved for TG data |

Matrices have time down rows and stations across columns.
The saved structures group related results; `O.M`, `O.L`, `O.N`,
`O.t` and `O.ALT` are the sea-level variables, and `R.MRD`, `R.LRD`, `R.NRD`,
`R.tRD` are the river variables. Gridded altimetry remains in `grids.mat`.

## Conventions

Sea-level records are in millimeters. ALT–TG is altimetry minus
tide-gauge sea level, with positive trends indicating uplift. Positive
tide-gauge trends indicate relative sea-level rise. Trend uncertainties in
the tables are one standard error calculated using calendar-month
Newey–West HAC24 covariance. Figure 9 uses 95% confidence intervals for VLM
and one ensemble standard deviation for GIA.

Full raw and corrected trends can use different available months. Paired raw
trends use exactly the corrected record's months. Figure 7 uses full records;
Figure 9 and S5 use paired records. Figure 7's correction effect is corrected
minus raw, whereas S5 shows raw minus corrected.

The correction uses centered, non-detrended predictors and therefore retains
their trends. The full model/S5 and S4 detrend each series before complete-case selection.
Figures 3/5 select common months first and then detrend using their actual times.
Final rates use actual-time OLS with calendar-month HAC errors.
Predictor-selection and other details are documented beside their calculations.
Wind PCA uses MATLAB `pca`, with centered, unscaled columns and enough components to reach
at least 95% cumulative variance. Model selection uses MATLAB `stepwiselm` with an
intercept-only start, linear terms, the SSE partial-F criterion,
`PEnter=0.05`, and `PRemove=0.10`.
