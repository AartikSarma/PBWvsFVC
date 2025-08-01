# PBW vs FVC - Predicted Body Weight versus Forced Vital Capacity in Mechanical Ventilation

**Initial Preprint**: [SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4898478)

**Published Paper**: [The Lancet Respiratory Medicine](https://www.thelancet.com/journals/lanres/article/PIIS2213-2600(25)00126-2/abstract)

## Repository Structure

### Core Analysis Files (Used for the preprint)

- **`1 - VTFVC - Simulated data.Rmd`** - Generates simulated patient data and calculates PBW-based and FVC-based tidal volumes across different demographics
- **`2 - VTFVC - RICU data import.Rmd`** - Imports and processes data from MIMIC-IV and eICU databases
- **`3 - VTFVC - Incusion criteria.Rmd`** - Applies inclusion/exclusion criteria to create the final analytical cohort
- **`4 - VTFVC - Table 1.Rmd`** - Creates baseline characteristics table (Table 1) for the study population
- **`5 - VTFVC - Regression and sensitivity analyses.R`** - Performs primary regression analyses and sensitivity analyses
- **`6 - VTFVC - Regression Tables.Rmd`** - Formats regression results into publication-ready tables
- **`7- VTFVC - Survival curves.R`** - Generates survival curves and related analyses

### Additional Analyses (For the published paper)

- **`PBW vs FVC analyses and figures.qmd`** - This file includes generates figures and analysis used in the final manuscript. Notable changes from the original preprint include a direct comparison of PBW and predicted FVC (rather than VT/PBW and VT/FVC) and additional investigation of driving pressure and normalized elastance.


## Requirements

### R Version
- R 4.0.0 or higher recommended

### Key R Packages

#### Data Manipulation & Core
- `tidyverse` - Core data science packages
- `data.table` - Fast data manipulation
- `lubridate` - Date/time handling

#### Statistical Analysis
- `lme4`, `lmerTest` - Mixed-effects models
- `survival`, `survminer` - Survival analysis
- `marginaleffects` - Marginal effects calculations
- `broom`, `broom.mixed` - Tidy model outputs

#### Visualization
- `ggplot2` (via tidyverse) - Core plotting
- `patchwork` - Combine plots
- `ggpubr`, `ggbeeswarm` - Publication-ready plots
- `ggdag`, `dagitty` - DAG visualization

#### Specialized
- `rspiro` - Spirometry calculations
- `furrr` - Parallel processing
- `officer`, `flextable` - Document generation


## Running the Analysis

### For the Preprint Analysis:
1. Run files 1-7 in numerical order
2. Each R Markdown file can be knitted to generate outputs
3. Results and figures will be saved to `results/` and `figures/` directories

### For the Published Paper Analysis:
1. Ensure the preprint analysis has been completed first
2. Run `PBW vs FVC analyses and figures.qmd` to generate additional analyses

## Output

- **Tables**: Generated as Word documents and CSV files in the `results/` directory
- **Figures**: Saved as high-resolution images in the `figures/` directory
- **Statistical Results**: Regression outputs and model summaries saved as RDS files

## Data Sources

This analysis uses data from:
- **MIMIC-IV**: Medical Information Mart for Intensive Care IV database
- **eICU**: eICU Collaborative Research Database

Both databases are available from [PhysioNet](https://physionet.org/) and are subject to PhysioNet's data use agreements.

## Reference 
Sarma, et al. Evaluating the generalisability of formulas used to set tidal volumes in mechanically ventilated patients: an observational, multicohort, retrospective study. Lancet Respiratory Medicine 2025
