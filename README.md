# tumor-metabolomics

A classifier that distinguishes lung adenocarcinoma tumor tissue from adjacent normal
tissue using GC-TOF metabolomics data (Metabolomics Workbench study
[ST000390](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST000390)).

**Status:** stage 1 (data loading and quality control) complete.

All data facts, methods and analysis decisions, with their reasons, are recorded in
[NOTES.md](NOTES.md).

## Data

Wikoff et al., *Cancer Prevention Research* 2015; UC Davis West Coast Metabolomics
Center. 39 patients (tumor + adjacent normal tissue each) plus 9 pooled QC samples,
183 named metabolites. Licensed CC BY-NC-ND 4.0, so the data are **not** stored in this
repository; the steps below download them from Metabolomics Workbench.

## Reproduce the analysis

Requirements: macOS or Linux, `bash`, `curl`, `git`, and Python 3.13. Run all commands
from the repository root.

1. **Get the code**
   ```bash
   git clone https://github.com/k3rrye/tumor-metabolomics.git
   cd tumor-metabolomics
   ```

2. **Create the Python environment** with the pinned package versions
   ```bash
   python3.13 -m venv .venv
   .venv/bin/python -m pip install -r requirements.txt
   ```
   Any Python 3.13 works (python.org, Homebrew, Anaconda); check with
   `python3.13 --version`.

3. **Download the data** (about 0.6 MB) into `data/`
   ```bash
   scripts/download.sh ST000390
   ```
   This saves the study summary, data, factors and mwTab file. ST000391 (the LC-MS
   companion study) is not used in the analysis; `scripts/download.sh ST000391`
   fetches it if wanted (about 5 MB).

4. **Stage 1: load and QC**
   ```bash
   .venv/bin/python scripts/01_load_qc.py
   ```
   The script stops with an error if the downloaded data differ from what the
   analysis expects (sample counts, groups, pairing, missing values), for example if
   the study is revised on Metabolomics Workbench. The analysis used study version 2,
   revision 1.

   Outputs:
   - `data/processed/samples.csv`, `intensities.csv`, `features.csv`: tidy tables
     used by later stages (not in git; rebuilt by this step)
   - `reports/01_qc.txt`: integrity checks and QC report
   - `figures/01_*.png`: QC figures

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Analysis scripts, numbered by stage |
| `reports/` | Text reports written by each stage |
| `figures/` | Figures written by each stage |
| `data/` | Downloads and processed tables (git-ignored, rebuilt by the steps above) |
| `NOTES.md` | Data facts, pairing method, decisions and stage results |
| `requirements.txt` | Pinned Python package versions |
