# tumor-metabolomics

A classifier that distinguishes lung adenocarcinoma tumor tissue from adjacent normal
tissue using GC-TOF metabolomics data (Metabolomics Workbench study
[ST000390](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?StudyID=ST000390)).

**Status:** stage 1 (data loading and quality control) complete.

## Data

Wikoff et al., *Cancer Prevention Research* 2015; UC Davis West Coast Metabolomics
Center. Malignant and adjacent nonmalignant lung tissue from 39 current or former
smokers with stage IA/IB adenocarcinoma, plus 9 pooled QC samples: 87 samples and 183
named metabolites measured by GC-TOF mass spectrometry. Licensed
[CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/), so the data are
**not** stored in this repository; the steps below download them from Metabolomics
Workbench.

Points about the deposited data that affect the analysis:

- **No patient identifier.** Every sample has its own subject ID, so which tumor and
  normal sample came from the same patient has to be inferred (see
  [Inferring patients](#inferring-patients)).
- **Already normalised.** Total signal is identical across all 87 samples to within
  0.06%, so the values are already total-signal normalised. That normalisation is not
  applied again.
- **The 9 extra samples are pooled QC samples, not blanks.** They are unlabelled
  (subject IDs RS01–RS09), but their total signal matches the tissue samples, their
  profile sits between the tumor and normal averages, and their metabolite variability
  is about 3× lower than across patients, as expected for repeated injections of one
  pool. They are used only for quality control, never in tests or the classifier.
- **Method metadata is unreliable.** The deposited sample-preparation text is a copied
  plasma protocol; the paper is the authority on methods.
- **Run order is inferred.** Sample IDs encode two run dates (`121129`, `121130`)
  and a position, which are used as run day and run order; neither is stated in the
  metadata.

A companion study, ST000391 (HILIC-QTOF LC-MS on the same samples), is not used: its
8,869 features are unannotated, which would make results hard to interpret and
overfitting likely with 78 samples.

## Methods and decisions

### Inferring patients

Tumor and normal samples are matched on all recorded clinical fields: sex, smoking
status, diagnosis, stage, age and race. A combination shared by exactly one tumor and
one normal sample is taken to be one patient. This gives **36 unique pairs**
(patients P01–P36, 72 samples). The other 6 samples (3 tumor, 3 normal) have identical
fields and cannot be paired; they form one group, `AMB`.

This assumes both samples from a patient carry identical recorded fields and that no
two patients share all six. Nothing in the metadata contradicts this, but it cannot be
verified, so the results are checked for sensitivity to it (decision 2).

### Analysis decisions

| # | Decision | Reason |
|---|---|---|
| 1 | **Cross-validation is grouped by inferred patient**: 36 two-sample groups plus the six-sample `AMB` group, so no patient's samples are split between training and test folds. | A tumor and normal sample from one patient share patient-level metabolism. Splitting them would let the model recognise the patient rather than the tissue, inflating accuracy. |
| 2 | **Paired tumor-vs-normal tests use the 36 unique pairs**; unpaired tests on all 78 tissue samples are run alongside as a sensitivity check. | Paired tests need to know which samples go together, so the ambiguous 6 are left out. Because the pairs are inferred, agreement with the unpaired tests shows the conclusions don't depend on the pairing. |
| 3 | **Metabolites with QC-pool CV > 30% are dropped** (33 of 183; 150 kept). | 30% is a common acceptance threshold for untargeted GC-MS. The filter uses only the unlabelled pools, so applying it before cross-validation leaks no tumor/normal information. |
| 4 | **Run day is not corrected for; it is included as a covariate in statistical tests.** | Tumor and normal are unevenly spread over the two run days (Fisher p = 0.068), and only 2 pools ran on day 2, too few for a reliable pool-based correction. |
| 5 | **Outlier sample `121130bbdsa23_1` is kept**; main results are repeated without it. | It shows no sign of technical failure (neighbouring injections are normal, and its extreme metabolites are ones measured precisely) and it is the most tumor-like tumor, so removing it would understate tumor/normal differences. |

Stage 1 changes no values; the filter in decision 3 and the analyses in decisions 1,
2, 4 and 5 are applied in later stages.

## Stage 1 results: data loading and QC

Full numbers are in [reports/01_qc.txt](reports/01_qc.txt).

- **Integrity:** all checks pass. 87 samples × 183 metabolites; 39 tumor, 39 normal, 9
  pools; no missing, zero or negative values; values identical in the two Metabolomics
  Workbench formats.
- **Measurement precision:** median CV per metabolite is 15% across QC pools versus 56%
  across tissue samples. 33 metabolites exceed 30% in the pools, and 7 vary as much
  across pools as across patients.
- **Run order and day:** no drift within runs (pools show no trend with run order).
  The median shift between days is small (−0.08 log2), but 14 metabolites shift more
  than 1.4-fold, and 18 of the 36 pairs have their two samples on different days.
- **Overall structure (PCA, log2, autoscaled):** pools cluster tightly at the centre;
  tumor and normal overlap substantially, with PC1 (18% of variance) moderately
  correlated with tumor status (r = 0.34). Samples don't cluster by run day.
- **Outliers:** one sample, `121130bbdsa23_1` (tumor), exceeds the 99% Hotelling T²
  limit; it is kept (decision 5).

| Figure | Shows |
|---|---|
| [01_pca_group.png](figures/01_pca_group.png) | PCA coloured by tumor / normal / pool |
| [01_pca_runday.png](figures/01_pca_runday.png) | PCA coloured by run day |
| [01_run_order.png](figures/01_run_order.png) | PC1 and PC2 scores in run order |
| [01_cv_pools.png](figures/01_cv_pools.png) | Metabolite CV in QC pools vs tissue samples |

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
     used by later stages (not in git; rebuilt by this step). `samples.csv` holds
     the inferred `patient_id` and the cross-validation group `cv_group`.
   - `reports/01_qc.txt`: integrity checks and QC report
   - `figures/01_*.png`: QC figures

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Analysis scripts, numbered by stage |
| `reports/` | Text reports written by each stage |
| `figures/` | Figures written by each stage |
| `data/` | Downloads and processed tables (git-ignored, rebuilt by the steps above) |
| `requirements.txt` | Pinned Python package versions |
