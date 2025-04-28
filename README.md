# Bayesian Methods to Predict CTCF Transcription Factor Binding Sites within the Human Genome

## Project Overview
This project uses a Bayesian statistical modeling approach to predict CTCF transcription factor binding sites in the human genome by integrating motif scores and ChIP-Seq signal data. We develop a generative model, perform Hamiltonian Monte Carlo inference, validate model diagnostics, and evaluate model performance.

---

## Repository Organization

###  `analysis/`
Contains main scripts and Jupyter notebooks for model development and evaluation:
- `model.py` — **[MAIN CODE]**: Runs Stan model fitting across different prior settings.
- `downsample.py` — Script to subsample large datasets for manageable modeling.
- `exploratory_analysis.ipynb` — Jupyter notebook for initial data exploration.
- `prior_posterior.ipynb` — Compares prior vs posterior distributions across models.
- `results_validation.ipynb` — **[MAIN NOTEBOOK]**: Performs model diagnostics, validation (trace plots, Rhat/ESS analysis, ROC curves).
- `posterior_outputs/` — Contains posterior draws, posterior predictive summaries, fitted models under different priors.

> Focus here on `model.py`, `prior_posterior.ipynb`, and `results_validation.ipynb` for core project logic and results.

---

###  `data/`
Contains raw and processed genomic data:
- `CTCF.meme` — Position Weight Matrix (PWM) for CTCF motif.
- `ENCFF278FNP.bed` — CTCF ChIP-Seq peaks from ENCODE.
- `hg38.fa` — Human genome sequence (large file).
- `ENCFF507CRU.bigWig` — ChIP-Seq coverage signal track (large file).
- BED files for positive and negative 500bp regions.

>  Large files like `hg38.fa` and `ENCFF507CRU.bigWig` were excluded from tar compression due to size.

---

###  `features/`
Contains processed datasets ready for modeling:
- `final_dataset.csv` — Full cleaned feature table (motif scores, ChIP signals, labels).
- `train_dataset.csv` and `test_dataset.csv` — Train/test splits for model training and validation.

---

###  `fimo_outputs/`
Contains motif scanning outputs from FIMO:
- `features/` — Motif scores (`fimo_neg_log_scores.tsv`, `fimo_pos_log_scores.tsv`) used for feature construction.
- Raw FIMO output directories for positives and negatives.

---

### `results/`
Contains visualization plots generated from exploratory data analysis and model outputs:
- Distribution plots, prior fits, likelihood fits, posterior diagnostics.

---

### `scripts/`
Helper scripts used for dataset preparation:
- `build_dataset.py` — **[IMPORTANT SCRIPT]**: Builds the final cleaned feature set used by the model.
- `extract_chip_signal.py` — Extracts ChIP-Seq signal overlaps.
- `process_fimo_scores.py` — Processes FIMO motif matches.
- `run_fimo.sh` — Bash script to automate FIMO motif scanning.

---

###  `sequences/`
Contains FASTA files for positive and negative sequences used for motif scanning:
- `positives.fa`
- `negatives.fa`

---

## Environment
The environment needed to run the code and reproduce results is specified in `environment.yml`. Install with:

```bash
conda env create -f environment.yml