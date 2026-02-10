# Spiky
A modular tool to predict fetal fraction from sequencing data using chromosome Y coverage.  
The tool processes BAM files, computes coverage statistics, predicts sample sex, identifies robust Y-chromosome regions, and builds/predicts fetal fraction (FFY) models.

---

## Table of Contents

- [Aims](#aims)  
- [Installation](#installation)  
- [Dependencies](#dependencies)  
- [Usage](#usage)  
  - [Modules](#modules)  
- [License](#license)  

---

## Aims

This tool is designed for **non-invasive prenatal testing (NIPT)** and related analyses.  
It allows users to:

1. Generate windowed coverage BED files from BAM sequencing data.  
2. Predict sample sex based on chromosome Y coverage.  
3. Identify robust regions on chromosome Y suitable for downstream analyses.  
4. Build linear regression models to estimate fetal fraction (FFY).  
5. Predict fetal fraction for new samples using the trained model.  

The modular design allows each step to be run independently or as part of a pipeline.

---

## Installation

1. Clone the repository:

Install dependencies:
pip install -r requirements.txt

Dependencies
Python 3.8+
pysam
numpy
scipy
pandas
scikit-learn

All dependencies can be installed via:
pip install pysam numpy scipy pandas scikit-learn

Usage
All modules are run through the main CLI script:
python ffy_calc/main.py --help
Modules are selected using double-dash flags, and each module has its own parameters.
Modules
--coverage-bed
Generate a windowed coverage BED file from a BAM.
python main.py --coverage-bed --bam sample.bam --bin-size 100000 --mapping-quality 20
--bam : Input BAM file
--bin-size : Window size in base pairs (default 1000)
--mapping-quality : Minimum mapping quality for reads (default 20)
Output: BED file with columns: chr, start, end, coverage, fractional coverage, lowq fraction
--predict-sex
Predict the sex of the sample using chromosome Y coverage.
python main.py --predict-sex --coverage-bed-file sample.coverage.bed --reference-fasta ref.fa
--coverage-bed-file : Input coverage BED file
--reference-fasta : Indexed reference FASTA
--y-frac-threshold : Threshold for calling male (default 0.02)
--lowq-threshold : Maximum fraction of low-quality reads per bin (default 0.1)
--max-n-frac : Maximum fraction of Ns per bin in reference (default 0.2)
Output: CSV table with predicted sex per sample.
--generate-regions
Generate robust Y-chromosome regions suitable for downstream analyses.
python main.py --generate-regions --sample-list predict_sex_outputs.txt --csv-out regions.csv
--sample-list : Text file listing predict-sex output files
--female-median-threshold : Maximum median chrY fraction for females (default 0.01)
--csv-out : Output CSV summarizing bins per sample
Output: BED and CSV files with median/percentile coverage statistics and status (ok or fail).
--model
Train a linear regression model to predict fetal fraction (FFY).
python main.py --model --training-csv training.csv --regions-bed regions.bed --model-out ffy_model.json
--training-csv : CSV file with BED path and fetal fraction
--regions-bed : BED of robust regions from generate-regions
--model-out : Output JSON file with linear regression model
Output: JSON model file and CSV table with predicted vs. observed FFY.
--predict
Predict FFY for a new sample using a trained model.
python main.py --predict --coverage-bed-file sample.coverage.bed --regions-bed regions.bed --model-file ffy_model.json
--coverage-bed-file : Input coverage BED file
--regions-bed : BED of robust regions
--model-file : Trained FFY model JSON
Output: Tab-separated file with columns: BED file, chrY median fraction, predicted FFY.
--version
Print the current version of the tool:
python main.py --version
License
This tool is released under the MIT License. See LICENSE for details.

---
