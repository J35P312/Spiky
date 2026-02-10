# Spiky

Spiky is a modular tool to predict fetal fraction (FFY) from sequencing data using chromosome Y coverage.
It complements the **Fluffy** NIPT pipeline and allows users to generate coverage statistics, predict sample sex, identify robust Y-chromosome regions, and estimate fetal fraction.

---

## Table of Contents

* [Aims](#aims)
* [Installation](#installation)
* [Dependencies](#dependencies)
* [Usage](#usage)

  * [Modules](#modules)
* [License](#license)

---

## Aims

Spiky is designed for **non-invasive prenatal testing (NIPT)** and related analyses.

It allows users to:

1. Generate windowed coverage BED files from BAM sequencing data.
2. Predict sample sex based on chromosome Y coverage.
3. Identify robust regions on chromosome Y suitable for downstream analyses.
4. Build linear regression models to estimate fetal fraction (FFY).
5. Predict fetal fraction for new samples using a trained model.

The modular design allows each step to be run independently or as part of a pipeline.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/J35P312/Spiky.git
cd Spiky
```

Install required packages:

```bash
pip install -r requirements.txt
```

Make sure the main script is executable:

```bash
chmod +x ffy_calc/main.py
```

---

## Dependencies

* Python 3.8+
* [pysam](https://pysam.readthedocs.io/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)
* [pandas](https://pandas.pydata.org/)
* [scikit-learn](https://scikit-learn.org/)

Install all dependencies via:

```bash
pip install pysam numpy scipy pandas scikit-learn
```

---

## Usage

All modules are run through the main CLI script:

```bash
python ffy_calc/main.py --help
```

Modules are selected using **double-dash flags**, and each module has its own parameters.

### Modules

#### `--coverage-bed`

Generate a windowed coverage BED file from a BAM:

```bash
python main.py --coverage-bed --bam sample.bam --bin-size 100000 --mapping-quality 20
```

* `--bam` : Input BAM file
* `--bin-size` : Window size in base pairs (default 1000)
* `--mapping-quality` : Minimum mapping quality for reads (default 20)

Output: BED file with columns: `chr, start, end, coverage, fractional coverage, lowq fraction`.

---

#### `--predict-sex`

Predict the sex of the sample using chromosome Y coverage:

```bash
python main.py --predict-sex --coverage-bed-file sample.coverage.bed --reference-fasta ref.fa
```

* `--coverage-bed-file` : Input coverage BED file
* `--reference-fasta` : Indexed reference FASTA
* `--y-frac-threshold` : Threshold for calling male (default 0.02)
* `--lowq-threshold` : Maximum fraction of low-quality reads per bin (default 0.1)
* `--max-n-frac` : Maximum fraction of Ns per bin in reference (default 0.2)

Output: CSV table with predicted sex per sample.

---

#### `--generate-regions`

Generate robust Y-chromosome regions suitable for downstream analyses:

```bash
python main.py --generate-regions --sample-list predict_sex_outputs.txt --csv-out regions.csv
```

* `--sample-list` : Text file listing predict-sex output files
* `--female-median-threshold` : Maximum median chrY fraction for females (default 0.01)
* `--csv-out` : Output CSV summarizing bins per sample

Output: BED and CSV files with median/percentile coverage statistics and status (`ok` or `fail`).

---

#### `--model`

Train a linear regression model to predict fetal fraction (FFY):

```bash
python main.py --model --training-csv training.csv --regions-bed regions.bed --model-out ffy_model.json
```

* `--training-csv` : CSV file with BED path and fetal fraction
* `--regions-bed` : BED of robust regions from generate-regions
* `--model-out` : Output JSON file with linear regression model

Output: JSON model file and CSV table with predicted vs. observed FFY.

---

#### `--predict`

Predict FFY for a new sample using a trained model:

```bash
python main.py --predict --coverage-bed-file sample.coverage.bed --regions-bed regions.bed --model-file ffy_model.json
```

* `--coverage-bed-file` : Input coverage BED file
* `--regions-bed` : BED of robust regions
* `--model-file` : Trained FFY model JSON

Output: Tab-separated file with columns: BED file, chrY median fraction, predicted FFY.

---

#### `--version`

Print the current version of the tool:

```bash
python main.py --version
```

Output:

```
0.0.0
```

---

## License

This tool is released under the MIT License. See [LICENSE](LICENSE) for details.

