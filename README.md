# SexDeterminer: An Automated Pipeline for Paleo-Proteomic Sex Estimation 

## Overview

SexDeterminer is a lightweight, open-source Nextflow pipeline designed for accurate, robust, and standardized sex determination from paleo-proteomic data. The pipeline uses a ratio-based method to calculate the proportion of male-specific amelogenin (AMELY) peptides relative to total amelogenin peptides detected, enabling reliable sex estimation even in degraded ancient samples where traditional morphological or DNA-based methods fail.

## Key Features

* Ratio-Based Sex Estimation: Uses R<sub>AMELY</sub> (n<sub>AMELY</sub>/(n<sub>AMELX</sub>+n<sub>AMELY</sub>)) ratio for classification instead of binary presence/absence detection of specific AMELY markers (e.g. 59M)

* Robust to False Positives: Statistical framework minimizes impact of false-positive AMELY signals

* Multi-Database Search Software Compatibility: Works with output from PEAKS, MaxQuant, pFind, and DIA-NN

* Pre-built Reference Databases: Includes amelogenin databases designed for different analysis contexts

* Automated Workflow: End-to-end processing from raw search results to final sex assessment report

## Usage 

We have provided an instruction on how to prepare input data and use the pipeline. Please refer to the `Manual.pdf` file in this repository.

## Example Data

We have provided the example raw data (`data/`), input files (`input/`), and output result files (`output/`) for reference under the `example` directory.

## License

Code released under GNU General Public License v3.0.

## Question & Bug

Please report questions, bugs, or any suggestions on **issues** page.