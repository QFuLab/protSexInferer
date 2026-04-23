# protSexInferer: An Automated Pipeline for Paleo-Proteomic Sex Estimation 

## Overview

protSexInferer is a lightweight, open-source Nextflow pipeline designed for accurate, robust, and standardized sex determination from paleo-proteomic data. The pipeline uses a ratio-based method to calculate the proportion of male-specific amelogenin (AMELY) peptides relative to total amelogenin peptides detected, enabling reliable sex estimation even in degraded ancient samples where traditional morphological or DNA-based methods fail.

## Citation

Fan Bai, Zhongyou Wu, Song Xing, Qiaomei Fu, Rapid and robust sex determination from ancient enamel proteomes using protSexInferer, _Journal of Genetics and Genomics_, 2026, ISSN 1673-8527, https://doi.org/10.1016/j.jgg.2026.04.012. 
(https://www.sciencedirect.com/science/article/pii/S167385272600144X)

## Key Features

* Ratio-Based Sex Estimation: Uses R<sub>AMELY</sub> (n<sub>AMELY</sub>/(n<sub>AMELX</sub>+n<sub>AMELY</sub>)) ratio for classification instead of binary presence/absence detection of specific AMELY markers (e.g. 59M)

* Robust to False Positives: Statistical framework minimizes impact of false-positive AMELY signals 

* Multi-Database Search Software Compatibility: Works with output from PEAKS, MaxQuant, pFind, and DIA-NN

* Generalizable Pre-built Reference Databases: Includes amelogenin databases designed for different analysis contexts

* Automated Workflow: End-to-end processing from raw LC-MS/MS data to final sex assessment report

## Usage 

We have provided an instruction on how to prepare input data and use the pipeline. Please refer to the `Manual.pdf` file in this repository.

## Example Data

We have provided the example raw data (`data/`), input files (`input/`), and output result files (`output/`) for reference under the `example` directory.

## License

Code released under GNU General Public License v3.0.

## Question & Bug

Please report questions, bugs, or any suggestions on [issues](https://github.com/QFuLab/protSexInferer/issues) page.