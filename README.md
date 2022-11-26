
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Genvej

### Read Trimmer for Next Generation Sequencing Data

<!-- badges: start -->
<!-- badges: end -->
<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

`Genvej` is Read Trimmer for Next Generation Sequencing Data program
written in Python 3.

`Genvej` provides a consistent set of verbs that help you trim reads
from FASTQ files.

These files may be compressed, uncompressed, encoded in with PHRED 64 or
PHRED 33 quality scores.

``` bash
python3 project_test_.py
```

------------------------------------------------------------------------

### Method 1

-   `--filename` is used before the name of the file is entered

-   `--shear` Trims 3’ prime end of a read by (x) amount of nucleotides

-   `--snip` Trims 5’ prime end of a read by (x) amount of nucleotides

``` bash
python3 project_test_.py --filename PHRED.gz --shear 5 --snip 5 
```

------------------------------------------------------------------------

### Method 2

-   `--filename` is used before the name of the file is entered

-   `--trim_minimum` Trims reads based on a minimum quality score using
    PHRED 33 only.

``` bash
python3 project_test_.py --filename BRISCOE.fastq --trim_minimum 20
```

------------------------------------------------------------------------

-   `--name` Receive a friendly welcome message. This is used before
    your name is entered.

``` bash
python3 project_test_.py --name 
```

-   `--version` Display the version of the program

``` bash
python3 project_test_.py --version 
```

-   `--cite` Display the citation for the paper associated with the
    program

``` bash
python3 project_test_.py --cite
```

-   Learn about the functions in the help section of the program
    `--help`.

``` bash
python3 project_test_.py --h 
```

If you are new to `Genvej`, the best place to start is this github
repository.

------------------------------------------------------------------------

## Installation

Download the python code from Github.

The link to the Github repository is :

<https://github.com/gavinakajiawen/Read-Trimmer-for-NGS-data>

------------------------------------------------------------------------

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on GitHub. For questions and other discussion,
please email `s172084@dtu.dk` or `s222372@student.dtu.dk`

Please note that this project is released with a Contributor Code of
Conduct. By participating in this project you agree to abide by its
terms.
