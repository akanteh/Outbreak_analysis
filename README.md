# Transmission dynamics of SARS-CoV-2 infections in a hospital in The Gambia

## Overview

This repository contains analysis scripts and manuscript outputs
associated with the study:

**Transmission dynamics of SARS-CoV-2 infections in a hospital in The
Gambia across two waves during the COVID-19 pandemic --- a retrospective
analysis of systematic surveillance data**

The study retrospectively investigated probable within-hospital
SARS-CoV-2 transmission at the Clinical Services Department (CSD) of the
Medical Research Council Unit The Gambia at the London School of Hygiene
& Tropical Medicine (MRCG at LSHTM). The analysed surveillance period
was July 2020 to March 2021 and covered the first epidemic wave
(July--September 2020) and the portion of the second national wave
captured by study surveillance (January--March 2021).

The primary transmission analysis used **Outbreaker2** to integrate:

-   symptom-onset dates;
-   pairwise genetic distances where genomic data were available; and
-   staff/patient location information within the hospital.

The inferred transmission links are probabilistic and should not be
interpreted as direct proof of person-to-person transmission.

## Study summary

A total of 1,414 RT-PCR tests were performed during the study period, of
which 130 were positive. Sixty-three hospital-associated cases met the
analysis criteria, and 20 (31.7%) generated high-quality, near-complete
SARS-CoV-2 genomes.

Three Outbreaker2 model specifications were explored:

1.  symptom-onset dates only;
2.  symptom-onset dates plus available genetic distances; and
3.  symptom-onset dates, available genetic distances, and staff/patient
    location information.

The combined model was used for the primary descriptive reconstruction.

For the first wave, the model reconstructed a major transmission chain
with substantial staff-associated transmission. For January--March 2021,
the reconstruction contained multiple introductions and a greater
relative contribution of patient-associated links. Interpretation is
limited by incomplete genomic sampling, uncertainty in symptom-onset
dates, unsampled infections, and differences in surveillance intensity
between waves.

## Repository structure

``` text
Outbreak_analysis/
├── README.md
├── .gitignore
├── scripts/
│   ├── outbreaker2AnalysisPipeline.R
│   └── TransmissionNetworkTree.R
├── results/
│   ├── figures/
│   │   ├── Figure1.pdf
│   │   ├── Figure2.pdf
│   │   ├── Figure3.pdf
│   │   ├── Figure4.pdf
│   │   ├── Figure5.pdf
│   │   └── Figure6.pdf
│   └── tables/
│       ├── TableS1.pdf
│       └── TableS2.pdf
└── docs/
    └── supplementaryFile.pdf
```

### `scripts/`

`outbreaker2AnalysisPipeline.R` contains the main Outbreaker2 workflow,
including data preparation, onset-date processing, sequence import,
pairwise SNP-distance calculation, generation-time and incubation-period
distributions, construction of the location/contact matrix, and
Outbreaker2 model runs.

`TransmissionNetworkTree.R` contains code for constructing and plotting
a directed transmission network from model-derived infector--infectee
links. Edge width/opacity represents posterior support and square nodes
identify inferred imports.

### `results/figures/`

Contains the six main manuscript figures in PDF format.

### `results/tables/`

Contains supplementary tables, including the Outbreaker2 configuration
and the descriptive distribution of reconstructed transmission-event
categories.

### `docs/`

Contains the combined supplementary material.

## Data and input requirements

The individual-level epidemiological and genomic input data are **not
included in this GitHub repository**. The analysis script currently
expects local input files including:

``` text
data/metaData.csv
analysis/oubreakerSamplesfirstwave.fasta
```

The metadata used by the current script include, at minimum, variables
corresponding to:

-   `lab_no` --- case/sample identifier;
-   `DateOfCollection` --- sample collection date;
-   `OnsetDate` --- symptom-onset date; and
-   `locations` --- hospital location/work area used to construct the
    contact matrix.

The analysis also requires helper R scripts referenced by
`outbreaker2AnalysisPipeline.R`:

``` text
src/functions_wards.R
src/functions.R
src/onset_distribution.R
src/functions_analyse.R
```

These helper files are referenced by the current analysis script but are
**not present in the repository snapshot used to prepare this README**.

## Genomic sequencing and primary analysis

Positive samples with Ct values ≤30 were selected for sequencing.
Libraries were prepared using the ARTIC Network V3 protocol and
sequenced using Oxford Nanopore GridION with R9.4.1 flow cells.

The manuscript reports the following primary genomic-analysis software:

-   Guppy 6.5.7 --- basecalling;
-   ARTIC bioinformatics pipeline 1.5.1 --- consensus genome generation;
-   Pangolin 4.3 (data version 1.22) --- lineage assignment;
-   Nextclade 3.10.0 --- sequence quality assessment;
-   MAFFT 7.490 --- multiple sequence alignment; and
-   AliView --- manual alignment inspection.

The SARS-CoV-2 reference genome used for mapping was **MN908947.3**.

Twenty high-quality (\>90% coverage), near-complete genomes were
included in the final genomic component of the analysis.

## Outbreaker2 model

The transmission reconstruction was performed in R using Outbreaker2,
which implements Bayesian MCMC inference of transmission chains.

### Epidemiological distributions

The current analysis script specifies:

``` r
incub <- distcrete::distcrete(
  "gamma",
  shape = 5.89,
  scale = 1.068555,
  interval = 1
)

gen <- distcrete::distcrete(
  "gamma",
  shape = 4.83,
  scale = 1.051851,
  interval = 1
)
```

These represent the incubation-period and generation-time distributions,
respectively.

### Main MCMC settings

The manuscript reports:

``` r
n_iter       = 1000000
sample_every = 500
```

The first 500 sampled iterations were discarded as burn-in, and
convergence was assessed using trace plots.

Supplementary Table 1 reports the following additional configuration:

  Parameter        Value
  ---------------- ---------------------------------------
  `n_iter`         1,000,000
  `sample_every`   500
  `move_mu`        TRUE; `init_mu = 0.2`, `sd_mu = 0.05`
  `move_pi`        TRUE; `prior_pi = c(5, 20)`
  `move_eps`       FALSE
  `move_tau`       FALSE
  `move_model`     FALSE
  `sd_t_inf`       3

Posterior probability ≥0.6 was used in the manuscript as the threshold
for primary interpretation of reconstructed transmission links.

## R dependencies

Packages explicitly loaded or used by the scripts include:

``` text
ape
chron
distcrete
epicontacts
epitrix
fitdistrplus
ggplot2
ggraph
glue
here
igraph
linelist
lubridate
magrittr
outbreaker2
pacman
pairsnp
remotes
rio
tidyverse
visNetwork
```

The main script also installs a development branch of Outbreaker2:

``` r
remotes::install_github("finlaycampbell/outbreaker2@timeline")
```

The manuscript reports that data management and statistical analyses
were performed using **R 4.2.1** and **RStudio 2022.7.1.554**.

## Running the analysis

> **Important:** the current repository is an archive of the analysis
> code and outputs and is not yet fully portable.

Before running `scripts/outbreaker2AnalysisPipeline.R`, the following
steps are required:

1.  Obtain authorised access to the study metadata and sequence data.
2.  Place the required input files in the paths expected by the script,
    or modify the paths accordingly.
3.  Add the helper scripts listed under **Data and input requirements**.
4.  Install the required R packages.
5.  Replace the hard-coded `setwd()` path in
    `outbreaker2AnalysisPipeline.R` with the local repository path or
    refactor the workflow to use relative paths.
6.  Confirm that metadata identifiers and FASTA sequence names are in
    the same order expected by the analysis.
7.  Run the Outbreaker2 models and assess MCMC convergence before
    interpreting reconstructed links.

The current script contains analysis objects and commands from iterative
model development. Users should therefore verify object names, model
settings, wave-specific inputs, and output paths before attempting a
complete rerun.

## Transmission-network visualisation

`TransmissionNetworkTree.R` uses `igraph` and `ggraph` to visualise
model-derived transmission links.

The script expects a table with the following fields:

``` text
from
to
time
support
generations
```

Rows with `from = NA` are treated as inferred imports. Directed edges
represent model-inferred infector-to-infectee relationships, and edge
width/opacity represents posterior support.

The script writes:

``` text
transmission_network.png
transmission_network.pdf
```

## Main outputs

The repository currently contains six manuscript figures under
`results/figures/`.

The supplementary material includes:

-   **Supplementary Table 1:** Outbreaker2 configuration and priors;
-   **Supplementary Table 2:** descriptive frequency of reconstructed
    transmission-event categories;
-   **Supplementary Figure 1:** distribution of samples tested and
    included in the analysis;
-   **Supplementary Figure 2:** symptom-onset-only Outbreaker2 model
    outputs; and
-   **Supplementary Figure 3:** Outbreaker2 model incorporating symptom
    onset, genetic data, and staff/patient location.

## Data availability and confidentiality

The manuscript states that the sequences, scripts, and metadata used to
generate the findings are deposited on Zenodo under **record 20680779**,
with access restricted to authorised persons.

The underlying individual-level epidemiological data should not be
committed to this public repository unless appropriate ethical,
institutional, and data-governance permissions explicitly allow public
release.

The study was approved by the MRCG Scientific Coordinating Committee and
The Gambia Government/MRC Joint Ethics Committee (reference **28200**).
The manuscript describes the work as a retrospective analysis using
anonymised clinical and genomic data.

## Interpretation and limitations

The reconstructed transmission networks should be interpreted as
**probabilistic model-based reconstructions**, not definitive evidence
of who infected whom.

Important limitations include:

-   high-quality genomes were available for only 20 of 63 analysed
    cases;
-   most links therefore relied substantially on symptom-onset and
    location information;
-   symptom-onset dates were imputed for some cases;
-   unsampled infections and shared exposures may produce alternative
    plausible networks;
-   surveillance intensity and case ascertainment differed between
    waves; and
-   study surveillance ended in March 2021 although the national second
    wave continued to May 2021.

These limitations should be considered when reproducing or reusing the
analysis.

## Citation

If using this repository, please cite the associated manuscript once the
final publication details are available:

> Kanteh A, Forrest K, Gomez ML, Campbell F, Mohammed N, Amambua-Ngwa A,
> Nadjm B, Roca A, Sesay AK, Hué S. *Transmission dynamics of SARS-CoV-2
> infections in a hospital in The Gambia across two waves during the
> COVID-19 pandemic --- a retrospective analysis of systematic
> surveillance data.*

The formal journal citation and DOI should be added here after
publication.

## Authors

**Abdoulie Kanteh** --- analysis and manuscript lead

The manuscript lists Karen Forrest, Marie Louise Gomez, Finlay Campbell,
Nuredin Mohammed, Alfred Amambua-Ngwa, Behzad Nadjm, Anna Roca, Abdul
Karim Sesay, and Stéphane Hué as co-authors.

## Licence

No software or data licence was specified in the repository materials
reviewed for this README. A licence should be selected before the
repository is made publicly reusable. The licence for the analysis code
should be considered separately from permissions governing the clinical,
epidemiological, and genomic data.
