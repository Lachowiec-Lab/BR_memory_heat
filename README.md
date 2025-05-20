# Data and Code from "Somatic memory of exogenous brassinosteroids alters histone expression and tempers heat responses"

This repository contains data and code supporting the findings of the study 
on the impacts of seed applied brassinosteroids on heat responses in spring wheat reproductive development.

We collected transcriptomic data on developing spikes post-BR application and heat stress.

## Description of Repository Contents

### Data

Data are in the `/data` directory:

* `metadata.xlsx`: Codebook providing detailed descriptions of RNA-seq samples. 
* `RNA-seq-counts.txt`: This is the raw sequencing counts data for all samples.
* `expDesign7.txt`: This includes the experimental treatment design for use with DESeq2
* `GOterms.csv.zip`: This contains a zipped file of the GO terms associated with the IWGSC Chinese Spring genome v1.1

### Code

Code used to process data and generate the manuscript's RNA-seq analysis and figures.

* `heatEffectsOnly.R`: R Script for analyzing the effects of heat of all samples and over time. This generates figures 3B, S3, S4.
* `BReffects.R`: R script for examining the BR effect prior to heat treatment. This generates figures 4.
* `16hTimepoint.R`: R script for examining the 16hr timepoint--specifically contrasting the mock at 29deg compared to three other treatments: mock 22deg, BR 22deg, and BR 29deg


## Citing this work

This repository contains data and code to support the submitted manuscript:

> Kothari A, Correr FH, Hinson C, King A, Astorga Bedoya A, Erwin D, Cook JP, Lachowiec J. (2025). Somatic memory of exogenous brassinosteroids alters histone expression and tempers heat responses.

And consider contributing cleaned data and code to this repository.

## Acknowlegements and Support

**Acknowledgments**

Computational efforts were performed on the Tempest High Performance Computing System, operated and supported by University Information Technology Research Cyberinfrastructure at Montana State University. 

**Funding**

This research was supported by the
U.S. Department of Agriculture, National Institute of Food and
Agriculture (2020-65114-30768) and the Montana Wheat and Barley Committee. 
