# DEPMAP24Q4 BAGEL Analysis

## Summary

This repository contains a workflow used to generate BAGEL essentiality scrores for the latest version of DepMap (24Q4) for comparing essentiality scores from Pan-Cancer samples with Acral melanomas. In brief, the DepMap data was reformatted in R and then corrected for copy number variation using [CRISPRcleanR](https://github.com/francescojm/CRISPRcleanR). The resulting copy-number corrected counts were then used as inputs for BAGEL. Finally, the BAGEL essentiality scores were aggregated to generate a matrix of essentiality scores for each cell line.

Running the scripts in this repository should retrieve and process the data from end-to-end.

## Dependencies
- R (4.4.0) - libraries recorded in `renv.lock` file
- Nextflow (>23.10.0)
- Singularity (3.14)

All R scripts were run on a Ubuntu 22.04.3 LTS system with R 4.4.0 installed via [rig](https://github.com/r-lib/rig). The Nextflow pipeline was run on Wellcome Trust Sanger Institutes's high-performance computing (HPC) cluster running Ubuntu 20.04 (Focal Fossa) with Singularity 3.14 installed and was configured to use the LSF scheduler. 

Runtime for each R script on a virtual machine with 8 cores is several minutes - with the execption of `scripts/03_apply_crispr_cleanr.R` which will take ~12 hours. Runs of the nextflow pipeline (`scripts/04_run_bagel_pipeline.sh`) requires approximately 70 hours wall-time using 3CPUs and 5GB of memory per BAGEL task. The pipeline is designed to run on a cluster with LSF scheduler, but can be modified to run anyhwhere by creating a profile in the `nextflow.config` file.

## Workflow steps

1) Counts data from DepMap 24Q4 was retrieved from the [DepMap Figshare repository](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248) (`scripts/01_fetch_depmap_data.R`).
2) Samples names in the counts matrix were updated according to DepMap model ID, replicate ID and Avana pDNA batch (`scripts/02_reformat_depmap_data.R`).
3) The counts matrix was split by cell line and CRISPRcleanR was run on each cell line (`scripts/03_apply_crispr_cleanr.R`). Guides with fewer than 30 counts in the plasmid sample that was used to transfect a cell line were removed at this stage. 
4) Steps from BAGEL (Fold change calculation, Bayes Factor estimation and Precision Recall calculation) were run on each cell line dataset (`scripts/04_run_bagel_pipeline.sh`) using a custom nextflow pipeline (`scripts/main.nf`)
5) The resulting BAGEL essentiality scores were aggregated into a single matrix (`scripts/05_aggregate_data.R`). Where there were instances of screens for a cell line having been perfomed with mutliple pDNA batches only the results from one pDNA batch were retained to avoid double-counting the same cell line in downstream analyses.

## Directory structure
```
└── scripts
    ├── 01_fetch_depmap_data.R
    ├── 02_reformat_depmap_data.R
    ├── 03_apply_crispr_cleanr.R
    ├── 04_run_bagel_pipeline.sh
    ├── 05_aggregate_data.R
    ├── bin
    ├── main.nf
    └── nextflow.config
```

## Contact
- Jamie Billington (jb63@sanger.ac.uk)