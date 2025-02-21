#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000
#BSUB -oo /lustre/scratch127/casm/team113da/users/jb63/bagel_depmap_24Q4/pipeline_%J.o
#BSUB -eo /lustre/scratch127/casm/team113da/users/jb63/bagel_depmap_24Q4/pipeline_%J.e
module load nextflow 
module load singularity
nextflow run main.nf -resume \
-c nextflow.config \
-profile farm22