#!/bin/sh
#SBATCH --job-name=Correlation-interTumor-R       # Job name
#SBATCH --mail-type=BEING,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=Correlation-interTumor_dump_%j.out
#SBATCH --get-user-env 


# Matriz de correlación de Spearman
Rscript inter_tumor_correlation_matrix.R melanoma &&
Rscript inter_tumor_correlation_matrix.R head_neck





