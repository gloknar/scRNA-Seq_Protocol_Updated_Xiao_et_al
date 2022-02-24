#!/bin/sh
#SBATCH --job-name=normalizacion-tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=normalization-R_dump_%j.out
#SBATCH --get-user-env 


# Normalizado de datos mediante distintos métodos
/opt/R/4.1.2/bin/Rscript normalization.R melanoma && 
/opt/R/4.1.2/bin/Rscript normalization.R head_neck




