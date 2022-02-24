#!/bin/sh
#SBATCH --job-name=imputacion-tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=imputation_dump_%j.out
#SBATCH --get-user-env 

# Imputamos
/opt/R/4.1.2/bin/Rscript impute_tpm.R melanoma 20
/opt/R/4.1.2/bin/Rscript impute_tpm.R head_neck 20


# Comprobamos los resultados de la imputación
/opt/R/4.1.2/bin/Rscript imputation_plots.R melanoma
/opt/R/4.1.2/bin/Rscript imputation_plots.R head_neck
