#!/bin/sh
#SBATCH --job-name=GSEA-Tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=GSEA_dump_%j.out
#SBATCH --get-user-env 

# Corremos GSEA en distintos tipos celulares no cancerosos
/opt/R/4.1.2/bin/Rscript non_malignant_subtype.R melanoma
/opt/R/4.1.2/bin/Rscript non_malignant_subtype.R head_neck
