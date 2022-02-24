#!/bin/sh
#SBATCH --job-name=UMAP-tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=R-UMAP_dump_%j.out
#SBATCH --get-user-env 


# UMAP (Opcional)
/opt/R/4.1.2/bin/Rscript umap_metabolic_genes.R melanoma
/opt/R/4.1.2/bin/Rscript umap_metabolic_genes.R head_neck



