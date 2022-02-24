#!/bin/sh
#SBATCH --job-name=scRNA_Pathway_Activity-tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=scRNA_Pathway_Activity-R_dump_%j.out
#SBATCH --get-user-env 

# Medidos actividad de rutas metabólicas en cada tipo celular presente en los datasets 
/opt/R/4.1.2/bin/Rscript scRNA_pathway_activity.R melanoma
/opt/R/4.1.2/bin/Rscript scRNA_pathway_activity.R head_neck
