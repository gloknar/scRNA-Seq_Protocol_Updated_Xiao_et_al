#!/bin/sh
#SBATCH --job-name=Pathway-heterogeneity-Tumor-R       # Job name
#SBATCH --mail-type=BEGIN,END              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=acm95@ugr.es     # Where to send mail 
#SBATCH --nodes=1                    # Use one node
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --output=Pathway-Heterogeneity_dump_%j.out
#SBATCH --get-user-env 


/opt/R/4.1.2/bin/Rscript intra_malignant_Heterogeneity.R melanoma
/opt/R/4.1.2/bin/Rscript intra_malignant_Heterogeneity.R head_neck

/opt/R/4.1.2/bin/Rscript intra_non_malignant_Heterogeneity.R melanoma
/opt/R/4.1.2/bin/Rscript intra_non_malignant_Heterogeneity.R head_neck

# Analizamos las rutas metabólicas OXPHOS, glicólisis y respuesta a hipoxia
/opt/R/4.1.2/bin/Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R melanoma
/opt/R/4.1.2/bin/Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R head_neck

# Obtenemos los genes diferencialmente expresados en células tumorales con 
# baja actividad en las 3 rutas metabólicas de interés
/opt/R/4.1.2/bin/Rscript GeneSignature-of-OXPHOS_Glycolysis_Hypoxia.R melanoma
/opt/R/4.1.2/bin/Rscript GeneSignature-of-OXPHOS_Glycolysis_Hypoxia.R head_neck
