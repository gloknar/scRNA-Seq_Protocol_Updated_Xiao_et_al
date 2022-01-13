#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(pheatmap)
library(gtools)
source("../utils.R")

# Opciones
options(stringsAsFactors = FALSE)
# argumento <- commandArgs()
# argumento <- argumento[6]
argumento <- "melanoma"
outDir <- file.path("./datasets",argumento,"low-OXPHOS-gly-hyp-activity-genes")
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}


# Leemos el dataset del head_neck/melanoma
filtered_sce <- readRDS(file.path("../1-ReadData/datasets/",argumento,"filtered_sce.rds"))


tumor_sce <- filtered_sce[,filtered_sce$cellType == "Malignant"]
tumor_metabolic_sce <- tumor_sce[rowData(tumor_sce)$metabolic,]

tumor_sce$tumor <- factor(tumor_sce$tumor)
tumores <- unique(tumor_sce$tumor)


# Leemos el archivo de las rutas en las que participan los 1566 genes
# metabólicos (este contiene 85 rutas metabólicas)
ruta_archivo_pathways <- "../Data/KEGG_metabolism.gmt"
pathways <- gmtPathways(ruta_archivo_pathways)

# Este contiene 50 rutas metabólicas (parece que contiene los genes que
# participan en respuesta a hipoxia)
hallmark_gmt <- '../Data/h.all.v6.1.symbols.gmt'
hallmarks <- gmtPathways(hallmark_gmt)



####################################################################################################

##############################################
###########     1. PLACEHOLDER     ###########
##############################################

#2.Tumor cells
all_low_cells <- c()
all_high_cells <- c()


for(t in tumores){
 # Usando los genes metabólicos
 each_metabolic_sce <- tumor_metabolic_sce[,tumor_metabolic_sce$tumor == t] # sce con células malignas que pertenecen al tumor t y los 1566 genes metabólicos
 each_metabolic_tpm <- assay(each_metabolic_sce, "exprs")
 each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm) > 0,] # Eliminamos los genes con dropout del 100%
 
 # Usando todos los genes
 each_tumor_sce <- tumor_sce[,tumor_sce$tumor == t]
 each_tumor_tpm <- assay(each_tumor_sce, "exprs")
 each_tumor_tpm <- each_tumor_tpm[rowSums(each_tumor_tpm) > 0,]
 
 oxphos_genes <- intersect(pathways[["Oxidative phosphorylation"]], rownames(each_tumor_tpm))
 glycolysis_genes <- intersect(pathways[["Glycolysis / Gluconeogenesis"]], rownames(each_tumor_tpm))
 hypoxia_genes <- intersect(hallmarks[["HALLMARK_HYPOXIA"]], rownames(each_tumor_tpm))
 
 three_all <- unique(c(oxphos_genes,glycolysis_genes,hypoxia_genes))

 oxphos_mean_exprs <- colMeans(each_tumor_tpm[three_all,], na.rm = T)
 
 oxphos_mean_exprs_quantile <- quantile(oxphos_mean_exprs,seq(0,1,0.2))
 low_cutoff <- oxphos_mean_exprs_quantile[["20%"]]
 high_cutoff <- oxphos_mean_exprs_quantile[["80%"]]
 oxphos_low <- which(oxphos_mean_exprs < low_cutoff)
 oxphos_high <- which(oxphos_mean_exprs > high_cutoff)
 
 all_low_cells <- c(all_low_cells, colnames(each_tumor_tpm)[oxphos_low])
 all_high_cells <- c(all_high_cells, colnames(each_tumor_tpm)[oxphos_high])

 if(length(oxphos_low) < 5 | length(oxphos_high) < 5){
   # Meter aquí un message
   next
 }
 condition <- factor(c(rep("oxphos_low",length(oxphos_low)),rep("oxphos_high",length(oxphos_high))),levels = c("oxphos_low","oxphos_high"))
 
 each_tumor_tpm_selected <- each_tumor_tpm[,c(oxphos_low,oxphos_high)] 
 
 pvalues <- sapply(X = 1:nrow(each_tumor_tpm_selected),
                   FUN = function(x) {
                     return(wilcox.test(each_tumor_tpm_selected[x,] ~ condition, alternative="greater")$p.value)
                   })
 pvalues_df <- data.frame(pvalues,row.names=rownames(each_tumor_tpm_selected))
 
 # Guardamos aquellos genes con un p-valor < 0.01 en un archivo
 write.table(rownames(pvalues_df[pvalues_df$pvalues<=0.01,,drop=F]),
             file=file.path(outDir,paste0(t,"_low_OXPHOS_glycolysis_hypoxia_signature.genes.txt")),
             quote=F,row.names=F,col.names=F)
}
##all together
selected_tumor_tpm <- assay(tumor_sce, "exprs")
selected_tumor_tpm <- selected_tumor_tpm[rowSums(selected_tumor_tpm) > 0,]


condition <- factor(c(rep("oxphos_low",length(all_low_cells)),rep("oxphos_high",length(all_high_cells))),levels = c("oxphos_low","oxphos_high"))

selected_tumor_tpm_selected <- selected_tumor_tpm[,c(all_low_cells,all_high_cells)] 

pvalues <- sapply(X = 1:nrow(selected_tumor_tpm_selected),
                  FUN = function(x) {
                    return(wilcox.test(selected_tumor_tpm_selected[x,] ~ condition, alternative="greater")$p.value)
                  })
pvalues_df <- data.frame(pvalues,row.names=rownames(selected_tumor_tpm_selected))


# Guardamos aquellos genes con un p-valor < 0.01 en un archivo
write.table(rownames(pvalues_df[pvalues_df$pvalues<=0.01,,drop=F]),
            file=file.path(outDir,paste0("ALL","_low_OXPHOS_glycolysis_hypoxia_signature.txt")),
            quote=F,row.names=F,col.names=F)
