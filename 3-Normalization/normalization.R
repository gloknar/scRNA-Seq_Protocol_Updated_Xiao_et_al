#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes
library(scater)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(viridis)
library(edgeR)
library(scran)


# Opciones
options(stringsAsFactors = FALSE)
# argumento <- commandArgs()
# argumento <- argumento[6]
argumento = "head_neck"
outDir <- file.path("./datasets", argumento)
if(!dir.exists(outDir)) {                 # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
  dir.create(outDir,recursive = TRUE)
}


# Leemos el dataset del head_neck/melanoma con la expresión génica imputada y el
# dataset con las células filtradas. Cargamos también los tipos celulares
# presentes en dichos dataset
imputed_sce <- readRDS(file.path("../2-Imputation/datasets", argumento, "imputed_sce.rds"))
imputed_sce$cellType <- droplevels(imputed_sce$cellType)
linajes_celulares <- unique(imputed_sce$cellType)



####################################################################################################

####################################################################
###########     1. Selección de genes de referencia      ###########
####################################################################

# Para normalizar la expresión de nuestros genes, usaremos como referencia
# aquellos que tengan una expresión elevada y sean detectados en 25% o más de
# las células (tasa de dropout < 75%)

limite_dropout <- 0.75

# Creamos una matriz de genes X tipos celulares llena de FALSE, que usaremos de
# plantilla para seleccionar, para cada linaje celular, los genes con un dropout
# < 75%
matriz_seleccion_genes <- matrix(FALSE,
                                 nrow = nrow(imputed_sce),
                                 ncol = length(linajes_celulares),
                                 dimnames = list(rownames(imputed_sce),linajes_celulares))

# Iteramos sobre cada tipo celular
for(c in linajes_celulares){
  sce_por_tipo_celular <- imputed_sce[,imputed_sce$cellType == c] # hacemos un subset para cada tipo celular del objeto sce con la expresión génica imputada
  matriz_genica_log_linaje_celular <- assay(sce_por_tipo_celular,"exprs")  # Calculamos la expresión en log2(TPM+1) de los genes de las células seleccionadas
  tasa_deteccion_gen <- apply(matriz_genica_log_linaje_celular,1, function(x) sum(x>0)/ncol(matriz_genica_log_linaje_celular))  # Calculamos para cada gen el % de células con expresión génica no nula
  genes_seleccionados <- tasa_deteccion_gen >= limite_dropout
  matriz_seleccion_genes[genes_seleccionados,c] <- TRUE
}

print("the number of genes selected:")
print(sum(rowSums(matriz_seleccion_genes) >= length(linajes_celulares)))
low_dropout_genes <- rownames(matriz_seleccion_genes)[rowSums(matriz_seleccion_genes) >= length(linajes_celulares)]

##########################################################################################
## different normalization methods:
##########################################################################################
#prepare the gene length file
#convert the TPM to count scale
all_gene_lengths <- read.table(file.path("../Data/","gene_length.txt"),sep="\t",header=F,row.names=1)
tmp <- intersect(rownames(all_gene_lengths),rownames(imputed_sce))
if (length(tmp) != nrow(imputed_sce)){
  warning("check the length file")
}
genelen <- all_gene_lengths[rownames(imputed_sce),]
genelen <- as.numeric(as.vector(genelen))

#1. up-quantile normalization
selected_impute_tpm <- tpm(imputed_sce)
selected_impute_counts <- sweep(selected_impute_tpm, 1, genelen, FUN = "*")

median.sf <- calcNormFactors(selected_impute_counts[low_dropout_genes,],method="upperquartile",p=0.75)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / median.sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path(outDir,"UpperQuartile_tpm.rds"))


#2.DESeq2 normalization
# get the function from "https://github.com/mikelove/DESeq2/blob/master/R/core.R"
source("../utils.R")
deseq2_sf <- estimateSizeFactorsForMatrix(selected_impute_counts[low_dropout_genes,])
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / deseq2_sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path(outDir,"RLE_tpm.rds"))

#3. edgeR (TMM)
edgeR_sf <- calcNormFactors(selected_impute_counts[low_dropout_genes,],method="TMM")
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / edgeR_sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm + 1)
#save
saveRDS(selected_impute_tpm_norm,file.path(outDir,"TMM_tpm.rds"))


# NOTA: Los normalizamos por TMM, RLE y up-quantile me salen iguales que en el
# paper, pero el de deconvolución NO!!!!! MIRAR QUÉ PASA AHÍ
#4. scran, deconvolution

# Parece que la solución era usar una funcion que sí admite matrices... o sea,
# usar "calculateSumFactors" en lugar de "computeSumFactors"
scran.sf <- scran::calculateSumFactors(selected_impute_counts[low_dropout_genes,],clusters=imputed_sce$cellType)
# scran.sf <- scran::computeSumFactors(selected_impute_counts[low_dropout_genes,],clusters=imputed_sce$cellType)
summary(scran.sf)
selected_impute_tpm_norm <- t(t(selected_impute_tpm) / scran.sf)
selected_impute_exp_norm <- log2(selected_impute_tpm_norm+1)
#save
saveRDS(selected_impute_tpm_norm,file.path(outDir,"Deconvolution_tpm.rds"))


###Evaluation of the normalization methods:
all_cell_type <- as.vector(imputed_sce$cellType)
for(m in c("RLE","TMM","UpperQuartile","Deconvolution"))
{
  rds_file <- file.path(outDir,paste0(m,"_tpm.rds"))
  norm_tpm <- readRDS(rds_file) 

  ##### check the ratio distribution
  low_dropout_genes_tpm <- norm_tpm[low_dropout_genes, ]
  low_dropout_genes_tpm_mean <- apply(low_dropout_genes_tpm, 1, function(x) by(x, all_cell_type, mean))
  low_dropout_genes_tpm_ratio <- t(low_dropout_genes_tpm_mean) / colMeans(low_dropout_genes_tpm_mean)
  dat <- reshape2::melt(low_dropout_genes_tpm_ratio)
  p <- ggplot(dat,aes(x=Var2,y=value)) +
    geom_boxplot(outlier.alpha=0.1)+ theme_classic() + 
    ylab("expression ratio") + xlab("") +    
    theme(axis.text.x = element_text(angle=45,hjust=1))
  ggsave(file.path(outDir,paste0(m,"_ratio_distribution.pdf")),p,width=3.5,height=2.5)
}