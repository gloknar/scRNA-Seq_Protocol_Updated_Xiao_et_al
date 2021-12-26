#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(Rtsne)
library(pheatmap)
library(RColorBrewer)
source("../utils.R")

# Opciones
options(stringsAsFactors=FALSE)
# argumento <- commandArgs()
# argumento <- args[6]
argumento <- "melanoma"

outDir <- file.path("./datasets",argumento)
if(!dir.exists(outDir)) {                    # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
  dir.create(outDir,recursive=TRUE)
}



# Leemos el dataset del head_neck/melanoma con la expresión génica imputada y de
# ahí generamos un objeto sce que contenga sólo las células tumorales y los
# genes de interés (en este caso, los metabólicos)
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
tumor_metabolic_sce <- imputed_sce[rowData(imputed_sce)$metabolic, imputed_sce$cellType == "Malignant"]

# Limpieza RAM
rm(imputed_sce)
gc(verbose = F)



####################################################################################################

########################################################################
###########     1. Calculo de la matriz de correlación       ###########
########################################################################

expr_genes_metab_tumor <- assay(tumor_metabolic_sce, "exprs")  # Cogemos la expresión génica en formato log2(TPM+1)
matriz_cor  <- cor(expr_genes_metab_tumor, method = "spearman") # Calculamos la matriz de correlación de spearman (no paramétrico)
??corplot
psych::cor.plot(matriz_cor[1:4,1:4])

hc <- hclust(as.dist(1-dist_dat),method="ward.D2")
mycolor = colorRampPalette(c("white","red"))(10)
mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)
col_row <- as.data.frame(colData(tumor_metabolic_sce)[,"argumento",drop=F])
col_row$argumento <- factor(col_row$argumento)
pdf(file.path(outDir,"malignant_metabolic_correlationMatrix.pdf"),width=5,height=4,onefile=T)
pheatmap(dist_dat,show_rownames = F,show_colnames = F,cluster_rows = hc,cluster_cols =hc,annotation_row = col_row,color=mycolor, annotation_legend = T)
dev.off()

####### For melanoma argumento, move MEL3 to after MEL12.
####### 
if(argumento == "melanoma"){
  order_names <- hc$labels[hc$order]
  MEL12_names <- rownames(col_row)[col_row$argumento=="MEL12"]
  order_names2 <- MEL12_names
  MEL3_names <- rownames(col_row)[col_row$argumento=="MEL3"]
  order_names2 <- c(order_names2,MEL3_names)
  other_names <- order_names
  other_names <- order_names[!(order_names %in% MEL3_names)]
  other_names <- other_names[!(other_names %in% MEL12_names)]
  order_names2 <- c(order_names2,other_names)
  dist_dat2 <- dist_dat[order_names2,order_names2]
  pdf(file.path(outDir,"malignant_metabolic_correlationMatrix2.pdf"),width=5,height=4,onefile=T)
  pheatmap(dist_dat2,show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols =F,annotation_row = col_row,color=mycolor, annotation_legend = T)
  dev.off()
  
}


