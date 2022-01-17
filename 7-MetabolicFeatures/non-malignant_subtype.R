#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(pheatmap)
source("../utils.R")
source("runGSEA.R")  # Pasar esta función a utils o cambiarle el nombre y revisar los nombres de los argumentos


# Opciones
options(stringsAsFactors = FALSE)
argumentos <- commandArgs()
tejido <- argumentos[6]
# tejido <- "melanoma"
outDir <- file.path("datasets",tejido)
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}


# Leemos el dataset del head_neck/melanoma y a partir de él creamos uno sólo
# para las células sanas y los genes metabólicos y otro para las células sanas y
# todos los genes
filtered_sce <- readRDS(file.path("../1-ReadData/datasets/",tejido,"filtered_sce.rds"))
filtered_sce$cellType <- factor(filtered_sce$cellType) # Casteamos a factor para eliminar niveles sobrantes

healthy_sce <- filtered_sce[, filtered_sce$cellType != "Malignant"] # Este ya viene con el factor arreglado, no hace falta tocarlo
healthy_metabolic_sce <- filtered_sce[rowData(filtered_sce)$metabolic, filtered_sce$cellType != "Malignant"]
tipos_celulares <- unique(healthy_metabolic_sce$cellType)

# Limpieza RAM
rm(filtered_sce)
invisible(gc())


# Leemos el archivo de las rutas en las que participan los 1566 genes
# metabólicos (este contiene 85 rutas metabólicas)
ruta_archivo_pathways <- "../Data/KEGG_metabolism.gmt"

# Este contiene 50 rutas metabólicas (parece que contiene los genes que
# participan en respuesta a hipoxia)
# hallmark_gmt <- '../Data/h.all.v6.1.symbols.gmt'



####################################################################################################

######################################################################
###########     1. Genes biomarcadores de linfocitos T     ###########
######################################################################


# Reunimos todos los linfocitos T y los guardamos en dos objetos sce: uno con
# todo el genoma y otro con sólo los 1566 genes metabólicos
linfoT_sce <- healthy_sce[,healthy_metabolic_sce$cellType == "T cell"]  # De 2887 células sanas, 2068 son linfos T
linfoT_metabolic_sce <- healthy_metabolic_sce[, healthy_metabolic_sce$cellType == "T cell"]

# Matrices de TPM de todo el genoma y de los 1566 genes metabólicos
linfoT_exprs <- assay(linfoT_sce, "exprs")
linfoT_metabolic_exprs <- assay(linfoT_metabolic_sce, "exprs")

length(linfoT_sce$cellType)  # 2068 Linfos T
length(linfoT_metabolic_sce$cellType) # 2068 Linfos T

dim(rowData(linfoT_sce))   # Todo el genoma de los linfos T
dim(rowData(linfoT_metabolic_sce))   # Sólo los 1566 genes metabólicos de los linfos T


# Separamos los linfocitos T CD4+ y los CD8+ ¡como cuando usamos el FACS!
# (citómetro de flujo)
pdf(file.path(outDir,"CD4_CD8_scatterplot.pdf"))
plot(t(linfoT_expr[c("CD4","CD8A"),]), main = "Scatterplot de linfocitos T")
dev.off()

# select the cells. Iniciamos un dataframe vacío de células x tipo
metadata_celulas <- data.frame(type = rep(NA, ncol(linfoT_sce)), 
                        row.names = colnames(linfoT_sce))

metadata_celulas[(linfoT_expr["CD4",] > 1) & (linfoT_expr["CD8A",] < 1),] <- "CD4"   # Catalogamos estas células como linfos T CD4+
metadata_celulas[(linfoT_expr["CD8A",] > 1) & (linfoT_expr["CD4",] < 1),] <- "CD8"   # Y estas como T CD8+


# Creamos un objeto sce para los linfos T CD4+ y los CD8+ 
linfosTCD4_CD8_sce <- linfoT_sce[, !is.na(metadata_celulas$type)]  # Cribamos los linfos NA (los que están en la región central del scatter plot)
linfosTCD4_CD8_metabolic_sce <- linfoT_metabolic_sce[, !is.na(metadata_celulas$type)]


length(linfosTCD4_CD8_sce$cellType)  # 1624 son Linfocitos T CD4+ o CD8+
length(linfosTCD4_CD8_metabolic_sce$cellType) # 1624 son Linfocitos T CD4+ o CD8+

metadata_CD4_CD8 <- metadata_celulas[!is.na(metadata_celulas$type), , drop = F]
linfosTCD4_CD8_metabolic_sce$type <- metadata_CD4_CD8$type

# Hacemos el GSEA con los genes metabólicos de los linfos T CD4+ y los CD8+
runGSEA(expr_data = linfosTCD4_CD8_metabolic_sce, 
        covariate = "type", 
        test = "CD4",
        control = "CD8",
        base.name = "t", 
        gmt.file = ruta_archivo_pathways, 
        outdir = file.path(outDir,"CD4_CD8_GSEA"))



####################################################################################################

######################################################################
###########     1. Genes biomarcadores de linfocitos T     ###########
######################################################################

# Select CD4 cells, and sort to Tregs and Th cells
cd4_tcell_sce <- linfosTCD4_CD8_sce[, linfosTCD4_CD8_sce$type == "CD4"]
cd4_tcell_expr <- assay(cd4_tcell_sce, "exprs")
tmp <- cd4_tcell_expr[c("FOXP3", "IL2RA"),]
cd4_col_annot <- data.frame(type = rep(NA, ncol(cd4_tcell_exp)), 
                            row.names = colnames(cd4_tcell_exp))


cd4_col_annot[colSums(tmp) >= 2,] <- "Tregs"
cd4_col_annot[colSums(tmp) == 0,] <- "Th"
select_cd4_tcell_sce <- cd4_tcell_sce[,!is.na(cd4_col_annot$type)]
select_cd4_tcell_metabolic_sce <- select_cd4_tcell_sce[rowData(select_cd4_tcell_sce)$metabolic,]
cd4_col_annot <- cd4_col_annot[!is.na(cd4_col_annot$type),,drop=F]
select_cd4_tcell_metabolic_sce$type <- cd4_col_annot$type
#plot selection
tmp <- cd4_tcell_exp[c("FOXP3","IL2RA"),]
pdf(file.path(outDir,"FOXP3_CD25_heatmap.pdf"),width=3,height=1)
pheatmap(tmp[,order(colSums(tmp))],show_rownames =T,show_colnames = F,cluster_rows = F,cluster_cols = F)
dev.off()
runGSEA(select_cd4_tcell_metabolic_sce,"type","Tregs","Th","t",ruta_archivo_pathways,file.path(outDir,"Tregs_Ths_GSEA"))



####################################################################################################

######################################################################
###########     2. Genes biomarcadores de fibroblastos     ###########
######################################################################

#3 Fibroblast cells: only for head and neck tumors
if (tumor == "head_neck"){
fib_sce <- selected_nontumor_sce[, selected_nontumor_sce$cellType == "Fibroblast"]
fib_exp <- assay(fib_sce, "exprs")

# Filter
select <- (fib_exp["FOS", ] >= 1) & (fib_exp["VIM",] >= 1)
select_fib_sce <- fib_sce[,select]
select_fib_exp <- assay(select_fib_sce, "exprs")

# Select CAF and Myofib, at least 2 of markers > 1
myofib_markers <- c("ACTA2", "MCAM", "MYLK", "MYL9", "PDGFA")
CAFs_markers <- c("FAP", "THY1", "PDPN", "PDGFRA", "PDGFRL", "MMP2")
select <- (apply(select_fib_exp[myofib_markers,], 2, function(x) sum(x >= 1)>=2)) | (apply(select_fib_exp[CAFs_markers,], 2, function(x) sum(x >= 1)>=2))
select_fib_sce2 <- select_fib_sce[,select]
select_fib_exp2 <- assay(select_fib_sce2, "exprs")

# Write the marker gene
dat <- select_fib_exp2[c("FOS", "VIM", myofib_markers, CAFs_markers),]
hr <- hclust(as.dist(1-cor(t(dat), method = "pearson")), method = "ward.D2")
hc <- hclust(as.dist(1-cor(dat, method = "pearson")), method = "ward.D2")

mybreaks <- c(seq(-2, 0, length.out = ceiling(200/2)+1),
              seq(2/200, 2, length.out = floor(200/2)))

library(RColorBrewer)

mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200)

pdf(file.path(outDir,"fibro_markers.pdf"), width = 3.5, height = 2, onefile = T)

pheatmap(dat, show_colnames = F, scale = "row", cluster_cols = hc, cluster_rows = hr, 
         breaks=mybreaks, legend = F, color = mycolor)

dev.off()

#cluster to CAFs or myofibroblasts
tmp = select_fib_exp2[c(CAFs_markers, myofib_markers),]
kmeans_res <- kmeans(t(tmp), centers = 2)
metadata_celulas <- data.frame(type = rep(NA, ncol(select_fib_exp2)), row.names = colnames(select_fib_exp2))

if(sum(tmp[CAFs_markers, kmeans_res$cluster == 1]) > sum(tmp[CAFs_markers, kmeans_res$cluster == 2])){
  metadata_celulas[kmeans_res$cluster==1,] <- "CAF"
  metadata_celulas[kmeans_res$cluster==2,] <- "Myofib"
}else{
  metadata_celulas[kmeans_res$cluster==1,] <- "Myofib"
  metadata_celulas[kmeans_res$cluster==2,] <- "CAF"
}
select_fib_metabolic_sce2 <- select_fib_sce2[rowData(select_fib_sce2)$metabolic,]
select_fib_metabolic_sce2$type <- metadata_celulas$type

runGSEA(select_fib_metabolic_sce2,"type","CAF","Myofib","t",ruta_archivo_pathways,file.path(outDir,"CAF_Myofib_GSEA"))
}



####################################################################################################

###################################################################
###########     3. Limpieza de archivos temporales      ###########
###################################################################

fecha <- Sys.Date()
fecha_split <- strsplit(as.character(fecha),"-")[[1]]   # Cortamos la string por los guiones
mes <- tolower(month.abb[as.numeric(date_split[2])])   # Para poner uns string en minúscula, como en Python
dia <- date_split[3]
unlink(paste0(mes, dia), recursive = T) # Eliminamos la carpeta vacía con la fecha del análisis GSEA