#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
source("../utils.R")

library(stringr)
library(reshape2)
library(scales)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)


# Opciones
options(stringsAsFactors = F)
# argumento <- commandArgs()
# argumento <- argumento[6]
argumento <- "melanoma"

outDir <- file.path("datasets",argumento)
if(!dir.exists(outDir) ) {dir.create(outDir, recursive = TRUE)}


# Leemos el dataset del head_neck/melanoma con la expresión génica imputada
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
all_cell_types <- as.vector(imputed_sce$cellType)
cell_types <- unique(all_cell_types)


# Leemos la matriz de TPM normalizada por deconvolución del dataset del
# head_neck/melanoma
set.seed(123)
metodo_normalizado <- "Deconvolution"
ruta_matriz_TPM_norm <- file.path("../3-Normalization/datasets",argumento,paste0(metodo_normalizado,"_tpm.rds"))
matriz_TPM_norm <- readRDS(ruta_matriz_TPM_norm)


# Para leer el archivo de las rutas metabólicas y los genes que participan en
# ellas (obtenido de KEGG), usamos la función auxiliar `gmtPathways()`
ruta_archivo_pathways <- "../Data/KEGG_metabolism.gmt"
pathways.gmt <- gmtPathways(ruta_archivo_pathways)
nombres_pathways <- names(pathways.gmt)


# Vamos a calcular con la función auxiliar `num_of_pathways()` el nº de rutas
# metabólicas en las que participan nuestros genes de interés (1566 genes
# metabólicos)
gene_pathway_number <- num_of_pathways(ruta_archivo_pathways, 
                                       rownames(imputed_sce)[rowData(imputed_sce)$metabolic])



####################################################################################################

###################################################################################
###########     1. Cálculo de la actividad de las rutas metabólicas     ###########
###################################################################################

# Inicializamos 3 matrices vacías de rutas metabólicas X tipos celulares:

# Una para almacenar los p-valores de la actividad metabólica (p-valores
# calculados con el método de shuffle)
matriz_pvalues <- matrix(NA, nrow = length(nombres_pathways), 
                         ncol = length(cell_types), dimnames = list(nombres_pathways, cell_types))

# Otra para #mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA, nrow = length(nombres_pathways), 
                           ncol = length(cell_types), dimnames = list(nombres_pathways, cell_types))

# Y otra para
mean_expression_noshuffle <- matrix(NA, nrow = length(nombres_pathways), 
                             ncol = length(cell_types), dimnames = list(nombres_pathways, cell_types))

# View(mean_expression_shuffle)
# View(mean_expression_noshuffle)
# View(matriz_TPM_norm[genes_compartidos, ])

for (ruta_metab in nombres_pathways) {
  
  # Comprobamos cuántos genes de la ruta metabólica están presentes tanto en el
  # archivo pathways.gmt como en la matriz TPM normalizada, y si dicha ruta
  # contiene menos de 5 genes, omitimos dicha ruta metabólica y pasamos a la
  # siguiente iteración del bucle
  genes_ruta <- pathways.gmt[[ruta_metab]]
  genes_compartidos <- intersect(genes_ruta, rownames(matriz_TPM_norm))
  if (length(genes_compartidos) < 5) {
    message('La ruta metabólica "',ruta_metab,'" presenta menos de 5 genes en la matriz de TPMs normalizada. Omitiendo dicha ruta...')
    next
  }

  # Generamos una matriz de TPM de los genes de la ruta en cuestión y
  # descartamos aquellos con expresión nula en 100% de las células del dataset,
  # i.e. tasa de dropout global de 100%
  matriz_TPM_pathway <- matriz_TPM_norm[genes_compartidos, ]
  matriz_TPM_pathway <- matriz_TPM_pathway[rowSums(matriz_TPM_pathway) > 0,] 
  
  # A partir de la matriz anterior creamos otra de dimensiones linaje celular X
  # gen de la ruta metabólica que contenga la expresión media de dicho gen en
  # cada tipo celular
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, all_cell_types, mean))

  # Eliminamos los genes que tengan expresión nula en cualquier linaje celular para
  # evitar ratios extremos
  genes_detectados_todos_tipos_celulares <- colnames(expr_media_por_linaje_cel)[colAlls(expr_media_por_linaje_cel > 0.001)]

  # Tras descartar los genes con expresión nula en cualquier linaje celular, si
  # la ruta metabólica conserva menos de 3 genes, la omitimos y pasamos a la
  # siguiente iteración del bucle
  if (length(genes_detectados_todos_tipos_celulares) < 3) {
    message('La ruta metabólica "',ruta_metab,'" contiene menos de 3 genes con expresión detectada en todos los tipos celulares. Omitiendo dicha ruta...')
    next
  }
  
  # Actualizamos la matriz de TPM con los genes que han pasado el cribado y,
  # para cada gen de la susodicha, imputamos las células con expresión nula con
  # el valor mínimo de dicho gen para evitar ratios extremos
  matriz_TPM_pathway <- matriz_TPM_pathway[genes_detectados_todos_tipos_celulares,]
  matriz_TPM_pathway <- t(apply(matriz_TPM_pathway, 1, function(x) {x[x <= 0] <- min(x[x > 0]); return(x)}))

  # Actualizamos también la matriz de linaje celular X genes
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, all_cell_types, mean))
  
  # Transformamos los valores de expresión TPM de la matriz
  # `expr_media_por_linaje_cel` a ratios (similar al fold change, expresión
  # relativa al resto de tipos celulares) y transponemos la matriz para que
  # tenga dimensiones genes X tipos celulares
  ratio_exp_por_linaje_cel <- t(expr_media_por_linaje_cel) / colMeans(expr_media_por_linaje_cel)
  
  
  
  # Calculamos para cada linaje celular los cuartiles de los ratios de expresión
  # génica y eliminamos los genes outliers (con valores de expresión anómalos)
  # para evitar ratios extremos. Si tras este cribado la ruta metabólica
  # contiene menos de 3 genes aptos, descartamos su análisis y procedemos a la
  # siguiente iteración del bucle
  cuartiles_ratios_exp_genica <- apply(ratio_exp_por_linaje_cel, 2, function(x) quantile(x, na.rm = T))
  cuartil1 <- cuartiles_ratios_exp_genica["25%",]
  limite_inferior <- cuartil1 / 3
  cuartil3 <- cuartiles_ratios_exp_genica["75%",]
  limite_superior <- cuartil3 * 3
  outliers <- apply(ratio_exp_por_linaje_cel, 1, function(x) {any((x < limite_inferior) | (x > limite_superior))})
  
  if(sum(!outliers) < 3) {
    message('Tras eliminar los genes con valores de expresión anormales (outliers), la ruta metabólica "',ruta_metab,'" conserva menos de 3 genes. Omitiendo dicha ruta...')
    next
  }
  
  # Cribamos los genes por 3a vez y actualizamos la matriz de TPM, la matriz de
  # linaje celular X genes y la matriz de ratios
  genes_detectados_todos_tipos_celulares <- names(outliers)[!outliers]
  matriz_TPM_pathway <- matriz_TPM_pathway[genes_detectados_todos_tipos_celulares,]
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, all_cell_types, mean))
  ratio_exp_por_linaje_cel <- t(expr_media_por_linaje_cel) / colMeans(expr_media_por_linaje_cel)
  
  
  # Procedemos a calcular la actividad de la ruta metabólica en cada linaje
  # celular. Para ello calcularemos la expresión media ponderada de los genes
  # que participan en dicha ruta (penalizamos los genes que participan en muchas
  # rutas metabólicas.
  pesos_genes_pathway = 1 / gene_pathway_number[genes_detectados_todos_tipos_celulares,]
  actividad_ponderada_pathway <- apply(ratio_exp_por_linaje_cel, 2, 
                                 function(x) weighted.mean(x, w = pesos_genes_pathway / sum(pesos_genes_pathway)))
  mean_expression_shuffle[ruta_metab, ] <-  actividad_ponderada_pathway[cell_types]
  mean_expression_noshuffle[ruta_metab, ] <-  actividad_ponderada_pathway[cell_types]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(matriz_TPM_pathway[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_por_linaje_cel_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pesos_genes_pathway/sum(pesos_genes_pathway)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_por_linaje_cel_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  actividad_ponderada_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(actividad_ponderada_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[ruta_metab,c])) next
    if(mean_expression_shuffle[ruta_metab,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[ruta_metab,c]) / 5000 
    }else if(mean_expression_shuffle[ruta_metab,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[ruta_metab,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[ruta_metab, c] <- NA  ### NA is  blank in heatmap
    matriz_pvalues[ruta_metab,c] <- pval
  }
}


all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]
#heatmap
dat <- mean_expression_shuffle

sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1
pdf(file.path(outDir,"KEGGpathway_activity_heatmap.pdf"),onefile=T,width=6,height=9)
mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("blue","white","red"))(100)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,breaks=mybreaks)
dev.off()

write.table(mean_expression_noshuffle,file=file.path(outDir,"KEGGpathway_activity_noshuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file=file.path(outDir,"KEGGpathway_activity_shuffle.txt"),row.names=T,col.names=T,quote=F,sep="\t")
write.table(matriz_pvalues,file=file.path(outDir,"KEGGpathway_activity_shuffle_pvalue.txt"),row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
scRNA_dat <- as.data.frame(mean_expression_noshuffle)
scRNA_dat$X <- NULL

scRNA_df <- melt(scRNA_dat)
scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
  scale_y_continuous(limits=c(0,3),breaks=0:3,labels=0:3)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue")+
  scale_fill_brewer(palette="Set2") +
  theme_classic() + 
  theme(legend.position="none",
  axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
  axis.text.y=element_text(colour="black", size = 6),
  axis.line=element_line(size=0.2,color="black"),
  axis.ticks = element_line(colour = "black",size=0.2),
  panel.border = element_blank(), panel.background = element_blank(),
  axis.ticks.length= unit(.5, "mm"))

ggsave(file.path(outDir,"pathway_activity_violinplot.pdf"),g,width = 2.5,height=1.5,units="in",device="pdf",useDingbats=FALSE)

