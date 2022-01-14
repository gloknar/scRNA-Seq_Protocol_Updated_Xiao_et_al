#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes
library(scater)
library(pheatmap)
library(RColorBrewer)

# Opciones
options(stringsAsFactors = FALSE)
argumento <- commandArgs()
argumento <- argumento[6]
# argumento <- "head_neck"

outDir <- file.path("./datasets",argumento)
if(!dir.exists(outDir)) {                    # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
  dir.create(outDir, recursive = TRUE)
}



# Leemos el dataset del head_neck/melanoma con la expresión génica imputada y de
# ahí generamos un objeto sce que contenga sólo las células tumorales y los
# genes de interés (en este caso, los metabólicos)
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
tumor_metabolic_sce <- imputed_sce[rowData(imputed_sce)$metabolic, imputed_sce$cellType == "Malignant"]

# Limpieza RAM
rm(imputed_sce)
invisible(gc(verbose = F))



####################################################################################################

########################################################################
###########     1. Cálculo de la matriz de correlación       ###########
########################################################################

# Para ver cuán similares o disimilares son los tumores de los distintos
# pacientes, vamos a generar y visualizar en un heatmap la matriz de correlación
# de las células de dichas neoplasias.
expr_genes_metab_tumor <- assay(tumor_metabolic_sce, "exprs")  # Cogemos la expresión génica en formato log2(TPM+1)
matriz_cor  <- cor(expr_genes_metab_tumor, method = "spearman") # Calculamos la matriz de correlación de spearman (no paramétrico)

# Nota: si quieres visualizar el dendrograma, quítale los nombres a las filas y
# las columnas del input de hclust:
# rownames(matriz_cor) <- NULL
# colnames(matriz_cor) <- NULL



# Generamos los clusters con el método de ward (el mejor, pero más lento) para
# el heatmap
clusters_tumores <- hclust(as.dist(1- matriz_cor), method = "ward.D2") # as.dist(1-matriz_cor) porque hclust trabaja con distancias, y la correlación es el opuesto de la distancia


# Creamos las paletas de colores del heatmap con el comando "colorRampPalette".
# Nota: Para visualizar los colores antes de usarlos:
#
# library(scales)
# show_col(my_colour)
colores_heatmap = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8) # creamos una paleta de color de 7 pasos, desde azul a rojo, pasando por el blanco. Esta es la paleta estándar de `pheatmap`
metadatos_heatmap <- as.data.frame(colData(tumor_metabolic_sce)[,"tumor", drop = F]) # Creamos un dataframe con cada célula y su tumor de procedencia
metadatos_heatmap$tumor <- factor(metadatos_heatmap$tumor) # Al castear de nuevo este factor, eliminamos niveles no usados. Sinónimo de usar `droplevels()`

# Creamos un pdf de 5x4 pulgadas donde guardaremos el heatmap
pdf(file.path(outDir,"malignant_metabolic_correlationMatrix.pdf"), 
    width = 7, height = 5, onefile = T) 

# Generamos la visualización de la matriz de correlación
pheatmap(matriz_cor,
         main = paste0("Matriz de correlación de spearman del dataset ",argumento),
         color=colores_heatmap,
         cluster_rows = clusters_tumores,     # Dendrograma creado con el método de ward.D2
         cluster_cols = clusters_tumores,
         annotation_row = metadatos_heatmap,  # Dataframe con el tumor al que pertenece cada célula maligna
         annotation_legend = T,               # Para mostrar los tumores
         show_rownames = F,                   # No mostramos los nombres de las 1167 células tumorales para poder visualizar el heatmap
         show_colnames = F)

# Guardamos el heatmap en el disco duro
dev.off() 



# En el dataset del melanoma, los pacientes MEL12 y MEL3 tienen un background
# genético parecido. Vamos a agruparlos en el heatmap
if (argumento == "melanoma"){
   order_names <- clusters_tumores$labels[clusters_tumores$order]
   MEL12_names <- rownames(metadatos_heatmap)[metadatos_heatmap$tumor == "MEL12"]
   order_names2 <- MEL12_names
   MEL3_names <- rownames(metadatos_heatmap)[metadatos_heatmap$tumor == "MEL3"]
   order_names2 <- c(order_names2,MEL3_names)
   other_names <- order_names
   other_names <- order_names[!(order_names %in% MEL3_names)]  # Quitamos las celulas que estan en MEL3 y MEL12
   other_names <- other_names[!(other_names %in% MEL12_names)]
   order_names2 <- c(order_names2,other_names)
   matriz_cor2 <- matriz_cor[order_names2,order_names2]
   
   # Creamos un pdf de 5x4 pulgadas donde guardaremos el heatmap
   pdf(file.path(outDir,"malignant_metabolic_correlationMatrix.pdf"),
       width = 7, height = 5, onefile = T)
   
   # Generamos la visualización de la matriz de correlación
   pheatmap(matriz_cor2,
            main = paste0("Matriz de correlación de spearman del dataset ",argumento),
            color = colores_heatmap,
            cluster_rows = F,
            cluster_cols = F,
            annotation_row = metadatos_heatmap,
            annotation_legend = T,
            show_rownames = F,
            show_colnames = F)
   
   # Guardamos el heatmap en el disco duro
   dev.off()
}
