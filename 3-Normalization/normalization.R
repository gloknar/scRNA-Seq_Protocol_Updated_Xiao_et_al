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
argumento <- commandArgs()
argumento <- argumento[6]
# argumento = "head_neck"
outDir <- file.path("./datasets",argumento)
if(!dir.exists(outDir)) {                 
  dir.create(outDir, recursive = TRUE)    # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
}


# Leemos el dataset del head_neck/melanoma con la expresión génica imputada y el
# dataset con las células filtradas. Cargamos también los tipos celulares
# presentes en dichos dataset
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
# imputed_sce$cellType <- droplevels(imputed_sce$cellType)
linajes_celulares <- unique(imputed_sce$cellType)



####################################################################################################

####################################################################
###########     1. Selección de genes de referencia      ###########
####################################################################

# Para normalizar la expresión de nuestros genes, usaremos como referencia
# aquellos que tengan un ratio de detección >= 75% en todas las estirpes
# celulares (tasa de dropout < 25%)

tasa_deteccion_minima <- 0.75

# Inicializamos una matriz de genes X tipos celulares llena de FALSE, que
# usaremos de plantilla para seleccionar, para cada linaje celular, los genes
# con un dropout < 25%
matriz_seleccion_genes <- matrix(FALSE,
                                 nrow = nrow(imputed_sce),
                                 ncol = length(linajes_celulares),
                                 dimnames = list(rownames(imputed_sce), linajes_celulares))


# Iteramos sobre cada estirpe celular
for (tipo in linajes_celulares){
  sce_linaje_celular <- imputed_sce[, imputed_sce$cellType == tipo] # hacemos un subset para cada tipo celular del objeto sce con la expresión génica imputada
  matriz_genica_log_linaje_celular <- assay(sce_linaje_celular, "exprs")  # Calculamos la expresión en log2(TPM+1) de los genes de las células seleccionadas
  tasa_deteccion_gen <- apply(matriz_genica_log_linaje_celular ,1, function(x) sum(x>0)/ncol(matriz_genica_log_linaje_celular))  # Calculamos para cada gen el % de células con expresión génica no nula
  genes_seleccionados <- tasa_deteccion_gen >= tasa_deteccion_minima # Seleccionamos como TRUE aquellos genes con expresion no nula en >= 75% de las células de la estirpe en cuestión (tasa dropout < 25%)
  matriz_seleccion_genes[genes_seleccionados, tipo] <- TRUE
}


print("Nº de genes seleccionados:")
print(sum(rowSums(matriz_seleccion_genes) >= length(linajes_celulares))) # Selecciona aquellos genes con una tasa de expresión no nula >= 75% en todas las estirpes celulares (dropout en cada estirpe celular < 25%) y los guardamos en el objeto "low_dropout_genes"
low_dropout_genes <- rownames(matriz_seleccion_genes)[rowSums(matriz_seleccion_genes) >= length(linajes_celulares)]


# Limpiamos RAM
rm(sce_linaje_celular, matriz_genica_log_linaje_celular, tasa_deteccion_gen,
   genes_seleccionados, tasa_deteccion_minima, tipo)
gc(verbose = F)



####################################################################################################

############################################################################
###########     2. Preparación de datos para el normalizado      ###########
############################################################################

# Cargamos la longitud de los genes en pares de bases, tal como hicimos en el
# apartado 2 del paso de imputado
all_gene_lengths <- read.table(file.path("../Data/","gene_length.txt"), 
                               sep = "\t", header = F, row.names = 1)


# Comprobamos que dicho archivo tenga los nombres correctos de los genes y
# contenga el mismo nº de genes que en nuestro objeto sce
comprobacion <- intersect(rownames(all_gene_lengths), rownames(imputed_sce))
if (length(comprobacion) != nrow(imputed_sce)){
  warning("check the length file")
}

# De todos los genes que hay en el archivo de longitudes, nos quedamos con los
# que aparecen en nuestro objeto sce y los casteamos a floating point para
# futuras operaciones matemáticas
genelen <- all_gene_lengths[rownames(imputed_sce),]
genelen <- as.numeric(genelen)


# A partir de los TPM del objeto sce imputado, calculamos los conteos de los
# transcritos de ARNm (i.e. tamaño de librería), paso necesario para calcular el
# factor de escalado de cada célula. Tamaño librería = TPM * longitud gen (en
# pares de bases)
matriz_conteos <- sweep(imputed_sce@assays@data$tpm, 1, genelen, FUN = "*") 

# Nota: sweep es un apply que recorre de manera recursiva el 3er argumento
# (STATS), o sea hace objeto[1] * genelen[1], objeto[2] * genelen[2]...



####################################################################################################

########################################################################
###########     3. Normalización Upper-Quartile (EdgeR)      ###########
########################################################################

# A partir de la matriz de conteos, calculamos el factor de escalado para cada
# célula usando como referencia los genes con poco dropout (apartado 1)
factores_escalado_UQ <- calcNormFactors(matriz_conteos[low_dropout_genes,],
                             method = "upperquartile", p = 0.75)

# Para normalizar la expresión génica, divide los TPM de cada célula por su
# factor de escalado correspondiente obtenido con el método "upperquartile".
matriz_tpm_norm <- sweep(imputed_sce@assays@data$tpm, 2, factores_escalado_UQ, "/")
# matriz_log_tpm_imputada_norm <- log2(matriz_tpm_norm + 1)


# Guardamos la matriz de TPM normalizada para más tarde
saveRDS(matriz_tpm_norm, file.path(outDir,"UpperQuartile_tpm.rds"))



####################################################################################################

#############################################################
###########     4. Normalización TMM (EdgeR)      ###########
#############################################################

# Calculamos factores de escalado celulares.
factores_escalado_TMM <- calcNormFactors(matriz_conteos[low_dropout_genes,], method = "TMM")

# Normalizamos con el método "TMM"
matriz_tpm_norm <- sweep(imputed_sce@assays@data$tpm, 2, factores_escalado_TMM, "/")
# matriz_log_tpm_imputada_norm <- log2(matriz_tpm_norm + 1)

# Guardamos la matriz de TPM normalizada para más tarde
saveRDS(matriz_tpm_norm,file.path(outDir,"TMM_tpm.rds"))



####################################################################################################

##############################################################
###########     5. Normalización RLE (DESeq2)      ###########
##############################################################

# Calculamos factores de escalado celulares. Aquí usamos la función
# "estimateSizeFactorsForMatrix" del paquete DESeq2 versión 1.35.0
# ("https://github.com/mikelove/DESeq2/blob/master/R/core.R")
source("../utils.R")
factores_escalado_RLE <- estimateSizeFactorsForMatrix(matriz_conteos[low_dropout_genes,])

# Normalizamos con el método "RLE"
matriz_tpm_norm <- sweep(imputed_sce@assays@data$tpm, 2, factores_escalado_RLE, "/")
# matriz_log_tpm_imputada_norm <- log2(matriz_tpm_norm + 1)

# Guardamos la matriz de TPM normalizada para más tarde
saveRDS(matriz_tpm_norm, file.path(outDir,"RLE_tpm.rds"))



####################################################################################################

############################################################################
###########     6. Normalización por deconvolución  (scran)      ###########
############################################################################

# Parece que la solución era usar una funcion que sí admite matrices... o sea,
# usar "calculateSumFactors" en lugar de "computeSumFactors"
factores_escalado_deconvolucion <- scran::calculateSumFactors(matriz_conteos[low_dropout_genes,], 
                                   clusters = imputed_sce$cellType)

# Normalizamos con el método de deconvolución
matriz_tpm_norm <- sweep(tpm(imputed_sce), 2, factores_escalado_deconvolucion, "/")
# matriz_log_tpm_imputada_norm <- log2(matriz_tpm_norm + 1)

# Guardamos la matriz de TPM normalizada para más tarde
saveRDS(matriz_tpm_norm, file.path(outDir,"Deconvolution_tpm.rds"))



####################################################################################################

############################################################################################
###########     7. Benchmark de los métodos de normalización para scRNA-seq      ###########
############################################################################################

all_cell_type <- as.vector(imputed_sce$cellType)

# Limpieza RAM
rm(imputed_sce, matriz_conteos)
gc(verbose = F)


# Creamos los boxplots de los distintos métodos de normalizado probados
for (metodo in c("RLE", "TMM", "UpperQuartile", "Deconvolution")) {
  
  # Cargamos la matriz de TPM normalizada con el método en cuestión
  ruta_archivo_tpm <- file.path(outDir,paste0(metodo,"_tpm.rds"))
  matriz_tpm_norm <- readRDS(ruta_archivo_tpm) 

  # Calculamos una especie de fold change para cada gen y cada estirpe celular
  low_dropout_genes_tpm <- matriz_tpm_norm[low_dropout_genes,]  # De los genes normalizados, cogemos solo los de referencia
  low_dropout_genes_tpm_mean_by_cell_type <- apply(low_dropout_genes_tpm, 1, function(x) by(x, all_cell_type, mean)) # Calculamos la expresión media de cada gen por cada tipo celular
  low_dropout_genes_tpm_ratio <- t(low_dropout_genes_tpm_mean_by_cell_type) / colMeans(low_dropout_genes_tpm_mean_by_cell_type) # Divide la expresión media de cada gen en cada tipo celular por la media global del gen para calcular una especie de fold change. 
  datos <- reshape2::melt(low_dropout_genes_tpm_ratio) # Aplanamos el dataframe para que ggplot2 pueda usarlo correctamente
  
  # Creamos los boxplots y los guardamos en un pdf
  TPM_norm_boxplot <- ggplot(datos, aes(x = Var2, y = value)) +
       geom_boxplot(outlier.alpha = 0.1) + theme_classic() + 
       ylab("Gene expression ratio") + xlab("") +    
       theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(outDir,paste0(metodo,"_ratio_distribution.pdf")), 
         TPM_norm_boxplot, width = 3.5, height = 2.5)
}
