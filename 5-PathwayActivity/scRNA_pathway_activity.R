#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(ggplot2)
library(pheatmap)
library(reshape2)
source("../utils.R")

# Opciones
options(stringsAsFactors = F)
argumentos <- commandArgs(trailingOnly = T)
argumento <- as.character(argumentos[1])    # "head_neck" o "melanoma"
outDir <- file.path("datasets",argumento)
if (!dir.exists(outDir)) {dir.create(outDir, recursive = TRUE)}


# Leemos el dataset del head_neck/melanoma con la expresión génica imputada
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
etiquetas_celulas <- as.vector(imputed_sce$cellType)
linajes_celulares <- unique(etiquetas_celulas)


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

# Para cada linaje celular de nuestro dataset, vamos a calcular la actividad
# ponderada de cada ruta y realizaremos un test de permutación aleatoria de
# tantas permutaciones como queramos (en este caso 5000). Para ello calcularemos
# la actividad de la ruta en cada linaje celular mediante la expresión media
# relativa ponderada de los genes que participan en dicha ruta (penalizamos los
# genes que participan en muchas rutas metabólicas) y compararemos esos valores
# con los de la distribución permutada, resultado de cambiar aleatoriamente las
# etiquetas de las células. (Más información sobre los tests de permutación
# aquí: https://www.jwilber.me/permutationtest/)


## Primero inicializamos 3 matrices vacías de dimensiones rutas metabólicas X
## tipos celulares:

# Una para almacenar los p-valores de la actividad metabólica (p-valores
# calculados con el test de permutación aleatoria)
matriz_pvalues <- matrix(NA, nrow = length(nombres_pathways), 
                         ncol = length(linajes_celulares), dimnames = list(nombres_pathways, linajes_celulares))

# Otra para almacenar la actividad de las rutas metabólicas permutadas 5000
# veces
distribucion_permutada_actividad_pathways <- matrix(NA, nrow = length(nombres_pathways), 
                                                    ncol = length(linajes_celulares), dimnames = list(nombres_pathways, linajes_celulares))

# Y otra para almacenar la actividad original de las rutas metabólicas en los
# distintos tipos celulares, sobre la que compararemos la distirbución permutada
# para calcular los p-valores
distribucion_empirica_actividad_pathways <- matrix(NA, nrow = length(nombres_pathways), 
                                                   ncol = length(linajes_celulares), dimnames = list(nombres_pathways, linajes_celulares))





## Lo segundo es definir las funciones que usaremos en el bucle a continuación
group_mean <- function(x){
  sapply(linajes_celulares, function(y) rowMeans(matriz_TPM_pathway[,etiquetas_celulas_dist_nula[[x]] == y, drop = F]))
}

column_weigth_mean <- function(x){
  apply(ratio_exp_por_linaje_cel_permutaciones[[x]], 2, function(y) weighted.mean(y, pesos_genes_pathway))
}



## Ahora iteramos sobre cada ruta metabólica del archivo pathways.gmt para
## calcular su actividad relativa (a lo "Fold Change") respecto a los distintos
## linajes celulares y su p-valor por linaje celular
for (ruta_metab in nombres_pathways) {
  
  ########################################################################
  # 1ª criba: Comprobamos cuántos genes de la ruta metabólica están presentes
  # tanto en el archivo pathways.gmt como en la matriz TPM normalizada, y si
  # dicha ruta contiene menos de 5 genes, omitimos dicha ruta metabólica y
  # pasamos a la siguiente iteración del bucle
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
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, etiquetas_celulas, mean))

  # Eliminamos los genes que tengan expresión nula en cualquier linaje celular para
  # evitar ratios extremos
  genes_detectados_todos_tipos_celulares <- colnames(expr_media_por_linaje_cel)[colAlls(expr_media_por_linaje_cel > 0.001)]
  
  
  ##########################################################################
  # 2ª criba: Tras descartar los genes con expresión nula en cualquier linaje
  # celular, si la ruta metabólica conserva menos de 3 genes, la omitimos y
  # pasamos a la siguiente iteración del bucle
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
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, etiquetas_celulas, mean))
  
  # Transformamos los valores de expresión TPM de la matriz
  # `expr_media_por_linaje_cel` a ratios (similar al fold change, expresión
  # relativa al resto de tipos celulares) y transponemos la matriz para que
  # tenga dimensiones genes X tipos celulares
  ratio_exp_por_linaje_cel <- t(expr_media_por_linaje_cel) / colMeans(expr_media_por_linaje_cel)
  
  
  ##############################################################################
  # 3ª criba: Calculamos para cada linaje celular los cuartiles de los ratios de
  # expresión génica y eliminamos los genes outliers (con valores de expresión
  # anómalos) para evitar ratios extremos. Si tras este cribado la ruta
  # metabólica contiene menos de 3 genes aptos, descartamos su análisis y
  # procedemos a la siguiente iteración del bucle
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
  
  # Filtramos los genes y actualizamos la matriz de TPM, la matriz de linaje
  # celular X genes y la matriz de ratios
  genes_detectados_todos_tipos_celulares <- names(outliers)[!outliers]
  matriz_TPM_pathway <- matriz_TPM_pathway[genes_detectados_todos_tipos_celulares,]
  expr_media_por_linaje_cel <- apply(matriz_TPM_pathway, 1, function(x) by(x, etiquetas_celulas, mean))
  ratio_exp_por_linaje_cel <- t(expr_media_por_linaje_cel) / colMeans(expr_media_por_linaje_cel)
  
  # Procedemos a calcular la actividad de la ruta metabólica en cada linaje
  # celular. Para ello calcularemos la expresión media ponderada de los genes
  # que participan en dicha ruta (penalizamos los genes que participan en muchas
  # rutas metabólicas.
  pesos_genes_pathway = 1 / gene_pathway_number[genes_detectados_todos_tipos_celulares,]
  pesos_genes_pathway <- pesos_genes_pathway/sum(pesos_genes_pathway)
  actividad_ponderada_pathway <- apply(ratio_exp_por_linaje_cel, 2, 
                                 function(x) weighted.mean(x, w = pesos_genes_pathway / sum(pesos_genes_pathway)))
  
  
  
  ###############################################################
  ###############################################################
  # Comenzamos el test de permutación aleatoria para calcular la significancia
  # estadística de la actividad de las rutas metabólicas en cada estirpe
  # celular. Rellenamos la fila correspondiente de las matrices de dimensiones
  # ruta metabólica X linaje celular
  distribucion_permutada_actividad_pathways[ruta_metab, ] <-  actividad_ponderada_pathway[linajes_celulares]
  distribucion_empirica_actividad_pathways[ruta_metab, ] <-  actividad_ponderada_pathway[linajes_celulares]
  
  # Generamos la distribución permutada (No sacarlo del bucle; los resultados cambian ligeramente respecto al paper original)
  permutaciones <- 1:5000
  etiquetas_celulas_dist_nula <- lapply(permutaciones, function(x) sample(etiquetas_celulas))
  names(etiquetas_celulas_dist_nula) <- permutaciones
  
  # Cada permutación consiste en repetir todos los pasos anteriores del presente
  # bucle con los genes cribados, i.e. 1º calcular expresión media cada gen en
  # cada linaje celular, 2º calcular su actividad relativa respecto a otros
  # tipos celulares y 3º hacer la media ponderada de los genes para obtener la
  # actividad de la ruta metabólica en cada estirpe celular
  expr_media_por_linaje_cel_permutaciones <- lapply(permutaciones, function(x) group_mean(x))
  ratio_exp_por_linaje_cel_permutaciones <- lapply(permutaciones,
      function(x) expr_media_por_linaje_cel_permutaciones[[x]] / rowMeans(expr_media_por_linaje_cel_permutaciones[[x]]))
  actividad_ponderada_pathway_permutaciones <- lapply(permutaciones, function(x) column_weigth_mean(x))
  
  # Agrupamos los resultados de las permutaciones en una matriz
  matriz_permutaciones <- matrix(unlist(actividad_ponderada_pathway_permutaciones), 
                            ncol = length(linajes_celulares), byrow = T, 
                            dimnames = list(permutaciones, linajes_celulares))

  
  for (c in linajes_celulares) {
    # Si esta ruta no pasó los cribados (ergo contiene NAs), no la procesamos
    if (is.na(distribucion_permutada_actividad_pathways[ruta_metab,c])) {  
      next
    }
    
    # Calculamos los p-values
    if (distribucion_permutada_actividad_pathways[ruta_metab, c] > 1) {
      pval <- sum(matriz_permutaciones[, c] > distribucion_permutada_actividad_pathways[ruta_metab, c]) / 5000 
    } else if (distribucion_permutada_actividad_pathways[ruta_metab,c] < 1) {
      pval <- sum(matriz_permutaciones[, c] < distribucion_permutada_actividad_pathways[ruta_metab, c]) / 5000
    }
    
    # Si la ruta no es estadísticamente significativa, le asignamos como p-valor
    # un NA, ya que se muestra blanco en el heatmap
    if (pval > 0.01) {
      distribucion_permutada_actividad_pathways[ruta_metab, c] <- NA
    }
    
    # Por último rellenamos la matriz de los p-valores
    matriz_pvalues[ruta_metab, c] <- pval
  }
}


# Eliminamos en la distribución permutada las rutas que no son significativas en
# ningún linaje celular (i.e. la fila sólo contiene NAs)
rutas_no_significativas <- rowAlls(is.na(distribucion_permutada_actividad_pathways))
distribucion_permutada_actividad_pathways <- distribucion_permutada_actividad_pathways[!rutas_no_significativas,]



####################################################################################################

################################################################################
###########     2. Generación y guardado de gráficos y resultados    ###########
################################################################################

# Renombramos la variable por comodidad
data <- distribucion_permutada_actividad_pathways # No deberíamos usar la distribución empírica?

sorted_rows <- c()
sorted_columns <- c()


# Ordenamos la actividad máxima de las rutas en cada tipo celular para que
# salgan escalonadas en el heatmap
for(i in colnames(data)){
  select_row <- which(rowMaxs(data, na.rm = T) == data[,i])  # Seleciona las rutas en las que el tipo celular i presenta la actividad máxima de entre todos los tipos celulares
  tmp <- rownames(data)[select_row][order(data[select_row,i], decreasing = T)] # Ordenamos dichas rutas, de mayor a menor actividad
  sorted_rows <- c(sorted_rows, tmp)
}


# Correr o no estas dos líneas no parece tener efectos en el heatmap ...
sorted_columns <- apply(data[sorted_rows,], 2, function(x) order(x)[nrow(data)])
sorted_columns <- names(sorted_columns)


# Las casillas del heatmap sin información se marcan como no significativas
# (fold change = 1)
data[is.na(data)] <- 1 

# Generamos el PDF donde guardaremos el heatmap
pdf(file.path(outDir,"KEGGpathway_activity_heatmap_permutation.pdf"), onefile = T,
    width = 9, height = 9)

# Generamos la paleta de colores: un gradiente de azul a rojo (pasando por el
# blanco) de 100 pasos
color <- colorRampPalette(c("blue", "white", "red"))(100)

# Le decimos al heatmap que use esa paleta de manera gradual, siendo el 2 rojo,
# el 1 blanco y el 0 azul
mybreaks <- c(
  seq(0, 0.5, length.out = 33),
  seq(0.51, 1.5, length.out = 33),
  seq(1.51, max(data), length.out = 34)
)

# Graficamos el heatmap
pheatmap(data[sorted_rows, sorted_columns], cluster_cols = F,
         cluster_rows = F, color = color, breaks = mybreaks)

# anyNA(data)  # No contiene NAs
# Guardamos heatmap en pdf a disco duro
dev.off()

##### PRUEBAS ####
saveRDS(data, "data.rds")
saveRDS(sorted_rows, "sorted_rows.rds")
saveRDS(sorted_columns, "sorted_columns.rds")
saveRDS(color, "color.rds")
saveRDS(mybreaks, "mybreaks.rds")

#####################

# Guardamos las matrices de actividad de las rutas metabólicas y las de sus
# p-valores
write.table(distribucion_empirica_actividad_pathways, 
            file = file.path(outDir,"KEGGpathway_activity_empirical_dist.tsv"),
            row.names = T, col.names = NA, quote = F, sep = "\t")   # con col.names = NA ponemos bien los nomrbes de las columnas, aunque no lo parezca a priori

write.table(distribucion_permutada_actividad_pathways,
            file = file.path(outDir,"KEGGpathway_activity_permutation_dist.tsv"),
            row.names = T, col.names = NA, quote = F, sep = "\t")

write.table(matriz_pvalues, file = file.path(outDir,"KEGGpathway_activity_pvalue.tsv"),
            row.names = T, col.names = NA, quote = F, sep = "\t") 



# Creamos un violinplot para comparar la actividad metabólica global entre tipos
# celulares
scRNA_data <- as.data.frame(distribucion_empirica_actividad_pathways)
scRNA_data_flattened <- melt(scRNA_data, na.rm = T)

graf_violin <- ggplot(scRNA_data_flattened, aes(x = variable, y = value, fill = variable)) +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3, labels = 0:3) +   # Establecemos el ylim en 0-3
  geom_violin(trim = F, size = 0.2, show.legend = F, width = 1.0) +    # Violines
  labs(y = NULL, x = NULL) +                                           # Eliminamos las etiquetas de los ejes X e Y
  geom_hline(yintercept = 1, color = "black", linetype = "dashed", size = .2) +
  stat_summary(fun = median, geom = "point", size = 1, color = "blue") +  # Marcamos con un punto azul la mediana en los violines
  scale_fill_brewer(palette = "Set2") +
  theme_classic() + 
  theme(legend.position = "none",                                    # Eliminamos la leyenda porque el gráfico es auto-explicatorio
  axis.text.x = element_text(colour = "black", size = 6, angle = 45, hjust = 1, vjust = 1),  # Propiedades del texto del eje X (tipos celulares)
  axis.text.y = element_text(colour = "black", size = 6),
  axis.line = element_line(size = 0.2, color = "black"),
  axis.ticks = element_line(colour = "black", size = 0.2),
  panel.border = element_blank(), panel.background = element_blank(),
  axis.ticks.length = unit(.5, "mm"))

# Guardamos a disco duro el violinplot
ggsave(file.path(outDir,"pathway_global_activity_violinplot.pdf"), graf_violin, 
       width = 2.5, height = 1.5, units = "in", device = "pdf",
       useDingbats = FALSE)



# Mensaje de fin
print("")
print("GRACIAS POR ASISTIR A MI CHARLA TED")
