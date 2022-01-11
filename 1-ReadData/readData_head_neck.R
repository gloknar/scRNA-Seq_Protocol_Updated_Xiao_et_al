######################################################################################
###########     0. Carga de paquetes, opciones y funciones de interés      ###########
######################################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(reshape2)
library(plyr)
source("../utils.R") # Cargamos funciones definidas en el archivo `utils.R` 

# Opciones
options(stringsAsFactors = FALSE)
outDir <- "datasets/head_neck"
if(!dir.exists(outDir)) {dir.create(outDir)}  # Si no existe la carpeta `head_neck`, la creamos



####################################################################################################

#############################################################################
###########     1. Leemos el dataset de head and neck cancer      ###########
#############################################################################

# Los valores de expresión génica vienen en formato log2(TPM+1), o sea
# cuasi-logaritmo binario de TPM (Transcripts Per Million). TPM se interpreta
# como "por cada 1,000,000 moléculas de ARNm en la muestra secuenciada, x vienen
# de este gen." Al ser una medida estandarizada por 10^6, los transcritos de
# cada célula/muestra deberían sumar siempre 10^6. Si no se cumple eso, es que
# nuestros datos están mal (fuentes: https://www.biostars.org/p/403916/ y
# https://github.com/gpertea/stringtie/issues/213)

dataset_crudo <- "datasets/GSE103322_HNSCC_all_data.txt"
temporary_data <- read.table(dataset_crudo, head = T, sep = "\t", row.names = 1,
                             quote = "\'", stringsAsFactors = F)

# Obtenemos nombres de los tumores
tumor_samples_names <- sapply(str_split(colnames(temporary_data),"_"), function(x) x[1], simplify = T) # Devuelve un vector con los nombres de las celulas (corta los nombres por las barrabajas y se queda con el primero de los cachitos, que tiene formas como "HN28" o "HNSCC26")
tumor_samples_names <- str_sub(tumor_samples_names, -2, -1)  # Elimina lo de HN y HNSCC... se queda sólo con los 2 últimos caracteres, que son 2 números
tumor_samples_names <- paste0("MEEI",str_replace(tumor_samples_names,"C","")) # Les añade el prefijo MEEI, de manera que se quedan MEEII26, MEEII28, etc... y si encuentra una C, la borra

# Obtenemos los tipos celulares
cell_type <- as.character(temporary_data[5,]) # Crea un vector de strings a partir de la fila 5, que contiene los tipos de célula (tiene mezclados números y strings, por eso lo castea todo a string)
malignant_cells <- as.character(temporary_data[3,]) == "1"   # La fila 3 codifica si la célula es normal ("0") o tumoral ("1")
cell_type[malignant_cells] <- "Malignant"  # Coge el vector de tipos celulares, y las que tengan en la fila 3 el valor 1, las recataloga como "Malignant"
cell_type[cell_type == 0] <- "Unknown" # Las células catalogadas como "0" se renombran a "Unknown" (tipo celular desconocido)

# `metadatos` es un dataframe con los metadatos de las células y es necesario
# para construir el objeto de tipo `sce` en el paso 3
metadatos <- data.frame(tumor = tumor_samples_names,
                       cellType = cell_type,
                       lymph = as.integer(temporary_data[2,]),
                       row.names = colnames(temporary_data))

# Con los metadatos de las células ya listos en `metadatos`, eliminamos las filas
# de las que provienen (de la 1 a la 5) para quedarnos con un dataframe de genes
# x células.
remove_rows <- c(1,2,3,4,5)
quasilog2_tpm <- temporary_data[-remove_rows,]
rm(temporary_data) # Eliminamos el dataframe temporal, pues no lo usaremos más



####################################################################################################

#############################################################################
###########     2. Marcamos los genes de las rutas de interés     ###########
#############################################################################

pathways <- gmtPathways("../Data/KEGG_metabolism.gmt") # Obtenemos una lista de rutas metabólicas a partir de un archivo .GMT
metabolics <- unique(as.vector(unname(unlist(pathways)))) # Nos quedamos con los nombres únicos (=no repetidos) de los genes que participan en dichas rutas metabólicas, que son 1667 (los guardamos en un vector de strings)
row_data <- data.frame(metabolic = rep(FALSE, nrow(quasilog2_tpm)), row.names = rownames(quasilog2_tpm)) # Creamos un dataframe que indica para cada gen si pertenece a las rutas de interés. Este objeto es necesario para crear el objeto de tipo `sce`
row_data[rownames(row_data) %in% metabolics, "metabolic"] = TRUE  # Marcamos como TRUE los 1667 genes de interés, de un total de 23686 genes secuenciados... nos quedaremos sólo con los genes de interés en el paso de generar el objeto de tipo `Single Cell Experiment`



####################################################################################################

#####################################################################################
###########     3. Creamos el objeto `Single Cell Experiment` de scater   ###########
#####################################################################################

# Casteamos el dataframe de genes x células a una matriz (= todos los valores
# son numéricos). Me aseguré de que no hay ningún NA ni Inf, los valores de
# log2(TPM+1) están acotados entre [0,6068], lo cual es imposible (TPM mide por
# millón de lecturas, no puede tener un valor superior a 10^6 = 2^19.93157). Voy
# a transformar todos los valores superiores a 16 en 16 para acotar los
# log2(TPM+1) entre [0,16]

quasilog2_tpm <- data.matrix(quasilog2_tpm)
quasilog2_tpm[quasilog2_tpm > 16] <- 16 # Con esta acotación deberíamos evitar crear valores infinitos
raw_tpm <- (2^quasilog2_tpm) - 1 # Ahora no deberían introducirse infinitos y el paso de la imputación debería funcionar
# hist(apply(raw_tpm, 2, sum)) # La mayoría de células deberían sumar 1 millón de TPMs. Si no, está mal
sce <- SingleCellExperiment(assays = list(tpm = raw_tpm, exprs = quasilog2_tpm),
                            colData = metadatos,
                            rowData = row_data)

16

####################################################################################################

##################################################
###########   4. Filtrado de células   ###########
##################################################

# Eliminaremos los tumores y las estirpes celulares que contengan < 50 células


# Creamos 2 objetos `sce`: Uno que contiene las células tumorales y otro que
# contiene las células sanas (las células de tipo Unknown no se usan)
tumor_sce <- sce[,sce$cellType == "Malignant"]
nontumor_sce <- sce[, !sce$cellType %in% c("Unknown", "Malignant")]

# Nos quedamos con los tumores de más de 50 células
tumor_sample_stats <- table(tumor_sce$tumor)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats >= 50] # Eliminaron los tipos celulares con <50 células
selected_tumor_sce <- tumor_sce[, tumor_sce$tumor %in% tumor_sample_select]

# Nos quedamos con los tipos celulares sanos con más de 50 células
nontumor_stats <- table(nontumor_sce$cellType)
nontumor_select <- names(nontumor_stats)[nontumor_stats >= 50]    # Eliminaron los tipos celulares con <50 células
selected_nontumor_sce <- nontumor_sce[, nontumor_sce$cellType %in% nontumor_select]

# Crea el objeto `filtered_sce`, con células de todos los tipos celulares
# excepto los `Unknown` y los que tienen < 50 células. Básicamente filtramos las
# células y pasamos de tener 5902 a 5502
selected_columns <- unique(c(colnames(selected_tumor_sce), colnames(selected_nontumor_sce)))
filtered_sce <- sce[, colnames(sce) %in% selected_columns]

# Renombramos los tumores de MEEIX a HNSX
filtered_sce$tumor <- factor(filtered_sce$tumor)
filtered_sce$tumor <- mapvalues(filtered_sce$tumor, 
                                from = c("MEEI5","MEEI6","MEEI7","MEEI8","MEEI10","MEEI12","MEEI13","MEEI16","MEEI17","MEEI18","MEEI20","MEEI22","MEEI23","MEEI24","MEEI25","MEEI26","MEEI28","MEEIC"),
                                to = paste0("HNS",seq(18))) # mapvalues sustituye los valores que tu le digas de un vector/factor por otros

# Nos aseguramos de que los factores contengan sólo los niveles presentes en el dataset
filtered_sce$tumor <- droplevels(filtered_sce$tumor)
filtered_sce$cellType <- factor(filtered_sce$cellType) 
filtered_sce$cellType <- droplevels(filtered_sce$cellType) 



####################################################################################################

####################################################################
###########   5. Guardamos los objetos sce resultantes   ###########
####################################################################

saveRDS(filtered_sce,file.path(outDir,"filtered_sce.rds"))
