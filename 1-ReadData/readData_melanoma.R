######################################################################################
###########     0. Carga de paquetes, opciones y funciones de interés      ###########
######################################################################################

source("../utils.R") # Cargamos funciones definidas en el archivo `utils.R` 
library(scater)
library(stringr)
library(reshape2)
options(stringsAsFactors=FALSE)
outdir <- "dataset/melanoma"
if(!dir.exists(outdir)) {dir.create(outdir)}  # Si no existe la carpeta `head_neck`, la creamos



####################################################################################################

#################################################################
###########     1. Leemos el dataset de melanoma      ###########
#################################################################

# Los valores de expresión génica vienen en formato log2(TPM+1), o sea
# cuasi-logaritmo binario de TPM (Transcripts Per Million). TPM se interpreta
# como "por cada 1,000,000 moléculas de ARNm en la muestra secuenciada, x vienen
# de este gen." Al ser una medida estandarizada por 10^6, los transcritos de
# cada célula/muestra deberían sumar siempre 10^6. Si no se cumple eso, es que
# nuestros datos están mal (fuentes: https://www.biostars.org/p/403916/ y
# https://github.com/gpertea/stringtie/issues/213)

dataset_crudo <- "dataset/GSE72056_melanoma_single_cell_corrected.txt"
temporary_data <- read.table(dataset_crudo, head = T, sep = "\t",
                             quote = NULL, stringsAsFactors = F)

# Hay dos genes duplicados: MARCH1 y MARCH2. Eliminaremos las copias que tengan
# la expresión media más baja
march1_rows <- grep("^MARCH1$",temporary_data$Cell)
march1_rows_mean <- apply(temporary_data[march1_rows,2:ncol(temporary_data)],1,mean)
remove_rows <- c(names(march1_rows_mean)[order(march1_rows_mean)][1])

march2_rows <- grep("^MARCH2$",temporary_data$Cell)
march2_rows_mean <- apply(temporary_data[march2_rows,2:ncol(temporary_data)],1,mean)
remove_rows <- c(remove_rows,names(march2_rows_mean)[order(march2_rows_mean)][1])
remove_rows <- as.integer(remove_rows)

# col_data` es un dataframe con los metadatos de las células y es necesario para
# construir el objeto de tipo `sce` en el paso 3
col_data = data.frame(tumor = t(temporary_data[1,2:ncol(temporary_data)]),
                      malignant = t(temporary_data[2,2:ncol(temporary_data)]),
                      cellType = t(temporary_data[3,2:ncol(temporary_data)]))

# Por algún motivo no se crean bien los nombres de las columnas de `col_data`,
# por lo que los asignamos manualmente
colnames(col_data) <- c("tumor","malignant","cellType")

# Renombramos los tumores al formato TXX
col_data$tumor <- factor(paste0("T",col_data$tumor)) 

# Recodificamos la columna `malignant` y la casteamos a factor 
col_data$malignant <- car::recode(col_data$malignant, '0="Unresolved";1="Non-malignant";2="Malignant"')
col_data$malignant <- factor(col_data$malignant)

# Ídem para la columna cellType
col_data$cellType <- car::recode(col_data$cellType, 
  recodes = '0="Unknown";1="T cell";2="B cell";3="Macrophage";4="Endothelial";5="CAF";6="NK"')

# Algunas células `Unknown` son en realidad tumorales (=`malignant`)
tumor_select <- (col_data$cellType=="Unknown") & (col_data$malignant=="Malignant")
col_data[tumor_select,"cellType"] <- "Malignant"
col_data$cellType <- factor(col_data$cellType)

# Con los metadatos de las células ya listos en `col_data`, eliminamos las filas
# de las que provienen (de la 1 a la 3, más las copias de los genes MARCH1 y 2)
# para quedarnos con un dataframe de genes x células. También eliminamos la
# primera columna, que contiene los nombres de las filas
remove_rows <- c(remove_rows,1,2,3)
quasilog2_tpm <- temporary_data[-remove_rows,]
rownames(quasilog2_tpm) <- quasilog2_tpm$Cell
quasilog2_tpm <- quasilog2_tpm[,-1]
rm(temporary_data)



####################################################################################################

#############################################################################
###########     2. Marcamos los genes de las rutas de interés     ###########
#############################################################################

pathways <- gmtPathways("../Data/KEGG_metabolism.gmt")  # Obtenemos una lista de rutas metabólicas a partir de un archivo .GMT
metabolics <- unique(as.vector(unname(unlist(pathways))))  # Nos quedamos con los nombres únicos (=no repetidos) de los genes que participan en dichas rutas metabólicas, que son 1667 (los guardamos en un vector de strings)
row_data <- data.frame(metabolic = rep(FALSE, nrow(quasilog2_tpm)), row.names = rownames(quasilog2_tpm))  # Creamos un dataframe que indica para cada gen si pertenece a las rutas de interés. Este objeto es necesario para crear el objeto de tipo `sce`
row_data[rownames(row_data)%in%metabolics,"metabolic"] = TRUE  # Marcamos como TRUE los 1667 genes de interés, de un total de 23686 genes secuenciados... nos quedaremos sólo con los genes de interés en el paso de generar el objeto de tipo `Single Cell Experiment`



####################################################################################################

#####################################################################################
###########     3. Creamos el objeto `Single Cell Experiment` de scater   ###########
#####################################################################################

# Casteamos el dataframe de genes x células a una matriz (= todos los valores
# son numéricos). Me aseguré de que no hay ningún NA ni Inf, los valores de
# log2(TPM+1) están acotados entre [0,16], por lo que no se crean valores Inf al
# calcular los TPM

quasilog2_tpm <- data.matrix(quasilog2_tpm)
raw_tpm <- (2^quasilog2_tpm) - 1
# hist(apply(raw_tpm, 2, sum)) # La mayoría de células deberían sumar 1 millón de TPMs. Si no, está mal
sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm,exprs=quasilog2_tpm),
  colData = col_data,
  rowData = row_data)



####################################################################################################

##################################################
###########   4. Filtrado de células   ###########
##################################################

# Eliminaremos los tumores y las estirpes celulares que contengan < 50 células


# Creamos 2 objetos `sce`: Uno que contiene las células tumorales y otro que
# contiene las células sanas (las células de tipo Unknown no se usan)
tumor_sce <- sce[,sce$cellType == "Malignant"]
nontumor_sce <- sce[,!sce$cellType%in%c("Unknown","Malignant")]

# Nos quedamos con los tumores de más de 50 células
tumor_sample_stats <- table(tumor_sce$tumor)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats>=50]
selected_tumor_sce <- tumor_sce[,tumor_sce$tumor %in% tumor_sample_select]

# Nos quedamos con los tipos celulares sanos con más de 50 células
nontumor_sample_stats <- table(nontumor_sce$cellType)
nontumor_sample_select <- names(nontumor_sample_stats)[nontumor_sample_stats>=50]
selected_nontumor_sce <- nontumor_sce[,nontumor_sce$cellType %in% nontumor_sample_select]

# Crea el objeto `filtered_sce`, con células de todos los tipos celulares
# excepto los `Unknown` y los que tienen < 50 células. Básicamente filtramos las
# células y pasamos de tener 5902 a 5502
selected_columns <- unique(c(colnames(selected_tumor_sce),colnames(selected_nontumor_sce)))
filtered_sce <- sce[,colnames(sce) %in% selected_columns]

# Renombramos los tumores de TXX a MELXX
unique_tumors <- unique(filtered_sce$tumor)
levels(filtered_sce$tumor) = paste0("MEL",seq(1,length(unique_tumors)))



####################################################################################################

######################################################################
###########   5. Guardamos los objetos `sce` resultantes   ###########
######################################################################

saveRDS(filtered_sce,file.path(outdir,"filtered_sce.rds"))
