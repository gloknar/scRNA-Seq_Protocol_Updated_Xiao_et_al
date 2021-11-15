# Carga de paquetes, opciones y funciones de interés
source("../utils.R") # Cargamos funciones definidas en el archivo `utils.R` 
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
library(plyr)

outdir <- "dataset/head_neck"
if(!dir.exists(outdir)) dir.create(outdir) # Si no existe la carpeta `head_neck`, la creamos



####################################################################################################

#############################################################################
###########     1. Leemos el dataset de head and neck cancer      ###########
#############################################################################

# Los valores de expresión génica vienen en formato log2(TPM+1), o sea
# cuasi-logaritmo binario de TPM (Transcripts Per Million). TPM se interpreta
# como "for every 1,000,000 RNA molecules in the RNA-seq sample, x came from
# this gene/transcript."
raw_tpm_file <- "dataset/GSE103322_HNSCC_all_data.txt"
tmp_data <- read.table(raw_tpm_file, head = T, sep = "\t", row.names = 1,
                       quote = "\'", stringsAsFactors = F)

# Obtenemos nombres de los tumores
tumor <- sapply(str_split(colnames(tmp_data),"_"),function(x) x[1], simplify = T) # Devuelve un vector con los nombres de las celulas (corta los nombres por las barrabajas y se queda con el primero de los cachitos, que tiene formas como "HN28" o "HNSCC26")
tumor <- str_sub(tumor,-2,-1)  # Elimina lo de HN y HNSCC... se queda sólo con los 2 últimos caracteres, que son 2 números
tumor <- paste0("MEEI",str_replace(tumor,"C","")) # Les añade el prefijo MEEI, de manera que se quedan MEEII26, MEEII28, etc... y si encuentra una C, la borra

# Obtenemos los tipos celulares
cell_type <- as.character(tmp_data[5,]) # Crea un vector de strings a partir de la fila 5, que contiene los tipos de célula (tiene mezclados números y strings, por eso lo castea todo a string)
malignant <- as.character(tmp_data[3,]) == "1"   # La fila 3 codifica si la célula es normal ("0") o tumoral ("1")
cell_type[malignant] <- "Malignant"  # Reformatea los valores de la fila 3 de "1" (tumoral) a "Malignant"
cell_type[cell_type==0] <- "Unknown" # Las células catalogadas como "0" se renombran a "Unknown" (tipo celular desconocido)


# column data (aquí creo que genera un dataframe con los metadatos de las células)
col_data <- data.frame(tumor=tumor,cellType=cell_type,
                       lymph=as.integer(tmp_data[2,]),
                       row.names=colnames(tmp_data))

# Eliminamos las filas que contienen los metadatos de las células (de la 1 a la
# 5) para quedarnos con un dataframe de genes x células. Los metadatos de las
# células están en `col_data`
remove_rows <- c(1,2,3,4,5)
all_data <- tmp_data[-remove_rows,]
rm(tmp_data) # Eliminamos el dataframe temporal, pues no lo usaremos más



####################################################################################################

#################################################################
###########     2. Marcamos los genes metabólicos     ###########
#################################################################

pathways <- gmtPathways("../Data/KEGG_metabolism.gmt") # Obtenemos una lista de rutas metabólicas a partir de un archivo .GMT
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE




####################################################################################################

#####################################################################################
###########     3. Creamos el objeto `Single Cell Experiment` de scater   ###########
#####################################################################################

all_data <- data.matrix(all_data)
raw_tpm <- 2^ all_data - 1
sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm, exprs=all_data),
  colData = col_data,
  rowData = row_data
)


####################################################################################################
#4. cell types
#malignant cells
tumor_sce <- sce[,sce$cellType == "Malignant"]
nontumor_sce <- sce[,!sce$cellType%in%c("Unknown","Malignant")]
#select tumor cells
tumor_sample_stats <- table(tumor_sce$tumor)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats>=50]
selected_tumor_sce <- tumor_sce[,tumor_sce$tumor %in% tumor_sample_select]
#select notumor
nontumor_stats <- table(nontumor_sce$cellType)
nontumor_select <- names(nontumor_stats)[nontumor_stats>=50]
selected_nontumor_sce <- nontumor_sce[,nontumor_sce$cellType %in% nontumor_select]
#select sce
selected_columns <- unique(c(colnames(selected_tumor_sce),colnames(selected_nontumor_sce)))
selected_sce <- sce[,colnames(sce) %in% selected_columns]

#rename the patients
selected_sce$tumor <- factor(selected_sce$tumor)
selected_sce$tumor <- mapvalues(selected_sce$tumor, from=c("MEEI5","MEEI6","MEEI7","MEEI8","MEEI10","MEEI12","MEEI13","MEEI16","MEEI17","MEEI18","MEEI20","MEEI22","MEEI23","MEEI24","MEEI25","MEEI26","MEEI28","MEEIC"),
                              to=paste0("HNS",seq(18))) 

#saveRDS(sce,file.path(outdir,"sce.rds"))
saveRDS(selected_sce,file.path(outdir,"selected_sce.rds"))
