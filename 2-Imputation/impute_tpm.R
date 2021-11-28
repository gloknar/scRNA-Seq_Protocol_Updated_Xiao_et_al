#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes
library(scImpute)
library(scater)

# Opciones
argumento <- commandArgs()[6]
outDir <- file.path("./datasets", argumento) # Carpeta donde guardaremos todos los archivos relacionados con la imputación del objeto `sce`
if(!dir.exists(outDir)) {                   # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
  dir.create(outDir,recursive = TRUE)
} 
num_cores <- 6                              # Usar 1 en Windows (scImpute usa mc.apply...)

# Leemos el dataset del head_neck/melanoma con las células filtradas
filtered_sce <- readRDS(file.path("../1-ReadData/datasets",argumento,"filtered_sce.rds"))



####################################################################################################

###########################################################################
###########     1. Preparado de datasets y matrices de TPM      ###########
###########################################################################

# Creamos 2 objetos `sce`: Uno que contiene las células tumorales y otro que
# contiene las células sanas (las imputaremos por separado)
filtered_sce_tumor <- filtered_sce[,filtered_sce$cellType=="Malignant"]
filtered_sce_nontumor <- filtered_sce[,filtered_sce$cellType!="Malignant"]

# Creamos las matrices de TPMs de ambos subconjuntos celulares
# NOTA: el bolsillo `filtered_sce_(non)tumor@assays@data$exprs` contiene la
# expresión génica de las células en formato log2(TPM+1);;
# El bolsillo `filtered_sce_(non)tumor@assays@data$tpm` contiene la expresión
# génica en formato TPM
filtered_sce_tumor_tpm <- tpm(filtered_sce_tumor)
filtered_sce_nontumor_tpm <- tpm(filtered_sce_nontumor) 
labels_tumor <- filtered_sce_tumor$tumor
labels_nontumor <- filtered_sce_nontumor$cellType

# Guardamos las matrices de TPMs en archivos
write.csv(filtered_sce_tumor_tpm,file.path(outDir,"tumor.tpm"))
write.csv(filtered_sce_nontumor_tpm,file.path(outDir,"nontumor.tpm"))



####################################################################################################

#############################################################################################
###########     2. Obtención de las longitudes de los genes en pares de bases     ###########
#############################################################################################

# Cargamos el archivo que contiene la longitud de todos los genes humanos en
# pares de bases (esto servía para pasar de TPM a read counts mediante TPM*gene
# length)
all_gene_lengths <- read.table("../Data/gene_length.txt", sep = "\t",
                               header = F, row.names = 1) # Longitudes en bp

# Comprobamos si tenemos el mismo nº de genes en nuestro dataset y en el archivo
# con las longitudes de los genes
temporary <- intersect(rownames(all_gene_lengths),rownames(filtered_sce_tumor_tpm)) # La intersección busca nombres de genes compartidos en ambos conjuntos de nombres de genes
if (length(temporary) != nrow(filtered_sce_tumor_tpm)){ # Compara si el nº de genes en filtered_sce_tumor_tpm y all_gene_lengths es igual. Si todo está bien, no pasa nada, pero si no, emite una advertencia y termina el programa
  warning("check the length file")
  print(setdiff(rownames(filtered_sce_tumor_tpm),rownames(all_gene_lengths)))
  q()
}

# Preparamos el vector de longitudes de genes necesario para scImpute
genelen <- all_gene_lengths[rownames(filtered_sce_tumor_tpm),]
genelen <- as.numeric(as.vector(genelen))



####################################################################################################

####################################################################################
###########     3.1 Imputación de la expresión génica (TPM) tumoral      ###########
####################################################################################

# Imputamos genes con dropout >= 0.5 para evitar sobre-imputación
scimpute(count_path = file.path(outDir, "tumor.tpm"), infile = "csv", 
         outfile = "csv", out_dir = paste(outDir,"malignant/", sep = "/"), 
         labeled = TRUE, labels = as.vector(labels_tumor), type = "TPM", 
         genelen = genelen, drop_thre = 0.5, ncores = num_cores)


# Con la imputación ya hecha, la cargamos en memoria y actualizamos los
# bolsillos `filtered_sce_tumor@assays@data$exprs` y
# `filtered_sce_tumor@assays@data$tpm`
imputed_tpm_tumor <- read.csv(file.path(outDir,"malignant/scimpute_count.csv"),
                              header = T, row.names = 1)
tpm(filtered_sce_tumor) <- data.matrix(imputed_tpm_tumor)
assay(filtered_sce_tumor,"exprs") <- data.matrix(log2(imputed_tpm_tumor + 1)) # Recuerda que `exprs` contiene la expresion génica en formato log2(TPM+1)



####################################################################################################

###################################################################################
###########     3.2 Imputación de la expresión génica (TPM) normal      ###########
###################################################################################

# Imputamos genes con dropout >= 0.5 para evitar sobre-imputación
scimpute(count_path = file.path(outDir, "nontumor.tpm"), infile = "csv", 
         outfile = "csv", out_dir = paste(outDir,"non-malignant/", sep = "/"), 
         labeled = TRUE, labels = as.vector(labels_nontumor),	type = "TPM", 
         genelen = genelen, drop_thre = 0.5, ncores = num_cores)


# Con la imputación ya hecha, la cargamos en memoria y actualizamos los
# bolsillos `filtered_sce_nontumor@assays@data$exprs` y
# `filtered_sce_nontumor@assays@data$tpm`
imputed_tpm_nontumor <- read.csv(file.path(outDir,"non-malignant/scimpute_count.csv"),
                                 header = T, row.names = 1)
tpm(filtered_sce_nontumor) <- data.matrix(imputed_tpm_nontumor) 
assay(filtered_sce_nontumor,"exprs") <- data.matrix(log2(imputed_tpm_nontumor + 1))



####################################################################################################

###################################################################################
###########     4 Construcción del nuevo objeto `sce` imputado      ###############
###################################################################################

# Concatenamos los TPMs de las células sanas y las tumorales y los ordenamos
# como estaban en el objeto original `filtered sce`. Repetimos el proceso para
# los log2(TPM+1)
imputed_tpm_total <- cbind(tpm(filtered_sce_tumor), tpm(filtered_sce_nontumor)) # Unimos en una sola matriz la expresión de células tumorales y sanas
imputed_tpm_total <- imputed_tpm_total[,colnames(filtered_sce)] # Ordenan las células para que se queden en el orden que estaban en el objeto original `filtered_sce`
imputed_exprs_total <- cbind(assay(filtered_sce_tumor,"exprs"),assay(filtered_sce_nontumor,"exprs"))
imputed_exprs_total <- imputed_exprs_total[,colnames(filtered_sce)] # Ordenan las células para que se queden en el orden que estaban en el objeto original `filtered_sce`

# Creamos el nuevo objeto de tipo `sce` con la expresión génica imputada
imputed_sce <- SingleCellExperiment(
  assays = list(tpm = imputed_tpm_total, exprs = imputed_exprs_total),
  colData = colData(filtered_sce),
  rowData = rowData(filtered_sce))

saveRDS(imputed_sce,file.path(outDir,"imputed_sce.rds"))