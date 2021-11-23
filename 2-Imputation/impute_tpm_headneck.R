library(scImpute)
library(scater)

args <- commandArgs()
tumor <- "head_neck"  # args[6]
num_cores <- 6 #for windows it should be 1

selected_sce <- readRDS(file.path("../1-ReadData/dataset/head_neck/selected_sce.rds"))


outDir <- file.path("dataset/head_neck") # "dataset",tumor
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE) # Crea la carpeta ./datasets/head_neck/  si no existe

#impute the tumor and non-tumor seperately
selected_tumor_sce <- selected_sce[,selected_sce$cellType=="Malignant"]
selected_nontumor_sce <- selected_sce[,selected_sce$cellType!="Malignant"]


# NOTA: el bolsillo `selected_tumor_sce@assays@data$exprs` contiene la
# expresión génica de las células en formato log2(TMP+1);;
# El bolsillo `selected_tumor_sce@assays@data$tpm` contiene la expresión génica
# en formato TPM

#write the tpm matrix
selected_tumor_tpm <- tpm(selected_tumor_sce)
selected_nontumor_tpm <- tpm(selected_nontumor_sce) 
labels_tumor <- selected_tumor_sce$tumor
labels_nontumor <- selected_nontumor_sce$cellType

write.csv(selected_tumor_tpm,file.path(outDir,"tumor.tpm"))
write.csv(selected_nontumor_tpm,file.path(outDir,"nontumor.tpm"))

##prepare the gene length file (esto servía para pasar de TPM a read counts mediante TPM*gene length)
all_gene_lengths <- read.table("../Data/gene_length.txt",sep="\t",header=F,row.names=1)
tmp <- intersect(rownames(all_gene_lengths),rownames(selected_tumor_tpm)) # La intersección busca nombres de genes compartidos en ambos conjuntos de nombres de genes
if (length(tmp) != nrow(selected_tumor_tpm)){ # Compara si el nº de genes en selected_tumor_tpm y all_gene_lengths es igual
  warning("check the length file")
	print(setdiff(rownames(selected_tumor_tpm),rownames(all_gene_lengths)))
	q()
}

genelen <- all_gene_lengths[rownames(selected_tumor_tpm),]
genelen <- as.numeric(as.vector(genelen)) # Crea una matriz numérica de 23686 * 1 con la longitud de cada gen


# Imputamos genes con dropout >= 0.5 para evitar sobre-imputación
scimpute(file.path(outDir, "tumor.tpm"), infile = "csv", outfile = "csv",
         out_dir = file.path(outDir, "malignant_"), labeled = TRUE, 
         labels = as.vector(labels_tumor), type = "TPM", genelen = genelen, 
         drop_thre = 0.5, ncores = num_cores)


imputed_tpm <- read.csv(file.path(outDir,"malignant_scimpute_count.csv"),header=T,row.names=1)
tpm(selected_tumor_sce) <- data.matrix(imputed_tpm) 
assay(selected_tumor_sce,"exprs") <- data.matrix(log2(imputed_tpm + 1)) # actualiza los bolsillos tpm y exprs con los valores imputados de los genes



# ME HE QUEDADO AQUÍ. AQUI ES DONDE DA ERROR EL SCRIPT. Hay valores infinitos, igual es eso
# Imputación del non tumor
getwd()
outDir
list.files("dataset/head_neck/nontumor.tpm")
min(selected_nontumor_tpm)

anyMissing(selected_nontumor_sce@assays@data$exprs)
anyNA(selected_nontumor_sce@assays@data$exprs)
anyMissing(selected_nontumor_sce@assays@data$tpm)
anyNA(selected_nontumor_sce@assays@data$tpm)
range(selected_nontumor_sce@assays@data$tpm)
range(selected_nontumor_sce@assays@data$exprs)

range(selected_tumor_sce@assays@data$tpm) # Este no tiene inf
range(selected_tumor_sce@assays@data$exprs)



# DA ERROR PORQUE SE CREAN TPMs infinitos
scimpute(file.path(outDir, "nontumor.tpm"), infile = "csv", outfile = "csv",
         out_dir = file.path(outDir, "non-malignant_"), labeled = TRUE,
         labels = as.vector(labels_nontumor),	type = "TPM", genelen = genelen,
         drop_thre = 0.5, ncores = num_cores)

imputed_tpm <- read.csv(file.path(outDir,"non-malignant_scimpute_count.csv"),header=T,row.names=1)
tpm(selected_nontumor_sce) <- data.matrix(imputed_tpm) 
assay(selected_nontumor_sce,"exprs") <- data.matrix(log2(imputed_tpm + 1))

#save as sce
impute_tpm <- cbind(tpm(selected_tumor_sce), tpm(selected_nontumor_sce))
impute_exprs <- cbind(assay(selected_tumor_sce,"exprs"),assay(selected_nontumor_sce,"exprs"))
impute_tpm <- impute_tpm[,colnames(selected_sce)]
impute_exprs <- impute_exprs[,colnames(selected_sce)]

selected_impute_sce <- SingleCellExperiment(
  assays = list(tpm = impute_tpm, exprs=impute_exprs),
  colData = colData(selected_sce),
  rowData = rowData(selected_sce)
)
saveRDS(selected_impute_sce,file.path(outDir,"selected_impute_sce.rds"))
