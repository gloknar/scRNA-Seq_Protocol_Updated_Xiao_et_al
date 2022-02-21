############
#### pruebas dropout   BORRAME LUEGO
# El dropout se mide por genes, un gen con dropout del 10%, del 20%... etc

library(scImpute)
library(scater)
options(stringsAsFactors = F)
argumento <- "melanoma"
outDir <- file.path("datasets",argumento)
if(!dir.exists(outDir) ) {dir.create(outDir, recursive = TRUE)}

before_imp <- readRDS(file.path("../1-ReadData/datasets",argumento,"filtered_sce.rds"))
before_imp <- before_imp[, before_imp$cellType == "B cell"]
dim(before_imp) # Tenemos 23684 genes y 515 linfos B
# View(filtered_sce_BCell@assays@data$tpm)

# De esas x celulas, calcular el porcentaje de dropout
before_imp <- tpm(before_imp)
gc()


# Vamos a calcular el porcentaje de droppout en el primer gen

# NOTA: Imputation was only applied to genes with dropout rates (i.e. the
# fraction of cells in which the corresponding gene has zero expression value)
# larger than 50% to avoid over-imputation

ratio_dropout_gen1 <- c()
aux <- c()

zeros_gen1 <- length(which(before_imp[1,] == 0))


for (i in c(1:nrow(bcell_tpm))) {  # Leemos la matriz por filas para calcular el % de dropout de cada gen
  aux <- length(which(bcell_tpm[i,] == 0))/2176
  # print(aux)
  tasas_dropout <- c(tasas_dropout, aux)
}


# LO TENGO
h = hist(tasas_dropout, breaks = 10) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE, col = "cyan", ylim = c(0,1), xlab = "Dropout rate", ylab = "Percentage of cell number", main = "B cell", las = 1)






# Despues del imputado
rm(h, bcell_tpm)
gc()
despues_imp <- read.csv(file.path(outDir,"non_malignant","scimpute_count.csv"), row.names = 1)
despues_imp[,filtered_sce$cellType == "B cell"]
cosas <- colnames(filtered_sce[, filtered_sce$cellType == "B cell"])
despues_imp <- despues_imp[,cosas]
# despues_imp <- readRDS(file.path(outDir,"imputed_sce.rds"))
View(despues_imp)

despues_imp$cellType <- factor(despues_imp$cellType)
despues_imp <- despues_imp[, despues_imp$cellType == "Malignant"]

despues_imp_tpm <- tpm(despues_imp)
rm(despues_imp)
gc()



dim(despues_imp_tpm)
tasas_dropout <- c()
aux <- c()

for (i in c(1:nrow(despues_imp))) {
  aux <- length(which(despues_imp[i,] == 0))/515
  # print(aux)
  tasas_dropout <- c(tasas_dropout, aux)
}

# LO TENGO
h = hist(tasas_dropout, breaks = 10) # or hist(x,plot=FALSE) to avoid the plot of the histogram
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE, col = "cyan", ylim = c(0,1), xlab = "Dropout rate", ylab = "Percentage of cell number", main = "B cell", las = 1)
