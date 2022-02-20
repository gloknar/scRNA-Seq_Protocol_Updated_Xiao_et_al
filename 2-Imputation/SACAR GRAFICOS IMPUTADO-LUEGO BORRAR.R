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


# De esas x celulas, calcular el porcentaje de dropout
before_imp <- tpm(before_imp)
gc()


# Imputation was only applied to genes with dropout rates (i.e. the fraction of
# cells in which the corresponding gene has zero expression value) larger than
# 50% to avoid over-imputation
resultados <- c()
# for (i in c(1:nrow(before_imp))) {  # Leemos la matriz por filas para calcular el % de dropout de cada gen
#   zeros_gen1 <- length(which(before_imp[i,] == 0))
#   dropout_gen1 <- zeros_gen1/dim(before_imp)[2]
#   resultados <- c(resultados,dropout_gen1)
# }

# Hacer lo del bucle sin bucle. MUCHO MAS EFICIENTE EN R
resultados <- rowCounts(before_imp,value = 0)/dim(before_imp)[2] # Refactorizado de las dos lineas anteriores

histograma <- hist(resultados, breaks = 10, freq = T)
histograma$density <-histograma$counts/sum(histograma$counts)
plot(histograma, freq = F, xlab = "Dropout rate", ylab = "Percentage of cell number", 
     main = "B cell", col = "cyan", las = 1, ylim = c(0,1))






after_imp <- readRDS(file.path("../2-Imputation/datasets/",argumento,"imputed_sce.rds"))
after_imp <- after_imp[, after_imp$cellType == "B cell"]
dim(after_imp) # Tenemos 23684 genes y 515 linfos B

# De esas x celulas, calcular el porcentaje de dropout
after_imp <- tpm(after_imp)
gc()


# Imputation was only applied to genes with dropout rates (i.e. the fraction of
# cells in which the corresponding gene has zero expression value) larger than
# 50% to avoid over-imputation
resultados <- c()

# Hacer lo del bucle sin bucle. MUCHO MAS EFICIENTE EN R
resultados <- rowCounts(after_imp, value = 0)/dim(after_imp)[2] # Refactorizado de las dos lineas anteriores

histograma <- hist(resultados, breaks = 10, freq = T)
histograma$density <-histograma$counts/sum(histograma$counts)
plot(histograma, freq = F, xlab = "Dropout rate", ylab = "Percentage of cell number", 
     main = "B cell", col = "cyan", las = 1, ylim = c(0,1))
