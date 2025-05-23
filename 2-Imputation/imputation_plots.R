############
#### Script para sacar las gráficas de los imputados
# El dropout se mide por genes, un gen con dropout del 10%, del 20%... etc

# Uso: imputation_plots.R head_neck  # "melanoma" o "head_neck"

library(scImpute)
library(scater)
options(stringsAsFactors = F)
argumentos <- commandArgs(trailingOnly = T)
argumento <- as.character(argumentos[1])   # "melanoma" o "head_neck"
outDir <- file.path("graficos_imputacion",argumento)
if(!dir.exists(outDir) ) {dir.create(outDir, recursive = TRUE)}


# Cargamos datos antes y despues del imputado
before_imp <- readRDS(file.path("../1-ReadData/datasets",argumento,"filtered_sce.rds"))
before_imp$cellType <- as.factor(before_imp$cellType)
before_imp@assays@data$tpm[before_imp@assays@data$tpm <= 1] = 0   # Ponemos como no detectado (0) genes con menos de 1 TPM

after_imp <- readRDS(file.path("datasets",argumento,"imputed_sce.rds"))
after_imp$cellType <- as.factor(after_imp$cellType)
after_imp@assays@data$tpm[after_imp@assays@data$tpm <= 1] = 0   # Ponemos como no detectado (0) genes con menos de 1 TPM




# Para calcular las gráficas ANTES del imputado
for (celulas in levels(before_imp$cellType)) {
  # Subset del sce
  celulas_before_imp <- before_imp[, before_imp$cellType == celulas]
  # Matriz de expresión génica de cada tipo celular
  tpm_celulas <- tpm(celulas_before_imp)
  # Tamaño de pobación celular
  n_celulas <- dim(celulas_before_imp)[2]
  
  # Calculamos para todos los genes su % de dropout
  zero <- as.numeric(0)
  ratio_dropout_genes <- matrixStats::rowCounts(x = tpm_celulas, value = zero)/n_celulas
  
  # Computamos el histograma interino
  histograma_tmp = hist(ratio_dropout_genes, breaks = 10, plot = F) # or hist(x,plot=FALSE) to avoid the plot of the histogram
  histograma_tmp$density = histograma_tmp$counts/sum(histograma_tmp$counts)
  
  # Computamos y guardamos el histograma final
  nombre_grafico <- paste0("before_imputation_",celulas,".png")
  png(filename=file.path(outDir,nombre_grafico))
  
  plot(histograma_tmp, freq=FALSE, col = "cyan", ylim = c(0,1), 
       xlab = "Dropout rate", ylab = "Percentage of cell number", 
       main = as.character(celulas), las = 1)
  
  dev.off()
  
}






# Para calcular las gráficas DESPUES del imputado
for (celulas in levels(after_imp$cellType)) {
  # Subset del sce
  celulas_after_imp <- after_imp[, after_imp$cellType == celulas]
  # Matriz de expresión génica de cada tipo celular
  tpm_celulas <- tpm(celulas_after_imp)
  # Tamaño de pobación celular
  n_celulas <- dim(celulas_after_imp)[2]
  
  # Calculamos para todos los genes su % de dropout
  zero <- as.numeric(0)
  ratio_dropout_genes <- matrixStats::rowCounts(x = tpm_celulas, value = zero)/n_celulas
  
  # Computamos el histograma interino
  histograma_tmp = hist(ratio_dropout_genes, breaks = 10, plot = F) # or hist(x,plot=FALSE) to avoid the plot of the histogram
  histograma_tmp$density = histograma_tmp$counts/sum(histograma_tmp$counts)
  
  # Computamos y guardamos el histograma final
  nombre_grafico <- paste0("after_imputation_",celulas,".png")
  png(filename=file.path(outDir,nombre_grafico))
  
  plot(histograma_tmp, freq=FALSE, col = "cyan", ylim = c(0,1), 
       xlab = "Dropout rate", ylab = "Percentage of cell number", 
       main = as.character(celulas), las = 1)
  
  dev.off()
  
}


print("")
print("GRACIAS POR ASISTIR A MI CHARLA TED")
