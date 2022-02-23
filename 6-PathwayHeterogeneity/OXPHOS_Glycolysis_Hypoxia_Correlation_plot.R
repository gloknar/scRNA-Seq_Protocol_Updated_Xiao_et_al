#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(pheatmap)
library(RColorBrewer)
source("../utils.R")


# Opciones
options(stringsAsFactors = FALSE)
argumentos <- commandArgs(trailingOnly = T)
argumento <- argumentos[1]  # "head_neck" o "melanoma"
outDir <- file.path("./datasets",argumento,"OXPHOS-gly-hyp-correlation")
if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}


# Leemos el dataset del head_neck/melanoma
filtered_sce <- readRDS(file.path("../1-ReadData/datasets/",argumento,"filtered_sce.rds"))


# Leemos el archivo de las rutas en las que participan los 1566 genes
# metabólicos (este contiene 85 rutas metabólicas)
ruta_archivo_pathways <- "../Data/KEGG_metabolism.gmt"
pathways <- gmtPathways(ruta_archivo_pathways)

# Este contiene 50 rutas metabólicas (parece que contiene los genes que
# participan en respuesta a hipoxia)
hallmark_gmt <- '../Data/h.all.v6.1.symbols.gmt'
hallmarks <- gmtPathways(hallmark_gmt)

# Cargamos las rutas metabólicas a estudiar
genes_OXPHOS <- pathways[["Oxidative phosphorylation"]]
genes_glicolisis <- pathways[["Glycolysis / Gluconeogenesis"]]
genes_hipoxia <- hallmarks[["HALLMARK_HYPOXIA"]]



####################################################################################################

######################################################################################################
###########     1. Correlación entre hipoxia, glicólisis y OXPHOS en células tumorales     ###########
######################################################################################################

# Hacemos un subset del objeto sce original para quedarnos sólo con las células
# tumorales
tumor_sce <- filtered_sce[, filtered_sce$cellType == "Malignant"]
tumor_sce$tumor <- factor(tumor_sce$tumor)  # El código es más robusto si hacemos volvemos a castear a factor en vez de usar droplevels() para eliminar niveles no usados del factor en cuestión
neoplasias <- unique(tumor_sce$tumor)


# Vamos a calcular con la función auxiliar `num_of_pathways()` el nº de rutas
# metabólicas en las que participan nuestros genes de interés (1566 genes
# metabólicos)
gene_pathway_num <- num_of_pathways(ruta_archivo_pathways, intersect(unlist(pathways), rownames(filtered_sce)))


# Calculamos la expresión media de los genes de cada proceso metabólico
all_genes_expr <- assay(tumor_sce, "exprs")

expr_oxphos <- all_genes_expr[rownames(all_genes_expr) %in% genes_OXPHOS,]         # Matriz TPM de los genes que participan en OXPHOS
expr_glicolisis <- all_genes_expr[rownames(all_genes_expr) %in% genes_glicolisis,] # Ídem para glicólisis
expr_hipoxia <- all_genes_expr[rownames(all_genes_expr) %in% genes_hipoxia,]       # Y para hipoxia

oxphos <- colMeans(as.matrix(expr_oxphos), na.rm = T)                          # Para cada célula, calculamos la expresion media de los genes que participan en la ruta metabólica OXHPOS
glicolisis <- colMeans(as.matrix(expr_glicolisis), na.rm = T)                  # Lo mismo para glicólisis
hipoxia <- colMeans(as.matrix(expr_hipoxia), na.rm = T)                        # Y para hipoxia
data <- data.frame(OXPHOS = oxphos, Glicolisis = glicolisis, Hipoxia = hipoxia)   # Agregamos todo en un dataframe


# Calculamos la matriz de correlación para las 3 rutas metabólicas
print("Correlación de las rutas metabólicas de interés:")
print(cor(data, method = "pearson"))  # Usa pearson porque tenemos más de 30 células, y podemos asumir que ambas distribuciones se aproximan a una normal... puedes comprobarlo si creas un objeto aov y lo ploteas, te saldra el qqplot


# Graficamos la regresión lineal para OXPHOS y glicólisis
data_min <- 0   # Rangos máximos y mínimos de dataframe
data_max <- 4

# Calculamos la ecuación de la regresión lineal para mostrarla en el gráfico
a <- lm(Glicolisis ~ OXPHOS, data)
b <- summary(a)
formula1 <- paste0("y = ",round(a$coefficients[[1]], 2), " + ", round(a$coefficients[[2]], 2), "x ;  R = ", round(sqrt(b$adj.r.squared), 2))
p_valor <- pf(q = b$fstatistic[[1]], df1= b$fstatistic[[2]], 
              df2= b$fstatistic[[3]], lower.tail = F)  # Ya que el summary de un objeto lm nos enseña el p-valor pero no nos los devuelve de manera programática, lo calculamos a partir del estadístico F
formula2 <- paste0("p-value = ", format(p_valor, digits = 3))


# Generamos el gráfico y lo guardamos en un pdf
p = ggplot(data, aes(x = OXPHOS, y = Glicolisis)) + 
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    xlim(data_min, data_max) + ylim(data_min, data_max) +
    theme_classic()  + theme(aspect.ratio = 0.8) +
    labs(x = "OXPHOS", y = "Glicolisis") +
    geom_text(x = 1, y = 3.5, label = formula1) +
    geom_text(x = .75, y = 3, label = formula2) +
    theme(axis.line = element_line(size = 0.3, colour = "black"),
          axis.ticks = element_line(size = 0.3, color = "black"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))

ggsave(filename = file.path(outDir,"malignant_OXPHOS_Glycolysis.pdf"), p, 
       device = "pdf", width = 6, height = 4, units = "in",
       useDingbats = FALSE)   # Evitamos usar la fuente Dingbats porque según la documentación de ggplot2, a veces da problemas



# Ídem para OXPHOS e hipoxia
a <- lm(Hipoxia ~ OXPHOS, data)
b <- summary(a)
formula1 <- paste0("y = ",round(a$coefficients[[1]], 2), " + ", round(a$coefficients[[2]], 2), "x ;  R = ", round(sqrt(b$adj.r.squared), 2))
p_valor <- pf(q = b$fstatistic[[1]], df1= b$fstatistic[[2]], 
              df2= b$fstatistic[[3]], lower.tail = F)
formula2 <- paste0("p-value = ", format(p_valor, digits = 3))


p = ggplot(data, aes(x = OXPHOS, y = Hipoxia)) + 
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    xlim(data_min, data_max) + ylim(data_min, data_max) +
    theme_classic()  + theme(aspect.ratio = 0.8) +
    labs(x = "OXPHOS", y = "Hipoxia") +
    geom_text(x = 1, y = 3.5, label = formula1) +
    geom_text(x = .75, y = 3, label = formula2) +
    theme(axis.line = element_line(size = 0.3, colour = "black"),
          axis.ticks = element_line(size = 0.3, color = "black"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))
  
ggsave(filename = file.path(outDir,"malignant_OXPHOS_Hypoxia.pdf"), p, 
       device = "pdf", width = 6, height = 4, units = "in",
       useDingbats = FALSE)


# Ídem para glicólisis e hipoxia
a <- lm(Hipoxia ~ Glicolisis, data)
b <- summary(a)
formula1 <- paste0("y = ",round(a$coefficients[[1]], 2), " + ", round(a$coefficients[[2]], 2), "x ;  R = ", round(sqrt(b$adj.r.squared), 2))
p_valor <- pf(q = b$fstatistic[[1]], df1= b$fstatistic[[2]], 
              df2= b$fstatistic[[3]], lower.tail = F)
formula2 <- paste0("p-value = ", format(p_valor, digits = 3))



p = ggplot(data, aes(x = Glicolisis, y = Hipoxia)) + 
    geom_point(size = 0.5) +
    geom_smooth(method = "lm",color = "red") +
    xlim(data_min, data_max) + ylim(data_min, data_max) +
    labs(x = "Glicolisis", y = "Hipoxia") +
    theme_classic()  + theme(aspect.ratio = 0.8) +
    geom_text(x = 1, y = 3.5, label = formula1) +
    geom_text(x = .75, y = 3, label = formula2) +
    theme(axis.line = element_line(size = 0.3, colour = "black"),
          axis.ticks = element_line(size = 0.3, color = "black"),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8))

ggsave(filename = file.path(outDir,"malignant_Glycolysis_Hypoxia.pdf"), p, 
       device = "pdf", width = 6, height = 4, units = "in",
       useDingbats = FALSE)


# # Correlación entre cada tumor (= paciente)
# for( i in neoplasias){
#   each_expr <- assay(tumor_sce[,tumor_sce$tumor == i], "exprs")
#   expr_oxphos <- each_expr[rownames(each_expr) %in% genes_OXPHOS,]
#   expr_glicolisis <- each_expr[rownames(each_expr) %in% genes_glicolisis,]
#   expr_hipoxia <- each_expr[rownames(each_expr) %in% genes_hipoxia,]
# 
#   oxphos <- colMeans(as.matrix(expr_oxphos), na.rm = T)
#   glicolisis <- colMeans(as.matrix(expr_glicolisis), na.rm = T)
#   hipoxia <- colMeans(as.matrix(expr_hipoxia), na.rm = T)
#   data <- data.frame(OXPHOS = oxphos, Glicolisis = glicolisis, Hipoxia = hipoxia)
#   print(paste0("correlacones entre pacientes:",i))
#   print(cor(data))
#  }


# Limpieza RAM
rm(tumor_sce)
invisible(gc(verbose = F))



####################################################################################################

##################################################################################################
###########     2. Correlación entre hipoxia, glicólisis y OXPHOS en células sanas     ###########
##################################################################################################

# Hacemos un subset del objeto sce original para quedarnos sólo con las células
# sanas
healthy_sce <- filtered_sce[,filtered_sce$cellType != "Malignant"]
healthy_sce$cellType <-  factor(healthy_sce$cellType)
cell_types <- unique(healthy_sce$cellType)

# Limpieza RAM, ya no necesitamos el objeto filtered_sce
rm(filtered_sce)
invisible(gc(verbose = F))


# Inicializamos una matriz de correlaciones vacía de tipos celulares X
# correlaciones y la rellenamos
matriz_corr <- matrix(NA, nrow = length(cell_types), ncol = 3,
                     dimnames = list(cell_types, c("oxphos_glicolosis","oxphos_hipoxia","glicolisis_hipoxia")))

for (c in cell_types) {
  each_expr <- assay(healthy_sce[,healthy_sce$cellType == c], "exprs")     # Obtenemos la matriz TPM de las células del linaje de interés
  expr_oxphos <- each_expr[rownames(each_expr) %in% genes_OXPHOS,]         # Obtenemos de ahi la matriz TPM de los genes de OXPHOS
  expr_glicolisis <- each_expr[rownames(each_expr) %in% genes_glicolisis,] # Ídem para glicólisis
  expr_hipoxia <- each_expr[rownames(each_expr) %in% genes_hipoxia,]       # Y para hipoxia
  
  # Para cada célula, calculamos la expresion media (i.e. actividad de la ruta)
  # de los genes que participan en la ruta metabólica OXHPOS, glicólisis e
  # hipoxia, respectivamente
  oxphos <- colMeans(as.matrix(expr_oxphos), na.rm = T)       
  glicolisis <- colMeans(as.matrix(expr_glicolisis), na.rm = T)
  hipoxia <- colMeans(as.matrix(expr_hipoxia), na.rm = T)
  data <- data.frame(OXPHOS = oxphos, Glicolisis = glicolisis, Hipoxia = hipoxia)  # Agregamos todo en un dataframe
  
  # Rellenamos la matriz de correlaciones previamente inicializada
  matriz_corr[c,1] <- cor(data[, c("OXPHOS", "Glicolisis")])[1,2]
  matriz_corr[c,2] <- cor(data[,c("OXPHOS", "Hipoxia")])[1,2]
  matriz_corr[c,3] <- cor(data[,c("Glicolisis", "Hipoxia")])[1,2]
}

# Guardamos las correlaciones entre las 3 rutas de interés en cada linaje
# celular en un archivo tsv
write.table(matriz_corr, 
    file.path(outDir,"non-malignant_OXPHOS_Glycolysis_Hypoxia_corr.tsv"),
    row.names = T, col.names = NA, sep = "\t")


# Generamos el archivo PDF donde guardaremos el heatmap
pdf(file.path(outDir,"non-malignant_pathways_corr_heatmap.pdf"), onefile = T,
    width = 6, height = 9)

# Generamos la paleta de colores: un gradiente de azul a rojo (pasando por el
# blanco) de 100 pasos
color <- colorRampPalette(c("white", "red"))(50)

# Le decimos al heatmap que use esa paleta de manera gradual, siendo el 2 rojo,
# el 1 blanco y el 0 azul
mybreaks <- c(
  seq(0, 0.5, length.out = 25),
  seq(0.51, max(matriz_corr), length.out = 25))

# Computamos el heatmap
pheatmap(matriz_corr, cluster_cols = F, cluster_rows = F, 
         color = color, breaks = mybreaks,
         main = "Persistencia de la correlación en células sanas")

# Lo guardamos en el disco duro
dev.off()


# Mensaje de fin
print("")
print("GRACIAS POR ASISTIR A MI CHARLA TED")
