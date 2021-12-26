#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes
library(scater)
library(Rtsne)

# Opciones
options(stringsAsFactors = FALSE)
argumento <- commandArgs()
argumento <- argumento[6]
# argumento <- "head_neck"

outDir <- file.path("./datasets",argumento)
if(!dir.exists(outDir)) {                       # Crea la carpeta ./datasets/<head_neck o melanoma>/  si no existe
  dir.create(outDir, recursive = TRUE)
}


# Leemos el dataset del head_neck/melanoma con la expresión génica imputada y de
# ahí generamos dos objetos sce: uno con células tumorales y los genes
# metabólicos, y otro para células sanas y genes metabólicos
imputed_sce <- readRDS(file.path("../2-Imputation/datasets",argumento,"imputed_sce.rds"))
tumor_metabolic_sce <- imputed_sce[rowData(imputed_sce)$metabolic, imputed_sce$cellType == "Malignant"]
healthy_metabolic_sce <- imputed_sce[rowData(imputed_sce)$metabolic, imputed_sce$cellType != "Malignant"]

# Limpieza de RAM
rm(imputed_sce)
gc(verbose = FALSE)



####################################################################################################

###############################################################
###########     1. t-SNE de células tumorales       ###########
###############################################################

# Establecemos semilla de aleatoriedad
set.seed(12345)

# Calculamos un t-SNE exacto (theta = 0)
tsne_tumor <- Rtsne(t(assay(tumor_metabolic_sce,"exprs")),
                    initial_dims = 20,theta = 0.0, perplexity = 30)

# Creamos un dataframe temporal con las coordenadas del t-SNE y el tumor de
# procedencia de cada célula
tmp <- data.frame(x = tsne_tumor$Y[,1], y = tsne_tumor$Y[,2], 
                  procedencia = colData(tumor_metabolic_sce)$tumor)

# Computamos y guardamos la visualización del t-SNE
visualizacion_tsne <- ggplot(tmp) + geom_point(aes(x, y, colour = procedencia), size = 1) +
                      labs(x = "t-SNE1",y = "t-SNE2") + theme_bw() + ggtitle("t-SNE de células tumorales")

ggsave(file.path(outDir,"tumor_metabolic_tsne.pdf"), visualizacion_tsne, 
       width = 5, height = 4)



####################################################################################################

##############################################################
###########     2. t-SNE de células normales       ###########
##############################################################

# Calculamos un t-SNE exacto (theta = 0)
tsne_no_tumor <- Rtsne(t(assay(healthy_metabolic_sce,"exprs")),
                    initial_dims = 20,theta = 0.0, perplexity = 30)

# Creamos un dataframe temporal con las coordenadas del t-SNE, el tumor de
# procedencia y la estirpe de cada célula
tmp <- data.frame(x = tsne_no_tumor$Y[,1], y = tsne_no_tumor$Y[,2], 
                  procedencia = colData(healthy_metabolic_sce)$tumor,
                  linaje = colData(healthy_metabolic_sce)$cellType)


# Computamos y guardamos la visualización del t-SNE, coloreando las células
# según su linaje
visualizacion_tsne <- ggplot(tmp) + geom_point(aes(x, y, colour = linaje), size = 1) +
  labs(x = "t-SNE1",y = "t-SNE2") + theme_bw() + ggtitle("t-SNE de células no tumorales (linaje)")

ggsave(file.path(outDir,"healthy_metabolic_tsne.pdf"), visualizacion_tsne, 
       width = 5, height = 4)

# Ídem, pero ahora coloreando en función del tumor de procedencia
visualizacion_tsne <- ggplot(tmp) + geom_point(aes(x, y, colour = procedencia), size = 1) +
  labs(x = "t-SNE1",y = "t-SNE2") + theme_bw() + ggtitle("t-SNE de células no tumorales (tumor de procedencia)")

ggsave(file.path(outDir,"healthy_metabolic_tsne2.pdf"), visualizacion_tsne, 
       width = 7, height = 5)




