#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes
library(scater)
library(umap)

# Opciones
options(stringsAsFactors = FALSE)
argumentos <- commandArgs(trailingOnly = T)
argumento <- argumentos[1]  # "head_neck" o "melanoma"
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
invisible(gc(verbose = FALSE))



####################################################################################################

##############################################################
###########     1. UMAP de células tumorales       ###########
##############################################################

# Establecemos semilla de aleatoriedad
set.seed(12345)

# Configuramos y calculamos el UMAP
umap_config <- umap.defaults
umap_config$n_components <- 2
umap_config$random_state <- 12345
umap_config$metric <- "euclidean"

umap_tumor <- umap(t(assay(tumor_metabolic_sce, "exprs")),
                   config = umap_config, 
                   method = "naive")

# Creamos un dataframe temporal con las coordenadas del UMAP y el tumor de
# procedencia de cada célula
tmp <- data.frame(x = umap_tumor$layout[,1], y = umap_tumor$layout[,2], 
                  procedencia = colData(tumor_metabolic_sce)$tumor)

# Computamos y guardamos la visualización del UMAP
visualizacion_umap <- ggplot(tmp) + geom_point(aes(x, y, colour = procedencia), size = 1) +
                      labs(x = "UMAP1",y = "UMAP2") + theme_bw() + ggtitle("UMAP de células tumorales")

ggsave(file.path(outDir,"tumor_metabolic_umap.pdf"), visualizacion_umap, 
       width = 5, height = 4)



####################################################################################################

#############################################################
###########     2. UMAP de células normales       ###########
#############################################################

# Computamos el UMAP
umap_no_tumor <- umap(t(assay(healthy_metabolic_sce, "exprs")),
                   config = umap_config, 
                   method = "naive")

# Creamos un dataframe temporal con las coordenadas del UMAP, el tumor de
# procedencia y la estirpe de cada célula
tmp <- data.frame(x = umap_no_tumor$layout[,1], y = umap_no_tumor$layout[,2], 
                  procedencia = colData(healthy_metabolic_sce)$tumor,
                  linaje = colData(healthy_metabolic_sce)$cellType)


# Computamos y guardamos la visualización del UMAP, coloreando las células
# según su linaje
visualizacion_umap <- ggplot(tmp) + geom_point(aes(x, y, colour = linaje), size = 1) +
  labs(x = "UMAP1",y = "UMAP2") + theme_bw() + ggtitle("UMAP de células no tumorales (linaje)")

ggsave(file.path(outDir,"healthy_metabolic_umap.pdf"), visualizacion_umap, 
       width = 5, height = 4)

# Ídem, pero ahora coloreando en función del tumor de procedencia
visualizacion_umap <- ggplot(tmp) + geom_point(aes(x, y, colour = procedencia), size = 1) +
  labs(x = "UMAP1",y = "UMAP2") + theme_bw() + ggtitle("UMAP de células no tumorales (tumor de procedencia)")

ggsave(file.path(outDir,"healthy_metabolic_umap2.pdf"), visualizacion_umap, 
       width = 7, height = 5)


# Mensaje de fin
print("")
print("GRACIAS POR ASISTIR A MI CHARLA TED")
