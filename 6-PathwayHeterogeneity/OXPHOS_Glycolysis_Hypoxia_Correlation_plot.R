#######################################################################
###########     0. Carga de paquetes, opciones y datos      ###########
#######################################################################

# Paquetes y funciones auxiliares
library(scater)
library(stringr)
library(RColorBrewer)
source("../utils.R")


# Opciones
options(stringsAsFactors = FALSE)
# argumento <- commandArgs()
# argumento <- argumento[6]
argumento <- "melanoma"
outDir <- file.path("./datasets",argumento,"oxphos-gly-hyp-corr")
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
all_pathways <- gmtPathways(hallmark_gmt)

# Cargamos las rutas metabólicas a estudiar
genes_OXPHOS <- pathways[["Oxidative phosphorylation"]]
genes_glicolisis <- pathways[["Glycolysis / Gluconeogenesis"]]
genes_hipoxia <- all_pathways[["HALLMARK_HYPOXIA"]]



####################################################################################################

######################################################################################################
###########     1. Correlación entre hipoxia, glicólisis y OXPHOS en células tumorales     ###########
######################################################################################################

# Hacemos un subset del objeto sce original para quedarnos sólo con las células
# tumorales
tumor_sce <- filtered_sce[, filtered_sce$cellType == "Malignant"]
tumor_sce$tumor <- droplevels(tumor_sce$tumor)
neoplasias <- unique(tumor_sce$tumor)


# Vamos a calcular con la función auxiliar `num_of_pathways()` el nº de rutas
# metabólicas en las que participan nuestros genes de interés (1566 genes
# metabólicos)
gene_pathway_num <- num_of_pathways(ruta_archivo_pathways, intersect(unlist(pathways), rownames(filtered_sce)))


# Calculamos la expresión media de los genes de cada proceso metabólico
all_gene_exp <- assay(tumor_sce, "exprs")

oxphos_exp <- all_gene_exp[rownames(all_gene_exp) %in% genes_OXPHOS,]
glycolysis_exp <- all_gene_exp[rownames(all_gene_exp) %in% genes_glicolisis,]
hypoxia_exp <- all_gene_exp[rownames(all_gene_exp) %in% genes_hipoxia,]

oxphos <- colMeans(as.matrix(oxphos_exp), na.rm = T)
glycolysis <- colMeans(as.matrix(glycolysis_exp), na.rm = T)
hypoxia <- colMeans(as.matrix(hypoxia_exp), na.rm = T)
data <- data.frame(OXPHOS = oxphos, Glycolysis = glycolysis, Hypoxia = hypoxia)


# Calculamos la matriz de correlación para las 3 rutas metabólicas
print("Correlación de las rutas metabólicas de interés:")
print(cor(data))

#correlation plot for each two of them
dat_min <- 0
dat_max <- 4
p=ggplot(data,aes(x=OXPHOS,y=Glycolysis)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  labs(x = "OXPHOS", y = "Glycolysis") +
  theme(axis.line=element_line(size=0.3,colour="black"),
       axis.ticks = element_line(size=0.3,color="black"),
       axis.text.x=element_text(size=6),
       axis.text.y=element_text(size=6),
       axis.title.x=element_text(size=8),
       axis.title.y=element_text(size=8))

ggsave(filename = file.path(outDir,"malignant_oxphos_glycolysis.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

p=ggplot(data,aes(x=OXPHOS,y=Hypoxia)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  labs(x = "OXPHOS", y = "Hypoxia") +
  theme(axis.line=element_line(size=0.3,colour="black"),
        axis.ticks = element_line(size=0.3,color="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))
  

ggsave(filename = file.path(outDir,"malignant_oxphos_hypoxia.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

p=ggplot(data,aes(x=Glycolysis,y=Hypoxia)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  labs(x = "Glycolysis", y = "Hypoxia") +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  theme(axis.line=element_line(size=0.3,colour="black"),
        axis.ticks = element_line(size=0.3,color="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))

ggsave(filename = file.path(outDir,"malignant_glycolysis_hypoxia.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

# #correlation in each tumor
# for( i in neoplasias){
#   each_exp <- assay(tumor_sce[,tumor_sce$tumor==i],"exprs")
#   oxphos_exp <- each_exp[rownames(each_exp)%in% genes_OXPHOS,]
#   glycolysis_exp <- each_exp[rownames(each_exp)%in% genes_glicolisis,]
#   hypoxia_exp <- each_exp[rownames(each_exp)%in% genes_hipoxia,]
  
#   oxphos <- colMeans(as.matrix(oxphos_exp),na.rm=T)
#   glycolysis <- colMeans(as.matrix(glycolysis_exp),na.rm=T)
#   hypoxia <- colMeans(as.matrix(hypoxia_exp),na.rm=T)
#   data <- data.frame(OXPHOS=oxphos,Glycolysis=glycolysis,Hypoxia=hypoxia)
#   print(paste0("correlations in patient:",i))
#   print(cor(data))
# }




####################################################################################################

##################################################################################################
###########     2. Correlación entre hipoxia, glicólisis y OXPHOS en células sanas     ###########
##################################################################################################

healthy_sce <- filtered_sce[,filtered_sce$cellType!="Malignant"]
cell_types <- unique(healthy_sce$cellType)
cor_matrix <- matrix(NA,nrow=length(cell_types),ncol=3,dimnames = list(cell_types,c("oxphos_glycolosis","oxphos_hypoxia","glycolysis_hypoxia")))
for(c in cell_types){
  each_exp <- assay(healthy_sce[,healthy_sce$cellType==c],"exprs")
  oxphos_exp <- each_exp[rownames(each_exp)%in% genes_OXPHOS,]
  glycolysis_exp <- each_exp[rownames(each_exp)%in% genes_glicolisis,]
  hypoxia_exp <- each_exp[rownames(each_exp)%in% genes_hipoxia,]
  
  oxphos <- colMeans(as.matrix(oxphos_exp),na.rm=T)
  glycolysis <- colMeans(as.matrix(glycolysis_exp),na.rm=T)
  hypoxia <- colMeans(as.matrix(hypoxia_exp),na.rm=T)
  data <- data.frame(OXPHOS=oxphos,Glycolysis=glycolysis,Hypoxia=hypoxia)
  #calculate correlation
  cor_matrix[c,1] <- cor(data[,1:2])[1,2]
  cor_matrix[c,2] <- cor(data[,c(1,3)])[1,2]
  cor_matrix[c,3] <- cor(data[,2:3])[1,2]
}
write.table(cor_matrix,file.path(outDir,"non-malignant_Oxphos_Glycolysis_Hypoxia_cor.txt"),row.names = T,col.names=T,sep="\t")
