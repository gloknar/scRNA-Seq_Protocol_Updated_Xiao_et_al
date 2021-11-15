#### Instalación paquetes en R 4.1

# Paralelizado instalación paquetes mediante `make` y `Ncpus` (4 hilos)
Ncpus <- 4L
Sys.setenv("MAKE" = "make -k -j 4")
# Sys.getenv("MAKE")

# Instalamos scater, scran, Rstne y biomaRt desde BiocManager (versiones de
# R>3.5)
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("Rtsne")
BiocManager::install("biomaRt")



#scimpute
if(library(devtools,logical.return=T)){
  install_github("Vivianstats/scImpute",force=TRUE)
}else{
  warning("please install devtools and try again")
  #install.packages("devtools")
}

#pheatmap
if(!library(pheatmap,logical.return=T)){
  warning("please install pheatmap and try again")
  #install.packages("pheatmap")
}

#ggrepel
if(!library(ggrepel,logical.return=T)){
  warning("please install pheatmap and try again")
  #install.packages("ggrepel")
}

# download gsea 4.1.0
if(!file.exists("GSEA_Linux_4.1.0.zip")){
  warning("please download gsea-4.1.0.zip from http://software.broadinstitute.org/gsea/downloads.jsp")
}

# Alternativamente, descarga el script actualizado de aquí (hace exactamente lo
# mismo que el archivo de arriba): https://github.com/GSEA-MSigDB/GSEA_R