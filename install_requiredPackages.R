#### Instalación paquetes en R 4.1

# Paralelizado instalación paquetes mediante `make` (Ubuntu) y `Ncpus`
# (cualquier SO)
Ncpus <- 4L
# Sys.setenv("MAKE" = "make -k -j 4")  # Activar si usar Ubuntu
# Sys.getenv("MAKE")



# Instalamos scater, scran, Rstne y biomaRt
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("Rtsne")
BiocManager::install("biomaRt")


# Instalamos scImpute
if(library(devtools,logical.return=T)){              # Si está instalado `devtools`, usa su comando para instalar scImpute desde su repo
  install_github("Vivianstats/scImpute",force=TRUE)
}else{
  warning("please install devtools and try again")   # Si no está instalado `devtools`, lanza una advertencia
  #install.packages("devtools")
}

# Instalamos pheatmap
if(!library(pheatmap,logical.return=T)){             # Si no está instalado pheatmap, lanza una advertencia
  warning("please install pheatmap and try again")
  #install.packages("pheatmap")
}

# Instalamos ggrepel
if(!library(ggrepel,logical.return=T)){              # Si no está instalado ggrepel, lanza una advertencia
  warning("please install pheatmap and try again")
  #install.packages("ggrepel")
}

# Instalamos car
if(!library(car,logical.return=T)){              # Si no está instalado car, lanza una advertencia
  warning("please install car and try again")
  #install.packages("car")
}


# Descargamos GSEA 4.1.0
if(!file.exists("GSEA_Linux_4.1.0.zip")){            # Si no hemos descargado el archivo comprimido de GSEA, lanza una advertencia
  warning("please download gsea-4.1.0.zip from http://software.broadinstitute.org/gsea/downloads.jsp")
}

# Alternativamente, descarga el script actualizado de aquí (hace exactamente lo
# mismo que el archivo de arriba): https://github.com/GSEA-MSigDB/GSEA_R