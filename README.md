# Proyecto scRNA-seq de ovarios

En su [publicación](https://www.nature.com/articles/s41467-019-11738-0) de Nature Communications (2019), Xiao _et al._ caracterizaron 2 tipos de tumores a nivel de sc-RNA-Seq y GSEA (gene expression pathways) mediante el empleo del lenguaje de programación estadístico `R 3.5` y `bash`, y se puede visitar su código original en el siguiente repositorio: https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape.

En este repositorio hemos clarificado, refactorizado, optimizado y traducido el protocolo de Xiao y compañeros a la versión más actualizada de R a fecha de redacción de este documento (`R 4.1`). Emplearemos este protocolo modificado con nuestras muestras del [Karolinska Institutet](https://ki.se/en), las cuales son de tipo bulk RNA-seq y scRNA-seq de ovarios, para caracterizar el transcriptoma a nivel de célula única de todo el tejido ovárico y sus distintos tipos celulares.

Este protocolo ha sido probado en una estación de trabajo con Ubuntu 20.04 LTS, 1 CPU de 12 núcleos, 24 hilos @ 3.8GHz, 16GB de RAM y 1TB de HDD.

<br>

Adam Casas, Investigador bioinformático de la UGR, 2021 | ~~acm95@ugr.es~~ adam.casas95@gmail.com

------------

## El protocolo a vista de pájaro

Este pipeline permite analizar perfiles de expresión de genes metabólicos en el ámbito de scRNA-seq de ovarios. El pipeline se divide en 7 pasos, visualizables en el esquema a continuación. Para mantener el repositorio organizado, cada paso tiene asignada su propia carpeta. Los scripts del pipeline están escritos en `R`, y están ordenados de tal manera que el workflow se puede ejecutar desde `bash` con los comandos de las siguientes secciones.

![Esquema del pipeline para análisis de datos scRNA-seq](pipeline.png)



## Requisitos

Además de emplear `R 4.1`, se requiere una serie de paquetes los cuales pueden ser instalados mediante la siguiente línea de `bash`:

``` bash
Rscript install_requiredPackages.R 
```

## Descarga y lectura de los datasets

``` bash
cd "1-ReadData"

# Descargamos los datasets sucios
bash download_dataset.sh

# Preprocesamos los datasets
/opt/R/4.1.2/bin/Rscript readData_head_neck.R
/opt/R/4.1.2/bin/Rscript readData_melanoma.R

cd ../
```

Con el script de bash `download_dataset.sh` descargamos los datasets desde el GEO, y con sendos scripts de R corregimos erratas presentes en dichos datasets, filtramos grupos celulares con menos de 50 células, creamos 2 objetos de tipo `Single Cell Experiment` y los guardamos en los archivos `<head_neck o melanoma>/filtered_sce.rds` 

## Imputación de NAs

``` bash
cd "2-Imputation"

# Imputamos
/opt/R/4.1.2/bin/Rscript impute_tpm.R melanoma 16  # el 2º parámetro es el nº de hilos a usar para el imputado. Por defecto es 1
/opt/R/4.1.2/bin/Rscript impute_tpm.R head_neck 16


# Graficamos los resultados de antes y despues del imputado
/opt/R/4.1.2/bin/Rscript imputation_plots.R melanoma
/opt/R/4.1.2/bin/Rscript imputation_plots.R head_neck


cd ../
```

Este paso utiliza el paquete ["scImpute"](https://github.com/Vivianstats/scImpute) para imputar los valores faltantes de aquellos genes que presentan una expresión génica nula en > 50% de las células estudiadas (_i.e. dropout rate_ > 50%). Nótese que los genes con un dropout del 100% no son imputados por falta de información.

## Normalizado y evaluación de distintos métodos de normalización

``` bash
cd "3-Normalization"

# Normalizado de datos mediante distintos métodos
/opt/R/4.1.2/bin/Rscript normalization.R melanoma
/opt/R/4.1.2/bin/Rscript normalization.R head_neck

cd ../
```

En este paso evaluamos en nuestros datasets la eficacia de 4 métodos de normalizado, a saber:
* Upper-Quartile (EdgeR)
* TMM (EdgeR)
* RLE (DESeq2)
* Deconvolution (scran)

Tras normalizar los datos, graficamos los resultados para poder compararlos y seguir el protocolo con el método que mejor nos funcione (en este caso, será el método de deconvolución de scran).

## Clustering y visualización de las células en función de los genes metabólicos

``` bash
cd "4-Clustering"

# t-SNE
/opt/R/4.1.2/bin/Rscript tsne_metabolic_genes.R melanoma 18  # El 2º parámetro es el nº de hilos para paralelizar el tSNE. Por defecto es 1
/opt/R/4.1.2/bin/Rscript tsne_metabolic_genes.R head_neck 18

# UMAP (Opcional)
/opt/R/4.1.2/bin/Rscript umap_metabolic_genes.R melanoma
/opt/R/4.1.2/bin/Rscript umap_metabolic_genes.R head_neck

# Matriz de correlación de Spearman
/opt/R/4.1.2/bin/Rscript inter_tumor_correlation_matrix.R melanoma
/opt/R/4.1.2/bin/Rscript inter_tumor_correlation_matrix.R head_neck

cd ../
```

En este paso nos fijaremos en los 1566 genes metabólicos para agrupar y visualizar las células con t-SNE y UMAP (las visualizaciones pueden variar ligeramente respecto al gráfico del artículo científico debido a la inicialización aleatoria del t-SNE y UMAP). En concreto, estudiaremos por un lado las células malignas y por otro las células sanas. 

También calcularemos la matriz de correlación de Spearman (no paramétrico) de los distintos tumores dentro de una misma patología para mostrar su heterogeneidad. Los distintos melanomas presentes en nuestro dataset están clasificados como melanomas, pero sus perfiles de expresión de genes metabólicos demuestran que son distintos entre sí (ídem para las neoplasias de HNSCC).


## Heatmap y violinplot de la actividad metabólica en todos los tipos celulares (expresión bulk vs scRNA-seq)

``` bash
cd 5-PathwayActivity

# Datasets de scRNA-seq
Rscript scRNA_pathway_activity.R melanoma
Rscript scRNA_pathway_activity.R head_neck

# Dataset de bulk RNA-seq
Rscript TCGA_pathway_activity.R
cd ..
```

En este paso calculamos la actividad de las distintas rutas metabólicas (1566 genes metabs) en cada linaje celular, tanto en el dataset de bulk RNA-seq como en el de scRNA-seq.
Con el script `scRNA_pathway_activity.R` creamos un heatmap donde desglosamos en detalle la actividad de cada ruta metabólica de interés en nuestros tipos celulares y un violinplot donde se muestra la actividad global de las rutas en cada tipo celular. Así mismo, se genera una matriz con los p-valores para la actividad de cada ruta en cada tipo celular. Los p-valores se obtuvieron remuestreando los datos 5000 veces (test de permutación aleatoria).

Para calcular la actividad de cada ruta, se obtuvo en cada linaje celular la expresión media de los genes que la constituyen, luego se dividieron estas medias por la actividad media del gen a lo largo de todos los linajes celulares para calcular la actividad relativa (a lo fold change), se penalizaron con pesos los genes que participan en más de una ruta y se sumó la actividad de todos los genes de dicha ruta. 



## Análisis de rutas metabólicas

``` bash
cd 6-PathwayHeterogeneity

Rscript intra_malignant_heterogeneity.R melanoma
Rscript intra_malignant_heterogeneity.R head_neck

Rscript intra_non-malignant_heterogeneity.R melanoma
Rscript intra_non-malignant_heterogeneity.R head_neck

# Analizamos las rutas metabólicas OXPHOS, glicólisis y respuesta a hipoxia
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R melanoma
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R head_neck

# Obtenemos los genes diferencialmente expresados en células tumorales con 
# baja actividad en las 3 rutas metabólicas de interés
Rscript GeneSignature-of-OXPHOS_Glycolysis_Hypoxia melanoma
Rscript GeneSignature-of-OXPHOS_Glycolysis_Hypoxia head_neck
cd ..
```

En este paso hacemos un PCA y un GSEA para investigar la expresión diferencial de las rutas metabólicas de las células de las poblaciones malignas y sanas.

También realizamos scatterplots para analizar la actividad y correlaciones entre las rutas de OXPHOS, glicólisis y respuesta a hipoxia en células tumorales. 

Las rutas metabólicas de las células con una baja actividad OXPHOS/glycolisis/hypoxia se identifican y son guardadas en archivos de texto, los cuales se pueden usar como input de un análisis GO en http://metascape.org




## Perfil metabólico de subtipos celulares sanos

``` bash
cd 7-MetabolicFeatures
Rscript non-malignant_subtype.R melanoma
Rscript non-malignant_subtype.R head_neck
cd ..
```

Se caracteriza la expresión de las rutas metabólicas de diversos subtipos de linfocitos T y fibroblastos
