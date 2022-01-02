# Proyecto scRNA-seq de ovarios

Xiao _et al._ publicaron en 2019 un [artículo](https://www.nature.com/articles/s41467-019-11738-0) en Nature Communications en el que caracterizan 2 tipos de tumores a nivel de transcriptoma de célula única y GSEA (Gene Set Enrichment Analysis). El producto de dicha investigación está disponible en forma de [repositorio](https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape) con todo el código para `R 3.5` + `bash`

Nosotros queremos replicar dicho protocolo en nuestras muestras del Karolinska Institutet, las cuales son bulk RNA-seq y scRNA-seq de ovarios. Para ello he creado el presente repositorio, en el cual se clarifica, refactoriza, optimiza y traduce el protocolo de Xiao y compañeros a la versión más actualizada de R a fecha de redacción de este documento, `R 4.1`.


<br>

Adam C, 2021 | acm95@ugr.es

------------

## El protocolo a vista de pájaro

Este pipeline permite analizar perfiles de expresión de genes metabólicos en el ámbito de scRNA-seq de ovarios. El pipeline se divide en 7 pasos, visualizables en el esquema a continuación. Para mantener el repositorio organizado, cada paso tiene asignada su propia carpeta. Los scripts del pipeline están escritos en `R`, y están ordenados de tal manera que el workflow se puede ejecutar desde `bash` con los comandos de las siguientes secciones.

![Esquema del pipeline para análisis de datos scRNA-seq](pipeline.png)



## Requisitos

Además de emplear `R 4.1`, se requiere una serie de paquetes, los cuales pueden ser instalados mediante la siguiente línea de `bash`:

``` bash
Rscript install_requiredPackages.R 
```

## Descarga y lectura de los datasets

``` bash
cd "1-ReadData"
bash download_dataset.sh
Rscript readData_head_neck.R
Rscript readData_melanoma.R
cd ../
```

Con el script de bash `download_dataset.sh` descargamos los datasets desde el GEO, y con sendos scripts de R corregimos erratas presentes en dichos datasets, filtramos grupos celulares con menos de 50 células, creamos 2 objetos de tipo `Single Cell Experiment` y los guardamos en los archivos `<head_neck o melanoma>/filtered_sce.rds` 

## Imputación de NAs

``` bash
cd "2-Imputation"
Rscript impute_tpm.R melanoma 
Rscript impute_tpm.R head_neck
cd ../
```

Este paso utiliza el paquete ["scImpute"](https://github.com/Vivianstats/scImpute) para imputar los valores faltantes de aquellos genes que presentan una expresión génica nula en > 50% de las células estudiadas (_i.e. dropout rate_ > 50%). Nótese que los genes con un dropout del 100% no son imputados por falta de información.

## Normalizado y evaluación de distintos métodos de normalización

``` bash
cd "3-Normalization"
Rscript normalization.R melanoma
Rscript normalization.R head_neck
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
Rscript tsne_metabolic_genes.R melanoma
Rscript tsne_metabolic_genes.R head_neck

# UMAP (Opcional)
Rscript umap_metabolic_genes.R melanoma
Rscript umap_metabolic_genes.R head_neck

# Matriz de correlación de Spearman
Rscript inter_tumor_correlation.R melanoma
Rscript inter_tumor_correlation.R head_neck
cd ../
```

En este paso nos fijaremos en los 1566 genes metabólicos para agrupar y visualizar las células con t-SNE y UMAP (las visualizaciones pueden variar ligeramente respecto al gráfico del artículo científico debido a la inicialización aleatoria del t-SNE y UMAP). En concreto, estudiaremos por un lado las células malignas y por otro las células sanas. 

También calcularemos la matriz de correlación de Spearman (no paramétrico) de los distintos tumores dentro de una misma patología para mostrar su heterogeneidad. Los distintos melanomas presentes en nuestro dataset están clasificados como melanomas, pero sus perfiles de expresión de genes metabólicos demuestran que son distintos entre sí (ídem para las neoplasias de HNSCC).

## Heatmap y violinplot de la actividad metabólica en los distintos tipos celulares


``` bash
cd 5-PathwayActivity
Rscript scRNA_pathway_activity.R melanoma
Rscript scRNA_pathway_activity.R head_neck
Rscript TCGA_pathway_activity.R
cd ..
```

En este paso calculamos la actividad de las distintas rutas metabólicas en cada linaje celular, tanto en el dataset de bulk RNA-seq como en el de scRNA-seq.
Con el script `scRNA_pathway_activity.R` creamos un heatmap donde desglosamos en detalle la actividad de cada ruta metabólica de interés en nuestros tipos celulares y un violinplot donde se muestra la actividad global de las rutas en cada tipo celular. Así mismo, se genera una matriz con los p-valores para la actividad de cada ruta en cada tipo celular.

Para calcular la actividad de cada ruta, se obtuvo en cada linaje celular la expresión media de los genes que la constituyen, luego se dividieron estas medias por la actividad media del gen a lo largo de todos los linajes celulares para calcular la actividad relativa (a lo fold change), se penalizaron con pesos los genes que participan en más de una ruta y se sumó la actividad de todos los genes de dicha ruta. 

*The bulk RNA-seq data was downloaded from TCGA website, please see the instruction of data downloading and preprocessing in Data/TCGA/README.md* 

Metabolic pathway heterogeneity
-------------------------------
``` bash
cd 6-PathwayHeterogeneity
Rscript intra_malignant_heterogeneity.R melanoma
Rscript intra_malignant_heterogeneity.R head_neck
Rscript intra_non-malignant_heterogeneity.R melanoma
Rscript intra_non-malignant_heterogeneity.R head_neck
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R melanoma
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot.R head_neck
Rscript OXPHOS_Glycolysis_Hypoxia_Correlation_plot-CCLE.R
Rscirpt GeneSignature-of-Low_OXPHOS_Glycolysis_Hypoxia.R melanoma
Rscript GeneSignature-of-Low_OXPHOS_Glycolysis_Hypoxia.R head_neck
cd ..
```
In this step, the PCA and GSEA analysis will be performed to investigate the metabolic pathway heterogeneity across single cells in malignant and non-malignant cell populations. The scatter plots will be performed to compare activities of OXPHOS, glycolysis and response to hypoxia in single malignant cells and cultured cell lines from CCLE database. The gene signatures in single cells with low OXPHOS/glycolysis/hypoxia activity will be identified and stored as the text files, which can be used as the input of GO analysis on the website: http://metascape.org

Metabolic features of nonmalignant cell subtypes
-----------------------------------
``` bash
cd 7-MetabolicFeatures
Rscript non-malignant_subtype.R melanoma
Rscript non-malignant_subtype.R head_neck
cd ..
```
The metabolic features in different T-cell subtypes and fibroblast subtypes will be identified in this step. 

