# Proyecto scRNA-seq de ovarios

Xiao _et al._ publicaron en 2019 un [artículo](https://www.nature.com/articles/s41467-019-11738-0) en Nature Communications en el que caracterizan 2 tipos de tumores a nivel de transcriptoma de célula única y GSEA (Gene Set Enrichment Analysis). El producto de dicha investigación está disponible en forma de [repositorio](https://github.com/LocasaleLab/Single-Cell-Metabolic-Landscape) con todo el código para `R 3.5` + `bash`

Nosotros queremos replicar dicho protocolo en nuestras muestras del Karolinska Institutet, las cuales son bulk RNA-seq y scRNA-seq de ovarios. Para ello he creado el presente repositorio, en el cual se traduce el protocolo de Xiao y compañeros a la versión más actualizada de R a fecha de redacción de este documento, `R 4.1`.


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

Con el script de bash descargamos los datasets desde el GEO, y con sendos scripts de R corregimos erratas presentes en dichos datasets, filtramos las células, creamos 2 objetos de tipo `Single Cell Experiment` y los guardamos en los archivos `<head_neck/melanoma>/filtered_sce.rds` 

## Imputación de NAs

``` bash
cd "2-Imputation"
Rscript impute_tpm.R melanoma 
Rscript impute_tpm.R head_neck
cd ../
```

Este paso utiliza el paquete ["scImpute"](https://github.com/Vivianstats/scImpute) para imputar los valores faltantes de aquellos genes que presentan una expresión génica nula en > 50% de las células estudiadas (_i.e. dropout rate_ > 50%).

## Normalizado y evaluación de distintos métodos de normalización

``` bash
cd "3-Normalization"
Rscript normalization.R melanoma
Rscript normalization.R head_neck
cd ../
```
Four commonly used data normalization methods are applied on each dataset. The distribution of relative gene expression of each cell type will be ploted to evaluate and select the best normalization method.

## Landscape of the metabolic gene expression profile

``` bash
cd "4-Clustering"
Rscript metabolic_landscape.R melanoma
Rscript metabolic_landscape.R head_neck
Rscript inter_tumor_distance.R melanoma
Rscript inter_tumor_distance.R head_neck
cd ../
```

En este paso se emplea el algoritmo t-SNE para visualizar la expresión de genes metabólicos en millones de células (El resultado puede variar ligeramente respecto al gráfico del artículo científico debido a la inicialización aleatoria del algoritmo). También se genera la matriz de correlación de spearman para mostrar la heterogeneidad inter-tumoral usando dichos genes metabólicos.

## Actividad de las rutas metabólicas en distintos tipos celulares


``` bash
cd 5-PathwayActivity
Rscript scRNA_pathway_activity.R melanoma
Rscript scRNA_pathway_activity.R head_neck
Rscript TCGA_pathway_activity.R
cd ..
```

This step will calculate the metabolic pathway activities for different single cell populations or bulk tumor/normal samples. The scatter plot will show the discrepancy of pathway activities between single malignant cells and bulk tumors. The violin plot will show the distribution of metabolic pathway activities in single cell populations or bulk tumor/normal samples.

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

