#!/usr/bin/env bash

# Si te da errores al ejecutar este archivo, prueba a hacer un backup del
# mismo y quitarle el retorno de carro (\r) con sed (o con dos2unix):
#
# sed -i 's/\r$//' download_dataset.sh


###### Descargamos el dataset del melanoma #####

url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz'  # Recuerda que las single quotes en bash no interpretan metacarácteres, preservan el valor literal de la string
wget $url || curl $url -o GSE72056_melanoma_single_cell_revised_v2.txt.gz     # Prueba a descargar el dataset mediante wget, y si falla, prueba con curl ( || significa OR en bash)
gunzip GSE72056_melanoma_single_cell_revised_v2.txt.gz  # Descomprime el archivo descargado, equivalente a gzip -d <nombre_comprimido>



###### Corregimos nombres de genes incorrectos (causado por una conversión incorrecta desde Excel) #####
# 1-Mar -> MARCH1
# 1-Dec -> DEC1
# 1-Sep -> SEPT1
# SEPT15 should be SEP15

command="sed"
for (( i = 15; i > 0; i-- )); do
	command="$command -e \"s/$i-Mar/MARCH$i/g\" -e \"s/$i-Sep/SEPT$i/g\" -e \"s/$i-Dec/DEC$i/g\" -e \"s/SEPT15/SEP15/g\""
done

eval "cat GSE72056_melanoma_single_cell_revised_v2.txt | $command > GSE72056_melanoma_single_cell_corrected.txt" # Ejecuta el comando sed tal como lo configuramos arriba y guarda el output en el archivo GSE72056_melanoma_single_cell_corrected.txt
rm GSE72056_melanoma_single_cell_revised_v2.txt       # Con el dataset ya corregido, eliminamos el dataset descargado (sin corregir)
# rm GSE72056_melanoma_single_cell_revised_v2.txt.gz  # El archivo comprimido se borra automáticamente al descomprimirlo con gzip/gunzip, por lo que esta línea es innecesaria



###### Descargamos el dataset del cancer de cabeza y cuello #####
url='ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz'
wget $url || curl $url -o GSE103322_HNSCC_all_data.txt.gz
gunzip GSE103322_HNSCC_all_data.txt.gz



###### Movemos los datasets a su carpeta #####
if [ ! -d "datasets" ]  # Si no existe la carpeta "datasets", la creamos con mkdir. Ten en cuenta que debes dejar un espacio entre cada cosa dentro de los corchetes, si pones if [-d "tal"], no funciona
then
  mkdir datasets
fi

mv {GSE103322_HNSCC_all_data.txt,GSE72056_melanoma_single_cell_corrected.txt} datasets/ # equivalente a mv GSE103322_HNSCC_all_data.txt dataset/ & mv GSE72056_melanoma_single_cell_corrected.txt dataset/
