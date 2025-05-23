##### REFORMATEAR EL HEATMAP PARA QUE QUEDE COMO EL PAPER ORIGINAL ###############
# Paquetes y funciones auxiliares
library(pheatmap)

# Opciones
options(stringsAsFactors = F)

data <- readRDS("data.rds")
sorted_rows <- readRDS("sorted_rows.rds")
sorted_columns <- readRDS("sorted_columns.rds")

color <- colorRampPalette(c("#0095ff", "#ffffff", "#ff6663"))(100)

# El máximo era 5, pero si los capamos a 3, la leyenda se queda igual que la imagen del paper. DEBES HACER ESTO ANTES DE LLAMAR A MYBREAKS
data[data >= 3] <- 3   

mybreaks <- c(
  seq(0, 0.5, length.out = 33),
  seq(0.51, 1.5, length.out = 33),
  seq(1.51, max(data), length.out = 34)
)


pdf("head_neck_V3.pdf", onefile = T,
    width = 9, height = 9)   # He aumentado la width de 6 a 9 para que se queden rectangulos en vez de cuadrados

pheatmap(data[sorted_rows, sorted_columns], cluster_cols = F,
         cluster_rows = F, color = color, breaks = mybreaks)

dev.off()

################


sorted_rows2 <-  read.csv("sorted_rows_prueba.csv")
sorted_rows2 <- unname(unlist(sorted_rows2))  # Para pasar de algo a vector...

sorted_rows2 <- sorted_rows2[-16] # Quitamos el Steroid hormone biosynthesis 
data[sorted_rows2[1:70],]
sorted_rows2[16]
# View(data)



# pheatmap(data[sorted_rows2,sorted_columns], cluster_cols = F, cluster_rows = F, color = color, breaks = mybreaks)



pdf("head_neck_reordenado.pdf", onefile = T,
    width = 9, height = 9)   # He aumentado la width de 6 a 9 para que se queden rectangulos en vez de cuadrados

pheatmap(data[sorted_rows2,sorted_columns], cluster_cols = F, cluster_rows = F, color = color, breaks = mybreaks)

dev.off()
