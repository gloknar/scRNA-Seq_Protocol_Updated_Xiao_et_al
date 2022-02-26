##### REFORMATEAR EL HEATMAP PARA QUE QUEDE COMO EL PAPER ORIGINAL ###############
# Paquetes y funciones auxiliares
library(pheatmap)

# Opciones
options(stringsAsFactors = F)

data <- readRDS("data.rds")
sorted_rows <- readRDS("sorted_rows.rds")
sorted_columns <- readRDS("sorted_columns.rds")

color <- colorRampPalette(c("blue", "white", "red"))(100)

mybreaks <- c(
  seq(0, 0.5, length.out = 33),
  seq(0.51, 1.5, length.out = 33),
  seq(1.51, max(data), length.out = 34)
)

data[data >= 3] <- 3   # EL mñaximoera 5, pero si los capamos a 3, la leyenda se queda igual que la imagen del paper

pdf("heatmap2.pdf", onefile = T,
    width = 9, height = 9)   # He aumentado la width de 6 a 9 para que se queden rectangulos en vez de cuadrados

pheatmap(data[sorted_rows, sorted_columns], cluster_cols = F,
         cluster_rows = F, color = color, breaks = mybreaks)

dev.off()

