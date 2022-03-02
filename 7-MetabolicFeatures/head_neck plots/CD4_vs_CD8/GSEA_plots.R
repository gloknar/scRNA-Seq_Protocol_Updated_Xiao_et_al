archivo_CD4 <- list.files(pattern = "gsea_report_for_CD4_*")
archivo_CD8 <- list.files(pattern = "gsea_report_for_CD8_*")

tabla_CD4 <- read.delim(archivo_CD4, header = T, row.names = 1)
tabla_CD8 <- read.delim(archivo_CD8, header = T, row.names = 1)

df_CD4 <- tabla_CD4[1:10,c(1,5:6)]
df_CD4$identidad <- rep("CD4+", 10)
  
df_CD8 <- tabla_CD8[1:10,c(1,5:6)]
df_CD8$identidad <- rep("CD8+", 10)

df_total <- rbind(df_CD4, df_CD8)
df_total$identidad <- as.factor(df_total$identidad)

# p-valores
for (i in c(1:nrow(df_total))) {
  if (df_total[i,3] <= 0.05) {
    df_total[i,1] <- paste0(df_total[i,1],"(Significant)")   
  }  
}


library(ggplot2)
p <- ggplot(df_total, aes(x = reorder(GS.br..follow.link.to.MSigDB, -NES), y = NES, fill = identidad)) + geom_bar(stat = "identity") + 
  coord_flip() + ggtitle("head_neck") + ylim(c(-2,2)) + scale_fill_manual(values=c("#00c0c5", "#f8776d")) + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p

ggsave(filename = "GSEA_CD4_vs_CD8.png", plot = p)

