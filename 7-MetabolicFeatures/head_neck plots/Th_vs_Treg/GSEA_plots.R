archivo_Th <- list.files(pattern = "gsea_report_for_Th_*")
archivo_Treg <- list.files(pattern = "gsea_report_for_Treg_*")

tabla_Th <- read.delim(archivo_Th, header = T, row.names = 1)
tabla_Treg <- read.delim(archivo_Treg, header = T, row.names = 1)

df_Th <- tabla_Th[1:10,c(1,5:6)]
df_Th$identidad <- rep("Th", 10)
  
df_Treg <- tabla_Treg[1:10,c(1,5:6)]
df_Treg$identidad <- rep("Treg", 10)

df_total <- rbind(df_Th, df_Treg)
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

ggsave(filename = "GSEA_Th_vs_Treg.png", plot = p)

