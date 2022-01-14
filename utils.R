##### Funciones relacionadas con rutas metabólicas #####

# Función para leer las rutas metabólicas contenidas en un archivo .gmt y
# devolverlas en forma de lista, obtenida de fgsea versión 1.21.0
# https://github.com/ctlab/fgsea/blob/master/R/pathways.R
gmtPathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
}


# Función casera para calcular el nº de rutas metabólicas en las que participa
# un gen
num_of_pathways <- function (gmtfile, overlapgenes){
  # Computes the sample covariance between two vectors.
  #
  # Argumentos:
  #   gmtfile: Ruta (relativa o absoluta) al archivo .gmt con las rutas metabólicas
  #   overlapgenes: vector de tipo "character" con los nombres de los genes a buscar en el archivo .gmt. 
  #       Deben ser genes presentes tanto en el objeto sce de interés como en el archivo gmtfile
  #
  # Output:
  #   Un dataframe de dimensiones overlapgenes X 1 que contiene el nº de rutas en las que participa
  #      cada gen
  pathways <- gmtPathways(gmtfile)
  pathway_names <- names(pathways)
  filter_pathways <- list()
  for (p in pathway_names){
    genes <- pathways[[p]]
    common_genes <- intersect(genes,overlapgenes)
    if(length(common_genes>=5)){
      filter_pathways[[p]] <- common_genes
    }
  }
  
  all_genes <- unique(as.vector(unlist(filter_pathways)))
  gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
  for(p in pathway_names){
    for(g in filter_pathways[[p]]){
      gene_times[g,"num"] = gene_times[g,"num"]+1
    }
  }
  gene_times
} 



# runGSEA_preRank.R
runGSEA_preRank<-function(preRank.matrix,gmt.file,outname){
  # descending numerical order
  # dump preRank into a tab-delimited txt file
  write.table(x = preRank.matrix,
              file = "prerank.rnk",
              quote = F,
              sep = "\t",
              col.names = F,
              row.names = T)
  
  # call java gsea version
  command <- paste('java -Xmx512m -cp ../gsea-3.0.jar xtools.gsea.GseaPreranked -gmx ', gmt.file, ' -norm meandiv -nperm 1000 -rnk prerank.rnk ',
                   ' -scoring_scheme weighted -make_sets true -rnd_seed 123456 -set_max 500 -set_min 15 -zip_report false ',
                   ' -out preRankResults -create_svgs true -gui false -rpt_label ',outname, sep='')
  
  if(get_os() == "win"){
    system(command, show.output.on.console = F)
  } else {
    system(command)
  }
  unlink(c('prerank.txt'))
}




##### Funciones de normalizado #####

# Función de normalizado de DESeq2 versión 1.35.0, obtenida de:
# https://github.com/mikelove/DESeq2/blob/master/R/core.R
estimateSizeFactorsForMatrix <- function(counts, locfunc = stats::median,
                                         geoMeans, controlGenes,
                                         type = c("ratio", "poscounts")) {
  type <- match.arg(type, c("ratio","poscounts"))
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    if (type == "ratio") {
      loggeomeans <- rowMeans(log(counts))
    } else if (type == "poscounts") {
      lc <- log(counts)
      lc[!is.finite(lc)] <- 0
      loggeomeans <- rowMeans(lc)
      allZero <- rowSums(counts) == 0
      loggeomeans[allZero] <- -Inf
    }
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}




##### Funciones misceláneas #####

# Función de rappdirs versión 0.3.3 para detectar el sistema operativo empleado,
# obtenida de: https://github.com/r-lib/rappdirs/blob/master/R/utils.r#L1
get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}

