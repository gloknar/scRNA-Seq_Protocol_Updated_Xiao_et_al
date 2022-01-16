options(stringsAsFactors=FALSE)

pathway_file <- "../Data/KEGG_metabolism.gmt"
hall_gmt <- '../Data/h.all.v6.1.symbols.gmt'

args <- commandArgs()

print(args)
class(args)

# args[6] es -n
num_nucleos <- as.integer(args[7])

print(paste("Nº núcleos activos:", num_nucleos))


