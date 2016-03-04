# genomic database
# danRer <- "/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7"
# genome <- "/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7"
genome <- "./data/genome"

gendb <- dir(genome, ".fa$", recursive = TRUE, full.names = TRUE)

# name of sequence
genlist.gz <- genlist <- basename(gendb)

# transcript database library
# this is inferred from genome name now
#txDB.dir <- dir("./data/txdb/", recursive = TRUE, full.names = TRUE)
#txDBl <- basename(txDB.dir)

################################################################################
# DEFAULT SETTING GUIDE
################################################################################

#d.start <- "23648474"
#d.length <- "5"
#d.chr <- "chr17"
#d.strand <-  "-"
g.seq <- "GGAGATCGCCACAAGATGTGAGG"

d.start <- NULL
d.length <- NULL
d.chr <- NULL
d.strand <-  "-"
#g.seq <- NULL

