# genomic database
# danRer <- "/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7"
# genome <- "/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7"
genome <- "./data/genome"

gendb <- dir(genome, recursive = TRUE, full.names = TRUE)

# name of sequence
genlist <- basename(gendb)
genlist.gz <- basename(gendb)

# transcript database library
txDB.dir <- dir("./data/txdb/", recursive = TRUE, full.names = TRUE)
txDBl <- basename(txDB.dir)
