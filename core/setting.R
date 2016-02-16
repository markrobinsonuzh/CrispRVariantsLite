# genomic database
# danRer <- "/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7"

genome <- "/data/genome"

gendb <- dir(genome,pattern = ".fa$", recursive = TRUE, full.names = TRUE)
gendb.gz <- dir(genome,pattern = "fa.gz$", recursive = TRUE, full.names = TRUE)
# name of sequence
# genlist <- gsub(".fa","",basename(gendb))
genlist <- basename(gendb)

genlist.gz <- basename(gendb.gz)

# transcript database library
txDB.dir <- dir("data/txdb/", recursive = TRUE, full.names = TRUE)
txDBl <- basename(txDB.dir)