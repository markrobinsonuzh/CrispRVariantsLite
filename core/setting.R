# genomic database
genome <- "./data/genome"

# List available genomes, ending in .fa
genlist <- dir(genome, ".fa$", recursive = TRUE)

# Remove extension for choosing genome
genlist_base <- tools::file_path_sans_ext(genlist)
names(genlist) <- genlist_base

# The transcript database is inferred from the genome name, these must match


################################################################################
# DEFAULT SETTING GUIDE
################################################################################

#d.start <- "23648474"
#d.length <- "5"
#d.chr <- "chr17"
#d.strand <-  "-"
#g.seq <- "GGAGATCGCCACAAGATGTGAGG"

#d.start <- NULL
#d.length <- NULL
#d.chr <- NULL
#d.strand <-  "-"
#g.seq <- NULL



# FOLDER NAMES AREN'T KEPT ANYMORE?
# BOX ON ALIGNMENT PLOT NOT POSITIONED CORRECTLY AB1?