# genomic database
genome <- "./data/genome"

# List available genomes, ending in .fa
gendb <- dir(genome, ".fa$", recursive = TRUE, full.names = TRUE)

# Genome names are the basename of the gendb files
genlist <- basename(gendb)
# Remove extension for choosing genome
genlist_base <- tools::file_path_sans_ext(genlist)
names(genlist) <- genlist_base

replace_first <- function(x, first_element){
    return(c(first_element, which(x != first_element)))
}


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
