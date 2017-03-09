# load packages
cat("library imports started", file = "crv_scratch.txt", append = TRUE)
cat(system("date", intern = TRUE), file = "crv_scratch.txt", append = TRUE)
library("shiny"); packageVersion("shiny")
library("shinydashboard"); packageVersion("shinydashboard")
library("shinyBS"); packageVersion("shinyBS")
library("ggplot2"); packageVersion("ggplot2")
library("CrispRVariants"); packageVersion("CrispRVariants")
#library("Rsamtools"); packageVersion("Rsamtools")
library("rhandsontable") ; packageVersion("rhandsontable")
#library("ShortRead"); packageVersion("ShortRead")
#library("rtracklayer"); packageVersion("rtracklayer")
#library("gdata"); packageVersion("gdata")
#library("gridExtra"); packageVersion("gridExtra")
#library("GenomicFeatures"); packageVersion("GenomicFeatures")
#library("AnnotationDbi"); packageVersion("AnnotationDbi")
#library("GenomicRanges"); packageVersion("GenomicRanges")
#library("IRanges"); packageVersion("IRanges")
#library("seqinr");packageVersion("seqinr")
cat("library imports finished", file = "crv_scratch.txt", append = TRUE)
cat(system("date", intern = TRUE), file = "crv_scratch.txt", append = TRUE)


source("./core/setting.R", local = TRUE)
source("./core/bsModal.R", local = TRUE)
#source("./core/misc.R", local = TRUE)

#_________________________________________________________________

# Utility functions

# For pasting times
simpletime = function(){gsub("\\D", "_", Sys.time())}

# Check errors of different types
is.error <- function(x) inherits(x, "try-error")


# Generate a random string.  Code from
# Bioconductor package Chimera
# https://bioconductor.org/packages/release/bioc/html/chimera.html
MHmakeRandomString <- function(n=1, length=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
      length, replace=TRUE),
      collapse = "")
  }
  return(randomString)
}
