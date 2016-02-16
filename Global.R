# load packages
library("shiny"); packageVersion("shiny")
library("shinydashboard"); packageVersion("shinydashboard")
library("shinyBS"); packageVersion("shinyBS")
library("ggplot2"); packageVersion("ggplot2")
library("CrispRVariants"); packageVersion("CrispRVariants")
library("Rsamtools"); packageVersion("Rsamtools")
library("rhandsontable") ; packageVersion("rhandsontable")
library("ShortRead"); packageVersion("ShortRead")
library("rtracklayer"); packageVersion("rtracklayer")
library("gdata"); packageVersion("gdata")
library("gridExtra"); packageVersion("gridExtra")
library("GenomicFeatures"); packageVersion("GenomicFeatures")
library("AnnotationDbi"); packageVersion("AnnotationDbi")
library("GenomicRanges"); packageVersion("GenomicRanges")
library("IRanges"); packageVersion("IRanges")


#Utilities
source("./core/setting.R", local = TRUE)

# For pasting times into things
simpletime = function(){gsub("\\D", "_", Sys.time())}


# MHmakeRandomString(n, length)
# function generates a random string random string of the
# length (length), made up of numbers, small and capital letters

MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
      lenght, replace=TRUE),
      collapse = "")
  }
  return(randomString)
}

# remove all the files in a directory
empty_dir <- function(dirname){
  ## Get vector of all file names
  ff <- dir(create_dir(dirname), recursive=TRUE, full.names=TRUE)
  ## Remove all files
  unlink(ff, recursive=TRUE, force=FALSE)
}

#create empty directory if not created
create_dir <- function(dirname){
  #Check existence of directory and create if doesn't exist
  #http://stackoverflow.com/questions/4216753/check-existence-of-directory-and-create-if-doesnt-exist
  ifelse(!dir.exists(file.path(dirname)), dir.create(file.path(dirname)), FALSE)
  outDir <- file.path(mainDir, subDir)
  return(outDir)
}


