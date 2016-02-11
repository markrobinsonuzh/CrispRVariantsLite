# load packages
library("shiny"); packageVersion("shiny")
#library("shinyAce"); packageVersion("shinyAce")
#library("shinyFiles"); packageVersion("shinyFiles")
library("shinydashboard"); packageVersion("shinydashboard")
#library("shinyjs"); packageVersion("shinyjs")
library("shinyBS"); packageVersion("shinyBS")
#library("shinythemes"); packageVersion("shinythemes")
library("ggplot2"); packageVersion("ggplot2")
library("CrispRVariants"); packageVersion("CrispRVariants")
library("Rsamtools"); packageVersion("Rsamtools")
library("rhandsontable") ; packageVersion("rhandsontable")
#library("Cairo"); packageVersion("Cairo")  # For nicer ggplot2 output when deployed on Linux
library("ShortRead"); packageVersion("ShortRead")
library("rtracklayer"); packageVersion("rtracklayer")
library("gdata"); packageVersion("gdata")
library("gridExtra"); packageVersion("gridExtra")
library("GenomicFeatures"); packageVersion("GenomicFeatures")
library("AnnotationDbi"); packageVersion("AnnotationDbi")
library("GenomicRanges"); packageVersion("GenomicRanges")
library("IRanges"); packageVersion("IRanges")

# Default options for app startup
#TODO: source("core/default-parameters.R", local = TRUE)

# Graphic-saving utilities
#source("core/ggsave.R", local = TRUE)

#Utilities
source("server/setting.R", local = TRUE)

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
      collapse="")
  }
  return(randomString)
}


load_files_from_dir <- function(dirname = "raw_data", uploaded = F) {
  inputDir <- create_dir(dirname)
  filenames <- list.files(inputDir, pattern = "*.fcs$")
  f <- c(f, filenames)
  return(f)
}

remove_files <- function(filename) {
  unlink(filename, recursive = T, force = F)
  #file.remove(filename)
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


# ui for selecting the channels
uichannel <- function(inputId, label, s = 10, ... ){
  selectInput(inputId, label,
    choices = list(),
    width = "100%",
    size = s, ... )
}


# ui for selecting user files (fcs files)
uifiles <- function(inputId, ... ){
  selectInput(inputId, label = "File/Sample",
    choices = load_files_from_dir(),
    selected = list(),  ...
  )
}

# ui for selecting user files (fcs files)
uifiles.norm <- function(inputId, dirname, ... ){
  selectInput(inputId, label = "File/Sample",
    choices = load_files_from_dir(dirname),
    selected = list(),  ...
  )
}

# ui for selecting user files (fcs files)
uifiles.xml <- function(inputId, ... ){
  selectInput(inputId, label = "Gate/Xml.file",
    choices = load_files_from_dir.xml(),
    selected = list(),  ...
  )
}


# ui for selecting data created by user (fcs files)
uidata <- function(inputId, label, folder, pattern,  ... ){
  selectInput(inputId, label,
    choices = list.files(paste0(folder), pattern),
    selected = list(),  ...
  )
}

#ui plot for creating plot function
uiplots <- function(i){
  id = paste0("nplot",i, sep="")
  plotOutput(outputId =  id,
    dblclick = paste0("dbclik", i),
    brush = brushOpts( id = paste0("plot_brush",i),
      resetOnNew = TRUE )
  )
}

