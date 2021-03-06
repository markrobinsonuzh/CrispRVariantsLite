## Introduction

This github repository contains the code underlying the *CrispRVariantsLite* Shiny-based web app, which accompanies the main [CrispRVariants](http://www.bioconductor.org/packages/CrispRVariants.html) package.  This web app is primarily available from the [Robinson lab](http://www.imls.uzh.ch/research/robinson.html) web server at [University of Zurich](http://www.uzh.ch/de.html) from the following link: [CrispRVariantsLite](http://imlspenticton.uzh.ch:3838/CrispRVariantsLite).  Instructions are given in the start up page of the app.

## Example data

Example ZIP files (and corresponding TXT files with description of genome location and/or guide+PAM sequence) for the different entry points can be found in the [/examples/](https://github.com/markrobinsonuzh/CrispRVariantsLite/tree/master/examples) directory.  If retrieving these files from the github repo, be sure to click on 'Raw' or 'View Raw'.

## Running this app locally

The CrispRVariantsLite web app can be run locally, after installation of the necessary R packages (see below) and well as the command line tools: [samtools](http://www.htslib.org/), [bwa mem](http://bio-bwa.sourceforge.net/) and the indices for the organism of interest. After all this is completed, clone a version of this repository, set up the symbolic links (see the [makelinks](https://github.com/markrobinsonuzh/CrispRVariantsLite/blob/master/makelinks) file) in the data/genome and data/txdb directories of the cloned repository and then the app can be run from the R console as:

    #install.packages("shiny")  # only needed if 'shiny' package is not already installed
    shiny::runApp( "/dir/to/CrispRVariantsLite/local/clone")

## Installing all the R packages

The following packages will be needed by the web app:

    source("https://bioconductor.org/biocLite.R")
    library(BiocInstaller)
    biocLite( c("CrispRVariants", "GenomicFeatures",
                "AnnotationDbi", "GenomicRanges",
                "IRanges", "Rsamtools", "Biostrings") )

    install.packages("ggplot2",quiet=TRUE)
    install.packages("shiny",quiet=TRUE)
    install.packages("shinydashboard",quiet=TRUE)
    install.packages("shinyBS",quiet=TRUE)
    install.packages("rhandsontable",quiet=TRUE)


## Other important details

The web application makes calls to bwa mem and samtools.  Therefore, for every organism of interest, a bwa index will need to be created from the genome sequence.  These need to be copied (or preferably, symbolically linked) to the data/genome directory; similarly, TxDb objects need to be saved to the data/txdb directory of the app.  The convention used is that all the files ending in .fa in the data/genome are presented to the user in the dropdown menu of the web application.  Importantly, the TxDB objects in the data/txdb directory must be named with the same prefix but a .sqlite extension.  That is, if you are working with the hg19 genome, the application is expecting hg19.fa (and hg19.fa.XXX for all the bwa indices) in the data/genome directory and hg19.sqlite in the data/txdb directory.
