## Introduction

This github repository contains the code underlying the *CrispRVariantsLite* Shiny-based web app, which accompanies the [CrispRVariants](http://www.bioconductor.org/packages/CrispRVariants.html) package.  This web app is primarily available from the [Robinson lab](http://www.imls.uzh.ch/research/robinson.html) web server at [University of Zurich](http://www.uzh.ch/de.html) from the following link: [CrispRVariantsLite](http://imlspenticton.uzh.ch:3838/CrispRVariantsLite).  Instructions are given in the start up page of the app.

## Running this app locally

The CrispRVariantsLite web app can be run locally, after installation of the necessary R packages (see below) and well as the command line tools: [samtools](http://www.htslib.org/), [bwa mem](http://bio-bwa.sourceforge.net/) and the indices for the organism of interest. After all this is completed, clone a version of this repository, set up the symbolic links (see the [makelinks](https://github.com/markrobinsonuzh/CrispRVariantsLite/blob/master/makelinks) file) in the data/genome and data/txdb directories of the cloned repository and then the app can be run from the R console as:

    #install.packages("shiny")  # only needed if 'shiny' package is not already installed
    shiny::runApp( "/dir/to/CrispRVariantsLite/local/clone")

## Installing all the R packages


