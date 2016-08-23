## Introduction

This github repository contains the *CrispRVariantsLite* Shiny-based web app to accompany the [CrispRVariants](http://www.bioconductor.org/packages/CrispRVariants.html) package.  This web app is primarily available from the [Robinson lab](http://www.imls.uzh.ch/research/robinson.html) web server at [University of Zurich](http://www.uzh.ch/de.html) from the following link: [CrispRVariantsLite](http://imlspenticton.uzh.ch:3838/CrispRVariantsLite).  Instructions are given in the start up page of the app.

## Running this app locally

The CrispRVariantsLite web app can be run locally, after installation of the necessary R packages (see below) and well as the command line tools and [bwa mem](http://bio-bwa.sourceforge.net/) and the indices for the organism of interest.  After all this is completed, the app can be run from the R console as:

    install.packages("shiny")
    shiny::runGitHub( "CrispRVariantsLite", "markrobinsonuzh") 

## Installing all the R packages, command line tools and indices

TODO
