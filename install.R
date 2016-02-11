################################################################################
# Check that the currently-installed version of R
# is at least the minimum required version.
################################################################################
R_min_version = "3.1.0"
R_version = paste0(R.Version()$major, ".", R.Version()$minor)
if(compareVersion(R_version, R_min_version) < 0){
  stop("You do not have the latest required version of R installed.\n",
       "Launch should fail.\n",
       "Go to http://cran.r-project.org/ and update your version of R.")
}
################################################################################
# Install basic required packages if not available/installed.
################################################################################
download_not_installed = function(x){
  availpacks = .packages(all.available = TRUE)
  source("http://bioconductor.org/biocLite.R")
  missingPackages = x[!(x %in% availpacks)]
  message("The following packages were missing. Installation attempted...")
  message(missingPackages)
  if(length(missingPackages) > 0){
    for(i in missingPackages){
      message("Installing", i, "package using biocLite... \n")
      biocLite(i)
    }
  }
}
vanilla_install_pkgs = c("ggplot2", "Rcpp", "metricsgraphics","RColorBrewer", "dplyr", "flowCore", "tidyr","Rtsne")


download_not_installed_cran = function(x){
  availpacks = .packages(all.available = TRUE)
  missingPackages = x[!(x %in% availpacks)]
  message("The following packages were missing. Installation attempted...")
  message(missingPackages)
  if(length(missingPackages) > 0){
    for(i in missingPackages){
      message("Installing", i, "package using biocLite... \n")
      install.packages(i)  }
  }
}

mango_install_pkgs = c("shiny","shinyFiles", "Cairo", "shinyBS", "shinydashboard", "shinytheme") 

download_not_installed(vanilla_install_pkgs)
download_not_installed_cran(mango_install_pkgs)
################################################################################
# Should use latest GitHub version of shiny
################################################################################
shiny_okay = FALSE
if("shiny" %in% .packages(all.available = TRUE)){
  shiny_min_version = "0.11"
  shiny_compare = compareVersion(as.character(packageVersion("shiny")), shiny_min_version)
  if( shiny_compare >= 0 ){
    shiny_okay <- TRUE
  }
}
if(!shiny_okay){
  install.packages("devtools")
  devtools::install_github("rstudio/shiny")
}
################################################################################
# Should use latest GitHub version of rmarkdown
# https://github.com/rstudio/rmarkdown
################################################################################
rmarkdown_okay = FALSE
if("rmarkdown" %in% .packages(all.available = TRUE)){
  rmarkdown_min_version = "0.5"
  rmarkdown_compare = compareVersion(as.character(packageVersion("rmarkdown")), rmarkdown_min_version)
  if( rmarkdown_compare >= 0 ){
    rmarkdown_okay <- TRUE
  }
}
if(!rmarkdown_okay){
  install.packages("devtools")
  devtools::install_github("rstudio/rmarkdown")
}

################################################################################
# Should use latest GitHub version of shinythemes
# https://github.com/rstudio/shinythemes
################################################################################
rshinythemes_okay = FALSE
if("shinythemes" %in% .packages(all.available = TRUE)){
  rshinythemes_min_version = "0.5"
  rshinythemes_compare = compareVersion(as.character(packageVersion("shinythemes")), rshinythemes_min_version)
  if( rshinythemescompare >= 0 ){
    rshinythemes_okay <- TRUE
  }
}
if(!rshinythemes_okay){
  install.packages("devtools")
  devtools::install_github("rstudio/shinythemes")
}

################################################################################
# flowCore existence/version test, and installation
################################################################################
flowCore_okay = FALSE
if("flowCore" %in% .packages(all.available = TRUE)){
  flowCore_min_version = "1.20.11"
  flowCore_compare = compareVersion(as.character(packageVersion("flowCore")), flowCore_min_version)
  if( flowCore_compare >= 0 ){
    flowCore_okay <- TRUE
  }
}
if(!flowCore_okay) {
  # Go through recommended phyloseq installation steps
  # (1) Load biocLite
  source("http://bioconductor.org/biocLite.R")
  # (2) Install latest devel version from BioC
  useDevel(devel = TRUE)
  biocLite("flowCore", suppressUpdates = TRUE)
  # (3) Restore biocLite to release status
  useDevel(devel = FALSE)
}
################################################################################
# Load packages that must be fully-loaded
################################################################################
shiny_phyloseq_full_load_packages = c("shiny", "flowCore", vanilla_install_pkgs)
for(i in shiny_phyloseq_full_load_packages){
  library(i, character.only = TRUE)
  message(packageVersion(i))
}
################################################################################
