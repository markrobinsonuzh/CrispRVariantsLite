################################################################################
# Check that the currently-installed version of R
# is at least the minimum required version.
################################################################################
R_min_version = "3.3.0"
R_version = paste0(R.Version()$major, ".", R.Version()$minor)
if(compareVersion(R_version, R_min_version) < 0){
  stop("You do not have the latest required version of R installed.\n",
       "Launch should fail.\n",
       "Go to http://cran.r-project.org/ and update your version of R.")
}
################################################################################
# Install basic required packages if not available/installed.
################################################################################

download_not_installed <- function(x, install.method = c("Bioconductor", "CRAN")){
    availpacks = .packages(all.available = TRUE)
    missingPackages = x[!(x %in% availpacks)]
    message("The following packages were missing. Installation attempted...")
    message(missingPackages)
  
    if(length(missingPackages) > 0){
        inst.method <- match.arg(install.method)
        inst.f <- ifelse(inst.method == "CRAN", install.packages, bioclite)
        if (inst.method  == "Bioconductor"){
            source("http://bioconductor.org/biocLite.R")
        }
        dummy <- lapply(missingPackages, function(pckg){
            message(sprintf("Installing package %s from %s \n", pckg, inst.method))
            inst.f(pckg)
        }) 
    }
}

cran_packages <- c("ggplot2", "shiny", "shinyBS", "shinydashboard",
                   "rhandsontable", "htmltools")
bioc_packages <- c("CrispRVariants", "GenomicRanges", "GenomicFeatures",
                   "AnnotationDbi", "Rsamtools", "Biostrings", "IRanges")
download_not_installed(cran_packages, "CRAN")
download_not_installed(bioc_packages, "Bioconductor")
