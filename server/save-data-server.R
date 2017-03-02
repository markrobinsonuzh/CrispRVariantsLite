# Modal for selecting input type
output$bams <- renderUI({
    
     if(!is.null(v$bm_fnames))
     {
         tags$div(
             p("BAM files stored")
         )
     }
     else
     {
         tags$div(
             bsButton("select_AB1", 'Upload ZIP of AB1s',
                      style = "primary", block = TRUE),
             bsButton("select_FastQ", 'Upload ZIP of FASTQs',
                      style = "primary", block = TRUE),
             p(),
             fileInput('upload_bams', 'Upload ZIP of BAMs',
                       multiple = F, width = "100%")
         )
     }
})

#___________________________________________________________________________________
# TO DO? - Make data download options depend on input data?
# Otherwise not necessary to have renderUI

# Modal for choosing data to save.
# Available output data will depend on input format
output$ui_options <- renderUI({
     base_choices <- c("save_plot", "save_fastq", "save_bam", "save_metadata",
                       "save_vc", "save_allele_seqs")
     names(base_choices) <- c("Allele plot", "FASTQ files", "Mapped BAM files",
                          "Metadata table", "Variant counts table",
                          "Consensus allele sequences")
     
     tags$div(style = "font-size: 16px;",
     checkboxGroupInput("data_to_save", label = NULL, choices = base_choices),
     downloadButton('download_zip', 'Download'))
})



output$download_zip <- downloadHandler(
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = sprintf("CrispRVariants_data_%s.zip", simpletime()),
  
    # the argument 'file' is a nonexistent temporary file name
    content = function(file){
        fns <- list("save_bam" = get_bam_fnames, "save_fastq" = get_fastq_fnames,
                    "save_metadata" = get_table_fname, "save_plot" = get_plot,
                    "save_vc" = get_vc, "save_allele_seqs" = get_alleles)
        cat(sprintf("functions to run %s \n", input$data_to_save))
        fs <- lapply(fns[input$data_to_save], function(f) f())
        cat(sprintf("filenames for zip %s \n", fs))      
        zip(zipfile = file, files = unlist(fs), flags = "-j")
    
    },
     contentType = "application/zip"
)  


#___________________________________________________________________________________
# These functions return names of files to be added to the zip.
# Files are created if necessary
get_bam_fnames <- function(){
    cat(sprintf("bam fnames %s \n", v$bam_fnames))
    result <- dir(v$bam_dir, pattern = ".bam$", full.names = TRUE, recursive = TRUE)
    cat(sprintf("result %s \n", result))
    result
}


get_fastq_fnames <- function(){
    cat(sprintf("fq fnames %s \n", v$fq_fnames))
    result <- dir(v$fq_dir, pattern = ".fastq$", full.names = TRUE, recursive = TRUE)
    cat(sprintf("result %s \n", result))
    result
}


get_table_fname <- function(){
    md.fname <- file.path(v$temp.dir, "metadata.csv")
    sprintf("md fname = %s \n", md.fname)
    write.table(t$DF, md.fname, sep = ",", row.names = FALSE, quote = FALSE)
    md.fname
}


get_plot <- function(){
    filename <- file.path(v$temp.dir, "CrispRVariants_alleles.pdf")
    
    pixelratio <- session$clientData$pixelratio
    cdata <- session$clientData
    br_wd <- cdata[["output_crispplots_width"]]

    br_ht <- 600 # MAKE THESE GLOBAL VARIABLES
    res <- 72

    ht <- (br_ht*pixelratio)/res
    wd <- (br_wd*pixelratio)/res
              
    pdf(filename, height = ht, width = wd, useDingbats = FALSE)
              
    r_ht <- as.numeric(strsplit(input$row.ht.ratio, ":")[[1]])
    c_wd <- as.numeric(strsplit(input$col.wdth.ratio, ":")[[1]])

    CrispRVariants:::arrangePlots(
      annotation_plot(), allele_plot(), frequency_heatmap(),
      col.wdth.ratio = c_wd, row.ht.ratio = r_ht,
      left.plot.margin = ggplot2::unit(c(1,0,5,2), "lines")) 
  
    dev.off()
    filename
}


get_vc <- function(){
    filename <- file.path(v$temp.dir, "CrispRVariants_allele_counts.csv")
    write.csv(variantCounts(d$cset), file = filename, row.names = FALSE)
    filename
}


get_alleles <- function(){
    filename <- file.path(v$temp.dir, "CrispRVariants_consensus_allele_seqs.fasta")
    sqs <- consensusSeqs(d$cset)
    cat(paste(sprintf(">%s\n%s\n", names(sqs), sqs), collapse = ""), file = filename)
    filename
}
