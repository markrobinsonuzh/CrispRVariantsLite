# ---------------
# Convert AB1 to FASTQ
# ---------------
convertAb1toFastq <- function(){

    temp <- unzip(input$ab1_input$datapath, exdir = v$ab1_dir)
    v$ab1_fnames <- dir(v$ab1_dir, "ab1$", recursive = TRUE)
   
    # Warn if no files found
    if (length(v$ab1_fnames) == 0){
        createAlert(session, "alertAB1", "prepAlertAB1",
        title = "WARNING", style = "warning", append = FALSE,
        content = "No files ending in '.ab1' found")
        return()
    }

    # get the sequence names
    v$sq_nms <- gsub(".ab1","", basename(v$ab1_fnames))
    # replace spaces and slashes in filename with underscores
    temp <- iconv(v$sq_nms, to = "ASCII//TRANSLIT")

    v$fq_fnames <- paste0(gsub("[[:punct:]]|\\s", "_", temp), ".fastq")
    v$fq_fnames <- file.path(v$fq_dir, v$fq_fnames)
    v$ab1_fnames <- file.path(v$ab1_dir, v$ab1_fnames)

    dummy <- mapply(function(u,v,w) {
       CrispRVariants::abifToFastq(u,v,w)
    }, v$sq_nms, v$ab1_fnames, v$fq_fnames)
     
    v$fq_fnames <- unique(v$fq_fnames)

}

# ---------------
# run mapping
# ---------------
mapFastQ <- function(){
  
    # If bwa cannot be run, stop the app
    if (system("bwa", ignore.stderr = TRUE) == 127){
        createAlert(session, "alertAB1", "alertAB1_1", title = "WARNING",
           content = "Please make sure BWA is available", 
           style = "warning", append = FALSE, dismiss = FALSE)
        stopApp()
    }
 
    # Update the selected genome in the guide panel to match the mapping genome
    slct <- input$select_genome
    updateSelectInput(session, "select_Refgenome", selected = slct, choices = genlist)
  
    #BWA indices were generated using bwa version 0.7.10
    ind <- paste0(genome,"/", genlist[input$select_genome])

    if (file.exists(ind) && !is.null(v$fq_fnames)) {
        bwa_index <- ind
        bm_fnames <- gsub(".fastq$",".bam",basename(v$fq_fnames))
        v$srt_bm_names <- file.path(v$bam_dir, gsub(".bam","_s",bm_fnames))
        
        bm_fnames <- file.path(v$bam_dir, bm_fnames)
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Creating BAM files", value = 0)
      
        n <- length(v$fq_fnames)
      
        #Map, sort and index the bam files, remove the unsorted bams
        for(i in 1:length(v$fq_fnames))
        {
          progress$inc(1/n, detail = paste0("Mapping FASTQs ", i, "/", n))
          cmd <- paste0("bwa mem -t 2 ", bwa_index, " ", 
                    v$fq_fnames[i]," | samtools view -Sb - > ", bm_fnames[i])
          cat(cmd, "\n"); system(cmd)
          Rsamtools::indexBam(Rsamtools::sortBam(bm_fnames[i],v$srt_bm_names[i]))
          unlink(bm_fnames[i])
        }
        
        v$bm_fnames <- dir(v$bam_dir, ".bam$", full.names = TRUE, recursive = T)
        
        cat(sprintf("bam fnames %s \n", v$bm_fnames))
               
    } else {
        createAlert(session, "alertAB1", "prepAlertFasQ", title = "WARNING",
          content = "BWA index doesn't exist", style = "warning", append = FALSE)
    }
}
