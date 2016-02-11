temp.dir <- tempdir()

# create the temp dir for the files
setDir <- reactive({
  
  bam_dir <- file.path(temp.dir, "bam")
  ifelse(!dir.exists(bam_dir), dir.create(bam_dir,  showWarnings = FALSE), FALSE)
  v$bam_dir <- bam_dir
  
  bamtemp_dir <- file.path(temp.dir, "bamtemp")
  ifelse(!dir.exists(bamtemp_dir), dir.create(bamtemp_dir,  showWarnings = FALSE), FALSE)
  v$bamtemp_dir <- bamtemp_dir
  
  fq_dir <- file.path(temp.dir, "fastq")
  ifelse(!dir.exists(fq_dir), dir.create(fq_dir,  showWarnings = FALSE), FALSE)
  v$fq_dir <- fq_dir
  
  ab1_dir <- file.path(temp.dir, "ab1")
  ifelse(!dir.exists(ab1_dir), dir.create(ab1_dir,  showWarnings = FALSE), FALSE)
  v$ab1_dir <- ab1_dir
  
})

# ---------------
# convert AB1 to FASTQ
# ---------------
convertAb1toFasq <- reactive({
  data_dir <- input$ab1_files
   if(!is.null(data_dir)){
     # list AB1 files
     temp <- unzip(data_dir$datapath, exdir = v$ab1_dir)
     
     v$ab1_fnames <- dir(v$ab1_dir, "ab1$", recursive = TRUE, full.names = TRUE)
     # get the sequence names
     v$sq_nms <- gsub(".ab1","",basename(v$ab1_fnames))
     
     # replace spaces ans slashes in filename with underscores
     v$fq_fnames  <- paste0(gsub("[\ |\\/]", "_", dirname(v$ab1_fnames)), ".fastq")
     print(dirname(v$ab1_fnames))
     v$fq_fnames  <- gsub("ab1_","",v$fq_fnames)
     v$fq_fnames <- file.path(v$fq_dir,v$fq_fnames)
     dummy <- mapply( function(u,v,w) {
       abifToFastq(u,v,w)
     }, v$sq_nms, v$ab1_fnames, v$fq_fnames)
     
   } else {
     #throw error : no file uploaded
   }
})


# ---------------
# run mapping
# ---------------
mapFastQ <- reactive({
  #BWA indices were generated using bwa version 0.7.10
  ind <- paste0("/home/Shared_taupo/data/annotation/Danio_rerio/genome_danRer7/", 
                input$select_genome)

    if(file.exists(ind))
    {
        bwa_index <- ind
        v$bm_fnames <- gsub(".fastq$",".bam",basename(v$fq_fnames))
        v$srt_bm_names <- file.path(v$bam_dir, gsub(".bam","_s",v$bm_fnames))
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        progress$set(message = "Preprocessing  part B ", value = 0)
      
        n <- length(v$fq_fnames)
      
        #Map, sort and index the bam files, remove the unsorted bams
        for(i in 1:length(v$fq_fnames))
        {
          progress$inc(1/n, detail = paste0("Mapping - part ", i))
          cmd <- paste0("bwa mem ", bwa_index, " ", v$fq_fnames[i]," | samtools view -Sb - > ", 
                         v$bamtemp_dir,"/",v$bm_fnames[i])
          message(cmd, "\n"); system(cmd)
          indexBam(sortBam(paste0(v$bamtemp_dir,"/",v$bm_fnames[i]),v$srt_bm_names[i]))
          unlink(paste0(v$bamtemp_dir,"/",v$bm_fnames[i]))
        }
    }
    else
    {
    createAlert(session, "alertAB1", "prepAlertFasQ", title = "WARNING",
      content = "BWA index doesn't exist", style = "warning", append = FALSE)
    }
})
