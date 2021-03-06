################################################################################
# Dialogue for uploading fastq files 
################################################################################

# List of available genomes
#select_genome  <- selectInput("select_Refgenome", "Reference Genome",
#                                 choices = genlist, width = "100%")


upload_FastQ <- fileInput("fq_input", 'Upload .fastq files in a zip file',
                          accept = ".zip", multiple = FALSE, width = "100%")

select_genome_fq <- selectInput("select_genome", "Select the genome",
                             choices = genlist_base, width = "100%")

back_fastq <- bsButton("back_fastq", 'Back', type = "action",
                       style = "default", block = TRUE)
run_fastq <- bsButton("run_fastq", 'Run', type = "action",
                      style = "success", block = TRUE)

modal_FASTQ <- .bsModal(
    "modal_FASTQ", "PREPROCESSING DATA | FASTQ FILES : ",
    "setting_btn", size = "small", close.button = FALSE,
    bsAlert("alertFASTQ"),
    fluidRow(
      column(width = 12,
        h6("Step 1"),
        p(paste("FASTQ file format (fastq) stores biological sequences and",
                "their quality scores in a simple plain text format ")),
        upload_FastQ,
        h6("Step 2"),
        p("Chose the genome for mapping the FastQ files"),
        select_genome_fq
      )
    ),
    footer = fluidRow(
      column(width = 6, run_fastq),
      column(width = 6, back_fastq)
    )
)
