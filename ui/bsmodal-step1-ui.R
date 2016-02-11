################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

upload_AB1 <- fileInput('ab1_files', 'Upload .ab1 files', multiple = TRUE, width = "100%")
upload_BAM <- fileInput('bam_files', 'Upload .bam files', multiple = TRUE, width = "100%")
upload_FastQ <- fileInput('fastq_files', 'Upload .fastq files', multiple = TRUE, width = "100%")
upload_Rdata <- fileInput('rdata_files', 'Upload .rdata files', multiple = TRUE, width = "100%")


select_AB1 <- bsButton("select_AB1", 'ab1 files', icon = icon("glyphicon glyphicon-plus"), style = "info", block = TRUE)
select_BAM <- bsButton("select_BAM", 'bam files', icon = icon("glyphicon glyphicon-plus"), style = "info", block = TRUE)
select_FastQ <- bsButton("select_FastQ", 'fastq files', icon = icon("glyphicon glyphicon-plus"), style = "info", block = TRUE)
select_Rdata <- bsButton("select_Rdata", 'Preprocess Data', style = "info", block = TRUE)

modal_1 <- bsModal(
  "modal_1", "CrispRVariants : ",
  "setting_btn", size = "small",
  fluidRow(column(width = 12,
    helpText("CrispRVariants is an R-based toolkit for counting, localising and plotting targeted insertion and deletion variant combinations."),
    p("Select the format of your data")
    )),
  fluidRow(
    column(width = 6,
      select_AB1
      ),
    column(width = 6,
      select_BAM
      )
  ),
  p(" "),
  fluidRow(
    column(width = 6,
      select_FastQ
      ),
    column(width = 6,
      select_Rdata
      )
  )
)
