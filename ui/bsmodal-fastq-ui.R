################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

upload_FastQ <- fileInput('fastq_files', 
  'Upload .fastq files in a zip file',
  accept = ".zip",
  multiple = TRUE,
  width = "100%")
select_genome <- selectInput("select_genome", "Select the genome", choices = genlist, width = "100%")

run_fastq <- bsButton("run_fastq", 'Run', icon =  icon("list-alt"), style = "primary", block = TRUE)

modal_FASTQ <- bsModal(
  "modal_FASTQ", "PREPROCESSING DATA | FASTQ FILES : ",
  "setting_btn", size = "small",
  bsAlert("alertFASTQ"),
  fluidRow(
    column(width = 12,
      h6("Step 1"),
      p("FASTQ file format (fastq) stores biologicalsequences and their quality scores in a simple plain text format "),
      upload_FastQ,
      h6("Step 2"),
      p("Chose the genome for mapping the FastQ files"),
      select_genome
    )
  ),
  fluidRow(
    column(width = 6),
    column(width = 6,
      run_fastq)
  )
  
)