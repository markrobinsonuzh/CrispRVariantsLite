################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

upload_AB1 <- fileInput('ab1_files', 'Upload .ab1 files in a .zip file', multiple = F, width = "100%")

select_genome <- selectInput("select_genome", "Select the genome", choices = genlist, width = "100%")

run_prep <- bsButton("run_prep", 'Run', icon =  icon("list-alt"), style = "primary", block = TRUE)

modal_AB1 <- bsModal(
  "modal_AB1", "PREPROCESSING DATA | AB1 FILES : ",
  "setting_btn", size = "small",
  bsAlert("alertAB1"),
  fluidRow(
    column(width = 12,
      h6("Step 1"),
      p("We convert AB1 files  ( format of the Applied Biosystems Sanger sequencer)to FASTQ format for mapping"),
      upload_AB1,
      h6("Step 2"),
      p("Chose the genome for mapping your FastQ file"),
      select_genome
      )
  ),
  fluidRow(
    column(width = 6),
    column(width = 6,
      run_prep)
  )
  
)
