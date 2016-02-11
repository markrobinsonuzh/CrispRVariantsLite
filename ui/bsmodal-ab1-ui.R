################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

upload_AB1 <- fileInput('ab1_files', 'Upload .ab1 files in a .zip file', multiple = F, width = "100%")

select_genome <- selectInput("select_genome", "Select the genome", choices = genlist, width = "100%")

run_prep <- bsButton("run_prep", 'Run', icon =  icon("list-alt"), style = "primary", block = TRUE)

modal_AB1 <- bsModal(
  "modal_AB1", "Preprocessiong The Data : ",
  "setting_btn", size = "small",
  bsAlert("alertAB1"),
  fluidRow(
    column(width = 12,
      h3("Step 1"),
      upload_AB1,
      h3("Step 2"),
      select_genome
      )
  ),
  fluidRow(
    column(width = 6),
    column(width = 6,
      run_prep)
  )
  
)
