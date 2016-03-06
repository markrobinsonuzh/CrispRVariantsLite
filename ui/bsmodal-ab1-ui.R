################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

upload_AB1 <- fileInput('ab1_files', 'Upload ZIP file with .AB1 files in directories', multiple = F, width = "100%")

select_genome <- selectInput("select_genome", "Select genome", choices = genlist, width = "100%")

run_prep <- bsButton("run_prep", 'Run',  style = "success", block = TRUE)
back_ab1 <- bsButton("back_ab1", 'Back',  style = "default", block = TRUE)
modal_AB1 <- .bsModal(
  "modal_AB1", "PREPROCESSING DATA | ZIP of directories with AB1 FILES: ",
  "setting_btn", size = "small",
  bsAlert("alertAB1"),
  fluidRow(
    column(width = 12,
      p("Step 1: convert AB1 files (Applied Biosystems Sanger sequencer) to FASTQ format for mapping"),
      upload_AB1,
      p("Step 2: chose genome for mapping FASTQ files"),
      select_genome
      )
  ),
  footer = fluidRow(
    column(width = 6,
    run_prep),
    column(width = 6,
    back_ab1  )
  ), close.button = F)
