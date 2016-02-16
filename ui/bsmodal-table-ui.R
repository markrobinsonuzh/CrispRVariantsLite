################################################################################
# Define  options modal Table 
################################################################################

#upload_table <- fileInput('upload_table', 'Apply table', multiple = TRUE, width = "100%")

upload_metadata <- fileInput('upload_metadata', 'Upload metadata', multiple = F, width = "100%")
select_genome <- selectInput("select_genome", "Select the genome", choices = genlist, width = "100%")

update_xls <- bsButton("update_xls", 'save data', icon =  icon("list-alt"), style = "primary", block = TRUE)

modal_table <- bsModal(
  "modal_table", "Edit Table : ",
  "nobtn", size = "large",
  fluidRow(
    column(width = 12,
      helpText("Changes to the table will be automatically saved to the source file."),
      htmlOutput("table")
    )
  ),
  p(),
  p(),
  hr(),
  fluidRow(
    column(width = 6,
      fluidRow(
        column(width = 6,
          radioButtons("filetype", "File type:",
            choices = c("csv", "tsv"), inline = T),
          downloadButton('downloadTable', 'Download')),
        column(width = 6)
        )),
    column(width = 6,
      fluidRow(
        column(width = 6),
        column(width = 6 )
      ))
  )
)
cd