################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

#--------------
# fileInput
#--------------
upload_bams <- htmlOutput("bams")
#reference <- fileInput('reference', 'Upload reference', multiple = F, width = "100%")


targert_seq <- selectInput("target_seq", "Target sequence", choices = list(), width = "100%")





info <- checkboxInput("info", "Show help box", value=FALSE, width=NULL)
#--------------
# button
#--------------

select_FastQ <- bsButton("select_FastQ", 'upload FastQ', style = "primary", block = TRUE)
select_AB1 <- bsButton("select_AB1", 'upload AB1', style = "primary", block = TRUE)

reset <- bsButton("reset", "Reset", icon = icon("fa fa-ellipsis-h"), style = "warning", block = TRUE)
run_plot <- bsButton("run_plot", 'Plot', icon =  icon("area-chart"), style = "success", block = TRUE)
select_data <- bsButton("select_data", 'instructions', icon =  icon("fa fa-pencil"), style = "primary", block = TRUE)


top.n = numericInput("top.n", "Top ranked variants", value = 50 )
min.freq = numericInput("min.freq", "Frequency cutoff", value = 0 )
min.count = numericInput("min.count", "No count cutoff", value = 0 )

ref_plot = plotOutput("ref_plot")

#sample_names <- textInput("sample_names", "Sample Names")
create_guides <- bsButton("create_guides", "Create guides", icon =  icon("fa fa-list"), style = "info", block = TRUE)
edit_xls <- bsButton("edit_xls", "metadata", type="action", icon =  icon("table"), style = "success", block = TRUE) 

modal_2 <- .bsModal(
  "modal_2", "Figures setting : ",
  "data_setting", size = "large",
  fluidRow(column(width = 12,
    bsAlert("alert2"),
    helpText("CrispRVariants narrows the alignments to the target region and renumbers variants with respect to a zero point, typically the Cas9 cut site three bases upstream of the protospacer adjacent motif (PAM). The target region and zero point can be flexibly specified, and the variants may be shown with respect to either DNA strand.")
    )),
  hr(),
  fluidRow(
    column(width = 3,
      h3("Load data : "),
      fluidRow(
        column(width = 12,
          upload_bams,
          create_guides,
          p(),
          edit_xls
        ))
    ),
    column(width = 3,
    h3("Filters: "),
      fluidRow(
          column(width = 12, helpText("Table of counts options"))),
          top.n,
          min.freq,
          min.count
      ),
      column(width = 6,
      h3("Guide Plot"),
      ref_plot
  ))
    , footer =   fluidRow(
    column(width = 6,
      fluidRow(
        column(width = 6),
        column(width = 6,
          info)
      )),
    column(width = 6,
      fluidRow(
        column(width = 6,
          select_data),
        column(width = 6, 
          run_plot)
        )
      )
    ), close.button = F)
