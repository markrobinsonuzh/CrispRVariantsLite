################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

#--------------
# fileInput
#--------------
upload_bams <- htmlOutput("bams")
reference <- fileInput('reference', 'Upload reference', multiple = F, width = "100%")


targert_seq <- selectInput("target_seq", "Target sequence", choices = list(), width = "100%")



target_loc <- sliderInput("target_loc", "Target location",
  min = 0, max = 30, value = 22, step= 1)


#--------------
# button
#--------------

reset <- bsButton("reset", "Reset", icon = icon("fa fa-ellipsis-h"), style = "warning", block = TRUE)
run_plot <- bsButton("run_plot", 'Plot', icon =  icon("area-chart"), style = "success", block = TRUE)
select_data <- bsButton("select_data", 'select data', icon =  icon("fa fa-pencil"), style = "primary", block = TRUE)


top.n = numericInput("top.n", "Top ranked variants", value = 50 )
min.freq = numericInput("min.freq", "Frequency cutoff", value = 0 )
min.count = numericInput("min.count", "No count cutoff", value = 0 )

ref_plot = plotOutput("ref_plot")

#sample_names <- textInput("sample_names", "Sample Names")
create_guides <- bsButton("create_guides", "Create guides", icon =  icon("fa fa-list"), style = "info", block = TRUE)
modal_2 <- bsModal(
  "modal_2", "Figures setting : ",
  "data_setting", size = "large",
  fluidRow(column(width = 12,
    bsAlert("alert2"),
    helpText("CRISPR-Cas9 mutagenesis experiment. Starting from a set of bam files, CrispRVariants narrows the alignments to the target region and renumbers variants with respect to a zero point, typically the Cas9 cut site three bases upstream of the protospacer adjacent motif (PAM). The target region and zero point can be flexibly specified, and the variants may be shown with respect to either DNA strand.")
    )),
  hr(),
  fluidRow(
    column(width = 3,
      h3("Load data : "),
      fluidRow(
        column(width = 12,
          upload_bams,
          p(),
          create_guides,
          p(),
          uiOutput("metadata")
        ))
    ),
    column(width = 3,
      h3("Plot setting : "),
      fluidRow(column(width = 12, helpText("Table of counts options"))),
          top.n,
          min.freq,
          min.count,
          target_loc
      ),
    column(width = 6,
      ref_plot)
  ),
  hr(),
  fluidRow(
    column(width = 6,
      fluidRow(
        column(width = 6),
        column(width = 6,
          reset)
      )),
    column(width = 6,
      fluidRow(
        column(width = 6,
          select_data),
        column(width = 6, 
          run_plot)
        )
      )
    )
  )
