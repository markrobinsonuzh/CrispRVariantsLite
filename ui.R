## ui.R ##

header <- dashboardHeader(title = "CrispRVariantsLite")

sidebar <- dashboardSidebar(disable = T)

#---------------------
# Source body
#---------------------
source("ui/bsmodal-step1-ui.R", local = T)
source("ui/bsmodal-fastq-ui.R", local = T)
source("ui/bsmodal-ab1-ui.R", local = T)
source("ui/bsmodal-step2-ui.R", local = T)
source("ui/bsmodal-table-ui.R", local = T)
source("ui/plotoptions-ui.R", local = T)
source("ui/bsmodal-guide-ui.R", local = T)
source("ui/bsmodal-help-ui.R", local = T)
source("ui/bsmodal-save-ui.R", local = T)

body <- dashboardBody(
  fluidRow(
    box(width = 8, height = "600px",
      solidHeader = T,
      htmlOutput("plots")),
    plotOptions,
    modal_1,
    modal_2,
    modal_AB1,
    modal_FASTQ,
    modal_table,
    modal_ref,
    modal_help,
    modal_save
  )
)

dashboardPage(header, sidebar, body)
