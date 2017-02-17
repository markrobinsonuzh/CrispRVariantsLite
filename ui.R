## ui.R ##

header <- dashboardHeader(title = "CrispRVariantsLite")

header$children[[3]]$children <-  tags$a(href='http://www.imls.uzh.ch/index.html',
    img(src='uzh_logo_d_pos_web_main_zone.jpg', height='60', align="right"))

# Institute of Molecular Life Sciences
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
#source("ui/bsmodal-help-ui.R", local = T)
source("ui/bsmodal-save-ui.R", local = T)
source("ui/bsmodal-reset-ui.R", local = T)

body <- dashboardBody(
  fluidRow(
    
    box(width = 8, height = "100%",
      solidHeader = T,
      bsAlert("AlertUI"),
      wellPanel(htmlOutput("plots"))
      ),
    plotOptions,
    modal_1,
    modal_2,
    modal_AB1,
    modal_FASTQ,
    modal_table,
    modal_ref,
    #modal_help,
    modal_save,
    modal_reset
  )
)

dashboardPage(header, sidebar, body, skin = "black")
