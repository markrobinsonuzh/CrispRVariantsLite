################################################################################
# Panel for entering and plotting the guide sequence
################################################################################

# List of available genomes


select_Refgenome  <- selectInput("select_Refgenome", "Reference Genome",
                                 choices = genlist, width = "100%")

# Option 1 for guide input: enter the sequence
ref_seqs <- textInput("ref_seqs", "Guide Sequence (inc PAM)", width = "100%",
                      value = NULL, placeholder = "ATGCTGCTGGTTATTAGATTAGT")

# Option 2 for guide input: enter the location
g.start <- textInput("g.start", "START", value = NULL)
g.chr <- textInput("g.chr", "CHR", value = NULL, placeholder = "chr12")
g.strand <- radioButtons("g.strand", "Strand", choices = c("+","-"),
                          selected = "-", inline = T, width = NULL)

# Location of cut site within target and padding around guide sequence
target_loc <- sliderInput("target_loc", "Target location", min = 0, max = 30,
                          value = 17, step= 1)
g.length <- sliderInput("g.length", "WINDOW AROUND GUIDE", min = 0, max = 20,
                        value = 5, step= 1)


# Action buttons - either create the guide or create plot when guide is okay
run_guide <- bsButton("run_guide", 'create guides', type = "action",
                      style = "success", block = TRUE)
run_plot_guide <- bsButton("run_plot_guide", 'create plot', type = "action",
                      style = "success", block = TRUE)

next_step <- uiOutput("next_step")

# Windows for plotting the guide and the transcripts
guide <- uiOutput("guide")
plot_anot <- plotOutput("plot_anot", width="100%", height="200px")

#___________________________________
# Panel arranging the above elements
modal_ref <- .bsModal(
  "modal_ref", "Specify the guide: ",
  "setting_btn", width = "800px", close.button = FALSE,
  fluidRow(
    column(width = 6,
      fluidRow(
        # Left column containing fields for user to enter the guide sequence
        column(width = 12,
          fluidRow(column(width=6, select_Refgenome), column(width=6)),
          fluidRow(column(width=12, helpText("Enter guide sequence"), ref_seqs)),
          h1("or", align = "center"),
          helpText(paste("Enter the chromosome, coordinates and strand of the guide",
                     "sequence and WIDTH (number of bases on each side of PAM+guide)")),
          fluidRow(
            column(width = 6,
              fluidRow(column(width=6, g.chr), column(width=6, g.strand))),
            column(width = 6, g.start)
          ),
          fluidRow(
            column(width=6, target_loc),
            column(width=6,g.length)
          )
        )
    )),
    # Right column containing plots showing plots of guide and transcript locations
    column(width = 6,
      wellPanel(bsAlert("alertRef"), guide, p(), plot_anot)
    )
  ),
  # Footer containing action buttons
  footer = fluidRow(
    column(width = 3),
    column(width = 3, run_guide),
    column(width = 3, run_plot_guide),
    column(width = 3, next_step)
  )
)
