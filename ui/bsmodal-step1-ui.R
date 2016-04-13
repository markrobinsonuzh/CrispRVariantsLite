################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################



start_btn <- bsButton("start_btn", "Start", style = "info", block = TRUE)


bscollapse_video <- bsCollapse(id = "bscollapse_video", open = "Overall workflow using BAM files", 
    bsCollapsePanel("Overall workflow using BAM files",
    #helpText(paste("To create CrispRVariants plot you need to create a guide.",
    #               "Click 'RePlot' to replot the data with the new options.")),
    #hr(),
    img(src="1.gif", width="100%"),
    style = "info"),
    

    bsCollapsePanel("Changing Plot Options",
    helpText(paste("These options control the layout, filtering, allele plot and the heatmapt.",
                   "Click 'RePlot' to replot the data with the new options.")),
    hr(),
    img(src="2.gif", width="100%"),
    style = "warning")
)

        

modal_1 <- .bsModal(
  "modal_1", "CrispRVariantsLite",
  "setting_btn", 
  fluidRow(column(width = 12,
    tabBox( width = 12,
      tabPanel("QUICK GUIDE",
      fluidRow(column(width=12,
        p("CrispRVariantsLite is an interactive web application that provides a graphical user interface to the R/Bioconductor CrispRVariants package. "),
        p("The starting point is a ZIP file that contains one of three possibilities: i) a set of already-mapped sequences in BAM format (one BAM file for each sample);  ii) a set of FASTQ files to be mapped (one FASTQ for each sample); iii) a set of directories that contain (Applied Biosytems) AB1 files (each directory represents a sample and contains 1 or more AB1 files).  users will need to know the genomic location of the target sequence (\"guide\").  Once the set of reads are loaded, the user will be given an opportunity to specify labels for each sample (by default, these are taken from the file or directory names) and a grouping variable for colouring the label.  In the allele variant plot, variant sequences will be renumbered using the target.loc as the zero point.  Note that for a 20 nucleotide guide RNA, the target location is base 17, where the double-strand cut occurs.")
        ))),
      tabPanel("QUICK INSTRUCTION",
        bscollapse_video),
      tabPanel("EXAMPLE ZIP FILES",
        includeMarkdown("examples/README.md")
      ),
      tabPanel("ABOUT",
        includeMarkdown("help_file/about.md")
      )
    )
    )),
    footer = tags$footer(
        fluidRow(
    column(width = 4),
    column(width = 4, 
           start_btn),
    column(width = 4)
    )), close.button = F)
