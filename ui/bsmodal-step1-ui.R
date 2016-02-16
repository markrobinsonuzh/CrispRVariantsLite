################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################



start_btn <- bsButton("start_btn", "Start CrisprRvariant light - 0.0.1", style = "info", block = TRUE)

modal_1 <- bsModal(
  "modal_1", "CrispRVariants : ",
  "setting_btn", size = "large",
  fluidRow(column(width = 12,
    tabBox( width = 12,
      tabPanel("Quick Guide",
        p("CrispRVariantsLite is an interactive web application that provides a graphical user interface to the microbiome analysis package for R, called CrispRVariants."),
        p("The usual starting point is a list of bam filenames (\"bams\"), the location of the target sequence (\"guide\") as a GenomicRanges::Granges object, and the reference sequence, as a Biostrings::DNAString (\"ref\"). You can optionally specify a vector of sample names (\"sample_names\") corresponding to the bam files, and the target location (target.loc) with respect to the reference. Variant sequences will be renumbered using the target.loc as the zero point. For a 20 nucleotide guide RNA, the target location is base 17, where the double-strand cut occurs. To count variants within the target region:")
        ),
      tabPanel("1. UPLOADING BAM FILES"#,
        #includeMarkdown("help_file/uploading_bams.md")
        ),
      tabPanel("2 EDITING THE METADATA",
        p()
        ),
      tabPanel("3 CREATING THE GUIDE",
        p()
        ),
      tabPanel("4 CREATE THE PLOTE",
        p()
        ),
      tabPanel("About"#,
        #includeMarkdown("help_file/about.md")
      )
    )
    )),
  fluidRow(
    column(width = 4),
    column(width = 4, start_btn),
    column(width = 4)
  )
  #,
  #p(" "),
  #fluidRow(
  #  column(width = 6,
  #    select_FastQ
  #    ),
  #  column(width = 6,
  #    select_Rdata
  #    )
  #)
)
