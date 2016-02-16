################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################



start_btn <- bsButton("start_btn", "Start", style = "info", block = TRUE)

modal_1 <- bsModal(
  "modal_1", "CrispRVariants Lite ",
  "setting_btn", size = "large",
  fluidRow(column(width = 12,
    tabBox( width = 12,
      tabPanel("QUICK GUIDE",
      fluidRow(column(width=12,
        p("CrispRVariantsLite is an interactive web application that provides a graphical user interface to the microbiome analysis package for R, called CrispRVariants."),
        p("The usual starting point is a list of bam filenames (\"bams\"), the location of the target sequence (\"guide\") as a GenomicRanges::Granges object, and the reference sequence, as a Biostrings::DNAString (\"ref\"). You can optionally specify a vector of sample names (\"sample_names\") corresponding to the bam files, and the target location (target.loc) with respect to the reference. Variant sequences will be renumbered using the target.loc as the zero point. For a 20 nucleotide guide RNA, the target location is base 17, where the double-strand cut occurs. To count variants within the target region.")
        ))),
      tabPanel("QUICK INSTRUCTION",
        img(src="2.gif", width="100%")
        ),
      tabPanel("ABOUT",
        includeMarkdown("help_file/about.md")
      )
    )
    )),
  fluidRow(
    column(width = 4),
    column(width = 4, 
           start_btn),
    column(width = 4)
  )
)
