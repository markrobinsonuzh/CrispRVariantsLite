################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

ref_seqs <- textInput("ref_seqs", "Guide Sequence", width = NULL, placeholder = "ATGCTGCTGGTTATTAGATTAGT")

select_Refgenome  <- selectInput("select_Refgenome", "Reference Genome", choices = genlist.gz, width = "100%")

txDb <- selectInput("txDb", "Transcript Database", choices = txDBl, selected = txDBl[1], width = "100%")

#run_guide <- bsButton("run_guide", 'create guides', type = "action", style = "success", block = TRUE)
run_guide <- actionButton("run_guide", 'create guides', width='100%')
next_step <- uiOutput("next_step")

g.start <- textInput("g.start", "START", value = d.start, placeholder = "20")
#g.length <- textInput("g.length", "WIDTH", value = d.length)
g.chr <- textInput("g.chr", "CHR", value = d.chr, placeholder = "chr12")
g.strand <- radioButtons("g.strand", "Strand", choices = c("+","-"), selected = d.strand, inline = T, width = NULL)

target_loc <- sliderInput("target_loc", "Target location", min = 0, max = 30, value = 17, step= 1)
g.length <- sliderInput("g.length", "WINDOW AROUND GUIDE", min = 0, max = 20, value = 5, step= 1)
guide = uiOutput("guide")

modal_ref <- bsModal(
  "modal_ref", "Specify the guide: ",
  "setting_btn", size = "large",
  bsAlert("alertRef"),
  fluidRow(
      column(width = 4,
       fluidRow(
    column(width = 12,
      select_Refgenome,
      txDb,
      helpText("Enter guide sequence"),
      ref_seqs,
      h1("or", align = "center"),
      helpText("Enter the chromosome, coordinates and strand of the guide sequence and WIDTH (number of bases on each side of PAM+guide)"),
      fluidRow(
        column(width = 6,
          g.chr,
          g.strand),
        column(width = 6,
          g.start
          )),
   fluidRow(
       column(width=12,
          p(),
          target_loc,
          g.length
          ))
          #)),
      #h1("or", align = "center"),
      #p(),
      #txDb
    )
  ),
  fluidRow(
    column(width = 6,
    next_step ),
    column(width = 6,
      run_guide)
  )),
      column(width = 8,
        guide
      )
  )
 
  
)
