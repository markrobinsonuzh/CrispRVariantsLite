################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

ref_seqs <- textInput("ref_seqs", "Reference", width = NULL, placeholder = "ATGCTGCTGGTTATTAGATTAGT")

select_Refgenome  <- selectInput("select_Refgenome", "Reference Genome", choices = genlist.gz, width = "100%")

txDb <- selectInput("txDb", "Transcript Database", choices = txDBl, selected = txDBl[1], width = "100%")

#run_guide <- bsButton("run_guide", 'create guides', type = "action", style = "success", block = TRUE)
run_guide <- actionButton("run_guide", 'create guides', width='100%')
next_step <- uiOutput("next_step")

g.start <- textInput("g.start", "START", value = d.start, placeholder = "20")
g.length <- textInput("g.length", "WIDTH", value = d.length)
g.chr <- textInput("g.chr", "CHR", value = d.chr, placeholder = "chr12")
g.strand <- radioButtons("g.strand", "Strands", choices = c("+","-"), selected = d.strand, inline = T, width = NULL)

target_loc <- sliderInput("target_loc", "Target location", min = 0, max = 30, value = 17, step= 1)
guide = uiOutput("guide")

modal_ref <- bsModal(
  "modal_ref", "Create the reference guide : ",
  "setting_btn", size = "large",
  bsAlert("alertRef"),
  fluidRow(
      column(width = 4,
       fluidRow(
    column(width = 12,
      helpText("Type de coordinate of the sequence of interest"),
      fluidRow(
        column(width = 6,
          g.start,
          g.length),
        column(width = 6,
          g.chr,
          g.strand
          )),
   fluidRow(
       column(width=12,
          target_loc
          )),
      h1("or", align = "center"),
      helpText("Type the sequence of interest"),
      ref_seqs,
      select_Refgenome,
      p(),
      txDb
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
