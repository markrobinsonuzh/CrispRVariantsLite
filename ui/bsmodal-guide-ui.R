################################################################################
# Define  options modal get AB1 files, sequence names, FASTQ files, 
################################################################################

ref_seqs <- textInput("ref_seqs", "Reference", width = NULL, placeholder = "ATGCTGCTGGTTATTAGATTAGT")

select_Refgenome  <- selectInput("select_Refgenome", "Reference Genome", choices = genlist.gz, width = "100%")

txDb <- selectInput("txDb", "Transcript Database", choices = txDBl, selected = txDBl[1], width = "100%")

run_guide <- bsButton("run_guide", 'create guides', icon =  icon("list-alt"), style = "primary", block = TRUE)

g.start <- textInput("g.start", "START", value = "38948341", placeholder = "20")
g.length <- textInput("g.length", "WIDTH", value = "0")
g.chr <- textInput("g.chr", "CHR", value = "chr12", placeholder = "chr12")
g.strand <- radioButtons("g.strand", "Strands", choices = c("+","-"), selected = "+", inline = T, width = NULL)

modal_ref <- bsModal(
  "modal_ref", "Create the reference guide : ",
  "setting_btn", size = "small",
  bsAlert("alertRef"),
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
      h1("or", align = "center"),
      helpText("Type the sequence of interest"),
      ref_seqs,
      select_Refgenome,
      p(),
      txDb,
      uiOutput("error1")
    )
  ),
  fluidRow(
    column(width = 6),
    column(width = 6,
      run_guide)
  )
  
)
