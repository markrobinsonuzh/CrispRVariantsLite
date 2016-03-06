data_setting <- bsButton("data_setting", 'Settings', type='action', icon =  icon("list-alt"), style = "primary", block = TRUE)
info_btn <- bsButton("info_btn", 'Help', type='action', icon =  icon("fa fa-pencil"), style = "primary", block = TRUE)
save_data <- bsButton("save_data", 'Save Data', type='action', icon =  icon("fa fa-pencil"), style = "info", block = TRUE)
replot <- bsButton("replot", 'Plot', type='action', icon =  icon("area-chart"),
            style = "info", block = TRUE)


#----------------
# Filter buttons
#----------------

top.n = numericInput("top.n", "Top ranked variants", value = 50, min = 0)
min.freq = numericInput("min.freq", "Frequency cutoff", value = 0, min = 0)
min.count = numericInput("min.count", "Count cutoff", value = 0, min = 0)

#----------------
# Select boxes
#----------------

# Row height ratio
row.ht.ratio <- selectInput("row.ht.ratio","Row height ratio", selected = "1:6",
    choices = c("1:6","1:5","1:4","1:3","1:2","1:1.5", "1:1", "2:1"), width = 75)

# Column width ratio
col.wdth.ratio <- selectInput("col.wdth.ratio", "Column width ratio",
    selected = "2:1", choices = c("4:1","3:1","2:1", "1.5:1", "1:1", "1:2"), width = 75)

reset <- bsButton("reset", 'Reset', type='action', style = "info", block = TRUE)
update <- bsButton("update", 'Update', type='action', style = "primary", block = TRUE)

replot <- bsButton("replot", 'Plot', type='action', icon =  icon("area-chart"),
            style = "info", block = TRUE)

#----------------
# Select boxes
#----------------

# Row height ratio
row.ht.ratio <- selectInput("row.ht.ratio","Row height ratio", selected = "1:6",
    choices = c("1:6","1:5","1:4","1:3","1:2","1:1","2:1"), width = 75)

# Column width ratio
col.wdth.ratio <- selectInput("col.wdth.ratio", "Column width ratio",
    selected = "2:1", choices = c("4:1","3:1","2:1","1:1", "1:2"), width = 75)


#----------------
# Slider
#----------------

# Simple integer interval
x.size <- sliderInput("x.size", "X-axis label text size:", min=0, max=20,
      value=10, step = 1)

# Decimal interval with step value
plot.text.size <- sliderInput("plot.text.size", "Text size of sequence letters / numbers in plots: ", 
    min = 0, max = 10, value = 2, step = 1)

# Specification of range within an interval
legend.text.size <- sliderInput("legend.text.size", "Legend text Size :", min = 1, max = 20, value = 8)

# Provide a custom currency format for value display, 
# with basic animation
x.angle <- sliderInput("x.angle", "Angle of sample label:", min = 1, max = 90, value = 90, step = 5, animate=TRUE)

axis.text.size  <- sliderInput("axis.text.size", "Text size of variant labels (e.g., 3:3D):", min = 1, max = 20, value = 10)

ins.size <- sliderInput("ins.size", "Size of insertion symbols: ", min = 1, max = 10, value = 3)
legend.symbol.size <- sliderInput("legend.symbol.size", "legendsymbol size :", min = 1, max = 20, value = 5)
legend.text.size <- sliderInput("legend.text.size", "Size of insertion legend:", min = 1, max = 20, value = 8)

plot.type <- radioButtons("plot.type", "Plot type", choices = c("proportions","counts"), selected = "counts", inline = TRUE) 


#----------------
# Sidebar panel
#----------------

bscollapse_1 <- bsCollapse(id = "bscollapse_1", open = "Plot Alignments Options",

  bsCollapsePanel("Layout options",
    helpText(paste("These options control the relative sizes of the plots.",
                   "Click 'Plot' to replot the data with the new options.")),
    hr(),
    row.ht.ratio,
    col.wdth.ratio,
    style = "info"),
  
  bsCollapsePanel("Filtering options",
    helpText(paste("These options control how many variant alleles are shown", 
                   "in the plot. Click 'Plot' to replot.")),
    hr(),
    top.n,
    min.freq,
    min.count,
    style = "danger"),
  
  bsCollapsePanel("Allele plot options",
    helpText(paste("These options control the appearance of the allele plot.",  
                   "Click 'Plot' to replot.")),

    hr(),
    helpText("Text sizes"),
    axis.text.size,
    plot.text.size,
    legend.text.size,
    hr(),
    ins.size,
    style = "warning"),
    
  bsCollapsePanel("Heatmap options", 
    helpText(paste("The size of the text within the plot is controlled",
                   "in the 'Allele plot options' panel")),
    hr(),
    x.size,
    x.angle,
    hr(),
    plot.type,
    style = "success" ))


#----------------
# Plot box
#----------------

plotOptions <- box(width = 4, 
  solidHeader = T,
  bscollapse_1,
  fluidRow(
    column(width=6,data_setting),
    column(width=6, save_data )
  ),
  hr(),
  fluidRow(
    column(width=6,replot),
    column(width=6,reset)
))
