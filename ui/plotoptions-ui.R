#----------------
# Action buttons
#----------------

data_setting <- bsButton("data_setting", 'Settings', type='action',
                         icon =  icon("list-alt"), style = "primary",
                         block = TRUE)
save_data <- bsButton("save_data", 'Save Data', type='action',
                      icon =  icon("fa fa-pencil"), style = "info", block = TRUE)

replot <- bsButton("replot", 'Plot', type='action', icon =  icon("area-chart"),
            style = "info", block = TRUE)
            
reset <- bsButton("reset", 'Reset', type='action', style = "info", block = TRUE)
update <- bsButton("update", 'Update', type='action', style = "primary", block = TRUE)


#----------------
# Filter buttons
#----------------

top.n = numericInput("top.n", "Top ranked variants", value = 50, min = 0)
min.freq = numericInput("min.freq", "Frequency cutoff", value = 0, min = 0)
min.count = numericInput("min.count", "Count cutoff", value = 0, min = 0)


#----------------
# Plot arrangement
#----------------

# Row height ratio
row.ht.ratio <- selectInput("row.ht.ratio","Row height ratio", selected = "1:6",
    choices = c("1:6","1:5","1:4","1:3","1:2","1:1","2:1"), width = 75)

# Column width ratio
col.wdth.ratio <- selectInput("col.wdth.ratio", "Column width ratio",
    selected = "2:1", choices = c("4:1","3:1","2:1","1:1", "1:2"), width = 75)

# Plot margins
plot.margins <- selectInput("plot.margins", "Space beneath plot",
    selected = "Normal", choices = c("Less", "Normal", "More"), width = 100)


#----------------
# Slider
#----------------

axis.text.size  <- sliderInput("axis.text.size",
                               "Text size of variant labels (e.g., 3:3D):",
                               min = 1, max = 20, value = 10)

# Decimal interval with step value
plot.text.size <- sliderInput("plot.text.size",
                             "Text size of sequence letters / numbers in plots: ", 
                              min = 0, max = 10, value = 4, step = 1)

# Specification of range within an interval
legend.text.size <- sliderInput("legend.text.size", "Legend text size :",
                                min = 1, max = 20, value = 8)

ins.size <- sliderInput("ins.size", "Size of insertion symbols: ",
                        min = 1, max = 10, value = 3)



# Simple integer interval
x.size <- sliderInput("x.size", "X-axis label text size:", min=0, max=20,
      value=10, step = 1)

# Provide a custom currency format for value display, 
# with basic animation
x.angle <- sliderInput("x.angle", "Angle of sample label:", min = 0, max = 90,
                       value = 90, step = 5, animate=TRUE)
                       
plot.type <- radioButtons("plot.type", "Plot type",
                          choices = c("proportions","counts"),
                          selected = "counts", inline = TRUE) 


legend.symbol.size <- sliderInput("legend.symbol.size", "legend symbol size :",
                                  min = 1, max = 20, value = 5)

legend.ncol <- sliderInput("legend.ncol", "Number of legend columns:",
                            min = 1, max = 20, value = 4)



#----------------
# Sidebar panel
#----------------

bscollapse_1 <- bsCollapse(id = "bscollapse_1", open = "Plot Alignments Options",

  bsCollapsePanel("Layout options",
    helpText(paste("These options control the relative sizes of the plots.",
                   "Click 'Replot' to replot the data with the new options.")),
    hr(),
    fluidRow(
        column(width = 4, row.ht.ratio),
        column(width = 4, col.wdth.ratio),
        column(width = 4, plot.margins)
        ),
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
    legend.symbol.size,
    legend.ncol,
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
    column(width=6, data_setting),
    column(width=6, replot)
  ),
  hr(),
  fluidRow(
    column(width=6, save_data ),
    column(width=6, reset)
))
