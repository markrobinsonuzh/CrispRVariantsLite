data_setting <- bsButton("data_setting", 'Settings', type='action', icon =  icon("list-alt"), style = "primary", block = TRUE)
info_btn <- bsButton("info_btn", 'Help', type='action', icon =  icon("fa fa-pencil"), style = "primary", block = TRUE)
save_data <- bsButton("save_data", 'Save Data', type='action', icon =  icon("fa fa-pencil"), style = "info", block = TRUE)
reset_btn <- bsButton("reset_btn", "Reset", type='action', style='primary', block = TRUE)

#----------------
# Slider
#----------------

# Simple integer interval
x.size <- sliderInput("x.size", "Line weight:", min=0, max=20, value=10, step = 1)

# Decimal interval with step value
plot.text.size <- sliderInput("plot.text.size", "Text size of sequence letters / numbers in heatmap: ", min = 0, max = 10, value = 2, step= 2)

# Specification of range within an interval
legend.text.size <- sliderInput("legend.text.size", "Legend text Size :", min = 1, max = 20, value = 8)

# Provide a custom currency format for value display, 
# with basic animation
x.angle <- sliderInput("x.angle", "Angle of sample label:", min = 1, max = 90, value = 90, step = 1, animate=TRUE)

axis.text.size  <- sliderInput("axis.text.size", "Text size of variant labels (e.g., 3:3D):", min = 1, max = 20, value = 10)

ins.size <- sliderInput("ins.size", "Size of insertion symbols: ", min = 1, max = 10, value = 3)
legend.symbol.size <- sliderInput("legend.symbol.size", "legendsymbol size :", min = 1, max = 20, value = 5)
legend.text.size <- sliderInput("legend.text.size", "Size of insertion legend:", min = 1, max = 20, value = 8)

plot.type <- radioButtons("plot.type", "plot type", choices = c("proportions","counts"), selected = "counts", inline = TRUE) 

bscollapse_1 <- bsCollapse(id = "bscollapse_1", open = "Plot Alignments Options",
  bsCollapsePanel("Plot Alignments Options",
    p("Plots update automatically when parameters are changed"),
    hr(),
    #helpText("Text sizes"),
    axis.text.size,
    plot.text.size,
    legend.text.size,
    #helpText("Insertion symbols paramerters"),
    hr(),
    ins.size,

    style = "info"),
  bsCollapsePanel("HeatMap Options", 
    helpText("HeatMap Options"),
    x.angle,
    plot.type,
    style = "success" ))


plotOptions <- box(width = 4, 
  solidHeader = T,
  bscollapse_1,
  fluidRow(
    column(width=4,
      data_setting),
    column(width=4,
      save_data ),
    column(width=4,
	reset_btn)
  )
)
