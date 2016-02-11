data_setting <- bsButton("data_setting", 'Setting', icon =  icon("list-alt"), style = "primary", block = TRUE)
info_btn <- bsButton("info_btn", 'Help', icon =  icon("fa fa-pencil"), style = "primary", block = TRUE)
save_data <- bsButton("save_data", 'Save Data', icon =  icon("fa fa-pencil"), style = "info", block = TRUE)

#----------------
# Slider
#----------------

# Simple integer interval
x.size <- sliderInput("x.size", "Line weight:", min=0, max=20, value=10, step = 1)

# Decimal interval with step value
plot.text.size <- sliderInput("plot.text.size", "Symbol size :", 
  min = 0, max = 10, value = 4, step= 2)

# Specification of range within an interval
legend.text.size <- sliderInput("legend.text.size", "Legend Symbol Size :",
  min = 1, max = 20, value = 8)

# Provide a custom currency format for value display, 
# with basic animation
x.angle <- sliderInput("x.angle", "X angle :", 
  min = 1, max = 180, value = 45, step = 1, animate=TRUE)

# Animation with custom interval (in ms) to control speed,
# plus looping
axis.text.size  <- sliderInput("axis.text.size", "x-axis size :", min = 1, max = 20, value = 4)

ins.size <- sliderInput("ins.size", "x-axis size :", min = 1, max = 20, value = 10)
legend.symbol.size <- sliderInput("legend.symbol.size", "legendsymbol size :", min = 1, max = 20, value = 5)
legend.text.size <- sliderInput("legend.text.size", "legend text size :", min = 1, max = 20, value = 10)


bscollapse_1 <- bsCollapse(id = "bscollapse_1", open = "Plot Options",
  bsCollapsePanel("Plot Options",
    helpText("Plots update automatically when parameters are changed"),
    x.size,
    plot.text.size,
    legend.text.size,
    x.angle,
    style = "info"),
  bsCollapsePanel("HeatMap Options", 
    p("HeatMap Options"),
    axis.text.size,
    ins.size,
    legend.symbol.size,
    legend.text.size, 
    style = "success" ))


plotOptions <- box(width = 4, 
  solidHeader = T,
  bscollapse_1,
  fluidRow(
    column(width=4,
      info_btn),
    column(width=4,
      data_setting),
    column(width=4,
      save_data )
  )
)
