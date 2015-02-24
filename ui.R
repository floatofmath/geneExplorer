shinyUI(fluidPage(
    titlePanel("Gene explorer"),
    sidebarLayout(
    sidebarPanel(
        selectInput("method",
                    label="Test Procedure",
                    choices=c("limma"="limma","edgeR"="edgeR"),
                    selected="e"),
        selectInput("sort",
                    label="Which genes to show:",
                    choices=c("top"="top","custom"="custom"),
                    selected="top"),
        selectInput("comparison",
                    label="Which comparison",
                    choices=c("global"="global","timepoint"="contrast"),
                    selected="global"),
        conditionalPanel(
            condition="input.comparison == 'contrast'",
            sliderInput("contrast",
                        label="Timepoint for comparison",
                        min=24,max=96,step=24,value=48)),
        conditionalPanel(
            condition="input.sort == 'custom' || input.gsort == 'custom'",
            HTML("Enter NM Identifiers separated by ', ':<br/><textarea id='custom' rows=10 cols=40>NM_001001144, NM_001001176, NM_001001179, NM_001001181, NM_001001182, NM_001001183, NM_001001184, NM_001001185, NM_001001187, NM_001001295</textarea><br/>"))
    ),
    mainPanel(
        tabsetPanel(
            tabPanel(title="Heatmap",
                     conditionalPanel(
                         condition="input.sort == 'top'",
                         sliderInput("toptags",
                                     label="Number of top tags to show",
                                     min=10,max=500,step=10,value=100)
                     ),
                     plotOutput("heatmap",height='1080px')),
            tabPanel("Correlation",
                     fluidRow(
                         column(4,
                                selectInput("gsort",
                                            label="Select genes to correlate",
                                            choices=c("top","all","custom"),
                                            selected="top"),
                                conditionalPanel(
                                    condition="input.gsort == 'top'",
                                    sliderInput("ngenes",
                                                label="Number of top tags to show",
                                                min=10,max=1000,step=10,value=100)),
                                textOutput("corsize")),
                         column(4,
                                selectInput("lsort",
                                            label="Select loci to correlate",
                                            choices=c("top","all"),
                                            selected="top"),
                                conditionalPanel(
                                    condition="input.lsort == 'top'",
                                    sliderInput("nloci",
                                                label="Number of top tags to show",
                                                min=10,max=1000,step=10,value=100)),
                                actionButton("compute","Apply and recompute correlation matrix")),
                         column(4,
                                sliderInput("cutoff",
                                            label="Cutoff above which an edge is drawn",
                                            min=0,
                                            max=1,
                                            value=0.7,step=0.01),
                                textOutput("graphsize"),
                                actionButton("plot","Apply cutoff and plot graph"))
                     ),
                     plotOutput("graph",height="1080px")),
            tabPanel("Table",
                     fluidRow(uiOutput("button")),
                     ## conditionalPanel(
                     ##     condition="input.method != 'c'",
                     ##     sliderInput("tabtags",
                     ##                 label="Number of top tags to show",
                     ##                 min=100,max=1500,step=200,value=100)
                     ## ),
                     dataTableOutput("table")
                     )
        )
    )
    )
))
        
    
    
