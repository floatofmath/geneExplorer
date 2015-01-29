shinyUI(fluidPage(
    titlePanel("Gene explorer"),
    sidebarLayout(
    sidebarPanel(
        selectInput("method",
                    label="Test Procedure",
                    choices=c("limma"="l","edgeR"="e","custom"="c"),
                    selected="e"),
        conditionalPanel(
            condition="input.method != 'c'",
            selectInput("contrast",
                        label="Comparison",
                        choices=c("global"="global","timepoint"="time","maxT"="maxT"),
                        selected="global")
            ),
        conditionalPanel(
            condition="input.contrast == 'time'",
            sliderInput("timepoint",
                        label="Timepoint for comparison",
                        min=24,max=96,step=24,value=48)),
        conditionalPanel(
            condition="input.method == 'c' || input.selectGenes == 'custom'",
            HTML("Enter NM Identifiers separated by ', ':<br/><textarea id='custom' rows=10 cols=40>NM_027498, NM_144885, NM_019979, NM_138314, NM_021468, NM_178589, NM_011766, NM_010846, NM_001164793, NM_001164155, NM_029645, NM_001109914, NM_175686, NM_001044384, NM_023503, NM_001146690, NM_001030307, NM_023449, NM_001285812, NM_145472, NM_001177856, NM_172122, NM_019406, NM_024434, NM_001276482, NM_001102410, NM_001081164, NM_016957, NM_026029, NM_176838</textarea><br/>"))
    ),
    mainPanel(
        tabsetPanel(
            tabPanel(title="Heatmap",
                     conditionalPanel(
                         condition="input.method != 'c'",
                         sliderInput("toptags",
                                     label="Number of top tags to show",
                                     min=10,max=500,step=10,value=100)
                     )), #plotOutput("heatmap")),
            tabPanel("Correlation",
                     fluidRow(
                         column(4,
                                selectInput("selectGenes",
                                            label="Select genes to correlate",
                                            choices=c("top","all","custom"),
                                            selected="top"),
                                conditionalPanel(
                                    condition="input.selectGenes == 'top'",
                                    sliderInput("toptags",
                                                label="Number of top tags to show",
                                                min=10,max=500,step=10,value=100))),
                         column(4,
                                selectInput("selectLoci",
                                            label="Select loci to correlate",
                                            choices=c("top","all"),
                                            selected="top"),
                                conditionalPanel(
                                    condition="input.selectLoci == 'top'",
                                    sliderInput("toptags",
                                                label="Number of top tags to show",
                                                min=10,max=500,step=10,value=100)))
                     )
                     ), #plotOutput("corgraph")),
            tabPanel("Table",
                     conditionalPanel(
                         condition="input.method != 'c'",
                         sliderInput("toptags",
                                     label="Number of top tags to show",
                                     min=100,max=1500,step=200,value=100)
                     ) ##tableOutput("table"))
                     )
        )
    )
    )
))
        
    
    
