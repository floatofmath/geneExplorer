method_changer <- function(...){
    selectInput("method",
                label="Test Procedure",
                choices=c("limma"="limma","edgeR"="edgeR"),
                selected="edgeR")
}


heatmap_output <- function(...){
    tabPanel(title="Heatmap",
                    conditionalPanel(
                        condition="input.sort == 'top'",
                        sliderInput("toptags",
                                    label="Number of top tags to show",
                                    min=10,max=500,step=10,value=100)
                    ),
                    plotOutput("heatmap",height='1080px'))
}
correlation_graph <- function(...){
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
             plotOutput("graph",height="1080px"))
}
gene_table <- function(...){
    tabPanel("Table",
                     ## conditionalPanel(
                     ##     condition="input.method != 'c'",
                     ##     sliderInput("tabtags",
                     ##                 label="Number of top tags to show",
                     ##                 min=100,max=1500,step=200,value=100)
                     ## ),
             fluidRow(uiOutput("button")),
             dataTableOutput("table")
             )
}
