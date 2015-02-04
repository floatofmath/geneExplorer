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
            HTML("Enter NM Identifiers separated by ', ':<br/><textarea id='custom' rows=10 cols=40>NM_026377, NM_001003719, NM_026454, locus.2725.3, NM_001040072, NM_007953, locus.1073.29, NM_146224, locus.5293.7, locus.1785.79, NM_001281829, NM_172310, locus.3971.1, NM_023363, locus.317.243, locus.6277.8, locus.10411.34, NM_001271500, locus.3938.3, NM_001167905</textarea><br/>"))
    ),
    mainPanel(
        tabsetPanel(
            heatmap_output(),
            correlation_graph(),
            gene_table()
        )
    )
    )
))
        
    
    
