shinyUI(fluidPage(
    titlePanel("Gene explorer"),
    sidebarLayout(
    sidebarPanel(
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
            HTML("Enter NM Identifiers separated by ', ':<br/><textarea id='custom' rows=10 cols=40>NM_026377, NM_001003719, NM_026454, NM_001040072, NM_007953</textarea><br/>"))
    ),
    mainPanel(
        tabsetPanel(
            heatmap_output(),
            gene_table()
        )
    )
    )
))
        
    
    
