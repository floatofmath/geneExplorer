shinyUI(fluidPage(
    titlePanel("Gene explorer"),
    sidebarLayout(
    sidebarPanel(
        selectInput("method",
                    label="Which test to compute:",
                    choices=c("limma"="limma","edgeR"="edgeR"),
                    selected="edgeR"),
        selectInput("sort",
                    label="Which genes to show:",
                    choices=c("top"="top","custom"="custom"),
                    selected="top"),
        selectInput("comparison",
                    label="Which comparison",
                    choices=c("global"="global","Diet"="diet","Tissue"="tissue"),
                    selected="global"),
        conditionalPanel(
            condition="input.comparison == 'tissue'",
            selectInput("contrast",
                        label="Which tissue to compare",
                        choices=c("BAT","interintestinal","perigonadal","retroperitoneal","sc"),
                        selected="BAT"),
            checkboxInput("interaction",
                          label="Interaction (i.e. low fat diet works different it Tissue)")
        ),
        conditionalPanel(
            condition="input.sort == 'custom' || input.gsort == 'custom'",
            HTML("Enter Ensemble Transcript Identifiers separated by ', ':<br/><textarea id='custom' rows=10 cols=40>ENSMUST00000023116, ENSMUST00000037324, ENSMUST00000001455, ENSMUST00000109986, ENSMUST00000151266, ENSMUST00000112498, ENSMUST00000152710, ENSMUST00000169282, ENSMUST00000069949, ENSMUST00000102743</textarea><br/>")),
        actionButton("redraw","Apply changes and update graph")
    ),
    mainPanel(
        tabsetPanel(
            heatmap_output(),
            gene_table()
        )
    )
    )
))
        
    
    
