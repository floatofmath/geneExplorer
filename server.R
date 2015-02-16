shinyServer(function(input,output){
    customTags <- reactive({
        strsplit(input$custom,',[[:blank:]]*\n*')[[1]]
    })
    fitted <- reactive({
        ## this is governed by a UI in other apps
        fit(input$method,input$comparison,input$contrast,input$interaction)
    })
    ## parcor <- reactive({
    ##     input$compute
    ##     pc <- isolate(TEparcor(fitted(),input$comparison,input$gsort,input$lsort,input$ngenes,input$nloci,input$custom))
    ##     pc
    ## })
    ## corgraph <- reactive({
    ##     g <- makeGraph(parcor(),input$cutoff)
    ## })
    output$heatmap <- renderPlot({
        input$redraw
        heatmap(isolate(fitted()),input$comparison,input$sort,input$toptags,strsplit(input$custom,',[[:blank:]]*')[[1]])
    })
    output$table <- renderDataTable({
        input$redraw
        tt <- TEtable(isolate(fitted()),input$comparison,input$sort,Inf,strsplit(input$custom,',[[:blank:]]*')[[1]])
        ##locus <- tt[,paste0(Chromosome,":",start,"-",end)]
        tt[,locus]
        tt[,nearest_ref_id:=makeUCSClink(locus,link=nearest_ref_id,project="&hgsid=202468101_ujaJjbSojoAWpGMMqKBM51YXk2QY")]
        tt
    },escape=F)
    output$downloadTable <- downloadHandler(
        filename=function() {           
            paste0('toptags-',Sys.Date(),'.csv')
        },
        content = function(con) {
            tt <- TEtable(isolate(fitted()),input$comparison,input$sort,Inf,customTags())
            write.csv(tt,con)
        },
        contentType="text/csv")
    output$button <- renderUI({
        downloadButton('downloadTable',paste("Download Spreadsheet for",            ifelse(input$sort=="custom",length(customTags()),"all"),"Top Tags"))})
    ## output$graph <- renderPlot({
    ##     input$plot
    ##     isolate(TEplotGraph(corgraph()))
    ## })
    ## output$graphsize <- renderText({
    ##     gs <- graphSize(corgraph())
    ##     paste("A graph with",gs[1],"nodes and",gs[2],"edges, will be drawn.")
    ## })
    ## output$corsize <- renderText({
    ##     ngenes <- ifelse(input$gsort == "all",length(grep("NM",rownames(resids))),input$ngenes)
    ##     if(input$gsort=="custom") ngenes <- length(strsplit(input$custom,',[[:blank:]]*')[[1]])
    ##     nloci <- ifelse(input$lsort ==  "all",length(grep("locus",rownames(resids))),input$nloci)
    ##     paste("A correlation matrix with",ngenes,"rows and",nloci,"columns, will be computed.")
    ## })
})
