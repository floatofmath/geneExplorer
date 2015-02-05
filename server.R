shinyServer(function(input,output){
    fitted <- reactive({
        time <- ifelse(input$comparison=="global","global",input$contrast)
        fit(input$method,time)
    })
    parcor <- reactive({
        input$compute
        pc <- isolate(TEparcor(fitted(),input$comparison,input$gsort,input$lsort,input$ngenes,input$nloci,input$custom))
        pc
    })
    corgraph <- reactive({
        g <- makeGraph(parcor(),input$cutoff)
    })
    output$heatmap <- renderPlot({
        heatmap(isolate(fitted()),input$comparison,input$sort,input$toptags,strsplit(input$custom,',[[:blank:]]*')[[1]])
    })
    output$table <- renderTable({
        tt <- TEtable(isolate(fitted()),input$comparison,input$sort,input$tabtags,strsplit(input$custom,',[[:blank:]]*')[[1]])
        locus <- tt[,paste0(Chromosome,":",start,"-",end)]
        tt[,NAMES:=makeUCSClink(locus,link=NAMES,project="&hgsid=202053665_L5FcliC2ZMly7BZLX6qdayCblg1D")]
        tt
    },sanitize.text.function = function(x) x)
    output$graph <- renderPlot({
        input$plot
        isolate(TEplotGraph(corgraph()))
    })
    output$graphsize <- renderText({
        gs <- graphSize(corgraph())
        paste("A graph with",gs[1],"nodes and",gs[2],"edges, will be drawn.")
    })
    output$corsize <- renderText({
        ngenes <- ifelse(input$gsort == "all",length(grep("NM",rownames(resids))),input$ngenes)
        if(input$gsort=="custom") ngenes <- length(strsplit(input$custom,',[[:blank:]]*')[[1]])
        nloci <- ifelse(input$lsort ==  "all",length(grep("locus",rownames(resids))),input$nloci)
        paste("A correlation matrix with",ngenes,"rows and",nloci,"columns, will be computed.")
    })
})
