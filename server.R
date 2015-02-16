shinyServer(function(input,output){
    customTags <- reactive({
        strsplit(input$custom,',[[:blank:]]*\n*')[[1]]
    })
    fitted <- reactive({
        time <- ifelse(input$comparison=="global","global",input$contrast)
        fit(input$method,time)
    })
    parcor <- reactive({
        input$compute
        pc <- isolate(TEparcor(fitted(),input$comparison,input$gsort,input$lsort,input$ngenes,input$nloci,customTags()))
        pc
    })
    corgraph <- reactive({
        g <- makeGraph(parcor(),input$cutoff)
    })
    output$heatmap <- renderPlot({
        heatmap(isolate(fitted()),input$comparison,input$sort,input$toptags,customTags())
    })
    output$table <- renderDataTable({
<<<<<<< HEAD
        tt <- TEtable(isolate(fitted()),input$comparison,input$sort,Inf,customTags())
        locus <- tt[,paste0(Chromosome,":",start,"-",end)]
        tt[,NAMES:=makeUCSClink(locus,link=NAMES,project="&hgsid=202053665_L5FcliC2ZMly7BZLX6qdayCblg1D")]
=======
        input$redraw
        tt <- TEtable(isolate(fitted()),input$comparison,input$sort,Inf,strsplit(input$custom,',[[:blank:]]*')[[1]])
        ##locus <- tt[,paste0(Chromosome,":",start,"-",end)]
        tt[,locus]
        tt[,nearest_ref_id:=makeUCSClink(locus,link=nearest_ref_id,project="&hgsid=202468101_ujaJjbSojoAWpGMMqKBM51YXk2QY")]
>>>>>>> 650514b... added ranks to toptable, simplified a few calls by autodetecting method
        tt
    },escape=FALSE)
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
        if(input$gsort=="custom") ngenes <- length(customTags())
        nloci <- ifelse(input$lsort ==  "all",length(grep("locus",rownames(resids))),input$nloci)
        paste("A correlation matrix with",ngenes,"rows and",nloci,"columns, will be computed.")
    })
})
