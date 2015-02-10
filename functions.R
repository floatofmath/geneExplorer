TEtopTags.limma.global <- function(fit,n,pattern){
    if(!is.null(pattern)) fit <- fit[grepl(pattern,rownames(fit)),]
    if(is.infinite(n)){ n <- nrow(fit) }
    rownames(fit)[order(fit$F.p.value)][1:n]
}
TEtopTags.limma.contrast <- function(fit,n,pattern){
    if(!is.null(pattern)) fit <- fit[grepl(pattern,rownames(fit)),]
    if(is.infinite(n)){ n <- nrow(fit) }
    rownames(fit)[order(fit$p.value)][1:n]
}
TEtopTags.edgeR <- function(comparison,fit,n,pattern){
    if(!is.null(pattern)) fit <- fit[grepl(pattern,rownames(fit)),]
    if(is.infinite(n)){ n <- nrow(fit) }
    rownames(fit)[order(fit$table$PValue)][1:n]
}
TEtopTags.limma <- function(comparison,fit,n,pattern){
    if(comparison=="global"){
        TEtopTags.limma.global(fit,n,pattern)
    } else {
        TEtopTags.limma.contrast(fit,n,pattern)
    }
}
TEtopTags <- function(method,comparison,fit,n,pattern=NULL){
    switch(method,
           "limma"=TEtopTags.limma,
           "edgeR"=TEtopTags.edgeR)(comparison,fit,n,pattern)
}
TEtopStats.limma.global <- function(fit,toptags){
    fit[toptags,]$F
}
TEtopStats.limma.contrast <- function(fit,toptags){
    as.numeric(coefficients(fit[toptags,]))
}
TEtopStats.edgeR.global <- function(fit,toptags){
    fit[toptags,]$table$LR
}
TEtopStats.edgeR.contrast <- function(fit,toptags){
    fit[toptags,]$table$logFC
}
TEtopStats.edgeR <- function(comparison,fit,toptags){
    if(comparison=="global"){
        TEtopStats.edgeR.global(fit,toptags)
    } else {
        TEtopStats.edgeR.contrast(fit,toptags)
    }
}
TEtopStats.limma <- function(comparison,fit,toptags){
    if(comparison=="global"){
        TEtopStats.limma.global(fit,toptags)
    } else {
        TEtopStats.limma.contrast(fit,toptags)
    }
}

TEtopStats <- function(method,comparison,fit,toptags){
    if(is.numeric(toptags) && length(toptags) == 1){
        toptags <- TEtopTags(method,comparison,fit,toptags)
    }
    switch(method,
           "limma"=TEtopStats.limma,
           "edgeR"=TEtopStats.edgeR)(comparison,fit,toptags)
}
topCounts <- function(toptags,method){
    switch(method,
           "limma"=dge.voom$E[toptags,],
           "edgeR"=log2(cpm(dge.edger)+1)[toptags,]
           )}
    
##' Selects tags for a given fit and displays a corresponding heatslide
##'
##' @title Print a heatslide for selected tags
##' @param fit model fit (either from edgeR or limma)
##' @param comparison which comparison to show stats for and select top genes from either "global" for a global test or "contrast" for a single contrast
##' @param sort use top tags or show a custom list of tags either "top" or "custom"
##' @param n number of to tags to show if \code{sort} is "top"
##' @param tags custom list of tags to show if \code{sort} is "custom"
##' @return 
##' @author float
heatmap <- function(fit,comparison,sort,n=NULL,tags=NULL){
    method <- c("DGELRT"="edgeR","MArrayLM"="limma")[class(fit)]
    caption <- ifelse(comparison=="global",
                      c("edgeR"="LR Statistic","limma"="F Statistik")[method],
                      "Log (Base 2) Foldchange")
    if(sort == "custom"){
        idx <- tags
    } else {
        idx <- TEtopTags(method,comparison,fit,n)
    }
    mat <- topCounts(idx,method)
    stat <- TEtopStats(method,comparison,fit,idx)
    gn <- annotation[idx,][["Gene Name"]]
    gn[is.na(gn)] <- idx[is.na(gn)]
    heatslide(mat,stat,as.factor(design$Group),genenames=gn,slidetitle=caption)
}


TEtopTable.edgeR <- function(fit,tags){
    fit <- fit[tags,]
    myTopTags(fit,Inf,annotation)
}

TEtopTable.limma <- function(fit,tags){
    fit <- fit[tags,]
    tt <- as.data.table(topTable(fit,number=Inf),keep.rownames=T)
    annotation[tt,]
}

    
TEtopTable <- function(fit,tags){
    method <- c("DGELRT"="edgeR","MArrayLM"="limma")[class(fit)]
    switch(method,
           "edgeR" = TEtopTable.edgeR,
           "limma" = TEtopTable.limma)(fit,tags)
}
           


TEtable <- function(fit,comparison,sort,n=NULL,tags=NULL){
    method <- c("DGELRT"="edgeR","MArrayLM"="limma")[class(fit)]
    statistic <- ifelse(comparison=="global",
                      c("edgeR"="LR-Statistic","limma"="F-Statistic")[method],
                      "Log2 Foldchange")
    if(sort == "custom"){
        idx <- tags
    } else {
        idx <- TEtopTags(method,comparison,fit,n)
    }
    counts <- as.data.table(topCounts(idx,method),keep.rownames=TRUE)
    stat <- TEtopTable(fit,idx)
    setkey(stat,nearest_ref_id)
    tab <- stat[counts,]
    tab
}


TEparcor <- function(fit,comparison,gsort,lsort,ngenes=NULL,nloci=NULL,genes=NULL){
    method <- c("DGELRT"="edgeR","MArrayLM"="limma")[class(fit)]
    if(lsort == "all"){
        nloci <- length(grep("locus",rownames(fit),perl=T))
        lsort <- "top"
    }
    loci <- TEtopTags(method,comparison,fit,nloci,pattern="locus")
    if(gsort == "all"){
        ngenes <- length(grep("^NM",rownames(fit),perl=T))
        gsort <- "top"
    }
    if(gsort == "custom") {
        genes <- genes
    } else {
        genes <- TEtopTags(method,comparison,fit,ngenes,pattern="^NM")
    }
    gres <- resids[genes,]
    lres <- resids[loci,]
    parCor <- bicor(t(gres),t(lres))
    parCor <- as.data.table(parCor,keep.rownames=TRUE)
    setnames(parCor,'rn','geneNAMES')
    parCor <- melt(parCor,id.vars='geneNAMES',value.name="parCor",variable.name="lncNAMES")
    setkey(parCor,'geneNAMES')
    setkey(annotation,"NAMES")
    parCor <- annotation[parCor,][,c("NAMES","Gene Name","lncNAMES","parCor"),with=F]
    setnames(parCor,'NAMES','geneNAMES')
    setnames(parCor,'Gene Name','gene_name')
    parCor[is.na(gene_name),gene_name:=geneNAMES]
#    parCor[,parCor:=abs(parCor)]
    return(parCor)
}

graphSize <- function(g){
  c(length(V(g)),length(E(g)))
}

makeGraph <- function(parCor,cut){
    parCor[,edge:=abs(parCor)>cut]
#    edgeList.nm <- as.matrix(parCor[edge==T,list(geneNAMES,lncNAMES)])
    edgeList <- as.matrix(parCor[edge>0,list(gene_name,lncNAMES,abs(parCor))])
    dups <- paste0(edgeList[,1],edgeList[,2])
    weights <- edgeList[!duplicated(dups),3]
    edgeList <- edgeList[!duplicated(dups),1:2]
##    out <- dcast(parCor[edge==T,],geneNAMES~lncNAMES,value.var='edge')
    g <- graph.edgelist(edgeList,directed=F)
    E(g)$weight <- weights
    return(g)
}

TEplotGraph <- function(g,...){
    plot(g,layout.kamada.kawai(g),
         vertex.size=2,
         edge.width=4*(as.numeric(E(g)$weight))^10,
         vertex.label.dist=.2,...)
}
