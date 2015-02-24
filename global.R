library(bt88.03.704)
library(limma)
library(edgeR)
library(WGCNA)
library(reshape2)
library(igraph)
source(functions.R)
load('design.rda')
load('targets.Rd')
load('normalized.Rd')
load('annotation.Rd')

MM <- model.matrix(~Group, data=targets)
colnames(MM) <- sub("Group","",colnames(MM))
colnames(MM)[1] <- "Intercept"
C <- as.list(with(targets,unique(paste(Group[Treatment=='phx'],Group[Treatment=="ehx"],sep="-"))))
C$levels <- MM

m <- lmFit(dge.voom,MM)
em <- eBayes(m)#,trend=T)
resids <- residuals(em,dge.voom)
annotation <- annotation[!duplicated(NAMES),]


contrasts <- function(contrast){
    comparison <- as.character(contrast)
    contrast <- c("global"=0,"24" = 1,"32"=2,"48" = 3,"72" = 4,"96" = 5)[comparison]
    if(contrast==0){
        return(do.call(makeContrasts,C))
    } else {
        return(makeContrasts(contrasts=C[[contrast]],levels=C$levels))
    }
}

fit <- function(method,contrast){
    fit <- switch(method,
                  "limma"= fit.limma(dge.voom,MM,contrasts(contrast)),
                  "edgeR"= fit.edger(dge.edger,MM,contrasts(contrast)))

    return(fit)
}
