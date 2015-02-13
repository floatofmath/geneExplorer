library(bt88.03.704)
library(limma)
library(edgeR)
library(WGCNA)
library(reshape2)
library(igraph)
#library(DT)
load('normalized.Rd')
load('annotation.Rd')
source('elements.R')
source('functions.R')


## download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ensGene.txt.gz","ensGene.mm10.txt.gz")
## system("gunzip ensGene.mm10.txt.gz -c > ensGene.mm10.txt")
## annotation <- fread("ensGene.mm10.txt")[,c(2:6,9,


MM.global <- model.matrix(~0+Group, data=design)
MM.diet <- model.matrix(~0+Diet+Tissue,data=design)
## annotation <- annotation[!duplicated(NAMES),]

for(i in levels(design$Tissue)){
    design[[i]] <- (design$Tissue == i)
}


fit <- function(method,comparison,contrast,interaction){
    fit <- switch(method,
                  "limma"= fit.limma(dge.voom,models(comparison,contrast),contrasts(comparison,contrast,interaction))[filter.mad,],
                  "edgeR"= fit.edger(dge.edger,models(comparison,contrast),contrasts(comparison,contrast,interaction)))
    return(fit)
}


contrasts <- function(comparison,tissue,interaction){
    model <- models(comparison,tissue)
    if(comparison == "global"){
        cm <- do.call(makeContrasts,c(list(paste(colnames(model)[c(T,F)],colnames(model)[c(F,T)],sep='-')),list(levels=model)))
    } else if(comparison=="diet"){
        cm <- makeContrasts('HFD-Lean',levels=model)
    } else if(comparison=="tissue"){
        if(interaction){
            cm <- makeContrasts(contrasts=tissue,levels=model)
        } else {
            cm <- makeContrasts(contrasts=paste0("Lean",tissue),levels=model)
        }
    }
    return(cm)
}

models <- function(comparison,tissue=NULL){
    if(!(comparison %in% c("global","diet","tissue"))) stop(paste("Invalid comparison:",comparison,"selected"))
    if(!(is.null(tissue) || (tissue %in% levels(design$Tissue)))) stop(paste("Invalid tissue:",tissue,"selected"))
    if(comparison=="global"){
        mm <- MM.global
    } else if(comparison=="diet"){
        mm <- MM.diet
    } else if(comparison=="tissue"){
        if(is.null(tissue)) stop("To compare tissues, tissue has to be selected")
        mm <- model.matrix(as.formula(paste0("~Diet*",tissue)),data=design)
    }
    colnames(mm) <- sub("\\(Intercept\\)","Intercept",colnames(mm))
    colnames(mm) <- gsub("Tissue|Diet|Group|TRUE|:","",colnames(mm))
    return(mm)
}

    




