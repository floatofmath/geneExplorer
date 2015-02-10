<!--pandoc
t: html
template: report.html
css: vignette.css
self-contained: 
-->
---
title: RNA Sequencing of lean and fat mice 
author: Florian Klinglmueller
date: 
client: 
project: 
number: 
review: 
email: florian.klinglmueller@meduniwien.ac.at 
...
```{r,Libraries}
library(data.table)
library(ggplot2)
library(plyr)
library(reshape2)
library(bt88.03.704)
library(edgeR)
library(rvest)
```

## Section 1

Read in data from website and save to data.table

```{r,DownloadReads,eval=F}

# http://www.biomedical-sequencing.at/projects/BSA_0013_Fatpads_84573b2f553f463087419b4270950417/mm10/rnaseq_report.html
# http://www.biomedical-sequencing.at/projects/BSA_0013_Fatpads_84573b2f553f463087419b4270950417/mm10/rnaseq_cufflinks_MUW_0000_B011BABXX_4_ArmpitHindleg_wtblk6lean_wt/rnaseq_cufflinks_MUW_0000_B011BABXX_4_ArmpitHindleg_wtblk6lean_wt_isoforms_fpkm_tracking.tsv
url <- "http://www.biomedical-sequencing.at/projects/BSA_0013_Fatpads_84573b2f553f463087419b4270950417/mm10/"
data_page <- html(paste0(url,"rnaseq_report.html"))

links <- data_page %>% html_nodes("table a") %>% html_attr(name="href")

reads <- list()
## trans <- list()
for(i in links[grep('isoforms_fpkm_tracking.tsv',links)]){
    reads[[i]] <- fread(paste0(url,i),sep='auto')
}
## for(i in links[grep('isoforms\\.fpkm_tracking',links)]){
##     trans[[i]] <- fread(paste0(url,i),sep='auto')
## }
## Set up design matrix
design <- do.call('rbind',strsplit(gsub("\\./\\w+/([[:alnum:]]+_){5}([[:digit:]])_([[:alnum:]]+)_([[:alnum:]]+)_\\w+.tsv","\\2 \\3 \\4",names(reads)),' '))
colnames(design) <- c('ID','FatType','Diet')
design <- as.data.frame(design)[,c(2,3,1)]
names(reads) <- apply(design,1,paste,collapse="_") # nice names
samples <- names(reads)
iso_counts <- fread("http://www.biomedical-sequencing.at/projects/BSA_0013_Fatpads_84573b2f553f463087419b4270950417/mm10/rnaseq_process_cuffdiff_global/rnaseq_process_cuffdiff_global_isoforms_counts_replicates.tsv")
gene_counts <- fread("http://www.biomedical-sequencing.at/projects/BSA_0013_Fatpads_84573b2f553f463087419b4270950417/mm10/rnaseq_process_cuffdiff_global/rnaseq_process_cuffdiff_global_genes_counts_replicates.tsv")
samples <- colnames(gene_counts)[4:23]
design <- do.call('rbind',lapply(strsplit(samples,"_"),tail,3))
save(reads,design,samples,file='data.Rd')
save(gene_counts,iso_counts,design,samples,file="counts.Rd")

```


## Preprocess reads

We put the FPKMs into a datatable, each sample in one column. We
extract FPKMs for those transcripts that have an *Ensemble Transcript
ID* and discard the rest.

```{r,PreprocessReads,eval=F}
load('data.Rd')
## Some summary stats
m <- length(reads) # number of samples
N <- sapply(reads,nrow) # annotated reads per sample
Nu <- sapply(reads,function(x) length(unique(x$ensembl_transcript_id)))
duplicated_reads <- any(N-Nu>0)

## Data cleaning for the moment let's go with transcript ids
reads <- lapply(reads,setkey,"ensembl_transcript_id")

#reads.long <- rbindlist(lapply(1:length(reads),function(i) reads[[i]][,c("FatType","Diet","ID","Sample"):=list(design$FatType[i],design$Diet[i],design$ID[i],samples[i])]))
                        
reads.clean <- lapply(reads,`[`,j=c("ensembl_transcript_id","FPKM"),with=FALSE)
for(i in 1:length(reads.clean)) setnames(reads.clean[[i]],2,names(reads.clean)[i])


## Make one big data.table with samples in columns
fpkms <- Reduce(function(a,b) merge(a,b,all=TRUE,by="ensembl_transcript_id"),reads.clean)
rm(reads.clean)


## remove everything that does not even have an ENSMUST id
fpkms.known <- fpkms[grepl("^ENSMUST",ensembl_transcript_id),]


                                        #reads.long.known <- reads.long[grepl("^ENSMUST",ensembl_transcript_id),]

## breaks it down to about 90k known transcripts
known <- nrow(fpkms.known)
#nrow(reads.long.known)
levels(design$ID) <- c(0,0,1,1)
save(fpkms.known,m,N,design,samples,file='preproc.rda')
```

```{r preprocess counts}
load('counts.Rd')
counts.known <- iso_counts[grep("ENSMUST",nearest_ref_id,perl=T),]
known <- nrow(counts.known)
m <- length(samples)
colnames(design) <- c('Tissue','Diet','Sample')
design <- as.data.frame(design)
design <- within(design,{ Lane=1:20
                          Label=samples
                          Group=paste(Tissue,Diet,sep='_')})
                 

counts.unique <- counts.known[,lapply(.SD,mean),by=nearest_ref_id,.SDcols=samples]
setkey(counts.known,nearest_ref_id)
setkey(counts.unique,nearest_ref_id)
## rejoin annotation data
counts.unique <- counts.known[,1:5,with=F][counts.unique,,mult='first']

save(counts.unique,file='unique.Rd')
DF <- as.data.frame(counts.unique[,samples,with=F])
rownames(DF) <- counts.unique[,nearest_ref_id]

dge <- DGEList(counts=DF,group=design$Group) #, group=with(targets,factor(Treatment):factor(Time)))w
rm(DF)

colnames(dge) <- design$Label
dim(dge)	



# We filter out very lowly expressed tags cpm>1 in at least 3 samples (as groups consists of 3 samples)
filter.thresh <- rowSums(cpm(dge)>1.5) >= 4
dge.thresh <- dge[filter.thresh,]
dim(dge.thresh)
filter.mad <- rank(-apply(cpm(dge.thresh),1,mad))<=8000
dge.mad <- dge.thresh[filter.mad,]

dge.mad$samples$lib.size <- colSums(dge.mad$counts)
dge.mad.tmm <- calcNormFactors(dge.mad)
dge.mad.tmm$samples
save(dge.mad.tmm,samples,design,file='tmm.Rd')
```

```{r plot distributions,eval=F}
long <- as.data.table(dge.mad.tmm$counts)
long[,NAMES:=rownames(dge.mad)]
long <- melt(long,id.vars='NAMES')
pdf('dist.pdf')
qplot(variable,log2(value),data=long,geom='boxplot')
dev.off()

pdf('mds.pdf')
plotMDS(dge.mad.tmm,top=Inf)
dev.off()

```
There are `r m` samples, with `r known` reads that have an *Ensemble Transcript ID*. 

## Filtering


```{r edger}
load('tmm.Rd')
MM.global <- model.matrix(~0+Group, data=design)
MM.diet <- model.matrix(~0+Diet+Tissue,data=design)

colnames(MM.diet) <- sub("Tissue|Diet|Group","",colnames(MM.diet))
colnames(MM.global) <- sub("Tissue|Diet|Group","",colnames(MM.global))

design <- within(design,BAT <- factor(Tissue=='BAT',labels=c('other','bat')))
MM.bat <- model.matrix(~Diet*BAT,data=design)
colnames(MM.bat) <- gsub("Tissue|Diet|Group|BAT","",colnames(MM.bat))

                                        #colnames(MM)[1] <- "Intercept"

dge.edger <- estimateDisp(dge.mad.tmm,prior.df=10) 
save(dge.edger,file="normalized.Rd")
save(design,file="design.Rd")

                                        # Disp = 0.15
dge.edger$common.dispersion
summary(dge.edger$tagwise.dispersion)

CM.global <- do.call(makeContrasts,c(list(paste(colnames(MM.global)[c(T,F)],colnames(MM.global)[c(F,T)],sep='-')),list(levels=MM.global)))
CM.diet <- makeContrasts('HFD-Lean',levels=MM.diet)




edger.fvl <- fit.edger(dge.edger,MM.diet,CM.diet)
edger.global <- fit.edger(dge.edger,MM.global,CM.global)
edger.bat <- fit.edger(dge.edger,MM.bat,NULL,coef=3)
edger.dietbat <- fit.edger(dge.edger,MM.bat,NULL,coef=4)

tt.fvl <- myTopTags(edger.fvl,8000,counts.unique)
tt.global <- myTopTags(edger.global,8000,counts.unique)
tt.bat <- myTopTags(edger.bat,8000,counts.unique)
tt.dietbat <- myTopTags(edger.dietbat,8000,counts.unique)

save(design,
    tt.fvl,
    tt.bat,
    tt.dietbat,
    tt.global,file='forApp.Rd')

tt <- myTopTags(edger.bat,30,counts.unique)
plotTopTags(tt,design$Tissue,11:30)
load('unique.Rd')
# Table with results

write.csv(myTopTags(edger.bat,1000,counts.unique),'toptable_hvl_BAT.csv')
write.csv(myTopTags(edger.fvl,1000,counts.unique),'toptable_hvl.csv')
write.csv(myTopTags(edger.global,1000,counts.unique),'toptable_global.csv')



makeUCSClink(tt,1)

```

























```{r, count Filtering}
## Additionally median normalization
long.raw <- melt(counts.known,id.vars=1:5)

## undo melt!
##wide.raw <- dcast.data.table(long.raw,ensembl_transcript_id~variable)

## Any NA
anyNA <- long.raw[,any(is.na(value)),by="nearest_ref_id"][,sum(V1)]
allNA <- long.raw[,all(is.na(value)),by="nearest_ref_id"][,sum(V1)]
## Any 0
any0 <- long.raw[!is.na(value),any(value==0),by="nearest_ref_id"][,sum(V1)]
## All > 0
allg0 <- long.raw[!is.na(value),all(value>0),by="nearest_ref_id"][,sum(V1)]
## Set NA to Zero
if(anyNA > 0) long.raw[is.na(value),value:=0]

## All > 0 !na
allpos <- long.raw[!is.na(value),all(value>0),by="nearest_ref_id"][,sum(V1)]

### Filtering
## At least 4 larger than too_low
too_low <- 1.5
at_least <- 4
long.fil <- long.raw[,if(sum(value>too_low,na.rm=T)>=at_least) .SD,by="nearest_ref_id"]


## MAD filtering
long.fil[,MAD:=mad(value),by="nearest_ref_id"]
pdf('dist.pdf')wide.fil <- dcast.data.table(long.fil,nearest_ref_id+MAD~variable,'mediadge.mad')
## 15000 transcripts above madold in at least 4 samples
filT <- nrow(wide.fil)

## select filT
setkeyv(wide.fil,c('MAD','nearest_ref_id'))
long.mad <- melt(wide.fil[filT:(filT-7999),],id.vars=c('nearest_ref_id','MAD'))
setkeyv(long.mad,c('MAD','ensembl_transcript_id','variable'))

```

```{r fpkm Filtering,eval=F}
load('preproc.rda')

## Additionally median normalization
long.raw <- melt(fpkms.known,id.vars="ensembl_transcript_id")

## undo melt!
##wide.raw <- dcast.data.table(long.raw,ensembl_transcript_id~variable)

## Any NA
anyNA <- long.raw[,any(is.na(value)),by="ensembl_transcript_id"][,sum(V1)]
allNA <- long.raw[,all(is.na(value)),by="ensembl_transcript_id"][,sum(V1)]
## Any 0
any0 <- long.raw[!is.na(value),any(value==0),by="ensembl_transcript_id"][,sum(V1)]
## All > 0
allg0 <- long.raw[!is.na(value),all(value>0),by="ensembl_transcript_id"][,sum(V1)]
## Set NA to Zero
long.raw[is.na(value),value:=0]

## All > 0 !na
allpos <- long.raw[!is.na(value),all(value>0),by="ensembl_transcript_id"][,sum(V1)]

### Filtering
## At least 4 larger than too_low
too_low <- 1.5
at_least <- 4
long.fil <- long.raw[,if(sum(value>too_low,na.rm=T)>=at_least) .SD,by="ensembl_transcript_id"]


## MAD filtering
long.fil[,MAD:=mad(value),by="ensembl_transcript_id"]
wide.fil <- dcast.data.table(long.fil,ensembl_transcript_id+MAD~variable)
## 15000 transcripts above threshold in at least 4 samples
filT <- nrow(wide.fil)

## select filT
setkeyv(wide.fil,c('MAD','ensembl_transcript_id'))
long.mad <- melt(wide.fil[filT:(filT-7999),],id.vars=c('ensembl_transcript_id','MAD'))
setkeyv(long.mad,c('MAD','ensembl_transcript_id','variable'))
```

```{r}


long.raw[,diet:=design[variable,'Diet']]
long.raw[,ID:=design[variable,'ID']]
long.raw[,fat:=design[variable,'FatType']]
long.fil[,diet:=design[variable,'Diet']]
long.fil[,ID:=design[variable,'ID']]
long.fil[,fat:=design[variable,'FatType']]
long.mad[,diet:=design[variable,'Diet']]
long.mad[,ID:=design[variable,'ID']]
long.mad[,fat:=design[variable,'FatType']]



a <- ggplot(long.raw,aes(fat:diet:ID,log2(value),fill=fat))+geom_boxplot()+opts(axis.text.x=theme_text(angle=-90))
b <- ggplot(long.fil,aes(fat:diet:ID,log2(value),fill=fat))+geom_boxplot()+opts(axis.text.x=theme_text(angle=-90))
c <- ggplot(long.mad,aes(fat:diet:ID,log2(value),fill=fat))+geom_boxplot()+opts(axis.text.x=theme_text(angle=-90))
pdf('distributions.pdf',width=20,height=9)
multiplot(a,b,c,cols=3,titles=c("All Transcripts","Above 1.5 in 4","Top 8000 MAD"))
dev.off()

```




```{r}
median.fpkm <-  long.raw[value>0,median(value)]
median.fpkm.sample <- long.raw[value>0,median(log2(value)),by=variable]
norm.fact <- median.fpkm.sample[,log2(V1)-median.fpkm]
setkey(median.fpkm.sample,Sample)


## median normalized
fpkms.known.mn <- fpkms[,(samples):=lapply(1:length(.SD),function(i) .SD[[i]]*2^norm.fact[i]),.SDcols=samples]

long <- melt(fpkms.known.mn,id.vars="ensembl_transcript_id")
wide <- dcast(long,variable~value)
#long[,c("fat","diet","ID"):=rbindlist(strsplit(as.character(variable),'_'))]

a <- ggplot(long,aes(x=variable,y=log2(value)),)+geom_boxplot()
b <- ggplot(long.raw,aes(x=variable,y=log2(value)),)+geom_boxplot()
multiplot(a,b)

```
