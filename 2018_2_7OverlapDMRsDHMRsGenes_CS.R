
library(GenomicRanges)
library(org.Hs.eg.db)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(annotatr)

load("parietal.mid.dmr.RData")
load("parietal.late.dmr.RData")
load("cingulate.mid.dmr.RData")
load("cingulate.late.dmr.RData")

########Overlap to known genes ############

# First create a GRanges with the locations of known genes. 
# Download gene locations from org.Hs.eg
gene.names = mappedkeys(org.Hs.egSYMBOL)
gene.names=as.data.frame(org.Hs.egSYMBOL[gene.names])

chrbegin = mappedkeys(org.Hs.egCHRLOC)
chrbegin=as.data.frame(org.Hs.egCHRLOC[chrbegin])

# Only keep the genes that are named
xx=match(chrbegin$gene_id,gene.names$gene_id)
gene.names=as.character(gene.names$symbol[xx])

# Download the chr coord for each gene
chrs=as.factor(chrbegin$Chromosome)
chrbegin=abs(as.numeric(chrbegin$start_location))

chrend = mappedkeys(org.Hs.egCHRLOCEND)
chrend=as.data.frame(org.Hs.egCHRLOCEND[chrend])
chrend=abs(as.numeric(chrend$end_location))

xx=grep("_",chrs,invert=TRUE)
chrs=as.factor(as.character(chrs[xx]))
chrbegin=chrbegin[xx]
chrend=chrend[xx]
gene.names.v=gene.names[xx]
gene.names=data.frame(gene.names[xx])

# Create an GRanges of gene locations 
starts.ends=IRanges(start=chrbegin,end=chrend)
gene.locs=GRanges(seqnames=Rle(chrs),ranges=starts.ends,mcols=gene.names)


# Now overlap these gene locations to our ranges
overlp.gene.parietal.mid=findOverlaps(parietal.mid.dmr.ranges,gene.locs)
overlp.gene.parietal.late=findOverlaps(parietal.late.dmr.ranges,gene.locs)
overlp.gene.cingulate.mid=findOverlaps(cingulate.mid.dmr.ranges,gene.locs)
overlp.gene.cingulate.late=findOverlaps(cingulate.late.dmr.ranges,gene.locs)

# Do the same for the control negatives
overlp.gene.parietal.mid.neg=findOverlaps(parietal.mid.dmr.neg.ranges,gene.locs)
overlp.gene.parietal.late.neg=findOverlaps(parietal.late.dmr.neg.ranges,gene.locs)
overlp.gene.cingulate.mid.neg=findOverlaps(cingulate.mid.dmr.neg.ranges,gene.locs)
overlp.gene.cingulate.late.neg=findOverlaps(cingulate.late.dmr.neg.ranges,gene.locs)


# Add gene names to its own column over to the modifcation meta data frames
parietal.mid.meta=data.frame(parietal.mid.dmr.meta[overlp.gene.parietal.mid@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.parietal.mid@to],stringsAsFactors=F)
parietal.late.meta=data.frame(parietal.late.dmr.meta[overlp.gene.parietal.late@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.parietal.late@to],stringsAsFactors=F)
cingulate.mid.meta=data.frame(cingulate.mid.dmr.meta[overlp.gene.cingulate.mid@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.cingulate.mid@to],stringsAsFactors=F)
cingulate.late.meta=data.frame(cingulate.late.dmr.meta[overlp.gene.cingulate.late@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.cingulate.late@to],stringsAsFactors=F)

colnames(cingulate.mid.meta)=colnames(cingulate.late.meta)=colnames(parietal.mid.meta)=colnames(parietal.late.meta)=c("chr.coord","noCpGs","minFDR","meanBeta","gene")
colnames(parietal.mid.meta)=c("chr.coord","noCpGs","minFDR","meanBeta","gene")

save(list=c("cingulate.mid.meta","cingulate.late.meta","parietal.mid.meta","parietal.late.meta"),file=paste0("genes.RData"))

write.table(parietal.mid.meta,file="parietal.mid.meta.txt",sep="\t")
write.table(parietal.late.meta,file="parietal.late.meta.txt",sep="\t")
write.table(cingulate.late.meta,file="cingulate.mid.meta.txt",sep="\t")
write.table(cingulate.late.meta,file="cingulate.late.meta.txt",sep="\t")


#Make data frames for neg data
parietal.mid.meta.neg=data.frame(parietal.mid.dmr.neg.meta[overlp.gene.parietal.mid.neg@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.parietal.mid.neg@to],stringsAsFactors=F)
parietal.late.meta.neg=data.frame(parietal.late.dmr.neg.meta[overlp.gene.parietal.late.neg@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.parietal.late.neg@to],stringsAsFactors=F)
cingulate.mid.meta.neg=data.frame(cingulate.mid.dmr.neg.meta[overlp.gene.cingulate.mid.neg@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.cingulate.mid.neg@to],stringsAsFactors=F)
cingulate.late.meta.neg=data.frame(cingulate.late.dmr.neg.meta[overlp.gene.cingulate.late.neg@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.cingulate.late.neg@to],stringsAsFactors=F)

colnames(cingulate.mid.meta.neg)=colnames(cingulate.late.meta.neg)=colnames(parietal.mid.meta.neg)=colnames(parietal.late.meta.neg)=c("chr.coord","noCpGs","medFDR","meanBeta","gene")

colnames(parietal.mid.meta.neg)=c("chr.coord","noCpGs","medFDR","meanBeta","gene")





pd.ranges=c("parietal.mid.dmr.ranges","parietal.late.dmr.ranges","cingulate.mid.dmr.ranges","cingulate.late.dmr.ranges","parietal.mid.dmr.neg.ranges","parietal.late.dmr.neg.ranges","cingulate.mid.dmr.neg.ranges","cingulate.late.dmr.neg.ranges")

# Use annotatr to annotate intergenic and/or regulatory regions
annot=grep("hg38",builtin_annotations(),value=TRUE)
annot=c("hg38_basicgenes","hg38_genes_3UTRs","hg38_genes_5UTRs","hg38_enhancers_fantom" ,"hg38_genes_introns","hg38_genes_exons","hg38_genes_exonintronboundaries","hg38_genes_intronexonboundaries","hg38_genes_promoters","hg38_lncrna_gencode"  )
annot=sort(annot)

for(ii in 1:length(pd.ranges)){
  pd.data=eval(parse(text=pd.ranges[ii]))
  
  for(jj in 1:length(annot)){
    #####parietal.mid.meta=data.frame(parietal.mid.dmr.meta[overlp.gene.parietal.mid@from,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.parietal.mid@to],stringsAsFactors=F)
    
    
    builtannot=build_annotations(genome ='hg38', annotations = annot[jj])
    builtannot@seqnames@values=gsub("chr","",builtannot@seqnames@values)
    builtannot@seqnames@values=gsub("_.{0,20}","", builtannot@seqnames@values)
    
    builtannot.ranges=IRanges(start=builtannot@ranges@start,width=builtannot@ranges@width)
    builtannot=GRanges(seqnames=Rle(builtannot@seqnames),ranges=builtannot.ranges,mcols=builtannot@elementMetadata@listData)
    
    xx=findOverlaps(builtannot,pd.data)
    df=data.frame(id= as.character(builtannot@elementMetadata@listData$mcols.id[xx@from]),
                  gene.id=as.numeric(builtannot@elementMetadata@listData$mcols.gene_id[xx@from]),
                  gene=builtannot@elementMetadata@listData$mcols.symbol[xx@from],
                  type=builtannot@elementMetadata@listData$mcols.type[xx@from],
                  stringsAsFactors = FALSE)
    
    pd.data
    xx=cbind(eval(parse(text=gsub("dmr.{1-4}ranges","meta",pd.ranges[jj])))[xx@to,],df)
    
  }
}





#Narrow down to one p val per gene
options(scipen=999)

pari.mid.netwas=tapply(as.numeric(parietal.mid.meta$minFDR),as.factor(as.character(parietal.mid.meta$gene)),min)
xx=match(parietal.mid.meta.neg$gene,names(pari.mid.netwas))
parietal.mid.meta.neg=parietal.mid.meta.neg[-which(!is.na(xx)),]
pari.mid.netwas.neg=tapply(as.numeric(parietal.mid.meta.neg$medFDR),as.factor(as.character(parietal.mid.meta.neg$gene)),max)
pari.mid.netwas=data.frame(c(names(pari.mid.netwas.neg),names(pari.mid.netwas)),c(pari.mid.netwas.neg,pari.mid.netwas))
pari.mid.netwas=pari.mid.netwas[,c(1,rep(2,8))]
write.table(pari.mid.netwas,"NetWasPariMidInput.txt",col.names=FALSE)

pari.late.netwas=tapply(as.numeric(parietal.late.meta$minFDR),as.factor(as.character(parietal.late.meta$gene)),min)
xx=match(parietal.late.meta.neg$gene,names(pari.late.netwas))
parietal.late.meta.neg=parietal.late.meta.neg[-which(!is.na(xx)),]
pari.late.netwas.neg=tapply(as.numeric(parietal.late.meta.neg$medFDR),as.factor(as.character(parietal.late.meta.neg$gene)),max)
pari.late.netwas=data.frame(c(names(pari.late.netwas.neg),names(pari.late.netwas)),c(pari.late.netwas.neg,pari.late.netwas))
pari.late.netwas=pari.late.netwas[,c(1,rep(2,8))]
write.table(pari.late.netwas,"NetWasPariLateInput.txt",col.names=FALSE)

dir.create(paste0(getwd(),"/CingMid"))
cing.mid.netwas=tapply(as.numeric(cingulate.mid.meta$minFDR),as.factor(as.character(cingulate.mid.meta$gene)),min)
xx=match(cingulate.mid.meta.neg$gene,names(cing.mid.netwas))
cingulate.mid.meta.neg=cingulate.mid.meta.neg[-which(!is.na(xx)),]
cing.mid.netwas.neg=tapply(as.numeric(cingulate.mid.meta.neg$medFDR),as.factor(as.character(cingulate.mid.meta.neg$gene)),max)
cing.mid.netwas=data.frame(c(names(cing.mid.netwas.neg),names(cing.mid.netwas)),c(cing.mid.netwas.neg,cing.mid.netwas))

write.table(cing.mid.netwas[,c(1,1,rep(2,8))],"CingMid/NetWasCingMidInput.txt",col.names=FALSE,row.names=FALSE)

for(ii in 1:10){
scramble=cbind(cing.mid.netwas[,c(1,1)],cing.mid.netwas[sample(1:nrow(cing.mid.netwas)),rep(2,8)])
write.table(scramble,file=paste0("CingMid/NetWasCingMidInputScrmbl",ii,".txt"),col.names=FALSE,row.names=FALSE)
}

cing.late.netwas=tapply(as.numeric(cingulate.late.meta$minFDR),as.factor(as.character(cingulate.late.meta$gene)),min)
xx=match(cingulate.late.meta.neg$gene,names(cing.late.netwas))
cingulate.late.meta.neg=cingulate.late.meta.neg[-which(!is.na(xx)),]
cing.late.netwas.neg=tapply(as.numeric(cingulate.late.meta.neg$medFDR),as.factor(as.character(cingulate.late.meta.neg$gene)),max)
cing.late.netwas=data.frame(c(names(cing.late.netwas.neg),names(cing.late.netwas)),c(cing.late.netwas.neg,cing.late.netwas))
cing.late.netwas=cing.late.netwas[,c(1,rep(2,8))]
write.table(cing.late.netwas,"NetWasCingLateInput.txt",col.names=FALSE)



### Submit these files to http://hb.flatironinstitute.org/netwas







