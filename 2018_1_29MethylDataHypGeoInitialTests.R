

pcut=.01
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggplot2)
library(GenomicRanges)
library(IRanges)
library(calibrate)
library(beepr)

binom=function(n,k){
  x=factorial(n)/(factorial(k)*factorial(n-k))
  return(x)
}

data=grep(".RData",dir(),value=TRUE)
load(data[2])


#### Download tissue specific functional data. 
setwd("/Users/cansav091/Desktop/Current Projects /P1MethylData/P1MethylData/Genome Annotation")
seg.annot=read.table("segway_encyclopedia.bed", sep="\t",skip=1)
colnames(seg.annot)=t(read.table("segway_encyclopedia.bed", sep="\t",nrows=1))
cing.annot=seg.annot$BRAIN_CINGULATE_GYRUS
sn.annot=seg.annot$BRAIN_SUBSTANTIA_NIGRA
seg.annot=GRanges(seqnames=Rle(gsub("chr","",seg.annot$chrom)),IRanges(start=seg.annot$start,end=seg.annot$end),mcols=data.frame(seg.annot$BRAIN_CINGULATE_GYRUS,seg.annot$BRAIN_SUBSTANTIA_NIGRA))

length(which(as.character(seg.annot@elementMetadata@listData$mcols.seg.annot.BRAIN_CINGULATE_GYRUS)==as.character(seg.annot@elementMetadata@listData$mcols.seg.annot.BRAIN_SUBSTANTIA_NIGRA)))
length(seg.annot@elementMetadata@listData$mcols.seg.annot.BRAIN_CINGULATE_GYRUS)


#Upload Methylation data
setwd("/Users/cansav091/Desktop/Current Projects /P1MethylData/P1MethylData")
lists=grep("five",ls()[-1])


# Upload annotation 
load("EPIC.manifest.rda")

#EPIC.manifest@seqnames@values=EPIC.manifest@elementMetadata@listData$chrmA
#EPIC.manifest@seqnames@values=as.factor(gsub("chr","",EPIC.manifest@seqnames@values))

gene.names = mappedkeys(org.Hs.egSYMBOL)
gene.names=as.data.frame(org.Hs.egSYMBOL[gene.names])

chrbegin = mappedkeys(org.Hs.egCHRLOC)
chrbegin=as.data.frame(org.Hs.egCHRLOC[chrbegin])

xx=match(chrbegin$gene_id,gene.names$gene_id)
gene.names=as.character(gene.names$symbol[xx])

chrs=as.factor(chrbegin$Chromosome)
chrbegin=abs(as.numeric(chrbegin$start_location))

### Download the chr coordinates for the ends of the genes
chrend = mappedkeys(org.Hs.egCHRLOCEND)
chrend=as.data.frame(org.Hs.egCHRLOCEND[chrend])
chrend=abs(as.numeric(chrend$end_location))

xx=grep("_",chrs,invert=TRUE)
chrs=as.factor(as.character(chrs[xx]))
chrbegin=chrbegin[xx]
chrend=chrend[xx]
gene.names.v=gene.names[xx]
gene.names=data.frame(gene.names[xx])
gene.names=data.frame(gene.names,rep("NA",nrow(gene.names)),rep("NA",nrow(gene.names)))
                                         
starts.ends=IRanges(start=chrbegin,end=chrend)
gene.locs=GRanges(seqnames=Rle(chrs),ranges=starts.ends,mcols=gene.names)

overlp.gene.func=findOverlaps(seg.annot,gene.locs,maxgap = 1000)

annot=data.frame(gene.names[overlp.gene.func@queryHits,1],seg.annot[overlp.gene.func@queryHits]@elementMetadata@listData$mcols.seg.annot.BRAIN_CINGULATE_GYRUS,seg.annot[overlp.gene.func@queryHits]@elementMetadata@listData$mcols.seg.annot.BRAIN_SUBSTANTIA_NIGRA)
colnames(annot)=c("Gene","CingFunc","SNFunc")

EPIC.manifest@ranges@width=as.integer(EPIC.manifest@ranges@width+4)


probe.ranges=GRanges(seqnames=Rle(gsub("chr","",EPIC.manifest@elementMetadata@listData$chrmA)),ranges =EPIC.manifest@ranges)

overlp.gene=findOverlaps(probe.ranges,gene.locs,maxgap = 1000) #865,918 elements
overlp.func=findOverlaps(probe.ranges,seg.annot,maxgap = 1000) #865,918 elements

length(unique(overlp.gene@queryHits)) #622418 overlap with a gene within 1 kb
length(unique(overlp.func@queryHits)) #496356 overlap with a function within 1 kb


CpG.gene.hits=probe.ranges[overlp.gene@queryHits]
CpG.gene.hits@elementMetadata@listData=annot[overlp.gene@subjectHits,]


CpGGenekey=data.frame(CpG.gene.hits@ranges,CpG.gene.hits@elementMetadata@listData)
write.csv(CpGGenekey,file="CpGtoGeneKeyOverlap.csv")

write.csv(table(CpGGenekey$CingFunc,CpGGenekey$SNFunc),file="CingtoSNFunctionalAnnot.csv")


length(which(CpGGenekey$CingFunc!=CpGGenekey$SNFunc))/length(CpGGenekey$SNFunc)
length(which(CpGGenekey$CingFunc==CpGGenekey$SNFunc))

regions=c(table(CpGGenekey$CingFunc),table(CpGGenekey$SNFunc))
regions[seq(from=1,to=length(regions),by=2)]=table(CpGGenekey$CingFunc)
regions[seq(from=2,to=length(regions),by=2)]=table(CpGGenekey$SNFunc)
names.regs=c(1:length(regions))
names.regs[seq(from=1,to=length(regions),by=2)]=paste(names(table(CpGGenekey$CingFunc)),"Cing Ctx")
names.regs[seq(from=2,to=length(regions),by=2)]=paste(names(table(CpGGenekey$CingFunc)),"SN")

jpeg("CingandSNannotationComparison.jpeg")
par(mar=c(8,5,1,5))
barplot(regions,las=2,names.arg=names.regs,col=rep(c("darkgreen","blue"),6))
dev.off()


datasets=grep("fivemc",ls(),value=TRUE)
datasets=grep(".fitted",datasets,invert=TRUE,value=TRUE)
output=c()




for(ii in 1:length(datasets)){
  
  dataset=eval(parse(text=datasets[ii]))
  xx=match(CpGGenekey$names,rownames(dataset))
  
  key=CpGGenekey[which(!is.na(xx)),]
  dataset=dataset[xx[!is.na(xx)],]

  which(key$SNFunc==
  
  probprgene=table(as.factor(as.character(key$Gene))) # 12,548 genes
  probprcingfunc=table(as.factor(as.character(key$CingFunc)))
  probprsnfunc=table(as.factor(as.character(key$SNFunc)))
  
  nn=length(unique(names(probprcingfunc)))
  hyp.geo.pval.cing.neo=1:nn
  hyp.geo.pval.sn.limb=1:nn
  hyp.geo.pval.cing.limb=1:nn
  hyp.geo.pval.sn.neo=1:nn
  Kneo=length(which(dataset$NeocorticalVSNoneFDR<pcut))
  Klimb=length(which(dataset$LimbicVSNoneFDR<pcut))
  ######### By functional annotation 
  for(jj in 1:nn){
  
    for(hh in 1:2){
    if(hh==1){
    xx=grep(names(probprcingfunc[jj]),key$CingFunc)
    n=probprcingfunc[jj]
    }else{
    xx=grep(names(probprsnfunc[jj]),key$SNFunc)
    n=probprsnfunc[jj]
    }
    
  Kneo=length(which(dataset$NeocorticalVSNoneFDR<pcut))
  Klimb=length(which(dataset$LimbicVSNoneFDR<pcut))
  N=nrow(dataset)    
  kneo=length(which(dataset$NeocorticalVSNoneFDR[xx]<pcut))
  klimb=length(which(dataset$LimbicVSNoneFDR[xx]<pcut))
  
  if(hh==1){
    if(kneo>0){
      hyp.geo.pval.cing.neo[jj]=dhyper(kneo,Kneo,N-Kneo,n)
      }
    if(klimb>0){
      hyp.geo.pval.cing.limb[jj]=dhyper(klimb,Klimb,N-Klimb,n)
      }
  }else{
    if(kneo>0){
      hyp.geo.pval.sn.neo[jj]=dhyper(kneo,Kneo,N-Kneo,n)
    }
    if(klimb>0){
      hyp.geo.pval.sn.limb[jj]=dhyper(klimb,Klimb,N-Klimb,n)
    }
      }
    }
  }
  reg.output=data.frame(hyp.geo.pval.cing.neo,hyp.geo.pval.cing.limb,hyp.geo.pval.sn.neo,hyp.geo.pval.sn.limb)
  reg.output=data.frame(reg.output,p.adjust(hyp.geo.pval.cing.neo,method="fdr"),p.adjust(hyp.geo.pval.cing.limb,method="fdr"),p.adjust(hyp.geo.pval.sn.neo,method="fdr"),p.adjust(hyp.geo.pval.sn.limb,method="fdr"))
  
  xx=c("EnrichCingNeo","EnrichCingLimb","EnrichSNNeo","EnrichSNLimb")
  colnames(reg.output)=c(xx,paste0("fdrqval",xx))
  rownames(reg.output)=unique(names(probprcingfunc))
  write.csv(reg.output,file=paste0("FuncRegionEnrichments",datasets[ii],".csv"),quote=F)
  ##### Plotting enrichment across brain regions and PD stages 
  jpeg(paste0("CingulateEnrichments",datasets[ii],".jpeg"))
  par(mar=c(8,5,1,5))
  barplot(-log2(c(reg.output$fdrqvalEnrichCingNeo,reg.output$fdrqvalEnrichCingLimb)), names.arg=c(paste("Neo",unique(names(probprcingfunc))),paste("Limbic",unique(names(probprcingfunc)))),las=2,cex.names = .7)
  abline(h=3,col="red")
  dev.off()
  
  jpeg(paste0("SNEnrichments",datasets[ii],".jpeg"))
  par(mar=c(8,5,1,5))
  barplot(-log2(c(reg.output$fdrqvalEnrichSNNeo,reg.output$fdrqvalEnrichSNLimb)), names.arg=c(paste("Neo",unique(names(probprcingfunc))),paste("Limbic",unique(names(probprcingfunc)))),las=2,cex.names = .7)
  abline(h=3,col="red")
  dev.off()
  
  #yy=gsub(".2$","",rownames(dataset))
  #yy=gsub(".1$","",yy)

  hyp.geo.pval.neo=hyp.geo.pval.limb=ratio.neo=ratio.limb=rep(NA,length(names(probprgene)))
  
  for(jj in 1:length(unique(names(probprgene)))){
  xx=grep(names(probprgene[jj]),key$Gene)
  
  kneo=length(which(dataset$NeocorticalVSNoneFDR[xx]<pcut))
  klimb=length(which(dataset$LimbicVSNoneFDR[xx]<pcut))
                                      
  n=probprgene[jj]
  ratio.neo[jj]=kneo/n
  ratio.limb[jj]=klimb/n
  
    if(kneo>0){
    hyp.geo.pval.neo[jj]=dhyper(kneo,Kneo,N-Kneo,n)
    }
    if(klimb>0){
    hyp.geo.pval.limb[jj]=dhyper(klimb,Klimb,N-Klimb,n)
      }
  }
  beep(sound=2)
  output=data.frame(probprgene,hyp.geo.pval.neo,p.adjust(hyp.geo.pval.neo,method="fdr"),ratio.neo,hyp.geo.pval.limb,p.adjust(hyp.geo.pval.limb,method="fdr"),ratio.limb)
  colnames(output)=c("GeneSymbol","NumberofAssocCpGs","Uncorrectedpvalneocortical","FDRqvalneocortical","RatioofSignif-Nonsignifprobesneo","Uncorrectedpvallimbic","FDRqvallimbic","RatioofSignif-Nonsignifprobeslimbic")
  write.csv(output,quote=FALSE,file=paste0("HypergeoTestGene",datasets[ii],".csv"))
  

  
  jpeg(paste0("NeocorticalvsLimbicCorrelation",datasets[ii],".jpeg"))
  plot(-log10(output$Uncorrectedpvalneocortical),-log10(output$Uncorrectedpvallimbic),pch=20)
  xx=which(is.na(output$Uncorrectedpvallimbic))
  xx=unique(c(xx,which(is.na(output$Uncorrectedpvalneocortical))))
  #textxy(-log10(output$Uncorrectedpvalneocortical[-xx]),-log10(output$Uncorrectedpvallimbic[-xx]),as.character(output$GeneSymbol)[-xx])
  dev.off()
  
}
beep(sound=4)


res=grep("Hypergeo",dir(),value=T)
cingres=read.csv(res[1])
parres=read.csv(res[2])

xx=which(cingres$FDRqvalneocortical<.01 & parres$FDRqvalneocortical<.01)
yy=which(cingres$FDRqvallimbic<.01 & parres$FDRqvallimbic<.01)

res=intersect(xx,yy)
final.res=data.frame(cingres$GeneSymbol[res],cingres$FDRqvalneocortical[res],parres$FDRqvalneocortical[res],cingres$FDRqvallimbic[res],parres$FDRqvallimbic[res])

write.csv(final.res,file="Finalresults.hmc.csv")

#Network Output
data=read.csv("NetWAS_Score.csv",skip=2)
network=as.character(data$X.NetWAS_Score..Distance.from.hyperplane)
network=unlist(strsplit(network,"\t"))
network=data.frame(network[seq(from=1,to=length(network),by=3)],network[seq(from=2,to=length(network),by=3)],network[seq(from=3,to=length(network),by=3)])
colnames(network)=c("Symbol","Direction","Strength")
network[,3]=as.numeric(as.character(network[,3]))

data=read.table("NetWasResultsFinalImage.csv",skip=1)
data=unlist(strsplit(as.character(data$V1),"\t"))
df.res=data.frame(data[seq(from=1, to=length(data),by=3)],data[seq(from=2, to=length(data),by=3)],data[seq(from=3, to=length(data),by=3)])
df.res=data.frame


col=read.table("NetWasResultsFinalImage.csv",nrow=1,sep="\t",colClasses="character")
col=unlist(strsplit(as.character(col),"\t"))
col=gsub(",","",col)

colnames(df.res)=col

df.res$NetWAS_Score=as.numeric(as.character(df.res$NetWAS_Score))

hist(df.res$NetWAS_Score)
barplot(table(df.res$Training_Set))

xx=match(df.res$Symbol[which(df.res$Training_Set==0)],output$GeneSymbol)
hist(output$FDRqvalneocortical[xx])


xx=df.res[which(df.res$Training_Set==1),] # 454 have 1
xx=xx[(order(xx$NetWAS_Score,decreasing=TRUE)),]

barplot(xx$NetWAS_Score[1:25],names.arg=xx$Symbol[1:25],las=2)









