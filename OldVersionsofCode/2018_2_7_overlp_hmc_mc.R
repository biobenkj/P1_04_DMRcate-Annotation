
setwd("/Users/cansav091/Desktop/CurrentProjects/P1MethylData")

library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(DMRcate)
library(limma)
library(GenomicRanges)
library(beepr)

files=grep(".dmr.Rdata",dir(),value=TRUE)
mc.files=grep("^methylcyto",files,value=TRUE)
hmc.files=grep("^hydroxymethylcyto",files,value=TRUE)

file.pairs=match(gsub("methylcyto","",mc.files),gsub("hydroxymethylcyto","", hmc.files))


pd.ranges=c()
# Run through each dataset and find overlaps between hmc and mc
for(ii in 1:length(file.pairs)){
  # Load Data
  brain.stage=substr(gsub(".dmr.Rdata","",mc.files[ii]),12,nchar(mc.files[ii])-10)
  load(mc.files[ii])
  load(hmc.files[file.pairs[ii]])
  # Make lists of mc and hmc objects
  in.obj=grep("dmr",tolower(grep(brain.stage,ls(),value=TRUE)),value=TRUE)
  mc.obj=grep("^methylcyto",in.obj,value=TRUE)
  hmc.obj=grep("hydroxymethylcyto",in.obj,value=TRUE)
  
  result.obj=c()
      for(jj in c("dmr.ranges","dmr.neg.ranges")){
      mc.range=eval(parse(text=grep(jj,mc.obj,value=TRUE)))
      hmc.range=eval(parse(text=grep(jj,hmc.obj,value=TRUE)))
      # Used the saved objects to find the unifiedranges of 5hmC and 5mC for both stages and both brain areas
      # Then save in order to create 4 files total brain area*stage.  
      
      ##############################
      # File for parietal
      # Mid stage parietal, find overlap between mc and hmc
      overlp=union(mc.range,hmc.range)
        # Create a meta data frame to store this info in 
        n=length(overlp)
        meta.dat=data.frame(rep("NA",n),as.numeric(rep("NA",n)),as.numeric(rep("NA",n)),as.numeric(rep("NA",n)),stringsAsFactors=F)
        colnames(meta.dat)=c("chr.coord","noCpGs","minFDR","meanBeta")
        # Pull stats from mc data
        xx=findOverlaps(mc.range,overlp)
        meta.dat[xx@from,1:4]=data.frame(mc.range@elementMetadata@listData[c(1:3,6)],stringsAsFactors=F)[xx@to,]
        # Pull stats from hmc data
        xx=findOverlaps(hmc.range,overlp)
        meta.dat[xx@from,1:4]=data.frame(hmc.range@elementMetadata@listData[c(1:3,6)],stringsAsFactors=F)[xx@to,]
  
      # Save objects to the global environment   
      assign(paste0(brain.stage,".",jj,".union.ranges"),overlp) 
      assign(paste0(brain.stage,".",jj,".union.meta"),meta.dat) 
      result.obj=c(result.obj,paste0(brain.stage,".",jj,".union.ranges"),paste0(brain.stage,".",jj,".union.meta"))
      pd.ranges=c(pd.ranges,paste0(brain.stage,".",jj,".union.ranges"),paste0(brain.stage,".",jj,".union.meta"))
      
      }
    
    save(list=result.obj,file=paste0(brain.stage,".union.RData"))
}


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






