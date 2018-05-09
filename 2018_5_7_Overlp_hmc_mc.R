
setwd("/Users/cansav091/Desktop/CurrentProjects/P1MethylData")
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(DMRcate)
library(limma)
library(GenomicRanges)
library(beepr)
library(annotatr)

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
      for(jj in c("dmr","dmr.neg")){
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
        meta.dat[xx@to,1:4]=data.frame(mc.range@elementMetadata@listData[c(1:3,6)],stringsAsFactors=F)[xx@from,]
        # Pull stats from hmc data
        xx=findOverlaps(hmc.range,overlp)
        meta.dat[xx@to,1:4]=data.frame(hmc.range@elementMetadata@listData[c(1:3,6)],stringsAsFactors=F)[xx@from,]
  
      # Save objects to the global environment   
      assign(paste0(brain.stage,".",jj,".ranges"),overlp) 
      assign(paste0(brain.stage,".",jj),meta.dat) 
      result.obj=c(result.obj,paste0(brain.stage,".",jj,".ranges"),paste0(brain.stage,".",jj))
      pd.ranges=c(pd.ranges,paste0(brain.stage,".",jj,".ranges"),paste0(brain.stage,".",jj))
      
      }
    
    save(list=result.obj,file=paste0(brain.stage,".RData"))
}






