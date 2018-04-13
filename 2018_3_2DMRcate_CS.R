
# Make sure you have these packages installed. 
# This main chunk of this code is modified from the cpg.annotate() function from the DMRcate package. 
setwd("/Users/cansav091/Desktop/Current Projects/P1MethylData/P1MethylData")
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(DMRcate)
library(limma)
library(GenomicRanges)
library(beepr)
library(plyr)

# Load in the beta value datasets
data=c("fivemc_glmmtmb_out.RData","fivehmc_gamlss_out.RData")
load(data[1])
load(data[2])


# Obtain a list of the beta value datasets
beta.data=grep(".fitted",ls(),value=TRUE)
rm(list=gsub(".fitted","",beta.data))

brain.reg=unlist(strsplit(beta.data,"\\."))
cytmod=brain.reg[seq(from=1,to=length(brain.reg),by=4)]
brain.reg=brain.reg[seq(from=3,to=length(brain.reg),by=4)]


# Do this for each data set. 
for(ii in 1:length(beta.data)){
samples=ifelse((length(grep("cingulate",beta.data[ii]))>0),1,10)
samples=c(samples:(samples+8))

beta.vals=apply(eval(parse(text=beta.data[ii])),2,as.numeric)
rownames(beta.vals)=rownames(eval(parse(text=beta.data[ii]))) # Set the CpG Ids as the rownames

# Make a design matrix for each late and mid stage comparisons to None
stage <- factor(rep(c("None","Mid"), each = 3),c("None","Mid"))
des.mat.mid=model.matrix(~stage)

stage <- factor(rep(c("None","Late"), each = 3),c("None","Late"))
des.mat.late=model.matrix(~stage)

# Settings for the analysis
for(jj in 1:2){
  if(jj==1){
    object= beta.vals[,samples][,c(1:6)]
    design=des.mat.mid
    stage="mid"
  }else{
    object= beta.vals[,samples][,c(1:3,7:9)]
    design=des.mat.late
    stage="late"
  }
na.rows=-which(is.na(object[,1]))
object=object[-which(is.na(object[,1])),]
datatype = "array"
annotation = "IlluminaHumanMethylationEPICanno.ilm10b2.hg19"
analysis.type = "differential"
contrasts = FALSE
cont.matrix = NULL
fdr = 0.05
coef=2

      
      # Download the EPIC annotation
      RSobject = RatioSet(object, annotation = annotation)
      RSanno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
      
      # Match the annotation to our dataset
      xx=match(rownames(object),RSanno@rownames)
      
      object=object[which(!is.na(xx)),]
      
      # Fit the design matrix
      fit = lmFit(object, design)
      fit = eBayes(fit)
      
      # Assess the number of significant CpG's
      tt = topTable(fit, coef = coef, number = nrow(object))
      nsig = sum(tt$adj.P.Val < fdr)
      
      # Fit the data to a Bayes model
      betafit = lmFit(ilogit2(object), design)
      betafit = eBayes(betafit)
      betatt = topTable(betafit, coef = coef, number = nrow(object))
      
      # match the CpGs to the fit and only keep those
      m = match(rownames(tt), rownames(betatt))
      tt$betafc = betatt$logFC[m]
      m = match(rownames(tt), rownames(object))
      object = object[m, ]
      # Combine stats and annotation into a dataframe
      stat = tt$t

      annotated = data.frame(ID = rownames(object), stat = stat[which(!is.na(xx))], CHR = RSanno$chr[xx[which(!is.na(xx))]],
                             pos = RSanno$pos[xx[which(!is.na(xx))]], betafc = tt$betafc,indfdr = tt$adj.P.Val)
      
      # Put in order of chr coord
      annotated = annotated[order(annotated$CHR, annotated$pos), ]
      
      # Coerce to object class "annotated"
      class(annotated) = "annot"
    
 
# Run the annotated object into the DMR analysis. Adjust the parameters as you see fit.    
DMR=dmrcate(annotated,
        lambda = 1000,
        C=NULL,
        p.adjust.method = "BH",
        pcutoff = "fdr",
        consec = FALSE,
        conseclambda = 10,
        betacutoff = NULL,
        min.cpgs = 2,
        mc.cores = 1)


# Pull the chromosomal coordinates out of the DMRcate results.
    coord=DMR$results$coord
    ends=as.numeric(unlist(strsplit(coord,"-"))[seq(from=2,to=length(coord)*2,by=2)])
    starts=unlist(strsplit(coord,":"))
    chrs=starts[seq(from=1,to=length(coord)*2,by=2)]
    chrs=gsub("chr","",chrs)
    starts=starts[seq(from=2,to=length(coord)*2,by=2)]
    starts=as.numeric(unlist(strsplit(starts,"-"))[seq(from=2,to=length(starts)*2,by=2)])

# Use the chr coordinates to make a GRanges and add the results matrix as mcols in the listData of the GRange
    DMR.ranges=IRanges(start=starts,end=ends)
    DMR.ranges=GRanges(seqnames=Rle(chrs),ranges=DMR.ranges,mcols=DMR$results)

    input=DMR$input
# Save GRange object and DMR results to an .RData file        
    save(list=c("DMR","DMR.ranges","input"),file=paste0(cytmod[ii],".",brain.reg[ii],".",stage,".dmr.Rdata"))

# Save objects to the global environment   
    assign(paste0(cytmod[ii],".",brain.reg[ii],".",stage,".range"),DMR.ranges,envir=.GlobalEnv) 
    assign(paste0(cytmod[ii],".",brain.reg[ii],".",stage,".dmr"),DMR,envir=.GlobalEnv) 
    }
}
beep(sound=2)



# Used the saved objects to find the unifiedranges of 5hmC and 5mC for both stages and both brain areas
# Then save in order to create 4 files total brain area*stage.  
      
      #######
      #File for parietal
      mod.pari.mid.range=union( fivemc.parietal.mid.range, fivehmc.parietal.mid.range)

      xx=findOverlaps(fivemc.parietal.mid.range,mod.pari.mid.range)
      df=data.frame(as.numeric(rep("NA",length(mod.pari.mid.range))),as.numeric(rep("NA",length(mod.pari.mid.range))),as.character(rep("NA",length(mod.pari.mid.range))))
      colnames(df)=c("minFDR","meanBeta","gene")
      
      mod.pari.mid.range@elementMetadata@listData=df
      mod.pari.mid.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivemc.parietal.mid.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.pari.mid.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivemc.parietal.mid.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
        
      xx=findOverlaps(fivehmc.parietal.mid.range,mod.pari.mid.range)
      mod.pari.mid.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivehmc.parietal.mid.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.pari.mid.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivehmc.parietal.mid.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      ################################
      
      mod.pari.late.range=union( fivemc.parietal.late.range, fivehmc.parietal.late.range)
      
      xx=findOverlaps(fivemc.parietal.late.range,mod.pari.late.range)
      df=data.frame(as.numeric(rep("NA",length(mod.pari.late.range))),as.numeric(rep("NA",length(mod.pari.late.range))),as.character(rep("NA",length(mod.pari.late.range))))
      colnames(df)=c("minFDR","meanBeta","gene")
      
      mod.pari.late.range@elementMetadata@listData=df
      mod.pari.late.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivemc.parietal.late.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.pari.late.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivemc.parietal.late.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      
      xx=findOverlaps(fivehmc.parietal.late.range,mod.pari.late.range)
      mod.pari.late.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivehmc.parietal.late.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.pari.late.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivehmc.parietal.late.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
    
      save(list=c("mod.pari.mid.range","mod.pari.late.range"),file=paste0("mod.parietal.dmr.range.RData"))
      ################################
      

      
      ################################
      #File for cingulate
      mod.cing.mid.range=union( fivemc.cingulate.mid.range, fivehmc.cingulate.mid.range)
      
      xx=findOverlaps(fivemc.cingulate.mid.range,mod.cing.mid.range)
      df=data.frame(as.numeric(rep("NA",length(mod.cing.mid.range))),as.numeric(rep("NA",length(mod.cing.mid.range))),as.character(rep("NA",length(mod.cing.mid.range))))
      colnames(df)=c("minFDR","meanBeta","gene")
      
      mod.cing.mid.range@elementMetadata@listData=df
      mod.cing.mid.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivemc.cingulate.mid.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.cing.mid.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivemc.cingulate.mid.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      
      xx=findOverlaps(fivehmc.cingulate.mid.range,mod.cing.mid.range)
      mod.cing.mid.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivehmc.cingulate.mid.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.cing.mid.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivehmc.cingulate.mid.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      ################################
      
      mod.cing.late.range=union( fivemc.cingulate.late.range, fivehmc.cingulate.late.range)
      
      xx=findOverlaps(fivemc.cingulate.late.range,mod.cing.late.range)
      df=data.frame(as.numeric(rep("NA",length(mod.cing.late.range))),as.numeric(rep("NA",length(mod.cing.late.range))),as.numeric(rep("NA",length(mod.cing.late.range))),as.character(rep("NA",length(mod.cing.late.range))))
      colnames(df)=c("chr.coord","minFDR","meanBeta","gene")
      
      mod.cing.late.range@elementMetadata@listData=df
      mod.cing.late.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivemc.cingulate.late.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.cing.late.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivemc.cingulate.late.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      
      xx=findOverlaps(fivehmc.cingulate.late.range,mod.cing.late.range)
      mod.cing.late.range@elementMetadata@listData$minFDR[xx@subjectHits]=as.numeric(fivehmc.cingulate.late.range@elementMetadata@listData$mcols.minfdr[xx@queryHits])
      mod.cing.late.range@elementMetadata@listData$meanBeta[xx@subjectHits]=fivehmc.cingulate.late.range@elementMetadata@listData$mcols.meanbetafc[xx@queryHits]
      
      save(list=c("mod.cing.mid.range","mod.cing.late.range"),file=paste0("mod.cingulate.dmr.range.RData"))
      ################################
      

    