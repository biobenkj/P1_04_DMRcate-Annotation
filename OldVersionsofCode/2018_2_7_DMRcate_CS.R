

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
  summary(annotated$indfdr)
  

sort(annotated$indfdr)[round(length(annotated$indfdr)*.0001)]  
      
      
# Run the annotated object into the dmr analysis. Adjust the parameters as you see fit.    
dmr=dmrcate(annotated,
        lambda = 1000,
        C=100,
        p.adjust.method = "BH",
        pcutoff = 1e-0300,
        consec = FALSE,
        conseclambda = 10,
        betacutoff = NULL,
        min.cpgs = 6,
        mc.cores = 1,
        neg=FALSE)

dmr.neg=dmrcate(annotated,
                lambda = 1000,
                C=100,
                p.adjust.method = "BH",
                pcutoff = 1e-0300,
                consec = FALSE,
                conseclambda = 10,
                betacutoff = NULL,
                min.cpgs = 6,
                mc.cores = 1,
                neg=TRUE)


# Pull the chromosomal coordinates out of the dmrcate results.
    pos.coord=dmr$results$coord
    neg.coord=dmr.neg$results$coord
    for(jj in 1:2){
      if(jj==1){
        coord=pos.coord
      }else{
        coord=neg.coord
      }
    ends=as.numeric(unlist(strsplit(coord,"-"))[seq(from=2,to=length(coord)*2,by=2)])
    starts=unlist(strsplit(coord,":"))
    chrs=starts[seq(from=1,to=length(coord)*2,by=2)]
    chrs=gsub("chr","",chrs)
    starts=starts[seq(from=2,to=length(coord)*2,by=2)]
    starts=as.numeric(unlist(strsplit(starts,"-"))[seq(from=1,to=length(starts)*2,by=2)])
    # Use the chr coordinates to make a GRanges and add the results matrix as mcols in the listData of the GRange
      if(jj==1){
      dmr.ranges=IRanges(start=starts,end=ends)
      dmr.ranges=GRanges(seqnames=Rle(chrs),ranges=dmr.ranges,mcols=dmr$results)
      input=dmr$input
        }else{
      dmr.neg.ranges=IRanges(start=starts,end=ends)
      dmr.neg.ranges=GRanges(seqnames=Rle(chrs),ranges=dmr.neg.ranges,mcols=dmr.neg$results)
      }
    
    }
    # Create names specific for each object
    obj=paste0(cytmod[ii],".",brain.reg[ii],".",stage,".")
    obj=paste0(obj,c("dmr.ranges","dmr","dmr.neg.ranges","dmr.neg"))
    # Save objects to the global environment   
    assign(obj[1],dmr.ranges) 
    assign(obj[2],dmr) 
    assign(obj[3],dmr.neg.ranges) 
    assign(obj[4],dmr.neg) 
  # Use the chr coordinates to make a GRanges and add the results matrix as mcols in the listData of the GRange
  # Save GRange object and DMR results to an .RData file        
    save(list=obj,file=paste0(cytmod[ii],".",brain.reg[ii],".",stage,".dmr.Rdata"))

  }
  beep(sound=2)
  }
beep(sound=8)







