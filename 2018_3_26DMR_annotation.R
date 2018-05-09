source("https://bioconductor.org/biocLite.R")
biocLite("annotatr")
library(org.Hs.eg.db)

source("https://bioconductor.org/biocLite.R")
biocLite("IRanges")
library(IRanges)
biocLite("GenomicRanges")
library(GenomicRanges)

########################################################
# Here is where you load in your info
chr.locs=read.csv('DMR_BED.csv')

# First you need to create a i ranges with your chromosome locations
i.ranges=IRanges(start=chr.locs$chromStart,end=chr.locs$chromEnd)

#now you can make a granges using the i.ranges object you just created and also putting in the corrresponding chromosomes.
# You can also add any other annotation info from your .csv file by putting it as an mcols in your granges
g.ranges=GRanges(seqnames=Rle(chr.locs$chrom),ranges=i.ranges,mcols=chr.locs$other)

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


overlp.gene.g.ranges=findOverlaps((g.ranges))

# Now overlap these gene locations to our ranges
g.ranges.meta=data.frame(g.ranges[overlp.gene.g.ranges@to,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.g.ranges@from],stringsAsFactors=F)

# Add gene names to its own column over to the modifcation meta data frames
g.ranges.meta=data.frame(g.ranges.meta[overlp.gene.g.ranges@to,],gene.locs@elementMetadata@listData$mcols.gene.names.xx.[overlp.gene.g.ranges@from],stringsAsFactors=F)

colnames(g.ranges.meta)=colnames(g.ranges.meta)=c("chr.coord","noCpGs","minFDR","meanBeta","gene")

save(list=c("g.ranges.meta"),file=paste0("genes.RData"))


write.table(g.ranges.meta,file="g.ranges.meta.txt",sep="\t")

