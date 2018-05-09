## Venn Diagrams for limma and glmm results.

library(gplots)
library(VennDiagram)
library(plotrix)
home="/Users/cansav091/Desktop/VennDiagrams"
setwd(home)
allfiles=dir(home)
stage=c("mid","late")
stage.labs=c("Mid","Late")
brain=c(".c.",".p.")
brain.labs=c("Cingulate","Parietal")
colrs=sapply(c("#CB5A28","#6E005F"),color.id)
  
for(ii in 1:length(stage)){
s.files=grep(stage[ii],allfiles,value=TRUE)
  for(jj in 1:length(brain)){
files=grep(brain[jj],s.files,value=TRUE)

glm=as.vector(t(read.table(grep("glmm",files,value=TRUE))))
lim=as.vector(t(read.table(grep("lim",files,value=TRUE))))

grid.newpage()
ven=draw.pairwise.venn(area1 = length(glm),
                   area2 = length(lim),
                   cat.pos=c(0,90),
                   cat.dist=c(-.1,.1),
                   cat.fontfamily="Arial",
                   cex=3,
                   cat.cex=3,
                   ext.percent=c(.5,.5,.5),
                   ext.pos=c(20,45),
                   ext.length=c(.7,.9),
                   ext.dist=c(-.08,.1),
                   fontfamily="Arial",
                   fill= c(adjustcolor(colrs[ii],alpha=.7),adjustcolor(colrs[ii],alpha=.7)),
                   cross.area=length(intersect(glm,lim)),
                   category = c("glmm-TMB","limma"),
                   main=paste0(stage.labs[ii],brain[jj]),
                   margin=.08)
tiff(filename=paste0("glmmvs.limmaVennDiagram",stage.labs[ii],brain[jj]))
grid.draw(ven)
dev.off()
  }
}

