
# Obtain average Beta values by sample groups

# Load in the 5mC Beta values and etc
load("fivemc_glmmtmb_out.RData")
load("fivehmc_gamlss_out.RData")

# Turn into numeric
mc.avg.pari=apply(fivemc.glmmtmb.parietal.fitted,2,as.numeric)
# Take averages by PD sample groups: None, Mid, Late
mc.avg.pari=data.frame(apply(mc.avg.pari[,10:12],1,mean),apply(mc.avg.pari[,13:15],1,mean),apply(mc.avg.pari[,16:18],1,mean))
# Add the FDR values as a column on the end
mc.avg.pari=data.frame(mc.avg.pari,fivemc.glmmtmb.parietal$LimbicVSNoneFDR,fivemc.glmmtmb.parietal$NeocorticalVSNoneFDR,mc.avg.pari[,2]-mc.avg.pari[,1],mc.avg.pari[,3]-mc.avg.pari[,1])

# Do the same for 5mc cingulate data
mc.avg.cing=apply(fivemc.glmmtmb.cingulate.fitted,2,as.numeric)
mc.avg.cing=data.frame(apply(mc.avg.cing[,1:3],1,mean),apply(mc.avg.cing[,4:6],1,mean),apply(mc.avg.cing[,7:9],1,mean))
mc.avg.cing=data.frame(mc.avg.cing,fivemc.glmmtmb.cingulate$LimbicVSNoneFDR,fivemc.glmmtmb.cingulate$NeocorticalVSNoneFDR,mc.avg.cing[,2]-mc.avg.cing[,1],mc.avg.cing[,3]-mc.avg.cing[,1])

# Do the same for 5hmc parietal data
hmc.avg.pari=apply(fivehmc.gamlss.parietal.fitted,2,as.numeric)
hmc.avg.pari=data.frame(apply(hmc.avg.pari[,10:12],1,mean),apply(hmc.avg.pari[,13:15],1,mean),apply(hmc.avg.pari[,16:18],1,mean))
hmc.avg.pari=data.frame(hmc.avg.pari,fivehmc.gamlss.parietal$LimbicVSNoneFDR,fivehmc.gamlss.parietal$NeocorticalVSNoneFDR,hmc.avg.pari[,2]-hmc.avg.pari[,1],hmc.avg.pari[,3]-hmc.avg.pari[,1])

# Do the same for 5hmc cingulate data
hmc.avg.cing=apply(fivehmc.gamlss.cingulate.fitted,2,as.numeric)
hmc.avg.cing=data.frame(apply(hmc.avg.cing[,1:3],1,mean),apply(hmc.avg.cing[,4:6],1,mean),apply(hmc.avg.cing[,7:9],1,mean))
hmc.avg.cing=data.frame(hmc.avg.cing,fivehmc.gamlss.cingulate$LimbicVSNoneFDR,fivehmc.gamlss.cingulate$NeocorticalVSNoneFDR,hmc.avg.cing[,2]-hmc.avg.cing[,1],hmc.avg.cing[,3]-hmc.avg.cing[,1])

# Carry over the Cpg IDs as the rownames
rownames(fivemc.glmmtmb.cingulate)=rownames(hmc.avg.pari)=rownames(mc.avg.pari)=rownames(hmc.avg.cing)=rownames(mc.avg.cing)

colnames(hmc.avg.pari)=colnames(mc.avg.pari)=colnames(hmc.avg.cing)=colnames(mc.avg.cing)=c("NoneAvg")

# Save everything to an .Rdata file
save(c(mc.avg.pari,mc.avg.cing,hmc.avg.pari,hmc.avg.cing),file="MethylAvgBetaResults.Rdata")   

