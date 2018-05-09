
length(which(fivemc.glmmtmb.cingulate$NeocorticalVSNoneFDR<.01))
length(which(fivemc.glmmtmb.cingulate$LimbicVSNoneFDR<.01))

length(which(fivehmc.gamlss.cingulate$NeocorticalVSNoneFDR<.01))
length(which(fivehmc.gamlss.cingulate$LimbicVSNoneFDR<.01))

length(which(fivemc.glmmtmb.parietal$NeocorticalVSNoneFDR<.01))
length(which(fivemc.glmmtmb.parietal$LimbicVSNoneFDR<.01))

length(which(fivehmc.gamlss.parietal$NeocorticalVSNoneFDR<.01))
length(which(fivehmc.gamlss.parietal$LimbicVSNoneFDR<.01))


length(intersect(which(fivemc.glmmtmb.cingulate$NeocorticalVSNoneFDR<.01),which(fivehmc.gamlss.cingulate$NeocorticalVSNoneFDR<.01)))
length(intersect(which(fivemc.glmmtmb.cingulate$LimbicVSNoneFDR<.01),which(fivehmc.gamlss.cingulate$LimbicVSNoneFDR<.01)))

length(intersect(which(fivemc.glmmtmb.parietal$NeocorticalVSNoneFDR<.01),which(fivehmc.gamlss.parietal$NeocorticalVSNoneFDR<.01)))
length(intersect(which(fivemc.glmmtmb.parietal$LimbicVSNoneFDR<.01),which(fivehmc.gamlss.parietal$LimbicVSNoneFDR<.01)))
