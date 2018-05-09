library(plyr)
library(dplyr)
library(gamlss)
library(parallel)

#Load in a chunk
all_bvals_hmc.t.chnk1 <- readRDS("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/all_bvals_hmc.t.chnk1.rda")

#New meta data that also includes the glial proportions
metadata.c <- read.table("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/meta_data_c.txt", header = T, check.names = F, stringsAsFactors = F, nrows = 18)
metadata.p <- read.table("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/meta_data_p.txt", header = T, check.names = F, stringsAsFactors = F, nrows = 18)

#Merge together
metadata <- rbind(metadata.c, metadata.p)

#Get the glial comps
glial <- metadata$glial[which(metadata$Assay == "B")]

#Build the design

brain_region <- as.factor(c(rep("Cingulate", 9), rep("Parietal", 9)))

C_None <- as.factor(c(rep(1, 3), rep(0, 15)))

C_Limbic <- as.factor(c(rep(0, 3), rep(1, 3), rep(0, 12)))

C_Neocortical <- as.factor(c(rep(0, 6), rep(1, 3), rep(0, 9)))

P_None <- as.factor(c(rep(0, 9), rep(1, 3), rep(0, 6)))

P_Limbic <- as.factor(c(rep(0, 12), rep(1, 3), rep(0, 3)))

P_Neocortical <- as.factor(c(rep(0, 15), rep(1, 3)))

#Set randeff
randeff <- all_bvals_hmc.t.chnk1$all_arrays

#Remove the all_arrays column so as not to run regression over it
all_bvals_hmc.t.chnk1$all_arrays <- NULL

#Filter columns with all zeros

keep <- colSums(all_bvals_hmc.t.chnk1)>0

all_bvals_hmc.t.chnk1 <- all_bvals_hmc.t.chnk1[,keep]

#recode values to 0 or 1 for filter
boolfilt <- t(all_bvals_hmc.t.chnk1)
boolfilt[boolfilt > 0] <- 1

#only keep probes that have at least 3 sample values greater than 0
keep <- rowSums(boolfilt) >= 3

all_bvals_hmc.t.chnk1 <- t(t(all_bvals_hmc.t.chnk1)[keep,])

#Create the coefs for plotting with se
createCoeftab <- function(GAM) {
    bGAM <- try(coef(GAM)[2:6])
    if (class(bGAM) == "try-error") {
      bGAM <- seq(1:5)
    }
    #seGAM <- try(diag(vcov(GAM))[2:6])
    #if (class(seGAM) == "try-error") {
    #  seGAM <- rep(0, 5)
    #}
    nms <- try(names(bGAM))
    if (class(nms) == "try-error") {
      nms <- paste0("Noname", 1:5)
    }
    df <- data.frame(model    = rep(c("gamlssBEZI"), each = 5),
                     term     = nms,
                     estimate = unname(c(bGAM)))
    #df <- transform(df,
    #                upper = estimate + sqrt(c(seGAM)),
    #                lower = estimate - sqrt(c(seGAM)))
    df
}

#ZIBeta regression with a logit link function
EPIC.beta.fit <- function(x, cellcomp, randeff) {
  x <- as.data.frame(x)
  colnames(x) <- "x"
  x$randeff <- randeff
  x$glial <- cellcomp

  #Fit models for all tissues
  #Fit with zibeta regression with C_None as the reference
  cnone.fit <- gamlss(x ~ C_Limbic + C_Neocortical + P_None + P_Limbic + P_Neocortical + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Refit with reduced model C_Limbic
  cnone.fit2 <- gamlss(x ~ C_Neocortical + P_None + P_Limbic + P_Neocortical + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Refit with reduced model C_Neocortical
  cnone.fit3 <- gamlss(x ~ C_Limbic + P_None + P_Limbic + P_Neocortical + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Calc LRT p-values for cingulate limbic tissues relative to C_none
  cnone.lrt.limbic <- try(LR.test(cnone.fit2, cnone.fit, print = F)$p.val)
  if (class(cnone.lrt.limbic) == "try-error") {
    cnone.lrt.limbic <- 1
  }
  #Calc LRT p-values for cingulate neocortical tissues relative to C_none
  cnone.lrt.neocortical <- try(LR.test(cnone.fit3, cnone.fit, print = F)$p.val)
  if (class(cnone.lrt.neocortical) == "try-error") {
    cnone.lrt.neocortical <- 1
  }
  #Extract fitted values
  cnone.fit.vals <- fitted(cnone.fit)
  #Fit with beta regression with C_Limbic as reference
  #climbic.fit <- glmmTMB(x ~ C_None + C_Neocortical + P_None + P_Limbic + P_Neocortical + (1 | randeff), data = x, family = list(family = "beta", link = "logit"), se = T)
  #Calc LRT p-values for cingulate tissues relative to C_none
  #climbic.lrt <- drop1(climbic.fit, scope = "C_Neocortical", test = "Chisq")
  #Extract fitted values
  #cnone.fit.vals2 <- fitted(climbic.fit)

  #Get the beta estimates for fold-change calcs
  cnone.fit.est <- createCoeftab(cnone.fit)[1:2,]
  #Build the output df
  cnone.fit.est$pval <- rbind(cnone.lrt.limbic, cnone.lrt.neocortical)

  #Fit models for parietal tissues
  #Fit with beta regression with P_None as the reference
  pnone.fit <- gamlss(x ~ C_None + C_Limbic + C_Neocortical + P_Limbic + P_Neocortical + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Fit reduced model with P_Limbic dropped
  pnone.fit2 <- gamlss(x ~ C_None + C_Limbic + C_Neocortical + P_Neocortical + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Fit reduced model with P_Neocortical dropped
  pnone.fit3 <- gamlss(x ~ C_None + C_Limbic + C_Neocortical + P_Limbic + glial, random = ~1 | randeff, data = x, family = BEZI, trace = F)
  #Calc LRT p-values for parietal tissues relative to P_none
  pnone.lrt.limbic <- try(LR.test(pnone.fit2, pnone.fit, print = F)$p.val)
  if (class(pnone.lrt.limbic) == "try-error") {
    pnone.lrt.limbic <- 1
  }
  #Calc LRT p-values for parietal neocortical tissues relative to P_none
  pnone.lrt.neocortical <- try(LR.test(pnone.fit3, pnone.fit, print = F)$p.val)
  if (class(pnone.lrt.neocortical) == "try-error") {
    pnone.lrt.neocortical <- 1
  }
  #Extract fitted values
  pnone.fit.vals <- fitted(pnone.fit)
  #Fit with beta regression with P_Limbic as reference
  #plimbic.fit <- glmmTMB(x ~ C_None + C_Limbic + C_Neocortical + P_None + P_Neocortical + (1 | randeff), data = x, family = list(family = "beta", link = "logit"), se = T)
  #Calc LRT p-values for cingulate tissues relative to P_none
  #plimbic.lrt <- drop1(plimbic.fit, scope = "P_Neocortical", test = "Chisq")

  #Get the beta estimates for fold-change calcs
  pnone.fit.est <- createCoeftab(pnone.fit)[4:5,]
  #Build the output df
  pnone.fit.est$pval <- rbind(pnone.lrt.limbic, pnone.lrt.neocortical)

  return(list(cingulate = cnone.fit.est, cingulate.fitted = cnone.fit.vals,
              parietal = pnone.fit.est, parietal.fitted = pnone.fit.vals))
}

#Use zibeta regression
#This takes about 2 days to run per assay type (e.g. BS and oxBS)
set.seed(12837)
hydroxymethylcyto.fit.1 <- mclapply(as.data.frame(all_bvals_hmc.t.chnk1), function(x) EPIC.beta.fit(x, glial, randeff), mc.preschedule = F, mc.cores = 32)

save(hydroxymethylcyto.fit.1, file = "/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/all_bvals_hmc.results.chnk1.RData")
