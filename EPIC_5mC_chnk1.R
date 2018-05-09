 library(plyr)
library(dplyr)
library(glmmTMB)
library(parallel)

#Load in a chunk
all_bvals.t.chnk1 <- readRDS("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/all_bvals.t.chnk1.rda")

#New meta data that also includes the glial proportions
metadata.c <- read.table("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/meta_data_c.txt", header = T, check.names = F, stringsAsFactors = F, nrows = 18)
metadata.p <- read.table("/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/meta_data_p.txt", header = T, check.names = F, stringsAsFactors = F, nrows = 18)

#Merge together
metadata <- rbind(metadata.c, metadata.p)

#Get the glial comps
glial <- metadata$glial[which(metadata$Assay == "A")]

#Build the design

brain_region <- as.factor(c(rep("Cingulate", 9), rep("Parietal", 9)))

C_None <- as.factor(c(rep(1, 3), rep(0, 15)))

C_Limbic <- as.factor(c(rep(0, 3), rep(1, 3), rep(0, 12)))

C_Neocortical <- as.factor(c(rep(0, 6), rep(1, 3), rep(0, 9)))

P_None <- as.factor(c(rep(0, 9), rep(1, 3), rep(0, 6)))

P_Limbic <- as.factor(c(rep(0, 12), rep(1, 3), rep(0, 3)))

P_Neocortical <- as.factor(c(rep(0, 15), rep(1, 3)))

#Set randeff
randeff <- all_bvals.t.chnk1$all_arrays

#Create the coefs for plotting with se
createCoeftab <- function(TMB) {
    bTMB <- fixef(TMB)$cond[-1]
    seTMB <- diag(vcov(TMB)$cond)[-1]
    nms <- names(bTMB)
    df <- data.frame(model    = rep(c("glmmTMB"), each = 6),
                     term     = nms,
                     estimate = unname(c(bTMB)))
    df <- transform(df,
                    upper = estimate + sqrt(c(seTMB)),
                    lower = estimate - sqrt(c(seTMB)))
    df
}

#Beta regression with a logit link function
EPIC.beta.fit <- function(x, cellcomp, randeff) {
  x <- as.data.frame(x)
  x$randeff <- randeff
  x$glial <- cellcomp
  
  #Fit models for all tissues
  #Fit with beta regression with C_None as the reference
  cnone.fit <- glmmTMB(x ~ C_Limbic + C_Neocortical + P_None + P_Limbic + P_Neocortical + glial + (1 | randeff), data = x, family = list(family = "beta", link = "logit"), se = T)
  #Calc LRT p-values for cingulate tissues relative to C_none
  cnone.lrt <- drop1(cnone.fit, scope = c("C_Limbic", "C_Neocortical"), test = "Chisq")
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
  #Extract the LRT p-values *UNADJUSTED*
  #None as the ref for cingulate
  cnone.fit.lrtpval <- cnone.lrt$`Pr(>Chi)`[2:3]
  #Build the output df
  cnone.fit.est$pval <- cnone.fit.lrtpval
  
  #Fit models for parietal tissues
  #Fit with beta regression with P_None as the reference
  pnone.fit <- glmmTMB(x ~ C_None + C_Limbic + C_Neocortical + P_Limbic + P_Neocortical + glial + (1 | randeff), data = x, family = list(family = "beta", link = "logit"), se = T)
  #Calc LRT p-values for cingulate tissues relative to P_none
  pnone.lrt <- drop1(pnone.fit, scope = c("P_Limbic", "P_Neocortical"), test = "Chisq")
  #Extract fitted values
  pnone.fit.vals <- fitted(pnone.fit)
  #Fit with beta regression with P_Limbic as reference
  #plimbic.fit <- glmmTMB(x ~ C_None + C_Limbic + C_Neocortical + P_None + P_Neocortical + (1 | randeff), data = x, family = list(family = "beta", link = "logit"), se = T)
  #Calc LRT p-values for cingulate tissues relative to P_none
  #plimbic.lrt <- drop1(plimbic.fit, scope = "P_Neocortical", test = "Chisq")
  
  #Get the beta estimates for fold-change calcs
  pnone.fit.est <- createCoeftab(pnone.fit)[4:5,]
  #Extract the LRT p-values *UNADJUSTED*
  #None as the ref for cingulate
  pnone.fit.lrtpval <- pnone.lrt$`Pr(>Chi)`[2:3]
  #Build the output df
  pnone.fit.est$pval <- pnone.fit.lrtpval

  return(list(cingulate = cnone.fit.est, cingulate.fitted = cnone.fit.vals,
              parietal = pnone.fit.est, parietal.fitted = pnone.fit.vals))
}

#Use beta regression
#This takes about 2 days to run per assay type (e.g. BS and oxBS)
set.seed(12837)
methylcyto.fit.1 <- mclapply(all_bvals.t.chnk1, function(x) EPIC.beta.fit(x, glial, randeff), mc.preschedule = F, mc.cores = 28)

save(methylcyto.fit.1, file = "/secondary/projects/bbc/research/BERA_20171017_EPIC/raw_data/all_bvals_mc.results.chnk1.RData")
