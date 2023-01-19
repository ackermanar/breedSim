
library(tidyverse)
library(magrittr)
library(data.table)
library(asreml)
library(rrBLUP)
library(gaston)
library(MASS)

# Upload gbs data

genoVCF <- read.vcf("IL_2022_all_regions_samp_filt_fullnames_dedup_imp.vcf.gz")
genoVCF@ped$id <- sub("^.*:", "", genoVCF@ped$id)

# Correct the line names

yrseries <- gsub("20", "", as.character(2000:2019))
for(i in 1:length(yrseries)){
  genoVCF@ped$id<- gsub(paste("IL", yrseries[i], sep=""), yrseries[i], genoVCF@ped$id)
}
genoVCF@ped$id <- gsub('IL20', "2020", genoVCF@ped$id)
genoVCF@ped$id <- gsub('IL21', "IL2021", genoVCF@ped$id)
genoVCF@ped$id <- gsub('16LCSDH', "IL16LCSDH", genoVCF@ped$id)
genoVCF@ped$id <- gsub("PIO-25R74", "Pio25R74", genoVCF@ped$id)
genoVCF@ped$id <- gsub("KASKASKIA", "Kaskaskia", genoVCF@ped$id)

dfGeno <- select.snps(genoVCF, maf > 0.01)
geno <- as.matrix(dfGeno)-1

# Manual Entry of variables ----------------------------------------------


# BreedSimX function -----------------------------------

breedSimX <- function(nloc,tcorr,ecorr, geno, s1, s2, s3, s4, H, deltaRes, deltaTrial, deltaLoc){
  
  nTrial <- 2 * nloc
  
  # Nested Fucntion: Convert correlation matrix to correlation matrix
  
  cor2cov_1 <- function(R,S){
    diag(S) %*% R %*% diag(S)
  }
  
  # Determine breeding values by location
  
  a <- matrix(rep(0,nloc*nloc), ncol=nloc)
  diag(a)<- 1
  b <- matrix(c(1,tcorr,tcorr,1), ncol=2)
  amat1 <- kronecker(a,b)
  
  # Make a block matrix for the environment correlations only
  
  amat1[which(amat1== 0)] <- ecorr
  amat1 <- cor2cov_1(amat1, rep(1, nrow(amat1)))   # Function 1: R = correlation matrix, S = Standard deviation
  cov <- mvrnorm(n=ncol(geno), mu =rep(0, ncol(amat1)), Sigma = amat1)
  
  breedSim1 <- data.table((geno %*% cov), keep.rownames = "germplasmName") %>%
    setnames(old = 1:ncol(cov)+1, rep(str_c("loc", rep(1:nloc, each = 2), "_", rep(c("prlm", "adv"), each= 1:nloc))))
  
  # Add cohort label for each year
  
  breedSim2 <- breedSim1[, plotDes := "check"][
    str_which(germplasmName, "^18-"), plotDes := "exp"][
      str_which(germplasmName, "^19-"), plotDes := "exp"][
        str_which(germplasmName, "^20-"), plotDes := "exp"][
          str_which(germplasmName, "^2020"), plotDes := "exp"][
            str_which(germplasmName, "^21-"), plotDes := "exp"][
              str_which(germplasmName, "IL2021-"), plotDes := "exp"]

  # Remove check lines that are not being used in the analysis
  
  checks <- breedSim2[plotDes == "check", ][
    !germplasmName %in% c("Kaskaskia", "PIONEER25R47", "02-18228", "07-19334")] # Insert checks being used for analysis here
  
  breedSim2 <- breedSim2[! germplasmName %in% checks$germplasmName]
  
  breedSim2 <- breedSim2[str_which(germplasmName, "^18-"), cohort := "S4"][
    str_which(germplasmName, "^19-"), cohort := "S3"][
      str_which(germplasmName, "^20-"), cohort := "S2"][
        str_which(germplasmName, "^2020"), cohort := "S2"][
          str_which(germplasmName, "^21-"), cohort := "S1"][
            str_which(germplasmName, "IL2021-"), cohort := "S1"]

  geno<-  geno[breedSim2$germplasmName,] # Filter for lines being used
 
  K <- A.mat(geno)
  K2 <- nearPD(K)
  K2 <- K2$mat
  
  saveRDS(K2, "K2.RData")
  
   # Replicate all entries
  
  checks <- rbind(
    rbindlist(replicate(s1, breedSim2[is.na(cohort),][,cohort := "S1"], simplify = FALSE)),
    rbindlist(replicate(s2, breedSim2[is.na(cohort),][,cohort := "S2"], simplify = FALSE)),
    rbindlist(replicate(s3, breedSim2[is.na(cohort),][,cohort := "S3"], simplify = FALSE)),
    rbindlist(replicate(s4, breedSim2[is.na(cohort),][,cohort := "S4"], simplify = FALSE))
  )
  
  checks <- checks["S1", unbalancedReps:= s1, on = "cohort"]["S1", balancedReps:= s1, on = "cohort"][
    "S2", unbalancedReps:= s2, on = "cohort"]["S2", balancedReps:= s2, on = "cohort"][
      "S3", unbalancedReps:= s3, on = "cohort"]["S3", balancedReps:= s3, on = "cohort"][
        "S4", unbalancedReps:= s4, on = "cohort"]["S4", balancedReps:= s4, on = "cohort"][
          ,design := "both"]
  
  exp <- rbind(
    rbindlist(replicate(s1, breedSim2[cohort == "S1" & plotDes != "check"], simplify = FALSE)),
    rbindlist(replicate(s2, breedSim2[cohort == "S2" & plotDes != "check"], simplify = FALSE)),
    rbindlist(replicate(s3, breedSim2[cohort == "S3" & plotDes != "check"], simplify = FALSE)),
    rbindlist(replicate(s4, breedSim2[cohort == "S4" & plotDes != "check"], simplify = FALSE)))
  
  exp <- exp[, unbalancedReps := s1][
    "S2", unbalancedReps := s2, on = "cohort"][
      "S3", unbalancedReps := s3, on = "cohort"][
        "S4", unbalancedReps := s4, on = "cohort"][
          , design := "unbalanced"][
            c("S1", "S2"), balancedReps := s2, on = "cohort"][
              c("S3", "S4"), balancedReps := s4, on = "cohort"]
  
  rcbd1 <- s2 - s1
  rcbd2 <- s4 - s3
  
  rcbd <- unique(exp, by = "germplasmName")
  
  rcbd <- rbind(rbindlist(replicate(rcbd1, exp[cohort == "S1"], simplify = FALSE)),
                rbindlist(replicate(rcbd2, exp[cohort == "S3"], simplify = FALSE)))
  
  
  rcbd <- rcbd[c("S1"), balancedReps:= s2, on = "cohort"][
    c("S3"), balancedReps:= s4, on = "cohort"][
      ,design := "balanced"]
  
  breedSim3 <- rbind(checks, exp, rcbd) %>% relocate(design, .after = cohort)
  
  # Seperate BVs by breeding and advanced then melt data frame 
  # Look into cleaning this chunk - messy 
  
  prelim <- breedSim3[cohort == "S1" | cohort == "S2", ][ ,.SD, .SDcols = ! patterns("adv")]
  colnames(prelim)<- gsub("_prlm","",colnames(prelim))
  adv <- breedSim3[cohort == "S3" | cohort == "S4", ][ ,.SD, .SDcols = ! patterns("prlm")]
  colnames(adv)<- gsub("_adv","",colnames(adv))
  breedSim4 <- rbind(prelim, adv)
  breedSim4 <- breedSim4[c("S1", "S2"), test := "prelim", on = "cohort"][
    c("S3", "S4"), test := "adv", on = "cohort"]
  breedSim4 <- melt.data.table(breedSim4, measure.vars = patterns("loc"), variable.name = "location", value.name = "bv")
  
  breedSim4<- breedSim4[, trial:= str_c(location, "_", test)] %>% relocate(trial, .before = bv)
  
  # Add trial effect
  
  nonGenVar <- (var(breedSim4$bv)/H) - var(breedSim4$bv)
  
  split <- split(breedSim4, list(breedSim4$location, breedSim4$test))
  
  breedSim5 <- data.table()
  
  trialVar <-  deltaTrial * nonGenVar
  
  for (i in 1:nTrial) {
    group <- split[[i]]
    group <- group[, trialEffect := rep(rnorm(1, mean = 0, sd = sqrt(trialVar)), times = nrow(group))]
    breedSim5 <- rbind(breedSim5, group)
  }
  
  breedSim5$trialEffect <-  (breedSim5$trialEffect - mean(breedSim5$trialEffect))/sqrt(var(breedSim5$trialEffect))
  breedSim5$trialEffect <- breedSim5$trialEffect * sqrt(trialVar)
  
  # Add main location effect and residual error
  
  split <- split(breedSim5, list(breedSim5$location))
  
  breedSim6 <- data.table()
  locVar <-  deltaLoc * nonGenVar
  
  for (i in 1:nloc) {
    group <- split[[i]]
    group <- group[, locEffect := rep(rnorm(1, 0, sd = sqrt(locVar)), times = nrow(group))]
    breedSim6 <- rbind(breedSim6, group)
  }
  
  breedSim6$locEffect <-  (breedSim6$locEffect - mean(breedSim6$locEffect))/sqrt(var(breedSim6$locEffect))
  breedSim6$locEffect <- breedSim6$locEffect * sqrt(locVar)
  
  # Add study effect for prep, portion of location effect variance + location effect variance
  
  deltaStudy <- deltaTrial + deltaLoc
  studyVar <- nonGenVar * deltaStudy
  
  split <- split(breedSim6, list(breedSim6$location))
  
  breedSim7 <- data.table()
  
  for (i in 1:nloc) {
    group <- split[[i]]
    group <- group[, studyEffect := rep(rnorm(1, 0, sd = sqrt(studyVar)), times = nrow(group))]
    breedSim7 <- rbind(breedSim7, group)
  }
  
  breedSim7$studyEffect <-  (breedSim7$studyEffect - mean(breedSim7$studyEffect))/sqrt(var(breedSim7$studyEffect))
  breedSim7$studyEffect <- breedSim7$studyEffect * sqrt(studyVar)
  
  # add residual error
  
  resVar <- deltaRes * nonGenVar
  breedSim8 <- breedSim7[, residual := rnorm(nrow(breedSim6), mean = 0, sd = sqrt(resVar))]
  breedSim8$residual <-  (breedSim8$residual - mean(breedSim8$residual))/sqrt(var(breedSim8$residual))
  breedSim8$residual <- breedSim8$residual * sqrt(resVar)
  
  # Create Pheno
  
  breedSim9 <- breedSim8[, rrPheno := bv + trialEffect + locEffect + residual][
    , rcbdPheno := bv + trialEffect + locEffect + residual][
      , prepPheno := bv + studyEffect + residual][
        ,trueBV := mean(bv), by = germplasmName]
  
  return(breedSim <- breedSim9 %>% relocate(trueBV, .before = bv) %>% 
           mutate(across(c(1:9), factor)))
}

breedSim <- breedSimX(geno = geno, nloc = 5,tcorr = .9,ecorr = .4, H = .4, 
                      s1 = 1 ,s2 = 2,s3 = 3,s4 = 4, 
                      deltaRes = (1/6), deltaTrial = (1/6), deltaLoc = (2/3)) 
version<- "1"

saveRDS(breedSim, str_c("breedSimInput_V", version, ".RData"))

# Perform two phase genomic selection -------------------------------------

# two phase model, BLUE by location then BLUP

asreml.options(workspace = "8gb", pworkspace = "8gb")

rm_dups <-  unique(breedSim[, c(1:3,10)], by = "germplasmName")

# predict for partial replication (unbalanced)

breedSim_pRep <-  breedSim[design ==  "unbalanced"] # Subset dataset using "design" column

BLUE_pRep <-  asreml(fixed = prepPheno ~ germplasmName,
                     random = ~ location + location:germplasmName,
                     residual = ~ idv(units),
                     data = breedSim_pRep,
                     na.action = na.method(y = "omit", x = "omit"))

pred_BLUE_pRep <-  as.data.table(predict.asreml(BLUE_pRep, classify = "germplasmName")$pvals)

BLUE_pRep<- rm_dups[pred_BLUE_pRep, on = .(germplasmName)]
BLUE_pRep[, locWeight :=  1/std.error^2]

gBLUP_pRep <- asreml(fixed = predicted.value ~ 1,
                    random = ~ vm(germplasmName, K2),
                    weights = locWeight,
                    family = asr_gaussian(dispersion = 1),
                    residual = ~ idv(units),
                    data = BLUE_pRep, 
                    na.action = na.method(y = "omit", x = "omit"))

pred_gBLUP_pRep<- as.data.table(predict.asreml(gBLUPE_pRep, classify = "germplasmName")$pvals)
saveRDS(pred_gBLUP_pRep, str_c("gBLUP_pRep_V", version, ".RData"))

gBLUP_pRep->  rm_dups[pred_gBLUP_pRep, on = .(germplasmName)]
cor_pRep<- gBLUP_pRep %>% cor(gBLUP_pREP$trueBV, gBLUP_pREP$predicted.value, method = "pearson")
setnames(gBLUP_pRep, "predicted.value", "pRep_gBLUP")



# predict for restricted randomization (unbalanced)

breedSim_RR_exp<- breedSim[design== "unbalanced"] # Filter dataset for restricted randomization
breedSim_RR_checks<- breedSim[plotDes==  "check" & cohort== "S2" | cohort==  "S4"]
breedSim_RR<- rbind(breedSim_RR_exp, breedSim_RR_checks)

BLUE_RR <- asreml(fixed = rrPheno ~ germplasmName,
                  random = ~ test + location + trial:location:germplasmName,
                  residual = ~ idv(units),
                  data = breedSim_RR,
                  na.action = na.method(y = "omit", x = "omit"))

pred_BLUE_RR <-  as.data.table(predict.asreml(BLUE_RR, classify = "germplasmName")$pvals)

BLUE_RR<-  rm_dups[pred_BLUE_RR, on = .(germplasmName)]
BLUE_RR[, locWeight :=  1/std.error^2]

gBLUP_RR <- asreml(fixed = predicted.value ~ 1,
                  random = ~ vm(germplasmName, K2),
                  weights = locWeight,
                  family = asr_gaussian(dispersion = 1),
                  residual = ~ idv(units),
                  data = BLUE_RR, 
                  na.action = na.method(y = "omit", x = "omit"))

pred_gBLUP_RR<- as.data.table(predict.asreml(gBLUP_RR, classify = "germplasmName")$pvals)
saveRDS(pred_gBLUP_RR, str_c("gBLUP_RR_V", version, ".RData"))

gBLUP_RR<- rm_dups[pred_gBLUP_RR, on = .(germplasmName)]
cor_RR<- cor(gBLUP_RR$trueBV, gBLUP_RR$predicted.value, method = "pearson")
setnames(gBLUP_RR, "predicted.value", "RR_gBLUP")

# Predict for randomized complete block design (balanced) 

breedSim_RCBD<- breedSim[cohort != "S1" | cohort != "S3" & design == "all"] # Filter dataset for RCBD, do not filter for balanced designs

BLUE_RCBD <- asreml(fixed = rcbdPheno ~ germplasmName,
                    random = ~ trial + location + trial:location:germplasmName,
                    residual = ~ idv(units),
                    data = breedSim_RCBD,
                    na.action = na.method(y = "omit", x = "omit"))

pred_BLUE_RCBD <-  as.data.table(predict.asreml(BLUE_RCBD, classify = "germplasmName")$pvals)

BLUE_RCBD<- rm_dups[pred_BLUE_RCBD, on = .(germplasmName)]
BLUE_RCBD[, locWeight :=  1/std.error^2]

gBLUP_RR <- asreml(fixed = predicted.value ~ 1,
                  random = ~ vm(germplasmName, K2),
                  weights = locWeight,
                  family = asr_gaussian(dispersion = 1),
                  residual = ~ idv(units),
                  data = BLUE_RCBD, 
                  na.action = na.method(y = "omit", x = "omit"))

pred_gBLUP_RCBD<- as.data.table(predict.asreml(gBLUP_RCBD, classify = "germplasmName")$pvals)
saveRDS(pred_gBLUP_RCBD, str_c("pred_RCBD_V", version, ".RData"))

gBLUP_RCBD<- rm_dups[pred_gBLUP_pRCBD, on = .(germplasmName)]
cor_RCBD<- cor(gBLUP_RCBD$trueBV, gBLUP_RCBD$predicted.value, method = "pearson")
setnames(gBLUP_RCBD, "predicted.value", "RCBD_gBLUP")

# Create csv for this version

pRep<- gBLUP_pRep[, .(germplasmName, pRep_gBLUP)]
RR<- gBLUP_RR[, .(germplasmName, RR_gBLUP)]
RCBD<- gBLUP_RCBD[, .(germplasmName, RCBD_gBLUP)]

results_df<-  rm_dups %>% left_join(pRep, by = "germplasmName") %>% 
  left_join(RR, by = "germplasmName") %>% 
  left_join(RCBD, by = "germplasmName")

write.csv(results_df, str_c("breedSimResults_V",version,".csv"))

results_cor<- data.table(pRep = cor_pRep, RR = cor_RR, RCBD = cor = cor_RCBD)

write.csv(results_cor, str_c("GS_ACC_Pearson_V",version,".csv"))


