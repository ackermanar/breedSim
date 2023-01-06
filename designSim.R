
setwd("/Users/ackerman/Documents/GitHub/Trial_Design_Project")

library(tidyverse)
library(magrittr)
library(data.table)
library(asreml)
library(asremlPlus)
library(rrBLUP)
library(gaston)
library(naniar)
library(FieldSimR)
library(MASS)
library(rbenchmark)

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

c <- 10
s1 <- 1
s2 <-  2
s3 <-  3
s4 <-  4

nloc <- 5 # number of studies
ntrial <- nloc * 2 # number of trials 
tcorr <- .9 # correlations between trials within locations
ecorr <- .4 # correlations between locations

H <- .4
deltaRes <- (1/6)
deltaTrial <-  (3/6)
deltaLoc <- (2/6)

# breedSim --------------------------------------------------------

breedSim <- function(nloc,tcorr,ecorr, geno, c,s1,s2,s3,s4, H, deltaRes, deltaTrial, deltaLoc){
  
  nTrial <- 2 * nloc
  
# Function 1: convert correlation matrix to correlation matrix

cor2cov_1 <- function(R,S){
  diag(S) %*% R %*% diag(S)
}

# Function 2: determine breeding values by location

  a <- matrix(rep(0,nloc*nloc), ncol=nloc)
  diag(a)<-1
  b <- matrix(c(1,tcorr,tcorr,1), ncol=2)
  amat1 <- kronecker(a,b)
  
# make a block matrix for the environment correlations only
  
  amat1[which(amat1==0)] <- ecorr
  amat1 <- cor2cov_1(amat1, rep(1, nrow(amat1)))   # Function 1: R = correlation matrix, S = Standard deviation
  cov <- mvrnorm(n=ncol(geno), mu =rep(0, ncol(amat1)), Sigma = amat1)
  
  breedVal <- data.table((geno %*% cov), keep.rownames = "germplasmName") %>%
           setnames(old = 1:ncol(cov)+1, rep(str_c("loc", rep(1:nloc, each = 2), "_", rep(c("prlm", "adv"), each= 1:nloc))))
         
# Add cohort label for each year

breedVal2 <- breedVal[, plotDes := "check"][
  str_which(germplasmName, "^18-"), plotDes := "exp"][
    str_which(germplasmName, "^19-"), plotDes := "exp"][
      str_which(germplasmName, "^20-"), plotDes := "exp"][
        str_which(germplasmName, "^2020"), plotDes := "exp"][
          str_which(germplasmName, "^21-"), plotDes := "exp"][
            str_which(germplasmName, "IL2021-"), plotDes := "exp"]

# remove check lines that are not being used in the analysis

checks <- breedVal2[plotDes == "check", ][
  !germplasmName %in% c("Kaskaskia", "PIONEER25R47", "02-18228", "07-19334")] # Insert checks being used for analysis here

breedVal2 <- breedVal2[! germplasmName %in% checks$germplasmName]

breedVal2 <- breedVal2[str_which(germplasmName, "^18-"), cohort := "S4"][
  str_which(germplasmName, "^19-"), cohort := "S3"][
    str_which(germplasmName, "^20-"), cohort := "S2"][
      str_which(germplasmName, "^2020"), cohort := "S2"][
        str_which(germplasmName, "^21-"), cohort := "S1"][
          str_which(germplasmName, "IL2021-"), cohort := "S1"]

# Replicate all entries

checks <- rbind(
  rep(breedVal2[is.na(cohort),][,cohort := "S1"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S2"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S3"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S4"], times = c)
  )

checks <- checks[, unbalancedReps := c][, balancedReps := c][
  ,design := "all"]

exp <- rbind(
  rep(breedVal2[cohort == "S1" & plotDes != "check"], each = s1),
  rep(breedVal2[cohort == "S2" & plotDes != "check"], each = s2),
  rep(breedVal2[cohort == "S3" & plotDes != "check"], each = s3),
  rep(breedVal2[cohort == "S4" & plotDes != "check"], each = s4))
 
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

rcbd <- rbind(rep(exp[cohort == "S1"], times = rcbd1),
              rep(exp[cohort == "S3"], times = rcbd2))
 

rcbd <- rcbd[c("S1"), balancedReps := s2, on = "cohort"][
           c("S3"), balancedReps := s4, on = "cohort"][
             ,design := "balanced"]

breedVal3 <- rbind(checks, exp, rcbd) %>% relocate(design, .after = cohort)

# Seperate BVs by breeding and advanced then melt data frame 
# Look into cleaning this chunk up - messy 

prelim <- breedVal3[cohort == "S1" | cohort == "S2", ][ ,.SD, .SDcols = ! patterns("adv")]
colnames(prelim)<-gsub("_prlm","",colnames(prelim))
adv <- breedVal3[cohort == "S3" | cohort == "S4", ][ ,.SD, .SDcols = ! patterns("prlm")]
colnames(adv)<-gsub("_adv","",colnames(adv))
breedVal4 <- rbind(prelim, adv)
breedVal4 <- breedVal4[c("S1", "S2"), test := "prelim", on = "cohort"][
                       c("S3", "S4"), test := "adv", on = "cohort"]
breedVal4 <- melt.data.table(breedVal4, measure.vars = patterns("loc"), variable.name = "location", value.name = "bv")

# Add trial effect

nonGenVar <- (var(breedVal4$bv)/H) - var(breedVal4$bv)

split <- split(breedVal4, list(breedVal4$location, breedVal4$test))

breedVal5 <- data.table()

trialVar <-  deltaTrial * nonGenVar

for (i in 1:ntrial) {
group <- split[[i]]
group <- group[, trialEffect := rep(rnorm(1, mean = 0, sd = sqrt(trialVar)), times = nrow(group))]
breedVal5 <- rbind(breedVal5, group)
}

breedVal5$trialEffect <-  (breedVal5$trialEffect - mean(breedVal5$trialEffect))/sqrt(var(breedVal5$trialEffect))
breedVal5$trialEffect <- breedVal5$trialEffect * sqrt(trialVar)

# add main location effect and residual error

split <- split(breedVal5, list(breedVal5$location))

breedVal6 <- data.table()
locVar <-  deltaLoc * nonGenVar

for (i in 1:nloc) {
  group <- split[[i]]
  group <- group[, locEffect := rep(rnorm(1, 0, sd = sqrt(locVar)), times = nrow(group))]
  breedVal6 <- rbind(breedVal6, group)
}

breedVal6$locEffect <-  (breedVal6$locEffect - mean(breedVal6$locEffect))/sqrt(var(breedVal6$locEffect))
breedVal6$locEffect <- breedVal6$locEffect * sqrt(locVar)

# add study effect for prep, portion of location effect variance + location effect variance

deltaStudy <- deltaTrial + deltaLoc
studyVar <- nonGenVar * deltaStudy

split <- split(breedVal6, list(breedVal6$location))

breedVal7 <- data.table()

for (i in 1:nloc) {
  group <- split[[i]]
  group <- group[, studyEffect := rep(rnorm(1, 0, sd = sqrt(studyVar)), times = nrow(group))]
  breedVal7 <- rbind(breedVal7, group)
}

breedVal7$studyEffect <-  (breedVal7$studyEffect - mean(breedVal7$studyEffect))/sqrt(var(breedVal7$studyEffect))
breedVal7$studyEffect <- breedVal7$studyEffect * sqrt(studyVar)

# add residual error

resVar <- deltaRes * nonGenVar
breedVal8 <- breedVal7[, residual := rnorm(nrow(breedVal6), mean = 0, sd = sqrt(resVar))]
breedVal8$residual <-  (breedVal8$residual - mean(breedVal8$residual))/sqrt(var(breedVal8$residual))
breedVal8$residual <- breedVal8$residual * sqrt(resVar)

# Create Pheno

breedVal9 <- breedVal8[, rrPheno := bv + trialEffect + locEffect + residual][
  , rcbdPheno := bv + trialEffect + locEffect + residual][
    , prepPheno := bv + studyEffect + residual][
      ,trueBV := mean(bv), by = germplasmName]

return(dfBreedSim <- as_tibble(breedVal9) %>%
  relocate(trueBV, .after = bv) %>% 
  mutate(across(c(1:6), factor)))
  
}

# End function ------------------------------------------------------------

dfBreedSim <- breedSim(geno = geno, nloc = 5,tcorr = .9,ecorr = .4, H = .4, 
                       c = 10,s1 = 1 ,s2 = 2,s3 = 3,s4 = 4, 
                       deltaRes = (1/6), deltaTrial = (1/6), deltaLoc = (2/3))

#### Test Function ####

testRRH <- var(dfBreedSim$bv)/(var(dfBreedSim$bv) + var(dfBreedSim$locEffect) + var(dfBreedSim$trialEffect) + var(dfBreedSim$residual))
print(testRRH)


testPrepH <- var(dfBreedSim$bv)/(var(dfBreedSim$bv) + var(dfBreedSim$residual) + var(dfBreedSim$studyEffect))
print(testPrepH)                               

print(testRCBDH)
                                 
# PHASE ONE: Obtain BLUEs and create df with weights--------------------------------------------------------

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

trialNames <- as.factor(c("T1", "T2", "T3", "T4", "T5"))
traitLevels <- c("pheno1")
blueVal1 <- tibble()
asreml.options(fixgammas=TRUE)

for (i in trialNames) {
  blueValTrait1 <- tibble()
  dfTrial <- breedSim %>% filter(trial == i)
  for (j in traitLevels) {
    dfTrait <- dfTrial %>% filter(trait == j)
    
    blueMod <- asreml(fixed = prepPheno ~ germplasmName,
                      random = ~ cohort,
                      residual = ~ idv(units),
                      data = dfTrait,
                      na.action = na.method(y = "include", x = "include"),
                      workspace="8gb")
    blueMod<- mkConv(blueMod)
    blueVal <- predict.asreml(blueMod, classify = "germplasmName", pworkspace = "8gb")$pvals
    
    blueVal %<>% mutate(trial = i, trait = j, weightTraitByLoc = 1/std.error^2)
    blueValTrait1 <- bind_rows(blueValTrait1, blueVal)
  }
  
  blueVal1 <- bind_rows(blueVal1, blueValTrait1)
  dfTrait <- tibble()
}


blueValAll <- filter(status != "Aliased") %>%
  left_join(breedSim, by = "germplasmName") %>%
  relocate(stage, .before = 6) %>%
  relocate(predicted.value, .after = 8) %>%
  relocate(std.error, .before = 7) %>%
  unite(locStage, c("locName", "stage"), remove = FALSE, sep = "_S") %>%
  mutate(across(c(1:6), factor)) %>%
  rename(traitLevel = trait) %>%
  print(width = Inf)

####end####

# Create relationship matrix ----------------------------------------------

K <- A.mat(geno)
K2 <- nearPD(K)
K2 <- K2$mat

####end####

# Create BLUPs and GEBVs using asreml -------------------------------------

# Create training and validation sets by masking lines at Neoga

# Create masked Neo set and unmasked Neo set

dfNA <- blueValAll %>% filter(locName == "YT_Neo_22") %>%
  filter(germplasmName %in% maskSet) %>%
  mutate(predicted.value = NA_real_) %>%
  mutate(set = "validation")

dfNeoMask <- blueValAll %>% filter(locName != "YT_Neo_22") %>%
  mutate(set = "training") %>%
  bind_rows(dfNA)

dfNeoAll <- pRepLong %>% filter(locName == "YT_Neo_22")

####Run asreml for prediction####

traitLevelIndex <- c("grainYield", "grainTestWeight")
blupValPredAll <- tibble()

for (i in traitLevelIndex) {
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i)
  
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i)
  
  #Training set without stage
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ germplasmName:locName + vm(germplasmName, K2),  #germplasm:locName?
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  
  #Training set with stage
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ germplasmName + stage:locName + vm(germplasmName, K2),
                                residual = ~ idv(units),
                                weights = weightTraitByLoc,
                                family = asr_gaussian(dispersion = 1),
                                data = dfNeoMaskFilter,
                                na.action = na.method(y = "include", x="include"),
                                workspace="8gb")
  blupModNeoMaskStage <-  mkConv(blupModNeoMaskStage)
  
  #Validation set - do not use stage as an effect
  blupModNeo <- asreml(fixed = value ~ 1,
                       random = ~ germplasmName + idv(units),
                       residual = ~ ar1v(colNumber):ar1(rowNumber),
                       data = dfNeoAllFilter,
                       na.action = na.method(y = "include", x="include"),
                       workspace="8gb")
  blupModNeo <-  mkConv(blupModNeo)
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", average = list(locName = NULL), pworkspace = "8gb")[[1]]
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "germplasmName", average = list(locName = NULL, stage = NULL), pworkspace = "8gb")[[1]]
  blupVal <- predict.asreml(blupModNeo, classify = "germplasmName",  pworkspace = "8gb")[[1]]
  
  blupValMask <- tibble(blupValMask) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "FALSE")
  blupVal <- tibble(blupVal) %>%
    mutate(maskingStatus = "completeSet", stageParam = "FALSE")
  blupValMaskStage <- tibble(blupValMaskStage) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "TRUE")
  
  blupValByTrait <- bind_rows(blupValMask, blupValMaskStage, blupVal) %>%
    mutate(trait = i)
  blupValPredAll <- bind_rows(blupValPredAll, blupValByTrait)
  
}

traitLevelIndex2 <- c("plantHeight", "julianDate")

UrbNeo <- c("YT_Urb_22", "YT_Neo_22")

for (i in traitLevelIndex2) {
  
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i, locName %in% UrbNeo)
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i, locName %in% UrbNeo)
  
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ germplasmName + vm(germplasmName, K2),
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ germplasmName + stage:locName + vm(germplasmName, K2),
                                residual = ~ idv(units),
                                weights = weightTraitByLoc,
                                family = asr_gaussian(dispersion = 1),
                                data = dfNeoMaskFilter,
                                na.action = na.method(y = "include", x="include"),
                                workspace="8gb")
  blupModNeoMaskStage <-  mkConv(blupModNeoMaskStage)
  
  #Validation - do not use stage as an effect
  blupModNeo <- asreml(fixed = value ~ 1,
                       random = ~ germplasmName + idv(units),
                       residual = ~ ar1v(colNumber):ar1(rowNumber),
                       data = dfNeoAllFilter,
                       na.action = na.method(y = "include", x="include"),
                       workspace="8gb")
  blupModNeo <-  mkConv(blupModNeo)
  
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", average = list(locName = NULL), pworkspace = "8gb")[[1]]
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "stage:germplasmName", present = list("stage"), average = list(locName = NULL), pworkspace = "8gb")[[1]]
  blupVal <- predict.asreml(blupModNeo, classify = "germplasmName",  pworkspace = "8gb")[[1]]
  
  blupValMask <- tibble(blupValMask) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "FALSE")
  blupVal <- tibble(blupVal) %>%
    mutate(maskingStatus = "completeSet", stageParam = "FALSE")
  blupValMaskStage <- tibble(blupValMaskStage) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "TRUE")
  
  
  blupValByTrait <- bind_rows(blupValMask, blupValMaskStage, blupVal) %>%
    mutate(trait = i)
  blupValPredAll <- bind_rows(blupValPredAll, blupValByTrait)
  
}

####end####



# Calculate accuracy ------------------------------------------------------

traitLevelIndex3 <- c("grainYield", "grainTestWeight", "plantHeight", "julianDate")
finalResultAll <- tibble()
blupValPredAll %<>% unite(maskStage, c("maskingStatus", "stageParam"))
noStage <- c("maskedSet_FALSE", "completeSet_FALSE")
stageVar <- c("maskedSet_TRUE", "completeSet_FALSE")

for (i in traitLevelIndex3) {
  
  blupNoStage <- blupValPredAll %>% filter(trait == i, maskStage %in% noStage & germplasmName %in% maskSet)
  
  completeNoStage <- blupNoStage %>% filter(maskStage == "completeSet_FALSE")
  maskNoStage <- blupNoStage %>% filter(maskStage == "maskedSet_FALSE")
  
  resultNoStage <- inner_join(completeNoStage, maskNoStage, by = "germplasmName") %>%
    rename(predictedValueGEBV = predicted.value.x, predictedValueBLUP = predicted.value.y) %>%
    select(1,2,7)
  
  corrNoStage <- cor(resultNoStage[[2]], resultNoStage[[3]])
  
  blupStage <- blupValPredAll %>% filter(trait == i, maskStage %in% stageVar & germplasmName %in% maskSet)
  
  completeStage <- blupStage %>% filter(maskStage == "completeSet_FALSE")
  maskStage <- blupStage %>% filter(maskStage == "maskedSet_TRUE")
  
  resultStage <- inner_join(completeStage, maskStage, by = "germplasmName") %>%
    rename(predictedValueGEBV = predicted.value.x, predictedValueBLUP = predicted.value.y) %>%
    select(1,2,7)
  
  corrStage <- cor(resultStage[[2]], resultStage[[3]])
  
  finalResult <- tibble(Trait = i, Stage = corrStage, noStage = corrNoStage)
  
  finalResultAll <- bind_rows(finalResultAll, finalResult)
  
}

####end####

stage <- c("1", "2", "3", "4")

setDT(pRep)

for (i in stage) {
  stage[i] <- pRep[stage %in% c(i, "0")]
  stage[i] <- melt(stage[i], c(1:25))
  stageJitter[i] <- jitter(stage[i]$value)
}


df <- read_tsv("YT_2022.tsv") %>%
  select(-c(replicate, `grainMoistureContent%`)) %>%
  replace_with_na(replace = list(25:28 < 0)) %>%
  print(width = Inf)

dfStage <- read_csv("FinalDesignPlan_Aug24.csv") %>%
  select(germplasmName, Stage) %>%
  rename(stage = Stage) %>%
  print(width = Inf)

# Add replication amount for conditional analysis

replicateNum <- df %>% group_by(germplasmName) %>%
  summarize(n()) %>%
  rename(replicateNum = 2)
replicateNum

pRep <- df %>% left_join(replicateNum, by = "germplasmName") %>%
  left_join(dfStage) %>%
  relocate(stage, .before= 3) %>%
  unite(locStage, c("locName", "stage"), remove = FALSE, sep = "_S") %>%
  relocate(replicateNum, .before = 4) %>%
  relocate(entryType, .before = 4) %>%
  arrange(colNumber, rowNumber) %>% #arrange to match spatial error structure col:row
  mutate(across(c(1:26), factor)) %>%
  select(-notes) %>%
  print(width = Inf)

pRepLong <- pivot_longer(pRep, c(27:30), names_to = "traitLevel"  )

# PHASE ONE: Obtain BLUEs and create df with weights--------------------------------------------------------

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

locNames <- c("YT_Neo_22", "YT_Urb_22", "YT_Stj_22", "YT_Stp_22", "YT_Blv_22")
traitLevels <- as.factor(c("grainYield", "grainTestWeight"))
blueVal1 <- tibble()

for (i in locNames) {
  blueValTrait1 <- tibble()
  dfLoc <- pRepLong %>% filter(locName == i)
  for (j in traitLevels) {
    dfTrait <- dfLoc %>% filter(traitLevel == j)
    blueMod <- asreml(fixed = value ~ germplasmName, #environment and genotype as fixed effect
                      random = ~ idv(units), #A nugget effect most often improved model fit.
                      residual = ~ ar1v(colNumber):ar1(rowNumber), #Two-dimensional spatial error models outperform one-dimensional spatial error models.
                      data = dfTrait,
                      na.action = na.method(y = "include", x = "include"),
                      workspace="8gb")
    blueMod<- mkConv(blueMod)
    blueVal <- predict.asreml(blueMod, classify = "germplasmName", pworkspace = "8gb")$pvals
    
    
    blueVal %<>% mutate(locName = i, trait = j, weightTraitByLoc = 1/std.error^2)
    blueValTrait1 <- bind_rows(blueValTrait1, blueVal)
  }
  
  blueVal1 <- bind_rows(blueVal1, blueValTrait1)
  
}

locNames <- c("YT_Neo_22", "YT_Urb_22")
traitLevels <- as.factor(c("plantHeight", "julianDate"))
blueVal2 <- tibble()

for (i in locNames) {
  blueValTrait2 <- tibble()
  dfLoc <- pRepLong %>% filter(locName == i)
  for (j in traitLevels) {
    dfTrait <- dfLoc %>% filter(traitLevel == j)
    blueMod <- asreml(fixed = value ~ germplasmName , #environment and genotype as fixed effect
                      random = ~ idv(units), #A nugget effect most often improved model fit.
                      residual = ~ ar1v(colNumber):ar1(rowNumber), #Two-dimensional spatial error models outperform one-dimensional spatial error models.
                      data = dfTrait,
                      na.action = na.method(y = "include", x = "include"),
                      workspace="8gb")
    blueMod<- mkConv(blueMod)
    blueVal <- predict.asreml(blueMod, classify = "germplasmName", pworkspace = "8gb")$pvals
    
    blueVal %<>% mutate(locName = i, trait = j, weightTraitByLoc = 1/std.error^2)
    blueValTrait2 <- bind_rows(blueValTrait2, blueVal)
  }
  
  blueVal2 <- bind_rows(blueVal2, blueValTrait2)
  
}

blueValAll <- bind_rows(blueVal1, blueVal2) %>%
  filter(status != "Aliased") %>%
  left_join(dfStage, by = "germplasmName") %>%
  relocate(stage, .before = 6) %>%
  relocate(predicted.value, .after = 8) %>%
  relocate(std.error, .before = 7) %>%
  mutate(across(c(1:5), factor)) %>%
  rename(traitLevel = trait) %>%
  print(width = Inf)

####end####

# Predict Neoga -----------------------------------------------------------

####Prep relationship matrix####

geno <- read.vcf("IL_2022_all_regions_samp_filt_fullnames_dedup_imp.vcf.gz")
geno@ped$id <- sub("^.*:", "", geno@ped$id)

#correct the line names

yrseries <- gsub("20", "", as.character(2000:2019))
for(i in 1:length(yrseries)){
  geno@ped$id<- gsub(paste("IL", yrseries[i], sep=""), yrseries[i], geno@ped$id)
}
geno@ped$id <- gsub('IL20', "2020", geno@ped$id)
geno@ped$id <- gsub('IL21', "IL2021", geno@ped$id)
geno@ped$id <- gsub('16LCSDH', "IL16LCSDH", geno@ped$id)
geno@ped$id <- gsub("PIO-25R74", "Pio25R74", geno@ped$id)
geno@ped$id <- gsub("KASKASKIA", "Kaskaskia", geno@ped$id)

dfGeno <- select.inds(geno, id %in% unique(blueValAll$germplasmName))
inter_lines <- intersect(unique(blueValAll$germplasmName), dfGeno@ped$id)
blueValAll %<>% filter(germplasmName %in% inter_lines)
dfGeno <- select.snps(dfGeno, maf > 0.01)
geno<- as.matrix(dfGeno)-1

K <- A.mat(geno)
K2 <- nearPD(K)
K2 <- K2$mat

####end####


# Create BLUPs and GEBVs using asreml -------------------------------------


# Create training and validation sets by masking lines at Neoga

#Find lines present in Neo that are not present in other locations
validNames <- blueValAll %>% filter(locName == "YT_Neo_22") %>%
  select(germplasmName)

validNames <- as.vector(unique(validNames[[1]]))

trainNames <- blueValAll %>% filter(locName != "YT_Neo_22") %>%
  select(germplasmName)

trainNames <- as.vector(unique(trainNames[[1]]))

maskSet <- validNames[!validNames %in% trainNames]

#Create masked Neo set and unmasked Neo set

dfNA <- blueValAll %>% filter(locName == "YT_Neo_22") %>%
  filter(germplasmName %in% maskSet) %>%
  mutate(predicted.value = NA_real_) %>%
  mutate(set = "validation")

dfNeoMask <- blueValAll %>% filter(locName != "YT_Neo_22") %>%
  mutate(set = "training") %>%
  bind_rows(dfNA)

dfNeoAll <- pRepLong %>% filter(locName == "YT_Neo_22")

####Run asreml for prediction####

traitLevelIndex <- c("grainYield", "grainTestWeight")
blupValPredAll <- tibble()

for (i in traitLevelIndex) {
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i)
  
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i)
  
  #Training set without stage
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ locName + vm(germplasmName, K2),
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  
  #Training set with stage
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ locStage + vm(germplasmName, K2),
                                residual = ~ idv(units),
                                weights = weightTraitByLoc,
                                family = asr_gaussian(dispersion = 1),
                                data = dfNeoMaskFilter,
                                na.action = na.method(y = "include", x="include"),
                                workspace="8gb")
  blupModNeoMaskStage <-  mkConv(blupModNeoMaskStage)
  
  #Validation set
  blupModNeo <- asreml(fixed = value ~ 1,
                       random = ~ germplasmName + idv(units),
                       residual = ~ ar1v(colNumber):ar1(rowNumber),
                       data = dfNeoAllFilter,
                       na.action = na.method(y = "include", x="include"),
                       workspace="8gb")
  blupModNeo <-  mkConv(blupModNeo)
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", pworkspace = "8gb")
  blupVal <- predict.asreml(blupModNeo, classify = "germplasmName", pworkspace = "8gb")
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "germplasmName", pworkspace = "8gb")
  
  blupValMask <- tibble(blupValMask) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "FALSE")
  blupVal <- tibble(blupVal) %>%
    mutate(maskingStatus = "completeSet", stageParam = "FALSE")
  blupValMaskStage <- tibble(blupValMaskStage) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "TRUE")
  
  blupValByTrait <- bind_rows(blupValMask, blupVal, blupValMaskStage) %>%
    mutate(trait = i)
  blupValPredAll <- bind_rows(blupValPredAll, blupValByTrait)
  
}

traitLevelIndex2 <- c("plantHeight", "julianDate")

UrbNeo <- c("YT_Urb_22", "YT_Neo_22")

for (i in traitLevelIndex2) {
  
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i, locName %in% UrbNeo)
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i, locName %in% UrbNeo)
  
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ vm(germplasmName, K2),
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ locStage + vm(germplasmName, K2),
                                residual = ~ idv(units),
                                weights = weightTraitByLoc,
                                family = asr_gaussian(dispersion = 1),
                                data = dfNeoMaskFilter,
                                na.action = na.method(y = "include", x="include"),
                                workspace="8gb")
  blupModNeoMaskStage <-  mkConv(blupModNeoMaskStage)
  
  
  blupModNeo <- asreml(fixed = value ~ locName,
                       random = ~ germplasmName + idv(units),
                       residual = ~ ar1v(colNumber):ar1(rowNumber),
                       data = dfNeoAllFilter,
                       na.action = na.method(y = "include", x="include"),
                       workspace="8gb")
  blupModNeo <-  mkConv(blupModNeo)
  
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", pworkspace = "8gb")[[1]]
  blupVal <- predict.asreml(blupModNeo, classify = "germplasmName", pworkspace = "8gb")[[1]]
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "germplasmName", pworkspace = "8gb")[[1]]
  
  blupValMask <- tibble(blupValMask) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "FALSE")
  blupVal <- tibble(blupVal) %>%
    mutate(maskingStatus = "completeSet", stageParam = "FALSE")
  blupValMaskStage <- tibble(blupValMaskStage) %>%
    mutate(maskingStatus = "maskedSet", stageParam = "TRUE")
  
  
  blupValByTrait <- bind_rows(blupValMask, blupVal, blupValMaskStage) %>%
    mutate(trait = i)
  blupValPredAll <- bind_rows(blupValPredAll, blupValByTrait)
  
}

####end####

# Calculate accuracy ------------------------------------------------------

traitLevelIndex3 <- c("grainYield", "grainTestWeight", "plantHeight", "julianDate")
finalResultAll <- tibble()

for (i in traitLevelIndex3) {
  
  blupNoStage <- blupValPredAll %>% filter(trait == i, stageParam == FALSE, germplasmName %in% maskSet)
  
  completeNoStage <- blupNoStage %>% filter(maskingStatus == "completeSet")
  maskNoStage <- blupNoStage %>% filter(maskingStatus == "maskedSet")
  
  resultNoStage <- inner_join(completeNoStage, maskNoStage, by = "germplasmName") %>%
    rename(predictedValueGEBV = predicted.value.x, predictedValueBLUP = predicted.value.y) %>%
    select(1,2,8)
  
  corrNoStage <- cor(resultNoStage[[2]], resultNoStage[[3]])
  
  blupStage <- blupValPredAll %>% filter(trait == i, stageParam == TRUE, germplasmName %in% maskSet)
  
  completeStage <- blupStage %>% filter(maskingStatus == "completeSet")
  maskStage <- blupStage %>% filter(maskingStatus == "maskedSet")
  
  resultStage <- inner_join(completeStage, maskStage, by = "germplasmName") %>%
    rename(predictedValueGEBV = predicted.value.x, predictedValueBLUP = predicted.value.y) %>%
    select(1,2,8)
  
  corrStage <- cor(resultStage[[2]], resultStage[[3]])
  
  finalResult <- tibble(Trait = i, Stage = corrStage, noStage = corrNoStage)
  
  finalResultAll <- bind_rows(finalResultAll, finalResult)
  
}
####end####

stage <- c("1", "2", "3", "4")

setDT(pRep)

for (i in stage) {
  stage[i] <- pRep[stage %in% c(i, "0")]
  stage[i] <- melt(stage[i], c(1:25))
  stageJitter[i] <- jitter(stage[i]$value)
}
