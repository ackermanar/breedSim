
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


# Upload Geno Data --------------------------------------------------------

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

dfGeno <- select.snps(geno, maf > 0.01)
dfGeno2<- as.matrix(dfGeno)-1
dfGenoT <- t(dfGeno2)

# Simulate Marker Matrix --------------------------------------------------

  # Function 1: convert correlation matrix to correlation matrix

cor2cov_1 <- function(R,S){
  diag(S) %*% R %*% diag(S)
}

  # Function 2: determine breeding values by location

breedValbyLoc <- function(n,t,l){
  
  nloc <- n
  tcorr <- t
  ecorr <- l
  
  a <- matrix(rep(0,nloc*nloc), ncol=nloc)
  diag(a)<-1
  b <- matrix(c(1,tcorr,tcorr,1), ncol=2)
  amat1 <- kronecker(a,b)
  
  #make a block matrix for the environment correlations only
  
  amat1[which(amat1==0)] <- ecorr
  amat1 <- cor2cov_1(amat1, rep(1, nrow(amat1)))   # Function 1: R = correlation matrix, S = Standard deviation
  cov <- mvrnorm(n=ncol(dfGeno2), mu =rep(0, ncol(amat1)), Sigma = amat1)
  
  c <- n*2
  
  return(data.table((dfGeno2 %*% cov), keep.rownames = "germplasmName") %>%
           setnames(old = 1:ncol(cov)+1, rep(str_c("study", rep(1:n, each = 2), "_", rep(c("prlm", "adv"), each= 1:n)))))
         }

breedVal <- breedValbyLoc(5,.9,.4) #Function 2: n = number of studies, t = correlation between trials, l = correlation between locations

# Add cohort label for each year ------------------------------------------

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

# Replicate all entries -------------------------------------------

replication <- function(c,s1,s2,s3,s4){

  checks <- rbindlist(list(
  rep(breedVal2[is.na(cohort),][,cohort := "S1"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S2"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S3"], times = c),
  rep(breedVal2[is.na(cohort),][,cohort := "S4"], times = c)
  ))


 checks <- checks[, reps := c]

 exp <- rbindlist(
   list(rep(breedVal2[cohort == "S1" & plotDes != "check"], each = s1),
        rep(breedVal2[cohort == "S2" & plotDes != "check"], each = s2),
        rep(breedVal2[cohort == "S3" & plotDes != "check"], each = s3),
        rep(breedVal2[cohort == "S4" & plotDes != "check"], each = s4)))
 
exp <- exp[, reps := s1][
   "S2", reps := s2, on = "cohort"][
     "S3", reps := s3, on = "cohort"][
       "S4", reps := s4, on = "cohort"]

 return(rbind(checks, exp))
 
}

# Enter amount of replication for checks per cohort, and S1, S2, S3, S4 entries
breedVal3 <- replication(10, 1, 2, 3, 4) 

# Seperate BVs by breeding and advanced and melt data frame ---------------


prelim <- breedVal4[cohort == "S1" | cohort == "S2", ][ ,.SD, .SDcols = ! patterns("adv")]
colnames(prelim)<-gsub("_prlm","",colnames(prelim))
adv <- breedVal4[cohort == "S3" | cohort == "S4", ][ ,.SD, .SDcols = ! patterns("prlm")]
colnames(adv)<-gsub("_adv","",colnames(adv))
breedVal5 <- rbind(prelim, adv)
breedVal5 <- breedVal5[c("S1", "S2"), test := "prelim", on = "cohort"][
                       c("S3", "S4"), test := "adv", on = "cohort"]
breedVal5 <- melt.data.table(breedVal5, measure.vars = patterns("study"), variable.name = "study", value.name = "bv")

# Add error for each entry ------------------------------------------------

breedVal6 <- breedVal4[ , prepResid := rnorm(length(breedVal4[[1]]),
                                             mean = 0, sd = sqrt(var(breedVal$bv)))]

cohortCount <- breedVal6 %>% count(cohort)
prelim <- sum(cohortCount[cohort == "S1" | cohort == "S2", 2])
adv <- sum(cohortCount[cohort == "S3" | cohort == "S4", 2])

breedVal6 <- breedVal6[  , rcbdResid := rnorm(length(breedVal6[[1]]), mean = 0, sd = sqrt(var(breedVal$bv)))][
  c("S1", "S2"), rcbdResid := rnorm(prelim,
                                    mean = 1, sd = sqrt(var(breedVal$bv))), on = "cohort"][
                                      c("S3", "S4"), rcbdResid := rnorm(adv,
                                                                        mean = 2, sd = sqrt(var(breedVal$bv))), on = "cohort"]

# Create Pheno ------------------------------------------------------------

breedVal7 <- breedVal6[, prepPheno := bv + trialResid + prepResid][, rcbdPheno := bv + trialResid + rcbdResid]

breedSim <- as_tibble(breedVal7) %>%
  mutate(across(c(1:6), factor))

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
  unite(studyStage, c("studyName", "stage"), remove = FALSE, sep = "_S") %>%
  mutate(across(c(1:6), factor)) %>%
  rename(traitLevel = trait) %>%
  print(width = Inf)

####end####

# Create relationship matrix ----------------------------------------------

K <- A.mat(dfGeno2)
K2 <- nearPD(K)
K2 <- K2$mat

####end####

# Create BLUPs and GEBVs using asreml -------------------------------------

# Create training and validation sets by masking lines at Neoga

# Create masked Neo set and unmasked Neo set

dfNA <- blueValAll %>% filter(studyName == "YT_Neo_22") %>%
  filter(germplasmName %in% maskSet) %>%
  mutate(predicted.value = NA_real_) %>%
  mutate(set = "validation")

dfNeoMask <- blueValAll %>% filter(studyName != "YT_Neo_22") %>%
  mutate(set = "training") %>%
  bind_rows(dfNA)

dfNeoAll <- pRepLong %>% filter(studyName == "YT_Neo_22")

####Run asreml for prediction####

traitLevelIndex <- c("grainYield", "grainTestWeight")
blupValPredAll <- tibble()

for (i in traitLevelIndex) {
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i)
  
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i)
  
  #Training set without stage
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ germplasmName:studyName + vm(germplasmName, K2),  #germplasm:studyName?
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  
  #Training set with stage
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ germplasmName + stage:studyName + vm(germplasmName, K2),
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
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", average = list(studyName = NULL), pworkspace = "8gb")[[1]]
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "germplasmName", average = list(studyName = NULL, stage = NULL), pworkspace = "8gb")[[1]]
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
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i, studyName %in% UrbNeo)
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i, studyName %in% UrbNeo)
  
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
                                random = ~ germplasmName + stage:studyName + vm(germplasmName, K2),
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
  
  
  blupValMask <- predict.asreml(blupModNeoMask, classify = "germplasmName", average = list(studyName = NULL), pworkspace = "8gb")[[1]]
  blupValMaskStage <- predict.asreml(blupModNeoMaskStage, classify = "stage:germplasmName", present = list("stage"), average = list(studyName = NULL), pworkspace = "8gb")[[1]]
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
  unite(studyStage, c("studyName", "stage"), remove = FALSE, sep = "_S") %>%
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

studyNames <- c("YT_Neo_22", "YT_Urb_22", "YT_Stj_22", "YT_Stp_22", "YT_Blv_22")
traitLevels <- as.factor(c("grainYield", "grainTestWeight"))
blueVal1 <- tibble()

for (i in studyNames) {
  blueValTrait1 <- tibble()
  dfLoc <- pRepLong %>% filter(studyName == i)
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
    
    
    blueVal %<>% mutate(studyName = i, trait = j, weightTraitByLoc = 1/std.error^2)
    blueValTrait1 <- bind_rows(blueValTrait1, blueVal)
  }
  
  blueVal1 <- bind_rows(blueVal1, blueValTrait1)
  
}

studyNames <- c("YT_Neo_22", "YT_Urb_22")
traitLevels <- as.factor(c("plantHeight", "julianDate"))
blueVal2 <- tibble()

for (i in studyNames) {
  blueValTrait2 <- tibble()
  dfLoc <- pRepLong %>% filter(studyName == i)
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
    
    blueVal %<>% mutate(studyName = i, trait = j, weightTraitByLoc = 1/std.error^2)
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
dfGeno2<- as.matrix(dfGeno)-1

K <- A.mat(dfGeno2)
K2 <- nearPD(K)
K2 <- K2$mat

####end####


# Create BLUPs and GEBVs using asreml -------------------------------------


# Create training and validation sets by masking lines at Neoga

#Find lines present in Neo that are not present in other locations
validNames <- blueValAll %>% filter(studyName == "YT_Neo_22") %>%
  select(germplasmName)

validNames <- as.vector(unique(validNames[[1]]))

trainNames <- blueValAll %>% filter(studyName != "YT_Neo_22") %>%
  select(germplasmName)

trainNames <- as.vector(unique(trainNames[[1]]))

maskSet <- validNames[!validNames %in% trainNames]

#Create masked Neo set and unmasked Neo set

dfNA <- blueValAll %>% filter(studyName == "YT_Neo_22") %>%
  filter(germplasmName %in% maskSet) %>%
  mutate(predicted.value = NA_real_) %>%
  mutate(set = "validation")

dfNeoMask <- blueValAll %>% filter(studyName != "YT_Neo_22") %>%
  mutate(set = "training") %>%
  bind_rows(dfNA)

dfNeoAll <- pRepLong %>% filter(studyName == "YT_Neo_22")

####Run asreml for prediction####

traitLevelIndex <- c("grainYield", "grainTestWeight")
blupValPredAll <- tibble()

for (i in traitLevelIndex) {
  blupValByTrait <- tibble()
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i)
  
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i)
  
  #Training set without stage
  blupModNeoMask <- asreml(fixed = predicted.value ~ 1,
                           random = ~ studyName + vm(germplasmName, K2),
                           weights = weightTraitByLoc,
                           family = asr_gaussian(dispersion = 1),
                           residual = ~ idv(units),
                           data = dfNeoMaskFilter,
                           na.action = na.method(y = "include", x="include"),
                           workspace="8gb")
  blupModNeoMask <-  mkConv(blupModNeoMask)
  
  
  #Training set with stage
  blupModNeoMaskStage <- asreml(fixed = predicted.value ~ 1,
                                random = ~ studyStage + vm(germplasmName, K2),
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
  
  dfNeoMaskFilter <- dfNeoMask %>%  filter(traitLevel == i, studyName %in% UrbNeo)
  dfNeoAllFilter <- dfNeoAll %>%  filter(traitLevel == i, studyName %in% UrbNeo)
  
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
                                random = ~ studyStage + vm(germplasmName, K2),
                                residual = ~ idv(units),
                                weights = weightTraitByLoc,
                                family = asr_gaussian(dispersion = 1),
                                data = dfNeoMaskFilter,
                                na.action = na.method(y = "include", x="include"),
                                workspace="8gb")
  blupModNeoMaskStage <-  mkConv(blupModNeoMaskStage)
  
  
  blupModNeo <- asreml(fixed = value ~ studyName,
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
