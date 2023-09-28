# This code generates RCBD and PREP trials, populated by lines with breeding values 
# generated from simulated marker effects from actual wheat GBS data, and tests
# asREML GS accuracies on resulting phenotypes. Wrote by AJ Ackerman. 
# Last edited on 

library(R.filesets)
library(MBESS)
library(MASS)
library(matrixcalc)
library(Matrix)
library(caret)
library(magrittr)
library(purrr)
library(furrr)
library(future)
library(asreml4)
library(data.table)
library(tidyverse)

finalResult <- tibble()

for(i in 1:50){
  
  # breedSim, runs simulation a single time, nested in breedSimX ------------
  
  # geno = gbs data, nLoc = number of unique environments, macroGxE = genetic correlation between environments, 
  # microGxE = genetic correlation between tests in RCBD and blocks in PREP, H2 = set heritabilities to test (designates residual error by default)
  # Usable macroGxE are currently .4, .5, and .6. Higher numbers deviate farther from expected correlations due
  # lack of positive definite matrices produced from combination of macro and micro GxE
  
  breedSimX <- function(geno, nLoc, macroGxE, microGxE, H2) {
    
    # Add cohort and plot designation to geno names
    
    uniqueLines <- data.table(germplasmName = row.names(geno))
    
    uniqueLines[, plotDes := "check"][
      str_which(germplasmName, "^18-"), plotDes := "exp"][
        str_which(germplasmName, "^19-"), plotDes := "exp"][
          str_which(germplasmName, "^20-"), plotDes := "exp"][
            str_which(germplasmName, "^2020"), plotDes := "exp"][
              str_which(germplasmName, "^21-"), plotDes := "exp"][
                str_which(germplasmName, "IL2021-"), plotDes := "exp"][
                  str_which(germplasmName, "^18-"), cohort := "S4"][
                    str_which(germplasmName, "^19-"), cohort := "S3"][
                      str_which(germplasmName, "^20-"), cohort := "S2"][
                        str_which(germplasmName, "^2020"), cohort := "S2"][
                          str_which(germplasmName, "^21-"), cohort := "S1"][
                            str_which(germplasmName, "IL2021-"), cohort := "S1"][
                              plotDes == "check", cohort := "check"
                            ]
    
    # Enter entries to design and assign test 
    
    breedSim1 <- rbind(data.table(uniqueLines, design = "RCBD"), data.table(uniqueLines, design = "PREP"))
    breedSim1[c("S1", "S2"), test := "prelim", on = "cohort"][c("S3","S4"), test := "adv", on = "cohort"]
    
    # Assign locations to each line
    
    breedSim2 <- data.table()
    for(i in seq_len(nLoc)){
      breedSim2 <- breedSim1[,location := as.factor(i)] |> 
        rbind(breedSim2)
    }
    
    # Replicate all entries
    
    breedSim3 <- replicate(length(replicateLevels), breedSim2, simplify = FALSE) # Replicate breedSim df for each level of corresponding replication
    
    PREPreps <- map(replicateLevels, ~{
      select(.x, -c(RCBDReplicates,design,test)) %>% 
        drop_na(PREPReplicates) %>% 
        rename(repTotal = PREPReplicates)
    }) # Extract PREP replication from replicateLevels
    
    breedSim3PREPexp <- map2(breedSim3, PREPreps, ~{
      filter(.x, design == "PREP" & plotDes == "exp") %>%
        left_join(.y, by = "cohort") %>% 
        uncount(repTotal, .remove = FALSE, .id = "repNum")
    }) # Replicate experimental entries for PREP
    
    breedSim3PREPcheck <- map2(breedSim3, PREPreps, \(breedSim, reps){
      
      breedSim <- breedSim %>% filter(design == "PREP")
      blocks <- reps %>% slice(-5) %>% pull(repTotal) %>% max()
      checksPerBlock <- reps %>% filter(cohort == "check") %>% pull(repTotal)
      
      breedSim %<>% filter(design == "PREP" & plotDes == "check") %>%
        mutate(blocks = blocks, repTotal = checksPerBlock * blocks, checksPerBlock = checksPerBlock, test = "check") %>% 
        uncount(checksPerBlock, .remove = TRUE, .id = "repNum") %>% 
        uncount(blocks, .remove = TRUE)
      
      
      return(breedSim)
      
    }) # Replicate check entries for PREP
    
    breedSim3PREP <- map2(breedSim3PREPexp, breedSim3PREPcheck, ~{
      bind_rows(.x, .y)
    }) # Combine PREP experimental and check entries into single df
    
    names(breedSim3PREP) <- map_chr(PREPreps, ~{
      slice(.x, -5) %>% 
        pull(repTotal) %>%
        str_c(., collapse = "-")
    }) # Name each list object for corresponding replication level
    
    RCBDReps <- map(replicateLevels, ~{
      select(.x, -PREPReplicates) %>% 
        drop_na(test) %>% 
        rename(repTotal = RCBDReplicates) %>% 
        select(-design)
    }) # Extract RCBD replication from replicateLevels
    
    breedSim3RCBDexp <- map2(breedSim3, RCBDReps, ~{
      filter(.x, design == "RCBD" & cohort != "check") %>%
        left_join(.y, by = c("cohort", "test")) %>% 
        uncount(repTotal, .remove = FALSE, .id = "repNum")
    }) # Replicate experimental entries for RCBD
    
    breedSim3RCBDprelimCheck <- map2(breedSim3, RCBDReps, \(breedSim, reps){
      
      breedSim <- breedSim %>% filter(design == "RCBD")
      blocks <- reps %>% filter(cohort == "S2" | cohort == "S1") %>% pull(repTotal) %>% max(.)
      checksPerBlock <- reps %>% filter(cohort == "check" & test == "prelim") %>% pull(repTotal)
      
      breedSim %<>% filter(design == "RCBD" & cohort == "check" ) %>%
        mutate(test = "prelim") %>% 
        left_join(reps) %>% 
        mutate(blocks = blocks, repTotal = blocks * checksPerBlock, checksPerBlock = checksPerBlock) %>% 
        uncount(blocks, .remove = TRUE, .id = "repNum") %>% 
        uncount(checksPerBlock, .remove = TRUE)
      
      return(breedSim)
      
    }) # Replicate check entries for RCBD
    
    breedSim3RCBDadvCheck <- map2(breedSim3, RCBDReps, \(breedSim, reps){
      
      breedSim <- breedSim %>% filter(design == "RCBD")
      blocks <- reps %>% filter(cohort == "S3" | cohort == "S4") %>% pull(repTotal) %>% max(.)
      prelimBlocks <- reps %>% filter(cohort == "S2" | cohort == "S1") %>% pull(repTotal) %>% max(.)
      checksPerBlock <- reps %>% filter(cohort == "check" & test == "adv") %>% pull(repTotal)
      
      
      breedSim %<>% filter(design == "RCBD" & cohort == "check") %>%
        mutate(test = "adv") %>% 
        left_join(reps) %>% 
        mutate(blocks = blocks, checksPerBlock = checksPerBlock, repTotal = blocks * checksPerBlock) %>% 
        uncount(blocks, .remove = TRUE, .id = "repNum") %>% 
        uncount(checksPerBlock, .remove = TRUE)
      
      return(breedSim)
      
    }) # Replicate check entries for RCBD
    
    
    breedSim3RCBDcheck <- map2(breedSim3RCBDadvCheck, breedSim3RCBDprelimCheck, ~{
      bind_rows(.x,.y)
    })
    
    breedSim3RCBD <- map2(breedSim3RCBDexp, breedSim3RCBDcheck, ~{
      bind_rows(.x,.y)
    })# Combine RCBD experimental and check entries into single df
    
    names(breedSim3RCBD) <- map_chr(RCBDReps, ~{
      slice(.x, 1:4) %>% 
        pull(repTotal) %>%
        str_c(., collapse = "-")
    }) # Name each list object for corresponding replication level
    
    breedSim4 <- imap(c(breedSim3RCBD, breedSim3PREP), ~mutate(.x, repCat = .y)) # Combine RCBD and PREP into single dataframe
    breedSim4 <- rbindlist(breedSim4)
    
    # Assign block to entries for PREP
    
    dfSplit <- split(breedSim4[design == "PREP"], by = c("repCat", "location"))
    
    breedSim5PREP <- map_dfr(dfSplit, \(breedSim) {
      
      breedSim <- setkey(breedSim, cohort)
      
      totalBlocks <- unique(as.integer(breedSim[cohort == "S4", repTotal]))
      S3 <- unique(as.integer(breedSim[cohort == "S3", repTotal]))
      S2 <- unique(as.integer(breedSim[cohort == "S2", repTotal]))
      S1 <- unique(as.integer(breedSim[cohort == "S1", repTotal]))
      
      breedSim["check", block := repNum]    
      
      breedSim["S4", block := repNum, on = "cohort"] 
      
      if (S3 == totalBlocks) {breedSim["S3", block := repNum, on = "cohort"]} else {
        folds <- createFolds(breedSim["S3", germplasmName], k = totalBlocks, list = FALSE)
        breedSim["S3", block := folds]
      }
      
      if (S2 == totalBlocks) {breedSim["S2", block := repNum, on = "cohort"]} else {
        folds <- createFolds(breedSim["S2", germplasmName], k = totalBlocks, list = FALSE)
        breedSim["S2", block := folds]
      }
      
      if (S1 == totalBlocks) {breedSim["S1", block := repNum, on = "cohort"]} else {
        folds <- createFolds(breedSim["S1", germplasmName], k = totalBlocks, list = FALSE)
        breedSim["S1", block := folds]
      }
      
      return(breedSim[, block := as.factor(block)])
      
    })
    
    #Assign blocks for each entry for RCBD
    
    dfSplit <- split(breedSim4[design == "RCBD"], by = "repCat")
    
    breedSim5RCBD <- map_dfr(dfSplit, \(breedSim){
      prelimBlock <- max(breedSim[cohort == "S2",repNum])
      breedSim["prelim", block := repNum, on = "test"][
        "adv", block := repNum + prelimBlock , on = "test"]
      return(breedSim)
    })
    
    # Combine RCBD and PREP
    
    breedSim5 <- rbind(breedSim5PREP, breedSim5RCBD) 
    
    # Create matrices for genetic correlations between environments (GxE)
    
    breedSim6 <- breedSim5[, nLoc := nLoc] 
    
    breedSim6 <- map2_dfr(replicate(length(macroGxE), breedSim5, simplify = FALSE), macroGxE,
                          ~mutate(.x, macroGxE = .y))
    
    breedSim6 <- map2_dfr(replicate(length(microGxE), breedSim6, simplify = FALSE), microGxE,
                          ~mutate(.x, microGxE = .y))
    
    # Generate breeding values with GxE
    
    dfSplit <- breedSim6 %>% group_split(macroGxE, microGxE) 
    
    bv <-  map2_dfr(dfSplit, replicate(length(dfSplit), geno, simplify = FALSE), \(breedSim, geno){
      
      locMat <- matrix(0, ncol = unique(breedSim$nLoc), nrow = unique(breedSim$nLoc))
      diag(locMat) <- 1
      
      testMat <- matrix(unique(breedSim$microGxE), ncol=max(breedSim$repTotal), nrow=max(breedSim$repTotal))
      diag(testMat) <- 1
      
      blockMat <- matrix(unique(breedSim$microGxE), ncol=max(breedSim$repTotal), nrow=max(breedSim$repTotal))
      diag(blockMat) <- 1
      
      amatRCBD <- kronecker(locMat,testMat)
      amatRCBD[which(amatRCBD==0)] <- unique(breedSim$macroGxE)
      if(is.positive.definite(amatRCBD)){print(str_c("RCBD Matrix developed from microGxE: ", microGxE, " and macroGxE: ", macroGxE," is positive definite"))} else {
        amatRCBD <- nearPD(amatRCBD, keepDiag = TRUE, corr =TRUE, do2eigen = FALSE)
        amatRCBD <- as.matrix(amatRCBD[["mat"]])
        print(str_c("RCBD Matrix developed from microGxE: ", microGxE, " and macroGxE: ", macroGxE," is NOT positive definite, corrected using nearPD"))
      }
      
      amatPREP <- kronecker(locMat,blockMat)
      amatPREP[which(amatPREP==0)] <- unique(breedSim$macroGxE)
      if(is.positive.definite(amatPREP)){print(str_c("PREP Matrix developed from microGxE: ", microGxE, " and macroGxE: ", macroGxE," is positive definite"))} else {
        amatPREP <- nearPD(amatPREP, keepDiag = TRUE, corr =TRUE, do2eigen = FALSE)
        amatPREP <- as.matrix(amatPREP[["mat"]])
        print(str_c("PREP Matrix developed from microGxE: ", microGxE, " and macroGxE: ", macroGxE," is NOT positive definite, corrected using nearPD"))
      }
      
      # Make a block matrix for the environment correlations only
      
      amatRCBD <- cor2cov(amatRCBD, rep(1, nrow(amatRCBD)), discrepancy = 1e-5) # Convert correlation matrix to covariance matrix
      amatPREP <- cor2cov(amatPREP, rep(1, nrow(amatPREP)), discrepancy = 1e-5) # Convert correlation matrix to covariance matrix
      
      covRCBD <- mvrnorm(n=ncol(geno), mu = rep(0, ncol(amatRCBD)), Sigma = amatRCBD) 
      covPREP <-  mvrnorm(n=ncol(geno), mu = rep(0, ncol(amatPREP)), Sigma = amatPREP) 
      
      bvRCBD <- data.table((geno %*% covRCBD), keep.rownames = "germplasmName") |>
        setnames(old = 1:ncol(covRCBD)+1, rep(str_c(rep(1:nLoc, each = max(breedSim$repTotal)), "_", rep(c("prelim", "adv", rep("matchBlock", times = (max(breedSim$repTotal) - 2))),  times = length(1:nLoc))))) |>
        melt(id.vars = "germplasmName", variable.name = "loc_test", variable.factor = TRUE, value.name = "bv")
      
      bvRCBD[, c("location","test") := tstrsplit(loc_test, "_")][,design := "RCBD"][, block := NA][, loc_test := NULL][, trueBV := mean(bv), by = germplasmName] 
      
      bvPREP <- data.table((geno %*% covPREP), keep.rownames = "germplasmName") |>
        setnames(old = 1:ncol(covPREP)+1, rep(str_c(rep(seq_len(nLoc), each = max(breedSim6$repTotal)), "_", seq_len(max(breedSim6$repTotal))))) |>
        melt(id.vars = "germplasmName", variable.name = "loc_block", variable.factor = TRUE, value.name = "bv")
      
      bvPREP[, c("location","block") := tstrsplit(loc_block, "_")][,design := "PREP"][, test := NA][, loc_block := NULL][, trueBV := mean(bv), by = germplasmName] 
      
      bv <- rbind(bvPREP, bvRCBD) %>% 
        mutate(nLoc = unique(breedSim$nLoc),                                                                
               macroGxE = unique(breedSim$macroGxE), 
               microGxE = unique(breedSim$microGxE)) %>%
        mutate(across(c(location, block), factor))
      
      return(bv)})
    
    # Use bv matrices to bind to breedSim
    
    breedSim7PREP <- bv[design == "PREP",][, test := NULL][
      breedSim6[design == "PREP",],  on = .(germplasmName, design, location, block, nLoc, microGxE, macroGxE)][
        , trial := str_c("Loc", location, "_", design)]
    breedSim7RCBD <- bv[design == "RCBD",][, block := NULL][
      breedSim6[design == "RCBD",], on = .(germplasmName, design, location, test, nLoc, microGxE, macroGxE)][
        , trial := str_c("Loc", location, "_", test)]
    
    breedSim7 <- rbind(breedSim7PREP, breedSim7RCBD)
    
    # Set heritability
    
    breedSim8 <- map2(replicate(length(H2), breedSim7, simplify = FALSE), H2,
                      ~mutate(.x, heritability = .y) %>%
                        relocate(heritability, .before = trueBV)) # .after within mutate not working?
    breedSim8 <- rbindlist(breedSim8) # rbindlist instead of map2_dfr to avoid internal self ref err from DT
    
    # Generate location and block effect
    
    dfSplit <- split(breedSim8, by = c("heritability", "repCat", "design", "location"))
    
    breedSim9 <- map_dfr(dfSplit, \(breedSim){
      breedSim[, genVarByLoc := var(bv)][
        ,locEff := rnorm(1, 0, sd = sqrt(unique(genVarByLoc)*2))][
          , blockEff := rnorm(1, 0, sd = sqrt(genVarByLoc/4))][
            , errorVar := (((1-heritability)/heritability)*(genVarByLoc))][
              , resErr := rnorm(nrow(breedSim), 0, sd = sqrt(errorVar))] |>
        (\(breedSim){
          breedSim$resErr <-  (breedSim$resErr - mean(breedSim$resErr))/sqrt(var(breedSim$resErr))
          breedSim$resErr <- breedSim$resErr * sqrt(breedSim$errorVar)
          return(breedSim)
        })()# Add block effect based on genetic variance at each location for RCBD
    })
    
    # Adjust vectors used in asREML to factor and nest data table for parallel analysis using furrr
    
    breedSimX <- breedSim9 %>% 
      mutate(phenotype = bv + locEff + blockEff + resErr) %>% 
      mutate(across(c("germplasmName", "location", "block", "test", "cohort"), factor)) %>% 
      group_split(design, repCat, nLoc, macroGxE, microGxE, heritability)
    
    return(breedSimX) # Return result
    
  }
  
  # gsAcc, applies ASReml across breedSimX output ---------------------------
  
  gsAcc <- function(breedSim){
    
    asreml.options(pworkspace = "8gb", ai.sing = TRUE, fail = "soft", threads = 6)
    
    breedSim1Loc <- breedSim %>% mutate(GSset = if_else(test == "prelim" & location != 1,"test", "train")) %>%
      mutate(phenotype = if_else(GSset == "test", NA_real_, phenotype))
    
      GEBV <-  asreml(fixed = phenotype ~ location,
                      random = ~ location:block + at(location):vm(germplasmName, K2),
                      residual = ~ dsum( ~ units | location),
                      workspace = "8gb",
                      data = breedSim1Loc,
                      na.action = na.method(y = "include", x = "include"))
      
      GEBVpred <-  as_tibble(predict.asreml(GEBV, classify = "germplasmName")$pvals)
          
      resultsByEntry <- breedSim1Loc %>% distinct(germplasmName, test, cohort, trueBV, GSset) %>% 
        left_join(GEBVpred, by = "germplasmName") %>%
        filter(GSset == "test" & cohort != "check")
      
      resultsOverall <- tibble(group = "testSet", cor = cor(resultsByEntry$predicted.value, resultsByEntry$trueBV, method = "pearson"), model = "GEBV") %>% 
        bind_cols(distinct(breedSim, design, heritability, nLoc, macroGxE, microGxE, repCat))

      resultsByCohort <- resultsByEntry %>% group_by(cohort) %>% 
        summarize(cor(predicted.value, trueBV, method = "pearson")) %>% 
        rename("cor" = 2, "group" = "cohort") %>% 
        mutate(model = "GEBV") %>% 
        bind_cols(distinct(breedSim, design, heritability, nLoc, macroGxE, microGxE, repCat))
      
      results <- bind_rows(resultsOverall, resultsByCohort)
    
    breedSimNoLoc <- breedSim %>% mutate(GSset = if_else(test == "prelim","test", "train")) %>%
      mutate(phenotype = if_else(GSset == "test", NA_real_, phenotype))
    
    GEBV <-  asreml(fixed = phenotype ~ location,
                    random = ~ location:block + at(location):vm(germplasmName, K2),
                    residual = ~ dsum( ~ units | location),
                    workspace = "8gb",
                    data = breedSimNoLoc,
                    na.action = na.method(y = "include", x = "include"))
    
    GEBVPred <-  as_tibble(predict.asreml(GEBV, classify = "germplasmName")$pvals)
    
    resultsByEntry <- breedSimNoLoc %>% distinct(germplasmName, test, cohort, trueBV, GSset) %>% 
      left_join(GEBVPred, by = "germplasmName") %>%
      filter(GSset == "test" & cohort != "check")
    
    resultsOverall <- tibble(group = "testSet", cor = cor(resultsByEntry$predicted.value, resultsByEntry$trueBV, method = "pearson"), model = "GEBVNoPrelim") %>% 
      bind_cols(distinct(breedSim, design, heritability, nLoc, macroGxE, microGxE, repCat))
    
    resultsByCohort <- resultsByEntry %>% group_by(cohort) %>% 
      summarize(cor(predicted.value, trueBV, method = "pearson")) %>% 
      rename("cor" = 2, "group" = "cohort") %>% 
      mutate(model = "GEBVNoPrelim") %>% 
      bind_cols(distinct(breedSim, design, heritability, nLoc, macroGxE, microGxE, repCat))
    
    results <- bind_rows(resultsOverall, resultsByCohort, results) %>% 
      mutate(time = Sys.time())
    
    return(results)
  }
  
  gsAccSafe <- safely(gsAcc)
  
  # Upload data for input ---------------------------------------------------
  
  geno <- loadRDS("geno.RData")
  cohort <- c("S1","S2","S3","S4", "check", "check", "check")
  design <- c(NA, NA, NA, NA, "PREP", "RCBD", "RCBD")
  test <- c("prelim", "prelim", "adv", "adv", NA, "prelim", "adv")
  checksPerBlockPREP <- 2
  checksPerBlockRCBD <- 2
  replicateLevels <- list(tibble(cohort = cohort, design = design, test = test, PREPReplicates = c(1,1,2,2,checksPerBlockPREP,NA,NA), RCBDReplicates = c(1,1,2,2,NA,checksPerBlockRCBD,checksPerBlockRCBD)), 
                          tibble(cohort = cohort, design = design, test = test, PREPReplicates = c(2,2,2,2,checksPerBlockPREP,NA,NA), RCBDReplicates = c(2,2,2,2,NA,checksPerBlockRCBD,checksPerBlockRCBD)),
                          tibble(cohort = cohort, design = design, test = test, PREPReplicates = c(2,2,3,3,checksPerBlockPREP,NA,NA), RCBDReplicates = c(2,2,3,3,NA,checksPerBlockRCBD,checksPerBlockRCBD)))
  
  # Run pipeline ------------------------------------------------------------
  
  breedSimResult <- breedSimX(geno = geno, nLoc = 5, macroGxE = .5, microGxE =  c(.25,.5,.75,1), H2 = .5)
  
  rm(geno, cohort, design, test, checksPerBlockPREP, checksPerBlockRCBD, replicateLevels, breedSimX)
  
  plan(multicore, workers = 4)
  furrr_options(scheduling = 2L)
  
  K2 <- loadRDS("K2.RData")
  
  allOutput <- breedSimResult %>% future_map(gsAccSafe)
  result <- allOutput %>% future_map("result") %>% compact() %>% rbindlist() %>% mutate(iteration = i)
  errors <- allOutput %>% future_map("error")
  
  write_csv(result, str_c("GEBVoutput/GEBVresults/GEBVresults_", i, ".csv"))
  saveRDS(errors, str_c("GEBVoutput/GEBVerrors/GEBVerrors_", i, ".Rdata"))
  
  finalResult <- bind_rows(finalResult,result)
  
  write_csv(finalResult, "GEBVoutput/GEBVresults_final/breedSimResult_GEBV_all_allIterations.csv")
  
}

