library(lavaan)

# Set the working directory
setwd("~/GitHub/OrdinalSim/Functions")

# Source the relevant functions
source("MMG-SEM.R")
source("E_Step.R")

setwd("~/GitHub/OrdinalSim")
source("DataGeneration.R")
source("evaluation.R")
source("evaluation_MM.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus            <- c(2, 4)                 # Number of clusters
ngroups          <- c(12)                   # Number of groups
coeff            <- c(0.2, 0.3, 0.4)        # Initial regression parameters
N_g              <- c(50, 100, 200)         # Sample size per groups
balance          <- c("bal", "unb")         # Cluster size
NonInvThreshSize <- c(0, 0.2, 0.4)          # Threshold non-invariance size
NonInvLoadSize   <- c(0.2, 0.6)             # Loading non-invariance size
c                <- c(2, 4, 5)              # Number of categories

model <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
    
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'
S1 <- '
    # factor loadings
    F1 =~ x1 + x2 + x3 + x4 + x5
    F2 =~ z1 + z2 + z3 + z4 + z5
    F3 =~ m1 + m2 + m3 + m4 + m5
    F4 =~ y1 + y2 + y3 + y4 + y5
'

S2 <- '
    # Regression parameters
    F4 ~ F1 + F3
    F3 ~ F1 + F2
'

# Get design matrix
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, NonInvThreshSize, c, NonInvLoadSize, model) 
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "NonInvThreshSize", "c", "NonInvLoadSize", "model")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, c, NonInvLoadSize, NonInvThreshSize) 

# Functions for the simulation
# First, to avoid stopping due to errors, create a function with data generation and MMGSEM
# Errors come from non positive definite cov matrices. This code allows the re-sample
genDat_analysis <- function(seed, RowDesign, k, NonInv){
  
  # Set seed per design condition (row) and replication (K)
  set.seed(seed)
  # Generate data
  #SimData <- do.call(what = DataGeneration, args = design[RowDesign, ])$SimData
  SimData <- DataGeneration(model            = design[RowDesign, "model"], 
                            nclus            = design[RowDesign, "nclus"], 
                            ngroups          = design[RowDesign, "ngroups"], 
                            reg_coeff        = design[RowDesign, "coeff"], 
                            N_g              = design[RowDesign, "N_g"], 
                            balance          = design[RowDesign, "balance"], 
                            NonInvThreshSize = design[RowDesign, "NonInvThreshSize"],
                            NonInvSize       = design[RowDesign, "NonInvLoadSize"],
                            c                = design[RowDesign, "c"])
  
  # Check that there are no empty categories
  SimData$SimData <- as.data.frame(SimData$SimData)
  for(g in 1:design[RowDesign, "ngroups"]){
    this_g <- SimData$SimData[SimData$SimData$group == g, ]
    categories <- apply(X = this_g, MARGIN = 2, FUN = unique)
    for (col in 1:(length(categories) - 1)) {
      if(!c(all(1:design[RowDesign, "c"] %in% categories[[col]]))){
        stop()
      }
    }
  }
  
  # Save data
  save(SimData, file = paste("Data/Data", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(seed, file = paste("Seed/Seed", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  # Get the non-invariant threshold parameters for the syntax
  fake_cfa <- cfa(S1, ordered = T, do.fit = F, data = SimData$SimData)
  PRT <- partable(fake_cfa)
  PRT$full <- paste0(PRT$lhs, PRT$op, PRT$rhs)
  NoninvThresh <- PRT$full[which(PRT$op == "|" & PRT$lhs %in% c("x2", "x3",
                                                                "m2", "m3",
                                                                "y2", "y3",
                                                                "z2", "z3"))]
  
  # Define non-invariances
  NonInv <- c("F1 =~ x2", "F1 =~ x3",
              "F2 =~ z2", "F2 =~ z3",
              "F3 =~ m2", "F3 =~ m3",
              "F4 =~ y2", "F4 =~ y3",
              NoninvThresh)
  # browser()
  # Fit model, including non-invariances
  ctime.cat <- system.time(
    fit.cat <- MMGSEM(dat = SimData$SimData, S1 = S1, S2 = S2, group = "group", 
                      nclus = design[RowDesign, "nclus"], seed = seed,
                      nstarts = 20, allG = T, est_method = "local", ordered = T, 
                      group.partial = NonInv, group.equal = c("loadings", "thresholds"))
  )
  
  ctime.con <- system.time(
    fit.con <- MMGSEM(dat = SimData$SimData, S1 = S1, S2 = S2, group = "group", 
                      nclus = design[RowDesign, "nclus"], seed = seed,
                      nstarts = 20, allG = T, est_method = "local", ordered = F, 
                      group.partial = NonInv, group.equal = c("loadings"))
  )
  
  # Fit model, ignoring non-invariances
  ctime.ign.cat <- system.time(
    fit.ign.cat <- MMGSEM(dat = SimData$SimData, S1 = S1, S2 = S2, group = "group", 
                          nclus = design[RowDesign, "nclus"], seed = seed,
                          nstarts = 20, allG = T, est_method = "local", ordered = T, 
                          group.partial = NULL, group.equal = c("loadings", "thresholds"))
  )
  
  ctime.ign.con <- system.time(
    fit.ign.con <- MMGSEM(dat = SimData$SimData, S1 = S1, S2 = S2, group = "group", 
                          nclus = design[RowDesign, "nclus"], seed = seed,
                          nstarts = 20, allG = T, est_method = "local", ordered = F, 
                          group.partial = NULL, group.equal = c("loadings"))
  )
  
  results <- list(fit.con = fit.con, fit.cat = fit.cat, 
                  fit.ign.con = fit.ign.con, fit.ign.cat = fit.ign.cat,
                  SimData = SimData)
  
  ### Save computation times
  save(ctime.con,     file = paste("Times/ConTime", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(ctime.ign.con, file = paste("Times/ConTimeIgn", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  save(ctime.cat,     file = paste("Times/CatTime", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(ctime.ign.cat, file = paste("Times/CatTimeIgn", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  # Save results if necessary
  save(fit.con,     file = paste("Fit/ConFit", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(fit.ign.con, file = paste("Fit/ConFitIgn", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  save(fit.cat,     file = paste("Fit/CatFit", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  save(fit.ign.cat, file = paste("Fit/CatFitIgn", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
  
  # If everything goes right, return results
  return(results)
}

# Create SAFE function (in case of errors)
genDat_analysis <- purrr::safely(.f = genDat_analysis, otherwise = NULL)

# Main simulation function
do_sim <- function(RowDesign){
  cat("\n", "Condition", RowDesign, "out of", nrow(design), "\n")
  # Create the original clustering matrix for comparison below
  original <- create_original(balance = design[RowDesign, "balance"], 
                              ngroups = design[RowDesign, "ngroups"], 
                              nclus   = design[RowDesign, "nclus"])
  # return(original)
  # Create matrices to store results
  # Including non-invariances
  # 10 columns for: ARI and RMSEA * 2 (continuous and ordinal) 
  ResultsRow <- matrix(data = NA, nrow = (K), ncol = 16)
  colnames(ResultsRow) <- c("ARI.con", "CC.con", "RMSE_B1.con", "RMSE_B2.con", "RMSE_B3.con", "RMSE_B4.con", "exo_mean.con", "cov_mean.con", 
                            "ARI.cat", "CC.cat", "RMSE_B1.cat", "RMSE_B2.cat", "RMSE_B3.cat", "RMSE_B4.cat", "exo_mean.cat", "cov_mean.cat")
  
  # 4 columns: RMSE (lambda and theta) * 2 (continuous and ordinal) 
  ResultsRow_MM <- matrix(data = NA, nrow = (K), ncol = 6)
  colnames(ResultsRow_MM) <- c("Tuck_lambda.con", "RMSE_lambda.con", "RMSE_theta.con", 
                               "Tuck_lambda.cat", "RMSE_lambda.cat", "RMSE_theta.cat")
  
  # Ignoring non-invariances
  # 10 columns for: ARI and RMSEA * 2 (continuous and ordinal) 
  ResultsRow.ign <- matrix(data = NA, nrow = (K), ncol = 16)
  colnames(ResultsRow.ign) <- c("ARI.con", "CC.con", "RMSE_B1.con", "RMSE_B2.con", "RMSE_B3.con", "RMSE_B4.con", "exo_mean.con", "cov_mean.con", 
                                "ARI.cat", "CC.cat", "RMSE_B1.cat", "RMSE_B2.cat", "RMSE_B3.cat", "RMSE_B4.cat", "exo_mean.cat", "cov_mean.cat")
  
  # 4 columns: RMSE (lambda and theta) * 2 (continuous and ordinal) 
  ResultsRow_MM.ign <- matrix(data = NA, nrow = (K), ncol = 6)
  colnames(ResultsRow_MM.ign) <- c("Tuck_lambda.con", "RMSE_lambda.con", "RMSE_theta.con",
                                   "Tuck_lambda.cat", "RMSE_lambda.cat", "RMSE_theta.cat")
  
  
  for(k in 1:K){
    # print(paste("Replication", k, "out of", K))
    # # browser()
    # Code to re-sample in case the covariance matrix is non positive definite
    attempts <- 10
    for(j in 1:attempts){
      # Seed will change if there is an error
      results <- genDat_analysis(seed = (RowDesign * k * j), RowDesign = RowDesign, k = k)
      test    <- results$result 
      if(!is.null(test)){
        # If there was no error, break the loop and continue
        break
      }
    }
    
    results <- results$result
    
    # ---------------------------------------------------------------
    # Evaluate the results
    if(is.null(test)){
      ResultsRow[k, 1:8]  <- NA  
      ResultsRow[k, 9:16] <- NA
      
      ResultsRow_MM[k, 1:3] <- NA  
      ResultsRow_MM[k, 4:6] <- NA
      
      ResultsRow.ign[k, 1:8]  <- NA  
      ResultsRow.ign[k, 9:16] <- NA
      
      ResultsRow_MM.ign[k, 1:3] <- NA  
      ResultsRow_MM.ign[k, 4:6] <- NA
      
    } else {
      # Evaluate main results ----------------------------------------------------------------------
      # Ignoring non-inv
      Evaluated.ign.con <- evaluation(z_gks    = results$fit.ign.con$posteriors,
                                      beta     = results$fit.ign.con$param$beta_ks,
                                      psi_gks  = results$fit.ign.con$param$psi_gks,
                                      original = original,
                                      nclus    = design[RowDesign, "nclus"],
                                      coeff    = design[RowDesign, "coeff"])
      Evaluated.ign.cat <- evaluation(z_gks    = results$fit.ign.cat$posteriors,
                                      beta     = results$fit.ign.cat$param$beta_ks,
                                      psi_gks  = results$fit.ign.cat$param$psi_gks,
                                      original = original,
                                      nclus    = design[RowDesign, "nclus"],
                                      coeff    = design[RowDesign, "coeff"])
      
      # Store the results
      ResultsRow.ign[k, 1:8]  <- unlist(Evaluated.ign.con)  
      ResultsRow.ign[k, 9:16] <- unlist(Evaluated.ign.cat)
      
      # Including non-inv
      Evaluated.con <- evaluation(z_gks    = results$fit.con$posteriors,
                                  beta     = results$fit.con$param$beta_ks,
                                  psi_gks  = results$fit.con$param$psi_gks,
                                  original = original,
                                  nclus    = design[RowDesign, "nclus"],
                                  coeff    = design[RowDesign, "coeff"])
      Evaluated.cat <- evaluation(z_gks    = results$fit.cat$posteriors,
                                  beta     = results$fit.cat$param$beta_ks,
                                  psi_gks  = results$fit.cat$param$psi_gks,
                                  original = original,
                                  nclus    = design[RowDesign, "nclus"],
                                  coeff    = design[RowDesign, "coeff"])
      
      # Store the results
      ResultsRow[k, 1:8]  <- unlist(Evaluated.con)  
      ResultsRow[k, 9:16] <- unlist(Evaluated.cat)
      
      # Evaluate MM parameters (secondary results) -------------------------------------------------
      # Including non-inv
      Evaluated.con_MM <- evaluation_MM(lambda_est = results$fit.con$param$lambda, 
                                        theta_est  = results$fit.con$param$theta, 
                                        lambda     = results$SimData$lambda, 
                                        theta      = results$SimData$theta, 
                                        ngroups    = design[RowDesign, "ngroups"])
      Evaluated.cat_MM <- evaluation_MM(lambda_est = results$fit.cat$param$lambda, 
                                        theta_est  = results$fit.cat$param$theta, 
                                        lambda     = results$SimData$lambda, 
                                        theta      = results$SimData$theta, 
                                        ngroups    = design[RowDesign, "ngroups"])
      
      # Store the results
      ResultsRow_MM[k, 1:3] <- unlist(Evaluated.con_MM)  
      ResultsRow_MM[k, 4:6] <- unlist(Evaluated.cat_MM)
      
      # Ignoring non-inv
      Evaluated.ign.con_MM <- evaluation_MM(lambda_est = results$fit.ign.con$param$lambda, 
                                            theta_est  = results$fit.ign.con$param$theta, 
                                            lambda     = results$SimData$lambda, 
                                            theta      = results$SimData$theta, 
                                            ngroups    = design[RowDesign, "ngroups"])
      Evaluated.ign.cat_MM <- evaluation_MM(lambda_est = results$fit.ign.cat$param$lambda, 
                                            theta_est  = results$fit.ign.cat$param$theta, 
                                            lambda     = results$SimData$lambda, 
                                            theta      = results$SimData$theta, 
                                            ngroups    = design[RowDesign, "ngroups"])
      
      # Store the results
      ResultsRow_MM.ign[k, 1:3] <- unlist(Evaluated.ign.con_MM)  
      ResultsRow_MM.ign[k, 4:6] <- unlist(Evaluated.ign.cat_MM)
    }
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Normal/Result", "Row", RowDesign,".Rdata" , sep =""))
  save(ResultsRow.ign, file = paste("Ignored/ResultIgn", "Row", RowDesign,".Rdata" , sep =""))
  
  save(ResultsRow_MM, file = paste("Normal/MM_Result", "Row", RowDesign,".Rdata" , sep =""))
  save(ResultsRow_MM.ign, file = paste("Ignored/MM_ResultIgn", "Row", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("C:/Users/perezalo/Documents/GitHub/OrdinalSim/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 1 # Number of replications per condition

Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 16))
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)









