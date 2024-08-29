library(lavaan)

# Set the working directory
setwd("D:/Andres/Functions")

# Source the relevant functions
source("MMG-SEM.R")
source("E_Step.R")
source("ModelSelection.R")

setwd("D:/Andres")
source("DataGeneration.R")
source("evaluation.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus            <- c(2, 4)                 # Number of clusters
ngroups          <- c(24)                   # Number of groups
coeff            <- c(0.3, 0.4)             # Initial regression parameters
N_g              <- c(100, 200)             # Sample size per groups
balance          <- c("bal", "unb")         # Cluster size
NonInvThreshSize <- c(0, 0.1)               # Threshold non-invariance size
c                <- c(2, 3, 4, 5)           # Number of categories

# What about the distribution of the categories?

#threshold <- c("equal", "unequal") # Equality of thresholds across items

# reliability <- c("low")
# NonInvSize <- c(0.6)
# ResRange <- 0.2
# NonInvItems <- 2
# NonInvG <- c(0.50)
# NonInvType <- c("fixed")

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
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, NonInvThreshSize, c, threshold, model) # , reliability, NonInvSize, ResRange,
                      # NonInvItems, NonInvG, NonInvType)
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "NonInvThreshSize", "c", "threshold", "model")
                      # "reliability", "NonInvSize", "ResRange", "NonInvItems", "NonInvG", "NonInvType")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, c, threshold, NonInvThreshSize) #, NonInvG, NonInvItems, NonInvSize, reliability, ResRange)

# Functions for the simulation
# First, to avoid stopping due to errors, create a function with data generation and MMGSEM
# Errors come from non positive definite cov matrices. This code allows the re-sample
genDat_analysis <- function(seed, Design, RowDesign, k, NonInv){
    # browser()
  tryCatch({
    # Set seed per design condition (row) and replication (K)
    set.seed(seed)
    
    # Generate data
    #SimData <- do.call(what = DataGeneration, args = Design[RowDesign, ])$SimData
    SimData <- DataGeneration(model            = Design[RowDesign, "model"], 
                              nclus            = Design[RowDesign, "nclus"], 
                              ngroups          = Design[RowDesign, "ngroups"], 
                              reg_coeff        = Design[RowDesign, "coeff"], 
                              N_g              = Design[RowDesign, "N_g"], 
                              balance          = Design[RowDesign, "balance"], 
                              NonInvThreshSize = Design[RowDesign, "NonInvThreshSize"],
                              c                = Design[RowDesign, "c"],
                              threshold        = Design[RowDesign, "threshold"])

    fit.con <- MMGSEM(dat = SimData$SimData, step1model = S1, step2model = S2, group = "group", 
                      nclus = Design[RowDesign, "nclus"], seed = seed,
                      nstarts = 20, allG = T, est_method = "local", ordered = F, 
                      NonInv = NonInv, constraints = c("loadings", "thresholds"))
    
    fit.cat <- MMGSEM(dat = SimData$SimData, step1model = S1, step2model = S2, group = "group", 
                      nclus = Design[RowDesign, "nclus"], seed = seed,
                      nstarts = 20, allG = T, est_method = "local", ordered = T, 
                      NonInv = NonInv, constraints = c("loadings"))#, "thresholds"))
    
    results <- list(fit.con = fit.con, fit.cat = fit.cat)
    
    # If everything goes right, return results
    return(results)
  }, error = function(e){
    return(NULL)
  })
}

# Main simulation function
do_sim <- function(Design, RowDesign, K){
  # Create the original clustering matrix for comparison below
  original <- create_original(balance = Design[RowDesign, "balance"], 
                              ngroups = Design[RowDesign, "ngroups"], 
                              nclus = Design[RowDesign, "nclus"])
  
  # Create matrix to store results
  # 10 columns for: ARI and RMSEA * 2 (continuous and ordinal) 
  ResultsRow <- matrix(data = NA, nrow = (K), ncol = 16)
  colnames(ResultsRow) <- c("ARI.con", "CC.con", "RMSE_B1.con", "RMSE_B2.con", "RMSE_B3.con", "RMSE_B4.con", "exo_mean.con", "cov_mean.con", 
                            "ARI.cat", "CC.cat", "RMSE_B1.cat", "RMSE_B2.cat", "RMSE_B3.cat", "RMSE_B4.cat", "exo_mean.cat", "cov_mean.cat")
  
  # Define non-invariances
  NonInv <- c("F1 =~ x2", "F1 =~ x3",
              "F2 =~ z2", "F2 =~ z3",
              "F3 =~ m2", "F3 =~ m3",
              "F4 =~ y2", "F4 =~ y3",
              "x2|t2", "z2|t2", "m2|t2", "y2|t2")
  
  for(k in 1:K){
    print(paste("Replication", k, "out of", K))
    # browser()
    # Code to re-sample in case the covariance matrix is non positive definite
    attempts <- 10
    for(j in 1:attempts){
      # Seed will change if there is an error
      ctimes <- system.time(results <- genDat_analysis(seed = (RowDesign * k * j), Design = Design, RowDesign = RowDesign, k = k, NonInv = NonInv))
      if(!is.null(results)){
        # If there was no error, break the loop and continue
        break
      }
    }
    
    # # Save computation times
    # save(ctimes, file = paste("Times/Time", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
    # save(ctimes_ignored, file = paste("Times/TimeIgn", "Row", RowDesign, "Rep", k, ".Rdata", sep = ""))
    # 
    # # Save results if necessary
    # results <- results$Overview
    # save(results, file = paste("Fit/Fit", "Row", RowDesign, "Rep", k, "-", j, ".Rdata" , sep = ""))
    # 
    # results_ignored <- results_ignored$Overview
    # save(results_ignored, file = paste("Fit/FitIgn", "Row", RowDesign, "Rep", k, "-", j, ".Rdata" , sep = ""))
    # browser()
    # ---------------------------------------------------------------
    # Evaluate the results
    Evaluated.con <- evaluation(z_gks    = results$fit.con$posteriors,
                                beta     = results$fit.con$param$beta_ks,
                                psi_gks  = results$fit.con$param$psi_gks,
                                original = original,
                                nclus    = Design[RowDesign, "nclus"],
                                coeff    = Design[RowDesign, "coeff"])
    Evaluated.cat <- evaluation(z_gks    = results$fit.cat$posteriors,
                                beta     = results$fit.cat$param$beta_ks,
                                psi_gks  = results$fit.cat$param$psi_gks,
                                original = original,
                                nclus    = Design[RowDesign, "nclus"],
                                coeff    = Design[RowDesign, "coeff"])
    
    # Store the results
    ResultsRow[k, 1:8] <- unlist(Evaluated.con);   
    ResultsRow[k, 9:16] <- unlist(Evaluated.cat)
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Result", "Row", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("D:/Andres/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 5 # Number of replications per condition

Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 16))
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)

system.time(for(i in 1:30){
  cat("\n", "Condition", i, "out of", nrow(design), "\n")
  Results <- do_sim(Design = design, RowDesign = i, K = K)
  Results_final[(K*(i-1)+1):(i*K), 1:16] <- Results
})

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")
