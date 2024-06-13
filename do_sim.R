library(lavaan)

# Set the working directory
setwd("~/GitHub/OrdinalSim/Functions")

# Source the relevant functions
source("MMG-SEM.R")
source("E_Step.R")
source("ModelSelection.R")

setwd("~/GitHub/OrdinalSim/")
source("DataGeneration.R")
source("evaluation.R")
source("evaluationARI.R")

# Simulation Design
# Which factors are going to be tested? For now:
nclus   <- c(2, 4)         # Number of clusters
ngroups <- c(24, 48)       # Number of groups
coeff   <- c(0.3, 0.4)     # Initial regression parameters
N_g     <- c(50, 100, 200) # Sample size per groups
balance <- c("bal", "unb") # Cluster size
sd      <- c(0, 0.05, 0.1) # Differences within a cluster (in betas)

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
design <- expand.grid(nclus, ngroups, coeff, N_g, balance, sd, model) # , reliability, NonInvSize, ResRange,
                      # NonInvItems, NonInvG, NonInvType)
colnames(design) <- c("nclus", "ngroups", "coeff", "N_g", "balance", "sd", "model")
                      # "reliability", "NonInvSize", "ResRange", "NonInvItems", "NonInvG", "NonInvType")

rownames(design) <- NULL
rm(balance, coeff, N_g, nclus, ngroups, sd) #, NonInvG, NonInvItems, NonInvSize, reliability, ResRange)

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
    SimData <- DataGeneration(model     = Design[RowDesign, "model"], 
                              nclus     = Design[RowDesign, "nclus"], 
                              ngroups   = Design[RowDesign, "ngroups"], 
                              reg_coeff = Design[RowDesign, "coeff"], 
                              N_g       = Design[RowDesign, "N_g"], 
                              balance   = Design[RowDesign, "balance"], 
                              sd        = Design[RowDesign, "sd"])

    
    # Run model selection from 1 to 6 clusters
    # 1. BOTH RES AND LOAD NON-INV ARE INCLUDED
    results <- ModelSelection(dat = SimData$SimData, step1model = S1, step2model = S2,
                              group = "group", clusters = c(1, 6), nstarts = 25, seed = (RowDesign * k), 
                              constraints = "loadings", allG = T, fit = "factors", NonInv = NonInv)
    
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
  # 12 columns for: BIC_G, BIC_N, AIC, AIC3, Chull, ICL 
  # There are 2 columns for each model selection measure (for factors and observed LL)
  ResultsRow <- matrix(data = NA, nrow = (K), ncol = 13)
  ResultsRow_ignored <- matrix(data = NA, nrow = (K), ncol = 13)
  
  # Create second matrix for the ARI
  ResultsRowARI <- matrix(data = NA, nrow = (K), ncol = 2)
  ResultsRowARI_ignored <- matrix(data = NA, nrow = (K), ncol = 2)
  
  # Define non-invariances
  NonInv <- c("F1 =~ x2", "F1 =~ x3",
              "F2 =~ z2", "F2 =~ z3",
              "F3 =~ m2", "F3 =~ m3",
              "F4 =~ y2", "F4 =~ y3")
  
  for(k in 1:K){
    print(paste("Replication", k, "out of", K))
    
    # Code to re-sample in case the covariance matrix is non positive definite
    attempts <- 10
    for(j in 1:attempts){
      # Seed will change if there is an error
      ctimes <- system.time(results <- genDat_analysis(seed = (RowDesign * k * j), Design = Design, RowDesign = RowDesign, k = k, NonInv = NonInv)) 
      ctimes_ignored <- system.time(results_ignored <- genDat_analysis(seed = (RowDesign * k * j), Design = Design, RowDesign = RowDesign, k = k, NonInv = NULL)) 
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
    
    # ---------------------------------------------------------------
    # Evaluate the results
    Evaluated            <- evaluation(res = results, clus = Design[RowDesign, "nclus"])
    Evaluated_ignored    <- evaluation(res = results, clus = Design[RowDesign, "nclus"])
    
    EvaluatedARI         <- evaluationARI(z_gks    = results$Models[[Design[RowDesign, "nclus"]]]$posteriors,
                                          original = original,
                                          nclus    = Design[RowDesign, "nclus"])
    EvaluatedARI_ignored <- evaluationARI(z_gks    = results$Models[[Design[RowDesign, "nclus"]]]$posteriors,
                                          original = original,
                                          nclus    = Design[RowDesign, "nclus"])
    
    # Store the results
    colnames(ResultsRow) <- colnames(Evaluated); colnames(ResultsRow_ignored) <- colnames(Evaluated_ignored)
    ResultsRow[k, ]      <- unlist(Evaluated);   ResultsRow_ignored[k, ]      <- unlist(Evaluated_ignored)
    
    colnames(ResultsRowARI) <- c("ARI", "CC");  colnames(ResultsRowARI_ignored) <- c("ARI", "CC")
    ResultsRowARI[k, ] <- unlist(EvaluatedARI); ResultsRowARI_ignored[k, ] <- unlist(EvaluatedARI_ignored)
  }
  
  # Save the results for each row
  save(ResultsRow, file = paste("Result", "Row", RowDesign,".Rdata" , sep =""))
  save(ResultsRow_ignored, file = paste("ResultIgn", "Row", RowDesign,".Rdata" , sep =""))
  
  save(ResultsRowARI, file = paste("Result", "Row", "ARI", RowDesign,".Rdata" , sep =""))
  save(ResultsRowARI_ignored, file = paste("Result_Ign", "Row", "ARI", RowDesign,".Rdata" , sep =""))
  
  # Return the final results
  return(ResultsRow)
}

# Set working directory for the results
# Post-IMPS
setwd("e:/Users/perezalo/Documents/ModelSelection_Simulation/Results")

# Create final results matrix 
# Everything is multiplied by 2 because we run the model twice (including and not including Non-Inv)
K <- 1 # Number of replications per condition

Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 13))
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)

system.time(for(i in 1:1){
  cat("\n", "Condition", i, "out of", nrow(design), "\n")
  Results <- do_sim(Design = design, RowDesign = i, K = K)
  Results_final[(K*(i-1)+1):(i*K), 1:13] <- Results
})

save(Results_final, file = "FinalResults.Rdata")
save(design, file = "design.Rdata")
