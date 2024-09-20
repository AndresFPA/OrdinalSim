library(MASS)
library(matrixcalc)

#' DataGeneration
#'
#' Generates simulated data according to the manipulated conditions. 
#'
#' INPUT: Arguments required by the function
#' @param model: full model in lavaan syntax (string).
#' @param nclus: number of clusters of the data. 
#' @param ngroups: number of groups for the generated data.
#' @param N_g: number of observations per group.
#' @param reg_coeff: regression coefficients of cluster 1. Must a numeric vector of length 3.
#' @param balance: can take values "balanced" and "unbalanced". Determines if the clusters will have the same
#'        number of groups.
#' @param reliability: reliability per item (manipulates ratio between loading and unique variance).
#' @param NonInvSize: size of the loading non-invariance (for the simulation: 0.2 and 0.4).
#' @param NonInvItems: number of items affected by the non-invariance (for now, only 2).
#' @param NonInvG: proportion of groups affected by the non-invariance (e.g., 0.25, 0.5).
#' @param randomVarX: define the variance of the exogenous variable. If TRUE, then random values will be defined
#'        for all groups. If FALSE, the variance will be 1 for all groups.
#' @param threshold: Equal or unequal thresholds for all items.
#' @param c: How many categories does the item have.
#'
#' OUTPUT
#' @return SimData: generated data
#' @return NonInvIdx: Index of the groups that presented the non-invariant loadings 
#' @return psi_g: generated psi matrix
#' @return OrTheta: generated theta matrix
#' @return cov_eta: generated cov_eta matrix (phi in the paper)

# Note:
# In the code, we use exog and endog to name the latent variables, while in the paper we use F1, F2, ..., etc.
# The corresponding matches are:
# F1 = exog2
# F2 = exog1
# F3 = endog1
# F4 = endog2

# Note 2:
# To ease the reading of the code, the names of the regression parameters can be seen below:
#   # exog1 -> endog2 = B1
#   # exog1 -> endog1 = B2
#   # endog1 -> endog2 = B3
#   # endog1 -> endog2 = B4


DataGeneration <- function(model, nclus, ngroups, N_g,
                           reg_coeff, balance, c,
                           reliability = "high", NonInvSize = 0.4, 
                           NonInvItems = 2, NonInvG = 0.5, NonInvType = "random",
                           ResRange = 0.2, NonInvGThresh = 0.5, NonInvThreshSize,
                           randomVarX = T){
  
  # Get number of variables
  par_table <- lavaanify(model) # Create parameter table
  lat_var <- lavNames(par_table, "lv")
  obs_var <- lavNames(par_table, "ov")
  m <- length(lat_var) # How many latent variables?
  p <- length(obs_var) # How many observed variables?
  
  # Identify type of latent variable
  endog1 <- lat_var[(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog2 <- lat_var[!c(lat_var %in% par_table$rhs[which(par_table$op == "~")]) &
                      (lat_var %in% par_table$lhs[which(par_table$op == "~")])]
  endog <- c(endog1, endog2)
  exog <- lat_var[!c(lat_var %in% endog)]
  
  # Reorder latent variables
  lat_var <- c(exog, endog)
  
  # Get a cluster label for each group (GperK)
  # E.g., If balanced: G = 6, K = 2. GperK = 111222
  if (balance == "bal"){
    GperK <- rep(x = 1:nclus, each = (ngroups/nclus))
  } else if (balance == "unb"){
    largest <- ngroups*.75; smaller <- ngroups - largest # largest and smaller clusters
    GperK <- c(rep(x = 1, times = largest), rep(x = 2:nclus, each = smaller/(nclus - 1)))
  }
  
  # Give names to the regression parameters (only for better understanding):
  B1 <- numeric(nclus)
  B2 <- numeric(nclus)
  B3 <- numeric(nclus)
  B4 <- numeric(nclus)
  # browser()
  # Get regression parameters for each cluster
  # Cluster 1
  B1[1] <- 0
  B2[1] <- reg_coeff
  B3[1] <- reg_coeff
  B4[1] <- reg_coeff
  
  # Cluster 2
  B1[2] <- reg_coeff
  B2[2] <- 0
  B3[2] <- reg_coeff
  B4[2] <- reg_coeff
  
  if (nclus == 4){
    # Cluster 3
    B1[3] <- reg_coeff
    B2[3] <- reg_coeff
    B3[3] <- 0
    B4[3] <- reg_coeff
    
    # Cluster 4
    B1[4] <- reg_coeff
    B2[4] <- reg_coeff
    B3[4] <- reg_coeff
    B4[4] <- 0
  }
  
  # Generate the structural parameters
  # BETA - Regression parameters
  # beta is cluster-specific. Only nclus matrices needed
  beta <- array(data = 0, dim = c(m, m, nclus), dimnames = list(lat_var, lat_var))
  colnames(beta) <- rownames(beta) <- lat_var
  
  # Fill in correct regression parameters
  # Initialize
  for(k in 1:nclus){
    # beta
    beta[endog2, exog[1], k] <- B1[k]
    beta[endog1, exog[1], k] <- B2[k]
    beta[endog1, exog[2], k] <- B3[k]
    beta[endog2, endog1, k] <- B4[k]
  }
  
  # psi is group- and cluster-specific. ngroups matrices are needed.
  # Generate enough psi matrices
  psi_g <- array(data = diag(m), dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))
  
  # Generate group-specific values for psi sampling from uniform distributions
  if (randomVarX == T){
    exog_var1 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_var2 <- runif(n = ngroups, min = 0.75, max = 1.25)
    exog_cov <- runif(n = ngroups, min = 0.30, max = 0.60)
    endo_var1 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
    endo_var2 <- runif(n = ngroups, min = 0.75, max = 1.25) # Total endogenous variance
  } else { # DEPRECATED
    exog_var1 <- rep(1, times = ngroups) 
    exog_var2 <- rep(1, times = ngroups) 
    exog_cov <- rep(1, times = ngroups) 
  }
  
  # Insert the corresponding group- and cluster-specific parts of psi
  for(g in 1:ngroups){
    # Insert group-specific parts
    psi_g[exog[1], exog[1], g] <- exog_var1[g] 
    psi_g[exog[2], exog[2], g] <- exog_var2[g] 
    psi_g[exog[1], exog[2], g] <- exog_cov[g] 
    psi_g[exog[2], exog[1], g] <- exog_cov[g] 
    
    # Insert the group-and-cluster-specific parts
    # For the endogenous variances, start from the total var (endog_var) and subtract the explained variance by the regression
    psi_g[endog1, endog1, g] <- endo_var1[g] - ((B2[GperK[g]]^2 * exog_var1[g]) + 
                                                  (B3[GperK[g]]^2 * exog_var2[g]) + 
                                                  (2 * B2[GperK[g]] * B3[GperK[g]] * exog_cov[g])) # Group-specific
    
    psi_g[endog2, endog2, g] <- endo_var2[g] - ((B1[GperK[g]]^2 * exog_var1[g]) + 
                                                  (B4[GperK[g]]^2 * endo_var1[g]) + 
                                                  (2 * B1[GperK[g]] * B4[GperK[g]] * ((B2[GperK[g]] * exog_var1[g]) + (B3[GperK[g]] * exog_cov[g])))) # Group-specific
  }
  
  # Create the covariance matrix of the factors (phi in the paper)
  I <- diag(m) # Identity matrix
  cov_eta <- array(data = 0, dim = c(m, m, ngroups), dimnames = list(lat_var, lat_var))

  for(g in 1:ngroups){
    cov_eta[, , g] <- solve(I - beta[, , GperK[g]]) %*% psi_g[, , g] %*% solve(t(I - beta[, , GperK[g]]))
  }
  
  # Generate the measurement model parameters
  # Lambda (depends on reliability lvl)
  if (reliability == "low"){load <- .4} else if (reliability == "high"){load <- .6}
  loadings <- sqrt(load)
  
  # Create Invariant Lambda
  Lambda <- matrix(data = rep(x = c(1, rep(loadings, 4), rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  
  # Create non-invariant Lambda (depending on the type of non-invariance)
  if(NonInvType == "fixed"){
    LambdaNonInv <- matrix(data = rep(x = c(1, 
                                            (loadings - NonInvSize), 
                                            (loadings + NonInvSize), 
                                            rep(loadings, (4 - NonInvItems)), 
                                            rep(0, p)), times = m)[1:(p*m)], nrow = p, ncol = m)
  } else if (NonInvType == "random"){
    # Sample random non-invariances
    NonInvariantLoadings <- sample(x = c(runif(100, min = (loadings - NonInvSize) - .1, max = (loadings - NonInvSize) + .1), 
                                         runif(100, min = (loadings + NonInvSize) - .1, max = (loadings + NonInvSize) + .1)),
                                   size = NonInvItems*4)
    
    # Create a non-invariant lambda matrix
    LambdaNonInv <- matrix(data = c(1, NonInvariantLoadings[1:NonInvItems], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[(NonInvItems + 1):(NonInvItems*2)], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[((NonInvItems*2) + 1):(NonInvItems*3)], rep(loadings, (4 - NonInvItems)), rep(0, p),
                                    1, NonInvariantLoadings[((NonInvItems*3) + 1):(NonInvItems*4)], rep(loadings, (4 - NonInvItems))),
                           nrow = p, ncol = m)
  }
  
  # Theta
  Theta <- array(data = 0, dim = c(p, p, ngroups))
  
  # Sample group-specific (diagonal of) theta values from an uniform distribution
  for (g in 1:ngroups){
    Theta[, , g] <- diag(runif(n = p, min = ((1 - load) - (ResRange/2)), max = ((1 - load) + (ResRange/2))))
  }
  
  # Generate sample covariance matrix (sigma) per groups
  Sigma <- array(data = 0, dim = c(p, p, ngroups))
  # browser()
  # Include non-invariance in the loadings according to NonInvG
  # Which groups are non-invariant? (random sampling)
  NonInvIdx <- sample(x = 1:ngroups, size = NonInvG*ngroups, replace = F)
  for(g in 1:ngroups){
    # Non-invariance - How many groups?
    if (g %in% NonInvIdx){
      Sigma[, , g] <- LambdaNonInv %*% cov_eta[, , g] %*% t(LambdaNonInv) + Theta[, , g]
    } else if (!c(g %in% NonInvIdx)){
      Sigma[, , g] <- Lambda %*% cov_eta[, , g] %*% t(Lambda) + Theta[, , g]
    }
  }

  
  # Data Generation final
  # Get non-invariant thresholds
  NonInvIdxThresh <- sample(x = 1:ngroups, size = NonInvGThresh*ngroups, replace = F)
  
  # Create object to save non-inv thresh
  NonInvthreshList <- vector(mode = "list", length = length(NonInvIdxThresh)) 
  g_counter <- 0
  
  # For now, mu would be 0 as we are only interested in centered variables
  SimData_ord <- c()
  SimData_con <- c()
  
  for(g in 1:ngroups){
    # Generate the data per group
    tmp_con <- mvrnorm(n = N_g, mu = rep(0, p), Sigma = Sigma[, , g], empirical = T) # Save continuous data, just in case
    tmp_ord <- scale(tmp_con) # For ordinal, first get z-scores to have an easier understanding of them
    
    # Get the group-data in ordinal format following the thresholds
    # Get the thresholds for all groups. 
    # Note: (g == 1) as everything is scaled, so the thresholds for the first group will mean the same for the remaining groups. 
    
    # Define the threshold by quantiles in the first item (anchor item).
    # Randomly shift the thresholds of the remaining items by 0.5 to have different thresholds per item
    if(g == 1){
      # Define thresholds
      thresh <- vector(mode = "list", length = 20)
      if (c == 2) {
        thresh[[1]] <- quantile(tmp_ord[, 1], probs = c(0.5))
        for(j in 2:20){
          thresh[[j]] <- thresh[[1]] + sample(x = c(-0.3, 0.3), size = 1)
        }
      } else if (c == 3) {
        thresh[[1]] <- quantile(tmp_ord[, 1], probs = c(0.33, 0.66))
        for(j in 1:20){
          thresh[[j]] <- thresh[[1]] + sample(x = c(-0.3, 0.3), size = 1)
        }
      } else if (c == 4){
        thresh[[1]] <- quantile(tmp_ord[, 1], probs = c(0.25, 0.50, 0.75))
        for(j in 1:20){
          thresh[[j]] <- thresh[[1]] + sample(x = c(-0.3, 0.3), size = 1)
        }
      } else if (c == 5){
        thresh[[1]] <- quantile(tmp_ord[, 1], probs = c(0.2, 0.4, 0.6, 0.8))
        for(j in 1:20){
          thresh[[j]] <- thresh[[1]] + sample(x = c(-0.3, 0.3), size = 1)
        }
      }
    }
    
    
    # Apply thresholds for each groups depending on whether they are in a non-invariant or invariant group
    # browser()
    if(g %in% NonInvIdxThresh){
      g_counter <- g_counter + 1
      # Name the non-inv thresh list
      names(NonInvthreshList)[g_counter] <- as.character(g)
      
      # First, do all invariant items normally
      for(j in c(1:20)[-c(2,3,7,8,12,13,17,18)]){
        tmp_ord[, j] <- as.factor(as.numeric(cut(tmp_ord[, j], breaks = c(-Inf, thresh[[j]], Inf))))
      }
      
      # Repeat it for the non-invariant items
      NonInvthresh <- thresh # done to make sure each group starts at the same point
      # For non-invariant, it depends on the number of categories
      for(j in c(2,3,7,8,12,13,17,18)){
        # Create the non-invariant thresholds by shifting ALL the thresholds randomly to one side 
        if (c == 2){
          NonInvthresh[[j]] <- NonInvthresh[[j]] + sample(x = c(-NonInvThreshSize, NonInvThreshSize), size = 1)
        } else if (c == 4){
          NonInvthresh[[j]] <- NonInvthresh[[j]] + sample(x = c(-NonInvThreshSize, NonInvThreshSize), size = 1)
        } else if (c == 5){
          NonInvthresh[[j]] <- NonInvthresh[[j]] + sample(x = c(-NonInvThreshSize, NonInvThreshSize), size = 1)
        }
        
        # Create the ordinal data using the non-inv thresholds
        tmp_ord[, j] <- as.factor(as.numeric(cut(tmp_ord[, j], breaks = c(-Inf, NonInvthresh[[j]], Inf))))
      }
      
      # Save the non-inv thresh
      NonInvthreshList[[g_counter]] <-  NonInvthresh
      
    } else {
      for(j in 1:20){ # If the group is invariant, simply use the 
        tmp_ord[, j] <- as.factor(as.numeric(cut(tmp_ord[, j], breaks = c(-Inf, thresh[[j]], Inf))))
      }
    }
      
    # browser()
    # Assemble the data with all groups
    SimData_ord <- rbind(SimData_ord, tmp_ord)
    SimData_con <- rbind(SimData_con, tmp_con)
  }
  
  # Add the final labels
  group <- rep(x = c(1:ngroups), each = N_g) # Group variable
  
  SimData_ord <- cbind(SimData_ord, group)
  SimData_con <- cbind(SimData_con, group)
  
  colnames(SimData_ord) <- c(obs_var, "group")
  colnames(SimData_con) <- c(obs_var, "group")
  
  # Return data
  return(list(SimData      = SimData_ord,
              SimData_con  = SimData_con,
              NonInv       = list(load   = NonInvIdx, 
                                  thresh = NonInvIdxThresh), 
              psi_g        = psi_g,
              Sigma        = Sigma, 
              cov_eta      = cov_eta, 
              threshValues = list(thresh       = thresh, 
                                  NonInvThresh = NonInvthreshList)
              )
         )
}
