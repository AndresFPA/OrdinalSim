#' Evaluation
#'
#' Evaluation function to assess the performance of the model according to the simulation results. 
#'
#' INPUT: Arguments required by the function
#' @param beta 
#' @param z_gks
#' @param original
#' @param global_max
#' @param runs_LL
#' @param nclus
#' @param reg_coef
#' @param reg_diff

library("combinat")

# evaluation_MM <- function(lambda_est, theta_est, lambda, theta, ngroups){
#   #browser()
#   # Evaluation - Measurement Parameters recovery --------------------------------------------------
#   # Original
#   lambdaVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambda[[x]])[c(lambda[[x]]) != 0]})
#   lambdaVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambdaVec[[x]])[c(lambdaVec[[x]]) != 1]}) # Remove fixed loading
#   lambdaVec <- unlist(lambdaVec)
#   
#   thetaVec <- lapply(X = 1:ngroups, FUN = function(x){c(theta[,,x])[c(theta[,,x]) != 0]})
#   thetaVec <- unlist(thetaVec)
#   
#   # Estimated
#   lambdaEstVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambda_est[[x]])[c(lambda_est[[x]]) != 0]})
#   lambdaEstVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambdaEstVec[[x]])[round(c(lambdaEstVec[[x]]), 2) != 1]})
#   lambdaEstVec <- unlist(lambdaEstVec)
#   
#   thetaEstVec <- lapply(X = 1:ngroups, FUN = function(x){c(theta_est[[x]])[c(theta_est[[x]]) != 0]})
#   thetaEstVec <- unlist(thetaEstVec)
#   
#   # Root Mean Squared Error (RMSE)
#   # SquaredErrorLambda <- numeric(length(lambdaEstVec[[1]])*ngroups)
#   # SquaredErrorTheta <- numeric(length(thetaEstVec[[1]])*ngroups)
#   
#   # Calculate corresponding measures
#   # Cluster 1
#   SquaredErrorLambda <- (lambdaVec - lambdaEstVec)^2
#   SquaredErrorTheta <- (thetaVec - thetaEstVec)^2
#   
#   # RMSE
#   RMSE_Lambda <- sqrt(mean(SquaredErrorLambda))
#   RMSE_Theta <- sqrt(mean(SquaredErrorTheta))
#   
#   # Return results ---------------------------------------------------------------------------------
#   return(list(ParamRecovery = list(RMSE_Lambda = RMSE_Lambda, RMSE_Theta = RMSE_Theta)))
# }

evaluation_MM <- function(lambda_est, theta_est, lambda, theta, ngroups){
  #browser()
  # Evaluation - Measurement Parameters recovery --------------------------------------------------
  # lambda - tuck congruence
  Tuck_Lambda_Group <- c()
  for(g in 1:ngroups){
    Tuck_Lambda_Group[g] <- mean(diag(psych::factor.congruence(x = lambda_est[[g]], y = lambda[[g]])))
  }
  
  # Re-scale lambda (will do nothing to continuous results)
  # Done due to the scaling used when using the categorical estimator (i.e., std.lv = T)
  for (g in 1:ngroups) {
    loadings <- apply(lambda_est[[g]], 2, \(x) {x[which(x != 0)]})[1, ]
    ratios <- 1/loadings # Ratio to rescale the remaining loadings
    lambda_est[[g]] <- sapply(1:length(ratios), FUN = \(x) {lambda_est[[g]][, x] * ratios[x]})
    colnames(lambda_est[[g]]) <- names(ratios)
  }
  
  # lambda - RMSE
  lambdaVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambda[[x]])[lambda[[x]] != 0 & round(lambda[[x]], 4) != 1]})
  lambdaVec <- unlist(lambdaVec)

  lambdaEstVec <- lapply(X = 1:ngroups, FUN = function(x){c(lambda_est[[x]])[lambda_est[[x]] != 0 & round(lambda[[x]], 4) != 1]})
  lambdaEstVec <- unlist(lambdaEstVec)
  
  # theta - RMSE
  thetaVec <- lapply(X = 1:ngroups, FUN = function(x){c(theta[,,x])[c(theta[,,x]) != 0]})
  thetaVec <- unlist(thetaVec)
  
  thetaEstVec <- lapply(X = 1:ngroups, FUN = function(x){c(theta_est[[x]])[c(theta_est[[x]]) != 0]})
  thetaEstVec <- unlist(thetaEstVec)
  
  # Root Mean Squared Error (RMSE)
  # Calculate corresponding measures
  SquaredErrorLambda <- (lambdaVec - lambdaEstVec)^2
  SquaredErrorTheta  <- (thetaVec - thetaEstVec)^2
  
  # RMSE
  # RMSE_Lambda <- sqrt(mean(SquaredErrorLambda))
  Tuck_Lambda <- mean(Tuck_Lambda_Group)
  RMSE_Lambda <- sqrt(mean(SquaredErrorLambda))
  RMSE_Theta  <- sqrt(mean(SquaredErrorTheta))
  
  # Relative bias
  RelativeBiasLambda <- (lambdaEstVec - lambdaVec)/lambdaVec
  RelativeBiasTheta <- (thetaVec - thetaEstVec)/thetaVec
  
  # Mean rel bias
  Mean_RelBiasLambda <- mean(RelativeBiasLambda)
  Mean_RelBiasTheta  <- mean(RelativeBiasTheta)
  
  # browser()
  # Return results ---------------------------------------------------------------------------------
  return(list(ParamRecovery = list(Tuck_Lambda    = Tuck_Lambda, 
                                   RMSE_Lambda    = RMSE_Lambda, 
                                   RMSE_Theta     = RMSE_Theta,
                                   RelBias_Lambda = Mean_RelBiasLambda,
                                   RelBias_Theta  = Mean_RelBiasTheta)))
}









