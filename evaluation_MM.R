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
  # browser()
  # Return results ---------------------------------------------------------------------------------
  return(list(ParamRecovery = list(Tuck_Lambda = Tuck_Lambda, 
                                   RMSE_Lambda = RMSE_Lambda, 
                                   RMSE_Theta  = RMSE_Theta)))
}









