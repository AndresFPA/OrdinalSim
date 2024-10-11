# __________________________________________________________________________________________________
# Start parallelization --------------------------------------------------------------------------
# __________________________________________________________________________________________________
library(foreach)
library(parallel)
library(doParallel)
library(doSNOW)

rm(results)

# Get object names
obj <- objects()

# Setup parallel cluster
cl <- makeCluster(2)

# registerDoParallel(cl)

# Divide rows per cluster
# rows_divided <- split(1:2, 1:2)

# Export the necessary function and variables to the cluster
clusterEvalQ(cl, {
  library(lavaan)
  library(MASS)
  library(combinat)
  })
clusterEvalQ(cl, setwd("C:/Users/perezalo/Documents/GitHub/OrdinalSim/Results"))
parallel::clusterExport(cl, varlist = obj)

# Parallel execution over RowDesign
# results <- parLapply(cl = cl, X = 1:2, fun = do_sim)
results <- foreach(RowDesign = 1:2, .options.snow = opts) %dopar% {
  do_sim(RowDesign)
}

stopCluster(cl)
