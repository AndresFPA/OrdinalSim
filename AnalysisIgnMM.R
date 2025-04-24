library(lavaan)
library(fpp3)
library(dplyr)
library(xtable)
library(ggpubr)
library(ggplot2)
library(qwraps2)
library(Cairo)

# CairoWin()

# Set wd
setwd("C:/Users/User/OneDrive - Tilburg University/1. Papers/Paper 3/R/24-11-06 Sim Results/Fit/CorrectedMM")
setwd("C:/Users/perezalo/OneDrive - Tilburg University/1. Papers/Paper 3/R/24-11-06 Sim Results/Fit/CorrectedMM")

# Load design matrix
load("design.Rdata")

# Create empty Final Results matrix
K <- 50 # Number of replications per condition
Results_final <- as.data.frame(matrix(data = NA, nrow = nrow(design)*K, ncol = 10))
Results_final$Replication <- rep(x = 1:K, times = nrow(design))
Results_final$Condition <- rep(x = 1:nrow(design), each = K)

# Organize the data frame
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")

# Add names to the unnamed cols
colnames(Results_final)[which(colnames(Results_final) == "V1"):which(colnames(Results_final) == "V10")] <- c("Tuck_lambda.con", "RMSE_lambda.con", "RMSE_theta.con", "RelBias_Lambda.con", "RelBias_Theta.con",
                                                                                                             "Tuck_lambda.cat", "RMSE_lambda.cat", "RMSE_theta.cat", "RelBias_Lambda.cat", "RelBias_Theta.cat")
# Re-order the cols
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g",
               "balance", "NonInvThreshSize", "NonInvLoadSize", "c", 
               "Tuck_lambda.con", "RMSE_lambda.con", "RMSE_theta.con", "RelBias_Lambda.con", "RelBias_Theta.con",
               "Tuck_lambda.cat", "RMSE_lambda.cat", "RMSE_theta.cat", "RelBias_Lambda.cat", "RelBias_Theta.cat")
Results_final <- Results_final[, col_order]
rm(col_order)

# Fill in the matrix with all results
ncond <- unique(Results_final$Condition) # How many conditions?
K <- length(unique(Results_final$Replication)) # How many replications?
# Indices for loading
start.idx <- which(colnames(Results_final) == "Tuck_lambda.con")
end.idx   <- ncol(Results_final)
for (i in ncond) {
  test <- NA
  # test <- try(load(paste0("Ignored/ResultIgnRow", i, ".Rdata")))
   test <- try(load(paste0("MM_ResultIgnRow", i, ".Rdata")))
  if(!c(class(test) == "try-error")){
   Results_final[(K*(i-1)+1):(i*K), start.idx:end.idx] <- ResultsRow_MM.ign
  }
}

# Remove NA results
Results_final <- Results_final %>% filter(!is.na(Tuck_lambda.con))

# TEMPORARY
# Identify the results where both the loading recovery AND the residual recovery looks impossible
# It mostly happened for categorical estimation
# How to define impossible is difficults. We set the threshold at:
# theta_rmse > 0.6
# loading_rmse > 2
# Remember that true theta and lambda are < 1. Such large rmse values imply very far away estimations.
imp.idx  <- which(Results_final$RMSE_lambda.cat > 2 )#| Results_final$RMSE_theta.cat > 0.6)
imp.idx2 <- which(Results_final$RMSE_lambda.con > 2 )#| Results_final$RMSE_theta.con > 0.6)
Results_final <- Results_final[-c(imp.idx, imp.idx2), ]

# General results
Results_final %>% dplyr::select(Tuck_lambda.con:RelBias_Theta.cat) %>% colMeans(., na.rm = T)
apply(Results_final[,c(3:10)], 2, table)
table(Results_final[, c(6, 10)])
table(Results_final[, c(9, 10)])

# Transform to numeric
# glimpse(Results_final)
# cl <- Results_final %>% select(MisClass:Exo_cov) %>% colnames()
# Results_final[, cl] <- Results_final %>% select(MisClass:Exo_cov) %>% mutate(across(where(is.character), ~ as.numeric(.x)))

# TRIAL GENERAL ------------------------------------------------------------------------------------
# Lambda recovery
for_plots <- Results_final %>% group_by(N_g, c, NonInvLoadSize) %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE_lambda.con, RMSE_lambda.cat), names_to = "Measure", 
                               values_to = "LambdaRecovery") %>% dplyr::select(Measure, LambdaRecovery, N_g, c, NonInvLoadSize)

# for_plots_long <- for_plots_long %>% filter(N_g != 50)

# CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)
# CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = LambdaRecovery, colour = Measure)) + 
  geom_line(size = 1) + 
  facet_grid(c ~ NonInvLoadSize) + labs(x = "Within-group Sample Size", y = "Lambda Recovery") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Loading Non-invariance size", breaks = NULL, labels = NULL)) + 
  # scale_x_discrete(name = "Loading Non-invariance size", position = "top") + 
  theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Number of categories", breaks = NULL, labels = NULL)) +
  #scale_color_hue(labels = c("Continuous", "Categorical")) + #scale_linetype_manual("Non-invariance threshold", values = c(0, 0.2, 0.4)) + #, labels = c("Equal", "Unequal")) +
  # scale_y_discrete(name = "Number of categories", position = "right") +
  ggtitle(label = "Lambda Recovery") + geom_point()

dev.off()

# Lambda recovery 2 - NONINV THRESH
for_plots <- Results_final %>% group_by(N_g, NonInvLoadSize, c, NonInvThreshSize) %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(RMSE_lambda.con, RMSE_lambda.cat), names_to = "Measure",
                               values_to = "LambdaRecovery") %>% dplyr::select(Measure, LambdaRecovery, N_g, NonInvLoadSize, c, NonInvThreshSize)

ggplot(for_plots_long, aes(x = N_g, y = LambdaRecovery, colour = Measure)) +
  geom_line(size = 1) +
  geom_point() +
  
  # Update facet grid with wrapped text
  facet_grid(
    NonInvThreshSize ~ NonInvLoadSize + c,
    labeller = labeller(
      NonInvThreshSize = function(x) paste("Non-Inv \n Threshold Size:\n", x),  # Add line break
      NonInvLoadSize = function(x) paste("Non-Inv \n Loading Size:\n", x),      # Add line break
      c = function(x) paste("Categories (c):\n", x)                   # Add line break
    )
  ) +
  
  # Axis and legend labels
  labs(
    x = "Within-group Sample Size", 
    y = "Lambda Recovery", 
    colour = "Measure"
  ) +
  
  # General theme settings
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),  # Centre title
    strip.text = element_text(size = 10),   # Adjust facet label text size
    legend.text = element_text(size = 10),  # Adjust legend size
    legend.title = element_text(size = 12)
  ) +
  
  # Add a clear plot title
  ggtitle("Effects of Non-Invariance and Sample Size on Lambda Recovery")

####################################################################################################
#################################### PLOTS - PARAMETER RECOVERY ####################################
####################################################################################################

# Parameter Recovery (Only one RMSE)
Results_final$RMSE.con <- Results_final %>% select(RMSE_B1.con:RMSE_B4.con) %>% rowMeans()
Results_final$RMSE.cat <- Results_final %>% select(RMSE_B1.cat:RMSE_B4.cat) %>% rowMeans()

for_plots <- Results_final %>% group_by(N_g, coeff, c) %>% summarise(across(ARI.con:RMSE.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = starts_with("RMSE."), # Specify the columns to pivot
                               names_to = c("Parameter", "Type"), # Use .value to split into multiple columns
                               names_pattern = "(RMSE)\\.(con|cat)",
                               values_to = "Value") %>% dplyr::select(Parameter, Type, Value, N_g, coeff, c)

# for_plots <- Results_final %>% group_by(N_g, coeff, c) %>% summarise(across(ARI.con:RMSE_B4.cat, mean))
# for_plots_long <- pivot_longer(data = for_plots, cols = starts_with("RMSE_"), # Specify the columns to pivot
#                                names_to = c("Parameter", "Type"), # Use .value to split into multiple columns
#                                names_pattern = "(RMSE_B[1-4])\\.(con|cat)",
#                                values_to = "Value") %>% dplyr::select(Parameter, Type, Value, N_g, coeff, c)
#for_plots_long <- for_plots_long %>% filter(coeff == 0.2)

# CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)
# CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = Value, linetype = Type, color = Parameter)) + 
  geom_line(size = 1) + 
  facet_grid(coeff ~ c) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
  # scale_x_continuous(sec.axis = sec_axis(~ . , name = "Model Specification", breaks = NULL, labels = NULL)) + 
  theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  # scale_linetype_manual("Non-inv Type", values = c(1, 2), labels = c("Fixed", "Random")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Regression Parameters", breaks = NULL, labels = NULL)) +
  ggtitle(label = "Parameter Recovery") + geom_point()

dev.off()

# Non-invariances plot
# Parameter Recovery (Only one RMSE)
Results_final$RMSE.con <- Results_final %>% select(RMSE_B1.con:RMSE_B4.con) %>% rowMeans()
Results_final$RMSE.cat <- Results_final %>% select(RMSE_B1.cat:RMSE_B4.cat) %>% rowMeans()

for_plots <- Results_final %>% group_by(N_g, NonInvLoadSize, c, NonInvThreshSize) %>% summarise(across(ARI.con:RMSE.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = starts_with("RMSE."), # Specify the columns to pivot
                               names_to = c("Parameter", "Type"), # Use .value to split into multiple columns
                               names_pattern = "(RMSE)\\.(con|cat)",
                               values_to = "Value") %>% dplyr::select(Parameter, Type, Value, N_g, NonInvLoadSize, c, NonInvThreshSize)

# for_plots <- Results_final %>% group_by(N_g, coeff, c) %>% summarise(across(ARI.con:RMSE_B4.cat, mean))
# for_plots_long <- pivot_longer(data = for_plots, cols = starts_with("RMSE_"), # Specify the columns to pivot
#                                names_to = c("Parameter", "Type"), # Use .value to split into multiple columns
#                                names_pattern = "(RMSE_B[1-4])\\.(con|cat)",
#                                values_to = "Value") %>% dplyr::select(Parameter, Type, Value, N_g, coeff, c)
#for_plots_long <- for_plots_long %>% filter(coeff == 0.2)

# CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)
# CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ParSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = Value, linetype = Type, color = Parameter)) + 
  geom_line(size = 1) + 
  # facet_wrap(~ c) +
  facet_grid(NonInvThreshSize ~ NonInvLoadSize + c) + labs(x = "Within-group Sample Size", y = "Root Mean Squared Error (RMSE)") + 
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "Non-Inv loadings + Number of categories", breaks = NULL, labels = NULL)) + 
  theme_bw() + theme(text=element_text(size=12), legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
  # scale_linetype_manual("Non-inv Type", values = c(1, 2), labels = c("Fixed", "Random")) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "Non-Inv Thresholds", breaks = NULL, labels = NULL)) +
  ggtitle(label = "Parameter Recovery") + geom_point()

dev.off()

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################

# Check mean results per simulation factor
a <- Results_final %>% group_by(nclus)            %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
b <- Results_final %>% group_by(ngroups)          %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
c <- Results_final %>% group_by(N_g)              %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
d <- Results_final %>% group_by(coeff)            %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
e <- Results_final %>% group_by(balance)          %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
f <- Results_final %>% group_by(c)                %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
g <- Results_final %>% group_by(NonInvLoadSize)   %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
h <- Results_final %>% group_by(NonInvThreshSize) %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
# i <- Results_final %>% group_by(ResRange)         %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
# j <- Results_final %>% group_by(NonInvG)          %>% summarise(across(Tuck_lambda.con:RelBias_Theta.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))

# list2 <- list(a, b, c, d, e, f, g, h, i, j)
list2 <- list(a, c, d, e, f, g, h)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[1]
  colnames(tmp)[1] <- c("Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(ncol(tmp), 1, c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g, h, i, j)
# current <- current[, c("Factor", "Level", "Non-inv Included", "ARI", "CorrectClus", "fARI", 
#                        "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "RMSE")]

for_paper <- current
# colnames(for_paper) <- c("Factor", "Level", "Non-inv Included", "ARI", "CorrectClus", "fARI", "RMSE_B1", "RMSE_B2", "RMSE_B3", "RMSE_B4", "RMSE")

# Re-organize table
# Get row indices of each factor
# factors <- unique(for_paper$Factor)
# for(i in 1:length(factors)){ 
#   assign(x = paste0("rn_", factors[i]), value = which(for_paper$Factor == factors[i]))
# }
# 
# for_paper <- for_paper[c(rn_coeff, rn_N_g, rn_nclus, rn_balance, rn_NonInvG,
#                          rn_NonInvSize, rn_NonInvType, rn_ResRange), ]
# 
# rm(rn_coeff, rn_ngroups, rn_N_g, rn_nclus, rn_balance, rn_reliability, rn_NonInvG, rn_NonInvSize, rn_NonInvType, rn_ResRange)

# Add total - Cluster
yes_tot <- Results_final %>% select(Tuck_lambda.con:RelBias_Theta.cat) %>%
  summarise(across(Tuck_lambda.con:RelBias_Theta.cat, mean_sd, denote_sd = "paren", digits = 3))

yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""

#Final - Cluster
for_paper_all <- rbind(for_paper, yes_tot)

# Final organization
for_paper_clus <- for_paper_clus[, c("Factor", "Level", 
                                "ARI_both", "fARI_both",
                                "ARI_none", "fARI_none",
                                "ARI_bothK", "fARI_bothK",
                                "ARI_noneK", "fARI_noneK")]

print(xtable(for_paper_clus, digits = 3), include.rownames = F)

rm(yes_tot, tmp)

# Add total - Parameter
yes_tot <- Results_final %>% filter(NonInvIncl %in% c("both", "bothK", "none", "noneK")) %>% group_by(NonInvIncl) %>%  
                             select(RMSE) %>% summarise(across(RMSE, mean_sd, denote_sd = "paren", digits = 3)) %>% 
                             pivot_wider(names_from = "NonInvIncl", values_from = c("RMSE"))
yes_tot <- as.data.frame(yes_tot)
yes_tot$Factor <- "Total"; yes_tot$Level <- ""
yes_tot <- yes_tot %>% select(Factor, Level, both:noneK)
colnames(yes_tot) <- c("Factor", "Level", "RMSE_both", "RMSE_bothK", "RMSE_none", "RMSE_noneK")

#Final - Parameter
for_paper_par <- for_paper[, c("Factor", "Level", 
                               "RMSE_both", "RMSE_bothK", "RMSE_none", "RMSE_noneK"
                               )]
for_paper_par <- rbind(for_paper_par, yes_tot)

# Final organization
for_paper_par <- for_paper_par[, c("Factor", "Level", 
                                   "RMSE_both", "RMSE_none", "RMSE_bothK", "RMSE_noneK")]

# for_paper_par <- data.frame(lapply(for_paper_par, gsub, pattern = " ", replacement = "")) # Remove unnecessary spaces
print(xtable(for_paper_par, digits = 3), include.rownames = F)

# Extra - Average parameter recovery per regression parameter across all models
Results_final %>% ungroup() %>% select(RMSE_B1:RMSE_B4) %>% colMeans()


rm(yes_tot) 

Results_final %>% filter(NonInvG == 0.5) %>% group_by(NonInvIncl) %>%  
  select(RMSE_B1:RMSE_B4) %>% summarise(across(RMSE_B1:RMSE_B4, mean_sd, denote_sd = "paren", digits = 3))


