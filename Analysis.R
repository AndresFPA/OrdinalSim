library(lavaan)
library(fpp3)
library(dplyr)
library(xtable)
library(ggpubr)
library(ggplot2)
library(qwraps2)
library(wesanderson)
# library(Cairo)
# CairoWin()

# Set wd
setwd("C:/Users/User/OneDrive - Tilburg University/1. Papers/Paper 3/R/24-11-06 Sim Results/1. Definitive Results")

# Load empty final results matrix
load("Ignored/FinalResults.Rdata")
load("Ignored/design.Rdata")
load("Ignored/NonConverged_idx.Rdata")

# Organize the data frame
design$Condition <- as.numeric(rownames(design))
Results_final <- merge(x = design, y = Results_final, by = "Condition")

# Re-order the cols
col_order <- c("Condition", "Replication", "nclus", "ngroups", "coeff", "N_g",
               "balance", "NonInvThreshSize", "NonInvLoadSize", "c", 
               "ARI.con", "CC.con", "RMSE_B1.con", "RMSE_B2.con", "RMSE_B3.con", "RMSE_B4.con", "exo_mean.con", "cov_mean.con", 
               "ARI.cat", "CC.cat", "RMSE_B1.cat", "RMSE_B2.cat", "RMSE_B3.cat", "RMSE_B4.cat", "exo_mean.cat", "cov_mean.cat")
Results_final <- Results_final[, col_order]
rm(col_order)

# Fill in the matrix with all results
ncond <- unique(Results_final$Condition) # How many conditions?
K <- length(unique(Results_final$Replication)) # How many replications?
# Indices for loading
start.idx <- which(colnames(Results_final) == "ARI.con")
end.idx   <- ncol(Results_final)
for (i in ncond) {
  test <- NA
  # test <- try(load(paste0("Ignored/ResultIgnRow", i, ".Rdata")))
  test <- try(load(paste0("Normal/ResultRow", i, ".Rdata")))
  if(!c(class(test) == "try-error")){
    # Results_final[(K*(i-1)+1):(i*K), start.idx:end.idx] <- ResultsRow.ign
     Results_final[(K*(i-1)+1):(i*K), start.idx:end.idx] <- ResultsRow
  }
}

# Remove non-run/non-converged rows
# Remove NA results (non-run cases and non-converged)
Results_final <- Results_final[-NA_idx,] # From the first step

# From the second step
Results_final <- Results_final %>% filter(!is.na(RMSE_B1.cat) & !is.na(RMSE_B1.con))

# General results
Results_final %>% dplyr::select(ARI.con:cov_mean.cat) %>% colMeans() %>% round(., 2)
apply(Results_final[,c(3:10)], 2, table)
table(Results_final[, c(6, 10)])
table(Results_final[, c(9, 10)])

# Transform to numeric
# glimpse(Results_final)
# cl <- Results_final %>% select(MisClass:Exo_cov) %>% colnames()
# Results_final[, cl] <- Results_final %>% select(MisClass:Exo_cov) %>% mutate(across(where(is.character), ~ as.numeric(.x)))

# TRIAL GENERAL ------------------------------------------------------------------------------------
# Cluster Recovery
for_plots <- Results_final %>% group_by(N_g, c, coeff) %>% summarise(across(ARI.con:RMSE_B4.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(ARI.con, ARI.cat), names_to = "Measure", 
                               values_to = "ClusRecovery") %>% dplyr::select(Measure, ClusRecovery, N_g, c, coeff)

# for_plots_long <- for_plots_long %>% filter(coeff == 0.2)

# CairoPNG(filename = "C:/Users/perezalo/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)
# CairoPNG(filename = "C:/Users/User/OneDrive - Tilburg University/First paper/ClusSim2.png", width = 650, height = 650)

ggplot(for_plots_long, aes(x = N_g, y = ClusRecovery, colour = Measure)) +
  geom_line(size = 1) +
  geom_point() +
  
  # Update facet grid with wrapped text
  facet_grid(
    c ~ coeff,
    labeller = labeller(
      coeff = function(x) paste("Beta:\n", x),      # Add line break
      c = function(x) paste("Categories (c):\n", x) # Add line break
    )
  ) +
  
  # Axis and legend labels
  labs(
    x = "Within-group Sample Size", 
    y = "Cluster Recovery", 
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
  
  scale_colour_manual(
    values = wes_palette(n = 2, name = "Darjeeling1"), # Keep your desired colours
    labels = c("ARI DWLS", "ARI ML")   # Change labels here
  ) # + 
  
  # Add a clear plot title
  # ggtitle("Effects of Regression Coefficient Size, Sample Size, and Number of Categories on Cluster Recovery")

dev.off()

# Cluster Recovery 2 - NONINV THRESH
for_plots <- Results_final %>% group_by(N_g, NonInvLoadSize, c, NonInvThreshSize) %>% summarise(across(ARI.con:RMSE_B4.cat, mean))
for_plots_long <- pivot_longer(data = for_plots, cols = c(ARI.con, ARI.cat), names_to = "Measure",
                               values_to = "ClusRecovery") %>% dplyr::select(Measure, ClusRecovery, N_g, NonInvLoadSize, c, NonInvThreshSize)

ggplot(for_plots_long, aes(x = N_g, y = ClusRecovery, colour = Measure)) +
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
    y = "Cluster Recovery", 
    colour = "Measure"
  ) +
  
  # General theme settings
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),  # Centre title
    strip.text = yelement_text(size = 10),   # Adjust facet label text size
    legend.text = element_text(size = 10),  # Adjust legend size
    legend.title = element_text(size = 12)
  ) +
  
  # Add a clear plot title
  ggtitle("Effects of Non-Invariance and Sample Size on Cluster Recovery")

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

ggplot(for_plots_long, aes(x = N_g, y = Value, colour = Type)) +
  geom_line(size = 1) +
  geom_point() +
  
  # Update facet grid with wrapped text
  facet_grid(
    c ~ coeff,
    labeller = labeller(
      coeff = function(x) paste("Beta:", x),      # Add line break
      c = function(x) paste("Categories (c):", x) # Add line break
    )
  ) +
  
  # Axis and legend labels
  labs(
    x = "Within-group Sample Size", 
    y = "Parameter Recovery", 
    colour = "Measure"
  ) +
  
  scale_colour_manual(
    values = wes_palette(n = 2, name = "Darjeeling1"), # Keep your desired colours
    labels = c(expression(RMSE~beta~DWLS), expression(RMSE~beta~ML))   # Change labels here
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
  ) # +
  
  # Add a clear plot title
  # ggtitle("Effects of Beta Coefficient Size, Sample Size, and Number of Categories on Parameter Recovery")

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

ggplot(for_plots_long, aes(x = N_g, y = Value, colour = Type)) +
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
    y = "Parameter Recovery", 
    colour = "Type"
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
  ggtitle("Effects of Non-Invariance and Sample Size on Parameter Recovery")

dev.off()

####################################################################################################
############################ TABLES - CLUSTER AND PARAMETER RECOVERY ###############################
####################################################################################################

# Add total RMSE
Results_final$RMSE.con <- Results_final %>% dplyr::select(RMSE_B1.con:RMSE_B4.con) %>% rowMeans()
Results_final$RMSE.cat <- Results_final %>% dplyr::select(RMSE_B1.cat:RMSE_B4.cat) %>% rowMeans()

# Check mean results per simulation factor
a <- Results_final %>% group_by(nclus)            %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
b <- Results_final %>% group_by(N_g)              %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
c <- Results_final %>% group_by(coeff)            %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
d <- Results_final %>% group_by(balance)          %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
e <- Results_final %>% group_by(c)                %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
f <- Results_final %>% group_by(NonInvLoadSize)   %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))
g <- Results_final %>% group_by(NonInvThreshSize) %>% summarise(across(ARI.con:RMSE.cat, qwraps2::mean_sd, denote_sd = "paren", digits = 3))

list2 <- list(a, b, c, d, e, f, g)
current <- c()

for(i in 1:length(list2)){
  tmp <- list2[[i]]
  tmp$Factor <- colnames(tmp)[1]
  colnames(tmp)[1] <- c("Level")
  tmp$Level <- as.factor(tmp$Level)
  tmp <- tmp[, c(1, ncol(tmp), c(2:(ncol(tmp) - 1)))]
  current <- rbind(current, tmp)
}

rm(a, b, c, d, e, f, g)
current <- current[, c(2, 1, 3:20)]
for_paper <- current %>% dplyr::select(-c(exo_mean.con, cov_mean.con, exo_mean.cat, cov_mean.cat))
for_paper <- for_paper %>% dplyr::select(Factor:RMSE_B4.con, RMSE.con, ARI.cat:RMSE_B4.cat, RMSE.cat)

# colnames(for_paper) <- c("Factor", "Level", "ARI.con", "PR.con", "RMSE_B1.con", "RMSE_B2.con", "RMSE_B3.con", "RMSE_B4.con", "RMSE.con",
#                          "ARI.cat", "PR.cat", "RMSE_B1.cat", "RMSE_B2.cat", "RMSE_B3.cat", "RMSE_B4.cat")

# Re-organize table
# Get row indices of each factor
factors <- unique(for_paper$Factor)
for(i in 1:length(factors)){ 
  assign(x = paste0("rn_", factors[i]), value = which(for_paper$Factor == factors[i]))
}

for_paper <- for_paper[c(rn_coeff, rn_N_g, rn_nclus, rn_balance, rn_NonInvLoadSize,
                         rn_NonInvThreshSize, rn_c), ]

rm(rn_coeff, rn_N_g, rn_nclus, rn_balance, rn_NonInvLoadSize, rn_NonInvThreshSize, rn_c)

# Add totals
Totals <- Results_final %>% dplyr::select(ARI.con:RMSE.cat) %>% apply(., 2, mean_sd, denote_sd = "paren", digits = 3) %>% as.data.frame() %>% t() %>% as.data.frame()
Totals$Factor <- "Total"; Totals$Level <- ""
Totals <- Totals %>% select(Factor, Level, ARI.con:RMSE_B4.con, RMSE.con, ARI.cat:RMSE_B4.cat, RMSE.cat)
for_paper <- rbind(for_paper, Totals)

#Final - Cluster
for_paper_clus <- for_paper %>% dplyr::select(Factor:CC.con, ARI.cat:CC.cat)
print(xtable(for_paper_clus, digits = 3), include.rownames = F)

# Final - Parameter
for_paper_par <- for_paper %>% dplyr::select(Factor:Level, RMSE.con, RMSE.cat)
print(xtable(for_paper_par, digits = 3), include.rownames = F)

Results_final %>% filter(NonInvG == 0.5) %>% group_by(NonInvIncl) %>%  
  select(RMSE_B1:RMSE_B4) %>% summarise(across(RMSE_B1:RMSE_B4, mean_sd, denote_sd = "paren", digits = 3))


