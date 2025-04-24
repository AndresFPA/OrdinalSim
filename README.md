Code used for the simulation study in the paper "Mixture Multigroup Structural Equation Modeling for Ordinal Data". A quick description of each file is listed below:

# Simulation Functions
**DataGeneration**: Script for the data generation.

**parallel_sim**: main script to run the simulation (parallel computation).

## Analyses & Evaluation Files
The analysis files process the results to bring a summary of the data in the form of a table or figure.

**Analysis**: Main analysis of the simulation. Contains the main results of the structural model and clustering. Compares both versions of MMG-SEM.

**AnalysisMM**: Analysis of the recovery of the measurement model (lambdas). Compares both version of MMG-SEM.

**AnalysisIgn**: Contains the appendix results of the structural model and clustering when *ignoring* the non-invariances. Compares both versions of MMG-SEM.

**AnalysisIgnMM**: Analysis of the recovery of the measurement model (lambdas) when *ignoring* the non-invariances. Compares both version of MMG-SEM.

The evaluation files evaluate the performance of the models per row of the design matrix. That is, they compute the root mean squared error, adjusted rand index, etc., for each condition in the simulation. 

## Functions
It contains the functions required for the simulation, such as the main function of MMG-SEM.

## Ignored & Normal Folders
These folders contain the results after using the evaluation function. In other words, there is one file per condition of the simulation in each folder. The "ignored" folder contains the results when we ignore the non-invariances, and the "normal" folder contains the results when we include the non-invariances.

## Extra Note
The simulation also included an intermediate step to identify some data sets for which MMG-SEM did not converge. Such intermediate "analysis" is not included here to avoid including all the code for re-runs, which would make the repository confusing. The results included here are *after* re-running the 'non-converged' data sets.
