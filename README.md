# MC-Simulation-Ensemble-Methods-Heterogeneous-Treatment-Effects

This repository features the accompanying code to the Master thesis 'Ensemble of Machine Learning Techniques to Detect Heterogeneous Treatment Effects in Selection on Observables Designs'.
It evolves around a Monte Carlo simulation that assesses the relative efficiency of Ensemble Methods. 

The code is split into three major files, 'functions.R', 'MCS.R' (Monte Carlo Simulation) and 'MCV.R' (Monte Carlo Visualization).

The first file 'functions.R' holds custom functions featured in the two latter files and can be ignored by the average user.

The file 'MCS.R', short for Monte Carlo Simulation, features the code for the actual Monte Carlo simulation, generates random variables, assembles them to DGPs and iterates through the estimation part. It returns and saves IATE estimates and informative conceptual data to a created folder called 'output', organized by sample sizes..
It can be modified to gather alternative information such as weights for Super Learning or similar.
Due to the complexity of the simulation, few special arrangements have to be setup before start.
First time users might need to install a Java Virtual Machine, required for certain packages used throughout the simulation.
Additionally, users need to establish a socket connection, implemented to receive intermediate feedback on the simulation's progress while running in parallel, that might require additional software, depending on the OS. Regardless of additional software requirements, each user must establish the connection by running "C:\nc64 -L -p 4000" in the command prompt of a Windows OS or "nc -l 4000" in the terminal under Mac OS and Linux.

The last file, 'MCV.R', short for Monte Carlo visualization, produces the analysis and visual components of the study.
It seamlessly takes on the output produced by 'MCS.R', converts the IATEs to GATEs and analysis the techniques' performances.
It reproduces all outputs, tables and figures identically to the paper and all versions of the simulation and saves them to the created folder 'analysis', organized by sample sizes.

Lastly, this repository features the two folders and data generated throughout the research to the underlying paper.


Abstract of the paper:
This paper investigates the value of Ensemble Methods when estimating heterogeneous treatment effects under confounding in finite samples, via a Monte Carlo simulation. It is the first simulation study that considers Ensemble Methods in Selection on Observables with heterogeneous treatment effects and thus closes the gap in the literature by extending Samii, Paler, and Daly (2016) to heterogeneous treatment effects, as well as extending Grimmer, Messing, and Westwood (2017) to settings with confounding. Ultimately, this research can guide practitioners as to whether Ensemble Method techniques, which proved very promising in forecasting and related fields, should be considered when estimating group average treatment effects (GATEs) in Selection on Observables. The study considers an intermediate aggregation level of causal effects since they are important as feasible action rules for practitioners.
I conduct two Monte Carlo simulations with different sample sizes to bench- mark a Naive Ensemble and a Super Learning Ensemble against their six respective component methods, being Elastic Net, Kernel Regularized Least Squares (KRLS), R-Learner, Random Forest, Bayesian generalized linear model (GLM) and Bayesian Additive Regression Tree (BART). The simulation is comprised of five diverse data generating processes (DGPs), differing in treatment-dimensionality, treatment interaction type, linearity and confound- ing strength. It features a diverse set of 385 covariates, addressing high- dimensional estimation problems.
I found that both Ensemble Methods systematically yield competitive mean squared errors (MSEs), smoothing the performance of their components along the diverse DGPs and different sample sizes. However, the results are less striking than in previous studies with different specifications and the application should be considered on a case by case basis in light of the computational demand. Additionally, this research confirmed certain theoretical findings on Ensemble Methods’ behaviour and contributed new insights to Super Learning with another benchmark study that future scholars can consider when assessing Super Learning within the field of Ensemble Methods.