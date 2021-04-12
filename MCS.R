rm(list = ls(all.names = TRUE))
set.seed(0815)

#Set working device
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import packages
pacman::p_load()

#Import custom functions
source('functions.R')

#Generate matrix of covariates X as outlined in Appendix B.


#Calculate treatment probability per unit based on joint distribution of con- founding variables.


#Draw treatment status Di per unit.


#Assign heterogeneous treatment effects per unit and DGP as outlined in Appendix C.


#Calculate outcomes Y per DGP.


#Randomly draw N pairs of outcomes and covariates plus treatment status.


#Randomly split sample of size N into folds of size ND , whereby I use D = 1 in-line with Grimmer et al. (2017).


#Generate out of sample predictions for all M component methods using all remaining D − 1-folds for training. 􏰎􏰎


#Obtain N × M matrix for Y per DGP, whereby Yi,m stands for unit i’s out of sample prediction from method m.


#Regress true response on out of sample predictions (Yi = 􏰋M wmYi,m + εi) i=1 per DGP with constraints 􏰋Mi=1 wm = 1 and wm ≥ 0, whereby ε is an error term.


#Obtain solutions to constrained regression w􏰎 as weights for final Ensemble per DGP.


#Obtain estimates for heterogeneous treatment effects per DGP using the full sample by all nine component methods, following the procedures outlined throughout chapter 2.3 and OLS as a benchmark.


#Calculate heterogeneous treatment effect estimates per DGP by a naive Ensemble and by Super Learning with weights from step 11.


#Iterate through steps 3 to 14 with different sample sizes N.


#Obtain performance and convergence measures for the Ensemble Methods, all components and OLS per DGP.