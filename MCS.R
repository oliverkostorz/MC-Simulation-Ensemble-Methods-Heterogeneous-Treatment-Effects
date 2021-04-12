rm(list = ls(all.names = TRUE))
set.seed(0815)

#Set working device
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import packages
pacman::p_load(purrr,extraDistr,poisbinom,actuar,circular,evd,rdetools)

#Import custom functions
source('functions.R')

#######################################################
############## Set simulation parameters ##############
#######################################################
#Simulation size
N <- 1000


#######################################################

#######################################################
##################### Covariates ######################
#######################################################
#Bernoulli
bernoulli <- as.numeric(rbernoulli(N, p = 0.5))

#Rademacher
rademacher <- rsign(N)

#Poisson binomial
poissonbinomial <- rpoisbinom(N, 0.4)

#Binomial
binomial <- rbinom(N, 6, 0.5) - 2

#Hypergeometric
hypergeometric <- rhyper(N, N/3, 2*N/3, N/2)

#Geometric distribution
geometric <- rgeom(N, 0.3)

#Logarithmic distribution
logarithmic <- rlogarithmic(N, 0.9)

#Poisson distribution
poisson <- rpois(N, 3)

#Beta distribution
beta <- rbeta(N, 2, 4, ncp = 0.5)

#Wrapped Cauchy distribution
wrappedcauchy <- rwrappedcauchy(N, mu = circular(4), rho = exp(-100))

#Wrapped Normal distribution
wrappednormal <- rwrappednormal(N, mu = circular(-4), sd = 9)

#Chi-squared distribution
chisquared <- rchisq(N, 4)

#Exponential distribution
exponential <- rexp(N, rate = 10)

#Gamma distribution
gamma <- rgamma(N, 0.8, scale = 1/1.5)

#Log-normal distribution 
lognormal <- rlnorm(N, meanlog = -2, sdlog = 4)

#Weibull distribution
weibull <- rweibull(N, 0.9, scale = 1)

#Gumbel distribution
gumbel <- rgumbel(N, loc = -5, scale = 10)

#Laplace distribution
laplace <- rlaplace(N, mu = 100, sigma = 5)

#Logistic distribution
logistic <- rlogis(N, location = 0, scale = 9)

#Normal distribution
normal <- rnorm(N, mean = -5, sd = 5)

#Student's t-distribution
studentst <- rt(N, 9)

#Bind randomly drawn covariates to dataframe
rvar <- cbind(bernoulli, rademacher, poissonbinomial, binomial,
              hypergeometric, geometric, logarithmic, poisson,
              beta, wrappedcauchy, wrappednormal, chisquared,
              exponential, gamma, lognormal, weibull, gumbel,
              laplace, logistic, normal, studentst)

rvar_nonbinary <- cbind(binomial, hypergeometric, geometric,
                        logarithmic, poisson, beta, wrappedcauchy,
                        wrappednormal, chisquared, exponential,
                        gamma, lognormal, weibull, gumbel,
                        laplace, logistic, normal, studentst)

#Sinc transformation
sincvar <- custom_sinc(rvar_nonbinary)

#Exponential transformation
expvar <- custom_exp(rvar_nonbinary)

#Bind transformed covariates to dataframe
transvar <- cbind(sincvar, expvar)

#Bind randomly drawn and transformed covariates to dataframe
vars <- cbind(rvar, sincvar, expvar)

#Crossproducts
crossvar <- custom_crossprod(vars)

#Bind all covariates to dataframe
X <- cbind(vars, crossvar)

#######################################################

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