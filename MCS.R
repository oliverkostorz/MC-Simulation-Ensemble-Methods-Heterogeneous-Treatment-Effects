rm(list = ls(all.names = TRUE))
set.seed(0815)

#Set working device
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Set path to Java
#Sys.setenv('JAVA_HOME' = '/Library/Java/JavaVirtualMachines/jdk1.8.0_291.jdk/') 
#dyn.load('/Library/Java/JavaVirtualMachines/1.8.0_291.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

#Import packages
pacman::p_load(purrr,extraDistr,poisbinom,actuar,circular,evd,rdetools,
               sets,glmnet,KRLS,devtools,stringr,randomForest)#,rJava,RWeka,SVMMatch,FindIt,mboost,GAMBoost)
install_github('xnie/rlearner')

#Import custom functions
source('functions.R')

#######################################################
############## Set simulation parameters ##############
#######################################################
#Simulation size
N <- 1000
sample_size <- 100 #Only choose sample sizes which are multiples of 10 or amend code for sampling folds


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

#Replace Inf from sinc function with 0
X <- data.frame(apply(X, MARGIN = 2, function(c) unlist(lapply(c, function(z) ifelse(is.finite(z), z, 0)))))

#Parameters for later stages in code
p <- ncol(X)
bXt <- which(colnames(X) == 'binomialXstudentst')
#Don't forget to add pseudo covariates later on
base_variables_name <- c('bernoulli', 'rademacher', 'poissonbinomial', 'binomial',
                         'hypergeometric', 'geometric', 'logarithmic', 'poisson',
                         'beta', 'wrappedcauchy', 'wrappednormal', 'chisquared',
                         'exponential', 'gamma', 'lognormal', 'weibull', 'gumbel',
                         'laplace', 'logistic', 'normal', 'studentst')
base_variables_name <- unlist(lapply(base_variables_name, function(x) c(x, paste(x, ':treated', sep = ''))))
base_variables_name <- c(base_variables_name, 'treated')

#Center and Standardize covariates for propensity score calculation
X_std <- data.frame(apply(X, MARGIN = 2, FUN = function(x) (x-mean(x))/sd(x)))

#######################################################
####################### Effects #######################
#######################################################

#DGPs 1 to 4 heterogeneous treatment effects
beta_d_1 <- c(rep(1, times = N))
beta_d_1[which(X[,1] == 0)] <- 0.5
beta_d_2 <- beta_d_3 <- beta_d_4 <- beta_d_1

#Propensity scores for linear DGPs (DGPs 1, 3, 5 and 7)
prop_score_linear <- 1/(1+exp(-rowSums(X_std[,1:21])))

#DGPs 5 to 8 heterogeneous treatment effects
beta_d_5 <- floor(X[,bXt]/sd(X[,bXt]))
beta_d_6 <- beta_d_7 <- beta_d_8 <- beta_d_5


#DGP 1 coefficients
beta_p_1 <- c(rep(0, times = p))
beta_p_1[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_1)-min(beta_d_1))/2)
treated_1 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 2 coefficients
beta_p_2 <- c(rep(0, times = p))
beta_p_2 <- as.numeric(rbernoulli(p, p = 0.7))*rnorm(p, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))/2)
prop_score_2 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_2 != 0)])))
treated_2 <- as.numeric(rbernoulli(N, p = prop_score_2))

#DGP 3 coefficients
beta_p_3 <- c(rep(0, times = p))
beta_p_3[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))*2)
treated_3 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 4 coefficients
beta_p_4 <- c(rep(0, times = p))
beta_p_4 <- as.numeric(rbernoulli(p, p = 0.7))*rnorm(p, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))*2)
prop_score_4 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_4 != 0)])))
treated_4 <- as.numeric(rbernoulli(N, p = prop_score_4))

#DGP 5 coefficients
beta_p_5 <- c(rep(0, times = p))
beta_p_5[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))/2)
treated_5 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 6 coefficients
beta_p_6 <- c(rep(0, times = p))
beta_p_6 <- as.numeric(rbernoulli(p, p = 0.7))*rnorm(p, mean = 0, sd = (max(beta_d_6)-min(beta_d_6))/2)
prop_score_6 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_6 != 0)])))
treated_6 <- as.numeric(rbernoulli(N, p = prop_score_6))

#DGP 7 coefficients
beta_p_7 <- c(rep(0, times = p))
beta_p_7[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_7)-min(beta_d_7))*2)
treated_7 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 8 coefficients
beta_p_8 <- c(rep(0, times = p))
beta_p_8 <- as.numeric(rbernoulli(p, p = 0.7))*rnorm(p, mean = 0, sd = (max(beta_d_8)-min(beta_d_8))*2)
prop_score_8 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_8 != 0)])))
treated_8 <- as.numeric(rbernoulli(N, p = prop_score_8))


#Calculate outcomes Y per DGP
Y_1 <- 0.5 + as.matrix(X)%*%beta_p_1 + treated_1 * beta_d_1 + rnorm(N, mean = 0, sd = (max(beta_d_1)-min(beta_d_1))/10)
Y_2 <- 0.5 + as.matrix(X)%*%beta_p_2 + treated_2 * beta_d_2 + rnorm(N, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))/10)
Y_3 <- 0.5 + as.matrix(X)%*%beta_p_3 + treated_3 * beta_d_3 + rnorm(N, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))/10)
Y_4 <- 0.5 + as.matrix(X)%*%beta_p_4 + treated_4 * beta_d_4 + rnorm(N, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))/10)
Y_5 <- 0.5 + as.matrix(X)%*%beta_p_5 + treated_5 * beta_d_5 + rnorm(N, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))/10)
Y_6 <- 0.5 + as.matrix(X)%*%beta_p_6 + treated_6 * beta_d_6 + rnorm(N, mean = 0, sd = (max(beta_d_6)-min(beta_d_6))/10)
Y_7 <- 0.5 + as.matrix(X)%*%beta_p_7 + treated_7 * beta_d_7 + rnorm(N, mean = 0, sd = (max(beta_d_7)-min(beta_d_7))/10)
Y_8 <- 0.5 + as.matrix(X)%*%beta_p_8 + treated_8 * beta_d_8 + rnorm(N, mean = 0, sd = (max(beta_d_8)-min(beta_d_8))/10)


#######################################################
###################### Iterations #####################
#######################################################

#Randomly draw N pairs of outcomes and covariates plus treatment status
sample <- sample(x = 1:N, size = sample_size)


#### Super Learner ####
#Randomly split sample into 10 folds
fold10 <- sample
fold1 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold1)]
fold2 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold2)]
fold3 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold3)]
fold4 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold4)]
fold5 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold5)]
fold6 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold6)]
fold7 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold7)]
fold8 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold8)]
fold9 <- sample(x = fold10, size = floor(sample_size/10))
fold10 <- fold10[-which(fold10 %in% fold9)]

folds <- list(fold1, fold2, fold3, fold4, fold5,
              fold6, fold7, fold8, fold9, fold10)


#Out of sample predictions for Ys
#Columns are component methods M and Rows units i
Y_1_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_1_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_1_hat) <- sort(sample)
Y_2_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_2_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_2_hat) <- sort(sample)
Y_3_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_3_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_3_hat) <- sort(sample)
Y_4_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_4_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_4_hat) <- sort(sample)
Y_5_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_5_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_5_hat) <- sort(sample)
Y_6_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_6_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_6_hat) <- sort(sample)
Y_7_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_7_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_7_hat) <- sort(sample)
Y_8_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 9))
colnames(Y_8_hat) <- c('Elastic Net', 'KRLS', 'R-learner', 'SVM', 'FindIt',
                       'Causal Boosting', 'Causal Forest', 'BGLM', 'BCF')
rownames(Y_8_hat) <- sort(sample)

a <- model.matrix(~matrix(data = 1, ncol = 3, nrow = 5)*c(1, 1, 0, 0, 1))

### Y_1 ###
for(i in 1:10){
  
  training_units <- sort(unlist(folds[c(1:10)[-i]]))
  training_sample <- model.matrix(~as.matrix(X[training_units,])*treated_1[training_units])
  
  colnames(training_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X\\[training_units, \\]\\)'), '_1\\[training_units\\]')
  base_variables <- which(colnames(training_sample) %in% base_variables_name)
  
  D_training_sample <- treated_1[training_units]
  Y_training_sample <- Y_1[training_units]
  
  test_units <- sort(folds[[i]])
  test_sample <- model.matrix(~as.matrix(X[test_units,])*treated_1[test_units])
  colnames(test_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X\\[test_units, \\]\\)'), '_1\\[test_units\\]')
  
  X_test_sample <- X[test_units,]
  D_test_sample <- treated_1[test_units]

  "
  #### Elastic-Net ####
  EN_fit <- cv.glmnet(as.matrix(training_sample), Y_training_sample, type.measure = 'mse', alpha = .5)
  Y_1_hat[which(rownames(Y_1_hat) %in% test_units),'Elastic Net'] <- predict(EN_fit, s = EN_fit$lambda.1se, newx = as.matrix(test_sample))

  
  #### KRLS ####
  non_constant <- which(!apply(training_sample[,base_variables], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
  KRLS_fit <- krls(X = training_sample[,base_variables][,non_constant], y = Y_training_sample, derivative = FALSE)
  Y_1_hat[which(rownames(Y_1_hat) %in% test_units),'KRLS'] <- predict(KRLS_fit, newdata = test_sample[,base_variables][,non_constant])$fit
  
  
  #### R-Learner ####
  r_fit <- rboost(as.matrix(training_sample[,base_variables]), D_training_sample, Y_training_sample)
  r_pred <- predict(r_fit, as.matrix(X_test_sample[,base_variables]), tau_only = FALSE)
  Y_1_hat[which(rownames(Y_1_hat) %in% test_units[which(D_test_sample==1)]),'R-learner'] <- r_pred$mu1[which(D_test_sample==1)]
  Y_1_hat[which(rownames(Y_1_hat) %in% test_units[which(D_test_sample==0)]),'R-learner'] <- r_pred$mu0[which(D_test_sample==0)]
  "
  
  #### SVMs #### 
  
  #Not done
  "
  SVM_fit <- SMO(Y ~ ., data = data.frame(Y = Y_training_sample, training_sample[,base_variables]), control = Weka_control(M = TRUE))
  predict(SVM_fit, newdata = data.frame(test_sample[,base_variables]), type = 'probability')[,2] 
  "
  
  #### FindIt ####
  
  #Skipped cause Y has to be binary?
  
  
  #### Causal Boosting ####
  
  #Not yet ready for publication / GAMBoost not for this R version
  
  
  #### Random Forest ####
  
  CF_fit <- randomForest(y = Y_training_sample, x = training_sample[,base_variables])
  Y_1_hat[which(rownames(Y_1_hat) %in% test_units),'Causal Forest'] <- predict(CF_fit, newdata = test_sample[,base_variables])
  
  
  #### BGLM ####
  
  
  
  #### BCF ####
  
  
  
}



#all M component methods using all remaining D − 1-folds for training.


#Obtain N × M matrix for Y per DGP, whereby Yi,m stands for unit i’s out of sample prediction from method m.


#Regress true response on out of sample predictions (Yi = 􏰋M wmYi,m + εi) i=1 per DGP with constraints 􏰋Mi=1 wm = 1 and wm ≥ 0, whereby ε is an error term.


#Obtain solutions to constrained regression w􏰎 as weights for final Ensemble per DGP.


#Obtain estimates for heterogeneous treatment effects per DGP using the full sample by all nine component methods, following the procedures outlined throughout chapter 2.3 and OLS as a benchmark.


#Calculate heterogeneous treatment effect estimates per DGP by a naive Ensemble and by Super Learning with weights from step 11.


#Iterate through steps 3 to 14 with different sample sizes N.


#Obtain performance and convergence measures for the Ensemble Methods, all components and OLS per DGP.