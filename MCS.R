#Next steps: 1. Include fake variables and add them to base variable description
# 2. Delete not used packages


rm(list = ls(all.names = TRUE))
set.seed(0815)

#Set working device
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import packages
pacman::p_load(purrr,extraDistr,poisbinom,actuar,circular,evd,rdetools,
               sets,glmnet,KRLS,mboost,devtools,stringr,randomForest,arm,
               BayesTree,bcf,fastDummies,pracma,quadprog,BBmisc)#,rJava,RWeka,SVMMatch,FindIt,GAMBoost)
#install_github('xnie/rlearner')
library(rlearner)

#Import custom functions
source('functions.R')


#######################################################
############## Set simulation parameters ##############
#######################################################
#Simulation size
N <- 1000L
sample_size <- 100L #Only choose sample sizes which are multiples of f or amend code for sampling folds
iterations <- 10L
f <- 10L #Folds for Super Learning

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
hypergeometric <- rhyper(N, N/3, 2*N/3, 25)

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
chisquared <- rchisq(N, 1)

#Exponential distribution
exponential <- rexp(N, rate = 10)

#Gamma distribution
gamma <- rgamma(N, 0.8, scale = 1/1.5)

#Log-normal distribution 
lognormal <- rlnorm(N, meanlog = -2, sdlog = 1)

#Weibull distribution
weibull <- rweibull(N, 0.9, scale = 1)

#Gumbel distribution
gumbel <- rgumbel(N, loc = -5, scale = 2)

#Laplace distribution
laplace <- rlaplace(N, mu = 10, sigma = 2)

#Logistic distribution
logistic <- rlogis(N, location = 0, scale = 2)

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

#Crossproducts
crossvar <- custom_crossprod(rvar)

#Sinc transformation
sincvar <- custom_sinc(rvar_nonbinary)

#Exponential transformation
expvar <- custom_exp(rvar_nonbinary)

#Bind all covariates to dataframe
X <- cbind(rvar, sincvar, expvar, crossvar)

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
beta_p_2 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))/2)
prop_score_2 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_2 != 0)])))
treated_2 <- as.numeric(rbernoulli(N, p = prop_score_2))

#DGP 3 coefficients
beta_p_3 <- c(rep(0, times = p))
beta_p_3[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))*2)
treated_3 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 4 coefficients
beta_p_4 <- c(rep(0, times = p))
beta_p_4 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))*2)
prop_score_4 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_4 != 0)])))
treated_4 <- as.numeric(rbernoulli(N, p = prop_score_4))

#DGP 5 coefficients
beta_p_5 <- c(rep(0, times = p))
beta_p_5[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))/2)
treated_5 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 6 coefficients
beta_p_6 <- c(rep(0, times = p))
beta_p_6 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_6)-min(beta_d_6))/2)
prop_score_6 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_6 != 0)])))
treated_6 <- as.numeric(rbernoulli(N, p = prop_score_6))

#DGP 7 coefficients
beta_p_7 <- c(rep(0, times = p))
beta_p_7[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_7)-min(beta_d_7))*2)
treated_7 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 8 coefficients
beta_p_8 <- c(rep(0, times = p))
beta_p_8 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_8)-min(beta_d_8))*2)
prop_score_8 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_8 != 0)])))
treated_8 <- as.numeric(rbernoulli(N, p = prop_score_8))

treateds <- list(treated_1, treated_2, treated_3, treated_4,
                 treated_5, treated_6, treated_7, treated_8)


#Calculate outcomes Y per DGP
Y_1 <- 0.5 + as.matrix(X)%*%beta_p_1 + treated_1 * beta_d_1 + rnorm(N, mean = 0, sd = (max(beta_d_1)-min(beta_d_1))/10)
Y_2 <- 0.5 + as.matrix(X)%*%beta_p_2 + treated_2 * beta_d_2 + rnorm(N, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))/10)
Y_3 <- 0.5 + as.matrix(X)%*%beta_p_3 + treated_3 * beta_d_3 + rnorm(N, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))/10)
Y_4 <- 0.5 + as.matrix(X)%*%beta_p_4 + treated_4 * beta_d_4 + rnorm(N, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))/10)
Y_5 <- 0.5 + as.matrix(X)%*%beta_p_5 + treated_5 * beta_d_5 + rnorm(N, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))/10)
Y_6 <- 0.5 + as.matrix(X)%*%beta_p_6 + treated_6 * beta_d_6 + rnorm(N, mean = 0, sd = (max(beta_d_6)-min(beta_d_6))/10)
Y_7 <- 0.5 + as.matrix(X)%*%beta_p_7 + treated_7 * beta_d_7 + rnorm(N, mean = 0, sd = (max(beta_d_7)-min(beta_d_7))/10)
Y_8 <- 0.5 + as.matrix(X)%*%beta_p_8 + treated_8 * beta_d_8 + rnorm(N, mean = 0, sd = (max(beta_d_8)-min(beta_d_8))/10)

Ys <- list(Y_1, Y_2, Y_3, Y_4,
           Y_5, Y_6, Y_7, Y_8)


#######################################################
###################### Simulation #####################
#######################################################

#Generate matrix to fill with output value
#Elements of list = DGPs, first column of elements = estimate, second = true value
output <- lapply(c(1:8), function(i) assign(paste('dgp_', i, sep = ''),
                                            list(EN = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 KRLS = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 RL = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 CF = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 BGLM = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 BCF = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 NE = matrix(data = NA, ncol = 2, nrow = sample_size * iterations),
                                                 SL = matrix(data = NA, ncol = 2, nrow = sample_size * iterations))))

#Take time of simulation process                                          
start_time <- Sys.time()

#Create folder for output
dir.create(file.path(getwd(), paste('/output/', start_time, sep = '')))

for(it in 1:iterations){
  
  #Print operational information
  now <- Sys.time()
  diff <- as.numeric(now) - as.numeric(start_time)
  
  days <- floor(diff/(60*60*24))
  hours <- floor(diff/(60*60)-days*24)
  minutes <- floor(diff/60-hours*60-days*24*60)
  seconds <- round(diff-minutes*60-hours*60*60-days*24*60*60)
  
  total_time_est <- diff/it*iterations
  end_time_est <- as.numeric(start_time) + total_time_est
  
  print(paste('Beginn iteration ', it, ' of ', iterations, '.', sep = ''))
  print(paste('Time elapsed: ',
              ifelse(days>1, paste(days, ' days, ', sep = ''), ''),
              ifelse(days==1, paste(days, ' day, ', sep = ''), ''),
              ifelse(hours>1, paste(hours, ' hours, ', sep = ''), ''),
              ifelse(hours==1, paste(hours, ' hour, ', sep = ''), ''),
              ifelse(minutes>1, paste(minutes, ' minutes, ', sep = ''), ''),
              ifelse(minutes==1, paste(minutes, ' minute, ', sep = ''), ''),
              seconds, ' seconds.',
              sep = ''))
  print(paste('Estimated time of completion: ',
              as.POSIXct(end_time_est, origin = "1970-01-01"),
              sep = ''))
  
  
  #Randomly draw N pairs of outcomes and covariates plus treatment status
  sample <- sample(x = 1:N, size = sample_size)
  
  
  #######################################################
  #################### Super Learner ####################
  #######################################################
  
  #Randomly split sample into f folds
  folds <- split(sample, rep(1:ceiling(length(sample)/f), each = length(sample)/f)[1:length(sample)])
  
  #DFs for out of sample predictions for Ys
  #Columns are component methods M and Rows units i
  Y_hats <- rep(list(data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                                       dimnames = list(sort(sample),
                                                       c('ElasticNet', 'KRLS', 'RLearner',
                                                         'CausalForest', 'BGLM', 'BCF'))))),
                times = 8)
  
  
  ### Y_1_hat to Y_4_hat ###
  for(dgps in 1:4){
    
    print(paste('Super Learning for DGP ', dgps,
                ' out of 8 of iteration ', it,
                ' out of ', iterations, '.', sep = ''))
    
    Y_hat <- Y_hats[[dgps]]
    treated <- treateds[[dgps]]
    
    for(fl in 1:f){
      
      print(paste('Fold ', fl, ' out of ', f,
                  ' of Super Learning for DGP ', dgps,
                  ' out of 8 of iteration ', it,
                  ' out of ', iterations, '.', sep = ''))
      
      training_units <- sort(unlist(folds[c(1:10)[-fl]]))
      training_sample <- model.matrix(~as.matrix(X[training_units,])*treated[training_units])
      
      colnames(training_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X\\[training_units, \\]\\)'), '_1\\[training_units\\]')
      base_variables <- which(colnames(training_sample) %in% base_variables_name)
      
      D_training_sample <- treated[training_units]
      Y_training_sample <- Ys[[dgps]][training_units]
      
      test_units <- sort(folds[[fl]])
      test_sample <- model.matrix(~as.matrix(X[test_units,])*treated[test_units])
      colnames(test_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X\\[test_units, \\]\\)'), '_1\\[test_units\\]')
      
      X_test_sample <- X[test_units,]
      D_test_sample <- treated[test_units]
      
      
      #### Elastic-Net ####
      EN_fit <- cv.glmnet(as.matrix(training_sample), Y_training_sample, type.measure = 'mse', alpha = .5)
      Y_hat[which(rownames(Y_hat) %in% test_units),'ElasticNet'] <- predict(EN_fit, s = EN_fit$lambda.1se, newx = as.matrix(test_sample))
      
      
      #### KRLS ####
      non_constant <- which(!apply(training_sample[,base_variables], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
      KRLS_fit <- krls(X = training_sample[,base_variables][,non_constant], y = Y_training_sample, derivative = FALSE)
      Y_hat[which(rownames(Y_hat) %in% test_units),'KRLS'] <- predict(KRLS_fit, newdata = test_sample[,base_variables][,non_constant])$fit
      
      
      #### R-Learner ####
      remove(r_fit)
      r_fit <- rboost(as.matrix(training_sample[,base_variables]), D_training_sample, Y_training_sample)
      r_pred <- predict(r_fit, as.matrix(test_sample[,base_variables]), tau_only = FALSE)
      Y_hat[which(rownames(Y_hat) %in% test_units[which(D_test_sample==1)]),'RLearner'] <- r_pred$mu1[which(D_test_sample==1)]
      Y_hat[which(rownames(Y_hat) %in% test_units[which(D_test_sample==0)]),'RLearner'] <- r_pred$mu0[which(D_test_sample==0)]
      
      
      #### Random Forest ####
      CF_fit <- randomForest(y = Y_training_sample, x = training_sample[,base_variables])
      Y_hat[which(rownames(Y_hat) %in% test_units),'CausalForest'] <- predict(CF_fit, newdata = test_sample[,base_variables])
      
      
      #### BGLM ####
      #Check if base variables are indeed the correct choice
      BGLM_fit <- bayesglm(Y_training_sample ~ training_sample[,base_variables])
      Y_hat[which(rownames(Y_hat) %in% test_units),'BGLM'] <- test_sample[,c(1, base_variables)] %*% BGLM_fit$coef
      
      
      #### BCF ####
      BART_fit <- bart(x.train = training_sample[,base_variables], y.train = Y_training_sample,
                       x.test = test_sample[,base_variables], ndpost = 1000, nskip = 500, usequants = T)
      Y_hat[which(rownames(Y_hat) %in% test_units),'BCF'] <- colMeans(BART_fit$yhat.test)
      
      gc()
      
    }
    
    Y_hats[[dgps]] <- Y_hat
    
  }
  
  
  ### Y_5_hat to Y_8_hat ###
  
  #Add dummy variables for heterogeneous group association
  colnams <- c(colnames(X), 'hetero_factor')
  X_dummy <- cbind(X, floor(X[,'binomialXstudentst']))
  colnames(X_dummy) <- colnams
  X_dummy <- dummy_cols(X_dummy, select_columns = 'hetero_factor')
  X_dummy <- X_dummy[,-which(colnames(X_dummy) == 'hetero_factor')]
  
  #Add dummies to base variable description
  dum_vars <- which(str_detect(colnames(X_dummy), 'hetero_factor', negate = FALSE))
  
  for(dgps in 5:8){
    
    print(paste('Super Learning for DGP ', dgps,
                ' out of 8 of iteration ', it,
                ' out of ', iterations, '.', sep = ''))
    
    Y_hat <- Y_hats[[dgps]]
    treated <- treateds[[dgps]]
    
    for(fl in 1:f){
      
      print(paste('Fold ', fl, ' out of ', f,
                  ' of Super Learning for DGP ', dgps,
                  ' out of 8 of iteration ', it,
                  ' out of ', iterations, '.', sep = ''))
      
      training_units <- sort(unlist(folds[c(1:10)[-fl]]))
      training_sample <- model.matrix(~as.matrix(X_dummy[training_units,])*treated[training_units])
      
      colnames(training_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X_dummy\\[training_units, \\]\\)'), '_1\\[training_units\\]')
      base_variables <- c(which(colnames(training_sample) %in% base_variables_name), dum_vars)
      
      D_training_sample <- treated[training_units]
      Y_training_sample <- Ys[[dgps]][training_units]
      
      test_units <- sort(folds[[fl]])
      test_sample <- model.matrix(~as.matrix(X_dummy[test_units,])*treated[test_units])
      colnames(test_sample) <- str_remove(str_remove(colnames(training_sample), 'as.matrix\\(X_dummy\\[test_units, \\]\\)'), '_1\\[test_units\\]')
      
      X_test_sample <- X_dummy[test_units,]
      D_test_sample <- treated[test_units]
      
      
      #### Elastic-Net ####
      EN_fit <- cv.glmnet(as.matrix(training_sample), Y_training_sample, type.measure = 'mse', alpha = .5)
      Y_hat[which(rownames(Y_hat) %in% test_units),'ElasticNet'] <- predict(EN_fit, s = EN_fit$lambda.1se, newx = as.matrix(test_sample))
      
      
      #### KRLS ####
      non_constant <- which(!apply(training_sample[,base_variables], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
      KRLS_fit <- krls(X = training_sample[,base_variables][,non_constant], y = Y_training_sample, derivative = FALSE)
      Y_hat[which(rownames(Y_hat) %in% test_units),'KRLS'] <- predict(KRLS_fit, newdata = test_sample[,base_variables][,non_constant])$fit
      
      
      #### R-Learner ####
      r_fit <- rboost(as.matrix(training_sample[,base_variables]), D_training_sample, Y_training_sample)
      r_pred <- predict(r_fit, as.matrix(test_sample[,base_variables]), tau_only = FALSE)
      Y_hat[which(rownames(Y_hat) %in% test_units[which(D_test_sample==1)]),'RLearner'] <- r_pred$mu1[which(D_test_sample==1)]
      Y_hat[which(rownames(Y_hat) %in% test_units[which(D_test_sample==0)]),'RLearner'] <- r_pred$mu0[which(D_test_sample==0)]
      
      
      #### Random Forest ####
      CF_fit <- randomForest(y = Y_training_sample, x = training_sample[,base_variables])
      Y_hat[which(rownames(Y_hat) %in% test_units),'CausalForest'] <- predict(CF_fit, newdata = test_sample[,base_variables])
      
      
      #### BGLM ####
      #Check if base variables are indeed the correct choice
      BGLM_fit <- bayesglm(Y_training_sample ~ training_sample[,base_variables])
      Y_hat[which(rownames(Y_hat) %in% test_units),'BGLM'] <- test_sample[,c(1, base_variables)] %*% BGLM_fit$coef
      
      
      #### BCF ####
      BART_fit <- bart(x.train = training_sample[,base_variables], y.train = Y_training_sample,
                       x.test = test_sample[,base_variables], ndpost = 1000, nskip = 500, usequants = T)
      Y_hat[which(rownames(Y_hat) %in% test_units),'BCF'] <- colMeans(BART_fit$yhat.test)
      
      gc()
      
    }
    
    Y_hats[[dgps]] <- Y_hat
    
  }
  
  #Obtain weights from Super Learning (needs to be scaled for technical reasons for large Y_i)
  wei <- try(mapply(function(x, y) lsqlincon(as.matrix(x), y[sort(sample)],
                                           Aeq = matrix(rep(1, ncol(x)), nrow = 1),
                                           beq = c(1),
                                           lb = rep(0, ncol(x)), ub = rep(1, ncol(x))),
                  Y_hats, Ys, SIMPLIFY = FALSE))  
  
  scale <- 1
  while(is.error(wei)){
    
    wei <- try(mapply(function(x, y) lsqlincon(as.matrix(x/(10^scale)), y[sort(sample)]/(10^scale),
                                               Aeq = matrix(rep(1, ncol(x)), nrow = 1),
                                               beq = c(1),
                                               lb = rep(0, ncol(x)), ub = rep(1, ncol(x))),
                      Y_hats, Ys, SIMPLIFY = FALSE))  
    
    scale <- scale + 1
    
  }
  
  
  #######################################################
  ############## Counterfactial Estimation ##############
  #######################################################
  
  Y_hats <- rep(list(data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                                       dimnames = list(sort(sample),
                                                       c('ElasticNet', 'KRLS', 'RLearner',
                                                         'CausalForest', 'BGLM', 'BCF'))))),
                times = 8)
  
  ### Y_1_hat to Y_4_hat ###
  for(dgps in 1:4){
    
    print(paste('Counterfactual estimation for DGP ', dgps,
                ' out of 8 of iteration ', it,
                ' out of ', iterations, '.', sep = ''))
    
    Y_hat <- Y_hats[[dgps]]
    treated <- treateds[[dgps]]
    
    
    X_sample <- model.matrix(~as.matrix(X[sort(sample),])*treated[sort(sample)])
    
    colnames(X_sample) <- str_remove(str_remove(colnames(X_sample), 'as.matrix\\(X\\[sort\\(sample\\), \\]\\)'), '_1\\[sort\\(sample\\)\\]')
    base_variables <- which(colnames(X_sample) %in% base_variables_name)
    
    D_sample <- treated[sort(sample)]
    Y_sample <- Ys[[dgps]][sort(sample)]
    
    counterfactual_sample <- model.matrix(~as.matrix(X[sort(sample),])*+(!treated[sort(sample)]))
    colnames(counterfactual_sample) <- colnames(X_sample)
    
    X_counterfactual_sample <- X[sort(sample),]
    D_counterfactual_sample <- +(!treated[sort(sample)])
    
    
    #### Elastic-Net ####
    EN_fit <- cv.glmnet(as.matrix(X_sample), Y_sample, type.measure = 'mse', alpha = .5)
    Y_hat[,'ElasticNet'] <- predict(EN_fit, s = EN_fit$lambda.1se, newx = as.matrix(counterfactual_sample))
    
    
    #### KRLS ####
    non_constant <- which(!apply(X_sample[,base_variables], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
    KRLS_fit <- krls(X = X_sample[,base_variables][,non_constant], y = Y_sample, derivative = FALSE)
    Y_hat[,'KRLS'] <- predict(KRLS_fit, newdata = counterfactual_sample[,base_variables][,non_constant])$fit
    
    
    #### R-Learner ####
    r_fit <- rboost(as.matrix(X_sample[,base_variables]), D_sample, Y_sample)
    r_pred <- predict(r_fit, as.matrix(X_sample[,base_variables]), tau_only = FALSE)
    Y_hat[which(rownames(Y_hat) %in% sample[which(D_sample==1)]),'RLearner'] <- r_pred$mu0[which(D_sample==1)]
    Y_hat[which(rownames(Y_hat) %in% sample[which(D_sample==0)]),'RLearner'] <- r_pred$mu1[which(D_sample==0)]
    
    
    #### Random Forest ####
    CF_fit <- randomForest(y = Y_sample, x = X_sample[,base_variables])
    Y_hat[,'CausalForest'] <- predict(CF_fit, newdata = counterfactual_sample[,base_variables])
    
    
    #### BGLM ####
    BGLM_fit <- bayesglm(Y_sample ~ X_sample[,base_variables])
    Y_hat[,'BGLM'] <- counterfactual_sample[,c(1, base_variables)] %*% BGLM_fit$coef
    
    
    #### BCF ####
    BART_fit <- bart(x.train = X_sample[,base_variables], y.train = Y_sample,
                     x.test = counterfactual_sample[,base_variables], ndpost = 1000, nskip = 500, usequants = T)
    Y_hat[,'BCF'] <- colMeans(BART_fit$yhat.test)
    
    Y_hats[[dgps]] <- Y_hat
    
    gc()
    
  }
  
  
  ### Y_5_hat to Y_8_hat ###
  
  #Add dummy variables for heterogeneous group association
  colnams <- c(colnames(X), 'hetero_factor')
  X_dummy <- cbind(X, floor(X[,'binomialXstudentst']))
  colnames(X_dummy) <- colnams
  X_dummy <- dummy_cols(X_dummy, select_columns = 'hetero_factor')
  X_dummy <- X_dummy[,-which(colnames(X_dummy) == 'hetero_factor')]
  
  #Add dummies to base variable description
  dum_vars <- which(str_detect(colnames(X_dummy), 'hetero_factor', negate = FALSE))
  
  ### Y_5_hat to Y_8_hat ###
  for(dgps in 5:8){
    
    print(paste('Counterfactual estimation for DGP ', dgps,
                ' out of 8 of iteration ', it,
                ' out of ', iterations, '.', sep = ''))
    
    Y_hat <- Y_hats[[dgps]]
    treated <- treateds[[dgps]]
    
    
    X_sample <- model.matrix(~as.matrix(X_dummy[sort(sample),])*treated[sort(sample)])
    
    colnames(X_sample) <- str_remove(str_remove(colnames(X_sample), 'as.matrix\\(X_dummy\\[sort\\(sample\\), \\]\\)'), '_1\\[sort\\(sample\\)\\]')
    base_variables <- c(which(colnames(X_sample) %in% base_variables_name), dum_vars)
    
    D_sample <- treated[sort(sample)]
    Y_sample <- Ys[[dgps]][sort(sample)]
    
    counterfactual_sample <- model.matrix(~as.matrix(X_dummy[sort(sample),])*+(!treated[sort(sample)]))
    colnames(counterfactual_sample) <- colnames(X_sample)
    
    X_counterfactual_sample <- X_dummy[sort(sample),]
    D_counterfactual_sample <- +(!treated[sort(sample)])
    
    
    #### Elastic-Net ####
    EN_fit <- cv.glmnet(as.matrix(X_sample), Y_sample, type.measure = 'mse', alpha = .5)
    Y_hat[,'ElasticNet'] <- predict(EN_fit, s = EN_fit$lambda.1se, newx = as.matrix(counterfactual_sample))
    
    
    #### KRLS ####
    non_constant <- which(!apply(X_sample[,base_variables], MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE)))
    KRLS_fit <- krls(X = X_sample[,base_variables][,non_constant], y = Y_sample, derivative = FALSE)
    Y_hat[,'KRLS'] <- predict(KRLS_fit, newdata = counterfactual_sample[,base_variables][,non_constant])$fit
    
    
    #### R-Learner ####
    r_fit <- rboost(as.matrix(X_sample[,base_variables]), D_sample, Y_sample)
    r_pred <- predict(r_fit, as.matrix(X_sample[,base_variables]), tau_only = FALSE)
    Y_hat[which(rownames(Y_hat) %in% sample[which(D_sample==1)]),'RLearner'] <- r_pred$mu0[which(D_sample==1)]
    Y_hat[which(rownames(Y_hat) %in% sample[which(D_sample==0)]),'RLearner'] <- r_pred$mu1[which(D_sample==0)]
    
    
    #### Random Forest ####
    CF_fit <- randomForest(y = Y_sample, x = X_sample[,base_variables])
    Y_hat[,'CausalForest'] <- predict(CF_fit, newdata = counterfactual_sample[,base_variables])
    
    
    #### BGLM ####
    BGLM_fit <- bayesglm(Y_sample ~ X_sample[,base_variables])
    Y_hat[,'BGLM'] <- counterfactual_sample[,c(1, base_variables)] %*% BGLM_fit$coef
    
    
    #### BCF ####
    BART_fit <- bart(x.train = X_sample[,base_variables], y.train = Y_sample,
                     x.test = counterfactual_sample[,base_variables], ndpost = 1000, nskip = 500, usequants = T)
    Y_hat[,'BCF'] <- colMeans(BART_fit$yhat.test)
    
    Y_hats[[dgps]] <- Y_hat
    
    gc()
    
  }
  
  save(Ys, Y_hats, wei, file = paste(getwd(), '/data.RData', sep = ''))
  
  ##
  #Load priliminarily data for development process
  load(paste(getwd(), '/data.RData', sep = ''))
  
  #######################################################
  ################# Effect Calculation ##################
  #######################################################
  
  #Obtain effect estimates per techique
  taus <- rep(list(data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                                     dimnames = list(sort(sample),
                                                     c('ElasticNet', 'KRLS', 'RLearner',
                                                       'CausalForest', 'BGLM', 'BCF'))))),
              times = 8)
  
  for(i in 1:8){
    
    observed <- Ys[[i]][sort(sample)]
    counter <- Y_hats[[i]]
    treatment <- treateds[[i]][sort(sample)]
    tau <- taus[[i]]
    
    for(j in 1:ncol(tau)){
      
      tau[,j] <- ifelse(treatment==1,
                        observed - counter[,j],
                        counter[,j] - observed)
      
    }
    
    taus[[i]] <- tau
    
  }
  
  #Obtain estimates of Naive Ensemble
  tau_EM_NE <- lapply(taus, function(x) as.matrix(x) %*% c(rep(1/ncol(x), times = ncol(x))))
  
  #Obtain Ensemble Estimates from Super Learner
  tau_EM_SL <- mapply(function(x, y) as.matrix(x) %*% y,
                      taus, wei, SIMPLIFY = FALSE)
  
  #Save estimates and real effects to matrix
  for(i in 1:length(output)){
    
    row_start <- iterations * sample_size - sample_size + 1
    row_end <- iterations * sample_size
    
    #Receive new information
    dgp <- output[[i]]
    tau <- taus[[i]]
    EM_NE <- tau_EM_NE[[i]]
    EM_SL <- tau_EM_SL[[i]]
    
    #Receive previous information
    EN <- dgp[[1]]
    KRLS <- dgp[[2]]
    RL <- dgp[[3]]
    CF <- dgp[[4]]
    BGLM <- dgp[[5]]
    BCF <- dgp[[6]]
    NE <- dgp[[7]]
    SL <- dgp[[8]]
    
    #Add new information
    EN[row_start:row_end,1] <- tau[,1]
    KRLS[row_start:row_end,1] <- tau[,2]
    RL[row_start:row_end,1] <- tau[,3]
    CF[row_start:row_end,1] <- tau[,4]
    BGLM[row_start:row_end,1] <- tau[,5]
    BCF[row_start:row_end,1] <- tau[,6]
    NE[row_start:row_end,1] <- EM_NE
    SL[row_start:row_end,1] <- EM_SL
    
    EN[row_start:row_end,2] <- beta_d_1[sort(sample)]
    KRLS[row_start:row_end,2] <- beta_d_2[sort(sample)]
    RL[row_start:row_end,2] <- beta_d_3[sort(sample)]
    CF[row_start:row_end,2] <- beta_d_4[sort(sample)]
    BGLM[row_start:row_end,2] <- beta_d_5[sort(sample)]
    BCF[row_start:row_end,2] <- beta_d_6[sort(sample)]
    NE[row_start:row_end,2] <- beta_d_7[sort(sample)]
    SL[row_start:row_end,2] <- beta_d_8[sort(sample)]
    
    #Save updated information
    dgp[[1]] <- EN
    dgp[[2]] <- KRLS
    dgp[[3]] <- RL
    dgp[[4]] <- CF
    dgp[[5]] <- BGLM
    dgp[[6]] <- BCF
    dgp[[7]] <- NE
    dgp[[8]] <- SL
    
    output[[i]] <- dgp
    
  }
  
  #Write current status to disk
  save(output, file = paste(getwd(), '/output/', start_time, '/output.RData', sep = '/'))
  
  gc()
  
}


#Print final operational information
now <- Sys.time()
diff <- as.numeric(now) - as.numeric(start_time)

days <- floor(diff/(60*60*24))
hours <- floor(diff/(60*60)-days*24)
minutes <- floor(diff/60-hours*60-days*24*60)
seconds <- round(diff-minutes*60-hours*60*60-days*24*60*60)

print('Simulation done.')
print(paste('Total number of iterations: ', iterations, sep = ''))
print(paste('Total time elapsed: ',
            ifelse(days>1, paste(days, ' days, ', sep = ''), ''),
            ifelse(days==1, paste(days, ' day, ', sep = ''), ''),
            ifelse(hours>1, paste(hours, ' hours, ', sep = ''), ''),
            ifelse(hours==1, paste(hours, ' hour, ', sep = ''), ''),
            ifelse(minutes>1, paste(minutes, ' minutes, ', sep = ''), ''),
            ifelse(minutes==1, paste(minutes, ' minute, ', sep = ''), ''),
            seconds, ' seconds.',
            sep = ''))