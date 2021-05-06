# To do before starting simulation
# Run "nc -l 4000" in terminal

rm(list = ls(all.names = TRUE))
set.seed(0815)

#Set working device
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Import packages
pacman::p_load(purrr,extraDistr,poisbinom,actuar,circular,evd,rdetools,
               sets,glmnet,KRLS,mboost,devtools,stringr,randomForest,arm,
               BayesTree,bcf,fastDummies,pracma,quadprog,BBmisc,doParallel,
               schoolmath,benchmarkme)
#install_github('xnie/rlearner')
library(rlearner)

#Import custom functions
source('functions.R')

#######################################################
############## Set simulation parameters ##############
#######################################################

#Simulation size
N <- 10000L
sample_size <- 1000L #Only choose sample sizes which are multiples of f or amend code for sampling folds
iterations <- 1000L
f <- 10L #Folds for Super Learning
n_pseudo <- 100 #Number of pseudo variables

sim_pars <- list(N, sample_size, iterations, f)


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

#Create pseudo variables
pseudovars <- as.data.frame(do.call(cbind, lapply(1:n_pseudo,
                                                  function(x) ifelse(is.even(rep(x, times = N)),
                                                                     rgamma(N, x/n_pseudo, rate = ceiling(sqrt(x))),
                                                                     rbeta(N, x/n_pseudo, ceiling(sqrt(x)))))))
                                                                     
colnames(pseudovars) <- unlist(lapply(1:n_pseudo, function(x) paste('pseudo', x, sep = '_')))

#Bind all covariates to dataframe
X <- cbind(rvar, sincvar, expvar, crossvar)

#######################################################

#Replace Inf from sinc function with 0
X <- data.frame(apply(X, MARGIN = 2, function(c) unlist(lapply(c, function(z) ifelse(is.finite(z), z, 0)))))

#Parameters for later stages in code
p <- ncol(X)
bXt <- which(colnames(X) == 'binomialXstudentst')
base_variables_name <- c('bernoulli', 'rademacher', 'poissonbinomial', 'binomial',
                         'hypergeometric', 'geometric', 'logarithmic', 'poisson',
                         'beta', 'wrappedcauchy', 'wrappednormal', 'chisquared',
                         'exponential', 'gamma', 'lognormal', 'weibull', 'gumbel',
                         'laplace', 'logistic', 'normal', 'studentst', colnames(pseudovars))
base_variables_name <- unlist(lapply(base_variables_name, function(x) c(x, paste(x, ':treated', sep = ''))))
base_variables_name <- c(base_variables_name, 'treated')

#Center and Standardize covariates for propensity score calculation
X_std <- data.frame(apply(X, MARGIN = 2, FUN = function(x) (x-mean(x))/sd(x)))

#######################################################
####################### Effects #######################
#######################################################

#Propensity scores for linear DGPs (DGPs 1, 2 and 3)
prop_score_linear <- 1/(1+exp(-rowSums(X_std[,1:21])))

#Heterogeneous treatment effects for low heterogeneous dimensionality (DGPs 1, 2 and 4)
beta_d_1 <- c(rep(1, times = N))
beta_d_1[which(X[,1] == 0)] <- 0.5
beta_d_2 <- beta_d_4 <- beta_d_1

#Heterogeneous treatment effects for high heterogeneous dimensionality (DGPs 3 and 5)
beta_d_5 <- floor(X[,bXt]/sd(X[,bXt]))
beta_d_3 <- beta_d_5


#DGP 1 coefficients
beta_p_1 <- c(rep(0, times = p))
beta_p_1[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_1)-min(beta_d_1))/2)
treated_1 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 2 coefficients
beta_p_2 <- c(rep(0, times = p))
beta_p_2[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))*2)
treated_2 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 3 coefficients
beta_p_3 <- c(rep(0, times = p))
beta_p_3 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))/2)
prop_score_3 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_3 != 0)])))
treated_3 <- as.numeric(rbernoulli(N, p = prop_score_3))

#DGP 4 coefficients
beta_p_4 <- c(rep(0, times = p))
beta_p_4[1:21] <- rnorm(21, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))/2)
treated_4 <- as.numeric(rbernoulli(N, p = prop_score_linear))

#DGP 5 coefficients
beta_p_5 <- c(rep(0, times = p))
beta_p_5 <- as.numeric(rbernoulli(p, p = 0.3))*rnorm(p, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))*2)
prop_score_5 <- 1/(1+exp(-rowSums(X_std[,which(beta_p_5 != 0)])))
treated_5 <- as.numeric(rbernoulli(N, p = prop_score_5))


treateds <- list(treated_1, treated_2, treated_3, treated_4, treated_5)

beta_ps <- list(beta_p_1, beta_p_2, beta_p_3, beta_p_4, beta_p_5)

beta_ds <- list(beta_d_1, beta_d_2, beta_d_3, beta_d_4, beta_d_5)

prop_scores <- list(prop_score_linear, prop_score_3, prop_score_5)


#Calculate outcomes Y per DGP
Y_1 <- 0.5 + as.matrix(X)%*%beta_p_1 + treated_1 * beta_d_1 + rnorm(N, mean = 0, sd = (max(beta_d_1)-min(beta_d_1))/10)
Y_2 <- 0.5 + as.matrix(X)%*%beta_p_2 + treated_2 * beta_d_2 + rnorm(N, mean = 0, sd = (max(beta_d_2)-min(beta_d_2))/10)
Y_3 <- 0.5 + as.matrix(X)%*%beta_p_3 + treated_3 * beta_d_3 + rnorm(N, mean = 0, sd = (max(beta_d_3)-min(beta_d_3))/10)
Y_4 <- 0.5 + as.matrix(X)%*%beta_p_4 + treated_4 * beta_d_4 + rnorm(N, mean = 0, sd = (max(beta_d_4)-min(beta_d_4))/10)
Y_5 <- 0.5 + as.matrix(X)%*%beta_p_5 + treated_5 * beta_d_5 + rnorm(N, mean = 0, sd = (max(beta_d_5)-min(beta_d_5))/10)

Ys <- list(Y_1, Y_2, Y_3, Y_4, Y_5)


#Add pseudo variables to X
X <- cbind(X, pseudovars)

#######################################################
###################### Simulation #####################
#######################################################

#Take time of simulation process                                          
start_time <- Sys.time()

#Create socket for progress output
log.socket <- make.socket(port = 4000)

#Create folder for output
dir.create(file.path(getwd(), paste('/output/', start_time, sep = '')))

#Save coefficients for descriptive evidence
save(sim_pars, X, Ys, treateds, beta_ps, beta_ds, prop_scores,
     file = paste(getwd(), '/output/', start_time, '/coefs.RData', sep = ''))


#Setup for parallel computing
n.cores <- min(detectCores() - 1, floor(as.numeric(get_ram())/10^9))
my.cluster <- makeCluster(n.cores, type = 'PSOCK')#, setup_strategy = 'sequential')
print(my.cluster)

### Simulation iterations
registerDoParallel(cl = my.cluster)


output <- foreach(it = 1:iterations, .inorder = FALSE,
                  .packages = c('sets', 'glmnet', 'KRLS', 'mboost',
                                'devtools', 'stringr', 'randomForest',
                                'arm', 'BayesTree', 'bcf', 'fastDummies',
                                'pracma', 'quadprog', 'rlearner', 'BBmisc',
                                'foreach')) %dopar% {
  
  #Randomly draw N pairs of outcomes and covariates plus treatment status
  sample <- sample(x = 1:N, size = sample_size)
  
  #######################################################
  #################### Super Learner ####################
  #######################################################
  #Randomly split sample into f folds
  folds <- split(sample, rep(1:ceiling(length(sample)/f), each = length(sample)/f)[1:length(sample)])
  
  ### Low dimensional heterogeneities
  ### Y_1_hat to Y_3_hat ###
  Y_hats_one <- foreach(dgps = c(1:3), .inorder = FALSE,
                        .packages = c('sets', 'glmnet', 'KRLS', 'mboost',
                                      'devtools', 'stringr', 'randomForest',
                                      'arm', 'BayesTree', 'bcf', 'fastDummies',
                                      'pracma', 'quadprog', 'rlearner', 'BBmisc')) %do% {

    Log('Started Super Learning for DGP %d of iteration %d.', dgps, it)
                                      
    Y_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                               dimnames = list(sort(sample),
                                               c('ElasticNet', 'KRLS', 'RLearner',
                                                 'CausalForest', 'BGLM', 'BCF'))))
    treated <- treateds[[dgps]]
    
    for(fl in 1:f){
      
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
      
    }
    
    return(Y_hat)
    
  }
  
  
  ### High dimensional heterogeneities
  ### Y_4_hat and Y_5_hat ###
  #Add dummy variables for heterogeneous group association
  colnams <- c(colnames(X), 'hetero_factor')
  X_dummy <- cbind(X, floor(X[,'binomialXstudentst']))
  colnames(X_dummy) <- colnams
  X_dummy <- dummy_cols(X_dummy, select_columns = 'hetero_factor')
  X_dummy <- X_dummy[,-which(colnames(X_dummy) == 'hetero_factor')]
  
  #Add dummies to base variable description
  dum_vars <- which(str_detect(colnames(X_dummy), 'hetero_factor', negate = FALSE))
  
  Y_hats_two <- foreach(dgps = c(4, 5), .inorder = FALSE,
                        .packages = c('sets', 'glmnet', 'KRLS', 'mboost',
                                      'devtools', 'stringr', 'randomForest',
                                      'arm', 'BayesTree', 'bcf', 'fastDummies',
                                      'pracma', 'quadprog', 'rlearner', 'BBmisc')) %do% {
                                        
    Log('Started Super Learning for DGP %d of iteration %d.', dgps, it)                 

    Y_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                               dimnames = list(sort(sample),
                                               c('ElasticNet', 'KRLS', 'RLearner',
                                                 'CausalForest', 'BGLM', 'BCF'))))
    treated <- treateds[[dgps]]
    
    for(fl in 1:f){
      
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
      
    }
    
    return(Y_hat)     
  
  }

  Y_hats <- c(Y_hats_one, Y_hats_two)
  
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
  
  ### Low dimensional heterogeneities
  ### Y_1_hat to Y_3_hat ###
  Y_hats_one <- foreach(dgps = c(1:3), .inorder = FALSE,
                        .packages = c('sets', 'glmnet', 'KRLS', 'mboost',
                                      'devtools', 'stringr', 'randomForest',
                                      'arm', 'BayesTree', 'bcf', 'fastDummies',
                                      'pracma', 'quadprog', 'rlearner', 'BBmisc')) %do% {
                                        
    Log('Started Counterfactial Estimation for DGP %d of iteration %d.', dgps, it)
  
    Y_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                               dimnames = list(sort(sample),
                                               c('ElasticNet', 'KRLS', 'RLearner',
                                                 'CausalForest', 'BGLM', 'BCF'))))
    
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
    
    
    return(Y_hat)
    
  }
  
  ### High dimensional heterogeneities
  #Add dummy variables for heterogeneous group association
  colnams <- c(colnames(X), 'hetero_factor')
  X_dummy <- cbind(X, floor(X[,'binomialXstudentst']))
  colnames(X_dummy) <- colnams
  X_dummy <- dummy_cols(X_dummy, select_columns = 'hetero_factor')
  X_dummy <- X_dummy[,-which(colnames(X_dummy) == 'hetero_factor')]
  
  #Add dummies to base variable description
  dum_vars <- which(str_detect(colnames(X_dummy), 'hetero_factor', negate = FALSE))
  
  ### Y_4_hat and Y_5_hat ###
  Y_hats_two <- foreach(dgps = c(4, 5), .inorder = FALSE,
                        .packages = c('sets', 'glmnet', 'KRLS', 'mboost',
                                      'devtools', 'stringr', 'randomForest',
                                      'arm', 'BayesTree', 'bcf', 'fastDummies',
                                      'pracma', 'quadprog', 'rlearner', 'BBmisc')) %do% {
                                        
    Log('Started Counterfactial Estimation for DGP %d of iteration %d.', dgps, it)
                                        
    Y_hat <- data.frame(matrix(data = NA, nrow = sample_size, ncol = 6,
                               dimnames = list(sort(sample),
                                               c('ElasticNet', 'KRLS', 'RLearner',
                                                 'CausalForest', 'BGLM', 'BCF'))))
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
    
    return(Y_hat)
  
  }
  
  Y_hats <- c(Y_hats_one, Y_hats_two)
  
  #######################################################
  ################# Effect Calculation ##################
  #######################################################
  
  #Obtain effect estimates per technique
  taus <- mapply(function(y, y_h, d) apply(y_h, MARGIN = 2, FUN = function(x)  ifelse(d[sort(sample)] == 1, y[sort(sample)] - x, x - y[sort(sample)])),
                 Ys, Y_hats, treateds, SIMPLIFY = FALSE)
  
  #Obtain estimates of Naive Ensemble
  tau_EM_NE <- lapply(taus, function(x) as.matrix(x) %*% c(rep(1/ncol(x), times = ncol(x))))
  
  #Obtain Ensemble Estimates from Super Learner
  tau_EM_SL <- mapply(function(x, y) as.matrix(x) %*% y,
                      taus, wei, SIMPLIFY = FALSE)
  
  #Save estimates and real effects to matrix
  #Elements of list = DGPs
  out_it <- mapply(function(b, t, ne, sl) data.frame(EN = t[,1],
                                                     KRLS = t[,2],
                                                     RL = t[,3],
                                                     CF = t[,4],
                                                     BGLM = t[,5],
                                                     BCF = t[,6],
                                                     NE = ne,
                                                     SL = sl,
                                                     true = b[sort(sample)]),
                   beta_ds, taus, tau_EM_NE, tau_EM_SL, SIMPLIFY = FALSE)

  #Print operational information
  now <- Sys.time()
  diff <- as.numeric(now) - as.numeric(start_time)
  
  days <- floor(diff/(60*60*24))
  hours <- floor(diff/(60*60)-days*24)
  minutes <- floor(diff/60-hours*60-days*24*60)
  seconds <- round(diff-minutes*60-hours*60*60-days*24*60*60)
  
  total_time_est <- diff/it*iterations
  end_time_est <- as.numeric(start_time) + total_time_est
  
  elapsed <- paste('Time elapsed: ',
                   ifelse(days>1, paste(days, ' days, ', sep = ''), ''),
                   ifelse(days==1, paste(days, ' day, ', sep = ''), ''),
                   ifelse(hours>1, paste(hours, ' hours, ', sep = ''), ''),
                   ifelse(hours==1, paste(hours, ' hour, ', sep = ''), ''),
                   ifelse(minutes>1, paste(minutes, ' minutes, ', sep = ''), ''),
                   ifelse(minutes==1, paste(minutes, ' minute, ', sep = ''), ''),
                   seconds, ' seconds.',
                   sep = '')
  complete <- paste('Estimated time of completion: ',
                    as.POSIXct(end_time_est, origin = "1970-01-01"),
                    sep = '')
  
  Log('Finished iteration %d of %d.', it, iterations)
  Log('%d', elapsed)
  Log('%d', complete)
  
  return(out_it)
  
}

#Stop parallel cluster
stopCluster(cl = my.cluster)

#Combine output of all iterations into handier list
out_gdp1 <- data.frame(matrix(data = NA, nrow = iterations*sample_size, ncol = 9))
colnames(out_gdp1) <- c('EN', 'KRLS', 'RL', 'CF',
                        'BGLM', 'BCF', 'NE', 'SL', 'True Value')
out_gdp2 <- data.frame(matrix(data = NA, nrow = iterations*sample_size, ncol = 9))
colnames(out_gdp2) <- c('EN', 'KRLS', 'RL', 'CF',
                        'BGLM', 'BCF', 'NE', 'SL', 'True Value')
out_gdp3 <- data.frame(matrix(data = NA, nrow = iterations*sample_size, ncol = 9))
colnames(out_gdp3) <- c('EN', 'KRLS', 'RL', 'CF',
                        'BGLM', 'BCF', 'NE', 'SL', 'True Value')
out_gdp4 <- data.frame(matrix(data = NA, nrow = iterations*sample_size, ncol = 9))
colnames(out_gdp4) <- c('EN', 'KRLS', 'RL', 'CF',
                        'BGLM', 'BCF', 'NE', 'SL', 'True Value')
out_gdp5 <- data.frame(matrix(data = NA, nrow = iterations*sample_size, ncol = 9))
colnames(out_gdp5) <- c('EN', 'KRLS', 'RL', 'CF',
                        'BGLM', 'BCF', 'NE', 'SL', 'True Value')


for(i in 1:length(output)){
  
  start <- i * sample_size - sample_size + 1
  end <- i * sample_size
  
  iter_out <- output[[i]]
  out_gdp1[start:end,] <- iter_out[[1]]
  out_gdp2[start:end,] <- iter_out[[2]]
  out_gdp3[start:end,] <- iter_out[[3]]
  out_gdp4[start:end,] <- iter_out[[4]]
  out_gdp5[start:end,] <- iter_out[[5]]
  
}

#Write output to disk
save(out_gdp1, out_gdp2, out_gdp3, out_gdp4, out_gdp5,
     file = paste(getwd(), '/output/', start_time, '/output.RData', sep = ''))

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

