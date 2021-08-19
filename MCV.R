#Monte Carlo Visualisation
rm(list = ls(all.names = TRUE))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

pacman::p_load(ggplot2,vtable,moments,reshape2,ggallin,xtable,matlib,viridis)

#Import custom functions
source('functions.R')

#Import data
dirs <- list.dirs(path = './output', full.names = TRUE, recursive = TRUE)
dirs <- dirs[-1]

weight_dir <- dirs[which(dirs == './output/weight sample')]
if('./output/weight sample' %in% dirs){
  
  dirs <- dirs[-which(dirs == './output/weight sample')]
  
}

#################### Variable analysis ####################
#Theme for graphs
nt_color <- '#1a9850'
t_color <- '#d73027'

own_theme <- theme_classic() +
  theme(text = element_text(size = 18,  family = 'Arial'))

for(i in dirs){
  
  load(paste(getwd(), sub('.', '', i), '/output.RData', sep = ''))
  load(paste(getwd(), sub('.', '', i), '/coefs.RData', sep = ''))
  
  samplesize <- sim_pars[[2]]
  
  #Create folder for output
  outputfolder <- paste(getwd(), '/analysis/', samplesize, 'n', sep = '')
  dir.create(outputfolder, recursive = TRUE)
  
  #Correct colnames of IATEs
  colnames(out_gdp1) <- c(colnames(out_gdp1)[1:3], 'RF', 'BGLM', 'BART', colnames(out_gdp1)[7:9])
  colnames(out_gdp2) <- c(colnames(out_gdp1)[1:3], 'RF', 'BGLM', 'BART', colnames(out_gdp1)[7:9])
  colnames(out_gdp3) <- c(colnames(out_gdp1)[1:3], 'RF', 'BGLM', 'BART', colnames(out_gdp1)[7:9])
  colnames(out_gdp4) <- c(colnames(out_gdp1)[1:3], 'RF', 'BGLM', 'BART', colnames(out_gdp1)[7:9])
  colnames(out_gdp5) <- c(colnames(out_gdp1)[1:3], 'RF', 'BGLM', 'BART', colnames(out_gdp1)[7:9])
  
  #################### Effect estimates analysis ####################
  ########## Aggregate GATEs from IATEs per iteration ##########
  out_gdp1 <- c.gate(out_gdp1)
  out_gdp2 <- c.gate(out_gdp2)
  out_gdp3 <- c.gate(out_gdp3)
  out_gdp4 <- c.gate(out_gdp4)
  out_gdp5 <- c.gate(out_gdp5)
  
  ########## Performance ##########
  perf_dgp1 <- perf(out_gdp1)
  perf_dgp2 <- perf(out_gdp2)
  perf_dgp3 <- perf(out_gdp3)
  perf_dgp4 <- perf(out_gdp4)
  perf_dgp5 <- perf(out_gdp5)
  
  if(max(abs(perf_dgp1))<999) print(xtable(perf_dgp1, digits = 2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp1', '.tex', sep = '')) else print(xtable(perf_dgp1, digits = -2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp1', '.tex', sep = ''))
  if(max(abs(perf_dgp2))<999) print(xtable(perf_dgp2, digits = 2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp2', '.tex', sep = '')) else print(xtable(perf_dgp2, digits = -2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp2', '.tex', sep = ''))
  if(max(abs(perf_dgp3))<999) print(xtable(perf_dgp3, digits = 2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp3', '.tex', sep = '')) else print(xtable(perf_dgp3, digits = -2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp3', '.tex', sep = ''))
  if(max(abs(perf_dgp4))<999) print(xtable(perf_dgp4, digits = 2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp4', '.tex', sep = '')) else print(xtable(perf_dgp4, digits = -2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp4', '.tex', sep = ''))
  if(max(abs(perf_dgp5))<999) print(xtable(perf_dgp5, digits = 2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp5', '.tex', sep = '')) else print(xtable(perf_dgp5, digits = -2, type = 'latex'), file = paste(outputfolder, '/performancetable_dgp5', '.tex', sep = ''))
  
  #Rank overview
  rank_overview <- cbind(colnames(perf_dgp1)[order(rank(perf_dgp1[1,]))],
                         colnames(perf_dgp2)[order(rank(perf_dgp2[1,]))],
                         colnames(perf_dgp3)[order(rank(perf_dgp3[1,]))],
                         colnames(perf_dgp4)[order(rank(perf_dgp4[1,]))],
                         colnames(perf_dgp5)[order(rank(perf_dgp5[1,]))])
                         
  colnames(rank_overview) <- c('DGP 1', 'DGP 2', 'DGP 3', 'DGP 4', 'DGP 5')
  
  print(xtable(rank_overview, type = 'latex'), file = paste(outputfolder, '/rank_overview', '.tex', sep = ''))

  ########## Boxplots ##########
  c.boxplot(out_gdp1, 1)
  c.sqrt_boxplot(out_gdp1, 1)
  c.log_boxplot(out_gdp1, 1)
  
  c.boxplot(out_gdp2, 2)
  c.sqrt_boxplot(out_gdp2, 2)
  c.log_boxplot(out_gdp2, 2)
  
  c.boxplot(out_gdp3, 3)
  c.sqrt_boxplot(out_gdp3, 3)
  c.log_boxplot(out_gdp3, 3)
  
  c.boxplot(out_gdp4, 4)
  c.sqrt_boxplot(out_gdp4, 4)
  c.log_boxplot(out_gdp4, 4)
  
  c.boxplot(out_gdp5, 5)
  c.sqrt_boxplot(out_gdp5, 5)
  c.log_boxplot(out_gdp5, 5)
  
  ########## Rank analysis ##########
  #Calculate rank per iteration
  rank_dgp1 <- c.rank(out_gdp1)
  rank_dgp2 <- c.rank(out_gdp2)
  rank_dgp3 <- c.rank(out_gdp3)
  rank_dgp4 <- c.rank(out_gdp4)
  rank_dgp5 <- c.rank(out_gdp5)
  
  #Rank table
  ranktable_dgp1 <- ranktable(rank_dgp1)
  ranktable_dgp2 <- ranktable(rank_dgp2)
  ranktable_dgp3 <- ranktable(rank_dgp3)
  ranktable_dgp4 <- ranktable(rank_dgp4)
  ranktable_dgp5 <- ranktable(rank_dgp5)
  
  print(xtable(ranktable_dgp1, type = 'latex'), file = paste(outputfolder, '/ranktable_dgp1', '.tex', sep = ''))
  print(xtable(ranktable_dgp2, type = 'latex'), file = paste(outputfolder, '/ranktable_dgp2', '.tex', sep = ''))
  print(xtable(ranktable_dgp3, type = 'latex'), file = paste(outputfolder, '/ranktable_dgp3', '.tex', sep = ''))
  print(xtable(ranktable_dgp4, type = 'latex'), file = paste(outputfolder, '/ranktable_dgp4', '.tex', sep = ''))
  print(xtable(ranktable_dgp5, type = 'latex'), file = paste(outputfolder, '/ranktable_dgp5', '.tex', sep = ''))
  
  #Rankcorrelation
  rankcortable_dgp1 <- c.rankcor(rank_dgp1)
  rankcortable_dgp2 <- c.rankcor(rank_dgp2)
  rankcortable_dgp3 <- c.rankcor(rank_dgp3)
  rankcortable_dgp4 <- c.rankcor(rank_dgp4)
  rankcortable_dgp5 <- c.rankcor(rank_dgp5)
  
  print(xtable(rankcortable_dgp1, type = 'latex'), file = paste(outputfolder, '/rankcortable_dgp1', '.tex', sep = ''))
  print(xtable(rankcortable_dgp2, type = 'latex'), file = paste(outputfolder, '/rankcortable_dgp2', '.tex', sep = ''))
  print(xtable(rankcortable_dgp3, type = 'latex'), file = paste(outputfolder, '/rankcortable_dgp3', '.tex', sep = ''))
  print(xtable(rankcortable_dgp4, type = 'latex'), file = paste(outputfolder, '/rankcortable_dgp4', '.tex', sep = ''))
  print(xtable(rankcortable_dgp5, type = 'latex'), file = paste(outputfolder, '/rankcortable_dgp5', '.tex', sep = ''))
  
  #Colored bar on rank
  rankchart <- function(df, dgp_num){
    
    df_long <- melt(df)
    df_long$iteration <- rep(1:nrow(df), times = ncol(df))
    
    plot <- ggplot(df_long, aes(x = iteration, y = variable)) +
      geom_tile(aes(fill = value)) +
      own_theme +
      scale_fill_distiller(palette = 'RdYlGn', direction = -1,
                           guide = 'legend', name = 'Rank',
                           breaks = c(1, 4, max(df_long$value))) +
      scale_y_discrete(limits = rev) +
      scale_x_discrete() +
      labs(title = paste('Relative rank of estimators in DGP ', dgp_num, sep = ''),
           subtitle = paste('Sample size: ', samplesize, sep = ''),
           x = 'Iterations', y = '')
    
    png(file = paste(outputfolder, '/rank_dgp', dgp_num, '.png', sep = ''),
        width = 1200, height = 600)
    print(plot)
    dev.off()
    
  }
  
  rankchart(rank_dgp1, 1)
  rankchart(rank_dgp2, 2)
  rankchart(rank_dgp3, 3)
  rankchart(rank_dgp4, 4)
  rankchart(rank_dgp5, 5)
  
}


################Covariates & Co

#Create folder for output
outputfolder <- paste(getwd(), '/analysis/descriptive', sep = '')
dir.create(outputfolder, recursive = TRUE)


########## Summary stats ##########
sumtable(X, out = 'latex',
         file = paste(outputfolder, '/sumtable_X', sep = ''))

Y <- as.data.frame(matrix(data = unlist(Ys), ncol = 5))
colnames(Y) <- as.character(1:5)
sumtable(Y, out = 'latex',
         file = paste(outputfolder, '/sumtable_Y', sep = ''))

p <- as.data.frame(matrix(data = unlist(prop_scores), ncol = 3))
colnames(p) <- as.character(1:3)
sumtable(p, out = 'latex',
         file = paste(outputfolder, '/sumtable_propscores', sep = ''))

########## Histograms ##########
#Create histograms of X variables

for(i in 1:ncol(X)){
  
  if(length(unique(X[,i])) < 5){
    
    title <- paste('Histogram of ', ifelse(grepl('sinc', colnames(X)[i], fixed = TRUE),
                                           paste(gsub('sinc_', 'sinc-transformed covariate originating from ', colnames(X)[i]), sep = ''),
                                           ifelse(grepl('exp', colnames(X)[i], fixed = TRUE),
                                                  paste(gsub('exp_', 'exponentially-transformed covariate originating from ', colnames(X)[i]), sep = ''),
                                                  ifelse(grepl('X', colnames(X)[i], fixed = TRUE),
                                                         paste('covariate originating from cross-product of ', unlist(strsplit(colnames(X)[i], 'X', fixed = TRUE))[1], ' and ', unlist(strsplit(colnames(X)[i], 'X', fixed = TRUE))[2], sep = ''),
                                                         paste('covariate originating from ', colnames(X)[i], sep = '')))),
                   ' distribution', sep = '')
    
    plot <- ggplot(X, aes(x = as.factor(X[,i]), fill = as.factor(treateds[[5]]), color = as.factor(treateds[[5]]))) +
      geom_bar(alpha = 0.8) +
      own_theme  +
      labs(title = title,
           x = 'Value', y = 'Count') +
      scale_fill_manual(name = ' ', labels = c('Not treated', 'Treated'),
                        values = c(nt_color, t_color)) +
      scale_color_manual(name = ' ', labels = c('Not treated', 'Treated'),
                         values = c(nt_color, t_color)) +
      theme(legend.position = 'bottom')
    
    png(file = paste(outputfolder, '/histogram_', colnames(X)[i], '.png', sep = ''),
        width = 600, height = 600)
    print(plot)
    dev.off()
    
  } else {
    
    title <- paste('Histogram of ', ifelse(grepl('sinc', colnames(X)[i], fixed = TRUE),
                                           paste(gsub('sinc_', 'sinc-transformed covariate originating from ', colnames(X)[i]), sep = ''),
                                           ifelse(grepl('exp', colnames(X)[i], fixed = TRUE),
                                                  paste(gsub('exp_', 'exponentially-transformed covariate originating from ', colnames(X)[i]), sep = ''),
                                                  ifelse(grepl('X', colnames(X)[i], fixed = TRUE),
                                                         paste('covariate originating from cross-product of ', unlist(strsplit(colnames(X)[i], 'X', fixed = TRUE))[1], ' and ', unlist(strsplit(colnames(X)[i], 'X', fixed = TRUE))[2], sep = ''),
                                                         paste('covariate originating from ', colnames(X)[i], sep = '')))),
                   ' distribution', sep = '')
    
    plot <- ggplot(X, aes(x = X[,i], fill = as.factor(treateds[[5]]), color = as.factor(treateds[[5]]))) +
      geom_histogram(bins = min(length(unique(X[,i])), 1000),
                     alpha = 0.8) +
      geom_vline(aes(xintercept = mean(X[which(treateds[[1]] == 0),i]), linetype = 'dashed'),
                 color = nt_color, alpha = .8, size = 1) +
      geom_vline(aes(xintercept = mean(X[which(treateds[[1]] == 1),i]), linetype = 'dashed'),
                 color = t_color, alpha = .8, size = 1) +
      geom_vline(aes(xintercept = median(X[which(treateds[[1]] == 0),i]), linetype = 'dotted'),
                 color = nt_color, alpha = .8, size = 1) +
      geom_vline(aes(xintercept = median(X[which(treateds[[1]] == 1),i]), linetype = 'dotted'),
                 color = t_color, alpha = .8, size = 1) +
      own_theme  +
      xlim(min(X[,i]) - sd(X[,i])/3, max(X[,i]) + sd(X[,i])/3) +
      labs(title = title,
           x = 'Value', y = 'Count') +
      scale_fill_manual(name = ' ', labels = c('Not treated', 'Treated'),
                        values = c(nt_color, t_color)) +
      scale_color_manual(name = ' ', labels = c('Not treated', 'Treated'),
                         values = c(nt_color, t_color)) +
      scale_linetype_discrete(name = ' ', labels = c('Mean', 'Median')) +
      theme(legend.position = 'bottom', legend.box = 'vertical')
    
    png(file = paste(outputfolder, '/histogram_', colnames(X)[i], '.png', sep = ''),
        width = 1200, height = 600)
    print(plot)
    dev.off()
    
  } 
  
}

#Create histograms of Y variables
for(i in 1:length(Ys)){
  
  title <- paste('Histogram of outcome variable of DGP ', i, sep = '')
  
  Y <- as.data.frame(Ys[[i]])
  
  plot <- ggplot(Y, aes(x = V1, fill = as.factor(treateds[[i]]), color = as.factor(treateds[[i]]))) +
    geom_histogram(bins = min(length(unique(Y$V1)), 1000),
                   alpha = 0.8) +
    geom_vline(aes(xintercept = mean(Y[which(treateds[[i]] == 0),1]), linetype = 'dashed'),
               color = nt_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = mean(Y[which(treateds[[i]] == 1),1]), linetype = 'dashed'),
               color = t_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = median(Y[which(treateds[[i]] == 0),1]), linetype = 'dotted'),
               color = nt_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = median(Y[which(treateds[[i]] == 1),1]), linetype = 'dotted'),
               color = t_color, alpha = .8, size = 1) +
    own_theme  +
    xlim(min(Y$V1) - sd(Y$V1)/3, max(Y$V1) + sd(Y$V1)/3) +
    labs(title = title,
         x = 'Value', y = 'Count') +
    scale_fill_manual(name = ' ', labels = c('Not treated', 'Treated'),
                      values = c(nt_color, t_color)) +
    scale_color_manual(name = ' ', labels = c('Not treated', 'Treated'),
                       values = c(nt_color, t_color)) +
    scale_linetype_discrete(name = ' ', labels = c('Mean', 'Median')) +
    theme(legend.position = 'bottom', legend.box = 'vertical')
  
  png(file = paste(outputfolder, '/histogram_Y', i, '.png', sep = ''),
      width = 1200, height = 600)
  print(plot)
  dev.off()
  
}

#Create histograms of propensity scores
for(i in 1:length(prop_scores)){
  
  title <- ifelse(i == 1, paste('Histogram of propensity scores of DGPs 1, 2 and 4', sep = ''),
                  ifelse(i == 2, paste('Histogram of propensity scores of DGP 3', sep = ''),
                         paste('Histogram of propensity scores of DGP 5', sep = '')))
  
  p <- as.data.frame(prop_scores[[i]])
  
  colnames(p) <- 'V1'
  
  if(i == 1){
    
    t <- treateds[[1]]
    
  }else if(i == 2){
    
    t <- treateds[[3]]
    
  } else {
    
    t <- treateds[[5]]
    
  }
  
  plot <- ggplot(p, aes(x = V1, fill = as.factor(t), color = as.factor(t))) +
    geom_histogram(bins = min(length(unique(p$V1)), 1000),
                   alpha = 0.8) +
    geom_vline(aes(xintercept = mean(Y[which(t == 0),1]), linetype = 'dashed'),
               color = nt_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = mean(Y[which(t == 1),1]), linetype = 'dashed'),
               color = t_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = median(Y[which(t == 0),1]), linetype = 'dotted'),
               color = nt_color, alpha = .8, size = 1) +
    geom_vline(aes(xintercept = median(Y[which(t == 1),1]), linetype = 'dotted'),
               color = t_color, alpha = .8, size = 1) +
    own_theme  +
    xlim(0, 1) +
    labs(title = title,
         x = 'Value', y = 'Count') +
    scale_fill_manual(name = ' ', labels = c('Not treated', 'Treated'),
                      values = c(nt_color, t_color)) +
    scale_color_manual(name = ' ', labels = c('Not treated', 'Treated'),
                       values = c(nt_color, t_color)) +
    scale_linetype_discrete(name = ' ', labels = c('Mean', 'Median')) +
    theme(legend.position = 'bottom', legend.box = 'vertical')
  
  png(file = paste(outputfolder, '/histogram_propscore_', i, '.png', sep = ''),
      width = 1200, height = 600)
  print(plot)
  dev.off()
  
}

bXt <- which(colnames(X) == 'binomialXstudentst')
heterogeneous_treatment <- floor(X[,bXt]/sd(X[,bXt]))

title <- 'Histogram of heterogeneous treatment value for DGPs 4 and 5'

plot <- ggplot(as.data.frame(heterogeneous_treatment), aes(x = heterogeneous_treatment, fill = as.factor(treateds[[5]]), color = as.factor(treateds[[5]]))) +
  geom_histogram(bins = min(length(unique(heterogeneous_treatment)), 1000),
                 alpha = 0.8) +
  geom_vline(aes(xintercept = mean(heterogeneous_treatment[which(treateds[[1]] == 0)]), linetype = 'dashed'),
             color = nt_color, alpha = .8, size = 1) +
  geom_vline(aes(xintercept = mean(heterogeneous_treatment[which(treateds[[1]] == 1)]), linetype = 'dashed'),
             color = t_color, alpha = .8, size = 1) +
  geom_vline(aes(xintercept = median(heterogeneous_treatment[which(treateds[[1]] == 0)]), linetype = 'dotted'),
             color = nt_color, alpha = .8, size = 1) +
  geom_vline(aes(xintercept = median(heterogeneous_treatment[which(treateds[[1]] == 1)]), linetype = 'dotted'),
             color = t_color, alpha = .8, size = 1) +
  own_theme  +
  xlim(min(heterogeneous_treatment) - sd(heterogeneous_treatment)/3, max(heterogeneous_treatment) + sd(heterogeneous_treatment)/3) +
  labs(title = title,
       x = 'Value', y = 'Count') +
  scale_fill_manual(name = ' ', labels = c('Not treated', 'Treated'),
                    values = c(nt_color, t_color)) +
  scale_color_manual(name = ' ', labels = c('Not treated', 'Treated'),
                     values = c(nt_color, t_color)) +
  scale_linetype_discrete(name = ' ', labels = c('Mean', 'Median')) +
  theme(legend.position = 'bottom', legend.box = 'vertical')

png(file = paste(outputfolder, '/histogram_heterogeneous_treatment.png', sep = ''),
    width = 1200, height = 600)
plot
dev.off()


#################### Analysis of Weight sample ####################

load(paste(getwd(), sub('.', '', weight_dir), '/output.RData', sep = ''))

#Create folder for output
outputfolder <- paste(getwd(), '/analysis/weight_sample', sep = '')
dir.create(outputfolder, recursive = TRUE)

#Combine into dfs
for(i in 1:5){
  
  out <- as.data.frame(matrix(data = NA, nrow = length(output), ncol = length(output[[1]][[1]])))
  colnames(out) <- c('EN', 'KRLS', 'RL', 'RF', 'BGLM', 'BART')
  
  for(j in 1:length(output)){
    
    it <- output[[j]]
    dgp <- it [[i]]
    
    out[j,] <- dgp
    
  }

  assign(paste('weights_dgp', i, sep = ''), out)
  
}


#Analyse weights for specific estimators
#Estimators
estimators <- c('EN', 'KRLS', 'RL', 'RF', 'BGLM', 'BART')

for(i in 1:length(estimators)){
  
  out <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 5))
  colnames(out) <- c('DGP 1', 'DGP 2', 'DGP 3', 'DGP 4', 'DGP 5')
  rownames(out) <- c('Max', 'Min', 'Mean', 'Median')
  
  for(j in 1:5){
    
    df <- get(paste('weights_dgp', j, sep = ''))
    
    df <- df[,i]
    
    out[1,j] <- max(df)
    out[2,j] <- min(df)
    out[3,j] <- mean(df)
    out[4,j] <- median(df)
    
  }
  
  print(xtable(out, type = 'latex'), file = paste(outputfolder, '/weights_', estimators[i], '.tex', sep = ''))
  
}

#Analyse weights for specific DGPs
dgp <- c('DGP 1', 'DGP 2', 'DGP 3', 'DGP 4', 'DGP 5')

for(i in 1:length(dgp)){
  
  out <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 6))
  colnames(out) <- c('EN', 'KRLS', 'RL', 'RF', 'BGLM', 'BART')
  rownames(out) <- c('Max', 'Min', 'Mean', 'Median')
  
  df <- get(paste('weights_dgp', i, sep = ''))
  out[1,] <- apply(df, MARGIN = 2, FUN = max)
  out[2,] <- apply(df, MARGIN = 2, FUN = min)
  out[3,] <- apply(df, MARGIN = 2, FUN = mean)
  out[4,] <- apply(df, MARGIN = 2, FUN = median)
  
  print(xtable(out, type = 'latex'), file = paste(outputfolder, '/weights_', dgp[i], '.tex', sep = ''))
  
}

#Analyse weight distribution per iterations
#Histogram of maximal / second biggest and median weight per dgp
for(i in 1:5){
  
  ext <- get(paste('weights_dgp', i, sep = ''))
  
  ext <- data.frame(max = apply(ext, MARGIN = 1, FUN = max),
                    secondmax = apply(ext, MARGIN = 1, FUN = function(x) max(x[-which(x == max(x))])),
                    median = apply(ext, MARGIN = 1, FUN = median),
                    min = apply(ext, MARGIN = 1, FUN = min))
  
  ext <- round(ext*10, 0)/10
  
  df <- data.frame(x = rep((0:10)/10, times = 4),
                   variable = unlist(lapply(colnames(ext), function(x) rep(x, times = 11))),
                   y = rep(0, times = 44))
  
  df$variable <- factor(df$variable, levels = c('max', 'secondmax', 'median', 'min'))

  for(j in 1:nrow(df)){
    
    df[j,'y'] <- sum(ext[,df[j,'variable']] == df[j,'x'])/19
    
  }
                  
  title <- paste('Share of iterations by weight-extrema for DGP ', i, sep = '')
  
  plot <- ggplot(df, aes(x = x, y = y, fill = variable)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    own_theme  +
    scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = c(0:10/10), labels = as.character(round(c(0:10/10), 1))) +
    scale_fill_viridis(discrete = T, name = '', labels = c('Max', 'Second max', 'Median', 'Min')) +
    labs(title = title,
         x = 'Weight', y = 'Share of iterations')
  
  for(j in 0:10){
    
    plot <- plot +
      geom_segment(aes(y = 0, yend = 0), x = j/10-4/90, xend = j/10+4/90,
                   color = '#3B3838')
    
  }
  
  png(file = paste(outputfolder, '/histogram_maxweight_', i, '.png', sep = ''),
      width = 1200, height = 600)
  print(plot)
  dev.off()
  
}