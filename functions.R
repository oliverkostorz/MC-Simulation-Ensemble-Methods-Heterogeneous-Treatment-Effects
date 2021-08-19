#pacman::p_load()

#######################################################
################## Custom Functions ###################
#######################################################

#Sinc function that returns appropriate column names
custom_sinc <- function(df){
  
  colnames <- colnames(df)
  
  colnames_sinc <- paste('sinc_', colnames, sep = '')
  
  df_sinc <- sinc(df)
  
  colnames(df_sinc) <- colnames_sinc
  
  return(df_sinc)
  
}

#Exponential function that returns appropriate column names
custom_exp <- function(df){
  
  colnames <- colnames(df)
  
  colnames_exp <- paste('exp_', colnames, sep = '')
  
  df_exp <- exp(df)
  
  colnames(df_exp) <- colnames_exp
  
  return(df_exp)
  
}

#Cross product that returns appropriate column names and excludes squared terms of binary variables
custom_crossprod <- function(df){
  
  max <- ncol(df)
  c <- 1
  ncol = sum(1:max)
  
  cross_df <- data.frame(matrix(data = NA, nrow = nrow(df), ncol = ncol))
  colnames <- list()
  
  for(i in 1:(max)){
    
    for(j in 1:(max)){
      
      if(i < j){
        
        cross_df[,c] <- df[,i] * df[,j]
        c <- c + 1
        colnames <- append(colnames, paste(colnames(df)[i], colnames(df)[j], sep = 'X'))
        
      }else if(i == j){
        
        if(length(unique(df[,i])) != 2){
          
          cross_df[,c] <- df[,i] * df[,j]
          c <- c + 1
          colnames <- append(colnames, paste(colnames(df)[i], colnames(df)[j], sep = 'X'))
          
        }
      
      }
      
    }
    
  }
  
  cross_df <-  cross_df[,which(colSums(is.na(cross_df))!=nrow(cross_df))]
  
  colnames(cross_df) <- colnames
  
  return(cross_df)
  
}


#Function to send text to socket
Log <- function(text, ...) {
  
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
  
}


#Rank per iteration
c.rank <- function(df){
  
  sq_error <- (df[,1:(ncol(df) - 1)] - df$`True Value`)^2
  out <- as.data.frame(matrix(data = NA, ncol = ncol(sq_error), nrow = nrow(sq_error)/samplesize))
  colnames(out) <- colnames(sq_error)
  
  for(i in 1:nrow(out)){
    
    out[i,] <- rank(colMeans(sq_error[((samplesize*(i-1)+1):(samplesize*i)),]),
                    ties.method = 'min')
    
  }
  
  return(out)
  
}


#Rank table
ranktable <- function(df){
  
  out <- as.data.frame(matrix(data = NA, ncol = ncol(df), nrow = ncol(df)))
  colnames(out) <- colnames(df)
  
  for(i in 1:nrow(out)){
    
    out[i,] <- apply(df, MARGIN = 2, FUN = function(x) sum(x == i))
    
  }
  
  return(out)
  
}


#Rankcorrelation
c.rankcor <- function(df){
  
  out <- as.data.frame(matrix(data = NA, nrow = ncol(df), ncol = ncol(df)))
  colnames(out) <- colnames(df)
  rownames(out) <- colnames(df)
  for(i in 1:ncol(out)){
    
    for(j in 1:nrow(out)){
      
      if(i > j){
        
        
      } else {
        
        out[j,i] <- cor(df[,j], df[,i], method = 'spearman')
        
      }
      
    }
    
  }
  
  return(out)
  
}

#Change IATEs to GATEs
c.gate_internal <- function(v, truevalue){
  
  gates <- c()
  
  for(j in unique(truevalue)){
    
    gates <- c(gates, mean(v[truevalue == j]))
    
  }
  
  out <- cbind(unique(truevalue), gates)
  
  return(out)
  
}

c.gate <- function(df){
  
  #Split into iteration
  for(i in 1:(nrow(df)/samplesize)){
    
    df_it <- df[((i-1)*samplesize+1):(i*samplesize),]
    truevalue <- df_it[,ncol(df_it)]
    
    #Replace IATEs wih GATEs
    for(j in 1:(ncol(df_it)-1)){
      
      gates <- as.data.frame(c.gate_internal(df_it[,j], truevalue))
      
      for(k in 1:nrow(gates)){
        
        df_it[which(truevalue == gates$V1[k]),j] <- gates$gates[k]
        
      }
      
    }
    
    df[((i-1)*samplesize+1):(i*samplesize),] <- df_it
    
  }
  
  return(df)
  
}


########## Performance ##########
perf <- function(df){
  
  true_value <- df$`True Value`
  df <- df[,1:(ncol(df) - 1)]
  error <- df - true_value
  absolute_bias <- colMeans(abs(error))
  max_abs_bias <- apply(error, MARGIN = 2, FUN = function(x) max(abs(x)))
  mean_bias <- colMeans(error)
  
  skewness <- skewness(error[1:samplesize,])
  kurtosis <- kurtosis(error[1:samplesize,])
  variance <- apply(error[1:samplesize,], MARGIN = 2, FUN = function(x) var(x))
  for(i in 2:(nrow(error)/samplesize)){
    
    sample <- error[((i-1)*samplesize+1):(samplesize*i),]
    
    skewness <- rbind(skewness, skewness(sample))
    kurtosis <- rbind(kurtosis, kurtosis(sample))
    variance <- rbind(variance, apply(sample, MARGIN = 2, FUN = function(x) var(x)))
    
  }
  mean_skewness <- colMeans(skewness)
  mean_kurtosis <- colMeans(kurtosis)
  mean_variance <- colMeans(variance)
  
  mse <- colMeans(error^2)
  standard_error <- apply(error, MARGIN = 2, FUN = function(x) sd(x))/sqrt(nrow(error))
  
  performance_measures <- rbind(mse, mean_variance, standard_error,
                                absolute_bias, mean_bias, max_abs_bias,
                                mean_skewness, mean_kurtosis)
  
  return(performance_measures)
  
}

#Boxplots of errors
c.boxplot <- function(df, dgp_num){
  
  error <- df[,1:(ncol(df) - 1)] - df$`True Value`
  error <- melt(error)
  
  plot <- ggplot(error, aes(x = variable, y = value)) + 
    stat_boxplot(geom = 'errorbar', width = .4) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = .8, color = 'darkgrey') +
    own_theme +
    scale_x_discrete() +
    #scale_y_continuous(trans = ssqrt_trans) +
    #scale_y_continuous(trans = pseudolog10_trans) +
    labs(title = paste('Biases of effect estimates of DGP ', dgp_num, sep = ''),
         subtitle = paste('Sample size: ', samplesize, sep = ''),
         y = 'Bias', x = '')
  
  png(file = paste(outputfolder, '/boxplot_dgp', dgp_num, '.png', sep = ''),
      width = 600, height = 600)
  print(plot)
  dev.off()
  
}

c.sqrt_boxplot <- function(df, dgp_num){
  
  error <- df[,1:(ncol(df) - 1)] - df$`True Value`
  error <- melt(error)
  
  plot <- ggplot(error, aes(x = variable, y = value)) + 
    stat_boxplot(geom = 'errorbar', width = .4) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = .8, color = 'darkgrey') +
    own_theme +
    scale_x_discrete() +
    scale_y_continuous(trans = ssqrt_trans) +
    #scale_y_continuous(trans = pseudolog10_trans) +
    labs(title = paste('Biases of effect estimates of DGP ', dgp_num, sep = ''),
         subtitle = paste('Sample size: ', samplesize, sep = ''),
         caption = 'Square-root-transformed scale',
         y = 'Bias', x = '')
  
  png(file = paste(outputfolder, '/sqrt_boxplot_dgp', dgp_num, '.png', sep = ''),
      width = 600, height = 600)
  print(plot)
  dev.off()
  
}

c.log_boxplot <- function(df, dgp_num){
  
  error <- df[,1:(ncol(df) - 1)] - df$`True Value`
  error <- melt(error)
  
  plot <- ggplot(error, aes(x = variable, y = value)) + 
    stat_boxplot(geom = 'errorbar', width = .4) +
    geom_boxplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = .8, color = 'darkgrey') +
    own_theme +
    scale_x_discrete() +
    #scale_y_continuous(trans = ssqrt_trans) +
    scale_y_continuous(trans = pseudolog10_trans) +
    labs(title = paste('Biases of effect estimates of DGP ', dgp_num, sep = ''),
         subtitle = paste('Sample size: ', samplesize, sep = ''),
         caption = 'logarithm-transformed scale',
         y = 'Bias', x = '')
  
  png(file = paste(outputfolder, '/log_boxplot_dgp', dgp_num, '.png', sep = ''),
      width = 600, height = 600)
  print(plot)
  dev.off()
  
}