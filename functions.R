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

