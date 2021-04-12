rm(list = ls())

crossprod <- function(df){
  
  max <- ncol(df)
  c <- 1
  ncol = 9999
  
  cross_df <- data.frame(matrix(data = NA, nrow = nrow(df), ncol = ncol))
  
  for(i in 1:(max)){
    
    for(j in 1:(max)){
      
      if(i>j){
        
      }else{
        
        cross_df[,c] <- df[,i] * df[,j]
        c <- c + 1
        
      }
      
    }
    
  }
  
  cross_df <-  cross_df[,which(colSums(is.na(cross_df))!=nrow(cross_df))]
  
  return(cross_df)
  
}



X <- data.frame(matrix(data = rep(1:10, times = 21), ncol = 21))

covariates <- ncol(X)

X <- cbind(X, X[,3:ncol(X)]*10, X[,3:21]*10)

cov_trans <- ncol(X)

X <- cbind(X, X[,3:ncol(X)]^2)

cov_trans_squared <- ncol(X)

X <- cbind(X, crossprod(X[,1:cov_trans]))

cov_trans_squared_cross <- ncol(X)-2


p <- 0
for(i in 1:59){
  
  p <- p + i
  
}

p + cov_trans_squared
  