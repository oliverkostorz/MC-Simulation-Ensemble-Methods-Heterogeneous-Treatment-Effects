pacman::p_load()

rm(list = ls())

crossprod <- function(df){
  
  max <- ncol(df)
  c <- 1
  ncol = sum(1:max)
  
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

