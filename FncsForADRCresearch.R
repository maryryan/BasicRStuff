#############
################# FUNCTION TO OMIT NAS IN SPECIFIC COLUMN
################# MARY RYAN
################# 8.29.2018
#############

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


#############
################# CONFIDENCE INTERVALS FOR GLM MODELS
################# MARY RYAN
################# 8.29.2018
#############

glmCI <- function( model, transform=TRUE, robust=FALSE ){
  link <- model$family$link
  coef <- summary( model )$coef[,1]
  se <- ifelse1( robust, robust.se.glm(model)[,2], summary( model )$coef[,2] )
  zvalue <- coef / se
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  
  if( transform & is.element(link, c("logit","log")) ){
    ci95.lo.exp <- exp( coef - qnorm(.975) * se )
    ci95.hi.exp <- exp( coef + qnorm(.975) * se )
    est.exp <- exp( coef )
  }
  
  ci95.lo <- as.numeric(coef - qnorm(.975) * se)
  ci95.hi <- as.numeric(coef + qnorm(.975) * se)
  est <- coef
  
  if(transform & is.element(link, c("logit","log"))){
    rslt <- cbind( round(est,4), paste0("(", round(ci95.lo,4), ",", round(ci95.hi,4), ")"),
                   round(est.exp, 4), 
                   paste0("(", round(ci95.lo.exp,4), ",", round(ci95.hi.exp,4), ")"))
    colnames(rslt) <- c("Est", "CI95", "exp( Est )", "exp( CI95 )")
    
  } else{
    rslt <- cbind( round(est,4), paste0("(", round(ci95.lo,4), ",", round(ci95.hi,4), ")"),
                   round(zvalue, 4), round(pvalue,4))
  }
  if(robust==T){
    colnames(rslt) <- c("Est", "robust ci95.lo", "robust ci95.hi", "robust z value", "robust Pr(>|z|)")
  }
  rslt
}


#############
################# CONFIDENCE INTERVALS FOR GEE MODELS
################# MARY RYAN
################# 8.29.2018
#############

geeCI <- function(model, robust=F, transform=F){
  est <- summary(model)$coef[,1]
  
  if(robust==F){
    lower <- summary(model)$coef[,1] - 1.96*summary(model)$coef[,2]
    upper <- summary(model)$coef[,1] + 1.96*summary(model)$coef[,2]
  } else{
    lower <- summary(model)$coef[,1] - 1.96*summary(model)$coef[,4]
    upper <- summary(model)$coef[,1] + 1.96*summary(model)$coef[,4]
  }
  if(transform==T){
    lower.exp <- exp(lower)
    upper.exp <- exp(upper)
    est.exp <- exp(est)
    cbind(names(summary(model)$coef[,1]),
          round(est,3),
          paste0("(", round(lower,3), ",", round(upper,3), ")"),
          round(est.exp, 3),
          paste0("(", round(lower.exp,3), ",", round(upper.exp,3), ")"))
  }else{
    
    cbind(names(summary(model)$coef[,1]),
          round(est,3),
          paste0("(", round(lower,3), ",", round(upper,3), ")") )
  }
}


#############
################# CONFIDENCE INTERVALS FOR MULTINOMIAL RESPONSE GEE MODELS
################# MARY RYAN
################# 8.29.2018
#############

geeMultinomCI <- function(model, robust=F, transform=T){
  
  if(robust==F){
    est <- model$coef
    lower <- model$coef - 1.96*summary(model)$coef[,2]
    upper <- model$coef + 1.96*summary(model)$coef[,2]
  } else{
    est <- model$coef
    lower <- model$coef - 1.96*sqrt(diag(model$robust.variance))
    upper <- model$coef + 1.96*sqrt(diag(model$robust.variance))
  }
  if(transform==T){
    lower <- exp(lower)
    upper <- exp(upper)
    est <- exp(est)
  }
  
  cbind(names(summary(model)$coef[,1]),
        round(est,3),
        paste0("(", round(lower,3), ",", round(upper,3), ")") )
}


#############
################# CONFIDENCE INTERVALS FOR LME MODELS
################# MARY RYAN
################# 8.29.2018
#############

lmeCI <- function(model, transform=F){
  est <- summary(model)$coef$fixed
  sd <- sqrt(diag(summary(model)$varFix))
  
  lower <- est - 1.96*sd
  upper <- est + 1.96*sd
  
  if(transform==T){
    lower.exp <- exp(lower)
    upper.exp <- exp(upper)
    est.exp <- exp(est)
    cbind(names(summary(model)$coef[,1]),
          round(est,3),
          paste0("(", round(lower,3), ",", round(upper,3), ")"),
          round(est.exp, 3),
          paste0("(", round(lower.exp,3), ",", round(upper.exp,3), ")"))
  }else{
    
    cbind(names(summary(model)$coef$fixed),
          round(est,3),
          paste0("(", round(lower,3), ",", round(upper,3), ")") )
  }
}



#############
################# CONFIDENCE INTERVALS FOR SURVIVAL MODELS
################# MARY RYAN
################# 8.29.2018
#############

survivalCI <- function(model){
  est <- summary(model)$coef[,1]
  se <- summary(model)$coef[,3]
  ci95.lo <- round(est - qnorm(.975) * se, 4)
  ci95.hi <- round(est + qnorm(.975) * se, 4)
  ci95 <- paste0("(", ci95.lo, ", ", ci95.hi,")")
  est.exp <- exp(est)
  ci95.lo.exp <- round(exp(est - qnorm(.975) * se), 4)
  ci95.hi.exp <- round(exp(est + qnorm(.975) * se), 4)
  ci95.exp <- paste0("(", ci95.lo.exp, ", ", ci95.hi.exp,")")
  rslt <- cbind(round(est,4), ci95, round(est.exp,4), ci95.exp)
  
  colnames(rslt) <- c("Est", "CI95", "exp( Est )", "exp( CI95 )")
  
  rslt
  
}




#############
################# K-FOLD CROSS VALIDATION FUNCTION
################# FOR RANDOM FOREST MSE
################# MARY RYAN
################# 8.29.2018
#############

## example data ##
# data(iris)
# bloop <- cbind(rep(seq(30),each=5),iris)
# colnames(bloop)<-c("ID",colnames(iris))
# B<- 100
# n_folds <- 4

rf.cv.kfold <- function(data, response, n_folds, B, cluster=NULL){
  ### data: dataframe containing all candidate predictors, response, and cluster IDs
  ### response: string of column name of data containing response
  ### n_folds: set how many folds in k-fold you want 
  ### B: how many bootstrap iterations you want
  ### cluster: optional argument if you want to sample on cluster IDs instead of row indices
    # strong of the column name of data containing cluster IDs
  
  require(randomForest)
  
  
   if( is.null(cluster)==TRUE ){
      # if you want to sample on row indices, set n_train to number of rows,
      # otherwise set n_train to number of unique clusters
      n_train <- nrow(data)
      # initialize cv storage matrix #
      cv_tmp <- matrix(NA, nrow = n_folds, ncol = B)

      # start the bootstrap #
      for(b in 1:B){
        # randomly order k-fold group numbers #
         folds_i <- sample( rep(1:n_folds, length.out=n_train) ) 
         
         # start the k-fold CV #
         for(k in 1:n_folds){
           # say the kth group will be the test set #
           test_i <- which(folds_i == k)
           # create a subset of the training dataset with all the observations not in the test IDs #
           train_data <-data[-test_i,]
           train_response <- data[-test_i,paste(response)]
           
           # create subset of the testing dataset with only the observations in the test IDs #
           test_data <- data[test_i,]
           test_response <- data[test_i,paste(response)]
            
           # run random forests on training set #
           fittedModel <- randomForest( train_response~.,train_data )
           
           # use model above to predict values from test set #
           pred <- predict( fittedModel, test_data )
           # find the MSE for this k-fold #
           cv_tmp[k,b] <- sapply(as.list(data.frame(pred)), function(pred) mean((test_response - pred)^2))
            
            
         }# repeat until each k-fold group has an MSE estimate #
         
      }# repeat process B times #
      # average the MSEs of each k-fold run together #
      cv.B <- colMeans(cv_tmp)
      # average the means of each bootstrap run together #
      cv <- mean(cv.B)
      
   } else{
     # if you want to sample on clusters, set n_train to unique number of clusters,
     # otherwise set n_train to number of rows in dataset
      n_train <- length( unique( data[,paste(cluster)] ) )
      # initialize cv storage matrix #
      cv_tmp <- matrix(NA, nrow = n_folds, ncol = B)
      
      # start the bootstrap #
      for(b in 1:B){
        # randomly order k-fold group numbers #
         folds_i <- sample( rep(1:n_folds, length.out=n_train) )
         # assign k-fold group numbers to each cluster ID (or row index) #
         cluster_folds <- cbind( unique( data[,paste(cluster)] ), folds_i )
         
         # start the k-fold CV #
         for(k in 1:n_folds){

           # say the kth group will be the test set #
           test_cluster <- data[which(cluster_folds[,2] == k),1]
           # create a subset of the training dataset with all the observations not in the test IDs #
           train_data <-data[which(!(data[,paste(cluster)] %in% test_cluster)),]
           train_response <- data[which(!(data[,paste(cluster)] %in% test_cluster)),paste(response)]
           
           # create subset of the testing dataset with only the observations in the test IDs #
           test_data <- data[which(data[,paste(cluster)] %in% test_cluster),]
           test_response <- data[which(data[,paste(cluster)] %in% test_cluster),paste(response)]
           
           # run random forests on training set #
           fittedModel <- randomForest( train_response~., train_data )
           
           # use model above to predict values from test set #
           pred <- predict( fittedModel, test_data )
           # find the MSE for this k-fold #
           cv_tmp[k,b] <- sapply(as.list(data.frame(pred)), function(pred) mean((test_response - pred)^2))
           

         }# repeat until each k-fold group has an MSE estimate #

      }# repeat process B times #
      # average the MSEs of each k-fold run together #
      cv.B <- colMeans(cv_tmp)
      # average the means of each bootstrap run together #
      cv <- mean(cv.B)
      
   }

  return(cv)
}
