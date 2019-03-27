#############
################# FUNCTION TO PUT RANDOM EFFECTS ON SBSS VIDEOS
################# MARY RYAN
################# 3.27.2019
#############
videoRandomEffects <- function(pi.true, sd.randEff, videos =1, ep= 1e-5, delta = 1e-5, iterPrint=FALSE){
  ## pi.true: vector of true probabilities for evaluating at each SBSS value
  ## sd.randEff: standard deviation you want for random effect
  ## correlation between videos in category
  ## videos: number of different videos want to produce
  ## ep: size of step for beta1 tuning
  ## delta: how close to 1 you would like sum of pi to get
  ## iterPrint: option to print how many iterations you've done on beta1 tuning
  ##             and current sum of pi
  
  # load MASS library for mvrnorm() #
  require( MASS )
  
  # intialize vectors and lists #
  beta0 <- rep(NA, videos)
  beta1 <- logit.pi <- pi.r <- matrix(nrow=6, ncol=videos)
  results <- rep( list(rep( NA, 4 )), videos )
  
  randomEffects <- rep(NA, length(sd.randEff))
  
  # generate random effect #
  randomEffect <- rnorm(videos, 0, sd.randEff)
  
  
  for( v in seq(videos) ){
    
    # calculate new beta0 #
    beta0[v] <- log( pi.true[1]/(1-pi.true[1]) ) + randomEffect[v]
    
    # calculate initial beta1s #
    beta1[,v] <-  log( (pi.true)/(1-pi.true) ) - beta0[v]
    
    # set beta1 for SBSS=0 to be 0 #
    beta1[1,v] <- 0
    
    # calculate pis from betas #
    logit.pi[,v] <- beta0[v] + beta1[,v]
    pi.r[,v] <- exp(logit.pi[,v])/(1+exp(logit.pi[,v]))
    
    # find beta1s necessary such that sum of pi is 1 #
    iter <- 0
    while( abs( 1 - sum(pi.r[,v]) ) > delta ){
      
      # take a little off the beta1s if we sum too high #
      if(sum(pi.r[,v]) > 1 + delta){
        
        beta1[2:6,v] <- beta1[2:6,v] - ep 
        logit.pi[,v] <- beta0[v] + beta1[,v]
        pi.r[,v] <- exp(logit.pi[,v])/(1+exp(logit.pi[,v]))
        
      }
      
      # add a little to the beta1s if we sum too low #
      if(sum(pi.r[,v]) < 1 - delta){
        
        beta1[2:6, v] <- beta1[2:6, v] + ep 
        logit.pi[,v] <- beta0[v] + beta1[,v]
        pi.r[,v] <- exp(logit.pi[,v])/(1+exp(logit.pi[,v]))
        
      }
      
      # print what iteration we're on and current pi sum if setting set to TRUE #
      if(iterPrint == TRUE){
        
        if( iter %% 1000 == 0 ) print( paste0( "Video", v, ", ", iter,", ", sum(pi.r[,v]) ) )
        iter <- iter + 1
        
      }
      
      
    } # end while
    
    
  } #end video for
  
  results <- list( pi = pi.r, beta0 = beta0,
                   beta1 = beta1, randomEffect = randomEffect )
  
  return( results )
}

#############
################# THEORETICAL KAPPA VARIANCE FUNCTION
################# MARY RYAN
################# 3.26.2019
#############
kappaSummary.theoretical <- function( cont.table, N, transform=T){
  ## cont.table: contiginecy table expressed as proportions of all observations
  ## transform: whether you want to output the Fisher-transformed kappa as well
  
  # calculating kappa statistic #
  po <- sum( diag(cont.table) )
  pe <- 0
  
  for( i in seq( nrow(cont.table) ) ){
    
    pe <- sum( cont.table[i,] )*sum( cont.table[,i] ) + pe
    
    
  }
  
  kappa <- (po - pe)/(1 - pe)
  
  # theoretical kappa variance #
  secondPart <- rep( NA, nrow(cont.table) )
  thirdPartA <- matrix( ncol=nrow(cont.table), nrow=nrow(cont.table) )
  
  p.i <- rowSums( cont.table )
  pj. <- colSums( cont.table )
  
  firstPart <- 1/(N*(1-pe)^4)
  
  for( i in seq( nrow(cont.table) ) ){
    
    secondPart[i] <- cont.table[i,i]*( (1-pe) - (p.i[i] + pj.[i])*(1-po) )^2
    
    for( j in seq( nrow(cont.table) ) ){
      
      if(i!=j){
        thirdPartA[i,j] <- cont.table[i,j]*(p.i[i] + pj.[j])^2
      } else{
        thirdPartA[i,j] <- NA
      }
      
    }
    
  }
  
  thirdPartA.sum <- sum( thirdPartA, na.rm=T )
  secondPart.sum <- sum( secondPart, na.rm=T )
  
  thirdPart <- (1-po)^2*thirdPartA.sum
  fourthPart <- (po*pe -2*pe + po)^2
  
  kappa.var <- firstPart*(secondPart.sum + thirdPart - fourthPart)
  
  if(transform == T){
    
    kappa.t <- log( (1+kappa)/(1-kappa) )
    kappa.var.t <- ( 4/( (1-kappa)^2 * (1 + kappa)^2 ) ) * kappa.var
    
    results <- cbind( c(kappa, kappa.t),
                      c(kappa.var, kappa.var.t) )
    colnames( results ) <- c("Statistic Est.", "Variance Est.")
    rownames( results ) <- c("Basic Kappa", "Fisher Transformed Kappa")
    
  } else{
    
    results <- cbind(kappa, kappa.var)
    colnames( results ) <- c("Statistic Est.", "Variance Est.")
    
  }
  
  return( results )
  
}


#############
################# FUNCTION TO SIMULATE SURGEON RATING DATA ON VIDEOS
################# MARY RYAN
################# 3.27.2019
#############
kappaSurgeonSimulations <- function(kappa = 0.8, videoRepeat = 3,
                                    nsims = 100000, surgeons = 500, seed = 12345,
                                    multipleVideos = 3, sdRandEff = rep(0, 6),
                                    alert = T){
  ### Function to generate surgeon SBSS rating data for kappa simulations
  ####
  ########################################
  ####
  ## kappa: what true kappa you want to simulate; defaults to 0.8
  ## videoRepeat: how many times you want to repeat each unique video;
  # defaults to 3
  ## nsims: number of simulations you want to perform; defaults to 100,000
  ## surgeons: number of surgeons you want rating videos in each simulation'
  # defaults to 500
  ## seed: specify number for set.seed; defaults to 12345
  ## multipleVideos: how many unique videos you want to generate
  # per SBSS category; defaults to 3
  ## sdRandEff: standard deviation settings for random effects when creating
  # multiple unique videos per SBSS category; defaults to rep(0, 6)
  ## alert: plays sound when function is done running if true; defaults to true
  ####
  ########################################
  ####
  
  ## initial warnings ##
  try( if( !(kappa %in% c(0.4, 0.6, 0.8)) ) stop("Invalid kappa input: specify 0.4, 0.6, or 0.8") )
  
  ## load libraries ##
  library(beepr)
  library(MASS)
  if( alert == T ) library(beepr)
  
  ## set seed ##
  set.seed(seed)
  
  ## set SBSS ##
  SBSS <- 0:5
  
  #### kappa = 0.8 ####
  if( kappa ==  0.8 ){
    
    ### SBSS = 0 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=0, #
    # and use coefficients to model this probability #
    pi.00t <- pbeta(0.3930378, 1, 6) #0.95
    pi.01t <- pbeta(0.55, 1, 6) - pbeta(0.3930378,1, 6) 
    pi.02t <- pbeta(0.7, 1, 6) - pbeta(0.55, 1, 6)
    pi.03t <- pbeta(0.9, 1, 6) - pbeta(0.7,1, 6)
    pi.04t <- pbeta(0.99, 1, 6) - pbeta(0.9,1, 6)
    pi.05t <- 1 - pbeta(0.99, 1, 6)
    pi.0t.true <- c(pi.00t, pi.01t, pi.02t, pi.03t, pi.04t, pi.05t)
    
    # intercept for true SBSS=0 model will be log-odds of getting 0 #
    beta0.0t <- log(pi.00t/(1-pi.00t))
    beta1.0t <- log( (pi.0t.true)/(1-pi.0t.true) ) - beta0.0t
    
    ## random videos! ##
    videos.0 <- videoRandomEffects( pi.true = pi.0t.true,
                                    sd.randEff = sdRandEff[1],
                                    video = multipleVideos )$pi
    
    ### SBSS = 1 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=1 #
    # and use coefficients to model this probability #
    pi.10t <- pbeta(0.1, 2,5)
    pi.11t <- pbeta(0.54, 2, 5) - pbeta(0.1, 2,5) #0.8095287
    pi.12t <- pbeta(0.65, 2, 5) - pbeta(0.54, 2,5)
    pi.13t <- pbeta(0.8, 2, 5) - pbeta(0.65, 2,5)
    pi.14t <- pbeta(0.95, 2, 5) - pbeta(0.8, 2,5)
    pi.15t <- 1 - pbeta(0.95, 2,5)
    pi.1t.true <- c(pi.10t, pi.11t, pi.12t, pi.13t, pi.14t, pi.15t)
    
    beta0.1t <- log(pi.10t/(1-pi.10t))#-2.047897
    beta1.1t <- log( (pi.1t.true)/(1-pi.1t.true) ) - beta0.1t
    
    ## random videos! ##
    videos.1 <- videoRandomEffects( pi.true = pi.1t.true,
                                    sd.randEff = sdRandEff[2],
                                    video = multipleVideos )$pi
    
    ### SBSS = 2 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=2 #
    # and use coefficients to model this probability #
    pi.20t <- pbeta(0.2, 4, 3)
    pi.21t <- pbeta(0.26, 4, 3) - pbeta(0.2, 4, 3)
    pi.22t <- pbeta(0.725, 4, 3) - pbeta(0.26, 4, 3) #0.7460092
    pi.23t <- pbeta(0.93, 4, 3) - pbeta(0.725, 4, 3)
    pi.24t <- pbeta(0.95, 4, 3) - pbeta(0.93, 4, 3)
    pi.25t <- 1 - pbeta(0.95, 4, 3)
    pi.2t.true <- c(pi.20t, pi.21t, pi.22t, pi.23t, pi.24t, pi.25t)
    
    beta0.2t <- log(pi.20t/(1-pi.20t))
    beta1.2t <- log( (pi.2t.true)/(1-pi.2t.true) ) - beta0.2t
    
    ## random videos! ##
    videos.2 <- videoRandomEffects( pi.true = pi.2t.true,
                                    sd.randEff = sdRandEff[3],
                                    video = multipleVideos )$pi
    
    ### SBSS = 3 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=3 #
    # and use coefficients to model this probability #
    pi.30t <- pbeta(0.2, 5, 3)
    pi.31t <- pbeta(0.4, 5, 3) - pbeta(0.2, 5, 3)
    pi.32t <- pbeta(0.47, 5, 3) - pbeta(0.4, 5, 3)
    pi.33t <- pbeta(0.85, 5, 3) - pbeta(0.47, 5, 3) #0.7458894
    pi.34t <- pbeta(0.95, 5, 3) - pbeta(0.85, 5, 3)
    pi.35t <- 1 - pbeta(0.95, 5, 3)
    pi.3t.true <- c(pi.30t, pi.31t, pi.32t, pi.33t, pi.34t, pi.35t)
    
    beta0.3t <- log(pi.30t/(1-pi.30t))
    beta1.3t <- log( (pi.3t.true)/(1-pi.3t.true) ) - beta0.3t
    
    ## random videos! ##
    videos.3 <- videoRandomEffects( pi.true = pi.3t.true,
                                    sd.randEff = sdRandEff[4],
                                    video = multipleVideos )$pi
    
    ### SBSS = 4 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=4 #
    # and use coefficients to model this probability #
    pi.40t <- pbeta(0.38, 12, 3.5)
    pi.41t <- pbeta(0.5, 12, 3.5) - pbeta(0.38, 12, 3.5)
    pi.42t <- pbeta(0.57, 12, 3.5) - pbeta(0.5, 12, 3.5)
    pi.43t <- pbeta(0.62, 12, 3.5) - pbeta(0.57, 12, 3.5)
    pi.44t <- pbeta(0.89, 12, 3.5) - pbeta(0.62, 12, 3.5) #0.7949393
    pi.45t <- 1-pbeta(0.89, 12, 3.5)
    pi.4t.true <- c(pi.40t, pi.41t, pi.42t, pi.43t, pi.44t, pi.45t)
    
    beta0.4t <- log(pi.40t/(1-pi.40t))
    beta1.4t <- log( (pi.4t.true)/(1-pi.4t.true) ) - beta0.4t
    
    ## random videos! ##
    videos.4 <- videoRandomEffects( pi.true = pi.4t.true,
                                    sd.randEff = sdRandEff[5],
                                    video = multipleVideos )$pi
    
    ### SBSS = 5 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=5 #
    # and use coefficients to model this probability #
    pi.50t <- pbeta(0.1, 18, 4)
    pi.51t <- pbeta(0.2, 18, 4) - pbeta(0.1, 18, 4)
    pi.52t <- pbeta(0.3, 18, 4) - pbeta(0.2, 18, 4)
    pi.53t <- pbeta(0.4, 18, 4) - pbeta(0.3, 18, 4)
    pi.54t <- pbeta(0.65, 18, 4) - pbeta(0.4, 18, 4)
    pi.55t <- 1-pbeta(0.65, 18, 4) #0.9669147
    pi.5t.true <- c(pi.50t, pi.51t, pi.52t, pi.53t, pi.54t, pi.55t)
    
    beta0.5t <- log(pi.50t/(1-pi.50t))
    beta1.5t <- log( (pi.5t.true)/(1-pi.5t.true) ) - beta0.5t
    
    ## random videos! ##
    videos.5 <- videoRandomEffects( pi.true = pi.5t.true,
                                    sd.randEff = sdRandEff[6],
                                    video = multipleVideos )$pi
    
  }
  
  #### kappa = 0.6 ####
  if( kappa == 0.6 ){
    
    ### SBSS = 0 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=0, #
    # and use coefficients to model this probability #
    pi.00t <- pbeta(0.3187079, 1, 6) #0.9
    pi.01t <- pbeta(0.55, 1, 6) - pbeta(0.3187079,1, 6)
    pi.02t <- pbeta(0.7, 1, 6) - pbeta(0.55, 1, 6)
    pi.03t <- pbeta(0.9, 1, 6) - pbeta(0.7,1, 6)
    pi.04t <- pbeta(0.99, 1, 6) - pbeta(0.9,1, 6)
    pi.05t <- 1 - pbeta(0.99, 1, 6)
    pi.0t.true <- c(pi.00t, pi.01t, pi.02t, pi.03t, pi.04t, pi.05t)
    
    # intercept for true SBSS=0 model will be log-odds of getting 0 #
    beta0.0t <- log(pi.00t/(1-pi.00t))
    beta1.0t <- log( (pi.0t.true)/(1-pi.0t.true) ) - beta0.0t
    
    ## random videos! ##
    videos.0 <- videoRandomEffects( pi.true = pi.0t.true,
                                    sd.randEff = sdRandEff[1],
                                    video = multipleVideos )$pi
    
    ### SBSS = 1 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=1 #
    # and use coefficients to model this probability #
    pi.10t <- pbeta(0.11, 2,5)
    pi.11t <- pbeta(0.45, 2, 5) - pbeta(0.11, 2,5) #0.7019618
    pi.12t <- pbeta(0.65, 2, 5) - pbeta(0.45, 2,5)
    pi.13t <- pbeta(0.8, 2, 5) - pbeta(0.65, 2,5)
    pi.14t <- pbeta(0.95, 2, 5) - pbeta(0.8, 2,5)
    pi.15t <- 1 - pbeta(0.95, 2,5)
    pi.1t.true <- c(pi.10t, pi.11t, pi.12t, pi.13t, pi.14t, pi.15t)
    
    beta0.1t <- log(pi.10t/(1-pi.10t))
    beta1.1t <- log( (pi.1t.true)/(1-pi.1t.true) ) - beta0.1t
    
    ## random videos! ##
    videos.1 <- videoRandomEffects( pi.true = pi.1t.true,
                                    sd.randEff = sdRandEff[2],
                                    video = multipleVideos )$pi
    
    ### SBSS = 2 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=2 #
    # and use coefficients to model this probability #
    pi.20t <- pbeta(0.2, 4, 3)
    pi.21t <- pbeta(0.4, 4, 3) - pbeta(0.2, 4, 3)
    pi.22t <- pbeta(0.615, 4, 3) - pbeta(0.4, 4, 3) #0.3961996
    pi.23t <- pbeta(0.93, 4, 3) - pbeta(0.615, 4, 3)
    pi.24t <- pbeta(0.95, 4, 3) - pbeta(0.93, 4, 3)
    pi.25t <- 1 - pbeta(0.95, 4, 3)
    pi.2t.true <- c(pi.20t, pi.21t, pi.22t, pi.23t, pi.24t, pi.25t)
    
    beta0.2t <- log(pi.20t/(1-pi.20t))#-4.059792
    beta1.2t <- log( (pi.2t.true)/(1-pi.2t.true) ) - beta0.2t
    
    ## random videos! ##
    videos.2 <- videoRandomEffects( pi.true = pi.2t.true,
                                    sd.randEff = sdRandEff[3],
                                    video = multipleVideos )$pi
    
    ### SBSS = 3 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=3 #
    # and use coefficients to model this probability #
    pi.30t <- pbeta(0.2, 5, 3)
    pi.31t <- pbeta(0.4, 5, 3) - pbeta(0.2, 5, 3)
    pi.32t <- pbeta(0.6, 5, 3) - pbeta(0.4, 5, 3)
    pi.33t <- pbeta(0.78, 5, 3) - pbeta(0.6, 5, 3) #0.3960115
    pi.34t <- pbeta(0.95, 5, 3) - pbeta(0.78, 5, 3)
    pi.35t <- 1 - pbeta(0.95, 5, 3)
    pi.3t.true <- c(pi.30t, pi.31t, pi.32t, pi.33t, pi.34t, pi.35t)
    
    beta0.3t <- log(pi.30t/(1-pi.30t))#-5.361485
    beta1.3t <- log( (pi.3t.true)/(1-pi.3t.true) ) - beta0.3t
    
    ## random videos! ##
    videos.3 <- videoRandomEffects( pi.true = pi.3t.true,
                                    sd.randEff = sdRandEff[4],
                                    video = multipleVideos )$pi
    
    ### SBSS = 4 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=4 #
    # and use coefficients to model this probability #
    pi.40t <- pbeta(0.38, 12, 3.5)
    pi.41t <- pbeta(0.55, 12, 3.5) - pbeta(0.38, 12, 3.5)
    pi.42t <- pbeta(0.63, 12, 3.5) - pbeta(0.55, 12, 3.5)
    pi.43t <- pbeta(0.69, 12, 3.5) - pbeta(0.63, 12, 3.5)
    pi.44t <- pbeta(0.9, 12, 3.5) - pbeta(0.69, 12, 3.5) # 0.702895
    pi.45t <- 1-pbeta(0.9, 12, 3.5)
    pi.4t.true <- c(pi.40t, pi.41t, pi.42t, pi.43t, pi.44t, pi.45t)
    
    beta0.4t <- log(pi.40t/(1-pi.40t))#-7.338262
    beta1.4t <- log( (pi.4t.true)/(1-pi.4t.true) ) - beta0.4t
    
    ## random videos! ##
    videos.4 <- videoRandomEffects( pi.true = pi.4t.true,
                                    sd.randEff = sdRandEff[5],
                                    video = multipleVideos )$pi
    
    ### SBSS = 5 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=5 #
    # and use coefficients to model this probability #
    pi.50t <- pbeta(0.1, 18, 4)
    pi.51t <- pbeta(0.2, 18, 4) - pbeta(0.1, 18, 4)
    pi.52t <- pbeta(0.3, 18, 4) - pbeta(0.2, 18, 4)
    pi.53t <- pbeta(0.4, 18, 4) - pbeta(0.3, 18, 4)
    pi.54t <- pbeta(0.695, 18, 4) - pbeta(0.4, 18, 4)
    pi.55t <- 1-pbeta(0.695, 18, 4) #0.9216484
    pi.5t.true <- c(pi.50t, pi.51t, pi.52t, pi.53t, pi.54t, pi.55t)
    
    beta0.5t <- log(pi.50t/(1-pi.50t))#-34.55209
    beta1.5t <- log( (pi.5t.true)/(1-pi.5t.true) ) - beta0.5t
    
    ## random videos! ##
    videos.5 <- videoRandomEffects( pi.true = pi.5t.true,
                                    sd.randEff = sdRandEff[6],
                                    video = multipleVideos )$pi
    
  }
  
  #### kappa = 0.4 ####
  if( kappa == 0.4 ){
    
    ### SBSS = 0 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=0, #
    # and use coefficients to model this probability #
    pi.00t <- pbeta(0.1818112, 1, 6) #0.7
    pi.01t <- pbeta(0.55, 1, 6) - pbeta(0.1818112,1, 6)
    pi.02t <- pbeta(0.7, 1, 6) - pbeta(0.55, 1, 6)
    pi.03t <- pbeta(0.9, 1, 6) - pbeta(0.7,1, 6)
    pi.04t <- pbeta(0.99, 1, 6) - pbeta(0.9,1, 6)
    pi.05t <- 1 - pbeta(0.99, 1, 6)
    pi.0t.true <- c(pi.00t, pi.01t, pi.02t, pi.03t, pi.04t, pi.05t)
    
    # intercept for true SBSS=0 model will be log-odds of getting 0 #
    beta0.0t <- log(pi.00t/(1-pi.00t))#2.944439
    beta1.0t <- log( (pi.0t.true)/(1-pi.0t.true) ) - beta0.0t
    
    ## random videos! ##
    videos.0 <- videoRandomEffects( pi.true = pi.0t.true,
                                    sd.randEff = sdRandEff[1],
                                    video = multipleVideos )$pi
    
    ### SBSS = 1 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=1 #
    # and use coefficients to model this probability #
    pi.10t <- pbeta(0.11, 2,5)
    pi.11t <- pbeta(0.33, 2, 5) - pbeta(0.11, 2,5) #0.5077461
    pi.12t <- pbeta(0.65, 2, 5) - pbeta(0.33, 2,5)
    pi.13t <- pbeta(0.8, 2, 5) - pbeta(0.65, 2,5)
    pi.14t <- pbeta(0.95, 2, 5) - pbeta(0.8, 2,5)
    pi.15t <- 1 - pbeta(0.95, 2,5)
    pi.1t.true <- c(pi.10t, pi.11t, pi.12t, pi.13t, pi.14t, pi.15t)
    
    beta0.1t <- log(pi.10t/(1-pi.10t))#-2.047897
    beta1.1t <- log( (pi.1t.true)/(1-pi.1t.true) ) - beta0.1t
    
    ## random videos! ##
    videos.1 <- videoRandomEffects( pi.true = pi.1t.true,
                                    sd.randEff = sdRandEff[2],
                                    video = multipleVideos )$pi
    
    ### SBSS = 2 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=2 #
    # and use coefficients to model this probability #
    pi.20t <- pbeta(0.2, 4, 3)
    pi.21t <- pbeta(0.47, 4, 3) - pbeta(0.2, 4, 3)
    pi.22t <- pbeta(0.62,4, 3) - pbeta(0.47, 4, 3) #0.2964177
    pi.23t <- pbeta(0.93, 4, 3) - pbeta(0.62, 4, 3)
    pi.24t <- pbeta(0.95, 4, 3) - pbeta(0.93, 4, 3)
    pi.25t <- 1 - pbeta(0.95, 4, 3)
    pi.2t.true <- c(pi.20t, pi.21t, pi.22t, pi.23t, pi.24t, pi.25t)
    
    beta0.2t <- log(pi.20t/(1-pi.20t))#-4.059792
    beta1.2t <- log( (pi.2t.true)/(1-pi.2t.true) ) - beta0.2t
    
    ## random videos! ##
    videos.2 <- videoRandomEffects( pi.true = pi.2t.true,
                                    sd.randEff = sdRandEff[3],
                                    video = multipleVideos )$pi
    
    ### SBSS = 3 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=3 #
    # and use coefficients to model this probability #
    pi.30t <- pbeta(0.2, 5, 3)
    pi.31t <- pbeta(0.4, 5, 3) - pbeta(0.2, 5, 3)
    pi.32t <- pbeta(0.625, 5, 3) - pbeta(0.4, 5, 3)
    pi.33t <- pbeta(0.76, 5, 3) - pbeta(0.625, 5, 3) #0.3015379
    pi.34t <- pbeta(0.95, 5, 3) - pbeta(0.76, 5, 3)
    pi.35t <- 1 - pbeta(0.95, 5, 3)
    pi.3t.true <- c(pi.30t, pi.31t, pi.32t, pi.33t, pi.34t, pi.35t)
    
    beta0.3t <- log(pi.30t/(1-pi.30t))#-5.361485
    beta1.3t <- log( (pi.3t.true)/(1-pi.3t.true) ) - beta0.3t
    
    ## random videos! ##
    videos.3 <- videoRandomEffects( pi.true = pi.3t.true,
                                    sd.randEff = sdRandEff[4],
                                    video = multipleVideos )$pi
    
    ### SBSS = 4 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=4 #
    # and use coefficients to model this probability #
    pi.40t <- pbeta(0.38, 12, 3.5)
    pi.41t <- pbeta(0.55, 12, 3.5) - pbeta(0.38, 12, 3.5)
    pi.42t <- pbeta(0.63, 12, 3.5) - pbeta(0.55, 12, 3.5)
    pi.43t <- pbeta(0.73, 12, 3.5) - pbeta(0.63, 12, 3.5)
    pi.44t <- pbeta(0.87, 12, 3.5) - pbeta(0.73, 12, 3.5) # 0.5090475
    pi.45t <- 1-pbeta(0.87, 12, 3.5)
    pi.4t.true <- c(pi.40t, pi.41t, pi.42t, pi.43t, pi.44t, pi.45t)
    
    beta0.4t <- log(pi.40t/(1-pi.40t))#-7.338262
    beta1.4t <- log( (pi.4t.true)/(1-pi.4t.true) ) - beta0.4t
    
    ## random videos! ##
    videos.4 <- videoRandomEffects( pi.true = pi.4t.true,
                                    sd.randEff = sdRandEff[5],
                                    video = multipleVideos )$pi
    
    ### SBSS = 5 ###
    ## TRUE VIDEO ##
    # put beta distribution on pi for true SBSS=5 #
    # and use coefficients to model this probability #
    pi.50t <- pbeta(0.1, 18, 4)
    pi.51t <- pbeta(0.2, 18, 4) - pbeta(0.1, 18, 4)
    pi.52t <- pbeta(0.3, 18, 4) - pbeta(0.2, 18, 4)
    pi.53t <- pbeta(0.4, 18, 4) - pbeta(0.3, 18, 4)
    pi.54t <- pbeta(0.784, 18, 4) - pbeta(0.4, 18, 4)
    pi.55t <- 1-pbeta(0.784, 18, 4) #0.695019
    pi.5t.true <- c(pi.50t, pi.51t, pi.52t, pi.53t, pi.54t, pi.55t)
    
    beta0.5t <- log(pi.50t/(1-pi.50t))#-34.55209
    beta1.5t <- log( (pi.5t.true)/(1-pi.5t.true) ) - beta0.5t
    
    ## random videos! ##
    videos.5 <- videoRandomEffects( pi.true = pi.5t.true,
                                    sd.randEff = sdRandEff[6],
                                    video = multipleVideos )$pi
    
  }
  
  
  #### find true kappa ####
  # progress message #
  print("Calculating true kappa...")
  
  screens <-videoRepeat*length(SBSS)
  rater.cont.true <- matrix( cbind( c(rowSums(videos.0)),
                                    c(rowSums(videos.1)),
                                    c(rowSums(videos.2)),
                                    c(rowSums(videos.3)),
                                    c(rowSums(videos.4)),
                                    c(rowSums(videos.5)) ),
                             ncol=6 )
  
  rater.cont.true <- rater.cont.true/screens
  
  po.true <- sum(diag(rater.cont.true))
  pe.true <- 0
  for(i in seq(length(SBSS))){
    
    pe.true <- sum(rater.cont.true[i,])*sum(rater.cont.true[,i]) + pe.true
    
    
  }
  
  kappa.true <- (po.true - pe.true)/(1 - pe.true)
  
  #### kappa simulations ####
  # progress message #
  print("Beginning simulations...")
  
  # generate contingency table #
  kappa <- kappa.var <- kappa.var.t <- rep(NA, nsims)
  
  for( j in seq(nsims) ){
    cont.table <- rep( list(rep( NA, 1 )), 3 )
    
    for( v in seq(multipleVideos) ){
      
      cont.table[[v]] <- matrix( cbind(matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.0[,v])), ncol=1),
                                       matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.1[,v])), ncol=1),
                                       matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.2[,v])), ncol=1),
                                       matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.3[,v])), ncol=1),
                                       matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.4[,v])), ncol=1),
                                       matrix(rowSums(rmultinom(surgeons, videoRepeat, videos.5[,v])), ncol=1)),
                                 ncol=6 )
      
    }# end multiple video for loop
    
    # combine contingency tables into one main table #
    cont.table <- Reduce( '+', cont.table )
    N <- sum( rowSums(cont.table) )
    cont.table <- cont.table/N
    
    # calculate kappa statistic and variance #
    kappa.summary <- kappaSummary.theoretical( cont.table, N )
    kappa[j] <- kappa.summary[1,1]
    if( kappa[j]==1 ) kappa[j]<- NA
    
    kappa.var[j] <- kappa.summary[1,2]
    kappa.var.t[j] <- kappa.summary[2,2]
    
    # progress message #
    if( j %in% c(nsims*seq(0, 1, by=0.25)) ){
      print( paste0( "Simulation progress: ",(j/nsims)*100, "%" ) )
    }
    
  }# end sim for loop 
  
  # progress message #
  print("Calculating results...")
  
  # transform true kappa and estimated kappa #
  kappa.true.t <- log( (1+kappa.true)/(1-kappa.true) )
  kappa.t <- log( (1+kappa)/(1-kappa) )
  
  # empirical variance #
  empVar.kappa.t <- var(kappa.t, na.rm=T)
  empVar.kappa <- var(kappa, na.rm=T)
  
  # check coverage #
  CI95.lo.t <- kappa.t- 1.96 * sqrt( kappa.var.t )
  CI95.up.t <- kappa.t+ 1.96 * sqrt( kappa.var.t )
  coverage.t <- length( which( kappa.true.t >= CI95.lo.t & kappa.true.t <= CI95.up.t ) )/(nsims-sum(is.na(kappa)))
  
  CI95.lo <- kappa- 1.96 * sqrt( kappa.var )
  CI95.up <- kappa+ 1.96 * sqrt( kappa.var )
  coverage <- length( which( kappa.true >= CI95.lo & kappa.true <= CI95.up ) )/(nsims-sum(is.na(kappa)))
  
  #### output results ####
  results <- cbind(c(kappa.true, kappa.true.t),
                   c(mean(kappa, na.rm=T), mean(kappa.t, na.rm=T)),
                   c(mean(kappa.var, na.rm=T), mean(kappa.var.t, na.rm=T)),
                   c(empVar.kappa, empVar.kappa.t),
                   c(coverage, coverage.t)
  )
  
  colnames( results ) <- c("True Kappa", "Mean Statistic Est.",
                           "Mean Theoretical Var.", "Empirical Var.",
                           "Coverage Prob.")
  rownames( results ) <- c("Basic Kappa", "Fisher Transformed Kappa")
  
  if( alert == T ) beep('coin')
  return( results )
  
}
