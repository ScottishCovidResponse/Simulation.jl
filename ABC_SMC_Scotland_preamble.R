run_model <- function(params){
  julia_assign("paramraw", params)
  julia_source('examples/Epidemiology/Scotland_run_inputs.jl')
  y <- julia_eval("abuns")
  #y<-y[41:48,,]
  return(apply(y,c(1,3),sum) )
}



calc_distance <- function(sind, D, D_star,ntimes){
  # Define sum of squared errors
  # Vectorised function
  xx <- dim(D)[1]
  switch(sind, sqrt(sum( (D[,ntimes] - D_star[,ntimes]) ^ 2)),
         sqrt( (sapply(1:xx,function(i)which.max(D[i,]))-sapply(xx,function(i)which.max(D_star[i,])))^2)
  )
  
}


calc_distance_final_size <- function(mu, mu_star){
  # Define sum of squared errors
  # Vectorised function
  dist <- sqrt( (mu - mu_star) ^ 2)
  return(dist)
}


summary_stat <- function(D){
  mu <- tail(D, n = 1)
  return(mu)
}



# Perturbation kernel 
rK <- function(mean, sigma){   
  return(rtmvnorm(1,mean=mean, sigma=sigma, lower=parrange[1,], upper=parrange[2,])) 
}

#  Identity function: H(x)= 1 if x=T
H <- function(x) as.numeric(x>0)

#  Test if prior is non zero
prior.non.zero<-function(par){
  prod(sapply(1:n_par, function(a) H(par[a]-parrange[1,a])* H(parrange[2,a]-par[a])))
}

Norm.Eucl.dist<-function(p1,p2){
  sqrt(sum(((p1-p2)/(parrange[2,]-parrange[1,]))^2)) }

#  Covariance based on M neighbours
getSigmaNeighbours<-function(M, theta, Theta){
  dist<- sapply(1:N, function(a) Norm.Eucl.dist(as.numeric(theta), as.numeric(Theta[a,])))
  temp<-data.frame(no=seq(1,N), dist)
  temp<-temp[order(temp$dist),]
  sigma<-cov(Theta[temp$no[1:(M+1)],])
  return(sigma)
}

# Function for parallel
runmodpar <- function(param){ #also needs [n,nsum,Dorig,epsilon,g]
  m<-0
  distances <- c()
  distanceall <- matrix(NA,ncol=nsumt,nrow=n)
  D_star<-run_model(param)
  for(j in 1:n){
    # Calculate distance for all sum stats untill any one is rejected  
    dthr <- 0
    sind <- 1
    while(dthr==0 & sind <= nsum){
      distance <- calc_distance(sind,Dorig, D_star,ntimes)
      if(any(distance > epsilon[sind,g])){ # If the distance is less than the tolerance
        # Store results
        dthr <- 1
      }else{
        distances <- c(distances,distance)
      }
      sind <- sind + 1
    }
    print(distances)
    if(dthr==0){ # If both distances are less than their tolerances
      distanceall[j,] <- distances
      m<-m+1
    }
  } 
  
  return(m)
}
