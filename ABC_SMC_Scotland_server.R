

cdargs <- commandArgs(trailingOnly = TRUE)

if(!require(tmvtnorm)){install.packages("tmvtnorms",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(tmvtnorm)}
if(!require(lhs)){install.packages("lhs",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(lhs)}
if(!require(JuliaCall)){install.packages("JuliaCall",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(JuliaCall)}
if(!require(foreach)){install.packages("foreach",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(foreach)}
if(!require(snow)){install.packages("snow",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(snow)}
if(!require(doSNOW)){install.packages("doSNOW",dependencies=TRUE,lib = Sys.getenv("R_LIBS_USER"), repos='http://cran.rstudio.com/');library(doSNOW)}

set.seed(12345)
parallel <- ifelse(as.integer(cdargs[2])>1,TRUE,FALSE) # Should be run in parallel? 
julia <- julia_setup()

setwd(cdargs[1]) #Set directory to be script directory

# Origianl parameter set for data simulation
parorig <- c(10.0, 10.0, 20, 100, 0.1, 1.0)
# Simulate some fake data
julia_assign("paramraw", parorig)

julia_source('examples/Epidemiology/Scotland_run_inputs.jl')
simdata <- julia_eval("sum(abuns[71:80,:,:],dims=2)")


# Load in model functions
source("ABC_SMC_Scotland_preamble_server.R")
Dorig <- drop(simdata)
rm(simdata)
#### ABC set up ####

G <- as.integer(cdargs[3]) # Number of particle generations
N <- as.integer(cdargs[4]) # Number of accepted particles
K <- as.integer(cdargs[5])
nreps <- 50
ntimes <- dim(Dorig)[2]
nsum <- 2 # Number of individual summary statistics
nsumt <- as.integer(cdargs[6]) # Number of individual summary statistics including compartment specific ones
epsilon <- matrix(c(sapply(1:G,function(i)60000*(0.85**(i-1))),sapply(1:G ,function(i)10*(0.9**(i-1)))),nrow=2,byrow=T) # Epsilon value(s) of length nsum
n_par <- 6 # How many parameters will be estimated
names.par <- c("Beta env","Beta force","Birth_rate","Death_rate","Virus growth","Virus decay")
res <- matrix(ncol = n_par + nsumt, nrow = N) # Empty matrix to store results


Asingle <- randomLHS(K, n_par)
A <- matrix(rep(t(Asingle), nreps), ncol = ncol(Asingle), byrow = TRUE)

parrange <- matrix(c(c(5,15),c(5,15),c(10,25),c(90,110),c(0,0.5),c(.5,2)),ncol=n_par)

Ascaled <- sapply(1:n_par,function(i)A[,i]*(parrange[2,i]-parrange[1,i])+parrange[1,i])


# Number of simulations for each parameter set
n <- 1

if(parallel){
  ncores <- as.integer(cdargs[2])#detectCores()-1
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  
}else{
  ncores <- 1
}

# Empty matrices to store results (5 model parameters)
res.old<-matrix(ncol=n_par,nrow=N)
res.new<-matrix(ncol=n_par,nrow=N)

# Empty vectors to store weights
w.old<-matrix(ncol=1,nrow=N)
w.new<-matrix(ncol=1,nrow=N)
currsamp <- c()
for(g in 1:G){
  
  #Initiate counter
  i<-1
  
  while(i <= N){ # While the number of accepted particles is less than N_particles
    
    if(parallel){
      if(g==1){
        # Sample from prior distributions
        # Sample from prior distributions
        ll <- min(max(N-(i-1),ncores),N)
        ll <- ncores-ll%%ncores
        L <- sample(setdiff(1:K,currsamp),ll)
        param <- lapply(1:ll,function(x)Ascaled[L[x],])
        
      } else {
        ll <- min(max(N-(i-1),ncores),N)
        ll <- ncores-ll%%ncores
        #  Select particle from previous generation
        p<-sample(seq(1,N),ll,prob=w.old,replace=T)
        param<- lapply(1:ll,function(x)rK(res.old[p[x],],sigma))
      }
      
      m <- foreach(i = (1:ll), .combine='c') %dopar% {
        library(JuliaCall)
        setwd(as.character(cdargs[1]))
        runmodpar(param[[i]])

      }
      
    }else{
      if(g==1){
        # Sample from prior distributions
        # Sample from prior distributions
        L <- sample(setdiff(1:K,currsamp),1)
        param <- lapply(1,function(x)Ascaled[L[x],])
        currsamp <- c(currsamp,L)
      } else {
        #  Select particle from previous generation
        p<-sample(seq(1,N),1,prob=w.old,replace=T)
        param<- lapply(1,function(x)rK(res.old[p[x],],sigma))
      }
      m <- runmodpar(param[[1]])
    }
    
    for(kk in 1:ll){
      if(i<=N){
        if(parallel)currsamp <- c(currsamp,L[kk])
        if (m[kk]>0){
          # Store results
          res.new[i,]<-param[[kk]]
          # Calculate weights
          w1<-prod(sapply(1:n_par, function(b) dunif(res.new[i,b], min=parrange[1,b], max=parrange[2,b])))
          if(g==1){
            w2<-1
          } else {
            w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=parrange[1,], upper=parrange[2,])))
          }
          w.new[i] <- (m[kk]/n)*w1/w2
          print(paste0('Generation: ', g, ", particle: ", i))
          # Update counter
          i <- i+1
  
          
        }
      }
    }
  }
  
  
  
  
  #While loop ends
  sigma <- cov(res.new)
  res.old<-res.new
  w.old<-w.new/sum(w.new)
  
  write.csv(res.new, file = paste("results_case_2_ABC_SMC_gen_",g,".csv",sep=""), row.names=FALSE)
  par(mfrow=c(2,3))
  sapply(1:n_par,function(x){hist(res.new[,x],xlim=c(parrange[1,x],parrange[2,x]),main=names.par[x],probability=T)
    abline(v=parorig[x],col=2);
    abline(v=parrange[1,x],col=3);
    abline(v=parrange[2,x],col=3)
  })
  

}

stopCluster(cl)

rm(cl)
# save(list=ls(),file="ABCmetadata.RData")
save(list=c(res.new,w.new,Dorig,parorig,currsamp,epsilon,names.par),file="ABCmetadata.RData")
ABC_SMC <- read.csv(paste0("results_case_2_ABC_SMC_gen_",G,".csv"))
par(mfrow=c(2,3))
sapply(1:n_par,function(x){hist(ABC_SMC[,x],xlim=c(parrange[1,x],parrange[2,x]),main=names.par[x],probability=T)
  abline(v=parorig[x],col=2);
  abline(v=parrange[1,x],col=3);
  abline(v=parrange[2,x],col=3)
})

xx<-data.frame(cbind(parorig,apply(ABC_SMC[,1:n_par],2,mean),apply(ABC_SMC[,1:n_par],2,median),apply(ABC_SMC[,1:n_par],2,sd),(parorig-apply(ABC_SMC[,1:n_par],2,median))/apply(ABC_SMC[,1:n_par],2,median)),row.names = names.par)
names(xx) <- c("Simulated","Mean","Median","SD","Relative error")

print(xx)
