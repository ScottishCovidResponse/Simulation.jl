#install.packages("JuliaConnectoR",dependencies=TRUE, repos='http://cran.rstudio.com/')

library(tmvtnorm)

library(lhs)
library(JuliaConnectoR)
library(parallel)

set.seed(12345)
parallel <- TRUE # Should be run in parallel? 

setwd(dirname(sys.frame(1)$ofile)) #Set directory to be script directory

# Origianl parameter set for data simulation
parorig <- c(10.0, 10.0, 20, 100, 0.1, 1.0)
# Simulate some fake data
juliaLet('paramraw = parorig',parorig=parorig)

datsim <- juliaEval("include(/Users/benswalloe/Simulation.jl/examples/Epidemiology/Scotland_run_inputs.jl)")
# simdata <- julia_eval("abuns")

# 
# # Load in model functions
# source("ABC_SMC_Scotland_preamble.R")
# Dorig <- simdata[41:48,,]
# Dorig <- apply(Dorig,c(1,3),sum) 
# rm(simdata)
# #### ABC set up ####
# 
# G <- 1 # Number of particle generations
# N <- 100 # Number of accepted particles
# K <- 1000
# nreps <- 50
# ntimes <- dim(Dorig)[2]
# nsum <- 2 # Number of individual summary statistics
# nsumt <- 9 # Number of individual summary statistics including compartment specific ones
# epsilon <- matrix(c(sapply(1:G,function(i)60000*(0.85**(i-1))),sapply(1:G ,function(i)10*(0.9**(i-1)))),nrow=2,byrow=T) # Epsilon value(s) of length nsum
# n_par <- 6 # How many parameters will be estimated
# names.par <- c("Beta env","Beta force","Birth_rate","Death_rate","Virus growth","Virus decay")
# res <- matrix(ncol = n_par + nsumt, nrow = N) # Empty matrix to store results
# 
# 
# Asingle <- randomLHS(K, n_par)
# A <- matrix(rep(t(Asingle), nreps), ncol = ncol(Asingle), byrow = TRUE)
# 
# parrange <- matrix(c(c(5,15),c(5,15),c(10,25),c(90,110),c(0,0.5),c(.5,2)),ncol=n_par)
# 
# Ascaled <- sapply(1:n_par,function(i)A[,i]*(parrange[2,i]-parrange[1,i])+parrange[1,i])
# 
# 
# # Number of simulations for each parameter set
# n <- 1
# 
# if(parallel){
#   ncores <- 2#detectCores()-1
#   myCluster <- makeCluster(ncores)
#   clusterEvalQ(myCluster, {
#     library(JuliaCall)
#     set.seed(12345)
#     
#     julia <- julia_setup()
#     
#     setwd("/Users/benswallow/Simulation.jl/")
#     JuliaCall:::.julia$cmd("using RCall")
#     
#     
#   })
#   clusterExport(myCluster, c("n","ntimes","nsum","nsumt","Dorig","epsilon","run_model","runmodpar","calc_distance"))
# }else{
#   ncores <- 1
# }
# 
# # Empty matrices to store results (5 model parameters)
# res.old<-matrix(ncol=n_par,nrow=N)
# res.new<-matrix(ncol=n_par,nrow=N)
# 
# # Empty vectors to store weights
# w.old<-matrix(ncol=1,nrow=N)
# w.new<-matrix(ncol=1,nrow=N)
# currsamp <- c()
# for(g in 1:G){
#   
#   #Initiate counter
#   i<-1
#   if(parallel)clusterExport(myCluster, c("g"))
#   while(i <= N){ # While the number of accepted particles is less than N_particles
#     
#     if(g==1){
#       # Sample from prior distributions
#       # Sample from prior distributions
#       ll <- min(N-(i-1),ncores)
#       L <- sample(setdiff(1:K,currsamp),ll)
#       param <- lapply(1:ll,function(x)Ascaled[L[x],])
#       currsamp <- c(currsamp,L)
#     } else {
#       #  Select particle from previous generation
#       p<-sample(seq(1,N),ll,prob=w.old,replace=T)
#       param<- lapply(1:ll,function(x)rK(res.old[p[x],],sigma))
#     }
#     if(parallel){
#       clusterExport(myCluster, c("param"))
#       m <- unlist(parLapply(myCluster,1:ll,function(ii){runmodpar(param[[ii]])}))
#     }else{
#       m <- runmodpar(param[[1]])
#     }
#     
#     for(kk in 1:ll){
#       if (m[kk]>0){
#         # Store results
#         res.new[i,]<-param[[kk]]
#         # Calculate weights
#         w1<-prod(sapply(1:n_par, function(b) dunif(res.new[i,b], min=parrange[1,b], max=parrange[2,b])))
#         if(g==1){
#           w2<-1
#         } else {
#           w2<-sum(sapply(1:N, function(a) w.old[a]* dtmvnorm(res.new[i,], mean=res.old[a,], sigma=sigma, lower=parrange[1,], upper=parrange[2,])))
#         }
#         w.new[i] <- (m[kk]/n)*w1/w2
#         print(paste0('Generation: ', g, ", particle: ", i))
#         # Update counter
#         i <- i+1
#         
#         
#         
#       }
#     }
#   }
#   
#   
#   
#   
#   #While loop ends
#   sigma <- cov(res.new)
#   res.old<-res.new
#   w.old<-w.new/sum(w.new)
#   
#   write.csv(res.new, file = paste("results_case_2_ABC_SMC_gen_",g,".csv",sep=""), row.names=FALSE)
#   par(mfrow=c(2,3))
#   sapply(1:n_par,function(x){hist(res.new[,x],xlim=c(parrange[1,x],parrange[2,x]),main=names.par[x],probability=T)
#     abline(v=parorig[x],col=2);
#     abline(v=parrange[1,x],col=3);
#     abline(v=parrange[2,x],col=3)
#   })
# }
# 
# stop(myCluster)
# # save(list=ls(),file="ABCmetadata.RData")
# ABC_SMC <- read.csv(paste0("~/Simulation.jl/results_case_2_ABC_SMC_gen_",G,".csv"))
# 
# sapply(1:n_par,function(x){hist(ABC_SMC[,x],xlim=c(parrange[1,x],parrange[2,x]),main=names.par[x],probability=T)
#   abline(v=parorig[x],col=2);
#   abline(v=parrange[1,x],col=3);
#   abline(v=parrange[2,x],col=3)
# })
# 
# xx<-data.frame(cbind(parorig,apply(ABC_SMC[,1:n_par],2,mean),apply(ABC_SMC[,1:n_par],2,median),(parorig-apply(ABC_SMC[,1:n_par],2,median))/apply(ABC_SMC[,1:n_par],2,median)),row.names = names.par)
# names(xx) <- c("Simulated","Mean","Median","Relative error")
# 
# print(xx)
# 
