# Parts of this code are based on Koop & Korobilis (2012)
# Florian Huber 2014 (fhuber7@gmail.com)
#Install R packages used later on
#install.packages("snowfall","compiler","matlab","coda","MASS","bayesm","psych","MCMCpack","mvtnorm","mnormt")
library(compiler)
require(snowfall)
enableJIT(2)
rm(list = ls()) #clears workspace
#Step I: Read in the data, weights and dummies
#-------------------------------------------------------------------------
wd <- getwd()

load(paste(wd,"data/Replication/dummies.RData",sep="/"))
load(paste(wd,"data/Replication/data.RData",sep="/"))
load(paste(wd,"data/Replication/weights.RData",sep="/"))

#Source helper scripts, BVAR estimation routines, IR/FEVD scripts etc.
source("Replication/auxiliary_functions.R")
source("Replication/bvarx_estimation.R")
source("Replication/solve_gvar.R")
source("Replication/impulse_response.R")
#-------------------------------------------------------------------------
#Step II: Creates a list object for the country-specific weights
#Creates a list object for the country-specific weights
gW <- list();Wstorage <- list()

cN <- names <- names(new.data)

xglobal <- getweights(1,1,Data=new.data,fin=W.list$finW0711,real=W.list$tradeW.12)
W1 <- xglobal$W
xglobal <- xglobal$bigx
#xglobal <- xglobal[64:nrow(xglobal),] #ONLY for long datasets

for (k in 1:length(names)){
  temp <- getweights(k,1,Data=new.data,fin=W.list$finW0711,real=W.list$tradeW.12)
  tempW <- apply(temp$W,2,function(x) x/(rowSums(temp$W))) #row standardize
  gW[[k]] <- tempW
}

#-------------------------------------------------------------------------
#Step III: Runs MCMC using snowfall for each country-specific model 
#Additional settings
horizon <- 0 #IRF, full sample is used
par <- 1  #IMPORTANT: snowfall required, implements parallel computing
diff <- 0 #if estimated in differences
predDraw <- 550 # 4 #draws from the global posterior - I think this number is set to 550 in the paper
nsave <- 10000 # 5 #saved gibbs draws at the local level - I think in the paper this is 10000
nburn <- 5000 # 15 #number of burn-ins at the local level - I think in the paper this is 5000
irhorizon <- 20
dum <- 1
Daten <- new.data
IRhorz <- 20
shockc <- c("US.y","US.Dp","US.stir","US.poil")
ssave <- 10 #how many rotation matrices do we have to save?

#Note: Prior settings can be chosen directly within the BVAR estimation routine
#-------------------------------------------------------------------------
globalG <- list()
predDens <- list()

sfInit(parallel=TRUE,cpus=4)
sfExport(list=list("cN","datahandling","Daten","xglobal","BVAR","horizon","dum","diff","rmse","gW","Dummy","nburn","nsave"))
predDens <- sfLapply(1:length(names),function(i) BVAR(NULL,i,NULL,horizon,dum=1,diff=diff,pdensity=1,pfore=0,gW=gW,bigx=xglobal,presmpl=0,Daten,Dummy,nburn=nburn,nsave=nsave))
sfStop()
#-------------------------------------------------------------------------
#Step IV: Perform IRF using the global posterior / Monte Carlo integration
#creates the global x matrix
VAR0 <- predDens[[1]]
x <- t(VAR0$End)
for (i in 2:length(names)){
  VAR1 <- predDens[[i]]
  x <- rbind(x,t(VAR1$End))
}

Liks <- matrix(0,predDraw,1)
IRF_post <- array(NA,dim=c(predDraw,nrow(x),length(shockc),IRhorz))

for (ksim in 1:predDraw){
  
  draw <- sample(1:dim(predDens[[1]]$ALPHA_draws)[1],1,replace=FALSE)
  
  for (i in 1:length(names)) globalG[[i]] <- predDens[[i]]$W
  
  VAR0 <- predDens[[1]]
  #stacking/solving like in pesaran
  #parameter stuff
  A0 <- cbind(diag(ncol(VAR0$End)),-t(VAR0$Lambda0post[draw,,]))
  W0 <- globalG[[1]]
  
  #creates matrices according to PSW 2004 for every lag
  for (kk in 1:length(VAR0$Theta)){
    assign(paste("B0",kk,sep=""),cbind(t(VAR0$Thetapost[[kk]][draw,,]),t(VAR0$Lambdapost[[kk]][draw,,])))
    assign(paste("H0",kk,sep=""),get(paste("B0",kk,sep=""))%*%W0)
  }
  
  a0 <- matrix(VAR0$a0post[draw,])
  a1 <- matrix(VAR0$a1post[draw,])
  G <- A0%*%W0
  #data stacking
  x <- t(VAR0$End)
  Ylag1 <- t(VAR0$Ylag)
  
  resids <- t(VAR0$End-VAR0$X%*%VAR0$ALPHA_draws[draw,,])
  #####COVTYPE WEG
  sigs <- list()
  sigs[[1]] <- VAR0$SIGMApost[draw,,]
  for (i in 2:length(names)){
    #N might create a list containing all country specific W matrices and use that in here (should be much faster)
    VAR1 <- predDens[[i]]
    sigs[[i]] <- VAR1$SIGMApost[draw,,]
    W1 <- globalG[[i]]
    A1 <- cbind(diag(ncol(VAR1$End)),-1*t(VAR1$Lambda0post[draw,,]))
    
    for (kk in 1:length(VAR0$Theta)){
      assign(paste("B1",kk,sep=""),cbind(t(VAR1$Thetapost[[kk]][draw,,]),t(VAR1$Lambdapost[[kk]][draw,,])))
      assign(paste("H0",kk,sep=""),rbind(get(paste("H0",kk,sep="")),get(paste("B1",kk,sep=""))%*%W1))
    }
    
    #B1 <- cbind(t(VAR1$Thetapost[draw,,]),t(VAR1$Lambdapost[draw,,]))
    G <- rbind(G,A1%*%W1)
    #H <- rbind(H,B1%*%W1)
    a0 <- rbind(a0,matrix(VAR1$a0post[draw,]))
    a1 <- rbind(a1,matrix(VAR1$a1post[draw,]))
    #Variables
    Ylag1 <- rbind(Ylag1,t(VAR1$Ylag))
    resids <- rbind(resids,t(VAR1$End-VAR1$X%*%VAR1$ALPHA_draws[draw,,]))
    x <- rbind(x,t(VAR1$End))
  }
  
  class(x) <- "numeric"
  
  b2 <- NULL
  b0 <- solve(G)%*%a0 #global constant
  b1 <- solve(G)%*%a1 #global trend
  
  for (kk in 1:length(VAR0$Theta)){
    assign(paste("F",kk,sep=""),solve(G)%*%get(paste("H0",kk,sep="")))
  }
  
  err <- solve(G)%*%resids
  F_post <- list()
  for (kk in 1:length(VAR0$Theta)){
    F_post[[kk]] <- get(paste("F",kk,sep=""))
  }
  
  #this stuff is needed to retrieve the companion form to check for stationarity
  p <- length(VAR0$Theta)
  nn <- nrow(x)
  nd <- 2
  nl <- p
  #----------------
  ALPHA <- NULL
  for (jj in 1:p){
    ALPHA <- cbind(ALPHA,F_post[[jj]])
  }
  ALPHA <- cbind(ALPHA,b0,b1)
  
  varndxv <- c(nn,nd,nl)
  nkk <- (nl*nn)+nd
  
  comp_parms <- get_companion(ALPHA,varndxv)
  eigs_parms <- eigen(comp_parms$MM)$values
  critval <- 1.05 #reject unstable draws, only necessary in a low percentage of the draws
  print(max(abs(Re(eigs_parms))))
  if (ksim==1){ 
    F_post1 <- F_post
    G_1 <- G
    ir_draw <- sign_zha(maxlag=p,G=G,F=F_post,err=err,x=x,horizon=20,ssave=ssave,covtype=0,shrinkage=0.9,sigs=sigs)  #covtype and shrinkage have no actual effect
    IRF_post[ksim,,,] <- ir_draw$IRF[,shockc,1:irhorizon]
    next
  }
  if (max(abs(Re(eigs_parms)))>=critval){
    IRF_post[ksim,,,] <- IRF_post[ksim-1,,,]
    print(ksim)
    next
  }else{
    F_post1 <- F_post
    G_1 <- G
    err_1 <- err
    sigs1 <- sigs
  }
  ir_draw <- sign_zha(maxlag=p,G=G,F=F_post,err=err,x=x,horizon=20,ssave=ssave,covtype=0,shrinkage=0.9,sigs=sigs)   
  
  IRF_post[ksim,,,] <- ir_draw$IRF[,shockc,1:irhorizon]
  print(ksim)
}
#IRF_post is the posterior of structural IRF, 1st dim: Draws, 2nd dim: Response, 3rd.dim: shock, 4rd. dim: horizon
dimnames(IRF_post)=list(NULL,rownames(x),shockc,NULL) 
saveRDS(IRF_post,"data/BGVAR_IRF_estimates.rds")





