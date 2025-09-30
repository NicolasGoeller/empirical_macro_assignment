#This function estimates the GVAR with SSVS prior and reports the posterior mean estimates
stacking <- function(hyper,prior=1,horizon=10,fhorz=1,levels=0,exo=0,dum=1,par=1,diff=0,pdensity=0,gW=gW,xglobal=xglobal,presmpl=1,Daten=Daten,Dummy,nsave,nburn){
  require(snowfall)
  require(forecast)
  require(bayesm)
  require(psych)
  require(MCMCpack)
  require(mvtnorm)
  require(mnormt)
  library(compiler)
  require(Matrix)
  
  names <- names(Daten)
  
  globalG <- list()
  allVARX <- list()
  
  if (par==0){
    for (i in 1:length(names)){
      allVARX[[i]]  <- BVAR(hyper,i,prior,horizon,dum,diff=diff,pdensity=0,pfore=0,gW=gW,bigx=xglobal,presmpl=0,Daten,Dummy)
    }
  }else{
    sfInit(parallel=TRUE,cpus=4)
    sfExport(list=list("nsave","nburn","cN","hyper","datahandling","Daten","xglobal","BVAR","prior","horizon","dum","diff","rmse","pdensity","remove_outliers","gW","Dummy"))
    allVARX <- sfLapply(1:length(names),function(i) BVAR(hyper,i,prior,horizon,dum,diff=diff,pdensity=0,pfore=0,gW=gW,bigx=xglobal,presmpl=0,Daten,Dummy,nsave=nsave,nburn=nburn))
    sfStop()
  }
  
  
  for (i in 1:length(names)) globalG[[i]] <- allVARX[[i]]$W
  
  VAR0 <- allVARX[[1]]
  #stacking/solving like in pesaran
  #parameter stuff
  A0 <- cbind(diag(ncol(VAR0$End)),-t(VAR0$Lambda0)) #includes contemp. wex
  W0 <- globalG[[1]]
  
  sigs <- list()
  sigs[[1]] <- VAR0$SIGMA
  
  #creates matrices according to PSW 2004 for every lag
  for (kk in 1:max(length(VAR0$Theta),length(VAR0$Lamda))){
    cc <- try(VAR0$Lambda[[kk]], silent=T);bb <-  try(VAR0$Theta[[kk]], silent=T);
    if(is(bb,"try-error")) {
      VAR0$Theta[[kk]] <- VAR0$Theta[[1]]*0
    }
    if (is(cc,"try-error")){
      VAR0$Lambda[[kk]] <- VAR0$Lambda[[1]]*0
    }
    assign(paste("B0",kk,sep=""),cbind(t(VAR0$Theta[[kk]]),t(VAR0$Lambda[[kk]])))
    assign(paste("H0",kk,sep=""),get(paste("B0",kk,sep=""))%*%W0)
  }
  
  
  G <- A0%*%W0
  #H01 <- B01%*%W0
  #H02 <- B01%*%W0
  a0 <- matrix(VAR0$a0)
  a1 <- matrix(VAR0$a1)
  #data stacking
  x <- t(VAR0$End)
  Ylag1 <- t(VAR0$Ylag)
  resids <- t(VAR0$Res)
  
  if(pdensity==1) PL <- log(mean(VAR0$PL))
  
  for (i in 2:length(names)){
    #N might create a list containing all country specific W matrices and use that in here (should be much faster)
    VAR1 <- allVARX[[i]]
    sigs[[i]] <- VAR1$SIGMA
    A1 <- cbind(diag(ncol(VAR1$End)),-1*t(VAR1$Lambda0))
    W1 <- globalG[[i]]
    
    for (kk in 1:length(VAR1$Theta)){
      
      cc <- try(VAR1$Lambda[[kk]], silent=T);bb <-  try(VAR1$Theta[[kk]], silent=T);
      if(is(bb,"try-error")) {
        VAR1$Theta[[kk]] <- VAR1$Theta[[1]]*0
      }
      if (is(cc,"try-error")){
        VAR1$Lambda[[kk]] <- VAR1$Lambda[[1]]*0
      }
      
      assign(paste("B1",kk,sep=""),cbind(t(VAR1$Theta[[kk]]),t(VAR1$Lambda[[kk]])))
      assign(paste("H0",kk,sep=""),rbind(get(paste("H0",kk,sep="")),get(paste("B1",kk,sep=""))%*%W1))
      
    }
    G <- rbind(G,A1%*%W1)
    
    a0 <- rbind(a0,matrix(VAR1$a0))
    a1 <- rbind(a1,matrix(VAR1$a1))
    #Variables
    Ylag1 <- rbind(Ylag1,t(VAR1$Ylag))
    resids <- rbind(resids,t(VAR1$Res))
    x <- rbind(x,t(VAR1$End))
    if (pdensity==1) PL <- rbind(PL,log(mean(VAR1$PL)))
  }
  
  class(x) <- "numeric"
  bigx <- t(x)
  
  b2 <- NULL
  b0 <- solve(G)%*%a0 #global constant
  b1 <- solve(G)%*%a1 #gloabl trend CHGCHG
  F_post <- list()
  
  for (kk in 1:length(VAR0$Theta)){
    F_post[[kk]] <- solve(G)%*%get(paste("H0",kk,sep=""))
  }
  
  err <- solve(G)%*%resids #residuals
  
  bigx <- cbind(1,seq(1:ncol(x)),t(Ylag1[order(rownames(Ylag1)),]))
  
  ALPHA <- cbind(b0,b1,F_post[[1]])
  
  fit <- ALPHA%*%t(bigx)
  gcov <- bdiag(sigs) #block diagonal global covariance matrix
  Sig <- solve(G)%*%as.matrix(gcov)%*%t(solve(G))
  
  PLS <- 0
  for (i in 2:(ncol(x)-horizon)){
    PLS <- PLS+sum(dmvnorm((x[,i]), bigx[i,]%*%t(ALPHA),Sig,log=TRUE))
  }
  
  ML <- -0.5*(-2*PLS+ncol(ALPHA)*nrow(ALPHA)*log(ncol(x)))
  
  return(list(b0=b0,b1=b1,b2=b2,G=G,F=F_post,err=err,x=x,ML=ML,allVARX=allVARX,sigs=sigs))
}