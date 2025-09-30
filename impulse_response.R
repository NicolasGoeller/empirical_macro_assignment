#Routines for impulses and sign restrictions - Note: Snowfall has to be installed to exploit parallel computing
sign_zha <- function(maxlag=length(gvarobj$F),G=gvarobj$G,F=gvarobj$F,err=gvarobj$err,x=gvarobj$x,horizon=40,ssave=10,covtype=0,shrinkage=0.9,sigs1=sigs){    
  require(abind)
  require(Matrix)
  covlist <- NULL
  names <- substr(rownames(x),1,2)
  names <- names[!duplicated(names)]
  N <- length(names)
  resids <- G%*%err
  rownames(resids) <- rownames(err) <- rownames(x)
  #creates a block diagonal global VC matrix
  for (i in 1:N){
    if (names[[i]]=="US"){
      covlist[[i]] <- t(chol(sigs1[[i]]))
    }else{
      covlist[[i]] <- sigs1[[i]]
    }
  }
  gcov <- bdiag(covlist) #block diagonal global covariance matrix
  Cm <- as.matrix(gcov)
  colnames(gcov) <- rownames(x); rownames(gcov) <- rownames(x)
  
  maxlag <- length(F)
  K <- nrow(x) #gvarobj$x
  
  lF <- array(0,c(K,K,maxlag))
  
  for (kk in 1:length(F)){
    lF[,,kk] <- F[[kk]]
  }
  
  PHIx <- array(0,c(K,K,maxlag+horizon+1));rownames(PHIx) <- colnames(PHIx) <- rownames(x);
  
  for (i in 1:maxlag){
    PHIx[,,i]  <-  matrix(0,K,K) #ok
  }
  PHIx[,,maxlag+1]  <-  diag(K)
  
  for (t in (maxlag+2):(maxlag+horizon+1)){
    acc = 0
    for (j in 1:maxlag){
      acc  <-  acc + lF[,,j]%*%PHIx[,,t-j]
    }
    PHIx[,,t]  <-  acc
  }
  
  #reindicize
  PHI  <-  PHIx[,,(maxlag+1):(maxlag+horizon+1)]
  invG <- solve(G)
  # preallocate matrix containing GIRFs
  eslct <- matrix(0,K,1)
  
  irfa  <- array(0,c(K,K,horizon+1))
  
  colnames(irfa) <- rownames(irfa) <- rownames(x)
  sirf <- NULL
  rotstore <- array(0,c(nrow(x),nrow(x),ssave))
  names <- substr(rownames(x),1,2)[!duplicated(substr(rownames(x),1,2))]
  h_sign <- 1:1
  sirf <- array(0,c(nrow(x),length(shockc),horizon+1,ssave))
  
  test <- list()
  sfInit(parallel=TRUE,cpus=4)
  sfExport(list=list("get_rot","x","irfa","horizon","x","K","h_sign","G","invG","Cm","PHI","restrictions","shockc"))
  test <- sfLapply(1:ssave,function(i) get_rot(i,irfa=irfa,horizon=horizon,K=K,x=x))
  sfStop()
  
  for (j in 1:ssave){
    sirf[,,,j] <- test[[j]]$IRF
    rotstore[,,j] <- test[[j]]$rot
  }
  
  irflev <- sirf
  rownames(irflev) <- rownames(x)
  if (ssave>1){
    signvars <- NULL
    signnames <- c(names(restrictions$D),names(restrictions$S),names(restrictions$OIL));signnames <- signnames[!duplicated(signnames)] # CHG
    
    for (jj in 1:(length(signnames))) signvars <- c(signvars,rownames(x)[which(substr(rownames(x),1,2)==signnames[[jj]])])
    
    irffp <- array(0,c(length(signvars),dim(irflev)[2],dim(irflev)[3],dim(irflev)[4]));rownames(irffp) <- (signvars)
    
    
    for (jj in 1:dim(irflev)[2]){ 
      for (k in 1:length(signvars)){
        for (i in 1:dim(irflev)[3]){
          for (j in 1:dim(irflev)[4]){
            irffp[signvars[[k]],jj,i,j] <- (irflev[signvars[[k]],jj,i,j] - median(irflev[signvars[[k]],jj,i,]))/sd(irflev[signvars[[k]],jj,i,])
          }
        }
      }
    }
    
    near <- NULL
    for (jj in 1:dim(irflev)[2]){
      for (k in 1:length(signvars)){
        for (j in 1:horizon){
          near <- c(near,which(abs(irffp[signvars[k],jj,j,]-0)==min(abs(irffp[signvars[k],jj,j,]-0))))
        }
      }
    }
    pickrot <- as.numeric(names(sort(-table(near)))[1])
    
    irfmed <- irflev[,,,pickrot];rownames(irfmed) <- rownames(x);colnames(irfmed) <- shockc
    rotpick <- rotstore[,,pickrot]
    
    st_impulses <- irfmed
  }else{
    st_impulses <- irflev
    rotpick <- rotS
  }
  return(list(IRF=st_impulses,R=rotpick))
}


get_rot <- function(nr,irfa=irfa,horizon=horizon,K=K,x=x){
  require(Matrix)
  sign_irf_oa <- irfa
  irfa1  <- array(0,c(K,K,horizon+1));colnames(irfa1) <- rownames(irfa1) <- rownames(x)
  rotS <- diag(nrow(x));colnames(rotS) <- rownames(rotS) <- substr(rownames(x),1,2)
  rotations <- list()
  PHIs <- irfa
  testsign <- array(0,c(length(which(substr(rownames(x),1,2)=="US")),nrow(x),3))
  
  colnames(testsign)  <- rownames(x)
  rownames(testsign) <- rownames(x)[grep("US",rownames(x))]
  icounter <- 0
  icheck <- 0
  a <- 0;b <- 0;c <- 0;d <- 0
  while(icheck<1){
    icounter <- icounter+1
    #print(c(a,b,c,d,icounter))
    sign_irf <- PHIs
    #icounter <- icounter+1
    #Step 1: Draw a rotation matrix
    k <- length(which(substr(rownames(x),1,2)=="US")) #where to impose the sign restrictions?
    A <- matrix(rnorm(k*k,0,1),k,k)
    qA <- qr(A)
    rotA <- qr.Q(qA)
    rotA <- rotA%*%diag(((diag(rotA)>0)-(diag(rotA)<0)))
    overall_log <- NULL
    colnames(rotS) <- rownames(rotS) <- substr(rownames(x),1,2)
    rotS[rownames(rotS)=="US",rownames(rotS)=="US"] <- rotA
    rownames(rotS) <- colnames(rotS) <- rownames(x)
    temp <- bdiag(rotS)%*%Cm
    
    temp  <-  invG%*%Cm%*%t(rotS)#%*%Cm
    slct <- diag(nrow(x))
    
    for (kk in 1:3){
      testsign[,,kk]  <-  as.matrix(PHI[grep("US",rownames(x)),,kk]%*%(temp)) #checks for the first three periods
    }
    nr <- 1
    # R>0, M<0, y<0, and P<0 for the irs periods.
    #a = (imf3hat(1:irs,3,3) > 0) .* (imf3hat(1:irs,4,3) < 0) .* (imf3hat(1:irs,1,3) < 0) .* (imf3hat(1:irs,2,3) < 0);
    #Impose AD shock
    a <- (testsign["US.y","US.y",1]>0)*(testsign["US.Dp","US.y",1]>0)*(testsign["US.stir","US.y",h_sign]>0)#*(testsign["US.poil","US.y",h_sign]>0)
    if (!(all(a)==1)){
      next
    }
    #impose restrictions on AS shock
    b <- (testsign["US.y","US.Dp",1]<0)*(testsign["US.Dp","US.Dp",1]>0)*(testsign["US.stir","US.Dp",1]>0)
    if (!(all(b)==1)){
      next
    }
    #impose restrictions on MP shock
    c <- (testsign["US.y","US.stir",1]<0)*(testsign["US.Dp","US.stir",1]<0)*(testsign["US.stir","US.stir",h_sign]>0)#*(testsign["US.poil","US.Dp",h_sign]<0)
    if (!(all(c)==1)){
      next
    }
    #impose Oil price shock restriction
    #d <- (testsign["US.y","US.poil",h_sign]<0)*(testsign["US.Dp","US.poil",h_sign]>0)*(testsign["US.stir","US.poil",h_sign]>0)#*(testsign["US.poil","US.poil",h_sign]>0)
    #if (max(d)==0){
    #  next
    #}
    icheck <- 1
  }
  temp  <- invG%*%Cm%*%t(rotS)
  for (ii in 1:(horizon+1)){
    sign_irf_oa[,,ii]  <-  PHI[,,ii]%*%as.matrix(temp) #checks for the first three periods
  }
  irflev <- sign_irf_oa[,shockc,]
  rownames(irflev) <- rownames(x)
  
  st_impulses <- irflev
  rotpick <- rotS
  return(list(IRF=st_impulses,rot=rotS))
}

#Restrictions, not needed actually, restrictions are imposed within the functions above in a static fashion
restrictions_supply <- list(US=list(var=c("US.y","US.Dp","US.stir"),sign=c(-1,1,1),horizon=list(h1=1,h2=1,h3=1)))

restrictions_demand <- list(US=list(var=c("US.y","US.Dp","US.stir","US.poil"),sign=c(1,1,1,1),horizon=list(h1=1,h2=1,h3=1,h4=1)))

restrictions_money <- list(US=list(var=c("US.y","US.Dp","US.stir","US.ltir","US.poil"),sign=c(-1,-1,1,-1,-1),horizon=list(h1=1,h2=1,h3=1,h4=1)))

restrictions_oil <- list(US=list(var=c("US.y","US.p","US.FFR","US.poil"),sign=c(-1,1,1,1),horizon=list(h1=1,h2=1,h3=1,h4=1)))

restrictions=list(D=restrictions_demand,S=restrictions_supply,OIL=restrictions_oil,M=restrictions_money)