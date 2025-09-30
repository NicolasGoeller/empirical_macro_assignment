#This function estimates the country-specific models -> prior=5 -> SSVS prior used in the paper
BVAR <- function(hyper,nr,prior,fperiods,dum,diff=0,pdensity=0,pfore=0,gW=gW,bigx=xglobal,presmpl=1,Daten,Dummy,nsave=5000,nburn=15000){
  require(matlab)
  require(coda)
  library(MASS) 
  require(bayesm)
  require(psych)
  require(MCMCpack)
  require(mvtnorm)
  require(mnormt)
  library(compiler)
  a_bar_1 <- hyper[1]
  a_bar_2 <- hyper[2]
  a_bar_3 <- hyper[3]
  #----------------------user input stuff---------------- 
  xglobal <- bigx
  gW <- gW
  nr <- nr
  p <- 1 #number of lags of the dependent variable
  pwex <- 1 #number of lags of weakly exogenous variables (will be included later) 
  pex <- 1 #number of lags on the exogenous variables
  cons <- 1 #1 includes constant, 0 exclude
  prior <- prior #1 1 diffuse, 2 litterman, 3 nat conjugate
  trend <- 1 #1 includes trend component i.e. alpha*t
  forecasting  <-  1     # 1: Compute h-step ahead predictions, 0: no prediction
  forecast_method <- 1 # 0: Direct forecasts 1: Iterated forecasts
  h  <- fperiods               # Number of forecast periods
  levels <- 0 #whether the model is estimated in levels or in differences (right now only differences are supported)
  exo <- 0 # no exogenous variables included
  dum <- dum # no dummies included
  diff <- diff #whether we use first differences for GDPPC
  pdensity <- pdensity #ONLY for analytical posteriors -> Simulating from the predictive density
  
  #--------------Gibbs-related prelims-------------------------
  ntot <- nsave+nburn # total draws
  #--------------------Reading in the Data----------------------
  #Data <- datahandling(Daten,nr)
  Names <- names(Daten)
  #Dummies <- Data[[2]]
  End <- xglobal[,substr(colnames(xglobal),1,2)==Names[[nr]]]
  Exogenous <- Dummy
  Dummies <- Dummy
  class(Dummies) <- "numeric"
  class(End) <- "numeric"
  class(Exogenous) <- "numeric"
  PIP <- NULL
  
  
  #xglobal <- globalV(Daten,diff=diff)
  class(xglobal) <- "numeric"
  W <- gW[[nr]]
  
  
  #-----------DIFFERENCED---------
  if (diff==1){
    End[2:nrow(End),1] <- diff(End[,1])
    End <- End[2:nrow(End),]
    Dummies <- Dummies[2:nrow(Dummies),]
    Exogenous <- Dummies
    all <- W%*%t(xglobal)
    Wex <- all[(ncol(End)+1):nrow(all),]
    Wex <- t(Wex)
    Wex[2:nrow(Wex),1] <- diff(Wex[,1])
    Wex <- Wex[2:nrow(Wex),]
    
  }else{
    all <- W%*%t(xglobal)
    Wex <- all[(ncol(End)+1):nrow(all),]
    Wex <- t(Wex)
    #Wex <- Wex[,!(colSums(abs(Wex)) == 0)]
  }
  
  #splitting the dataset -> PRE-Sample Data 
  class(Wex) <- "numeric"  
  #-----------------------Dummy selection routine--------------
  #----choose suitable dummies for each country model----------
  #----------------gets the last n chars of a string-----------
  if (dum==1){
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    
    actname <- Names[[nr]]
    dnames <- substr(colnames(Dummies),5,6)
    xnames <- substrRight(colnames(Exogenous),2)
    dbig <- cbind(Exogenous,Dummies)
    colnames(dbig) <- c((xnames),(dnames))
    IND <- matrix(NA,1,length(colnames(dbig)))
    IND <- as.logical(IND)
    
    for (i in 1:length(colnames(dbig))) if (grepl(actname,colnames(dbig)[[i]],ignore.case=TRUE)) IND[i] <- TRUE
    
    dummynames <- subset(colnames(cbind(Exogenous,Dummies)),IND)
    
    if (all(is.na(IND))) Dummies <- matrix(0,0,0) else Dummies <- t(subset(t(dbig),IND,drop=FALSE))
    
    #check if Dummy is the zero vector & delete it (could be the case if we cut the dummy periods (by))
    Dummies <- Dummies[,!(colSums(abs(Dummies)) == 0),drop=FALSE] 
    
    if (exo==1) Exogenous <- matrix(Exogenous[,1]) else Exogenous <- matrix(0,0,0)
    if (dum==1) Dummies <- Dummies else Dummies <- matrix(0,0,0)
    
    X <- End
    Yraw <- X
    
    if (ncol(Dummies)!=0){
      tempXs <- cbind(1,1:nrow(Yraw),Yraw,Dummies)
      tempXs <- tempXs[(p+1):nrow(tempXs),]
      cc <- try(solve(t(tempXs)%*%tempXs), silent=T) 
      if(is(cc,"try-error")) {  
        detc <- substr(colnames(Dummies),7,11)
        inc <- which(detc=="Dummy")
        Dummies <- Dummies[,inc,drop=FALSE]
        dummynames <- dummynames[inc]
      }
    }
  }
  
  X <- End
  Yraw <- X
  
  Yraw <- as.matrix(Yraw)
  class(Yraw) <- "numeric"
  #--------------------Data Preparation-------------------------
  dimYraw <- dim(Yraw) #get dimensions of the dependent variable
  Traw <- dimYraw[1]
  N <- dimYraw[2]
  
  #--------------------Model Specification----------------------
  Y1 <- Yraw
  Y2 <- Yraw
  #preliminary functions -------------------------------------
  wish <- function(h,n){
    A = t(chol(h))%*%matrix(rnorm(size(h,1)*n,0,1),size(h,1),n)
    A = A%*%t(A);
  }
  
  mlag <- function(X,lag)
  {
    p <- lag
    X <- as.matrix(X)
    Traw <- nrow(X)
    N <- ncol(X)
    Xlag <- matrix(0,Traw,p*N)
    for (ii in 1:p){
      Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
    }
    return(Xlag)  
  }
  
  if (exo==1) Exogenous <- matrix(Exogenous[,1]) else Exogenous <- matrix(0,0,0)
  if (dum==1) Dummies <- Dummies else Dummies <- matrix(0,0,0)
  #-------------------------------------------------------------
  #uses the function above
  Ylag <- mlag(Y2,p)
  Wexlag <- mlag(Wex,pwex)
  
  colnames(Wex) <- rep("Wex",ncol(Wex))
  
  #creates nametags for the lags, makes it easier to select the corresponding coefficient matrices later on
  nameslags <- NULL
  wexnameslags <- NULL
  for (ii in 1:p){
    
    nameslags <- c(nameslags,rep(paste("Ylag",ii,sep=""),ncol(Y2)))
    wexnameslags <- c(wexnameslags,rep(paste("Wexlag",ii,sep=""),ncol(Wex)))
  }
  colnames(Ylag) <- nameslags
  colnames(Wexlag) <- wexnameslags
  
  X1 <- cbind(Wex,Wexlag,Ylag)
  Exlag <- matrix(0,0,0)
  
  if (ncol(Dummies)!=0){ 
    colnames(Dummies) <- rep("Dummies",ncol(Dummies))
    X1 <- cbind(Dummies,X1)
    
  }else X1 <- X1
  
  X1 <- X1[(p+1):Traw,]
  #---------we need to include this for the trend component------------#
  if (trend==1){
    X1 <- cbind(1:nrow(X1),X1)
    colnames(X1)[[1]] <- "trend"
  } else {
    X1 <- X1
  }
  
  #----------------includes constant if chosen by the user-------------#
  if(cons==1) {
    X1 <- cbind(1,X1)
    colnames(X1)[[1]] <- "constant" 
  } else {
    X1 <- X1
  }
  
  #-----------------get size of the final matrix X---------------------#
  dimX <- dim(X1)
  Traw3 <- dimX[1]
  K <- dimX[2]
  T <- Traw-p
  
  #------------------create block diagonal matrix z-------------------#
  Z1 <- kronecker(diag(N),X1)
  Y1 <- Yraw[(p+1):Traw,]
  
  #----------------Forecasting SET-UP---------------------------------#
  #keeps the last "h" or 1 observations to evaluate forecasts
  Yfull <- Y1 #saves the full data
  Xfull <- X1 
  
  Y <- Y1[1:(nrow(Y1)-h),]
  X <- X1[1:(nrow(X1)-h),]
  Z <- kronecker(diag(N),X)
  T <- T-h 
  
  
  #----------------------PRIOR FUN------------------------------
  #----------------------getting OLS estimates------------------ 
  A_OLS  <-  ginv(t(X)%*%X)%*%(t(X)%*%Y)
  a_OLS <- as.vector(A_OLS) #vectorizes A_OLS, i.e. a_OLS=vec(A_OLS)
  SSE  <-  t((Y - X%*%A_OLS))%*%(Y - X%*%A_OLS)
  SIGMA_OLS  <-  SSE/(T-K+1)
  
  #-------------Initialize Bayesian Posteriors using OLS values-
  alpha <- a_OLS # single draw from the posterior of alpha
  ALPHA <- A_OLS # -,,- ALPHA
  SSE_Gibbs <- SSE # of SSE
  SIGMA <- SIGMA_OLS #of SIGMA
  #------------PRESAMPLE FUN-----------------------------------
  if (presmpl==1){
    if (nrow(Xpre)<=10 || is.null(dim(Xpre))){
      uniY <- Yraw 
    }else{
      uniY <- Xpre
    }
  }else{
    uniY <- Yraw
  }
  
  #store the draws:
  alpha_draws <- matrix(0,nsave,K*N)
  ALPHA_draws <- array(0,dim=c(nsave,K,N))
  SIGMA_draws <- array(0,dim=c(nsave,N,N)) 
  PL <- matrix(0,nsave,1)
  #-------------------PRIOR Hyperparameters--------------------
  n <- K*N
  
  gamma_draws <- matrix(0,nsave,n)
  M <- N
  #---------------------PRIOR SETUP FOR SSVS PRIOR----------------#
  # Priors on SIGMA
  # SSVS variances for non-diagonal elements of SIGMA     
  kappa_0  <-  10 #Set kappa_[0ij], kappa_[1ij]
  kappa_1  <-  10
  
  # Hyperparameters for diagonal elements of SIGMA (Gamma density)
  a_i  <-  0.01
  b_i  <-  0.01
  
  # Hyperparameters for Gamma ~ BERNOULLI(m,p_i), see eq. (14)
  p_i  <-  .5
  
  # Hyperparameters for Omega_[j] ~ BERNOULLI(j,q_ij), see eq. (17)
  q_ij  <-  .5
  
  # Initialize Gamma and Omega vectors
  gammas  <-  matrix(1,n,1)       # vector of Gamma
  omega <- list() #!!!!!!!
  for (kk_1 in 1:(M-1)){
    omega[[kk_1]]  <-  matrix(1,kk_1,1)  # Omega_j
  }
  
  bernoulli <- function(p){
    u <- runif(1)
    if (u<p){
      x=0
    }else{
      x=1
    } 
    return(x)
  }
  
  #-------------------POSTERIORS-----------------------------------------
  Xt <- X
  Yt <- Y
  T <- nrow(Xt); 
  #----------------------PRIOR FUN------------------------------
  #----------------------getting OLS estimates------------------ 
  XtXinv <- ginv(t(Xt)%*%Xt)
  A_OLS  <-  XtXinv%*%(t(Xt)%*%Yt)
  a_OLS <- as.vector(A_OLS) #vectorizes A_OLS, i.e. a_OLS=vec(A_OLS)
  SSE  <-  t((Yt - Xt%*%A_OLS))%*%(Yt - Xt%*%A_OLS)
  SIGMA_OLS  <-  SSE/(T-K+1)
  IXY  <-   kronecker(diag(M),(t(Xt)%*%Yt))
  #-------------Initialize Bayesian Posteriors using OLS values-
  alpha <- a_OLS # single draw from the posterior of alpha
  ALPHA <- A_OLS # -,,- ALPHA
  SSE_Gibbs <- SSE # of SSE
  SIGMA <- SIGMA_OLS #of SIGMA
  
  a_prior  <-  matrix(0,n,1)
  
  # This is the std of the OLS estimate ofalpha. You can use this to 
  # scale tau_0 and tau_1 (see below) if you want.
  # SSVS variances for alpha
  sigma_alpha  <-  sqrt(diag(kronecker(SIGMA_OLS,XtXinv)))
  
  tau_0  <-  0.1*sigma_alpha   # Set tau_[0i], tau_[1i]
  tau_1  <-  3*sigma_alpha
  
  
  for (k in 1:ntot){
    # Draw psi|alpha,gamma,omega,DATA from the GAMMA dist.
    # Get S_[j] - upper-left [j x j] submatrices of SSE
    # The following loop creates a cell array with elements S_1,
    # S_2,...,S_j with respective dimensions 1x1, 2x2,...,jxj
    S <- list() #!!!!!!!
    for (kk_2 in 1:M){                                      
      S[[kk_2]]  <-  SSE_Gibbs[1:kk_2,1:kk_2]  
    }
    # Set also SSE =(s_[i,j]) & get vectors s_[j]=(s_[1,j] , ... , s_[j-1,j])
    s <- list()
    for (kk_3 in 2:M){
      s[[kk_3 - 1]]  <-  SSE_Gibbs[1:(kk_3 - 1),kk_3]
    }
    # Parameters for Eta|omega ~ N_[j-1](0,D_[j]*R_[j]*D_[j]), see eq. (15)
    # Create and update h_[j] matrix
    # If omega_[ij] = 0 => h_[ij] = kappa0, else...
    hh <- vector('list',M-1)
    
    for (kk_4 in 1:(M-1)){                 
      omeg  <-  omega[[kk_4]]
      het  <-  hh[[kk_4]]
      for (kkk in 1:dim(omeg)[1]){           
        if (omeg[kkk,1] == 0){
          het[kkk]  <-  kappa_0
        }else{                        
          het[kkk]  <-  kappa_1              
        }
      }
      hh[[kk_4]]  <-  het
    }            
    # D_j = diag(hh_[1j],...,hh_[j-1,j])
    D_j <- vector('list',M-1)
    for (kk_5 in 1:(M-1)) {          
      D_j[[kk_5]]  <-  diag(hh[[kk_5]])
    }
    D_j[[1]] <- hh[[kk_5]][1]
    # Now create covariance matrix D_[j]*R_[j]*D_[j], see eq. (15)
    DD_j <- vector('list',M-1)
    for (kk_6 in 1:(M-1)){
      DD  <-  D_j[[kk_6]]
      DD_j[[kk_6]] = (DD%*%DD);
    }
    # Create B_[i] matrix
    B <- vector('list',M)
    for (rr in 1:M) {          
      if (rr == 1){
        B[[rr]]  <-  b_i + 0.5*(SSE[rr,rr]);
      }else{if (rr > 1){
        s_i  <-  s[[rr-1]]
        S_i  <-  S[[rr-1]]
        DiDi  <-  DD_j[[rr-1]]
        B[[rr]]  <-  b_i + 0.5*(SSE_Gibbs[rr,rr] - t(s_i)%*%solve(S_i + solve(DiDi))%*%s_i)
      }
      }
    }
    #Now get B_i from cell array B, and generate (psi_[ii])^2
    B_i  <-  as.matrix(unlist(B));
    psi_ii_sq  <- matrix(0,M,1)
    for (kk_7 in 1:M){                
      psi_ii_sq[kk_7,1]  <-  rgamma(1,a_i + 0.5*T,B_i[kk_7,1]) #gamm_rnd(1,1,(a_i + 0.5*T),B_i(1,kk_7));
    }
    
    # Draw eta|psi,phi,gamma,omega,DATA from the [j-1]-variate
    # NORMAL dist.
    eta  <-  vector('list',M-1)
    for (kk_8 in 1:(M-1)){       
      s_i  <-  s[[kk_8]]
      S_i  <-  S[[kk_8]]
      DiDi  <-  DD_j[[kk_8]]
      miu_j  <-  - sqrt(psi_ii_sq[kk_8+1])*(solve(S_i + solve(DiDi))%*%s_i);
      Delta_j  <-  solve(S_i + solve(DiDi))
      
      eta[[kk_8]]  <-  miu_j + t(chol(Delta_j))%*%rnorm(kk_8,0,1);
    }
    
    # Draw omega|eta,psi,phi,gamma,omega,DATA from BERNOULLI dist.
    omega_vec  <-  NULL #temporary vector to store draws of omega
    for (kk_9 in 1:(M-1)){       
      omeg_g  <-  omega[[kk_9]]
      eta_g  <-  eta[[kk_9]]
      for (nn in 1:length(omeg_g)){  # u_[ij1], u_[ij2], see eqs. (32 - 33)                          
        u_ij1  <-  (1/kappa_0)*exp(-0.5*((eta_g[nn])^2)/((kappa_0)^2))*q_ij
        u_ij2  <-  (1/kappa_1)*exp(-0.5*((eta_g[nn])^2)/((kappa_1)^2))*(1-q_ij)
        ost  <-  u_ij1/(u_ij1 + u_ij2)
        omeg_g[nn,1]  <-  bernoulli(ost)
        omega_vec  <-  rbind(omega_vec , omeg_g[nn,1])
      }
      omega[[kk_9]]  <-  omeg_g
    }
    # Create PSI matrix from individual elements of "psi_ii_sq" and "eta"
    PSI_ALL  <-  matrix(0,M,M)
    for (nn_1 in 1:M){ # first diagonal elements
      PSI_ALL[nn_1,nn_1] <-  sqrt(psi_ii_sq[nn_1,1]);   
    }
    for (nn_2 in 1:(M-1)){ # Now non-diagonal elements
      eta_gg = eta[[nn_2]]
      for (nnn in 1:dim(eta_gg)[1]){
        PSI_ALL[nnn,nn_2+1]  <-  eta_gg[nnn]
      }
    }
    # Create SIGMA
    SIGMA  <- solve(PSI_ALL%*%t(PSI_ALL))                
    # Draw alpha              
    # Hyperparameters for alpha|gamma ~ N_[m](0,D*D)
    h_i = matrix(0,n,1)   # h_i is tau_0 if gamma=0 and tau_1 if gamma=1
    for (nn_3 in 1:n){
      if (gammas[nn_3,1] == 0){               
        h_i[nn_3,1]  <-  tau_0[nn_3]
      }else{if (gammas[nn_3,1] == 1){        
        h_i[nn_3,1]  <-  tau_1[nn_3]       
      }
      }
    }
    D  <-  diag(as.vector(t(h_i)%*%diag(n)))# Create D. Here D=diag(h_i) will also do
    DD  <-  D%*%D   # Prior covariance matrix for Phi_m
    isig <- solve(SIGMA)
    psi_xx  <-  kronecker(solve(SIGMA),(t(Xt)%*%Xt))
    V_post  <-  solve(psi_xx + ginv(DD))
    
    visig <- as.vector(isig)
    a_post  <-  V_post%*%(IXY%*%visig + ginv(DD)%*%a_prior);
    alpha  <-  a_post + t(chol(V_post))%*%rnorm(n,0,1) # Draw alpha
    
    ALPHA  <-  matrix(alpha,K,M) # Draw of ALPHA
    
    # Draw gamma|phi,psi,eta,omega,DATA from BERNOULLI dist.    
    for (nn_6 in 1:n){
      u_i1  <-  (1/tau_0[nn_6])*exp(-0.5*(alpha[nn_6]/(tau_0[nn_6]))^2)*p_i        
      u_i2  <-  (1/tau_1[nn_6])*exp(-0.5*(alpha[nn_6]/(tau_1[nn_6]))^2)*(1-p_i)
      gst  <-  u_i1/(u_i1 + u_i2)
      if (gst=="NaN") gst <- 0 
      gammas[nn_6,1]  <-  bernoulli(gst) 
      gammas[is.na(gammas)] <- 1
    }
    # Save new Sum of Squared Errors (SSE) based on draw of ALPHA  
    SSE_Gibbs  <-  t(Yt - Xt%*%ALPHA)%*%(Yt - Xt%*%ALPHA)
    colnames(ALPHA) <- colnames(Y) 
    
    if (k>nburn){
      ALPHA_draws[k-nburn,,] <- ALPHA
      SIGMA_draws[k-nburn,,] <- SIGMA
      gamma_draws[k-nburn,] <- gammas
      
    }
    print(paste("Gibbs Sampling started:", round(k/ntot*100),"% finished"))
    print(paste("Burn-In-Phase:",if(k/nburn>=1){c("FINISHED!")}else{paste(round(k/nburn*100),"%")}))
  }
  dimnames(ALPHA_draws)=list(NULL,colnames(X),colnames(A_OLS))
  
  A_post <- apply(ALPHA_draws, c(2,3), mean)
  S_post <- apply(SIGMA_draws,c(2,3),mean)
  PIP <- matrix(apply(gamma_draws,2,mean),K,N)
  colnames(PIP) <- colnames(A_post);rownames(PIP) <- rownames(A_post)
  
  Res <- Y-X%*%A_post
  
  dims <- dimnames(ALPHA_draws)[[2]] 
  #splitting the output of the gibbs sampler for predictive density later on
  a0post <- ALPHA_draws[,which(dims=="constant"),]
  a1post <- ALPHA_draws[,which(dims=="trend"),] #coefficients on the trend
  a2post <- NULL
  alpha_var <- NULL
  dummiespost <- ALPHA_draws[,which(dims=="Dummies"),] #coefficients for the dummies
  postExpost <- ALPHA_draws[,which(dims=="Exogenous"),] #coefficients for the cont.exogenous
  postExlpost <- ALPHA_draws[,which(dims=="ExogenousLag"),] #coefficients for the lagged. ex
  Lambda0post <- ALPHA_draws[,which(dims=="Wex"),]#coefficients on WEX
  SIGMApost <- SIGMA_draws
  
  Lambdapost <- NULL
  Thetapost <- NULL
  
  for (jj in 1:p) {
    Lambdapost[[jj]] <- ALPHA_draws[,which(dims==paste("Wexlag",jj,sep="")),]
    Thetapost[[jj]]<- ALPHA_draws[,which(dims==paste("Ylag",jj,sep="")),]
  }
  
  
  
  if (pdensity==0){
    a0post <- NULL
    a1post <- NULL 
    a2post <- NULL
    alpha_var <- NULL
    dummiespost <- NULL 
    postExpost <- NULL 
    postExlpost <- NULL 
    Lambda0post <- NULL
    Lambdapost <- NULL 
    Thetapost <- NULL
    SIGMApost <- NULL 
  }
  
  #--------------Splits the posterior matrix & creates residuals---------------------
  dims <- dimnames(A_post)[[1]] 
  
  fit <- X%*%A_post
  
  #Res <- Y-X%*%A_post
  Lik <- 0
  for (jj in 1:T){
    Lik <- Lik+dmvnorm(Y[jj,],(X%*%A_post)[jj,],S_post,log=TRUE)
  }
  ML <- -0.5*(-2*Lik+ncol(ALPHA)*nrow(ALPHA)*log(T))
  
  
  a0 <- A_post[which(dims=="constant"),] #coefficients on the constant
  a1 <- A_post[which(dims=="trend"),] #coefficients on the trend
  a2 <- NULL
  dummies <- A_post[which(dims=="Dummies"),] #coefficients for the dummies
  postEx <- A_post[which(dims=="Exogenous"),] #coefficients for the cont.exogenous
  postExl <- A_post[which(dims=="ExogenousLag"),] #coefficients for the lagged. ex
  Lambda0 <- A_post[which(dims=="Wex"),]#coefficients on WEX
  Lambda <- A_post[which(dims=="Wexlag"),] #coefficients on Wexlag
  Theta <- A_post[which(dims=="Ylag"),]#coefficients on endogenous lagged
  
  Lambda <- NULL
  Theta <- NULL
  
  for (jj in 1:p) {
    Lambda[[jj]] <- A_post[which(dims==paste("Wexlag",jj,sep="")),]
    Theta[[jj]]<- A_post[which(dims==paste("Ylag",jj,sep="")),]
  }
  
  
  #Data Splitting cbind(Dummies,Exogenous,Exlag,Wex,Wexlag,Ylag)
  
  DDummies <- Xfull[,which(dims=="Dummies")]
  DExogenous <- Xfull[,which(dims=="Exogenous")]
  DExogenousL <- Xfull[,which(dims=="ExogenousLag")]
  DWex <- Xfull[,which(dims=="Wex")]
  DWexlag <- Xfull[,which(dims=="Wexlag")]
  DEndlag <- Ylag[(p+1):nrow(Ylag),]
  
  
  return(list(Dummies=DDummies,End=Yfull,Exogenous=DExogenous,Exlag=DExogenousL,Wex=DWex,Wexlag=DWexlag,Ylag=DEndlag,A_post=A_post,
              Res=Res,a0=a0,a1=a1,a2=a2,dummies=dummies,postEx=postEx,postExl=postExl,Lambda0=Lambda0,Lambda=Lambda,Theta=Theta,p=p,bigX=X,X=X,W=W,a0post=a0post,a1post=a1post,a2post=a2post,
              dummiespost=dummiespost,postExpost=postExpost,postExlpost=postExlpost,Lambda0post=Lambda0post,Lambdapost=Lambdapost,Thetapost=Thetapost,SIGMApost=SIGMApost,alpha_Var=alpha_var,PL=PL,ALPHA_draws=ALPHA_draws,
              ML=ML,SIGMA=S_post,PIP=PIP))
}
