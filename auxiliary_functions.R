#This function computes the weighting matrices used in the local estimation, nr = country index, mixed=1 if mixed weights are used, real and fin correspond to the N times N 
#weighting matrices (trade and financial based) and Data corresponds to the datalist new.data
getweights <- function(nr,mixed,real,fin,Data){
  data1 <- datahandling(Data,nr)
  cnames <- data1[[1]]
  xglobal <- Data[[1]]
  class(xglobal) <- "numeric"
  colnames(xglobal) <- paste(cnames[1],colnames(xglobal),sep=".")
  xglobal <- xglobal[,colSums(is.na(xglobal))<nrow(xglobal)]
  
  for (jj in 2:length(Data)){
    temp <- Data[[jj]];class(temp) <- "numeric"
    temp <- temp[,colSums(is.na(temp))<nrow(temp)]
    colnames(temp) <- paste(cnames[jj],colnames(temp),sep=".")
    xglobal <- cbind(xglobal,temp)
  }
  
  #xglobal <- xglobal[,!substr(colnames(xglobal),4,10)=="epepsus"];xglobal <- xglobal[,!substr(colnames(xglobal),4,10)=="p"];#xglobal <- xglobal[,!substr(colnames(xglobal),4,10)=="eq"] #deletes variables not needed
  colnames(xglobal) <- substr(colnames(xglobal),1,7); xglobal <- xglobal[,!duplicated(colnames(xglobal))];#xglobal <- xglobal[,-which(colnames(xglobal)=="EG.stir")]
  varnames <- substr(colnames(xglobal),4,7); varnames <- varnames[!duplicated(varnames)] #creates a dynamic list of variables
  varnames <- varnames[!varnames=="tb"]
  endnames <- substr(colnames(xglobal[,substr(colnames(xglobal),1,2)==cnames[[nr]]]),4,7);endnames <- endnames[!duplicated(endnames)]
  
  W <- matrix(0,length(varnames),ncol(xglobal));colnames(W) <- colnames(xglobal);rownames(W) <- varnames
  WeightsREAL <- real ; WeightsFIN <- fin 
  
  if (mixed==1){
    financials <- c("ltir.a","stir.a","eq","tc")
    wghtsreal1 <- WeightsREAL[nr,]; wghtsfin1 <- WeightsFIN[nr,]
    
    for (k in 1:(nrow(W)-1)){
      if (rownames(W)[[k]] %in% financials){
        W[k,substr(colnames(W),4,7)==rownames(W)[k]] <- wghtsfin1[substr(names(W[k,substr(colnames(W),4,7)==rownames(W)[k]] ),1,2)]
      }else{
        W[k,substr(colnames(W),4,7)==rownames(W)[k]] <- wghtsreal1[substr(names(W[k,substr(colnames(W),4,7)==rownames(W)[k]] ),1,2)]
      }
    }
    if (cnames[[nr]]=="US") W["poil", "US.poil"] <- wghtsreal1["US"] else W["poil", "US.poil"] <- wghtsreal1["US"]/wghtsreal1["US"]
    
  }else{
    wghtsreal1 <- WeightsREAL[nr,]
    for (k in 1:(nrow(W)-1)){
      W[k,substr(colnames(W),4,7)==rownames(W)[k]] <- wghtsreal1[substr(names(W[k,substr(colnames(W),4,7)==rownames(W)[k]] ),1,2)]
    }
  }
  
  if (cnames[[nr]]=="US") W["poil", "US.poil"] <- wghtsreal1["US"] else W["poil", "US.poil"] <- wghtsreal1["US"]/wghtsreal1["US"]
  #cross check for "NaN"
  W[which(W=="NaN")] <- 1
  #----------------------here we specify the part for the endogenous variabes-----------------------------------------#
  endoW <- matrix(0,length(endnames),ncol(W))
  endonr <- xglobal[,substr(colnames(xglobal),1,2)==cnames[nr]];rownames(endoW) <- colnames(endonr)
  colnames(endoW) <- colnames(W)
  namesW <- colnames(endoW)
  namesNr <- colnames(endonr)
  
  for (j in 1:nrow(endoW)){
    for (i in 1:length(namesW)){
      
      if (grepl(namesNr[[j]],namesW[[i]])) endoW[j,i]=1
      
    }
  }
  WfinNR <- rbind(endoW,W)
  
  WfinNR <- WfinNR[!(rowSums(abs(WfinNR)) == 0),]
  return(list(W=WfinNR,bigx=xglobal))
}
datahandling <- function(Data,nr){
  #---------------------------------------------------------------------------------------------------------------#  
  #Input:Data and the number of the country in the dataset which is the current home country
  #Output: Data_Country , a list, where Data_Country[[2]] refers to the true exogenous variables, Data_Country[[3]] 
  #are the weak exogenous variables and Data_Country[[4]] are the endogenous variables and Country[[1]] are the mn
  #This function feeds in DATA and creates a 1 x K vector which indicates which variables are exogenous, own lags 
  #lags of other variables
  #indicator 0 means home country, 1 denotes the rest of the world and 2 are truly exogenous variables
  #WARNING: original data needs labels for the countries in order to create the indicator vector
  #---------------------------------------------------------------------------------------------------------------# 
  Daten <- as.data.frame(Data) #reads in  data and converts it to a data.frame
  names  <-  colnames(Daten)
  cols <- ncol(Daten)
  indicator  <-  matrix("NA",1,ncol(Daten))
  ind  <- substr(names,1,2)
  
  number  <-nr #indicates which country is selected
  
  IND <- ind[!duplicated(ind)]
  
  for (i in 1:cols){
    
    if (grepl(IND[number],names[i]))  {
      
      indicator[i]=0 
      
    }else indicator[i]=1
    if(grepl(IND[length(IND)],names[i])) indicator[i]=2
    
  }
  
  sortingM <- rbind(indicator,as.matrix(Daten))
  
  #optional: would put the choosen country's columns at the start of the matrix, using indicator would be more efficient i guess
  sort1.data <- sortingM[,order(sortingM[1,],decreasing="TRUE")]
  
  #splitting the data into (true) exogenous, weak exogenous and system variables
  #could be more elegant, but it works
  ex <- as.matrix(which(sort1.data[1,]==2, arr.ind=T))
  wex <- as.matrix(which(sort1.data[1,]==1, arr.ind=T))
  end <- as.matrix(which(sort1.data[1,]==0, arr.ind=T))
  #Exogenous variables need extra treatment to extract dummys
  Exogenous <- sort1.data[2:nrow(sort1.data),ex]
  Exnames <- colnames(Exogenous)
  inddummy <- substr(Exnames,7,10)
  dummies <- matrix("NA",1,ncol(Exogenous))
  
  for (j in 1:ncol(Exogenous)){
    
    if (grepl("Dumm",Exnames[j])) dummies[j]=1 else dummies[j]=0
  }
  class(dummies) <- "numeric"
  
  Exogenous <- rbind(dummies,Exogenous)
  sort.Ex <- Exogenous[,order(Exogenous[1,],decreasing="TRUE")]
  
  tex <- as.matrix(which(sort.Ex[1,]==0),arr.ind=T)
  dum <- as.matrix(which(sort.Ex[1,]==1),arr.ind=T)
  
  Dummies <- sort.Ex[2:nrow(sort.Ex),dum]
  Exogenous <- sort.Ex[2:nrow(sort.Ex),tex]
  
  
  WExogenous <- sort1.data[2:nrow(sort1.data),wex]
  Endogenous <- sort1.data[2:nrow(sort1.data),end]
  
  assign(paste("WEx_",IND[number],sep=""),WExogenous) #creates a matrix WEx_Countryname
  assign(paste("END_",IND[number],sep=""),Endogenous) # same, we dont need this for Exogenous because it stays the same for all countries
  assign(paste("DataFULL_",IND[number],sep=""),rbind(indicator,as.matrix(Daten))) #returns the full matrix with a indicator row
  
  assign(paste("",IND[number],sep=""),list(IND,Dummies,Exogenous,assign(paste("WEx_",IND[number],sep=""),WExogenous),assign(paste("END_",IND[number],sep=""),Endogenous)))
  
  return(get(IND[nr]))
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

get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  if (nl==1) MM <- t(Beta_)
  
  MM <- rbind((Beta_),cbind(diag((nl-1)*nn), matrix(0,(nl-1)*nn,nn+2)),c(matrix(0,1,(nn*nl)),1,0),c(matrix(0,1,(nn*nl)+1),1))
  
  return(list(MM=MM,Jm=Jm))
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
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}
rmse <- function(x,y) sqrt(mean((x-y)^2))

wish <- function(h,n){
  A = t(chol(h))%*%matrix(rnorm(size(h,1)*n,0,1),size(h,1),n)
  A = A%*%t(A);
}