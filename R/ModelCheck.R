
ModelCheck=function(y, Backpay,Ind.Var=Ind.Var, Expe.Var=NULL,ch=0.5){
  if (length(y)!=length(Ind.Var)){
    stop("The number of data points in y must be the same as the length of vector Ind.Var")
  }
  
  if (!is.factor(Ind.Var)) {
    stop("Argument 'Ind.Var' must be of type 'factor'.")
  }
  
  uniqInd.Var=sort(unique(Ind.Var))
  Indvar1=rep(0,length(Ind.Var));
  for (x in 1:length(Ind.Var)) {Indvar1[which(Ind.Var==uniqInd.Var[x])]=x-1;}
  NbrGps=length(unique(Expe.Var));
  n=length(Ind.Var)
  ExpVar=rep(0,n)
  if (is.null(Expe.Var)){
    print("The experimental variable Expe.Var  is not specified.")
    NbrGps=1;
    Expe.Var=as.factor(rep(NbrGps,n))
  } else {
    if (!is.factor(Expe.Var)) {
      stop("Argument 'Expe.Var' must be of type 'factor'.")
    }
    if (length(y)!=length(Expe.Var)){
      stop("The number of data points in y  must be equal to the length of Expe.Var")
    }
    uniqExp.Var=sort(unique(Expe.Var))
    for (x in 1:length(uniqExp.Var)) {
      ExpVar[which(Expe.Var==uniqExp.Var[x])]=x-1;
    }
    #ExpVar=unclass(Expe.Var)-1
  }
  
  IndVar=Indvar1[order(Expe.Var,Ind.Var)];
  y=y[order(Expe.Var,Ind.Var)]
  nbrcov=length(unique(Ind.Var))-1
  
  if (nbrcov<1){
    stop("The independent variable must have at least two modalities or conditions")
  }
  H=3^nbrcov;
  Nk=table(ExpVar)
  Beta=Backpay$BetaPost;mu=Backpay$muPost;pi=Backpay$Pih;sigma2=Backpay$sigmaPost
  ystand=y;    Xcov=list();    C=matrix(1,nbrcov,nbrcov);
  C=C*lower.tri(C,diag = T)
  for (k in 1:NbrGps){
    IdK=Ind.Var[ExpVar==k-1]
    Xcov[[k]]=matrix(0,Nk[k],nbrcov)
    for (i in 1:Nk[k]){
      for (l in 1:nbrcov){
        Xcov[[k]][i,l]=(IdK[i]==l-1)
      }
    }
    Xcov[[k]]=Xcov[[k]]%*%C; 
    yk=y[ExpVar==k-1]
    prob=rep(0,H);
    for (h in 1:H){
      prob[h]=log(pi[h])+dmvnorm(yk, mean=as.vector(Xcov[[k]]%*%Beta[h,]+mu[h]),sigma=sigma2[h]*(ch+diag(Nk[k])),log = T)
    }
    h=seq(1,H)[prob==max(prob)]
    ystand[Expe.Var==k-1]=solve(t(chol(sigma2[h]*(ch+diag(Nk[k])))))%*%(yk-as.vector(Xcov[[k]]%*%Beta[h,]+mu[h]))
  }
  
  test <- ks.test(ystand, pnorm)
  test$p.value
}

pmixnorm <- function(y,Ind.Var=Ind.Var, Expe.Var=NULL,k, mu=0,Beta=0,sigma2=1){
  ## K is the number of experimental conditions
  ## H is the number of clusters
  cdf=1;
  #library("mvtnorm")
  NbrGrps=length(unique(Expe.Var))
  Nk=table(Expe.Var)
  Xcov=list();
  nbrcov=length(unique(Ind.Var))-1;  #H=3^nbrcov;
  C=matrix(1,nbrcov,nbrcov);
  C=C*lower.tri(C,diag = T)
  #for (k in 1:NbrGrps){
  IdK=Ind.Var[Expe.Var==k-1]
  Xcov=matrix(0,Nk[k],nbrcov)
  for (i in 1:Nk[k]){
    for (l in 1:nbrcov){
      Xcov[i,l]=(IdK[i]==l-1)
    }
  }
  Xcov=Xcov%*%C;
  yk=y[Expe.Var==k-1]
  pnorm(as.numeric(yk),mean=as.vector(Xcov%*%Beta+mu),sd=sqrt(sigma2))
  #pmvnorm(upper=as.numeric(yk), mean=as.vector(Xcov%*%Beta+mu), 
     #     sigma=sigma2*(ch+diag(Nk[k])))
  
}

