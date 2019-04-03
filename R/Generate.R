#' @export
#' @useDynLib backpay generedata
Generate <-
function(NbrModCov=2,NbrGps=3,p=p,nbrDuplic=2,varyslope="features",
         prob=rep(1,3^(NbrModCov-1))/3^(NbrModCov-1),bmin=0.5,bmax=1,smin=0.2,smax=0.2,seed=seed){
  if (NbrModCov<=1){
    stop("The number of modalities of the independent variable must be at least 2.")
  }
  if ((NbrGps<=0) | (p<=0) | (nbrDuplic<=0) | (bmin<=0) |(bmax<=0) |(smin<=0) | (smax<=0)){
   stop("Arguments must be all positive") 
  }
  if (length(which(prob<=0))>=1){
    stop("Values of the argument prob must be all positive")
  }
  var=1
  if (sum(prob)!=1){
  prob=prob/sum(prob)
  }
  if (smax<smin){
    stop("Argument smax must not be smaller than smax.")
  }
  if (varyslope=="features"){
beta=runif(p,bmin,bmax)
  }
  else {
    beta=c(runif(floor(NbrGps/2),bmin,(bmin+bmax)/2),runif(NbrGps-floor(NbrGps/2),(bmin+bmax)/2,bmax));
    #beta=runif(NbrGps,bmin,bmax);
    print(beta)
    var=0;
  }
sig=runif(p,smin,smax)
n=nbrDuplic*NbrModCov*NbrGps;
Ind.Var=rep(0,n);Expe.Var=rep(0,n);
H=3^(NbrModCov-1)
Data <- .C("generedata",nbrTimeMod=as.integer(NbrModCov),inbrResistMod1=as.integer(NbrGps),
           TimePoint=as.integer(Ind.Var),Resist=as.integer(Expe.Var),p1=as.integer(p),
           yv=as.double(rep(0,n*p)),nbrsampl1=as.integer(nbrDuplic), vary1=as.integer(var),
           probInit=as.double(prob),beta=as.double(beta),seed=as.integer(seed),Rho=as.integer(rep(0,NbrGps*p*H)),sig=as.double(sig));
y=matrix(Data$yv,p,n,byrow=T);
RhoKnown<-aperm(array(Data$Rho,dim=c(H,p,NbrGps)));
list(data=y,Ind.Var=Data$TimePoint,Expe.Var=Data$Resist,RhoKnown=RhoKnown);
}
