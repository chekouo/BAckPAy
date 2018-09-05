#' @export
#' @useDynLib backpay generedata
Generate <-
function(NbrModCov=2,NbrGps=3,p=p,nbrDuplic=2,bmin=0.5,bmax=1,smin=0.2,smax=0.2,seed=seed){
beta=runif(p,bmin,bmax)
sig=runif(p,smin,smax)
n=nbrDuplic*NbrModCov*NbrGps;
Ind.Var=rep(0,n);Expe.Var=rep(0,n);
H=3^(NbrModCov-1)
Data <- .C("generedata",nbrTimeMod=as.integer(NbrModCov),inbrResistMod1=as.integer(NbrGps),TimePoint=as.integer(Ind.Var),Resist=as.integer(Expe.Var),p1=as.integer(p),yv=as.double(rep(0,n*p)),nbrsampl1=as.integer(nbrDuplic),beta=as.double(beta),seed=as.integer(seed),Rho=as.integer(rep(0,NbrGps*p*H)),sig=as.double(sig));
y=matrix(Data$yv,p,n,byrow=T);
RhoKnown<-aperm(array(Data$Rho,dim=c(H,p,NbrGps)));
list(data=y,Ind.Var=Data$TimePoint,Expe.Var=Data$Resist,RhoKnown=RhoKnown);
}
