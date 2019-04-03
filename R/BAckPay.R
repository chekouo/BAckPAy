#' @useDynLib backpay mainfunction
#' @export
BAckPay <-
function(data=data, Ind.Var=Ind.Var, Expe.Var=NULL, sample=10000, burnin=1000,a.tau=a.tau,b.tau=b.tau,c.h=0.5,b.beta=0.1,alpha.h=c(20,20)) {
  #function(data=data, Ind.Var=Ind.Var, Expe.Var=NULL, sample=10000, burnin=1000,
           #a.tau=a.tau,b.tau=b.tau,c.h=0.5,b.beta=0.1)
if (ncol(data)!=length(Ind.Var)){
  stop("The number of columns in the data matrix must be the same as the length of vector Ind.Var")
}
  
  if (missing(a.tau)|| missing(b.tau)){
    stop("Arguments a.tau and b.tau must be specified.")
  }  
  if (a.tau<=0 || b.tau<=0){
    stop("Arguments a.tau and b.tau must be positive.")
  }
  
if (sample<=burnin){
stop("Argument burnin must be smaller than sample: the number of MCMC iterations.")
}
if (sample<=20){
stop("Please specify a larger number of MCMC iterations")
}

 if (!is.factor(Ind.Var)) {
    stop("Argument 'Ind.Var' must be of type 'factor'.")
  }

uniqInd.Var=sort(unique(Ind.Var))
Indvar1=rep(0,length(Ind.Var));
for (x in 1:length(Ind.Var)) {Indvar1[which(Ind.Var==uniqInd.Var[x])]=x-1;}

data=as.matrix(data)
p=nrow(data);n=ncol(data);
NbrGps=length(unique(Expe.Var));
ExpVar=rep(0,n)
if (is.null(Expe.Var)){
print("The experimental variable Expe.Var  is not specified.")
NbrGps=1;
Expe.Var=as.factor(rep(NbrGps,n))
} else {
  if (!is.factor(Expe.Var)) {
    stop("Argument 'Expe.Var' must be of type 'factor'.")
  }
  if (ncol(data)!=length(Expe.Var)){
    stop("The number of columns in the data matrix must be equal to the length of Expe.Var")
  }
uniqExp.Var=sort(unique(Expe.Var))
for (x in 1:length(uniqExp.Var)) {
  ExpVar[which(Expe.Var==uniqExp.Var[x])]=x-1;
  }
#ExpVar=unclass(Expe.Var)-1
}

IndVar=Indvar1[order(Expe.Var,Ind.Var)];
ExpVar=ExpVar[order(Expe.Var,Ind.Var)]
data=data[,order(Expe.Var,Ind.Var)]

data1=as.vector(t(data));
NbrCov=length(unique(IndVar))-1;

if (NbrCov<1){
stop("The independent variable must have at least two modalities or conditions")
}

H=3^NbrCov;
comb=H^NbrGps;
rhoMean1=rep(0,p*NbrGps*H)
ProbProt1=rep(0,p*comb);
probDiff=rep(0,p);
alphah=rep(alpha.h[1],H);alphah[(H+1)/2]=alpha.h[2];
result <- .C("mainfunction",
               p1=as.integer(p),n1=as.integer(n),NbrGps1=as.integer(NbrGps),ProtExp1=as.double(data1),
               Time=as.integer(IndVar),Resist=as.integer(ExpVar),sample=as.integer(sample), 
		burnin=as.integer(burnin),atau=as.double(a.tau),btau=as.double(b.tau),c_h=as.double(c.h),
		b_beta=as.double(b.beta),alpha=as.double(alphah),rhoMean1=as.double(rhoMean1), ProbProt1=as.double(ProbProt1),
		probDiff=as.double(probDiff),BetaMean1=as.double(rep(0,NbrCov*H)), sigma2Mean=as.double(rep(0,H)),
		muMean=as.double(rep(0,H)),Pih=as.double(rep(0,H)))
#,PACKAGE = "BAckPay");
ProbProtGrp=matrix(result$ProbProt1,comb,byrow=T);
BetaMean=matrix(result$BetaMean1,H,byrow=T)
#ProbProtGrp=matrix(result$ProbProt1,p,byrow=T);
NameClust=Nameclust(NbrCov+1,NbrGps);
rownames(ProbProtGrp)=NameClust$namegroup
probDiff=result$probDiff;
rhoMean<-aperm(array(result$rhoMean1,dim=c(H,p,NbrGps)))
dimnames(rhoMean)[[3]] <- NameClust$namecl;
if (length(rownames(data))==0){ 
dimnames(rhoMean)[[2]]<-1:nrow(data);
names(probDiff)=1:nrow(data);
colnames(ProbProtGrp)=1:nrow(data);
} else {
colnames(ProbProtGrp)=rownames(data);
dimnames(rhoMean)[[2]]<-rownames(data);
names(probDiff)=rownames(data);
}
q.valDiff=sapply(probDiff, function(t) sum(1-probDiff[probDiff>=t])/sum(probDiff>=t))
names(q.valDiff)=rownames(data);
return(list(data=data,ExpVar=ExpVar,IndVar=IndVar,probDiff=probDiff,
            q.valueDiff=q.valDiff,rhoMean=rhoMean,ProbProtGrp=ProbProtGrp, 
            BetaPost=BetaMean, muPost=result$muMean,sigmaPost=result$sigma2Mean,
            Pih=result$Pih/sum(result$Pih)))
}
