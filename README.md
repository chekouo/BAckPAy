## BAckPAy 

This is the backpay package

License: GPL

For more information please contact thierry.chekouotekou@ucalgary.ca

Installation

1. SystemRequirements: GSL (GNU Scientific Library) C library.
2. Using the 'devtools' package:

> install.packages("devtools")
> library(devtools)
> install_github('chekouo/BAckPAy')

Usage and Example, see the pdf manual for more details. 

> library(backpay)
> Gen=Generate(NbrModCov=2,NbrGps=3,p=2000,nbrDuplic=1,seed=1)
> data=Gen$data
> IndVar=as.factor(Gen$Ind.Var);
> ExpVar=as.factor(Gen$Expe.Var);
> Result=BAckPay(data=data,Ind.Var=IndVar,Expe.Var=ExpVar, sample=10000,burnin=5000, a.tau=8,b.tau=10);
> round(Result$probDiff[1:10],digits=4)
> dim(Result$rhoMean)
> round(head(Result$rhoMean[1,,]),digits=2)
> library(AUC)
> KnownRho=Gen$RhoKnown ## true clustering memberships
> H=3^(length(unique(IndVar))-1);
> rhodiffTrue=1-apply(KnownRho[,,(H+1)/2],2,prod)
## AUC for detecting differential features
> AUC=auc(roc(Result$probDiff, as.factor(rhodiffTrue)))
> AUC
> Names=Nameclust(NbrModCov=2, NbrGps=3)
> Names$namegroup[1]
# [1] "Up-Up-Up"
### We represent the pattern Up-Up-Up.
> PlotResult= PlotThreeD(data=data,Ind.Var = IndVar, Expe.Var = ExpVar,ProbProtGrp=Result$ProbProtGrp,patternname=Names$namegroup[1],thres=0.5)


Reference: Thierry Chekouo et al (2018), Investigating Protein Patterns in Human Leukemia Cell Line Experiments:
A Bayesian Approach for Extremely Small Sample Sizes, submitted.

