\name{BAckPay}
\alias{BAckPay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BAyesian mixture model for identifying Clusters of features (e.g. proteins) with similar ``pre-defined" expression PAtterns.
}
\description{
Calculate (i) the marginal posterior probabilities of inclusions for each feature, and (ii) estimated false discovery rate for detecting differential features.  
}
\usage{
BAckPay(data = data, Ind.Var = Ind.Var, Expe.Var = NULL, sample = 10000, 
burnin = 1000, a.tau = a.tau, b.tau = b.tau, c.h = 0.5, b.beta = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Continuous (\eqn{\log_2}-) expression data of size \eqn{p\times n}, with \eqn{p} the number of features (e.g proteins or genes) and \eqn{n} the number of samples  
}
  \item{Ind.Var}{
%%     ~~Describe \code{Ind.Var} here~~
Independent explanatory (categorial) variable. We aim to group features based on the change (or no change) in expression between the modalities of this variable. This argument must be of type "factor".
}
  \item{Expe.Var}{
%%     ~~Describe \code{Expe.Var} here~~
Experimental explanatory (categorial) variable. It's the "confounding" variable. The patterns are compared for every modality of Expe.Var. This argument must be of type "factor".
}
  \item{sample}{
%%     ~~Describe \code{sample} here~~
Total number of MCMC draws. It must be larger than burnin.
}
  \item{burnin }{ Number of draws to discard for burn-in
}
  \item{a.tau }{Shape  of \eqn{\tau_{hl}\sim \textrm{Gamma}(a_{\tau},b_{\tau})}, the truncated parameter of  the independent variable  effects \eqn{\beta_{hl}}. 
}
 \item{b.tau }{Rate hyperparameters of \eqn{\tau_{hl}\sim \textrm{Gamma}(a_{\tau},b_{\tau})}, the truncated parameter of  the independent variable  effects \eqn{\beta_{hl}}. We have \eqn{E(\beta_{hl})\geq a_{\tau}/b_{\tau}}. The choice of these hyperparameters is guided by the experimenter, who would need to decide on a threshold for biological significance. Note that \eqn{E(\beta_{hl})} can be interpreted as the \eqn{\log_2} fold change if the feature expression  is \eqn{\log_2} transformed. 
}
%  \item{b.tau}{
%%     ~~Describe \code{b.tau} here~~
%}
  \item{c.h }{Variance hyperparameter of the random effect \eqn{a_{jh}}, that captures positive correlations between samples of the same type.
}
  \item{b.beta }{Variance hyperparameter of the effects \eqn{\beta_{hl}}'s. This parameter is involved in capturing positive correlations between features within the same cluster.
}
}
\details{
The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}
\value{
 \item{probDiff}{ a numeric vector of the probability of differential expression for each feature in the data set.}
	\item{q.valueDiff}{ Estimated ``q values" of detecting differential expression.} %,  the smallest FDR (False discovery rate) at which a particular feature would just stay on the list of positives.
	\item{rhoMean}{ Estimated marginal posterior probabilities of cluster memberships, \eqn{\rho_{jk}}, \eqn{P(\rho_{jk}=h)}, where \eqn{k=1,...,T}, \eqn{h=1,....,H},  \eqn{H=3^{q-1}}  the number of clusters, \eqn{q} the number of modalities of \emph{Ind.Var}, and  \eqn{T} is the number of modalities of the experimental variable (\emph{Expe.Var}) . It's an array of size \eqn{T\times p\times H}.
	}
	\item{ProbProtGrp}{ Estimated (joint) marginal  posterior probabilities of proteins for patterns or groups (jMPP). It's a matrix of size \eqn{q\times p} where \eqn{q=H^{T}=3^{(q-1)*T}} is the total number of patterns or groups. For instance, if \emph{Ind.Var} has two modalities: Resistant and Sensitive, and the experimental variable has also two modalities: Time 0h and 48h. Hence group \emph{Up-Up} is the pattern or group of features that go up from Resistance to Sensitive for both times 0h and 48h.} 
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Generate}}
}
\examples{
library(backpay)
##---- Simulate data
Gen=Generate(NbrModCov=2,NbrGps=3,p=1000,nbrDuplic=1,seed=1)
data=Gen$data
IndVar=as.factor(Gen$Ind.Var);
ExpVar=as.factor(Gen$Expe.Var);
Result=BAckPay(data=data,Ind.Var=IndVar,Expe.Var=ExpVar, sample=10000,burnin=5000, a.tau=8,b.tau=10);
round(Result$probDiff[1:10],digits=4)
dim(Result$rhoMean)
round(head(Result$rhoMean[1,,]),digits=2)
dim(Result$ProbProtGrp)
print(round(head(Result$ProbProtGrp[,1:10]),digits=2))
library(AUC)
KnownRho=Gen$RhoKnown ## true clustering memberships
H=3^(length(unique(IndVar))-1);
rhodiffTrue=1-apply(KnownRho[,,(H+1)/2],2,prod)
AUC=auc(roc(Result$probDiff, as.factor(rhodiffTrue)))
AUC
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
