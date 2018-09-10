\name{PlotThreeD}
\alias{PlotThreeD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot 3D of patterns/groups
}
\description{
In addition to the 3D plot, the function also returns (i) the list of the top features with their respective joint marginal probabilities, and (ii) the expression data used to make the plot.
}

\usage{
PlotThreeD(data=data,Ind.Var = Ind.Var, Expe.Var = NULL,ProbProtGrp=ProbProtGrp,
patternname=patternname,thres=thres)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data }{Continuous (\eqn{\log_2}-) expression data of size \eqn{p\times n}, with \eqn{p} the number of features (e.g proteins or genes) and \eqn{n} the number of samples. It should be the same dataset used in the function \emph{BAckPAy}.  
}
  \item{Ind.Var }{Independent explanatory (categorial) variable. We aim to group features based on the change (or no change) in expression between the modalities of this variable. It should be the same \emph{Ind.Var} used in the function \emph{BAckPAy}.
}
  \item{Expe.Var }{Experimental explanatory (categorial) variable. It's the "confounding" variable. The patterns are compared for every modality of Expe.Var. It should be the same \emph{Expe.Var} used in the function \emph{BAckPAy}. 
}
	\item{ProbProtGrp }{Estimated (joint) marginal  posterior probabilities of proteins for patterns or groups (jMPP). It's a matrix of size \eqn{q\times p} where \eqn{q=H^{T}} is the total number of patterns or groups. This is obtained from the BAckPay function.}
  \item{patternname }{Pattern name to plot
}
   \item{thres }{Threshold used on joint marginal probabilities (see paper) to plot features
}

}

\value{
 \item{ProbSort }{a numeric vector of the higher joint marginal posterior probabilities obtained with a threshold of \emph{thres}.}
	\item{data }{data used to make the plot; centered data of features belonging to pattern \emph{patternname}}
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
\code{\link{BAckPay}}
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
Names=Nameclust(NbrModCov=2, NbrGps=3)
Names$namegroup[1]
PlotResult= PlotThreeD(data=data,Ind.Var = IndVar, Expe.Var = ExpVar,ProbProtGrp=Result$ProbProtGrp,
patternname=Names$namegroup[1],thres=0.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
