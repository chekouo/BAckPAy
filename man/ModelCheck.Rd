\name{ModelCheck}
\alias{ModelCheck}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Model checking using one-sample kolmogorov text.
}
\description{
Model checking  as explained in Section A.3 of the Supplementary material of the reference paper.
}
\usage{
ModelCheck(y, Backpay,Ind.Var=Ind.Var, Expe.Var=NULL,ch=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{y }{
Observed data expression for a specific feature (e.g protein). It's a vector.
}
\item{Backpay}{
A Backpay object obtained using the function \emph{BAckPay}.
}

\item{Ind.Var}{
%%     ~~Describe \code{Ind.Var} here~~
Independent explanatory (categorial) variable. We aim to group features based on the change (or no change) in expression between the modalities of this variable. This argument must be of type "factor".
}
\item{Expe.Var}{
%%     ~~Describe \code{Expe.Var} here~~
Experimental explanatory (categorial) variable. It's the "confounding" variable. The patterns are compared for every modality of Expe.Var. This argument must be of type "factor".
}
\item{ch}{Variance hyperparameter of the random effect \eqn{a_{jh}}, that captures positive correlations between samples of the same type.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{}{The function returns pvalues for each feature from the one-sample Kolmogorov Smirnov test.}
}
\references{
Thierry Chekouo et al (2019), \emph{Investigating Protein Patterns in Human Leukemia Cell Line Experiments: A Bayesian Approach for Extremely Small Sample Sizes}, \emph{submitted}.
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
library(backpay)
data("Isoform")
SurfaceTF1=Isoform$SurfaceTF1
SampleInfo=Isoform$SampleInfo
NameSurf=names(SurfaceTF1)
ExpVar=as.factor(SampleInfo$time[SampleInfo$Name \%in\% NameSurf])
IndVar=as.factor(SampleInfo$resistant[SampleInfo$Name \%in\% NameSurf])
Result=BAckPay(data=SurfaceTF1,Ind.Var=IndVar,Expe.Var=ExpVar, 
sample=10000,burnin=5000, a.tau=5,b.tau=10,c.h=0.01,b.beta=1);
Pval=NULL
for (j in 1:nrow(SurfaceTF1)){
YY=datNorm(data=SurfaceTF1[j,],Expe.Var=ExpVar,Ind.Var=IndVar)
Pval[j]=ModelCheck(as.vector(YY$y), Backpay=Result,Ind.Var=IndVar, 
Expe.Var=ExpVar,ch=0.01)
}
min(Pval)
boxplot(Pval)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
