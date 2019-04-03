\name{ClusterPerformance}
\alias{ClusterPerformance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compute three criteria for clustering performance after fitting BACkPAy.
}
\description{
This function computes three clustering performance criteria: \eqn{A_{wc}} is the area over the proportion of well-classified (WC)
feature profiles, \eqn{A_{mc}} is the area under the proportion of misclassifed (MC)
feature profiles, and \eqn{A_{nc}} is the area under the proportion of non-classified (NC) feature profiles
}
\usage{
ClusterPerformance(BackPay=BackPay,Generate=Generate,thres=seq(0,1,0.02))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{BackPay}{
A Backpay object obtained using the function \emph{BAckPay}.
}
\item{Generate}{
A Generate object obtained using the function \emph{Generate}.
}
\item{thres}{
A sequence of cut-off values  on the marginal probabilities for clustering.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{}{The function returns \eqn{A_{wc}}, \eqn{A_{mc}} and \eqn{A_{nc}} .}
}
\references{
Thierry Chekouo et al (2019), \emph{Investigating Protein Patterns in Human Leukemia Cell Line Experiments: A Bayesian Approach for Extremely Small Sample Sizes}, \emph{submitted}.
}

\examples{
	Gen=Generate(NbrModCov=2,NbrGps=3,p=500,nbrDuplic=2,varyslope="features",seed=1)
		data=Gen$data
		IndVar=as.factor(Gen$Ind.Var);
		ExpVar=as.factor(Gen$Expe.Var);
		Result=BAckPay(data=data,Ind.Var=IndVar,Expe.Var=ExpVar, sample=10000,burnin=5000, 
		a.tau=7,b.tau=10,c.h=0.5,b.beta=100);
	 ClusterPerformance(Result,Gen)
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
