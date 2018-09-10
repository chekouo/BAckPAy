\name{Generate}
\alias{Generate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate independent features (e.g proteins or genes) with \emph{NbrModCov} and \emph{NbrGps} modalities of the independent and  experimental variables respectively (see the reference for more details).
}
\description{
Generate simulated data as esplained in the reference.
}
\usage{
Generate(NbrModCov = 2, NbrGps = 3, p = p, nbrDuplic = 2, bmin = 0.5, bmax = 1, 
smin = 0.2, smax = 0.2, seed = seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NbrModCov }{This is the number of modalities of the independent variable (\eqn{q=}NbrModCov)
}
  \item{NbrGps }{This is the number of modalities of the experimental variable (\eqn{T=}NbrGps)
}
  \item{p }{number of features (e.g. proteins or genes).
}
  \item{nbrDuplic }{number of duplicates for each combination \emph{Ind.Var}/\emph{Expe.Var}. For instance, if \emph{nbrDuplic=1}, only one sample is available for each combination \emph{Ind.Var}/\emph{Expe.Var}.
}
  \item{bmin }{Minimum of the coefficients effects \eqn{\beta}, which   are generated as Uniform(bmin,bmax).  
}
  \item{bmax }{Maximum of the coefficients effects \eqn{\beta}, which   are generated as Uniform(bmin,bmax).  
}
  \item{smin }{Minimum and maximum of the error variances of proteins in cluster,  \eqn{\sigma}, which are generated as Uniform(smin,smax).  
}
  \item{smax }{Maximum of the error variances of proteins in cluster,  \eqn{\sigma}, which are generated as Uniform(smin,smax).  
}
  \item{seed }{seed number for generating random numbers.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{data }{Expression data of size \eqn{p\times n}, with \eqn{p} the number of features (e.g proteins or genes) and \eqn{n} the number of samples  (\eqn{n=nbrDuplic*NbrModCov*NbrGps}).}
\item{Ind.Var }{A vector of size \eqn{n} with modalities \eqn{0,1,...,NbrModCov-1.}}
\item{Expe.Var }{A vector of size \eqn{n} with modalities \eqn{0,1,...,NbrGps-1.}}
\item{KnownRho }{A binary array of dimension \eqn{T\times p\times H}, of known (true) cluster memberships where KnownRho\eqn{[j,k,h]=1} if \eqn{\rho_{jk}=h} and 0 otherwise.} 
}
\references{
Thierry Chekouo et al (2018), \emph{Investigating Protein Patterns in Human Leukemia Cell Line Experiments: A Bayesian Approach for Extremely Small Sample Sizes}, \emph{submitted}.
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
Gen=Generate(NbrModCov=2,NbrGps=3,p=1000,nbrDuplic=1,seed=1)
dim(Gen$data)
round(head(Gen$data),digits=2)
IndVar=as.factor(Gen$Ind.Var);
ExpVar=as.factor(Gen$Expe.Var);
IndVar
ExpVar
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
