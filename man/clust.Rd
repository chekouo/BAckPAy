\name{clust}
\alias{clust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Internal function: List possible clusters along with the signs of the effect of each covariate.
}
\description{
Provide all possible clusters 
}
\usage{
clust(VV=c(1,0,-1), nbrcov=nbrcov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VV}{
Possible signs of  the effects in each cluster.
}
  \item{nbrcov}{
Number of covariates included in the model. It's NbrModCov-1 where NbrModCov is the number of modalities of the independent variable. 
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
It gives a matrix of dimension \eqn{H\times 3}, where \eqn{H=3^{nbrcov}}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
clust(nbrcov=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
