\name{Nameclust}
\alias{Nameclust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Internal function: Define the names of the obtained cluster patterns from BAckPAy.
}
\description{
List all the cluster pattern names (or groups) with respect to the number of modalities of both the independent and  experimental variable. 
}
\usage{
Nameclust(NbrModCov, NbrGps)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NbrModCov }{Number of modalities of the independent variable.
}
  \item{NbrGps }{Number of modalities of the experimental variable (\eqn{T=}NbrGps)}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{namecl }{Cluster names with respect of the independent variable. For instance, if the independent variable has 3 modalities, then cluster \emph{UpDown} contains features that are up from modality 1 from 2, and down from 2 to 3.}
  \item{namegroup }{Pattern (or groups) names obtained with combinations of both  the independent and  experimental variables. If both variables have 2 modalities, then the group  pattern \emph{DownFlat-UpUp}  contains features that are DownFlat and UpUp for the first and second modality of the experimental variable respectively.}
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
Res=Nameclust(NbrModCov=3, NbrGps=2)
#List of cluster names with 3 modalities from the indep. variable
Res$namecl
#List of all cluster pattern (or groups) names with 3 and 2 modalities 
#from the indep. and experimental variables respectively. 
Res$namegroup
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
