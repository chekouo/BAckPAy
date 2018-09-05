\name{datNorm}
\alias{datNorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function: Mean-centering data   for each modality of the experimental variable
}
\description{
This function mean-centers features for each modality of the experimental variable.    
}
\usage{
datNorm(data=data,Expe.Var=NULL,Ind.Var=Ind.Var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Continuous (\eqn{\log_2}-) expression data of size \eqn{p\times n}, with \eqn{p} the number of features (e.g proteins or genes) and \eqn{n} the number of samples  
}
  \item{Ind.Var}{
%%     ~~Describe \code{Ind.Var} here~~
Independent explanatory (categorial) variable. We aim to group features based on the change (or no change) in expression between the modalities of this variable. 
}
  \item{Expe.Var}{
%%     ~~Describe \code{Expe.Var} here~~
Experimental explanatory (categorial) variable. It's the "confounding" variable. The patterns are compared for every modality of Expe.Var. 
}
 
}

\value{
 \item{y }{Mean-centered data matrix \eqn{p\times n}.}
 \item{covlist }{Independent variable for each  modality of the experimental variable.}  
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
%\code{\link{Generate}}
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
