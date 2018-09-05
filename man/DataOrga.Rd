\name{DataOrga}
\alias{DataOrga}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal Function: Organize the data to plot in 3D. 
}
\usage{
  DataOrga(Proba=Proba,thres=thres,data=data, Expe.Var =Expe.Var,varcovlist=varcovlist)
}
\description{
It returns averaged data over the replicates.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{Proba }{Joint marginal probabilities (see paper) to plot features
}
   \item{thres }{Threshold used on joint marginal probabilities (see paper) to plot features
}
   \item{data }{Continuous  expression data of size \eqn{p\times n}, with \eqn{p} the number of features (e.g proteins or genes) and \eqn{n} the number of samples.  
}
  \item{Expe.Var }{Experimental explanatory (categorial) variable. It's the "confounding" variable. The patterns are compared for every modality of Expe.Var.}
\item{varcovlist }{Independent variable for each  modality of the experimental variable. This is obtained using function \emph{datNorm}.}  
}
\value{
 \item{DataAveraged }{Data averaged over the replicates.}
\item{names }{Feature names of the new data set.}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
