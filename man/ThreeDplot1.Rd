\name{ThreeDplot1}
\alias{ThreeDplot1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function: Preparing to represent data in 3D
}
\description{
Intermediate function to represent patterns.
}
\usage{
ThreeDplot1(meanF=meanF,nbrMark=nbrMark,NbrGps=NbrGps,NbrModCov=NbrModCov,
LabelsLegend=LabelsLegend,LabelConf=LabelConf,patternName=patternName,miny=miny,
maxy=maxy,thres=thres,maxprob=maxprob)
}
\arguments{
  \item{meanF }{Average expression data over the number of replicates for each modality of the experimental variable.}
\item{nbrMark }{Number of features to plot}
\item{NbrGps }{This is the number of modalities of the experimental variable (\eqn{T=}NbrGps)}
  \item{NbrModCov }{This is the number of modalities of the independent variable (\eqn{q=}NbrModCov)}
  \item{LabelsLegend }{Label names of the independent variable}
  \item{LabelConf }{Label names of the experimental  variable}
  \item{patternName }{Pattern name to plot}
  \item{miny }{Minimium value of \emph{meanF}}
  \item{maxy }{Maximum value of \emph{meanF}}
   \item{thres }{Threshold used on joint marginal probabilities (see reference) to plot features}
  \item{maxprob }{The maximum of  joint marginal probabilities (see reference) for pattern name \emph{patternName}.}
}

%\value{
%}
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
