
\docType{data}
\name{Isoform}
\alias{Isoform}
\title{Data sets for the Isoform protein data}
\format{A list of 10 data frames:
\describe{
  \item{SampleInfo}{Sample, cell line and  time point  information}
  \item{NuclearHL60}{Nuclear proteins for cell line HL60.}
  \item{NuclearTF1}{Nuclear proteins for cell line TF1.}
  \item{Nuclearu937}{Nuclear proteins for cell line u937.}
  \item{SurfaceHL60}{Surface proteins for cell line HL60.}
  \item{SurfaceTF1}{Surface proteins for cell line TF1}
  \item{Surfaceu937}{Surface proteins for cell line u937.}
  \item{TCEHL60}{TCE proteins for cell line HL60.}
  \item{TCETF1}{TCE proteins for cell line TF1.}
  \item{TCEu937}{TCE proteins for cell line u937.}
}}

\usage{
Isoform
}
\description{
A list of 10 data frames that contains the isoform protein expressions and sample information 
}
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
}
\references{
Thierry Chekouo et al (2019), \emph{Investigating Protein Patterns in Human Leukemia Cell Line Experiments: A Bayesian Approach for Extremely Small Sample Sizes}, \emph{submitted}.
}
\keyword{datasets}
