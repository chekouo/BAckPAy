#' @export
clust <-
function(VV=c(1,0,-1),nbrcov=nbrcov){
  p=nbrcov;n=length(VV);k=n^p;
  clustpattern=matrix(0,k,p)
  for (i in c(0:(k-1))){
    b=1+floor(i*n^-c(0:(p-1)))%%n
    clustpattern[i+1,]=VV[b];
  }
  clustpattern
}
