#' @export
Nameclust <-function(NbrModCov,NbrGps){
  nbrcov=NbrModCov-1; nbrgrp=NbrGps;
  H=3^nbrcov;
signClus=clust(c(1,0,-1),nbrcov)
#signClus[,1]=-signClus[,1];
Indx=matrix(1,H,nbrcov)
Indx=replace(signClus,signClus==-1 ,3)
Indx=replace(Indx,signClus==0 ,2)

labb=rep(0,H)
VV=c("Up","Flat","Down")
labb=apply(Indx,1,function(x) paste(VV[x],collapse=""))
nbrc=H^nbrgrp
nameclust=rep("0",nbrc);
for (xx in 1:(nbrc)){
  nameclust[xx]=paste(labb[(1+floor((xx-1)*H^-c(0:(nbrgrp-1)))%%H)],collapse="-");
}
list(namecl=labb,namegroup=nameclust)
}
