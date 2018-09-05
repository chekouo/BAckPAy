#datNorm=function(yy,confound.var,cova){
 
  datNorm=function(data=data,Expe.Var=NULL,Ind.Var=Ind.Var){
    yy=data;cova=Ind.Var;confound.var=Expe.Var;
  nbrgrp=length(unique(confound.var))
  covlist=list();
  y=matrix(0,nrow(yy),length(Ind.Var));
  if (nbrgrp!=0){
 
  for (l in 1:nbrgrp){
    vv=sort(unique(confound.var))[l];
    XX=yy[,confound.var==vv];
    #XX=log2(XX)
    XX=(XX-apply(XX,1,mean));
    y[,confound.var==vv]=as.matrix(XX);
    covlist[[l]]=cova[confound.var==vv]
  }
  rownames(y)=rownames(yy)
  } else {
    covlist[[1]]=rep(1,length(Expe.Var))
    y=(y-apply(y,1,mean));
  }
  return(list("y"=y,"covlist"=covlist))
}