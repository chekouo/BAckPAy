#DataOrga=function(ProbProtein,thres,y, vargroup, varcovlist){
  DataOrga=function(Proba=Proba,thres=thres,data=data, Expe.Var =Expe.Var, varcovlist=varcovlist){
    ProbProtein=Proba; vargroup=Expe.Var;y=data;
  Nam=NULL;
  listp=(ProbProtein>=thres)
  wlist=which(listp>0)
  if (is.null(rownames(y))){
    rownames(y)=1:nrow(y) 
  }
  names(listp)=rownames(y);
  summ=sum(listp)
  uniqgrp=sort(unique(vargroup));
  n=length(uniqgrp);
  agreg=NULL
  if (summ>0){
    for (l in 1:n){
      y1=as.matrix(y[wlist,vargroup==uniqgrp[l]])
      Meany1=matrix(0,summ,length(unique(varcovlist[[l]])))
      for (i in 1:summ){
        if (summ>1){
          yy=as.vector(as.matrix(y1[i,]));
        } else if (summ==1) {yy=y1;}
        Meany1[i,]=as.vector(as.matrix(aggregate(yy, by=list(varcovlist[[l]]), FUN=mean)[2]))
      }
      agreg=rbind(agreg,Meany1);
    }
    Nam=names(which(listp>0));
  }
  return(list("DataAveraged"=agreg,"names"=Nam))
}