#PlotThreeD=function(data=data,expe.var,ind.var,ProbRes,patterntype,thres){
  PlotThreeD=function(data=data,Ind.Var = Ind.Var, Expe.Var = NULL,ProbProtGrp=ProbProtGrp,patternname=patternname,thres=thres){
    expe.var= Expe.Var;ind.var=Ind.Var;ProbRes=ProbProtGrp;patterntype=patternname;
    data=as.matrix(data);
    if (is.null(rownames(data))){
      rownames(data)=1:nrow(data);
    }
    if (is.null(colnames(data))){
      colnames(data)=1:ncol(data);
    }
  Pb=ProbRes[patterntype,]
  listP=(Pb>=thres)
  if (sum(listP)>0){
    Probsort=sort(Pb[listP],decreasing = T);
  } else {Probsort=NULL
  stop(paste("Note: We do not have any features in the pattern",patterntype,
             "with the threshold of", thres))
  }
  if (is.null(Expe.Var)){
    expe.var=rep(1,length(Ind.Var))
  }
  ord=order(expe.var,ind.var)
  confound.var=expe.var[ord];
  cova=ind.var[ord];
  yy=data[,ord]
  datNN=datNorm(data=yy,Expe.Var=confound.var,Ind.Var=cova)
  #datNorm=function(data=data,Expe.Var=NULL,Ind.Var=Ind.Var)
  nbrgrp=length(unique(confound.var))
  y=datNN$y; 
  #y=yy;
  covlist=datNN$covlist
  miny=min(y);maxy=max(y);
  LIST=DataOrga(Proba=ProbRes[patterntype,],thres=thres,data=y, Expe.Var=confound.var, varcovlist=covlist);
  MeanyF=LIST$DataAveraged;
  #print(LIST)
  listMarkers=list();
  ll=0;
  if (is.null(MeanyF)==F){
    ll=ll+1;
    listMarkers[[ll]]=LIST$names;
    coll=NULL
    nbrcov=length(unique(cova))-1;
    for (l1 in 1:nbrgrp){
      coll=c(coll,rep(rainbow(20)[2*l1],length(listMarkers[[ll]])));
    }
    lle=length(LIST$names)
    LabelConf=as.character(sort(unique(confound.var)))
    Labels=sort(as.character(unique(ind.var)))
    #requireNamespace(scatterplot3d)
    ThreeDplot1(meanF=MeanyF,nbrMark=lle,NbrGps=nbrgrp,NbrModCov=nbrcov+1,LabelsLegend=Labels,LabelConf=LabelConf,
               patternName=patterntype,miny=miny,maxy=maxy,thres=thres,maxprob=Probsort[1]);
    #ThreeDplot1=function(meanF=meanF,nbrMark=nbrMark,NbrGps=NbrGps,NbrModCov=NbrModCov, LabelsLegend=LabelsLegend,LabelConf=LabelConf,
                        # patternName=patternName,miny=miny,maxy=maxy,thres,maxprob=maxprob)
    lle1=min(lle,36)
  } 
  return (list(ProbSort=Probsort, dataplot=y[LIST$names,]))
  #return (list(prob=sort(ProbRes[patterntype,][ProbRes[patterntype,]>0.5],decreasing = T),plott=s3d))
}
