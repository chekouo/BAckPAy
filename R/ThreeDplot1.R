  ThreeDplot1=function(meanF=meanF,nbrMark=nbrMark,NbrGps=NbrGps,NbrModCov=NbrModCov, LabelsLegend=LabelsLegend,LabelConf=LabelConf,
                       patternName=patternName,miny=miny,maxy=maxy,thres=thres,maxprob=maxprob){
    Nbrgrp=NbrGps;k=NbrModCov;group=patternName;th=thres; 
  y=as.vector(t(meanF)) #Matrix nbrMark*Nbrgroup times k
  x=rep(0:(k-1),nbrMark*Nbrgrp);
  z=NULL;
  for (l in 1:Nbrgrp){
    z=c(z,rep(l,nbrMark*k));
  }
  label1=LabelConf[1];
  for (l in 2:k){
    label1=c(label1," ", LabelConf[l]);
  }
  label2=LabelsLegend[1];
  for (l in 2:k){
    label2=c(label2," ", LabelsLegend[l]);
  }
  if (length(LabelsLegend)==2){
    label2=rep("",6)
    label2[1]=LabelsLegend[1];label2[6]=LabelsLegend[2]
  }
  YLAB="";
  s3d=scatterplot3d(x,z,y,pch=20,color=z+1, xlab="",zlab="Expression",ylab="",yaxt="n",
                    x.ticklabs = label2,y.ticklabs = YLAB,cex.lab = 1,cex.axis = 1.1,
                    box = T,#zlim=c(min(y,-1.5),max(y,1.5)),
                    zlim=c(miny,maxy),
                    main=paste(group, ", pb=",th, ", pma=",round(maxprob, digits = 2), ", N=",nbrMark,sep=""))
  ngg=k*nbrMark
  
  ng=0;
  coll=2;
  while (ng<k*nbrMark*Nbrgrp){
    k1=0;
    while(k1<k*nbrMark){
      for (l in 1:(k-1)){
        s3d$points3d(x[(k1+l+ng):(k1+l+ng+1)],z[(k1+l+ng):(k1+l+ng+1)],
                     y[(k1+l+ng):(k1+l+ng+1)],type="l",col=coll)
      }
      k1=k1+k;
    }
    ng=ng+k*nbrMark;
    coll=coll+1;
  }
  legend("topleft",lwd=3,bty = "n",cex=0.7,
         legend=LabelConf,col = 2:(length(LabelConf)+1), pch=20, horiz = F)
  #return (s3d)
}