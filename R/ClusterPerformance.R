
## Function to compute clustering performance for BACkPAY

ClusterPerformance=function(BackPay=BackPay,Generate=Generate,thres=seq(0,1,0.02)){
#  library(flux)
  rhoK=Generate$RhoKnown
  Dim=dim(rhoK)
  TT=Dim[1];p=Dim[2];H=Dim[3]
  rhoMean=NULL;rho=NULL;
  for (k in 1:TT){
    rhoMean=rbind(rhoMean, BackPay$rhoMean[k,,])
    rho=rbind(rho,rhoK[k,,])
  }
  #print(head(rhoMean))
  #rhoEst=apply(rhoMean,1,function (v) (v>=max(v))*1)
  #print(rhoEst%*%as.matrix(rho))
  #print(ari(rhoEst%*%as.matrix(rho)))
  mceTh=class=nclass=mc=rep(0,length(thres))
  for (j in 1:length(thres)){
    th=thres[j];
    rhoEst=apply(rhoMean,1,function (v) (v>=max(th,max(v)))*1)
    idx=which(apply(rhoEst,2,sum)>=2)
    if (length(idx)!=0){
    for (l in 1:length(idx)){
    idx1=which(rhoEst[,idx[l]]==1);
    rhoEst[idx1,idx[l]]=((1:length(idx1))==1)*1;
    }
  }
    XX=rhoEst%*%as.matrix(rho)
    sumX=sum(XX)
    nclass[j]=TT*p-sumX
    diag(XX)=NA
    mc[j]=sum(XX,na.rm=T)
    class[j]=sumX- mc[j]
  }
  A1=nclass/(TT*p);A2=(nclass+mc)/(TT*p);
  Anc=round(flux::auc(thres,A1),digits = 3);
  Awc=round(1-flux::auc(thres,A2),digits = 3);
  Amc=round(1-(Anc+Awc),digits = 3)
  list(Anc=Anc,Amc=Amc,Awc=Awc)
}






