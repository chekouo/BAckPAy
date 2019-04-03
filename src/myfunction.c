#include <stddef.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdbool.h>
#include <math.h>
#include <Rmath.h>

#include "utils.h"
#include "myfunction.h"
int Rho(int h_old,double * data, _Bool ** Time, double **beta,double *mu, double *sigma2,double *alpha,double sumalpha,double *C,_Bool* rho, int *nbrprot,int p,int n,int nbrcov,int H, int **signCoef, const gsl_rng * r){
// This function sample rho's, cluster memberships....
int h;
// p is the number of proteins or features, i is the protein index
// n is the size of covariates, H is the number of clusters
// nbrcov is the number of covariates
double* logpl=malloc(H*sizeof(double));
int l,j; double pll=0;
for (h=0;h<H;h++){
 double sumy2=0;double sum2y=0;
double yesti;
for (j=0;j<n;j++){
yesti=0;
for (l=0;l<nbrcov;l++){
yesti+=(beta[h][l]*signCoef[h][l]*Time[j][l]);
}
yesti+=mu[h];
sumy2+=pow(data[j]-yesti,2);
sum2y+=(data[j]-yesti);
}
sum2y=pow(sum2y,2);
if (rho[h]==1){
pll=nbrprot[h]-1;
} else {pll=nbrprot[h];}
logpl[h]=log(pll+alpha[h])-log(p-1+sumalpha)-(n/2.0)*log(sigma2[h])-0.5*log(1+n*C[h])-(1/(2*sigma2[h]))*(sumy2-(C[h]/(1+n*C[h]))*sum2y);
}

double maxlogpl=max(H,logpl);
double  prob[H];double sumprob=0;
  for (h=0;h<H;h++){
    prob[h]=exp(logpl[h]-maxlogpl);
    sumprob+=prob[h];
}
free(logpl);
unsigned int rho1[H];
gsl_ran_multinomial (r, H,1,prob,rho1);
// multinomial(H,1,prob,rho1);
int h1=h_old;
for (h=0;h<H;h++){
rho[h]=rho1[h];
if (rho[h]==1)
h1=h;
}
rho[h1]=1;
return h1;
}

double Mu(double *** data,_Bool *** Time,int h, double *beta,double sigma2,double C,double c_mu,_Bool*** rho, int *nbrprot,int *P,int* N,int nbrcov, int nbrgrp, int * signCoef,const gsl_rng * r){
// n is the size of covariates,
int re;  
double invsig=0;
for (re=0;re<nbrgrp;re++){
invsig+=nbrprot[re]*N[re]/(sigma2*(N[re]*C+1));
}


double  sigma2mu=1/((1/c_mu)+invsig);
int i,j,l;double yy,yestj;
double muh=0;
for (re=0;re<nbrgrp;re++){
for (j=0;j<N[re];j++){
yy=0;
yestj=0;
for (l=0;l<nbrcov;l++){
yestj+=(beta[l]*signCoef[l]*Time[re][j][l]);
}
for (i=0;i<P[re];i++){
if (rho[re][i][h]==1){
yy+=data[re][i][j]-yestj;
}
}
muh+=(1-N[re]*(C/(N[re]*C+1)))*yy;
}
}
muh=sigma2mu*muh/sigma2;

return muh+gsl_ran_gaussian (r, sqrt(sigma2mu)); 
//return muh+rnorm(0,sqrt(sigma2mu));

}


double Sigma2(double *** data,_Bool *** Time, int h, double * beta,double mu,double C,double a,double b,_Bool*** rho, int* nbrprot,int *P,int* N,int nbrcov,int nbrgrp, int* signCoef,const gsl_rng * r){
double asig=a;
int i,j,l,re; 
for (re=0;re<nbrgrp;re++){
asig+=N[re]*nbrprot[re]/2.0;
}
double sum2yi,sumyi2,yestj;
double bsig=0;
for (re=0;re<nbrgrp;re++){
for (i=0;i<P[re];i++){
if (rho[re][i][h]==1){
sum2yi=0; sumyi2=0;
for (j=0;j<N[re];j++){
yestj=0;
for (l=0;l<nbrcov;l++){
yestj+=(beta[l]*signCoef[l]*Time[re][j][l]);
}
yestj+=mu;
double err=data[re][i][j]-yestj;
sumyi2+=pow(err,2);sum2yi+=err;
}
sum2yi=pow(sum2yi,2);
bsig+=sumyi2-(C/(N[re]*C+1))*sum2yi;
}
}
}
bsig=b+0.5*bsig;

return 1/gsl_ran_gamma (r, asig, 1/bsig);
//return 1/rgamma (asig, 1/bsig);
}

double Beta(double *** data,_Bool *** Time,int h,int l, double *beta,double sigma2,
            double mu,double C,double c_beta,_Bool*** rho,
            int *nbrprot,int *P,int* N,int nbrcov, int nbrgrp,double lower,int *signCoef,const gsl_rng * r){
    //P is the number of proteins per cell type (sensitive and resistant)
    //nbrprot number of proteins in group per cell type (sensitive and resistant)
    //n is the size of covariates,
    
    int i,l1;
    double tim=0;double t2=0;
    int j,re;
    for (re=0;re<nbrgrp;re++){
        double t22=0;
        for (j=0;j<N[re];j++){
            tim+=(nbrprot[re])*pow(signCoef[l]*Time[re][j][l],2);
            t22+=sqrt(nbrprot[re])*sqrt(C/(N[re]*C+1))*signCoef[l]*Time[re][j][l];
        }
        t2+=pow(t22,2);
    }
    tim=tim-t2;
    double  sigma2beta=sigma2/(1/c_beta+tim);
    double yy,yestj;
    double timy=0;
    double sumyx=0;
    for (re=0;re<nbrgrp;re++){
        double sumy=0;double sumt=0;
        for (j=0;j<N[re];j++){
            yy=0;yestj=0;
            for (l1=0;l1<nbrcov;l1++){
                if (l1!=l){
                    yestj+=(beta[l1]*signCoef[l1]*Time[re][j][l1]);
                }
            }
            yestj+=mu;
            for (i=0;i<P[re];i++){
                if (rho[re][i][h]==1){
                    yy+=data[re][i][j]-yestj;
                }
            }
            timy+=signCoef[l]*Time[re][j][l]*yy;
            sumt+=(C/(N[re]*C+1))*signCoef[l]*Time[re][j][l];
            sumy+=yy;
        }
        sumyx+=sumt*sumy;
        
    }
    double mubeta=sigma2beta*(lower/c_beta+(timy-sumyx))/sigma2;
   // printf("\nLower=%lf= ",lower);
    //printf("Sigma2=%lf= ",sigma2);
    //printf("Prod=%lf=",timy-sumyx);
    //printf("SigBeta=%lf= ",sqrt(sigma2beta));
    //printf("muBeta=%lf= ",mubeta);
    // printf("H=%d= \n",h);
    //return Truncate(mubeta,sqrt(sigma2beta), lower, r);
    
    return r_lefttruncnorm(lower, mubeta,sqrt(sigma2beta),r);
    
}



double Beta1(double *** data,_Bool *** Time,int h,int l, double *beta,double sigma2,
            double mu,double C,double c_beta,_Bool*** rho, 
            int *nbrprot,int *P,int* N,int nbrcov, int nbrgrp,double lower,int *signCoef,const gsl_rng * r){
//P is the number of proteins per cell type (sensitive and resistant)
//nbrprot number of proteins in group per cell type (sensitive and resistant)
//n is the size of covariates,  

int i,l1;
double tim=0;double t2=0;
int j,re;
for (re=0;re<nbrgrp;re++){
double t22=0;
for (j=0;j<N[re];j++){
tim+=(nbrprot[re])*pow(signCoef[l]*Time[re][j][l],2);
t22+=sqrt(nbrprot[re])*sqrt(C/(N[re]*C+1))*signCoef[l]*Time[re][j][l];
}
t2+=pow(t22,2);
}
tim=tim-t2;
double  sigma2beta=sigma2/(1/c_beta+tim);
double yy,yestj;
double timy=0;
double sumt=0;double sumy=0;
for (re=0;re<nbrgrp;re++){
for (j=0;j<N[re];j++){
yy=0;yestj=0;
for (l1=0;l1<nbrcov;l1++){
if (l1!=l){
yestj+=(beta[l1]*signCoef[l1]*Time[re][j][l1]);
}
}
yestj+=mu;
for (i=0;i<P[re];i++){
if (rho[re][i][h]==1){
yy+=data[re][i][j]-yestj;
}
}
timy+=signCoef[l]*Time[re][j][l]*yy;
sumt+=(C/(N[re]*C+1))*signCoef[l]*Time[re][j][l];
sumy+=yy;
}
}
double mubeta=sigma2beta*(lower/c_beta+(timy-sumt*sumy))/sigma2;
    printf("\nLower=%lf= ",lower);
    printf("Sigma2=%lf= ",sigma2);
    printf("Prod=%lf=",timy-sumt*sumy);
printf("SigBeta=%lf= ",sqrt(sigma2beta));
    printf("muBeta=%lf= \n",mubeta);
//return Truncate(mubeta,sqrt(sigma2beta), lower, r);

return r_lefttruncnorm(lower, mubeta,sqrt(sigma2beta),r);

}


double logpost(double *** data, _Bool *** Time, double **beta,double *mu, double *sigma2,double *alpha,double sumalpha,double *C,_Bool*** rho, int *nbrprot,int* P,int *N,int nbrcov,int nbrgrp,int H, int **signCoef, double a, double b, double al, double bl,double c_mu, double c_beta, double** truncbeta){
double logpost=0;
int i,l,h,j,re; 
for (h=0;h<H;h++){
double xx,x;
if (nbrprot[h]>1){
xx=x=0;double loggamma=0;
loggamma= gsl_sf_lngamma (nbrprot[h]+alpha[h]);
 //Rprintf("LogGamma1 is %lf\n",loggamma);
//loggamma=lgammafn(nbrprot[h]+alpha[h]);
 // Rprintf("LogGamma2 is %lf\n",loggamma);
for (re=0;re<nbrgrp;re++){
for (i=0;i<P[re];i++){
if (rho[re][i][h]==1){
double sumy2=0;double sum2y=0;
double yesti;
for (j=0;j<N[re];j++){
yesti=0;
for (l=0;l<nbrcov;l++){
yesti+=(beta[h][l]*signCoef[h][l]*Time[re][j][l]);
}
yesti+=mu[h];
sumy2+=pow(data[re][i][j]-yesti,2);
sum2y+=(data[re][i][j]-yesti);
}
sum2y=pow(sum2y,2);
xx+=sumy2-(C[h]/(N[re]*C[h]+1))*sum2y;
x+=-0.5*log(1+N[re]*C[h])+(-N[re]/2.0)*(log(sigma2[h]));
}
}
}
double lltrunc=0;
double llbeta=0;
for (l=0;l<nbrcov;l++){
if (abs(signCoef[h][l])>0){
llbeta+=-0.5*log(sigma2[h])-0.5*pow(beta[h][l]-truncbeta[h][l],2)/(pow(c_beta,2)*sigma2[h]);
lltrunc+=(al-1)*log(truncbeta[h][l])-bl*truncbeta[h][l];
}
}
double llsig=0; 
llsig=(-a-1)*log(sigma2[h])-b/sigma2[h];
double llmu=0;
llmu=-0.5*(pow(mu[h],2))/pow(c_mu,2);
logpost+=llbeta+lltrunc+llmu+llsig+loggamma+x-(0.5/sigma2[h])*xx;
}
}
return logpost;
}

double Lower1(double  lower, double beta,double  sigma2, double c_beta,double a,double b,double * accept,const gsl_rng * r){
//double var=0.01;
double var=0.1;
double betaold=lower/var+0.0001;
double alphaold=betaold*lower+0.0001;
//printf("Betaold=%lf\n",betaold);
//printf("Alphaold=%lf\n",alphaold);
//double lowerprop=MAX(gsl_ran_gamma (r, alphaold, 1/betaold),gsl_ran_flat (r, 0, 0.05));
double lowerprop=gsl_ran_gamma (r, alphaold, 1/betaold)+gsl_ran_flat (r, 0, 0.05);
      // printf("lowerprop=%lf ",lowerprop);
//double lowerprop=rgamma (alphaold, 1/betaold);


//printf("Prop=%lf\n",lowerprop);
double betanew=lowerprop/var;
double alphanew=lowerprop*betanew;

double logliknew=0;double loglikold=0;
    //loglikold=-0.5*pow(beta-lower,2)/(pow(c_beta,2)*sigma2);
    //logliknew=-0.5*pow(beta-lowerprop,2)/(pow(c_beta,2)*sigma2);
loglikold=-0.5*pow(beta-lower,2)/(c_beta*sigma2);
logliknew=-0.5*pow(beta-lowerprop,2)/(c_beta*sigma2);

//double difflogprop=log(gsl_ran_gamma_pdf(lower,alphanew,1/betanew))-log(gsl_ran_gamma_pdf(lowerprop,alphaold,1/betaold));
//double difflogprop=dgamma(lower,alphanew,1/betanew,1)-dgamma(lowerprop,alphaold,1/betaold,1);
double difflogprop=alphanew*log(betanew)-gsl_sf_lngamma (alphanew)+
  (alphanew-1)*log(lower)-betanew*lower-((alphaold-1)*log(lowerprop)-betaold*lowerprop+alphaold*log(betaold)-gsl_sf_lngamma (alphaold));
//x^{a-1} e^{-x/b} dx
//printf("Difflog1=%lf\n",difflogprop);
double diff=difflogprop+logliknew-loglikold+(a-1)*(log(lowerprop)-log(lower))-b*(lowerprop-lower);
//printf("LogLower=%lf\n",log(lowerprop));
double  uni=gsl_ran_flat (r, 0, 1);
//double  uni=runif(0,1);
if (log(uni)<diff){
lower=lowerprop;
*accept+=1;
}
   // printf("Bet=%lf ",beta);
    //printf("Lower=%lf\n",lower);
return lower;
}







