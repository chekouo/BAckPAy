/*BACkPAy:  C main function for identifying expression patterns
 * */

#include <time.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "utils.h"
#include "myfunction.h"
#include <R.h>
void  mainfunction(int *p1, int *n1,int *NbrGps1,double * ProtExp1, int * Time, int  * Resist,int *sample1,int *burnin1,double *atau, double *btau, double *c_h,double *b_beta,  double * rhoMean1, double * ProbProt1,double * probDiff) {
//Time is the independent variable
//resist is the experimental variable (confounding variable)
int i,h,l,j,j1,l1;
//clock_t t = clock();
double c_mu=100; 
double c_beta=b_beta[0];
int p=p1[0];
Rprintf("The number of features (i.e genes or proteins) is %d\n",p);
int n=n1[0];
Rprintf("The number of samples is %d\n",n);
int modal[n];int nbrmodal=n;
for (j=0;j<n;j++){
modal[j]=Time[j];
}
uniqvalues(modal,&nbrmodal);
Rprintf("The number of modalities of the indep. var. is %d\n",nbrmodal);
int idx[nbrmodal];
sort(nbrmodal,modal,idx);
int NbrCov=nbrmodal-1;
int H=pow(3,NbrCov); //nber of clusters
int NbrGps=NbrGps1[0];

Rprintf("The number of modalities of the Experimental var.  is %d\n",NbrGps);

double a=(n+NbrCov)/2.0;double b=1/pow(H,NbrCov/n);
long seed=100;
gsl_rng * r = gsl_rng_alloc (gsl_rng_rand48);
gsl_rng_set (r, seed);
//setvbuf(stdout, NULL, _IONBF, 0);
//These are the nine patterns
//int SignCoef[9][2]={{-1,1},{-1,0},{0,0},{1,0},{1,-1},{0,-1},{0,1},{1,1},{-1,-1}};
int ** SignCoef=malloc(H*sizeof(int *));
for (h=0;h<H;h++)
SignCoef[h]=malloc(NbrCov*sizeof(int));
int v[3]={1,0,-1};
Permutations(SignCoef,v, 3,NbrCov);
double ** ProtExp=dmatrix(0, p-1,0,n-1);

for (i=0;i<p;i++){
for (j=0;j<n;j++){
ProtExp[i][j]=ProtExp1[i*n+j];
}}
_Bool **TimeCov=bmatrix(0, n-1, 0, NbrCov-1);
for (j=0;j<n;j++){
for (l1=0;l1<NbrCov;l1++){
TimeCov[j][l1]=0;
int ss=0;
int te=0;
for (ss=l1+1;ss<nbrmodal;ss++){
if (Time[j]==modal[ss])
te+=1; 
}
if (te>0) TimeCov[j][l1]=1;
}
}


int *N=malloc(NbrGps*sizeof(int));
int ** idxRes=malloc(NbrGps*sizeof(int *));
for (l1=0;l1<NbrGps;l1++){
idxRes[l1]=malloc(n*sizeof(int));
N[l1]=0;
}
//Unique values of Resist
int nbrdist=n;
int * res=malloc(n*sizeof(int));
for (j=0;j<n;j++)
res[j]=Resist[j];
uniqvalues(res,&nbrdist); 
for (j=0;j<nbrdist;j++)
for (l1=0;l1<nbrdist;l1++){
j1=0;
for (j=0;j<n;j++){
if (Resist[j]==res[l1]){
N[l1]+=1; 
idxRes[l1][j1]=j;j1++;
} 
}
}
int *P=malloc(NbrGps*sizeof(int));
_Bool *** TimeCovRes=malloc(NbrGps*sizeof(_Bool**));
for (l1=0;l1<NbrGps;l1++){
TimeCovRes[l1]=bmatrix(0, N[l1]-1, 0, NbrCov-1);
P[l1]=p;
}
for (l1=0;l1<NbrGps;l1++){
for (j=0;j<N[l1];j++){
for (l=0;l<NbrCov;l++){
TimeCovRes[l1][j][l]=TimeCov[idxRes[l1][j]][l];
}
}
}

double ***ProtExpreResist=malloc(NbrGps*sizeof(double**));
for (l1=0;l1<NbrGps;l1++){
ProtExpreResist[l1]=malloc(P[l1]*sizeof(double*));
for (i=0;i<P[l1];i++){
ProtExpreResist[l1][i]=malloc(N[l1]*sizeof(double));
for (j=0;j<N[l1];j++){
ProtExpreResist[l1][i][j]=ProtExp[i][idxRes[l1][j]];
}
}
Center(P[l1],N[l1],ProtExpreResist[l1]);
}

int nbrprot[H];
int **nbrprotRes=malloc(H*sizeof(int*));
for (h=0;h<H;h++){
nbrprotRes[h]=malloc(NbrGps*sizeof(int));
}
int **clust=malloc(NbrGps*sizeof(int*));
for (l1=0;l1<NbrGps;l1++){
clust[l1]=  malloc(P[l1]*sizeof(int));
}
double * alpha=malloc(H*sizeof(double));
double * C=malloc(H*sizeof(double));
double ** beta=dmatrix(0,H-1,0, NbrCov-1);
double ** lower=dmatrix(0,H-1,0, NbrCov-1);
double ** betaMean=dmatrix(0,H-1,0, NbrCov-1);
double * mu=malloc(H*sizeof(double));
double * muMean=malloc(H*sizeof(double));
double * sigma2=malloc(H*sizeof(double));
double * sigma2Mean=malloc(H*sizeof(double));
double sumalpha=0;
double al=atau[0];
double bl=btau[0];
double lower1=al/bl;

_Bool *** rho=malloc(NbrGps*sizeof(_Bool **));
for (l1=0;l1<NbrGps;l1++){
rho[l1]=bmatrix(0, P[l1]-1, 0, H-1);
}
double probInit[H];
for (h=0;h<H;h++){
mu[h]=0;muMean[h]=0;sigma2Mean[h]=0;
C[h]=c_h[0];
sigma2[h]=1;
probInit[h]=1.0/H;
double alphah=20;
alpha[h]=alphah;
sumalpha+=alpha[h];
for (l=0;l<NbrCov;l++){
lower[h][l]=lower1*abs(SignCoef[h][l]);
beta[h][l]=lower[h][l]*abs(SignCoef[h][l]);
betaMean[h][l]=0;
}
}
for (l1=0;l1<NbrGps;l1++){
for (i=0;i<P[l1];i++){
 unsigned int rho1[H];
gsl_ran_multinomial (r, H,1,probInit,rho1);
for (h=0;h<H;h++){
rho[l1][i][h]=rho1[h];
}
}
}
for (h=0;h<H;h++)
nbrprot[h]=0;


int comb=pow(H,NbrGps);
double ** ProbProt=malloc(comb*sizeof(double*));
int h1;
for (h1=0;h1<comb;h1++){
ProbProt[h1]=malloc(p*sizeof(double));
for (i=0;i<p;i++)
ProbProt[h1][i]=0;
}
for (i=0;i<p;i++){
probDiff[i]=0;
}
double *** rhoMean=malloc(NbrGps*sizeof(double**));
for (l1=0;l1<NbrGps;l1++){
rhoMean[l1]=dmatrix(0,P[l1]-1,0,H-1);
for (h=0;h<H;h++){
nbrprotRes[h][l1]=0;
for (i=0;i<P[l1];i++){
rhoMean[l1][i][h]=0;
nbrprot[h]+=rho[l1][i][h];
nbrprotRes[h][l1]+=rho[l1][i][h];
if (rho[l1][i][h]==1)
clust[l1][i]=h;
}
}
}

int sample=sample1[0]; 
int burnin=burnin1[0];
int s;
double acceptLower=0;
double *Logpost=malloc(sample*sizeof(double));
double *** beta_sample=malloc((sample-burnin)*sizeof(double **));
for (s=0;s<sample-burnin;s++){
beta_sample[s]=malloc(H*sizeof(double*));
for (h=0;h<H;h++){
beta_sample[s][h]=malloc(NbrCov*sizeof(double));
for (l=0;l<NbrCov;l++)
beta_sample[s][h][l]=beta[h][l];
}
}
int h_old;
for (s=0;s<sample;s++){
for (l1=0;l1<NbrGps;l1++){
for (i=0;i<P[l1];i++){
h_old=clust[l1][i];
int h_new=Rho(h_old,ProtExpreResist[l1][i], TimeCovRes[l1], beta,mu,sigma2,alpha,sumalpha,C,rho[l1][i], nbrprot,p,N[l1],NbrCov,H,SignCoef,r);
nbrprot[h_old]=nbrprot[h_old]-1;
nbrprot[h_new]=nbrprot[h_new]+1;
clust[l1][i]=h_new;
nbrprotRes[h_old][l1]=nbrprotRes[h_old][l1]-1;
nbrprotRes[h_new][l1]=nbrprotRes[h_new][l1]+1;
if (s>=burnin){
rhoMean[l1][i][h_new]+=1;
}
}
}
if (s>=burnin){
for (i=0;i<p;i++){
int ti=1;
for (l1=0;l1<NbrGps;l1++){ 
ti*=rho[l1][i][(H+1)/2-1];
}
if (ti==0) probDiff[i]+=1;

for (h1=0;h1<comb;h1++){
int prd=1;
for (l1=0;l1<NbrGps;l1++){
int bj;
bj=((int)floor(h1*pow(H,-l1)))%H;
prd*=rho[l1][i][bj];
}
if (prd==1)
ProbProt[h1][i]+=1;
}
}
}
for (h=0;h<H;h++){

mu[h]=Mu(ProtExpreResist,TimeCovRes,h, beta[h],sigma2[h],C[h],c_mu,rho, nbrprotRes[h],P,N,NbrCov,NbrGps,SignCoef[h],r);
sigma2[h]=Sigma2(ProtExpreResist,TimeCovRes, h, beta[h],mu[h],C[h],a,b,rho,nbrprotRes[h],P,N,NbrCov,NbrGps,SignCoef[h],r);

for (l=0;l<NbrCov;l++){
if (SignCoef[h][l]!=0){
if (nbrprot[h]>=2){
lower[h][l]=Lower1(lower[h][l], beta[h][l],sigma2[h],c_beta, al,bl,&acceptLower,r);
} else {
lower[h][l]=gsl_ran_gamma (r, al, 1/bl);
}
beta[h][l]=Beta(ProtExpreResist,TimeCovRes,h,l, beta[h],sigma2[h], mu[h],C[h], c_beta, rho, nbrprotRes[h], P, N,NbrCov,NbrGps, lower[h][l],SignCoef[h], r);
if (s>=burnin){
beta_sample[s-burnin][h][l]=beta[h][l];
}
}
}
if (s>=burnin){
muMean[h]+=mu[h];
sigma2Mean[h]+=sigma2[h];
for (l=0;l<NbrCov;l++){
if (SignCoef[h][l]!=0){
betaMean[h][l]+=beta[h][l];
}}
}
}
Logpost[s]=logpost(ProtExpreResist, TimeCovRes, beta,mu, sigma2,alpha,sumalpha,C,rho, nbrprot,P,N,NbrCov,NbrGps,H, SignCoef,a,b, al,bl,c_mu, c_beta, lower);
if (s%(sample/5)==1){
Rprintf("The number of MCMC iterations is %d\n",s);
Rprintf("\n\n");
R_CheckUserInterrupt();
} // end if 

}

//process the results
for (h=0;h<H;h++){
for (l=0;l<NbrCov;l++){
betaMean[h][l]=betaMean[h][l]/(sample-burnin);
}
sigma2Mean[h]=sigma2Mean[h]/(sample-burnin);
muMean[h]=muMean[h]/(sample-burnin);
for (l1=0;l1<NbrGps;l1++){
for (i=0;i<p;i++){
rhoMean[l1][i][h]=rhoMean[l1][i][h]/(sample-burnin);
}
}
}

for (i=0;i<p;i++){
probDiff[i]=probDiff[i]/(sample-burnin);
for (h1=0;h1<comb;h1++){
ProbProt[h1][i]=ProbProt[h1][i]/(sample-burnin);
}
}

free_dmatrix(ProtExp,0, p-1,0,n-1);
free_dmatrix(beta,0,H-1, 0,NbrCov-1);
free(mu);free(sigma2);free(alpha);
free_dmatrix(betaMean,0,H-1, 0,NbrCov-1);
free_dmatrix(lower,0,H-1, 0,NbrCov-1);

for (h=0;h<comb;h++){
for (j=0;j<p;j++){
ProbProt1[h * p + j]=ProbProt[h][j];
}
}
for (l1=0;l1<NbrGps;l1++){
for (i=0;i<p;i++){
for (h=0;h<H;h++){
rhoMean1[l1*p*H+i*H+h]=rhoMean[l1][i][h];
}
}
}



free_dmatrix(ProbProt,0,comb-1,0,p-1);

for (l1=0;l1<NbrGps;l1++){
free(idxRes[l1]);
free_bmatrix(TimeCovRes[l1],0, N[l1]-1, 0, NbrCov-1);
free_bmatrix(rho[l1],0, P[l1]-1, 0, H-1);
free_dmatrix(rhoMean[l1],0, P[l1]-1, 0, H-1);
free_dmatrix(ProtExpreResist[l1],0, P[l1]-1,0,N[l1]-1);
free(clust[l1]);
}
free(ProtExpreResist);
gsl_rng_free (r);free(muMean);free(sigma2Mean);
free(C);
free(rhoMean);
free(rho);
free(clust);
free(Logpost);
for (s=0;s<sample-burnin;s++)
free_dmatrix(beta_sample[s],0,H-1, 0,NbrCov-1);
free(beta_sample);
for (h=0;h<H;h++){
free(nbrprotRes[h]);free(SignCoef[h]);}
free(TimeCovRes);
free(nbrprotRes);
free(idxRes);
free(SignCoef);free(P);
free(N);free(res);
//t = clock() - t;
//double  time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
  //  printf("\n\nTime taken in seconds is %f\n",time_taken);
   // printf("\nTime taken in minutes is %f\n",time_taken/60);
}
