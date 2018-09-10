#include <stddef.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdbool.h>
#include <math.h>
#include "utils.h"
#include "myfunction.h"
#include <R.h>
#include <Rmath.h>

void generedata(int *nbrTimeMod, int *nbrResistMod1,int *TimePoint,int * Resist, int *p1, double * yv,int *nbrsampl1,double* beta,int *seed,int * Rho,double* sig){
 GetRNGstate();
int l,h,i,j,l1;
int p=p1[0];
int nbrsampl=nbrsampl1[0];
int nbrResistMod=nbrResistMod1[0];
int n=nbrResistMod*nbrsampl*nbrTimeMod[0];
double ** y=dmatrix(0, p-1,0,n-1);
gsl_rng * r = gsl_rng_alloc (gsl_rng_rand48);
long seedd=*seed;
gsl_rng_set (r, seedd);
l=0;
for (i=0;i<nbrResistMod*nbrsampl;i++){
for (j=0;j<nbrTimeMod[0];j++){
TimePoint[l]=j;
l++;
}
}
l=0;
for (j=0;j<nbrResistMod;j++){
for (i=0;i<nbrTimeMod[0]*nbrsampl;i++){
Resist[l]=j;
l++;
}
}
int NbrCov=nbrTimeMod[0]-1; //length(Ts)-1
int NbrGps=nbrResistMod; //length(R)
int Ts[nbrTimeMod[0]];
for (j=0;j<nbrTimeMod[0];j++){
Ts[j]=j;
}
_Bool **TimeCov=bmatrix(0, n-1, 0, NbrCov-1);

ConvertToBinaryMatP(TimePoint,Ts,n,nbrTimeMod[0],TimeCov);

_Bool *** TimeCovRes=malloc(NbrGps*sizeof(_Bool**));

int ** idxRes=malloc(NbrGps*sizeof(int *));
for (l1=0;l1<NbrGps;l1++){
idxRes[l1]=malloc(nbrTimeMod[0]*nbrsampl*sizeof(int));
}
int j1;
for (l1=0;l1<NbrGps;l1++){
j1=0;
for (j=0;j<n;j++){
if (Resist[j]==l1){
idxRes[l1][j1]=j;
j1++;
} 
}
}
for (l1=0;l1<NbrGps;l1++){
TimeCovRes[l1]=bmatrix(0, nbrTimeMod[0]*nbrsampl-1, 0, NbrCov-1);
}
for (l1=0;l1<NbrGps;l1++){
for (j=0;j<nbrTimeMod[0]*nbrsampl;j++){
for (l=0;l<NbrCov;l++){
TimeCovRes[l1][j][l]=TimeCov[idxRes[l1][j]][l];
}
}
}


int H=pow(3,NbrCov);
double probInit[H];
int ** SignCoef=malloc(H*sizeof(int *));
for (h=0;h<H;h++){
SignCoef[h]=malloc(NbrCov*sizeof(int));
probInit[h]=1.0/H;
}
int v[3]={1,0,-1};
Permutations(SignCoef,v, 3,NbrCov);


_Bool *** rho=malloc(NbrGps*sizeof(_Bool **));
gsl_rng * r1 = gsl_rng_alloc (gsl_rng_rand48);
gsl_rng_set (r1, 1);
for (l1=0;l1<NbrGps;l1++){
rho[l1]=bmatrix(0, p-1, 0, H-1);
for (i=0;i<p;i++){
 unsigned int rho1[H];
gsl_ran_multinomial (r1, H,1,probInit,rho1);
// multinomial(H,1,probInit,rho1);
for (h=0;h<H;h++){
rho[l1][i][h]=rho1[h];
Rho[l1*p*H+i*H+h]=rho[l1][i][h];
}
}
}
//save3dBool("KnownRho1.txt",NbrGps,p,H, rho);
//double sig=0.2;
double mu[H][p];
for (h=0;h<H;h++){
for (i=0;i<p;i++){
mu[h][i]=gsl_ran_gaussian (r, 0.1);
// mu[h][i]=rnorm (0, 0.1);
}
}
double err=0;
double xb;
for (i=0;i<p;i++){
j1=0;
for (l1=0;l1<NbrGps;l1++){
for (j=0;j<nbrTimeMod[0]*nbrsampl;j++){
for (h=0;h<H;h++){
if (rho[l1][i][h]==1){ 
xb=0;
for (l=0;l<NbrCov;l++){
xb+=TimeCovRes[l1][j][l]*SignCoef[h][l]*beta[i];
}
y[i][j1]=mu[h][i]+xb;
err=gsl_ran_gaussian (r, sig[i]);
// err=rnorm (0, sig[i]);
}
}
y[i][j1]+=err;
yv[i*n+j1]=y[i][j1];
j1++;
}
}
}

for (l1=0;l1<NbrGps;l1++){
free(idxRes[l1]);
free_bmatrix(TimeCovRes[l1],0, nbrTimeMod[0]*nbrsampl-1, 0, NbrCov-1);
}
free(TimeCovRes);
free(idxRes);
for (h=0;h<H;h++)
free(SignCoef[h]);
free(SignCoef);
for (i=0;i<n;i++)
free(TimeCov[i]);
free(TimeCov);
for (l1=0;l1<NbrGps;l1++)
free_bmatrix(rho[l1],0, p-1, 0, H-1);//free(nbrprot);
gsl_rng_free (r1);free(rho);
gsl_rng_free (r);
free_dmatrix(y,0, p-1,0,n-1);
PutRNGstate();
}
