#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include <R.h>
#include <Rmath.h>
static const double t4 = 0.45;

void multinomial (size_t K,unsigned int  N, const double p[], unsigned int n[])
{
   GetRNGstate();
  size_t k;
  double norm = 0.0;
  double sum_p = 0.0;
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */

  for (k = 0; k < K; k++)
    {
      norm += p[k];
    }

  for (k = 0; k < K; k++)
    {
      if (p[k] > 0.0)
        {
          n[k] = rbinom (N - sum_n,p[k] / (norm - sum_p));
        }
      else
        {
          n[k] = 0;
        }

      sum_p += p[k];
      sum_n += n[k];
    }
 PutRNGstate();
}




void sortd(int n,double *x,int *idx)
{
int i,j;
double a;
int id;
for (i = 0; i < n; i++)
idx[i]=i;
for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (x[i] <= x[j])
            {
                a =  x[i];
             id=idx[i];
                idx[i]=idx[j];
                x[i] = x[j];
                idx[j]=id;
                x[j] = a;
            }
        }
    }

}

double auc(int n, double * esti,_Bool class[n]){
double fpr[n+2],tpr[n+2];
double auc1=0;
int P=0;//P=positive instances
int i,j;
double esti1[n];
for (i=0;i<n;i++){
esti1[i]=esti[i];
if (class[i]==1)
P+=1;
}
int idx[n];
sortd(n,esti1,idx);
fpr[n+1]=1;tpr[n+1]=1;
fpr[0]=0;tpr[0]=0;
for (i=n;i>=1;--i){
double af=0;double at=0;
for (j=0;j<n;j++){
if (esti[j]>esti1[i-1]){
if (class[j]==0){
af+=1;
}
else {
at+=1;
} } }
tpr[i]=at/P;
fpr[i]=af/(n-P);
auc1+=(fpr[i+1]-fpr[i])*(tpr[i+1]+tpr[i]);
}
auc1+=(fpr[1]-fpr[0])*(tpr[1]+tpr[0]);
auc1=0.5*(auc1);
return auc1;
}

void uniqvalues(int * v,int *n){
int i,j;
int v1[*n];
int countd=0;
v1[countd]=v[countd];
for (i=0;i<*n;i++){
for (j=0;j<countd;j++){
if (v[i]==v1[j])
break;
}
if (j==countd){
v1[countd]=v[i];
countd++;
}
}
for (j=0;j<countd;j++)
v[j]=v1[j];
*n=countd;
}

//_Bool **TimeCov=bmatrix(0, n-1, 0, NbrCov-1);
void ConvertToBinaryMat(int * Time,int* uniqmodTime,int n,int nbrmodal,_Bool **TimeCov){
int j,l1;
int NbrCov=nbrmodal-1;
for (j=0;j<n;j++){
for (l1=0;l1<NbrCov;l1++){
TimeCov[j][l1]=0;
if (l1==0){
if (Time[j]==uniqmodTime[l1])
TimeCov[j][l1]=1;
}
else if (l1==NbrCov-1){
if (Time[j]==uniqmodTime[l1+1])
TimeCov[j][l1]=1;
}
else if (nbrmodal>3){
if ((Time[j]==uniqmodTime[l1+1])||(Time[j]==uniqmodTime[l1+2]))
TimeCov[j][l1]=1;
}
}
}
}

void ConvertToBinaryMatP(int * Time,int* uniqmodTime,int n,int nbrmodal,_Bool **TimeCov){
int j,l1;
int NbrCov=nbrmodal-1;
for (j=0;j<n;j++){
for (l1=0;l1<NbrCov;l1++){
TimeCov[j][l1]=0;
int ss=0;
int te=0;
for (ss=l1+1;ss<nbrmodal;ss++){
if (Time[j]==uniqmodTime[ss])
te+=1;
}
if (te>0) TimeCov[j][l1]=1;
}
}
}








void nrerror(char error_text[])
{
        Rprintf("Utils run-time error...\n");
        Rprintf("%s\n",error_text);
//        Rprintf("...now exiting to system...\n");
  //      exit(1);
}
_Bool *bvector(int nl, int nh)
{
        _Bool *v;

        v=(_Bool *)malloc( (nh-nl+1)*sizeof(_Bool));
        if (!v) nrerror("allocation failure in dvector()");
        return v-nl;
}


double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));

        free((char*) (m+nrl));
}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
        int i;
        _Bool **m;

        m=(_Bool **) malloc( (nrh-nrl+1)*sizeof(_Bool*));
        if (!m) nrerror("allocation failure 1 in dmatrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++)
   {
                m[i]=(_Bool *) malloc((nch-ncl+1)*sizeof(_Bool));
                if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
                m[i] -= ncl;
        }
        return m;
}
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));

        free((char*) (m+nrl));
}


double max(int n, double * x){
double xmax=x[0];
int i;
for (i=0;i<n;i++){
if (x[i]>=xmax)
xmax=x[i];
}
return xmax;
}
double min(int n, double * x){
double xmin=x[0];
int i;
for (i=0;i<n;i++){
if (x[i]<=xmin)
xmin=x[i];
}
return xmin;
}

double Truncate(double mu,double sd, double lower, const gsl_rng * r){
double lowern=(lower-mu)/sd;
double alphaopt=(lowern+sqrt(pow(lowern,2)+4))/2;
double z=lowern+gsl_ran_exponential (r, 1/alphaopt);
//double z=lowern+rexp (1/alphaopt);
double qz=exp(-pow(z-alphaopt,2)/2);
double u=gsl_ran_flat (r, 0,1);
//double u=runif(0,1);
while (u>qz){
z=lowern+gsl_ran_exponential (r, 1/alphaopt);
// z=lowern+rexp (1/alphaopt);
qz=exp(-pow(z-alphaopt,2)/2);
u=gsl_ran_flat (r, 0,1);
//u=runif(0,1);
}
return z*sd+mu;
}

    
void NormalizeRow(int nR,int nC,double ** x){
double *Rowmea =Rowmean(nR,nC,x);
double *RowVar =Rowvar(nR,nC,x);
int i,j;
double th=pow(pow(10,-8),2);
for (i=0;i<nR;i++){
if (RowVar[i]<th){
RowVar[i]=th;
}
for (j=0;j<nC;j++){
x[i][j]=(x[i][j]-Rowmea[i])/sqrt(RowVar[i]);
}
}
free(Rowmea);free(RowVar);
}

void Center(int nR,int nC,double ** x){
double *Rowmea =Rowmean(nR,nC,x);
int i,j;
for (i=0;i<nR;i++){
for (j=0;j<nC;j++)
x[i][j]=(x[i][j]-Rowmea[i]);
}
free(Rowmea);
}


double * Rowmean(int nR,int nC,double ** x){
int i,j;
double *Mean=malloc(nR*sizeof(double));
for (i=0;i<nR;i++){
double me=0;
for (j=0;j<nC;j++)
me+=x[i][j];
Mean[i]=me/nC;
}
return Mean;
}
double * Rowvar(int nR,int nC,double ** x){
int i,j;
double *Rowmea1 =Rowmean(nR,nC,x);
double *RowVar=malloc(nR*sizeof(double));
for (i=0;i<nR;i++){
double va=0;
for (j=0;j<nC;j++)
va+=pow(x[i][j]-Rowmea1[i],2);
RowVar[i]=va/(nC-1);
}
free(Rowmea1);
return RowVar;
}

void Permutations(int ** Permut,int*v, int n,int p ){
//we count n-tuples
int k=pow(n,p);
int i,j;
for (i=0;i<k;i++){
for (j=p-1;j>=0;j--){
int bj;
bj=((int)floor(i*pow(n,-j)))%n; //gives all the binary combinations
Permut[i][j]=v[bj];
}
}
}


void sort(int n,int *x,int *idx)
{
int i,j;
double a;
int id;
for (i = 0; i < n; i++)
idx[i]=i;
for (i = 0; i < n; ++i)
    {
        for (j = i + 1; j < n; ++j)
        {
            if (x[i] >= x[j])
            {
                a =  x[i];
 	     id=idx[i];
		idx[i]=idx[j];
                x[i] = x[j];
		idx[j]=id;
                x[j] = a;
            }
        }
    }

}

/* Exponential rejection sampling (a,inf) */
double ers_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("ers_a_inf", a, R_PosInf);
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    //x = rexp(ainv) + a; /* rexp works with 1/lambda */
    x=gsl_ran_exponential (r, ainv)+a;
    rho = exp(-0.5 * pow((x - a), 2));
  } while (gsl_ran_flat (r, 0,1) > rho);
  return x;
}

double nrs_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("nrs_a_inf", a, R_PosInf);
  //double x = -DBL_AX;
  double x = gsl_ran_ugaussian(r);
  while (x < a) {
    x = gsl_ran_ugaussian(r);
  }
  return x;
}


double r_lefttruncnorm(double a, double mean, double sd,const gsl_rng * r) {
  const double alpha = (a - mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha,r);
  } else {
    return mean + sd * ers_a_inf(alpha,r);
  }
}
double r_righttruncnorm(double b, double mean, double sd,const gsl_rng * r) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0,r);
}


