#include <stddef.h>
double Lower1(double  lower, double beta,double  sigma2, double c_beta,double a,double b,double * accept,const gsl_rng * r);
int Rho(int h_old,double * data, _Bool ** Time, double **beta,double *mu, double *sigma2,double *alpha,double sumalpha,double *C,_Bool* rho, int *nbrprot,int p,int n,int nbrcov,int H, int **signCoef,const gsl_rng * r);
double Mu(double *** data,_Bool *** Time,int h, double *beta,double sigma2,double C,double c_mu,_Bool*** rho, int *nbrprot,int *P,int* N,int nbrcov, int nbrgrp,int *signCoef,const gsl_rng * r);

double Sigma2(double *** data,_Bool *** Time, int h, double * beta,double mu,double C,double a,double b,_Bool*** rho, int* nbrprot,int *P,int* N,int nbrcov, int nbrgrp,int* signCoef,const gsl_rng * r);


double Beta(double *** data,_Bool *** Time,int h,int l, double *beta,double sigma2,double mu,double C,double c_beta,_Bool*** rho, int *nbrprot,int *P,int* N,int nbrcov,int nbrgrp, double lower,int *signCoef,const gsl_rng * r);

double logpost(double *** data, _Bool *** Time, double **beta,double *mu, double *sigma2,double *alpha,double sumalpha,double *C,_Bool*** rho, int *nbrprot,int* P,int *N,int nbrcov,int nbrgrp,int H, int ** signCoef, double a, double b, double al, double bl,double c_mu, double c_beta, double** truncbeta);

double Lower(double  lower, double ** beta,double * sigma2, double c_beta,int nbrcov,int H, int ** signCoef,double a,double b,double * accept,const gsl_rng * r);


