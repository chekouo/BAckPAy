#include <stddef.h>
#include <gsl/gsl_rng.h>
double nrs_a_inf(double a,const gsl_rng * r);
double ers_a_inf(double a,const gsl_rng * r);
double r_lefttruncnorm(double a, double mean, double sd,const gsl_rng * r) ;
double r_righttruncnorm(double b, double mean, double sd,const gsl_rng * r);
void multinomial (size_t K,unsigned int  N, const double p[], unsigned int n[]); 
void sortd(int n,double *x,int *idx);
double auc(int n, double * esti,_Bool class[n]);
void Center(int nR,int nC,double ** x);
void sort(int n,int *x,int *idx);
void readIntVector(char *filename, int nRows,  int * data );
void uniqvalues(int * v,int *n);
void Permutations(int ** Permut,int*v, int n,int p );
void NormalizeRow(int nR,int nC,double ** x);
double * Rowmean(int nR,int nC,double ** x);
double * Rowvar(int nR,int nC,double ** x);
void nrerror(char error_text[]);
void readBoolArray(char *filename, int nRows, int nCols, _Bool ** data );
void readBoolVector(char *filename, int nRows,  _Bool * data );
void readDoubleArray(char *filename, int nRows, int nCols, double ** data );
double min(int n, double * x);
double max(int n, double * x);
_Bool *bvector(int nl, int nh);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch);
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch);
double Truncate(double mu,double sd, double lower, const gsl_rng * r);
 void save3d(char *filename,int n,int p,int q, double *** data);
void save2d(char *filename,int n,int p,double ** data);
void save1ic(char *filename,int n,int * data);
void save1d(char *filename,int n,double * data);
void save3dBool(char *filename,int n,int p,int q, _Bool *** data);
void ConvertToBinaryMat(int * Time,int * uniqmodTime,int n,int nbrmodal,_Bool **TimeCov);
void ConvertToBinaryMatP(int * Time,int * uniqmodTime,int n,int nbrmodal,_Bool **TimeCov);
