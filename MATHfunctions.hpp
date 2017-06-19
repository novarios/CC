#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <complex>
#include <vector>

#define PI 3.141592653589793
#define e 2.718281828459045

#define hbarc_MeVfm 197.3269788 // MeV fm
#define hbarc_eVum 0.1973269788 // eV um
#define hbarc_HartA 72.5163303659 // Hart A
#define m_neutronc2 939.5654133 // MeV
#define m_protonc2 938.2720813 // MeV
//#define m_protonc2 939.565378 // MeV
#define m_electronc2 0.5109989461 // MeV
#define m_electronc2_Hart 18778.865727 // Hart
#define eVs_in_Hartree 27.21138505 // eV
#define fine_struct 0.007297352566355

//LAPACK functions
extern "C" void dgemm_(char* ta,char* tb,int* m,int* n,int* k,double* al,double* a,int* la,double* b,int* lb,double* be,double* c,int* lc);
extern "C" void dgetrf_(int* M,int* N,double* A,int* lda,int* ipiv,int* info);
extern "C" void dgetri_(int* N,double* A,int* lda,int* ipiv,double* work,int* lwork,int* info);
extern "C" void dgeev_(char* jobvl,char* jobvr,int* N,double* A,int* lda,double* wr,double* wi,double* vl,int* ldvl,double* vr,int* ldvr,double* work,int* lwork,int* info);

extern "C" void dnaupd_(int* ido,char* bmat,int* N,char* which,int* nev,double* tol,double* resid,int* ncv,double* v,int* ldv,int* iparam,int* ipntr,double* workd,double* workl,int* lworkl,int* info);
extern "C" void dneupd_(bool* rvec,char* howmny,int* select,double* dr,double* di,double* z,int* ldz,double* sigmar,double* sigmai,double* workev,char* bmat,int* N,char* which,int* nev,double* tol,double* resid,int* ncv,double* v,int* ldv,int* iparam,int* ipntr,double* workd,double* workl,int* lworkl,int* info);

#define dgemm_NN(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, n, A, k, beta, C, n)
#define dgemm_NT(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, k, A, k, beta, C, n)
#define dgemm_TN(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, n, A, m, beta, C, n)
#define dgemm_TT(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, k, A, m, beta, C, n)

long long factorial(const int &n);
long long factorial(const double &n);
double logfac(const int &n);
long long factorial2(const int &n);
long long factorial2(const double &n);
double CGC(double j1, double m1, double j2, double m2, double jtot, double mtot);
double CGC3(const double &j1, const double &m1, const double &j2, const double &m2, const double &jtot, const double &mtot);
double CGC6(const double &j1, const double &j2, const double &j3, const double &j4, const double &j5, const double &j6);
double Chi_J(const int &a, const int &b, const int &c, const int &d, const int &J1, const int &J2);
double Legendre(const double &x, const int &l, const int &m);
std::complex<double> SphericalY_C(const double &theta, const double &phi, const int &l, const int &m);
double SphericalY(const double &theta, const double &phi, const int &l, const int &m);
double SphericalYTens(const double &theta, const double &phi, const double &j, const int &l, const int &s, const int &ms);
double Erf(const double &z);
double Laguerre(const int &k, const double &alpha, const double &x);
double HOfunction(const double &hw, const int &k, const int &l, const int &m, const double &r, const double &theta, const double &phi);
int choose(const int &int1, const int &int2);
double logratio1(const int &int1, const int &int2, const int &int3, const int &int4);
double logratio2(const int &G);
double product1(const int &n1, const int &m1, const int &n2, const int &m2, const int &n3, const int &m3, const int &n4, const int &m4);
double logproduct2(const int &n1, const int &m1, const int &n2, const int &m2, const int &n3, const int &m3, const int &n4, const int &m4, const int &j1, const int &j2, const int &j3, const int &j4);
double logproduct3(const int &l1, const int &l2, const int &l3, const int &l4, const int &g1, const int &g2, const int &g3, const int &g4);
double loggamma(const double &x);
void projection(double *u, double *v, double *proj, const int &size);
void GramSchmidt(double *Vectors, const int &size);
void GramSchmidt(double **Vectors, const int &size);
double rand_normal(double mean, double stddev);
void Asym_Diagonalize1(double *Ham, int &N, double &eigenvalue, double &norm1p, int &np0);
void Asym_Diagonalize2(double *Ham, int &N, double &eigenvalue, double &norm1p, int &np0);

#endif
