#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <complex>
#include <vector>

//#define PI 3.141592653589793 // real value
#define PI 3.141592741012573 // Morten's value
#define e 2.718281828459045

//#define hbarc_MeVfm 197.3269788 // MeV fm
#define hbarc_MeVfm 197.326968 // MeV fm ( Morten's value )
#define hbarc_eVum 0.1973269788 // eV um
#define hbarc_HartA 72.5163303659 // Hart A
//#define m_neutronc2 939.5654133 // MeV
//#define m_protonc2 938.2720813 // MeV
////#define m_protonc2 939.565378 // MeV
#define m_neutronc2 938.918725 // MeV ( Morten's value )
#define m_protonc2 938.918725 // MeV ( Mortens's value )
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

struct Input_Parameters;
struct Model_Space;
struct HF_Matrix_Elements;

#define dgemm_NN(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, n, A, k, beta, C, n)
#define dgemm_NT(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, k, A, k, beta, C, n)
#define dgemm_TN(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, n, A, m, beta, C, n)
#define dgemm_TT(A, B, C, m, n, k, alpha, beta, transA, transB) dgemm_(transB, transA, n, m, k, alpha, B, k, A, m, beta, C, n)

double phase(int arg);
double phase2(int arg);
double d_ij(int i, int j);
double delta_ij(int i, int j);
long long fac(int n);
//long long factorial(double n);
double logfac(int n);
double logfac2(int n);
long long fac2(int n);
//long long factorial2(double n);
//double CGC(double j1, double m1, double j2, double m2, double jtot, double mtot);
//double CGC3(double j1, double m1, double j2, double m2, double jtot, double mtot);
//double CGC6(double j1, double j2, double j3, double j4, double j5, double j6);
double CGC(int j1, int m1, int j2, int m2, int jtot, int mtot);
double CGC3(int j1, int m1, int j2, int m2, int jtot, int mtot);
double CGC6(int j1, int j2, int j3, int j4, int j5, int j6);
double Chi_J(int a, int b, int c, int d, int J1, int J2);
double Legendre(double x, int l, int m);
std::complex<double> SphericalY_C(double theta, double phi, int l, int m);
double SphericalY(double theta, double phi, int l, int m);
double SphericalYTens(double theta, double phi, double j, int l, int s, int ms);
double Erf(double z);
double Laguerre(int k, double alpha, double x);
double HOfunction(double hw, int k, int l, int m, double r, double theta, double phi);
int choose(int int1, int int2);
double logratio1(int int1, int int2, int int3, int int4);
double logratio2(int G);
double product1(int n1, int m1, int n2, int m2, int n3, int m3, int n4, int m4);
double logproduct2(int n1, int m1, int n2, int m2, int n3, int m3, int n4, int m4, int j1, int j2, int j3, int j4);
double logproduct3(int l1, int l2, int l3, int l4, int g1, int g2, int g3, int g4);
double loggamma(double x);
void projection(double *u, double *v, double *proj, int size);
void GramSchmidt(double *Vectors, int size);
//void GramSchmidt(double **Vectors, int size);
double rand_normal(double mean, double stddev);
void Asym_Diagonalize1(double *Ham, int &N, double *eigenvalues, double *eigenvectors_L, double *eigenvectors_R, int num);
//void Asym_Diagonalize2(double *Ham, int &N, double &eigenvalue, double &norm1p, int &np0);
//void Asym_Diagonalize1(double *Ham, int &N, double *eigenvalues, double *eigenvectors_R, double *eigenvectors_L, int num);
void Asym_Diagonalize2(double *Ham, int &N, double *eigenvalues, double *eigenvectors_L, double *eigenvectors_R, int num);
void Asym_Diagonalize2_0(double *Ham, int &N, double *eigenvalues, double *eigenvectors, int num, char type);

void gauss_legendre(double x1, double x2, double *x, double *w, int n);
void laguerre_general(int n, double alpha, double x, double *cx);
void setup_ho_cutoff(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME);
void ho_wfunction(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME);
double kinetic_energy(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME, int a, int c);

double CGC(int j1, int m1, int j2, int m2, int jtot, int mtot); // for 2*j
double CGC3(int j1, int m1, int j2, int m2, int jtot, int mtot); // for 2*j
double CGC6(int j1, int j2, int j3, int j4, int j5, int j6); // for 2*j

#endif
