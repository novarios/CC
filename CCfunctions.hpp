#ifndef COUPLEDCLUSTER_H
#define COUPLEDCLUSTER_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <bitset>
#include <iomanip>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <time.h>
#include "omp.h"

//LAPACK functions
extern "C" void dgemm_(char* ta,char* tb,int* m,int* n,int* k,double* al,double* a,int* la,double* b,int* lb,double* be,double* c,int* lc);
extern "C" void dgetri_(int *n, double*a, int *lda, int *ipiv, double *work, int *lwork, int *info);
extern "C" void dgetrf_(int *m, int *n, double*a, int *lda, int *ipiv, int *info);
#define RM_dgemm(a, b, c, d, e, f, g, h, i, j, k, l, m) dgemm_(b, a, d, c, e, f, i, d, g, j, k, l, d)

const std::string PATH = "inputs/";

struct Input_Parameters;
struct Model_Space;
struct Channels;
struct CCD;
struct CC_Matrix_Elements;

int Index(const std::vector<std::vector<int> > &vec1, const std::vector<std::vector<int> > &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index2(const std::vector<int> &vec1, const std::vector<std::vector<int> > &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
//Channels HO_Setup_Channels(const Model_Space &Space);
Channels CART_Setup_Channels(const Model_Space &Space);
Input_Parameters Get_Input_Parameters(std::string &infile);
//Model_Space Build_Model_Space(const Input_Parameters &Parameters);
Model_Space CART_Build_Model_Space(const Input_Parameters &Parameters);
//CC_Matrix_Elements Read_Matrix_Elements(const std::string &MEfile, const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
CCD Perform_CCD(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, const Channels &Chan);
void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME);
double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME);
void Doubles_Step(const Model_Space &Space, const Channels &Chan, CC_Matrix_Elements &ME, CCD &CC, CCD &CC2);
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
CC_Matrix_Elements Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
int CART_tbInd1(const Model_Space &Space, const int &Nx, const int &Ny, const int &Nz, const double &M, const double &T);
int CART_tbInd2(const Model_Space &Space, const int &Nx2, const int &Ny2, const int &Nz2, const double &M2, const double &T2);
//int HO_tbInd1(const Model_Space &Space, const double &P, const double &M, const double &T);
//int HO_tbInd2(const Model_Space &Space, const double &P2, const double &M2, const double &T2);
void Print_Parameters(const Input_Parameters &Parameters);

//Structure for holding Input parameters
struct Input_Parameters{
  std::string basis; //HO or CART
  double obstrength; //one-body multiplication strength
  double tbstrength; //two-body multiplication strength
  int N; //number of valence neutrons
  int P; //number of valence protons
  double density;
  int Nmax;
  std::string LevelScheme; //level scheme path
  std::string MatrixElements; //matrix elements path
};

//Structure for holding all model space info
struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits
  std::vector<int> levelsind; //list of single particle state indicies (1,2...)
  std::vector<double> levelsm; //list of single particle state total angular momentum projection
  std::vector<double> levelst; //list of single particle state isospins
  std::vector<std::string> levelstype; //list of single particle state types
  std::vector<double> levelsen; //list of single particle energies
  //for HO
  std::vector<int> levelsn; //list of single particle state principal quantum numbers
  std::vector<int> levelsl; //list of single particle state orbital angular momentum
  //for CART
  std::vector<int> levelsnx; //list of single particle state x-momeuntum quantum number
  std::vector<int> levelsny; //list of single particle state x-momeuntum quantum number
  std::vector<int> levelsnz; //list of single particle state x-momeuntum quantum number

  int nmax;
  int Nxmin;
  int tb1Indsize;
  std::vector<int> Nymin;
  std::vector<std::vector<int> > Nzmin;
  std::vector<std::vector<std::vector<int> > > tb1Indvec;
  int Nx2min;
  int tb2Indsize;
  std::vector<int> Ny2min;
  std::vector<std::vector<int> > Nz2min;
  std::vector<std::vector<std::vector<int> > > tb2Indvec;

  //int Nxmin, Nxmax, Nxsize;
  //int Nymin, Nymax, Nysize;
  //int Nzmin, Nzmax, Nzsize;
  int Mmin, Mmax, Msize;
  int Tmin, Tmax, Tsize;
  int Pmin, Pmax, Psize;
  //int Nx2min, Nx2max, Nx2size;
  //int Ny2min, Ny2max, Ny2size;
  //int Nz2min, Nz2max, Nz2size;
  int M2min, M2max, M2size;
  int T2min, T2max, T2size;
  int P2min, P2max, P2size;
};

//Structure for holding channel information
struct Channels{
  int size1;
  int size2;
  int size3;

  std::vector<int> indvec;

  std::vector<int> hh;
  std::vector<int> pp;
  std::vector<int> hp1;
  std::vector<int> hp2;
  std::vector<int> ph1;
  std::vector<int> ph2;
  std::vector<int> h;
  std::vector<int> hpp;
  std::vector<int> p;
  std::vector<int> hhp;

  std::vector<std::vector<int> > hhvec1;
  std::vector<std::vector<int> > ppvec1;
  std::vector<std::vector<int> > hp1vec1;
  std::vector<std::vector<int> > hp2vec1;
  std::vector<std::vector<int> > pvec1;
  std::vector<std::vector<int> > hhpvec1;
  std::vector<std::vector<int> > hvec1;
  std::vector<std::vector<int> > hppvec1;
};

struct CCD{
  std::vector<std::vector<int> > Tmap;
  std::vector<std::vector<double> > Evec;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > T2T;
  std::vector<std::vector<double> > T3;
  std::vector<std::vector<double> > T4;
  std::vector<std::vector<double> > T5;
  std::vector<std::vector<double> > T6;
  std::vector<std::vector<double> > T7;
  std::vector<std::vector<double> > S1;
  std::vector<std::vector<double> > S2;
  std::vector<std::vector<double> > S3;
  std::vector<std::vector<double> > S4;
  std::vector<std::vector<double> > S4T;
  double CCDE;
  CCD(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  double get_T(int, int) const;
};

struct CC_Matrix_Elements{
  std::vector<std::vector<double> > PPPP;
  std::vector<std::vector<double> > HHHH;
  std::vector<std::vector<double> > HPHP1;
  std::vector<std::vector<double> > HPHP2;
  std::vector<std::vector<double> > HHPP1;
  std::vector<std::vector<double> > HHPP2;
  std::vector<std::vector<double> > HHPP3;
  std::vector<std::vector<double> > HHPP4;
  std::vector<std::vector<double> > HHPP4T;
  CC_Matrix_Elements(Channels);
};

#endif
