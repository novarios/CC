#ifndef CCFUNCTIONS_H
#define CCFUNCTIONS_H

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
#include <numeric>
#include <omp.h>
#include <complex>

//LAPACK functions
extern "C" void dgemm_(char* ta,char* tb,int* m,int* n,int* k,double* al,double* a,int* la,double* b,int* lb,double* be,double* c,int* lc);
extern "C" void dgetrf_(int* M,int* N,double* A,int* lda,int* ipiv,int* info);
extern "C" void dgetri_(int* N,double* A,int* lda,int* ipiv,double* work,int* lwork,int* info);
extern "C" void dgeev_(char* jobvl,char* jobvr,int* N,double* A,int* lda,double* wr,double* wi,double* vl,int* ldvl,double* vr,int* ldvr,double* work,int* lwork,int* info);
#define RM_dgemm(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_B, transf_A, n, m, k, alpha, B, n, A, k, beta, C, n)
#define RM_dgemm2(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_A, transf_B, m, n, k, alpha, A, k, B, n, beta, C, n)


#define min(a,b) (a <= b ? a : b)
#define max(a,b) (a >= b ? a : b)


const std::string PATH = "inputs/";

struct Input_Parameters;
struct Model_Space;
struct Model_Space_J;
struct Channels;
struct CCD;
struct CC_Matrix_Elements;
struct CC_Eff;
struct V_Conv;

int Index(const std::vector<std::vector<int> > &vec1, const std::vector<std::vector<int> > &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index2(const std::vector<int> &vec1, const std::vector<std::vector<int> > &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
Channels HO_Setup_Channels(const Model_Space &Space);
Channels CART_Setup_Channels(const Model_Space &Space);
Input_Parameters Get_Input_Parameters(std::string &infile);
Model_Space Build_Model_Space(Input_Parameters &Parameters);
Model_Space_J Build_Model_Space_J1(Input_Parameters &Parameters);
Model_Space Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space_J &Space_J);
Model_Space CART_Build_Model_Space(Input_Parameters &Parameters);
Model_Space CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &tx, const double &ty, const double &tz);
CC_Matrix_Elements Read_Matrix_Elements(const std::string &MEfile, const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
CC_Matrix_Elements Read_Matrix_Elements_J(const std::string &MEfile, const Input_Parameters &Parameters, const Model_Space &Space, const Model_Space_J &Space_J, const Channels &Chan);
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
int HO_tbInd1(const Model_Space &Space, const double &P, const double &M, const double &T);
int HO_tbInd2(const Model_Space &Space, const double &P2, const double &M2, const double &T2);
int CART_tbInd2_rel(const Model_Space &Space, const int &Nx2, const int &Ny2, const int &Nz2, const double &m1, const double &m2);
void Print_Parameters(const Input_Parameters &Parameters);
CC_Eff Build_CC_Eff(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, CCD &CC, const Channels &Chan);
void EE_EOM(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan);
void EE_EOM_HO(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan);
std::vector<unsigned long long> bitsetup(const int &i, const int &j, const int &a, const int &b, const int &Nbit);
double matrixe1(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const double &ME, const int &test);
double matrixe2(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &r, const int &s, const double &ME);
double matrixe3(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &r, const int &t, const int &u, const int &v, const double &ME);
double vint_Conversion(const Model_Space &Space, const std::vector<std::vector<double> > &ME, const int &qi, const int &qj, const int &qk, const int &ql);

//Structure for holding Input parameters
struct Input_Parameters{
  std::string basis; //HO or CART
  double obstrength; //one-body multiplication strength
  double tbstrength; //two-body multiplication strength
  int Nshells; //number of neutrons shells
  int Pshells; //number of protons shells
  int N; //number of neutrons
  int P; //number of protons
  double density;
  int Nmax;
  std::string LevelScheme; //level scheme path
  std::string MatrixElements; //matrix elements path
  //For Excited States
  int Nx, Ny, Nz;
  double M, T, Par;
};

//Structure for holding all model space info
struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits
  std::vector<int> levelsind; //list of single particle state indicies (1,2...)
  std::vector<int> levelsm; //list of single particle state total angular momentum projection
  std::vector<int> levelst; //list of single particle state isospins
  std::vector<std::string> levelstype; //list of single particle state types
  std::vector<double> levelsen; //list of single particle energies
  //for HO
  std::vector<int> levelsn; //list of single particle state principal quantum numbers
  std::vector<int> levelsl; //list of single particle state orbital angular momentum
  std::vector<double> levelsj; //list of single particle state orbital angular momentum
  int Pmin;
  int HO_tb1Indsize;
  std::vector<int> HO_tb1Indvec;
  int P2min;
  int HO_tb2Indsize;
  std::vector<int> HO_tb2Indvec;
  
  //for CART
  std::vector<int> levelsnx; //list of single particle state x-momeuntum quantum number
  std::vector<int> levelsny; //list of single particle state y-momeuntum quantum number
  std::vector<int> levelsnz; //list of single particle state z-momeuntum quantum number
  int nmax;
  int Nxmin;
  int CART_tb1Indsize;
  std::vector<int> Nymin;
  std::vector<std::vector<int> > Nzmin;
  std::vector<std::vector<std::vector<int> > > CART_tb1Indvec;
  int Nx2min;
  int CART_tb2Indsize;
  std::vector<int> Ny2min;
  std::vector<std::vector<int> > Nz2min;
  std::vector<std::vector<std::vector<int> > > CART_tb2Indvec;

  int Mmin, Msize, M2min, M2size;
  int Tmin, Tsize, T2min, T2size;
};

//Structure for holding all model space info
struct Model_Space_J{
  std::vector<int> indvec;
  std::vector<int> nvec;
  std::vector<int> lvec;
  std::vector<double> jvec;
  std::vector<double> tzvec;
  std::vector<double> envec;
  std::vector<std::vector<int> > shellsm;
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

struct CC_Eff{
  std::vector<std::vector<double> > V1;
  std::vector<std::vector<double> > V2;
  std::vector<std::vector<double> > V3;
  std::vector<std::vector<double> > V4;
  std::vector<std::vector<double> > V5;
  std::vector<std::vector<double> > V6;
  std::vector<std::vector<double> > V7;
  CC_Eff(Channels);
};

struct V_Conv{
  /*std::vector<double> JME;
  std::vector<std::vector<int> > LmlS;
  std::vector<std::vector<int> > JJz;
  std::vector<std::vector<double> > sigmak;
  std::vector<std::vector<std::complex<double> > > C1;
  std::vector<std::vector<std::complex<double> > > C2;
  std::vector<std::vector<std::complex<double> > > Y1;
  std::vector<std::vector<std::complex<double> > > Y2;*/
  
  std::vector<double> JME;
  std::vector<std::vector<int> > JL;
  std::vector<std::vector<int> > SzTz;
  std::vector<std::vector<int> > Skhat0;
  std::vector<std::vector<int> > Skhat1;

  std::vector<std::vector<double> > Y1_JLSk;
  std::vector<std::vector<double> > Y2_JLSk;
  std::vector<std::vector<double> > V_JL;

  V_Conv(Channels, Input_Parameters, Model_Space, int);
};

#endif
