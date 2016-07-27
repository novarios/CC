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
struct State;
struct Model_Space;
struct Channels;
struct Amplitudes;
struct Interactions;

struct Model_Space2;
struct Channels2;

struct Doubles_1;
struct Singles_1;
struct Doubles_ME1;
struct Singles_ME1;
struct CC_Eff;

//struct V_Conv;

int Index1(const std::vector<int> &vec1, const int &p);
int Index11(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q);
int Index2(const std::vector<int> &vec1, const int &p, const int &q);
int Index22(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index13(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index31(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);

int ChanInd_2b_dir2(const Model_Space2 &Space, const State &State);
int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State);
int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State);
void plus(State &S, const State &S1, const State &S2);
void minus(State &S, const State &S1, const State &S2);

void Setup_Channels(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan);

void Get_Input_Parameters(std::string &infile, Input_Parameters &Parameters);
void Print_Parameters(const Input_Parameters &Parameters);

void Setup_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps);
void Setup_Ints(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints);

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space);
void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void CART_Build_Model_Space2(Input_Parameters &Parameters, Model_Space2 &Space);
void EG_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
//void CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &tx, const double &ty, const double &tz);

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
//void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const Model_Space &Space_J, const Channels &Chan, Interactions &Ints);

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps);
void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);

void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, Interactions &Int);
double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Int);

double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
double vint_Minnesota_Momentum2(const Model_Space2 &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql);

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff);
void EE_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff);
void bitsetup(const std::vector<int> &vec, std::vector<unsigned long long> &state);
double matrixe(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const std::vector<int> &p, const std::vector<int> &q, const double &ME);

void Setup_Channels_MBPT(const Input_Parameters &Parameters, const Model_Space2 &Space, Channels2 &Chan);
void Perform_MBPT_0(const Input_Parameters &Parameters, const Model_Space2 &Space);
void Perform_MBPT_1(const Input_Parameters &Parameters, const Model_Space2 &Space);
void Perform_MBPT_2(const Input_Parameters &Parameters, const Model_Space2 &Space, const Channels2 &Chan);
void Perform_MBPT_3(const Input_Parameters &Parameters, const Model_Space2 &Space, const Channels2 &Chan);
void Perform_MBPT_4(const Input_Parameters &Parameters, const Model_Space2 &Space, const Channels2 &Chan);
void Perform_MBPT_5(const Input_Parameters &Parameters, const Model_Space2 &Space, const Channels2 &Chan);

//Structure for holding Input parameters
struct Input_Parameters{
  std::string calc_case; //nuclear or electronic
  std::string basis; //infinite, finite (HO), or finite_J (HO_J)
  std::string approx; //doubles, singles, or triples
  int HF; //1 to perform HF, anything else for bypass
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
  int extra; // -1 for pr, 0 for es, 1 for pa
  int Nx, Ny, Nz;
  double M, T, Par;
  int MBPT; // 0 for serial, 1 for parallel, 2 for block serial, 3 for block parallel, 4 for block M-M
};

struct State{
  int t; // x2
  int m; // x2
  int nx;
  int ny;
  int nz;
  int ml;
  int n;
  int j; // x2
  int par; // -1,+1
  double energy;
  std::string type;
};

//Structure for holding all model space info
struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits

  std::vector<struct State> qnums;
  struct State qmins;
  struct State qmaxs;
  struct State qsizes;
  std::vector<std::vector<int> > shellsm; // for j

  int nmax;
  std::vector<int> map_2b;
  int Chansize_2b;
};

//Structure for holding all model space info
struct Model_Space2{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits

  State *qnums;
  struct State qmins;
  struct State qmaxs;
  struct State qsizes;

  int nmax;
  int *map_2b;
  int Chansize_2b;
};

//Structure for holding channel information
struct Channels{
  int size1;
  int size2;
  int size3;

  std::vector<struct State> qnums1;
  std::vector<struct State> qnums2;
  std::vector<struct State> qnums3;
  
  std::vector<int> indvec;

  std::vector<int> hh;
  std::vector<int> pp;
  std::vector<int> hp;
  std::vector<int> hp1;
  std::vector<int> hp2;
  std::vector<int> ph1;
  std::vector<int> ph2;
  std::vector<int> h;
  std::vector<int> p;
  std::vector<int> hpp;
  std::vector<int> hhp;
  std::vector<int> hpp2;
  std::vector<int> hhp2;
  std::vector<int> hh1;
  std::vector<int> pp1;
  std::vector<int> hhh;
  std::vector<int> ppp;

  std::vector<std::vector<int> > hhvec1;
  std::vector<std::vector<int> > ppvec1;
  std::vector<std::vector<int> > hpvec1;
  std::vector<std::vector<int> > hp1vec1;
  std::vector<std::vector<int> > hp2vec1;
  std::vector<std::vector<int> > pvec1;
  std::vector<std::vector<int> > hvec1;
  std::vector<std::vector<int> > hhpvec1;
  std::vector<std::vector<int> > hppvec1;
  std::vector<std::vector<int> > hhp2vec1;
  std::vector<std::vector<int> > hpp2vec1;
  std::vector<std::vector<int> > hh1vec1;
  std::vector<std::vector<int> > pp1vec1;
  std::vector<std::vector<int> > hhhvec1;
  std::vector<std::vector<int> > pppvec1;

  int ind0; // index of i-i cross channel for singles
  std::vector<int> hhppsize;
};

//Structure for holding channel information
struct Channels2{
  int size1;
  int size3;

  int* hh;
  int* pp;
  int** hhvec1;
  int** ppvec1;
};

//Structure for holding channel information
/*struct HF_Channels_M{
  int size1;
  int size3;
  std::vector<int> indvec;
  std::vector<int> tb;
  std::vector<std::vector<int> > tbvec1;
  std::vector<int> ob;
  std::vector<std::vector<int> > obvec1;
  };*/

struct Doubles_1{
  std::vector<std::vector<int> > Tmap;
  std::vector<std::vector<double> > Evec;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > T3;
  std::vector<std::vector<double> > T4;
  std::vector<std::vector<double> > T5;
  std::vector<std::vector<double> > T6;
  std::vector<std::vector<double> > T7;
  std::vector<std::vector<double> > T8;
  std::vector<std::vector<double> > T9;
  std::vector<std::vector<double> > S1;
  std::vector<std::vector<double> > S2;
  std::vector<std::vector<double> > S3;
  std::vector<std::vector<double> > S4;
  std::vector<std::vector<double> > S5;
  std::vector<std::vector<double> > S6;
  std::vector<std::vector<double> > S7;
  std::vector<std::vector<double> > Q11;
  std::vector<std::vector<double> > Q21;
  std::vector<std::vector<double> > Q12;
  std::vector<std::vector<double> > Q22;
  std::vector<std::vector<int> > Qmap1;
  std::vector<std::vector<int> > Qmap2;

  Doubles_1(); //default constructor
  Doubles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int, int) const;
};

struct Singles_1{
  std::vector<int> Tmap;
  std::vector<int> Tmap2;
  std::vector<double> Evec;
  std::vector<double> T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > T3;
  std::vector<std::vector<double> > S1;
  std::vector<std::vector<double> > S2;
  std::vector<double> S3;
  std::vector<double> S4;
  std::vector<std::vector<double> > E1;
  std::vector<std::vector<double> > E2;
  std::vector<std::vector<double> > E3;
  std::vector<std::vector<double> > E4;
  std::vector<std::vector<double> > E5;
  std::vector<std::vector<double> > E6;
  std::vector<std::vector<double> > E7;
  std::vector<std::vector<double> > E8;
  std::vector<std::vector<double> > E9;
  std::vector<std::vector<double> > Q11;
  std::vector<std::vector<double> > Q12;
  std::vector<std::vector<double> > Q21;
  std::vector<std::vector<double> > Q22;
  std::vector<double> Q31;
  std::vector<std::vector<double> > Q32;
  std::vector<double> Q41;
  std::vector<std::vector<double> > Q42;
  std::vector<std::vector<double> > Q51;
  std::vector<std::vector<double> > Q52;
  std::vector<std::vector<double> > Q61;
  std::vector<std::vector<double> > Q62;
  std::vector<std::vector<int> > Qmap1;
  std::vector<std::vector<int> > Qmap2;
  std::vector<int> Qmap3;
  std::vector<int> Qmap4;
  std::vector<std::vector<int> > Qmap5;
  std::vector<std::vector<int> > Qmap6;
  Singles_1(); //default constructor
  Singles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int) const;
};

struct Amplitudes{
  Doubles_1 D1; // for doubles only
  Singles_1 S1; // for singles part of singles
  double get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints);
  void zero(const Channels &Chan, const Input_Parameters &Parameters);
};

struct Doubles_ME1{
  std::vector<std::vector<double> > V1;
  std::vector<std::vector<double> > V2;
  std::vector<std::vector<double> > V3;
  std::vector<std::vector<double> > V4;
  std::vector<std::vector<double> > V5;
  std::vector<std::vector<double> > V6;
  std::vector<std::vector<double> > V7;
  std::vector<std::vector<double> > V8;
  std::vector<std::vector<double> > V9;
  std::vector<std::vector<double> > V10;
  Doubles_ME1();
  Doubles_ME1(const Channels &Chan);
};

struct Singles_ME1{
  std::vector<std::vector<double> > V11;
  std::vector<std::vector<double> > V12;
  std::vector<std::vector<double> > V15;
  std::vector<std::vector<double> > V16;
  std::vector<std::vector<double> > V13;
  std::vector<std::vector<double> > V14;
  std::vector<std::vector<double> > V19;
  std::vector<std::vector<double> > V20;
  std::vector<std::vector<double> > V17;
  std::vector<std::vector<double> > V18;
  Singles_ME1();
  Singles_ME1(const Channels &Chan);
};

struct Interactions{
  Doubles_ME1 D_ME1; // for doubles only
  Singles_ME1 S_ME1; // for singles part of singles
};

struct CC_Eff{
  std::vector<double> X_ia1;
  std::vector<std::vector<double> > X_ia2;
  std::vector<std::vector<double> > X_ia3;
  std::vector<int> Map_ia;

  std::vector<double> X_ab1;
  std::vector<std::vector<double> > X_ab2;
  std::vector<std::vector<double> > X_ab3;
  std::vector<int> Map_ab;

  std::vector<double> X_ij1;
  std::vector<std::vector<double> > X_ij2;
  std::vector<std::vector<double> > X_ij3;
  std::vector<double> X1_ij1;
  std::vector<std::vector<double> > X1_ij2;
  std::vector<std::vector<double> > X1_ij3;
  std::vector<int> Map_ij;

  std::vector<double> X_ai1;
  std::vector<std::vector<double> > X_ai2;
  std::vector<std::vector<double> > X_ai3;
  std::vector<int> Map_ai;

  std::vector<std::vector<double> > X_ijab1;

  std::vector<std::vector<double> > X1_iabc1;
  std::vector<std::vector<double> > X1_iabc2;
  std::vector<std::vector<double> > X1_iabc3;
  std::vector<std::vector<double> > X_iabc1;
  std::vector<std::vector<double> > X_iabc3;
  std::vector<std::vector<double> > X_iabc4;
  std::vector<std::vector<double> > X_iabc5;
  std::vector<std::vector<int> > Map_iabc;

  std::vector<std::vector<double> > X1_ijka1;
  std::vector<std::vector<double> > X1_ijka2;
  std::vector<std::vector<double> > X_ijka1;
  std::vector<std::vector<double> > X_ijka4;
  std::vector<std::vector<double> > X_ijka5;
  std::vector<std::vector<int> > Map_ijka;

  std::vector<std::vector<double> > X1_abcd1;
  std::vector<std::vector<double> > X1_abcd2;
  std::vector<std::vector<double> > X1_abcd3;
  std::vector<std::vector<double> > X_abcd1;
  std::vector<std::vector<double> > V_abcd;
  std::vector<std::vector<int> > Map_abcd;

  std::vector<std::vector<double> > X_ijkl1;
  std::vector<std::vector<double> > X_ijkl2;
  std::vector<std::vector<double> > X_ijkl3;
  std::vector<std::vector<double> > X_ijkl4;
  std::vector<std::vector<double> > V_ijkl;
  std::vector<std::vector<int> > Map_ijkl;

  std::vector<std::vector<double> > X1_iajb1;
  std::vector<std::vector<double> > X1_iajb2;
  std::vector<std::vector<double> > X1_iajb3;
  std::vector<std::vector<double> > X1_iajb4;
  std::vector<std::vector<double> > X3_iajb1;
  std::vector<std::vector<double> > X3_iajb2;
  std::vector<std::vector<double> > X3_iajb3;
  std::vector<std::vector<double> > X3_iajb5;
  std::vector<std::vector<double> > X_iajb1;
  std::vector<std::vector<double> > X_iajb3;
  std::vector<std::vector<int> > Map_iajb;

  std::vector<std::vector<double> > X_abic1;
  std::vector<std::vector<double> > X_abic2;
  std::vector<std::vector<double> > X_abic3;
  std::vector<std::vector<double> > X_abic4;
  std::vector<std::vector<double> > X_abic5;
  std::vector<std::vector<double> > X_abic6;
  std::vector<std::vector<double> > X_abic7;
  std::vector<std::vector<int> > Map_abic;

  std::vector<std::vector<double> > X2_iajk1;
  std::vector<std::vector<double> > X2_iajk2;
  std::vector<std::vector<double> > X2_iajk3;
  std::vector<std::vector<double> > X2_iajk4;
  std::vector<std::vector<double> > X2_iajk5;
  std::vector<std::vector<double> > X2_iajk6;
  std::vector<std::vector<double> > X2_iajk7;
  std::vector<std::vector<double> > X_iajk1;
  std::vector<std::vector<int> > Map_iajk;

  CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
  void set_X_ia(const Channels &Chan);
  void set_X_ab(const Channels &Chan);
  void set_X_ij(const Channels &Chan);
  void set_X1_ij(const Channels &Chan);
  void set_X_ai(const Channels &Chan);
  void set_X1_iabc(const Channels &Chan);
  void set_X_iabc(const Channels &Chan);
  void set_X1_ijka(const Channels &Chan);
  void set_X_ijka(const Channels &Chan);
  void set_X1_abcd(const Channels &Chan);
  void set_X_ijkl(const Channels &Chan);
  void set_X1_iajb(const Channels &Chan);
  void set_X3_iajb(const Channels &Chan);
  void set_X_iajb(const Channels &Chan);
  void set_X_abic(const Channels &Chan);
  void set_X2_iajk(const Channels &Chan);
};

/*struct V_Conv{  
  std::vector<double> JME;
  std::vector<std::vector<int> > JL;
  std::vector<std::vector<int> > SzTz;
  std::vector<std::vector<int> > Skhat0;
  std::vector<std::vector<int> > Skhat1;

  std::vector<std::vector<double> > Y1_JLSk;
  std::vector<std::vector<double> > Y2_JLSk;
  std::vector<std::vector<double> > V_JL;

  V_Conv(Channels, Input_Parameters, Model_Space, int);
};*/

#endif
