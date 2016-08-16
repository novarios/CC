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

struct Doubles_1;
struct Singles_1;
struct Doubles_ME1;
struct Singles_ME1;
struct CC_Eff;

//struct V_Conv;

int Index11(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q);
int Index2(const int *vec1, const int &num1, const int &p, const int &q);
int Index1(const int *vec1, const int &num1, const int &p);
int Index22(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index13(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index31(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);

/*int Index1(const std::vector<int> &vec1, const int &p);
int Index11(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q);
int Index2(const std::vector<int> &vec1, const int &p, const int &q);
int Index22(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index13(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index31(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);*/

int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State);
int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State);
void plus(State &S, const State &S1, const State &S2);
void minus(State &S, const State &S1, const State &S2);
bool equal(const State &S1, const State &S2);

//void Setup_Channels(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan);
//void Delete_Channels(Channels &Chan);

void Get_Input_Parameters(std::string &infile, Input_Parameters &Parameters);
void Print_Parameters(const Input_Parameters &Parameters);

//void Setup_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps);
//void Delete_Amps(const Input_Parameters &Parameters, const Channels &Chan, Amplitudes &Amps);
//void Setup_Ints(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints);
//void Delete_Ints(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints);

//void Delete_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space);
void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
//void CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &tx, const double &ty, const double &tz);

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
//void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const Model_Space &Space_J, const Channels &Chan, Interactions &Ints);

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps);
void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);

void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, Interactions &Int);
double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Int);

double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql);

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff);
void EE_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff);
//void bitsetup(int *vec, unsigned long long *state, const int &size);
void bitsetup(int *vec, unsigned long long *state, const int &begin, const int &size);
//void bitsetup(int *vec, const int &begin, const int &size0, unsigned long long *state, const int &size);
//void bitsetup(const std::vector<int> &vec, std::vector<unsigned long long> &state);
double matrixe(unsigned long long *bra, unsigned long long *ket, const int &size, int *p, const int &psize, int *q, const int &qsize, const double &ME);
//double matrixe(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const std::vector<int> &p, const std::vector<int> &q, const double &ME);

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
  int Shells; //Nmax -> Shells
  std::string LevelScheme; //level scheme path
  std::string MatrixElements; //matrix elements path

  //For Excited States
  int extra; // -1 for pr, 0 for es, 1 for pa
  int Nx, Ny, Nz;
  double M, T, Par;

  Input_Parameters(){};
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
 
  State(){
    t = 0;
    m = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    ml = 0;
    n = 0;
    j = 0;
    par = 1;
    energy = -100;
    type = "none";
  };
};

//Structure for holding all model space info
struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits

  State *qnums;
  State qmins;
  State qmaxs;
  State qsizes;
  int **shellsm; // for j

  int Nmax;
  int nmax;
  int *map_2b;
  int size_2b;
  Model_Space(){};
  void delete_struct(Input_Parameters &Parameters);
};

//Structure for holding channel information
struct Channels{
  int size1;
  int size2;
  int size3;

  State *qnums1;
  State *qnums2;
  State *qnums3;
  
  int *indvec;

  int *nhh;
  int *npp;
  int *nhp;
  int *nhp1;
  int *nhp2;
  int *nh;
  int *np;
  int *nhhp;
  int *nhpp;
  int *nhhp1;
  int *nhpp1;
  int *nhh1;
  int *npp1;
  int *nhhh;
  int *nppp;

  int **hhvec;
  int **ppvec;
  int **hpvec;
  int **hp1vec;
  int **hp2vec;
  int **pvec;
  int **hvec;
  int **hhpvec;
  int **hppvec;
  int **hhp1vec;
  int **hpp1vec;
  int **hh1vec;
  int **pp1vec;
  int **hhhvec;
  int **pppvec;

  int ind0; // index of i-i cross channel for singles
  Channels(){};
  Channels(const Input_Parameters &Parameters, const Model_Space &Space);
  void delete_struct();
};

struct Doubles_1{
  int **Tmap;
  double **Evec;
  double **T1;
  double **T2;
  double **T3;
  double **T4;
  double **T5;
  double **T6;
  double **T7;
  double **T8;
  double **T9;
  double **S1;
  double **S2;
  double **S3;
  double **S4;
  double **S5;
  double **S6;
  double **S7;
  double **Q11;
  double **Q21;
  double **Q12;
  double **Q22;
  int **Qmap1;
  int **Qmap2;

  Doubles_1(){}; //default constructor
  Doubles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan); //constructor
  void delete_struct(const Channels &Chan);
  void zero(const Channels &Chan);
  void set_T(int, int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int, int) const;
};

struct Singles_1{
  int *Tmap;
  int *Tmap2;
  double *Evec;
  double *T1;
  double **T2;
  double **T3;
  double **S1;
  double **S2;
  double *S3;
  double *S4;
  double **E1;
  double **E2;
  double **E3;
  double **E4;
  double **E5;
  double **E6;
  double **E7;
  double **E8;
  double **E9;
  double **Q11;
  double **Q12;
  double **Q21;
  double **Q22;
  double *Q31;
  double **Q32;
  double *Q41;
  double **Q42;
  double **Q51;
  double **Q52;
  double **Q61;
  double **Q62;
  int **Qmap1;
  int **Qmap2;
  int *Qmap3;
  int *Qmap4;
  int **Qmap5;
  int **Qmap6;
  Singles_1(){}; //default constructor
  Singles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan); //constructor
  void delete_struct(const Channels &Chan);
  void zero(const Channels &Chan);
  void set_T(int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int) const;
};

struct Amplitudes{
  Doubles_1 D1; // for doubles only
  Singles_1 S1; // for singles part of singles
  Amplitudes(){};
  Amplitudes(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
  void delete_struct(const Input_Parameters &Parameters, const Channels &Chan);
  void zero(const Input_Parameters &Parameters, const Channels &Chan);
  double get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints);
};

struct Doubles_ME1{
  double **V1;
  double **V2;
  double **V3;
  double **V4;
  double **V5;
  double **V6;
  double **V7;
  double **V8;
  double **V9;
  double **V10;
  Doubles_ME1(){};
  Doubles_ME1(const Channels &Chan);
  void delete_struct(const Channels &Chan);
};

struct Singles_ME1{
  double **V11;
  double **V12;
  double **V15;
  double **V16;
  double **V13;
  double **V14;
  double **V19;
  double **V20;
  double **V17;
  double **V18;
  Singles_ME1(){};
  Singles_ME1(const Channels &Chan);
  void delete_struct(const Channels &Chan);
};

struct Interactions{
  Doubles_ME1 D_ME1; // for doubles only
  Singles_ME1 S_ME1; // for singles part of singles
  Interactions(){};
  Interactions(const Input_Parameters &Parameters, const Channels &Chan);
  void delete_struct(const Input_Parameters &Parameters, const Channels &Chan);
};

struct CC_Eff{
  double *X_ia1;
  double **X_ia2;
  double **X_ia3;
  int *Map_ia;

  double *X_ab1;
  double **X_ab2;
  double **X_ab3;
  int *Map_ab;

  double *X_ij1;
  double **X_ij2;
  double **X_ij3;
  double *X1_ij1;
  double **X1_ij2;
  double **X1_ij3;
  int *Map_ij;

  double *X_ai1;
  double **X_ai2;
  double **X_ai3;
  int *Map_ai;

  double **X_ijab1;

  double **X1_iabc1;
  double **X1_iabc2;
  double **X1_iabc3;
  double **X_iabc1;
  double **X_iabc3;
  double **X_iabc4;
  double **X_iabc5;
  int **Map_iabc;

  double **X1_ijka1;
  double **X1_ijka2;
  double **X_ijka1;
  double **X_ijka4;
  double **X_ijka5;
  int **Map_ijka;

  double **X1_abcd1;
  double **X1_abcd2;
  double **X1_abcd3;
  double **X_abcd1;
  double **V_abcd;
  int **Map_abcd;

  double **X_ijkl1;
  double **X_ijkl2;
  double **X_ijkl3;
  double **X_ijkl4;
  double **V_ijkl;
  int **Map_ijkl;

  double **X1_iajb1;
  double **X1_iajb2;
  double **X1_iajb3;
  double **X1_iajb4;
  double **X3_iajb1;
  double **X3_iajb2;
  double **X3_iajb3;
  double **X3_iajb5;
  double **X_iajb1;
  double **X_iajb3;
  int **Map_iajb;

  double **X_abic1;
  double **X_abic2;
  double **X_abic3;
  double **X_abic4;
  double **X_abic5;
  double **X_abic6;
  double **X_abic7;
  int **Map_abic;

  double **X2_iajk1;
  double **X2_iajk2;
  double **X2_iajk3;
  double **X2_iajk4;
  double **X2_iajk5;
  double **X2_iajk6;
  double **X2_iajk7;
  double **X_iajk1;
  int **Map_iajk;

  CC_Eff(){};
  CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
  void delete_struct(const Channels &Chan);
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
