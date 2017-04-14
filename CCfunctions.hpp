#ifndef CCFUNCTIONS_H
#define CCFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <cstring>
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
#include <unordered_map>

const std::string PATH = "inputs/";

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
//#define RM_dgemm(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_B, transf_A, n, m, k, alpha, B, n, A, k, beta, C, n)
//#define RMT_dgemm(A, B, C, m, n, k, alpha, beta, transf_A, transf_B) dgemm_(transf_B, transf_A, n, m, k, alpha, B, k, A, k, beta, C, n)

#define min(a,b) (a <= b ? a : b)
#define max(a,b) (a >= b ? a : b)

struct Input_Parameters;
struct State;
struct Model_Space;
struct Channels;

struct Amplitudes;
struct Doubles_1;
struct Singles_1;
struct CC_Eff;

struct HF_Channels;
struct HF_Matrix_Elements;

struct Interactions;

int Hash2(const int &p, const int &q, const int &size);
int Hash3(const int &p, const int &q, const int &r, const int &size);

int Index11(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q);
int Index2(const int *vec1, const int &num1, const int &p, const int &q);
int Index1(const int *vec1, const int &num1, const int &p);
int Index22(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index13(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index31(const int *vec1, const int *vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);

void Get_Input_Parameters(std::string &infile, Input_Parameters &Parameters);
void Print_Parameters(const Input_Parameters &Parameters, const Model_Space &Space);

//void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME);
void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps);
void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_J(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_2_J(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step_J(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);

void Random_Step(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, double &width, double &error2);
void Randomize_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, double &width);
void Gather_Amps0(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix);
void Gather_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix);
void Gather_Amps2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix, double &checkdot, int &N, double ***delp);
void Gather_Amps3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix, int &N, double ***p, double ***delp, double *B);
void CC_Error(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error);
void Print_Amps(const Input_Parameters &Parameters, const Channels &Chan, Amplitudes &Amps);
void Doubles_Step_explicit(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2, double &error);
void Doubles_Step_explicit2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps1, Amplitudes &Amps2, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, double &error);

void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, Interactions &Int);
double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Int);

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff);
void EE_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff);
void PA_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums);
void PR_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums);
void CC_compare_JM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, std::string &inputfile);

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

struct Doubles_1{
  //int **Tmap;
  //int **TJnum;
  //int ****TJmap;
  int **Tnum;
  int ****Tmap;
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
  Doubles_1(const Channels &Chan, const Doubles_1 &D1); //constructor
  Doubles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan); //constructor
  void copy_Doubles_1(const Channels &Chan, const Doubles_1 &D1);
  void delete_struct(const Channels &Chan);
  void zero(const Channels &Chan);
  void zero1(const Channels &Chan);
  void set_T(int, int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int, int) const;
  void set_TJ(const Model_Space &Space, const Channels &Chan, int &chan, int &hhpp, int &i, int &j, int &a, int &b, double T);
  double get_TJ(const Model_Space &Space, const Channels &Chan, int &chan, int &hhpp, int &i, int &j, int &a, int &b) const;
  void set_T_2J(const Model_Space &Space, const Channels &Chan, Interactions &Ints);
};

struct Singles_1{
  int *Tmap;
  //int *Tmap2;
  //int *TJnum;
  //int **TJmap;
  int *Tnum2;
  int **Tmap2;
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
  int **Qnum1;
  int **Qnum2;
  int ***Qmap1;
  int ***Qmap2;
  //int **QJnum1;
  //int **QJnum2;
  //int ***QJmap1;
  //int ***QJmap2;
  int *Qmap3;
  int *Qmap4;
  int **Qmap5;
  int **Qmap6;
  Singles_1(){}; //default constructor
  Singles_1(const Channels &Chan, const Singles_1 &S1); //constructor
  Singles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan); //constructor
  void copy_Singles_1(const Channels &Chan, const Singles_1 &S1);
  void delete_struct(const Channels &Chan);
  void zero(const Channels &Chan);
  void zero1(const Channels &Chan);
  void set_T(int, double);
  void set_T_2(const Channels &Chan, Interactions &Ints);
  double get_T(int) const;
  void set_TJ(const Model_Space &Space, int &hp, int &i, int &a, double T);
  void set_T_2J(const Model_Space &Space, const Channels &Chan, Interactions &Ints);
  double get_TJ(const Model_Space &Space, int &hp, int &i, int &a) const;
};

struct Amplitudes{
  Doubles_1 D1; // for doubles only
  Singles_1 S1; // for singles part of singles
  Amplitudes(){};
  Amplitudes(const Input_Parameters &Parameters, const Channels &Chan, const Amplitudes &Amps);
  Amplitudes(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan);
  void copy_Amplitudes(const Input_Parameters &Parameters, const Channels &Chan, const Amplitudes &Amps);
  void delete_struct(const Input_Parameters &Parameters, const Channels &Chan);
  void zero(const Input_Parameters &Parameters, const Channels &Chan);
  void zero1(const Input_Parameters &Parameters, const Channels &Chan);
  double get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints);
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
