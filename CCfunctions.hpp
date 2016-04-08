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
struct Amplitudes;
struct Interactions;

struct  Doubles_1;
struct  Singles_1;
struct  Doubles_2;
struct  Singles_2;
struct  Doubles_3;

struct  Doubles_ME1;
struct  Singles_ME1;
struct  Singles_ME2;

//struct CCD;
//struct CCD_2; // for singles
//struct CCS;
//struct CC_Matrix_Elements;
//struct CCS_Matrix_Elements; // for singles

//struct CC_Eff;
//struct V_Conv;

int Index0(const std::vector<int> &vec1, const int &p, const int &q);
int Index(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index2(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int Index3(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s);
int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const int &ind1, const int &ind2);
int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const int &ind1, const int &ind2);

void Setup_Channels(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan);
void Get_Input_Parameters(std::string &infile, Input_Parameters &Parameters);
void Print_Parameters(const Input_Parameters &Parameters);
void Setup_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps);
void Setup_Ints(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints);

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void Build_Model_Space_J1(Input_Parameters &Parameters, Model_Space_J &Space_J);
void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space_J &Space_J, Model_Space &Space);
void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void EG_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
//void CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &tx, const double &ty, const double &tz);

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Coulomb_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const Model_Space_J &Space_J, const Channels &Chan, Interactions &Ints);

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps);
void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);

void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Int);
double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Int);

double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
double Coulomb(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);

/*CC_Eff Build_CC_Eff(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, CCD &CC, const Channels &Chan);
void EE_EOM(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan);
void EE_EOM_HO(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan);
std::vector<unsigned long long> bitsetup(const int &i, const int &j, const int &a, const int &b, const int &Nbit);
double matrixe1(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const double &ME, const int &test);
double matrixe2(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &r, const int &s, const double &ME);
double matrixe3(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &r, const int &t, const int &u, const int &v, const double &ME);
double vint_Conversion(const Model_Space &Space, const std::vector<std::vector<double> > &ME, const int &qi, const int &qj, const int &qk, const int &ql);*/


//Structure for holding Input parameters
struct Input_Parameters{
  std::string calc_case; //nuclear or electronic
  std::string basis; //infinite, finite (HO), or finite_J (HO_J)
  std::string approx; //doubles, singles, or triples
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
  int Pmin, Pmax, Psize, P2size;
  
  //for CART
  std::vector<int> levelsnx; //list of single particle state x-momeuntum quantum number
  std::vector<int> levelsny; //list of single particle state y-momeuntum quantum number
  std::vector<int> levelsnz; //list of single particle state z-momeuntum quantum number
  int nmax;
  int Nxmin, Nxmax, Nxsize, Nx2size;
  int Nymin, Nymax, Nysize, Ny2size;
  int Nzmin, Nzmax, Nzsize, Nz2size;

  int Mmin, Mmax, Msize, M2size;
  int Tmin, Tmax, Tsize, T2size;
  std::vector<int> map_2b_dir;
  std::vector<int> map_2b_cross;
  int Chansize_2b_dir;
  int Chansize_2b_cross;
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
  std::vector<int> hp;
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
  std::vector<std::vector<int> > hpvec1;
  std::vector<std::vector<int> > hp1vec1;
  std::vector<std::vector<int> > hp2vec1;
  std::vector<std::vector<int> > pvec1;
  std::vector<std::vector<int> > hhpvec1;
  std::vector<std::vector<int> > hvec1;
  std::vector<std::vector<int> > hppvec1;

  int ind0; // index of i-i cross channel for singles
};

struct Doubles_1{
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
  Doubles_1(); //default constructor
  Doubles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  double get_T(int, int) const;
};

struct Singles_1{
  std::vector<std::vector<int> > Tmap;
  std::vector<int> Tmap2;
  std::vector<std::vector<int> > Tmap3;
  std::vector<std::vector<int> > Tmap4;
  std::vector<std::vector<double> > Evec;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > S1;
  std::vector<std::vector<double> > S2;
  std::vector<std::vector<double> > E1;
  std::vector<std::vector<double> > E2;
  std::vector<double> T3;
  std::vector<double> S3;
  std::vector<double> S4;
  std::vector<double> E3;
  Singles_1(); //default constructor
  Singles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  void set_T_2(const Channels &Chan);
  double get_T(int, int) const;
};

struct Doubles_2{
  std::vector<std::vector<int> > Tmap;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<double> T3;
  Doubles_2(); //default constructor
  Doubles_2(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(const Channels &Chan, const int &i, const int &j, const double &T);
  double get_T(int, int) const;
};

struct Singles_2{
  std::vector<std::vector<int> > Tmap;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > T3;
  Singles_2(); //default constructor
  Singles_2(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  double get_T(int, int) const;
};

struct Doubles_3{
  std::vector<std::vector<int> > Tmap;
  std::vector<std::vector<double> > T1;
  std::vector<std::vector<double> > T2;
  std::vector<std::vector<double> > T3;
  Doubles_3(); //default constructor
  Doubles_3(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space); //constructor
  void set_T(int, int, double);
  double get_T(int, int) const;
};

struct Amplitudes{
  Doubles_1 D1; // for doubles only
  Singles_1 S1; // for singles part of singles
  Doubles_2 D2; // for doubles part of singles
  //Singles_2 S2; // for singles part of doubles
  //Doubles_3 D3; // for second doubles of doubles
  double get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints);
};

struct Doubles_ME1{
  std::vector<std::vector<double> > PPPP;
  std::vector<std::vector<double> > HHHH;
  std::vector<std::vector<double> > HPHP1;
  std::vector<std::vector<double> > HPHP2;
  std::vector<std::vector<double> > HHPP1;
  std::vector<std::vector<double> > HHPP2;
  std::vector<std::vector<double> > HHPP3;
  std::vector<std::vector<double> > HHPP4;
  std::vector<std::vector<double> > HHPP4T;
  Doubles_ME1();
  Doubles_ME1(const Channels &Chan);
};

struct Singles_ME1{
  std::vector<std::vector<double> > HHPP1;
  std::vector<std::vector<double> > HHPP2;
  std::vector<double> HHPP3;
  std::vector<std::vector<double> > HPPP;
  std::vector<std::vector<double> > HHHP;
  std::vector<double> HPHP;
  Singles_ME1();
  Singles_ME1(const Channels &Chan);
};

struct Singles_ME2{
  std::vector<std::vector<double> > HHPP1;
  std::vector<std::vector<double> > HHPP2;
  std::vector<double> HHPP3;
  std::vector<std::vector<double> > HPPP;
  std::vector<std::vector<double> > HHHP;
  Singles_ME2();
  Singles_ME2(const Channels &Chan);
};

struct Interactions{
  Doubles_ME1 D_ME1; // for doubles only
  Singles_ME1 S_ME1; // for singles part of singles
  Singles_ME2 S_ME2; // for singles part of doubles
};

/*struct CC_Eff{
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
