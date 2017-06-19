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

struct Input_Parameters;
struct State;
struct Model_Space;
struct Channels;

struct Amplitudes;
struct Doubles_1;
struct Singles_1;

struct HF_Channels;
struct HF_Matrix_Elements;

struct Interactions;
struct Eff_Interactions;

int Hash2(int &p, int &q, int &size);
int Hash3(int &p, int &q, int &r, int &size);

int Index11(int *vec1, int *vec2, int &num1, int &num2, int &p, int &q);
int Index2(int *vec1, int &num1, int &p, int &q);
int Index1(int *vec1, int &num1, int &p);
int Index22(int *vec1, int *vec2, int &num1, int &num2, int &p, int &q, int &r, int &s);
int Index13(int *vec1, int *vec2, int &num1, int &num2, int &p, int &q, int &r, int &s);
int Index31(int *vec1, int *vec2, int &num1, int &num2, int &p, int &q, int &r, int &s);

void Get_Input_Parameters(std::string &infile, Input_Parameters &Parameters);
void Print_Parameters(Input_Parameters &Parameters, Model_Space &Space);

//void Perform_CC(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME);
void Update_Heff_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);
void Update_Heff_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);
void Update_Heff_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);
void Doubles_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Doubles_Step_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);

void Perform_CC(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);
/*void Doubles_Step_2(Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step(Model_Space &Space, Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_J(Model_Space &Space, Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);
void Doubles_Step_2_J(Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step_J(Model_Space &Space, Channels &Chan, Interactions &Int, Amplitudes &Amp1, Amplitudes &Amp2);*/

void Random_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, double &width, double &error, double &error2);
void Randomize_Amps(Input_Parameters &Parameters, Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, double &width);
void Print_Amps(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps);
void Gather_Amps(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix);
void CC_Error(Input_Parameters &Parameters, Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error);

void HF(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Int);
double E_Ref(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Int);

void Initialize_DIIS(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl);
void Delete_DIIS(Input_Parameters &Parameters, Channels &Chan, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl);
void Perform_DIIS(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, Amplitudes &Amps0, double &mix, double *&p, double *&delp, double *&tempdelp, double *&B, int &N, int &maxl, int &DIIS_count);
void Update_B1(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);
void Update_B2(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);

void PA_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, State *states, double *nums);
void PR_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, State *states, double *nums);
void CC_compare_JM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, std::string &inputfile);


void Map_4_count(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int type, int offset, int &length, int &count, int &p, int &q, int &r, int &s);
void Map_4(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int type, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, int &p, int &q, int &r, int &s, double &J);
void Map_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, int &p, int &q);
void direct_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &q, int &r, int &s, int &jmin1, State &tb1);
void cross_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &s, int &r, int &q, int &jmin2, State &tb2);

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
  int *Evec_chan;
  int *Evec_ind;
  double *T1;
  double *T2_1;
  double *T2_2;
  double *T2_3;
  double *T2_4;
  double *T3_1;
  double *T3_2;
  double *T3_3;
  double *T3_4;

  int *T1_index;
  int *T2_1_index;
  int *T2_2_index;
  int *T2_3_index;
  int *T2_4_index;
  int *T3_1_index;
  int *T3_2_index;
  int *T3_3_index;
  int *T3_4_index;
  int T1_length;
  int T2_1_length;
  int T2_2_length;
  int T2_3_length;
  int T2_4_length;
  int T3_1_length;
  int T3_2_length;
  int T3_3_length;
  int T3_4_length;

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  Doubles_1(){}; //default constructor
  Doubles_1(Channels &Chan, Doubles_1 &D1); //constructor
  Doubles_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan); //constructor
  void copy_Doubles_1(Channels &Chan, Doubles_1 &D1);
  void delete_struct(Channels &Chan);
  void zero(Channels &Chan, bool flag);
  void set_T(int &ind, double &T);
  double get_T(int &ind);
};

struct Singles_1{
  int *evec_chan;
  int *evec_ind;
  double *t2;
  double *t3;
  
  int *t3_index;
  int t2_length;
  int t3_length;

  int *map_chan;
  int *map_ind;
  double *map_fac1;
  double *map_fac2;

  Singles_1(){}; //default constructor
  Singles_1(Channels &Chan, Singles_1 &S1); //constructor
  Singles_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan); //constructor
  void copy_Singles_1(Channels &Chan, Singles_1 &S1);
  void delete_struct(Channels &Chan);
  void zero(Channels &Chan, bool flag);
  void set_T(int &ind, double &t);
  double get_T(int &ind);
};

struct Amplitudes{
  Doubles_1 D1; // for doubles only
  Singles_1 S1; // for singles part of singles
  Amplitudes(){};
  Amplitudes(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps);
  Amplitudes(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void copy_Amplitudes(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void zero(Input_Parameters &Parameters, Channels &Chan, bool t1flag);
  double get_energy(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints);
};

#endif
