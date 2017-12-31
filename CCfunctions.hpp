#ifndef CCFUNCTIONS_H
#define CCFUNCTIONS_H

#include <string>

struct Channels;
struct Amplitudes;
struct Doubles;
struct Singles;
struct Interactions;
struct Eff_Interactions;

void Update_Amps(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &tempAmps);
void Doubles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Singles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);
void Doubles_Step_2(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2);

void CC_Algorithm(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);

void Random_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, double &width, double &error, double &error2);
void Randomize_Amps(Amplitudes &Amps0, Amplitudes &Amps, double &width);
void Print_Amps(Channels &Chan, Amplitudes &Amps);
void Gather_Amps(Channels &Chan, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix);
void CC_Error(Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error);

void Initialize_DIIS(Amplitudes &Amps, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl);
void Delete_DIIS();
void Perform_DIIS(Amplitudes &Amps, Amplitudes &Amps0, double &mix, double *&p, double *&delp, double *&tempdelp, double *&B, int &N, int &maxl, int &DIIS_count);
void Update_B1(Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);
void Update_B2(Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);

struct Doubles{
  //int *Evec_chan;
  int *Evec_ind;
  double *Evec_fac;
  double *T1;

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

  int *map21_index;
  int *map21_num;
  int *map21_chan;
  int *map21_ind;
  double *map21_fac1;
  double *map21_fac2;

  int *map22_index;
  int *map22_num;
  int *map22_ind;
  double *map22_fac1;

  int *map23_index;
  int *map23_num;
  int *map23_ind;
  double *map23_fac1;
  double *map23_fac2;

  int *map24_index;
  int *map24_num;
  int *map24_ind;
  double *map24_fac1;

  int *map31_ind;
  double *map31_fac1;
  double *map31_fac2;

  int *map32_ind;
  double *map32_fac1;
  double *map32_fac2;

  int *map33_ind;
  double *map33_fac1;
  double *map33_fac2;

  int *map34_ind;
  double *map34_fac1;
  double *map34_fac2;

  Doubles(){}; //default constructor
  void Build(Channels &Chan); //constructor
  void Build(Channels &Chan, Doubles &D1); //copy constructor
  void Copy(Doubles &D1);
  void Copy(Doubles &D1, double *vec);
  void Delete();
  void Zero(bool flag);
  void Set_T(int ind, double T);
  void Set_T(double *vec);
  double Get_T(int ind);

  void Set_T2_1(double *T2, int length, int offset);
  void Set_T2_3(double *T2, int length, int offset);
  void Gather_T2_1(double *T2, int length, int offset);
  void Gather_T2_2(double *T2, int length, int offset);
  void Gather_T2_3(double *T2, int length, int offset);
  void Gather_T2_4(double *T2, int length, int offset);
  void Set_T3_1(double *T3, int length, int offset);
  void Set_T3_2(double *T3, int length, int offset);
  void Set_T3_3(double *T3, int length, int offset);
  void Set_T3_4(double *T3, int length, int offset);
  void Gather_T3_1(double *T3, int length, int offset);
  void Gather_T3_2(double *T3, int length, int offset);
  void Gather_T3_3(double *T3, int length, int offset);
  void Gather_T3_4(double *T3, int length, int offset);
};

struct Singles{
  //int *evec_chan;
  int *evec_ind;
  double *evec_fac;
  double *t2;
  //double *t3;
  
  int *t3_index;
  int t2_length;
  int t3_length;

  //int *map_chan;
  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;

  Singles(){}; //default constructor
  void Build(Channels &Chan); //constructor
  void Build(Channels &Chan, Singles &S1); //copy constructor
  void Copy(Singles &S1);
  void Copy(Singles &S1, double *vec, int offset);
  void Delete();
  void Zero(bool flag);
  void Set_T(int ind, double t);
  void Set_T(double *vec, int offset);
  double Get_T(int ind);
  void Gather_t3(double *t3, int length, int offset);
  void Set_t3(double *t3, int length, int offset);
};

struct Amplitudes{
  Doubles D1; // for doubles only
  Singles S1; // for singles part of singles
  Amplitudes(){};
  void Build(Channels &Chan);
  void Build(Channels &Chan, Amplitudes &Amps);
  void Copy(Amplitudes &Amps);
  void Copy(Amplitudes &Amps, double *vec);
  void Delete();
  void Zero(bool t1flag);
  void Set_T(double *vec);
  
  double dE;
  void get_dE(Channels &Chan, Interactions &Ints);
};

#endif
