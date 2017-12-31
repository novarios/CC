#ifndef LAMBDAFUNCTIONS_H
#define LAMBDAFUNCTIONS_H

#include <string>

struct Channels;
struct LAmplitudes;
struct Amplitudes;
struct Doubles;
struct Singles;
struct Interactions;
struct Eff_Interactions;

void Update_LAmps(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps, LAmplitudes &tempLAmps);
void LDoubles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps1, LAmplitudes &LAmps2);
void LSingles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps1, LAmplitudes &LAmps2);

void LCC_Algorithm(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps);

void Random_LStep(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps0, LAmplitudes &LAmps, LAmplitudes &LAmps2, LAmplitudes &tempLAmps, double &mix, double &width, double &error, double &error2);
void Randomize_LAmps(LAmplitudes &LAmps0, LAmplitudes &LAmps, double &width);
void Print_LAmps(Channels &Chan, LAmplitudes &LAmps);
void Gather_LAmps(Channels &Chan, Eff_Interactions &Eff_Ints, LAmplitudes &LAmps, LAmplitudes &LAmps2, double &mix);
void LCC_Error(Interactions &Ints, LAmplitudes &LAmps, LAmplitudes &LAmps2, double &error);

struct LAmplitudes{
  int *Evec_ind;
  double *Evec_fac;
  double *L1;

  int *L1_index;
  int *L2_1_index;
  int *L2_2_index;
  int *L2_3_index;
  int *L2_4_index;
  int *L3_1_index;
  int *L3_2_index;
  int *L3_3_index;
  int *L3_4_index;

  int L1_length;
  int L2_1_length;
  int L2_2_length;
  int L2_3_length;
  int L2_4_length;
  int L3_1_length;
  int L3_2_length;
  int L3_3_length;
  int L3_4_length;

  int *map21_index;
  int *map21_num;
  int *map21_chan;
  int *map21_ind;
  double *map21_fac1;
  double *map21_fac2;

  int *map22_index;
  int *map22_num;
  int *map22_chan;
  int *map22_ind;
  double *map22_fac1;

  int *map23_index;
  int *map23_num;
  int *map23_chan;
  int *map23_ind;
  double *map23_fac1;

  int *map24_index;
  int *map24_num;
  int *map24_chan;
  int *map24_ind;
  double *map24_fac1;

  int *map31_chan;
  int *map31_ind;
  double *map31_fac1;
  double *map31_fac2;

  int *map32_chan;
  int *map32_ind;
  double *map32_fac1;
  double *map32_fac2;

  int *map33_chan;
  int *map33_ind;
  double *map33_fac1;
  double *map33_fac2;

  int *map34_chan;
  int *map34_ind;
  double *map34_fac1;
  double *map34_fac2;

  /*LDoubles(){}; //default constructor
  void Build(Channels &Chan); //constructor
  void Build(Channels &Chan, LDoubles &LD1); //copy constructor
  void Copy(LDoubles &LD1);
  void Copy(LDoubles &LD1, double *vec);
  void Delete();
  void Zero(bool flag);
  void Set_L(int ind, double L);
  void Set_L(double *vec);
  double Get_L(int ind);*/
  void Set_L2_1(double *L2, int length, int offset);
  void Gather_L2_1(double *L2, int length, int offset);
  void Gather_L2_2(double *L2, int length, int offset);
  void Gather_L2_3(double *L2, int length, int offset);
  void Gather_L2_4(double *L2, int length, int offset);
  void Set_L3_1(double *L3, int length, int offset);
  void Set_L3_2(double *L3, int length, int offset);
  void Set_L3_3(double *L3, int length, int offset);
  void Set_L3_4(double *L3, int length, int offset);
  void Gather_L3_1(double *L3, int length, int offset);
  void Gather_L3_2(double *L3, int length, int offset);
  void Gather_L3_3(double *L3, int length, int offset);
  void Gather_L3_4(double *L3, int length, int offset);

  int *evec_ind;
  double *evec_fac;
  double *l2;
  
  int *l3_index;
  int l2_length;
  int l3_length;

  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;

  /*LSingles(){}; //default constructor
  void Build(Channels &Chan); //constructor
  void Build(Channels &Chan, LSingles &LS1); //copy constructor
  void Copy(LSingles &LS1);
  void Copy(LSingles &LS1, double *vec, int offset);
  void Delete();
  void Zero(bool flag);
  void Set_L(int ind, double t);
  void Set_L(double *vec, int offset);
  double Get_L(int ind);*/
  void Gather_l3(double *l3, int length, int offset);
  void Set_l3(double *l3, int length, int offset);

  LAmplitudes(){};
  void Build(Channels &Chan);
  void Build(Channels &Chan, LAmplitudes &LAmps);
  void Copy(LAmplitudes &LAmps);
  void Copy(LAmplitudes &LAmps, double *vec);
  void Delete();
  //void Zero(bool t1flag);
  void Zero();
  void Set_T(double *vec);
};

/*struct LSingles{
  int *evec_ind;
  double *evec_fac;
  double *l2;
  
  int *l3_index;
  int l2_length;
  int l3_length;

  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;

  LSingles(){}; //default constructor
  void Build(Channels &Chan); //constructor
  void Build(Channels &Chan, LSingles &LS1); //copy constructor
  void Copy(LSingles &LS1);
  void Copy(LSingles &LS1, double *vec, int offset);
  void Delete();
  void Zero(bool flag);
  void Set_L(int ind, double t);
  void Set_L(double *vec, int offset);
  double Get_L(int ind);
  void Gather_l3(double *l3, int length, int offset);
  void Set_l3(double *l3, int length, int offset);
};

struct LAmplitudes{
  LDoubles LD1; // for doubles only
  LSingles LS1; // for singles part of singles
  LAmplitudes(){};
  void Build(Channels &Chan);
  void Build(Channels &Chan, LAmplitudes &LAmps);
  void Copy(LAmplitudes &LAmps);
  void Copy(LAmplitudes &LAmps, double *vec);
  void Delete();
  void Zero(bool t1flag);
  void Set_T(double *vec);
  };*/

#endif
