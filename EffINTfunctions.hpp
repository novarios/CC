#ifndef EFFINTFUNCTIONS_H
#define EFFINTFUNCTIONS_H

#include <iostream>
#include <time.h>

struct Channels;
struct Amplitudes;
struct Interactions;
struct Eff_Interactions;
struct X_hh;
struct X_pp;
struct X_hp;
struct X_hhhh;
struct X_pppp;
struct X_hhpp;
struct X_pphh;
struct X_hphp;
struct X_hhhp;
struct X_hppp;
struct X_hphh;
struct X_pphp;

struct X_hh{
  double *X1_2;
  double *X_2;
  //double *X_3;
  //double *X1_3;
  //double *X_3d;
  //double *X1_3d;
  //double *X_3od;
  //double *X1_3od;

  int *X_3_index;
  //int *X_3d_index;
  int X_2_length;
  int X_3_length;
  //int X_3d_length;

  //int *map_chan;
  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;
  double *map3_fac3;

  X_hh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  //void X1_Gather();
  //void X_Gather();
  void X1_Zero(bool flag);
  void X_Zero(bool flag);
  //void Diagonalize(Channels &Chan);
  void Gather_X1_3(double *X3, int length, int offset);
  void Gather_X_3(double *X3, int length, int offset);
  void Set_X1_3(double *X3, int length, int offset);
  void Set_X_3(double *X3, int length, int offset);
  void Set_X1_3_OD(double *X3, int length, int offset);
  void Set_X_3_OD(double *X3, int length, int offset);
};

struct X_pp{
  double *X_2;
  //double *X_3;
  //double *X_3d;
  //double *X_3od;

  int *X_3_index;
  //int *X_3d_index;
  int X_2_length;
  int X_3_length;
  //int X_3d_length;

  int *map3_chan;
  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;
  double *map3_fac3;

  X_pp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  //void X_Gather();
  void X_Zero(bool flag);
  //void Diagonalize(Channels &Chan);
  void Gather_X_3(double *X3, int length, int offset);
  void Set_X_3(double *X3, int length, int offset);
  void Set_X_3_OD(double *X3, int length, int offset);
};

struct X_hp{
  double *X_2;
  //double *X_3;

  int *X_3_index;
  int X_2_length;
  int X_3_length;

  //int *map_chan;
  int *map3_ind;
  double *map3_fac1;
  double *map3_fac2;

  X_hp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  //void X_Gather();
  void X_Zero(bool flag);
  void Gather_X_3(double *X3, int length, int offset);
  void Set_X_3(double *X3, int length, int offset);
};

struct X_hhhh{
  double *X_1;

  int *X_1_index;
  int *X_3_3_index;
  int *X_3_4_index;

  int X_1_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map33_ind;
  double *map33_fac1;

  int *map34_ind;
  double *map34_fac1;

  X_hhhh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);

  void Gather_X_3_3(double *X3, int length, int offset);
  void Gather_X_3_4(double *X3, int length, int offset);
};

struct X_pppp{
  double *X1_1;
  double *X_1;

  int *X_1_index;
  int *X_3_1_index;
  int *X_3_2_index;

  int X_1_length;
  int X_3_1_length;
  int X_3_2_length;

  int *map31_ind;
  double *map31_fac1;

  int *map32_ind;
  double *map32_fac1;

  X_pppp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);

  void Gather_X1_3_1(double *X3, int length, int offset);
  void Gather_X1_3_2(double *X3, int length, int offset);
};

struct X_hphp{
  double *X_2_1;
  double *X1_2_1;
  double *X2_2_1;
  double *X3_2_1;

  int *X_2_1_index;
  int *X_2_2_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int *X_3_3_index;
  int *X_3_4_index;

  int X_2_1_length;
  int X_2_2_length;
  int X_3_1_length;
  int X_3_2_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map22_ind;
  double *map22_fac2;

  int *map31_index;
  int *map31_num;
  int *map31_ind;
  double *map31_fac2;

  int *map32_index;
  int *map32_num;
  int *map32_ind;
  double *map32_fac1;

  int *map33_index;
  int *map33_num;
  int *map33_ind;
  double *map33_fac1;

  int *map34_index;
  int *map34_num;
  int *map34_ind;
  double *map34_fac2;

  X_hphp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);
  void X2_Zero(bool flag);
  void X3_Zero(bool flag);

  void Set_X_2_2(double *X2, int length, int offset);
  void Gather_X1_3_2(double *X3, int length, int offset);
  void Gather_X1_3_3(double *X3, int length, int offset);
  void Set_X1_3_1(double *X3, int length, int offset);
  void Gather_X2_3_2(double *X3, int length, int offset);
  void Gather_X2_3_3(double *X3, int length, int offset);
  void Set_X2_3_1(double *X3, int length, int offset);
  void Gather_X3_3_2(double *X3, int length, int offset);
  void Gather_X3_3_3(double *X3, int length, int offset);
  void Set_X3_3_4(double *X3, int length, int offset);
  void Gather_X_3_2(double *X3, int length, int offset);
  void Gather_X_3_3(double *X3, int length, int offset);
};

struct X_hhhp{
  double *X_3_3;
  double *X1_3_3;

  int *X_1_index;
  int *X_2_1_index;
  int *X_2_3_index;
  int *X_3_3_index;
  int *X_3_4_index;

  int X_1_length;
  int X_2_1_length;
  int X_2_3_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map1_ind;
  double *map1_fac2;

  int *map21_index;
  int *map21_num;
  int *map21_ind;
  double *map21_fac2;

  int *map23_index;
  int *map23_num;
  int *map23_ind;
  double *map23_fac2;

  int *map34_ind;
  double *map34_fac2;

  X_hhhp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);

  void Set_X1_3_4(double *X3, int length, int offset);
  void Set_X_2_1(double *X2, int length, int offset);
  void Set_X_2_3(double *X2, int length, int offset);
  void Set_X_1(double *X1, int length, int offset);
};

struct X_hppp{
  double *X_3_2;
  double *X1_3_2;

  int *X_1_index;
  int *X_2_3_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int *X_3_3_index;

  int X_1_length;
  int X_2_3_length;
  int X_3_1_length;
  int X_3_2_length;
  int X_3_3_length;

  int *map1_ind;
  double *map1_fac2;

  int *map23_index;
  int *map23_num;
  int *map23_ind;
  double *map23_fac2;

  int *map31_ind;
  double *map31_fac2;

  int *map33_ind;
  double *map33_fac2;

  X_hppp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);

  void Set_X1_3_1(double *X3, int length, int offset);
  void Set_X1_3_3(double *X3, int length, int offset);
  void Set_X_3_3(double *X3, int length, int offset);
  void Set_X_2_3(double *X2, int length, int offset);
  void Set_X_1(double *X1, int length, int offset);
};

struct X_hphh{
  double *X_3_2;
  double *X1_3_2;

  int *X_1_index;
  int *X_2_1_index;
  int *X_2_3_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int *X_3_3_index;
  int *X_3_4_index;

  int X_1_length;
  int X_2_1_length;
  int X_2_3_length;
  int X_3_1_length;
  int X_3_2_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map1_ind;
  double *map1_fac1;

  int *map21_index;
  int *map21_num;
  int *map21_ind;
  double *map21_fac1;

  int *map23_index;
  int *map23_ind;
  int *map23_num;
  double *map23_fac1;

  int *map31_ind;
  double *map31_fac1;
  double *map31_fac2;

  int *map33_ind;
  double *map33_fac1;

  int *map34_ind;
  double *map34_fac1;

  X_hphh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);

  void Set_X1_3_1(double *X3, int length, int offset);
  void Set_X_3_1(double *X3, int length, int offset);
  void Gather_X_3_1(double *X3, int length, int offset);
  void Gather_X_3_3(double *X3, int length, int offset);
  void Gather_X_3_4(double *X3, int length, int offset);
  void Gather_X_2_1(double *X2, int length, int offset);
  void Gather_X_2_3(double *X2, int length, int offset);
  void Gather_X_1(double *X1, int length, int offset);
};

struct X_pphp{
  double *X_3_3;
  double *X1_3_3;

  int *X_1_index;
  int *X_2_2_index;
  int *X_2_3_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int *X_3_3_index;
  int *X_3_4_index;

  int X_1_length;
  int X_2_2_length;
  int X_2_3_length;
  int X_3_1_length;
  int X_3_2_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map1_ind;
  double *map1_fac1;

  int *map22_index;
  int *map22_num;
  int *map22_ind;
  double *map22_fac1;

  int *map23_index;
  int *map23_num;
  int *map23_ind;
  double *map23_fac1;

  int *map31_ind;
  double *map31_fac1;

  int *map32_ind;
  double *map32_fac1;

  int *map34_ind;
  double *map34_fac1;
  double *map34_fac2;

  X_pphp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
  void X_Zero(bool flag);
  void X1_Zero(bool flag);

  void Set_X1_3_4(double *X3, int length, int offset);
  void Set_X_3_4(double *X3, int length, int offset);
  void Gather_X1_3_1(double *X3, int length, int offset);
  void Gather_X1_3_2(double *X3, int length, int offset);
  void Gather_X_3_1(double *X3, int length, int offset);
  void Gather_X_3_2(double *X3, int length, int offset);
  void Gather_X_3_4(double *X3, int length, int offset);
  void Gather_X_2_2(double *X2, int length, int offset);
  void Gather_X_2_3(double *X2, int length, int offset);
  void Gather_X_1(double *X1, int length, int offset);
};

struct Eff_Interactions{
  X_hh Xhh;
  X_pp Xpp;
  X_hp Xhp;
  X_hhhh Xhhhh;
  X_pppp Xpppp;
  X_hphp Xhphp;
  X_hhhp Xhhhp;
  X_hppp Xhppp;
  X_hphh Xhphh;
  X_pphp Xpphp;
  Eff_Interactions(){}; //default constructor
  void Build(Channels &Chan);
  void Update_1(Channels &Chan, Interactions &Ints, Amplitudes &Amps);
  void Update_2(Channels &Chan, Interactions &Ints, Amplitudes &Amps);
  void Update_3(Channels &Chan, Interactions &Ints, Amplitudes &Amps);
  void Delete();
  void Zero(bool flag);
  void Diagonalize(Channels &Chan);
};

#endif


