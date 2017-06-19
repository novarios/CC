#ifndef INTFUNCTIONS_H
#define INTFUNCTIONS_H

struct Input_Parameters;
struct Model_Space;
struct Channels;
struct HF_Channels;
struct HF_Matrix_Elements;

struct Interactions;
struct F_matrix;
struct V_hhhh;
struct V_pppp;
struct V_hhpp;
struct V_pphh;
struct V_hphp;
struct V_hhhp;
struct V_hppp;
struct V_hphh;
struct V_pphp;
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

void Read_Matrix_Elements_J(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_M(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_QD(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME);
void Read_QD_ME_From_File(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME);

void Coulomb_Inf_Matrix_Elements(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints);
double Coulomb_Inf(Model_Space &Space, int &qi, int &qj, int &qk, int &ql, double &L);
double Coulomb_HO(Input_Parameters &Parameters, Model_Space &Space, int &qi, int &qj, int &qk, int &ql);

void Minnesota_Matrix_Elements(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints);
int kron_del(int &i, int &j);
int spinExchangeMtxEle(int &i, int &j, int &k, int &l);
double vint_Minnesota_Momentum(Model_Space &Space, int &qi, int &qj, int &qk, int &ql, double &L);

void Get_Fock_Matrix(Input_Parameters &Paramters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &HF, Model_Space &Space, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_J(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_JM(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints);

struct F_matrix{
  double *hh_3;
  int *hh_3_index;
  int hh_3_length;
  double *hp_3;
  int *hp_3_index;
  int hp_3_length;
  double *ph_3;
  int *ph_3_index;
  int ph_3_length;
  double *pp_3;
  int *pp_3_index;
  int pp_3_length;
  double *hh_2;
  int hh_2_length;
  double *hp_2;
  int hp_2_length;
  double *ph_2;
  int ph_2_length;
  double *pp_2;
  int pp_2_length;
  F_matrix(){}; //default constructor
  F_matrix(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hhhh{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hhhh(){}; //default constructor
  V_hhhh(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_pppp{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  double *V_3_3;
  int *V_3_3_index;
  int V_3_3_length;
  V_pppp(){}; //default constructor
  V_pppp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hhpp{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  double *V_2_1;
  int *V_2_1_index;
  int V_2_1_length;
  double *V_3_1;
  int *V_3_1_index;
  int V_3_1_length;
  double *V_3_3;
  int *V_3_3_index;
  int V_3_3_length;
  double *V_2_3;
  int *V_2_3_index;
  int V_2_3_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hhpp(){}; //default constructor
  V_hhpp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_pphh{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  V_pphh(){}; //default constructor
  V_pphh(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hphp{
  double *V_2_1;
  int *V_2_1_index;
  int V_2_1_length;
  double *V_2_2;
  int *V_2_2_index;
  int V_2_2_length;
  V_hphp(){}; //default constructor
  V_hphp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hhhp{
  double *V_2_3;
  int *V_2_3_index;
  int V_2_3_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  double *V_3_3;
  int *V_3_3_index;
  int V_3_3_length;
  V_hhhp(){}; //default constructor
  V_hhhp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hppp{
  double *V_2_4;
  int *V_2_4_index;
  int V_2_4_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hppp(){}; //default constructor
  V_hppp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_hphh{
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hphh(){}; //default constructor
  V_hphh(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct V_pphp{
  double *V_3_3;
  int *V_3_3_index;
  int V_3_3_length;
  V_pphp(){}; //default constructor
  V_pphp(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct Interactions{
  F_matrix Fmatrix;
  V_hhhh Vhhhh;
  V_pppp Vpppp;
  V_hhpp Vhhpp;
  V_pphh Vpphh;
  V_hphp Vhphp;
  V_hhhp Vhhhp;
  V_hppp Vhppp;
  V_hphh Vhphh;
  V_pphp Vpphp;
  Interactions(){}; //default constructor
  Interactions(Input_Parameters &Parameters, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
};

struct X_hh{
  double *X1_2;
  double *X_2;
  double *X_3;
  double *X1_3;
  double *X_3d;
  double *X1_3d;
  double *X_3od;
  double *X1_3od;

  int *X_3_index;
  int *X_3d_index;
  int X_2_length;
  int X_3_length;
  int X_3d_length;

  int *map_chan;
  int *map_ind;
  double *map_fac1;
  double *map_fac2;

  X_hh(){}; //default constructor
  X_hh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void diagonalize(Input_Parameters &Parameters, Channels &Chan);
};

struct X_pp{
  double *X_2;
  double *X_3;
  double *X_3d;
  double *X_3od;

  int *X_3_index;
  int *X_3d_index;
  int X_2_length;
  int X_3_length;
  int X_3d_length;

  int *map_chan;
  int *map_ind;
  double *map_fac1;
  double *map_fac2;

  X_pp(){}; //default constructor
  X_pp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void diagonalize(Input_Parameters &Parameters, Channels &Chan);
};

struct X_hp{
  double *X_2;
  double *X_3;

  int *X_3_index;
  int X_2_length;
  int X_3_length;

  int *map_chan;
  int *map_ind;
  double *map_fac1;
  double *map_fac2;

  X_hp(){}; //default constructor
  X_hp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_hhhh{
  double *X_1;
  double *X_3_3;
  double *X_3_4;

  int *X_1_index;
  int *X_3_3_index;
  int *X_3_4_index;
  int X_1_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_hhhh(){}; //default constructor
  X_hhhh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_pppp{
  double *X1_1;
  double *X_1;
  double *X1_3_1;
  double *X1_3_2;

  int *X_1_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int X_1_length;
  int X_3_1_length;
  int X_3_2_length;

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_pppp(){}; //default constructor
  X_pppp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void zero(Input_Parameters &Parameters, Channels &Chan);
};

struct X_hphp{
  double *X_2_1;
  double *X1_2_1;
  double *X2_2_1;
  double *X3_2_1;
  double *X1_3_1;
  double *X2_3_1;
  double *X_3_2;
  double *X1_3_2;
  double *X2_3_2;
  double *X3_3_2;
  double *X_3_3;
  double *X1_3_3;
  double *X2_3_3;
  double *X3_3_3;
  double *X3_3_4;

  int *X_2_1_index;
  int *X_3_1_index;
  int *X_3_2_index;
  int *X_3_3_index;
  int *X_3_4_index;
  int X_1_length;
  int X_2_1_length;
  int X_3_1_length;
  int X_3_2_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_hphp(){}; //default constructor
  X_hphp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X2_gather(Input_Parameters &Parameters, Channels &Chan);
  void X3_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X2_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X3_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_hhhp{
  double *X_1;
  double *X_2_3;
  double *X_3_3;
  double *X1_3_3;
  double *X1_3_4;

  int *X_1_index;
  int *X_2_3_index;
  int *X_3_3_index;
  int *X_3_4_index;
  int X_1_length;
  int X_2_3_length;
  int X_3_3_length;
  int X_3_4_length;

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_hhhp(){}; //default constructor
  X_hhhp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_hppp{
  double *X_1;
  double *X_2_3;
  double *X1_3_1;
  double *X_3_2;
  double *X1_3_2;
  double *X_3_3;
  double *X1_3_3;

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

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_hppp(){}; //default constructor
  X_hppp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_hphh{
  double *X_1;
  double *X_2_1;
  double *X_2_3;
  double *X_3_1;
  double *X1_3_1;
  double *X_3_2;
  double *X1_3_2;
  double *X_3_3;
  double *X_3_4;

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

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_hphh(){}; //default constructor
  X_hphh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
};

struct X_pphp{
  double *X_1;
  double *X_2_2;
  double *X_2_3;
  double *X_3_1;
  double *X1_3_1;
  double *X_3_2;
  double *X1_3_2;
  double *X_3_3;
  double *X1_3_3;
  double *X_3_4;
  double *X1_3_4;

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

  int *map_index;
  int *map_chan;
  int *map_ind;
  int *map_num;
  double *map_fac1;
  double *map_fac2;

  X_pphp(){}; //default constructor
  X_pphp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void X_gather(Input_Parameters &Parameters, Channels &Chan);
  void X1_gather(Input_Parameters &Parameters, Channels &Chan);
  void X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
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
  Eff_Interactions(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan);
  void delete_struct(Input_Parameters &Parameters, Channels &Chan);
  void zero(Input_Parameters &Parameters, Channels &Chan, bool flag);
  void diagonalize(Input_Parameters &Parameters, Channels &Chan);
};

#endif


