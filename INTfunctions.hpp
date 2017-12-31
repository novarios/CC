#ifndef INTFUNCTIONS_H
#define INTFUNCTIONS_H

#include <string>
#include <zlib.h>
#include <cstring>

struct Channels;
struct HF_Channels;
struct HF_Matrix_Elements;
struct Single_Particle_States;
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

void Darmstadt_Setup(HF_Channels &Chan, HF_Matrix_Elements &ME);
void Darmstadt_Read(std::string filepath, double *ME0, int count, int type);
void Darmstadt_Fill(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, double *ME0, double factor);
void Darmstadt3_Setup(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME);
void Darmstadt3_Read(std::string filepath, float *ME0, long long count, int type);
double ME3J_GetME_pn(HF_Matrix_Elements &HF_ME, int p, int q, int r, int s, int t, int u, int j, int jj, int J);
double ME3J_GetME(HF_Matrix_Elements &HF_ME, int nlj1, int nlj2, int nlj3, int nnlj1, int nnlj2, int nnlj3, int j12, int jj12, int twoJ, int t12, int tt12, int twoT);
double ME3J_Overlap(HF_Matrix_Elements &ME, int nlj1, int nlj2, int nlj3, int j12, int t12, int nlja, int nljb, int nljc, int jab, int tab, int twoJ, int twoT);
float ME3J_GetME_Ordered(HF_Matrix_Elements &HF_ME, int nlj1, int nlj2, int nlj3, int nnlj1, int nnlj2, int nnlj3, int jab0, int jjab0, int twoJ0, int tab0, int ttab0, int twoT0);

double Kinetic_Energy(int p, int q);
double Hcom_1_Body(int p, int q, double hwBar);

void Coulomb_Inf_Matrix_Elements(Channels &Chan, Interactions &Ints);
double Coulomb_Inf(int &qi, int &qj, int &qk, int &ql, double &L);
double Coulomb_HO(int &qi, int &qj, int &qk, int &ql);

void Minnesota_Matrix_Elements(Channels &Chan, Interactions &Ints);
int kron_del(int &i, int &j);
int spinExchangeMtxEle(int &i, int &j, int &k, int &l);
double vint_Minnesota_Momentum(int &qi, int &qj, int &qk, int &ql, double &L);

void Get_Fock_Matrix(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &HF, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_J(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_JM(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints);

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
  void Build(Channels &Chan);
  void Delete();
};

struct V_hhhh{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hhhh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
};

struct V_pppp{
  double *V_1;
  int *V_1_index;
  int V_1_length;
  //double *V_3_3;
  //int *V_3_3_index;
  //int V_3_3_length;
  V_pppp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
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
  void Build(Channels &Chan);
  void Delete();
};

struct V_pphh{
  //double *V_1;
  //int *V_1_index;
  //int V_1_length;
  V_pphh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
};

struct V_hphp{
  double *V_2_1;
  int *V_2_1_index;
  int V_2_1_length;
  double *V_2_2;
  int *V_2_2_index;
  int V_2_2_length;
  V_hphp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
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
  void Build(Channels &Chan);
  void Delete();
};

struct V_hppp{
  double *V_2_4;
  int *V_2_4_index;
  int V_2_4_length;
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hppp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
};

struct V_hphh{
  double *V_3_2;
  int *V_3_2_index;
  int V_3_2_length;
  V_hphh(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
};

struct V_pphp{
  double *V_3_3;
  int *V_3_3_index;
  int V_3_3_length;
  V_pphp(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
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

  double Eref;
  void get_Eref(Channels &Chan);
  Interactions(){}; //default constructor
  void Build(Channels &Chan);
  void Delete();
};

#endif


