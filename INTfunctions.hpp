#ifndef INTFUNCTIONS_H
#define INTFUNCTIONS_H

struct Input_Parameters;
struct Model_Space;
struct Channels;
struct HF_Channels;
struct HF_Matrix_Elements;

struct Interactions;
struct Doubles_ME1;
struct Singles_ME1;

void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_M(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_QD(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME);
void Read_QD_ME_From_File(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME);

void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql);

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
int kron_del(const int &i, const int &j);
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l);
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L);

void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_JM(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);

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

#endif
