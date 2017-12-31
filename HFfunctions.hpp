#ifndef HFFUNCTIONS_H
#define HFFUNCTIONS_H

#include <unordered_map>

struct Channels;
struct State;
struct one_body;
struct two_body;
struct Single_Particle_States;
struct HF_Channels;
struct HF_Matrix_Elements;

void Hartree_Fock_States(HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME);
void Hartree_Fock_Step(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &Bshift0, double &error);
void Randomize_HF(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, double &width);
void Convert_To_HF_Matrix_Elements(HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME);
//void Setup_HF_Space(Model_Space &Space, Single_Particle_States &States, HF_Channels &Chan);
void Random_Step_HF(HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &width, double &error2, double &Bshift);
void Randomize_HF(HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, double &width);
void Three_Body_Mono(HF_Channels &Chan, HF_Matrix_Elements &ME);
double Eref3(HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME);

struct Single_Particle_States{
  int* hole_size;
  int* vector_size;
  int* vector_index;
  int* energy_index;
  double* vectors; 
  double* energies;
  int vector_length;
  int energy_length;

  double E3; // 3-body correction

  Single_Particle_States(){};
  void Build(HF_Channels &Chan);
  void Delete();
};

//Structure for holding channel information
struct HF_Channels{
  int size3;
  State *qnums3;
  int *nob;
  one_body *ob_vec;
  int *ob_index;
  std::unordered_map<int, int> *ob_map;

  int size1;
  State *qnums1;
  int *ntb;
  two_body *tb_vec;
  int *tb_index;
  std::unordered_map<int, int> *tb_map;

  HF_Channels(){};
  void Build();
  void Delete();
  one_body ob_state(int &chan3, int &ind3);
  two_body tb_state(int &chan1, int &ind1);
};

struct HF_Matrix_Elements{
  double *V;
  int *Index;

  float *V3;
  long long ******ME3Idx;

  long long V3mon_count;
  long long *V3mon_index;
  float *V3mon;

  /*int n_rel;
  double cutoff;
  double *ra;
  double *wra;
  double *hol;*/

  HF_Matrix_Elements(){};
  void Build(HF_Channels &HF_Chan);
  void Delete();
  double get_sixJ(int j1, int j2, int J, int jj1, int jj2, int JJ);
  void Read_J(HF_Channels &Chan);
  void Read_M(HF_Channels &Chan);
  void Read_QD(HF_Channels &Chan);
  void Read_QDFile(HF_Channels &Chan);
};

#endif
