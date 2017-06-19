#ifndef HFFUNCTIONS_H
#define HFFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <bitset>
#include <iomanip>
#include <omp.h>
#include <cstdarg>
#include <cstdlib>
//#include <unordered_map>

extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

struct Input_Parameters;
struct Model_Space;
struct Channels;
struct State;
struct one_body;
struct two_body;

struct Single_Particle_States;
struct HF_Channels;
struct HF_Matrix_Elements;

void Hartree_Fock_States(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME);
void Hartree_Fock_Step(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &Bshift0, double &error);
void Randomize_HF(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, double &width);
void Convert_To_HF_Matrix_Elements(Input_Parameters &Parameters, HF_Channels &Chan, Model_Space &Space, Single_Particle_States &States, HF_Matrix_Elements &ME);
void Setup_HF_Space(Model_Space &Space, Single_Particle_States &States, HF_Channels &Chan);

void Initialize_DIIS_HF(HF_Channels &Chan, Single_Particle_States &HF, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl);
void Delete_DIIS_HF(HF_Channels &Chan, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl);
void Perform_DIIS_HF(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF0, double *&p, double *&delp, double *&tempdelp, double *&B, int &N, int &maxl, int &DIIS_count);
void Update_B1_HF(HF_Channels &Chan, Single_Particle_States &HF, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);
void Update_B2_HF(HF_Channels &Chan, Single_Particle_States &HF, int &N, double *&p, double *&delp, double *&tempdelp, double *&B);
void Random_Step_HF(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &width, double &error2, double &Bshift);
void Randomize_HF(HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, double &width);

struct Single_Particle_States{
  int* hole_size;
  int* vector_size;
  int* vector_index;
  int* energy_index;
  double* vectors; 
  double* energies;
  int vector_length;
  int energy_length;

  Single_Particle_States(){};
  Single_Particle_States(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan);
  void delete_struct(HF_Channels &Chan);
  void Separate(HF_Channels &Chan);
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
  HF_Channels(Input_Parameters &Parameters, Model_Space &Space);
  void delete_struct();
  one_body ob_state(int &chan3, int &ind3);
  two_body tb_state(int &chan1, int &ind1);
};

struct HF_Matrix_Elements{
  double *V;
  int *Index;
  HF_Matrix_Elements(HF_Channels &Chan);
  HF_Matrix_Elements(){};
  void delete_struct(HF_Channels &Chan);
};

#endif
