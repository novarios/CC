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

struct Single_Particle_States;
struct HF_Channels;
struct HF_Matrix_Elements;

void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States, const HF_Matrix_Elements &ME);
void Hartree_Fock_Step(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, const HF_Matrix_Elements &ME, const double &Bshift, double &error);
void Randomize_HF(const HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, double &width);
void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME);
void Setup_HF_Space(Model_Space &Space, const Single_Particle_States &States, const HF_Channels &Chan);

struct Single_Particle_States{
  int hp;
  int hn;
  double*** holes; //list of sp hole states given as vectors of coefficients
  double*** particles; //list of sp particle states given as vectors of coefficients
  double** h_energies; //list of sp-hole energies
  double** pt_energies; //list of sp-particle energies
  double*** vectors; 
  double** energies;
  int* h; //number of holes in OBchan
  int* p; //number of particles in OBchan

  Single_Particle_States(){};
  Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan);
  void delete_struct(const HF_Channels &Chan);
  void Separate(const HF_Channels &Chan);
};

//Structure for holding channel information
struct HF_Channels{
  int size1;
  int size3;
  State *qnums1;
  State *qnums3;
  int *indvec;
  int *ntb;
  int **tbvec;
  std::unordered_map<int, int> *tb_map;
  int *nob;
  int **obvec;
  std::unordered_map<int, int> *ob_map;
  HF_Channels(){};
  HF_Channels(const Input_Parameters &Parameters, const Model_Space &Space);
  void delete_struct();
};

struct HF_Matrix_Elements{
  double **V;
  HF_Matrix_Elements(const HF_Channels &Chan);
  HF_Matrix_Elements(){};
  void delete_struct(const HF_Channels &Chan);
};

#endif
