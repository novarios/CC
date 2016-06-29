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

extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

struct Single_Particle_States;
struct HF_Channels;
struct HF_Matrix_Elements;

void Separate_Particles_Holes(Single_Particle_States &States, const HF_Channels &Chan);
void Build_Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States);
void Setup_Channels_HF(const Input_Parameters &Parameters, const Model_Space &Space, HF_Channels &Chan);
void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Read_Matrix_Elements_M(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &HF_ME);
void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States, const HF_Matrix_Elements &ME);
void Hartree_Fock_States_J(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States, const HF_Matrix_Elements &ME);
void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME);
void Setup_HF_Space(Model_Space &Space, const Single_Particle_States &States, const HF_Channels &Chan);
void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints);
void Read_Matrix_Elements_HO(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME);

struct Single_Particle_States{
  int hp;
  int hn;
  std::vector<std::vector<std::vector<double> > > holes; //list of sp hole states given as vectors of coefficients
  std::vector<std::vector<std::vector<double> > > particles; //list of sp particle states given as vectors of coefficients
  std::vector<std::vector<double> > h_energies; //list of sp-hole energies
  std::vector<std::vector<double> > pt_energies; //list of sp-particle energies
  std::vector<std::vector<std::vector<double> > > vectors; 
  std::vector<std::vector<double> > energies;
  std::vector<int> h; //number of holes in OBchan
  std::vector<int> p; //number of particles in OBchan
};

//Structure for holding channel information
struct HF_Channels{
  int size1;
  int size3;
  std::vector<struct State> qnums1;
  std::vector<struct State> qnums3;
  std::vector<int> indvec;
  std::vector<int> tb;
  std::vector<std::vector<int> > tbvec1;
  std::vector<int> ob;
  std::vector<std::vector<int> > obvec1;
};

struct HF_Matrix_Elements{
  std::vector<std::vector<double> > V;
  HF_Matrix_Elements(const HF_Channels &Chan);
  HF_Matrix_Elements();
};

#endif
