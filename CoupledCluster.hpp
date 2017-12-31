#ifndef COUPLEDCLUSTER_H
#define COUPLEDCLUSTER_H

#include <string>
#include <iomanip>

const std::string PATH = "inputs/";
#define PI 3.141592653589793 // real value
//#define PI 3.141592741012573 // Morten's value
#define e 2.718281828459045
#define hbarc_MeVfm 197.3269788 // MeV fm
//#define hbarc_MeVfm 197.326968 // MeV fm ( Morten's value )
#define hbarc_eVum 0.1973269788 // eV um
#define hbarc_HartA 72.5163303659 // Hart A
#define m_neutronc2 939.5654133 // MeV
#define m_protonc2 938.2720813 // MeV
////#define m_protonc2 939.565378 // MeV
//#define m_neutronc2 938.918725 // MeV ( Morten's value )
//#define m_protonc2 938.918725 // MeV ( Mortens's value )
#define m_electronc2 0.5109989461 // MeV
#define m_electronc2_Hart 18778.865727 // Hart
#define eVs_in_Hartree 27.21138505 // eV
#define fine_struct 0.007297352566355
#define INVM 20.7355285386 // hbarc^2/(2 MN)  [MeV fm^2]

struct HF_Channels;
struct HF_Matrix_Elements;
struct Single_Particle_States;
struct Channels;
struct Amplitudes;
struct Interactions;
struct Eff_Interactions;

void Initialize_From_File(std::string &inputfile, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME);
void Initialize_Darmstadt(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME);
void Setup_From_File(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &States, Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints);
void Setup_Inf_Nuclear(Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints);
void Setup_Inf_Electronic(Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints);
void Setup_QD(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &States, Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints);

void CC_Ground_State(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps);
void Print_Parameters();
void Print_CC_Result(Interactions &Ints, Amplitudes &Amps);

//Structure for holding Input parameters
struct Input_Parameters{
  std::string calc_case; //nuclear or electronic
  std::string basis; //infinite, finite (HO), or finite_J (HO_J)
  std::string approx; //doubles, singles, or triples

  int CoM; //1 for CoM correction
  double CoM_hw;
  int ME3; //1 for 3N-ME
  double CoM_fac; //+1, -1 for E_CoM
  int HF; //1 to perform HF, anything else for bypass
  int JM; //for J->M runs

  double obstrength; //one-body multiplication strength
  double tbstrength; //two-body multiplication strength

  int Nshells; //number of neutrons shells
  int Pshells; //number of protons shells
  int Shells; //Nmax -> Shells
  int N; //number of neutrons
  int P; //number of protons
  int A; //P + N

  double density;
  double ho_energy;
  double ho_length;

  std::string LevelScheme;      // level scheme path
  std::string MatrixElements;   // matrix elements path
  std::string MatrixElements3;  // 3N-matrix elements path

  //For Excited States
  int extra; // -1 for pr, 0 for es, 1 for pa

  Input_Parameters(){};
  void Initialize(std::string &infile);
};
extern struct Input_Parameters PAR;

#endif
