#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "HFfunctions.hpp"
#include "CCfunctions.hpp"
#include "LAMBDAfunctions.hpp"
#include "INTfunctions.hpp"
#include "EffINTfunctions.hpp"
#include "EOMfunctions.hpp"
#include "AngMom.hpp"

struct Input_Parameters PAR;
struct Model_Space SPB;
struct AngMom JCOUP;

int main(int argc, char * argv[])
{
  // basis = "infinite", "finite", "finite_J"
  // case = "nuclear", "electronic"
  // approx = "doubles", "singles", "triples"
  std::srand(time(NULL));
  struct timespec time1, time2;
  double elapsed0 = 0.0;

  Channels Chan;
  Amplitudes Amps;
  Interactions Ints;
  Eff_Interactions Eff_Ints;
  HF_Channels HF_Chan;
  HF_Matrix_Elements HF_ME;
  Single_Particle_States States;
  std::string calctype;
  std::string inputfile;

  PAR.HF = 1;
  PAR.CoM = 1; // 0 = no COM correction, 1 = regular COM correction (intrinsic H), 2 = calc E_com
  PAR.CoM_fac = 0.0;
  PAR.obstrength = 1.0;
  PAR.tbstrength = 1.0;
  PAR.ME3 = 0; // No 3N matrix elements default
  PAR.JM = 0;

  clock_gettime(CLOCK_MONOTONIC, &time1);
  if(argc == 1 || argc == 2){
    if(argc == 1){ inputfile = "input.dat"; }
    else{ inputfile = argv[1]; }
    Initialize_From_File(inputfile, HF_Chan, HF_ME);
    Setup_From_File(HF_Chan, HF_ME, States, Chan, Amps, Ints, Eff_Ints);
    CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
  }
  else{
    PAR.calc_case = argv[1]; // infinite_e, infinite_n, quantum_dot, finite_J, finite_M, darmstadt
    if( PAR.calc_case == "inifinite_e" ){
      PAR.density = atof(argv[2]);
      PAR.Shells = atoi(argv[3]);
      PAR.Pshells = atoi(argv[4]);
      PAR.Nshells = 0;
      Setup_Inf_Electronic(Chan, Amps, Ints, Eff_Ints);
      CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
    }
    else if( PAR.calc_case == "inifinite_n" ){
      PAR.density = atof(argv[2]);
      PAR.Shells = atoi(argv[3]);
      PAR.Pshells = atoi(argv[4]);
      PAR.Nshells = atoi(argv[5]);
      Setup_Inf_Nuclear(Chan, Amps, Ints, Eff_Ints);
      CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
    }
    else if( PAR.calc_case == "quantum_dot" ){
      PAR.density = atof(argv[2]);
      PAR.Shells = atoi(argv[3]);
      PAR.Pshells = atoi(argv[4]);
      PAR.Nshells = 0;
      Setup_QD(HF_Chan, HF_ME, States, Chan, Amps, Ints, Eff_Ints);
      CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
    }
    else if( PAR.calc_case == "darmstadt" ){
      PAR.MatrixElements = argv[2];
      PAR.Pshells = atoi(argv[3]);
      PAR.Nshells = atoi(argv[4]);
      if( argc == 6 ){
	PAR.MatrixElements3 = argv[5];
	PAR.ME3 = 1;
      }
      Initialize_Darmstadt(HF_Chan, HF_ME);
      Setup_From_File(HF_Chan, HF_ME, States, Chan, Amps, Ints, Eff_Ints);

      if(PAR.CoM < 2){
	CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
      }
      else{
	//PAR.CoM_fac = 0.001;
	PAR.CoM_fac = 1.0;
	Initialize_Darmstadt(HF_Chan, HF_ME);
	Setup_From_File(HF_Chan, HF_ME, States, Chan, Amps, Ints, Eff_Ints);
	CC_Ground_State(Chan, Ints, Eff_Ints, Amps);
	/*PAR.CoM_fac = -0.001;
	Initialize_Darmstadt(HF_Chan, HF_ME);
	Setup_From_File(HF_Chan, HF_ME, States, Chan, Amps, Ints, Eff_Ints, Energy0);
	CC_Ground_State(Chan, Ints, Eff_Ints, Amps, Energy0, dEnergy, Energy2);
	CoM_Energy = (Energy1 - Energy2)/0.002;
	std::cout << "!!!   CoM_Energy = " << CoM_Energy - 1.5*PAR.ho_energy << std::endl;*/
      }
    }
  }

  Eff_Ints.Update_3(Chan, Ints, Amps);

  //LAmplitudes LAmps;
  //LAmps.Build(Chan);
  //LCC_Algorithm(Chan, Ints, Eff_Ints, Amps, LAmps);

  double Energy = 100.0;
  EOM EOM_1PA, EOM_1PA1;
  EOM_1PA.PA1_EOM(Chan, Ints, Amps, Eff_Ints);
  EOM_1PA.Print_EOM_1P(Energy);
  EOM_1PA1.EOM_PA(Chan, Ints, Amps, Eff_Ints);
  EOM_1PA1.Print_EOM_1P(Energy);

  EOM_1PA.Delete();
  EOM_1PA1.Delete();

  EOM EOM_1PR, EOM_1PR1;
  EOM_1PR.PR1_EOM(Chan, Ints, Amps, Eff_Ints);
  EOM_1PR.Print_EOM_1P(Energy);
  EOM_1PR1.EOM_PR(Chan, Ints, Amps, Eff_Ints);
  EOM_1PR1.Print_EOM_1P(Energy);

  EOM_1PR.Delete();
  EOM_1PR1.Delete();
    
  Ints.Delete();
  Amps.Delete();
  Eff_Ints.Delete();
  if( PAR.basis == "finite_J" || PAR.basis == "finite_JM" ){ JCOUP.Delete(); }
  Chan.Delete();
  SPB.Delete();

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "Total Runtime = " << elapsed0 << " sec. " << std::endl;

  return 0;
}

void Initialize_From_File(std::string &inputfile, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME)
{
  PAR.Initialize(inputfile);
  if(PAR.basis == "infinite"){ PAR.approx == "doubles"; }
  SPB.Build();
  Print_Parameters();
  HF_Chan.Build();
  HF_ME.Build(HF_Chan);
  if(PAR.basis == "finite_J" || PAR.basis == "finite_JM"){ HF_ME.Read_J(HF_Chan); }
  else{ HF_ME.Read_M(HF_Chan); }
}

void Initialize_Darmstadt(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME)
{
  PAR.calc_case = "nuclear";
  PAR.basis = "finite_J";
  PAR.approx = "singles";
  Darmstadt_Setup(HF_Chan, HF_ME);
  if( PAR.ME3 == 1 ){ Darmstadt3_Setup(HF_Chan, HF_ME); }
}

void Setup_From_File(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &States, Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints)
{
  if( PAR.basis == "finite_J" || PAR.basis == "finite_JM" ){ JCOUP.Build(); }
  States.Build(HF_Chan);
  if(PAR.HF == 1){
    Hartree_Fock_States(HF_Chan, States, HF_ME);
    Convert_To_HF_Matrix_Elements(HF_Chan, States, HF_ME);
  }
  if(PAR.basis == "finite_JM"){
    PAR.basis = "finite_M";
    PAR.JM = 1;
    SPB.Build_J2();
  }
  Chan.Build();
  Amps.Build(Chan);
  Ints.Build(Chan);
  if( PAR.ME3 == 1 ){ Ints.Eref += Eref3(HF_Chan, States, HF_ME); }
  
  if( PAR.HF == 0 ){ Get_Fock_Matrix(HF_Chan, HF_ME, States, Chan, Ints); }
  if( PAR.basis == "finite_J" ){ Get_Matrix_Elements_J(HF_Chan, HF_ME, Chan, Ints); }
  else if( PAR.JM == 1 ){ Get_Matrix_Elements_JM(HF_Chan, HF_ME, Chan, Ints); }
  else{ Get_Matrix_Elements(HF_Chan, HF_ME, Chan, Ints); }
  States.Delete();
  HF_ME.Delete();
  HF_Chan.Delete();
  Eff_Ints.Build(Chan);
}

void Setup_Inf_Nuclear(Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints)
{
  PAR.basis = "infinite";
  PAR.approx = "doubles";
  SPB.Build_CART();
  Print_Parameters();
  Chan.Build();
  Amps.Build(Chan);
  Ints.Build(Chan);
  Minnesota_Matrix_Elements(Chan, Ints);
  Eff_Ints.Build(Chan);
}

void Setup_Inf_Electronic(Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints)
{
  PAR.basis = "infinite";
  PAR.approx = "doubles";
  SPB.Build_CART();
  Print_Parameters();
  Chan.Build();
  Amps.Build(Chan);
  Ints.Build(Chan);
  Coulomb_Inf_Matrix_Elements(Chan, Ints);
  Eff_Ints.Build(Chan);
}

void Setup_QD(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &States, Channels &Chan, Amplitudes &Amps, Interactions &Ints, Eff_Interactions &Eff_Ints)
{
  PAR.HF = 1;
  PAR.basis = "finite_HO";
  PAR.approx = "singles";
  SPB.Build_QD();
  Print_Parameters();
  HF_Chan.Build();
  HF_ME.Build(HF_Chan);
  States.Build(HF_Chan);
  //HF_ME.Read_QD(HF_Chan);
  HF_ME.Read_QDFile(HF_Chan);  
  if(PAR.HF == 1){
    Hartree_Fock_States(HF_Chan, States, HF_ME);
    Convert_To_HF_Matrix_Elements(HF_Chan, States, HF_ME);
  }
  Chan.Build();
  Amps.Build(Chan);
  Ints.Build(Chan);
  if(PAR.HF == 0){ Get_Fock_Matrix(HF_Chan, HF_ME, States, Chan, Ints); }
  Get_Matrix_Elements(HF_Chan, HF_ME, Chan, Ints);
  States.Delete();
  HF_ME.Delete();
  HF_Chan.Delete();
  Eff_Ints.Build(Chan);
}
  
void CC_Ground_State(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  CC_Algorithm(Chan, Ints, Eff_Ints, Amps);
  Ints.get_Eref(Chan);
  Amps.get_dE(Chan, Ints);
  Print_CC_Result(Ints, Amps);
}

void Input_Parameters::Initialize(std::string &infile)
{ 
  std::string path; // PAR File Path
  std::string line; // Line from input file
  std::ifstream filestream; // File stream
  int index; // Count input line
  size_t colon; // Size for finding ':' which proceeds every input
  std::string substr; // Substd::string for extracting input

  path = PATH + infile;
  filestream.open(path.c_str());
  if (!filestream.is_open()){ std::cerr << "PAR file, " << path << ", does not exist" << std::endl; exit(1); }

  //find lines that start with '\*'
  index = 0;
  while (getline(filestream, line)){ 
    if(line[0] == '\\' && line[1] == '*'){  
      ++index;
      colon = line.find(':');
      if( colon == line.size() - 1 ){ continue; };
      substr = line.substr(colon + 2, line.size());
      switch(index){
      case 1:
	PAR.calc_case = substr;
	break;
      case 2:
	PAR.basis = substr;
	break;
      case 3:
	PAR.approx = substr;
	break;
      case 4:
	PAR.obstrength = atof(substr.c_str());
	break;
      case 5:
	PAR.tbstrength = atof(substr.c_str());
	break;
      case 6:
	PAR.Pshells = atoi(substr.c_str());
	break;
      case 7:
	PAR.Nshells = atoi(substr.c_str());
	break;
      case 8:
	PAR.LevelScheme = substr;
	break;
      case 9:
	PAR.MatrixElements = substr;
	break;
      case 10:
	PAR.extra = atoi(substr.c_str());
	break;
      } 
    }
    else{ continue; };
  }
}

void Print_Parameters()
{
  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << "Case = " << PAR.calc_case << ", Basis = " << PAR.basis << ", Approximation = " << PAR.approx << std::endl;
  if(PAR.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << PAR.LevelScheme << std::endl;
    if(PAR.MatrixElements.size() > 0){ std::cout << "Interaction = " << PAR.MatrixElements << std::endl; }
  }
  if(PAR.calc_case == "nuclear"){
    std::cout << "Number of Shells = " << PAR.Shells << ", Total States = " << SPB.num_states << std::endl;
    std::cout << "Proton Shells = " << PAR.Pshells << ", Neutron Shells = " << PAR.Nshells << std::endl;
    std::cout << "Protons = " << PAR.P << ", Neutrons = " << PAR.N << std::endl;
    if(PAR.calc_case == "infinite"){ std::cout << "Density = " << PAR.density << std::endl; }
  }
  else if(PAR.calc_case == "electronic"){
    std::cout << "Number of Shells = " << PAR.Shells << ", Total States = " << SPB.num_states << std::endl;
    std::cout << "Electron Shells = " << PAR.Pshells << ", Electrons = " << PAR.P << std::endl;
    std::cout << "Density = " << PAR.density << std::endl;
  }
  else if(PAR.calc_case == "quantum_dot"){
    std::cout << "Number of Shells = " << PAR.Shells << ", Total States = " << SPB.num_states << std::endl;
    std::cout << "Electron Shells = " << PAR.Pshells << ", Electrons = " << PAR.P << std::endl;
    std::cout << "Oscillator Energy = " << PAR.density << std::endl;
  }
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void Print_CC_Result(Interactions &Ints, Amplitudes &Amps)
{
  std::cout << std::endl;
  std::cout << std::setprecision(10);
  std::cout << "Eref = " << Ints.Eref << ", dCC = " << Amps.dE << std::endl;
  std::cout << "Eref/A = " << Ints.Eref/(PAR.P + PAR.N) << ", dCC/A = " << Amps.dE/(PAR.P + PAR.N) << std::endl;
  std::cout << "E = " << (Ints.Eref + Amps.dE) << ", E/A = " << (Ints.Eref + Amps.dE)/(PAR.P + PAR.N) << std::endl << std::endl;
}
