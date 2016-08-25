#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"
#include "HFfunctions.hpp"

int main(int argc, char * argv[])
{
  // basis = "infinite", "finite", "finite_J"
  // case = "nuclear", "electronic"
  // approx = "doubles", "singles", "triples"

  Input_Parameters Parameters;
  Model_Space Space;
  Channels Chan;
  Amplitudes Amps;
  Interactions Ints;
  HF_Channels HF_Chan;
  HF_Matrix_Elements HF_ME;
  Single_Particle_States States;
  double EperA;
  //omp_set_num_threads(8);
  std::cout << std::setprecision(12);

  if(argc == 1 || argc == 2){
    std::string inputfile;
    if(argc == 1){ inputfile = "input.dat"; }
    else{ inputfile = argv[1]; }

    Get_Input_Parameters(inputfile, Parameters);
    if(Parameters.basis == "infinite"){ Parameters.approx == "doubles"; }

    Build_Model_Space(Parameters, Space);
    Print_Parameters(Parameters, Space);

    std::cout << "Performing Hartree-Fock Diagonalization ..." << std::endl;
    HF_Chan = HF_Channels(Parameters, Space);
    States = Single_Particle_States(Parameters, Space, HF_Chan);
    Read_Matrix_Elements_M(Parameters, Space, HF_Chan, HF_ME);
    Hartree_Fock_States(Parameters, Space, HF_Chan, States, HF_ME);
    Convert_To_HF_Matrix_Elements(HF_Chan, States, HF_ME);
    States.delete_struct(HF_Chan);
    Chan = Channels(Parameters, Space);
    Amps = Amplitudes(Parameters, Space, Chan);
    Ints = Interactions(Parameters, Chan);
    Get_Matrix_Elements(Parameters, HF_Chan, HF_ME, Space, Chan, Ints);
    HF_ME.delete_struct(HF_Chan);
    HF_Chan.delete_struct();

    Perform_CC(Parameters, Space, Chan, Ints, Amps);
    EperA = E_Ref(Parameters, Space, Chan, Ints);

    if(Parameters.extra == 1 || Parameters.extra == 0){
      CC_Eff V_Eff(Parameters, Space, Chan);
      Build_CC_Eff(Parameters, Space, Chan, Ints, Amps, V_Eff);
      if(Parameters.extra == 0){ EE_EOM(Parameters, Space, Chan, V_Eff); }
      else if(Parameters.extra == 1){ PA_EOM(Parameters, Space, Chan, V_Eff); }
      V_Eff.delete_struct(Chan);
    }

    double en = Amps.get_energy(Parameters, Chan, Ints);
    std::cout << "Eref = " << EperA << ", dCCD = " << en << std::endl;
    EperA += en;
    std::cout << "E = " << EperA << ", E/A = " << EperA/(Parameters.P + Parameters.N) << std::endl << std::endl;
  }
  else if(argc == 6 || argc == 7){
    Parameters.calc_case = argv[1];
    Parameters.density = atof(argv[2]);
    Parameters.Shells = atoi(argv[3]);
    Parameters.Pshells = atoi(argv[4]);
    Parameters.Nshells = atoi(argv[5]);
    Parameters.obstrength = 1.0;
    Parameters.tbstrength = 1.0;
    
    if(Parameters.calc_case == "nuclear"){
      Parameters.basis = "infinite";
      Parameters.approx = "doubles";
      CART_Build_Model_Space(Parameters, Space);
      Print_Parameters(Parameters, Space);
      Chan = Channels(Parameters, Space);
      Amps = Amplitudes(Parameters, Space, Chan);
      Ints = Interactions(Parameters, Chan);
      Minnesota_Matrix_Elements(Parameters, Space, Chan, Ints);
      Perform_CC(Parameters, Space, Chan, Ints, Amps);
      EperA = E_Ref(Parameters, Space, Chan, Ints);
    }
    else if(Parameters.calc_case == "electronic"){
      Parameters.basis = "infinite";
      Parameters.approx = "doubles";
      CART_Build_Model_Space(Parameters, Space);
      Print_Parameters(Parameters, Space);
      Chan = Channels(Parameters, Space);
      Amps = Amplitudes(Parameters, Space, Chan);
      Ints = Interactions(Parameters, Chan);
      Coulomb_Inf_Matrix_Elements(Parameters, Space, Chan, Ints);
      Perform_CC(Parameters, Space, Chan, Ints, Amps);
      EperA = E_Ref(Parameters, Space, Chan, Ints);
    }
    else if(Parameters.calc_case == "quantum_dot"){
      Parameters.basis = "finite_HO";
      Parameters.approx = "singles";
      QD_Build_Model_Space(Parameters, Space);
      Print_Parameters(Parameters, Space);

      std::cout << "Performing Hartree-Fock Diagonalization ..." << std::endl;
      HF_Chan = HF_Channels(Parameters, Space);
      States = Single_Particle_States(Parameters, Space, HF_Chan);
      Read_Matrix_Elements_HO(Parameters, Space, HF_Chan, HF_ME);
      Hartree_Fock_States(Parameters, Space, HF_Chan, States, HF_ME);
      Convert_To_HF_Matrix_Elements(HF_Chan, States, HF_ME);
      States.delete_struct(HF_Chan);
      Chan = Channels(Parameters, Space);
      Amps = Amplitudes(Parameters, Space, Chan);
      Ints = Interactions(Parameters, Chan);
      Get_Matrix_Elements(Parameters, HF_Chan, HF_ME, Space, Chan, Ints);
      HF_ME.delete_struct(HF_Chan);
      HF_Chan.delete_struct();

      Perform_CC(Parameters, Space, Chan, Ints, Amps);
      EperA = E_Ref(Parameters, Space, Chan, Ints);
    }
    if(argc == 7){
      Parameters.extra = atoi(argv[6]);
      CC_Eff V_Eff(Parameters, Space, Chan);
      Build_CC_Eff(Parameters, Space, Chan, Ints, Amps, V_Eff);
      if(Parameters.extra == 0){ EE_EOM(Parameters, Space, Chan, V_Eff); }
      else if(Parameters.extra == 1){ PA_EOM(Parameters, Space, Chan, V_Eff); }
      else if(Parameters.extra == -1){ PR_EOM(Parameters, Space, Chan, V_Eff); }
      V_Eff.delete_struct(Chan);
    }
    double en = Amps.get_energy(Parameters, Chan, Ints);
    std::cout << "Eref = " << EperA << ", dCCD = " << en << std::endl;
    EperA += en;
    std::cout << "E = " << EperA << ", E/A = " << EperA/(Parameters.P + Parameters.N) << std::endl << std::endl;
  }
  Ints.delete_struct(Parameters, Chan);
  Amps.delete_struct(Parameters, Chan);
  Chan.delete_struct();
  Space.delete_struct(Parameters);

  /*Model_Space Space = CART_Build_Model_Space(Parameters);
  Print_Parameters(Parameters);    
  
  Channels Chan = CART_Setup_Channels(Space);
  
  CC_Matrix_Elements CCME = Minnesota_Matrix_Elements(Parameters, Space, Chan);
  
  HF(Parameters, Space, Chan, CCME);
  
  CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
  
  double EperA = E_Ref(Parameters, Space, Chan, CCME);
  EperA += CC.CCDE;
  EperA /= (Parameters.P + Parameters.N);
  std::cout << "E/A = " << EperA << std::endl << std::endl;*/
  
  //CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);*/
  
  //double hbarc = 197.3269788; // MeVfm
  //double m_neutronc2 = 939.565378; // MeV
  //double m_protonc2 = 938.272046; // MeV
  //double m_protonc2 = 939.565378; // MeV
  //double neutron_prefac = hbarc*hbarc/(2.0*m_neutronc2);
  //double proton_prefac = hbarc*hbarc/(2.0*m_protonc2);
  //double kf = pow(3.0 * M_PI * M_PI * Parameters.density, 1.0/3.0);
  //double L = pow((Parameters.P + Parameters.N) / Parameters.density, 1.0/3.0);
  
  //Print_Parameters(Parameters);
  
  /*std::vector<double> wi(5);
    wi[0] = 0.2369268851;
    wi[1] = 0.4786286705;
    wi[2] = 0.5688888889;
    wi[3] = 0.4786286705;
    wi[4] = 0.2369268851;
    std::vector<double> xi(5);
    xi[0] = -0.9061798459;
    xi[1] = -0.5384693101;
    xi[2] = 0.0;
    xi[3] = 0.5384693101;
    xi[4] = 0.9061798459;

    double energy, Etot = 0.0;
    int N = Parameters.P + Parameters.N;
    for(int tx = 0; tx < 5; ++tx){
      for(int ty = 0; ty < 5; ++ty){
	for(int tz = 0; tz < 5; ++tz){
	  std::cout << 25*tx+5*ty+tz << " ";
	  Model_Space Space = CART_Build_Model_Space_Twist(Parameters, 0.5*M_PI*(xi[tx]+1), 0.5*M_PI*(xi[ty]+1), 0.5*M_PI*(xi[tz]+1));
	  double L = pow((Parameters.P + Parameters.N) / Parameters.density, 1.0/3.0);
	  energy = 0.0;
	  for(int i = 0; i < Space.indtot; ++i){
	    if(Space.levelstype[i] == "hole"){
	      for(int j = 0; j < Space.indtot; ++j){
		if(Space.levelstype[j] == "hole" && i != j){
		  energy += 0.5 * vint_Minnesota_Momentum(Space, i, j, i, j, L);
		}
	      }
	    }
	  }
	  Etot += 0.125*wi[tx]*wi[ty]*wi[tz]*energy;
	}
      }
    }
    Etot = Etot/N;*/

    /*double Etot = 0.0;
    Model_Space Space = CART_Build_Model_Space(Parameters);
    int N = Parameters.P + Parameters.N;
    double L = pow((Parameters.P + Parameters.N) / Parameters.density, 1.0/3.0);
    for(int i = 0; i < Space.indtot; ++i){
      if(Space.levelstype[i] == "hole"){
	for(int j = 0; j < Space.indtot; ++j){
	  if(Space.levelstype[j] == "hole" && i != j){
	    Etot += 0.5 * vint_Minnesota_Momentum(Space, i, j, i, j, L);
	  }
	}
      }
    }
    Etot = Etot/N;*/

    //double energyi = 2.0*neutron_prefac*(3.0/10.0)*kf*kf;
    //double energyi = -21.4791;
    //std::cout << Etot << " " << energyi << std::endl;

    /*std::ofstream results;
    results.open("output_numholes2.txt", std::ios_base::app);
    results << Parameters.Nmax << "\t" << Parameters.P << "\t" << Parameters.N << "\t";
    results << Parameters.density << "\t" << EperA << "\n";*/

    /*std::ofstream results;
    results.open("output_particlesHF.txt", std::ios_base::app);
    results << Parameters.Nmax << "\t" << Parameters.P << "\t" << Parameters.N << "\t";
    results << Parameters.density << "\t" << 1 - Etot/energyi << "\n";*/
  
  return 0;
}
