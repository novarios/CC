#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

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
  double EperA;
  omp_set_num_threads(4);

  if(argc == 1 || argc == 2){
    std::string inputfile;
    if(argc == 1){ inputfile = "input.dat"; }
    else{ inputfile = argv[1]; }
    Get_Input_Parameters(inputfile, Parameters);
    if(Parameters.basis != "finite_J"){
      if(Parameters.basis == "infinite"){ Parameters.approx == "doubles"; }
      Build_Model_Space(Parameters, Space);
      Print_Parameters(Parameters);
      Setup_Channels(Parameters, Space, Chan);
      Setup_Amps(Parameters, Space, Chan, Amps);
      Setup_Ints(Parameters, Chan, Ints);
      Read_Matrix_Elements(Parameters, Space, Chan, Ints);
      HF(Parameters, Space, Chan, Ints);
      Perform_CC(Parameters, Space, Chan, Ints, Amps);
      EperA = E_Ref(Parameters, Space, Chan, Ints);
    }
    else if(Parameters.basis == "finite_J"){
      Model_Space_J Space_J;
      Build_Model_Space_J1(Parameters, Space_J);
      Build_Model_Space_J2(Parameters, Space_J, Space);
      Print_Parameters(Parameters);
      Setup_Channels(Parameters, Space, Chan);
      Setup_Amps(Parameters, Space, Chan, Amps);
      Setup_Ints(Parameters, Chan, Ints);
      Read_Matrix_Elements_J(Parameters, Space, Space_J, Chan, Ints);
      Perform_CC(Parameters, Space, Chan, Ints, Amps);
      EperA = E_Ref(Parameters, Space, Chan, Ints);
    }
    /*CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);
      EE_EOM_HO(Space, Parameters, H_Eff, CC, Chan);*/
  }
  else if(argc == 6){
    Parameters.calc_case = argv[1];
    Parameters.density = atof(argv[2]);
    Parameters.Nmax = atoi(argv[3]);
    Parameters.Pshells = atoi(argv[4]);
    Parameters.Nshells = atoi(argv[5]);
    Parameters.basis = "infinite";
    Parameters.approx = "doubles";
    Parameters.obstrength = 1.0;
    Parameters.tbstrength = 1.0;
    double hbarc = 0.1973269788; // eV um
    double eVs_in_Hartree = 27.21138505; // eV
    hbarc *= 10000/eVs_in_Hartree; // Hartree Angstrom
    double massc2 = 0.5109989461; // MeV
    massc2 *= 1000000/eVs_in_Hartree; // Hartree    
    double fine_struct = 1.0/137.035999139;
    double r_b = hbarc/(massc2*fine_struct);
    
    if(Parameters.calc_case == "nuclear"){
      CART_Build_Model_Space(Parameters, Space);
    }
    else if(Parameters.calc_case == "electronic"){
      Parameters.density = 3.0/(4.0*M_PI*pow(Parameters.density*r_b, 3));
      EG_Build_Model_Space(Parameters, Space);
    }

    Print_Parameters(Parameters);
    Setup_Channels(Parameters, Space, Chan);
    Setup_Amps(Parameters, Space, Chan, Amps);
    Setup_Ints(Parameters, Chan, Ints);
    if(Parameters.calc_case == "nuclear"){ Minnesota_Matrix_Elements(Parameters, Space, Chan, Ints); }
    else if(Parameters.calc_case == "electronic"){ Coulomb_Matrix_Elements(Parameters, Space, Chan, Ints); }
    HF(Parameters, Space, Chan, Ints);
    Perform_CC(Parameters, Space, Chan, Ints, Amps);
    EperA = E_Ref(Parameters, Space, Chan, Ints);
  }

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

  double en = Amps.get_energy(Parameters, Chan, Ints);
  std::cout << "Eref = " << EperA << ", dCCD = " << en << std::endl;
  EperA += en;
  EperA /= (Parameters.P + Parameters.N);
  std::cout << std::setprecision(10) << "E/A = " << EperA << std::endl << std::endl;
  
  return 0;
}
