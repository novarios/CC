#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

int main(int argc, char * argv[])
{
  omp_set_num_threads(4);
  if(argc == 1){
    std::string inputfile = "input.dat";
    Input_Parameters Parameters = Get_Input_Parameters(inputfile);
    Parameters.Par = 1.0;
    Parameters.M = 0;
    Parameters.T = 0;
    Model_Space_J Space_J;
    Model_Space Space;
    if(Parameters.basis == "HO_J"){
      Space_J = Build_Model_Space_J1(Parameters);
      Space = Build_Model_Space_J2(Parameters, Space_J);
    }
    else{
      Space = Build_Model_Space(Parameters);
    }
    Print_Parameters(Parameters);

    Channels Chan;
    if(Parameters.basis == "HO" || Parameters.basis == "HO_J"){ Chan = HO_Setup_Channels(Space); }
    else{ Chan = CART_Setup_Channels(Space); }
    
    CC_Matrix_Elements CCME(Chan);
    if(Parameters.basis == "HO_J"){
      CCME = Read_Matrix_Elements_J(Parameters.MatrixElements, Parameters, Space, Space_J, Chan);
    }
    else{
      CCME = Read_Matrix_Elements(Parameters.MatrixElements, Parameters, Space, Chan);
    }
    
    if(Parameters.basis != "HO_J"){ HF(Parameters, Space, Chan, CCME); }
    
    CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
    
    double EperA = E_Ref(Parameters, Space, Chan, CCME);
    //std::cout << "Eref = " << EperA << std::endl;
    std::cout << "Eref = " << EperA << ", dCCD = " << CC.CCDE << std::endl;
    EperA += CC.CCDE;
    EperA /= (Parameters.P + Parameters.N);
    std::cout << "E/A = " << EperA << std::endl << std::endl;

    /*CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);
    
    EE_EOM_HO(Space, Parameters, H_Eff, CC, Chan);*/
  }
  else if(argc == 2){
    std::string inputfile = argv[1];
    Input_Parameters Parameters = Get_Input_Parameters(inputfile);
    Print_Parameters(Parameters);
    Model_Space Space = Build_Model_Space(Parameters);
    Channels Chan;
    if(Parameters.basis == "HO"){ Chan = HO_Setup_Channels(Space); }
    else{ Chan = CART_Setup_Channels(Space); }
    
    CC_Matrix_Elements CCME = Read_Matrix_Elements(Parameters.MatrixElements, Parameters, Space, Chan);
    
    HF(Parameters, Space, Chan, CCME);
    CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
    
    double EperA = E_Ref(Parameters, Space, Chan, CCME);
    EperA += CC.CCDE;
    EperA /= (Parameters.P + Parameters.N);
    std::cout << "E/A = " << EperA << std::endl << std::endl;

    CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);
  }
  else if(argc == 5){  
    Input_Parameters Parameters;
    Parameters.density = atof(argv[1]);
    Parameters.Nmax = atoi(argv[2]);
    Parameters.Pshells = atoi(argv[3]);
    Parameters.Nshells = atoi(argv[4]);
    Parameters.basis = "CART";
    Parameters.obstrength = 1.0;
    Parameters.tbstrength = 1.0;
    
    Model_Space Space = CART_Build_Model_Space(Parameters);
    Print_Parameters(Parameters);    

    Channels Chan = CART_Setup_Channels(Space);
    
    CC_Matrix_Elements CCME = Minnesota_Matrix_Elements(Parameters, Space, Chan);
    
    HF(Parameters, Space, Chan, CCME);
    
    CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
    
    double EperA = E_Ref(Parameters, Space, Chan, CCME);
    EperA += CC.CCDE;
    EperA /= (Parameters.P + Parameters.N);
    std::cout << "E/A = " << EperA << std::endl << std::endl;
    
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
  }
 else if(argc == 10){  
    Input_Parameters Parameters;
    Parameters.density = atof(argv[1]);
    Parameters.Nmax = atoi(argv[2]);
    Parameters.Pshells = atoi(argv[3]);
    Parameters.Nshells = atoi(argv[4]);
    //For Excited States
    Parameters.Nx = atoi(argv[5]);
    Parameters.Ny = atoi(argv[6]);
    Parameters.Nz = atoi(argv[7]);
    Parameters.M = atof(argv[8]);
    Parameters.T = atof(argv[9]);
    //////
    Parameters.basis = "CART";
    Parameters.obstrength = 1.0;
    Parameters.tbstrength = 1.0;
    
    Model_Space Space = CART_Build_Model_Space(Parameters);
    Print_Parameters(Parameters);
    Channels Chan = CART_Setup_Channels(Space);
    
    CC_Matrix_Elements CCME = Minnesota_Matrix_Elements(Parameters, Space, Chan);

    HF(Parameters, Space, Chan, CCME);

    CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
    
    double EperA = E_Ref(Parameters, Space, Chan, CCME);
    EperA += CC.CCDE;
    EperA /= (Parameters.P + Parameters.N);
    std::cout << "E/A = " << EperA << std::endl << std::endl;
    
    CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);
    
    EE_EOM(Space, Parameters, H_Eff, CC, Chan);
  }
  
  return 0;

}
