#include "CCfunctions.hpp"

int main(int argc, char * argv[])
{
  omp_set_num_threads(4);
  if(argc == 1){
    std::string inputfile = "input.dat";
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
    Parameters.P = atoi(argv[3]);
    Parameters.N = atoi(argv[4]);
    Parameters.basis = "CART";
    Parameters.obstrength = 1.0;
    Parameters.tbstrength = 1.0;
    Print_Parameters(Parameters);
    
    Model_Space Space = CART_Build_Model_Space(Parameters);
    
    Channels Chan = CART_Setup_Channels(Space);
    
    CC_Matrix_Elements CCME = Minnesota_Matrix_Elements(Parameters, Space, Chan);
    
    HF(Parameters, Space, Chan, CCME);
    
    CCD CC = Perform_CCD(Space, Parameters, CCME, Chan);
    
    double EperA = E_Ref(Parameters, Space, Chan, CCME);
    EperA += CC.CCDE;
    EperA /= (Parameters.P + Parameters.N);
    std::cout << "E/A = " << EperA << std::endl << std::endl;
    
    CC_Eff H_Eff = Build_CC_Eff(Space, Parameters, CCME, CC, Chan);
    
    /*std::ofstream results;
    results.open("run3.txt", std::ios_base::app);
    results << Parameters.Nmax << "\t" << Parameters.P << "\t" << Parameters.N << "\t";
    results << Parameters.density << "\t" << EperA << "\n";*/
  }
  
  return 0;

}
