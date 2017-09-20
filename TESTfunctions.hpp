#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

void CC_Test_Full(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);

void CC_Test_T1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
double Diagram_D1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_D_2c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_D_5e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_7b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_D_2d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_D_3a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_D_5f(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_7a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_X_2c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);
double Diagram_X_2d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print);

void CC_Test_T2_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
void CC_Test_T2_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
void CC_Test_T2_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
void CC_Test_T2_4(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
double Diagram_D_2e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_3b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_5c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_5d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_7c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_X_2e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);

void CC_Test_T3_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
void CC_Test_T3_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
double Diagram_D_3d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_5g(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_7e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_4b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_6b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_X_2a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_X_4b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);

void CC_Test_T3_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
void CC_Test_T3_4(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print);
double Diagram_D_3c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_5h(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_7d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_4a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_6a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_6c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_8a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_8b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_D_9(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_X_2b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);
double Diagram_X_4a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute);

void CC_Test_Full_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T2_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T2_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T2_3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T2_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T3_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T3_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T3_3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_T3_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_t2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);
void CC_Test_t3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps);

void CC_Test_Xpp2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xhp2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xhh2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xpphp3_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xhppp3_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xhphp2_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);
void CC_Test_Xpppp1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints);

#endif
