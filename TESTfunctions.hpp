#ifndef TESTFUNCTIONS_H
#define TESTFUNCTIONS_H

void Doubles_Step_explicit(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2, double &error);
void Doubles_Step_explicit2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps1, Amplitudes &Amps2, const HF_Channels &HF_Chan,
			    const HF_Matrix_Elements &HF_ME, double &error);
void CC_compare_JM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, std::string &inputfile);

#endif
