#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

void CC_Test_Full(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  CC_Test_T1(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T2_1(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T2_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T2_3(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T2_4(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T3_1(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T3_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T3_3(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
  CC_Test_T3_4(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, Amps2, print);
}

void CC_Test_T1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind1, chan1, chanind_T;
  int nhh, npp;
  int a, b, i, j;
  two_body pp, hh;
  double tempt, D1, D2, D3;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    chanind_T = Amps1.D1.T1_index[chan1];
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      pp = Chan.pp_state(chan1, pp_ind);
      a = pp.v1;
      b = pp.v2;
      if(a == b){ continue; }
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	i = hh.v1;
	j = hh.v2;
	if(i == j){ continue; }
	ind1 = pp_ind * nhh + hh_ind;

	D1 = 0.0;
	D2 = 0.0;
	D3 = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T1 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T1(ab|ij)  <-  V(ab|ij)  (D1)
	D1 += Diagram_D1(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);

	// T1(ab|ij)  <-  (1/2).X'(ab|cd).T(cd|ij)
	// --------------------------------------------------
	// T1(ab|ij)  <-  (1/2).V(ab|cd).T(cd|ij)  (D_2c)
	D2 += Diagram_D_2c(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);
	// T1(ab|ij)  <-  -(1/2).P(ab).V(kb|cd).t(a|k).T(cd|ij)  (D_5e)
	D2 += Diagram_D_5e(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D2 += Diagram_D_5e(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T1(ab|ij)  <-  (1/2).V(kl|cd).T(cd|ij).t(a|k).t(b|l)  (D_7b)
	D2 += Diagram_D_7b(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);

	// T1(ab|ij)  <-  (1/2).X(kl|ij).T(ab|kl)
	// --------------------------------------------------
	// T1(ab|ij)  <-  (1/2).V(kl|ij).T(ab|kl)  (D_2d)
	D3 += Diagram_D_2d(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);
	// T1(ab|ij)  <-  (1/4).V(kl|cd).T(cd|ij).T(ab|kl)  (D_3a)
	D3 += Diagram_D_3a(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);
	// T1(ab|ij)  <-  -(1/2).P(ij).V(kl|jc).t(c|i).T(ab|kl)  (D_5f)
	D3 += Diagram_D_5f(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D3 += Diagram_D_5f(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T1(ab|ij)  <-  (1/2).V(kl|cd).T(ab|kl).t(c|i).t(d|j)  (D_7a)
	D3 += Diagram_D_7a(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print);

	tempt = D1 + D2 + D3;
	if(std::fabs(Amps2.D1.T1[chanind_T + ind1] - tempt) > 1e-12){
	  double X2, X3;
	  std::cout << "!! T1 < " << a << "," << b << " |t| " << i << "," << j << " > = X: " << Amps2.D1.T1[chanind_T + ind1] << ", ";
	  std::cout << "D: " << tempt << std::endl;
	  X2 = Diagram_X_2c(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print);
	  X3 = Diagram_X_2d(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print);
	  std::cout << " X_2c: " << X2 << ", D: " << D2 << std::endl;
	  std::cout << " X_2d: " << X3 << ", D: " << D3 << std::endl;
	}
      }
    }
  }
}

// T1(ab|ij)  <-  V(ab|ij)  (D1)
double Diagram_D1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int ind_abij = key_ab * Chan.nhh[chan_ij] + key_ij;
  int chanind_V = Ints.Vpphh.V_1_index[chan_ab];
  double term = Ints.Vpphh.V_1[chanind_V + ind_abij];
  if(print == 2){
    std::cout << "D_1 = < '" << a << "','" << b << "' |V| '" << i << "','" << j << "' > = " << term << std::endl;
  }
  if(print != 0){ std::cout << "D_1 = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  (1/2).V(ab|cd).T(cd|ij)  (D_2c)
double Diagram_D_2c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij, cd;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_cd;
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_cd;
  int ind_abcd, ind_cdij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ab){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_abcd = key_ab * Chan.npp[chan_cd] + key_cd;
      ind_cdij = key_cd * Chan.nhh[chan_ij] + key_ij;
      chanind_V = Ints.Vpppp.V_1_index[chan_ab];
      chanind_T = Amps.D1.T1_index[chan_cd];
      V = Ints.Vpppp.V_1[chanind_V + ind_abcd];
      T = Amps.D1.T1[chanind_T + ind_cdij];
      term0 = 0.5 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_2c += (1/2) * < '"<< a <<"','"<< b <<"' |V| "<< c <<","<< d <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > ";
	std::cout <<"= (1/2) * "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
      term += term0;
    }
  }
  if(print != 0){ std::cout << "D_2c = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  (1/2).X'(ab|cd).T(cd|ij)  (X_2c)
double Diagram_X_2c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij, cd;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_cd;
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_cd;
  int ind_abcd, ind_cdij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ab){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_abcd = key_ab * Chan.npp[chan_cd] + key_cd;
      ind_cdij = key_cd * Chan.nhh[chan_ij] + key_ij;
      chanind_V = Eff_Ints.Xpppp.X_1_index[chan_ab];
      chanind_T = Amps.D1.T1_index[chan_cd];
      V = Eff_Ints.Xpppp.X1_1[chanind_V + ind_abcd];
      T = Amps.D1.T1[chanind_T + ind_cdij];
      term0 = 0.5 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2c += (1/2) * < '"<< a <<"','"<< b <<"' |X'| "<< c <<","<< d <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > ";
	std::cout <<"= (1/2) * "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
      term += term0;
    }
  }
  if(print != 0){ std::cout << "X_2c = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  -(1/2).P(ab).V(kb|cd).t(a|k).T(cd|ij)  (D_5e)
double Diagram_D_5e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ij, cd, cdk;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_cd, chan_cdk, chan_k;
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_b = Chan.p_map[chan_b][b];
  int key_cd, key_cdk, key_ak;
  int ind_kbcd, ind_cdij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ij){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_cdij = key_cd * Chan.nhh[chan_ij] + key_ij;
      for(int k = 0; k < Space.num_hol; ++k){
	minus(cdk, cd, Space.qnums[k]);
	chan_cdk = Space.ind_1b(Parameters.basis, cdk);
	if(chan_cdk != chan_b){ continue; }
	key_cdk = Chan.pph_map[chan_cdk][Space.hash3(c, d, k, cd.j)];
	ind_kbcd = key_b * Chan.npph[chan_cdk] + key_cdk;

	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_a){ continue; }
	key_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];

	chanind_V = Ints.Vhppp.V_3_2_index[chan_b];
	chanind_T = Amps.D1.T1_index[chan_cd];
	V = Ints.Vhppp.V_3_2[chanind_V + ind_kbcd];
	T1 = Amps.S1.t2[key_ak];
	T2 = Amps.D1.T1[chanind_T + ind_cdij];
	if(permute == 1){
	  term0 = -0.5 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5e(1) += -(1/2) * < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > ";
	    std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = 0.5 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5e(2) += (1/2) * < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5e(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5e(2) = " << term << std::endl; }
  }
  return term;
}

// T1(ab|ij)  <-  (1/2).V(kl|cd).T(cd|ij).t(a|k).t(b|l)  (D_7b)
double Diagram_D_7b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ij, cd, kl;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_cd, chan_kl, chan_k, chan_l;
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_cd, key_kl, key_ak, key_bl;
  int ind_klcd, ind_cdij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2, T3;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ij){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_cdij = key_cd * Chan.nhh[chan_ij] + key_ij;
      for(int k = 0; k < Space.num_hol; ++k){
	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_a){ continue; }
	key_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];
	for(int l = 0; l < Space.num_hol; ++l){
	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_cd){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
	  if(chan_l != chan_b){ continue; }
	  key_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T = Amps.D1.T1_index[chan_cd];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T + ind_cdij];
	  T2 = Amps.S1.t2[key_ak];
	  T3 = Amps.S1.t2[key_bl];

	  term0 = 0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_7b += (1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > ";
	    std::cout <<"* < '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	  term += term0;
	}
      }
    }
  }
  if(print != 0){ std::cout << "D_7b = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  (1/2).V(kl|ij).T(ab|kl)  (D_2d)
double Diagram_D_2d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij, kl;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_kl;
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl;
  int ind_abkl, ind_klij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int k = 0; k < Space.num_hol; ++k){
    for(int l = 0; l < Space.num_hol; ++l){
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      if(chan_kl != chan_ij){ continue; }
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      ind_abkl = key_ab * Chan.nhh[chan_kl] + key_kl;
      ind_klij = key_kl * Chan.nhh[chan_ij] + key_ij;

      chanind_V = Ints.Vhhhh.V_1_index[chan_kl];
      chanind_T = Amps.D1.T1_index[chan_ab];
      V = Ints.Vhhhh.V_1[chanind_V + ind_klij];
      T = Amps.D1.T1[chanind_T + ind_abkl];
      term0 = 0.5 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_2d += (1/2) * < "<< k <<","<< l <<" |V| '"<< i <<"','"<< j <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	std::cout <<"= (1/2) * "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
      term += term0;
    }
  }
  if(print != 0){ std::cout << "D_2d = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  (1/2).X(kl|ij).T(ab|kl)  (D_2d)
double Diagram_X_2d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij, kl;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_kl;
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl;
  int ind_abkl, ind_klij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int k = 0; k < Space.num_hol; ++k){
    for(int l = 0; l < Space.num_hol; ++l){
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      if(chan_kl != chan_ij){ continue; }
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      ind_abkl = key_ab * Chan.nhh[chan_kl] + key_kl;
      ind_klij = key_kl * Chan.nhh[chan_ij] + key_ij;

      chanind_V = Eff_Ints.Xhhhh.X_1_index[chan_kl];
      chanind_T = Amps.D1.T1_index[chan_ab];
      V = Eff_Ints.Xhhhh.X_1[chanind_V + ind_klij];
      T = Amps.D1.T1[chanind_T + ind_abkl];
      term0 = 0.5 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2d += (1/2) * < "<< k <<","<< l <<" |X| '"<< i <<"','"<< j <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	std::cout <<"= (1/2) * "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
      term += term0;
    }
  }
  if(print != 0){ std::cout << "X_2d = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  (1/4).V(kl|cd).T(cd|ij).T(ab|kl)  (D_3a)
double Diagram_D_3a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, ij, kl, cd;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_kl, chan_cd;
  if(chan_ab != chan_ij){ return 0.0; }
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl, key_cd;
  int ind_klcd, ind_cdij, ind_abkl;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T1, chanind_T2;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ij){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_cdij = key_cd * Chan.nhh[chan_ij] + key_ij;
      for(int k = 0; k < Space.num_hol; ++k){
	for(int l = 0; l < Space.num_hol; ++l){
	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_ab){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_abkl = key_ab * Chan.nhh[chan_kl] + key_kl;
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T1 = Amps.D1.T1_index[chan_cd];
	  chanind_T2 = Amps.D1.T1_index[chan_ab];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T1 + ind_cdij];
	  T2 = Amps.D1.T1[chanind_T2 + ind_abkl];
	  term0 = 0.25 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_3a += (1/4) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * < "<< c <<","<< d <<" |t| '"<< i <<"','"<< j <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	    std::cout <<"= (1/4) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	  term += term0;
	}
      }
    }
  }
  if(print != 0){ std::cout << "D_3a = " << term << std::endl; }
  return term;
}

// T1(ab|ij)  <-  -(1/2).P(ij).V(kl|jc).t(c|i).T(ab|kl)  (D_5f)
double Diagram_D_5f(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ab, kl, klc;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_kl, chan_klc, chan_c;
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_j = Chan.h_map[chan_j][j];
  int key_kl, key_klc, key_ci;
  int ind_kljc, ind_abkl;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int k = 0; k < Space.num_hol; ++k){
    for(int l = 0; l < Space.num_hol; ++l){
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      if(chan_kl != chan_ab){ continue; }
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      ind_abkl = key_ab * Chan.nhh[chan_kl] + key_kl;
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	minus(klc, kl, Space.qnums[c]);
	chan_klc = Space.ind_1b(Parameters.basis, klc);
	if(chan_klc != chan_j){ continue; }
	key_klc = Chan.hhp_map[chan_klc][Space.hash3(k, l, c, kl.j)];
	ind_kljc = key_klc * Chan.nh[chan_j] + key_j;

	chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
	if(chan_c != chan_i){ continue; }
	key_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];

	chanind_V = Ints.Vhhhp.V_3_3_index[chan_j];
	chanind_T = Amps.D1.T1_index[chan_ab];
	V = Ints.Vhhhp.V_3_3[chanind_V + ind_kljc];
	T1 = Amps.S1.t2[key_ci];
	T2 = Amps.D1.T1[chanind_T + ind_abkl];
	if(permute == 1){
	  term0 = -0.5 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5f(1) += -(1/2) * < "<< k <<","<< l <<" |V| '"<< j <<"',"<< c <<" > * ";
	    std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	    std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = 0.5 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5f(2) += (1/2) * < "<< k <<","<< l <<" |V| '"<< j <<"',"<< c <<" > * ";
	    std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5f(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5f(2) = " << term << std::endl; }
  }
  return term;
}

// T1(ab|ij)  <-  (1/2).V(kl|cd).T(ab|kl).t(c|i).t(d|j)  (D_7a)
double Diagram_D_7a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print)
{
  State ab, kl, cd;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_kl, chan_cd, chan_c, chan_d;
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_kl, key_cd, key_ci, key_dj;
  int ind_klcd, ind_abkl;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2, T3;
  for(int k = 0; k < Space.num_hol; ++k){
    for(int l = 0; l < Space.num_hol; ++l){
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      if(chan_kl != chan_ab){ continue; }
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      ind_abkl = key_ab * Chan.nhh[chan_kl] + key_kl;
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	for(int d = Space.num_hol; d < Space.num_states; ++d){
	  plus(cd, Space.qnums[c], Space.qnums[d]);
	  chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
	  if(chan_cd != chan_kl){ continue; }
	  key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
	  if(chan_c != chan_i){ continue; }
	  key_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];

	  chan_d = Space.ind_1b(Parameters.basis, Space.qnums[d]);
	  if(chan_d != chan_j){ continue; }
	  key_dj = Chan.ph1_map[Chan.ind0][Space.hash2(d, j, Chan.qnums2[Chan.ind0].j)];

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T = Amps.D1.T1_index[chan_ab];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T + ind_abkl];
	  T2 = Amps.S1.t2[key_ci];
	  T3 = Amps.S1.t2[key_dj];

	  term0 = 0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_7a += (1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * < '"<< a <<"','"<< b <<"' |t| "<< k <<","<< l <<" > ";
	    std::cout <<"* < "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	  term += term0;
	}
      }
    }
  }
  if(print != 0){ std::cout << "D_7a = " << term << std::endl; }
  return term;
}

void CC_Test_T2_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind2, chan2, chanind_T;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tempt;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind_T = Amps1.D1.T2_1_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      a = ph1.v1;
      j = ph1.v2;
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	i = hp1.v1;
	b = hp1.v2;
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	tempt = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T2_1 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T2_1(ab|ij)  <-  -X(kb|ic).T(ac|kj)
	// --------------------------------------------------
	// T2_1(ab|ij)  <-  -V(kb|ic).T(ac|kj)  (D_2e)
	tempt += Diagram_D_2e(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T2_1(ab|ij)  <-  (1/2).V(kl|cd).T(ad|kj).T(cb|il)  (D_3b)
	tempt += Diagram_D_3b(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T2_1(ab|ij)  <-  -V(kb|cd).T(ad|kj).t(c|i)  (D_5c)
	tempt += Diagram_D_5c(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T2_1(ab|ij)  <-  V(kl|ic).T(ac|kj).t(b|l)  (D_5d)
	tempt += Diagram_D_5d(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T2_1(ab|ij)  <-  V(kl|cd).T(ad|kj).t(c|i).t(b|l)  (D_7c)
	tempt += Diagram_D_7c(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);

	if(std::fabs(Amps2.D1.T2_1[chanind_T + ind2] - tempt) > 1e-12){
	  std::cout <<"!! T2_1 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T2_1[chanind_T + ind2] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  double X = Diagram_X_2e(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print, 1);
	  std::cout << " X_2e(1): " << X << ", D: " << tempt << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind2, chan2, chanind_T;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tempt;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind_T = Amps1.D1.T2_2_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      b = ph1.v1;
      i = ph1.v2;
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	j = hp1.v1;
	a = hp1.v2;
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	tempt = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T2_2 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T2_2(ab|ij)  <-  -X(ka|jc).T(bc|ki)
	// --------------------------------------------------
	// T2_2(ab|ij)  <-  -V(ka|jc).T(bc|ki)  (D_2e)
	tempt += Diagram_D_2e(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T2_2(ab|ij)  <-  (1/2).V(kl|cd).T(bd|ki).T(ca|jl)  (D_3b)
	tempt += Diagram_D_3b(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T2_2(ab|ij)  <-  -V(ka|cd).T(bd|ki).t(c|j)  (D_5c)
	tempt += Diagram_D_5c(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T2_2(ab|ij)  <-  V(kl|jc).T(bc|ki).t(a|l)  (D_5d)
	tempt += Diagram_D_5d(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T2_2(ab|ij)  <-  V(kl|cd).T(bd|ki).t(c|j).t(a|l)  (D_7c)
	tempt += Diagram_D_7c(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);

	if(std::fabs(Amps2.D1.T2_2[chanind_T + ind2] - tempt) > 1e-12){
	  std::cout <<"!! T2_2 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T2_2[chanind_T + ind2] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  double X = Diagram_X_2e(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, b, a, j, i, print, 1);
	  std::cout << " X_2e(2): " << X << ", D: " << tempt << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind2, chan2, chanind_T;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tempt;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind_T = Amps1.D1.T2_3_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      a = ph1.v1;
      i = ph1.v2;
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	j = hp1.v1;
	b = hp1.v2;
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	tempt = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T2_3 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T2_3(ab|ij)  <-  -X(kb|jc).T(ac|ki)
	// --------------------------------------------------
	// T2_3(ab|ij)  <-  -V(kb|jc).T(ac|ki)  (D_2e)
	tempt += Diagram_D_2e(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T2_3(ab|ij)  <-  (1/2).V(kl|cd).T(ad|ki).T(cb|jl)  (D_3b)
	tempt += Diagram_D_3b(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T2_3(ab|ij)  <-  -V(kb|cd).T(ad|ki).t(c|j)  (D_5c)
	tempt += Diagram_D_5c(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T2_3(ab|ij)  <-  V(kl|jc).T(ac|ki).t(b|l)  (D_5d)
	tempt += Diagram_D_5d(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T2_3(ab|ij)  <-  V(kl|cd).T(ad|ki).t(c|j).t(b|l)  (D_7c)
	tempt += Diagram_D_7c(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);

	if(std::fabs(Amps2.D1.T2_3[chanind_T + ind2] - tempt) > 1e-12){
	  std::cout <<"!! T2_3 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T2_3[chanind_T + ind2] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  double X = Diagram_X_2e(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, j, i, print, 2);
	  std::cout << " X_2e(3): " << X << ", D: " << tempt << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_4(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind2, chan2, chanind_T;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tempt;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind_T = Amps1.D1.T2_4_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      b = ph1.v1;
      j = ph1.v2;
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	i = hp1.v1;
	a = hp1.v2;
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	tempt = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T2_4 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T2_4(ab|ij)  <-  -X(ka|ic).T(bc|kj)
	// --------------------------------------------------
	// T2_4(ab|ij)  <-  -V(ka|ic).T(bc|kj)  (D_2e)
	tempt += Diagram_D_2e(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T2_4(ab|ij)  <-  (1/2).V(kl|cd).T(bd|kj).T(ca|il)  (D_3b)
	tempt += Diagram_D_3b(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T2_4(ab|ij)  <-  -V(ka|cd).T(bd|kj).t(c|i)  (D_5c)
	tempt += Diagram_D_5c(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T2_4(ab|ij)  <-  V(kl|ic).T(bc|kj).t(a|l)  (D_5d)
	tempt += Diagram_D_5d(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T2_4(ab|ij)  <-  V(kl|cd).T(bd|kj).t(c|i).t(a|l)  (D_7c)
	tempt += Diagram_D_7c(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);

	if(std::fabs(Amps2.D1.T2_4[chanind_T + ind2] - tempt) > 1e-12){
	  std::cout <<"!! T2_4 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T2_4[chanind_T + ind2] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  double X = Diagram_X_2e(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, b, a, i, j, print, 2);
	  std::cout << " X_2e(4): " << X << ", D: " << tempt << std::endl;
	}
      }
    }
  }
}

// T2_1(ab|ij)  <-  -V(kb|ic).T(ac|kj)  (D_2e)
double Diagram_D_2e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State aj, kc, ib;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  minus(ib, Space.qnums[i], Space.qnums[b]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_ib = Space.ind_2b_cross(Parameters.basis, ib);
  int chan_kc;
  if(chan_aj != chan_ib){ return 0.0; }
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_ib = Chan.hp1_map[chan_ib][Space.hash2(i, b, ib.j)];
  int key_kc;
  int ind_kcib, ind_ajkc;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int k = 0; k < Space.num_hol; ++k){
      minus(kc, Space.qnums[k], Space.qnums[c]);
      chan_kc = Space.ind_2b_cross(Parameters.basis, kc);
      if(chan_kc != chan_aj){ continue; }
      key_kc = Chan.hp1_map[chan_kc][Space.hash2(k, c, kc.j)];
      ind_kcib = key_kc * Chan.nhp1[chan_ib] + key_ib;
      ind_ajkc = key_aj * Chan.nhp1[chan_kc] + key_kc;

      chanind_V = Ints.Vhphp.V_2_1_index[chan_kc];
      chanind_T = Amps.D1.T2_1_index[chan_aj];
      V = Ints.Vhphp.V_2_1[chanind_V + ind_kcib];
      T = Amps.D1.T2_1[chanind_T + ind_ajkc];
      if(permute == 1){
	term0 = -1.0 * V * T;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_2e(1) += - < "<< k <<",'"<< b <<"' |V| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<","<< j << " > ";
	  std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
	}
      }
      else if(permute == 2){
	term0 = V * T;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_2e(2) += < "<< k <<",'"<< b <<"' |V| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<","<< j << " > ";
	  std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
	}
      }
      term += term0;
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_2e(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_2e(2) = " << term << std::endl; }
  }
  return term;
}

// T2_1(ab|ij)  <-  -X(kb|ic).T(ac|kj)  (X_2e)
double Diagram_X_2e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State aj, kc, ib;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  minus(ib, Space.qnums[i], Space.qnums[b]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_ib = Space.ind_2b_cross(Parameters.basis, ib);
  int chan_kc;
  if(chan_aj != chan_ib){ return 0.0; }
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_ib = Chan.hp1_map[chan_ib][Space.hash2(i, b, ib.j)];
  int key_kc;
  int ind_kcib, ind_ajkc;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;

  for(int c = Space.num_hol; c < Space.num_states; ++c){
    for(int k = 0; k < Space.num_hol; ++k){
      minus(kc, Space.qnums[k], Space.qnums[c]);
      chan_kc = Space.ind_2b_cross(Parameters.basis, kc);
      if(chan_kc != chan_aj){ continue; }
      key_kc = Chan.hp1_map[chan_kc][Space.hash2(k, c, kc.j)];
      ind_kcib = key_kc * Chan.nhp1[chan_ib] + key_ib;
      ind_ajkc = key_aj * Chan.nhp1[chan_kc] + key_kc;

      chanind_V = Eff_Ints.Xhphp.X_2_1_index[chan_kc];
      chanind_T = Amps.D1.T2_1_index[chan_aj];
      V = Eff_Ints.Xhphp.X_2_1[chanind_V + ind_kcib];
      T = Amps.D1.T2_1[chanind_T + ind_ajkc];
      if(permute == 1){
	term0 = -1.0 * V * T;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"X_2e(1) += - < "<< k <<",'"<< b <<"' |X| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<","<< j << " > ";
	  std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
	}
      }
      else if(permute == 2){
	term0 = V * T;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"X_2e(2) += < "<< k <<",'"<< b <<"' |X| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<","<< j << " > ";
	  std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
	}
      }
      term += term0;
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "X_2e(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "X_2e(2) = " << term << std::endl; }
  }
  return term;
}

// T2_1(ab|ij)  <-  (1/2).V(kl|cd).T(ad|kj).T(cb|il)  (D_3b)
double Diagram_D_3b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State aj, kd, cl, ib;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  minus(ib, Space.qnums[i], Space.qnums[b]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_ib = Space.ind_2b_cross(Parameters.basis, ib);
  int chan_kd, chan_cl;
  if(chan_aj != chan_ib){ return 0.0; }
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_ib = Chan.hp1_map[chan_ib][Space.hash2(i, b, ib.j)];
  int key_kd, key_cl;
  int ind_kdcl, ind_ajkd, ind_clib;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T1, chanind_T2;
  double V, T1, T2;
  for(int d = Space.num_hol; d < Space.num_states; ++d){
    for(int k = 0; k < Space.num_hol; ++k){
      minus(kd, Space.qnums[k], Space.qnums[d]);
      chan_kd = Space.ind_2b_cross(Parameters.basis, kd);
      if(chan_kd != chan_aj){ continue; }
      key_kd = Chan.hp1_map[chan_kd][Space.hash2(k, d, kd.j)];
      ind_ajkd = key_aj * Chan.nhp1[chan_kd] + key_kd;
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	for(int l = 0; l < Space.num_hol; ++l){
	  minus(cl, Space.qnums[c], Space.qnums[l]);
	  chan_cl = Space.ind_2b_cross(Parameters.basis, cl);
	  if(chan_cl != chan_ib){ continue; }
	  key_cl = Chan.ph1_map[chan_cl][Space.hash2(c, l, cl.j)];
	  ind_clib = key_cl * Chan.nhp1[chan_ib] + key_ib;
	  ind_kdcl = key_kd * Chan.nph1[chan_cl] + key_cl;

	  chanind_V = Ints.Vhhpp.V_2_1_index[chan_kd];
	  chanind_T1 = Amps.D1.T2_1_index[chan_aj];
	  chanind_T2 = Amps.D1.T2_1_index[chan_cl];
	  V = Ints.Vhhpp.V_2_1[chanind_V + ind_kdcl];
	  T1 = Amps.D1.T2_1[chanind_T1 + ind_ajkd];
	  T2 = Amps.D1.T2_1[chanind_T2 + ind_clib];
	  if(permute == 1){
	    term0 = 0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3b(1) += (1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"',"<< d <<" |t| "<< k <<","<< j << " > * < "<< c <<",'"<< b <<"' |t| '"<< i <<"',"<< l << " > ";
	      std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = -0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3b(2) += -(1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"',"<< d <<" |t| "<< k <<","<< j << " > * < "<< c <<",'"<< b <<"' |t| '"<< i <<"',"<< l << " > ";
	      std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_3b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_3b(2) = " << term << std::endl; }
  }
  return term;
}

// T2_1(ab|ij)  <-  -V(kb|cd).T(ad|kj).t(c|i)  (D_5c)
double Diagram_D_5c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State bd, ck, aj, kd, ci;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_bd, chan_ck, chan_kd, chan_c;
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_bd, key_ck, key_kd, key_ci;
  int ind_bdck, ind_ajkd;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_i){ continue; }
    key_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      minus(bd, Space.qnums[b], Space.qnums[d]);
      chan_bd = Space.ind_2b_cross(Parameters.basis, bd);
      key_bd = Chan.pp1_map[chan_bd][Space.hash2(b, d, bd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	minus(ck, Space.qnums[c], Space.qnums[k]);
	chan_ck = Space.ind_2b_cross(Parameters.basis, ck);
	if(chan_ck != chan_bd){ continue; }
	key_ck = Chan.ph1_map[chan_ck][Space.hash2(c, k, ck.j)];
	ind_bdck = key_bd * Chan.nph1[chan_ck] + key_ck;
	minus(kd, Space.qnums[k], Space.qnums[d]);
	chan_kd = Space.ind_2b_cross(Parameters.basis, kd);
	if(chan_kd != chan_aj){ continue; }
	key_kd = Chan.hp1_map[chan_kd][Space.hash2(k, d, kd.j)];
	ind_ajkd = key_aj * Chan.nhp1[chan_kd] + key_kd;

	chanind_V = Ints.Vhppp.V_2_4_index[chan_bd];
	chanind_T = Amps.D1.T2_1_index[chan_aj];
	V = Ints.Vhppp.V_2_4[chanind_V + ind_bdck];
	T1 = Amps.D1.T2_1[chanind_T + ind_ajkd];
	T2 = Amps.S1.t2[key_ci];
	if(permute == 1){
	  term0 = -1.0 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5c(1) += - < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * < '"<< a <<"',"<< d <<" |t| "<< k <<",'"<< j <<"' > * < "<< c <<" |t| '"<< i <<"' > ";
	    std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5c(2) += < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * < '"<< a <<"',"<< d <<" |t| "<< k <<",'"<< j <<"' > * < "<< c <<" |t| '"<< i <<"' > ";
	    std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5c(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5c(2) = " << term << std::endl; }
  }
  return term;
}

// T2_1(ab|ij)  <-  V(kl|ic).T(ac|kj).t(b|l)  (D_5d)
double Diagram_D_5d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ki, cl, aj, kc, bl;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_ki, chan_cl, chan_kc, chan_l;
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_ki, key_cl, key_kc, key_bl;
  int ind_kicl, ind_ajkc;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int l = 0; l < Space.num_hol; ++l){
    chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
    if(chan_l != chan_b){ continue; }
    key_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];
    for(int k = 0; k < Space.num_hol; ++k){
      minus(ki, Space.qnums[k], Space.qnums[i]);
      chan_ki = Space.ind_2b_cross(Parameters.basis, ki);
      key_ki = Chan.hh1_map[chan_ki][Space.hash2(k, i, ki.j)];
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	minus(cl, Space.qnums[c], Space.qnums[l]);
	chan_cl = Space.ind_2b_cross(Parameters.basis, cl);
	if(chan_cl != chan_ki){ continue; }
	key_cl = Chan.ph1_map[chan_cl][Space.hash2(c, l, cl.j)];
	ind_kicl = key_ki * Chan.nph1[chan_cl] + key_cl;
	minus(kc, Space.qnums[k], Space.qnums[c]);
	chan_kc = Space.ind_2b_cross(Parameters.basis, kc);
	if(chan_kc != chan_aj){ continue; }
	key_kc = Chan.hp1_map[chan_kc][Space.hash2(k, c, kc.j)];
	ind_ajkc = key_aj * Chan.nhp1[chan_kc] + key_kc;

	chanind_V = Ints.Vhhhp.V_2_3_index[chan_ki];
	chanind_T = Amps.D1.T2_1_index[chan_aj];
	V = Ints.Vhhhp.V_2_3[chanind_V + ind_kicl];
	T1 = Amps.D1.T2_1[chanind_T + ind_ajkc];
	T2 = Amps.S1.t2[key_bl];
	if(permute == 1){
	  term0 = V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5d(1) += < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<",'"<< j <<"' > * < '"<< b <<"' |t| "<< l <<" > ";
	    std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = -1.0 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5d(2) += - < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * < '"<< a <<"',"<< c <<" |t| "<< k <<",'"<< j <<"' > * < '"<< b <<"' |t| "<< l <<" > ";
	    std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5d(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5d(2) = " << term << std::endl; }
  }
  return term;
}

// T2_1(ab|ij)  <-  V(kl|cd).T(ad|kj).t(c|i).t(b|l)  (D_7c)
double Diagram_D_7c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kd, cl, aj;
  minus(aj, Space.qnums[a], Space.qnums[j]);
  int chan_aj = Space.ind_2b_cross(Parameters.basis, aj);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_kd, chan_cl, chan_c, chan_l;
  int key_aj = Chan.ph1_map[chan_aj][Space.hash2(a, j, aj.j)];
  int key_kd, key_cl, key_ci, key_bl;
  int ind_kdcl, ind_ajkd;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2, T3;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_i){ continue; }
    key_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];
    for(int l = 0; l < Space.num_hol; ++l){
      chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
      if(chan_l != chan_b){ continue; }
      key_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];
      minus(cl, Space.qnums[c], Space.qnums[l]);
      chan_cl = Space.ind_2b_cross(Parameters.basis, cl);
      key_cl = Chan.ph1_map[chan_cl][Space.hash2(c, l, cl.j)];
      for(int d = Space.num_hol; d < Space.num_states; ++d){
	for(int k = 0; k < Space.num_hol; ++k){
	  minus(kd, Space.qnums[k], Space.qnums[d]);
	  chan_kd = Space.ind_2b_cross(Parameters.basis, kd);
	  if(chan_kd != chan_cl || chan_kd != chan_aj){ continue; }
	  key_kd = Chan.hp1_map[chan_kd][Space.hash2(k, d, kd.j)];
	  ind_kdcl = key_kd * Chan.nph1[chan_cl] + key_cl;
	  ind_ajkd = key_aj * Chan.nhp1[chan_kd] + key_kd;

	  chanind_V = Ints.Vhhpp.V_2_1_index[chan_kd];
	  chanind_T = Amps.D1.T2_1_index[chan_aj];
	  V = Ints.Vhhpp.V_2_1[chanind_V + ind_kdcl];
	  T1 = Amps.D1.T2_1[chanind_T + ind_ajkd];
	  T2 = Amps.S1.t2[key_ci];
	  T3 = Amps.S1.t2[key_bl];
	  if(permute == 1){
	    term0 = V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7c(1) += < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * < '"<< a <<"',"<< d <<" |t| "<< k <<",'"<< j <<"' > * ";
	      std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < '"<< b <<"' |t| "<< l <<" > ";
	      std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = -1.0 * V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7c(2) += - < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * < '"<< a <<"',"<< d <<" |t| "<< k <<",'"<< j <<"' > * ";
	      std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < '"<< b <<"' |t| "<< l <<" > ";
	      std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_7c(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_7c(2) = " << term << std::endl; }
  }
  return term;
}

void CC_Test_T3_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind3, chan3, chanind_T;
  int np, nhhp;
  int a, b, i, j;
  one_body p;
  three_body hhp;
  double tempt, D1, D2;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    chanind_T = Amps1.D1.T3_1_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
	hhp = Chan.hhp_state(chan3, hhp_ind);
	i = hhp.v1;
	j = hhp.v2;
	b = hhp.v3;
	if(a == b || i == j){ continue; }
	ind3 = p_ind * nhhp + hhp_ind;

	D1 = 0.0;
	D2 = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T3_1 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T3_1(ab|ij)  <-  X(a|c).T(cb|ij)  (X_2a)
	// --------------------------------------------------
	// T3_1(ab|ij)  <-  -(1/2).V(kl|cd).T(ac|lk).T(db|ij)  (D_3d)
	D1 += Diagram_D_3d(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_1(ab|ij)  <-  V(ka|cd).T(db|ij).t(c|k)  (D_5g)
	D1 += Diagram_D_5g(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_1(ab|ij)  <-  -V(kl|cd).T(db|ij).t(a|l).t(c|k)  (D_7e)
	D1 += Diagram_D_7e(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);


	// T3_1(ab|ij)  <-  -X'(kb|ij).t(a|k)  (X_4b)
	// --------------------------------------------------
	// T3_1(ab|ij)  <-  -V(kb|ij).t(a|k)  (D_4b)
	D2 += Diagram_D_4b(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_1(ab|ij)  <-  (1/2).V(kl|ij).t(a|k).t(b|l)  (D_6b)
	D2 += Diagram_D_6b(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);

	tempt = D1 + D2;
	if(std::fabs(Amps2.D1.T3_1[chanind_T + ind3] - tempt) > 1e-12){
	  double X1, X2;
	  std::cout <<"!! T3_1 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T3_1[chanind_T + ind3] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  X1 = Diagram_X_2a(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print, 1);
	  X2 = Diagram_X_4b(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print, 1);
	  std::cout << " X_2a: " << X1 << ", D: " << D1 << std::endl;
	  std::cout << " X_4b: " << X2 << ", D: " << D2 << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind3, chan3, chanind_T;
  int np, nhhp;
  int a, b, i, j;
  one_body p;
  three_body hhp;
  double tempt, D1, D2;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    chanind_T = Amps1.D1.T3_2_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      b = p.v1;
      for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
	hhp = Chan.hhp_state(chan3, hhp_ind);
	i = hhp.v1;
	j = hhp.v2;
	a = hhp.v3;
	if(a == b || i == j){ continue; }
	ind3 = p_ind * nhhp + hhp_ind;

	D1 = 0.0;
	D2 = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T3_2 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T3_2(ab|ij)  <-  X(b|c).T(ca|ij)  (X_2a)
	// --------------------------------------------------
	// T3_2(ab|ij)  <-  -(1/2).V(kl|cd).T(bc|lk).T(da|ij)  (D_3d)
	D1 += Diagram_D_3d(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_2(ab|ij)  <-  V(kb|cd).T(da|ij).t(c|k)  (D_5g)
	D1 += Diagram_D_5g(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_2(ab|ij)  <-  -V(kl|cd).T(da|ij).t(b|l).t(c|k)  (D_7e)
	D1 += Diagram_D_7e(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);


	// T3_2(ab|ij)  <-  -X'(ka|ij).t(b|k)  (X_4b)
	// --------------------------------------------------
	// T3_2(ab|ij)  <-  -V(ka|ij).t(b|k)  (D_4b)
	D2 += Diagram_D_4b(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_2(ab|ij)  <-  (1/2).V(kl|ij).t(b|k).t(a|l)  (D_6b)
	D2 += Diagram_D_6b(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);

	tempt = D1 + D2;
	if(std::fabs(Amps2.D1.T3_2[chanind_T + ind3] - tempt) > 1e-12){
	  double X1, X2;
	  std::cout <<"!! T3_2 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T3_2[chanind_T + ind3] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  X1 = Diagram_X_2a(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, b, a, i, j, print, 2);
	  X2 = Diagram_X_4b(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, b, a, i, j, print, 2);
	  std::cout << " X_2a: " << X1 << ", D: " << D1 << std::endl;
	  std::cout << " X_4b: " << X2 << ", D: " << D2 << std::endl;
	}
      }
    }
  }
}

// T3_1(ab|ij)  <-  X(a|c).T(cb|ij)  (X_2a)
double Diagram_X_2a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State cb, ij;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_cb, chan_c;
  int key_a = Chan.p_map[chan_a][a];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_c, key_cb;
  int ind_cbij, ind_ac;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_a){ continue; }
    key_c = Chan.p_map[chan_c][c];
    ind_ac = key_a * Chan.np[chan_c] + key_c;

    plus(cb, Space.qnums[c], Space.qnums[b]);
    chan_cb = Space.ind_2b_dir(Parameters.basis, cb);
    if(chan_cb != chan_ij){ continue; }
    key_cb = Chan.pp_map[chan_cb][Space.hash2(c, b, cb.j)];
    ind_cbij = key_cb * Chan.nhh[chan_ij] + key_ij;

    chanind_V = Eff_Ints.Xpp.X_3_index[chan_a];
    chanind_T = Amps.D1.T1_index[chan_cb];
    V = Eff_Ints.Xpp.X_3od[chanind_V + ind_ac];
    T = Amps.D1.T1[chanind_T + ind_cbij];
    if(permute == 1){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2a(1) += < '"<< a <<"' |X| "<< c <<" > * < "<< c <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2a(2) += - < '"<< a <<"' |X| "<< c <<" > * < "<< c <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "X_2a(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "X_2a(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  (1/2).V(kl|cd).T(ac|kl).T(db|ij)  (D_3d)
double Diagram_D_3d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, cd, ac, db, ij;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_kl, chan_cd, chan_ac, chan_db;
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl, key_cd, key_ac, key_db;
  int ind_klcd, ind_ackl, ind_dbij;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T1, chanind_T2;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    plus(ac, Space.qnums[a], Space.qnums[c]);
    chan_ac = Space.ind_2b_dir(Parameters.basis, ac);
    key_ac = Chan.pp_map[chan_ac][Space.hash2(a, c, ac.j)];
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      if(d == a){ continue; }
      plus(db, Space.qnums[d], Space.qnums[b]);
      chan_db = Space.ind_2b_dir(Parameters.basis, db);
      if(chan_db != chan_ij){ continue; }
      key_db = Chan.pp_map[chan_db][Space.hash2(d, b, db.j)];
      ind_dbij = key_db * Chan.nhh[chan_ij] + key_ij;
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	for(int l = 0; l < Space.num_hol; ++l){
	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_cd || chan_kl != chan_ac){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;
	  ind_ackl = key_ac * Chan.nhh[chan_kl] + key_kl;

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T1 = Amps.D1.T1_index[chan_ac];
	  chanind_T2 = Amps.D1.T1_index[chan_db];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T1 + ind_ackl];
	  T2 = Amps.D1.T1[chanind_T2 + ind_dbij];
	  if(permute == 1){
	    term0 = 0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3d(1) += (1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"',"<< c <<" |t| "<< k <<","<< l <<" > * < "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > ";
	      std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = -0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3d(2) += -(1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"',"<< c <<" |t| "<< k <<","<< l <<" > * < "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > ";
	      std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_3d(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_3d(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  V(ka|cd).T(db|ij).t(c|k)  (D_5g)
double Diagram_D_5g(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State cd, cdk, db, ij;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_cdk, chan_db, chan_c, chan_k;
  int key_a = Chan.p_map[chan_a][a];
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_cdk, key_db;
  int ind_acdk, ind_dbij, ind_ck;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int d = Space.num_hol; d < Space.num_states; ++d){
    if(d == a){ continue; }
    plus(db, Space.qnums[d], Space.qnums[b]);
    chan_db = Space.ind_2b_dir(Parameters.basis, db);
    if(chan_db != chan_ij){ continue; }
    key_db = Chan.pp_map[chan_db][Space.hash2(d, b, db.j)];
    ind_dbij = key_db * Chan.nhh[chan_ij] + key_ij;
    for(int c = Space.num_hol; c < Space.num_states; ++c){
      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
      for(int k = 0; k < Space.num_hol; ++k){
	minus(cdk, cd, Space.qnums[k]);
	chan_cdk = Space.ind_1b(Parameters.basis, cdk);
	if(chan_cdk != chan_a){ continue; }
	key_cdk = Chan.pph_map[chan_cdk][Space.hash3(c, d, k, cd.j)];
	ind_acdk = key_a * Chan.npph[chan_cdk] + key_cdk;

	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_c){ continue; }
	ind_ck = Chan.ph1_map[Chan.ind0][Space.hash2(c, k, Chan.qnums2[Chan.ind0].j)];

	chanind_V = Ints.Vhppp.V_3_2_index[chan_a];
	chanind_T = Amps.D1.T1_index[chan_db];
	V = Ints.Vhppp.V_3_2[chanind_V + ind_acdk];
	T1 = Amps.D1.T1[chanind_T + ind_dbij];
	T2 = Amps.S1.t2[ind_ck];
	if(permute == 1){
	  term0 = V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5g(1) += < "<< k <<",'"<< a <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > * < "<< c <<" |t| "<< k <<" > ";
	    std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = -1.0 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5g(2) += - < "<< k <<",'"<< a <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > * < "<< c <<" |t| "<< k <<" > ";
	    std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5g(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5g(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  -V(kl|cd).T(db|ij).t(a|l).t(c|k)  (D_7e)
double Diagram_D_7e(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, cd, db, ij;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_kl, chan_cd, chan_db, chan_l, chan_c, chan_k;
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl, key_cd, key_db;
  int ind_klcd, ind_dbij, ind_al, ind_ck;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2, T3;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      if(d == a){ continue; }
      plus(db, Space.qnums[d], Space.qnums[b]);
      chan_db = Space.ind_2b_dir(Parameters.basis, db);
      if(chan_db != chan_ij){ continue; }
      key_db = Chan.pp_map[chan_db][Space.hash2(d, b, db.j)];
      ind_dbij = key_db * Chan.nhh[chan_ij] + key_ij;

      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_c){ continue; }
	ind_ck = Chan.ph1_map[Chan.ind0][Space.hash2(c, k, Chan.qnums2[Chan.ind0].j)];
	for(int l = 0; l < Space.num_hol; ++l){
	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_cd){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
	  if(chan_l != chan_a){ continue; }
	  ind_al = Chan.ph1_map[Chan.ind0][Space.hash2(a, l, Chan.qnums2[Chan.ind0].j)];

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T = Amps.D1.T1_index[chan_db];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T + ind_dbij];
	  T2 = Amps.S1.t2[ind_al];
	  T3 = Amps.S1.t2[ind_ck];
	  if(permute == 1){
	    term0 = -1.0 * V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7e(1) += - < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< l <<" > * < "<< c <<" |t| "<< k <<" > ";
	      std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7e(2) += < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< "<< d <<",'"<< b <<"' |t| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< l <<" > * < "<< c <<" |t| "<< k <<" > ";
	      std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_7e(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_7e(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  -V(kb|ij).t(a|k)  (D_4b)
double Diagram_D_4b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ij, ijk;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_ijk, chan_k;
  int key_b = Chan.p_map[chan_b][b];
  int key_ijk;
  int ind_bijk, ind_ak;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T;
  for(int k = 0; k < Space.num_hol; ++k){
    minus(ijk, ij, Space.qnums[k]);
    chan_ijk = Space.ind_1b(Parameters.basis, ijk);
    if(chan_ijk != chan_b){ continue; }
    key_ijk = Chan.hhh_map[chan_ijk][Space.hash3(i, j, k, ij.j)];
    ind_bijk = key_b * Chan.nhhh[chan_ijk] + key_ijk;
    
    chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
    if(chan_k != chan_a){ continue; }
    ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];

    chanind_V = Ints.Vhphh.V_3_2_index[chan_b];    
    V = Ints.Vhphh.V_3_2[chanind_V + ind_bijk];
    T = Amps.S1.t2[ind_ak];
    if(permute == 1){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_4b(1) += - < "<< k <<",'"<< b <<"' |V| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< k <<" > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_4b(2) += < "<< k <<",'"<< b <<"' |V| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< k <<" > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_4b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_4b(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  -X'(kb|ij).t(a|k)  (X_4b)
double Diagram_X_4b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ij, ijk;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_ijk, chan_k;
  int key_b = Chan.p_map[chan_b][b];
  int key_ijk;
  int ind_bijk, ind_ak;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T;
  for(int k = 0; k < Space.num_hol; ++k){
    minus(ijk, ij, Space.qnums[k]);
    chan_ijk = Space.ind_1b(Parameters.basis, ijk);
    if(chan_ijk != chan_b){ continue; }
    key_ijk = Chan.hhh_map[chan_ijk][Space.hash3(i, j, k, ij.j)];
    ind_bijk = key_b * Chan.nhhh[chan_ijk] + key_ijk;
    
    chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
    if(chan_k != chan_a){ continue; }
    ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];
    
    chanind_V = Eff_Ints.Xhphh.X_3_2[chan_b];
    V = Eff_Ints.Xhphh.X1_3_2[chanind_V + ind_bijk];
    T = Amps.S1.t2[ind_ak];
    if(permute == 1){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_4b(1) += - < "<< k <<",'"<< b <<"' |X'| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< k <<" > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_4b(2) += < "<< k <<",'"<< b <<"' |X'| '"<< i <<"','"<< j <<"' > * < '"<< a <<"' |t| "<< k <<" > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "X_4b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "X_4b(2) = " << term << std::endl; }
  }
  return term;
}

// T3_1(ab|ij)  <-  (1/2).V(kl|ij).t(a|k).t(b|l)  (D_6b)
double Diagram_D_6b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, ij;
  plus(ij, Space.qnums[i], Space.qnums[j]);
  int chan_ij = Space.ind_2b_dir(Parameters.basis, ij);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_kl, chan_k, chan_l;
  int key_ij = Chan.hh_map[chan_ij][Space.hash2(i, j, ij.j)];
  int key_kl;
  int ind_klij, ind_ak, ind_bl;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2;
  for(int k = 0; k < Space.num_hol; ++k){
    chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
    if(chan_k != chan_a){ continue; }
    ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];
    for(int l = 0; l < Space.num_hol; ++l){
      chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
      if(chan_l != chan_b){ continue; }
      ind_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];

      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      if(chan_kl != chan_ij){ continue; }
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      ind_klij = key_kl * Chan.nhh[chan_ij] + key_ij;

      chanind_V = Ints.Vhhhh.V_1_index[chan_kl];
      V = Ints.Vhhhh.V_1[chanind_V + ind_klij];
      T1 = Amps.S1.t2[ind_ak];
      T2 = Amps.S1.t2[ind_bl];
      if(permute == 1){
	term0 = 0.5 * V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6b(1) += (1/2) * < "<< k <<","<< l <<" |V| '"<< i <<"','"<< j <<"' > * ";
	  std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > ";
	  std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      else if(permute == 2){
	term0 = -0.5 * V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6b(2) += -(1/2) * < "<< k <<","<< l <<" |V| '"<< i <<"','"<< j <<"' > * ";
	  std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > ";
	  std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      term += term0;
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_6b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_6b(2) = " << term << std::endl; }
  }
  return term;
}

void CC_Test_T3_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind3, chan3, chanind_T;
  int npph, nh;
  int a, b, i, j;
  three_body pph;
  one_body h;
  double tempt, D1, D2;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    npph = Chan.npph[chan3];
    nh = Chan.nh[chan3];
    chanind_T = Amps1.D1.T3_3_index[chan3];
    for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
      pph = Chan.pph_state(chan3, pph_ind);
      a = pph.v1;
      b = pph.v2;
      j = pph.v3;
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	i = h.v1;
	if(a == b || i == j){ continue; }
	ind3 = pph_ind * nh + h_ind;

	D1 = 0.0;
	D2 = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T3_3 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T3_3(ab|ij)  <-  -X(k|i).T(ab|kj)  (X_2b)
	// --------------------------------------------------
	// T3_3(ab|ij)  <-  (1/2).V(kl|cd).T(cd|ik).T(ab|lj)  (D_3c)
	D1 += Diagram_D_3c(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_3(ab|ij)  <-  -V(kl|ci).T(ab|lj).t(c|k)  (D_5h)
	D1 += Diagram_D_5h(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_3(ab|ij)  <-  -V(kl|cd).T(ab|lj).t(d|i).t(c|k)  (D_7d)
	D1 += Diagram_D_7d(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);

	// T3_3(ab|ij)  <-  X'(ab|jc).t(c|i)  (X_4a)
	// --------------------------------------------------
	// T3_3(ab|ij)  <-  V(ab|jc).t(c|i)  (D_4a)
	D2 += Diagram_D_4a(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T3_3(ab|ij)  <-  (1/2).V(ab|cd).t(c|j).t(d|i)  (D_6a)
	D2 += Diagram_D_6a(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T3_3(ab|ij)  <-  -P(ab).V(kb|jc).t(a|k).t(c|i)  (D_6c)
	D2 += Diagram_D_6c(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	D2 += Diagram_D_6c(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T3_3(ab|ij)  <-  -(1/2).P(ab).V(kb|cd).t(a|k).t(c|j).t(d|i)  (D_8a)
	D2 += Diagram_D_8a(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	D2 += Diagram_D_8a(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T3_3(ab|ij)  <-  (1/2).P(ab).V(kl|jc).t(a|k).t(b|l).t(c|i)  (D_8b)
	D2 += Diagram_D_8b(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	D2 += Diagram_D_8b(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);
	// T3_3(ab|ij)  <-  (1/4).P(ab).V(kl|cd).t(a|k).t(b|l).t(c|j).t(d|i)  (D_9)
	D2 += Diagram_D_9(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	D2 += Diagram_D_9(Parameters, Space, Chan, Ints, Amps1, b, a, j, i, print, 1);

	tempt = D1 + D2;
	if(std::fabs(Amps2.D1.T3_3[chanind_T + ind3] - tempt) > 1e-12){
	  double X1, X2;
	  std::cout <<"!! T3_3 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T3_3[chanind_T + ind3] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  X1 = Diagram_X_2b(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print, 1);
	  X2 = Diagram_X_4a(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, i, j, print, 1);
	  std::cout << " X_2b(1): " << X1 << ", D: " << D1 << std::endl;
	  std::cout << " X_4a(1): " << X2 << ", D: " << D2 << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_4(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2, int &print)
{
  int ind3, chan3, chanind_T;
  int npph, nh;
  int a, b, i, j;
  three_body pph;
  one_body h;
  double tempt, D1, D2;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    npph = Chan.npph[chan3];
    nh = Chan.nh[chan3];
    chanind_T = Amps1.D1.T3_4_index[chan3];
    for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
      pph = Chan.pph_state(chan3, pph_ind);
      a = pph.v1;
      b = pph.v2;
      i = pph.v3;
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	j = h.v1;
	if(a == b || i == j){ continue; }
	ind3 = pph_ind * nh + h_ind;

	D1 = 0.0;
	D2 = 0.0;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! T3_4 < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	}

	// T3_4(ab|ij)  <-  -X(k|j).T(ab|ki)  (X_2b)
	// --------------------------------------------------
	// T3_4(ab|ij)  <-  (1/2).V(kl|cd).T(cd|jk).T(ab|li)  (D_3c)
	D1 += Diagram_D_3c(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T3_4(ab|ij)  <-  -V(kl|cj).T(ab|li).t(c|k)  (D_5h)
	D1 += Diagram_D_5h(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);
	// T3_4(ab|ij)  <-  -V(kl|cd).T(ab|li).t(d|j).t(c|k)  (D_7d)
	D1 += Diagram_D_7d(Parameters, Space, Chan, Ints, Amps1, a, b, j, i, print, 2);

	// T3_4(ab|ij)  <-  X'(ab|ic).t(c|j)  (X_4a)
	// --------------------------------------------------
	// T3_4(ab|ij)  <-  V(ab|ic).t(c|j)  (D_4a)
	D2 += Diagram_D_4a(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_4(ab|ij)  <-  (1/2).V(ab|cd).t(c|i).t(d|j)  (D_6a)
	D2 += Diagram_D_6a(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	// T3_4(ab|ij)  <-  -P(ab).V(kb|ic).t(a|k).t(c|j)  (D_6c)
	D2 += Diagram_D_6c(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D2 += Diagram_D_6c(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_4(ab|ij)  <-  -(1/2).P(ab).V(kb|cd).t(a|k).t(c|i).t(d|j)  (D_8a)
	D2 += Diagram_D_8a(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D2 += Diagram_D_8a(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_4(ab|ij)  <-  (1/2).P(ab).V(kl|ic).t(a|k).t(b|l).t(c|j)  (D_8b)
	D2 += Diagram_D_8b(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D2 += Diagram_D_8b(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);
	// T3_4(ab|ij)  <-  (1/4).P(ab).V(kl|cd).t(a|k).t(b|l).t(c|i).t(d|j)  (D_9)
	D2 += Diagram_D_9(Parameters, Space, Chan, Ints, Amps1, a, b, i, j, print, 1);
	D2 += Diagram_D_9(Parameters, Space, Chan, Ints, Amps1, b, a, i, j, print, 2);

	tempt = D1 + D2;
	if(std::fabs(Amps2.D1.T3_4[chanind_T + ind3] - tempt) > 1e-12){
	  double X1, X2;
	  std::cout <<"!! T3_4 < '"<< a <<"','"<< b <<"' |t| '"<< i <<"','"<< j <<"' > = X: "<< Amps2.D1.T3_4[chanind_T + ind3] <<", ";
	  std::cout <<"D: "<< tempt << std::endl;
	  X1 = Diagram_X_2b(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, j, i, print, 2);
	  X2 = Diagram_X_4a(Parameters, Space, Chan, Ints, Eff_Ints, Amps1, a, b, j, i, print, 2);
	  std::cout << " X_2b(2): " << X1 << ", D: " << D1 << std::endl;
	  std::cout << " X_4a(2): " << X2 << ", D: " << D2 << std::endl;
	}
      }
    }
  }
}

// T3_3(ab|ij)  <-  -X(k|i).T(ab|kj)  (X_2b)
double Diagram_X_2b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ab, kj;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_kj, chan_k;
  int key_i = Chan.h_map[chan_i][i];
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_k, key_kj;
  int ind_abkj, ind_ki;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T;
  for(int k = 0; k < Space.num_hol; ++k){
    chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
    if(chan_k != chan_i){ continue; }
    key_k = Chan.h_map[chan_k][k];
    ind_ki = key_k * Chan.nh[chan_i] + key_i;

    plus(kj, Space.qnums[k], Space.qnums[j]);
    chan_kj = Space.ind_2b_dir(Parameters.basis, kj);
    if(chan_kj != chan_ab){ continue; }
    key_kj = Chan.hh_map[chan_kj][Space.hash2(k, j, kj.j)];
    ind_abkj = key_ab * Chan.nhh[chan_kj] + key_kj;

    chanind_V = Eff_Ints.Xhh.X_3_index[chan_k];
    chanind_T = Amps.D1.T1_index[chan_ab];
    V = Eff_Ints.Xhh.X_3od[chanind_V + ind_ki];
    T = Amps.D1.T1[chanind_T + ind_abkj];
    if(permute == 1){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2b(1) += - < "<< k <<" |X| '"<< i <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<",'"<< j <<"' > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_2b(2) += < "<< k <<" |X| '"<< i <<"' > * < '"<< a <<"','"<< b <<"' |t| "<< k <<",'"<< j <<"' > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "X_2b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "X_2b(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  (1/2).V(kl|cd).T(cd|ik).T(ab|lj)  (D_3c)
double Diagram_D_3c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, cd, ik, ab, lj;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_kl, chan_cd, chan_ik, chan_lj;
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_kl, key_cd, key_ik, key_lj;
  int ind_klcd, ind_cdik, ind_ablj;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T1, chanind_T2;
  double V, T1, T2;
  for(int k = 0; k < Space.num_hol; ++k){
    plus(ik, Space.qnums[i], Space.qnums[k]);
    chan_ik = Space.ind_2b_dir(Parameters.basis, ik);
    key_ik = Chan.hh_map[chan_ik][Space.hash2(i, k, ik.j)];
    for(int l = 0; l < Space.num_hol; ++l){
      if(l == i){ continue; }
      plus(lj, Space.qnums[l], Space.qnums[j]);
      chan_lj = Space.ind_2b_dir(Parameters.basis, lj);
      if(chan_lj != chan_ab){ continue; }
      key_lj = Chan.hh_map[chan_lj][Space.hash2(l, j, lj.j)];
      ind_ablj = key_ab * Chan.nhh[chan_lj] + key_lj;
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
      key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	for(int d = Space.num_hol; d < Space.num_states; ++d){
	  plus(cd, Space.qnums[c], Space.qnums[d]);
	  chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
	  if(chan_cd != chan_kl || chan_cd != chan_ik){ continue; }
	  key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;
	  ind_cdik = key_cd * Chan.nhh[chan_ik] + key_ik;

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T1 = Amps.D1.T1_index[chan_cd];
	  chanind_T2 = Amps.D1.T1_index[chan_ab];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T1 + ind_cdik];
	  T2 = Amps.D1.T1[chanind_T2 + ind_ablj];
	  if(permute == 1){
	    term0 = 0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3c(1) += (1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< "<< c <<","<< d <<" |t| '"<< i <<"',"<< k <<" > * < '"<< a <<"','"<< b <<"' |t| "<< l <<",'"<< j <<"' > ";
	      std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = -0.5 * V * T1 * T2;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_3c(2) += -(1/2) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< "<< c <<","<< d <<" |t| '"<< i <<"',"<< k <<" > * < '"<< a <<"','"<< b <<"' |t| "<< l <<",'"<< j <<"' > ";
	      std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_3c(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_3c(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  -V(kl|ic).T(ab|kj).t(c|l)  (D_5h)
double Diagram_D_5h(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, klc, ab, kj;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_klc, chan_kj, chan_c, chan_l;
  int key_i = Chan.h_map[chan_i][i];
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_klc, key_kj;
  int ind_klci, ind_abkj, ind_cl;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2;
  for(int k = 0; k < Space.num_hol; ++k){
    if(k == i){ continue; }
    plus(kj, Space.qnums[k], Space.qnums[j]);
    chan_kj = Space.ind_2b_dir(Parameters.basis, kj);
    if(chan_kj != chan_ab){ continue; }
    key_kj = Chan.hh_map[chan_kj][Space.hash2(k, j, kj.j)];
    ind_abkj = key_ab * Chan.nhh[chan_kj] + key_kj;
    for(int l = 0; l < Space.num_hol; ++l){
      plus(kl, Space.qnums[k], Space.qnums[l]);
      chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	minus(klc, kl, Space.qnums[c]);
	chan_klc = Space.ind_1b(Parameters.basis, klc);
	if(chan_klc != chan_i){ continue; }
	key_klc = Chan.hhp_map[chan_klc][Space.hash3(k, l, c, kl.j)];
	ind_klci = key_klc * Chan.nh[chan_i] + key_i;

	chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
	if(chan_c != chan_l){ continue; }
	ind_cl = Chan.ph1_map[Chan.ind0][Space.hash2(c, l, Chan.qnums2[Chan.ind0].j)];

	chanind_V = Ints.Vhhhp.V_3_3_index[chan_klc];
	chanind_T = Amps.D1.T1_index[chan_ab];
	V = Ints.Vhhhp.V_3_3[chanind_V + ind_klci];
	T1 = Amps.D1.T1[chanind_T + ind_abkj];
	T2 = Amps.S1.t2[ind_cl];
	if(permute == 1){
	  term0 = -1.0 * V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5h(1) += - < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * ";
	    std::cout <<"< '"<< a <<"','"<< b <<"' |t| "<< k <<",'"<< j <<"' > * < "<< c <<" |t| "<< l <<" > ";
	    std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = V * T1 * T2;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_5h(2) += < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * ";
	    std::cout <<"< '"<< a <<"','"<< b <<"' |t| "<< k <<",'"<< j <<"' > * < "<< c <<" |t| "<< l <<" > ";
	    std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_5h(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_5h(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  -V(kl|cd).T(ab|lj).t(d|i).t(c|k)  (D_7d)
double Diagram_D_7d(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, cd, ab, lj;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_kl, chan_cd, chan_lj, chan_d, chan_c, chan_k;
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_kl, key_cd, key_lj;
  int ind_klcd, ind_ablj, ind_di, ind_ck;
  double term = 0.0;
  double term0;
  int chanind_V, chanind_T;
  double V, T1, T2, T3;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      chan_d = Space.ind_1b(Parameters.basis, Space.qnums[d]);
      if(chan_d != chan_i){ continue; }
      ind_di = Chan.ph1_map[Chan.ind0][Space.hash2(d, i, Chan.qnums2[Chan.ind0].j)];

      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_c){ continue; }
	ind_ck = Chan.ph1_map[Chan.ind0][Space.hash2(c, k, Chan.qnums2[Chan.ind0].j)];
	for(int l = 0; l < Space.num_hol; ++l){
	  if(l == i){ continue; }
	  plus(lj, Space.qnums[l], Space.qnums[j]);
	  chan_lj = Space.ind_2b_dir(Parameters.basis, lj);
	  if(chan_lj != chan_ab){ continue; }
	  key_lj = Chan.hh_map[chan_lj][Space.hash2(l, j, lj.j)];
	  ind_ablj = key_ab * Chan.nhh[chan_lj] + key_lj;

	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_cd){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  chanind_T = Amps.D1.T1_index[chan_ab];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.D1.T1[chanind_T + ind_ablj];
	  T2 = Amps.S1.t2[ind_di];
	  T3 = Amps.S1.t2[ind_ck];
	  if(permute == 1){
	    term0 = -1.0 * V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7d(1) += - < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"','"<< b <<"' |t| "<< l <<",'"<< j <<"' > * < "<< d <<" |t| '"<< i <<"' > * < "<< c <<" |t| "<< k <<" > ";
	      std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = V * T1 * T2 * T3;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_7d(2) += < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"','"<< b <<"' |t| "<< l <<",'"<< j <<"' > * < "<< d <<" |t| '"<< i <<"' > * < "<< c <<" |t| "<< k <<" > ";
	      std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_7d(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_7d(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  V(ab|ic).t(c|j)  (D_4a)
double Diagram_D_4a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ab, abc;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_abc, chan_c;
  int key_i = Chan.h_map[chan_i][i];
  int key_abc;
  int ind_abci, ind_cj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    minus(abc, ab, Space.qnums[c]);
    chan_abc = Space.ind_1b(Parameters.basis, abc);
    if(chan_abc != chan_i){ continue; }
    key_abc = Chan.ppp_map[chan_abc][Space.hash3(a, b, c, ab.j)];
    ind_abci = key_abc * Chan.nh[chan_i] + key_i;
    
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_j){ continue; }
    ind_cj = Chan.ph1_map[Chan.ind0][Space.hash2(c, j, Chan.qnums2[Chan.ind0].j)];

    chanind_V = Ints.Vpphp.V_3_3_index[chan_abc];    
    V = Ints.Vpphp.V_3_3[chanind_V + ind_abci];
    T = Amps.S1.t2[ind_cj];
    if(permute == 1){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_4a(1) += < '"<< a <<"','"<< b <<"' |V| '"<< i <<"',"<< c <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"D_4a(2) += - < '"<< a <<"','"<< b <<"' |V| '"<< i <<"',"<< c <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_4a(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_4a(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  X'(ab|ic).t(c|j)  (X_4a)
double Diagram_X_4a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ab, abc;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_abc, chan_c;
  int key_i = Chan.h_map[chan_i][i];
  int key_abc;
  int ind_abci, ind_cj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    minus(abc, ab, Space.qnums[c]);
    chan_abc = Space.ind_1b(Parameters.basis, abc);
    if(chan_abc != chan_i){ continue; }
    key_abc = Chan.ppp_map[chan_abc][Space.hash3(a, b, c, ab.j)];
    ind_abci = key_abc * Chan.nh[chan_i] + key_i;
    
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_j){ continue; }
    ind_cj = Chan.ph1_map[Chan.ind0][Space.hash2(c, j, Chan.qnums2[Chan.ind0].j)];

    chanind_V = Eff_Ints.Xpphp.X_3_3_index[chan_abc];
    V = Eff_Ints.Xpphp.X1_3_3[chanind_V + ind_abci];
    T = Amps.S1.t2[ind_cj];
    if(permute == 1){
      term0 = V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_4a(1) += < '"<< a <<"','"<< b <<"' |X'| '"<< i <<"',"<< c <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	std::cout <<"= "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    else if(permute == 2){
      term0 = -1.0 * V * T;
      if(print == 2 && term0 != 0.0){
	std::cout <<"X_4a(2) += - < '"<< a <<"','"<< b <<"' |X'| '"<< i <<"',"<< c <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	std::cout <<"= - "<< V <<" * "<< T <<" = "<< term0 << std::endl;
      }
    }
    term += term0;
  }
  if(permute == 1){
    if(print != 0){ std::cout << "X_4a(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "X_4a(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  (1/2).V(ab|cd).t(c|i).t(d|j)  (D_6a)
double Diagram_D_6a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ab, cd;
  plus(ab, Space.qnums[a], Space.qnums[b]);
  int chan_ab = Space.ind_2b_dir(Parameters.basis, ab);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_cd, chan_c, chan_d;
  int key_ab = Chan.pp_map[chan_ab][Space.hash2(a, b, ab.j)];
  int key_cd;
  int ind_abcd, ind_ci, ind_dj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_i){ continue; }
    ind_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      chan_d = Space.ind_1b(Parameters.basis, Space.qnums[d]);
      if(chan_d != chan_j){ continue; }
      ind_dj = Chan.ph1_map[Chan.ind0][Space.hash2(d, j, Chan.qnums2[Chan.ind0].j)];

      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      if(chan_cd != chan_ab){ continue; }
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      ind_abcd = key_ab * Chan.npp[chan_cd] + key_cd;

      chanind_V = Ints.Vpppp.V_1_index[chan_ab];    
      V = Ints.Vpppp.V_1[chanind_V + ind_abcd];
      T1 = Amps.S1.t2[ind_ci];
      T2 = Amps.S1.t2[ind_dj];
      if(permute == 1){
	term0 = 0.5 * V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6a(1) += (1/2) * < '"<< a <<"','"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	  std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	  std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      else if(permute == 2){
	term0 = -0.5 * V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6a(2) += -(1/2) * < '"<< a <<"','"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	  std::cout <<"< "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	  std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      term += term0;
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_6a(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_6a(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  -V(kb|ic).t(a|k).t(c|j)  (D_6c)
double Diagram_D_6c(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kc, ib;
  minus(ib, Space.qnums[i], Space.qnums[b]);
  int chan_ib = Space.ind_2b_cross(Parameters.basis, ib);
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_kc, chan_k, chan_c;
  int key_ib = Chan.hp1_map[chan_ib][Space.hash2(i, b, ib.j)];
  int key_kc;
  int ind_kcib, ind_ak, ind_cj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_j){ continue; }
    ind_cj = Chan.ph1_map[Chan.ind0][Space.hash2(c, j, Chan.qnums2[Chan.ind0].j)];
    for(int k = 0; k < Space.num_hol; ++k){
      chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
      if(chan_k != chan_a){ continue; }
      ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];

      minus(kc, Space.qnums[k], Space.qnums[c]);
      chan_kc = Space.ind_2b_cross(Parameters.basis, kc);
      if(chan_kc != chan_ib){ continue; }
      key_kc = Chan.hp1_map[chan_kc][Space.hash2(k, c, kc.j)];
      ind_kcib = key_kc * Chan.nhp1[chan_ib] + key_ib;

      chanind_V = Ints.Vhphp.V_2_1_index[chan_kc];    
      V = Ints.Vhphp.V_2_1[chanind_V + ind_kcib];
      T1 = Amps.S1.t2[ind_ak];
      T2 = Amps.S1.t2[ind_cj];
      if(permute == 1){
	term0 = -1.0 * V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6c(1) += - < "<< k <<",'"<< b <<"' |V| '"<< i <<"',"<< c <<" > * ";
	  std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	  std::cout <<"= - "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      else if(permute == 2){
	term0 = V * T1 * T2;
	if(print == 2 && term0 != 0.0){
	  std::cout <<"D_6c(2) += < "<< k <<",'"<< b <<"' |V| '"<< i <<"',"<< c <<" > * ";
	  std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	  std::cout <<"= "<< V <<" * "<< T1 <<" * "<< T2 <<" = "<< term0 << std::endl;
	}
      }
      term += term0;
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_6c(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_6c(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  -(1/2).V(kb|cd).t(a|k).t(c|i).t(d|j)  (D_8a)
double Diagram_D_8a(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State bd, ck;
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_bd, chan_ck, chan_k, chan_c, chan_d;
  int key_bd, key_ck;
  int ind_bdck, ind_ak, ind_ci, ind_dj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2, T3;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_i){ continue; }
    ind_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      chan_d = Space.ind_1b(Parameters.basis, Space.qnums[d]);
      if(chan_d != chan_j){ continue; }
      ind_dj = Chan.ph1_map[Chan.ind0][Space.hash2(d, j, Chan.qnums2[Chan.ind0].j)];

      minus(bd, Space.qnums[b], Space.qnums[d]);
      chan_bd = Space.ind_2b_cross(Parameters.basis, bd);
      key_bd = Chan.pp1_map[chan_bd][Space.hash2(b, d, bd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_a){ continue; }
	ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];

	minus(ck, Space.qnums[c], Space.qnums[k]);
	chan_ck = Space.ind_2b_cross(Parameters.basis, ck);
	if(chan_ck != chan_bd){ continue; }
	key_ck = Chan.ph1_map[chan_ck][Space.hash2(c, k, ck.j)];
	ind_bdck = key_bd * Chan.nph1[chan_ck] + key_ck;

	chanind_V = Ints.Vhppp.V_2_4_index[chan_bd];
	V = Ints.Vhppp.V_2_4[chanind_V + ind_bdck];
	T1 = Amps.S1.t2[ind_ak];
	T2 = Amps.S1.t2[ind_ci];
	T3 = Amps.S1.t2[ind_dj];
	if(permute == 1){
	  term0 = -0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_8a(1) += -(1/2) * < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	    std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = 0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_8a(2) += (1/2) * < "<< k <<",'"<< b <<"' |V| "<< c <<","<< d <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_8a(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_8a(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  (1/2).V(kl|ic).t(a|k).t(b|l).t(c|j)  (D_8b)
double Diagram_D_8b(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State ki, cl;
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_ki, chan_cl, chan_k, chan_l, chan_c;
  int key_ki, key_cl;
  int ind_kicl, ind_ak, ind_bl, ind_cj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2, T3;
  for(int l = 0; l < Space.num_hol; ++l){
    chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
    if(chan_l != chan_b){ continue; }
    ind_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];
    for(int k = 0; k < Space.num_hol; ++k){
      chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
      if(chan_k != chan_a){ continue; }
      ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];

      minus(ki, Space.qnums[k], Space.qnums[i]);
      chan_ki = Space.ind_2b_cross(Parameters.basis, ki);
      key_ki = Chan.hh1_map[chan_ki][Space.hash2(k, i, ki.j)];
      for(int c = Space.num_hol; c < Space.num_states; ++c){
	chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
	if(chan_c != chan_j){ continue; }
	ind_cj = Chan.ph1_map[Chan.ind0][Space.hash2(c, j, Chan.qnums2[Chan.ind0].j)];

	minus(cl, Space.qnums[c], Space.qnums[l]);
	chan_cl = Space.ind_2b_cross(Parameters.basis, cl);
	if(chan_cl != chan_ki){ continue; }
	key_cl = Chan.ph1_map[chan_cl][Space.hash2(c, l, cl.j)];
	ind_kicl = key_ki * Chan.nph1[chan_cl] + key_cl;

	chanind_V = Ints.Vhhhp.V_2_3_index[chan_ki];
	V = Ints.Vhhhp.V_2_3[chanind_V + ind_kicl];
	T1 = Amps.S1.t2[ind_ak];
	T2 = Amps.S1.t2[ind_bl];
	T3 = Amps.S1.t2[ind_cj];
	if(permute == 1){
	  term0 = 0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_8b(1) += (1/2) * < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	    std::cout <<"= (1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	}
	else if(permute == 2){
	  term0 = -0.5 * V * T1 * T2 * T3;
	  if(print == 2 && term0 != 0.0){
	    std::cout <<"D_8b(2) += -(1/2) * < "<< k <<","<< l <<" |V| '"<< i <<"',"<< c <<" > * ";
	    std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > * < "<< c <<" |t| '"<< j <<"' > ";
	    std::cout <<"= -(1/2) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" = "<< term0 << std::endl;
	  }
	}
	term += term0;
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_8b(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_8b(2) = " << term << std::endl; }
  }
  return term;
}

// T3_3(ab|ij)  <-  (1/4).V(kl|cd).t(a|k).t(b|l).t(c|i).t(d|j)  (D_9)
double Diagram_D_9(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Amplitudes &Amps, int &a, int &b, int &i, int &j, int &print, int permute)
{
  State kl, cd;
  int chan_a = Space.ind_1b(Parameters.basis, Space.qnums[a]);
  int chan_b = Space.ind_1b(Parameters.basis, Space.qnums[b]);
  int chan_i = Space.ind_1b(Parameters.basis, Space.qnums[i]);
  int chan_j = Space.ind_1b(Parameters.basis, Space.qnums[j]);
  int chan_kl, chan_cd, chan_k, chan_l, chan_c, chan_d;
  int key_kl, key_cd;
  int ind_klcd, ind_ak, ind_bl, ind_ci, ind_dj;
  double term = 0.0;
  double term0;
  int chanind_V;
  double V, T1, T2, T3, T4;
  for(int c = Space.num_hol; c < Space.num_states; ++c){
    chan_c = Space.ind_1b(Parameters.basis, Space.qnums[c]);
    if(chan_c != chan_i){ continue; }
    ind_ci = Chan.ph1_map[Chan.ind0][Space.hash2(c, i, Chan.qnums2[Chan.ind0].j)];
    for(int d = Space.num_hol; d < Space.num_states; ++d){
      chan_d = Space.ind_1b(Parameters.basis, Space.qnums[d]);
      if(chan_d != chan_j){ continue; }
      ind_dj = Chan.ph1_map[Chan.ind0][Space.hash2(d, j, Chan.qnums2[Chan.ind0].j)];

      plus(cd, Space.qnums[c], Space.qnums[d]);
      chan_cd = Space.ind_2b_dir(Parameters.basis, cd);
      key_cd = Chan.pp_map[chan_cd][Space.hash2(c, d, cd.j)];
      for(int k = 0; k < Space.num_hol; ++k){
	chan_k = Space.ind_1b(Parameters.basis, Space.qnums[k]);
	if(chan_k != chan_a){ continue; }
	ind_ak = Chan.ph1_map[Chan.ind0][Space.hash2(a, k, Chan.qnums2[Chan.ind0].j)];
	for(int l = 0; l < Space.num_hol; ++l){
	  chan_l = Space.ind_1b(Parameters.basis, Space.qnums[l]);
	  if(chan_l != chan_b){ continue; }
	  ind_bl = Chan.ph1_map[Chan.ind0][Space.hash2(b, l, Chan.qnums2[Chan.ind0].j)];

	  plus(kl, Space.qnums[k], Space.qnums[l]);
	  chan_kl = Space.ind_2b_dir(Parameters.basis, kl);
	  if(chan_kl != chan_cd){ continue; }
	  key_kl = Chan.hh_map[chan_kl][Space.hash2(k, l, kl.j)];
	  ind_klcd = key_kl * Chan.npp[chan_cd] + key_cd;

	  chanind_V = Ints.Vhhpp.V_1_index[chan_kl];
	  V = Ints.Vhhpp.V_1[chanind_V + ind_klcd];
	  T1 = Amps.S1.t2[ind_ak];
	  T2 = Amps.S1.t2[ind_bl];
	  T3 = Amps.S1.t2[ind_ci];
	  T4 = Amps.S1.t2[ind_dj];
	  if(permute == 1){
	    term0 = 0.25 * V * T1 * T2 * T3 * T4;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_9(1) += (1/4) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > * < "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	      std::cout <<"= (1/4) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" * "<< T4 <<" = "<< term0 << std::endl;
	    }
	  }
	  else if(permute == 2){
	    term0 = -0.25 * V * T1 * T2 * T3 * T4;
	    if(print == 2 && term0 != 0.0){
	      std::cout <<"D_9(2) += -(1/4) * < "<< k <<","<< l <<" |V| "<< c <<","<< d <<" > * ";
	      std::cout <<"< '"<< a <<"' |t| "<< k <<" > * < '"<< b <<"' |t| "<< l <<" > * < "<< c <<" |t| '"<< i <<"' > * < "<< d <<" |t| '"<< j <<"' > ";
	      std::cout <<"= -(1/4) * "<< V <<" * "<< T1 <<" * "<< T2 <<" * "<< T3 <<" * "<< T4 <<" = "<< term0 << std::endl;
	    }
	  }
	  term += term0;
	}
      }
    }
  }
  if(permute == 1){
    if(print != 0){ std::cout << "D_9(1) = " << term << std::endl; }
  }
  else if(permute == 2){
    if(print != 0){ std::cout << "D_9(2) = " << term << std::endl; }
  }
  return term;
}


void CC_Test_Full_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  CC_Test_T1_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T2_1_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T2_2_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T2_3_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T2_4_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T3_1_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T3_2_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T3_3_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_T3_4_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_t2_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
  CC_Test_t3_J(ParametersJ, SpaceJ, ChanJ, AmpsJ, Parameters, Space, Chan, Amps);
}

void CC_Test_T1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind1, chan1;
  int nhh, npp;
  int a, b, i, j;
  two_body pp, hh;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    chanind = Amps.D1.T1_index[chan1];
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      pp = Chan.pp_state(chan1, pp_ind);
      a = pp.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      b = pp.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      if(a == b){ continue; }
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	i = hh.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	j = hh.v2;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	if(i == j){ continue; }
	ind1 = pp_ind * nhh + hh_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_j) > minJ ){ minJ = abs(j_i - j_j); }
	maxJ = j_a + j_b;
	if( j_i + j_j < maxJ ){ maxJ = j_i + j_j; }
	M = m_a + m_b;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  plus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_b]);
	  plus(stateJ2, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_j]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  if( shell_a == shell_b && J%4 != 0 ){ continue; }
	  if( shell_i == shell_j && J%4 != 0 ){ continue; }
	  chanJ1 = SpaceJ.ind_2b_dir(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_dir(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T1_index[chanJ1];
	  keyJ1 = ChanJ.pp_map[chanJ1][SpaceJ.hash2(shell_a, shell_b, J)];
	  keyJ2 = ChanJ.hh_map[chanJ2][SpaceJ.hash2(shell_i, shell_j, J)];
	  tJ += CGC(j_a, m_a, j_b, m_b, J, M) * CGC(j_i, m_i, j_j, m_j, J, M) * AmpsJ.D1.T1[chanind_J + (keyJ1*ChanJ.nhh[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T1[chanind + ind1])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T1 < " << a << "," << b << " |t| " << i << "," << j << " > !! -- " << Amps.D1.T1[chanind + ind1] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind2, chan2;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind = Amps.D1.T2_1_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      a = ph1.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      j = ph1.v2;
      state_j = Space.qnums[j];
      state_j.j = Space.shellsj[j];
      key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
      shell_j = SpaceJ.map_state[key_j];
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	i = hp1.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	b = hp1.v2;
	state_b = Space.qnums[b];
	state_b.j = Space.shellsj[b];
	key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
	shell_b = SpaceJ.map_state[key_b];
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = -1 * Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = -1 * Space.qnums[j].m;

	minJ = abs(j_a - j_j);
	if( abs(j_i - j_b) > minJ ){ minJ = abs(j_i - j_b); }
	maxJ = j_a + j_j;
	if( j_i + j_b < maxJ ){ maxJ = j_i + j_b; }
	M = m_a + m_j;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  minus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_j]);
	  minus(stateJ2, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_b]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  chanJ1 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T2_1_index[chanJ1];
	  keyJ1 = ChanJ.ph1_map[chanJ1][SpaceJ.hash2(shell_a, shell_j, J)];
	  keyJ2 = ChanJ.hp1_map[chanJ2][SpaceJ.hash2(shell_i, shell_b, J)];
	  tJ += -1.0 * phase2(j_b + m_b + j_j + m_j) * CGC(j_a, m_a, j_j, m_j, J, M) * CGC(j_i, m_i, j_b, m_b, J, M) * AmpsJ.D1.T2_1[chanind_J + (keyJ1*ChanJ.nhp1[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T2_1[chanind + ind2])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T2_1 < " << a << "," << b << " |t| " << i << "," << j << " > -- < " << shell_a << "," << shell_b << " |t| " << shell_i << "," << shell_j << " > (" << minJ << "," << maxJ << ")  :  " << Amps.D1.T2_1[chanind + ind2] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind2, chan2;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind = Amps.D1.T2_2_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      b = ph1.v1;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      i = ph1.v2;
      state_i = Space.qnums[i];
      state_i.j = Space.shellsj[i];
      key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
      shell_i = SpaceJ.map_state[key_i];
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	j = hp1.v1;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	a = hp1.v2;
	state_a = Space.qnums[a];
	state_a.j = Space.shellsj[a];
	key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
	shell_a = SpaceJ.map_state[key_a];
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	j_a = Space.shellsj[a];
	m_a = -1 * Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = -1 * Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_b - j_i);
	if( abs(j_j - j_a) > minJ ){ minJ = abs(j_j - j_a); }
	maxJ = j_b + j_i;
	if( j_j + j_a < maxJ ){ maxJ = j_j + j_a; }
	M = m_b + m_i;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  minus(stateJ1, SpaceJ.qnums[shell_b], SpaceJ.qnums[shell_i]);
	  minus(stateJ2, SpaceJ.qnums[shell_j], SpaceJ.qnums[shell_a]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  chanJ1 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T2_2_index[chanJ1];
	  keyJ1 = ChanJ.ph1_map[chanJ1][SpaceJ.hash2(shell_b, shell_i, J)];
	  keyJ2 = ChanJ.hp1_map[chanJ2][SpaceJ.hash2(shell_j, shell_a, J)];
	  tJ += -1.0 * phase2(j_a + m_a + j_i + m_i) * CGC(j_b, m_b, j_i, m_i, J, M) * CGC(j_j, m_j, j_a, m_a, J, M) * AmpsJ.D1.T2_2[chanind_J + (keyJ1*ChanJ.nhp1[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T2_2[chanind + ind2])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T2_2 < " << a << "," << b << " |t| " << i << "," << j << " > -- < " << shell_a << "," << shell_b << " |t| " << shell_i << "," << shell_j << " > (" << minJ << "," << maxJ << ")  :  " << Amps.D1.T2_2[chanind + ind2] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind2, chan2;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind = Amps.D1.T2_3_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      a = ph1.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      i = ph1.v2;
      state_i = Space.qnums[i];
      state_i.j = Space.shellsj[i];
      key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
      shell_i = SpaceJ.map_state[key_i];
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	j = hp1.v1;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	b = hp1.v2;
	state_b = Space.qnums[b];
	state_b.j = Space.shellsj[b];
	key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
	shell_b = SpaceJ.map_state[key_b];
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = -1 * Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = -1 * Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_a - j_i);
	if( abs(j_j - j_b) > minJ ){ minJ = abs(j_j - j_b); }
	maxJ = j_a + j_i;
	if( j_j + j_b < maxJ ){ maxJ = j_j + j_b; }
	M = m_a + m_i;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  minus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_i]);
	  minus(stateJ2, SpaceJ.qnums[shell_j], SpaceJ.qnums[shell_b]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  chanJ1 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T2_3_index[chanJ1];
	  keyJ1 = ChanJ.ph1_map[chanJ1][SpaceJ.hash2(shell_a, shell_i, J)];
	  keyJ2 = ChanJ.hp1_map[chanJ2][SpaceJ.hash2(shell_j, shell_b, J)];
	  tJ += phase2(j_b + m_b + j_i + m_i) * CGC(j_a, m_a, j_i, m_i, J, M) * CGC(j_j, m_j, j_b, m_b, J, M) * AmpsJ.D1.T2_3[chanind_J + (keyJ1*ChanJ.nhp1[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T2_3[chanind + ind2])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T2_3 < " << a << "," << b << " |t| " << i << "," << j << " > -- < " << shell_a << "," << shell_b << " |t| " << shell_i << "," << shell_j << " > (" << minJ << "," << maxJ << ")  :  " << Amps.D1.T2_3[chanind + ind2] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T2_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind2, chan2;
  int nph1, nhp1;
  int a, b, i, j;
  two_body ph1, hp1;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    chanind = Amps.D1.T2_4_index[chan2];
    for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
      ph1 = Chan.ph1_state(chan2, ph1_ind);
      b = ph1.v1;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      j = ph1.v2;
      state_j = Space.qnums[j];
      state_j.j = Space.shellsj[j];
      key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
      shell_j = SpaceJ.map_state[key_j];
      for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
	hp1 = Chan.hp1_state(chan2, hp1_ind);
	i = hp1.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	a = hp1.v2;
	state_a = Space.qnums[a];
	state_a.j = Space.shellsj[a];
	key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
	shell_a = SpaceJ.map_state[key_a];
	if(a == b || i == j){ continue; }
	ind2 = ph1_ind * nhp1 + hp1_ind;

	j_a = Space.shellsj[a];
	m_a = -1 * Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = -1 * Space.qnums[j].m;

	minJ = abs(j_b - j_j);
	if( abs(j_i - j_a) > minJ ){ minJ = abs(j_i - j_a); }
	maxJ = j_b + j_j;
	if( j_i + j_a < maxJ ){ maxJ = j_i + j_a; }
	M = m_b + m_j;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  minus(stateJ1, SpaceJ.qnums[shell_b], SpaceJ.qnums[shell_j]);
	  minus(stateJ2, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_a]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  chanJ1 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T2_4_index[chanJ1];
	  keyJ1 = ChanJ.ph1_map[chanJ1][SpaceJ.hash2(shell_b, shell_j, J)];
	  keyJ2 = ChanJ.hp1_map[chanJ2][SpaceJ.hash2(shell_i, shell_a, J)];
	  tJ += phase2(j_a + m_a + j_j + m_j) * CGC(j_b, m_b, j_j, m_j, J, M) * CGC(j_i, m_i, j_a, m_a, J, M) * AmpsJ.D1.T2_4[chanind_J + (keyJ1*ChanJ.nhp1[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T2_4[chanind + ind2])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T2_4 < " << a << "," << b << " |t| " << i << "," << j << " > -- < " << shell_a << "," << shell_b << " |t| " << shell_i << "," << shell_j << " > (" << minJ << "," << maxJ << ")  :  " << Amps.D1.T2_4[chanind + ind2] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind3, chan3;
  int np, nhhp;
  int a, b, i, j;
  one_body p;
  three_body hhp;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    chanind = Amps.D1.T3_1_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
	hhp = Chan.hhp_state(chan3, hhp_ind);
	i = hhp.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	j = hhp.v2;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	b = hhp.v3;
	state_b = Space.qnums[b];
	state_b.j = Space.shellsj[b];
	key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
	shell_b = SpaceJ.map_state[key_b];
	if(a == b || i == j){ continue; }
	ind3 = p_ind * nhhp + hhp_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = -1 * Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_j) > minJ ){ minJ = abs(j_i - j_j); }
	maxJ = j_a + j_b;
	if( j_i + j_j < maxJ ){ maxJ = j_i + j_j; }
	M = m_i + m_j;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_a]);
	  plus(stateJ1, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_j]);
	  stateJ1.j = J;
	  if(shell_i == shell_j && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_b]);
	  stateJ2.j = SpaceJ.qnums[shell_a].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T3_1_index[chanJ1];
	  keyJ1 = ChanJ.p_map[chanJ1][shell_a];
	  keyJ2 = ChanJ.hhp_map[chanJ2][SpaceJ.hash3(shell_i, shell_j, shell_b, J)];
	  tJ += phase2(j_b + m_b) * CGC(J, M, j_b, m_b, j_a, m_a) * CGC(j_i, m_i, j_j, m_j, J, M) * AmpsJ.D1.T3_1[chanind_J + (keyJ1*ChanJ.nhhp[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T3_1[chanind + ind3])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T3_1 < " << a << "," << b << " |t| " << i << "," << j << " > !! -- " << Amps.D1.T3_1[chanind + ind3] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind3, chan3;
  int np, nhhp;
  int a, b, i, j;
  one_body p;
  three_body hhp;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    chanind = Amps.D1.T3_2_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      b = p.v1;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
	hhp = Chan.hhp_state(chan3, hhp_ind);
	i = hhp.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	j = hhp.v2;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	a = hhp.v3;
	state_a = Space.qnums[a];
	state_a.j = Space.shellsj[a];
	key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
	shell_a = SpaceJ.map_state[key_a];
	if(a == b || i == j){ continue; }
	ind3 = p_ind * nhhp + hhp_ind;

	j_a = Space.shellsj[a];
	m_a = -1 * Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_j) > minJ ){ minJ = abs(j_i - j_j); }
	maxJ = j_a + j_b;
	if( j_i + j_j < maxJ ){ maxJ = j_i + j_j; }
	M = m_i + m_j;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_b]);
	  plus(stateJ1, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_j]);
	  stateJ1.j = J;
	  if(shell_i == shell_j && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_a]);
	  stateJ2.j = SpaceJ.qnums[shell_b].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T3_2_index[chanJ1];
	  keyJ1 = ChanJ.p_map[chanJ1][shell_b];
	  keyJ2 = ChanJ.hhp_map[chanJ2][SpaceJ.hash3(shell_i, shell_j, shell_a, J)];
	  tJ += -1.0 * phase2(j_a + m_a) * CGC(J, M, j_a, m_a, j_b, m_b) * CGC(j_i, m_i, j_j, m_j, J, M) * AmpsJ.D1.T3_2[chanind_J + (keyJ1*ChanJ.nhhp[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T3_2[chanind + ind3])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T3_2 < " << a << "," << b << " |t| " << i << "," << j << " > !! -- " << Amps.D1.T3_2[chanind + ind3] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind3, chan3;
  int npph, nh;
  int a, b, i, j;
  three_body pph;
  one_body h;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    npph = Chan.npph[chan3];
    nh = Chan.nh[chan3];
    chanind = Amps.D1.T3_3_index[chan3];
    for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
      pph = Chan.pph_state(chan3, pph_ind);
      a = pph.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      b = pph.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      j = pph.v3;
      state_j = Space.qnums[j];
      state_j.j = Space.shellsj[j];
      key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
      shell_j = SpaceJ.map_state[key_j];
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	i = h.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	if(a == b || i == j){ continue; }
	ind3 = pph_ind * nh + h_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = -1 * Space.qnums[j].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_j) > minJ ){ minJ = abs(j_i - j_j); }
	maxJ = j_a + j_b;
	if( j_i + j_j < maxJ ){ maxJ = j_i + j_j; }
	M = m_a + m_b;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_i]);
	  plus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_b]);
	  stateJ1.j = J;
	  if(shell_a == shell_b && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_j]);
	  stateJ2.j = SpaceJ.qnums[shell_i].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T3_3_index[chanJ1];
	  keyJ1 = ChanJ.pph_map[chanJ1][SpaceJ.hash3(shell_a, shell_b, shell_j, J)];
	  keyJ2 = ChanJ.h_map[chanJ2][shell_i];
	  tJ += phase2(j_j + m_j) * CGC(J, M, j_j, m_j, j_i, m_i) * CGC(j_a, m_a, j_b, m_b, J, M) * AmpsJ.D1.T3_3[chanind_J + (keyJ1*ChanJ.nh[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T3_3[chanind + ind3])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T3_3 < " << a << "," << b << " |t| " << i << "," << j << " > !! -- " << Amps.D1.T3_3[chanind + ind3] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_T3_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind3, chan3;
  int npph, nh;
  int a, b, i, j;
  three_body pph;
  one_body h;
  double tJ;
  int key_a, key_b, key_i, key_j;
  int shell_a, shell_b, shell_i, shell_j;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_j;
  int j_a, j_b, j_i, j_j;
  int m_a, m_b, m_i, m_j, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    npph = Chan.npph[chan3];
    nh = Chan.nh[chan3];
    chanind = Amps.D1.T3_4_index[chan3];
    for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
      pph = Chan.pph_state(chan3, pph_ind);
      a = pph.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      b = pph.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      i = pph.v3;
      state_i = Space.qnums[i];
      state_i.j = Space.shellsj[i];
      key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
      shell_i = SpaceJ.map_state[key_i];
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	j = h.v1;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	if(a == b || i == j){ continue; }
	ind3 = pph_ind * nh + h_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = -1 * Space.qnums[i].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_j) > minJ ){ minJ = abs(j_i - j_j); }
	maxJ = j_a + j_b;
	if( j_i + j_j < maxJ ){ maxJ = j_i + j_j; }
	M = m_a + m_b;

	tJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_j]);
	  plus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_b]);
	  stateJ1.j = J;
	  if(shell_a == shell_b && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_i]);
	  stateJ2.j = SpaceJ.qnums[shell_j].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = AmpsJ.D1.T3_4_index[chanJ1];
	  keyJ1 = ChanJ.pph_map[chanJ1][SpaceJ.hash3(shell_a, shell_b, shell_i, J)];
	  keyJ2 = ChanJ.h_map[chanJ2][shell_j];
	  tJ += -1.0 * phase2(j_i + m_i) * CGC(J, M, j_i, m_i, j_j, m_j) * CGC(j_a, m_a, j_b, m_b, J, M) * AmpsJ.D1.T3_4[chanind_J + (keyJ1*ChanJ.nh[chanJ1] + keyJ2)];
	}
	if( std::fabs((tJ - Amps.D1.T3_4[chanind + ind3])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! T3_4 < " << a << "," << b << " |t| " << i << "," << j << " > !! -- " << Amps.D1.T3_4[chanind + ind3] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_t2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind2, chan2;
  int nph1;
  int a, i;
  two_body ph1;
  double tJ;
  int key_a, key_i;
  int shell_a, shell_i;
  int chanJ1, keyJ1;
  State state_a, state_i;
  int j_a, j_i, J;
  int m_a, m_i, M;
  chan2 = Chan.ind0;
  nph1 = Chan.nph1[chan2];
  for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
    ph1 = Chan.ph1_state(chan2, ph1_ind);
    a = ph1.v1;
    state_a = Space.qnums[a];
    state_a.j = Space.shellsj[a];
    key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
    shell_a = SpaceJ.map_state[key_a];
    i = ph1.v2;
    state_i = Space.qnums[i];
    state_i.j = Space.shellsj[i];
    key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
    shell_i = SpaceJ.map_state[key_i];
    ind2 = ph1_ind;

    j_a = Space.shellsj[a];
    m_a = Space.qnums[a].m;
    j_i = Space.shellsj[i];
    m_i = -1 * Space.qnums[i].m;
    if( j_a != j_i || m_a != m_i ){ continue; }
    J = 0;
    M = 0;

    chanJ1 = ChanJ.ind0;
    keyJ1 = ChanJ.ph1_map[chanJ1][SpaceJ.hash2(shell_a, shell_i, J)];
    tJ = phase2(j_i + m_i) * CGC(j_a, m_a, j_i, m_i, J, M) * AmpsJ.S1.t2[keyJ1];
    if( std::fabs((tJ - Amps.S1.t2[ph1_ind])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
      std::cout << "!! t2 < " << a << " |t| " << i << " > !! -- " << Amps.S1.t2[ph1_ind] << " " << tJ << std::endl;
    }
  }
}

void CC_Test_t3_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Amplitudes &AmpsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Amplitudes &Amps)
{
  int ind3, chan3;
  int nh, np;
  int a, i;
  one_body p, h;
  double tJ;
  int key_a, key_i;
  int shell_a, shell_i;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_i;
  int j_a, j_i;
  int m_a, m_i;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nh = Chan.nh[chan3];
    chanind = Amps.S1.t3_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	i = h.v1;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	ind3 = p_ind * nh + h_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	if( j_a != j_i || m_a != m_i ){ continue; }

	chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_a]);
	chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_i]);
	if( chanJ1 != chanJ2 ){ continue; }
	chanind_J = AmpsJ.S1.t3_index[chanJ1];
	keyJ1 = ChanJ.p_map[chanJ1][shell_a];
	keyJ2 = ChanJ.h_map[chanJ2][shell_i];
	tJ = AmpsJ.S1.t3[chanind_J + (keyJ1*ChanJ.nh[chanJ1] + keyJ2)];
	if( std::fabs((tJ - Amps.S1.t3[chanind + ind3])/(std::fabs(tJ) + 1e-12)) > 1e-8 && std::fabs(tJ) > 1e-12 ){
	  std::cout << "!! t3 < " << a << " |t| " << i << " > !! -- " << Amps.S1.t3[chanind + ind3] << " " << tJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_Xpp2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind2, chan2;
  int npp1;
  int a, b;
  two_body pp1;
  double xJ;
  int key_a, key_b;
  int shell_a, shell_b;
  int chanJ1, keyJ1;
  State state_a, state_b;
  int j_a, j_b, J;
  int m_a, m_b, M;
  chan2 = Chan.ind0;
  npp1 = Chan.npp1[chan2];
  for(int pp1_ind = 0; pp1_ind < npp1; ++pp1_ind){
    pp1 = Chan.pp1_state(chan2, pp1_ind);
    a = pp1.v1;
    state_a = Space.qnums[a];
    state_a.j = Space.shellsj[a];
    key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
    shell_a = SpaceJ.map_state[key_a];
    b = pp1.v2;
    state_b = Space.qnums[b];
    state_b.j = Space.shellsj[b];
    key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
    shell_b = SpaceJ.map_state[key_b];
    ind2 = pp1_ind;

    j_a = Space.shellsj[a];
    m_a = Space.qnums[a].m;
    j_b = Space.shellsj[b];
    m_b = -1 * Space.qnums[b].m;
    if( j_a != j_b || m_a != m_b ){ continue; }
    J = 0;
    M = 0;

    chanJ1 = ChanJ.ind0;
    keyJ1 = ChanJ.pp1_map[chanJ1][SpaceJ.hash2(shell_a, shell_b, J)];
    xJ = phase2(j_b + m_b) * CGC(j_a, m_a, j_b, m_b, J, M) * Eff_IntsJ.Xpp.X_2[keyJ1];
    if( std::fabs((xJ - Eff_Ints.Xpp.X_2[pp1_ind])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
      std::cout << "!! xpp2 < " << a << " |X| " << b << " > !! -- " << Eff_Ints.Xpp.X_2[pp1_ind] << " " << xJ << std::endl;
    }
  }
}

void CC_Test_Xhp2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind2, chan2;
  int nhp1;
  int i, a;
  two_body hp1;
  double xJ;
  int key_i, key_a;
  int shell_i, shell_a;
  int chanJ1, keyJ1;
  State state_i, state_a;
  int j_i, j_a, J;
  int m_i, m_a, M;
  chan2 = Chan.ind0;
  nhp1 = Chan.nhp1[chan2];
  for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
    hp1 = Chan.hp1_state(chan2, hp1_ind);
    i = hp1.v1;
    state_i = Space.qnums[i];
    state_i.j = Space.shellsj[i];
    key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
    shell_i = SpaceJ.map_state[key_i];
    a = hp1.v2;
    state_a = Space.qnums[a];
    state_a.j = Space.shellsj[a];
    key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
    shell_a = SpaceJ.map_state[key_a];
    ind2 = hp1_ind;

    j_i = Space.shellsj[i];
    m_i = Space.qnums[i].m;
    j_a = Space.shellsj[a];
    m_a = -1 * Space.qnums[a].m;
    if( j_i != j_a || m_i != m_a ){ continue; }
    J = 0;
    M = 0;

    chanJ1 = ChanJ.ind0;
    keyJ1 = ChanJ.hp1_map[chanJ1][SpaceJ.hash2(shell_i, shell_a, J)];
    xJ = phase2(j_a + m_a) * CGC(j_i, m_i, j_a, m_a, J, M) * Eff_IntsJ.Xpp.X_2[keyJ1];
    if( std::fabs((xJ - Eff_Ints.Xpp.X_2[hp1_ind])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
      std::cout << "!! Xhp2 < " << i << " |X| " << a << " > !! -- " << Eff_Ints.Xpp.X_2[hp1_ind] << " " << xJ << std::endl;
    }
  }
}

void CC_Test_Xhh2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind2, chan2;
  int nhh1;
  int i, j;
  two_body hh1;
  double xJ;
  int key_i, key_j;
  int shell_i, shell_j;
  int chanJ1, keyJ1;
  State state_i, state_j;
  int j_i, j_j, J;
  int m_i, m_j, M;
  chan2 = Chan.ind0;
  nhh1 = Chan.nhh1[chan2];
  for(int hh1_ind = 0; hh1_ind < nhh1; ++hh1_ind){
    hh1 = Chan.hh1_state(chan2, hh1_ind);
    i = hh1.v1;
    state_i = Space.qnums[i];
    state_i.j = Space.shellsj[i];
    key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
    shell_i = SpaceJ.map_state[key_i];
    j = hh1.v2;
    state_j = Space.qnums[j];
    state_j.j = Space.shellsj[j];
    key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
    shell_j = SpaceJ.map_state[key_j];
    ind2 = hh1_ind;

    j_i = Space.shellsj[i];
    m_i = Space.qnums[i].m;
    j_j = Space.shellsj[j];
    m_j = -1 * Space.qnums[j].m;
    if( j_i != j_j || m_i != m_j ){ continue; }
    J = 0;
    M = 0;

    chanJ1 = ChanJ.ind0;
    keyJ1 = ChanJ.hh1_map[chanJ1][SpaceJ.hash2(shell_i, shell_j, J)];
    xJ = phase2(j_j + m_j) * CGC(j_i, m_i, j_j, m_j, J, M) * Eff_IntsJ.Xhh.X_2[keyJ1];
    if( std::fabs((xJ - Eff_Ints.Xhh.X_2[hh1_ind])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
      std::cout << "!! Xhh2 < " << i << " |X| " << j << " > !! -- " << Eff_Ints.Xhh.X_2[hh1_ind] << " " << xJ << std::endl;
    }
  }
}

void CC_Test_Xpphp3_4_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind3, chan3;
  int npph, np;
  int a, b, i, c;
  three_body pph;
  one_body p;
  double xJ;
  int key_a, key_b, key_i, key_c;
  int shell_a, shell_b, shell_i, shell_c;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_i, state_c;
  int j_a, j_b, j_i, j_c;
  int m_a, m_b, m_i, m_c, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    npph = Chan.npph[chan3];
    np = Chan.np[chan3];
    chanind = Eff_Ints.Xpphp.X_3_4_index[chan3];
    for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
      pph = Chan.pph_state(chan3, pph_ind);
      a = pph.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      b = pph.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      i = pph.v3;
      state_i = Space.qnums[i];
      state_i.j = Space.shellsj[i];
      key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
      shell_i = SpaceJ.map_state[key_i];
      for(int p_ind = 0; p_ind < np; ++p_ind){
	p = Chan.p_state(chan3, p_ind);
	c = p.v1;
	state_c = Space.qnums[c];
	state_c.j = Space.shellsj[c];
	key_c = SpaceJ.ind_state(ParametersJ.basis, state_c);
	shell_c = SpaceJ.map_state[key_c];
	if( a == b ){ continue; }
	ind3 = pph_ind * np + p_ind;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_i = Space.shellsj[i];
	m_i = -1 * Space.qnums[i].m;
	j_c = Space.shellsj[c];
	m_c = Space.qnums[c].m;

	minJ = abs(j_a - j_b);
	if( abs(j_i - j_c) > minJ ){ minJ = abs(j_i - j_c); }
	maxJ = j_a + j_b;
	if( j_i + j_c < maxJ ){ maxJ = j_i + j_c; }
	M = m_a + m_b;

	xJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_c]);
	  plus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_b]);
	  stateJ1.j = J;
	  if(shell_a == shell_b && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_i]);
	  stateJ2.j = SpaceJ.qnums[shell_c].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = Eff_IntsJ.Xpphp.X_3_4_index[chanJ1];
	  keyJ1 = ChanJ.pph_map[chanJ1][SpaceJ.hash3(shell_a, shell_b, shell_i, J)];
	  keyJ2 = ChanJ.p_map[chanJ2][shell_c];
	  xJ += -1.0 * phase2(j_i + m_i) * CGC(J, M, j_i, m_i, j_c, m_c) * CGC(j_a, m_a, j_b, m_b, J, M) * Eff_IntsJ.Xpphp.X_3_4[chanind_J + (keyJ1*ChanJ.np[chanJ1] + keyJ2)];
	}
	if( std::fabs((xJ - Eff_Ints.Xpphp.X_3_4[chanind + ind3])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
	  std::cout << "!! Xpphp3_4 < " << a << "," << b << " |X| " << i << "," << c << " > !! -- " << Eff_Ints.Xpphp.X_3_4[chanind + ind3] << " " << xJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_Xhppp3_2_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind3, chan3;
  int np, npph;
  int i, a, b, c;
  one_body p;
  three_body pph;
  double xJ;
  int key_i, key_a, key_b, key_c;
  int shell_i, shell_a, shell_b, shell_c;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_i, state_a, state_b, state_c;
  int j_i, j_a, j_b, j_c;
  int m_i, m_a, m_b, m_c, M;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    chanind = Eff_Ints.Xhppp.X_3_2_index[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
	pph = Chan.pph_state(chan3, pph_ind);
	b = pph.v1;
	state_b = Space.qnums[b];
	state_b.j = Space.shellsj[b];
	key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
	shell_b = SpaceJ.map_state[key_b];
	c = pph.v2;
	state_c = Space.qnums[c];
	state_c.j = Space.shellsj[c];
	key_c = SpaceJ.ind_state(ParametersJ.basis, state_c);
	shell_c = SpaceJ.map_state[key_c];
	i = pph.v3;
	state_i = Space.qnums[i];
	state_i.j = Space.shellsj[i];
	key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
	shell_i = SpaceJ.map_state[key_i];
	if( b == c ){ continue; }
	ind3 = p_ind * npph + pph_ind;

	j_i = Space.shellsj[i];
	m_i = -1 * Space.qnums[i].m;
	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_c = Space.shellsj[c];
	m_c = Space.qnums[c].m;

	minJ = abs(j_i - j_a);
	if( abs(j_b - j_c) > minJ ){ minJ = abs(j_b - j_c); }
	maxJ = j_i + j_a;
	if( j_b + j_c < maxJ ){ maxJ = j_b + j_c; }
	M = m_b + m_c;

	xJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  chanJ1 = SpaceJ.ind_1b(ParametersJ.basis, SpaceJ.qnums[shell_a]);
	  plus(stateJ1, SpaceJ.qnums[shell_b], SpaceJ.qnums[shell_c]);
	  stateJ1.j = J;
	  if(shell_b == shell_c && J%4 != 0){ continue; }
	  minus(stateJ2, stateJ1, SpaceJ.qnums[shell_i]);
	  stateJ2.j = SpaceJ.qnums[shell_a].j;
	  chanJ2 = SpaceJ.ind_1b(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = Eff_IntsJ.Xhppp.X_3_2_index[chanJ1];
	  keyJ1 = ChanJ.p_map[chanJ1][shell_a];
	  keyJ2 = ChanJ.pph_map[chanJ2][SpaceJ.hash3(shell_b, shell_c, shell_i, J)];
	  xJ += -1.0 * phase2(j_i + m_i) * CGC(J, M, j_i, m_i, j_a, m_a) * CGC(j_b, m_b, j_c, m_c, J, M) * Eff_IntsJ.Xhppp.X_3_2[chanind_J + (keyJ1*ChanJ.npph[chanJ1] + keyJ2)];
	}
	if( std::fabs((xJ - Eff_Ints.Xhppp.X_3_2[chanind + ind3])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
	  std::cout << "!! Xhppp3_2 < " << i << "," << a << " |X| " << b << "," << c << " > !! -- " << Eff_Ints.Xhppp.X_3_2[chanind + ind3] << " " << xJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_Xhphp2_1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind2, chan2;
  int nhp1;
  int i, a, j, b;
  two_body hp1, hp2;
  double xJ;
  int key_i, key_a, key_j, key_b;
  int shell_i, shell_a, shell_j, shell_b;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_i, state_a, state_j, state_b;
  int j_i, j_a, j_j, j_b;
  int m_i, m_a, m_j, m_b, M;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    chanind = Eff_Ints.Xhphp.X_2_1_index[chan2];
    for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
      hp1 = Chan.hp1_state(chan2, hp1_ind);
      i = hp1.v1;
      state_i = Space.qnums[i];
      state_i.j = Space.shellsj[i];
      key_i = SpaceJ.ind_state(ParametersJ.basis, state_i);
      shell_i = SpaceJ.map_state[key_i];
      b = hp1.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      for(int hp2_ind = 0; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	j = hp2.v1;
	state_j = Space.qnums[j];
	state_j.j = Space.shellsj[j];
	key_j = SpaceJ.ind_state(ParametersJ.basis, state_j);
	shell_j = SpaceJ.map_state[key_j];
	a = hp2.v2;
	state_a = Space.qnums[a];
	state_a.j = Space.shellsj[a];
	key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
	shell_a = SpaceJ.map_state[key_a];
	ind2 = hp1_ind * nhp1 + hp2_ind;

	j_i = Space.shellsj[i];
	m_i = Space.qnums[i].m;
	j_a = Space.shellsj[a];
	m_a = -1 * Space.qnums[a].m;
	j_j = Space.shellsj[j];
	m_j = Space.qnums[j].m;
	j_b = Space.shellsj[b];
	m_b = -1 * Space.qnums[b].m;

	minJ = abs(j_i - j_b);
	if( abs(j_j - j_a) > minJ ){ minJ = abs(j_j - j_a); }
	maxJ = j_i + j_b;
	if( j_j + j_a < maxJ ){ maxJ = j_j + j_a; }
	M = m_i + m_b;

	xJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  minus(stateJ1, SpaceJ.qnums[shell_i], SpaceJ.qnums[shell_b]);
	  minus(stateJ2, SpaceJ.qnums[shell_j], SpaceJ.qnums[shell_a]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  chanJ1 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_cross(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = Eff_IntsJ.Xhphp.X_2_1_index[chanJ1];
	  keyJ1 = ChanJ.hp1_map[chanJ1][SpaceJ.hash2(shell_i, shell_b, J)];
	  keyJ2 = ChanJ.hp1_map[chanJ2][SpaceJ.hash2(shell_j, shell_a, J)];
	  xJ += -1.0 * phase2(j_a + m_a + j_b + m_b) * CGC(j_i, m_i, j_b, m_b, J, M) * CGC(j_j, m_j, j_a, m_a, J, M) * Eff_IntsJ.Xhphp.X_2_1[chanind_J + (keyJ1*ChanJ.nhp1[chanJ1] + keyJ2)];
	}
	if( std::fabs((xJ - Eff_Ints.Xhphp.X_2_1[chanind + ind2])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
	  std::cout << "!! Xhphp2_1 < " << i << "," << a << " |X| " << j << "," << b << " > -- " << Eff_Ints.Xhphp.X_2_1[chanind + ind2] << " " << xJ << std::endl;
	}
      }
    }
  }
}

void CC_Test_Xpppp1_J(Input_Parameters &ParametersJ, Model_Space &SpaceJ, Channels &ChanJ, Eff_Interactions &Eff_IntsJ, Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  int ind1, chan1;
  int npp;
  int a, b, c, d;
  two_body pp1, pp2;
  double xJ;
  int key_a, key_b, key_c, key_d;
  int shell_a, shell_b, shell_c, shell_d;
  int minJ, maxJ;
  int chanind, chanJ1, chanJ2, chanind_J, keyJ1, keyJ2;
  State stateJ1, stateJ2, state_a, state_b, state_c, state_d;
  int j_a, j_b, j_c, j_d;
  int m_a, m_b, m_c, m_d, M;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    chanind = Eff_Ints.Xpppp.X_1_index[chan1];
    for(int pp_ind1 = 0; pp_ind1 < npp; ++pp_ind1){
      pp1 = Chan.pp_state(chan1, pp_ind1);
      a = pp1.v1;
      state_a = Space.qnums[a];
      state_a.j = Space.shellsj[a];
      key_a = SpaceJ.ind_state(ParametersJ.basis, state_a);
      shell_a = SpaceJ.map_state[key_a];
      b = pp1.v2;
      state_b = Space.qnums[b];
      state_b.j = Space.shellsj[b];
      key_b = SpaceJ.ind_state(ParametersJ.basis, state_b);
      shell_b = SpaceJ.map_state[key_b];
      if(a == b){ continue; }
      for(int pp_ind2 = 0; pp_ind2 < npp; ++pp_ind2){
	pp2 = Chan.pp_state(chan1, pp_ind2);
	c = pp2.v1;
	state_c = Space.qnums[c];
	state_c.j = Space.shellsj[c];
	key_c = SpaceJ.ind_state(ParametersJ.basis, state_c);
	shell_c = SpaceJ.map_state[key_c];
	d = pp2.v2;
	state_d = Space.qnums[d];
	state_d.j = Space.shellsj[d];
	key_d = SpaceJ.ind_state(ParametersJ.basis, state_d);
	shell_d = SpaceJ.map_state[key_d];
	if(c == d){ continue; }
	ind1 = pp_ind1 * npp + pp_ind2;

	j_a = Space.shellsj[a];
	m_a = Space.qnums[a].m;
	j_b = Space.shellsj[b];
	m_b = Space.qnums[b].m;
	j_c = Space.shellsj[c];
	m_c = Space.qnums[c].m;
	j_d = Space.shellsj[d];
	m_d = Space.qnums[d].m;

	minJ = abs(j_a - j_b);
	if( abs(j_c - j_d) > minJ ){ minJ = abs(j_c - j_d); }
	maxJ = j_a + j_b;
	if( j_c + j_d < maxJ ){ maxJ = j_c + j_d; }
	M = m_a + m_b;

	xJ = 0.0;
	for(int J = minJ; J <= maxJ; J+=2){
	  if( std::abs(M) > J ){ continue; }
	  plus(stateJ1, SpaceJ.qnums[shell_a], SpaceJ.qnums[shell_b]);
	  plus(stateJ2, SpaceJ.qnums[shell_c], SpaceJ.qnums[shell_d]);
	  stateJ1.j = J;
	  stateJ2.j = J;
	  if((shell_a == shell_b || shell_c == shell_d) && J%4 != 0){ continue; }
	  chanJ1 = SpaceJ.ind_2b_dir(ParametersJ.basis, stateJ1);
	  chanJ2 = SpaceJ.ind_2b_dir(ParametersJ.basis, stateJ2);
	  if( chanJ1 != chanJ2 ){ continue; }
	  chanind_J = Eff_IntsJ.Xpppp.X_1_index[chanJ1];
	  keyJ1 = ChanJ.pp_map[chanJ1][SpaceJ.hash2(shell_a, shell_b, J)];
	  keyJ2 = ChanJ.pp_map[chanJ2][SpaceJ.hash2(shell_c, shell_d, J)];
	  xJ += CGC(j_a, m_a, j_b, m_b, J, M) * CGC(j_c, m_c, j_d, m_d, J, M) * Eff_IntsJ.Xpppp.X_1[chanind_J + (keyJ1*ChanJ.npp[chanJ1] + keyJ2)];
	}
	if( std::fabs((xJ - Eff_Ints.Xpppp.X_1[chanind + ind1])/(std::fabs(xJ) + 1e-12)) > 1e-8 && std::fabs(xJ) > 1e-12 ){
	  std::cout << "!! Xpppp1 < " << a << "," << b << " |X| " << c << "," << d << " > !! -- " << Eff_Ints.Xpppp.X_1[chanind + ind1] << " " << xJ << std::endl;
	}
      }
    }
  }
}
