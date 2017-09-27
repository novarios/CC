#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

double E_Ref(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  double energy = 0.0;
  int nh, nhh, chan_ind, ind;
  State tb;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    chan_ind = Ints.Vhhhh.V_1_index[chan];
    for(int hh = 0; hh < nhh; ++hh){
      ind = hh*nhh + hh;
      if(Parameters.basis == "finite_J"){ energy += (Chan.qnums1[chan].j + 1) * Ints.Vhhhh.V_1[chan_ind + ind]; }
      else{ energy += Ints.Vhhhh.V_1[chan_ind + ind]; }
    }
  }
  energy *= -0.5;

  if(Parameters.HF == 1){
    for(int i = 0; i < Space.num_hol; ++i){
      if(Parameters.basis == "finite_J"){ energy += (Space.qnums[i].j + 1) * Space.qnums[i].energy; }
      else{ energy += Space.qnums[i].energy; }
    }
  }
  else if(Parameters.HF == 0){
    for(int chan = 0; chan < Chan.size3; ++chan){
      nh = Chan.nh[chan];
      chan_ind = Ints.Fmatrix.hh_3_index[chan];
      for(int h = 0; h < nh; ++h){
	ind = h*nh + h;
	if(Parameters.basis == "finite_J"){ energy += (Chan.qnums3[chan].j + 1) * Ints.Fmatrix.hh_3[chan_ind + ind]; }
	else{ energy += Ints.Fmatrix.hh_3[chan_ind + ind]; }
      }
    }
  }
  return energy;
}

double Amplitudes::get_energy(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  double energy = 0.0;
  int nhh, npp, ind1, ind2, chan_ind1, chan_ind2;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    chan_ind1 = D1.T1_index[chan];
    chan_ind2 = Ints.Vhhpp.V_1_index[chan];
    for(int pp = 0; pp < npp; ++pp){
      for(int hh = 0; hh < nhh; ++hh){
	ind1 = pp * nhh + hh;
	ind2 = hh * npp + pp;
	if(Parameters.basis == "finite_J"){ energy += (Chan.qnums1[chan].j + 1) * D1.T1[chan_ind1 + ind1] * Ints.Vhhpp.V_1[chan_ind2 + ind2]; }
	else{ energy += D1.T1[chan_ind1 + ind1] * Ints.Vhhpp.V_1[chan_ind2 + ind2]; }
      }
    }
  }
  energy *= 0.25;
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    double energy1 = 0.0;
    int a, i, hp_ind;
    int nph0 = Chan.nph1[Chan.ind0];
    two_body ph;
    chan_ind1 = Ints.Vhhpp.V_2_3_index[Chan.ind0];
    for(int ph1 = 0; ph1 < nph0; ++ph1){
      ph = Chan.ph1_state(Chan.ind0, ph1);
      a = ph.v1;
      i = ph.v2;
      hp_ind = Chan.hp1_map[Chan.ind0][Space.hash2(i, a, Chan.qnums2[Chan.ind0].j)];
      for(int ph2 = 0; ph2 < nph0; ++ph2){
	ind1 = hp_ind * nph0 + ph2;
	energy1 += S1.t2[ph1] * S1.t2[ph2] * Ints.Vhhpp.V_2_3[chan_ind1 + ind1];
      }
    }
    energy += 0.5 * energy1;
  }
  if(Parameters.HF == 0){
    int np, nh;
    for(int chan = 0; chan < Chan.size3; ++chan){
      np = Chan.np[chan];
      nh = Chan.nh[chan];
      chan_ind1 = Ints.Fmatrix.hp_3_index[chan];
      chan_ind2 = S1.t3_index[chan];
      for(int p = 0; p < np; ++p){
	for(int h = 0; h < nh; ++h){
	  ind1 = h * np + p;
	  ind2 = p * nh + h;
	  energy += Ints.Fmatrix.hp_3[chan_ind1 + ind1] * S1.t3[chan_ind2 + ind2];
	}
      }
    }
  }
  return energy;
}

void Amplitudes::copy_Amplitudes(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps)
{
  if(Parameters.approx == "doubles"){
    D1.copy_Doubles_1(Chan, Amps.D1);
  }
  else if(Parameters.approx == "singles"){
    D1.copy_Doubles_1(Chan, Amps.D1);
    S1.copy_Singles_1(Chan, Amps.S1);
  }
}

Amplitudes::Amplitudes(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1 = Doubles_1(Parameters, Space, Chan);
  }
  else if(Parameters.approx == "singles"){
    D1 = Doubles_1(Parameters, Space, Chan);
    S1 = Singles_1(Parameters, Space, Chan);
  }
}

void Amplitudes::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.delete_struct(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.delete_struct(Chan);
    S1.delete_struct(Chan);
  }
}

void Amplitudes::zero(Input_Parameters &Parameters, Channels &Chan, bool t1flag)
{
  if(Parameters.approx == "doubles"){
    D1.zero(Chan, t1flag);
  }
  else if(Parameters.approx == "singles"){
    D1.zero(Chan, t1flag);
    S1.zero(Chan, t1flag);
  }
}

void Doubles_1::copy_Doubles_1(Channels &Chan, Doubles_1 &D)
{
  for(int ind = 0; ind < T1_length; ++ind){ T1[ind] = D.T1[ind]; }
  for(int ind = 0; ind < T2_1_length; ++ind){ T2_1[ind] = D.T2_1[ind]; }
  for(int ind = 0; ind < T2_2_length; ++ind){ T2_2[ind] = D.T2_2[ind]; }
  for(int ind = 0; ind < T2_3_length; ++ind){ T2_3[ind] = D.T2_3[ind]; }
  for(int ind = 0; ind < T2_4_length; ++ind){ T2_4[ind] = D.T2_4[ind]; }
  for(int ind = 0; ind < T3_1_length; ++ind){ T3_1[ind] = D.T3_1[ind]; }
  for(int ind = 0; ind < T3_2_length; ++ind){ T3_2[ind] = D.T3_2[ind]; }
  for(int ind = 0; ind < T3_3_length; ++ind){ T3_3[ind] = D.T3_3[ind]; }
  for(int ind = 0; ind < T3_4_length; ++ind){ T3_4[ind] = D.T3_4[ind]; }
}

Doubles_1::Doubles_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count, chan;
  int a, b, i, j, npp, nhh;
  two_body pp, hh;
  four_body abij, abij_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;

  length = 0;
  T1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    T1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhh[chan1];
  }
  T1 = new double[length];
  map_num = new int[8 * length]; // T2_1, T2_2, T2_3, T2_4, T3_1, T3_2, T3_3, T3_4
  map_index = new int[8 * length];
  Evec_chan = new int[4 * length];
  Evec_ind = new int[4 * length];
  fb_ind = new four_body[length];
  fb_j = new four_body[length];
  J = new int[length];
  for(ind = 0; ind < length; ++ind){
    T1[ind] = 0.0;
  }
  T1_length = length;

  length = 0;
  T2_1_index = new int[Chan.size2];
  T2_2_index = new int[Chan.size2];
  T2_3_index = new int[Chan.size2];
  T2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    T2_1_index[chan2] = length;
    T2_2_index[chan2] = length;
    T2_3_index[chan2] = length;
    T2_4_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.nhp1[chan2];
  }
  T2_1 = new double[length];
  T2_2 = new double[length];
  T2_3 = new double[length];
  T2_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    T2_1[ind] = 0.0;
    T2_2[ind] = 0.0;
    T2_3[ind] = 0.0;
    T2_4[ind] = 0.0;
  }
  T2_1_length = length;
  T2_2_length = length;
  T2_3_length = length;
  T2_4_length = length;

  length = 0;
  T3_1_index = new int[Chan.size3];
  T3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    T3_1_index[chan3] = length;
    T3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhp[chan3];
  }
  T3_1 = new double[length];
  T3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    T3_1[ind] = 0.0;
    T3_2[ind] = 0.0;
  }
  T3_1_length = length;
  T3_2_length = length;

  length = 0;
  T3_3_index = new int[Chan.size3];
  T3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    T3_3_index[chan3] = length;
    T3_4_index[chan3] = length;
    length += Chan.npph[chan3] * Chan.nh[chan3];
  }
  T3_3 = new double[length];
  T3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    T3_3[ind] = 0.0;
    T3_4[ind] = 0.0;
  }
  T3_3_length = length;
  T3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      pp = Chan.pp_state(chan1, pp_ind);
      a = pp.v1;
      b = pp.v2;
      abij.v1 = a;
      abij.v2 = b;
      abij_j.v1 = Space.qnums[a].j;
      abij_j.v2 = Space.qnums[b].j;
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	i = hh.v1;
	j = hh.v2;
	abij.v3 = i;
	abij.v4 = j;
	abij_j.v3 = Space.qnums[i].j;
	abij_j.v4 = Space.qnums[j].j;
	fb_ind[count] = abij;
	fb_j[count] = abij_j;
	J[count] = Chan.qnums1[chan1].j;
	// T2_1, T2_2, T2_3, T2_4, T3_1, T3_2, T3_3, T3_4
	Map_4_count_1(Parameters,Space, map_index,map_num, 8,0,length,count, abij);
	Map_4_count_2(Parameters,Space, map_index,map_num, 8,1,length,count, abij);
	Map_4_count_3(Parameters,Space, map_index,map_num, 8,2,length,count, abij);
	Map_4_count_4(Parameters,Space, map_index,map_num, 8,3,length,count, abij);
	Map_4_count_5678(Parameters,Space, map_index,map_num, 8,4,length,count);
	Map_4_count_5678(Parameters,Space, map_index,map_num, 8,5,length,count);
	Map_4_count_5678(Parameters,Space, map_index,map_num, 8,6,length,count);
	Map_4_count_5678(Parameters,Space, map_index,map_num, 8,7,length,count);
	++count;
      }
    }
  }
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];

  #pragma omp parallel for schedule(static) private(a, b, i, j, chan)
  for(int pphh = 0; pphh < T1_length; ++pphh){
    // T2_1, T2_2, T2_3, T2_4, T3_1, T3_2, T3_3, T3_4
    Map_4_1(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,0,pphh, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, fb_ind,fb_j,J);
    Map_4_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,1,pphh, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, fb_ind,fb_j,J);
    Map_4_3(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,2,pphh, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, fb_ind,fb_j,J);
    Map_4_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,3,pphh, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, fb_ind,fb_j,J);
    Map_4_5(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,4,pphh, Chan.p_map,Chan.hhp_map, Chan.nhhp, fb_ind,fb_j,J);
    Map_4_6(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,5,pphh, Chan.p_map,Chan.hhp_map, Chan.nhhp, fb_ind,fb_j,J);
    Map_4_7(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,6,pphh, Chan.pph_map,Chan.h_map, Chan.nh, fb_ind,fb_j,J);
    Map_4_8(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,7,pphh, Chan.pph_map,Chan.h_map, Chan.nh, fb_ind,fb_j,J);

    a = fb_ind[pphh].v1;
    b = fb_ind[pphh].v2;
    i = fb_ind[pphh].v3;
    j = fb_ind[pphh].v4;
    chan = Space.ind_1b(Parameters.basis, Space.qnums[i]);
    Evec_chan[4*pphh] = chan;
    Evec_ind[4*pphh] = Chan.h_map[chan][i];
    chan = Space.ind_1b(Parameters.basis, Space.qnums[j]);
    Evec_chan[4*pphh + 1] = chan;
    Evec_ind[4*pphh + 1] = Chan.h_map[chan][j];
    chan = Space.ind_1b(Parameters.basis, Space.qnums[a]);
    Evec_chan[4*pphh + 2] = chan;
    Evec_ind[4*pphh + 2] = Chan.p_map[chan][a];
    chan = Space.ind_1b(Parameters.basis, Space.qnums[b]);
    Evec_chan[4*pphh + 3] = chan;
    Evec_ind[4*pphh + 3] = Chan.p_map[chan][b];
  }
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
}

void Doubles_1::delete_struct(Channels &Chan)
{
  delete[] T1;
  delete[] T2_1;
  delete[] T2_2;
  delete[] T2_3;
  delete[] T2_4;
  delete[] T3_1;
  delete[] T3_2;
  delete[] T3_3;
  delete[] T3_4;
  delete[] Evec_chan;
  delete[] Evec_ind;

  delete[] T1_index;
  delete[] T2_1_index;
  delete[] T2_2_index;
  delete[] T2_3_index;
  delete[] T2_4_index;
  delete[] T3_1_index;
  delete[] T3_2_index;
  delete[] T3_3_index;
  delete[] T3_4_index;
  
  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

void Doubles_1::zero(Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < T1_length; ++ind){
      T1[ind] = 0.0;
    }
  }
  for(ind = 0; ind < T2_1_length; ++ind){
    T2_1[ind] = 0.0;
  }
  for(ind = 0; ind < T2_2_length; ++ind){
    T2_2[ind] = 0.0;
  }
  for(ind = 0; ind < T2_3_length; ++ind){
    T2_3[ind] = 0.0;
  }
  for(ind = 0; ind < T2_4_length; ++ind){
    T2_4[ind] = 0.0;
  }
  for(ind = 0; ind < T3_1_length; ++ind){
    T3_1[ind] = 0.0;
  }
  for(ind = 0; ind < T3_2_length; ++ind){
    T3_2[ind] = 0.0;
  }
  for(ind = 0; ind < T3_3_length; ++ind){
    T3_3[ind] = 0.0;
  }
  for(ind = 0; ind < T3_4_length; ++ind){
    T3_4[ind] = 0.0;
  }
}

void Doubles_1::set_T(int &ind, double &T)
{
  int index;
  T1[ind] = T;
  for(int n = 0; n < map_num[8*ind]; ++n){
    index = map_index[8*ind] + n;
    T2_1[T2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 1]; ++n){
    index = map_index[8*ind + 1] + n;
    T2_2[T2_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 2]; ++n){
    index = map_index[8*ind + 2] + n;
    T2_3[T2_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 3]; ++n){
    index = map_index[8*ind + 3] + n;
    T2_4[T2_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 4]; ++n){
    index = map_index[8*ind + 4] + n;
    T3_1[T3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 5]; ++n){
    index = map_index[8*ind + 5] + n;
    T3_2[T3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 6]; ++n){
    index = map_index[8*ind + 6] + n;
    T3_3[T3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
  for(int n = 0; n < map_num[8*ind + 7]; ++n){
    index = map_index[8*ind + 7] + n;
    T3_4[T3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * T;
  }
}

double Doubles_1::get_T(int &ind)
{
  int index;
  double T = 0.0;
  T += T1[ind];
  for(int n = 0; n < map_num[8*ind]; ++n){
    index = map_index[8*ind] + n;
    T += map_fac2[index] * T2_1[T2_1_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 1]; ++n){
    index = map_index[8*ind + 1] + n;
    T += map_fac2[index] * T2_2[T2_2_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 2]; ++n){
    index = map_index[8*ind + 2] + n;
    T += map_fac2[index] * T2_3[T2_3_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 3]; ++n){
    index = map_index[8*ind + 3] + n;
    T += map_fac2[index] * T2_4[T2_4_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 4]; ++n){
    index = map_index[8*ind + 4] + n;
    T += map_fac2[index] * T3_1[T3_1_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 5]; ++n){
    index = map_index[8*ind + 5] + n;
    T += map_fac2[index] * T3_2[T3_2_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 6]; ++n){
    index = map_index[8*ind + 6] + n;
    T += map_fac2[index] * T3_3[T3_3_index[map_chan[index]] + map_ind[index]];
  }
  for(int n = 0; n < map_num[8*ind + 7]; ++n){
    index = map_index[8*ind + 7] + n;
    T += map_fac2[index] * T3_4[T3_4_index[map_chan[index]] + map_ind[index]];
  }
  return T;
}

void Singles_1::copy_Singles_1(Channels &Chan, Singles_1 &S)
{
  for(int ind = 0; ind < t2_length; ++ind){ t2[ind] = S.t2[ind]; }
  for(int ind = 0; ind < t3_length; ++ind){ t3[ind] = S.t3[ind]; }
}

Singles_1::Singles_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan3, ind, length;
  int a, i;
  two_body ph1;

  length = Chan.nph1[Chan.ind0];
  t2 = new double[length];
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];
  evec_chan = new int[2 * length];
  evec_ind = new int[2 * length];
  for(ind = 0; ind < length; ++ind){
    t2[ind] = 0.0;
  }
  t2_length = length;

  length = 0;
  t3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    t3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nh[chan3];
  }
  t3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    t3[ind] = 0.0;
  }
  t3_length = length;

  for(ind = 0; ind < t2_length; ++ind){
    ph1 = Chan.ph1_state(Chan.ind0, ind);
    a = ph1.v1;
    i = ph1.v2;
    Map_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, ind, Chan.p_map,Chan.h_map,Chan.nh, ph1, Space.qnums[a].j);

    chan3 = Space.ind_1b(Parameters.basis, Space.qnums[i]);
    evec_chan[2*ind] = chan3;
    evec_ind[2*ind] = Chan.h_map[chan3][i];
    chan3 = Space.ind_1b(Parameters.basis, Space.qnums[a]);
    evec_chan[2*ind + 1] = chan3;
    evec_ind[2*ind + 1] = Chan.p_map[chan3][a];
  }
}

void Singles_1::delete_struct(Channels &Chan)
{
  delete[] t2;
  delete[] t3;
  delete[] evec_chan;
  delete[] evec_ind;

  delete[] t3_index;

  delete[] map_chan;
  delete[] map_ind;
  delete[] map_fac1;
  delete[] map_fac2;
}

void Singles_1::zero(Channels &Chan, bool flag)
{
  if( !flag ){
    for(int ind = 0; ind < t2_length; ++ind){ t2[ind] = 0.0; }
  }
  for(int ind = 0; ind < t3_length; ++ind){ t3[ind] = 0.0; }
}

void Singles_1::set_T(int &ind, double &t)
{
  t2[ind] = t;
  t3[t3_index[map_chan[ind]] + map_ind[ind]] += map_fac1[ind] * t;
}

double Singles_1::get_T(int &ind)
{
  double t = 0.0;
  t += t2[ind];
  t += map_fac2[ind] * t3[t3_index[map_chan[ind]] + map_ind[ind]];
  return t;
}

//Initialize program from input file
void Get_Input_Parameters(std::string &infile, Input_Parameters &Input)
{ 
  std::string path; // Input File Path
  std::string line; // Line from input file
  std::ifstream filestream; // File stream
  int index; // Count input line
  size_t colon; // Size for finding ':' which proceeds every input
  std::string substr; // Substd::string for extracting input

  path = PATH + infile;
  filestream.open(path.c_str());
  if (!filestream.is_open()){ std::cerr << "Input file, " << path << ", does not exist" << std::endl; exit(1); }

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
	Input.calc_case = substr;
	break;
      case 2:
	Input.basis = substr;
	break;
      case 3:
	Input.approx = substr;
	break;
      case 4:
	Input.obstrength = atof(substr.c_str());
	break;
      case 5:
	Input.tbstrength = atof(substr.c_str());
	break;
      case 6:
	Input.Pshells = atoi(substr.c_str());
	break;
      case 7:
	Input.Nshells = atoi(substr.c_str());
	break;
      case 8:
	Input.LevelScheme = substr;
	break;
      case 9:
	Input.MatrixElements = substr;
	break;
      case 10:
	Input.extra = atoi(substr.c_str());
	break;
      } 
    }
    else{ continue; };
  }
}

void Print_Parameters(Input_Parameters &Parameters, Model_Space &Space)
{
  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << "Case = " << Parameters.calc_case << ", Basis = " << Parameters.basis << ", Approximation = " << Parameters.approx << std::endl;
  if(Parameters.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << Parameters.LevelScheme << std::endl;
    if(Parameters.MatrixElements.size() > 0){ std::cout << "Interaction = " << Parameters.MatrixElements << std::endl; }
  }
  if(Parameters.calc_case == "nuclear"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.num_states << std::endl;
    std::cout << "Proton Shells = " << Parameters.Pshells << ", Neutron Shells = " << Parameters.Nshells << std::endl;
    std::cout << "Protons = " << Parameters.P << ", Neutrons = " << Parameters.N << std::endl;
    if(Parameters.calc_case == "infinite"){ std::cout << "Density = " << Parameters.density << std::endl; }
  }
  else if(Parameters.calc_case == "electronic"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.num_states << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Density = " << Parameters.density << std::endl;
  }
  else if(Parameters.calc_case == "quantum_dot"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.num_states << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Oscillator Energy = " << Parameters.density << std::endl;
  }
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void Perform_CC(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  double mixmax = 0.75;
  double mixmin = 0.01;
  double mix;
  double width = 1.0;
  int Rand_count = 0;

  double error = 1e20, error2 = 1e20;
  int ind = 0;
  double CCoutE;
  Amplitudes Amps0 = Amplitudes(Parameters, Space, Chan);
  Amplitudes Amps2 = Amplitudes(Parameters, Space, Chan);
  Amplitudes tempAmps = Amplitudes(Parameters, Space, Chan);

  /// For DIIS ///
  int DIIS_count = 0;
  double DIISstart = 0.01;
  int maxl = 5;
  int N = 0;
  double *p = NULL;
  double *delp = NULL;
  double *tempdelp = NULL;
  double *B = NULL;
  Initialize_DIIS(Parameters, Chan, Amps, p, delp, tempdelp, B, maxl);
  ///////////////

  // Initialize Amplitudes //
  mix = mixmin;
  while(error > 1e-12 && ind < 2000){
    Update_CC(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);

    Amps2.copy_Amplitudes(Parameters, Chan, Amps);
    Amps2.zero(Parameters, Chan, true);
    Gather_Amps(Parameters, Space, Chan, Eff_Ints, Amps2, tempAmps, mix);
    
    CC_Error(Parameters, Chan, Ints, Amps, Amps2, error);
    error /= mix;
    if( !std::isfinite(error) ){ std::cerr << std::endl << ind << ", error = " << error << " : CC Solution Diverged!! " << std::endl; exit(1); }
    if( ind > 1000 && double(Rand_count)/double(ind) > 0.9 ){ std::cerr << std::endl << ind << " : CC Solution Not Converged!!" << std::endl; exit(1); }
    
    if( error < error2 || error > 1.0 ){
      if( error < error2 ){ mix = std::pow(mixmax, 0.035) * std::pow(mix, 0.965); }
      Amps0.copy_Amplitudes(Parameters, Chan, Amps);
      Amps.copy_Amplitudes(Parameters, Chan, Amps2);
      error2 = error;
      if( error < DIISstart && mix > 0.2 * mixmax ){ Perform_DIIS(Parameters, Chan, Amps, Amps0, mix, p, delp, tempdelp, B, N, maxl, DIIS_count); }
    }
    else{
      ++Rand_count;
      if(error2 > 1.0){ width = 0.001; }
      else{ width = 0.001 * error2; }
      Random_Step(Parameters, Space, Chan, Ints, Eff_Ints, Amps0, Amps, Amps2, tempAmps, mix, width, error, error2);
    }

    CCoutE = Amps.get_energy(Parameters, Space, Chan, Ints);
    if( !std::isfinite(CCoutE) ){ std::cerr << std::endl << ind << ", Energy = " << CCoutE << " : CC Solution Diverged!! " << std::endl; exit(1); }
    std::cout << "Iteration Number = " << ind << ", Energy = " << CCoutE << ", error = " << error << ", mix = " << mix << ", ";
    std::cout << "Random = " << Rand_count << ", DIIS = " << DIIS_count << std::endl;
    ++ind;
  }

  std::cout << std::endl << std::endl;
  if( error > 1e-12 ){
    std::cout << Parameters.Shells << ", " << Parameters.Pshells << ", " << Parameters.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". CC Solution Not Converged!!" << std::endl;
  }
  Amps0.delete_struct(Parameters, Chan);
  Amps2.delete_struct(Parameters, Chan);
  tempAmps.delete_struct(Parameters, Chan);
  Delete_DIIS(Parameters, Chan, p, delp, tempdelp, B, maxl);
}

void Perform_CC_Test(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  double mix = 1.0;
  int ind = 0;
  double CCoutE, CCoutEM;
  Amplitudes Amps2 = Amplitudes(Parameters, Space, Chan);
  Amplitudes tempAmps = Amplitudes(Parameters, Space, Chan);

  /// For J-Testing ///
  Input_Parameters ParametersM;
  Model_Space SpaceM;
  Channels ChanM;
  Amplitudes AmpsM;
  Interactions IntsM;
  Eff_Interactions Eff_IntsM;
  HF_Channels HF_ChanM;
  HF_Matrix_Elements HF_MEM;
  Single_Particle_States StatesM;
  ParametersM = Parameters;
  ParametersM.basis = "finite_JM";
  Build_Model_Space(ParametersM, SpaceM);
  HF_ChanM = HF_Channels(ParametersM, SpaceM);
  StatesM = Single_Particle_States(ParametersM, SpaceM, HF_ChanM);
  Read_Matrix_Elements_J(ParametersM, SpaceM, HF_ChanM, HF_MEM);
  if(Parameters.HF == 1){
    Hartree_Fock_States(ParametersM, SpaceM, HF_ChanM, StatesM, HF_MEM);
    Convert_To_HF_Matrix_Elements(ParametersM, HF_ChanM, SpaceM, StatesM, HF_MEM);
  }
  ParametersM.basis = "finite_M";
  Build_Model_Space_J2(ParametersM, SpaceM);
  ChanM = Channels(ParametersM, SpaceM);
  AmpsM = Amplitudes(ParametersM, SpaceM, ChanM);
  IntsM = Interactions(ParametersM, ChanM);
  if(Parameters.HF == 0){ Get_Fock_Matrix(ParametersM, HF_ChanM, HF_MEM, StatesM, SpaceM, ChanM, IntsM); }
  Eff_IntsM = Eff_Interactions(ParametersM, SpaceM, ChanM);
  Get_Matrix_Elements_JM(ParametersM, HF_ChanM, HF_MEM, SpaceM, ChanM, IntsM);
  Amplitudes Amps2M = Amplitudes(ParametersM, SpaceM, ChanM);
  Amplitudes tempAmpsM = Amplitudes(ParametersM, SpaceM, ChanM);
  /////////////////////

  while( ind < 20 ){
    Update_CC(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);           // J
    Update_CC(ParametersM, SpaceM, ChanM, IntsM, Eff_IntsM, AmpsM, tempAmpsM);    // M

    //Print_Amps(ParametersM, ChanM, tempAmpsM);

    std::cout << " ****  tempAmps  **** " << std::endl;
    CC_Test_Full_J(Parameters, Space, Chan, tempAmps, ParametersM, SpaceM, ChanM, tempAmpsM);

    Amps2.copy_Amplitudes(Parameters, Chan, Amps);                         // J
    Amps2.zero(Parameters, Chan, true);
    Gather_Amps(Parameters, Space, Chan, Eff_Ints, Amps2, tempAmps, mix);
    Amps2M.copy_Amplitudes(ParametersM, ChanM, AmpsM);                     // M
    Amps2M.zero(ParametersM, ChanM, true);
    Gather_Amps(ParametersM, SpaceM, ChanM, Eff_IntsM, Amps2M, tempAmpsM, mix);

    Amps.copy_Amplitudes(Parameters, Chan, Amps2);      // J
    AmpsM.copy_Amplitudes(ParametersM, ChanM, Amps2M);  // M

    //Print_Amps(Parameters, Chan, Amps);

    std::cout << " ****  Amps  **** " << std::endl;
    CC_Test_Full_J(Parameters, Space, Chan, Amps, ParametersM, SpaceM, ChanM, AmpsM);

    CCoutE = Amps.get_energy(Parameters, Space, Chan, Ints);        // J
    CCoutEM = AmpsM.get_energy(ParametersM, SpaceM, ChanM, IntsM);  // M

    std::cout << "Iteration Number = " << ind << ", EnergyJ = " << CCoutE << ", EnergyM = " << CCoutEM << std::endl;
    ++ind;
  }

  Update_Heff_3(Parameters, Space, Chan, Ints, Eff_Ints, Amps);       // J
  Update_Heff_3(ParametersM, SpaceM, ChanM, IntsM, Eff_IntsM, AmpsM); // M

  /*EOM_1P(Parameters, Space, Chan, Eff_Ints, CCoutE);
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  EOM_1P(ParametersM, SpaceM, ChanM, Eff_IntsM, CCoutEM);*/

  EOM EOM_J1, EOM_J2, EOM_M1, EOM_M2;
  EOM_J1.PA1_EOM(Parameters, Space, Chan, Eff_Ints);
  EOM_J1.Print_EOM_1P(Parameters, CCoutE);
  EOM_J2.PR1_EOM(Parameters, Space, Chan, Eff_Ints);
  EOM_J2.Print_EOM_1P(Parameters, CCoutE);
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  EOM_M1.PA1_EOM(ParametersM, SpaceM, ChanM, Eff_IntsM);
  EOM_M1.Print_EOM_1P(ParametersM, CCoutEM);
  EOM_M2.PR1_EOM(ParametersM, SpaceM, ChanM, Eff_IntsM);
  EOM_M2.Print_EOM_1P(ParametersM, CCoutEM);

  EOM_J1.delete_struct();
  EOM_J2.delete_struct();
  EOM_M1.delete_struct();
  EOM_M2.delete_struct();
  Amps2.delete_struct(Parameters, Chan);      // J
  tempAmps.delete_struct(Parameters, Chan);
  StatesM.delete_struct(HF_ChanM);            // M
  HF_MEM.delete_struct(HF_ChanM);
  HF_ChanM.delete_struct();
  IntsM.delete_struct(ParametersM, ChanM);
  AmpsM.delete_struct(ParametersM, ChanM);
  Amps2M.delete_struct(ParametersM, ChanM);
  tempAmpsM.delete_struct(ParametersM, ChanM);
  ChanM.delete_struct(ParametersM);
  SpaceM.delete_struct(ParametersM);
}

void Update_CC(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &tempAmps)
{
  Eff_Ints.zero(Parameters, Chan, false);
  tempAmps.zero(Parameters, Chan, false);
  //  Build Effective Hamiltonian
  Update_Heff_1(Parameters, Space, Chan, Ints, Eff_Ints, Amps);
  if(Parameters.approx == "singles"){
    Update_Heff_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps);
  }
  //  Solve Amplitude Equations
  Eff_Ints.diagonalize(Parameters, Chan);
  Doubles_Step(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
  if(Parameters.approx == "singles"){
    Singles_Step(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
    Doubles_Step_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
  }
}

void Gather_Amps(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix)
{
  double tempt, tempen;
  int length = Amps.D1.T1_length;
  for(int ind = 0; ind < length; ++ind){
    tempt = Amps2.D1.get_T(ind);
    tempen = 0.0;
    tempen += Eff_Ints.Xhh.X_3d[Eff_Ints.Xhh.X_3d_index[Amps.D1.Evec_chan[4*ind]] + Amps.D1.Evec_ind[4*ind]];
    tempen += Eff_Ints.Xhh.X_3d[Eff_Ints.Xhh.X_3d_index[Amps.D1.Evec_chan[4*ind + 1]] + Amps.D1.Evec_ind[4*ind + 1]];
    tempen -= Eff_Ints.Xpp.X_3d[Eff_Ints.Xpp.X_3d_index[Amps.D1.Evec_chan[4*ind + 2]] + Amps.D1.Evec_ind[4*ind + 2]];
    tempen -= Eff_Ints.Xpp.X_3d[Eff_Ints.Xpp.X_3d_index[Amps.D1.Evec_chan[4*ind + 3]] + Amps.D1.Evec_ind[4*ind + 3]];
    tempt /= tempen;
    tempt = mix*tempt + (1.0-mix)*Amps.D1.T1[ind];
    Amps.D1.set_T(ind, tempt);
  }
  if(Parameters.approx == "singles"){
    length = Amps.S1.t2_length;
    for(int ind = 0; ind < length; ++ind){
      tempt = Amps2.S1.get_T(ind);
      tempen = 0.0;
      tempen += Eff_Ints.Xhh.X1_3d[Eff_Ints.Xhh.X_3d_index[Amps.S1.evec_chan[2*ind]] + Amps.S1.evec_ind[2*ind]];
      tempen -= Eff_Ints.Xpp.X_3d[Eff_Ints.Xpp.X_3d_index[Amps.S1.evec_chan[2*ind + 1]] + Amps.S1.evec_ind[2*ind + 1]];
      tempt /= tempen;
      tempt = mix*tempt + (1.0-mix)*Amps.S1.t2[ind];
      Amps.S1.set_T(ind, tempt);
    }
  }
}

void CC_Error(Input_Parameters &Parameters, Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error)
{
  error = 0.0;
  int ind;
  int length = Amps.D1.T1_length;
  double norm = 1.0e-16;
  for(ind = 0; ind < length; ++ind){
    error += (Amps2.D1.T1[ind] - Amps.D1.T1[ind])*(Amps2.D1.T1[ind] - Amps.D1.T1[ind]);
    norm += Amps2.D1.T1[ind] * Amps2.D1.T1[ind];
  }
  if(Parameters.approx == "singles"){  
    length = Amps.S1.t2_length;
    for(ind = 0; ind < length; ++ind){
      error += (Amps2.S1.t2[ind] - Amps.S1.t2[ind])*(Amps2.S1.t2[ind] - Amps.S1.t2[ind]);
      norm += Amps2.S1.t2[ind] * Amps2.S1.t2[ind];
    }
  }
  error = std::sqrt(error/norm);
}

void Print_Amps(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps)
{
  int length1 = Chan.size1;
  int length2 = Chan.nph1[Chan.ind0];
  int npp, nhh;
  int chanind, ind, index;
  int a, b, i, j;
  two_body pp, hh, ph;
  ind = 0;

  std::cout << std::endl;
  for(int chan1 = 0; chan1 < length1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    chanind = Amps.D1.T1_index[chan1];
    for(int pp1 = 0; pp1 < npp; ++pp1){
      pp = Chan.pp_state(chan1, pp1);
      a = pp.v1;
      b = pp.v2;
      for(int hh1 = 0; hh1 < nhh; ++hh1){
	hh = Chan.hh_state(chan1, hh1);
	i = hh.v1;
	j = hh.v2;
	std::cout << "! < " << a << "," << b << " |t| " << i << "," << j << " >^ " << Chan.qnums1[chan1].j << " = " << chanind << " " << (pp1*nhh + hh1) << ", " << Amps.D1.T1[chanind + (pp1*nhh + hh1)] << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind]; ++n){
	  index = Amps.D1.map_index[8*ind] + n;
	  std::cout << " " << Amps.D1.T2_1[Amps.D1.T2_1_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 1]; ++n){
	  index = Amps.D1.map_index[8*ind + 1] + n;
	  std::cout << " " << Amps.D1.T2_2[Amps.D1.T2_2_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 2]; ++n){
	  index = Amps.D1.map_index[8*ind + 2] + n;
	  std::cout << " " << Amps.D1.T2_3[Amps.D1.T2_3_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 3]; ++n){
	  index = Amps.D1.map_index[8*ind + 3] + n;
	  std::cout << " " << Amps.D1.T2_4[Amps.D1.T2_4_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 4]; ++n){
	  index = Amps.D1.map_index[8*ind + 4] + n;
	  std::cout << " " << Amps.D1.T3_1[Amps.D1.T3_1_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 5]; ++n){
	  index = Amps.D1.map_index[8*ind + 5] + n;
	  std::cout << " " << Amps.D1.T3_2[Amps.D1.T3_2_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 6]; ++n){
	  index = Amps.D1.map_index[8*ind + 6] + n;
	  std::cout << " " << Amps.D1.T3_3[Amps.D1.T3_3_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << ",";
	for(int n = 0; n < Amps.D1.map_num[8*ind + 7]; ++n){
	  index = Amps.D1.map_index[8*ind + 7] + n;
	  std::cout << " " << Amps.D1.T3_4[Amps.D1.T3_4_index[Amps.D1.map_chan[index]] + Amps.D1.map_ind[index]];
	}
	std::cout << std::endl;
	++ind;
      }
    }
  }
  if(Parameters.approx == "singles"){
    for(int ph1 = 0; ph1 < length2; ++ph1){
      ph = Chan.ph1_state(Chan.ind0, ph1);
      a = ph.v1;
      i = ph.v2;
      std::cout << "! < " << a << " |t| " << i << " > = " << Amps.S1.t2[ph1] << " " << Amps.S1.t3[Amps.S1.t3_index[Amps.S1.map_chan[ph1]] + Amps.S1.map_ind[ph1]] << std::endl;
    }
  }
  std::cout << std::endl;
}

void Random_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, double &width, double &error, double &error2)
{
  int go = 0;
  int ind = 0;
  while(go == 0){
    ++ind;
    Amps.zero(Parameters, Chan, false);
    Randomize_Amps(Parameters, Chan, Ints, Amps0, Amps, width);

    Eff_Ints.zero(Parameters, Chan, false);
    tempAmps.zero(Parameters, Chan, false);
    //  Build Effective Hamiltonian
    Update_Heff_1(Parameters, Space, Chan, Ints, Eff_Ints, Amps);
    if(Parameters.approx == "singles"){
      Update_Heff_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps);
    }
    //  Solve Amplitude Equations
    Eff_Ints.diagonalize(Parameters, Chan);
    Doubles_Step(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
    if(Parameters.approx == "singles"){
      Singles_Step(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
      Doubles_Step_2(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
    }
    Amps2.copy_Amplitudes(Parameters, Chan, Amps);
    Amps2.zero(Parameters, Chan, true);
    Gather_Amps(Parameters, Space, Chan, Eff_Ints, Amps2, tempAmps, mix);
    CC_Error(Parameters, Chan, Ints, Amps, Amps2, error);
    error /= mix;

    std::cout << "!!  " << ind << ", " << error << " " << error2 << ", " << width << std::endl;
    if(error < error2){
      ++go;
      Amps0.copy_Amplitudes(Parameters, Chan, Amps);
      Amps.copy_Amplitudes(Parameters, Chan, Amps2);
      error2 = error;
    }
    width *= 0.975;
    if(width < 10-12){ width = 10e-12; }
    if(ind >= 1000){
      std::cout << Parameters.Shells << ", " << Parameters.Pshells << ", " << Parameters.density << std::endl;
      std::cerr << std::endl << "Random Step Unsuccessful!!" << std::endl; exit(1);
    }
  }
}

void Randomize_Amps(Input_Parameters &Parameters, Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, double &width)
{
  int ind;
  double tempt;
  double rand;
  size_t key;
  std::unordered_map<size_t,double> t_map;
  int length = Amps.D1.T1_length;
  for(ind = 0; ind < length; ++ind){
    tempt = Amps0.D1.T1[ind];
    if(fabs(tempt) > 1.0e-12){
      rand = rand_normal(0.0, width * fabs(tempt));
      key = std::hash<float>{}(float(fabs(tempt)));
      t_map[key] = rand;
    }
  }
  if(Parameters.approx == "singles"){  
    length = Amps.S1.t2_length;
    for(ind = 0; ind < length; ++ind){
      tempt = Amps0.S1.t2[ind];
      if(fabs(tempt) > 1.0e-12){
	rand = rand_normal(0.0, width * fabs(tempt));
	key = std::hash<float>{}(float(fabs(tempt)));
	t_map[key] = rand;
      }
    }
  }

  length = Amps.D1.T1_length;
  for(ind = 0; ind < length; ++ind){
    tempt = Amps0.D1.T1[ind];
    key = std::hash<float>{}(float(fabs(tempt)));
    if(tempt > 1.0e-12){ tempt += t_map[key]; }
    else if(tempt < -1.0e-12){ tempt -= t_map[key]; }
    Amps.D1.set_T(ind, tempt);
  }
  if(Parameters.approx == "singles"){  
    length = Amps.S1.t2_length;
    for(ind = 0; ind < length; ++ind){
      tempt = Amps0.S1.t2[ind];
      key = std::hash<float>{}(float(fabs(tempt)));
      if(tempt > 1.0e-12){ tempt += t_map[key]; }
      else if(tempt < -1.0e-12){ tempt -= t_map[key]; }
      Amps.S1.set_T(ind, tempt);
    }
  }
}

void Initialize_DIIS(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl)
{
  int size = Amps.D1.T1_length;
  if(Parameters.approx == "singles"){ size += Amps.S1.t2_length; }
  B = new double[1];
  B[0] = 0.0;

  p = new double[maxl * size];
  delp = new double[maxl * size];
  tempdelp = new double[size];
  for(int l = 0; l < maxl; ++l){
    for(int ind = 0; ind < size; ++ind){
      p[l * ind] = 0.0;
      delp[l * ind] = 0.0;
    }
  }
  for(int ind = 0; ind < size; ++ind){
    tempdelp[ind] = 0.0;
  }
}

void Delete_DIIS(Input_Parameters &Parameters, Channels &Chan, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl)
{
  delete[] B;
  delete[] p;
  delete[] delp;
  delete[] tempdelp;
}

void Perform_DIIS(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, Amplitudes &Amps0, double &mix, double *&p, double *&delp, double *&tempdelp, double *&B, int &N, int &maxl, int &DIIS_count)
{
  double norm = 0.0;
  double norm0 = 0.0;
  double tempt;
  int length1, length2, size, ind;
  int ortho;
  int P = N + 1;
  int lwork;
  int *ipiv;
  double *work;
  int info;
  double *B2;
  double dot_prod, norm1, norm2;  
  length1 = Amps.D1.T1_length;
  length2 = Amps.S1.t2_length;
  size = length1 + length2;
  
  // Fill temp_delp
  for(ind = 0; ind < length1; ++ind){
    tempdelp[ind] = (Amps.D1.T1[ind] - Amps0.D1.T1[ind]) / mix;
    norm0 += Amps0.D1.T1[ind] * Amps0.D1.T1[ind];
  }
  if(Parameters.approx == "singles"){
    for(ind = 0; ind < length2; ++ind){
      tempdelp[length1 + ind] = (Amps.S1.t2[ind] - Amps0.S1.t2[ind]) / mix;
      norm0 += Amps0.S1.t2[ind] * Amps0.S1.t2[ind];
    }	
  }
  norm0 = std::sqrt(norm0);
  for(ind = 0; ind < length1; ++ind){
    tempdelp[ind] /= norm0;
    norm += tempdelp[ind] * tempdelp[ind];
  }
  if(Parameters.approx == "singles"){
    for(ind = 0; ind < length2; ++ind){
      tempdelp[length1 + ind] /= norm0;
      norm += tempdelp[length1 + ind] * tempdelp[length1 + ind];
    }	
  }

  // check orthogonality of tempdelp
  ortho = 1;
  for(int l = 0; l < N; ++l){
    dot_prod = 0.0;
    norm1 = 0.0;
    norm2 = 0.0;
    for(ind = 0; ind < length1; ++ind){
      dot_prod += Amps.D1.T1[ind] * p[l * size + ind];
      norm1 += Amps.D1.T1[ind] * Amps.D1.T1[ind];
      norm2 += p[l * size + ind] * p[l * size + ind];
    }
    if(Parameters.approx == "singles"){
      for(ind = 0; ind < length2; ++ind){
	dot_prod += Amps.S1.t2[ind] * p[l * size + (length1 + ind)];
	norm1 += Amps.S1.t2[ind] * Amps.S1.t2[ind];
	norm2 += p[l * size + (length1 + ind)] * p[l * size + (length1 + ind)];
      }	
    }
    dot_prod /= std::sqrt(norm1 * norm2);
    //std::cout << "!! norm, dot_prod, B[P*l+l] " << norm << " " << dot_prod << " " << B[P*l + l] << std::endl;
    if(norm > B[P * l + l] || dot_prod > (1.0 - std::pow(norm, 1.5))){ ortho = 0; break; }
  }
  if(ortho != 1){ return; }

  ++DIIS_count;
  if(N < maxl){ Update_B1(Parameters, Chan, Amps, N, p, delp, tempdelp, B); }
  else{ Update_B2(Parameters, Chan, Amps, N, p, delp, tempdelp, B); }
  P = N + 1;

  // copy B into B2 for inversion
  B2 = new double[P * P];
  for(int j = 0; j < P*P; ++j){ B2[j] = B[j]; }
  
  ipiv = new int[P];
  work = new double[sizeof(double) * P];
  lwork = sizeof(double) * P;
  info = 0;
  dgetrf_(&P, &P, B2, &P, ipiv, &info);
  dgetri_(&P, B2, &P, ipiv, work, &lwork, &info);
  Amps.zero(Parameters, Chan, false);
  for(ind = 0; ind < length1; ++ind){
    tempt = 0.0;
    for(int l = 0; l < N; ++l){ tempt += -1.0 * B2[P * l + N] * p[l * size + ind]; }
    Amps.D1.set_T(ind, tempt);
  }
  if(Parameters.approx == "singles"){
    for(ind = 0; ind < length2; ++ind){
      tempt = 0.0;
      for(int l = 0; l < N; ++l){ tempt += -1.0 * B2[P * l + N] * p[l * size + (length1 + ind)]; }
      Amps.S1.set_T(ind, tempt);
    }
  }
  delete[] B2;
  delete[] ipiv;
  delete[] work;
}

void Update_B1(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B)
{
  double *B2; // used to copy B into updated B
  int P = N+1;
  int ind, length1, length2, size;
  length1 = Amps.D1.T1_length;
  length2 = Amps.S1.t2_length;
  size = length1 + length2;

  for(ind = 0; ind < length1; ++ind){
    p[N * size + ind] = Amps.D1.T1[ind];
    delp[N * size + ind] = tempdelp[ind];
  }
  if(Parameters.approx == "singles"){
    for(ind = 0; ind < length2; ++ind){
      p[N * size + (length1 + ind)] = Amps.S1.t2[ind];
      delp[N * size + (length1 + ind)] = tempdelp[length1 + ind];
    }
  }

  B2 = new double[(P+1) * (P+1)];
  for(int j = 0; j < N; ++j){
    for(int k = 0; k < N; ++k){
      B2[(P+1) * j + k] = B[P * j + k];
    }
  }
  for(int l = 0; l < P; ++l){
    B2[(P+1) * N + l] = 0.0;
    if(l != N){ B2[(P+1) * l + N] = 0.0; }
    for(ind = 0; ind < length1; ++ind){
      B2[(P+1) * N + l] += delp[N * size + ind] * delp[l * size + ind];
      if(l != N){ B2[(P+1) * l + N] += delp[l * size + ind] * delp[N * size + ind]; }
    }
    if(Parameters.approx == "singles"){
      for(ind = 0; ind < length2; ++ind){
	B2[(P+1) * N + l] += delp[N * size + (length1 + ind)] * delp[l * size + (length1 + ind)];
	if(l != N){ B2[(P+1) * l + N] += delp[l * size + (length1 + ind)] * delp[N * size + (length1 + ind)]; }
      }
    }
  }
  for(int l = 0; l < P; ++l){
    B2[(P+1) * P + l] = -1.0;
    B2[(P+1) * l + P] = -1.0;
  }
  B2[(P+1) * P + P] = 0.0;
  ++N;
  ++P;
  delete[] B;
  B = new double[P * P];
  for(int ind = 0; ind < P * P; ++ind){ B[ind] = B2[ind]; }
  delete[] B2;
}

void Update_B2(Input_Parameters &Parameters, Channels &Chan, Amplitudes &Amps, int &N, double *&p, double *&delp, double *&tempdelp, double *&B)
{
  int maxind, ind, length1, length2, size;
  double maxnorm;
  int P = N + 1;
  length1 = Amps.D1.T1_length;
  length2 = Amps.S1.t2_length;
  size = length1 + length2;

  // Find largest norm of B to remove that vector
  maxind = -1;
  maxnorm = 0.0;
  for(int j = 0; j < N; ++j){
    if(fabs(B[P * j + j]) > maxnorm){
      maxind = j;
      maxnorm = fabs(B[P * j + j]);
    }
  }
  // Replace maxnorm vector with new vector
  for(ind = 0; ind < length1; ++ind){
    p[maxind * size + ind] = Amps.D1.T1[ind];
    delp[maxind * size + ind] = tempdelp[ind];
  }
  if(Parameters.approx == "singles"){
    for(ind = 0; ind < length2; ++ind){
      p[maxind * size + (length1 + ind)] = Amps.S1.t2[ind];
      delp[maxind * size + (length1 + ind)] = tempdelp[length1 + ind];
    }
  }
  for(int l = 0; l < N; ++l){
    B[P * maxind + l] = 0.0;
    if(l != maxind){ B[P * l + maxind] = 0.0; }	
    for(ind = 0; ind < length1; ++ind){
      B[P * maxind + l] += delp[maxind * size + ind] * delp[l * size + ind];
      if(l != maxind){ B[P * l + maxind] += delp[l * size + ind] * delp[maxind * size + ind]; }
    }
    if(Parameters.approx == "singles"){
      for(ind = 0; ind < length2; ++ind){
	B[P * maxind + l] += delp[maxind * size + (length1 + ind)] * delp[l * size + (length1 + ind)];
	if(l != maxind){ B[P * l + maxind] += delp[l * size + (length1 + ind)] * delp[maxind * size + (length1 + ind)]; }
      }
    }
  }
}

void Update_Heff_1(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m12 = -0.5, fac;
  int nh, np, npph, nhhp, nhh, npp, nhp1, nph1, nppp, nhhh, length;
  int p_ind, h_ind;
  int chan_ind1, chan_ind2, chan_ind3;
  char N = 'N';
  
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];
    //////////////   Xpp   ///////////////
    if(Parameters.HF == 1){
      chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
      // Xpp3(a|b){a,b}  <-  fpp(a|a) del_(ab)
      for(int p = 0; p < np; ++p){
 	p_ind = Chan.p_state(chan3, p).v1;
 	Eff_Ints.Xpp.X_3[chan_ind1 + (p*np + p)] += Space.qnums[p_ind].energy;
      }
    }
    else if(Parameters.HF == 0){
      chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
      chan_ind2 = Ints.Fmatrix.pp_3_index[chan3];
      // Xpp3(a|b){a,b}  <-  fpp(a|b){a,b}
      for(int p1 = 0; p1 < np; ++p1){
	for(int p2 = 0; p2 < np; ++p2){
	  Eff_Ints.Xpp.X_3[chan_ind1 + (p1*np + p2)] += Ints.Fmatrix.pp_3[chan_ind2 + (p1*np + p2)];
	}
      }
    }
    if(np * nhhp != 0){
      chan_ind1 = Amps.D1.T3_1_index[chan3];
      chan_ind2 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind3 = Eff_Ints.Xpp.X_3_index[chan3];
      // Xpp3(a|b){a,b}  =  -(1/2) T3_1(ac|kl){a,klc'}.Vhhpp3_3(kl|bc){klc',b}
      dgemm_NN((Amps.D1.T3_1 + chan_ind1), (Ints.Vhhpp.V_3_3 + chan_ind2), (Eff_Ints.Xpp.X_3 + chan_ind3), &np, &np, &nhhp, &m12, &p1, &N, &N);
    }

    //////////////   X1hh,Xhh   ///////////////
    if(Parameters.HF == 1){
      chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
      // X1hh3(i|j){i,j}  <-  fhh(i|i) del_(ij)
      for(int h = 0; h < nh; ++h){
	h_ind = Chan.h_state(chan3, h).v1;
	Eff_Ints.Xhh.X1_3[chan_ind1 + (h*nh + h)] += Space.qnums[h_ind].energy;
	Eff_Ints.Xhh.X_3[chan_ind1 + (h*nh + h)] += Space.qnums[h_ind].energy;
      }
    }
    else if(Parameters.HF == 0){
      chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
      chan_ind2 = Ints.Fmatrix.hh_3_index[chan3];
      // X1hh3(i|j){i,j}  <-  fhh(i|j){i,j
      for(int h1 = 0; h1 < nh; ++h1){
	for(int h2 = 0; h2 < nh; ++h2){
	  Eff_Ints.Xhh.X1_3[chan_ind1 + (h1*nh + h2)] += Ints.Fmatrix.hh_3[chan_ind2 + (h1*nh + h2)];
	  Eff_Ints.Xhh.X_3[chan_ind1 + (h1*nh + h2)] += Ints.Fmatrix.hh_3[chan_ind2 + (h1*nh + h2)];
	}
      }
    }
    if(nh * npph != 0){
      chan_ind1 = Ints.Vhhpp.V_3_1_index[chan3];
      chan_ind2 = Amps.D1.T3_3_index[chan3];
      chan_ind3 = Eff_Ints.Xhh.X_3_index[chan3];
      // X1hh3(i|j){i,j}  =  +(1/2) Vhhpp3_1(ik|cd){i,cdk'}.T3_3(cd|jk){cdk',j}
      dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), (Amps.D1.T3_3 + chan_ind2), (Eff_Ints.Xhh.X1_3 + chan_ind3), &nh, &nh, &npph, &p12, &p1, &N, &N);
      // Xhh3(i|j){i,j}  =  +(1/2) Vhhpp3_1(ik|cd){i,cdk'}.T3_3(cd|jk){cdk',j}
      dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), (Amps.D1.T3_3 + chan_ind2), (Eff_Ints.Xhh.X_3 + chan_ind3), &nh, &nh, &npph, &p12, &p1, &N, &N);
    }
  }

  //////////////   X1hhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    nhhp = Chan.nhhp[chan3];
    length = nhhp * nh;
    chan_ind1 = Eff_Ints.Xhhhp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vhhhp.V_3_3_index[chan3];
    // X1hhhp3_3 = Vhhhp3_3
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhhhp.X1_3_3[chan_ind1 + ind] += Ints.Vhhhp.V_3_3[chan_ind2 + ind];
    }
  }

  //////////////   X1hppp,Xhppp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    length = np * npph;
    chan_ind1 = Eff_Ints.Xhppp.X_3_2_index[chan3];
    chan_ind2 = Ints.Vhppp.V_3_2_index[chan3];
    // X1hppp3_2 = Vhppp3_2, Xhppp3_2 = Vhppp3_2
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhppp.X1_3_2[chan_ind1 + ind] += Ints.Vhppp.V_3_2[chan_ind2 + ind];
      Eff_Ints.Xhppp.X_3_2[chan_ind1 + ind] += Ints.Vhppp.V_3_2[chan_ind2 + ind];
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];

    //////////////   Xhhhh   ///////////////
    length = nhh * nhh;
    chan_ind1 = Eff_Ints.Xhhhh.X_1_index[chan1];
    chan_ind2 = Ints.Vhhhh.V_1_index[chan1];
    // Xhhhh1(ij|kl){ij,kl}  =  Vhhhh1(ij|kl){ij,kl}
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhhhh.X_1[chan_ind1 + ind] += Ints.Vhhhh.V_1[chan_ind2 + ind];
    }
    if(nhh * npp != 0){
      chan_ind1 = Ints.Vhhpp.V_1_index[chan1];
      chan_ind2 = Amps.D1.T1_index[chan1];
      chan_ind3 = Eff_Ints.Xhhhh.X_1_index[chan1];
      // Xhhhh1(ij|kl){ij,kl}  <-  +(1/2) Vhhpp1(ij|cd){ij,cd}.T1(cd|kl){cd,kl}
      dgemm_NN((Ints.Vhhpp.V_1 + chan_ind1), (Amps.D1.T1 + chan_ind2), (Eff_Ints.Xhhhh.X_1 + chan_ind3), &nhh, &nhh, &npp, &p12, &p1, &N, &N);
    }

    //////////////   X1pppp   ///////////////
    length = npp * npp;
    chan_ind1 = Eff_Ints.Xpppp.X_1_index[chan1];
    chan_ind2 = Ints.Vpppp.V_1_index[chan1];
    // X1pppp1(ab|cd){ab,cd}  =  Vpppp1(ab|cd){ab,cd}
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xpppp.X1_1[chan_ind1 + ind] += Ints.Vpppp.V_1[chan_ind2 + ind];
    }
  }

  //////////////   X1hphp,X2hphp,Xhphp   ///////////////
  if(Parameters.basis != "finite_J"){ fac = m12; }
  else{ fac = p12; }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    length = nhp1 * nhp1;
    chan_ind1 = Eff_Ints.Xhphp.X_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    // X(1)(2)hphp2_1(ia|jb){ib',ja'}  =  Vhphp2_1(ia|jb){ib',ja'}
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhphp.X1_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
      Eff_Ints.Xhphp.X2_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
      Eff_Ints.Xhphp.X_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
    }
    if(nhp1 * nph1 != 0){
      chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
      chan_ind2 = Amps.D1.T2_1_index[chan2];
      chan_ind3 = Eff_Ints.Xhphp.X_2_1_index[chan2];
      // Xhphp2_1(ia|jb){ib',ja'}  <-  -(+)(1/2) Vhhpp2_1(ik|cb){ib',ck'}.T2_1(ca|jk){ck',ja'}
      dgemm_NN((Ints.Vhhpp.V_2_1 + chan_ind1), (Amps.D1.T2_1 + chan_ind2), (Eff_Ints.Xhphp.X_2_1 + chan_ind3), &nhp1, &nhp1, &nph1, &fac, &p1, &N, &N);
    }
  }

  //////////// Xpphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    nppp = Chan.nppp[chan3];
    length = nppp * nh;
    chan_ind1 = Eff_Ints.Xpphp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vpphp.V_3_3_index[chan3];
    // X1pphp3_3 = Vpphp3_3
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xpphp.X1_3_3[chan_ind1 + ind] += Ints.Vpphp.V_3_3[chan_ind2 + ind];
    }
  }

  //////////// Xhphh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    length = np * nhhh;
    chan_ind1 = Eff_Ints.Xhphh.X_3_2_index[chan3];
    chan_ind2 = Ints.Vhphh.V_3_2_index[chan3];
    // X1hphh3_2 = Vhphh3_2
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhphh.X1_3_2[chan_ind1 + ind] += Ints.Vhphh.V_3_2[chan_ind2 + ind];
    }
  }
}

void Update_Heff_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5, fac;
  int nh, np, npph, nhhp, nppp, nhhh, nhpp, nhph, one = 1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  int nhp0 = Chan.nhp1[Chan.ind0];
  int nph0 = Chan.nph1[Chan.ind0];
  int npp0 = Chan.npp1[Chan.ind0];
  int nhh0 = Chan.nhh1[Chan.ind0];
  char N = 'N';

  ////////////   Xhp   ///////////////
  if(nhp0 != 0){
    if(Parameters.HF == 0){
      // Xhp2(i|a){ia'}  <-  fhp(i|a){ia'}
      for(int hp1 = 0; hp1 < nhp0; ++hp1){
	Eff_Ints.Xhp.X_2[hp1] += Ints.Fmatrix.hp_2[hp1];
      }
    }
    chan_ind1 = Ints.Vhhpp.V_2_3_index[Chan.ind0];
    // Xhp2(i|a){ia'}  =  Vhhpp2_3(ik|ac){ia',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhpp.V_2_3 + chan_ind1), Amps.S1.t2, Eff_Ints.Xhp.X_2, &nhp0, &one, &nph0, &p1, &p1, &N, &N);
  }
  Eff_Ints.Xhp.X_gather(Parameters, Chan);
 
  //////////////   Xpp   ///////////////
  if(nph0 * npp0 != 0){
    chan_ind1 = Ints.Vhppp.V_2_4_index[Chan.ind0];
    // Xpp2(a|b){ab'}  =  Vhppp2_4(ka|cb){ab',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhppp.V_2_4 + chan_ind1), Amps.S1.t2, Eff_Ints.Xpp.X_2, &npp0, &one, &nph0, &p1, &p1, &N, &N);
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    if(nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Eff_Ints.Xhp.X_3_index[chan3];
      chan_ind3 = Eff_Ints.Xpp.X_3_index[chan3];
      // Xpp3(a|b){a,b}  <-  - t3(a|k){a,k}.Xhp3(k|b){k,b}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhp.X_3 + chan_ind2), (Eff_Ints.Xpp.X_3 + chan_ind3), &np, &np, &nh, &m1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xpp.X_gather(Parameters, Chan);

  //////////////   Xhh   ///////////////
  if(nph0 * nhh0 != 0){
    chan_ind1 = Ints.Vhhhp.V_2_3_index[Chan.ind0];
    // X1hh2(i|j){ij'}  =  Vhhhp2_3(ik|jc){ij',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhhp.V_2_3 + chan_ind1), Amps.S1.t2, Eff_Ints.Xhh.X1_2, &nhh0, &one, &nph0, &p1, &p1, &N, &N);
    // Xhh2(i|j){ij'}  =  Vhhhp2_3(ik|jc){ij',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhhp.V_2_3 + chan_ind1), Amps.S1.t2, Eff_Ints.Xhh.X_2, &nhh0, &one, &nph0, &p1, &p1, &N, &N);
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    if(nh * np != 0){
      chan_ind1 = Eff_Ints.Xhp.X_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhh.X_3_index[chan3];
      // Xhh3(i|j){i,j}  <-  Xhp3(i|c){i,c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhp.X_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhh.X_3 + chan_ind3), &nh, &nh, &np, &p1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhh.X1_gather(Parameters, Chan);
  Eff_Ints.Xhh.X_gather(Parameters, Chan);

  //////////////   X1hppp,Xhppp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    if(nh * np * npph != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhpp.V_3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhppp.X_3_2_index[chan3];
      // X1hppp3_2(ia|bc){a,bci'}  <-  - (1/2) t3(a|k){a,k}.Vhhpp3_2(ik|bc){k,bci'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhpp.V_3_2 + chan_ind2), (Eff_Ints.Xhppp.X1_3_2 + chan_ind3), &np, &npph, &nh, &m12, &p1, &N, &N);
      // Xhppp3_2(ia|bc){a,bci'}  <-  - t3(a|k){a,k}.Vhhpp3_2(ik|bc){k,bci'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhpp.V_3_2 + chan_ind2), (Eff_Ints.Xhppp.X_3_2 + chan_ind3), &np, &npph, &nh, &m1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhppp.X1_gather(Parameters, Chan);
  Eff_Ints.Xhppp.X_gather(Parameters, Chan);

  //////////////   Xhhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    if(nh * np * nhhp != 0){
      chan_ind1 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhhhp.X_3_3_index[chan3];
      // X1hhhp3_3(ik|ja){ika',j}  <-  +(1/2) Vhhpp3_3(ik|ca){ika',c}.t3(c|j){c,j}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhhhp.X1_3_3 + chan_ind3), &nhhp, &nh, &np, &p12, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhhhp.X1_gather(Parameters, Chan);

  //////////////   X1pppp   ///////////////
  if(Parameters.basis != "finite_J"){ fac = p1; }
  else{ fac = m1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    if(nh * np * nppp){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Eff_Ints.Xhppp.X_3_1_index[chan3];
      chan_ind3 = Eff_Ints.Xpppp.X_3_1_index[chan3];
      chan_ind4 = Eff_Ints.Xpppp.X_3_2_index[chan3];
      // X1pppp3_1(ab|cd){a,cdb'}  =  - t3(a|k){a,k}.X1hppp3_1(kb|cd){k,cdb'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhppp.X1_3_1 + chan_ind2), (Eff_Ints.Xpppp.X1_3_1 + chan_ind3), &np, &nppp, &nh, &m1, &p1, &N, &N);
      // X1pppp3_2(ab|cd){b,cda'}  =  +(-) t3(b|k){b,k}.X1hppp3_1(ka|cd){k,cda'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhppp.X1_3_1 + chan_ind2), (Eff_Ints.Xpppp.X1_3_2 + chan_ind4), &np, &nppp, &nh, &fac, &p1, &N, &N);
    }
  }
  Eff_Ints.Xpppp.X1_gather(Parameters, Chan);

  //////////////   Xhhhh   ///////////////
  if(Parameters.basis != "finite_J"){ fac = m1; }
  else{ fac = p1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    if(nh * np * nhhh){
      chan_ind1 = Eff_Ints.Xhhhp.X_3_4_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhhhh.X_3_4_index[chan3];
      chan_ind4 = Eff_Ints.Xhhhh.X_3_3_index[chan3];
      // Xhhhh3_4(ij|kl){ijk',l}  =  + X1hhhp3_4(ij|kc){ijk',c}.t3(c|l){c,l}
      dgemm_NN((Eff_Ints.Xhhhp.X1_3_4 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhhhh.X_3_4 + chan_ind3), &nhhh, &nh, &np, &p1, &p1, &N, &N);
      // Xhhhh3_3(ij|kl){ijk',l}  =  -(+) X1hhhp3_4(ij|kc){ijk',c}.t3(c|l){c,l}
      dgemm_NN((Eff_Ints.Xhhhp.X1_3_4 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhhhh.X_3_3 + chan_ind4), &nhhh, &nh, &np, &fac, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhhhh.X_gather(Parameters, Chan);

  //////////////   X1hphp,X2hphp,Xhphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhpp = Chan.nhpp[chan3];
    nhph = Chan.nhph[chan3];
    if(nhpp * nh * np != 0){
      chan_ind1 = Eff_Ints.Xhppp.X_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhphp.X_3_3_index[chan3];
      // X1hphp3_3(ia|jb){iab',j}  =  X1hppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhppp.X1_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphp.X1_3_3 + chan_ind3), &nhpp, &nh, &np, &p1, &p1, &N, &N);
      // X2hphp3_3(ia|jb){iab',j}  =  +(1/2) X1hppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhppp.X1_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphp.X2_3_3 + chan_ind3), &nhpp, &nh, &np, &p12, &p1, &N, &N);
      // Xhphp3_3(ia|jb){iab',j}  =  Xhppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhppp.X_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphp.X_3_3 + chan_ind3), &nhpp, &nh, &np, &p1, &p1, &N, &N);
    }
    if(nhph * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhhp.V_3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhphp.X_3_2_index[chan3];
      // X1hphp3_2(ia|jb){a,jbi'}  =  -(1/2) t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhp.V_3_2 + chan_ind2), (Eff_Ints.Xhphp.X1_3_2 + chan_ind3), &np, &nhph, &nh, &m12, &p1, &N, &N);
      // X2hphp3_2(ia|jb){a,jbi'}  =  -(1/2) t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhp.V_3_2 + chan_ind2), (Eff_Ints.Xhphp.X2_3_2 + chan_ind3), &np, &nhph, &nh, &m12, &p1, &N, &N);
      // Xhphp3_2(ia|jb){a,jbi'}  =  - t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhp.V_3_2 + chan_ind2), (Eff_Ints.Xhphp.X_3_2 + chan_ind3), &np, &nhph, &nh, &m1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhphp.X1_gather(Parameters, Chan);
  Eff_Ints.Xhphp.X2_gather(Parameters, Chan);
  Eff_Ints.Xhphp.X_gather(Parameters, Chan);

  //////////////   X1pphp   ///////////////
  if(Parameters.basis != "finite_J"){ fac = p1; }
  else{ fac = m1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    nhpp = Chan.nhpp[chan3];
    if(nppp * nh * np != 0){
      chan_ind1 = Ints.Vpppp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_3_index[chan3];
      // X1pphp3_3(ab|ic){abc',i}  <-  +(1/2) Vpppp3_3(ab|dc){abc',d}.t3(d|i){d,i}
      dgemm_NN((Ints.Vpppp.V_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xpphp.X1_3_3 + chan_ind3), &nppp, &nh, &np, &p12, &p1, &N, &N);
    }
    if(nhpp * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Eff_Ints.Xhphp.X_3_1_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_1_index[chan3];
      chan_ind4 = Eff_Ints.Xpphp.X_3_2_index[chan3];
      // X1pphp3_1(ab|ic){a,icb'}  =  - t3(a|k){a,l}.X2hphp3_1(kb|ic){k,icb'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhphp.X2_3_1 + chan_ind2), (Eff_Ints.Xpphp.X1_3_1 + chan_ind3), &np, &nhpp, &nh, &m1, &p1, &N, &N);
      // X1pphp3_2(ab|ic){b,ica'}  =  +(-) t3(b|k){b,l}.X2hphp3_1(ka|ic){k,ica'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhphp.X2_3_1 + chan_ind2), (Eff_Ints.Xpphp.X1_3_2 + chan_ind4), &np, &nhpp, &nh, &fac, &p1, &N, &N);
    }
  }
  Eff_Ints.Xpphp.X1_gather(Parameters, Chan);

  //////////////   X1hphh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    if(nhhh * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhhh.V_3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhphh.X_3_2_index[chan3];
      // X1hphh3_2(ia|jk){a,jki'}  <-  -(1/2) t3(a|l){a,l}.Vhhhh3_2(il|jk){l,jki'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhh.V_3_2 + chan_ind2), (Eff_Ints.Xhphh.X1_3_2 + chan_ind3), &np, &nhhh, &nh, &m12, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhphh.X1_gather(Parameters, Chan);
}

void Doubles_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, fac1, fac2;
  int nh, np, npph, nhhp, nhh, npp, nhp1, nph1, length;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4, chan_ind5;
  char N = 'N';
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    length = npp * nhh;
    chan_ind1 = Amps1.D1.T1_index[chan1];
    chan_ind2 = Ints.Vpphh.V_1_index[chan1];
    //T1(ab|ij){ab,ij}  =  Vpphh1(ab|ij){ab,ij}
    for(int ind = 0; ind < length; ++ind){ Amps2.D1.T1[chan_ind1 + ind] += Ints.Vpphh.V_1[chan_ind2 + ind]; }

    if(length != 0){
      chan_ind2 = Eff_Ints.Xpppp.X_1_index[chan1];
      chan_ind3 = Eff_Ints.Xhhhh.X_1_index[chan1];
      //T1(ab|ij){ab,ij}  <-  (1/2).X1pppp1(ab|cd){ab,cd}.T1(cd|ij){cd,ij}
      dgemm_NN((Eff_Ints.Xpppp.X1_1 + chan_ind2), (Amps1.D1.T1 + chan_ind1), (Amps2.D1.T1 + chan_ind1), &npp, &nhh, &npp, &p12, &p1, &N, &N);
      //T1(ab|ij){ab,ij}  <-  (1/2).T1(ab|kl){ab,kl}.Xhhhh1(kl|ij){kl,ij}
      dgemm_NN((Amps1.D1.T1 + chan_ind1), (Eff_Ints.Xhhhh.X_1 + chan_ind3), (Amps2.D1.T1 + chan_ind1), &npp, &nhh, &nhh, &p12, &p1, &N, &N);
    }
  }
  if(Parameters.basis != "finite_J"){ fac1 = m1; }
  else{ fac1 = p1; }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    nhp1 = Chan.nhp1[chan2];
    length = nph1 * nhp1;
    if(length != 0){
      chan_ind1 = Amps1.D1.T2_1_index[chan2];
      chan_ind2 = Amps1.D1.T2_2_index[chan2];
      chan_ind3 = Amps1.D1.T2_3_index[chan2];
      chan_ind4 = Amps1.D1.T2_4_index[chan2];
      chan_ind5 = Eff_Ints.Xhphp.X_2_1_index[chan2];
      //T2_1(ab|ij){aj',ib'}  =  -(+) T2_1(ac|kj){aj',kc'}.Xhphp2_1(kb|ic){kc',ib'}
      dgemm_NN((Amps1.D1.T2_1 + chan_ind1), (Eff_Ints.Xhphp.X_2_1 + chan_ind5), (Amps2.D1.T2_1 + chan_ind1), &nph1, &nhp1, &nhp1, &fac1, &p1, &N, &N);
      //T2_2(ab|ij){bi',ja'}  =  -(+) T2_1(bc|ki){bi',kc'}.Xhphp2_1(ka|jc){kc',ja'}
      dgemm_NN((Amps1.D1.T2_1 + chan_ind1), (Eff_Ints.Xhphp.X_2_1 + chan_ind5), (Amps2.D1.T2_2 + chan_ind2), &nph1, &nhp1, &nhp1, &fac1, &p1, &N, &N);
      //T2_3(ab|ij){ai',jb'}  =   T2_1(ac|ki){ai',kc'}.Xhphp2_1(kb|jc){kc',jb'}
      dgemm_NN((Amps1.D1.T2_1 + chan_ind1), (Eff_Ints.Xhphp.X_2_1 + chan_ind5), (Amps2.D1.T2_3 + chan_ind3), &nph1, &nhp1, &nhp1, &p1, &p1, &N, &N);
      //T2_4(ab|ij){bj',ia'}  =   T2_1(bc|kj){bj',kc'}.Xhphp2_1(ka|ic){kc',ia'}
      dgemm_NN((Amps1.D1.T2_1 + chan_ind1), (Eff_Ints.Xhphp.X_2_1 + chan_ind5), (Amps2.D1.T2_4 + chan_ind4), &nph1, &nhp1, &nhp1, &p1, &p1, &N, &N);
    }
  }

  if(Parameters.basis != "finite_J"){ fac1 = m1; fac2 = p1; }
  else{ fac1 = p1; fac2 = m1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];
    length = np * nhhp;
    if(length != 0){
      chan_ind1 = Amps1.D1.T3_1_index[chan3];
      chan_ind2 = Amps1.D1.T3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xpp.X_3_index[chan3];
      //T3_1(ab|ij){a,ijb'}  =  +(+) Xpp3_od(a|c){a,c}.T3_1(cb|ij){c,ijb'}
      dgemm_NN((Eff_Ints.Xpp.X_3od + chan_ind3), (Amps1.D1.T3_1 + chan_ind1), (Amps2.D1.T3_1 + chan_ind1), &np, &nhhp, &np, &p1, &p1, &N, &N);
      //T3_2(ab|ij){b,ija'}  =  -(+) Xpp3_od(b|c){b,c}.T3_1(ca|ij){c,ija'}
      dgemm_NN((Eff_Ints.Xpp.X_3od + chan_ind3), (Amps1.D1.T3_1 + chan_ind1), (Amps2.D1.T3_2 + chan_ind2), &np, &nhhp, &np, &fac1, &p1, &N, &N);
    }
    length = nh * npph;
    if(length != 0){
      chan_ind1 = Amps1.D1.T3_3_index[chan3];
      chan_ind2 = Amps1.D1.T3_4_index[chan3];
      chan_ind3 = Eff_Ints.Xhh.X_3_index[chan3];
      //T3_3(ab|ij){abj',i}  =  -(-) T3_3(ab|kj){abj',k}.Xhh3_od(k|i){k,i}
      dgemm_NN((Amps1.D1.T3_3 + chan_ind1), (Eff_Ints.Xhh.X_3od + chan_ind3), (Amps2.D1.T3_3 + chan_ind1), &npph, &nh, &nh, &m1, &p1, &N, &N);
      //T3_4(ab|ij){abi',j}  =  +(-) T3_3(ab|ki){abi',k}.Xhh3_od(k|j){k,j}
      dgemm_NN((Amps1.D1.T3_3 + chan_ind1), (Eff_Ints.Xhh.X_3od + chan_ind3), (Amps2.D1.T3_4 + chan_ind2), &npph, &nh, &nh, &fac2, &p1, &N, &N);
    }
  }
}

void Singles_Step(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5, fac;
  int nh, np, npph, nhhp, one = 1;
  int nhp0 = Chan.nhp1[Chan.ind0];
  int nph0 = Chan.nph1[Chan.ind0];
  int chan_ind1, chan_ind2, chan_ind3;
  char N = 'N';

  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];
    chan_ind1 = Amps1.S1.t3_index[chan3];
    if(nh * np != 0){
      if(Parameters.HF == 0){
	chan_ind2 = Ints.Fmatrix.ph_3_index[chan3];
	// t3(a|i){a,i}  <-  fph(a|i){a,i}
	for(int p = 0; p < np; ++p){
	  for(int h = 0; h < nh; ++h){
	    Amps2.S1.t3[chan_ind1 + (p*nh + h)] += Ints.Fmatrix.ph_3[chan_ind2 + (p*nh + h)];
	  }
	}
      }

      chan_ind2 = Eff_Ints.Xpp.X_3_index[chan3];
      //t3(a|i){a,i}  <-  Xpp3_od(a|c){a,c}.t3(c|i){c,i}
      dgemm_NN((Eff_Ints.Xpp.X_3od + chan_ind2), (Amps1.S1.t3 + chan_ind1), (Amps2.S1.t3 + chan_ind1), &np, &nh, &np, &p1, &p1, &N, &N);

      chan_ind2 = Eff_Ints.Xhh.X_3_index[chan3];
      //t3(a|i){a,i}  <-  - t3(a|k){a,k}.X1hh3_od(k|i){k,i}
      dgemm_NN((Amps1.S1.t3 + chan_ind1), (Eff_Ints.Xhh.X1_3od + chan_ind2), (Amps2.S1.t3 + chan_ind1), &np, &nh, &nh, &m1, &p1, &N, &N);
    }
    if(npph * np * nh != 0){
      chan_ind2 = Ints.Vhppp.V_3_2_index[chan3];
      chan_ind3 = Amps1.D1.T3_4_index[chan3];
      //t3(a|i){a,i}  <-  (1/2) Vhppp3_2(ka|cd){a,cdk'}.T3_4(cd|ki){cdk',i}
      dgemm_NN((Ints.Vhppp.V_3_2 + chan_ind2), (Amps1.D1.T3_4 + chan_ind3), (Amps2.S1.t3 + chan_ind1), &np, &nh, &npph, &p12, &p1, &N, &N);
    }
    if(nhhp * np * nh != 0){
      chan_ind2 = Amps1.D1.T3_1_index[chan3];
      chan_ind3 = Ints.Vhhhp.V_3_3_index[chan3];
      //t3(a|i){a,i}  <-  -(1/2) T3_1(ac|kl){a,klc'}.Vhhhp3_3(kl|ic){klc',i}
      dgemm_NN((Amps1.D1.T3_1 + chan_ind2), (Ints.Vhhhp.V_3_3 + chan_ind3), (Amps2.S1.t3 + chan_ind1), &np, &nh, &nhhp, &m12, &p1, &N, &N);
    }
  }
  if(Parameters.basis != "finite_J"){ fac = m1; }
  else{ fac = p1; }
  if(nhp0 != 0){
    chan_ind1 = Ints.Vhphp.V_2_2_index[Chan.ind0];
    //t2(a|i){ai'}  =  -(+) Vhphp2_2(ka|ic){ai',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhphp.V_2_2 + chan_ind1), Amps1.S1.t2, Amps2.S1.t2, &nph0, &one, &nph0, &fac, &p1, &N, &N);
  }
  if(nph0 * nhp0 != 0){
    chan_ind1 = Amps1.D1.T2_3_index[Chan.ind0];
    //t2(a|i){ai'}  <-  + T2_3(ac|ik){ai',kc'}.Xhp2(k|c){kc'}
    dgemm_NN((Amps1.D1.T2_3 + chan_ind1), Eff_Ints.Xhp.X_2, Amps2.S1.t2, &nph0, &one, &nhp0, &p1, &p1, &N, &N);
  }
}


void Doubles_Step_2(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, m1 = -1.0, fac1, fac2;
  int nh, np, npph, nhhp;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';

  if(Parameters.basis != "finite_J"){ fac1 = p1; fac2 = m1; }
  else{ fac1 = m1; fac2 = p1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];
    if(nh * np * nhhp != 0){
      chan_ind1 = Amps1.S1.t3_index[chan3];
      chan_ind2 = Eff_Ints.Xhphh.X_3_1_index[chan3];
      chan_ind3 = Amps1.D1.T3_1_index[chan3];
      chan_ind4 = Amps1.D1.T3_2_index[chan3];
      //T3_1(ab|ij){a,ijb'}  =  - t3(a|k){a,k}.X1hphh3_1(kb|ij){k,ijb'}
      dgemm_NN((Amps1.S1.t3 + chan_ind1), (Eff_Ints.Xhphh.X1_3_1 + chan_ind2), (Amps2.D1.T3_1 + chan_ind3), &np, &nhhp, &nh, &m1, &p1, &N, &N);
      //T3_2(ab|ij){b,ija'}  =  +(-) t3(b|k){b,k}.X1hphh3_1(ka|ij){k,ija'}
      dgemm_NN((Amps1.S1.t3 + chan_ind1), (Eff_Ints.Xhphh.X1_3_1 + chan_ind2), (Amps2.D1.T3_2 + chan_ind4), &np, &nhhp, &nh, &fac1, &p1, &N, &N);
    }
    if(nh * np * npph != 0){
      chan_ind1 = Eff_Ints.Xpphp.X_3_4_index[chan3];
      chan_ind2 = Amps1.S1.t3_index[chan3];
      chan_ind3 = Amps1.D1.T3_4_index[chan3];
      chan_ind4 = Amps1.D1.T3_3_index[chan3];
      //T3_4(ab|ij){abi',j}  =  + X1pphp3_4(ab|ic){abi',c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xpphp.X1_3_4 + chan_ind1), (Amps1.S1.t3 + chan_ind2), (Amps2.D1.T3_4 + chan_ind3), &npph, &nh, &np, &p1, &p1, &N, &N);
      //T3_3(ab|ij){abj',i}  =  -(+) X1pphp3_4(ab|jc){abj',c}.t3(c|i){c,i}
      dgemm_NN((Eff_Ints.Xpphp.X1_3_4 + chan_ind1), (Amps1.S1.t3 + chan_ind2), (Amps2.D1.T3_3 + chan_ind4), &npph, &nh, &np, &fac2, &p1, &N, &N);
    }
  }
}

void Update_Heff_3(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5, fac;
  int nh, np, npph, nhhp, nhph, nhpp, nhh, npp, nhp, nhp1, nph1, npp1, nhh1, nhhh, nppp, length;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';

  //////////////   Xhhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    length = nhhp * nh;
    chan_ind1 = Eff_Ints.Xhhhp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vhhhp.V_3_3_index[chan3];
    // Xhhhp3_3 = Vhhhp3_3
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhhhp.X_3_3[chan_ind1 + ind] += Ints.Vhhhp.V_3_3[chan_ind2 + ind];
    }
    if(nh * np * nhhp != 0){
      chan_ind1 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhhhp.X_3_3_index[chan3];
      // Xhhhp3_3(ik|ja){ika',j}  <-  + Vhhpp3_3(ik|ca){ika',c}.t3(c|j){c,j}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhhhp.X_3_3 + chan_ind3), &nhhp, &nh, &np, &p1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhhhp.X_gather(Parameters, Chan);

  //////////////   Xpppp   ///////////////
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    length = npp * npp;
    chan_ind1 = Eff_Ints.Xpppp.X_1_index[chan1];
    // Xpppp1(ab|cd){ab,cd}  <-  X1pppp1(ab|cd){ab,cd}
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xpppp.X_1[chan_ind1 + ind] += Eff_Ints.Xpppp.X1_1[chan_ind1 + ind];
    }
    if(nhh * npp != 0){
      chan_ind1 = Amps.D1.T1_index[chan1];
      chan_ind2 = Ints.Vhhpp.V_1_index[chan1];
      chan_ind3 = Eff_Ints.Xpppp.X_1_index[chan1];
      // Xpppp1(ab|cd){ab,cd}  <-  +(1/2) T1(ab|kl){ab,kl}.Vhhpp1(kl|cd){kl,cd}
      dgemm_NN((Amps.D1.T1 + chan_ind1), (Ints.Vhhpp.V_1 + chan_ind2), (Eff_Ints.Xpppp.X_1 + chan_ind3), &npp, &npp, &nhh, &p12, &p1, &N, &N);
    }
  }

  //////////////   X3hphp,Xhphp   ///////////////
  if(Parameters.basis != "finite_J"){ fac = m12; }
  else{ fac = p12; }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    length = nhp1 * nhp1;
    chan_ind1 = Eff_Ints.Xhphp.X_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    // X3hphp2_1 = Vhphp2_1
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhphp.X3_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
    }
    if(nhp1 * nph1 != 0){
      chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
      chan_ind2 = Amps.D1.T2_1_index[chan2];
      chan_ind3 = Eff_Ints.Xhphp.X_2_1_index[chan2];
      // Xhphp2_1(ia|jb){ib',ja'}  <-  -(+)(1/2) Vhhpp2_1(ik|cb){ib',ck'}.T2_1(ca|jk){ck',ja'}
      dgemm_NN((Ints.Vhhpp.V_2_1 + chan_ind1), (Amps.D1.T2_1 + chan_ind2), (Eff_Ints.Xhphp.X_2_1 + chan_ind3), &nhp1, &nhp1, &nph1, &fac, &p1, &N, &N);
    }
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhpp = Chan.nhpp[chan3];
    nhph = Chan.nhph[chan3];
    if(nhpp * nh * np != 0){
      chan_ind1 = Eff_Ints.Xhppp.X_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhphp.X_3_3_index[chan3];
      // X3hphp3_3(ia|jb){iab',j}  =  +(1/2) Xhppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhppp.X_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphp.X3_3_3 + chan_ind3), &nhpp, &nh, &np, &p12, &p1, &N, &N);
    }
    if(nhph * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhhp.V_3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhphp.X_3_2_index[chan3];
      // X3hphp3_2(ia|jb){a,jbi'}  =  - t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhp.V_3_2 + chan_ind2), (Eff_Ints.Xhphp.X3_3_2 + chan_ind3), &np, &nhph, &nh, &m1, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhphp.X3_gather(Parameters, Chan);

  //////////////   Xpphp   ///////////////
  if(Parameters.basis != "finite_J"){ fac = p1; }
  else{ fac = m1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    nhpp = Chan.nhpp[chan3];
    npph = Chan.npph[chan3];

    length = nppp * nh;
    chan_ind1 = Eff_Ints.Xpphp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vpphp.V_3_3_index[chan3];
    // Xpphp3_3 = Vpphp3_3
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xpphp.X_3_3[chan_ind1 + ind] += Ints.Vpphp.V_3_3[chan_ind2 + ind];
    }
    if(nppp * nh * np != 0){
      chan_ind1 = Ints.Vpppp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_3_index[chan3];
      // Xpphp3_3(ab|ic){abc',i}  <-  Vpppp3_3(ab|dc){abc',d}.t3(d|i){d,i}
      dgemm_NN((Ints.Vpppp.V_3_3 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xpphp.X_3_3 + chan_ind3), &nppp, &nh, &np, &p1, &p1, &N, &N);
    }
    if(nhpp * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Eff_Ints.Xhphp.X_3_1_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_1_index[chan3];
      chan_ind4 = Eff_Ints.Xpphp.X_3_2_index[chan3];
      // Xpphp3_1(ab|ic){a,icb'}  =  - t3(a|k){a,l}.X1hphp3_1(kb|ic){k,icb'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhphp.X1_3_1 + chan_ind2), (Eff_Ints.Xpphp.X_3_1 + chan_ind3), &np, &nhpp, &nh, &m1, &p1, &N, &N);
      // Xpphp3_2(ab|ic){a,icb'}  =  +(-) t3(a|k){a,l}.X1hphp3_1(kb|ic){k,icb'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Eff_Ints.Xhphp.X1_3_1 + chan_ind2), (Eff_Ints.Xpphp.X_3_2 + chan_ind4), &np, &nhpp, &nh, &fac, &p1, &N, &N);
    }
    if(npph * np * nh != 0){
      chan_ind1 = Amps.D1.T3_4_index[chan3];
      chan_ind2 = Eff_Ints.Xhp.X_3_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_4_index[chan3];
      // Xpphp3_4(ab|ic){abi,c}  <-  - T3_4(ab|ik){abi,k}.Xhp3(k|c){k,c}
      dgemm_NN((Amps.D1.T3_4 + chan_ind1), (Eff_Ints.Xhp.X_3 + chan_ind2), (Eff_Ints.Xpphp.X_3_4 + chan_ind3), &npph, &np, &nh, &m1, &p1, &N, &N);
    }
  }
  if(Parameters.basis != "finite_J"){ fac = m1; }
  else{ fac = p1; }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    npp1 = Chan.npp1[chan2];
    if(nph1 * npp1 * nhp1 != 0){
      chan_ind1 = Amps.D1.T2_3_index[chan2];
      chan_ind2 = Eff_Ints.Xhppp.X_2_3_index[chan2];
      chan_ind3 = Eff_Ints.Xpphp.X_2_3_index[chan2];
      chan_ind4 = Eff_Ints.Xpphp.X_2_2_index[chan2];
      // Xpphp2_3(ab|ic){ai',cb'}  <-  T2_3(ad|ik){ai',kd'}.Xhppp2_3(kb|dc){kd',cb'}
      dgemm_NN((Amps.D1.T2_3 + chan_ind1), (Eff_Ints.Xhppp.X_2_3 + chan_ind2), (Eff_Ints.Xpphp.X_2_3 + chan_ind3), &nph1, &npp1, &nhp1, &p1, &p1, &N, &N);
      // Xpphp2_2(ab|ic){bi',ca'}  <-  -(+) T2_3(bd|ik){bi',kd'}.Xhppp2_3(ka|dc){kd',ca'}
      dgemm_NN((Amps.D1.T2_3 + chan_ind1), (Eff_Ints.Xhppp.X_2_3 + chan_ind2), (Eff_Ints.Xpphp.X_2_2 + chan_ind4), &nph1, &npp1, &nhp1, &fac, &p1, &N, &N);
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    nhp = Chan.nhp[chan1];
    if(nhh * npp * nhp != 0){
      chan_ind1 = Amps.D1.T1_index[chan1];
      chan_ind2 = Eff_Ints.Xhhhp.X_1_index[chan1];
      chan_ind3 = Eff_Ints.Xpphp.X_1_index[chan1];
      // Xpphp1(ab|ic){ab,ic}  <-  +(1/2) T1(ab|kl){ab,kl}.Xhhhp1(kl|ic){kl,ic}
      dgemm_NN((Amps.D1.T1 + chan_ind1), (Eff_Ints.Xhhhp.X_1 + chan_ind2), (Eff_Ints.Xpphp.X_1 + chan_ind3), &npp, &nhp, &nhh, &p12, &p1, &N, &N);
    }
  }
  Eff_Ints.Xpphp.X_gather(Parameters, Chan);

  //////////////   Xhphh   ///////////////
  if(Parameters.basis != "finite_J"){ fac = m1; }
  else{ fac = p1; }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhph = Chan.nhph[chan3];
    nhhp = Chan.nhhp[chan3];
    nhhh = Chan.nhhh[chan3];

    length = np * nhhh;
    chan_ind1 = Eff_Ints.Xhphh.X_3_2_index[chan3];
    chan_ind2 = Ints.Vhphh.V_3_2_index[chan3];
    // Xhphh3_2  <-  Vhphh3_2
    for(int ind = 0; ind < length; ++ind){
      Eff_Ints.Xhphh.X_3_2[chan_ind1 + ind] += Ints.Vhphh.V_3_2[chan_ind2 + ind];
    }
    if(nhhh * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhhh.V_3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhphh.X_3_2_index[chan3];
      // Xhphh3_2(ia|jk){a,jki'}  <-  - t3(a|l){a,l}.Vhhhh3_2(il|jk){l,jki'}
      dgemm_NN((Amps.S1.t3 + chan_ind1), (Ints.Vhhhh.V_3_2 + chan_ind2), (Eff_Ints.Xhphh.X_3_2 + chan_ind3), &np, &nhhh, &nh, &m1, &p1, &N, &N);
    }
    if(nhph * np * nh != 0){
      chan_ind1 = Eff_Ints.Xhphp.X_3_4_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = Eff_Ints.Xhphh.X_3_4_index[chan3];
      chan_ind4 = Eff_Ints.Xhphh.X_3_3_index[chan3];
      // Xhphh3_4(ia|jk){iaj,k}  <-  X3hphp3_4(ia|jc){iaj,c}.t3(c|k){c,k}
      dgemm_NN((Eff_Ints.Xhphp.X3_3_4 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphh.X_3_4 + chan_ind3), &nhph, &nh, &np, &p1, &p1, &N, &N);
      // Xhphh3_3(ia|jk){iak,j}  <-  -(+) X3hphp3_4(ia|kc){iak,c}.t3(c|j){c,j}
      dgemm_NN((Eff_Ints.Xhphp.X3_3_4 + chan_ind1), (Amps.S1.t3 + chan_ind2), (Eff_Ints.Xhphh.X_3_3 + chan_ind4), &nhph, &nh, &np, &fac, &p1, &N, &N);
    }
    if(nh * nhhp * np != 0){
      chan_ind1 = Eff_Ints.Xhp.X_3_index[chan3];
      chan_ind2 = Amps.D1.T3_1_index[chan3];
      chan_ind3 = Eff_Ints.Xhphh.X_3_1_index[chan3];
      // Xhphh3_1(ia|jk){i,jka}  <-  Xhp3(i|c){i,c}.T3_1(ca|jk){c,jka}
      dgemm_NN((Eff_Ints.Xhp.X_3 + chan_ind1), (Amps.D1.T3_1 + chan_ind2), (Eff_Ints.Xhphh.X_3_1 + chan_ind3), &nh, &nhhp, &np, &p1, &p1, &N, &N);
    }
  }
  if(Parameters.basis != "finite_J"){ fac = m1; }
  else{ fac = p1; }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    nhh1 = Chan.nhh1[chan2];
    if(nhh1 * nhp1 * nph1 != 0){
      chan_ind1 = Eff_Ints.Xhhhp.X_2_3_index[chan2];
      chan_ind2 = Amps.D1.T2_3_index[chan2];
      chan_ind3 = Eff_Ints.Xhphh.X_2_3_index[chan2];
      chan_ind4 = Eff_Ints.Xhphh.X_2_1_index[chan2];
      // Xhphh2_3(ia|jk){ij',ka'}  <-  Xhhhp2_3(il|jc){ij',cl'}.T2_3(ca|lk){cl',ka'}
      dgemm_NN((Eff_Ints.Xhhhp.X_2_3 + chan_ind1), (Amps.D1.T2_3 + chan_ind2), (Eff_Ints.Xhphh.X_2_3 + chan_ind3), &nhh1, &nhp1, &nph1, &p1, &p1, &N, &N);
      // Xhphh2_1(ia|jk){ik',ja'}  <-  -(+) Xhhhp2_3(il|kc){ik',cl'}.T2_3(ca|lj){cl',ja'}
      dgemm_NN((Eff_Ints.Xhhhp.X_2_3 + chan_ind1), (Amps.D1.T2_3 + chan_ind2), (Eff_Ints.Xhphh.X_2_1 + chan_ind4), &nhh1, &nhp1, &nph1, &fac, &p1, &N, &N);
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    nhp = Chan.nhp[chan1];
    if(nhh * npp * nhp != 0){
      chan_ind1 = Eff_Ints.Xhppp.X_1_index[chan1];
      chan_ind2 = Amps.D1.T1_index[chan1];
      chan_ind3 = Eff_Ints.Xhphh.X_1_index[chan1];
      // Xhphh1(ia|jk){ia,jk}  <-  +(1/2) Xhppp1(ia|cd){ia,cd}.T1(cd|jk){cd,jk}
      dgemm_NN((Eff_Ints.Xhppp.X_1 + chan_ind1), (Amps.D1.T1 + chan_ind2), (Eff_Ints.Xhphh.X_1 + chan_ind3), &nhp, &nhh, &npp, &p12, &p1, &N, &N);
    }
  }
  Eff_Ints.Xhphh.X_gather(Parameters, Chan);
}

void EOM::PA1_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  double *Ham;
  State state;
  int length;
  int np0, npph0, npph, np_tot, npph_tot;
  three_body *pph_vec;
  State *J_state;
  int *J_vec;
  int *a_vec;
  int *b_vec;
  int *i_vec;
  int chan_j;
  int chan, N;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int p1, p2, h1;

  // Count N_states and setup EOM structures //
  count_states(Parameters, Space, Chan, 1);
  np_tot = 0;
  npph_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan]; 
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (Parameters.basis == "finite_J" && p1 == p2) ){ ++npph0; }
    }
    nob[n] = np0;
    nthb[n] = npph0;
    nstate[n] = np0 + npph0;
    ob_index[n] = np_tot;
    thb_index[n] = npph_tot;
    state_index[n] = np_tot + npph_tot;
    np_tot += np0;
    npph_tot += npph0;
  }
  ob_vec = new one_body[np_tot];
  thb_vec = new three_body[npph_tot];
  thb_qnums = new State[npph_tot];
  state_vec_R = new double[np_tot + npph_tot];
  state_vec_L = new double[np_tot + npph_tot];
  /////////////////////////////////////////////

  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan];
    chan_j = Chan.qnums3[chan].j;

    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (Parameters.basis == "finite_J" && p1 == p2) ){ ++npph0; }
    }
    pph_vec = new three_body[npph0];
    J_state = new State[npph0];
    J_vec = new int[npph0];
    a_vec = new int[npph0];
    b_vec = new int[npph0];
    i_vec = new int[npph0];

    // set EOM array
    npph0 = 0;
    for(int p = 0; p < np0; ++p){ ob_vec[ob_index[n] + p] = Chan.p_state(chan, p); }
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      h1 = Chan.pph_state(chan, pph).v3;
      if( (p1 < p2) || (Parameters.basis == "finite_J" && p1 == p2) ){
	pph_vec[npph0] = Chan.pph_state(chan, pph);
	J_state[npph0] = Chan.pph_j[Chan.pph_index[chan] + pph];
	J_vec[npph0] = Chan.pph_j[Chan.pph_index[chan] + pph].j;
	a_vec[npph0] = Space.qnums[p1].j;
	b_vec[npph0] = Space.qnums[p2].j;
	i_vec[npph0] = Space.qnums[h1].j;
	// set EOM arrays
	thb_vec[thb_index[n] + npph0] = Chan.pph_state(chan, pph);
	thb_qnums[thb_index[n] + npph0] = Chan.pph_j[Chan.pph_index[chan] + pph];
	++npph0;
      }
    }

    N = np0 + npph0;
    Ham = new double[N*N];
    eigenvalues = new double[1];
    eigenvectors_R = new double[N];
    eigenvectors_L = new double[N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }

    #pragma omp parallel
    {
      double ME, X;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind, jmin;
      int bra, ket, key1, key2;
      State tb1, tb2;
      double fac1, fac2;

      length = N * N;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%N);
	bra = int((braket - ket)/N);
	fac2 = 1.0;
	ME = 0.0;

	if( bra >= np0 ){
	  bra -= np0;
	  b1 = pph_vec[bra].v1;
	  b2 = pph_vec[bra].v2;
	  b3 = pph_vec[bra].v3;
	  if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	  if( ket >= np0 ){  //  < pph | pph >
	    ket -= np0;
	    k1 = pph_vec[ket].v1;
	    k2 = pph_vec[ket].v2;
	    k3 = pph_vec[ket].v3;
	    if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      if(b1 == k1 && equal(Space.qnums[b2], Space.qnums[k2])){
		minus(tb1, Space.qnums[b2], Space.qnums[k2]);
		ind = Chan.pp1_map[Chan.ind0][Space.hash2(b2, k2, 0)];
		ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(b_vec[bra] + 1.0);
	      }
	      if(b1 == k2 && equal(Space.qnums[b2], Space.qnums[k1])){
		minus(tb1, Space.qnums[b2], Space.qnums[k1]);
		ind = Chan.pp1_map[Chan.ind0][Space.hash2(b2, k1, 0)];
		ME -= Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(b_vec[bra] + 1.0);
	      }
	      if(b2 == k2 && equal(Space.qnums[b1], Space.qnums[k1])){
		minus(tb1, Space.qnums[b1], Space.qnums[k1]);
		ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k1, 0)];
		ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[bra] + 1.0);
	      }
	      if(b2 == k1 && equal(Space.qnums[b1], Space.qnums[k2])){
		minus(tb1, Space.qnums[b1], Space.qnums[k2]);
		ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k2, 0)];
		ME -= Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(a_vec[bra] + 1.0);
	      }
	    }

	    if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) && equal(Space.qnums[b3], Space.qnums[k3])){
	      minus(tb1, Space.qnums[k3], Space.qnums[b3]);
	      ind = Chan.hh1_map[Chan.ind0][Space.hash2(k3, b3, 0)];
	      ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0);
	      if(b1 == b2){ ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0); }
	    }
	    
	    fac1 = -1.0;
	    if(Parameters.basis == "finite_J"){ fac1 *= -1.0; }
	    if(b1 == k1){
	      minus(tb1, Space.qnums[k3], Space.qnums[k2]);
	      minus(tb2, Space.qnums[b3], Space.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k3].j - Space.qnums[k2].j);
		if(abs(Space.qnums[b3].j - Space.qnums[b2].j) > jmin){ jmin = abs(Space.qnums[b3].j - Space.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b1 == k2){
	      minus(tb1, Space.qnums[k3], Space.qnums[k1]);
	      minus(tb2, Space.qnums[b3], Space.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k3].j - Space.qnums[k1].j);
		if(abs(Space.qnums[b3].j - Space.qnums[b2].j) > jmin){ jmin = abs(Space.qnums[b3].j - Space.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k2){
	      minus(tb1, Space.qnums[k3], Space.qnums[k1]);
	      minus(tb2, Space.qnums[b3], Space.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k3].j - Space.qnums[k1].j);
		if(abs(Space.qnums[b3].j - Space.qnums[b1].j) > jmin){ jmin = abs(Space.qnums[b3].j - Space.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(a_vec[bra] + b_vec[bra] - J_vec[bra] + a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k1){
	      minus(tb1, Space.qnums[k3], Space.qnums[k2]);
	      minus(tb2, Space.qnums[b3], Space.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k3].j - Space.qnums[k2].j);
		if(abs(Space.qnums[b3].j - Space.qnums[b1].j) > jmin){ jmin = abs(Space.qnums[b3].j - Space.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    
	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      chan0 = Space.ind_2b_dir(Parameters.basis, J_state[bra]);
	      key1 = Chan.pp_map[chan0][Space.hash2(b1, b2, J_state[bra].j)];
	      key2 = Chan.pp_map[chan0][Space.hash2(k1, k2, J_state[bra].j)];
	      chan_ind = Eff_Ints.Xpppp.X_1_index[chan0];
	      ind = key1 * Chan.npp[chan0] + key2;
	      ME += Eff_Ints.Xpppp.X_1[chan_ind + ind];
	    }
	    ket += np0;
	  }
	  else{  //  < pph | p >
	    k1 = Chan.p_state(chan, ket).v1;
	    
	    fac1 = -1.0;
	    if(Parameters.basis == "finite_J"){ fac1 *= -1.0; }
	    key1 = Chan.pph_map[chan][Space.hash3(b1, b2, b3, J_vec[bra])];
	    key2 = Chan.p_map[chan][k1];
	    chan_ind = Eff_Ints.Xpphp.X_3_4_index[chan];
	    ind = key1 * np0 + key2;
	    ME += fac1 * Eff_Ints.Xpphp.X_3_4[chan_ind + ind];
	  }
	  bra += np0;
	}
	else if( ket >= np0 ){  //  < p | pph >
	  ket -= np0;
	  k1 = pph_vec[ket].v1;
	  k2 = pph_vec[ket].v2;
	  k3 = pph_vec[ket].v3;
	  b1 = Chan.p_state(chan, bra).v1;
	  if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	  if(k1 == b1 && equal(Space.qnums[k3], Space.qnums[k2])){
	    minus(tb1, Space.qnums[k3], Space.qnums[k2]);
	    ind = Chan.hp1_map[Chan.ind0][Space.hash2(k3, k2, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (b_vec[ket] + 1.0))) * phase2(a_vec[ket] + b_vec[ket] - J_vec[ket]);
	    ME += Eff_Ints.Xhp.X_2[ind] * X;
	  }
	  if(k2 == b1 && equal(Space.qnums[k3], Space.qnums[k1])){
	    minus(tb1, Space.qnums[k3], Space.qnums[k1]);
	    ind = Chan.hp1_map[Chan.ind0][Space.hash2(k3, k1, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (a_vec[ket] + 1.0)));
	    ME -= Eff_Ints.Xhp.X_2[ind] * X;
	  }
	  
	  fac1 = -1.0;
	  if(Parameters.basis == "finite_J"){ fac1 *= -1.0; }
	  plus(tb1, Space.qnums[k1], Space.qnums[k2]);
	  key1 = Chan.p_map[chan][b1];
	  key2 = Chan.pph_map[chan][Space.hash3(k1, k2, k3, J_vec[ket])];
	  chan_ind = Eff_Ints.Xhppp.X_3_2_index[chan];
	  ind = key1 * npph + key2;
	  ME += fac1 * Eff_Ints.Xhppp.X_3_2[chan_ind + ind];
	  ket += np0;
	}
	else{  //  < p | p >
	  b1 = Chan.p_state(chan, bra).v1;
	  k1 = Chan.p_state(chan, ket).v1;

	  minus(tb1, Space.qnums[b1], Space.qnums[k1]);
	  ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k1, 0)];
	  ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(chan_j + 1.0);
	}
	Ham[N*ket + bra] += fac2 * ME;
      }
    }

    if(N <= 1000){ Asym_Diagonalize1(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }

    norm1p = 0.0;
    for(int i = 0; i < np0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);

    del_E[n] = eigenvalues[0];
    norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
    delete[] pph_vec;
    delete[] J_state;
    delete[] J_vec;
    delete[] a_vec;
    delete[] b_vec;
    delete[] i_vec;
    delete[] Ham;
  }
}

void EOM::PR1_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints)
{
  double *Ham;
  State state;
  int length;
  int nh0, nhhp0, nhhp, nh_tot, nhhp_tot;
  three_body *hhp_vec;
  State *J_state;
  int *J_vec;
  int *i_vec;
  int *j_vec;
  int *a_vec;
  int chan_j;
  int chan, N;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int h1, h2, p1;

  // Count N_states and setup EOM structures //
  count_states(Parameters, Space, Chan, 0);
  nh_tot = 0;
  nhhp_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan]; 
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (Parameters.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    nob[n] = nh0;
    nthb[n] = nhhp0;
    nstate[n] = nh0 + nhhp0;
    ob_index[n] = nh_tot;
    thb_index[n] = nhhp_tot;
    state_index[n] = nh_tot + nhhp_tot;
    nh_tot += nh0;
    nhhp_tot += nhhp0;
  }
  ob_vec = new one_body[nh_tot];
  thb_vec = new three_body[nhhp_tot];
  thb_qnums = new State[nhhp_tot];
  state_vec_R = new double[nh_tot + nhhp_tot];
  state_vec_L = new double[nh_tot + nhhp_tot];
  /////////////////////////////////////////////

  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan];
    chan_j = Chan.qnums3[chan].j;

    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (Parameters.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    hhp_vec = new three_body[nhhp0];
    J_state = new State[nhhp0];
    J_vec = new int[nhhp0];
    i_vec = new int[nhhp0];
    j_vec = new int[nhhp0];
    a_vec = new int[nhhp0];

    // set EOM array    
    nhhp0 = 0;
    for(int h = 0; h < nh0; ++h){ ob_vec[ob_index[n] + h] = Chan.h_state(chan, h); }
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      p1 = Chan.hhp_state(chan, hhp).v3;
      if( (h1 < h2) || (Parameters.basis == "finite_J" && h1 == h2) ){
	hhp_vec[nhhp0] = Chan.hhp_state(chan, hhp);
	J_state[nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp];
	J_vec[nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp].j;
	i_vec[nhhp0] = Space.qnums[h1].j;
	j_vec[nhhp0] = Space.qnums[h2].j;
	a_vec[nhhp0] = Space.qnums[p1].j;
	// set EOM arrays
	thb_vec[thb_index[n] + nhhp0] = Chan.hhp_state(chan, hhp);
	thb_qnums[thb_index[n] + nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp];
	++nhhp0;
      }
    }

    N = nh0 + nhhp0;
    Ham = new double[N*N];
    eigenvalues = new double[1];
    eigenvectors_R = new double[N];
    eigenvectors_L = new double[N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }
    
    #pragma omp parallel
    {
      double ME, X;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind, jmin;
      int bra, ket, key1, key2;
      State tb1, tb2;
      double fac1, fac2;

      length = N * N;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%N);
	bra = int((braket - ket)/N);
	fac2 = 1.0;
	ME = 0.0;
	  
	if( bra >= nh0 ){
	  bra -= nh0;
	  b1 = hhp_vec[bra].v1;
	  b2 = hhp_vec[bra].v2;
	  b3 = hhp_vec[bra].v3;
	  if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	  if( ket >= nh0 ){  //  < hhp | hhp >
	    ket -= nh0;
	    k1 = hhp_vec[ket].v1;
	    k2 = hhp_vec[ket].v2;
	    k3 = hhp_vec[ket].v3;
	    if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      if(b1 == k1 && equal(Space.qnums[b2], Space.qnums[k2])){
		minus(tb1, Space.qnums[k2], Space.qnums[b2]);
		ind = Chan.hh1_map[Chan.ind0][Space.hash2(k2, b2, 0)];
		ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(j_vec[bra] + 1.0);
	      }
	      if(b1 == k2 && equal(Space.qnums[b2], Space.qnums[k1])){
		minus(tb1, Space.qnums[k1], Space.qnums[k2]);
		ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b2, 0)];
		ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(j_vec[bra] + 1.0);
	      }
	      if(b2 == k2 && equal(Space.qnums[b1], Space.qnums[k1])){
		minus(tb1, Space.qnums[k1], Space.qnums[b1]);
		ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b1, 0)];
		ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[bra] + 1.0);
	      }
	      if(b2 == k1 && equal(Space.qnums[b1], Space.qnums[k2])){
		minus(tb1, Space.qnums[k2], Space.qnums[b1]);
		ind = Chan.hh1_map[Chan.ind0][Space.hash2(k2, b1, 0)];
		ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(i_vec[bra] + 1.0);
	      }
	    }

	    if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) && equal(Space.qnums[b3], Space.qnums[k3])){
	      minus(tb1, Space.qnums[b3], Space.qnums[k3]);
	      ind = Chan.pp1_map[Chan.ind0][Space.hash2(b3, k3, 0)];
	      ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0);
	      if(b1 == b2){ ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0); }
	    }

	    fac1 = -1.0;
	    if(Parameters.basis == "finite_J"){ fac1 *= -1.0; }
	    if(b1 == k1){
	      minus(tb1, Space.qnums[k2], Space.qnums[k3]);
	      minus(tb2, Space.qnums[b2], Space.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k2].j - Space.qnums[k3].j);
		if(abs(Space.qnums[b2].j - Space.qnums[b3].j) > jmin){ jmin = abs(Space.qnums[b2].j - Space.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(j_vec[bra] + j_vec[ket] + a_vec[bra] + a_vec[ket]);
		  X *= CGC6(i_vec[bra],j_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(i_vec[ket],j_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k2, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b2, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b1 == k2){
	      minus(tb1, Space.qnums[k1], Space.qnums[k3]);
	      minus(tb2, Space.qnums[b2], Space.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k1].j - Space.qnums[k3].j);
		if(abs(Space.qnums[b2].j - Space.qnums[b3].j) > jmin){ jmin = abs(Space.qnums[b2].j - Space.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);;
		  X *= phase2(j_vec[bra] + i_vec[ket] + a_vec[bra] + a_vec[ket] + i_vec[ket] + j_vec[ket] - J_vec[ket]);
		  X *= CGC6(i_vec[bra],j_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(j_vec[ket],i_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k1, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b2, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k2){
	      minus(tb1, Space.qnums[k1], Space.qnums[k3]);
	      minus(tb2, Space.qnums[b1], Space.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k1].j - Space.qnums[k3].j);
		if(abs(Space.qnums[b1].j - Space.qnums[b3].j) > jmin){ jmin = abs(Space.qnums[b1].j - Space.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(i_vec[bra] + i_vec[ket] + a_vec[bra] + a_vec[ket]);
		  X *= phase2(i_vec[bra] + j_vec[bra] - J_vec[bra] + i_vec[ket] + j_vec[ket] - J_vec[ket]);
		  X *= CGC6(j_vec[bra],i_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(j_vec[ket],i_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k1, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b1, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k1){
	      minus(tb1, Space.qnums[k2], Space.qnums[k3]);
	      minus(tb2, Space.qnums[b1], Space.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(Space.qnums[k2].j - Space.qnums[k3].j);
		if(abs(Space.qnums[b1].j - Space.qnums[b3].j) > jmin){ jmin = abs(Space.qnums[b1].j - Space.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(i_vec[bra] + j_vec[ket] + a_vec[bra] + a_vec[ket] + i_vec[bra] + j_vec[bra] - J_vec[bra]);
		  X *= CGC6(j_vec[bra],i_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(i_vec[ket],j_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Space.ind_2b_cross(Parameters.basis, tb1);
		  key1 = Chan.hp1_map[chan0][Space.hash2(k2, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Space.hash2(b1, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= fac1 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      chan0 = Space.ind_2b_dir(Parameters.basis, J_state[bra]);
	      key1 = Chan.hh_map[chan0][Space.hash2(k1, k2, J_state[bra].j)];
	      key2 = Chan.hh_map[chan0][Space.hash2(b1, b2, J_state[bra].j)];
	      chan_ind = Eff_Ints.Xhhhh.X_1_index[chan0];
	      ind = key1 * Chan.nhh[chan0] + key2;
	      ME += Eff_Ints.Xhhhh.X_1[chan_ind + ind];
	    }
	    ket += nh0;
	  }
	  else{  //  < hhp | h >
	    k1 = Chan.h_state(chan, ket).v1;

	    key1 = Chan.h_map[chan][k1];
	    key2 = Chan.hhp_map[chan][Space.hash3(b1, b2, b3, J_vec[bra])];
	    chan_ind = Eff_Ints.Xhphh.X_3_1_index[chan];
	    ind = key1 * nhhp + key2;
	    ME += -1.0 * Eff_Ints.Xhphh.X_3_1[chan_ind + ind];
	  }
	  bra += nh0;
	}
	else if( ket >= nh0 ){  //  < h | hhp >
	  ket -= nh0;
	  k1 = hhp_vec[ket].v1;
	  k2 = hhp_vec[ket].v2;
	  k3 = hhp_vec[ket].v3;
	  b1 = Chan.h_state(chan, bra).v1;
	  if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	    
	  if(k1 == b1 && equal(Space.qnums[k2], Space.qnums[k3])){
	    minus(tb1, Space.qnums[k2], Space.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Space.hash2(k2, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (j_vec[ket] + 1.0))) * phase2(i_vec[ket] + j_vec[ket] - J_vec[ket]);
	    ME += Eff_Ints.Xhp.X_2[ind] * X;
	  }
	  if(k2 == b1 && equal(Space.qnums[k1], Space.qnums[k3])){
	    minus(tb1, Space.qnums[k1], Space.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Space.hash2(k1, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (i_vec[ket] + 1.0)));
	    ME -= Eff_Ints.Xhp.X_2[ind] * X;
	  }

	  key1 = Chan.hhp_map[chan][Space.hash3(k1, k2, k3, J_vec[ket])];
	  key2 = Chan.h_map[chan][b1];
	  chan_ind = Eff_Ints.Xhhhp.X_3_3_index[chan];
	  ind = key1 * nh0 + key2;
	  ME += -1.0 * Eff_Ints.Xhhhp.X_3_3[chan_ind + ind];
	  ket += nh0;
	}
	else{  //  < h | h >
	  b1 = Chan.h_state(chan, bra).v1;
	  k1 = Chan.h_state(chan, ket).v1;

	  minus(tb1, Space.qnums[k1], Space.qnums[b1]);
	  ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b1, 0)];
	  ME += -1.0 * Eff_Ints.Xhh.X_2[ind] / std::sqrt(chan_j + 1.0);
	}
	Ham[N*ket + bra] += fac2 * ME;
      }
    }
    
    if(N <= 100){ Asym_Diagonalize1(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }

    norm1p = 0.0;
    for(int i = 0; i < nh0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);

    del_E[n] = eigenvalues[0];
    norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
    delete[] hhp_vec;
    delete[] J_state;
    delete[] J_vec;
    delete[] i_vec;
    delete[] j_vec;
    delete[] a_vec;
    delete[] Ham;
  }
}

void EOM::count_states(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, int type)
{
  int count1 = 0, count2 = 0, count_limit1 = 0, count_limit2 = 0;
  int num, ind, index;
  double en, tempen;
  double limit;
  int p_limit = 3;
  int n_limit = 3;
  if( type == 1 ){ limit = 1.0e10; }
  else if( type == 0 ){ limit = -1.0e10; }

  if( Parameters.Pshells != 0 ){
    for(int chan = 0; chan < Chan.size3; ++chan){
      if( type == 1 ){ num = Chan.np[chan]; }
      else if( type == 0 ){ num = Chan.nh[chan]; }
      if( num == 0 || Chan.qnums3[chan].t == 1 ){ continue; } // skip if no particles/holes in this channel
      if(Parameters.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
      else if(Parameters.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
      ++count_limit1;
    }
    if( count_limit1 > p_limit ){ count_limit1 = p_limit; }
  }
  if( Parameters.Nshells != 0 ){
    for(int chan = 0; chan < Chan.size3; ++chan){
      if( type == 1 ){ num = Chan.np[chan]; }
      else if( type == 0 ){ num = Chan.nh[chan]; }
      if( num == 0 || Chan.qnums3[chan].t == -1 ){ continue; } // skip if no particles/holes in this channel
      if(Parameters.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
      else if(Parameters.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
      ++count_limit2;
    }
    if( count_limit2 > n_limit ){ count_limit2 = n_limit; }
  }
  N_states = count_limit1 + count_limit2;   // EOM variable
  chan_vec = new int[N_states];             // EOM variable
  qnums = new State[N_states];              // EOM variable
  del_E = new double[N_states];             // EOM variable
  norm_1p = new double[N_states];           // EOM variable
  nob = new int[N_states];                  // EOM variable
  ob_index = new int[N_states];             // EOM variable
  nthb = new int[N_states];                 // EOM variable
  thb_index = new int[N_states];            // EOM variable
  nstate = new int[N_states];               // EOM variable
  state_index = new int[N_states];          // EOM variable

  if( Parameters.Pshells != 0 ){
    while( count1 < count_limit1 ){
      en = limit;
      for(int chan = 0; chan < Chan.size3; ++chan){
	if( type == 1 ){ num = Chan.np[chan]; index = 0; }
	else if( type == 0 ){ num = Chan.nh[chan]; index = num-1; }
	if( num == 0 || Chan.qnums3[chan].t == 1 ){ continue; } // skip if no particles/holes in this channel
	if(Parameters.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
        else if(Parameters.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
	for(int n = 0; n < count1; ++n){ if( chan_vec[n] == chan ){ goto stop1; } }
	if( type == 1 ){ tempen = Space.qnums[Chan.p_state(chan, index).v1].energy; }
	else if( type == 0 ){ tempen = Space.qnums[Chan.h_state(chan, index).v1].energy; }
	if( (type == 1 && tempen < en) || (type == 0 && tempen > en) ){
	  en = tempen;
	  ind = chan;
	}
      stop1:;
      }
      chan_vec[count1] = ind;
      qnums[count1] = Chan.qnums3[ind];
      ++count1;
    }
  }
  if( Parameters.Nshells != 0 ){
    while( count2 < count_limit2 ){
      en = limit;
      for(int chan = 0; chan < Chan.size3; ++chan){
	if( type == 1 ){ num = Chan.np[chan]; index = 0; }
	else if( type == 0 ){ num = Chan.nh[chan]; index = num-1; }
	if( num == 0 || Chan.qnums3[chan].t == -1 ){ continue; } // skip if no particles/holes in this channel
	if(Parameters.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
        else if(Parameters.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
	for(int n = count1; n < count1+count2; ++n){ if( chan_vec[n] == chan ){ goto stop2; } }
	if( type == 1 ){ tempen = Space.qnums[Chan.p_state(chan, index).v1].energy; }
	else if( type == 0 ){ tempen = Space.qnums[Chan.h_state(chan, index).v1].energy; }
	if( (type == 1 && tempen < en) || (type == 0 && tempen > en) ){
	  en = tempen;
	  ind = chan;
	}
      stop2:;
      }
      chan_vec[count1 + count2] = ind;
      qnums[count1 + count2] = Chan.qnums3[ind];
      ++count2;
    }
  }
}

void EOM::Print_EOM_1P(Input_Parameters &Parameters, double Energy)
{
  std::cout << std::fixed;
  if(Parameters.basis == "finite_HO"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << Parameters.Shells << std::setw(5) << Parameters.Pshells << std::setw(5) << qnums[i].ml;
      std::cout << std::setw(5) << qnums[i].m << std::setprecision(2) << std::setw(8) << Parameters.density;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
  else if(Parameters.basis == "finite_J"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << Parameters.Shells << std::setw(5) << Parameters.Pshells << std::setw(5) << Parameters.Nshells;
      std::cout << std::setw(5) << qnums[i].j << std::setw(5) << qnums[i].par << std::setw(5) << qnums[i].t;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
  else if(Parameters.basis == "finite_M"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << Parameters.Shells << std::setw(5) << Parameters.Pshells << std::setw(5) << Parameters.Nshells;
      std::cout << std::setw(5) << qnums[i].m << std::setw(5) << qnums[i].par << std::setw(5) << qnums[i].t;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
}

void EOM::delete_struct()
{
  if(N_states != 0){
    delete[] qnums;
    delete[] chan_vec;
    delete[] del_E;
    delete[] norm_1p;
    delete[] nob;
    delete[] ob_vec;
    delete[] ob_index;
    delete[] nthb;
    delete[] thb_vec;
    delete[] thb_index;
    delete[] thb_qnums;
    delete[] nstate;
    delete[] state_vec_R;
    delete[] state_vec_L;
    delete[] state_index;
  }
}
