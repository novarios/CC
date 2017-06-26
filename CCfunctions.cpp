#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

//   Function to return Hash index for 2 indices
int Hash2(int &p, int &q, int &size){ return size*p + q; }
//   Function to return Hash index for 3 indices
int Hash3(int &p, int &q, int &r, int &size){ return size*size*p + q*size + r; }

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
      energy += Space.qnums[i].energy;
    }
  }
  else if(Parameters.HF == 0){
    for(int chan = 0; chan < Chan.size3; ++chan){
      nh = Chan.nh[chan];
      chan_ind = Ints.Fmatrix.hh_3_index[chan];
      for(int h = 0; h < nh; ++h){
	ind = h*nh + h;
	energy += Ints.Fmatrix.hh_3[chan_ind + ind];
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
  int chan1, chan2, chan3, ind, length, count;
  int a, b, i, j, npp, nhh;
  two_body pp, hh;
  double J;
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
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	pp = Chan.pp_state(chan1, pp_ind);
	hh = Chan.hh_state(chan1, hh_ind);
	a = pp.v1;
	b = pp.v2;
	i = hh.v1;
	j = hh.v2;
	// T2_1, T2_2, T2_3, T2_4, T3_1, T3_2, T3_3, T3_4
	Map_4_count(Parameters,Space, map_index,map_num, 8,1,0,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,2,1,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,3,2,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,4,3,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,5,4,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,6,5,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,7,6,length,count, a,b,i,j);
	Map_4_count(Parameters,Space, map_index,map_num, 8,8,7,length,count, a,b,i,j);
	++count;
      }
    }
  }
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];
  
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	pp = Chan.pp_state(chan1, pp_ind);
	hh = Chan.hh_state(chan1, hh_ind);
	a = pp.v1;
	b = pp.v2;
	i = hh.v1;
	j = hh.v2;
	// T2_1, T2_2, T2_3, T2_4, T3_1, T3_2, T3_3, T3_4
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,1,0,count, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,2,1,count, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,3,2,count, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,4,3,count, Chan.ph1_map,Chan.hp1_map, Chan.nhp1, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,5,4,count, Chan.p_map,Chan.hhp_map, Chan.nhhp, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,6,5,count, Chan.p_map,Chan.hhp_map, Chan.nhhp, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,7,6,count, Chan.pph_map,Chan.h_map, Chan.nh, a,b,i,j,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 8,8,7,count, Chan.pph_map,Chan.h_map, Chan.nh, a,b,i,j,J);

	chan3 = Space.ind_1b(Parameters.basis, Space.qnums[i]);
	Evec_chan[4*count] = chan3;
	Evec_ind[4*count] = Chan.h_map[chan3][i];
	chan3 = Space.ind_1b(Parameters.basis, Space.qnums[j]);
	Evec_chan[4*count+1] = chan3;
	Evec_ind[4*count+1] = Chan.h_map[chan3][j];
	chan3 = Space.ind_1b(Parameters.basis, Space.qnums[a]);
	Evec_chan[4*count+2] = chan3;
	Evec_ind[4*count+2] = Chan.p_map[chan3][a];
	chan3 = Space.ind_1b(Parameters.basis, Space.qnums[b]);
	Evec_chan[4*count+3] = chan3;
	Evec_ind[4*count+3] = Chan.p_map[chan3][b];
	++count;
      }
    }
  }
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

  int a, i;
  two_body ph1;
  for(ind = 0; ind < t2_length; ++ind){
    ph1 = Chan.ph1_state(Chan.ind0, ind);
    a = ph1.v1;
    i = ph1.v2;
    Map_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, ind, Chan.p_map,Chan.h_map,Chan.nh, a,i);

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
  while(error > 1e-12 && ind < 10000){
    Update_CC(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);

    //Print_Amps(Parameters, Chan, tempAmps);
    Amps2.copy_Amplitudes(Parameters, Chan, Amps);
    Amps2.zero(Parameters, Chan, true);
    Gather_Amps(Parameters, Space, Chan, Eff_Ints, Amps2, tempAmps, mix);
    //Print_Amps(Parameters, Chan, Amps2);
    
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

  while( ind < 10 ){
    // J
    Update_CC(Parameters, Space, Chan, Ints, Eff_Ints, Amps, tempAmps);
    // M
    Update_CC(ParametersM, SpaceM, ChanM, IntsM, Eff_IntsM, AmpsM, tempAmpsM);


    std::cout << " ****  tempAmps  **** " << std::endl;
    CC_Test_Full_J(Parameters, Space, Chan, Ints, tempAmps, ParametersM, SpaceM, ChanM, IntsM, tempAmpsM);


    // J
    Amps2.copy_Amplitudes(Parameters, Chan, Amps);
    Amps2.zero(Parameters, Chan, true);
    Gather_Amps(Parameters, Space, Chan, Eff_Ints, Amps2, tempAmps, mix);
    // M
    Amps2M.copy_Amplitudes(ParametersM, ChanM, AmpsM);
    Amps2M.zero(ParametersM, ChanM, true);
    Gather_Amps(ParametersM, SpaceM, ChanM, Eff_IntsM, Amps2M, tempAmpsM, mix);


    // J
    Amps.copy_Amplitudes(Parameters, Chan, Amps2);
    // M
    AmpsM.copy_Amplitudes(ParametersM, ChanM, Amps2M);


    std::cout << " ****  Amps  **** " << std::endl;
    CC_Test_Full_J(Parameters, Space, Chan, Ints, Amps, ParametersM, SpaceM, ChanM, IntsM, AmpsM);


    // J
    CCoutE = Amps.get_energy(Parameters, Space, Chan, Ints);
    // M
    CCoutEM = AmpsM.get_energy(ParametersM, SpaceM, ChanM, IntsM);
    std::cout << "Iteration Number = " << ind << ", EnergyJ = " << CCoutE << ", EnergyM = " << CCoutEM << std::endl;
    ++ind;
  }

  // J
  Amps2.delete_struct(Parameters, Chan);
  tempAmps.delete_struct(Parameters, Chan);
  // M
  StatesM.delete_struct(HF_ChanM);
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
  int ind;
  int length = Amps.D1.T1_length;
  for(ind = 0; ind < length; ++ind){
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
    for(ind = 0; ind < length; ++ind){
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
      //dgemm_NN((Ints.Vhhpp.V_2_1 + chan_ind1), (Amps.D1.T2_1 + chan_ind2), (Eff_Ints.Xhphp.X_2_1 + chan_ind3), &nhp1, &nhp1, &nph1, &m12, &p1, &N, &N);
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
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    length = nhh * npp;
    chan_ind1 = Amps1.D1.T1_index[chan1];
    chan_ind2 = Ints.Vpphh.V_1_index[chan1];
    //T1(ab|ij){ab,ij}  =  Vpphh1(ab|ij){ab,ij}
    for(int ind = 0; ind < length; ++ind){
      Amps2.D1.T1[chan_ind1 + ind] += Ints.Vpphh.V_1[chan_ind2 + ind];
    }

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


void PA_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, State *states, double *nums)
{
  double *Ham;
  State state;
  int length;
  int np0, npph0, npph;
  int *pph_vec;
  int N;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }
    
    np0 = Chan.np[chan];
    npph0 = Chan.npph[chan];
    npph = 0;
    for(int pph = 0; pph < npph0; ++pph){
      if( Chan.pph_state(chan, pph).v1 < Chan.pph_state(chan, pph).v2 ){
	++npph;
      }
    }
    pph_vec = new int[3*npph];
    npph = 0;
    for(int pph = 0; pph < npph0; ++pph){
      if( Chan.pph_state(chan, pph).v1 < Chan.pph_state(chan, pph).v2 ){
	pph_vec[3*npph] = Chan.pph_state(chan, pph).v1;
	pph_vec[3*npph + 1] = Chan.pph_state(chan, pph).v2;
	pph_vec[3*npph + 2] = Chan.pph_state(chan, pph).v3;
	++npph;
      }
    }

    N = np0 + npph;
    if(N == 0){ continue; }
    Ham = new double[N*N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind;
      int bra, ket, key1, key2;
      State tb;
      length = np0 * np0;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%np0);
	bra = int((braket - ket)/np0);
	b1 = Chan.p_state(chan, bra).v1;
	k1 = Chan.p_state(chan, ket).v1;
	minus(tb, Space.qnums[b1], Space.qnums[k1]);
	ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k1, tb.j)];
	ME = Eff_Ints.Xpp.X_2[ind];
	//std::cout << "1: " << bra << " " << ket << " += Xpp(" << b1 << "," << k1 << ") = " << Eff_Ints.Xpp.X_2[ind] << std::endl;
	Ham[N*ket + bra] += ME;
      }
      length = npph * np0;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%np0);
	bra = int((braket - ket)/np0);
	b1 = pph_vec[3*bra];
	b2 = pph_vec[3*bra + 1];
	b3 = pph_vec[3*bra + 2];
	k1 = Chan.p_state(chan, ket).v1;
	plus(tb, Space.qnums[b1], Space.qnums[b2]);
	key1 = Chan.pph_map[chan][Space.hash3(b1, b2, b3, tb.j)];
	key2 = Chan.p_map[chan][k1];
	chan_ind = Eff_Ints.Xpphp.X_3_4_index[chan];
	ind = key1 * np0 + key2;
	ME = -1.0 * Eff_Ints.Xpphp.X_3_4[chan_ind + ind];
	//std::cout << "2: " << (bra+np0) << " " << ket << " += Xpphp(" << b1 << "," << b2 << "," << b3 << "," << k1 << ") = " << -1.0 * Eff_Ints.Xpphp.X_3_4[chan_ind + ind] << std::endl;
	Ham[N*ket + (bra + np0)] += ME;
      }
      length = np0 * npph;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%np0);
	ket = int((braket - bra)/np0);
	k1 = pph_vec[3*ket];
	k2 = pph_vec[3*ket + 1];
	k3 = pph_vec[3*ket + 2];
	b1 = Chan.p_state(chan, bra).v1;
	ME = 0.0;

	if(k1 == b1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  ind = Chan.hp1_map[Chan.ind0][Space.hash2(k3, k2, tb.j)];
	  ME += Eff_Ints.Xhp.X_2[ind];
	  //std::cout << "3_1: " << bra << " " << (ket+np0) << " += Xhp(" << k3 << "," << k2 << ") = " << Eff_Ints.Xhp.X_2[ind] << std::endl;
	}
	else if(k2 == b1){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind = Chan.hp1_map[Chan.ind0][Space.hash2(k3, k1, tb.j)];
	  ME -= Eff_Ints.Xhp.X_2[ind];
	  //std::cout << "3_2: " << bra << " " << (ket+np0) << " += Xhp(" << k3 << "," << k1 << ") = " << -1.0 * Eff_Ints.Xhp.X_2[ind] << std::endl;
	}

 	plus(tb, Space.qnums[k1], Space.qnums[k2]);
	key1 = Chan.p_map[chan][b1];
	key2 = Chan.pph_map[chan][Space.hash3(k1, k2, k3, tb.j)];
	chan_ind = Eff_Ints.Xhppp.X_3_2_index[chan];
	ind = key1 * npph0 + key2;
	ME -= Eff_Ints.Xhppp.X_3_2[chan_ind + ind];
	//std::cout << "3_3: " << bra << " " << (ket+np0) << " += Xhppp(" << b1 << "," << k1 << "," << k2 << "," << k3 << ") = " << -1.0 * Eff_Ints.Xhppp.X_3_2[chan_ind + ind] << std::endl;
	Ham[N*(ket + np0) + bra] += ME;
      }
      length = npph * npph;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%npph);
	bra = int((braket - ket)/npph);
	b1 = pph_vec[3*bra];
	b2 = pph_vec[3*bra + 1];
	b3 = pph_vec[3*bra + 2];
	k1 = pph_vec[3*ket];
	k2 = pph_vec[3*ket + 1];
	k3 = pph_vec[3*ket + 2];
	ME = 0.0;

	if(b3 == k3){
	  if(b1 == k1){
	    minus(tb, Space.qnums[b2], Space.qnums[k2]);
	    ind = Chan.pp1_map[Chan.ind0][Space.hash2(b2, k2, tb.j)];
	    ME += Eff_Ints.Xpp.X_2[ind];
	    //std::cout << "4_1: " << (bra+np0) << " " << (ket+np0) << " += Xpp(" << b2 << "," << k2 << ") = " << Eff_Ints.Xpp.X_2[ind] << std::endl;
	  }
	  else if(b1 == k2){
	    minus(tb, Space.qnums[b2], Space.qnums[k1]);
	    ind = Chan.pp1_map[Chan.ind0][Space.hash2(b2, k1, tb.j)];
	    ME -= Eff_Ints.Xpp.X_2[ind];
	    //std::cout << "4_2: " << (bra+np0) << " " << (ket+np0) << " += Xpp(" << b2 << "," << k1 << ") = " << -1.0 * Eff_Ints.Xpp.X_2[ind] << std::endl;
	  }

	  if(b2 == k2){
	    minus(tb, Space.qnums[b1], Space.qnums[k1]);
	    ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k1, tb.j)];
	    ME += Eff_Ints.Xpp.X_2[ind];
	    //std::cout << "4_3: " << (bra+np0) << " " << (ket+np0) << " += Xpp(" << b1 << "," << k1 << ") = " << Eff_Ints.Xpp.X_2[ind] << std::endl;
	  }
	  else if(b2 == k1){
	    minus(tb, Space.qnums[b1], Space.qnums[k2]);
	    ind = Chan.pp1_map[Chan.ind0][Space.hash2(b1, k2, tb.j)];
	    ME -= Eff_Ints.Xpp.X_2[ind];
	    //std::cout << "4_4: " << (bra+np0) << " " << (ket+np0) << " += Xpp(" << b1 << "," << k2 << ") = " << -1.0 * Eff_Ints.Xpp.X_2[ind] << std::endl;
	  }
	}

	if(b1 == k1 && b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[b3]);
	  ind = Chan.hh1_map[Chan.ind0][Space.hash2(k3, b3, tb.j)];
	  ME -= Eff_Ints.Xhh.X_2[ind];
	  //std::cout << "4_5: " << (bra+np0) << " " << (ket+np0) << " += Xhh(" << k3 << "," << b3 << ") = " << -1.0 * Eff_Ints.Xhh.X_2[ind] << std::endl;
	}

	if(b1 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k3, k2, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b3, b2, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_6: " << (bra+np0) << " " << (ket+np0) << " += Xhphp(" << k3 << "," << k2 << "," << b3 << "," << b2 << ") = " << -1.0 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}
	else if(b1 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k3, k1, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b3, b2, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_7: " << (bra+np0) << " " << (ket+np0) << " += Xhphp(" << k3 << "," << k1 << "," << b3 << "," << b2 << ") = " << Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}

	if(b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k3, k1, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b3, b1, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_8: " << (bra+np0) << " " << (ket+np0) << " += Xhphp(" << k3 << "," << k1 << "," << b3 << "," << b1 << ") = " << -1.0 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}
	else if(b2 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k3, k2, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b3, b1, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_9: " << (bra+np0) << " " << (ket+np0) << " += Xhphp(" << k3 << "," << k2 << "," << b3 << "," << b1 << ") = " << Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}

	if(b3 == k3){
	  plus(tb, Space.qnums[b1], Space.qnums[b2]);
	  chan0 = Space.ind_2b_dir(Parameters.basis, tb);
	  plus(tb, Space.qnums[k1], Space.qnums[k2]);
	  if( chan0 == Space.ind_2b_dir(Parameters.basis, tb) ){
	    key1 = Chan.pp_map[chan0][Space.hash2(b1, b2, tb.j)];
	    key2 = Chan.pp_map[chan0][Space.hash2(k1, k2, tb.j)];
	    chan_ind = Eff_Ints.Xpppp.X_1_index[chan0];
	    ind = key1 * Chan.npp[chan0] + key2;
	    ME += Eff_Ints.Xpppp.X_1[chan_ind + ind];
	    //std::cout << "4_10: " << (bra+np0) << " " << (ket+np0) << " += Xpppp(" << b1 << "," << b2 << "," << k1 << "," << k2 << ") = " << Eff_Ints.Xpppp.X_1[chan_ind + ind] << std::endl;
	  }
	}
	Ham[N*(ket + np0) + (bra + np0)] += ME;
      }
    }

    double norm1p;
    double eigenvalue;
    if(N <= 50){ Asym_Diagonalize1(Ham, N, eigenvalue, norm1p, np0); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalue, norm1p, np0); }
    states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
    nums[2*Chan.qnums3[chan].ml] = eigenvalue;
    nums[2*Chan.qnums3[chan].ml + 1] = norm1p;

    delete[] pph_vec;
    delete[] Ham;
  }
}

void PR_EOM(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Eff_Interactions &Eff_Ints, State *states, double *nums)
{
  double *Ham;
  State state;
  int length;
  int nh0, nhhp0, nhhp;
  int *hhp_vec;
  int N;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }

    nh0 = Chan.nh[chan];
    nhhp0 = Chan.nhhp[chan];
    nhhp = 0;
    for(int hhp = 0; hhp < nhhp0; ++hhp){
      if( Chan.hhp_state(chan, hhp).v1 < Chan.hhp_state(chan, hhp).v2 ){
	++nhhp;
      }
    }
    hhp_vec = new int[3*nhhp];
    nhhp = 0;
    for(int hhp = 0; hhp < Chan.nhhp[chan]; ++hhp){
      if( Chan.hhp_state(chan, hhp).v1 < Chan.hhp_state(chan, hhp).v2 ){
	hhp_vec[3*nhhp] = Chan.hhp_state(chan, hhp).v1;
	hhp_vec[3*nhhp + 1] = Chan.hhp_state(chan, hhp).v2;
	hhp_vec[3*nhhp + 2] = Chan.hhp_state(chan, hhp).v3;
	++nhhp;
      }
    }

    N = nh0 + nhhp;
    if(N == 0){ continue; }
    Ham = new double[N*N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind;
      int bra, ket, key1, key2;
      State tb;
      length = nh0 * nh0;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%nh0);
	bra = int((braket - ket)/nh0);
	Ham[N*ket + bra] = 0.0;
	b1 = Chan.h_state(chan, bra).v1;
	k1 = Chan.h_state(chan, ket).v1;
	minus(tb, Space.qnums[k1], Space.qnums[b1]);
	ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b1, tb.j)];
	ME = -1.0 * Eff_Ints.Xhh.X_2[ind];
	//std::cout << "1: " << bra << " " << ket << " += Xhh(" << k1 << "," << b1 << ") = " << -1.0 * Eff_Ints.Xhh.X_2[ind] << std::endl;
	Ham[N*ket + bra] = ME;
      }
      length = nhhp * nh0;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%nh0);
	bra = int((braket - ket)/nh0);
	b1 = hhp_vec[3*bra];
	b2 = hhp_vec[3*bra + 1];
	b3 = hhp_vec[3*bra + 2];
	k1 = Chan.h_state(chan, ket).v1;
	plus(tb, Space.qnums[b1], Space.qnums[b2]);
	key1 = Chan.h_map[chan][k1];
	key2 = Chan.hhp_map[chan][Space.hash3(b1, b2, b3, tb.j)];
	chan_ind = Eff_Ints.Xhphh.X_3_1_index[chan];
	ind = key1 * nhhp0 + key2;
	ME = -1.0 * Eff_Ints.Xhphh.X_3_1[chan_ind + ind];
	//std::cout << "2: " << (bra+nh0) << " " << ket << " += Xhphh(" << k1 << "," << b1 << "," << b2 << "," << b3 << ") = " << -1.0 * Eff_Ints.Xhphh.X_3_1[chan_ind + ind] << std::endl;
	Ham[N*ket + (bra + nh0)] = ME;
      }
      length = nh0 * nhhp;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%nh0);
	ket = int((braket - bra)/nh0);
	k1 = hhp_vec[3*ket];
	k2 = hhp_vec[3*ket + 1];
	k3 = hhp_vec[3*ket + 2];
	b1 = Chan.h_state(chan, bra).v1;
	ME = 0.0;

	if(k1 == b1){
	  minus(tb, Space.qnums[k2], Space.qnums[k3]);
	  ind = Chan.hp1_map[Chan.ind0][Space.hash2(k2, k3, tb.j)];
	  ME += Eff_Ints.Xhp.X_2[ind];
	  //std::cout << "3_1: " << bra << " " << (ket+nh0) << " += Xhp(" << k2 << "," << k3 << ") = " << Eff_Ints.Xhp.X_2[ind] << std::endl;
	}
	else if(k2 == b1){
	  minus(tb, Space.qnums[k1], Space.qnums[k3]);
	  ind = Chan.hp1_map[Chan.ind0][Space.hash2(k1, k3, tb.j)];
	  ME -= Eff_Ints.Xhp.X_2[ind];
	  //std::cout << "3_2: " << bra << " " << (ket+nh0) << " += Xhp(" << k1 << "," << k3 << ") = " << -1.0 * Eff_Ints.Xhp.X_2[ind] << std::endl;
	}

	plus(tb, Space.qnums[k1], Space.qnums[k2]);
	key1 = Chan.hhp_map[chan][Space.hash3(k1, k2, k3, tb.j)];
	key2 = Chan.h_map[chan][b1];
	chan_ind = Eff_Ints.Xhhhp.X_3_3_index[chan];
	ind = key1 * nh0 + key2;
	ME -= Eff_Ints.Xhhhp.X_3_3[chan_ind + ind];
	//std::cout << "3_3: " << bra << " " << (ket+nh0) << " += Xhhhp(" << k1 << "," << k2 << "," << k3 << "," << b1 << ") = " << -1.0 * Eff_Ints.Xhhhp.X_3_3[chan_ind + ind] << std::endl;
	Ham[N*(ket + nh0) + bra] = ME;
      }
      length = nhhp * nhhp;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%nhhp);
	bra = int((braket - ket)/nhhp);
	b1 = hhp_vec[3*bra];
	b2 = hhp_vec[3*bra + 1];
	b3 = hhp_vec[3*bra + 2];
	k1 = hhp_vec[3*ket];
	k2 = hhp_vec[3*ket + 1];
	k3 = hhp_vec[3*ket + 2];
	ME = 0.0;

	if(b3 == k3){
	  if(b1 == k1){
	    minus(tb, Space.qnums[k2], Space.qnums[b2]);
	    ind = Chan.hh1_map[Chan.ind0][Space.hash2(k2, b2, tb.j)];
	    ME -= Eff_Ints.Xhh.X_2[ind];
	    //std::cout << "4_1: " << (bra+nh0) << " " << (ket+nh0) << " += Xhh(" << k2 << "," << b2 << ") = " << -1.0 * Eff_Ints.Xhh.X_2[ind] << std::endl;
	  }
	  else if(b1 == k2){
	    minus(tb, Space.qnums[k1], Space.qnums[k2]);
	    ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b2, tb.j)];
	    ME += Eff_Ints.Xhh.X_2[ind];
	    //std::cout << "4_2: " << (bra+nh0) << " " << (ket+nh0) << " += Xhh(" << k1 << "," << b2 << ") = " << Eff_Ints.Xhh.X_2[ind] << std::endl;
	  }

	  if(b2 == k2){
	    minus(tb, Space.qnums[k1], Space.qnums[b1]);
	    ind = Chan.hh1_map[Chan.ind0][Space.hash2(k1, b1, tb.j)];
	    ME -= Eff_Ints.Xhh.X_2[ind];
	    //std::cout << "4_3: " << (bra+nh0) << " " << (ket+nh0) << " += Xhh(" << k1 << "," << b1 << ") = " << -1.0 * Eff_Ints.Xhh.X_2[ind] << std::endl;
	  }
	  else if(b2 == k1){
	    minus(tb, Space.qnums[k2], Space.qnums[b1]);
	    ind = Chan.hh1_map[Chan.ind0][Space.hash2(k2, b1, tb.j)];
	    ME += Eff_Ints.Xhh.X_2[ind];
	    //std::cout << "4_4: " << (bra+nh0) << " " << (ket+nh0) << " += Xhh(" << k2 << "," << b1 << ") = " << Eff_Ints.Xhh.X_2[ind] << std::endl;
	  }
	}

	if(b1 == k1 && b2 == k2){
	  minus(tb, Space.qnums[b3], Space.qnums[k3]);
	  ind = Chan.pp1_map[Chan.ind0][Space.hash2(b3, k3, tb.j)];
	  ME += Eff_Ints.Xpp.X_2[ind];
	  //std::cout << "4_5: " << (bra+nh0) << " " << (ket+nh0) << " += Xpp(" << b3 << "," << k3 << ") = " << Eff_Ints.Xpp.X_2[ind] << std::endl;
	}

	if(b1 == k1){
	  minus(tb, Space.qnums[k2], Space.qnums[k3]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b3]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k2, k3, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b2, b3, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_6: " << (bra+nh0) << " " << (ket+nh0) << " += Xhphp(" << k2 << "," << k3 << "," << b2 << "," << b3 << ") = " << -1.0 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}
	else if(b1 == k2){
	  minus(tb, Space.qnums[k1], Space.qnums[k3]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b3]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k1, k3, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b2, b3, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_7: " << (bra+nh0) << " " << (ket+nh0) << " += Xhphp(" << k1 << "," << k3 << "," << b2 << "," << b3 << ") = " << Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}

	if(b2 == k2){
	  minus(tb, Space.qnums[k1], Space.qnums[k3]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b1], Space.qnums[b3]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k1, k3, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b1, b3, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_8: " << (bra+nh0) << " " << (ket+nh0) << " += Xhphp(" << k1 << "," << k3 << "," << b1 << "," << b3 << ") = " << -1.0 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}
	if(b2 == k1){
	  minus(tb, Space.qnums[k2], Space.qnums[k3]);
	  chan0 = Space.ind_2b_cross(Parameters.basis, tb);
	  minus(tb, Space.qnums[b1], Space.qnums[b3]);
	  if( chan0 == Space.ind_2b_cross(Parameters.basis, tb) ){
	    key1 = Chan.hp1_map[chan0][Space.hash2(k2, k3, tb.j)];
	    key2 = Chan.hp1_map[chan0][Space.hash2(b1, b3, tb.j)];
	    chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
	    ind = key1 * Chan.nhp1[chan0] + key2;
	    ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind];
	    //std::cout << "4_9: " << (bra+nh0) << " " << (ket+nh0) << " += Xhphp(" << k2 << "," << k3 << "," << b1 << "," << b3 << ") = " << Eff_Ints.Xhphp.X_2_1[chan_ind + ind] << std::endl;
	  }
	}

	if(b3 == k3){
	  plus(tb, Space.qnums[b1], Space.qnums[b2]);
	  chan0 = Space.ind_2b_dir(Parameters.basis, tb);
	  plus(tb, Space.qnums[k1], Space.qnums[k2]);
	  if( chan0 == Space.ind_2b_dir(Parameters.basis, tb)){
	    key1 = Chan.hh_map[chan0][Space.hash2(k1, k2, tb.j)];
	    key2 = Chan.hh_map[chan0][Space.hash2(b1, b2, tb.j)];
	    chan_ind = Eff_Ints.Xhhhh.X_1_index[chan0];
	    ind = key1 * Chan.nhh[chan0] + key2;
	    ME += Eff_Ints.Xhhhh.X_1[chan_ind + ind];
	    //std::cout << "4_10: " << (bra+nh0) << " " << (ket+nh0) << " += Xhhhh(" << k1 << "," << k2 << "," << b1 << "," << b2 << ") = " << Eff_Ints.Xhhhh.X_1[chan_ind + ind] << std::endl;
	  }
	}
	Ham[N*(ket + nh0) + (bra + nh0)] += ME;
      }
    }

    double norm1p;
    double eigenvalue;
    if(N <= 50){ Asym_Diagonalize1(Ham, N, eigenvalue, norm1p, nh0); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalue, norm1p, nh0); }
    states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
    nums[2*Chan.qnums3[chan].ml] = eigenvalue;
    nums[2*Chan.qnums3[chan].ml + 1] = norm1p;
    
    delete[] hhp_vec;
    delete[] Ham;
  }
}

void Map_4_count(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int type, int offset, int &length, int &count, int &p, int &q, int &r, int &s)
{
  State tb;
  int jmin;
  int num;
  if(type == 1){
    cross_state(Parameters, Space, p, s, r, q, jmin, tb); // 2_1
    num = int(0.5*(tb.j - jmin) + 1);
    map_num[size * count + offset] = num;
    map_index[size * count + offset] = length;
    length += num;
  }
  else if(type == 2){
    cross_state(Parameters, Space, q, r, s, p, jmin, tb); // 2_2
    num = int(0.5*(tb.j - jmin) + 1);
    map_num[size * count + offset] = num;
    map_index[size * count + offset] = length;
    length += num;
  }
  else if(type == 3){
    cross_state(Parameters, Space, p, r, s, q, jmin, tb); // 2_3
    num = int(0.5*(tb.j - jmin) + 1);
    map_num[size * count + offset] = num;
    map_index[size * count + offset] = length;
    length += num;
  }
  else if(type == 4){
    cross_state(Parameters, Space, q, s, r, p, jmin, tb); // 2_4
    num = int(0.5*(tb.j - jmin) + 1);
    map_num[size * count + offset] = num;
    map_index[size * count + offset] = length;
    length += num;
  }
  else if(type <= 8){
    map_num[size * count + offset] = 1;	                  // 3_1, 3_2, 3_3, 3_4
    map_index[size * count + offset] = length;
    length += 1;
  }
}

void Map_4(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int type, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, int &p, int &q, int &r, int &s, double &J)
{
  State tb;
  int jmin, key1, key2;
  int chan, ind;
  double jp, jq, jr, js, J1, X;
  int J2 = int(2 * J);
  jp = 0.5 * Space.qnums[p].j;
  jq = 0.5 * Space.qnums[q].j;
  jr = 0.5 * Space.qnums[r].j;
  js = 0.5 * Space.qnums[s].j;
  
  int index = map_index[size * count + offset];
  int num = map_num[size * count + offset];
  for(int n = 0; n < num; ++n){
    if(type == 1){
      cross_state(Parameters, Space, p, s, r, q, jmin, tb);
      tb.j -= 2*n;
      chan = Space.ind_2b_cross(Parameters.basis, tb);
      key1 = map1[chan][Space.hash2(p, s, tb.j)];
      key2 = map2[chan][Space.hash2(r, q, tb.j)];
      J1 = 0.5 * tb.j;
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = -1.0 * CGC6(jp, jq, J, jr, js, J1); }
      map_fac1[index + n] = (2.0 * J + 1) * X;
      map_fac2[index + n] = (2.0 * J1 + 1) * X;
    }
    else if(type == 2){
      cross_state(Parameters, Space, q, r, s, p, jmin, tb);
      tb.j -= 2*n;
      chan = Space.ind_2b_cross(Parameters.basis, tb);
      key1 = map1[chan][Space.hash2(q, r, tb.j)];
      key2 = map2[chan][Space.hash2(s, p, tb.j)];
      J1 = 0.5 * tb.j;
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = -1.0 * std::pow(-1.0, jp + jq + jr + js) * CGC6(jq, jp, J, js, jr, J1); }
      map_fac1[index + n] = (2.0 * J + 1) * X;
      map_fac2[index + n] = (2.0 * J1 + 1) * X;
    }
    else if(type == 3){
      cross_state(Parameters, Space, p, r, s, q, jmin, tb);
      tb.j -= 2*n;
      chan = Space.ind_2b_cross(Parameters.basis, tb);
      key1 = map1[chan][Space.hash2(p, r, tb.j)];
      key2 = map2[chan][Space.hash2(s, q, tb.j)];
      J1 = 0.5 * tb.j;
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = std::pow(-1.0, jr + js - J) * CGC6(jp, jq, J, js, jr, J1); }
      map_fac1[index + n] = (2.0 * J + 1) * X;
      map_fac2[index + n] = (2.0 * J1 + 1) * X;
    }
    else if(type == 4){
      cross_state(Parameters, Space, q, s, r, p, jmin, tb);
      tb.j -= 2*n;
      chan = Space.ind_2b_cross(Parameters.basis, tb);
      key1 = map1[chan][Space.hash2(q, s, tb.j)];
      key2 = map2[chan][Space.hash2(r, p, tb.j)];
      J1 = 0.5 * tb.j;
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = std::pow(-1.0, jp + jq - J) * CGC6(jq, jp, J, jr, js, J1); }
      map_fac1[index + n] = (2.0 * J + 1) * X;
      map_fac2[index + n] = (2.0 * J1 + 1) * X;
    }
    else if(type == 5){
      chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
      key1 = map1[chan][p];
      key2 = map2[chan][Space.hash3(r, s, q, J2)];
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = std::pow(-1.0, jp + jq - J) * std::sqrt((2.0 * J + 1)/(2.0 * jp + 1)); }
      map_fac1[index + n] = X;
      map_fac2[index + n] = 1.0 / X;
    }
    else if(type == 6){
      chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
      key1 = map1[chan][q];
      key2 = map2[chan][Space.hash3(r, s, p, J2)];
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = -1.0 * std::sqrt((2.0 * J + 1)/(2.0 * jq + 1)); }
      map_fac1[index + n] = X;
      map_fac2[index + n] = 1.0 / X;
    }
    else if(type == 7){
      chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
      key2 = map2[chan][r];
      key1 = map1[chan][Space.hash3(p, q, s, J2)];
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = std::pow(-1.0, jr + js - J) * std::sqrt((2.0 * J + 1)/(2.0 * jr + 1)); }
      map_fac1[index + n] = X;
      map_fac2[index + n] = 1.0 / X;
    }
    else if(type == 8){
      chan = Space.ind_1b(Parameters.basis, Space.qnums[s]);
      key2 = map2[chan][s];
      key1 = map1[chan][Space.hash3(p, q, r, J2)];
      if( jp == 0.0 ){ X = 1.0; }
      else{ X = -1.0 * std::sqrt((2.0 * J + 1)/(2.0 * js + 1)); }
      map_fac1[index + n] = X;
      map_fac2[index + n] = 1.0 / X;
    }
    ind = key1 * num2[chan] + key2;
    map_chan[index + n] = chan;
    map_ind[index + n] = ind;
  }
}

void Map_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, int &p, int &q)
{
  int chan, ind, key1, key2;
  double X, jp = 0.5 * Space.qnums[p].j;
  chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
  key1 = map1[chan][p];
  key2 = map2[chan][q];
  ind = key1 * num2[chan] + key2;
  X = std::sqrt(2.0 * jp + 1);
  map_chan[count] = chan;
  map_ind[count] = ind;
  map_fac1[count] = 1.0/X;
  map_fac2[count] = X;
}

void direct_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &q, int &r, int &s, int &jmin1, State &tb1)
{
  plus(tb1, Space.qnums[p], Space.qnums[q]);
  jmin1 = abs(Space.qnums[p].j - Space.qnums[q].j);
  if(abs(Space.qnums[r].j - Space.qnums[s].j) > jmin1){ jmin1 = abs(Space.qnums[r].j - Space.qnums[s].j); }
  if(Space.qnums[r].j + Space.qnums[s].j < tb1.j){ tb1.j = Space.qnums[r].j + Space.qnums[s].j; }
}

void cross_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &s, int &r, int &q, int &jmin2, State &tb2)
{
  minus(tb2, Space.qnums[p], Space.qnums[s]);
  jmin2 = abs(Space.qnums[p].j - Space.qnums[s].j);
  if(abs(Space.qnums[r].j - Space.qnums[q].j) > jmin2){ jmin2 = abs(Space.qnums[r].j - Space.qnums[q].j); }
  if(Space.qnums[r].j + Space.qnums[q].j < tb2.j){ tb2.j = Space.qnums[r].j + Space.qnums[q].j; }
}
