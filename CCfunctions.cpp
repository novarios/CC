#include "CCfunctions.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "INTfunctions.hpp"
#include "EffINTfunctions.hpp"
#include "DIISfunctions.hpp"
#include "MATHfunctions.hpp"

void Amplitudes::get_dE(Channels &Chan, Interactions &Ints)
{
  double energy = 0.0;
  double energy1 = 0.0; // for singles part
  two_body ph;
  int nhh, npp, nh, np, nph0;
  int ind1, ind2, chan_ind1, chan_ind2, hp_ind, a, i;

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    chan_ind1 = D1.T1_index[chan];
    chan_ind2 = Ints.Vhhpp.V_1_index[chan];
    for(int pp = 0; pp < npp; ++pp){
      for(int hh = 0; hh < nhh; ++hh){
	ind1 = pp * nhh + hh;
	ind2 = hh * npp + pp;
	if(PAR.basis == "finite_J"){ energy += (Chan.qnums1[chan].j + 1) * D1.T1[chan_ind1 + ind1] * Ints.Vhhpp.V_1[chan_ind2 + ind2]; }
	else{ energy += D1.T1[chan_ind1 + ind1] * Ints.Vhhpp.V_1[chan_ind2 + ind2]; }
      }
    }
  }
  energy *= 0.25;
  if(PAR.approx == "singles" || PAR.approx == "triples"){
    nph0 = Chan.nph1[Chan.ind0];
    chan_ind1 = Ints.Vhhpp.V_2_3_index[Chan.ind0];
    for(int ph1 = 0; ph1 < nph0; ++ph1){
      ph = Chan.ph1_state(Chan.ind0, ph1);
      a = ph.v1;
      i = ph.v2;
      hp_ind = Chan.hp1_map[Chan.ind0][Hash(i, a, 0)];
      for(int ph2 = 0; ph2 < nph0; ++ph2){
	ind1 = hp_ind * nph0 + ph2;
	energy1 += S1.t2[ph1] * S1.t2[ph2] * Ints.Vhhpp.V_2_3[chan_ind1 + ind1];
      }
    }
    energy += 0.5 * energy1;
  }
  if(PAR.HF == 0){
    for(int chan = 0; chan < Chan.size3; ++chan){
      np = Chan.np[chan];
      nh = Chan.nh[chan];
      for(int p = 0; p < np; ++p){
	a = Chan.p_state(chan, p).v1;
	for(int h = 0; h < nh; ++h){
	  i = Chan.h_state(chan, h).v1;
	  ind1 = Chan.hp1_map[Chan.ind0][Hash(i, a, 0)];
	  ind2 = Chan.ph1_map[Chan.ind0][Hash(a, i, 0)];
	  energy += Ints.Fmatrix.hp_2[ind1] * S1.t2[ind2];
	}
      }
    }
  }

  this->dE = energy;
}

void Amplitudes::Copy(Amplitudes &Amps)
{
  if(PAR.approx == "doubles"){
    this->D1.Copy(Amps.D1);
  }
  else if(PAR.approx == "singles"){
    this->D1.Copy(Amps.D1);
    this->S1.Copy(Amps.S1);
  }
}

void Amplitudes::Copy(Amplitudes &Amps, double *vec)
{
  if(PAR.approx == "doubles"){
    this->D1.Copy(Amps.D1, vec);
  }
  else if(PAR.approx == "singles"){
    this->D1.Copy(Amps.D1, vec);
    this->S1.Copy(Amps.S1, vec, Amps.D1.T1_length);
  }
}

void Amplitudes::Build(Channels &Chan)
{
  if(PAR.approx == "doubles"){
    this->D1.Build(Chan);
  }
  else if(PAR.approx == "singles"){
    this->D1.Build(Chan);
    this->S1.Build(Chan);
  }
}

void Amplitudes::Build(Channels &Chan, Amplitudes &Amps)
{
  if(PAR.approx == "doubles"){
    this->D1.Build(Chan, Amps.D1);
  }
  else if(PAR.approx == "singles"){
    this->D1.Build(Chan, Amps.D1);
    this->S1.Build(Chan, Amps.S1);
  }
}

void Amplitudes::Delete()
{
  if(PAR.approx == "doubles"){
    this->D1.Delete();
  }
  else if(PAR.approx == "singles"){
    this->D1.Delete();
    this->S1.Delete();
  }
}

void Amplitudes::Zero(bool t1flag)
{
  if(PAR.approx == "doubles"){
    this->D1.Zero(t1flag);
  }
  else if(PAR.approx == "singles"){
    this->D1.Zero(t1flag);
    this->S1.Zero(t1flag);
  }
}

void Amplitudes::Set_T(double *vec)
{
  if(PAR.approx == "doubles"){
    this->D1.Set_T(vec);
  }
  else if(PAR.approx == "singles"){
    this->D1.Set_T(vec);
    this->S1.Set_T(vec, D1.T1_length);
  }
}

void Doubles::Copy(Doubles &D1)
{
  for(int ind = 0; ind < this->T1_length; ++ind){ this->T1[ind] = D1.T1[ind]; }
}

void Doubles::Copy(Doubles &D1, double *vec)
{
  for(int ind = 0; ind < this->T1_length; ++ind){
    this->T1[ind] = D1.T1[ind];
    vec[ind] = D1.T1[ind];
  }
}

void Doubles::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int a, b, i, j, npp, nhh;
  two_body pp, hh;
  four_body abij, abij_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J21_map, *J22_map, *J23_map, *J24_map;

  length = 0;
  this->T1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->T1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhh[chan1];
  }
  this->T1 = new double[length];
  //this->Evec_chan = new int[4 * length];
  this->Evec_ind = new int[4 * length];
  this->Evec_fac = new double[4 * length];
  for(ind = 0; ind < length; ++ind){
    this->T1[ind] = 0.0;
  }
  this->T1_length = length;

  length = 0;
  this->T2_1_index = new int[Chan.size2];
  this->T2_2_index = new int[Chan.size2];
  this->T2_3_index = new int[Chan.size2];
  this->T2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->T2_1_index[chan2] = length;
    this->T2_2_index[chan2] = length;
    this->T2_3_index[chan2] = length;
    this->T2_4_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.nhp1[chan2];
  }
  this->T2_1_length = length;
  this->T2_2_length = length;
  this->T2_3_length = length;
  this->T2_4_length = length;
  this->map21_num = new int[length];
  this->map22_num = new int[length];
  this->map23_num = new int[length];
  this->map24_num = new int[length];
  this->map21_index = new int[length];
  this->map22_index = new int[length];
  this->map23_index = new int[length];
  this->map24_index = new int[length];
  J21_map = new std::unordered_map<int,int>[length];
  J22_map = new std::unordered_map<int,int>[length];
  J23_map = new std::unordered_map<int,int>[length];
  J24_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){
    this->map21_num[ind] = 0;
    this->map22_num[ind] = 0;
    this->map23_num[ind] = 0;
    this->map24_num[ind] = 0;
  }

  length = 0;
  this->T3_1_index = new int[Chan.size3];
  this->T3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->T3_1_index[chan3] = length;
    this->T3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhp[chan3];
  }
  this->T3_1_length = length;
  this->T3_2_length = length;

  length = 0;
  this->T3_3_index = new int[Chan.size3];
  this->T3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->T3_3_index[chan3] = length;
    this->T3_4_index[chan3] = length;
    length += Chan.npph[chan3] * Chan.nh[chan3];
  }
  this->T3_3_length = length;
  this->T3_4_length = length;

  fb_ind = new four_body[this->T1_length];
  fb_j = new four_body[this->T1_length];
  J = new int[this->T1_length];
  length = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      pp = Chan.pp_state(chan1, pp_ind);
      a = pp.v1;
      b = pp.v2;
      abij.v1 = a;
      abij.v2 = b;
      abij_j.v1 = SPB.qnums[a].j;
      abij_j.v2 = SPB.qnums[b].j;
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	i = hh.v1;
	j = hh.v2;
	abij.v3 = i;
	abij.v4 = j;
	abij_j.v3 = SPB.qnums[i].j;
	abij_j.v4 = SPB.qnums[j].j;
	fb_ind[length] = abij;
	fb_j[length] = abij_j;
	J[length] = Chan.qnums1[chan1].j;
	Map_Count_1_to_21(this->map21_num, this->T2_1_index, Chan.ph1_map, Chan.hp1_map, Chan.nhp1, abij, J[length], J21_map);
	Map_Count_1_to_22(this->map22_num, this->T2_2_index, Chan.ph1_map, Chan.hp1_map, Chan.nhp1, abij, J[length], J22_map);
	Map_Count_1_to_23(this->map23_num, this->T2_3_index, Chan.ph1_map, Chan.hp1_map, Chan.nhp1, abij, J[length], J23_map);
	Map_Count_1_to_24(this->map24_num, this->T2_4_index, Chan.ph1_map, Chan.hp1_map, Chan.nhp1, abij, J[length], J24_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < this->T2_1_length; ++ind){
    this->map21_index[ind] = length;
    length += this->map21_num[ind];
  }
  this->map21_ind = new int[length];
  this->map21_fac1 = new double[length];
  this->map21_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map21_ind[ind] = 0;
    this->map21_fac1[ind] = 0.0;
    this->map21_fac2[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->T2_2_length; ++ind){
    this->map22_index[ind] = length;
    length += this->map22_num[ind];
  }
  this->map22_ind = new int[length];
  this->map22_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map22_ind[ind] = 0;
    this->map22_fac1[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->T2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  this->map23_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac1[ind] = 0.0;
    this->map23_fac2[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->T2_4_length; ++ind){
    this->map24_index[ind] = length;
    length += this->map24_num[ind];
  }
  this->map24_ind = new int[length];
  this->map24_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map24_ind[ind] = 0;
    this->map24_fac1[ind] = 0.0;
  }

  this->map31_ind = new int[T3_1_length];
  this->map31_fac1 = new double[T3_1_length];
  this->map31_fac2 = new double[T3_1_length];
  for(int ind = 0; ind < this->T3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac1[ind] = 0.0;
    this->map31_fac2[ind] = 0.0;
  }

  this->map32_ind = new int[T3_2_length];
  this->map32_fac1 = new double[T3_2_length];
  this->map32_fac2 = new double[T3_2_length];
  for(int ind = 0; ind < this->T3_2_length; ++ind){
    this->map32_ind[ind] = 0;
    this->map32_fac1[ind] = 0.0;
    this->map32_fac2[ind] = 0.0;
  }

  this->map33_ind = new int[T3_3_length];
  this->map33_fac1 = new double[T3_3_length];
  this->map33_fac2 = new double[T3_3_length];
  for(int ind = 0; ind < this->T3_3_length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac1[ind] = 0.0;
    this->map33_fac2[ind] = 0.0;
  }

  this->map34_ind = new int[T3_4_length];
  this->map34_fac1 = new double[T3_4_length];
  this->map34_fac2 = new double[T3_4_length];
  for(int ind = 0; ind < this->T3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac1[ind] = 0.0;
    this->map34_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static) private(a, b, i, j)
  for(int pphh = 0; pphh < this->T1_length; ++pphh){
    Map_1_to_21(this->map21_ind,this->map21_fac1,this->map21_index, this->T2_1_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J21_map,1);
    Map_1_to_21(this->map21_ind,this->map21_fac2,this->map21_index, this->T2_1_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J21_map,2);
    Map_1_to_22(this->map22_ind,this->map22_fac1,this->map22_index, this->T2_2_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J22_map,1);
    Map_1_to_23(this->map23_ind,this->map23_fac1,this->map23_index, this->T2_3_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J23_map,1);
    Map_1_to_23(this->map23_ind,this->map23_fac2,this->map23_index, this->T2_3_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J23_map,2);
    Map_1_to_24(this->map24_ind,this->map24_fac1,this->map24_index, this->T2_4_index,pphh, Chan.ph1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J24_map,1);
    Map_1_to_31(this->map31_ind,this->map31_fac1,this->T3_1_index,pphh, Chan.p_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 1);
    Map_1_to_31(this->map31_ind,this->map31_fac2,this->T3_1_index,pphh, Chan.p_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 2);
    Map_1_to_32(this->map32_ind,this->map32_fac1,this->T3_2_index,pphh, Chan.p_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 1);
    Map_1_to_32(this->map32_ind,this->map32_fac2,this->T3_2_index,pphh, Chan.p_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 2);
    Map_1_to_33(this->map33_ind,this->map33_fac1,this->T3_3_index,pphh, Chan.pph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
    Map_1_to_33(this->map33_ind,this->map33_fac2,this->T3_3_index,pphh, Chan.pph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 2);
    Map_1_to_34(this->map34_ind,this->map34_fac1,this->T3_4_index,pphh, Chan.pph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
    Map_1_to_34(this->map34_ind,this->map34_fac2,this->T3_4_index,pphh, Chan.pph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 2);

    a = fb_ind[pphh].v1;
    b = fb_ind[pphh].v2;
    i = fb_ind[pphh].v3;
    j = fb_ind[pphh].v4;
    this->Evec_ind[4*pphh] = Chan.hh1_map[Chan.ind0][Hash(i, i, 0)];
    this->Evec_fac[4*pphh] = 1.0 / std::sqrt(fb_j[pphh].v3 + 1.0);
    this->Evec_ind[4*pphh + 1] = Chan.hh1_map[Chan.ind0][Hash(j, j, 0)];
    this->Evec_fac[4*pphh + 1] = 1.0 / std::sqrt(fb_j[pphh].v4 + 1.0);
    this->Evec_ind[4*pphh + 2] = Chan.pp1_map[Chan.ind0][Hash(a, a, 0)];
    this->Evec_fac[4*pphh + 2] = 1.0 / std::sqrt(fb_j[pphh].v1 + 1.0);
    this->Evec_ind[4*pphh + 3] = Chan.pp1_map[Chan.ind0][Hash(b, b, 0)];
    this->Evec_fac[4*pphh + 3] = 1.0 / std::sqrt(fb_j[pphh].v2 + 1.0);
  }
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J21_map;
  delete[] J22_map;
  delete[] J23_map;
  delete[] J24_map;
}

void Doubles::Build(Channels &Chan, Doubles &D1)
{
  int chan1, chan2, chan3, ind, length;
  this->T1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){ this->T1_index[chan1] = D1.T1_index[chan1]; }
  this->T1 = new double[D1.T1_length];
  this->Evec_ind = new int[4 * D1.T1_length];
  this->Evec_fac = new double[4 * D1.T1_length];
  for(ind = 0; ind < D1.T1_length; ++ind){ this->T1[ind] = D1.T1[ind]; }
  for(ind = 0; ind < 4 * D1.T1_length; ++ind){
    this->Evec_ind[ind] = D1.Evec_ind[ind];
    this->Evec_fac[ind] = D1.Evec_fac[ind];
  }
  this->T1_length = D1.T1_length;

  this->T2_1_index = new int[Chan.size2];
  this->T2_2_index = new int[Chan.size2];
  this->T2_3_index = new int[Chan.size2];
  this->T2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->T2_1_index[chan2] = D1.T2_1_index[chan2];
    this->T2_2_index[chan2] = D1.T2_2_index[chan2];
    this->T2_3_index[chan2] = D1.T2_3_index[chan2];
    this->T2_4_index[chan2] = D1.T2_4_index[chan2];
  }
  this->T2_1_length = D1.T2_1_length;
  this->T2_2_length = D1.T2_2_length;
  this->T2_3_length = D1.T2_3_length;
  this->T2_4_length = D1.T2_4_length;

  this->T3_1_index = new int[Chan.size3];
  this->T3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->T3_1_index[chan3] = D1.T3_1_index[chan3];
    this->T3_2_index[chan3] = D1.T3_2_index[chan3];
  }
  this->T3_1_length = D1.T3_1_length;
  this->T3_2_length = D1.T3_2_length;

  this->T3_3_index = new int[Chan.size3];
  this->T3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->T3_3_index[chan3] = D1.T3_3_index[chan3];
    this->T3_4_index[chan3] = D1.T3_4_index[chan3];
  }
  this->T3_3_length = D1.T3_3_length;
  this->T3_4_length = D1.T3_4_length;

  this->map21_num = new int[D1.T2_1_length];
  this->map21_index = new int[D1.T2_1_length];
  this->map22_num = new int[D1.T2_2_length];
  this->map22_index = new int[D1.T2_2_length];
  this->map23_num = new int[D1.T2_3_length];
  this->map23_index = new int[D1.T2_3_length];
  this->map24_num = new int[D1.T2_4_length];
  this->map24_index = new int[D1.T2_4_length];

  length = 0;
  for(ind = 0; ind < D1.T2_1_length; ++ind){
    this->map21_num[ind] = D1.map21_num[ind];
    this->map21_index[ind] = D1.map21_index[ind];
    length += D1.map21_num[ind];
  }
  this->map21_ind = new int[length];
  this->map21_fac1 = new double[length];
  this->map21_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map21_ind[ind] = D1.map21_ind[ind];
    this->map21_fac1[ind] = D1.map21_fac1[ind];
    this->map21_fac2[ind] = D1.map21_fac2[ind];
  }

  length = 0;
  for(ind = 0; ind < D1.T2_2_length; ++ind){
    this->map22_num[ind] = D1.map22_num[ind];
    this->map22_index[ind] = D1.map22_index[ind];
    length += D1.map22_num[ind];
  }
  this->map22_ind = new int[length];
  this->map22_fac1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map22_ind[ind] = D1.map22_ind[ind];
    this->map22_fac1[ind] = D1.map22_fac1[ind];
  }

  length = 0;
  for(ind = 0; ind < D1.T2_3_length; ++ind){
    this->map23_num[ind] = D1.map23_num[ind];
    this->map23_index[ind] = D1.map23_index[ind];
    length += D1.map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  this->map23_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = D1.map23_ind[ind];
    this->map23_fac1[ind] = D1.map23_fac1[ind];
    this->map23_fac2[ind] = D1.map23_fac2[ind];
  }

  length = 0;
  for(ind = 0; ind < D1.T2_4_length; ++ind){
    this->map24_num[ind] = D1.map24_num[ind];
    this->map24_index[ind] = D1.map24_index[ind];
    length += D1.map24_num[ind];
  }
  this->map24_ind = new int[length];
  this->map24_fac1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map24_ind[ind] = D1.map24_ind[ind];
    this->map24_fac1[ind] = D1.map24_fac1[ind];
  }

  this->map31_ind = new int[D1.T3_1_length];
  this->map31_fac1 = new double[D1.T3_1_length];
  this->map31_fac2 = new double[D1.T3_1_length];
  for(ind = 0; ind < D1.T3_1_length; ++ind){
    this->map31_ind[ind] = D1.map31_ind[ind];
    this->map31_fac1[ind] = D1.map31_fac1[ind];
    this->map31_fac2[ind] = D1.map31_fac2[ind];
  }

  this->map32_ind = new int[D1.T3_2_length];
  this->map32_fac1 = new double[D1.T3_2_length];
  this->map32_fac2 = new double[D1.T3_2_length];
  for(ind = 0; ind < D1.T3_2_length; ++ind){
    this->map32_ind[ind] = D1.map32_ind[ind];
    this->map32_fac1[ind] = D1.map32_fac1[ind];
    this->map32_fac2[ind] = D1.map32_fac2[ind];
  }

  this->map33_ind = new int[D1.T3_3_length];
  this->map33_fac1 = new double[D1.T3_3_length];
  this->map33_fac2 = new double[D1.T3_3_length];
  for(ind = 0; ind < D1.T3_3_length; ++ind){
    this->map33_ind[ind] = D1.map33_ind[ind];
    this->map33_fac1[ind] = D1.map33_fac1[ind];
    this->map33_fac2[ind] = D1.map33_fac2[ind];
  }

  this->map34_ind = new int[D1.T3_4_length];
  this->map34_fac1 = new double[D1.T3_4_length];
  this->map34_fac2 = new double[D1.T3_4_length];
  for(ind = 0; ind < D1.T3_3_length; ++ind){
    this->map34_ind[ind] = D1.map34_ind[ind];
    this->map34_fac1[ind] = D1.map34_fac1[ind];
    this->map34_fac2[ind] = D1.map34_fac2[ind];
  }
}

void Doubles::Delete()
{
  delete[] this->T1;
  //delete[] this->Evec_chan;
  delete[] this->Evec_ind;
  delete[] this->Evec_fac;

  delete[] this->T1_index;
  delete[] this->T2_1_index;
  delete[] this->T2_2_index;
  delete[] this->T2_3_index;
  delete[] this->T2_4_index;
  delete[] this->T3_1_index;
  delete[] this->T3_2_index;
  delete[] this->T3_3_index;
  delete[] this->T3_4_index;

  delete[] this->map21_index;
  delete[] this->map21_num;
  delete[] this->map21_ind;
  delete[] this->map21_fac1;
  delete[] this->map21_fac2;

  delete[] this->map22_index;
  delete[] this->map22_num;
  delete[] this->map22_ind;
  delete[] this->map22_fac1;

  delete[] this->map23_index;
  delete[] this->map23_ind;
  delete[] this->map23_num;
  delete[] this->map23_fac1;
  delete[] this->map23_fac2;

  delete[] this->map24_index;
  delete[] this->map24_num;
  delete[] this->map24_ind;
  delete[] this->map24_fac1;

  delete[] this->map31_ind;
  delete[] this->map31_fac1;
  delete[] this->map31_fac2;

  delete[] this->map32_ind;
  delete[] this->map32_fac1;
  delete[] this->map32_fac2;

  delete[] this->map33_ind;
  delete[] this->map33_fac1;
  delete[] this->map33_fac2;

  delete[] this->map34_ind;
  delete[] this->map34_fac1;
  delete[] this->map34_fac2;
}

void Doubles::Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->T1_length; ++ind){ this->T1[ind] = 0.0; } }
}

void Doubles::Set_T2_1(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    T2[ind] = 0.0;
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      T2[ind] += this->map21_fac2[index] * this->T1[index1];
    }
  }
}

void Doubles::Set_T2_3(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    T2[ind] = 0.0;
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      T2[ind] += this->map23_fac2[index] * this->T1[index1];
    }
  }
}

void Doubles::Gather_T2_1(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      this->T1[index1] += this->map21_fac1[index] * T2[ind];
      /*if(index1 == 113){
	std::cout << "T2_1[" << index1 << "]: += " << this->map21_fac1[index] << " * " << T2[ind] << std::endl;
	}*/
    }
    T2[ind] = 0.0;
  }
}

void Doubles::Gather_T2_2(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map22_num[ind + offset]; ++n){
      index = this->map22_index[ind + offset] + n;
      index1 = this->map22_ind[index];
      this->T1[index1] += this->map22_fac1[index] * T2[ind];
      /*if(index1 == 113){
	std::cout << "T2_2[" << index1 << "]: += " << this->map22_fac1[index] << " * " << T2[ind] << std::endl;
	}*/
    }
    T2[ind] = 0.0;
  }
}

void Doubles::Gather_T2_3(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      this->T1[index1] += this->map23_fac1[index] * T2[ind];
      /*if(index1 == 113){
	std::cout << "T2_3[" << index1 << "]: += " << this->map23_fac1[index] << " * " << T2[ind] << std::endl;
	}*/
    }
    T2[ind] = 0.0;
  }
}

void Doubles::Gather_T2_4(double *T2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map24_num[ind + offset]; ++n){
      index = this->map24_index[ind + offset] + n;
      index1 =this->map24_ind[index];
      this->T1[index1] += this->map24_fac1[index] * T2[ind];
      /*if(index1 == 113){
	std::cout << "T2_4[" << index1 << "]: += " << this->map24_fac1[index] << " * " << T2[ind] << std::endl;
	}*/
    }
    T2[ind] = 0.0;
  }
}

void Doubles::Set_T3_1(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    T3[ind] = this->map31_fac2[ind + offset] * this->T1[index];
  }
}

void Doubles::Set_T3_2(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    T3[ind] = this->map32_fac2[ind + offset] * this->T1[index];
  }
}

void Doubles::Set_T3_3(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    T3[ind] = this->map33_fac2[ind + offset] * this->T1[index];
  }
}

void Doubles::Set_T3_4(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    T3[ind] = this->map34_fac2[ind + offset] * this->T1[index];
  }
}

void Doubles::Gather_T3_1(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->T1[index] += this->map31_fac1[ind + offset] * T3[ind];
    /*if(index == 113){
      std::cout << "T3_1[" << index << "]: += " << this->map31_fac1[ind + offset] << " * " << T3[ind] << std::endl;
      }*/
    T3[ind] = 0.0;
  }
}

void Doubles::Gather_T3_2(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    this->T1[index] += this->map32_fac1[ind + offset] * T3[ind];
    /*if(index == 113){
      std::cout << "T3_2[" << index << "]: += " << this->map32_fac1[ind + offset] << " * " << T3[ind] << std::endl;
      }*/
    T3[ind] = 0.0;
  }
}

void Doubles::Gather_T3_3(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    this->T1[index] += this->map33_fac1[ind + offset] * T3[ind];
    /*if(index == 113){
      std::cout << "T3_3[" << index << "]: += " << this->map33_fac1[ind + offset] << " * " << T3[ind] << std::endl;
      }*/
    T3[ind] = 0.0;
  }
}

void Doubles::Gather_T3_4(double *T3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    this->T1[index] += this->map34_fac1[ind + offset] * T3[ind];
    /*if(index == 113){
      std::cout << "T3_4[" << index << "]: += " << this->map34_fac1[ind + offset] << " * " << T3[ind] << std::endl;
      }*/
    T3[ind] = 0.0;
  }
}

void Doubles::Set_T(double *vec)
{
  for(int ind = 0; ind < this->T1_length; ++ind){ T1[ind] = vec[ind]; }
}

void Singles::Copy(Singles &S1)
{
  for(int ind = 0; ind < this->t2_length; ++ind){ this->t2[ind] = S1.t2[ind]; }
}

void Singles::Copy(Singles &S1, double *vec, int offset)
{
  for(int ind = 0; ind < this->t2_length; ++ind){
    this->t2[ind] = S1.t2[ind];
    vec[offset + ind] = S1.t2[ind];
  }
}

void Singles::Build(Channels &Chan)
{
  int chan3, ind, length;
  int a, i, nph1;
  two_body ph1;
  two_body ai, ai_j;
  two_body *tb_ind;
  two_body *tb_j;

  length = Chan.nph1[Chan.ind0];
  this->t2 = new double[length];
  this->evec_ind = new int[2 * length];
  this->evec_fac = new double[2 * length];
  for(ind = 0; ind < length; ++ind){
    this->t2[ind] = 0.0;
  }
  this->t2_length = length;

  length = 0;
  this->t3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->t3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nh[chan3];
  }
  this->t3_length = length;

  tb_ind = new two_body[this->t2_length];
  tb_j = new two_body[this->t2_length];
  length = 0;
  nph1 = Chan.nph1[Chan.ind0];
  for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
    ph1 = Chan.ph1_state(Chan.ind0, ph1_ind);
    a = ph1.v1;
    i = ph1.v2;
    ai.v1 = a;
    ai.v2 = i;
    ai_j.v1 = SPB.qnums[a].j;
    ai_j.v2 = SPB.qnums[i].j;
    tb_ind[length] = ai;
    tb_j[length] = ai_j;
    ++length;
  }

  this->map3_ind = new int[t3_length];
  this->map3_fac1 = new double[t3_length];
  this->map3_fac2 = new double[t3_length];
  for(int ind = 0; ind < t3_length; ++ind){
    this->map3_ind[ind] = 0;
    this->map3_fac1[ind] = 0.0;
    this->map3_fac2[ind] = 0.0;
  }

  for(int ph = 0; ph < this->t2_length; ++ph){
    Map_2_to_3(this->map3_ind, this->map3_fac1, t3_index, ph, Chan.p_map, Chan.h_map, Chan.nh, tb_ind, tb_j, 1);
    Map_2_to_3(this->map3_ind, this->map3_fac2, t3_index, ph, Chan.p_map, Chan.h_map, Chan.nh, tb_ind, tb_j, 2);

    a = tb_ind[ph].v1;
    i = tb_ind[ph].v2;
    this->evec_ind[2*ph] = Chan.hh1_map[Chan.ind0][Hash(i, i, 0)];
    this->evec_fac[2*ph] = 1.0 / std::sqrt(tb_j[ph].v2 + 1.0);
    this->evec_ind[2*ph + 1] = Chan.pp1_map[Chan.ind0][Hash(a, a, 0)];
    this->evec_fac[2*ph + 1] = 1.0 / std::sqrt(tb_j[ph].v1 + 1.0);
  }
  delete[] tb_ind;
  delete[] tb_j;
}

void Singles::Build(Channels &Chan, Singles &S1)
{
  int chan3, ind, length;
  length = Chan.nph1[Chan.ind0];
  this->t2 = new double[length];
  this->evec_ind = new int[2 * length];
  this->evec_fac = new double[2 * length];
  for(ind = 0; ind < length; ++ind){
    this->t2[ind] = S1.t2[ind];
  }
  this->t2_length = length;

  this->t3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->t3_index[chan3] = S1.t3_index[chan3];
  }
  this->t3_length = S1.t3_length;
  this->map3_ind = new int[S1.t3_length];
  this->map3_fac1 = new double[S1.t3_length];
  this->map3_fac2 = new double[S1.t3_length];
  for(ind = 0; ind < S1.t3_length; ++ind){
    this->map3_ind[ind] = S1.map3_ind[ind];
    this->map3_fac1[ind] = S1.map3_fac1[ind];
    this->map3_fac2[ind] = S1.map3_fac2[ind];
  }
  for(ind = 0; ind < 2 * S1.t2_length; ++ind){
    this->evec_ind[ind] = S1.evec_ind[ind];
    this->evec_fac[ind] = S1.evec_fac[ind];
  }
}

void Singles::Delete()
{
  delete[] this->t2;
  delete[] this->evec_ind;
  delete[] this->evec_fac;
  delete[] this->t3_index;

  delete[] this->map3_ind;
  delete[] this->map3_fac1;
  delete[] this->map3_fac2;
}

void Singles::Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->t2_length; ++ind){ this->t2[ind] = 0.0; } }
}

void Singles::Gather_t3(double *t3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->t2[index] += this->map3_fac1[ind + offset] * t3[ind];
    t3[ind] = 0.0;
  }
}

void Singles::Set_t3(double *t3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    t3[ind] = this->map3_fac2[ind + offset] * this->t2[index];
  }
}

void Singles::Set_T(double *vec, int offset)
{
  for(int ind = 0; ind < this->t2_length; ++ind){
    t2[ind] = vec[ind + offset];
  }
}

void CC_Algorithm(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps)
{
  //double mixmax = 0.75;
  //double mixmin = 0.01;
  double mixmax = 1.0;
  double mixmin = 0.5;
  double mix;
  double width = 1.0;
  int Rand_count = 0;

  double error = 1e20, error2 = 1e20;
  int ind = 0;
  double CCoutE;
  Amplitudes Amps0;
  Amplitudes Amps2;
  Amplitudes tempAmps;
  Amps0.Build(Chan, Amps);
  Amps2.Build(Chan, Amps);
  tempAmps.Build(Chan, Amps);

  std::cout << std::setprecision(10);

  /// For DIIS ///
  DIIS DIIS;
  int maxl = 5;
  double DIIS_start = 0.01;
  int size = Amps.D1.T1_length;
  if(PAR.approx == "singles"){ size += Amps.S1.t2_length; }
  double *vec = new double[size];
  double *vec0 = new double[size];
  DIIS.Build(size, maxl);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*four_body *Pairlist0 = new four_body[1];
  four_body *Pairlist1 = new four_body[8];
  four_body *Pairlist2 = new four_body[6];
  four_body state4, state4_1, state4_2;
  double ME1, ME2, tempen, Eref;
  int index1, index2, index, chan1, ind0 = 8;
  State tb1, tb2;

  state4.v1 = 0, state4.v2 = 1, state4.v3 = 2, state4.v4 = 3;
  Pairlist0[0] = state4;
  state4.v1 = 2, state4.v2 = 3, state4.v3 = 4, state4.v4 = 5;
  Pairlist1[0] = state4;
  state4.v1 = 2, state4.v2 = 3, state4.v3 = 6, state4.v4 = 7;
  Pairlist1[1] = state4;
  state4.v1 = 0, state4.v2 = 1, state4.v3 = 4, state4.v4 = 5;
  Pairlist1[2] = state4;
  state4.v1 = 2, state4.v2 = 3, state4.v3 = 8, state4.v4 = 9;
  Pairlist1[3] = state4;
  state4.v1 = 0, state4.v2 = 1, state4.v3 = 6, state4.v4 = 7;
  Pairlist1[4] = state4;
  state4.v1 = 2, state4.v2 = 3, state4.v3 = 10, state4.v4 = 11;
  Pairlist1[5] = state4;
  state4.v1 = 0, state4.v2 = 1, state4.v3 = 8, state4.v4 = 9;
  Pairlist1[6] = state4;
  state4.v1 = 0, state4.v2 = 1, state4.v3 = 10, state4.v4 = 11;
  Pairlist1[7] = state4;
  state4.v1 = 4, state4.v2 = 5, state4.v3 = 6, state4.v4 = 7;
  Pairlist2[0] = state4;
  state4.v1 = 4, state4.v2 = 5, state4.v3 = 8, state4.v4 = 9;
  Pairlist2[1] = state4;
  state4.v1 = 4, state4.v2 = 5, state4.v3 = 10, state4.v4 = 11;
  Pairlist2[2] = state4;
  state4.v1 = 6, state4.v2 = 7, state4.v3 = 8, state4.v4 = 9;
  Pairlist2[3] = state4;
  state4.v1 = 6, state4.v2 = 7, state4.v3 = 10, state4.v4 = 11;
  Pairlist2[4] = state4;
  state4.v1 = 8, state4.v2 = 9, state4.v3 = 10, state4.v4 = 11;
  Pairlist2[5] = state4;*/
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Initialize Amplitudes //
  mix = mixmin;
  while(error > 1e-12 && ind < 2000){
    Update_Amps(Chan, Ints, Eff_Ints, Amps, tempAmps);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*if( ind == 0 || ind == 1 || ind == 9 ){
      Eff_Ints.Update_3(Chan, Ints, Amps);
      ///////////
      Eref = 0.0;
      state4_1 = Pairlist0[0];
      //////
      index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v1, state4_1.v1, 0)];
      Eref += Eff_Ints.Xhh.X_2[index];
      index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v2, state4_1.v2, 0)];
      Eref += Eff_Ints.Xhh.X_2[index];
      index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v3, state4_1.v3, 0)];
      Eref += Eff_Ints.Xhh.X_2[index];
      index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v4, state4_1.v4, 0)];
      Eref += Eff_Ints.Xhh.X_2[index];
      //////
      plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
      chan1 = Ind_Chan1(tb1);
      plus(tb2, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
      if( Ind_Chan1(tb2) == chan1 ){
	index1 = Chan.hh_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
	index2 = Chan.hh_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
	index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
	Eref -= Eff_Ints.Xhhhh.X_1[index];
	index1 = Chan.hh_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
	index2 = Chan.hh_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
	index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
        Eref -= Eff_Ints.Xhhhh.X_1[index];
      }
      ///////////

      //std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      for(int ind1 = 0; ind1 < 15; ++ind1){
	//state4_1 = Pairlist[ind1];
	for(int ind2 = 0; ind2 < 15; ++ind2){
	  //state4_2 = Pairlist[ind2];

	  ME1 = 0.0;
	  ME2 = 0.0;

	  if( ind1 == ind2 ){
	    ME1 += Eref;
	    if( ind1 > 0 ){
	      if( ind1 > ind0 ){  // <4p4h|4p4h>
		state4_1 = Pairlist0[0];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v1, state4_1.v1, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v2, state4_1.v2, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v3, state4_1.v3, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v4, state4_1.v4, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		state4_1 = Pairlist2[ind1 - ind0 - 1];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v1, state4_1.v1, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v2, state4_1.v2, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v3, state4_1.v3, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v4, state4_1.v4, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
	      }
	      else{             // <2p2h|2p2h>
		state4_1 = Pairlist1[ind1 - 1];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v1, state4_1.v1, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		index = Chan.hh1_map[Chan.ind0][Hash(state4_1.v2, state4_1.v2, 0)];
		ME1 -= Eff_Ints.Xhh.X_2[index];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v3, state4_1.v3, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
		index = Chan.pp1_map[Chan.ind0][Hash(state4_1.v4, state4_1.v4, 0)];
		ME1 += Eff_Ints.Xpp.X_2[index];
	      }
	    }
	  }

	  

	  if( ind1 > 0 ){
	    if( ind1 > ind0 ){
	      if( ind2 > 0 ){
		if( ind2 > ind0 ){ // <4p4h|4p4h>

		  state4_1 = Pairlist2[ind1 - ind0 - 1];
		  state4_2 = Pairlist2[ind2 - ind0 - 1];
		  plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
		  chan1 = Ind_Chan1(tb1);
		  index1 = Chan.pp_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
		  plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		  if( Ind_Chan1(tb2) == chan1 && state4_1.v3 == state4_2.v3 && state4_1.v4 == state4_2.v4 ){
		    index2 = Chan.pp_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
		    index = Eff_Ints.Xpppp.X_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		    ME2 += Eff_Ints.Xpppp.X_1[index];
		  }
		  plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		  if( Ind_Chan1(tb2) == chan1 && state4_1.v3 == state4_2.v1 && state4_1.v4 == state4_2.v2 ){
		    index2 = Chan.pp_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
		    index = Eff_Ints.Xpppp.X_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		    ME2 += Eff_Ints.Xpppp.X_1[index];
		  }
		  plus(tb1, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
		  chan1 = Ind_Chan1(tb1);
		  index1 = Chan.pp_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
		  plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		  if( Ind_Chan1(tb2) == chan1 && state4_1.v1 == state4_2.v3 && state4_1.v2 == state4_2.v4 ){
		    index2 = Chan.pp_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
		    index = Eff_Ints.Xpppp.X_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		    ME2 += Eff_Ints.Xpppp.X_1[index];
		  }
		  plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		  if( Ind_Chan1(tb2) == chan1 && state4_1.v1 == state4_2.v1 && state4_1.v2 == state4_2.v2 ){
		    index2 = Chan.pp_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
		    index = Eff_Ints.Xpppp.X_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		    ME2 += Eff_Ints.Xpppp.X_1[index];
		  }

		  if( ind1 == ind2 ){
		    state4_1 = Pairlist0[0];
		    state4_2 = Pairlist0[0];
		    plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
		    chan1 = Ind_Chan1(tb1);
		    index1 = Chan.hh_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
		    plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		    if( Ind_Chan1(tb2) == chan1 && state4_1.v3 == state4_2.v3 && state4_1.v4 == state4_2.v4 ){
		      index2 = Chan.hh_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
		      index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		      ME2 += Eff_Ints.Xhhhh.X_1[index];
		    }
		    plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		    if( Ind_Chan1(tb2) == chan1 && state4_1.v3 == state4_2.v1 && state4_1.v4 == state4_2.v2 ){
		      index2 = Chan.hh_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
		      index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		      ME2 += Eff_Ints.Xhhhh.X_1[index];
		    }
		    plus(tb1, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
		    chan1 = Ind_Chan1(tb1);
		    index1 = Chan.hh_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
		    plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		    if( Ind_Chan1(tb2) == chan1 && state4_1.v1 == state4_2.v3 && state4_1.v2 == state4_2.v4 ){
		      index2 = Chan.hh_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
		      index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		      ME2 += Eff_Ints.Xhhhh.X_1[index];
		    }
		    plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		    if( Ind_Chan1(tb2) == chan1 && state4_1.v1 == state4_2.v1 && state4_1.v2 == state4_2.v2 ){
		      index2 = Chan.hh_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
		      index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		      ME2 += Eff_Ints.Xhhhh.X_1[index];
		    }
		  }

		}
		else{              // <4p4h|2p2h>

		  state4 = Pairlist0[0];
		  state4_1 = Pairlist2[ind1 - ind0 - 1];
		  state4_2 = Pairlist1[ind2 - 1];

		  if( state4.v1 == state4_2.v1 && state4.v2 == state4_2.v2 ){
		    if( state4_1.v1 == state4_2.v3 && state4_1.v2 == state4_2.v4 ){

		      plus(tb1, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.pp_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
		      plus(tb2, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.hh_map[chan1][Hash(state4.v3, state4.v4, 0)];
			index = Amps.D1.T1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
			tempen = 0.0;
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index]] * Amps.D1.Evec_fac[4*index];
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index + 1]] * Amps.D1.Evec_fac[4*index + 1];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 2]] * Amps.D1.Evec_fac[4*index + 2];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 3]] * Amps.D1.Evec_fac[4*index + 3];
			ME2 += tempAmps.D1.T1[index] - tempen*Amps.D1.T1[index];
		      }

		    }
		    if( state4_1.v3 == state4_2.v3 && state4_1.v4 == state4_2.v4 ){

		      plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.pp_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
		      plus(tb2, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.hh_map[chan1][Hash(state4.v3, state4.v4, 0)];
			index = Amps.D1.T1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
			tempen = 0.0;
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index]] * Amps.D1.Evec_fac[4*index];
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index + 1]] * Amps.D1.Evec_fac[4*index + 1];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 2]] * Amps.D1.Evec_fac[4*index + 2];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 3]] * Amps.D1.Evec_fac[4*index + 3];
			ME2 += tempAmps.D1.T1[index] - tempen*Amps.D1.T1[index];
		      }

		    }
		  }
		  if( state4.v3 == state4_2.v1 && state4.v4 == state4_2.v2 ){
		    if( state4_1.v1 == state4_2.v3 && state4_1.v2 == state4_2.v4 ){

		      plus(tb1, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.pp_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
		      plus(tb2, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
			index = Amps.D1.T1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
			tempen = 0.0;
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index]] * Amps.D1.Evec_fac[4*index];
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index + 1]] * Amps.D1.Evec_fac[4*index + 1];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 2]] * Amps.D1.Evec_fac[4*index + 2];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 3]] * Amps.D1.Evec_fac[4*index + 3];
			ME2 += tempAmps.D1.T1[index] - tempen*Amps.D1.T1[index];
		      }

		    }
		    if( state4_1.v3 == state4_2.v3 && state4_1.v4 == state4_2.v4 ){

		      plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.pp_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
		      plus(tb2, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
			index = Amps.D1.T1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
			tempen = 0.0;
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index]] * Amps.D1.Evec_fac[4*index];
			tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index + 1]] * Amps.D1.Evec_fac[4*index + 1];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 2]] * Amps.D1.Evec_fac[4*index + 2];
			tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 3]] * Amps.D1.Evec_fac[4*index + 3];
			ME2 += tempAmps.D1.T1[index] - tempen*Amps.D1.T1[index];
		      }

		    }
		  }

		}
	      }
	    }
	    else{
	      if( ind2 > 0 ){
		if( ind2 > ind0 ){ // <2p2h|4p4h>

		  state4 = Pairlist0[0];
		  state4_1 = Pairlist1[ind1 - 1];
		  state4_2 = Pairlist2[ind2 - ind0 - 1];

		  if( state4.v1 == state4_1.v1 && state4.v2 == state4_1.v2 ){
		    if( state4_2.v1 == state4_1.v3 && state4_2.v2 == state4_1.v4 ){

		      plus(tb1, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.hh_map[chan1][Hash(state4.v3, state4.v4, 0)];
		      plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.pp_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
			index = Ints.Vhhpp.V_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
			ME2 += Ints.Vhhpp.V_1[index];
		      }

		    }
		    if( state4_2.v3 == state4_1.v3 && state4_2.v4 == state4_1.v4 ){

		      plus(tb1, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.hh_map[chan1][Hash(state4.v3, state4.v4, 0)];
		      plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.pp_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
			index = Ints.Vhhpp.V_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
			ME2 += Ints.Vhhpp.V_1[index];
		      }

		    }
		  }
		  if( state4.v3 == state4_1.v1 && state4.v4 == state4_1.v2 ){
		    if( state4_2.v1 == state4_1.v3 && state4_2.v2 == state4_1.v4 ){

		      plus(tb1, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
		      plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.pp_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
			index = Ints.Vhhpp.V_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
			ME2 += Ints.Vhhpp.V_1[index];
		      }
		      
		    }
		    if( state4_2.v3 == state4_1.v3 && state4_2.v4 == state4_1.v4 ){

		      plus(tb1, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		      chan1 = Ind_Chan1(tb1);
		      index1 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
		      plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		      if( Ind_Chan1(tb2) == chan1 ){
			index2 = Chan.pp_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
			index = Ints.Vhhpp.V_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
			ME2 += Ints.Vhhpp.V_1[index];
		      }
		      
		    }
		  }

		}
		else{              // <2p2h|2p2h>		  
		  state4_1 = Pairlist1[ind1 - 1];
		  state4_2 = Pairlist1[ind2 - 1];
		  
		  if( state4_1.v1 == state4_2.v1 && state4_1.v2 == state4_2.v2 ){

		    plus(tb1, SPB.qnums[state4_1.v3], SPB.qnums[state4_1.v4]);
		    chan1 = Ind_Chan1(tb1);
		    index1 = Chan.pp_map[chan1][Hash(state4_1.v3, state4_1.v4, 0)];
		    plus(tb2, SPB.qnums[state4_2.v3], SPB.qnums[state4_2.v4]);
		    if( Ind_Chan1(tb2) == chan1 ){
		      index2 = Chan.pp_map[chan1][Hash(state4_2.v3, state4_2.v4, 0)];
		      index = Eff_Ints.Xpppp.X_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		      ME2 += Eff_Ints.Xpppp.X_1[index];
		    }
			  
		  }
		  if( state4_1.v3 == state4_2.v3 && state4_1.v4 == state4_2.v4 ){

		    plus(tb1, SPB.qnums[state4_1.v1], SPB.qnums[state4_1.v2]);
		    chan1 = Ind_Chan1(tb1);
		    index1 = Chan.hh_map[chan1][Hash(state4_1.v1, state4_1.v2, 0)];
		    plus(tb2, SPB.qnums[state4_2.v1], SPB.qnums[state4_2.v2]);
		    if( Ind_Chan1(tb2) == chan1 ){
		      index2 = Chan.hh_map[chan1][Hash(state4_2.v1, state4_2.v2, 0)];
		      index = Eff_Ints.Xhhhh.X_1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		      ME2 += Eff_Ints.Xhhhh.X_1[index];
		    }

		  }

		}
	      }
	      else{                // <2p2h|0p0h>

		state4 = Pairlist1[ind1 - 1];

		plus(tb1, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		chan1 = Ind_Chan1(tb1);
		index1 = Chan.pp_map[chan1][Hash(state4.v3, state4.v4, 0)];
		plus(tb2, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		if( Ind_Chan1(tb2) == chan1 ){
		  index2 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
		  index = Amps.D1.T1_index[chan1] + (index1 * Chan.nhh[chan1] + index2);
		  tempen = 0.0;
		  tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index]] * Amps.D1.Evec_fac[4*index];
		  tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*index + 1]] * Amps.D1.Evec_fac[4*index + 1];
		  tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 2]] * Amps.D1.Evec_fac[4*index + 2];
		  tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*index + 3]] * Amps.D1.Evec_fac[4*index + 3];
		  ME2 += tempAmps.D1.T1[index] - tempen*Amps.D1.T1[index];
		}
		
	      }
	    }
	  }
	  else{
	    if( ind2 > 0 ){
	      if( ind2 <= ind0 ){ // <0p0h|2p2h>

		state4 = Pairlist1[ind2 - 1];

		plus(tb1, SPB.qnums[state4.v1], SPB.qnums[state4.v2]);
		chan1 = Ind_Chan1(tb1);
		index1 = Chan.hh_map[chan1][Hash(state4.v1, state4.v2, 0)];
		plus(tb2, SPB.qnums[state4.v3], SPB.qnums[state4.v4]);
		if( Ind_Chan1(tb2) == chan1 ){
		  index2 = Chan.pp_map[chan1][Hash(state4.v3, state4.v4, 0)];
		  index = Ints.Vhhpp.V_1_index[chan1] + (index1 * Chan.npp[chan1] + index2);
		  ME2 += Ints.Vhhpp.V_1[index];
		}
		
	      }
	    }
	  }
	  //std::cout << std::setprecision(5) << ME1 + ME2 << ", ";
	}
	//std::cout << std::endl;
      }
      //std::cout << std::endl;
    }*/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Print_Amps(Chan, tempAmps);
    Amps2.Copy(Amps);
    Amps2.Zero(true);
    Gather_Amps(Chan, Eff_Ints, Amps2, tempAmps, mix);
    //if( ind == 2 ){ Print_Amps(Chan, Amps2); }

    CC_Error(Ints, Amps, Amps2, error);
    error /= mix;
    if( !std::isfinite(error) ){ std::cerr << std::endl << ind << ", error = " << error << " : CC Solution Diverged!! " << std::endl; exit(1); }
    if( ind > 1000 && double(Rand_count)/double(ind) > 0.9 ){ std::cerr << std::endl << ind << " : CC Solution Not Converged!!" << std::endl; exit(1); }

    if( error < error2 || error > 1.0 ){
      if( error < error2 ){ mix = std::pow(mixmax, 0.035) * std::pow(mix, 0.965); }
      Amps0.Copy(Amps, vec0);
      Amps.Copy(Amps2, vec);
      error2 = error;
      if( error < DIIS_start && mix > 0.2 * mixmax ){
	DIIS.Perform(vec, vec0, mix);
	Amps.Zero(false);
	Amps.Set_T(vec);
      }
    }
    else{
      ++Rand_count;
      if(error2 > 1.0){ width = 0.001; }
      else{ width = 0.001 * error2; }
      Random_Step(Chan, Ints, Eff_Ints, Amps0, Amps, Amps2, tempAmps, mix, width, error, error2);
    }

    Amps.get_dE(Chan, Ints);
    if( !std::isfinite(CCoutE) ){ std::cerr << std::endl << ind << ", Energy = " << Amps.dE << " : CC Solution Diverged!! " << std::endl; exit(1); }
    std::cout << "Iteration Number = " << ind << ", Energy = " << Amps.dE << ", error = " << error << ", mix = " << mix << ", ";
    std::cout << "Random = " << Rand_count << ", DIIS = " << DIIS.count << std::endl;
    ++ind;
  }

  std::cout << std::endl << std::endl;
  if( error > 1e-12 ){
    std::cout << PAR.Shells << ", " << PAR.Pshells << ", " << PAR.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". CC Solution Not Converged!!" << std::endl;
  }
  Amps0.Delete();
  Amps2.Delete();
  tempAmps.Delete();
  DIIS.Delete();
  delete[] vec;
  delete[] vec0;
}

void Update_Amps(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &tempAmps)
{
  Eff_Ints.Zero(false);
  tempAmps.Zero(false);
  //  Build Effective Hamiltonian
  Eff_Ints.Update_1(Chan, Ints, Amps);
  if(PAR.approx == "singles"){
    Eff_Ints.Update_2(Chan, Ints, Amps);
  }
  //  Solve Amplitude Equations
  Doubles_Step(Chan, Ints, Eff_Ints, Amps, tempAmps);
  if(PAR.approx == "singles"){
    Singles_Step(Chan, Ints, Eff_Ints, Amps, tempAmps);
    Doubles_Step_2(Chan, Ints, Eff_Ints, Amps, tempAmps);
  }
}

void Gather_Amps(Channels &Chan, Eff_Interactions &Eff_Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix)
{
  double tempt, tempen;
  int length = Amps.D1.T1_length;
  for(int ind = 0; ind < length; ++ind){
    tempt = Amps2.D1.T1[ind];
    tempen = 0.0;
    tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*ind]] * Amps.D1.Evec_fac[4*ind];
    tempen += Eff_Ints.Xhh.X_2[Amps.D1.Evec_ind[4*ind + 1]] * Amps.D1.Evec_fac[4*ind + 1];
    tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*ind + 2]] * Amps.D1.Evec_fac[4*ind + 2];
    tempen -= Eff_Ints.Xpp.X_2[Amps.D1.Evec_ind[4*ind + 3]] * Amps.D1.Evec_fac[4*ind + 3];
    tempt /= tempen;
    tempt = mix*tempt + (1.0-mix)*Amps.D1.T1[ind];
    Amps.D1.T1[ind] = tempt;
  }
  if(PAR.approx == "singles"){
    length = Amps.S1.t2_length;
    for(int ind = 0; ind < length; ++ind){
      tempt = Amps2.S1.t2[ind];
      tempen = 0.0;
      tempen += Eff_Ints.Xhh.X1_2[Amps.S1.evec_ind[2*ind]] * Amps.S1.evec_fac[2*ind];
      tempen -= Eff_Ints.Xpp.X_2[Amps.S1.evec_ind[2*ind + 1]] * Amps.S1.evec_fac[2*ind + 1];
      tempt /= tempen;
      tempt = mix*tempt + (1.0-mix)*Amps.S1.t2[ind];
      Amps.S1.t2[ind] = tempt;
    }
  }
}

void CC_Error(Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error)
{
  error = 0.0;
  int ind;
  int length = Amps.D1.T1_length;
  double norm = 1.0e-16;
  for(ind = 0; ind < length; ++ind){
    error += (Amps2.D1.T1[ind] - Amps.D1.T1[ind])*(Amps2.D1.T1[ind] - Amps.D1.T1[ind]);
    norm += Amps2.D1.T1[ind] * Amps2.D1.T1[ind];
  }
  if(PAR.approx == "singles"){  
    length = Amps.S1.t2_length;
    for(ind = 0; ind < length; ++ind){
      error += (Amps2.S1.t2[ind] - Amps.S1.t2[ind])*(Amps2.S1.t2[ind] - Amps.S1.t2[ind]);
      norm += Amps2.S1.t2[ind] * Amps2.S1.t2[ind];
    }
  }
  error = std::sqrt(error/norm);
}

void Print_Amps(Channels &Chan, Amplitudes &Amps)
{
  int length1 = Chan.size1;
  int length2 = Chan.nph1[Chan.ind0];
  int npp, nhh;
  int chanind, ind;//, index;
  int a, b, i, j;
  two_body pp, hh, ph;
  ind = 0;

  std::cout << std::setprecision(10) << std::endl;
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
	if( std::fabs(Amps.D1.T1[chanind + (pp1*nhh + hh1)]) > 1.0e-15 ){
	  std::cout << "! " << chanind + (pp1*nhh + hh1) << ": <" << a << "," << b << " |t| " << i << "," << j << " >^ " << Chan.qnums1[chan1].j << " = " << chanind << " " << (pp1*nhh + hh1) << ", " << Amps.D1.T1[chanind + (pp1*nhh + hh1)] << std::endl;
	}
	++ind;
      }
    }
  }
  if(PAR.approx == "singles"){
    for(int ph1 = 0; ph1 < length2; ++ph1){
      ph = Chan.ph1_state(Chan.ind0, ph1);
      a = ph.v1;
      i = ph.v2;
      if( std::fabs(Amps.S1.t2[ph1]) > 1.0e-15 ){
	std::cout << "! < " << a << " |t| " << i << " > = " << Amps.S1.t2[ph1] << std::endl;
      }
    }
  }
  std::cout << std::endl;
}

void Random_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, double &width, double &error, double &error2)
{
  int go = 0;
  int ind = 0;
  while(go == 0){
    ++ind;
    Amps.Zero(false);
    Randomize_Amps(Amps0, Amps, width);

    Eff_Ints.Zero(false);
    tempAmps.Zero(false);
    //  Build Effective Hamiltonian
    Eff_Ints.Update_1(Chan, Ints, Amps);
    if(PAR.approx == "singles"){
      Eff_Ints.Update_2(Chan, Ints, Amps);
    }
    //  Solve Amplitude Equations
    Doubles_Step(Chan, Ints, Eff_Ints, Amps, tempAmps);
    if(PAR.approx == "singles"){
      Singles_Step(Chan, Ints, Eff_Ints, Amps, tempAmps);
      Doubles_Step_2(Chan, Ints, Eff_Ints, Amps, tempAmps);
    }
    Amps2.Copy(Amps);
    Amps2.Zero(true);
    Gather_Amps(Chan, Eff_Ints, Amps2, tempAmps, mix);
    CC_Error(Ints, Amps, Amps2, error);
    error /= mix;

    std::cout << "!!  " << ind << ", " << error << " " << error2 << ", " << width << std::endl;
    if(error < error2){
      ++go;
      Amps0.Copy(Amps);
      Amps.Copy(Amps2);
      error2 = error;
    }
    width *= 0.975;
    if(width < 10-12){ width = 10e-12; }
    if(ind >= 1000){
      std::cout << PAR.Shells << ", " << PAR.Pshells << ", " << PAR.density << std::endl;
      std::cerr << std::endl << "Random Step Unsuccessful!!" << std::endl; exit(1);
    }
  }
}

void Randomize_Amps(Amplitudes &Amps0, Amplitudes &Amps, double &width)
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
  if(PAR.approx == "singles"){  
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
    Amps.D1.T1[ind] = tempt;
    //Amps.D1.Set_T(ind, tempt);
  }
  if(PAR.approx == "singles"){  
    length = Amps.S1.t2_length;
    for(ind = 0; ind < length; ++ind){
      tempt = Amps0.S1.t2[ind];
      key = std::hash<float>{}(float(fabs(tempt)));
      if(tempt > 1.0e-12){ tempt += t_map[key]; }
      else if(tempt < -1.0e-12){ tempt -= t_map[key]; }
      Amps.S1.t2[ind] = tempt;
      //Amps.S1.Set_T(ind, tempt);
    }
  }
}

void Doubles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0;
  int nh, np, npph, nhhp, nhh, npp, nhp1, nph1, length;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4, chan_ind5;
  double *T2, *TT2, *T3, *TT3, *X3;
  char N = 'N';
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    length = npp * nhh;
    chan_ind1 = Amps1.D1.T1_index[chan1];
    chan_ind2 = Ints.Vhhpp.V_1_index[chan1];
    //T1(ab|ij){ab,ij}  =  Vpphh1(ab|ij){ab,ij}
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	Amps2.D1.T1[chan_ind1 + (nhh * pp_ind + hh_ind)] += Ints.Vhhpp.V_1[chan_ind2 + (npp * hh_ind + pp_ind)]; // Vhhpp -> Vpphh
      }
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

      TT2 = new double[nph1 * nhp1];
      Amps1.D1.Set_T2_1(TT2, nph1*nhp1, chan_ind1);
      T2 = new double[nph1 * nhp1];
      for(int ind = 0; ind < nph1 * nhp1; ++ind){ T2[ind] = 0.0; }
      //T2_1(ab|ij){aj',ib'}  =  -T2_1(ac|kj){aj',kc'}.Xhphp2_1(kb|ic){kc',ib'}
      dgemm_NN(TT2, (Eff_Ints.Xhphp.X_2_1 + chan_ind5), T2, &nph1, &nhp1, &nhp1, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T2_1(T2, nph1*nhp1, chan_ind1);
      //T2_2(ab|ij){bi',ja'}  =  -T2_1(bc|ki){bi',kc'}.Xhphp2_1(ka|jc){kc',ja'}
      dgemm_NN(TT2, (Eff_Ints.Xhphp.X_2_1 + chan_ind5), T2, &nph1, &nhp1, &nhp1, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T2_2(T2, nph1*nhp1, chan_ind2);
      //T2_3(ab|ij){ai',jb'}  =   T2_1(ac|ki){ai',kc'}.Xhphp2_1(kb|jc){kc',jb'}
      dgemm_NN(TT2, (Eff_Ints.Xhphp.X_2_1 + chan_ind5), T2, &nph1, &nhp1, &nhp1, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T2_3(T2, nph1*nhp1, chan_ind3);
      //T2_4(ab|ij){bj',ia'}  =   T2_1(bc|kj){bj',kc'}.Xhphp2_1(ka|ic){kc',ia'}
      dgemm_NN(TT2, (Eff_Ints.Xhphp.X_2_1 + chan_ind5), T2, &nph1, &nhp1, &nhp1, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T2_4(T2, nph1*nhp1, chan_ind4);

      delete[] TT2;
      delete[] T2;
    }
  }

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

      X3 = new double[np * np];
      Eff_Ints.Xpp.Set_X_3_OD(X3, np*np, chan_ind3);
      TT3 = new double[np * nhhp];
      Amps1.D1.Set_T3_1(TT3, np*nhhp, chan_ind1);
      T3 = new double[np * nhhp];
      for(int ind = 0; ind < np*nhhp; ++ind){ T3[ind] = 0.0; }
      //T3_1(ab|ij){a,ijb'}  =  Xpp3_od(a|c){a,c}.T3_1(cb|ij){c,ijb'}
      dgemm_NN(X3, TT3, T3, &np, &nhhp, &np, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T3_1(T3, np*nhhp, chan_ind1);
      //T3_2(ab|ij){b,ija'}  =  -Xpp3_od(b|c){b,c}.T3_1(ca|ij){c,ija'}
      dgemm_NN(X3, TT3, T3, &np, &nhhp, &np, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T3_2(T3, np*nhhp, chan_ind2);

      delete[] X3;
      delete[] TT3;
      delete[] T3;
    }
    length = nh * npph;
    //std::cout << "!!!!!  Length[" << chan3 << "] = " << length << std::endl;
    if(length != 0){
      chan_ind1 = Amps1.D1.T3_3_index[chan3];
      chan_ind2 = Amps1.D1.T3_4_index[chan3];
      chan_ind3 = Eff_Ints.Xhh.X_3_index[chan3];

      X3 = new double[nh * nh];
      Eff_Ints.Xhh.Set_X_3_OD(X3, nh*nh, chan_ind3);
      TT3 = new double[npph * nh];
      Amps1.D1.Set_T3_3(TT3, npph*nh, chan_ind1);
      T3 = new double[npph * nh];
      for(int ind = 0; ind < npph*nh; ++ind){ T3[ind] = 0.0; }

      /*if(chan3 == 16 || chan3 == 17){
	std::cout << "Chan3 = " << chan3 << ", T3_3 = " << std::endl;
	for(int pph = 0; pph < npph; ++pph){
	  for(int h = 0; h < nh; ++h){
	    std::cout << TT3[pph*nh + h] << " ";
	  }
	  //std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "Xhh3 = " << std::endl;
	for(int h1 = 0; h1 < nh; ++h1){
	  for(int h2 = 0; h2 < nh; ++h2){
	    std::cout << X3[h1*nh + h2] << " ";
	  }
	  //std::cout << std::endl;
	}
	std::cout << std::endl;
	}*/

      //T3_3(ab|ij){abj',i}  =  -T3_3(ab|kj){abj',k}.Xhh3_od(k|i){k,i}
      dgemm_NN(TT3, X3, T3, &npph, &nh, &nh, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T3_3(T3, npph*nh, chan_ind1);
      //T3_4(ab|ij){abi',j}  =  T3_3(ab|ki){abi',k}.Xhh3_od(k|j){k,j}
      dgemm_NN(TT3, X3, T3, &npph, &nh, &nh, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T3_4(T3, npph*nh, chan_ind2);

      delete[] X3;
      delete[] TT3;
      delete[] T3;
    }
  }
}

void Singles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int nh, np, npph, nhhp, one = 1;
  int nhp0 = Chan.nhp1[Chan.ind0];
  int nph0 = Chan.nph1[Chan.ind0];
  int chan_ind1, chan_ind2, chan_ind3;
  double *t3_1, *t3_2, *T2, *T3, *Xpp3, *Xhh3;
  char N = 'N';

  if(PAR.HF == 0){
    // t2(a|i){ai}  <-  fph(a|i){ai}
    for(int ph = 0; ph < nph0; ++ph){ Amps2.S1.t2[ph] += Ints.Fmatrix.ph_2[ph]; }
  }

  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];

    chan_ind1 = Amps1.S1.t3_index[chan3];
    t3_1 = new double[np * nh];
    Amps1.S1.Set_t3(t3_1, np*nh, chan_ind1);
    t3_2 = new double[np * nh];
    for(int ind = 0; ind < np*nh; ++ind){ t3_2[ind] = 0.0; }

    if(nh * np != 0){
      chan_ind2 = Eff_Ints.Xpp.X_3_index[chan3];
      Xpp3 = new double[np * np];
      Eff_Ints.Xpp.Set_X_3_OD(Xpp3, np*np, chan_ind2);
      //t3(a|i){a,i}  <-  Xpp3_od(a|c){a,c}.t3(c|i){c,i}
      dgemm_NN(Xpp3, t3_1, t3_2, &np, &nh, &np, &p1, &p1, &N, &N);

      chan_ind2 = Eff_Ints.Xhh.X_3_index[chan3];
      Xhh3 = new double[nh * nh];
      Eff_Ints.Xhh.Set_X1_3_OD(Xhh3, nh*nh, chan_ind2);
      //t3(a|i){a,i}  <-  - t3(a|k){a,k}.X1hh3_od(k|i){k,i}
      dgemm_NN(t3_1, Xhh3, t3_2, &np, &nh, &nh, &m1, &p1, &N, &N);

      delete[] Xpp3;
      delete[] Xhh3;
    }
    if(npph * np * nh != 0){
      chan_ind2 = Ints.Vhppp.V_3_2_index[chan3];
      chan_ind3 = Amps1.D1.T3_4_index[chan3];

      T3 = new double[npph * nh];
      Amps1.D1.Set_T3_4(T3, npph*nh, chan_ind3);
      //t3(a|i){a,i}  <-  (1/2) Vhppp3_2(ka|cd){a,cdk'}.T3_4(cd|ki){cdk',i}
      dgemm_NN((Ints.Vhppp.V_3_2 + chan_ind2), T3, t3_2, &np, &nh, &npph, &p12, &p1, &N, &N);

      delete[] T3;
    }
    if(nhhp * np * nh != 0){
      chan_ind2 = Amps1.D1.T3_1_index[chan3];
      chan_ind3 = Ints.Vhhhp.V_3_3_index[chan3];

      T3 = new double[np * nhhp];
      Amps1.D1.Set_T3_1(T3, np*nhhp, chan_ind2);
      //t3(a|i){a,i}  <-  -(1/2) T3_1(ac|kl){a,klc'}.Vhhhp3_3(kl|ic){klc',i}
      dgemm_NN(T3, (Ints.Vhhhp.V_3_3 + chan_ind3), t3_2, &np, &nh, &nhhp, &m12, &p1, &N, &N);

      delete[] T3;
    }
    Amps2.S1.Gather_t3(t3_2, np*nh, chan_ind1);

    delete[] t3_1;
    delete[] t3_2;
  }
  if(nhp0 != 0){
    chan_ind1 = Ints.Vhphp.V_2_2_index[Chan.ind0];
    //t2(a|i){ai'}  =  -Vhphp2_2(ka|ic){ai',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhphp.V_2_2 + chan_ind1), Amps1.S1.t2, Amps2.S1.t2, &nph0, &one, &nph0, &m1, &p1, &N, &N);
  }
  if(nph0 * nhp0 != 0){
    chan_ind1 = Amps1.D1.T2_3_index[Chan.ind0];

    T2 = new double[nph0 * nhp0];
    Amps1.D1.Set_T2_3(T2, nph0*nhp0, chan_ind1);
    //t2(a|i){ai'}  <-  T2_3(ac|ik){ai',kc'}.Xhp2(k|c){kc'}
    dgemm_NN(T2, Eff_Ints.Xhp.X_2, Amps2.S1.t2, &nph0, &one, &nhp0, &p1, &p1, &N, &N);

    delete[] T2;
  }
}

void Doubles_Step_2(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, m1 = -1.0;
  int nh, np, npph, nhhp;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  double *Xhphh, *Xpphp, *T3, *t3;
  char N = 'N';

  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];

    chan_ind1 = Amps1.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps1.S1.Set_t3(t3, np*nh, chan_ind1);

    if(nh * np * nhhp != 0){
      chan_ind2 = Eff_Ints.Xhphh.X_3_1_index[chan3];
      chan_ind3 = Amps1.D1.T3_1_index[chan3];
      chan_ind4 = Amps1.D1.T3_2_index[chan3];

      Xhphh = new double[nh * nhhp];
      Eff_Ints.Xhphh.Set_X1_3_1(Xhphh, nh*nhhp, chan_ind2);
      T3 = new double[np * nhhp];
      for(int ind = 0; ind < np*nhhp; ++ind){ T3[ind] = 0.0; }
      //T3_1(ab|ij){a,ijb'}  =  -t3(a|k){a,k}.X1hphh3_1(kb|ij){k,ijb'}
      dgemm_NN(t3, Xhphh, T3, &np, &nhhp, &nh, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T3_1(T3, np*nhhp, chan_ind3);
      //T3_2(ab|ij){b,ija'}  =  t3(b|k){b,k}.X1hphh3_1(ka|ij){k,ija'}
      dgemm_NN(t3, Xhphh, T3, &np, &nhhp, &nh, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T3_2(T3, np*nhhp, chan_ind4);

      delete[] T3;
      delete[] Xhphh;
    }
    //std::cout << "?????  Length[" << chan3 << "] = " << nh * np * npph << std::endl;
    if(nh * np * npph != 0){
      chan_ind2 = Eff_Ints.Xpphp.X_3_4_index[chan3];
      chan_ind3 = Amps1.D1.T3_4_index[chan3];
      chan_ind4 = Amps1.D1.T3_3_index[chan3];

      Xpphp = new double[npph * np];
      Eff_Ints.Xpphp.Set_X1_3_4(Xpphp, npph*np, chan_ind2);
      T3 = new double[npph * nh];
      for(int ind = 0; ind < npph*nh; ++ind){ T3[ind] = 0.0; }

      //if( chan3 == 17 ){ std::cout << "!! " << npph*np << " " << chan_ind2 << std::endl; }
      /*if( chan3 == 17 ){
	std::cout << "Chan3 = " << chan3 << std::endl << std::endl;
	for(int pph = 0; pph < npph; ++pph){
	  std::cout << "(" << Chan.pph_state(chan3, pph).v1 << "," << Chan.pph_state(chan3, pph).v2 << "," << Chan.pph_state(chan3, pph).v3 << ")  ";
	}
	std::cout << std::endl;
	for(int p = 0; p < np; ++p){
	  std::cout << "(" << Chan.p_state(chan3, p).v1 << ")  ";
	}
	std::cout << std::endl;
	std::cout << "X1pphp3_4 (" << chan_ind2 << ") = " << std::endl;
	for(int pph = 0; pph < npph; ++pph){
	  for(int p = 0; p < np; ++p){
	    std::cout << Xpphp[pph*np + p] << " ";
	  }
	  //std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
	std::cout << "t3 = " << std::endl;
	for(int p = 0; p < np; ++p){
	  for(int h = 0; h < nh; ++h){
	    std::cout << t3[p*nh + h] << " ";
	  }
	  //std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
	}*/

      //T3_4(ab|ij){abi',j}  =  X1pphp3_4(ab|ic){abi',c}.t3(c|j){c,j}
      dgemm_NN(Xpphp, t3, T3, &npph, &nh, &np, &p1, &p1, &N, &N);
      Amps2.D1.Gather_T3_4(T3, npph*nh, chan_ind3);
      //T3_3(ab|ij){abj',i}  =  -X1pphp3_4(ab|jc){abj',c}.t3(c|i){c,i}
      dgemm_NN(Xpphp, t3, T3, &npph, &nh, &np, &m1, &p1, &N, &N);
      Amps2.D1.Gather_T3_3(T3, npph*nh, chan_ind4);

      delete[] T3;
      delete[] Xpphp;
    }
    delete[] t3;
  }
}
