#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

Interactions::Interactions(Input_Parameters &Parameters, Channels &Chan)
{
  #pragma omp parallel sections
  {
    #pragma omp section
    { Fmatrix = F_matrix(Parameters, Chan); }
    #pragma omp section
    { Vhhhh = V_hhhh(Parameters, Chan); }
    #pragma omp section
    { Vpppp = V_pppp(Parameters, Chan); }
    #pragma omp section
    { Vhhpp = V_hhpp(Parameters, Chan); }
    #pragma omp section
    { Vpphh = V_pphh(Parameters, Chan); }
    #pragma omp section
    { Vhphp = V_hphp(Parameters, Chan); }
    #pragma omp section
    { Vhhhp = V_hhhp(Parameters, Chan); }
    #pragma omp section
    { Vhppp = V_hppp(Parameters, Chan); }
    #pragma omp section
    { Vhphh = V_hphh(Parameters, Chan); }
    #pragma omp section
    { Vpphp = V_pphp(Parameters, Chan); }
  }
}

void Interactions::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  Fmatrix.delete_struct(Parameters, Chan);
  Vhhhh.delete_struct(Parameters, Chan);
  Vpppp.delete_struct(Parameters, Chan);
  Vhhpp.delete_struct(Parameters, Chan);
  Vpphh.delete_struct(Parameters, Chan);
  Vhphp.delete_struct(Parameters, Chan);
  Vhhhp.delete_struct(Parameters, Chan);
  Vhppp.delete_struct(Parameters, Chan);
  Vhphh.delete_struct(Parameters, Chan);
  Vpphp.delete_struct(Parameters, Chan);
}

F_matrix::F_matrix(Input_Parameters &Parameters, Channels &Chan)
{
  int length, ind, chan3;
  hh_2_length = Chan.nhh1[Chan.ind0];
  hp_2_length = Chan.nhp1[Chan.ind0];
  ph_2_length = Chan.nph1[Chan.ind0];
  pp_2_length = Chan.npp1[Chan.ind0];
  hh_2 = new double[hh_2_length];
  hp_2 = new double[hp_2_length];
  ph_2 = new double[ph_2_length];
  pp_2 = new double[pp_2_length];
  for(ind = 0; ind < hh_2_length; ++ind){ hh_2[ind] = 0.0; }
  for(ind = 0; ind < hp_2_length; ++ind){ hp_2[ind] = 0.0; }
  for(ind = 0; ind < ph_2_length; ++ind){ ph_2[ind] = 0.0; }
  for(ind = 0; ind < pp_2_length; ++ind){ pp_2[ind] = 0.0; }

  length = 0;
  hh_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    hh_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nh[chan3];
  }
  hh_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ hh_3[ind] = 0.0; }
  hh_3_length = length;

  length = 0;
  hp_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    hp_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.np[chan3];
  }
  hp_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ hp_3[ind] = 0.0; }
  hp_3_length = length;

  length = 0;
  ph_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    ph_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nh[chan3];
  }
  ph_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ ph_3[ind] = 0.0; }
  ph_3_length = length;

  length = 0;
  pp_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    pp_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.np[chan3];
  }
  pp_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ pp_3[ind] = 0.0; }
  pp_3_length = length;
}

void F_matrix::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] hh_2;
  delete[] hp_2;
  delete[] ph_2;
  delete[] pp_2;
  delete[] hh_3;
  delete[] hp_3;
  delete[] ph_3;
  delete[] pp_3;
  delete[] hh_3_index;
  delete[] hp_3_index;
  delete[] ph_3_index;
  delete[] pp_3_index;
}

V_hhhh::V_hhhh(Input_Parameters &Parameters, Channels &Chan)
{
  int chan1, chan3, ind, length;
  length = 0;
  V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    V_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhh[chan1];
  }
  V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_1[ind] = 0.0; }
  V_1_length = length;

  length = 0;
  V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhhh[chan3];
  }
  V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_2[ind] = 0.0; }
  V_3_2_length = length;
}

void V_hhhh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_1;
  delete[] V_3_2;
  delete[] V_1_index;
  delete[] V_3_2_index;
}

V_pppp::V_pppp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan1, chan3, ind, length;
  length = 0;
  V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    V_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.npp[chan1];
  }
  V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_1[ind] = 0.0; }
  V_1_length = length;

  length = 0;
  V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_3_index[chan3] = length;;
    length += Chan.nppp[chan3] * Chan.np[chan3];
  }
  V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_3[ind] = 0.0; }
  V_3_3_length = length;
}

void V_pppp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_1;
  delete[] V_3_3;
  delete[] V_1_index;
  delete[] V_3_3_index;
}

V_hhpp::V_hhpp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  length = 0;
  V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    V_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.npp[chan1];
  }
  V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_1[ind] = 0.0; }
  V_1_length = length;

  length = 0;
  V_2_1_index = new int[Chan.size2];
  V_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    V_2_1_index[chan2] = length;
    V_2_3_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nph1[chan2];
  }
  V_2_1 = new double[length];
  V_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    V_2_1[ind] = 0.0;
    V_2_3[ind] = 0.0;
  }
  V_2_1_length = length;
  V_2_3_length = length;

  length = 0;
  V_3_1_index = new int[Chan.size3];
  V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_1_index[chan3] = length;
    V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.npph[chan3];
  }
  V_3_1 = new double[length];
  V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    V_3_1[ind] = 0.0;
    V_3_2[ind] = 0.0;
  }
  V_3_1_length = length;
  V_3_2_length = length;

  length = 0;
  V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.np[chan3];
  }
  V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_3[ind] = 0.0; }
  V_3_3_length = length;
}

void V_hhpp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_1;
  delete[] V_2_1;
  delete[] V_2_3;
  delete[] V_3_1;
  delete[] V_3_2;
  delete[] V_3_3;
  delete[] V_1_index;
  delete[] V_2_1_index;
  delete[] V_2_3_index;
  delete[] V_3_1_index;
  delete[] V_3_2_index;
  delete[] V_3_3_index;
}

V_pphh::V_pphh(Input_Parameters &Parameters, Channels &Chan)
{
  int chan1, ind, length;
  length = 0;
  V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    V_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhh[chan1];
  }
  V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_1[ind] = 0.0; }
  V_1_length = length;
}

void V_pphh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_1;
  delete[] V_1_index;
}

V_hphp::V_hphp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan2, ind, length;
  length = 0;
  V_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    V_2_1_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nhp1[chan2];
  }
  V_2_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_2_1[ind] = 0.0; }
  V_2_1_length = length;

  length = 0;
  V_2_2_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    V_2_2_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.nph1[chan2];
  }
  V_2_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_2_2[ind] = 0.0; }
  V_2_2_length = length;
}

void V_hphp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_2_1;
  delete[] V_2_2;
  delete[] V_2_1_index;
  delete[] V_2_2_index;
}

V_hhhp::V_hhhp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan2, chan3, ind, length;
  length = 0;
  V_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    V_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nph1[chan2];
  }
  V_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_2_3[ind] = 0.0; }
  V_2_3_length = length;

  length = 0;
  V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhph[chan3];
  }
  V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_2[ind] = 0.0; }
  V_3_2_length = length;

  length = 0;
  V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.nh[chan3];
  }
  V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_3[ind] = 0.0; }
  V_3_3_length = length;
}

void V_hhhp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_2_3;
  delete[] V_3_2;
  delete[] V_3_3;
  delete[] V_2_3_index;
  delete[] V_3_2_index;
  delete[] V_3_3_index;
}

V_hppp::V_hppp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan2, chan3, ind, length;
  length = 0;
  V_2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    V_2_4_index[chan2] = length;
    length += Chan.npp1[chan2] * Chan.nph1[chan2];
  }
  V_2_4 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_2_4[ind] = 0.0; }
  V_2_4_length = length;

  length = 0;
  V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.npph[chan3];
  }
  V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_2[ind] = 0.0; }
  V_3_2_length = length;
}

void V_hppp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_2_4;
  delete[] V_3_2;
  delete[] V_2_4_index;
  delete[] V_3_2_index;
}

V_hphh::V_hphh(Input_Parameters &Parameters, Channels &Chan)
{
  int chan3, ind, length;
  length = 0;
  V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhh[chan3];
  }
  V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_2[ind] = 0.0; }
  V_3_2_length = length;
}

void V_hphh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_3_2;
  delete[] V_3_2_index;
}

V_pphp::V_pphp(Input_Parameters &Parameters, Channels &Chan)
{
  int chan3, ind, length;
  length = 0;
  V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    V_3_3_index[chan3] = length;
    length += Chan.nppp[chan3] * Chan.nh[chan3];
  }
  V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ V_3_3[ind] = 0.0; }
  V_3_3_length = length;
}

void V_pphp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] V_3_3;
  delete[] V_3_3_index;
}

Eff_Interactions::Eff_Interactions(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  #pragma omp parallel sections
  {
    #pragma omp section
    { Xhh = X_hh(Parameters, Space, Chan); }
    #pragma omp section
    { Xpp = X_pp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhp = X_hp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhhhh = X_hhhh(Parameters, Space, Chan); }
    #pragma omp section
    { Xpppp = X_pppp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhphp = X_hphp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhhhp = X_hhhp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhppp = X_hppp(Parameters, Space, Chan); }
    #pragma omp section
    { Xhphh = X_hphh(Parameters, Space, Chan); }
    #pragma omp section
    { Xpphp = X_pphp(Parameters, Space, Chan); }
  }
}

void Eff_Interactions::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  Xhh.delete_struct(Parameters, Chan);
  Xpp.delete_struct(Parameters, Chan);
  Xhp.delete_struct(Parameters, Chan);
  Xhhhh.delete_struct(Parameters, Chan);
  Xpppp.delete_struct(Parameters, Chan);
  Xhphp.delete_struct(Parameters, Chan);
  Xhhhp.delete_struct(Parameters, Chan);
  Xhppp.delete_struct(Parameters, Chan);
  Xhphh.delete_struct(Parameters, Chan);
  Xpphp.delete_struct(Parameters, Chan);
}

void Eff_Interactions::zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  Xhh.X1_zero(Parameters, Chan, flag);
  Xhh.X_zero(Parameters, Chan, flag);
  Xpp.X_zero(Parameters, Chan, flag);
  Xhp.X_zero(Parameters, Chan, flag);
  Xhhhh.X_zero(Parameters, Chan, flag);
  Xpppp.X1_zero(Parameters, Chan, flag);
  Xpppp.X_zero(Parameters, Chan, flag);
  Xhphp.X_zero(Parameters, Chan, flag);
  Xhphp.X1_zero(Parameters, Chan, flag);
  Xhphp.X2_zero(Parameters, Chan, flag);
  Xhphp.X3_zero(Parameters, Chan, flag);
  Xhhhp.X1_zero(Parameters, Chan, flag);
  Xhhhp.X_zero(Parameters, Chan, flag);
  Xhppp.X1_zero(Parameters, Chan, flag);
  Xhppp.X_zero(Parameters, Chan, flag);
  Xhphh.X1_zero(Parameters, Chan, flag);
  Xhphh.X_zero(Parameters, Chan, flag);
  Xpphp.X1_zero(Parameters, Chan, flag);
  Xpphp.X_zero(Parameters, Chan, flag);
}

void Eff_Interactions::diagonalize(Input_Parameters &Parameters, Channels &Chan)
{
  Xhh.diagonalize(Parameters, Chan);
  Xpp.diagonalize(Parameters, Chan);
}

X_hh::X_hh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan3, ind, length;
  length = Chan.nhh1[Chan.ind0];
  X1_2 = new double[length];
  X_2 = new double[length];
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_2[ind] = 0.0;
    X_2[ind] = 0.0;
  }
  X_2_length = length;

  length = 0;
  X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nh[chan3];
  }
  X1_3 = new double[length];
  X1_3od = new double[length];
  X_3 = new double[length];
  X_3od = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3[ind] = 0.0;
    X1_3od[ind] = 0.0;
    X_3[ind] = 0.0;
    X_3od[ind] = 0.0;
  }
  X_3_length = length;

  length = 0;
  X_3d_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3d_index[chan3] = length;
    length += Chan.nh[chan3];
  }
  X1_3d = new double[length];
  X_3d = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3d[ind] = 0.0;
    X_3d[ind] = 0.0;
  }
  X_3d_length = length;

  int i, j;
  two_body hh1;
  for(ind = 0; ind < X_2_length; ++ind){
    hh1 = Chan.hh1_state(Chan.ind0, ind);
    i = hh1.v1;
    j = hh1.v2;
    Map_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, ind, Chan.h_map,Chan.h_map,Chan.nh, i,j);
  }
}

void X_hh::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_2_length; ++ind){
      X1_2[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_length; ++ind){
    X1_3[ind] = 0.0;
    X1_3od[ind] = 0.0;
  }
  for(ind = 0; ind < X_3d_length; ++ind){
    X1_3d[ind] = 0.0;
  }
}

void X_hh::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_2_length; ++ind){
      X_2[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_length; ++ind){
    X_3[ind] = 0.0;
    X_3od[ind] = 0.0;
  }
  for(ind = 0; ind < X_3d_length; ++ind){
    X_3d[ind] = 0.0;
  }
}

void X_hh::diagonalize(Input_Parameters &Parameters, Channels &Chan)
{
  int nh, chan3, h, ind;
  for(ind = 0; ind < X_3_length; ++ind){
    X1_3od[ind] = X1_3[ind];
    X_3od[ind] = X_3[ind];
  }
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    for(h = 0; h < nh; ++h){
      ind = h * nh + h;
      X1_3od[X_3_index[chan3] + ind] = 0.0;
      X_3od[X_3_index[chan3] + ind] = 0.0;
      X1_3d[X_3d_index[chan3] + h] = X1_3[X_3_index[chan3] + ind];
      X_3d[X_3d_index[chan3] + h] = X_3[X_3_index[chan3] + ind];
    }
  }
}

void X_hh::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  double x;
  for(int ind = 0; ind < X_2_length; ++ind){
    x = 0.0;
    x += X1_2[ind];
    x += map_fac2[ind] * X1_3[X_3_index[map_chan[ind]] + map_ind[ind]];
    X1_2[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_2_length; ++ind){
    x = X1_2[ind];
    X1_3[X_3_index[map_chan[ind]] + map_ind[ind]] += map_fac1[ind] * x;
  }
}

void X_hh::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  double x;
  for(int ind = 0; ind < X_2_length; ++ind){
    x = 0.0;
    x += X_2[ind];
    x += map_fac2[ind] * X_3[X_3_index[map_chan[ind]] + map_ind[ind]];
    X_2[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_2_length; ++ind){
    x = X_2[ind];
    X_3[X_3_index[map_chan[ind]] + map_ind[ind]] += map_fac1[ind] * x;
  }
}

void X_hh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_2;
  delete[] X1_2;
  delete[] X_3;
  delete[] X1_3;
  delete[] X_3d;
  delete[] X1_3d;
  delete[] X_3od;
  delete[] X1_3od;

  delete[] X_3_index;
  delete[] X_3d_index;

  delete[] map_chan;
  delete[] map_ind;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_pp::X_pp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan3, ind, length;
  length = Chan.npp1[Chan.ind0];
  X_2 = new double[length];
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2[ind] = 0.0;
  }
  X_2_length = length;

  length = 0;
  X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.np[chan3];
  }
  X_3 = new double[length];
  X_3od = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3[ind] = 0.0;
    X_3od[ind] = 0.0;
  }
  X_3_length = length;

  length = 0;
  X_3d_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3d_index[chan3] = length;
    length += Chan.np[chan3];
  }
  X_3d = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3d[ind] = 0.0;
  }
  X_3d_length = length;

  int a, b;
  two_body pp1;
  for(ind = 0; ind < X_2_length; ++ind){
    pp1 = Chan.pp1_state(Chan.ind0, ind);
    a = pp1.v1;
    b = pp1.v2;
    Map_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, ind, Chan.p_map,Chan.p_map,Chan.np, a,b);
  }
}

void X_pp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_2_length; ++ind){
      X_2[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_length; ++ind){
    X_3[ind] = 0.0;
    X_3od[ind] = 0.0;
  }
  for(ind = 0; ind < X_3d_length; ++ind){
    X_3d[ind] = 0.0;
  }
}

void X_pp::diagonalize(Input_Parameters &Parameters, Channels &Chan)
{
  int np, chan3, p, ind;
  for(ind = 0; ind < X_3_length; ++ind){
    X_3od[ind] = X_3[ind];
  }
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    for(p = 0; p < np; ++p){
      ind = p * np + p;
      X_3od[X_3_index[chan3] + ind] = 0.0;
      X_3d[X_3d_index[chan3] + p] = X_3[X_3_index[chan3] + ind];
    }
  }
}

void X_pp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  double x;
  for(int ind = 0; ind < X_2_length; ++ind){
    x = 0.0;
    x += X_2[ind];
    x += map_fac2[ind] * X_3[X_3_index[map_chan[ind]] + map_ind[ind]];
    X_2[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_2_length; ++ind){
    x = X_2[ind];
    X_3[X_3_index[map_chan[ind]] + map_ind[ind]] += map_fac1[ind] * x;
  }
}

void X_pp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_2;
  delete[] X_3;
  delete[] X_3d;
  delete[] X_3od;

  delete[] X_3_index;
  delete[] X_3d_index;

  delete[] map_chan;
  delete[] map_ind;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hp::X_hp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan3, ind, length;
  length = Chan.nhp1[Chan.ind0];
  X_2 = new double[length];
  map_chan = new int[length];
  map_ind = new int[length];
  map_fac1 = new double[length];
  map_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2[ind] = 0.0;
  }
  X_2_length = length;

  length = 0;
  X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.np[chan3];
  }
  X_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3[ind] = 0.0;
  }
  X_3_length = length;

  int i, a;
  two_body hp1;
  for(ind = 0; ind < X_2_length; ++ind){
    hp1 = Chan.hp1_state(Chan.ind0, ind);
    i = hp1.v1;
    a = hp1.v2;
    Map_2(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, ind, Chan.h_map,Chan.p_map,Chan.np, i,a);
  }
}

void X_hp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_2_length; ++ind){
      X_2[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_length; ++ind){
    X_3[ind] = 0.0;
  }
}

void X_hp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  double x;
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_2_length; ++ind){
    x = X_2[ind];
    X_3[X_3_index[map_chan[ind]] + map_ind[ind]] += map_fac1[ind] * x;
  }
}

void X_hp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_2;
  delete[] X_3;

  delete[] X_3_index;

  delete[] map_chan;
  delete[] map_ind;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hhhh::X_hhhh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan3, ind, length, count;
  int i, j, k, l, nhh;
  two_body hh1, hh2;
  double J;

  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhh[chan1];
  }
  X_1 = new double[length];
  map_num = new int[2 * length]; // X_3_3, X_3_4
  map_index = new int[2 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    X_3_4_index[chan3] = length;
    length += Chan.nhhh[chan3] * Chan.nh[chan3];
  }
  X_3_3 = new double[length];
  X_3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_3[ind] = 0.0;
    X_3_4[ind] = 0.0;
  }
  X_3_3_length = length;
  X_3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
      for(int hh2_ind = 0; hh2_ind < nhh; ++hh2_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	hh2 = Chan.hh_state(chan1, hh2_ind);
	i = hh1.v1;
	j = hh1.v2;
	k = hh2.v1;
	l = hh2.v2;
	Map_4_count(Parameters,Space, map_index,map_num, 2,7,0,length,count, i,j,k,l);
	Map_4_count(Parameters,Space, map_index,map_num, 2,8,1,length,count, i,j,k,l);
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
    nhh = Chan.nhh[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
      for(int hh2_ind = 0; hh2_ind < nhh; ++hh2_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	hh2 = Chan.hh_state(chan1, hh2_ind);
	i = hh1.v1;
	j = hh1.v2;
	k = hh2.v1;
	l = hh2.v2;
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 2,7,0,count, Chan.hhh_map,Chan.h_map, Chan.nh, i,j,k,l,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 2,8,1,count, Chan.hhh_map,Chan.h_map, Chan.nh, i,j,k,l,J);
	++count;
      }
    }
  }
}

void X_hhhh::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){
      X_1[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_3_length; ++ind){
    X_3_3[ind] = 0.0;
  }
  for(ind = 0; ind < X_3_4_length; ++ind){
    X_3_4[ind] = 0.0;
  }
}

void X_hhhh::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    x += X_1[ind];
    for(int n = 0; n < map_num[2*ind]; ++n){
      index = map_index[2*ind] + n;
      x += map_fac2[index] * X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[2*ind + 1]; ++n){
      index = map_index[2*ind + 1] + n;
      x += map_fac2[index] * X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]];
    }
    X_1[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1[ind];
    for(int n = 0; n < map_num[2*ind]; ++n){
      index = map_index[2*ind] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[2*ind + 1]; ++n){
      index = map_index[2*ind + 1] + n;
      X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_hhhh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X_3_3;
  delete[] X_3_4;

  delete[] X_1_index;
  delete[] X_3_3_index;
  delete[] X_3_4_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_pppp::X_pppp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan3, ind, length, count;
  int a, b, c, d, npp;
  two_body pp1, pp2;
  double J;
  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.npp[chan1];
  }
  X_1 = new double[length];
  X1_1 = new double[length];
  map_num = new int[2 * length]; // X1_3_1, X1_3_2
  map_index = new int[2 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
    X1_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_3_1_index = new int[Chan.size3];
  X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_1_index[chan3] = length;
    X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nppp[chan3];
  }
  X1_3_1 = new double[length];
  X1_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_1[ind] = 0.0;
    X1_3_2[ind] = 0.0;
  }
  X_3_1_length = length;
  X_3_2_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    for(int pp1_ind = 0; pp1_ind < npp; ++pp1_ind){
      for(int pp2_ind = 0; pp2_ind < npp; ++pp2_ind){
	pp1 = Chan.pp_state(chan1, pp1_ind);
	pp2 = Chan.pp_state(chan1, pp2_ind);
	a = pp1.v1;
	b = pp1.v2;
	c = pp2.v1;
	d = pp2.v2;
	Map_4_count(Parameters,Space, map_index,map_num, 2,5,0,length,count, a,b,c,d);
	Map_4_count(Parameters,Space, map_index,map_num, 2,6,1,length,count, a,b,c,d);
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
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int pp1_ind = 0; pp1_ind < npp; ++pp1_ind){
      for(int pp2_ind = 0; pp2_ind < npp; ++pp2_ind){
	pp1 = Chan.pp_state(chan1, pp1_ind);
	pp2 = Chan.pp_state(chan1, pp2_ind);
	a = pp1.v1;
	b = pp1.v2;
	c = pp2.v1;
	d = pp2.v2;
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 2,5,0,count, Chan.p_map,Chan.ppp_map, Chan.nppp, a,b,c,d,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 2,6,1,count, Chan.p_map,Chan.ppp_map, Chan.nppp, a,b,c,d,J);
	++count;
      }
    }
  }
}

void X_pppp::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){
      X1_1[ind] = 0.0;
    }
  }
  for(ind = 0; ind < X_3_1_length; ++ind){
    X1_3_1[ind] = 0.0;
  }
  for(ind = 0; ind < X_3_2_length; ++ind){
    X1_3_2[ind] = 0.0;
  }
}

void X_pppp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){
      X1_1[ind] = 0.0;
    }
  }
}

void X_pppp::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    x += X1_1[ind];
    for(int n = 0; n < map_num[2*ind]; ++n){
      index = map_index[2*ind] + n;
      x += map_fac2[index] * X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[2*ind + 1]; ++n){
      index = map_index[2*ind + 1] + n;
      x += map_fac2[index] * X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    X1_1[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X1_1[ind];
    for(int n = 0; n < map_num[2*ind]; ++n){
      index = map_index[2*ind] + n;
      X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[2*ind + 1]; ++n){
      index = map_index[2*ind + 1] + n;
      X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_pppp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X1_1;
  delete[] X1_3_1;
  delete[] X1_3_2;

  delete[] X_1_index;
  delete[] X_3_1_index;
  delete[] X_3_2_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hphp::X_hphp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count;
  int i, a, j, b, nhp;
  two_body hp1, hp2;
  double J;
  length = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    length += Chan.nhp[chan1] * Chan.nhp[chan1];
  }
  map_num = new int[5 * length]; // X_2_1, X_3_1, X_3_2, X_3_3, X_3_4
  map_index = new int[5 * length];
  X_1_length = length;

  length = 0;
  X_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_1_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nhp1[chan2];
  }
  X_2_1 = new double[length];
  X1_2_1 = new double[length];
  X2_2_1 = new double[length];
  X3_2_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_1[ind] = 0.0;
    X1_2_1[ind] = 0.0;
    X2_2_1[ind] = 0.0;
    X3_2_1[ind] = 0.0;
  }
  X_2_1_length = length;

  length = 0;
  X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhpp[chan3];
  }
  X1_3_1 = new double[length];
  X2_3_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_1[ind] = 0.0;
    X2_3_1[ind] = 0.0;
  }
  X_3_1_length = length;

  length = 0;
  X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhph[chan3];
  }
  X_3_2 = new double[length];
  X1_3_2 = new double[length];
  X2_3_2 = new double[length];
  X3_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_2[ind] = 0.0;
    X1_3_2[ind] = 0.0;
    X2_3_2[ind] = 0.0;
    X3_3_2[ind] = 0.0;
  }
  X_3_2_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    length += Chan.nhpp[chan3] * Chan.nh[chan3];
  }
  X_3_3 = new double[length];
  X1_3_3 = new double[length];
  X2_3_3 = new double[length];
  X3_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_3[ind] = 0.0;
    X1_3_3[ind] = 0.0;
    X2_3_3[ind] = 0.0;
    X3_3_3[ind] = 0.0;
  }
  X_3_3_length = length;

  length = 0;
  X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_4_index[chan3] = length;
    length += Chan.nhph[chan3] * Chan.np[chan3];
  }
  X3_3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X3_3_4[ind] = 0.0;
  }
  X_3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhp = Chan.nhp[chan1];
    for(int hp1_ind = 0; hp1_ind < nhp; ++hp1_ind){
      for(int hp2_ind = 0; hp2_ind < nhp; ++hp2_ind){
	hp1 = Chan.hp_state(chan1, hp1_ind);
	hp2 = Chan.hp_state(chan1, hp2_ind);
	i = hp1.v1;
	a = hp1.v2;
	j = hp2.v1;
	b = hp2.v2;
	// X_2_1, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4_count(Parameters,Space, map_index,map_num, 5,1,0,length,count, i,a,j,b);
	Map_4_count(Parameters,Space, map_index,map_num, 5,5,1,length,count, i,a,j,b);
	Map_4_count(Parameters,Space, map_index,map_num, 5,6,2,length,count, i,a,j,b);
	Map_4_count(Parameters,Space, map_index,map_num, 5,7,3,length,count, i,a,j,b);
	Map_4_count(Parameters,Space, map_index,map_num, 5,8,4,length,count, i,a,j,b);
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
    nhp = Chan.nhp[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int hp1_ind = 0; hp1_ind < nhp; ++hp1_ind){
      for(int hp2_ind = 0; hp2_ind < nhp; ++hp2_ind){
	hp1 = Chan.hp_state(chan1, hp1_ind);
	hp2 = Chan.hp_state(chan1, hp2_ind);
	i = hp1.v1;
	a = hp1.v2;
	j = hp2.v1;
	b = hp2.v2;
	// X_2_1, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 5,1,0,count, Chan.hp1_map,Chan.hp1_map, Chan.nhp1, i,a,j,b,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 5,5,1,count, Chan.h_map,Chan.hpp_map, Chan.nhpp, i,a,j,b,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 5,6,2,count, Chan.p_map,Chan.hph_map, Chan.nhph, i,a,j,b,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 5,7,3,count, Chan.hpp_map,Chan.h_map, Chan.nh, i,a,j,b,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 5,8,4,count, Chan.hph_map,Chan.p_map, Chan.np, i,a,j,b,J);
	++count;
      }
    }
  }
}

void X_hphp::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_2_1_length; ++ind){ X1_2_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_1_length; ++ind){ X1_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X1_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X1_3_3[ind] = 0.0; }
}

void X_hphp::X2_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_2_1_length; ++ind){ X2_2_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_1_length; ++ind){ X2_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X2_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X2_3_3[ind] = 0.0; }
}

void X_hphp::X3_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_2_1_length; ++ind){ X3_2_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X3_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X3_3_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_4_length; ++ind){ X3_3_4[ind] = 0.0; }
}

void X_hphp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_2_1_length; ++ind){ X_2_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X_3_3[ind] = 0.0; }
}

void X_hphp::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      x += map_fac2[index] * X1_2_1[X_2_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      x += map_fac2[index] * X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      x += map_fac2[index] * X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X_1temp[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1temp[ind];
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      X1_2_1[X_2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 1]; ++n){
      index = map_index[5*ind + 1] + n;
      X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X_1temp;
}

void X_hphp::X2_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      x += map_fac2[index] * X2_2_1[X_2_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      x += map_fac2[index] * X2_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      x += map_fac2[index] * X2_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X_1temp[ind] = x;
  }
  X2_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1temp[ind];
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      X2_2_1[X_2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 1]; ++n){
      index = map_index[5*ind + 1] + n;
      X2_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      X2_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      X2_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X_1temp;
}

void X_hphp::X3_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      x += map_fac2[index] * X3_2_1[X_2_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      x += map_fac2[index] * X3_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      x += map_fac2[index] * X3_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X_1temp[ind] = x;
  }
  X3_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1temp[ind];
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      X3_2_1[X_2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      X3_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      X3_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 4]; ++n){
      index = map_index[5*ind + 4] + n;
      X3_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X_1temp;
}

void X_hphp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      x += map_fac2[index] * X_2_1[X_2_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      x += map_fac2[index] * X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      x += map_fac2[index] * X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X_1temp[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1temp[ind];
    for(int n = 0; n < map_num[5*ind]; ++n){
      index = map_index[5*ind] + n;
      X_2_1[X_2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 2]; ++n){
      index = map_index[5*ind + 2] + n;
      X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[5*ind + 3]; ++n){
      index = map_index[5*ind + 3] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X_1temp;
}

void X_hphp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X1_2_1;
  delete[] X2_2_1;
  delete[] X3_2_1;
  delete[] X_2_1;
  delete[] X1_3_1;
  delete[] X2_3_1;
  delete[] X1_3_2;
  delete[] X2_3_2;
  delete[] X3_3_2;
  delete[] X_3_2;
  delete[] X1_3_3;
  delete[] X2_3_3;
  delete[] X3_3_3;
  delete[] X_3_3;
  delete[] X3_3_4;

  delete[] X_2_1_index;
  delete[] X_3_1_index;
  delete[] X_3_2_index;
  delete[] X_3_3_index;
  delete[] X_3_4_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hppp::X_hppp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count;
  int i, a, b, c, nhp, npp;
  two_body hp, pp;
  double J;
  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.nhp[chan1] * Chan.npp[chan1];
  }
  X_1 = new double[length];
  map_num = new int[4 * length]; // X_2_3, X_3_1, X_3_2, X_3_3
  map_index = new int[4 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_3_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.npp1[chan2];
  }
  X_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_3[ind] = 0.0;
  }
  X_2_3_length = length;

  length = 0;
  X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nppp[chan3];
  }
  X1_3_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_1[ind] = 0.0;
  }
  X_3_1_length = length;

  length = 0;
  X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.npph[chan3];
  }
  X_3_2 = new double[length];
  X1_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_2[ind] = 0.0;
    X1_3_2[ind] = 0.0;
  }
  X_3_2_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    length += Chan.nhpp[chan3] * Chan.np[chan3];
  }
  X_3_3 = new double[length];
  X1_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_3[ind] = 0.0;
    X1_3_3[ind] = 0.0;
  }
  X_3_3_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhp = Chan.nhp[chan1];
    npp = Chan.npp[chan1];
    for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
      for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
	hp = Chan.hp_state(chan1, hp_ind);
	pp = Chan.pp_state(chan1, pp_ind);
	i = hp.v1;
	a = hp.v2;
	b = pp.v1;
	c = pp.v2;
	// X_2_3, X_3_1, X_3_2, X_3_3
	Map_4_count(Parameters,Space, map_index,map_num, 4,3,0,length,count, i,a,b,c);
	Map_4_count(Parameters,Space, map_index,map_num, 4,5,1,length,count, i,a,b,c);
	Map_4_count(Parameters,Space, map_index,map_num, 4,6,2,length,count, i,a,b,c);
	Map_4_count(Parameters,Space, map_index,map_num, 4,7,3,length,count, i,a,b,c);
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
    nhp = Chan.nhp[chan1];
    npp = Chan.npp[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
      for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
	hp = Chan.hp_state(chan1, hp_ind);
	pp = Chan.pp_state(chan1, pp_ind);
	i = hp.v1;
	a = hp.v2;
	b = pp.v1;
	c = pp.v2;
	// X_2_3, X_3_1, X_3_2, X_3_3
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 4,3,0,count, Chan.hp1_map,Chan.pp1_map, Chan.npp1, i,a,b,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 4,5,1,count, Chan.h_map,Chan.ppp_map, Chan.nppp, i,a,b,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 4,6,2,count, Chan.p_map,Chan.pph_map, Chan.npph, i,a,b,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 4,7,3,count, Chan.hpp_map,Chan.p_map, Chan.np, i,a,b,c,J);
	++count;
      }
    }
  }
}

void X_hppp::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_3_1_length; ++ind){ X1_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X1_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X1_3_3[ind] = 0.0; }
}

void X_hppp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){ X_1[ind] = 0.0; }
  }
  for(ind = 0; ind < X_2_3_length; ++ind){ X_2_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X_3_3[ind] = 0.0; }
}

void X_hppp::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X1_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[4*ind + 2]; ++n){
      index = map_index[4*ind + 2] + n;
      x += map_fac2[index] * X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    X1_1temp[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X1_1temp[ind];
    for(int n = 0; n < map_num[4*ind + 1]; ++n){
      index = map_index[4*ind + 1] + n;
      X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[4*ind + 2]; ++n){
      index = map_index[4*ind + 2] + n;
      X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[4*ind + 3]; ++n){
      index = map_index[4*ind + 3] + n;
      X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X1_1temp;
}

void X_hppp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[4*ind + 2]; ++n){
      index = map_index[4*ind + 2] + n;
      x += map_fac2[index] * X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    X_1[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1[ind];
    for(int n = 0; n < map_num[4*ind]; ++n){
      index = map_index[4*ind] + n;
      X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[4*ind + 2]; ++n){
      index = map_index[4*ind + 2] + n;
      X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[4*ind + 3]; ++n){
      index = map_index[4*ind + 3] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_hppp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X_2_3;
  delete[] X1_3_1;
  delete[] X1_3_2;
  delete[] X_3_2;
  delete[] X1_3_3;
  delete[] X_3_3;

  delete[] X_1_index;
  delete[] X_2_3_index;
  delete[] X_3_1_index;
  delete[] X_3_2_index;
  delete[] X_3_3_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hhhp::X_hhhp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count;
  int i, j, k, a, nhh, nhp;
  two_body hh, hp;
  double J;
  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhp[chan1];
  }
  X_1 = new double[length];
  map_num = new int[3 * length]; // X_2_3, X_3_3, X_3_4
  map_index = new int[3 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nph1[chan2];
  }
  X_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_3[ind] = 0.0;
  }
  X_2_3_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.nh[chan3];
  }
  X1_3_3 = new double[length];
  X_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_3[ind] = 0.0;
    X_3_3[ind] = 0.0;
  }
  X_3_3_length = length;

  length = 0;
  X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_4_index[chan3] = length;
    length += Chan.nhhh[chan3] * Chan.np[chan3];
  }
  X1_3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_4[ind] = 0.0;
  }
  X_3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    nhp = Chan.nhp[chan1];
    for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
      for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	hp = Chan.hp_state(chan1, hp_ind);
	i = hh.v1;
	j = hh.v2;
	k = hp.v1;
	a = hp.v2;
	// X_2_3, X_3_3, X_3_4
	Map_4_count(Parameters,Space, map_index,map_num, 3,3,0,length,count, i,j,k,a);
	Map_4_count(Parameters,Space, map_index,map_num, 3,7,1,length,count, i,j,k,a);
	Map_4_count(Parameters,Space, map_index,map_num, 3,8,2,length,count, i,j,k,a);
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
    nhh = Chan.nhh[chan1];
    nhp = Chan.nhp[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
      for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
	hh = Chan.hh_state(chan1, hh_ind);
	hp = Chan.hp_state(chan1, hp_ind);
	i = hh.v1;
	j = hh.v2;
	k = hp.v1;
	a = hp.v2;
	// X_2_3, X_3_3, X_3_4
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 3,3,0,count, Chan.hh1_map,Chan.ph1_map, Chan.nph1, i,j,k,a,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 3,7,1,count, Chan.hhp_map,Chan.h_map, Chan.nh, i,j,k,a,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 3,8,2,count, Chan.hhh_map,Chan.p_map, Chan.np, i,j,k,a,J);
	++count;
      }
    }
  }
}

void X_hhhp::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_3_3_length; ++ind){ X1_3_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_4_length; ++ind){ X1_3_4[ind] = 0.0; }
}

void X_hhhp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){ X_1[ind] = 0.0; }
  }
  for(ind = 0; ind < X_2_3_length; ++ind){ X_2_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X_3_3[ind] = 0.0; }
}

void X_hhhp::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X1_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[3*ind + 1]; ++n){
      index = map_index[3*ind + 1] + n;
      x += map_fac2[index] * X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X1_1temp[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X1_1temp[ind];
    for(int n = 0; n < map_num[3*ind + 1]; ++n){
      index = map_index[3*ind + 1] + n;
      X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[3*ind + 2]; ++n){
      index = map_index[3*ind + 2] + n;
      X1_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X1_1temp;
}

void X_hhhp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[3*ind + 1]; ++n){
      index = map_index[3*ind + 1] + n;
      x += map_fac2[index] * X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X_1[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1[ind];
    for(int n = 0; n < map_num[3*ind]; ++n){
      index = map_index[3*ind] + n;
      X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[3*ind + 1]; ++n){
      index = map_index[3*ind + 1] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_hhhp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X_2_3;
  delete[] X1_3_3;
  delete[] X_3_3;
  delete[] X1_3_4;

  delete[] X_1_index;
  delete[] X_2_3_index;
  delete[] X_3_3_index;
  delete[] X_3_4_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_hphh::X_hphh(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count;
  int i, a, j, k, nhp, nhh;
  two_body hp, hh;
  double J;
  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.nhp[chan1] * Chan.nhh[chan1];
  }
  X_1 = new double[length];
  map_num = new int[6 * length]; // X_2_1, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
  map_index = new int[6 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_1_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nhp1[chan2];
  }
  X_2_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_1[ind] = 0.0;
  }
  X_2_1_length = length;

  length = 0;
  X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nhp1[chan2];
  }
  X_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_3[ind] = 0.0;
  }
  X_2_3_length = length;

  length = 0;
  X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhhp[chan3];
  }
  X1_3_1 = new double[length];
  X_3_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_1[ind] = 0.0;
    X_3_1[ind] = 0.0;
  }
  X_3_1_length = length;

  length = 0;
  X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhh[chan3];
  }
  X1_3_2 = new double[length];
  X_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_2[ind] = 0.0;
    X_3_2[ind] = 0.0;
  }
  X_3_2_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    length += Chan.nhph[chan3] * Chan.nh[chan3];
  }
  X_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_3[ind] = 0.0;
  }
  X_3_3_length = length;

  length = 0;
  X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_4_index[chan3] = length;
    length += Chan.nhph[chan3] * Chan.nh[chan3];
  }
  X_3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_3_4[ind] = 0.0;
  }
  X_3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhp = Chan.nhp[chan1];
    nhh = Chan.nhh[chan1];
    for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hp = Chan.hp_state(chan1, hp_ind);
	hh = Chan.hh_state(chan1, hh_ind);
	i = hp.v1;
	a = hp.v2;
	j = hh.v1;
	k = hh.v2;
	// X_2_1, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4_count(Parameters,Space, map_index,map_num, 6,1,0,length,count, i,a,j,k);
	Map_4_count(Parameters,Space, map_index,map_num, 6,3,1,length,count, i,a,j,k);
	Map_4_count(Parameters,Space, map_index,map_num, 6,5,2,length,count, i,a,j,k);
	Map_4_count(Parameters,Space, map_index,map_num, 6,6,3,length,count, i,a,j,k);
	Map_4_count(Parameters,Space, map_index,map_num, 6,7,4,length,count, i,a,j,k);
	Map_4_count(Parameters,Space, map_index,map_num, 6,8,5,length,count, i,a,j,k);
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
    nhp = Chan.nhp[chan1];
    nhh = Chan.nhh[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
      for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
	hp = Chan.hp_state(chan1, hp_ind);
	hh = Chan.hh_state(chan1, hh_ind);
	i = hp.v1;
	a = hp.v2;
	j = hh.v1;
	k = hh.v2;
	// X_2_1, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,1,0,count, Chan.hh1_map,Chan.hp1_map, Chan.nhp1, i,a,j,k,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,3,1,count, Chan.hh1_map,Chan.hp1_map, Chan.nhp1, i,a,j,k,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,5,2,count, Chan.h_map,Chan.hhp_map, Chan.nhhp, i,a,j,k,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,6,3,count, Chan.p_map,Chan.hhh_map, Chan.nhhh, i,a,j,k,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,7,4,count, Chan.hph_map,Chan.h_map, Chan.nh, i,a,j,k,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,8,5,count, Chan.hph_map,Chan.h_map, Chan.nh, i,a,j,k,J);
	++count;
      }
    }
  }
}

void X_hphh::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_3_1_length; ++ind){ X1_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X1_3_2[ind] = 0.0; }
}

void X_hphh::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){ X_1[ind] = 0.0; }
  }
  for(ind = 0; ind < X_2_1_length; ++ind){ X_2_1[ind] = 0.0; }
  for(ind = 0; ind < X_2_3_length; ++ind){ X_2_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_1_length; ++ind){ X_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X_3_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_4_length; ++ind){ X_3_4[ind] = 0.0; }
}

void X_hphh::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X1_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      x += map_fac2[index] * X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    X1_1temp[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X1_1temp[ind];
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X1_1temp;
}

void X_hphh::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    x += X_1[ind];
    for(int n = 0; n < map_num[6*ind]; ++n){
      index = map_index[6*ind] + n;
      x += map_fac2[index] * X_2_1[X_2_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 1]; ++n){
      index = map_index[6*ind + 1] + n;
      x += map_fac2[index] * X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      x += map_fac2[index] * X_3_1[X_3_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      x += map_fac2[index] * X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      x += map_fac2[index] * X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 5]; ++n){
      index = map_index[6*ind + 5] + n;
      x += map_fac2[index] * X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]];
    }
    X_1[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1[ind];
    for(int n = 0; n < map_num[6*ind]; ++n){
      index = map_index[6*ind] + n;
      X_2_1[X_2_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 1]; ++n){
      index = map_index[6*ind + 1] + n;
      X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      X_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 5]; ++n){
      index = map_index[6*ind + 5] + n;
      X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_hphh::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X_2_1;
  delete[] X_2_3;
  delete[] X1_3_1;
  delete[] X_3_1;
  delete[] X1_3_2;
  delete[] X_3_2;
  delete[] X_3_3;
  delete[] X_3_4;

  delete[] X_1_index;
  delete[] X_2_1_index;
  delete[] X_2_3_index;
  delete[] X_3_1_index;
  delete[] X_3_2_index;
  delete[] X_3_3_index;
  delete[] X_3_4_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

X_pphp::X_pphp(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan)
{
  int chan1, chan2, chan3, ind, length, count;
  int a, b, i, c, npp, nhp;
  two_body pp, hp;
  double J;
  length = 0;
  X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    X_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhp[chan1];
  }
  X_1 = new double[length];
  map_num = new int[6 * length]; // X_2_2, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
  map_index = new int[6 * length];
  for(ind = 0; ind < length; ++ind){
    X_1[ind] = 0.0;
  }
  X_1_length = length;

  length = 0;
  X_2_2_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_2_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.npp1[chan2];
  }
  X_2_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_2[ind] = 0.0;
  }
  X_2_2_length = length;

  length = 0;
  X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    X_2_3_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.npp1[chan2];
  }
  X_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X_2_3[ind] = 0.0;
  }
  X_2_3_length = length;

  length = 0;
  X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_1_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhpp[chan3];
  }
  X1_3_1 = new double[length];
  X_3_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_1[ind] = 0.0;
    X_3_1[ind] = 0.0;
  }
  X_3_1_length = length;

  length = 0;
  X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhpp[chan3];
  }
  X1_3_2 = new double[length];
  X_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_2[ind] = 0.0;
    X_3_2[ind] = 0.0;
  }
  X_3_2_length = length;

  length = 0;
  X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_3_index[chan3] = length;
    length += Chan.nppp[chan3] * Chan.nh[chan3];
  }
  X1_3_3 = new double[length];
  X_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_3[ind] = 0.0;
    X_3_3[ind] = 0.0;
  }
  X_3_3_length = length;

  length = 0;
  X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    X_3_4_index[chan3] = length;
    length += Chan.npph[chan3] * Chan.np[chan3];
  }
  X1_3_4 = new double[length];
  X_3_4 = new double[length];
  for(ind = 0; ind < length; ++ind){
    X1_3_4[ind] = 0.0;
    X_3_4[ind] = 0.0;
  }
  X_3_4_length = length;

  length = 0;
  count = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhp = Chan.nhp[chan1];
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
	pp = Chan.pp_state(chan1, pp_ind);
	hp = Chan.hp_state(chan1, hp_ind);
	a = pp.v1;
	b = pp.v2;
	i = hp.v1;
	c = hp.v2;
	// X_2_2, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4_count(Parameters,Space, map_index,map_num, 6,2,0,length,count, a,b,i,c);
	Map_4_count(Parameters,Space, map_index,map_num, 6,3,1,length,count, a,b,i,c);
	Map_4_count(Parameters,Space, map_index,map_num, 6,5,2,length,count, a,b,i,c);
	Map_4_count(Parameters,Space, map_index,map_num, 6,6,3,length,count, a,b,i,c);
	Map_4_count(Parameters,Space, map_index,map_num, 6,7,4,length,count, a,b,i,c);
	Map_4_count(Parameters,Space, map_index,map_num, 6,8,5,length,count, a,b,i,c);
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
    nhp = Chan.nhp[chan1];
    if(Parameters.basis != "finite_J"){ J = 0.0; }
    else{ J = 0.5 * Chan.qnums1[chan1].j; }
    for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
      for(int hp_ind = 0; hp_ind < nhp; ++hp_ind){
	pp = Chan.pp_state(chan1, pp_ind);
	hp = Chan.hp_state(chan1, hp_ind);
	a = pp.v1;
	b = pp.v2;
	i = hp.v1;
	c = hp.v2;
	// X_2_2, X_2_3, X_3_1, X_3_2, X_3_3, X_3_4
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,2,0,count, Chan.ph1_map,Chan.pp1_map, Chan.npp1, a,b,i,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,3,1,count, Chan.ph1_map,Chan.pp1_map, Chan.npp1, a,b,i,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,5,2,count, Chan.p_map,Chan.hpp_map, Chan.nhpp, a,b,i,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,6,3,count, Chan.p_map,Chan.hpp_map, Chan.nhpp, a,b,i,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,7,4,count, Chan.ppp_map,Chan.h_map, Chan.nh, a,b,i,c,J);
	Map_4(Parameters,Space, map_chan,map_ind,map_fac1,map_fac2, map_index,map_num, 6,8,5,count, Chan.pph_map,Chan.p_map, Chan.np, a,b,i,c,J);
	++count;
      }
    }
  }
}

void X_pphp::X1_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  for(ind = 0; ind < X_3_1_length; ++ind){ X1_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X1_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X1_3_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_4_length; ++ind){ X1_3_4[ind] = 0.0; }
}

void X_pphp::X_zero(Input_Parameters &Parameters, Channels &Chan, bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < X_1_length; ++ind){ X_1[ind] = 0.0; }
  }
  for(ind = 0; ind < X_2_2_length; ++ind){ X_2_2[ind] = 0.0; }
  for(ind = 0; ind < X_2_3_length; ++ind){ X_2_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_1_length; ++ind){ X_3_1[ind] = 0.0; }
  for(ind = 0; ind < X_3_2_length; ++ind){ X_3_2[ind] = 0.0; }
  for(ind = 0; ind < X_3_3_length; ++ind){ X_3_3[ind] = 0.0; }
  for(ind = 0; ind < X_3_4_length; ++ind){ X_3_4[ind] = 0.0; }
}

void X_pphp::X1_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  double *X1_1temp = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      x += map_fac2[index] * X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      x += map_fac2[index] * X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      x += map_fac2[index] * X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    X1_1temp[ind] = x;
  }
  X1_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X1_1temp[ind];
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      X1_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      X1_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      X1_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 5]; ++n){
      index = map_index[6*ind + 5] + n;
      X1_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
  delete[] X1_1temp;
}

void X_pphp::X_gather(Input_Parameters &Parameters, Channels &Chan)
{
  int index;
  double x;
  for(int ind = 0; ind < X_1_length; ++ind){
    x = 0.0;
    x += X_1[ind];
    for(int n = 0; n < map_num[6*ind]; ++n){
      index = map_index[6*ind] + n;
      x += map_fac2[index] * X_2_2[X_2_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 1]; ++n){
      index = map_index[6*ind + 1] + n;
      x += map_fac2[index] * X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      x += map_fac2[index] * X_3_1[X_3_1_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      x += map_fac2[index] * X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      x += map_fac2[index] * X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]];
    }
    for(int n = 0; n < map_num[6*ind + 5]; ++n){
      index = map_index[6*ind + 5] + n;
      x += map_fac2[index] * X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]];
    }
    X_1[ind] = x;
  }
  X_zero(Parameters, Chan, true);

  for(int ind = 0; ind < X_1_length; ++ind){
    x = X_1[ind];
    for(int n = 0; n < map_num[6*ind]; ++n){
      index = map_index[6*ind] + n;
      X_2_2[X_2_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 1]; ++n){
      index = map_index[6*ind + 1] + n;
      X_2_3[X_2_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 2]; ++n){
      index = map_index[6*ind + 2] + n;
      X_3_1[X_3_1_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 3]; ++n){
      index = map_index[6*ind + 3] + n;
      X_3_2[X_3_2_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 4]; ++n){
      index = map_index[6*ind + 4] + n;
      X_3_3[X_3_3_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
    for(int n = 0; n < map_num[6*ind + 5]; ++n){
      index = map_index[6*ind + 5] + n;
      X_3_4[X_3_4_index[map_chan[index]] + map_ind[index]] += map_fac1[index] * x;
    }
  }
}

void X_pphp::delete_struct(Input_Parameters &Parameters, Channels &Chan)
{
  delete[] X_1;
  delete[] X_2_2;
  delete[] X_2_3;
  delete[] X1_3_1;
  delete[] X_3_1;
  delete[] X1_3_2;
  delete[] X_3_2;
  delete[] X1_3_3;
  delete[] X_3_3;
  delete[] X1_3_4;
  delete[] X_3_4;

  delete[] X_1_index;
  delete[] X_2_2_index;
  delete[] X_2_3_index;
  delete[] X_3_1_index;
  delete[] X_3_2_index;
  delete[] X_3_3_index;
  delete[] X_3_4_index;

  delete[] map_index;
  delete[] map_chan;
  delete[] map_ind;
  delete[] map_num;
  delete[] map_fac1;
  delete[] map_fac2;
}

void Read_Matrix_Elements_J(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1 = PATH + Parameters.MatrixElements + ".int"; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME, hom, r2, p2; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, coupJ, coupT, par; // interaction file contents
  int chan1, ind, key1, key2, key3, key4;
  State tb;

  ME = HF_Matrix_Elements(Chan);

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while (number != "Total"){ 
    getline(interaction, interactionline);
    interactionstream.str(interactionline);
    interactionstream >> number;
  }
  index1 = interactionline.find_first_of("0123456789");
  index2 = interactionline.find_last_of("0123456789");
  NumElements = std::atoi( interactionline.substr(index1, index2 - index1 + 1).c_str() );

  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while(number != "Tz"){
    getline(interaction, interactionline);
    interactionstream.str(interactionline);
    interactionstream >> number;
  }

  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> coupT >> par >> coupJ >> shell1 >> shell2 >> shell3 >> shell4 >> TBME >> hom >> r2 >> p2;
    //TBME *= Parameters.tbstrength;
    //std::cout << coupT << " " << par << " " << coupJ << " " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if(shell1 == shell2){ TBME *= std::sqrt(2.0); } // !! check
    if(shell3 == shell4){ TBME *= std::sqrt(2.0); } // !! check
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    tb.j = coupJ;
    chan1 = Space.ind_2b_dir(Parameters.basis, tb);
    key1 = Chan.tb_map[chan1][Space.hash2(shell1, shell2, tb.j)];
    key2 = Chan.tb_map[chan1][Space.hash2(shell3, shell4, tb.j)];
    key3 = Chan.tb_map[chan1][Space.hash2(shell2, shell1, tb.j)];
    key4 = Chan.tb_map[chan1][Space.hash2(shell4, shell3, tb.j)];
    ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
    ME.V[ind] = TBME;
    ind = ME.Index[chan1] + (key3 * Chan.ntb[chan1] + key2);
    ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
    ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key4);
    ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
    ind = ME.Index[chan1] + (key3 * Chan.ntb[chan1] + key4);
    ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = ME.Index[chan1] + (key2 * Chan.ntb[chan1] + key1);
      ME.V[ind] = TBME;
      ind = ME.Index[chan1] + (key4 * Chan.ntb[chan1] + key1);
      ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
      ind = ME.Index[chan1] + (key2 * Chan.ntb[chan1] + key3);
      ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
      ind = ME.Index[chan1] + (key4 * Chan.ntb[chan1] + key3);
      ME.V[ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    }
  }
  interaction.close();
}

void Read_Matrix_Elements_M(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int chan1, ind, key1, key2, key3, key4;
  State tb;
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;
  std::cout << "!! " << Parameters.tbstrength << std::endl;
  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
    TBME *= Parameters.tbstrength;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    std::cout << "? " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << "  " << TBME << std::endl;
    if(shell1 == shell2 || shell3 == shell4){ continue; }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    chan1 = Space.ind_2b_dir(Parameters.basis, tb);
    key1 = Chan.tb_map[chan1][Space.hash2(shell1, shell2, tb.j)];
    key2 = Chan.tb_map[chan1][Space.hash2(shell3, shell4, tb.j)];
    key3 = Chan.tb_map[chan1][Space.hash2(shell2, shell1, tb.j)];
    key4 = Chan.tb_map[chan1][Space.hash2(shell4, shell3, tb.j)];
    ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
    ME.V[ind] = TBME;
    ind = ME.Index[chan1] + (key3 * Chan.ntb[chan1] + key2);
    ME.V[ind] = -1.0 * TBME;
    ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key4);
    ME.V[ind] = -1.0 * TBME;
    ind = ME.Index[chan1] + (key3 * Chan.ntb[chan1] + key4);
    ME.V[ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = ME.Index[chan1] + (key2 * Chan.ntb[chan1] + key1);
      ME.V[ind] = TBME;
      ind = ME.Index[chan1] + (key4 * Chan.ntb[chan1] + key1);
      ME.V[ind] = -1.0 * TBME;
      ind = ME.Index[chan1] + (key2 * Chan.ntb[chan1] + key3);
      ME.V[ind] = -1.0 * TBME;
      ind = ME.Index[chan1] + (key4 * Chan.ntb[chan1] + key3);
      ME.V[ind] = TBME;
    }
  }
  interaction.close();
}

void Read_Matrix_Elements_QD(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::cout << "Computing Coulomb Matrix Elements for QD..." << std::endl;
  struct timespec time1, time2;
  double elapsed0 = 0.0;
  ME = HF_Matrix_Elements(Chan);
  clock_gettime(CLOCK_MONOTONIC, &time1);
  
  int pq, rs, p, q, r, s;
  int ntb, ntb0;
  int ind1, ind2, ind3, ind4;
  int length, length0;
  two_body *tbvec0;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tb_vec[Chan.tb_index[chan] + i].v1 < Chan.tb_vec[Chan.tb_index[chan] + i].v2){
	++ntb0;
      }
    }
    tbvec0 = new two_body[ntb0];
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tb_vec[Chan.tb_index[chan] + i].v1 < Chan.tb_vec[Chan.tb_index[chan] + i].v2){
	tbvec0[ntb0].v1 = Chan.tb_vec[Chan.tb_index[chan] + i].v1;
	tbvec0[ntb0].v2 = Chan.tb_vec[Chan.tb_index[chan] + i].v2;
	++ntb0;
      }
    }
    
    length = int(0.5 * ntb0 * (ntb0 + 1));
    #pragma omp parallel private(pq, rs, p, q, r, s, length0, ind1, ind2, ind3, ind4, TBME)
    {
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < length; ++pqrs){
	pq = std::floor((2*ntb0 - 1 - std::sqrt(1 + 4*ntb0 + 4*ntb0*ntb0 - 8*pqrs))/2) + 1;
	length0 = int(0.5 * pq * (2*ntb0 - pq + 1));
	rs = int(pq + pqrs - length0);
	p = tbvec0[pq].v1;
	q = tbvec0[pq].v2;
	r = tbvec0[rs].v1;
	s = tbvec0[rs].v2;
	ind1 = Chan.tb_map[chan][Space.hash2(p, q, Chan.qnums1[chan].j)];
	ind2 = Chan.tb_map[chan][Space.hash2(r, s, Chan.qnums1[chan].j)];
	ind3 = Chan.tb_map[chan][Space.hash2(q, p, Chan.qnums1[chan].j)];
	ind4 = Chan.tb_map[chan][Space.hash2(s, r, Chan.qnums1[chan].j)];
	TBME = Coulomb_HO(Parameters, Space, p, q, r, s);
	ME.V[ME.Index[chan] + (ind1 * ntb + ind2)] = TBME;
	ME.V[ME.Index[chan] + (ind3 * ntb + ind2)] = -1.0 * TBME;
	ME.V[ME.Index[chan] + (ind1 * ntb + ind4)] = -1.0 * TBME;
	ME.V[ME.Index[chan] + (ind3 * ntb + ind4)] = TBME;
	if(ind1 != ind2){
	  ME.V[ME.Index[chan] + (ind2 * ntb + ind1)] = TBME;
	  ME.V[ME.Index[chan] + (ind4 * ntb + ind1)] = -1.0 * TBME;
	  ME.V[ME.Index[chan] + (ind2 * ntb + ind3)] = -1.0 * TBME;
	  ME.V[ME.Index[chan] + (ind4 * ntb + ind3)] = TBME;
	}
      }
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "!! Runtime = " << elapsed0 << " sec. " << std::endl;
}

void Read_QD_ME_From_File(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  std::ifstream interaction;	// interaction file
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + "coulomb-ho2d-elements-20-shells.dat";

  interaction.open(fullpath1.c_str(), std::ios::binary);
  if(interaction.is_open()){
    //get length of file
    interaction.seekg(0, interaction.end);
    int length = interaction.tellg();
    interaction.seekg(0, interaction.beg);

    char *buffer = new char[length];
    interaction.read(buffer, length);
    interaction.close();
    length /= 16;

    unsigned int neg1 = 128;
    unsigned int neg2 = 4294967040;
    int nmax = Parameters.Shells;
    #pragma omp parallel
    {
      double TBME;
      unsigned int n1, n2, n3, n4;
      int ml1, ml2, ml3, ml4;
      int key, key1, key2;
      State tb;
      int ind, chan1;
      State statep, stateq, stater, states;
      int p, q, r, s;
      statep.t = -1, stateq.t = -1, stater.t = -1, states.t = -1;

      #pragma omp for schedule(static)
      for(int i = 0; i < length; ++i){
	n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	ml1 = 0, ml2 = 0, ml3 = 0, ml4 = 0;
	
	n1 = *(buffer + 16*i);
	ml1 = *(buffer + 16*i + 1);
	if((neg1 & ml1) != 0){ ml1 = (ml1 | neg2); }// ml1 < 0
	n2 = *(buffer + 16*i + 2);
	ml2 = *(buffer + 16*i + 3);
	if((neg1 & ml2) != 0){ ml2 = (ml2 | neg2); }// ml2 < 0
	n3 = *(buffer + 16*i + 4);
	ml3 = *(buffer + 16*i + 5);
	if((neg1 & ml3) != 0){ ml3 = (ml3 | neg2); }// ml3 < 0
	n4 = *(buffer + 16*i + 6);
	ml4 = *(buffer + 16*i + 7);
	if((neg1 & ml4) != 0){ ml4 = (ml4 | neg2); }// ml4 < 0
	TBME = *(double*)(buffer + 16*i + 8);
	if(int(2*n1 + abs(ml1)) >= nmax || int(2*n2 + abs(ml2)) >= nmax || int(2*n3 + abs(ml3)) >= nmax || int(2*n4 + abs(ml4)) >= nmax){ continue; }
	TBME *= std::sqrt(Parameters.density);
	statep.n = n1;
	statep.ml = ml1;
	stateq.n = n2;
	stateq.ml = ml2;
	stater.n = n3;
	stater.ml = ml3;
	states.n = n4;
	states.ml = ml4;
	for(int s1 = -1; s1 <= 1; s1 += 2){
	  statep.m = s1;
	  key = Space.ind_state(Parameters.basis, statep);
	  p = Space.map_state[key];
	  for(int s2 = -1; s2 <= 1; s2 += 2){
	    stateq.m = s2;
	    key = Space.ind_state(Parameters.basis, stateq);
	    q = Space.map_state[key];
	    if(p == q){ continue; }
	    for(int s3 = -1; s3 <= 1; s3 += 2){
	      stater.m = s3;
	      key = Space.ind_state(Parameters.basis, stater);
	      r = Space.map_state[key];
	      if(s3 != s1){ continue; }
	      for(int s4 = -1; s4 <= 1; s4 += 2){
		states.m = s4;
		key = Space.ind_state(Parameters.basis, states);
		s = Space.map_state[key];
		if(r == s || s4 != s2){ continue; }
		
		plus(tb, Space.qnums[p], Space.qnums[q]);
		chan1 = Space.ind_2b_dir(Parameters.basis, tb);

		// C(p1q2r3s4) -> <p1q2 || r3s4>
		key1 = Chan.tb_map[chan1][Space.hash2(p, q, tb.j)];
		key2 = Chan.tb_map[chan1][Space.hash2(r, s, tb.j)];
		ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		ME.V[ind] += TBME;
		// C(p1q2r3s4) -> -<p1q2 || s4r3>
		key2 = Chan.tb_map[chan1][Space.hash2(s, r, tb.j)];
		ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		ME.V[ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(q2p1s4r3) -> <q2p1 || s4r3>
		  key1 = Chan.tb_map[chan1][Space.hash2(q, p, tb.j)];
		  key2 = Chan.tb_map[chan1][Space.hash2(s, r, tb.j)];
		  ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  ME.V[ind] += TBME;
		  // C(p1q2r3s4) = C(q2p1s4r3) -> -<q2p1 || r3s4>
		  key2 = Chan.tb_map[chan1][Space.hash2(r, s, tb.j)];
		  ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  ME.V[ind] -= TBME;
		}
		if(((n1 == n3 && ml1 == ml3) && (n2 == n4 && ml2 == ml4)) || ((n1 == n4 && ml1 == ml4) && (n2 == n3 && ml2 == ml3))){ continue; }
		// C(p1q2r3s4) = C(r3s4p1q2) -> <r3s4 || p1q2>
		key1 = Chan.tb_map[chan1][Space.hash2(r, s, tb.j)];
		key2 = Chan.tb_map[chan1][Space.hash2(p, q, tb.j)];
		ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		ME.V[ind] += TBME;
		// C(p1q2r3s4) = C(r3s4p1q2) -> -<r3s4 || q2p1>
		key2 = Chan.tb_map[chan1][Space.hash2(q, p, tb.j)];
		ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		ME.V[ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(s4r3q2p1) -> <s4r3 || q2p1>
		  key1 = Chan.tb_map[chan1][Space.hash2(s, r, tb.j)];
		  key2 = Chan.tb_map[chan1][Space.hash2(q, p, tb.j)];
		  ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  ME.V[ind] += TBME;
		  // C(p1q2r3s4) = C(s4r3q2p1) -> -<s4r3 || p1q2>
		  key2 = Chan.tb_map[chan1][Space.hash2(p, q, tb.j)];
		  ind = ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  ME.V[ind] -= TBME;
		}
	      }
	    }
	  }
	}
      }
    }
    delete[] buffer;
  }
}

void Coulomb_Inf_Matrix_Elements(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  double L = pow(Space.num_hol/Parameters.density, 1.0/3.0);
  int h1, h2, h3, h4, p1, p2, p3, p4;
  int nhh, npp, nh, np, npph, nhhp, nhp1, nph1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  one_body h, p;
  two_body pp1, pp2, hh1, hh2, hp1, hp2, ph1;
  three_body pph, hhp;
  double TBME;

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    chan_ind1 = Ints.Vpppp.V_1_index[chan1];
    chan_ind2 = Ints.Vpphh.V_1_index[chan1];
    chan_ind3 = Ints.Vhhpp.V_1_index[chan1];
    chan_ind4 = Ints.Vhhhh.V_1_index[chan1];
    for(int pp1_ind = 0; pp1_ind < npp; ++pp1_ind){
      pp1 = Chan.pp_state(chan1, pp1_ind);
      p1 = pp1.v1;
      p2 = pp1.v2;
      if(p1 == p2){ continue; }
      for(int pp2_ind = pp1_ind; pp2_ind < npp; ++pp2_ind){
	pp2 = Chan.pp_state(chan1, pp2_ind);
	p3 = pp2.v1;
	p4 = pp2.v2;
	if(p3 == p4){ continue; }
	TBME = Coulomb_Inf(Space, p1, p2, p3, p4, L);
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp1_ind + pp2_ind)] = TBME;
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp2_ind + pp1_ind)] = TBME;
      }
      for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	h1 = hh1.v1;
	h2 = hh1.v2;
	if(h1 == h2){ continue; }
	TBME = Coulomb_Inf(Space, p1, p2, h1, h2, L);
	Ints.Vpphh.V_1[chan_ind2 + (nhh*pp1_ind + hh1_ind)] = TBME;
	Ints.Vhhpp.V_1[chan_ind3 + (npp*hh1_ind + pp1_ind)] = TBME;
      }
    }
    for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
      hh1 = Chan.hh_state(chan1, hh1_ind);
      h1 = hh1.v1;
      h2 = hh1.v2;
      if(h1 == h2){ continue; }
      for(int hh2_ind = hh1_ind; hh2_ind < nhh; ++hh2_ind){
	hh2 = Chan.hh_state(chan1, hh2_ind);
	h3 = hh1.v1;
	h4 = hh1.v2;
	if(h3 == h4){ continue; }
	TBME = Coulomb_Inf(Space, h1, h2, h3, h4, L);
	Ints.Vhhhh.V_1[chan_ind4 + (nhh*hh1_ind + hh2_ind)] = TBME;
	Ints.Vhhhh.V_1[chan_ind4 + (nhh*hh2_ind + hh1_ind)] = TBME;
      }
    }
  }

  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
      hp1 = Chan.hp1_state(chan2, hp1_ind);
      h1 = hp1.v1;
      p2 = hp1.v2;
      for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
	ph1 = Chan.ph1_state(chan2, ph1_ind);
	p1 = ph1.v1;
	h2 = ph1.v2;
	if(h1 == h2 || p1 == p2){ continue; }
	TBME = Coulomb_Inf(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_2_1[chan_ind1 + (nph1*hp1_ind + ph1_ind)] = TBME;
      }
      for(int hp2_ind = hp1_ind; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	h2 = hp2.v1;
	p1 = hp2.v2;
	TBME = Coulomb_Inf(Space, h1, p1, h2, p2, L);
	Ints.Vhphp.V_2_1[chan_ind2 + (nhp1*hp1_ind + hp2_ind)] = TBME;
	Ints.Vhphp.V_2_1[chan_ind2 + (nhp1*hp2_ind + hp1_ind)] = TBME;
      }
    }
  }
  
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nh = Chan.nh[chan3];
    nhhp = Chan.nhhp[chan3];
    npph = Chan.npph[chan3];
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan3];
    chan_ind2 = Ints.Vhhpp.V_3_3_index[chan3];
    for(int h_ind = 0; h_ind < nh; ++h_ind){
      h = Chan.h_state(chan3, h_ind);
      h1 = h.v1;
      for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
	pph = Chan.pph_state(chan3, pph_ind);
	p1 = pph.v1;
	p2 = pph.v2;
	h2 = pph.v3;
	if(h1 == h2 || p1 == p2){ continue; }
	TBME = Coulomb_Inf(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_3_1[chan_ind1 + (npph*h_ind + pph_ind)] = TBME;
      }
    }
    for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
      hhp = Chan.hhp_state(chan3, hhp_ind);
      h1 = hhp.v1;
      h2 = hhp.v2;
      p2 = hhp.v3;
      if(h1 == h2){ continue; }
      for(int p_ind = 0; p_ind < np; ++p_ind){
	p = Chan.p_state(chan3, p_ind);
	p1 = p.v1;
	if(p1 == p2){ continue; }
	TBME = Coulomb_Inf(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_3_3[chan_ind2 + (np*hhp_ind + p_ind)] = TBME;
      }
    }
  }
}

double Coulomb_Inf(Model_Space &Space, int &qi, int &qj, int &qk, int &ql, double &L)
{
  double term = 0.0;
  double e_sq = hbarc_HartA * fine_struct;
  double prefactor = e_sq/(L*L*L);
  double qSquared1;
  double kX1, kY1, kZ1;
  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx || 
     Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny || 
     Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qk].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qk].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qk].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-15){ term += 0.0; }
    else{ term += 4.0 * prefactor * M_PI / qSquared1; }
  }
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-15){ term -= 0.0; }
    else{ term -= 4.0 * prefactor * M_PI / qSquared1; }
  }
  return term;
}

void Minnesota_Matrix_Elements(Input_Parameters &Parameters, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  int h1, h2, h3, h4, p1, p2, p3, p4;
  int nhh, npp, nh, np, npph, nhhp, nhp1, nph1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  one_body h, p;
  two_body pp1, pp2, hh1, hh2, hp1, hp2, ph1;
  three_body pph, hhp;
  double TBME;

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    chan_ind1 = Ints.Vpppp.V_1_index[chan1];
    chan_ind2 = Ints.Vpphh.V_1_index[chan1];
    chan_ind3 = Ints.Vhhpp.V_1_index[chan1];
    chan_ind4 = Ints.Vhhhh.V_1_index[chan1];
    for(int pp1_ind = 0; pp1_ind < npp; ++pp1_ind){
      pp1 = Chan.pp_state(chan1, pp1_ind);
      p1 = pp1.v1;
      p2 = pp1.v2;
      if(p1 == p2){ continue; }
      for(int pp2_ind = pp1_ind; pp2_ind < npp; ++pp2_ind){
	pp2 = Chan.pp_state(chan1, pp2_ind);
	p3 = pp2.v1;
	p4 = pp2.v2;
	if(p3 == p4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, p1, p2, p3, p4, L);
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp1_ind + pp2_ind)] = TBME;
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp2_ind + pp1_ind)] = TBME;
      }
      for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	h1 = hh1.v1;
	h2 = hh1.v2;
	if(h1 == h2){ continue; }
	TBME = vint_Minnesota_Momentum(Space, p1, p2, h1, h2, L);
	Ints.Vpphh.V_1[chan_ind2 + (nhh*pp1_ind + hh1_ind)] = TBME;
	Ints.Vhhpp.V_1[chan_ind3 + (npp*hh1_ind + pp1_ind)] = TBME;
      }
    }
    for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
      hh1 = Chan.hh_state(chan1, hh1_ind);
      h1 = hh1.v1;
      h2 = hh1.v2;
      if(h1 == h2){ continue; }
      for(int hh2_ind = hh1_ind; hh2_ind < nhh; ++hh2_ind){
	hh2 = Chan.hh_state(chan1, hh2_ind);
	h3 = hh1.v1;
	h4 = hh1.v2;
	if(h3 == h4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, h1, h2, h3, h4, L);
	Ints.Vhhhh.V_1[chan_ind4 + (nhh*hh1_ind + hh2_ind)] = TBME;
	Ints.Vhhhh.V_1[chan_ind4 + (nhh*hh2_ind + hh1_ind)] = TBME;
      }
    }
  }
  
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
      hp1 = Chan.hp1_state(chan2, hp1_ind);
      h1 = hp1.v1;
      p2 = hp1.v2;
      for(int ph1_ind = 0; ph1_ind < nph1; ++ph1_ind){
	ph1 = Chan.ph1_state(chan2, ph1_ind);
	p1 = ph1.v1;
	h2 = ph1.v2;
	if(h1 == h2 || p1 == p2){ continue; }
	TBME = vint_Minnesota_Momentum(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_2_1[chan_ind1 + (nph1*hp1_ind + ph1_ind)] = TBME;
      }
      for(int hp2_ind = hp1_ind; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	h2 = hp2.v1;
	p1 = hp2.v2;
	TBME = vint_Minnesota_Momentum(Space, h1, p1, h2, p2, L);
	Ints.Vhphp.V_2_1[chan_ind2 + (nhp1*hp1_ind + hp2_ind)] = TBME;
	Ints.Vhphp.V_2_1[chan_ind2 + (nhp1*hp2_ind + hp1_ind)] = TBME;
      }
    }
  }
  
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nh = Chan.nh[chan3];
    nhhp = Chan.nhhp[chan3];
    npph = Chan.npph[chan3];
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan3];
    chan_ind2 = Ints.Vhhpp.V_3_3_index[chan3];
    for(int h_ind = 0; h_ind < nh; ++h_ind){
      h = Chan.h_state(chan3, h_ind);
      h1 = h.v1;
      for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
	pph = Chan.pph_state(chan3, pph_ind);
	p1 = pph.v1;
	p2 = pph.v2;
	h2 = pph.v3;
	if(h1 == h2 || p1 == p2){ continue; }
	TBME = vint_Minnesota_Momentum(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_3_1[chan_ind1 + (npph*h_ind + pph_ind)] = TBME;
      }
    }
    for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
      hhp = Chan.hhp_state(chan3, hhp_ind);
      h1 = hhp.v1;
      h2 = hhp.v2;
      p2 = hhp.v3;
      if(h1 == h2){ continue; }
      for(int p_ind = 0; p_ind < np; ++p_ind){
	p = Chan.p_state(chan3, p_ind);
	p1 = p.v1;
	if(p1 == p2){ continue; }
	TBME = vint_Minnesota_Momentum(Space, h1, h2, p1, p2, L);
	Ints.Vhhpp.V_3_3[chan_ind2 + (np*hhp_ind + p_ind)] = TBME;
      }
    }
  }
}

int kron_del(int &i, int &j)
{
  if(i != j){ return 0; }
  return 1;
}
int spinExchangeMtxEle(int &i, int &j, int &k, int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}

// Minnesota Potential for momentum basis
double vint_Minnesota_Momentum(Model_Space &Space, int &qi, int &qj, int &qk, int &ql, double &L)
{
  double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
  double V_0R, V_0T, V_0S;
  double kappa_R, kappa_T, kappa_S;
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx){ return 0.0; }
  if(Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny){ return 0.0; }
  if(Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m + Space.qnums[qj].m != Space.qnums[qk].m + Space.qnums[ql].m){ return 0.0; }
  if(Space.qnums[qi].t + Space.qnums[qj].t != Space.qnums[qk].t + Space.qnums[ql].t){ return 0.0; }

  kX1 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[qk].nx + Space.qnums[ql].nx);
  kY1 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[qk].ny + Space.qnums[ql].ny);
  kZ1 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[qk].nz + Space.qnums[ql].nz);

  kX2 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[ql].nx + Space.qnums[qk].nx);
  kY2 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[ql].ny + Space.qnums[qk].ny);
  kZ2 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[ql].nz + Space.qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[qk].m, Space.qnums[ql].m);
  isoSpinEx1 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[qk].t, Space.qnums[ql].t);

  spinEx2 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[ql].m, Space.qnums[qk].m);
  isoSpinEx2 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[ql].t, Space.qnums[qk].t);
  
  IsIt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m) * kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsIt1 = spinEx1 * kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m)*kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsIt2 = spinEx2 * kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 + 
    0.25 * (V_T1 - V_S1) * PsIt1 - 
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 - 
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 - 
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 + 
    0.25 * (V_T2 - V_S2) * IsPt2;
}


// Minnesota Potential for momentum basis
double Coulomb_HO(Input_Parameters &Parameters, Model_Space &Space, int &qi, int &qj, int &qk, int &ql)
{
  int g1, g2, g3, g4, G, L;
  double dir = 0.0;
  double exch = 0.0;
  int n1, m1, n2, m2, n3, m3, n4, m4;
  n1 = Space.qnums[qi].n; //1
  m1 = Space.qnums[qi].ml;
  n2 = Space.qnums[qj].n; //2
  m2 = Space.qnums[qj].ml;
  n3 = Space.qnums[ql].n; //4
  m3 = Space.qnums[ql].ml;
  n4 = Space.qnums[qk].n; //3
  m4 = Space.qnums[qk].ml;
  if((m1 + m2 == m3 + m4) && Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = int(j1 + j4 + 0.5*(abs(m1) + m1) + 0.5*(abs(m4) - m4));
	    g2 = int(j2 + j3 + 0.5*(abs(m2) + m2) + 0.5*(abs(m3) - m3));
	    g3 = int(j3 + j2 + 0.5*(abs(m3) + m3) + 0.5*(abs(m2) - m2));
	    g4 = int(j4 + j1 + 0.5*(abs(m4) + m4) + 0.5*(abs(m1) - m1));
	    G = g1 + g2 + g3 + g4;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    double LogRatio2 = logratio2(G);
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    dir += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    dir *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }

  n1 = Space.qnums[qi].n; //1
  m1 = Space.qnums[qi].ml;
  n2 = Space.qnums[qj].n; //2
  m2 = Space.qnums[qj].ml;
  n3 = Space.qnums[qk].n; //3
  m3 = Space.qnums[qk].ml;
  n4 = Space.qnums[ql].n; //4
  m4 = Space.qnums[ql].ml;
  if((m1 + m2 == m3 + m4) && Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = int(j1 + j4 + 0.5*(abs(m1) + m1) + 0.5*(abs(m4) - m4));
	    g2 = int(j2 + j3 + 0.5*(abs(m2) + m2) + 0.5*(abs(m3) - m3));
	    g3 = int(j3 + j2 + 0.5*(abs(m3) + m3) + 0.5*(abs(m2) - m2));
	    g4 = int(j4 + j1 + 0.5*(abs(m4) + m4) + 0.5*(abs(m1) - m1));
	    G = g1 + g2 + g3 + g4;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    double LogRatio2 = logratio2(G);
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    exch += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    exch *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }
  return std::sqrt(Parameters.density)*(dir - exch);
}

void Get_Fock_Matrix(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &HF, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  int p_ind, q_ind, t_ind;
  int key1, key2, chan1, ind1, minj;
  int nob0, nob;
  double fock, fock0;
  State tb;
  int ind3, ind2;
  int chan_ind;
  int vector_ind;

  for(int chan0 = 0; chan0 < HF_Chan.size3; ++chan0){
    vector_ind = HF.vector_index[chan0];
    nob0 = HF_Chan.nob[chan0];
    for(int p = 0; p < nob0; ++p){
      p_ind = HF_Chan.ob_state(chan0, p).v1;
      for(int q = 0; q < nob0; ++q){
	q_ind = HF_Chan.ob_state(chan0, q).v1;
	fock = 0.0;
	fock0 = 0.0;
	for(int t = 0; t < nob0; ++t){
	  t_ind = HF_Chan.ob_state(chan0, t).v1;
	  fock0 += Space.qnums[t_ind].energy * HF.vectors[vector_ind + (p*nob0 + t)] * HF.vectors[vector_ind + (q*nob0 + t)];
	}
	fock += fock0;
	for(int chan = 0; chan < HF_Chan.size3; ++chan){
	  nob = HF_Chan.nob[chan];
	  for(int t = 0; t < nob; ++t){ // Sum over occupied levels
	    t_ind = HF_Chan.ob_state(chan, t).v1;
	    if(t_ind >= Space.num_hol){ continue; } 
	    plus(tb, Space.qnums[p_ind], Space.qnums[t_ind]);
	    minj = abs(HF_Chan.qnums3[chan0].j - HF_Chan.qnums3[chan].j);
	    if(Parameters.basis != "finite_J" && Parameters.basis != "finite_JM"){
	      tb.j = 0;
	      minj = 0;
	    }
	    while(tb.j >= minj){
	      chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	      key1 = HF_Chan.tb_map[chan1][Space.hash2(p_ind, t_ind, tb.j)];
	      key2 = HF_Chan.tb_map[chan1][Space.hash2(q_ind, t_ind, tb.j)];
	      ind1 = HF_ME.Index[chan1] + (key1 * HF_Chan.ntb[chan1] + key2);
	      fock += (tb.j + 1.0)/(HF_Chan.qnums3[chan0].j + 1.0) * HF_ME.V[ind1];
	      tb.j -= 2;
	    }
	  }
	}
	if(p_ind < Space.num_hol){
	  key1 = Chan.h_map[chan0][p_ind];
	  if(q_ind < Space.num_hol){ // hh
	    key2 = Chan.h_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.hh_3_index[chan0];
	    ind3 = key1 * Chan.nh[chan0] + key2;
	    minus(tb, Space.qnums[p_ind], Space.qnums[q_ind]);
	    ind2 = Chan.hh1_map[Chan.ind0][Space.hash2(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.hh_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.hh_2[ind2] = fock;
	  }
	  else{ // hp
	    key2 = Chan.p_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.hp_3_index[chan0];
	    ind3 = key1 * Chan.np[chan0] + key2;
 	    minus(tb, Space.qnums[p_ind], Space.qnums[q_ind]);
	    ind2 = Chan.hp1_map[Chan.ind0][Space.hash2(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.hp_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.hp_2[ind2] = fock;
	  }
	}
	else{
	  key1 = Chan.p_map[chan0][p_ind];
	  if(q_ind < Space.num_hol){ // ph
	    key2 = Chan.h_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.ph_3_index[chan0];
	    ind3 = key1 * Chan.nh[chan0] + key2;
 	    minus(tb, Space.qnums[p_ind], Space.qnums[q_ind]);
	    ind2 = Chan.ph1_map[Chan.ind0][Space.hash2(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.ph_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.ph_2[ind2] = fock;
	  }
	  else{ // pp
	    key2 = Chan.p_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.pp_3_index[chan0];
	    ind3 = key1 * Chan.np[chan0] + key2;
 	    minus(tb, Space.qnums[p_ind], Space.qnums[q_ind]);
	    ind2 = Chan.pp1_map[Chan.ind0][Space.hash2(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.pp_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.pp_2[ind2] = fock;
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  int ntb, length;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    length = ntb * ntb;
    #pragma omp parallel
    {
      int p, q, r, s, tb1_ind, tb2_ind, chan_ind;
      two_body tb1, tb2;
      int ind, chan, key1, key2;
      std::string ptype, qtype, rtype, stype;
      double TBME;
      State tb;
      #pragma omp for schedule(static)
      for(int ind1 = 0; ind1 < length; ++ind1){
	tb2_ind = ind1%ntb;
	tb1_ind = (ind1 - tb2_ind)/ntb;
	tb1 = HF_Chan.tb_state(chan1, tb1_ind);
	tb2 = HF_Chan.tb_state(chan1, tb2_ind);
	p = tb1.v1;
	q = tb1.v2;
	ptype = Space.qnums[p].type;
	qtype = Space.qnums[q].type;
	if(p == q){ continue; }
	r = tb2.v1;
	s = tb2.v2;
	rtype = Space.qnums[r].type;
	stype = Space.qnums[s].type;
	if(r == s){ continue; }
	TBME = HF_ME.V[HF_ME.Index[chan1] + (tb1_ind*ntb + tb2_ind)];
	/*if(p < q && r < s && p <= r){
	  std::cout << std::setprecision(12) << "V_hf: " << p << " " << q << " " << r << " " << s << " = " << TBME << std::endl;
	  }*/
	if(rtype == "particle" && stype == "particle"){
	  if(ptype == "particle" && qtype == "particle"){
	    key1 = Chan.pp_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	    key2 = Chan.pp_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vpppp.V_1_index[chan1];
	    ind = key1 * Chan.npp[chan1] + key2;
	    Ints.Vpppp.V_1[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	    key1 = Chan.ppp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	    key2 = Chan.p_map[chan][r];
	    chan_ind = Ints.Vpppp.V_3_3_index[chan];
	    ind = key1 * Chan.np[chan] + key2;
	    Ints.Vpppp.V_3_3[chan_ind + ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "particle"){
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	    key1 = Chan.p_map[chan][q];
	    key2 = Chan.pph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhppp.V_3_2_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhppp.V_3_2[chan_ind + ind] = TBME;
	    
	    minus(tb, Space.qnums[q], Space.qnums[s]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.pp1_map[chan][Space.hash2(q, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(r, p, tb.j)];
	    chan_ind = Ints.Vhppp.V_2_4_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhppp.V_2_4[chan_ind + ind] = TBME;
	    
	    // Vpphp -> rspq
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	    key1 = Chan.ppp_map[chan][Space.hash3(r, s, q, Chan.qnums1[chan1].j)];
	    key2 = Chan.h_map[chan][p];
	    chan_ind = Ints.Vpphp.V_3_3_index[chan];
	    ind = key1 * Chan.nh[chan] + key2;
	    Ints.Vpphp.V_3_3[chan_ind + ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole"){
	    key1 = Chan.hh_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	    key2 = Chan.pp_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhpp.V_1_index[chan1];
	    ind = key1 * Chan.npp[chan1] + key2;
	    Ints.Vhhpp.V_1[chan_ind + ind] = TBME;
	    chan_ind = Ints.Vpphh.V_1_index[chan1];
	    ind = key2 * Chan.nhh[chan1] + key1;
	    Ints.Vpphh.V_1[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	    key1 = Chan.h_map[chan][p];
	    key2 = Chan.pph_map[chan][Space.hash3(r, s, q, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhpp.V_3_1_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhhpp.V_3_1[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.pph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhpp.V_3_2_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhhpp.V_3_2[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	    key1 = Chan.hhp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	    key2 = Chan.p_map[chan][r];
	    chan_ind = Ints.Vhhpp.V_3_3_index[chan];
	    ind = key1 * Chan.np[chan] + key2;
	    Ints.Vhhpp.V_3_3[chan_ind + ind] = TBME;
	    
	    minus(tb, Space.qnums[p], Space.qnums[s]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(r, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_1_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhpp.V_2_1[chan_ind + ind] = TBME;
	    
	    minus(tb, Space.qnums[p], Space.qnums[r]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhpp.V_2_3[chan_ind + ind] = TBME;
	  } 
	}
	else if(rtype == "hole" && stype == "particle"){
	  if(ptype == "hole" && qtype == "particle"){
	    minus(tb, Space.qnums[p], Space.qnums[s]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, s, tb.j)];
	    key2 = Chan.hp1_map[chan][Space.hash2(r, q, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_1_index[chan];
	    ind = key1 * Chan.nhp1[chan] + key2;
	    Ints.Vhphp.V_2_1[chan_ind + ind] = TBME;
	    
	    minus(tb, Space.qnums[q], Space.qnums[r]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.ph1_map[chan][Space.hash2(q, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, p, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_2_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhphp.V_2_2[chan_ind + ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole"){
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.hph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhhp.V_3_2_index[chan];
	    ind = key1 * Chan.nhph[chan] + key2;
	    Ints.Vhhhp.V_3_2[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	    key1 = Chan.hhp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	    key2 = Chan.h_map[chan][r];
	    chan_ind = Ints.Vhhhp.V_3_3_index[chan];
	    ind = key1 * Chan.nh[chan] + key2;
	    Ints.Vhhhp.V_3_3[chan_ind + ind] = TBME;
	    
	    minus(tb, Space.qnums[p], Space.qnums[r]);
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hh1_map[chan][Space.hash2(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, q, tb.j)];
	    chan_ind = Ints.Vhhhp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhhp.V_2_3[chan_ind + ind] = TBME;
	    
	    // Vhphh -> rspq
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[s]);
	    key1 = Chan.p_map[chan][s];
	    key2 = Chan.hhh_map[chan][Space.hash3(p, q, r, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhphh.V_3_2_index[chan];
	    ind = key1 * Chan.nhhh[chan] + key2;
	    Ints.Vhphh.V_3_2[chan_ind + ind] = TBME;
	  }
	}
	else if(rtype == "hole" && stype == "hole"){
	  if(ptype == "hole" && qtype == "hole"){
	    key1 = Chan.hh_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	    key2 = Chan.hh_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhhh.V_1_index[chan1];
	    ind = key1 * Chan.nhh[chan1] + key2;
	    Ints.Vhhhh.V_1[chan_ind + ind] = TBME;
	    
	    chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.hhh_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	    chan_ind = Ints.Vhhhh.V_3_2_index[chan];
	    ind = key1 * Chan.nhhh[chan] + key2;
	    Ints.Vhhhh.V_3_2[chan_ind + ind] = TBME;
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements_J(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  int ntb, length;
  double J;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    length = ntb * ntb;
    J = 0.5 * HF_Chan.qnums1[chan1].j;
    #pragma omp parallel
    {
      int p, q, r, s, tb1_ind, tb2_ind, chan_ind;
      two_body tb1, tb2;
      double pj, qj, rj, sj, tbj, X;
      int ind, chan, key1, key2, jmin;
      std::string ptype, qtype, rtype, stype;
      double TBME;
      State tb;
      #pragma omp for schedule(static)
      for(int ind1 = 0; ind1 < length; ++ind1){
	tb2_ind = ind1%ntb;
	tb1_ind = (ind1 - tb2_ind)/ntb;
	tb1 = HF_Chan.tb_state(chan1, tb1_ind);
	tb2 = HF_Chan.tb_state(chan1, tb2_ind);
	p = tb1.v1;
	q = tb1.v2;
	ptype = Space.qnums[p].type;
	qtype = Space.qnums[q].type;
	pj = 0.5 * Space.qnums[p].j;
	qj = 0.5 * Space.qnums[q].j;
	r = tb2.v1;
	s = tb2.v2;
	rtype = Space.qnums[r].type;
	stype = Space.qnums[s].type;
	rj = 0.5 * Space.qnums[r].j;
	sj = 0.5 * Space.qnums[s].j;
	TBME = HF_ME.V[HF_ME.Index[chan1] + (tb1_ind*ntb + tb2_ind)];
	/*if(p < q && r < s && p <= r){
	  std::cout << std::setprecision(12) << "V_hf: " << p << " " << q << " " << r << " " << s << ", " << 0.5*HF_Chan.qnums1[chan].j << " = " << TBME << std::endl;
	  }*/
	if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	  key1 = Chan.pp_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	  key2 = Chan.pp_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vpppp.V_1_index[chan1];
	  ind = key1 * Chan.npp[chan1] + key2;
	  Ints.Vpppp.V_1[chan_ind + ind] = TBME;

	  chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	  key1 = Chan.ppp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	  key2 = Chan.p_map[chan][r];
	  chan_ind = Ints.Vpppp.V_3_3_index[chan];
	  ind = key1 * Chan.np[chan] + key2;
	  Ints.Vpppp.V_3_3[chan_ind + ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	  key1 = Chan.hh_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	  key2 = Chan.hh_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhhh.V_1_index[chan1];
	  ind = key1 * Chan.nhh[chan1] + key2;
	  Ints.Vhhhh.V_1[chan_ind + ind] = TBME;

	  chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.hhh_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhhh.V_3_2_index[chan];
	  ind = key1 * Chan.nhhh[chan] + key2;
	  Ints.Vhhhh.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	}
	else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	  minus(tb, Space.qnums[p], Space.qnums[s]);
	  if(Space.qnums[r].j + Space.qnums[q].j < tb.j){ tb.j = Space.qnums[r].j + Space.qnums[q].j; }
	  jmin = abs(Space.qnums[p].j - Space.qnums[s].j);
	  if(abs(Space.qnums[r].j - Space.qnums[q].j) > jmin){ jmin = abs(Space.qnums[r].j - Space.qnums[q].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, s, tb.j)];
	    key2 = Chan.hp1_map[chan][Space.hash2(r, q, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_1_index[chan];
	    ind = key1 * Chan.nhp1[chan] + key2;
	    X = -1.0 * (2.0 * J + 1) * CGC6(pj,qj,J,rj,sj,tbj);
	    Ints.Vhphp.V_2_1[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }

	  minus(tb, Space.qnums[q], Space.qnums[r]);
	  if(Space.qnums[s].j + Space.qnums[p].j < tb.j){ tb.j = Space.qnums[s].j + Space.qnums[p].j; }
	  jmin = abs(Space.qnums[q].j - Space.qnums[r].j);
	  if(abs(Space.qnums[s].j - Space.qnums[p].j) > jmin){ jmin = abs(Space.qnums[s].j - Space.qnums[p].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.ph1_map[chan][Space.hash2(q, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, p, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_2_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = -1.0 * std::pow(-1.0, pj + qj + rj + sj) * (2.0 * J + 1) * CGC6(qj,pj,J,sj,rj,tbj);
	    Ints.Vhphp.V_2_2[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	  key1 = Chan.hh_map[chan1][Space.hash2(p, q, Chan.qnums1[chan1].j)];
	  key2 = Chan.pp_map[chan1][Space.hash2(r, s, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhpp.V_1_index[chan1];
	  ind = key1 * Chan.npp[chan1] + key2;
	  Ints.Vhhpp.V_1[chan_ind + ind] = TBME;
	  chan_ind = Ints.Vpphh.V_1_index[chan1];
	  ind = key2 * Chan.nhh[chan1] + key1;
	  Ints.Vpphh.V_1[chan_ind + ind] = TBME;

	  chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	  key1 = Chan.h_map[chan][p];
	  key2 = Chan.pph_map[chan][Space.hash3(r, s, q, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhpp.V_3_1_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhhpp.V_3_1[chan_ind + ind] = std::pow(-1.0, pj + qj - J) * std::sqrt((2.0*J + 1)/(2.0*pj + 1)) * TBME;

	  chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
 	  key1 = Chan.hhp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	  key2 = Chan.p_map[chan][r];
	  chan_ind = Ints.Vhhpp.V_3_3_index[chan];
	  ind = key1 * Chan.np[chan] + key2;
	  Ints.Vhhpp.V_3_3[chan_ind + ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	   
	  minus(tb, Space.qnums[p], Space.qnums[s]);
	  if(Space.qnums[r].j + Space.qnums[q].j < tb.j){ tb.j = Space.qnums[r].j + Space.qnums[q].j; }
	  jmin = abs(Space.qnums[p].j - Space.qnums[s].j);
	  if(abs(Space.qnums[r].j - Space.qnums[q].j) > jmin){ jmin = abs(Space.qnums[r].j - Space.qnums[q].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(r, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_1_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = -1.0 * (2.0 * J + 1) * CGC6(pj,qj,J,rj,sj,tbj);
	    Ints.Vhhpp.V_2_1[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }

	  chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.pph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhpp.V_3_2_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhhpp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	  
	  minus(tb, Space.qnums[p], Space.qnums[r]);
	  if(Space.qnums[s].j + Space.qnums[q].j < tb.j){ tb.j = Space.qnums[s].j + Space.qnums[q].j; }
	  jmin = abs(Space.qnums[p].j - Space.qnums[r].j);
	  if(abs(Space.qnums[s].j - Space.qnums[q].j) > jmin){ jmin = abs(Space.qnums[s].j - Space.qnums[q].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hp1_map[chan][Space.hash2(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = std::pow(-1.0, rj + sj - J) * (2.0 * J + 1) * CGC6(pj,qj,J,sj,rj,tbj);
	    Ints.Vhhpp.V_2_3[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	}
	else if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	  chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	  key1 = Chan.p_map[chan][q];
	  key2 = Chan.pph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhppp.V_3_2_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhppp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	  
	  minus(tb, Space.qnums[q], Space.qnums[s]);
	  if(Space.qnums[r].j + Space.qnums[p].j < tb.j){ tb.j = Space.qnums[r].j + Space.qnums[p].j; }
	  jmin = abs(Space.qnums[q].j - Space.qnums[s].j);
	  if(abs(Space.qnums[r].j - Space.qnums[p].j) > jmin){ jmin = abs(Space.qnums[r].j - Space.qnums[p].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.pp1_map[chan][Space.hash2(q, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(r, p, tb.j)];
	    chan_ind = Ints.Vhppp.V_2_4_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = std::pow(-1.0, pj + qj - J) * (2.0 * J + 1) * CGC6(qj,pj,J,rj,sj,tbj);
	    Ints.Vhppp.V_2_4[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }	  

	  // Vpphp -> rspq
	  chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	  key1 = Chan.ppp_map[chan][Space.hash3(r, s, q, Chan.qnums1[chan1].j)];
	  key2 = Chan.h_map[chan][p];
	  chan_ind = Ints.Vpphp.V_3_3_index[chan];
	  ind = key1 * Chan.nh[chan] + key2;
	  Ints.Vpphp.V_3_3[chan_ind + ind] = std::pow(-1.0, pj + qj - J) * std::sqrt((2.0*J + 1)/(2.0*pj + 1)) * TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	  chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	  key1 = Chan.hhp_map[chan][Space.hash3(p, q, s, Chan.qnums1[chan1].j)];
	  key2 = Chan.h_map[chan][r];
	  chan_ind = Ints.Vhhhp.V_3_3_index[chan];
	  ind = key1 * Chan.nh[chan] + key2;
	  Ints.Vhhhp.V_3_3[chan_ind + ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	  
	  chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.hph_map[chan][Space.hash3(r, s, p, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhhhp.V_3_2_index[chan];
	  ind = key1 * Chan.nhph[chan] + key2;
	  Ints.Vhhhp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	  
	  minus(tb, Space.qnums[p], Space.qnums[r]);
	  if(Space.qnums[s].j + Space.qnums[q].j < tb.j){ tb.j = Space.qnums[s].j + Space.qnums[q].j; }
	  jmin = abs(Space.qnums[p].j - Space.qnums[r].j);
	  if(abs(Space.qnums[s].j - Space.qnums[q].j) > jmin){ jmin = abs(Space.qnums[s].j - Space.qnums[q].j); }
	  while(tb.j >= jmin){
	    tbj = 0.5 * tb.j;
	    chan = Space.ind_2b_cross(Parameters.basis, tb);
	    key1 = Chan.hh1_map[chan][Space.hash2(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Space.hash2(s, q, tb.j)];
	    chan_ind = Ints.Vhhhp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = std::pow(-1.0, rj + sj - J) * (2.0 * J + 1) * CGC6(pj,qj,J,sj,rj,tbj);
	    Ints.Vhhhp.V_2_3[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }

	  // Vhphh -> rspq
	  chan = Space.ind_1b(Parameters.basis, Space.qnums[s]);
	  key1 = Chan.p_map[chan][s];
	  key2 = Chan.hhh_map[chan][Space.hash3(p, q, r, Chan.qnums1[chan1].j)];
	  chan_ind = Ints.Vhphh.V_3_2_index[chan];
	  ind = key1 * Chan.nhhh[chan] + key2;
	  Ints.Vhphh.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*sj + 1)) * TBME;
	}
      }
    }
  }
}

void Get_Matrix_Elements_JM(Input_Parameters &Parameters, HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Model_Space &Space, Channels &Chan, Interactions &Ints)
{
  double J;
  int ntb, length;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    length = ntb * ntb;
    J = 0.5 * HF_Chan.qnums1[chan1].j;
    #pragma omp parallel
    {
      int p, q, r, s, pind, qind, rind, sind, tb1_ind, tb2_ind, chan_ind;
      two_body tb1, tb2;
      double pj, qj, rj, sj, pm, qm, rm, sm, M;
      double m1, m2;
      int t1, t2;
      double CGC1, CGC2;
      int ind, chan, key1, key2;
      std::string ptype, qtype, rtype, stype;
      double TBME, TBME0;
      State tb;
      #pragma omp for schedule(static)
      for(int ind1 = 0; ind1 < length; ++ind1){
	tb2_ind = ind1%ntb;
	tb1_ind = (ind1 - tb2_ind)/ntb;
	tb1 = HF_Chan.tb_state(chan1, tb1_ind);
	tb2 = HF_Chan.tb_state(chan1, tb2_ind);
	p = tb1.v1;
	q = tb1.v2;
	pj = 0.5 * Space.qnums[Space.shellsm[p][0]].j;
	qj = 0.5 * Space.qnums[Space.shellsm[q][0]].j;
	ptype = Space.qnums[Space.shellsm[p][0]].type;
	qtype = Space.qnums[Space.shellsm[q][0]].type;
	r = tb2.v1;
	s = tb2.v2;
	rj = 0.5 * Space.qnums[Space.shellsm[r][0]].j;
	sj = 0.5 * Space.qnums[Space.shellsm[s][0]].j;
	rtype = Space.qnums[Space.shellsm[r][0]].type;
	stype = Space.qnums[Space.shellsm[s][0]].type;
	TBME0 = HF_ME.V[HF_ME.Index[chan1] + (tb1_ind*ntb + tb2_ind)];

	for(int jz = -HF_Chan.qnums1[chan1].j; jz <= HF_Chan.qnums1[chan1].j; jz+=2){
	  M = 0.5 * jz;
	  for(int p1 = 0; p1 < -1*Space.qnums[Space.shellsm[p][0]].m + 1; ++p1){
	    pind = Space.shellsm[p][p1];
	    pm = 0.5 * Space.qnums[pind].m;
	    for(int q1 = 0; q1 < -1*Space.qnums[Space.shellsm[q][0]].m + 1; ++q1){
	      qind = Space.shellsm[q][q1];
	      qm = 0.5 * Space.qnums[qind].m;
	      m1 = pm + qm;
	      t1 = Space.qnums[pind].t + Space.qnums[qind].t;
	      if(pind == qind){ continue; }
	      for(int r1 = 0; r1 < -1*Space.qnums[Space.shellsm[r][0]].m + 1; ++r1){
		rind = Space.shellsm[r][r1];
		rm = 0.5 * Space.qnums[rind].m;
		for(int s1 = 0; s1 < -1*Space.qnums[Space.shellsm[s][0]].m + 1; ++s1){
		  sind = Space.shellsm[s][s1];
		  sm = 0.5 * Space.qnums[sind].m;
		  qm = 0.5 * Space.qnums[qind].m;
		  m2 = rm + sm;
		  t2 = Space.qnums[rind].t + Space.qnums[sind].t;
		  if(rind == sind){ continue; }
		  
		  if(t1 != t2 || m1 != M || m2 != M){ continue; }
		  CGC1 = CGC(pj, pm, qj, qm, J, M);
		  CGC2 = CGC(rj, rm, sj, sm, J, M);
		  TBME = TBME0 * CGC1 * CGC2;

		  if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		    key1 = Chan.pp_map[chan1][Space.hash2(pind, qind, Chan.qnums1[chan1].j)];
		    key2 = Chan.pp_map[chan1][Space.hash2(rind, sind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vpppp.V_1_index[chan1];
		    ind = key1 * Chan.npp[chan1] + key2;
		    Ints.Vpppp.V_1[chan_ind + ind] = TBME;

		    chan = Space.ind_1b(Parameters.basis, Space.qnums[rind]);
		    key1 = Chan.ppp_map[chan][Space.hash3(pind, qind, sind, Chan.qnums1[chan1].j)];
		    key2 = Chan.p_map[chan][rind];
		    chan_ind = Ints.Vpppp.V_3_3_index[chan];
		    ind = key1 * Chan.np[chan] + key2;
		    Ints.Vpppp.V_3_3[chan_ind + ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
		    key1 = Chan.hh_map[chan1][Space.hash2(pind, qind, Chan.qnums1[chan1].j)];
		    key2 = Chan.hh_map[chan1][Space.hash2(rind, sind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhhh.V_1_index[chan1];
		    ind = key1 * Chan.nhh[chan1] + key2;
		    Ints.Vhhhh.V_1[chan_ind + ind] = TBME;

		    chan = Space.ind_1b(Parameters.basis, Space.qnums[qind]);
		    key1 = Chan.h_map[chan][qind];
		    key2 = Chan.hhh_map[chan][Space.hash3(rind, sind, pind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhhh.V_3_2_index[chan];
		    ind = key1 * Chan.nhhh[chan] + key2;
		    Ints.Vhhhh.V_3_2[chan_ind + ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
		    minus(tb, Space.qnums[pind], Space.qnums[sind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.hp1_map[chan][Space.hash2(pind, sind, tb.j)];
		    key2 = Chan.hp1_map[chan][Space.hash2(rind, qind, tb.j)];
		    chan_ind = Ints.Vhphp.V_2_1_index[chan];
		    ind = key1 * Chan.nhp1[chan] + key2;
		    Ints.Vhphp.V_2_1[chan_ind + ind] = TBME;

		    minus(tb, Space.qnums[qind], Space.qnums[rind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.ph1_map[chan][Space.hash2(qind, rind, tb.j)];
		    key2 = Chan.ph1_map[chan][Space.hash2(sind, pind, tb.j)];
		    chan_ind = Ints.Vhphp.V_2_2_index[chan];
		    ind = key1 * Chan.nph1[chan] + key2;
		    Ints.Vhphp.V_2_2[chan_ind + ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
		    key1 = Chan.hh_map[chan1][Space.hash2(pind, qind, Chan.qnums1[chan1].j)];
		    key2 = Chan.pp_map[chan1][Space.hash2(rind, sind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhpp.V_1_index[chan1];
		    ind = key1 * Chan.npp[chan1] + key2;
		    Ints.Vhhpp.V_1[chan_ind + ind] = TBME;
		    chan_ind = Ints.Vpphh.V_1_index[chan1];
		    ind = key2 * Chan.nhh[chan1] + key1;
		    Ints.Vpphh.V_1[chan_ind + ind] = TBME;

		    chan = Space.ind_1b(Parameters.basis, Space.qnums[pind]);
		    key1 = Chan.h_map[chan][pind];
		    key2 = Chan.pph_map[chan][Space.hash3(rind, sind, qind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhpp.V_3_1_index[chan];
		    ind = key1 * Chan.npph[chan] + key2;
		    Ints.Vhhpp.V_3_1[chan_ind + ind] = TBME;
	    
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[rind]);
		    key1 = Chan.hhp_map[chan][Space.hash3(pind, qind, sind, Chan.qnums1[chan1].j)];
		    key2 = Chan.p_map[chan][rind];
		    chan_ind = Ints.Vhhpp.V_3_3_index[chan];
		    ind = key1 * Chan.np[chan] + key2;
		    Ints.Vhhpp.V_3_3[chan_ind + ind] = TBME;
	    
		    minus(tb, Space.qnums[pind], Space.qnums[sind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.hp1_map[chan][Space.hash2(pind, sind, tb.j)];
		    key2 = Chan.ph1_map[chan][Space.hash2(rind, qind, tb.j)];
		    chan_ind = Ints.Vhhpp.V_2_1_index[chan];
		    ind = key1 * Chan.nph1[chan] + key2;
		    Ints.Vhhpp.V_2_1[chan_ind + ind] = TBME;

		    chan = Space.ind_1b(Parameters.basis, Space.qnums[qind]);
		    key1 = Chan.h_map[chan][qind];
		    key2 = Chan.pph_map[chan][Space.hash3(rind, sind, pind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhpp.V_3_2_index[chan];
		    ind = key1 * Chan.npph[chan] + key2;
		    Ints.Vhhpp.V_3_2[chan_ind + ind] = TBME;
		    
		    minus(tb, Space.qnums[pind], Space.qnums[rind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.hp1_map[chan][Space.hash2(pind, rind, tb.j)];
		    key2 = Chan.ph1_map[chan][Space.hash2(sind, qind, tb.j)];
		    chan_ind = Ints.Vhhpp.V_2_3_index[chan];
		    ind = key1 * Chan.nph1[chan] + key2;
		    Ints.Vhhpp.V_2_3[chan_ind + ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[qind]);
		    key1 = Chan.p_map[chan][qind];
		    key2 = Chan.pph_map[chan][Space.hash3(rind, sind, pind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhppp.V_3_2_index[chan];
		    ind = key1 * Chan.npph[chan] + key2;
		    Ints.Vhppp.V_3_2[chan_ind + ind] = TBME;
		    
		    minus(tb, Space.qnums[qind], Space.qnums[sind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.pp1_map[chan][Space.hash2(qind, sind, tb.j)];
		    key2 = Chan.ph1_map[chan][Space.hash2(rind, pind, tb.j)];
		    chan_ind = Ints.Vhppp.V_2_4_index[chan];
		    ind = key1 * Chan.nph1[chan] + key2;
		    Ints.Vhppp.V_2_4[chan_ind + ind] = TBME;

		    // Vpphp -> rspq
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[pind]);
		    key1 = Chan.ppp_map[chan][Space.hash3(rind, sind, qind, Chan.qnums1[chan1].j)];
		    key2 = Chan.h_map[chan][pind];
		    chan_ind = Ints.Vpphp.V_3_3_index[chan];
		    ind = key1 * Chan.nh[chan] + key2;
		    Ints.Vpphp.V_3_3[chan_ind + ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[rind]);
		    key1 = Chan.hhp_map[chan][Space.hash3(pind, qind, sind, Chan.qnums1[chan1].j)];
		    key2 = Chan.h_map[chan][rind];
		    chan_ind = Ints.Vhhhp.V_3_3_index[chan];
		    ind = key1 * Chan.nh[chan] + key2;
		    Ints.Vhhhp.V_3_3[chan_ind + ind] = TBME;
	      
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[qind]);
		    key1 = Chan.h_map[chan][qind];
		    key2 = Chan.hph_map[chan][Space.hash3(rind, sind, pind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhhhp.V_3_2_index[chan];
		    ind = key1 * Chan.nhph[chan] + key2;
		    Ints.Vhhhp.V_3_2[chan_ind + ind] = TBME;

		    minus(tb, Space.qnums[pind], Space.qnums[rind]);
		    chan = Space.ind_2b_cross(Parameters.basis, tb);
		    key1 = Chan.hh1_map[chan][Space.hash2(pind, rind, tb.j)];
		    key2 = Chan.ph1_map[chan][Space.hash2(sind, qind, tb.j)];
		    chan_ind = Ints.Vhhhp.V_2_3_index[chan];
		    ind = key1 * Chan.nph1[chan] + key2;
		    Ints.Vhhhp.V_2_3[chan_ind + ind] = TBME;

		    // Vhphh -> rspq
		    chan = Space.ind_1b(Parameters.basis, Space.qnums[sind]);
		    key1 = Chan.p_map[chan][sind];
		    key2 = Chan.hhh_map[chan][Space.hash3(pind, qind, rind, Chan.qnums1[chan1].j)];
		    chan_ind = Ints.Vhphh.V_3_2_index[chan];
		    ind = key1 * Chan.nhhh[chan] + key2;
		    Ints.Vhphh.V_3_2[chan_ind + ind] = TBME;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
