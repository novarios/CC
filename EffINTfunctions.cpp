#include "EffINTfunctions.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "CCfunctions.hpp"
#include "INTfunctions.hpp"
#include "MATHfunctions.hpp"

void Eff_Interactions::Build(Channels &Chan)
{
  struct timespec time1, time2;
  double elapsed0 = 0.0;
  clock_gettime(CLOCK_MONOTONIC, &time1);

  this->Xhh.Build(Chan);
  this->Xpp.Build(Chan);
  this->Xhp.Build(Chan);
  this->Xhhhh.Build(Chan);
  this->Xpppp.Build(Chan);
  this->Xhphp.Build(Chan);
  this->Xhhhp.Build(Chan);
  this->Xhppp.Build(Chan);
  this->Xhphh.Build(Chan);
  this->Xpphp.Build(Chan);

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "Eff_INT Runtime = " << elapsed0 << " sec. " << std::endl;
}

void Eff_Interactions::Delete()
{
  this->Xhh.Delete();
  this->Xpp.Delete();
  this->Xhp.Delete();
  this->Xhhhh.Delete();
  this->Xpppp.Delete();
  this->Xhphp.Delete();
  this->Xhhhp.Delete();
  this->Xhppp.Delete();
  this->Xhphh.Delete();
  this->Xpphp.Delete();
}

void Eff_Interactions::Zero(bool flag)
{
  this->Xhh.X1_Zero(flag);
  this->Xhh.X_Zero(flag);
  this->Xpp.X_Zero(flag);
  this->Xhp.X_Zero(flag);
  this->Xhhhh.X_Zero(flag);
  this->Xpppp.X1_Zero(flag);
  this->Xhphp.X_Zero(flag);
  this->Xhphp.X1_Zero(flag);
  this->Xhphp.X2_Zero(flag);
  this->Xhphp.X3_Zero(flag);
  this->Xhhhp.X1_Zero(flag);
  this->Xhhhp.X_Zero(flag);
  this->Xhppp.X1_Zero(flag);
  this->Xhppp.X_Zero(flag);
  this->Xhphh.X1_Zero(flag);
  this->Xhphh.X_Zero(flag);
  this->Xpphp.X1_Zero(flag);
  this->Xpphp.X_Zero(flag);
}

void X_hh::Build(Channels &Chan)
{
  int chan3, ind, length;
  int i, j, nhh1;
  two_body hh1;
  two_body ij, ij_j;
  two_body *tb_ind;
  two_body *tb_j;

  length = Chan.nhh1[Chan.ind0];
  this->X1_2 = new double[length];
  this->X_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X1_2[ind] = 0.0;
    this->X_2[ind] = 0.0;
  }
  this->X_2_length = length;

  length = 0;
  this->X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nh[chan3];
  }
  this->X_3_length = length;

  tb_ind = new two_body[this->X_2_length];
  tb_j = new two_body[this->X_2_length];
  length = 0;
  nhh1 = Chan.nhh1[Chan.ind0];
  for(int hh1_ind = 0; hh1_ind < nhh1; ++hh1_ind){
    hh1 = Chan.hh1_state(Chan.ind0, hh1_ind);
    i = hh1.v1;
    j = hh1.v2;
    ij.v1 = i;
    ij.v2 = j;
    ij_j.v1 = SPB.qnums[i].j;
    ij_j.v2 = SPB.qnums[j].j;
    tb_ind[length] = ij;
    tb_j[length] = ij_j;
    ++length;
  }

  this->map3_ind = new int[X_3_length];
  this->map3_fac1 = new double[X_3_length];
  this->map3_fac2 = new double[X_3_length];
  this->map3_fac3 = new double[X_3_length];
  for(int ind = 0; ind < X_3_length; ++ind){
    this->map3_ind[ind] = 0;
    this->map3_fac1[ind] = 0.0;
    this->map3_fac2[ind] = 0.0;
    this->map3_fac3[ind] = 0.0;
  }

  for(int hh = 0; hh < this->X_2_length; ++hh){
    Map_2_to_3(this->map3_ind, this->map3_fac1, X_3_index, hh, Chan.h_map, Chan.h_map, Chan.nh, tb_ind, tb_j, 1);
    Map_2_to_3(this->map3_ind, this->map3_fac2, X_3_index, hh, Chan.h_map, Chan.h_map, Chan.nh, tb_ind, tb_j, 2);
    Map_2_to_3(this->map3_ind, this->map3_fac3, X_3_index, hh, Chan.h_map, Chan.h_map, Chan.nh, tb_ind, tb_j, 3);
  }
  delete[] tb_ind;
  delete[] tb_j;
}

void X_hh::X1_Zero(bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < this->X_2_length; ++ind){ this->X1_2[ind] = 0.0; }
  }
}

void X_hh::X_Zero(bool flag)
{
  int ind;
  if( !flag ){
    for(ind = 0; ind < this->X_2_length; ++ind){ this->X_2[ind] = 0.0; }
  }
}

void X_hh::Gather_X1_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->X1_2[index] += this->map3_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hh::Gather_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->X_2[index] += this->map3_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hh::Set_X1_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac2[ind + offset] * this->X1_2[index];
  }
}

void X_hh::Set_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac2[ind + offset] * this->X_2[index];
  }
}

void X_hh::Set_X1_3_OD(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac3[ind + offset] * this->X1_2[index];
  }
}

void X_hh::Set_X_3_OD(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac3[ind + offset] * this->X_2[index];
  }
}

void X_hh::Delete()
{
  delete[] this->X_2;
  delete[] this->X1_2;
  delete[] this->X_3_index;

  delete[] this->map3_ind;
  delete[] this->map3_fac1;
  delete[] this->map3_fac2;
  delete[] this->map3_fac3;
}

void X_pp::Build(Channels &Chan)
{
  int chan3, ind, length;
  int a, b, npp1;
  two_body pp1;
  two_body ab, ab_j;
  two_body *tb_ind;
  two_body *tb_j;

  length = Chan.npp1[Chan.ind0];
  this->X_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_2[ind] = 0.0;
  }
  this->X_2_length = length;

  length = 0;
  this->X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.np[chan3];
  }
  this->X_3_length = length;

  tb_ind = new two_body[this->X_2_length];
  tb_j = new two_body[this->X_2_length];
  length = 0;
  npp1 = Chan.npp1[Chan.ind0];
  for(int pp1_ind = 0; pp1_ind < npp1; ++pp1_ind){
    pp1 = Chan.pp1_state(Chan.ind0, pp1_ind);
    a = pp1.v1;
    b = pp1.v2;
    ab.v1 = a;
    ab.v2 = b;
    ab_j.v1 = SPB.qnums[a].j;
    ab_j.v2 = SPB.qnums[b].j;
    tb_ind[length] = ab;
    tb_j[length] = ab_j;
    ++length;
  }

  this->map3_ind = new int[X_3_length];
  this->map3_fac1 = new double[X_3_length];
  this->map3_fac2 = new double[X_3_length];
  this->map3_fac3 = new double[X_3_length];
  for(int ind = 0; ind < X_3_length; ++ind){
    this->map3_ind[ind] = 0;
    this->map3_fac1[ind] = 0.0;
    this->map3_fac2[ind] = 0.0;
    this->map3_fac3[ind] = 0.0;
  }

  for(int pp = 0; pp < this->X_2_length; ++pp){
    Map_2_to_3(this->map3_ind, this->map3_fac1, X_3_index, pp, Chan.p_map, Chan.p_map, Chan.np, tb_ind, tb_j, 1);
    Map_2_to_3(this->map3_ind, this->map3_fac2, X_3_index, pp, Chan.p_map, Chan.p_map, Chan.np, tb_ind, tb_j, 2);
    Map_2_to_3(this->map3_ind, this->map3_fac3, X_3_index, pp, Chan.p_map, Chan.p_map, Chan.np, tb_ind, tb_j, 3);
  }
  delete[] tb_ind;
  delete[] tb_j;
}

void X_pp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_length; ++ind){ this->X_2[ind] = 0.0; } }
}

void X_pp::Gather_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->X_2[index] += this->map3_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pp::Set_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac2[ind + offset] * this->X_2[index];
  }
}

void X_pp::Set_X_3_OD(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac3[ind + offset] * this->X_2[index];
  }
}

void X_pp::Delete()
{
  delete[] this->X_2;
  delete[] this->X_3_index;

  delete[] this->map3_ind;
  delete[] this->map3_fac1;
  delete[] this->map3_fac2;
  delete[] this->map3_fac3;
}

void X_hp::Build(Channels &Chan)
{
  int chan3, ind, length;
  int i, a, nhp1;
  two_body hp1;
  two_body ia, ia_j;
  two_body *tb_ind;
  two_body *tb_j;

  length = Chan.nhp1[Chan.ind0];
  this->X_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_2[ind] = 0.0;
  }
  this->X_2_length = length;

  length = 0;
  this->X_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.np[chan3];
  }
  this->X_3_length = length;

  tb_ind = new two_body[this->X_2_length];
  tb_j = new two_body[this->X_2_length];
  length = 0;
  nhp1 = Chan.nhp1[Chan.ind0];
  for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
    hp1 = Chan.hp1_state(Chan.ind0, hp1_ind);
    i = hp1.v1;
    a = hp1.v2;
    ia.v1 = i;
    ia.v2 = a;
    ia_j.v1 = SPB.qnums[i].j;
    ia_j.v2 = SPB.qnums[a].j;
    tb_ind[length] = ia;
    tb_j[length] = ia_j;
    ++length;
  }

  this->map3_ind = new int[X_3_length];
  this->map3_fac1 = new double[X_3_length];
  this->map3_fac2 = new double[X_3_length];
  for(int ind = 0; ind < X_3_length; ++ind){
    this->map3_ind[ind] = 0;
    this->map3_fac1[ind] = 0.0;
    this->map3_fac2[ind] = 0.0;
  }

  for(int hp = 0; hp < this->X_2_length; ++hp){
    Map_2_to_3(this->map3_ind, this->map3_fac1, X_3_index, hp, Chan.h_map, Chan.p_map, Chan.np, tb_ind, tb_j, 1);
    Map_2_to_3(this->map3_ind, this->map3_fac2, X_3_index, hp, Chan.h_map, Chan.p_map, Chan.np, tb_ind, tb_j, 2);
  }
  delete[] tb_ind;
  delete[] tb_j;
}

void X_hp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_length; ++ind){ this->X_2[ind] = 0.0; } }
}

void X_hp::Gather_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->X_2[index] += this->map3_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hp::Set_X_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    X3[ind] = this->map3_fac2[ind + offset] * this->X_2[index];
  }
}

void X_hp::Delete()
{
  delete[] this->X_2;
  delete[] this->X_3_index;

  delete[] this->map3_ind;
  delete[] this->map3_fac1;
  delete[] this->map3_fac2;
}

void X_hhhh::Build(Channels &Chan)
{
  int chan1, chan3, ind, length;
  int i, j, k, l, nhh;
  two_body hh1, hh2;
  four_body ijkl, ijkl_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhh[chan1];
  }
  this->X_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_1[ind] = 0.0;
  }
  this->X_1_length = length;

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  this->X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    this->X_3_4_index[chan3] = length;
    length += Chan.nhhh[chan3] * Chan.nh[chan3];
  }
  this->X_3_3_length = length;
  this->X_3_4_length = length;

  fb_ind = new four_body[this->X_1_length];
  fb_j = new four_body[this->X_1_length];
  J = new int[this->X_1_length];
  length = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
      hh1 = Chan.hh_state(chan1, hh1_ind);
      i = hh1.v1;
      j = hh1.v2;
      ijkl.v1 = i;
      ijkl.v2 = j;
      ijkl_j.v1 = SPB.qnums[i].j;
      ijkl_j.v2 = SPB.qnums[j].j;
      for(int hh2_ind = 0; hh2_ind < nhh; ++hh2_ind){
	hh2 = Chan.hh_state(chan1, hh2_ind);
	k = hh2.v1;
	l = hh2.v2;
	ijkl.v3 = k;
	ijkl.v4 = l;
	ijkl_j.v3 = SPB.qnums[k].j;
	ijkl_j.v4 = SPB.qnums[l].j;
	fb_ind[length] = ijkl;
	fb_j[length] = ijkl_j;
	J[length] = Chan.qnums1[chan1].j;
	++length;
      }
    }
  }

  this->map33_ind = new int[X_3_3_length];
  this->map33_fac1 = new double[X_3_3_length];
  for(int ind = 0; ind < X_3_3_length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac1[ind] = 0.0;
  }
  this->map34_ind = new int[X_3_4_length];
  this->map34_fac1 = new double[X_3_4_length];  
  for(int ind = 0; ind < X_3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac1[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hhhh = 0; hhhh < this->X_1_length; ++hhhh){
    Map_1_to_33(this->map33_ind,this->map33_fac1,this->X_3_3_index,hhhh, Chan.hhh_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
    Map_1_to_34(this->map34_ind,this->map34_fac1,this->X_3_4_index,hhhh, Chan.hhh_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
}

void X_hhhh::Gather_X_3_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    this->X_1[index] += this->map33_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hhhh::Gather_X_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    this->X_1[index] += this->map34_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hhhh::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_1_length; ++ind){ this->X_1[ind] = 0.0; } }
}

void X_hhhh::Delete()
{
  delete[] this->X_1;

  delete[] this->X_1_index;
  delete[] this->X_3_3_index;
  delete[] this->X_3_4_index;

  delete[] this->map33_ind;
  delete[] this->map33_fac1;

  delete[] this->map34_ind;
  delete[] this->map34_fac1;
}

void X_pppp::Build(Channels &Chan)
{
  int chan1, chan3, ind, length;
  int a, b, c, d, npp;
  two_body pp1, pp2;
  four_body abcd, abcd_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.npp[chan1];
  }
  this->X1_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X1_1[ind] = 0.0;
  }
  this->X_1_length = length;

  length = 0;
  this->X_3_1_index = new int[Chan.size3];
  this->X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_1_index[chan3] = length;
    this->X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nppp[chan3];
  }
  this->X_3_1_length = length;
  this->X_3_2_length = length;

  fb_ind = new four_body[this->X_1_length];
  fb_j = new four_body[this->X_1_length];
  J = new int[this->X_1_length];
  length = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    for(int pp1_ind = 0; pp1_ind < npp; ++pp1_ind){
      pp1 = Chan.pp_state(chan1, pp1_ind);
      a = pp1.v1;
      b = pp1.v2;
      abcd.v1 = a;
      abcd.v2 = b;
      abcd_j.v1 = SPB.qnums[a].j;
      abcd_j.v2 = SPB.qnums[b].j;
      for(int pp2_ind = 0; pp2_ind < npp; ++pp2_ind){
	pp2 = Chan.pp_state(chan1, pp2_ind);
	c = pp2.v1;
	d = pp2.v2;
	abcd.v3 = c;
	abcd.v4 = d;
	abcd_j.v3 = SPB.qnums[c].j;
	abcd_j.v4 = SPB.qnums[d].j;
	fb_ind[length] = abcd;
	fb_j[length] = abcd_j;
	J[length] = Chan.qnums1[chan1].j;
	++length;
      }
    }
  }

  this->map31_ind = new int[X_3_1_length];
  this->map31_fac1 = new double[X_3_1_length];
  for(int ind = 0; ind < X_3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac1[ind] = 0.0;
  }
  this->map32_ind = new int[X_3_2_length];
  this->map32_fac1 = new double[X_3_2_length];  
  for(int ind = 0; ind < X_3_2_length; ++ind){
    this->map32_ind[ind] = 0;
    this->map32_fac1[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int pppp = 0; pppp < this->X_1_length; ++pppp){
    Map_1_to_31(this->map31_ind,this->map31_fac1,this->X_3_1_index,pppp, Chan.p_map,Chan.ppp_map,Chan.nppp, fb_ind,fb_j,J, 1);
    Map_1_to_32(this->map32_ind,this->map32_fac1,this->X_3_2_index,pppp, Chan.p_map,Chan.ppp_map,Chan.nppp, fb_ind,fb_j,J, 1);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
}

void X_pppp::Gather_X1_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->X1_1[index] += this->map31_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pppp::Gather_X1_3_2(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    this->X1_1[index] += this->map32_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pppp::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_1_length; ++ind){ this->X1_1[ind] = 0.0; } }
}

void X_pppp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_1_length; ++ind){ this->X_1[ind] = 0.0; } }
}

void X_pppp::Delete()
{
  delete[] this->X1_1;

  delete[] this->X_1_index;
  delete[] this->X_3_1_index;
  delete[] this->X_3_2_index;

  delete[] this->map31_ind;
  delete[] this->map31_fac1;

  delete[] this->map32_ind;
  delete[] this->map32_fac1;
}

void X_hphp::Build(Channels &Chan)
{
  int chan2, chan3, ind, length;
  int i, a, j, b, nhp1;
  two_body hp1, hp2;
  four_body iajb, iajb_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J31_map, *J32_map, *J33_map, *J34_map;

  length = 0;
  this->X_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_1_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nhp1[chan2];
  }
  this->X_2_1 = new double[length];
  this->X1_2_1 = new double[length];
  this->X2_2_1 = new double[length];
  this->X3_2_1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_2_1[ind] = 0.0;
    this->X1_2_1[ind] = 0.0;
    this->X2_2_1[ind] = 0.0;
    this->X3_2_1[ind] = 0.0;
  }
  this->X_2_1_length = length;

  length = 0;
  this->X_2_2_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_2_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.nph1[chan2];
  }
  this->X_2_2_length = length;

  length = 0;
  this->X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhpp[chan3];
  }
  this->X_3_1_length = length;
  this->map31_num = new int[length];
  this->map31_index = new int[length];
  J31_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map31_num[ind] = 0; }

  length = 0;
  this->X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhph[chan3];
  }
  this->X_3_2_length = length;
  this->map32_num = new int[length];
  this->map32_index = new int[length];
  J32_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map32_num[ind] = 0; }

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    length += Chan.nhpp[chan3] * Chan.nh[chan3];
  }
  this->X_3_3_length = length;
  this->map33_num = new int[length];
  this->map33_index = new int[length];
  J33_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map33_num[ind] = 0; }

  length = 0;
  this->X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_4_index[chan3] = length;
    length += Chan.nhph[chan3] * Chan.np[chan3];
  }
  this->X_3_4_length = length;
  this->map34_num = new int[length];
  this->map34_index = new int[length];
  J34_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map34_num[ind] = 0; }

  fb_ind = new four_body[this->X_2_1_length];
  fb_j = new four_body[this->X_2_1_length];
  J = new int[this->X_2_1_length];
  length = 0;
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    for(int hp1_ind = 0; hp1_ind < nhp1; ++hp1_ind){
      hp1 = Chan.hp1_state(chan2, hp1_ind);
      i = hp1.v1;
      b = hp1.v2;
      iajb.v1 = i;
      iajb.v4 = b;
      iajb_j.v1 = SPB.qnums[i].j;
      iajb_j.v4 = SPB.qnums[b].j;
      for(int hp2_ind = 0; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	j = hp2.v1;
	a = hp2.v2;
	iajb.v3 = j;
	iajb.v2 = a;
	iajb_j.v3 = SPB.qnums[j].j;
	iajb_j.v2 = SPB.qnums[a].j;
	fb_ind[length] = iajb;
	fb_j[length] = iajb_j;
	J[length] = Chan.qnums2[chan2].j;
	Map_Count_2_to_31(this->map31_num, this->X_3_1_index, Chan.h_map, Chan.hpp_map, Chan.nhpp, iajb, J[length], J31_map);
	Map_Count_2_to_32(this->map32_num, this->X_3_2_index, Chan.p_map, Chan.hph_map, Chan.nhph, iajb, J[length], J32_map);
	Map_Count_2_to_33(this->map33_num, this->X_3_3_index, Chan.hpp_map, Chan.h_map, Chan.nh, iajb, J[length], J33_map);
	Map_Count_2_to_34(this->map34_num, this->X_3_4_index, Chan.hph_map, Chan.p_map, Chan.np, iajb, J[length], J34_map);
	++length;
      }
    }
  }

  this->map22_ind = new int[X_2_2_length];
  this->map22_fac2 = new double[X_2_2_length];
  for(int ind = 0; ind < X_2_2_length; ++ind){
    this->map22_ind[ind] = 0;
    this->map22_fac2[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < X_3_1_length; ++ind){
    this->map31_index[ind] = length;
    length += this->map31_num[ind];
  }
  this->map31_ind = new int[length];
  this->map31_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac2[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < X_3_2_length; ++ind){
    this->map32_index[ind] = length;
    length += this->map32_num[ind];
  }
  this->map32_ind = new int[length];
  this->map32_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map32_ind[ind] = 0;
    this->map32_fac1[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < X_3_3_length; ++ind){
    this->map33_index[ind] = length;
    length += this->map33_num[ind];
  }
  this->map33_ind = new int[length];
  this->map33_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac1[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < X_3_4_length; ++ind){
    this->map34_index[ind] = length;
    length += this->map34_num[ind];
  }
  this->map34_ind = new int[length];
  this->map34_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hphp = 0; hphp < this->X_2_1_length; ++hphp){
    Map_21_to_22(this->map22_ind,this->map22_fac2,this->X_2_2_index,hphp, Chan.ph1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, 2);
    Map_21_to_31(this->map31_ind,this->map31_fac2,this->map31_index, this->X_3_1_index,hphp, Chan.h_map,Chan.hpp_map,Chan.nhpp, fb_ind,fb_j,J, J31_map,2);
    Map_21_to_32(this->map32_ind,this->map32_fac1,this->map32_index, this->X_3_2_index,hphp, Chan.p_map,Chan.hph_map,Chan.nhph, fb_ind,fb_j,J, J32_map,1);
    Map_21_to_33(this->map33_ind,this->map33_fac1,this->map33_index, this->X_3_3_index,hphp, Chan.hpp_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, J33_map,1);
    Map_21_to_34(this->map34_ind,this->map34_fac2,this->map34_index, this->X_3_4_index,hphp, Chan.hph_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, J34_map,2);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J31_map;
  delete[] J32_map;
  delete[] J33_map;
  delete[] J34_map;
}

void X_hphp::Set_X_2_2(double *X2, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map22_ind[ind + offset];
    X2[ind] = this->map22_fac2[ind + offset] * this->X_2_1[index];
  }
}

void X_hphp::Gather_X1_3_2(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map32_num[ind + offset]; ++n){
      index = this->map32_index[ind + offset] + n;
      index1 = this->map32_ind[index];
      this->X1_2_1[index1] += this->map32_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Gather_X1_3_3(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map33_num[ind + offset]; ++n){
      index = this->map33_index[ind + offset] + n;
      index1 = this->map33_ind[index];
      this->X1_2_1[index1] += this->map33_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Set_X1_3_1(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X3[ind] = 0.0;
    for(int n = 0; n < this->map31_num[ind + offset]; ++n){
      index = this->map31_index[ind + offset] + n;
      index1 = this->map31_ind[index];
      X3[ind] += this->map31_fac2[index] * this->X1_2_1[index1];
    }
  }
}

void X_hphp::Gather_X2_3_2(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map32_num[ind + offset]; ++n){
      index = this->map32_index[ind + offset] + n;
      index1 = this->map32_ind[index];
      this->X2_2_1[index1] += this->map32_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Gather_X2_3_3(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map33_num[ind + offset]; ++n){
      index = this->map33_index[ind + offset] + n;
      index1 = this->map33_ind[index];
      this->X2_2_1[index1] += this->map33_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Set_X2_3_1(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X3[ind] = 0.0;
    for(int n = 0; n < this->map31_num[ind + offset]; ++n){
      index = this->map31_index[ind + offset] + n;
      index1 = this->map31_ind[index];
      X3[ind] += this->map31_fac2[index] * this->X2_2_1[index1];
    }
  }
}

void X_hphp::Gather_X3_3_2(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map32_num[ind + offset]; ++n){
      index = this->map32_index[ind + offset] + n;
      index1 = this->map32_ind[index];
      this->X3_2_1[index1] += this->map32_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Gather_X3_3_3(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map33_num[ind + offset]; ++n){
      index = this->map33_index[ind + offset] + n;
      index1 = this->map33_ind[index];
      this->X3_2_1[index1] += this->map33_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Set_X3_3_4(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X3[ind] = 0.0;
    for(int n = 0; n < this->map34_num[ind + offset]; ++n){
      index = this->map34_index[ind + offset] + n;
      index1 = this->map34_ind[index];
      X3[ind] += this->map34_fac2[index] * this->X3_2_1[index1];
    }
  }
}

void X_hphp::Gather_X_3_2(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map32_num[ind + offset]; ++n){
      index = this->map32_index[ind + offset] + n;
      index1 = this->map32_ind[index];
      this->X_2_1[index1] += this->map32_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::Gather_X_3_3(double *X3, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map33_num[ind + offset]; ++n){
      index = this->map33_index[ind + offset] + n;
      index1 = this->map33_ind[index];
      this->X_2_1[index1] += this->map33_fac1[index] * X3[ind];
    }
    X3[ind] = 0.0;
  }
}

void X_hphp::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_1_length; ++ind){ this->X1_2_1[ind] = 0.0; } }
}

void X_hphp::X2_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_1_length; ++ind){ this->X2_2_1[ind] = 0.0; } }
}

void X_hphp::X3_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_1_length; ++ind){ this->X3_2_1[ind] = 0.0; } }
}

void X_hphp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_2_1_length; ++ind){ this->X_2_1[ind] = 0.0; } }
}

void X_hphp::Delete()
{
  delete[] this->X1_2_1;
  delete[] this->X2_2_1;
  delete[] this->X3_2_1;
  delete[] this->X_2_1;

  delete[] this->X_2_1_index;
  delete[] this->X_3_1_index;
  delete[] this->X_3_2_index;
  delete[] this->X_3_3_index;
  delete[] this->X_3_4_index;

  delete[] this->map31_index;
  delete[] this->map31_num;
  delete[] this->map31_ind;
  delete[] this->map31_fac2;

  delete[] this->map32_index;
  delete[] this->map32_num;
  delete[] this->map32_ind;
  delete[] this->map32_fac1;

  delete[] this->map33_index;
  delete[] this->map33_num;
  delete[] this->map33_ind;
  delete[] this->map33_fac1;

  delete[] this->map34_index;
  delete[] this->map34_num;
  delete[] this->map34_ind;
  delete[] this->map34_fac2;
}

void X_hppp::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int i, a, b, c, np, npph;
  three_body pph;
  one_body p;
  four_body iabc, iabc_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J23_map;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.nhp[chan1] * Chan.npp[chan1];
  }
  this->X_1_length = length;

  length = 0;
  this->X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_3_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.npp1[chan2];
  }
  this->X_2_3_length = length;
  this->map23_num = new int[length];
  this->map23_index = new int[length];
  J23_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map23_num[ind] = 0; }

  length = 0;
  this->X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nppp[chan3];
  }
  this->X_3_1_length = length;

  length = 0;
  this->X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.npph[chan3];
  }
  this->X_3_2 = new double[length];
  this->X1_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_3_2[ind] = 0.0;
    this->X1_3_2[ind] = 0.0;
  }
  this->X_3_2_length = length;

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    length += Chan.nhpp[chan3] * Chan.np[chan3];
  }
  this->X_3_3_length = length;

  fb_ind = new four_body[this->X_3_2_length];
  fb_j = new four_body[this->X_3_2_length];
  J = new int[this->X_3_2_length];
  length = 0;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      iabc.v2 = a;
      iabc_j.v2 = SPB.qnums[a].j;
      for(int pph_ind = 0; pph_ind < npph; ++pph_ind){
	pph = Chan.pph_state(chan3, pph_ind);
	b = pph.v1;
	c = pph.v2;
	i = pph.v3;
	iabc.v1 = i;
	iabc.v3 = b;
	iabc.v4 = c;
	iabc_j.v1 = SPB.qnums[i].j;
	iabc_j.v3 = SPB.qnums[b].j;
	iabc_j.v4 = SPB.qnums[c].j;
	fb_ind[length] = iabc;
	fb_j[length] = iabc_j;
	J[length] = Chan.pph_j[Chan.pph_index[chan3] + pph_ind].j;
	Map_Count_1_to_23(this->map23_num, this->X_2_3_index, Chan.hp1_map, Chan.pp1_map, Chan.npp1, iabc, J[length], J23_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < X_2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac2[ind] = 0.0;
  }
  this->map1_ind = new int[X_1_length];
  this->map1_fac2 = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    this->map1_ind[ind] = 0;
    this->map1_fac2[ind] = 0.0;
  }
  this->map31_ind = new int[X_3_1_length];
  this->map31_fac2 = new double[X_3_1_length];
  for(int ind = 0; ind < X_3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac2[ind] = 0.0;
  }
  this->map33_ind = new int[X_3_3_length];
  this->map33_fac2 = new double[X_3_3_length];  
  for(int ind = 0; ind < X_3_3_length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hppp = 0; hppp < this->X_3_2_length; ++hppp){
    Map_32_to_1(this->map1_ind,this->map1_fac2,this->X_1_index,hppp, Chan.hp_map,Chan.pp_map,Chan.npp, fb_ind,fb_j,J, 2);
    Map_32_to_23(this->map23_ind,this->map23_fac2,this->map23_index, this->X_2_3_index,hppp, Chan.hp1_map,Chan.pp1_map,Chan.npp1, fb_ind,fb_j,J, J23_map,2);
    Map_32_to_31(this->map31_ind,this->map31_fac2,this->X_3_1_index,hppp, Chan.h_map,Chan.ppp_map,Chan.nppp, fb_ind,fb_j,J, 2);
    Map_32_to_33(this->map33_ind,this->map33_fac2,this->X_3_3_index,hppp, Chan.hpp_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 2);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J23_map;
}

void X_hppp::Set_X1_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    X3[ind] = this->map31_fac2[ind + offset] * this->X1_3_2[index];
  }
}

void X_hppp::Set_X1_3_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    X3[ind] = this->map33_fac2[ind + offset] * this->X1_3_2[index];
  }
}

void X_hppp::Set_X_3_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    X3[ind] = this->map33_fac2[ind + offset] * this->X_3_2[index];
  }
}

void X_hppp::Set_X_2_3(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X2[ind] = 0.0;
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      X2[ind] += this->map23_fac2[index] * this->X_3_2[index1];
    }
  }
}

void X_hppp::Set_X_1(double *X1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map1_ind[ind + offset];
    X1[ind] = this->map1_fac2[ind + offset] * this->X_3_2[index];
  }
}

void X_hppp::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_2_length; ++ind){ this->X1_3_2[ind] = 0.0; } }
}

void X_hppp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_2_length; ++ind){ this->X_3_2[ind] = 0.0; } }
}

void X_hppp::Delete()
{
  delete[] this->X1_3_2;
  delete[] this->X_3_2;

  delete[] this->X_1_index;
  delete[] this->X_2_3_index;
  delete[] this->X_3_1_index;
  delete[] this->X_3_2_index;
  delete[] this->X_3_3_index;

  delete[] this->map1_ind;
  delete[] this->map1_fac2;

  delete[] this->map23_index;
  delete[] this->map23_num;
  delete[] this->map23_ind;
  delete[] this->map23_fac2;

  delete[] this->map31_ind;
  delete[] this->map31_fac2;

  delete[] this->map33_ind;
  delete[] this->map33_fac2;
}

void X_hhhp::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int i, j, k, a, nh, nhhp;
  three_body hhp;
  one_body h;
  four_body ijka, ijka_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J21_map, *J23_map;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhp[chan1];
  }
  this->X_1_length = length;

  length = 0;
  this->X_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_1_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nhh1[chan2];
  }
  this->X_2_1_length = length;
  this->map21_num = new int[length];
  this->map21_index = new int[length];
  J21_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map21_num[ind] = 0; }

  length = 0;
  this->X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nph1[chan2];
  }
  this->X_2_3_length = length;
  this->map23_num = new int[length];
  this->map23_index = new int[length];
  J23_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){ this->map23_num[ind] = 0; }

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.nh[chan3];
  }
  this->X_3_3 = new double[length];
  this->X1_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_3_3[ind] = 0.0;
    this->X1_3_3[ind] = 0.0;
  }
  this->X_3_3_length = length;

  length = 0;
  this->X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_4_index[chan3] = length;
    length += Chan.nhhh[chan3] * Chan.np[chan3];
  }
  this->X_3_4_length = length;

  fb_ind = new four_body[this->X_3_3_length];
  fb_j = new four_body[this->X_3_3_length];
  J = new int[this->X_3_3_length];
  length = 0;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    nhhp = Chan.nhhp[chan3];
    nh = Chan.nh[chan3];
    for(int hhp_ind = 0; hhp_ind < nhhp; ++hhp_ind){
      hhp = Chan.hhp_state(chan3, hhp_ind);
      i = hhp.v1;
      j = hhp.v2;
      a = hhp.v3;
      ijka.v1 = i;
      ijka.v2 = j;
      ijka.v4 = a;
      ijka_j.v1 = SPB.qnums[i].j;
      ijka_j.v2 = SPB.qnums[j].j;
      ijka_j.v4 = SPB.qnums[a].j;
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	k = h.v1;
	ijka.v3 = k;
	ijka_j.v3 = SPB.qnums[k].j;
	fb_ind[length] = ijka;
	fb_j[length] = ijka_j;
	J[length] = Chan.hhp_j[Chan.hhp_index[chan3] + hhp_ind].j;
	Map_Count_1_to_21(this->map21_num, this->X_2_1_index, Chan.hp1_map, Chan.hh1_map, Chan.nhh1, ijka, J[length], J21_map);
	Map_Count_1_to_23(this->map23_num, this->X_2_3_index, Chan.hh1_map, Chan.ph1_map, Chan.nph1, ijka, J[length], J23_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < X_2_1_length; ++ind){
    this->map21_index[ind] = length;
    length += this->map21_num[ind];
  }
  this->map21_ind = new int[length];
  this->map21_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map21_ind[ind] = 0;
    this->map21_fac2[ind] = 0.0;
  }
  length = 0;
  for(int ind = 0; ind < X_2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac2 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac2[ind] = 0.0;
  }
  this->map1_ind = new int[X_1_length];
  this->map1_fac2 = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    this->map1_ind[ind] = 0;
    this->map1_fac2[ind] = 0.0;
  }
  this->map34_ind = new int[X_3_4_length];
  this->map34_fac2 = new double[X_3_4_length];
  for(int ind = 0; ind < X_3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hhhp = 0; hhhp < this->X_3_3_length; ++hhhp){
    Map_33_to_1(this->map1_ind,this->map1_fac2, this->X_1_index,hhhp, Chan.hh_map,Chan.hp_map,Chan.nhp, fb_ind,fb_j,J, 2);
    Map_33_to_21(this->map21_ind,this->map21_fac2, this->map21_index,this->X_2_1_index,hhhp, Chan.hp1_map,Chan.hh1_map,Chan.nhh1, fb_ind,fb_j,J,J21_map,2);
    Map_33_to_23(this->map23_ind,this->map23_fac2, this->map23_index,this->X_2_3_index,hhhp, Chan.hh1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J,J23_map,2);
    Map_33_to_34(this->map34_ind,this->map34_fac2, this->X_3_4_index,hhhp, Chan.hhh_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 2);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J21_map;
  delete[] J23_map;
}

void X_hhhp::Set_X1_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    X3[ind] = this->map34_fac2[ind + offset] * this->X1_3_3[index];
  }
}

void X_hhhp::Set_X_2_1(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X2[ind] = 0.0;
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      X2[ind] += this->map21_fac2[index] * this->X_3_3[index1];
    }
  }
}

void X_hhhp::Set_X_2_3(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    X2[ind] = 0.0;
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      X2[ind] += this->map23_fac2[index] * this->X_3_3[index1];
    }
  }
}

void X_hhhp::Set_X_1(double *X1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map1_ind[ind + offset];
    X1[ind] = this->map1_fac2[ind + offset] * this->X_3_3[index];
  }
}

void X_hhhp::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_3_length; ++ind){ this->X1_3_3[ind] = 0.0; } }
}

void X_hhhp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_3_length; ++ind){ this->X_3_3[ind] = 0.0; } }
}

void X_hhhp::Delete()
{
  delete[] this->X1_3_3;
  delete[] this->X_3_3;

  delete[] this->X_1_index;
  delete[] this->X_2_1_index;
  delete[] this->X_2_3_index;
  delete[] this->X_3_3_index;
  delete[] this->X_3_4_index;

  delete[] this->map1_ind;
  delete[] this->map1_fac2;

  delete[] this->map21_index;
  delete[] this->map21_num;
  delete[] this->map21_ind;
  delete[] this->map21_fac2;

  delete[] this->map23_index;
  delete[] this->map23_num;
  delete[] this->map23_ind;
  delete[] this->map23_fac2;

  delete[] this->map34_ind;
  delete[] this->map34_fac2;
}

void X_hphh::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int i, a, j, k, np, nhhh;
  three_body hhh;
  one_body p;
  four_body iajk, iajk_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J21_map, *J23_map;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.nhp[chan1] * Chan.npp[chan1];
  }
  this->X_1_length = length;

  length = 0;
  this->X_2_1_index = new int[Chan.size2];
  this->X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_1_index[chan2] = length;
    this->X_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nhp1[chan2];
  }
  this->X_2_1_length = length;
  this->X_2_3_length = length;
  this->map21_num = new int[length];
  this->map23_num = new int[length];
  this->map21_index = new int[length];
  this->map23_index = new int[length];
  J21_map = new std::unordered_map<int,int>[length];
  J23_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){
    this->map21_num[ind] = 0;
    this->map23_num[ind] = 0;
  }

  length = 0;
  this->X_3_1_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_1_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhhp[chan3];
  }
  this->X_3_1_length = length;

  length = 0;
  this->X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhh[chan3];
  }
  this->X_3_2 = new double[length];
  this->X1_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_3_2[ind] = 0.0;
    this->X1_3_2[ind] = 0.0;
  }
  this->X_3_2_length = length;

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  this->X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    this->X_3_4_index[chan3] = length;
    length += Chan.nhph[chan3] * Chan.nh[chan3];
  }
  this->X_3_3_length = length;
  this->X_3_4_length = length;

  fb_ind = new four_body[this->X_3_2_length];
  fb_j = new four_body[this->X_3_2_length];
  J = new int[this->X_3_2_length];
  length = 0;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    for(int p_ind = 0; p_ind < np; ++p_ind){
      p = Chan.p_state(chan3, p_ind);
      a = p.v1;
      iajk.v2 = a;
      iajk_j.v2 = SPB.qnums[a].j;
      for(int hhh_ind = 0; hhh_ind < nhhh; ++hhh_ind){
	hhh = Chan.hhh_state(chan3, hhh_ind);
	j = hhh.v1;
	k = hhh.v2;
	i = hhh.v3;
	iajk.v1 = i;
	iajk.v3 = j;
	iajk.v4 = k;
	iajk_j.v1 = SPB.qnums[i].j;
	iajk_j.v3 = SPB.qnums[j].j;
	iajk_j.v4 = SPB.qnums[k].j;
	fb_ind[length] = iajk;
	fb_j[length] = iajk_j;
	J[length] = Chan.hhh_j[Chan.hhh_index[chan3] + hhh_ind].j;
	Map_Count_1_to_23(this->map21_num, this->X_2_1_index, Chan.hh1_map, Chan.hp1_map, Chan.nhp1, iajk, J[length], J21_map);
	Map_Count_1_to_23(this->map23_num, this->X_2_3_index, Chan.hh1_map, Chan.hp1_map, Chan.nhp1, iajk, J[length], J23_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < X_2_1_length; ++ind){
    this->map21_index[ind] = length;
    length += this->map21_num[ind];
  }
  this->map21_ind = new int[length];
  this->map21_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map21_ind[ind] = 0;
    this->map21_fac1[ind] = 0.0;
  }
  length = 0;
  for(int ind = 0; ind < X_2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac1[ind] = 0.0;
  }
  this->map1_ind = new int[X_1_length];
  this->map1_fac1 = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    this->map1_ind[ind] = 0;
    this->map1_fac1[ind] = 0.0;
  }
  this->map31_ind = new int[X_3_1_length];
  this->map31_fac1 = new double[X_3_1_length];
  this->map31_fac2 = new double[X_3_1_length];
  for(int ind = 0; ind < X_3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac1[ind] = 0.0;
    this->map31_fac2[ind] = 0.0;
  }
  this->map33_ind = new int[X_3_3_length];
  this->map33_fac1 = new double[X_3_3_length];  
  for(int ind = 0; ind < X_3_3_length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac1[ind] = 0.0;
  }
  this->map34_ind = new int[X_3_4_length];
  this->map34_fac1 = new double[X_3_4_length];  
  for(int ind = 0; ind < X_3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac1[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hphh = 0; hphh < this->X_3_2_length; ++hphh){
    Map_32_to_1(this->map1_ind,this->map1_fac1,this->X_1_index,hphh, Chan.hp_map,Chan.hh_map,Chan.nhh, fb_ind,fb_j,J, 1);
    Map_32_to_21(this->map21_ind,this->map21_fac1,this->map21_index, this->X_2_1_index,hphh, Chan.hh1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J21_map,1);
    Map_32_to_23(this->map23_ind,this->map23_fac1,this->map23_index, this->X_2_3_index,hphh, Chan.hh1_map,Chan.hp1_map,Chan.nhp1, fb_ind,fb_j,J, J23_map,1);
    Map_32_to_31(this->map31_ind,this->map31_fac1,this->X_3_1_index,hphh, Chan.h_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 1);
    Map_32_to_31(this->map31_ind,this->map31_fac2,this->X_3_1_index,hphh, Chan.h_map,Chan.hhp_map,Chan.nhhp, fb_ind,fb_j,J, 2);
    Map_32_to_33(this->map33_ind,this->map33_fac1,this->X_3_3_index,hphh, Chan.hph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
    Map_32_to_34(this->map34_ind,this->map34_fac1,this->X_3_4_index,hphh, Chan.hph_map,Chan.h_map,Chan.nh, fb_ind,fb_j,J, 1);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J21_map;
  delete[] J23_map;
}

void X_hphh::Set_X1_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    X3[ind] = this->map31_fac2[ind + offset] * this->X1_3_2[index];
  }
}

void X_hphh::Set_X_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    X3[ind] = this->map31_fac2[ind + offset] * this->X_3_2[index];
  }
}

void X_hphh::Gather_X_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->X_3_2[index] += this->map31_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hphh::Gather_X_3_3(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    this->X_3_2[index] += this->map33_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hphh::Gather_X_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    this->X_3_2[index] += this->map34_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_hphh::Gather_X_2_1(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      this->X_3_2[index1] += this->map21_fac1[index] * X2[ind];
    }
    X2[ind] = 0.0;
  }
}

void X_hphh::Gather_X_2_3(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      this->X_3_2[index1] += this->map23_fac1[index] * X2[ind];
    }
    X2[ind] = 0.0;
  }
}

void X_hphh::Gather_X_1(double *X1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map1_ind[ind + offset];
    this->X_3_2[index] += this->map1_fac1[ind + offset] * X1[ind];
    X1[ind] = 0.0;
  }
}

void X_hphh::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_2_length; ++ind){ this->X1_3_2[ind] = 0.0; } }
}

void X_hphh::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_2_length; ++ind){ this->X_3_2[ind] = 0.0; } }
}

void X_hphh::Delete()
{
  delete[] this->X1_3_2;
  delete[] this->X_3_2;

  delete[] this->X_1_index;
  delete[] this->X_2_1_index;
  delete[] this->X_2_3_index;
  delete[] this->X_3_1_index;
  delete[] this->X_3_2_index;
  delete[] this->X_3_3_index;
  delete[] this->X_3_4_index;

  delete[] this->map1_ind;
  delete[] this->map1_fac1;

  delete[] this->map21_index;
  delete[] this->map21_num;
  delete[] this->map21_ind;
  delete[] this->map21_fac1;

  delete[] this->map23_index;
  delete[] this->map23_num;
  delete[] this->map23_ind;
  delete[] this->map23_fac1;

  delete[] this->map31_ind;
  delete[] this->map31_fac1;
  delete[] this->map31_fac2;

  delete[] this->map33_ind;
  delete[] this->map33_fac1;

  delete[] this->map34_ind;
  delete[] this->map34_fac1;
}

void X_pphp::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int a, b, i, c, nh, nppp;
  three_body ppp;
  one_body h;
  four_body abic, abic_j;
  four_body *fb_ind;
  four_body *fb_j;
  int *J;
  std::unordered_map<int,int> *J22_map, *J23_map;

  length = 0;
  this->X_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->X_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhp[chan1];
  }
  this->X_1_length = length;

  length = 0;
  this->X_2_2_index = new int[Chan.size2];
  this->X_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->X_2_2_index[chan2] = length;
    this->X_2_3_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.npp1[chan2];
  }
  this->X_2_2_length = length;
  this->X_2_3_length = length;
  this->map22_num = new int[length];
  this->map23_num = new int[length];
  this->map22_index = new int[length];
  this->map23_index = new int[length];
  J22_map = new std::unordered_map<int,int>[length];
  J23_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){
    this->map22_num[ind] = 0;
    this->map23_num[ind] = 0;
  }

  length = 0;
  this->X_3_1_index = new int[Chan.size3];
  this->X_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_1_index[chan3] = length;
    this->X_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhpp[chan3];
  }
  this->X_3_1_length = length;
  this->X_3_2_length = length;

  length = 0;
  this->X_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_3_index[chan3] = length;
    length += Chan.nppp[chan3] * Chan.nh[chan3];
  }
  this->X_3_3 = new double[length];
  this->X1_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->X_3_3[ind] = 0.0;
    this->X1_3_3[ind] = 0.0;
  }
  this->X_3_3_length = length;

  length = 0;
  this->X_3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->X_3_4_index[chan3] = length;
    length += Chan.npph[chan3] * Chan.np[chan3];
  }
  this->X_3_4_length = length;

  /*for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    std::cout << "(" << chan3 << "," << this->X_3_1_index[chan3] << "," << this->X_3_2_index[chan3] << "," << this->X_3_3_index[chan3] << "," << this->X_3_4_index[chan3] << ")" << std::endl;
    }*/

  fb_ind = new four_body[this->X_3_3_length];
  fb_j = new four_body[this->X_3_3_length];
  J = new int[this->X_3_3_length];
  length = 0;
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    nppp = Chan.nppp[chan3];
    nh = Chan.nh[chan3];
    for(int ppp_ind = 0; ppp_ind < nppp; ++ppp_ind){
      ppp = Chan.ppp_state(chan3, ppp_ind);
      a = ppp.v1;
      b = ppp.v2;
      c = ppp.v3;
      abic.v1 = a;
      abic.v2 = b;
      abic.v4 = c;
      abic_j.v1 = SPB.qnums[a].j;
      abic_j.v2 = SPB.qnums[b].j;
      abic_j.v4 = SPB.qnums[c].j;
      for(int h_ind = 0; h_ind < nh; ++h_ind){
	h = Chan.h_state(chan3, h_ind);
	i = h.v1;
	abic.v3 = i;
	abic_j.v3 = SPB.qnums[i].j;
	fb_ind[length] = abic;
	fb_j[length] = abic_j;
	J[length] = Chan.ppp_j[Chan.ppp_index[chan3] + ppp_ind].j;
	Map_Count_1_to_22(this->map22_num, this->X_2_2_index, Chan.ph1_map, Chan.pp1_map, Chan.npp1, abic, J[length], J22_map);
	Map_Count_1_to_23(this->map23_num, this->X_2_3_index, Chan.ph1_map, Chan.pp1_map, Chan.npp1, abic, J[length], J23_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < X_2_2_length; ++ind){
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
  for(int ind = 0; ind < X_2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac1[ind] = 0.0;
  }
  this->map1_ind = new int[X_1_length];
  this->map1_fac1 = new double[X_1_length];
  for(int ind = 0; ind < X_1_length; ++ind){
    this->map1_ind[ind] = 0;
    this->map1_fac1[ind] = 0.0;
  }
  this->map31_ind = new int[X_3_1_length];
  this->map31_fac1 = new double[X_3_1_length];
  for(int ind = 0; ind < X_3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac1[ind] = 0.0;
  }
  this->map32_ind = new int[X_3_2_length];
  this->map32_fac1 = new double[X_3_2_length];  
  for(int ind = 0; ind < X_3_2_length; ++ind){
    this->map32_ind[ind] = 0;
    this->map32_fac1[ind] = 0.0;
  }
  this->map34_ind = new int[X_3_4_length];
  this->map34_fac1 = new double[X_3_4_length];
  this->map34_fac2 = new double[X_3_4_length];
  for(int ind = 0; ind < X_3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac1[ind] = 0.0;
    this->map34_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int pphp = 0; pphp < this->X_3_3_length; ++pphp){
    Map_33_to_1(this->map1_ind,this->map1_fac1,this->X_1_index,pphp, Chan.pp_map,Chan.hp_map,Chan.nhp, fb_ind,fb_j,J, 1);
    Map_33_to_22(this->map22_ind,this->map22_fac1,this->map22_index, this->X_2_2_index,pphp, Chan.ph1_map,Chan.pp1_map,Chan.npp1, fb_ind,fb_j,J,J22_map,1);
    Map_33_to_23(this->map23_ind,this->map23_fac1,this->map23_index, this->X_2_3_index,pphp, Chan.ph1_map,Chan.pp1_map,Chan.npp1, fb_ind,fb_j,J,J23_map,1);
    Map_33_to_31(this->map31_ind,this->map31_fac1,this->X_3_1_index,pphp, Chan.p_map,Chan.hpp_map,Chan.nhpp, fb_ind,fb_j,J, 1);
    Map_33_to_32(this->map32_ind,this->map32_fac1,this->X_3_2_index,pphp, Chan.p_map,Chan.hpp_map,Chan.nhpp, fb_ind,fb_j,J, 1);
    Map_33_to_34(this->map34_ind,this->map34_fac1,this->X_3_4_index,pphp, Chan.pph_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 1);
    Map_33_to_34(this->map34_ind,this->map34_fac2,this->X_3_4_index,pphp, Chan.pph_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 2);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J22_map;
  delete[] J23_map;
}

void X_pphp::Set_X1_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    X3[ind] = this->map34_fac2[ind + offset] * this->X1_3_3[index];
    //if( ind + offset == 3489 ){ std::cout << "3_4_set: " << index << ", " << this->map34_fac2[ind + offset] << ", " << this->X1_3_3[index] << std::endl; }
  }
}

void X_pphp::Set_X_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    X3[ind] = this->map34_fac2[ind + offset] * this->X_3_3[index];
  }
}

void X_pphp::Gather_X1_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->X1_3_3[index] += this->map31_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
    //if( ind + offset == 31 ){ std::cout << "3_1_gather: " << index << ", " << this->map31_fac1[ind + offset] << ", " << X3[ind] << std::endl; }
  }
}

void X_pphp::Gather_X1_3_2(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    this->X1_3_3[index] += this->map32_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
    //if( ind + offset == 1356 ){ std::cout << "3_2_gather: " << index << ", " << this->map32_fac1[ind + offset] << ", " << X3[ind] << std::endl; }
  }
}

void X_pphp::Gather_X_3_1(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->X_3_3[index] += this->map31_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pphp::Gather_X_3_2(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    this->X_3_3[index] += this->map32_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pphp::Gather_X_3_4(double *X3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    this->X_3_3[index] += this->map34_fac1[ind + offset] * X3[ind];
    X3[ind] = 0.0;
  }
}

void X_pphp::Gather_X_2_2(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map22_num[ind + offset]; ++n){
      index = this->map22_index[ind + offset] + n;
      index1 = this->map22_ind[index];
      this->X_3_3[index1] += this->map22_fac1[index] * X2[ind];
    }
    X2[ind] = 0.0;
  }
}

void X_pphp::Gather_X_2_3(double *X2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      this->X_3_3[index1] += this->map23_fac1[index] * X2[ind];
    }
    X2[ind] = 0.0;
  }
}

void X_pphp::Gather_X_1(double *X1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map1_ind[ind + offset];
    this->X_3_3[index] += this->map1_fac1[ind + offset] * X1[ind];
    X1[ind] = 0.0;
  }
}

void X_pphp::X1_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_3_length; ++ind){ this->X1_3_3[ind] = 0.0; } }
}

void X_pphp::X_Zero(bool flag)
{
  if( !flag ){ for(int ind = 0; ind < this->X_3_3_length; ++ind){ this->X_3_3[ind] = 0.0; } }
}

void X_pphp::Delete()
{
  delete[] this->X1_3_3;
  delete[] this->X_3_3;

  delete[] this->X_1_index;
  delete[] this->X_2_2_index;
  delete[] this->X_2_3_index;
  delete[] this->X_3_1_index;
  delete[] this->X_3_2_index;
  delete[] this->X_3_3_index;
  delete[] this->X_3_4_index;

  delete[] this->map1_ind;
  delete[] this->map1_fac1;

  delete[] this->map22_index;
  delete[] this->map22_num;
  delete[] this->map22_ind;
  delete[] this->map22_fac1;

  delete[] this->map23_index;
  delete[] this->map23_num;
  delete[] this->map23_ind;
  delete[] this->map23_fac1;

  delete[] this->map31_ind;
  delete[] this->map31_fac1;

  delete[] this->map32_ind;
  delete[] this->map32_fac1;

  delete[] this->map34_ind;
  delete[] this->map34_fac1;
  delete[] this->map34_fac2;
}

void Eff_Interactions::Update_1(Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m12 = -0.5;
  int nh, np, npph, nhhp, nhh, npp, nhp1, nph1, nppp, nhhh, length;
  int npp0 = Chan.npp1[Chan.ind0];
  int nhh0 = Chan.nhh1[Chan.ind0];
  int p_ind, h_ind;
  int chan_ind1, chan_ind2, chan_ind3;
  double *Xpp, *Xhh, *T2, *T3;
  char N = 'N';

  //////////////   Xpp   ///////////////
  if(PAR.HF == 1){
    // Xpp2(a|b){ab}  <-  fpp(a|a) del_(ab)
    for(int p = SPB.num_hol; p < SPB.num_states; ++p){
      p_ind = Chan.pp1_map[Chan.ind0][Hash(p, p, 0)];
      this->Xpp.X_2[p_ind] += SPB.qnums[p].energy * std::sqrt(SPB.qnums[p].j + 1.0);
    }
  }
  else if(PAR.HF == 0){
    // Xpp2(a|b){ab}  <-  fpp(a|b){ab}
    for(int pp1 = 0; pp1 < npp0; ++pp1){
      this->Xpp.X_2[pp1] += Ints.Fmatrix.pp_2[pp1];
    }
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    if(np * nhhp != 0){
      chan_ind1 = Amps.D1.T3_1_index[chan3];
      chan_ind2 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind3 = this->Xpp.X_3_index[chan3];

      T3 = new double[np * nhhp];
      Amps.D1.Set_T3_1(T3, np*nhhp, chan_ind1);
      Xpp = new double[np * np];
      for(int ind = 0; ind < np*np; ++ind){ Xpp[ind] = 0.0; }
      // Xpp3(a|b){a,b}  =  -(1/2) T3_1(ac|kl){a,klc'}.Vhhpp3_3(kl|bc){klc',b}
      dgemm_NN(T3, (Ints.Vhhpp.V_3_3 + chan_ind2), Xpp, &np, &np, &nhhp, &m12, &p1, &N, &N);
      this->Xpp.Gather_X_3(Xpp, np*np, chan_ind3);

      delete[] Xpp;
      delete[] T3;
    }
  }

  //////////////   X1hh,Xhh   ///////////////
  if(PAR.HF == 1){
    // Xhh2(i|j){ij}  <-  fhh(i|i) del_(ij)
    for(int h = 0; h < SPB.num_hol; ++h){
      h_ind = Chan.hh1_map[Chan.ind0][Hash(h, h, 0)];
      this->Xhh.X1_2[h_ind] += SPB.qnums[h].energy * std::sqrt(SPB.qnums[h].j + 1.0);
      this->Xhh.X_2[h_ind] += SPB.qnums[h].energy * std::sqrt(SPB.qnums[h].j + 1.0);
    }
  }
  else if(PAR.HF == 0){
    // Xhh2(i|j){ij}  <-  fhh(i|j){ij}
    for(int hh1 = 0; hh1 < nhh0; ++hh1){
      this->Xhh.X1_2[hh1] += Ints.Fmatrix.hh_2[hh1];
      this->Xhh.X_2[hh1] += Ints.Fmatrix.hh_2[hh1];
    }
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    npph = Chan.npph[chan3];
    if(nh * npph != 0){
      chan_ind1 = Ints.Vhhpp.V_3_1_index[chan3];
      chan_ind2 = Amps.D1.T3_3_index[chan3];
      chan_ind3 = this->Xhh.X_3_index[chan3];

      T3 = new double[npph * nh];
      Amps.D1.Set_T3_3(T3, npph*nh, chan_ind2);
      Xhh = new double[nh * nh];
      for(int ind = 0; ind < nh*nh; ++ind){ Xhh[ind] = 0.0; }
      // X1hh3(i|j){i,j}  =  +(1/2) Vhhpp3_1(ik|cd){i,cdk'}.T3_3(cd|jk){cdk',j}
      dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), T3, Xhh, &nh, &nh, &npph, &p12, &p1, &N, &N);
      this->Xhh.Gather_X1_3(Xhh, nh*nh, chan_ind3);
      // Xhh3(i|j){i,j}  =  +(1/2) Vhhpp3_1(ik|cd){i,cdk'}.T3_3(cd|jk){cdk',j}
      dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), T3, Xhh, &nh, &nh, &npph, &p12, &p1, &N, &N);
      this->Xhh.Gather_X_3(Xhh, nh*nh, chan_ind3);

      delete[] Xhh;
      delete[] T3;
    }
  }

  //////////////   X1hhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    nhhp = Chan.nhhp[chan3];
    length = nhhp * nh;
    chan_ind1 = this->Xhhhp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vhhhp.V_3_3_index[chan3];
    // X1hhhp3_3 = Vhhhp3_3
    for(int ind = 0; ind < length; ++ind){
      this->Xhhhp.X1_3_3[chan_ind1 + ind] += Ints.Vhhhp.V_3_3[chan_ind2 + ind];
    }
  }

  //////////////   X1hppp,Xhppp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    length = np * npph;
    chan_ind1 = this->Xhppp.X_3_2_index[chan3];
    chan_ind2 = Ints.Vhppp.V_3_2_index[chan3];
    // X1hppp3_2 = Vhppp3_2, Xhppp3_2 = Vhppp3_2
    for(int ind = 0; ind < length; ++ind){
      this->Xhppp.X1_3_2[chan_ind1 + ind] += Ints.Vhppp.V_3_2[chan_ind2 + ind];
      this->Xhppp.X_3_2[chan_ind1 + ind] += Ints.Vhppp.V_3_2[chan_ind2 + ind];
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];

    //////////////   Xhhhh   ///////////////
    length = nhh * nhh;
    chan_ind1 = this->Xhhhh.X_1_index[chan1];
    chan_ind2 = Ints.Vhhhh.V_1_index[chan1];
    // Xhhhh1(ij|kl){ij,kl}  =  Vhhhh1(ij|kl){ij,kl}
    for(int ind = 0; ind < length; ++ind){
      this->Xhhhh.X_1[chan_ind1 + ind] += Ints.Vhhhh.V_1[chan_ind2 + ind];
    }
    if(nhh * npp != 0){
      chan_ind1 = Ints.Vhhpp.V_1_index[chan1];
      chan_ind2 = Amps.D1.T1_index[chan1];
      chan_ind3 = this->Xhhhh.X_1_index[chan1];
      // Xhhhh1(ij|kl){ij,kl}  <-  +(1/2) Vhhpp1(ij|cd){ij,cd}.T1(cd|kl){cd,kl}
      dgemm_NN((Ints.Vhhpp.V_1 + chan_ind1), (Amps.D1.T1 + chan_ind2), (this->Xhhhh.X_1 + chan_ind3), &nhh, &nhh, &npp, &p12, &p1, &N, &N);
    }

    //////////////   X1pppp   ///////////////
    length = npp * npp;
    chan_ind1 = this->Xpppp.X_1_index[chan1];
    chan_ind2 = Ints.Vpppp.V_1_index[chan1];
    // X1pppp1(ab|cd){ab,cd}  =  Vpppp1(ab|cd){ab,cd}
    for(int ind = 0; ind < length; ++ind){
      this->Xpppp.X1_1[chan_ind1 + ind] += Ints.Vpppp.V_1[chan_ind2 + ind];
    }
  }

  //////////////   X1hphp,X2hphp,Xhphp   ///////////////
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    length = nhp1 * nhp1;
    chan_ind1 = this->Xhphp.X_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    // X(1)(2)hphp2_1(ia|jb){ib',ja'}  =  Vhphp2_1(ia|jb){ib',ja'}
    for(int ind = 0; ind < length; ++ind){
      this->Xhphp.X1_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
      this->Xhphp.X2_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
      this->Xhphp.X_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
    }
    if(nhp1 * nph1 != 0){
      chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
      chan_ind2 = Amps.D1.T2_1_index[chan2];
      chan_ind3 = this->Xhphp.X_2_1_index[chan2];

      T2 = new double[nph1 * nhp1];
      Amps.D1.Set_T2_1(T2, nph1*nhp1, chan_ind2);
      // Xhphp2_1(ia|jb){ib',ja'}  <-  -(1/2) Vhhpp2_1(ik|cb){ib',ck'}.T2_1(ca|jk){ck',ja'}
      dgemm_NN((Ints.Vhhpp.V_2_1 + chan_ind1), T2, (this->Xhphp.X_2_1 + chan_ind3), &nhp1, &nhp1, &nph1, &m12, &p1, &N, &N);

      delete[] T2;
    }
  }

  //////////// Xpphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    nppp = Chan.nppp[chan3];
    length = nppp * nh;
    chan_ind1 = this->Xpphp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vpphp.V_3_3_index[chan3];
    /*if( chan3 == 17 ){
      std::cout << "Xpphp.X1_3_3(0) = " << this->Xpphp.X1_3_3[762] << std::endl;
      }*/

    for(int ind = 0; ind < length; ++ind){
      this->Xpphp.X1_3_3[chan_ind1 + ind] += Ints.Vpphp.V_3_3[chan_ind2 + ind];
    }

    /*if( chan3 == 17 ){
      std::cout << "Xpphp.X1_3_3(1) = " << this->Xpphp.X1_3_3[762] << std::endl;
      }*/
  }

  //////////// Xhphh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    length = np * nhhh;
    chan_ind1 = this->Xhphh.X_3_2_index[chan3];
    chan_ind2 = Ints.Vhphh.V_3_2_index[chan3];
    // X1hphh3_2 = Vhphh3_2
    for(int ind = 0; ind < length; ++ind){
      this->Xhphh.X1_3_2[chan_ind1 + ind] += Ints.Vhphh.V_3_2[chan_ind2 + ind];
    }
  }
}

void Eff_Interactions::Update_2(Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int nh, np, npph, nhhp, nppp, nhhh, nhpp, nhph, one = 1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4, index;
  double *t3, *Xpp, *Xhh, *Xhp, *Xhppp, *Xhhhp, *Xpppp, *Xhhhh, *Xhphp, *Xpphp, *Vpppp;
  int nhp0 = Chan.nhp1[Chan.ind0];
  int nph0 = Chan.nph1[Chan.ind0];
  int npp0 = Chan.npp1[Chan.ind0];
  int nhh0 = Chan.nhh1[Chan.ind0];
  char N = 'N';

  ////////////   Xhp   ///////////////
  if(nhp0 != 0){
    if(PAR.HF == 0){
      // Xhp2(i|a){ia'}  <-  fhp(i|a){ia'}
      for(int hp1 = 0; hp1 < nhp0; ++hp1){ this->Xhp.X_2[hp1] += Ints.Fmatrix.hp_2[hp1]; }
    }
    chan_ind1 = Ints.Vhhpp.V_2_3_index[Chan.ind0];
    // Xhp2(i|a){ia'}  =  Vhhpp2_3(ik|ac){ia',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhpp.V_2_3 + chan_ind1), Amps.S1.t2, this->Xhp.X_2, &nhp0, &one, &nph0, &p1, &p1, &N, &N);
  }
 
  //////////////   Xpp   ///////////////
  if(nph0 * npp0 != 0){
    chan_ind1 = Ints.Vhppp.V_2_4_index[Chan.ind0];
    // Xpp2(a|b){ab'}  =  Vhppp2_4(ka|cb){ab',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhppp.V_2_4 + chan_ind1), Amps.S1.t2, this->Xpp.X_2, &npp0, &one, &nph0, &p1, &p1, &N, &N);
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    if(nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = this->Xhp.X_3_index[chan3];
      chan_ind3 = this->Xpp.X_3_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind1);
      Xhp = new double[nh * np];
      this->Xhp.Set_X_3(Xhp, nh*np, chan_ind2);
      Xpp = new double[np * np];
      for(int ind = 0; ind < np*np; ++ind){ Xpp[ind] = 0.0; }
      // Xpp3(a|b){a,b}  <-  - t3(a|k){a,k}.Xhp3(k|b){k,b}
      dgemm_NN(t3, Xhp, Xpp, &np, &np, &nh, &m1, &p1, &N, &N);
      this->Xpp.Gather_X_3(Xpp, np*np, chan_ind3);

      delete[] t3;
      delete[] Xhp;
      delete[] Xpp;
    }
  }

  //////////////   Xhh   ///////////////
  if(nph0 * nhh0 != 0){
    chan_ind1 = Ints.Vhhhp.V_2_3_index[Chan.ind0];
    // X1hh2(i|j){ij'}  =  Vhhhp2_3(ik|jc){ij',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhhp.V_2_3 + chan_ind1), Amps.S1.t2, this->Xhh.X1_2, &nhh0, &one, &nph0, &p1, &p1, &N, &N);
    // Xhh2(i|j){ij'}  =  Vhhhp2_3(ik|jc){ij',ck'}.t2(c|k){ck'}
    dgemm_NN((Ints.Vhhhp.V_2_3 + chan_ind1), Amps.S1.t2, this->Xhh.X_2, &nhh0, &one, &nph0, &p1, &p1, &N, &N);
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    if(nh * np != 0){
      chan_ind1 = this->Xhp.X_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = this->Xhh.X_3_index[chan3];

      Xhp = new double[nh * np];
      this->Xhp.Set_X_3(Xhp, nh*np, chan_ind1);
      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind2);
      Xhh = new double[nh * nh];
      for(int ind = 0; ind < nh*nh; ++ind){ Xhh[ind] = 0.0; }
      // Xhh3(i|j){i,j}  <-  Xhp3(i|c){i,c}.t3(c|j){c,j}
      dgemm_NN(Xhp, t3, Xhh, &nh, &nh, &np, &p1, &p1, &N, &N);
      this->Xhh.Gather_X_3(Xhh, nh*nh, chan_ind3);

      delete[] t3;
      delete[] Xhp;
      delete[] Xhh;
    }
  }

  //////////////   X1hppp,Xhppp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    if(nh * np * npph != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhpp.V_3_2_index[chan3];
      chan_ind3 = this->Xhppp.X_3_2_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind1);
      // X1hppp3_2(ia|bc){a,bci'}  <-  - (1/2) t3(a|k){a,k}.Vhhpp3_2(ik|bc){k,bci'}
      dgemm_NN(t3, (Ints.Vhhpp.V_3_2 + chan_ind2), (this->Xhppp.X1_3_2 + chan_ind3), &np, &npph, &nh, &m12, &p1, &N, &N);
      // Xhppp3_2(ia|bc){a,bci'}  <-  - t3(a|k){a,k}.Vhhpp3_2(ik|bc){k,bci'}
      dgemm_NN(t3, (Ints.Vhhpp.V_3_2 + chan_ind2), (this->Xhppp.X_3_2 + chan_ind3), &np, &npph, &nh, &m1, &p1, &N, &N);

      delete[] t3;
    }
  }

  //////////////   Xhhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    if(nh * np * nhhp != 0){
      chan_ind1 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = this->Xhhhp.X_3_3_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind2);
      // X1hhhp3_3(ik|ja){ika',j}  <-  +(1/2) Vhhpp3_3(ik|ca){ika',c}.t3(c|j){c,j}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind1), t3, (this->Xhhhp.X1_3_3 + chan_ind3), &nhhp, &nh, &np, &p12, &p1, &N, &N);

      delete[] t3;
    }
  }

  //////////////   X1pppp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    if(nh * np * nppp){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = this->Xhppp.X_3_1_index[chan3];
      chan_ind3 = this->Xpppp.X_3_1_index[chan3];
      chan_ind4 = this->Xpppp.X_3_2_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind1);
      Xhppp = new double[nh * nppp];
      this->Xhppp.Set_X1_3_1(Xhppp, nh*nppp, chan_ind2);
      Xpppp = new double[np * nppp];
      for(int ind = 0; ind < np*nppp; ++ind){ Xpppp[ind] = 0.0; }
      // X1pppp3_1(ab|cd){a,cdb'}  =  -t3(a|k){a,k}.X1hppp3_1(kb|cd){k,cdb'}
      dgemm_NN(t3, Xhppp, Xpppp, &np, &nppp, &nh, &m1, &p1, &N, &N);
      this->Xpppp.Gather_X1_3_1(Xpppp, np*nppp, chan_ind3);
      // X1pppp3_2(ab|cd){b,cda'}  =  t3(b|k){b,k}.X1hppp3_1(ka|cd){k,cda'}
      dgemm_NN(t3, Xhppp, Xpppp, &np, &nppp, &nh, &p1, &p1, &N, &N);
      this->Xpppp.Gather_X1_3_2(Xpppp, np*nppp, chan_ind4);

      delete[] t3;
      delete[] Xpppp;
      delete[] Xhppp;
    }
  }

  //////////////   Xhhhh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];
    if(nh * np * nhhh){
      chan_ind1 = this->Xhhhp.X_3_4_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = this->Xhhhh.X_3_4_index[chan3];
      chan_ind4 = this->Xhhhh.X_3_3_index[chan3];

      Xhhhp = new double[nhhh * np];
      this->Xhhhp.Set_X1_3_4(Xhhhp, nhhh*np, chan_ind1);
      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind2);
      Xhhhh = new double[nhhh * nh];
      for(int ind = 0; ind < nhhh*nh; ++ind){ Xhhhh[ind] = 0.0; }
      // Xhhhh3_4(ij|kl){ijk',l}  =  X1hhhp3_4(ij|kc){ijk',c}.t3(c|l){c,l}
      dgemm_NN(Xhhhp, t3, Xhhhh, &nhhh, &nh, &np, &p1, &p1, &N, &N);
      this->Xhhhh.Gather_X_3_4(Xhhhh, nhhh*nh, chan_ind3);
      // Xhhhh3_3(ij|kl){ijk',l}  =  -X1hhhp3_4(ij|kc){ijk',c}.t3(c|l){c,l}
      dgemm_NN(Xhhhp, t3, Xhhhh, &nhhh, &nh, &np, &m1, &p1, &N, &N);
      this->Xhhhh.Gather_X_3_3(Xhhhh, nhhh*nh, chan_ind4);

      delete[] t3;
      delete[] Xhhhh;
      delete[] Xhhhp;
    }
  }

  //////////////   X1hphp,X2hphp,Xhphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhpp = Chan.nhpp[chan3];
    nhph = Chan.nhph[chan3];
    
    chan_ind1 = Amps.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps.S1.Set_t3(t3, np*nh, chan_ind1);
    if(nhpp * nh * np != 0){
      chan_ind2 = this->Xhppp.X_3_3_index[chan3];
      chan_ind3 = this->Xhphp.X_3_3_index[chan3];

      Xhppp = new double[nhpp * np];
      this->Xhppp.Set_X1_3_3(Xhppp, nhpp*np, chan_ind2);
      Xhphp = new double[nhpp * nh];
      for(int ind = 0; ind < nhpp*nh; ++ind){ Xhphp[ind] = 0.0; }
      // X1hphp3_3(ia|jb){iab',j}  =  X1hppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN(Xhppp, t3, Xhphp, &nhpp, &nh, &np, &p1, &p1, &N, &N);
      this->Xhphp.Gather_X1_3_3(Xhphp, nhpp*nh, chan_ind3);
      // X2hphp3_3(ia|jb){iab',j}  =  +(1/2) X1hppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN(Xhppp, t3, Xhphp, &nhpp, &nh, &np, &p12, &p1, &N, &N);
      this->Xhphp.Gather_X2_3_3(Xhphp, nhpp*nh, chan_ind3);

      this->Xhppp.Set_X_3_3(Xhppp, nhpp*np, chan_ind2);
      // Xhphp3_3(ia|jb){iab',j}  =  Xhppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN(Xhppp, t3, Xhphp, &nhpp, &nh, &np, &p1, &p1, &N, &N);
      this->Xhphp.Gather_X_3_3(Xhphp, nhpp*nh, chan_ind3);

      delete[] Xhphp;
      delete[] Xhppp;
    }
    if(nhph * nh * np != 0){
      chan_ind2 = Ints.Vhhhp.V_3_2_index[chan3];
      chan_ind3 = this->Xhphp.X_3_2_index[chan3];

      Xhphp = new double[np * nhph];
      for(int ind = 0; ind < np*nhph; ++ind){ Xhphp[ind] = 0.0; }
      // X1hphp3_2(ia|jb){a,jbi'}  =  -(1/2) t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN(t3, (Ints.Vhhhp.V_3_2 + chan_ind2), Xhphp, &np, &nhph, &nh, &m12, &p1, &N, &N);
      this->Xhphp.Gather_X1_3_2(Xhphp, np*nhph, chan_ind3);
      // X2hphp3_2(ia|jb){a,jbi'}  =  -(1/2) t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN(t3, (Ints.Vhhhp.V_3_2 + chan_ind2), Xhphp, &np, &nhph, &nh, &m12, &p1, &N, &N);
      this->Xhphp.Gather_X2_3_2(Xhphp, np*nhph, chan_ind3);
      // Xhphp3_2(ia|jb){a,jbi'}  =  - t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN(t3, (Ints.Vhhhp.V_3_2 + chan_ind2), Xhphp, &np, &nhph, &nh, &m1, &p1, &N, &N);
      this->Xhphp.Gather_X_3_2(Xhphp, np*nhph, chan_ind3);

      delete[] Xhphp;
    }
    delete[] t3;
  }

  //////////////   X1pphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    nhpp = Chan.nhpp[chan3];

    chan_ind1 = Amps.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps.S1.Set_t3(t3, np*nh, chan_ind1);
    if(nppp * nh * np != 0){
      chan_ind2 = this->Xpppp.X_3_1_index[chan3]; // for V3_3 -> V3_1
      chan_ind3 = this->Xpphp.X_3_3_index[chan3];

      Vpppp = new double[np*nppp];
      for(int p = 0; p < np; ++p){
	for(int ppp = 0; ppp < nppp; ++ppp){
	  if( this->Xpppp.map31_fac1[chan_ind2 + (p*nppp + ppp)] == 0.0){ Vpppp[ppp*np + p] = 0.0; continue; }
	  index = this->Xpppp.map31_ind[chan_ind2 + (p*nppp + ppp)];
	  Vpppp[ppp*np + p] = Ints.Vpppp.V_1[index] / this->Xpppp.map31_fac1[chan_ind2 + (p*nppp + ppp)];
	}
      }

      /*if( chan3 == 17 ){
	std::cout << "Xpphp.X1_3_3(2) = " << this->Xpphp.X1_3_3[762] << std::endl;
	}*/

      /*if( chan3 == 16 ){
	std::cout << "Vpppp = " << std::endl;
	for(int p = 0; p < np; ++p){
	  for(int h = 0; h < nh; ++h){
	    std::cout << t3[p*nh + h] << " ";
	  }
	}
	std::cout << std::endl << std::endl;
	std::cout << "t3 = " << std::endl;
	for(int p = 0; p < np; ++p){
	  for(int h = 0; h < nh; ++h){
	    std::cout << t3[p*nh + h] << " ";
	  }
	}
	std::cout << std::endl << std::endl;
	}*/

      // X1pphp3_3(ab|ic){abc',i}  <-  +(1/2) Vpppp3_3(ab|dc){abc',d}.t3(d|i){d,i}
      dgemm_NN(Vpppp, t3, (this->Xpphp.X1_3_3 + chan_ind3), &nppp, &nh, &np, &p12, &p1, &N, &N);
      /*if( chan3 == 17 ){
	std::cout << "Xpphp.X1_3_3(3) = " << this->Xpphp.X1_3_3[762] << std::endl;
	}*/
      
      delete[] Vpppp;
    }
    //std::cout << "Xpphp.X1_3_3(3.5) = " << this->Xpphp.X1_3_3[762] << std::endl;
    if(nhpp * nh * np != 0){
      chan_ind2 = this->Xhphp.X_3_1_index[chan3];
      chan_ind3 = this->Xpphp.X_3_1_index[chan3];
      chan_ind4 = this->Xpphp.X_3_2_index[chan3];

      Xhphp = new double[nh * nhpp];
      this->Xhphp.Set_X2_3_1(Xhphp, nh*nhpp, chan_ind2);
      Xpphp = new double[np * nhpp];
      for(int ind = 0; ind < np * nhpp; ++ind){ Xpphp[ind] = 0.0; }
      // X1pphp3_1(ab|ic){a,icb'}  =  -t3(a|k){a,l}.X2hphp3_1(kb|ic){k,icb'}
      dgemm_NN(t3, Xhphp, Xpphp, &np, &nhpp, &nh, &m1, &p1, &N, &N);
      /*if( chan3 == 2 ){
	std::cout << "Xpphp.X1_3_1(4) = " << Xpphp[31] << std::endl;
	}*/
      this->Xpphp.Gather_X1_3_1(Xpphp, np*nhpp, chan_ind3);

      /*if( chan3 == 2 ){
	std::cout << "Xpphp.X1_3_3(5) = " << this->Xpphp.X1_3_3[762] << std::endl;
	}*/

      // X1pphp3_2(ab|ic){b,ica'}  =  t3(b|k){b,l}.X2hphp3_1(ka|ic){k,ica'}
      dgemm_NN(t3, Xhphp, Xpphp, &np, &nhpp, &nh, &p1, &p1, &N, &N);
      /*if( chan3 == 7 ){
	std::cout << "Xpphp.X1_3_2(6) = " << Xpphp[0] << std::endl;
	}*/
      this->Xpphp.Gather_X1_3_2(Xpphp, np*nhpp, chan_ind4);

      /*if( chan3 == 7 ){
	std::cout << "Xpphp.X1_3_3(7) = " << this->Xpphp.X1_3_3[762] << std::endl;
	}*/

      delete[] Xhphp;
      delete[] Xpphp;
    }
    delete[] t3;
  }

  //std::cout << "Xpphp.X1_3_3(8) = " << this->Xpphp.X1_3_3[762] << std::endl;

  //////////////   X1hphh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhh = Chan.nhhh[chan3];

    if(nhhh * nh * np != 0){
      chan_ind1 = Amps.S1.t3_index[chan3];
      chan_ind2 = Ints.Vhhhh.V_3_2_index[chan3];
      chan_ind3 = this->Xhphh.X_3_2_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind1);
      // X1hphh3_2(ia|jk){a,jki'}  <-  -(1/2) t3(a|l){a,l}.Vhhhh3_2(il|jk){l,jki'}
      dgemm_NN(t3, (Ints.Vhhhh.V_3_2 + chan_ind2), (this->Xhphh.X1_3_2 + chan_ind3), &np, &nhhh, &nh, &m12, &p1, &N, &N);

      delete[] t3;
    }
  }
}


void Eff_Interactions::Update_3(Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int nh, np, npph, nhhp, nhph, nhpp, nhh, npp, nhp, nhp1, nph1, npp1, nhh1, nhhh, nppp, length;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4, index;
  double *t3, *T2, *T3, *Xhppp, *Xhhhp, *Xhphp, *Xhphh, *Xpphp, *Vpppp, *Xhp;
  char N = 'N';

  /*int a, b, i, j;
  for(int ind = 0; ind < Chan.nhh1[Chan.ind0]; ++ind){
    i = Chan.hh1_state(Chan.ind0, ind).v1;
    j = Chan.hh1_state(Chan.ind0, ind).v2;
    if(i == j){ std::cout << "Xhh_" << i << " = " << this->Xhh.X_2[ind] << std::endl; }
  }
  std::cout << std::endl;
  for(int ind = 0; ind < Chan.npp1[Chan.ind0]; ++ind){
    a = Chan.pp1_state(Chan.ind0, ind).v1;
    b = Chan.pp1_state(Chan.ind0, ind).v2;
    if(a == b){ std::cout << "Xpp_" << a << " = " << this->Xpp.X_2[ind] << std::endl; }
  }
  std::cout << std::endl;*/

  this->Xpppp.X_1 = Xpppp.X1_1;  // reassign X1_1 -> X_1

  //////////////   Xhhhp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhhp = Chan.nhhp[chan3];
    length = nhhp * nh;
    chan_ind1 = this->Xhhhp.X_3_3_index[chan3];
    chan_ind2 = Ints.Vhhhp.V_3_3_index[chan3];
    // Xhhhp3_3 = Vhhhp3_3
    for(int ind = 0; ind < length; ++ind){
      this->Xhhhp.X_3_3[chan_ind1 + ind] += Ints.Vhhhp.V_3_3[chan_ind2 + ind];
    }
    if(nh * np * nhhp != 0){
      chan_ind1 = Ints.Vhhpp.V_3_3_index[chan3];
      chan_ind2 = Amps.S1.t3_index[chan3];
      chan_ind3 = this->Xhhhp.X_3_3_index[chan3];

      t3 = new double[np * nh];
      Amps.S1.Set_t3(t3, np*nh, chan_ind2);
      // Xhhhp3_3(ik|ja){ika',j}  <-  + Vhhpp3_3(ik|ca){ika',c}.t3(c|j){c,j}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind1), t3, (this->Xhhhp.X_3_3 + chan_ind3), &nhhp, &nh, &np, &p1, &p1, &N, &N);

      delete[] t3;
    }
  }

  //////////////   Xpppp   ///////////////
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    length = npp * npp;
    chan_ind1 = this->Xpppp.X_1_index[chan1];
    if(nhh * npp != 0){
      chan_ind1 = Amps.D1.T1_index[chan1];
      chan_ind2 = Ints.Vhhpp.V_1_index[chan1];
      chan_ind3 = this->Xpppp.X_1_index[chan1];
      // Xpppp1(ab|cd){ab,cd}  <-  +(1/2) T1(ab|kl){ab,kl}.Vhhpp1(kl|cd){kl,cd}
      dgemm_NN((Amps.D1.T1 + chan_ind1), (Ints.Vhhpp.V_1 + chan_ind2), (this->Xpppp.X_1 + chan_ind3), &npp, &npp, &nhh, &p12, &p1, &N, &N);
    }
  }

  //////////////   X3hphp,Xhphp   ///////////////
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    length = nhp1 * nhp1;
    chan_ind1 = this->Xhphp.X_2_1_index[chan2];
    chan_ind2 = Ints.Vhphp.V_2_1_index[chan2];
    // X3hphp2_1 = Vhphp2_1
    for(int ind = 0; ind < length; ++ind){
      this->Xhphp.X3_2_1[chan_ind1 + ind] += Ints.Vhphp.V_2_1[chan_ind2 + ind];
    }
    if(nhp1 * nph1 != 0){
      chan_ind1 = Ints.Vhhpp.V_2_1_index[chan2];
      chan_ind2 = Amps.D1.T2_1_index[chan2];
      chan_ind3 = this->Xhphp.X_2_1_index[chan2];

      T2 = new double[nph1 * nhp1];
      Amps.D1.Set_T2_1(T2, nph1*nhp1, chan_ind2);
      // Xhphp2_1(ia|jb){ib',ja'}  <-  -(1/2) Vhhpp2_1(ik|cb){ib',ck'}.T2_1(ca|jk){ck',ja'}
      dgemm_NN((Ints.Vhhpp.V_2_1 + chan_ind1), T2, (this->Xhphp.X_2_1 + chan_ind3), &nhp1, &nhp1, &nph1, &m12, &p1, &N, &N);

      delete[] T2;
    }
  }
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhpp = Chan.nhpp[chan3];
    nhph = Chan.nhph[chan3];
    chan_ind1 = Amps.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps.S1.Set_t3(t3, np*nh, chan_ind1);
    if(nhpp * nh * np != 0){
      chan_ind2 = this->Xhppp.X_3_3_index[chan3];
      chan_ind3 = this->Xhphp.X_3_3_index[chan3];

      Xhppp = new double[nhpp * np];
      this->Xhppp.Set_X_3_3(Xhppp, nhpp*np, chan_ind2);
      Xhphp = new double[nhpp * nh];
      for(int ind = 0; ind < nhpp*nh; ++ind){ Xhphp[ind] = 0.0; }
      // X3hphp3_3(ia|jb){iab',j}  =  +(1/2) Xhppp3_3(ia|cb){iab',c}.t3(c|j){c,j}
      dgemm_NN(Xhppp, t3, Xhphp, &nhpp, &nh, &np, &p12, &p1, &N, &N);
      this->Xhphp.Gather_X3_3_3(Xhphp, nhpp*nh, chan_ind3);

      delete[] Xhppp;
      delete[] Xhphp;
    }
    if(nhph * nh * np != 0){
      chan_ind2 = Ints.Vhhhp.V_3_2_index[chan3];
      chan_ind3 = this->Xhphp.X_3_2_index[chan3];

      Xhphp = new double[np * nhph];
      for(int ind = 0; ind < np*nhph; ++ind){ Xhphp[ind] = 0.0; }
      // X3hphp3_2(ia|jb){a,jbi'}  =  -t3(a|k){a,k}.Vhhhp3_2(ik|jb){k,jbi'}
      dgemm_NN(t3, (Ints.Vhhhp.V_3_2 + chan_ind2), Xhphp, &np, &nhph, &nh, &m1, &p1, &N, &N);
      this->Xhphp.Gather_X3_3_2(Xhphp, np*nhph, chan_ind3);

      delete[] Xhphp;
    }
    delete[] t3;
  }

  //////////////   Xpphp   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nppp = Chan.nppp[chan3];
    nhpp = Chan.nhpp[chan3];
    npph = Chan.npph[chan3];

    chan_ind1 = Amps.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps.S1.Set_t3(t3, np*nh, chan_ind1);

    length = nppp * nh;
    chan_ind2 = this->Xpphp.X_3_3_index[chan3];
    chan_ind3 = Ints.Vpphp.V_3_3_index[chan3];
    // Xpphp3_3 = Vpphp3_3
    for(int ind = 0; ind < length; ++ind){
      this->Xpphp.X_3_3[chan_ind2 + ind] += Ints.Vpphp.V_3_3[chan_ind3 + ind];
    }
    if(nppp * nh * np != 0){
      chan_ind2 = this->Xpppp.X_3_1_index[chan3];
      chan_ind3 = this->Xpphp.X_3_3_index[chan3];

      Vpppp = new double[np*nppp];
      for(int p = 0; p < np; ++p){
	for(int ppp = 0; ppp < nppp; ++ppp){
	  if( this->Xpppp.map31_fac1[chan_ind2 + (p*nppp + ppp)] == 0.0){ Vpppp[ppp*np + p] = 0.0; continue; }
	  index = this->Xpppp.map31_ind[chan_ind2 + (p*nppp + ppp)];
	  Vpppp[ppp*np + p] = Ints.Vpppp.V_1[index] / this->Xpppp.map31_fac1[chan_ind2 + (p*nppp + ppp)];
	}
      }

      // Xpphp3_3(ab|ic){abc',i}  <-  Vpppp3_3(ab|dc){abc',d}.t3(d|i){d,i}
      dgemm_NN(Vpppp, t3, (this->Xpphp.X_3_3 + chan_ind3), &nppp, &nh, &np, &p1, &p1, &N, &N);

      delete[] Vpppp;
    }
    if(nhpp * nh * np != 0){
      chan_ind2 = this->Xhphp.X_3_1_index[chan3];
      chan_ind3 = this->Xpphp.X_3_1_index[chan3];
      chan_ind4 = this->Xpphp.X_3_2_index[chan3];

      Xhphp = new double[nh * nhpp];
      this->Xhphp.Set_X1_3_1(Xhphp, nh*nhpp, chan_ind2);
      Xpphp = new double[np * nhpp];
      for(int ind = 0; ind < np * nhpp; ++ind){ Xpphp[ind] = 0.0; }
      // Xpphp3_1(ab|ic){a,icb'}  =  -t3(a|k){a,l}.X1hphp3_1(kb|ic){k,icb'}
      dgemm_NN(t3, Xhphp, Xpphp, &np, &nhpp, &nh, &m1, &p1, &N, &N);
      this->Xpphp.Gather_X_3_1(Xpphp, np*nhpp, chan_ind3);
      // Xpphp3_2(ab|ic){a,icb'}  =  t3(a|k){a,l}.X1hphp3_1(kb|ic){k,icb'}
      dgemm_NN(t3, Xhphp, Xpphp, &np, &nhpp, &nh, &p1, &p1, &N, &N);
      this->Xpphp.Gather_X_3_2(Xpphp, np*nhpp, chan_ind4);

      delete[] Xhphp;
      delete[] Xpphp;
    }
    if(npph * np * nh != 0){
      chan_ind1 = Amps.D1.T3_4_index[chan3];
      chan_ind2 = this->Xhp.X_3_index[chan3];
      chan_ind3 = this->Xpphp.X_3_4_index[chan3];

      T3 = new double[npph * nh];
      Amps.D1.Set_T3_4(T3, npph*nh, chan_ind1);
      Xhp = new double[nh * np];
      this->Xhp.Set_X_3(Xhp, nh*np, chan_ind2);
      Xpphp = new double[npph * np];
      for(int ind = 0; ind < npph*np; ++ind){ Xpphp[ind] = 0.0; }
      // Xpphp3_4(ab|ic){abi,c}  <-  -T3_4(ab|ik){abi,k}.Xhp3(k|c){k,c}
      dgemm_NN(T3, Xhp, Xpphp, &npph, &np, &nh, &m1, &p1, &N, &N);
      this->Xpphp.Gather_X_3_4(Xpphp, npph*np, chan_ind3);

      delete[] Xhp;
      delete[] T3;
      delete[] Xpphp;
    }
    delete[] t3;
  }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    npp1 = Chan.npp1[chan2];
    if(nph1 * npp1 * nhp1 != 0){
      chan_ind1 = Amps.D1.T2_3_index[chan2];
      chan_ind2 = this->Xhppp.X_2_3_index[chan2];
      chan_ind3 = this->Xpphp.X_2_3_index[chan2];
      chan_ind4 = this->Xpphp.X_2_2_index[chan2];

      T2 = new double[nph1 * nhp1];
      Amps.D1.Set_T2_3(T2, nph1*nhp1, chan_ind1);
      Xhppp = new double[nhp1 * npp1];
      this->Xhppp.Set_X_2_3(Xhppp, nhp1*npp1, chan_ind2);
      Xpphp = new double[nph1 * npp1];
      for(int ind = 0; ind < nph1*npp1; ++ind){ Xpphp[ind] = 0.0; }
      // Xpphp2_3(ab|ic){ai',cb'}  <-  T2_3(ad|ik){ai',kd'}.Xhppp2_3(kb|dc){kd',cb'}
      dgemm_NN(T2, Xhppp, Xpphp, &nph1, &npp1, &nhp1, &p1, &p1, &N, &N);
      this->Xpphp.Gather_X_2_3(Xpphp, nph1*npp1, chan_ind3);
      // Xpphp2_2(ab|ic){bi',ca'}  <-  -T2_3(bd|ik){bi',kd'}.Xhppp2_3(ka|dc){kd',ca'}
      dgemm_NN(T2, Xhppp, Xpphp, &nph1, &npp1, &nhp1, &m1, &p1, &N, &N);
      this->Xpphp.Gather_X_2_2(Xpphp, nph1*npp1, chan_ind4);
      
      delete[] T2;
      delete[] Xhppp;
      delete[] Xpphp;
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    nhp = Chan.nhp[chan1];
    if(nhh * npp * nhp != 0){
      chan_ind1 = Amps.D1.T1_index[chan1];
      chan_ind2 = this->Xhhhp.X_1_index[chan1];
      chan_ind3 = this->Xpphp.X_1_index[chan1];

      Xhhhp = new double[nhh * nhp];
      this->Xhhhp.Set_X_1(Xhhhp, nhh*nhp, chan_ind2);
      Xpphp = new double[npp * nhp];
      for(int ind = 0; ind < npp*nhp; ++ind){ Xpphp[ind] = 0.0; }
      // Xpphp1(ab|ic){ab,ic}  <-  +(1/2) T1(ab|kl){ab,kl}.Xhhhp1(kl|ic){kl,ic}
      dgemm_NN((Amps.D1.T1 + chan_ind1), Xhhhp, Xpphp, &npp, &nhp, &nhh, &p12, &p1, &N, &N);
      this->Xpphp.Gather_X_1(Xpphp, npp*nhp, chan_ind3);

      delete[] Xhhhp;
      delete[] Xpphp;
    }
  }

  //////////////   Xhphh   ///////////////
  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    nhph = Chan.nhph[chan3];
    nhhp = Chan.nhhp[chan3];
    nhhh = Chan.nhhh[chan3];

    chan_ind1 = Amps.S1.t3_index[chan3];
    t3 = new double[np * nh];
    Amps.S1.Set_t3(t3, np*nh, chan_ind1);

    length = np * nhhh;
    chan_ind2 = this->Xhphh.X_3_2_index[chan3];
    chan_ind3 = Ints.Vhphh.V_3_2_index[chan3];
    // Xhphh3_2  <-  Vhphh3_2
    for(int ind = 0; ind < length; ++ind){
      this->Xhphh.X_3_2[chan_ind2 + ind] += Ints.Vhphh.V_3_2[chan_ind3 + ind];
    }
    if(nhhh * nh * np != 0){
      chan_ind2 = Ints.Vhhhh.V_3_2_index[chan3];
      chan_ind3 = this->Xhphh.X_3_2_index[chan3];
      // Xhphh3_2(ia|jk){a,jki'}  <-  - t3(a|l){a,l}.Vhhhh3_2(il|jk){l,jki'}
      dgemm_NN(t3, (Ints.Vhhhh.V_3_2 + chan_ind2), (this->Xhphh.X_3_2 + chan_ind3), &np, &nhhh, &nh, &m1, &p1, &N, &N);
    }
    if(nhph * np * nh != 0){
      chan_ind2 = this->Xhphp.X_3_4_index[chan3];
      chan_ind3 = this->Xhphh.X_3_4_index[chan3];
      chan_ind4 = this->Xhphh.X_3_3_index[chan3];

      Xhphp = new double[nhph * np];
      this->Xhphp.Set_X3_3_4(Xhphp, nhph*np, chan_ind2);
      Xhphh = new double[nhph * nh];
      for(int ind = 0; ind < nhph*nh; ++ind){ Xhphh[ind] = 0.0; }
      // Xhphh3_4(ia|jk){iaj,k}  <-  X3hphp3_4(ia|jc){iaj,c}.t3(c|k){c,k}
      dgemm_NN(Xhphp, t3, Xhphh, &nhph, &nh, &np, &p1, &p1, &N, &N);
      this->Xhphh.Gather_X_3_4(Xhphh, nhph*nh, chan_ind3);
      // Xhphh3_3(ia|jk){iak,j}  <-  -X3hphp3_4(ia|kc){iak,c}.t3(c|j){c,j}
      dgemm_NN(Xhphp, t3, Xhphh, &nhph, &nh, &np, &m1, &p1, &N, &N);
      this->Xhphh.Gather_X_3_3(Xhphh, nhph*nh, chan_ind4);

      delete[] Xhphp;
      delete[] Xhphh;
    }
    if(nh * nhhp * np != 0){
      chan_ind1 = this->Xhp.X_3_index[chan3];
      chan_ind2 = Amps.D1.T3_1_index[chan3];
      chan_ind3 = this->Xhphh.X_3_1_index[chan3];

      Xhp = new double[nh * np];
      this->Xhp.Set_X_3(Xhp, nh*np, chan_ind1);
      T3 = new double[np * nhhp];
      Amps.D1.Set_T3_1(T3, np*nhhp, chan_ind2);
      Xhphh = new double[nh * nhhp];
      for(int ind = 0; ind < nh*nhhp; ++ind){ Xhphh[ind] = 0.0; }
      // Xhphh3_1(ia|jk){i,jka}  <-  Xhp3(i|c){i,c}.T3_1(ca|jk){c,jka}
      dgemm_NN(Xhp, T3, Xhphh, &nh, &nhhp, &np, &p1, &p1, &N, &N);
      this->Xhphh.Gather_X_3_1(Xhphh, nh*nhhp, chan_ind3);

      delete[] Xhp;
      delete[] T3;
      delete[] Xhphh;
    }
    delete[] t3;
  }
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    nhh1 = Chan.nhh1[chan2];
    if(nhh1 * nhp1 * nph1 != 0){
      chan_ind1 = this->Xhhhp.X_2_3_index[chan2];
      chan_ind2 = Amps.D1.T2_3_index[chan2];
      chan_ind3 = this->Xhphh.X_2_3_index[chan2];
      chan_ind4 = this->Xhphh.X_2_1_index[chan2];

      T2 = new double[nph1 * nhp1];
      Amps.D1.Set_T2_3(T2, nph1*nhp1, chan_ind2);
      Xhhhp = new double[nhh1 * nph1];
      this->Xhhhp.Set_X_2_3(Xhhhp, nhh1*nph1, chan_ind1);
      Xhphh = new double[nhh1 * nhp1];
      for(int ind = 0; ind < nhh1*nhp1; ++ind){ Xhphh[ind] = 0.0; }
      // Xhphh2_3(ia|jk){ij',ka'}  <-  Xhhhp2_3(il|jc){ij',cl'}.T2_3(ca|lk){cl',ka'}
      dgemm_NN(Xhhhp, T2, Xhphh, &nhh1, &nhp1, &nph1, &p1, &p1, &N, &N);
      this->Xhphh.Gather_X_2_3(Xhphh, nhh1*nhp1, chan_ind3);
      // Xhphh2_1(ia|jk){ik',ja'}  <-  -Xhhhp2_3(il|kc){ik',cl'}.T2_3(ca|lj){cl',ja'}
      dgemm_NN(Xhhhp, T2, Xhphh, &nhh1, &nhp1, &nph1, &m1, &p1, &N, &N);
      this->Xhphh.Gather_X_2_1(Xhphh, nhh1*nhp1, chan_ind4);

      delete[] T2;
      delete[] Xhhhp;
      delete[] Xhphh;
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    nhh = Chan.nhh[chan1];
    nhp = Chan.nhp[chan1];
    if(nhh * npp * nhp != 0){
      chan_ind1 = this->Xhppp.X_1_index[chan1];
      chan_ind2 = Amps.D1.T1_index[chan1];
      chan_ind3 = this->Xhphh.X_1_index[chan1];

      Xhppp = new double[nhp * npp];
      this->Xhppp.Set_X_1(Xhppp, nhp*npp, chan_ind1);
      Xhphh = new double[nhp * nhh];
      for(int ind = 0; ind < nhp*nhh; ++ind){ Xhphh[ind] = 0.0; }
      // Xhphh1(ia|jk){ia,jk}  <-  +(1/2) Xhppp1(ia|cd){ia,cd}.T1(cd|jk){cd,jk}
      dgemm_NN(Xhppp, (Amps.D1.T1 + chan_ind2), Xhphh, &nhp, &nhh, &npp, &p12, &p1, &N, &N);
      this->Xhphh.Gather_X_1(Xhphh, nhp*nhh, chan_ind3);

      delete[] Xhppp;
      delete[] Xhphh;
    }
  }
}
