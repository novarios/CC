#include "CCfunctions.hpp"
#include "LAMBDAfunctions.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "INTfunctions.hpp"
#include "EffINTfunctions.hpp"
#include "DIISfunctions.hpp"
#include "MATHfunctions.hpp"

void LAmplitudes::Copy(LAmplitudes &LAmps)
{
  for(int ind = 0; ind < this->L1_length; ++ind){ this->L1[ind] = LAmps.L1[ind]; }
  for(int ind = 0; ind < this->l2_length; ++ind){ this->l2[ind] = LAmps.l2[ind]; }  
}

void LAmplitudes::Copy(LAmplitudes &LAmps, double *vec)
{
  for(int ind = 0; ind < this->L1_length; ++ind){
    this->L1[ind] = LAmps.L1[ind];
    vec[ind] = LAmps.L1[ind];
  }
  for(int ind = 0; ind < this->l2_length; ++ind){
    this->l2[ind] = LAmps.l2[ind];
    vec[this->L1_length + ind] = LAmps.l2[ind];
  }
}

void LAmplitudes::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  int i, j, a, b, nhh, npp, nhp1;
  two_body hh, pp, hp1, ia, ia_j;
  four_body ijab, ijab_j;
  four_body *fb_ind;
  four_body *fb_j;
  two_body *tb_ind;
  two_body *tb_j;
  int *J;
  std::unordered_map<int,int> *J21_map, *J22_map, *J23_map, *J24_map;

  /////////////////// l2 ///////////////////////////////////
  length = Chan.nhp1[Chan.ind0];
  this->l2 = new double[length];
  this->evec_ind = new int[2 * length];
  this->evec_fac = new double[2 * length];
  for(ind = 0; ind < length; ++ind){
    this->l2[ind] = 0.0;
  }
  this->l2_length = length;

  length = 0;
  this->l3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->l3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.np[chan3];
  }
  this->l3_length = length;

  tb_ind = new two_body[this->l2_length];
  tb_j = new two_body[this->l2_length];
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

  this->map3_ind = new int[l3_length];
  this->map3_fac1 = new double[l3_length];
  this->map3_fac2 = new double[l3_length];
  for(int ind = 0; ind < l3_length; ++ind){
    this->map3_ind[ind] = 0;
    this->map3_fac1[ind] = 0.0;
    this->map3_fac2[ind] = 0.0;
  }

  for(int hp = 0; hp < this->l2_length; ++hp){
    Map_2_to_3(this->map3_ind, this->map3_fac1, l3_index, hp, Chan.h_map, Chan.p_map, Chan.np, tb_ind, tb_j, 1);
    Map_2_to_3(this->map3_ind, this->map3_fac2, l3_index, hp, Chan.h_map, Chan.p_map, Chan.np, tb_ind, tb_j, 2);

    i = tb_ind[hp].v1;
    a = tb_ind[hp].v2;
    this->evec_ind[2*hp] = Chan.hh1_map[Chan.ind0][Hash(i, i, 0)];
    this->evec_fac[2*hp] = 1.0 / std::sqrt(tb_j[hp].v2 + 1.0);
    this->evec_ind[2*hp + 1] = Chan.pp1_map[Chan.ind0][Hash(a, a, 0)];
    this->evec_fac[2*hp + 1] = 1.0 / std::sqrt(tb_j[hp].v1 + 1.0);
  }
  delete[] tb_ind;
  delete[] tb_j;

  /////////////////// L1 ///////////////////////////////////
  length = 0;
  this->L1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->L1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.npp[chan1];
  }
  this->L1 = new double[length];
  this->Evec_ind = new int[4 * length];
  this->Evec_fac = new double[4 * length];
  for(ind = 0; ind < length; ++ind){
    this->L1[ind] = 0.0;
  }
  this->L1_length = length;

  length = 0;
  this->L2_1_index = new int[Chan.size2];
  this->L2_2_index = new int[Chan.size2];
  this->L2_3_index = new int[Chan.size2];
  this->L2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->L2_1_index[chan2] = length;
    this->L2_2_index[chan2] = length;
    this->L2_3_index[chan2] = length;
    this->L2_4_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nph1[chan2];
  }
  this->L2_1_length = length;
  this->L2_2_length = length;
  this->L2_3_length = length;
  this->L2_4_length = length;
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
  this->L3_1_index = new int[Chan.size3];
  this->L3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->L3_1_index[chan3] = length;
    this->L3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.npph[chan3];
  }
  this->L3_1_length = length;
  this->L3_2_length = length;

  length = 0;
  this->L3_3_index = new int[Chan.size3];
  this->L3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->L3_3_index[chan3] = length;
    this->L3_4_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.np[chan3];
  }
  this->L3_3_length = length;
  this->L3_4_length = length;

  fb_ind = new four_body[this->L1_length];
  fb_j = new four_body[this->L1_length];
  J = new int[this->L1_length];
  length = 0;
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    for(int hh_ind = 0; hh_ind < nhh; ++hh_ind){
      hh = Chan.hh_state(chan1, hh_ind);
      i = hh.v1;
      j = hh.v2;
      ijab.v1 = i;
      ijab.v2 = j;
      ijab_j.v1 = SPB.qnums[i].j;
      ijab_j.v2 = SPB.qnums[j].j;
      for(int pp_ind = 0; pp_ind < npp; ++pp_ind){
	pp = Chan.pp_state(chan1, pp_ind);
	a = pp.v1;
	b = pp.v2;
	ijab.v3 = a;
	ijab.v4 = b;
	ijab_j.v3 = SPB.qnums[a].j;
	ijab_j.v4 = SPB.qnums[b].j;
	fb_ind[length] = ijab;
	fb_j[length] = ijab_j;
	J[length] = Chan.qnums1[chan1].j;
	Map_Count_1_to_21(this->map21_num, this->L2_1_index, Chan.hp1_map, Chan.ph1_map, Chan.nph1, ijab, J[length], J21_map);
	Map_Count_1_to_22(this->map22_num, this->L2_2_index, Chan.hp1_map, Chan.ph1_map, Chan.nph1, ijab, J[length], J22_map);
	Map_Count_1_to_23(this->map23_num, this->L2_3_index, Chan.hp1_map, Chan.ph1_map, Chan.nph1, ijab, J[length], J23_map);
	Map_Count_1_to_24(this->map24_num, this->L2_4_index, Chan.hp1_map, Chan.ph1_map, Chan.nph1, ijab, J[length], J24_map);
	++length;
      }
    }
  }

  length = 0;
  for(int ind = 0; ind < L2_1_length; ++ind){
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
  for(int ind = 0; ind < L2_2_length; ++ind){
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
  for(int ind = 0; ind < L2_3_length; ++ind){
    this->map23_index[ind] = length;
    length += this->map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = 0;
    this->map23_fac1[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < L2_4_length; ++ind){
    this->map24_index[ind] = length;
    length += this->map24_num[ind];
  }
  this->map24_ind = new int[length];
  this->map24_fac1 = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->map24_ind[ind] = 0;
    this->map24_fac1[ind] = 0.0;
  }

  this->map31_ind = new int[L3_1_length];
  this->map31_fac1 = new double[L3_1_length];
  this->map31_fac2 = new double[L3_1_length];
  for(int ind = 0; ind < L3_1_length; ++ind){
    this->map31_ind[ind] = 0;
    this->map31_fac1[ind] = 0.0;
    this->map31_fac2[ind] = 0.0;
  }

  this->map32_ind = new int[L3_2_length];
  this->map32_fac1 = new double[L3_2_length];
  this->map32_fac2 = new double[L3_2_length];
  for(int ind = 0; ind < L3_2_length; ++ind){
    this->map32_ind[ind] = 0;
    this->map32_fac1[ind] = 0.0;
    this->map32_fac2[ind] = 0.0;
  }

  this->map33_ind = new int[L3_3_length];
  this->map33_fac1 = new double[L3_3_length];
  this->map33_fac2 = new double[L3_3_length];
  for(int ind = 0; ind < L3_3_length; ++ind){
    this->map33_ind[ind] = 0;
    this->map33_fac1[ind] = 0.0;
    this->map33_fac2[ind] = 0.0;
  }

  this->map34_ind = new int[L3_4_length];
  this->map34_fac1 = new double[L3_4_length];
  this->map34_fac2 = new double[L3_4_length];
  for(int ind = 0; ind < L3_4_length; ++ind){
    this->map34_ind[ind] = 0;
    this->map34_fac1[ind] = 0.0;
    this->map34_fac2[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static) private(i, j, a, b)
  for(int hhpp = 0; hhpp < this->L1_length; ++hhpp){
    Map_1_to_21(this->map21_ind,this->map21_fac1,this->map21_index, this->L2_1_index,hhpp, Chan.hp1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, J21_map,1);
    Map_1_to_21(this->map21_ind,this->map21_fac2,this->map21_index, this->L2_1_index,hhpp, Chan.hp1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, J21_map,2);
    Map_1_to_22(this->map22_ind,this->map22_fac1,this->map22_index, this->L2_2_index,hhpp, Chan.hp1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, J22_map,1);
    Map_1_to_23(this->map23_ind,this->map23_fac1,this->map23_index, this->L2_3_index,hhpp, Chan.hp1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, J23_map,1);
    Map_1_to_24(this->map24_ind,this->map24_fac1,this->map24_index, this->L2_4_index,hhpp, Chan.hp1_map,Chan.ph1_map,Chan.nph1, fb_ind,fb_j,J, J24_map,1);
    Map_1_to_31(this->map31_ind,this->map31_fac1,this->L3_1_index,hhpp, Chan.h_map,Chan.pph_map,Chan.npph, fb_ind,fb_j,J, 1);
    Map_1_to_31(this->map31_ind,this->map31_fac2,this->L3_1_index,hhpp, Chan.h_map,Chan.pph_map,Chan.npph, fb_ind,fb_j,J, 2);
    Map_1_to_32(this->map32_ind,this->map32_fac1,this->L3_2_index,hhpp, Chan.h_map,Chan.pph_map,Chan.npph, fb_ind,fb_j,J, 1);
    Map_1_to_32(this->map32_ind,this->map32_fac2,this->L3_2_index,hhpp, Chan.h_map,Chan.pph_map,Chan.npph, fb_ind,fb_j,J, 2);
    Map_1_to_33(this->map33_ind,this->map33_fac1,this->L3_3_index,hhpp, Chan.hhp_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 1);
    Map_1_to_33(this->map33_ind,this->map33_fac2,this->L3_3_index,hhpp, Chan.hhp_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 2);
    Map_1_to_34(this->map34_ind,this->map34_fac1,this->L3_4_index,hhpp, Chan.hhp_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 1);
    Map_1_to_34(this->map34_ind,this->map34_fac2,this->L3_4_index,hhpp, Chan.hhp_map,Chan.p_map,Chan.np, fb_ind,fb_j,J, 2);

    i = fb_ind[hhpp].v1;
    j = fb_ind[hhpp].v2;
    a = fb_ind[hhpp].v3;
    b = fb_ind[hhpp].v4;
    this->Evec_ind[4*hhpp] = Chan.hh1_map[Chan.ind0][Hash(i, i, 0)];
    this->Evec_fac[4*hhpp] = 1.0 / std::sqrt(fb_j[hhpp].v3 + 1.0);
    this->Evec_ind[4*hhpp + 1] = Chan.hh1_map[Chan.ind0][Hash(j, j, 0)];
    this->Evec_fac[4*hhpp + 1] = 1.0 / std::sqrt(fb_j[hhpp].v4 + 1.0);
    this->Evec_ind[4*hhpp + 2] = Chan.pp1_map[Chan.ind0][Hash(a, a, 0)];
    this->Evec_fac[4*hhpp + 2] = 1.0 / std::sqrt(fb_j[hhpp].v1 + 1.0);
    this->Evec_ind[4*hhpp + 3] = Chan.pp1_map[Chan.ind0][Hash(b, b, 0)];
    this->Evec_fac[4*hhpp + 3] = 1.0 / std::sqrt(fb_j[hhpp].v2 + 1.0);
  }
  
  delete[] fb_ind;
  delete[] fb_j;
  delete[] J;
  delete[] J21_map;
  delete[] J22_map;
  delete[] J23_map;
  delete[] J24_map;
}

void LAmplitudes::Build(Channels &Chan, LAmplitudes &LAmps)
{
  int chan1, chan2, chan3, ind, length;
  /////////////////// l2 ///////////////////////////////////
  length = Chan.nhp1[Chan.ind0];
  this->l2 = new double[length];
  this->evec_ind = new int[2 * length];
  this->evec_fac = new double[2 * length];
  for(ind = 0; ind < length; ++ind){
    this->l2[ind] = LAmps.l2[ind];
  }
  this->l2_length = length;

  this->l3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->l3_index[chan3] = LAmps.l3_index[chan3];
  }
  this->l3_length = LAmps.l3_length;
  this->map3_ind = new int[LAmps.l3_length];
  this->map3_fac1 = new double[LAmps.l3_length];
  this->map3_fac2 = new double[LAmps.l3_length];
  for(ind = 0; ind < LAmps.l3_length; ++ind){
    this->map3_ind[ind] = LAmps.map3_ind[ind];
    this->map3_fac1[ind] = LAmps.map3_fac1[ind];
    this->map3_fac2[ind] = LAmps.map3_fac2[ind];
  }
  for(ind = 0; ind < 2 * LAmps.l2_length; ++ind){
    this->evec_ind[ind] = LAmps.evec_ind[ind];
    this->evec_fac[ind] = LAmps.evec_fac[ind];
  }

  /////////////////// L1 ///////////////////////////////////
  this->L1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){ this->L1_index[chan1] = LAmps.L1_index[chan1]; }
  this->L1 = new double[LAmps.L1_length];
  this->Evec_ind = new int[4 * LAmps.L1_length];
  this->Evec_fac = new double[4 * LAmps.L1_length];
  for(ind = 0; ind < LAmps.L1_length; ++ind){ this->L1[ind] = LAmps.L1[ind]; }
  for(ind = 0; ind < 4 * LAmps.L1_length; ++ind){
    this->Evec_ind[ind] = LAmps.Evec_ind[ind];
    this->Evec_fac[ind] = LAmps.Evec_fac[ind];
  }
  this->L1_length = LAmps.L1_length;

  this->L2_1_index = new int[Chan.size2];
  this->L2_2_index = new int[Chan.size2];
  this->L2_3_index = new int[Chan.size2];
  this->L2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->L2_1_index[chan2] = LAmps.L2_1_index[chan2];
    this->L2_2_index[chan2] = LAmps.L2_2_index[chan2];
    this->L2_3_index[chan2] = LAmps.L2_3_index[chan2];
    this->L2_4_index[chan2] = LAmps.L2_4_index[chan2];
  }
  this->L2_1_length = LAmps.L2_1_length;
  this->L2_2_length = LAmps.L2_2_length;
  this->L2_3_length = LAmps.L2_3_length;
  this->L2_4_length = LAmps.L2_4_length;

  this->L3_1_index = new int[Chan.size3];
  this->L3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->L3_1_index[chan3] = LAmps.L3_1_index[chan3];
    this->L3_2_index[chan3] = LAmps.L3_2_index[chan3];
  }
  this->L3_1_length = LAmps.L3_1_length;
  this->L3_2_length = LAmps.L3_2_length;

  this->L3_3_index = new int[Chan.size3];
  this->L3_4_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->L3_3_index[chan3] = LAmps.L3_3_index[chan3];
    this->L3_4_index[chan3] = LAmps.L3_4_index[chan3];
  }
  this->L3_3_length = LAmps.L3_3_length;
  this->L3_4_length = LAmps.L3_4_length;

  this->map21_num = new int[LAmps.L2_1_length];
  this->map21_index = new int[LAmps.L2_1_length];
  this->map22_num = new int[LAmps.L2_2_length];
  this->map22_index = new int[LAmps.L2_2_length];
  this->map23_num = new int[LAmps.L2_3_length];
  this->map23_index = new int[LAmps.L2_3_length];
  this->map24_num = new int[LAmps.L2_4_length];
  this->map24_index = new int[LAmps.L2_4_length];

  length = 0;
  for(ind = 0; ind < LAmps.L2_1_length; ++ind){
    this->map21_num[ind] = LAmps.map21_num[ind];
    this->map21_index[ind] = LAmps.map21_index[ind];
    length += LAmps.map21_num[ind];
  }
  this->map21_ind = new int[length];
  this->map21_fac1 = new double[length];
  this->map21_fac2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map21_ind[ind] = LAmps.map21_ind[ind];
    this->map21_fac1[ind] = LAmps.map21_fac1[ind];
    this->map21_fac2[ind] = LAmps.map21_fac2[ind];
  }

  length = 0;
  for(ind = 0; ind < LAmps.L2_2_length; ++ind){
    this->map22_num[ind] = LAmps.map22_num[ind];
    this->map22_index[ind] = LAmps.map22_index[ind];
    length += LAmps.map22_num[ind];
  }
  this->map22_ind = new int[length];
  this->map22_fac1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map22_ind[ind] = LAmps.map22_ind[ind];
    this->map22_fac1[ind] = LAmps.map22_fac1[ind];
  }

  length = 0;
  for(ind = 0; ind < LAmps.L2_3_length; ++ind){
    this->map23_num[ind] = LAmps.map23_num[ind];
    this->map23_index[ind] = LAmps.map23_index[ind];
    length += LAmps.map23_num[ind];
  }
  this->map23_ind = new int[length];
  this->map23_fac1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map23_ind[ind] = LAmps.map23_ind[ind];
    this->map23_fac1[ind] = LAmps.map23_fac1[ind];
  }

  length = 0;
  for(ind = 0; ind < LAmps.L2_4_length; ++ind){
    this->map24_num[ind] = LAmps.map24_num[ind];
    this->map24_index[ind] = LAmps.map24_index[ind];
    length += LAmps.map24_num[ind];
  }
  this->map24_ind = new int[length];
  this->map24_fac1 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->map24_ind[ind] = LAmps.map24_ind[ind];
    this->map24_fac1[ind] = LAmps.map24_fac1[ind];
  }

  this->map31_ind = new int[LAmps.L3_1_length];
  this->map31_fac1 = new double[LAmps.L3_1_length];
  this->map31_fac2 = new double[LAmps.L3_1_length];
  for(ind = 0; ind < LAmps.L3_1_length; ++ind){
    this->map31_ind[ind] = LAmps.map31_ind[ind];
    this->map31_fac1[ind] = LAmps.map31_fac1[ind];
    this->map31_fac2[ind] = LAmps.map31_fac2[ind];
  }

  this->map32_ind = new int[LAmps.L3_2_length];
  this->map32_fac1 = new double[LAmps.L3_2_length];
  this->map32_fac2 = new double[LAmps.L3_2_length];
  for(ind = 0; ind < LAmps.L3_2_length; ++ind){
    this->map32_ind[ind] = LAmps.map32_ind[ind];
    this->map32_fac1[ind] = LAmps.map32_fac1[ind];
    this->map32_fac2[ind] = LAmps.map32_fac2[ind];
  }

  this->map33_ind = new int[LAmps.L3_3_length];
  this->map33_fac1 = new double[LAmps.L3_3_length];
  this->map33_fac2 = new double[LAmps.L3_3_length];
  for(ind = 0; ind < LAmps.L3_3_length; ++ind){
    this->map33_ind[ind] = LAmps.map33_ind[ind];
    this->map33_fac1[ind] = LAmps.map33_fac1[ind];
    this->map33_fac2[ind] = LAmps.map33_fac2[ind];
  }

  this->map34_ind = new int[LAmps.L3_4_length];
  this->map34_fac1 = new double[LAmps.L3_4_length];
  this->map34_fac2 = new double[LAmps.L3_4_length];
  for(ind = 0; ind < LAmps.L3_4_length; ++ind){
    this->map34_ind[ind] = LAmps.map34_ind[ind];
    this->map34_fac1[ind] = LAmps.map34_fac1[ind];
    this->map34_fac2[ind] = LAmps.map34_fac2[ind];
  }
}

void LAmplitudes::Delete()
{
  /////////////////// l2 ///////////////////////////////////
  delete[] this->l2;
  delete[] this->evec_ind;
  delete[] this->evec_fac;
  delete[] this->l3_index;

  delete[] this->map3_ind;
  delete[] this->map3_fac1;
  delete[] this->map3_fac2;

  /////////////////// L1 ///////////////////////////////////
  delete[] this->L1;
  delete[] this->Evec_ind;
  delete[] this->Evec_fac;

  delete[] this->L1_index;
  delete[] this->L2_1_index;
  delete[] this->L2_2_index;
  delete[] this->L2_3_index;
  delete[] this->L2_4_index;
  delete[] this->L3_1_index;
  delete[] this->L3_2_index;
  delete[] this->L3_3_index;
  delete[] this->L3_4_index;

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

void LAmplitudes::Zero()
{
  for(int ind = 0; ind < this->L1_length; ++ind){ this->L1[ind] = 0.0; }
  for(int ind = 0; ind < this->l2_length; ++ind){ this->l2[ind] = 0.0; }
}

void LAmplitudes::Set_T(double *vec)
{
  for(int ind = 0; ind < this->L1_length; ++ind){ L1[ind] = vec[ind]; }
  for(int ind = 0; ind < this->l2_length; ++ind){ l2[ind] = vec[ind + this->L1_length]; }
}

void LAmplitudes::Set_L2_1(double *L2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    L2[ind] = 0.0;
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      L2[ind] += this->map21_fac2[index] * this->L1[index1];
    }
  }
}

void LAmplitudes::Gather_L2_1(double *L2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map21_num[ind + offset]; ++n){
      index = this->map21_index[ind + offset] + n;
      index1 = this->map21_ind[index];
      this->L1[index1] += this->map21_fac1[index] * L2[ind];
    }
    L2[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L2_2(double *L2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map22_num[ind + offset]; ++n){
      index = this->map22_index[ind + offset] + n;
      index1 = this->map22_ind[index];
      this->L1[index1] += this->map22_fac1[index] * L2[ind];
    }
    L2[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L2_3(double *L2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map23_num[ind + offset]; ++n){
      index = this->map23_index[ind + offset] + n;
      index1 = this->map23_ind[index];
      this->L1[index1] += this->map23_fac1[index] * L2[ind];
    }
    L2[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L2_4(double *L2, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->map24_num[ind + offset]; ++n){
      index = this->map24_index[ind + offset] + n;
      index1 =this->map24_ind[index];
      this->L1[index1] += this->map24_fac1[index] * L2[ind];
    }
    L2[ind] = 0.0;
  }
}

void LAmplitudes::Set_L3_1(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    L3[ind] = this->map31_fac2[ind + offset] * this->L1[index];
  }
}

void LAmplitudes::Set_L3_2(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    L3[ind] = this->map32_fac2[ind + offset] * this->L1[index];
  }
}

void LAmplitudes::Set_L3_3(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    L3[ind] = this->map33_fac2[ind + offset] * this->L1[index];
  }
}

void LAmplitudes::Set_L3_4(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    L3[ind] = this->map34_fac2[ind + offset] * this->L1[index];
  }
}

void LAmplitudes::Gather_L3_1(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map31_ind[ind + offset];
    this->L1[index] += this->map31_fac1[ind + offset] * L3[ind];
    L3[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L3_2(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map32_ind[ind + offset];
    this->L1[index] += this->map32_fac1[ind + offset] * L3[ind];
    L3[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L3_3(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map33_ind[ind + offset];
    this->L1[index] += this->map33_fac1[ind + offset] * L3[ind];
    L3[ind] = 0.0;
  }
}

void LAmplitudes::Gather_L3_4(double *L3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map34_ind[ind + offset];
    this->L1[index] += this->map34_fac1[ind + offset] * L3[ind];
    L3[ind] = 0.0;
  }
}

void LAmplitudes::Gather_l3(double *l3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    this->l2[index] += this->map3_fac1[ind + offset] * l3[ind];
    l3[ind] = 0.0;
  }
}

void LAmplitudes::Set_l3(double *l3, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->map3_ind[ind + offset];
    l3[ind] = this->map3_fac2[ind + offset] * this->l2[index];
  }
}

void LCC_Algorithm(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps)
{
  double mixmax = 0.75;
  double mixmin = 0.01;
  double mix;
  double width = 1.0;
  int Rand_count = 0;

  double error = 1e20, error2 = 1e20;
  int ind = 0;
  LAmplitudes LAmps0;
  LAmplitudes LAmps2;
  LAmplitudes tempLAmps;
  LAmps0.Build(Chan, LAmps);
  LAmps2.Build(Chan, LAmps);
  tempLAmps.Build(Chan, LAmps);

  /// For DIIS ///
  DIIS DIIS;
  int maxl = 5;
  double DIIS_start = 0.01;
  int size = LAmps.L1_length;
  if(PAR.approx == "singles"){ size += LAmps.l2_length; }
  double *vec = new double[size];
  double *vec0 = new double[size];
  DIIS.Build(size, maxl);

  // Initialize LAmplitudes //
  mix = mixmin;
  while(error > 1e-12 && ind < 2000){
    Update_LAmps(Chan, Ints, Eff_Ints, Amps, LAmps, tempLAmps);

    //Print_LAmps(Chan, tempLAmps);
    LAmps2.Copy(LAmps);
    Gather_LAmps(Chan, Eff_Ints, LAmps2, tempLAmps, mix);

    LCC_Error(Ints, LAmps, LAmps2, error);
    error /= mix;
    if( !std::isfinite(error) ){ std::cerr << std::endl << ind << ", error = " << error << " : LCC Solution Diverged!! " << std::endl; exit(1); }
    if( ind > 1000 && double(Rand_count)/double(ind) > 0.9 ){ std::cerr << std::endl << ind << " : LCC Solution Not Converged!!" << std::endl; exit(1); }

    if( error < error2 || error > 1.0 ){
      if( error < error2 ){ mix = std::pow(mixmax, 0.035) * std::pow(mix, 0.965); }
      LAmps0.Copy(LAmps, vec0);
      LAmps.Copy(LAmps2, vec);
      error2 = error;
      if( error < DIIS_start && mix > 0.2 * mixmax ){
	DIIS.Perform(vec, vec0, mix);
	LAmps.Zero();
	LAmps.Set_T(vec);
      }
    }
    else{
      ++Rand_count;
      if(error2 > 1.0){ width = 0.001; }
      else{ width = 0.001 * error2; }
      Random_LStep(Chan, Ints, Eff_Ints, Amps, LAmps0, LAmps, LAmps2, tempLAmps, mix, width, error, error2);
    }

    //CCoutE = LAmps.get_energy(Chan, Ints);
    //if( !std::isfinite(CCoutE) ){ std::cerr << std::endl << ind << ", Energy = " << CCoutE << " : CC Solution Diverged!! " << std::endl; exit(1); }
    //std::cout << "Iteration Number = " << ind << ", Energy = " << CCoutE << ", error = " << error << ", mix = " << mix << ", ";
    //std::cout << "Random = " << Rand_count << ", DIIS = " << DIIS.count << std::endl;
    std::cout << "Lambda Iteration Number = " << ind << ", error = " << error << ", mix = " << mix << ", ";
    std::cout << "Random = " << Rand_count << ", DIIS = " << DIIS.count << std::endl;
    ++ind;
  }

  std::cout << std::endl << std::endl;
  if( error > 1e-12 ){
    std::cout << PAR.Shells << ", " << PAR.Pshells << ", " << PAR.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". LCC Solution Not Converged!!" << std::endl;
  }
  LAmps0.Delete();
  LAmps2.Delete();
  tempLAmps.Delete();
  DIIS.Delete();
  delete[] vec;
  delete[] vec0;
}

void Update_LAmps(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps, LAmplitudes &tempLAmps)
{
  tempLAmps.Zero();
  //  Solve LAmplitude Equations
  LSingles_Step(Chan, Ints, Eff_Ints, Amps, LAmps, tempLAmps);
  LDoubles_Step(Chan, Ints, Eff_Ints, Amps, LAmps, tempLAmps);
}

void Gather_LAmps(Channels &Chan, Eff_Interactions &Eff_Ints, LAmplitudes &LAmps, LAmplitudes &LAmps2, double &mix)
{
  double tempt, tempen;
  int length = LAmps.L1_length;
  for(int ind = 0; ind < length; ++ind){
    tempt = LAmps2.L1[ind];
    tempen = 0.0;
    tempen += Eff_Ints.Xhh.X_2[LAmps.Evec_ind[4*ind]] * LAmps.Evec_fac[4*ind];
    tempen += Eff_Ints.Xhh.X_2[LAmps.Evec_ind[4*ind + 1]] * LAmps.Evec_fac[4*ind + 1];
    tempen -= Eff_Ints.Xpp.X_2[LAmps.Evec_ind[4*ind + 2]] * LAmps.Evec_fac[4*ind + 2];
    tempen -= Eff_Ints.Xpp.X_2[LAmps.Evec_ind[4*ind + 3]] * LAmps.Evec_fac[4*ind + 3];
    tempt /= tempen;
    tempt = mix*tempt + (1.0-mix)*LAmps.L1[ind];
    LAmps.L1[ind] = tempt;
  }
  if(PAR.approx == "singles"){
    length = LAmps.l2_length;
    for(int ind = 0; ind < length; ++ind){
      tempt = LAmps2.l2[ind];
      tempen = 0.0;
      tempen += Eff_Ints.Xhh.X1_2[LAmps.evec_ind[2*ind]] * LAmps.evec_fac[2*ind];
      tempen -= Eff_Ints.Xpp.X_2[LAmps.evec_ind[2*ind + 1]] * LAmps.evec_fac[2*ind + 1];
      tempt /= tempen;
      tempt = mix*tempt + (1.0-mix)*LAmps.l2[ind];
      LAmps.l2[ind] = tempt;
    }
  }
}

void LCC_Error(Interactions &Ints, LAmplitudes &LAmps, LAmplitudes &LAmps2, double &error)
{
  error = 0.0;
  int ind;
  int length = LAmps.L1_length;
  double norm = 1.0e-16;
  for(ind = 0; ind < length; ++ind){
    error += (LAmps2.L1[ind] - LAmps.L1[ind])*(LAmps2.L1[ind] - LAmps.L1[ind]);
    norm += LAmps2.L1[ind] * LAmps2.L1[ind];
  }
  if(PAR.approx == "singles"){  
    length = LAmps.l2_length;
    for(ind = 0; ind < length; ++ind){
      error += (LAmps2.l2[ind] - LAmps.l2[ind])*(LAmps2.l2[ind] - LAmps.l2[ind]);
      norm += LAmps2.l2[ind] * LAmps2.l2[ind];
    }
  }
  error = std::sqrt(error/norm);
}

void Print_LAmps(Channels &Chan, LAmplitudes &LAmps)
{
  int length1 = Chan.size1;
  int length2 = Chan.nph1[Chan.ind0];
  int npp, nhh;
  int chanind, ind;//, index;
  int a, b, i, j;
  two_body pp, hh, hp;
  ind = 0;

  std::cout << std::setprecision(5) << std::endl;
  for(int chan1 = 0; chan1 < length1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    chanind = LAmps.L1_index[chan1];
    for(int hh1 = 0; hh1 < nhh; ++hh1){
      hh = Chan.hh_state(chan1, hh1);
      i = hh.v1;
      j = hh.v2;
      for(int pp1 = 0; pp1 < npp; ++pp1){
	pp = Chan.pp_state(chan1, pp1);
	a = pp.v1;
	b = pp.v2;
	std::cout << "! < " << i << "," << j << " |l| " << a << "," << b << " >^ " << Chan.qnums1[chan1].j << " = " << chanind << " " << (hh1*npp + pp1) << ", " << LAmps.L1[chanind + (hh1*npp + pp1)] << std::endl;
	++ind;
      }
    }
  }
  if(PAR.approx == "singles"){
    for(int hp1 = 0; hp1 < length2; ++hp1){
      hp = Chan.hp1_state(Chan.ind0, hp1);
      i = hp.v1;
      a = hp.v2;
      std::cout << "! < " << i << " |l| " << a << " > = " << LAmps.l2[hp1] << std::endl;
    }
  }
  std::cout << std::endl;
}

void Random_LStep(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps0, LAmplitudes &LAmps, LAmplitudes &LAmps2, LAmplitudes &tempLAmps, double &mix, double &width, double &error, double &error2)
{
  int go = 0;
  int ind = 0;
  while(go == 0){
    ++ind;
    LAmps.Zero();
    Randomize_LAmps(LAmps0, LAmps, width);

    tempLAmps.Zero();
    //  Solve LAmplitude Equations  //
    LSingles_Step(Chan, Ints, Eff_Ints, Amps, LAmps, tempLAmps);
    LDoubles_Step(Chan, Ints, Eff_Ints, Amps, LAmps, tempLAmps);
    //////////////////////////////////
    LAmps2.Copy(LAmps);
    Gather_LAmps(Chan, Eff_Ints, LAmps2, tempLAmps, mix);
    LCC_Error(Ints, LAmps, LAmps2, error);
    error /= mix;

    std::cout << "!!  " << ind << ", " << error << " " << error2 << ", " << width << std::endl;
    if(error < error2){
      ++go;
      LAmps0.Copy(LAmps);
      LAmps.Copy(LAmps2);
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

void Randomize_LAmps(LAmplitudes &LAmps0, LAmplitudes &LAmps, double &width)
{
  int ind;
  double tempt;
  double rand;
  size_t key;
  std::unordered_map<size_t,double> t_map;
  int length = LAmps.L1_length;
  for(ind = 0; ind < length; ++ind){
    tempt = LAmps0.L1[ind];
    if(fabs(tempt) > 1.0e-12){
      rand = rand_normal(0.0, width * fabs(tempt));
      key = std::hash<float>{}(float(fabs(tempt)));
      t_map[key] = rand;
    }
  }
  if(PAR.approx == "singles"){  
    length = LAmps.l2_length;
    for(ind = 0; ind < length; ++ind){
      tempt = LAmps0.l2[ind];
      if(fabs(tempt) > 1.0e-12){
	rand = rand_normal(0.0, width * fabs(tempt));
	key = std::hash<float>{}(float(fabs(tempt)));
	t_map[key] = rand;
      }
    }
  }

  length = LAmps.L1_length;
  for(ind = 0; ind < length; ++ind){
    tempt = LAmps0.L1[ind];
    key = std::hash<float>{}(float(fabs(tempt)));
    if(tempt > 1.0e-12){ tempt += t_map[key]; }
    else if(tempt < -1.0e-12){ tempt -= t_map[key]; }
    LAmps.L1[ind] = tempt;
  }
  if(PAR.approx == "singles"){  
    length = LAmps.l2_length;
    for(ind = 0; ind < length; ++ind){
      tempt = LAmps0.l2[ind];
      key = std::hash<float>{}(float(fabs(tempt)));
      if(tempt > 1.0e-12){ tempt += t_map[key]; }
      else if(tempt < -1.0e-12){ tempt -= t_map[key]; }
      LAmps.l2[ind] = tempt;
    }
  }
}

void LDoubles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps1, LAmplitudes &LAmps2)
{
  double p1 = 1.0, p12 = 0.5, m12 = -0.5, m1 = -1.0;
  int nh, np, npph, nhhp, nhh, npp, nhp1, nph1, hp2, length;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4, chan_ind5;
  int a, b, i, j;
  double *L2, *LL2, *L3, *LL3, *T3, *X3, *Q3;
  char N = 'N';
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    npp = Chan.npp[chan1];
    length = nhh * npp;
    chan_ind1 = LAmps1.L1_index[chan1];
    chan_ind2 = Ints.Vhhpp.V_1_index[chan1];
    //L1(ij|ab){ij,ab}  =  Vhhpp1(ij|ab){ij,ab}
    for(int ind = 0; ind < length; ++ind){ LAmps2.L1[chan_ind1 + ind] += Ints.Vhhpp.V_1[chan_ind2 + ind]; }

    if(length != 0){
      chan_ind2 = Eff_Ints.Xpppp.X_1_index[chan1];
      chan_ind3 = Eff_Ints.Xhhhh.X_1_index[chan1];
      //L1(ij|ab){ij,ab}  <-  (1/2).L1(ij|cd){ij,cd}.Xpppp1(cd|ab){cd,ab}
      dgemm_NN((LAmps1.L1 + chan_ind1), (Eff_Ints.Xpppp.X_1 + chan_ind2), (LAmps2.L1 + chan_ind1), &nhh, &npp, &npp, &p12, &p1, &N, &N);
      //L1(ij|ab){ij,ab}  <-  (1/2).Xhhhh1(ij|kl){ij,kl}.L1(kl|ab){kl,ab}
      dgemm_NN((Eff_Ints.Xhhhh.X_1 + chan_ind3), (LAmps1.L1 + chan_ind1), (LAmps2.L1 + chan_ind1), &nhh, &npp, &nhh, &p12, &p1, &N, &N);
    }
  }

  ///////////////////////////////////////////////////
  nhp1 = Chan.nhp1[Chan.ind0];
  nph1 = Chan.nph1[Chan.ind0];
  length = nhp1 * nph1;
  chan_ind1 = LAmps1.L2_1_index[Chan.ind0];
  chan_ind2 = LAmps1.L2_2_index[Chan.ind0];
  chan_ind3 = LAmps1.L2_3_index[Chan.ind0];
  chan_ind4 = LAmps1.L2_4_index[Chan.ind0];
  L2 = new double[length];
  for(int ind = 0; ind < length; ++ind){ L2[ind] = 0.0; }
  // L2_1(ij|ab){ib',aj'} = - L3(i|b){ib'}.Xhp3(j|a){ja'}
  for(int hp1 = 0; hp1 < nhp1; ++hp1){
    for(int ph1 = 0; ph1 < nph1; ++ph1){
      a = Chan.ph1_state(Chan.ind0, ph1).v1;
      j = Chan.ph1_state(Chan.ind0, ph1).v2;
      hp2 = Chan.hp1_map[Chan.ind0][Hash(j, a, 0)];
      L2[hp1 * nph1 + ph1] += -1.0 * LAmps1.l2[hp1] * Eff_Ints.Xhp.X_2[hp2];
    }
  }
  LAmps2.Gather_L2_1(L2, nhp1*nph1, chan_ind1);
  // L2_2(ij|ab){ja',bi'} = - L3(j|a){ja'}.Xhp3(i|b){ib'}
  for(int hp1 = 0; hp1 < nhp1; ++hp1){
    for(int ph1 = 0; ph1 < nph1; ++ph1){
      b = Chan.ph1_state(Chan.ind0, ph1).v1;
      i = Chan.ph1_state(Chan.ind0, ph1).v2;
      hp2 = Chan.hp1_map[Chan.ind0][Hash(i, b, 0)];
      L2[hp1 * nph1 + ph1] += -1.0 * LAmps1.l2[hp1] * Eff_Ints.Xhp.X_2[hp2];
    }
  }
  LAmps2.Gather_L2_2(L2, nhp1*nph1, chan_ind2);
  // L2_3(ij|ab){ia',bj'} = + L3(i|a){ia'}.Xhp3(j|b){jb'}
  for(int hp1 = 0; hp1 < nhp1; ++hp1){
    for(int ph1 = 0; ph1 < nph1; ++ph1){
      b = Chan.ph1_state(Chan.ind0, ph1).v1;
      j = Chan.ph1_state(Chan.ind0, ph1).v2;
      hp2 = Chan.hp1_map[Chan.ind0][Hash(j, b, 0)];
      L2[hp1 * nph1 + ph1] += LAmps1.l2[hp1] * Eff_Ints.Xhp.X_2[hp2];
    }
  }
  LAmps2.Gather_L2_3(L2, nhp1*nph1, chan_ind3);
  // L2_4(ij|ab){jb',ai'} = + L3(j|b){jb'}.Xhp3(i|a){ia'}
  for(int hp1 = 0; hp1 < nhp1; ++hp1){
    for(int ph1 = 0; ph1 < nph1; ++ph1){
      a = Chan.ph1_state(Chan.ind0, ph1).v1;
      i = Chan.ph1_state(Chan.ind0, ph1).v2;
      hp2 = Chan.hp1_map[Chan.ind0][Hash(i, a, 0)];
      L2[hp1 * nph1 + ph1] += LAmps1.l2[hp1] * Eff_Ints.Xhp.X_2[hp2];
    }
  }
  LAmps2.Gather_L2_4(L2, nhp1*nph1, chan_ind4);
  ///////////////////////////////////////////////////

  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    nph1 = Chan.nph1[chan2];
    length = nhp1 * nph1;
    if(length != 0){
      chan_ind1 = LAmps1.L2_1_index[chan2];
      chan_ind2 = LAmps1.L2_2_index[chan2];
      chan_ind3 = LAmps1.L2_3_index[chan2];
      chan_ind4 = LAmps1.L2_4_index[chan2];
      chan_ind5 = Eff_Ints.Xhphp.X_2_1_index[chan2];

      LL2 = new double[nhp1 * nph1];
      LAmps1.Set_L2_1(LL2, nhp1*nph1, chan_ind1);
      L2 = new double[nhp1 * nph1];
      for(int ind = 0; ind < nhp1 * nph1; ++ind){ L2[ind] = 0.0; }
      //L2_1(ij|ab){ib',aj'}  =  -Xhphp2_1(ic|kb){ib',kc'}.L2_1(kj|ac){kc',aj'}
      dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind5), LL2, L2, &nhp1, &nph1, &nhp1, &m1, &p1, &N, &N);
      LAmps2.Gather_L2_1(L2, nhp1*nph1, chan_ind1);
      //L2_2(ij|ab){ja',bi'}  =  -Xhphp2_1(jc|ka){ja',kc'}.L2_1(ki|bc){kc',bi'}
      dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind5), LL2, L2, &nhp1, &nph1, &nhp1, &m1, &p1, &N, &N);
      LAmps2.Gather_L2_2(L2, nhp1*nph1, chan_ind2);
      //L2_3(ij|ab){ia',bj'}  =  Xhphp2_1(ic|ka){ia',kc'}.L2_1(kj|bc){kc',bj'}
      dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind5), LL2, L2, &nhp1, &nph1, &nhp1, &p1, &p1, &N, &N);
      LAmps2.Gather_L2_3(L2, nhp1*nph1, chan_ind3);
      //L2_4(ij|ab){jb',ai'}  =  Xhphp2_1(jc|kb){jb',kc'}.L2_1(ki|ac){kc',ai'}
      dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind5), LL2, L2, &nhp1, &nph1, &nhp1, &p1, &p1, &N, &N);
      LAmps2.Gather_L2_4(L2, nhp1*nph1, chan_ind4);

      delete[] LL2;
      delete[] L2;
    }
  }

  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];
    length = nh * npph;
    if(length != 0){
      chan_ind1 = LAmps1.L3_1_index[chan3];
      chan_ind2 = LAmps1.L3_2_index[chan3];
      chan_ind3 = Eff_Ints.Xhh.X_3_index[chan3];
      chan_ind4 = Amps.D1.T3_3_index[chan3];
      chan_ind5 = Ints.Vhhpp.V_3_1_index[chan3];

      Q3 = new double[nh * nh];
      for(int ind = 0; ind < nh*nh; ++ind){ Q3[ind] = 0.0; }
      T3 = new double[npph * nh];
      Amps.D1.Set_T3_3(T3, npph*nh, chan_ind4);
      X3 = new double[nh * nh];
      Eff_Ints.Xhh.Set_X_3_OD(X3, nh*nh, chan_ind3);
      LL3 = new double[nh * npph];
      LAmps1.Set_L3_1(LL3, nh*npph, chan_ind1);
      //Q_3(i|k){i,k}  =  L3_1(il|cd){i,cdl'}.T3_3(cd|kl){cdl',k}
      dgemm_NN(LL3, T3, Q3, &nh, &nh, &npph, &p1, &p1, &N, &N);

      L3 = new double[nh * npph];
      for(int ind = 0; ind < nh*npph; ++ind){ L3[ind] = 0.0; }
      //L3_1(ij|ab){i,abj'}  =  -Xhh3_od(i|k){i,k}.L3_1(kj|ab){k,abj'}
      dgemm_NN(X3, LL3, L3, &nh, &npph, &nh, &m1, &p1, &N, &N);
      //L3_1(ij|ab){i,abj'}  =  -(1/2).L3_1(il|cd){i,cdl'}.T3_3(cd|kl){cdl',k}.Vhhpp3_1(kj|ab){k,abj'}
      dgemm_NN(Q3, (Ints.Vhhpp.V_3_1 + chan_ind5), L3, &nh, &npph, &nh, &m12, &p1, &N, &N);
      LAmps2.Gather_L3_1(L3, nh*npph, chan_ind1);

      //L3_2(ij|ab){j,abi'}  =  Xhh3_od(j|k){j,k}.L3_1(ki|ab){k,abi'}
      dgemm_NN(X3, LL3, L3, &nh, &npph, &nh, &p1, &p1, &N, &N);
      //L3_2(ij|ab){j,abi'}  =  (1/2).L3_1(jl|cd){j,cdl'}.T3_3(cd|kl){cdl',k}.Vhhpp3_1(ki|ab){k,abi'}
      dgemm_NN(Q3, (Ints.Vhhpp.V_3_1 + chan_ind5), L3, &nh, &npph, &nh, &p12, &p1, &N, &N);
      LAmps2.Gather_L3_2(L3, nh*npph, chan_ind2);

      delete[] Q3;
      delete[] X3;
      delete[] LL3;
      delete[] L3;
    }
    length = np * nhhp;
    if(length != 0){
      chan_ind1 = LAmps1.L3_3_index[chan3];
      chan_ind2 = LAmps1.L3_4_index[chan3];
      chan_ind3 = Eff_Ints.Xpp.X_3_index[chan3];
      chan_ind4 = Amps.D1.T3_1_index[chan3];
      chan_ind5 = Ints.Vhhpp.V_3_3_index[chan3];

      Q3 = new double[np * np];
      for(int ind = 0; ind < np*np; ++ind){ Q3[ind] = 0.0; }
      T3 = new double[np * nhhp];
      Amps.D1.Set_T3_1(T3, np*nhhp, chan_ind4);
      X3 = new double[np * np];
      Eff_Ints.Xpp.Set_X_3_OD(X3, np*np, chan_ind3);
      LL3 = new double[nhhp * np];
      LAmps1.Set_L3_3(LL3, nhhp*np, chan_ind1);
      //Q_3(c|a){c,a}  =  T3_1(cd|kl){c,kld'}.L3_3(kl|ad){kld',a}
      dgemm_NN(T3, LL3, Q3, &np, &np, &nhhp, &p1, &p1, &N, &N);

      L3 = new double[nhhp * np];
      for(int ind = 0; ind < nhhp*np; ++ind){ L3[ind] = 0.0; }
      //L3_3(ij|ab){ijb',a}  =  L3_3(ij|ca){ija',c}.Xpp3_od(c|a){c,a}
      dgemm_NN(LL3, X3, L3, &nhhp, &np, &np, &p1, &p1, &N, &N);
      //L3_3(ij|ab){ijb',a}  =  -(1/2).Vhhpp3_3(ij|cb){ijb',c}.T3_1(cd|kl){c,kld'}.L3_3(kl|ad){kld',a}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind5), Q3, L3, &nhhp, &np, &np, &m12, &p1, &N, &N);
      LAmps2.Gather_L3_3(L3, nhhp*np, chan_ind1);

      //L3_4(ij|ab){ija',b}  =  -L3_3(ij|cb){ijb',c}.Xpp3_od(c|b){c,b}
      dgemm_NN(LL3, X3, L3, &nhhp, &np, &np, &m1, &p1, &N, &N);
      //L3_3(ij|ab){ija',b}  =  (1/2).Vhhpp3_3(ij|ca){ija',c}.T3_1(cd|kl){c,kld'}.L3_3(kl|bd){kld',b}
      dgemm_NN((Ints.Vhhpp.V_3_3 + chan_ind5), Q3, L3, &nhhp, &np, &np, &p12, &p1, &N, &N);
      LAmps2.Gather_L3_4(L3, nhhp*np, chan_ind2);

      delete[] Q3;
      delete[] X3;
      delete[] LL3;
      delete[] L3;
    }
  }
}

void LSingles_Step(Channels &Chan, Interactions &Ints, Eff_Interactions &Eff_Ints, Amplitudes &Amps, LAmplitudes &LAmps1, LAmplitudes &LAmps2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int nh, np, npph, nhhp, one = 1;
  int nhp0 = Chan.nhp1[Chan.ind0];
  int nhh0 = Chan.nhh1[Chan.ind0];
  int npp0 = Chan.npp1[Chan.ind0];
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  int d, c, j, k, pp1, hh1;
  double *l3_1, *l3_2, *L3, *T3, *Xpp3, *Xhh3, *Qhh, *Qpp, *Q3, *Xhphh, *Xpphp, *Xhhhp, *Xhppp;
  char N = 'N';

  Qhh = new double[nhh0];
  for(int hh = 0; hh < nhh0; ++hh){ Qhh[hh] = 0.0; }
  Qpp = new double[npp0];
  for(int pp = 0; pp < npp0; ++pp){ Qpp[pp] = 0.0; }

  // l2(i|a){ia'}  <-  Xhp(i|a){ia'}
  for(int hp = 0; hp < nhp0; ++hp){ LAmps2.l2[hp] += Eff_Ints.Xhp.X_2[hp]; }

  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
    nh = Chan.nh[chan3];
    np = Chan.np[chan3];
    npph = Chan.npph[chan3];
    nhhp = Chan.nhhp[chan3];

    chan_ind1 = LAmps1.l3_index[chan3];
    l3_1 = new double[nh * np];
    LAmps1.Set_l3(l3_1, nh*np, chan_ind1);
    l3_2 = new double[nh * np];
    for(int ind = 0; ind < nh*np; ++ind){ l3_2[ind] = 0.0; }

    if(nh * np != 0){
      chan_ind2 = Eff_Ints.Xpp.X_3_index[chan3];
      Xpp3 = new double[np * np];
      Eff_Ints.Xpp.Set_X_3_OD(Xpp3, np*np, chan_ind2);
      //l3(i|a){i,a}  <-  l3(i|c){i,c}.Xpp3_od(c|a){c,a}
      dgemm_NN(l3_1, Xpp3, l3_2, &nh, &np, &np, &p1, &p1, &N, &N);

      chan_ind2 = Eff_Ints.Xhh.X_3_index[chan3];
      Xhh3 = new double[nh * nh];
      Eff_Ints.Xhh.Set_X_3_OD(Xhh3, nh*nh, chan_ind2);
      //l3(i|a){i,a}  <-  - Xhh3_od(i|k){i,k}.l3(k|a){k,a}
      dgemm_NN(Xhh3, l3_1, l3_2, &nh, &np, &nh, &m1, &p1, &N, &N);

      delete[] Xpp3;
      delete[] Xhh3;
    }
    if(nhhp * np * nh != 0){
      chan_ind2 = Eff_Ints.Xhphh.X_3_2_index[chan3];
      chan_ind3 = LAmps1.L3_3_index[chan3];
      chan_ind4 = Amps.D1.T3_2_index[chan3];

      Xhphh = new double[nh * nhhp];
      Eff_Ints.Xhphh.Set_X_3_1(Xhphh, nh*nhhp, chan_ind2);
      L3 = new double[nhhp * np];
      LAmps1.Set_L3_3(L3, nhhp*np, chan_ind3);
      //l3(i|a){i,a}  <-  -(1/2) Xhphh3_1(ic|kl){i,klc'}.L3_3(kl|ac){klc',a}
      dgemm_NN(Xhphh, L3, l3_2, &nh, &np, &nhhp, &m12, &p1, &N, &N);
      
      T3 = new double[np * nhhp];
      Amps.D1.Set_T3_2(T3, np*nhhp, chan_ind4);
      Q3 = new double[np * np];
      for(int ind = 0; ind < np*np; ++ind){ Q3[ind] = 0.0; }
      //Q3(d|c){d,c}  <-  T3_2(bd|kl){d,klb'}.L3_3(kl|cb){klb',c}
      dgemm_NN(T3, L3, Q3, &np, &np, &nhhp, &p1, &p1, &N, &N);
      for(int p1 = 0; p1 < np; ++p1){
	d = Chan.p_state(chan3, p1).v1;
	for(int p2 = 0; p2 < np; ++p2){
	  c = Chan.p_state(chan3, p2).v1;
	  pp1 = Chan.pp1_map[Chan.ind0][Hash(d, c, 0)];
	  Qpp[pp1] = Q3[p1 * np + p2] * std::sqrt(SPB.qnums[d].j + 1.0);
	}
      }

      delete[] T3;
      delete[] Q3;
      delete[] Xhphh;
      delete[] L3;
    }
    if(npph * np * nh != 0){
      chan_ind2 = LAmps1.L3_1_index[chan3];
      chan_ind3 = Eff_Ints.Xpphp.X_3_4_index[chan3];
      chan_ind4 = Amps.D1.T3_3_index[chan3];

      L3 = new double[nh * npph];
      LAmps1.Set_L3_1(L3, nh*npph, chan_ind2);
      Xpphp = new double[npph * np];
      Eff_Ints.Xpphp.Set_X_3_4(Xpphp, npph*np, chan_ind3);
      //l3(i|a){i,a}  <-  -(1/2) L3_1(ik|cd){i,cdk'}.Xpphp3_4(cd|ka){cdk',a}
      dgemm_NN(L3, Xpphp, l3_2, &nh, &np, &npph, &m12, &p1, &N, &N);

      T3 = new double[npph * nh];
      Amps.D1.Set_T3_3(T3, npph*nh, chan_ind4);
      Q3 = new double[nh * nh];
      for(int ind = 0; ind < nh*nh; ++ind){ Q3[ind] = 0.0; }
      //Q3(j|k){j,k}  <-  L3_1(jl|cd){j,cdl'}.T3_3(cd|kl){cdl',k}
      dgemm_NN(L3, T3, Q3, &nh, &nh, &npph, &p1, &p1, &N, &N);
      for(int h1 = 0; h1 < nh; ++h1){
	j = Chan.h_state(chan3, h1).v1;
	for(int h2 = 0; h2 < nh; ++h2){
	  k = Chan.h_state(chan3, h2).v1;
	  hh1 = Chan.hh1_map[Chan.ind0][Hash(j, k, 0)];
	  Qhh[hh1] = Q3[h1 * nh + h2] * std::sqrt(SPB.qnums[j].j + 1.0);
	}
      }

      delete[] T3;
      delete[] Q3;
      delete[] L3;
      delete[] Xpphp;
    }
    LAmps2.Gather_l3(l3_2, np*nh, chan_ind1);

    delete[] l3_1;
    delete[] l3_2;
  }
  if(nhp0 != 0){
    chan_ind1 = Eff_Ints.Xhphp.X_2_1_index[Chan.ind0];
    //l2(i|a){ia'}  =  -Xhphp2_1(ic|ka){ia',kc'}.l2(k|c){kc'}
    dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind1), LAmps1.l2, LAmps2.l2, &nhp0, &one, &nhp0, &m1, &p1, &N, &N);
  }
  if(nhp0 * nhh0 != 0){
    chan_ind1 = Eff_Ints.Xhhhp.X_2_1_index[Chan.ind0];

    Xhhhp = new double[nhp0 * nhh0];
    Eff_Ints.Xhhhp.Set_X_2_1(Xhhhp, nhp0*nhh0, chan_ind1);
    //l2(i|a){ia'}  <-  (1/2).Xhhhp2_1(ik|ja){ia',jk'}.Qhh(j|k){jk'}
    dgemm_NN(Xhhhp, Qhh, LAmps2.l2, &nhp0, &one, &nhh0, &p12, &p1, &N, &N);

    delete[] Xhhhp;
  }
  if(nhp0 * npp0 != 0){
    chan_ind1 = Eff_Ints.Xhppp.X_2_3_index[Chan.ind0];

    Xhppp = new double[nhp0 * npp0];
    Eff_Ints.Xhppp.Set_X_2_3(Xhppp, nhp0*npp0, chan_ind1);
    //l2(i|a){ia'}  <-  -(1/2).Xhppp2_3(ic|ad){ia',dc'}.Qpp(d|c){dc'}
    dgemm_NN(Xhppp, Qpp, LAmps2.l2, &nhp0, &one, &npp0, &m12, &p1, &N, &N);

    delete[] Xhppp;
  }
  delete[] Qhh;
  delete[] Qpp;
}
