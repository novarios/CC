#include "INTfunctions.hpp"
#include "CoupledCluster.hpp"
#include "HFfunctions.hpp"
#include "BASISfunctions.hpp"
#include "MATHfunctions.hpp"
#include "AngMom.hpp"

void Interactions::Build(Channels &Chan)
{
  this->Eref = 0.0;
  this->Fmatrix.Build(Chan);
  this->Vhhhh.Build(Chan);
  this->Vpppp.Build(Chan);
  this->Vhhpp.Build(Chan);
  this->Vpphh.Build(Chan);
  this->Vhphp.Build(Chan);
  this->Vhhhp.Build(Chan);
  this->Vhppp.Build(Chan);
  this->Vhphh.Build(Chan);
  this->Vpphp.Build(Chan);
}

void Interactions::Delete()
{
  this->Fmatrix.Delete();
  this->Vhhhh.Delete();
  this->Vpppp.Delete();
  this->Vhhpp.Delete();
  this->Vpphh.Delete();
  this->Vhphp.Delete();
  this->Vhhhp.Delete();
  this->Vhppp.Delete();
  this->Vhphh.Delete();
  this->Vpphp.Delete();
}

void F_matrix::Build(Channels &Chan)
{
  int length, ind, chan3;
  this->hh_2_length = Chan.nhh1[Chan.ind0];
  this->hp_2_length = Chan.nhp1[Chan.ind0];
  this->ph_2_length = Chan.nph1[Chan.ind0];
  this->pp_2_length = Chan.npp1[Chan.ind0];
  this->hh_2 = new double[this->hh_2_length];
  this->hp_2 = new double[this->hp_2_length];
  this->ph_2 = new double[this->ph_2_length];
  this->pp_2 = new double[this->pp_2_length];
  for(ind = 0; ind < this->hh_2_length; ++ind){ this->hh_2[ind] = 0.0; }
  for(ind = 0; ind < this->hp_2_length; ++ind){ this->hp_2[ind] = 0.0; }
  for(ind = 0; ind < this->ph_2_length; ++ind){ this->ph_2[ind] = 0.0; }
  for(ind = 0; ind < this->pp_2_length; ++ind){ this->pp_2[ind] = 0.0; }

  length = 0;
  this->hh_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->hh_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nh[chan3];
  }
  this->hh_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->hh_3[ind] = 0.0; }
  this->hh_3_length = length;

  length = 0;
  this->hp_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->hp_3_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.np[chan3];
  }
  this->hp_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->hp_3[ind] = 0.0; }
  this->hp_3_length = length;

  length = 0;
  this->ph_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->ph_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nh[chan3];
  }
  this->ph_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->ph_3[ind] = 0.0; }
  this->ph_3_length = length;

  length = 0;
  this->pp_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->pp_3_index[chan3] = length;
    length += Chan.np[chan3] * Chan.np[chan3];
  }
  this->pp_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->pp_3[ind] = 0.0; }
  this->pp_3_length = length;
}

void F_matrix::Delete()
{
  delete[] this->hh_2;
  delete[] this->hp_2;
  delete[] this->ph_2;
  delete[] this->pp_2;
  delete[] this->hh_3;
  delete[] this->hp_3;
  delete[] this->ph_3;
  delete[] this->pp_3;
  delete[] this->hh_3_index;
  delete[] this->hp_3_index;
  delete[] this->ph_3_index;
  delete[] this->pp_3_index;
}

void V_hhhh::Build(Channels &Chan)
{
  int chan1, chan3, ind, length;
  length = 0;
  this->V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->V_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.nhh[chan1];
  }
  this->V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_1[ind] = 0.0; }
  this->V_1_length = length;

  length = 0;
  this->V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhhh[chan3];
  }
  this->V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_2[ind] = 0.0; }
  this->V_3_2_length = length;
}

void V_hhhh::Delete()
{
  delete[] this->V_1;
  delete[] this->V_3_2;
  delete[] this->V_1_index;
  delete[] this->V_3_2_index;
}

void V_pppp::Build(Channels &Chan)
{
  int chan1, ind, length;
  length = 0;
  this->V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->V_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.npp[chan1];
  }
  this->V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_1[ind] = 0.0; }
  this->V_1_length = length;
}

void V_pppp::Delete()
{
  delete[] this->V_1;
  delete[] this->V_1_index;
}

void V_hhpp::Build(Channels &Chan)
{
  int chan1, chan2, chan3, ind, length;
  length = 0;
  this->V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->V_1_index[chan1] = length;
    length += Chan.nhh[chan1] * Chan.npp[chan1];
  }
  this->V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_1[ind] = 0.0; }
  this->V_1_length = length;

  length = 0;
  this->V_2_1_index = new int[Chan.size2];
  this->V_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->V_2_1_index[chan2] = length;
    this->V_2_3_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nph1[chan2];
  }
  this->V_2_1 = new double[length];
  this->V_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->V_2_1[ind] = 0.0;
    this->V_2_3[ind] = 0.0;
  }
  this->V_2_1_length = length;
  this->V_2_3_length = length;

  length = 0;
  this->V_3_1_index = new int[Chan.size3];
  this->V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_1_index[chan3] = length;
    this->V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.npph[chan3];
  }
  this->V_3_1 = new double[length];
  this->V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->V_3_1[ind] = 0.0;
    this->V_3_2[ind] = 0.0;
  }
  this->V_3_1_length = length;
  this->V_3_2_length = length;

  length = 0;
  this->V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.np[chan3];
  }
  this->V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_3[ind] = 0.0; }
  this->V_3_3_length = length;
}

void V_hhpp::Delete()
{
  delete[] this->V_1;
  delete[] this->V_2_1;
  delete[] this->V_2_3;
  delete[] this->V_3_1;
  delete[] this->V_3_2;
  delete[] this->V_3_3;
  delete[] this->V_1_index;
  delete[] this->V_2_1_index;
  delete[] this->V_2_3_index;
  delete[] this->V_3_1_index;
  delete[] this->V_3_2_index;
  delete[] this->V_3_3_index;
}

void V_pphh::Build(Channels &Chan)
{
  /*int chan1, ind, length;
  length = 0;
  this->V_1_index = new int[Chan.size1];
  for(chan1 = 0; chan1 < Chan.size1; ++chan1){
    this->V_1_index[chan1] = length;
    length += Chan.npp[chan1] * Chan.nhh[chan1];
  }
  this->V_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_1[ind] = 0.0; }
  this->V_1_length = length;*/
}

void V_pphh::Delete()
{
  //delete[] this->V_1;
  //delete[] this->V_1_index;
}

void V_hphp::Build(Channels &Chan)
{
  int chan2, ind, length;
  length = 0;
  this->V_2_1_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->V_2_1_index[chan2] = length;
    length += Chan.nhp1[chan2] * Chan.nhp1[chan2];
  }
  this->V_2_1 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_2_1[ind] = 0.0; }
  this->V_2_1_length = length;

  length = 0;
  this->V_2_2_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->V_2_2_index[chan2] = length;
    length += Chan.nph1[chan2] * Chan.nph1[chan2];
  }
  this->V_2_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_2_2[ind] = 0.0; }
  this->V_2_2_length = length;
}

void V_hphp::Delete()
{
  delete[] this->V_2_1;
  delete[] this->V_2_2;
  delete[] this->V_2_1_index;
  delete[] this->V_2_2_index;
}

void V_hhhp::Build(Channels &Chan)
{
  int chan2, chan3, ind, length;
  length = 0;
  this->V_2_3_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->V_2_3_index[chan2] = length;
    length += Chan.nhh1[chan2] * Chan.nph1[chan2];
  }
  this->V_2_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_2_3[ind] = 0.0; }
  this->V_2_3_length = length;

  length = 0;
  this->V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_2_index[chan3] = length;
    length += Chan.nh[chan3] * Chan.nhph[chan3];
  }
  this->V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_2[ind] = 0.0; }
  this->V_3_2_length = length;

  length = 0;
  this->V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_3_index[chan3] = length;
    length += Chan.nhhp[chan3] * Chan.nh[chan3];
  }
  this->V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_3[ind] = 0.0; }
  this->V_3_3_length = length;
}

void V_hhhp::Delete()
{
  delete[] this->V_2_3;
  delete[] this->V_3_2;
  delete[] this->V_3_3;
  delete[] this->V_2_3_index;
  delete[] this->V_3_2_index;
  delete[] this->V_3_3_index;
}

void V_hppp::Build(Channels &Chan)
{
  int chan2, chan3, ind, length;
  length = 0;
  this->V_2_4_index = new int[Chan.size2];
  for(chan2 = 0; chan2 < Chan.size2; ++chan2){
    this->V_2_4_index[chan2] = length;
    length += Chan.npp1[chan2] * Chan.nph1[chan2];
  }
  this->V_2_4 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_2_4[ind] = 0.0; }
  this->V_2_4_length = length;

  length = 0;
  this->V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.npph[chan3];
  }
  this->V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_2[ind] = 0.0; }
  this->V_3_2_length = length;
}

void V_hppp::Delete()
{
  delete[] this->V_2_4;
  delete[] this->V_3_2;
  delete[] this->V_2_4_index;
  delete[] this->V_3_2_index;
}

void V_hphh::Build(Channels &Chan)
{
  int chan3, ind, length;
  length = 0;
  this->V_3_2_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_2_index[chan3] = length;
    length += Chan.np[chan3] * Chan.nhhh[chan3];
  }
  this->V_3_2 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_2[ind] = 0.0; }
  this->V_3_2_length = length;
}

void V_hphh::Delete()
{
  delete[] this->V_3_2;
  delete[] this->V_3_2_index;
}

void V_pphp::Build(Channels &Chan)
{
  int chan3, ind, length;
  length = 0;
  this->V_3_3_index = new int[Chan.size3];
  for(chan3 = 0; chan3 < Chan.size3; ++chan3){
    this->V_3_3_index[chan3] = length;
    length += Chan.nppp[chan3] * Chan.nh[chan3];
  }
  this->V_3_3 = new double[length];
  for(ind = 0; ind < length; ++ind){ this->V_3_3[ind] = 0.0; }
  this->V_3_3_length = length;
}

void V_pphp::Delete()
{
  delete[] this->V_3_3;
  delete[] this->V_3_3_index;
}


void Interactions::get_Eref(Channels &Chan)
{
  double energy = 0.0;
  State tb;
  int nh, nhh, chan_ind, ind;

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    chan_ind = this->Vhhhh.V_1_index[chan];
    for(int hh = 0; hh < nhh; ++hh){
      ind = hh*nhh + hh;
      if(PAR.basis == "finite_J"){ energy += (Chan.qnums1[chan].j + 1) * this->Vhhhh.V_1[chan_ind + ind]; }
      else{ energy += this->Vhhhh.V_1[chan_ind + ind]; }
      //std::cout << std::setprecision(10) << "Eref1: " << this->Vhhhh.V_1[chan_ind + ind] << std::endl;
    }
  }
  energy *= -0.5;

  if(PAR.HF == 1){
    for(int i = 0; i < SPB.num_hol; ++i){
      if(PAR.basis == "finite_J"){ energy += (SPB.qnums[i].j + 1) * SPB.qnums[i].energy; }
      else{ energy += SPB.qnums[i].energy; }
      //std::cout << std::setprecision(10) << "Eref2: " << SPB.qnums[i].energy << std::endl;
    }
  }
  else if(PAR.HF == 0){
    for(int chan = 0; chan < Chan.size3; ++chan){
      nh = Chan.nh[chan];
      chan_ind = this->Fmatrix.hh_3_index[chan];
      for(int h = 0; h < nh; ++h){
	ind = h*nh + h;
	if(PAR.basis == "finite_J"){ energy += (Chan.qnums3[chan].j + 1) * this->Fmatrix.hh_3[chan_ind + ind]; }
	else{ energy += this->Fmatrix.hh_3[chan_ind + ind]; }
      }
    }
  }
  this->Eref += energy;
}


void Coulomb_Inf_Matrix_Elements(Channels &Chan, Interactions &Ints)
{
  double L = pow(SPB.num_hol/PAR.density, 1.0/3.0);
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
    //chan_ind2 = Ints.Vpphh.V_1_index[chan1];
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
	TBME = Coulomb_Inf(p1, p2, p3, p4, L);
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp1_ind + pp2_ind)] = TBME;
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp2_ind + pp1_ind)] = TBME;
      }
      for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	h1 = hh1.v1;
	h2 = hh1.v2;
	if(h1 == h2){ continue; }
	TBME = Coulomb_Inf(p1, p2, h1, h2, L);
	//Ints.Vpphh.V_1[chan_ind2 + (nhh*pp1_ind + hh1_ind)] = TBME;
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
	TBME = Coulomb_Inf(h1, h2, h3, h4, L);
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
	TBME = Coulomb_Inf(h1, h2, p1, p2, L);
	Ints.Vhhpp.V_2_1[chan_ind1 + (nph1*hp1_ind + ph1_ind)] = TBME;
      }
      for(int hp2_ind = hp1_ind; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	h2 = hp2.v1;
	p1 = hp2.v2;
	TBME = Coulomb_Inf(h1, p1, h2, p2, L);
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
	TBME = Coulomb_Inf(h1, h2, p1, p2, L);
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
	TBME = Coulomb_Inf(h1, h2, p1, p2, L);
	Ints.Vhhpp.V_3_3[chan_ind2 + (np*hhp_ind + p_ind)] = TBME;
      }
    }
  }
}

double Coulomb_Inf(int &qi, int &qj, int &qk, int &ql, double &L)
{
  double term = 0.0;
  double e_sq = hbarc_HartA * fine_struct;
  double prefactor = e_sq/(L*L*L);
  double qSquared1;
  double kX1, kY1, kZ1;
  if(SPB.qnums[qi].nx + SPB.qnums[qj].nx != SPB.qnums[qk].nx + SPB.qnums[ql].nx || 
     SPB.qnums[qi].ny + SPB.qnums[qj].ny != SPB.qnums[qk].ny + SPB.qnums[ql].ny || 
     SPB.qnums[qi].nz + SPB.qnums[qj].nz != SPB.qnums[qk].nz + SPB.qnums[ql].nz){ return 0.0; }
  if(SPB.qnums[qi].m == SPB.qnums[qk].m && SPB.qnums[qj].m == SPB.qnums[ql].m){
    kX1 = (2.0*M_PI/L) * (SPB.qnums[qi].nx - SPB.qnums[qk].nx);
    kY1 = (2.0*M_PI/L) * (SPB.qnums[qi].ny - SPB.qnums[qk].ny);
    kZ1 = (2.0*M_PI/L) * (SPB.qnums[qi].nz - SPB.qnums[qk].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-15){ term += 0.0; }
    else{ term += 4.0 * prefactor * M_PI / qSquared1; }
  }
  if(SPB.qnums[qi].m == SPB.qnums[ql].m && SPB.qnums[qj].m == SPB.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (SPB.qnums[qi].nx - SPB.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (SPB.qnums[qi].ny - SPB.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (SPB.qnums[qi].nz - SPB.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-15){ term -= 0.0; }
    else{ term -= 4.0 * prefactor * M_PI / qSquared1; }
  }
  return term;
}

void Minnesota_Matrix_Elements(Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;  
  double L = pow((PAR.P + PAR.N)/PAR.density, 1./3.);
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
    //chan_ind2 = Ints.Vpphh.V_1_index[chan1];
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
	TBME = vint_Minnesota_Momentum(p1, p2, p3, p4, L);
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp1_ind + pp2_ind)] = TBME;
	Ints.Vpppp.V_1[chan_ind1 + (npp*pp2_ind + pp1_ind)] = TBME;
      }
      for(int hh1_ind = 0; hh1_ind < nhh; ++hh1_ind){
	hh1 = Chan.hh_state(chan1, hh1_ind);
	h1 = hh1.v1;
	h2 = hh1.v2;
	if(h1 == h2){ continue; }
	TBME = vint_Minnesota_Momentum(p1, p2, h1, h2, L);
	//Ints.Vpphh.V_1[chan_ind2 + (nhh*pp1_ind + hh1_ind)] = TBME;
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
	TBME = vint_Minnesota_Momentum(h1, h2, h3, h4, L);
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
	TBME = vint_Minnesota_Momentum(h1, h2, p1, p2, L);
	Ints.Vhhpp.V_2_1[chan_ind1 + (nph1*hp1_ind + ph1_ind)] = TBME;
      }
      for(int hp2_ind = hp1_ind; hp2_ind < nhp1; ++hp2_ind){
	hp2 = Chan.hp1_state(chan2, hp2_ind);
	h2 = hp2.v1;
	p1 = hp2.v2;
	TBME = vint_Minnesota_Momentum(h1, p1, h2, p2, L);
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
	TBME = vint_Minnesota_Momentum(h1, h2, p1, p2, L);
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
	TBME = vint_Minnesota_Momentum(h1, h2, p1, p2, L);
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
double vint_Minnesota_Momentum(int &qi, int &qj, int &qk, int &ql, double &L)
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

  if(SPB.qnums[qi].nx + SPB.qnums[qj].nx != SPB.qnums[qk].nx + SPB.qnums[ql].nx){ return 0.0; }
  if(SPB.qnums[qi].ny + SPB.qnums[qj].ny != SPB.qnums[qk].ny + SPB.qnums[ql].ny){ return 0.0; }
  if(SPB.qnums[qi].nz + SPB.qnums[qj].nz != SPB.qnums[qk].nz + SPB.qnums[ql].nz){ return 0.0; }
  if(SPB.qnums[qi].m + SPB.qnums[qj].m != SPB.qnums[qk].m + SPB.qnums[ql].m){ return 0.0; }
  if(SPB.qnums[qi].t + SPB.qnums[qj].t != SPB.qnums[qk].t + SPB.qnums[ql].t){ return 0.0; }

  kX1 = (M_PI/L) * (SPB.qnums[qi].nx - SPB.qnums[qj].nx - SPB.qnums[qk].nx + SPB.qnums[ql].nx);
  kY1 = (M_PI/L) * (SPB.qnums[qi].ny - SPB.qnums[qj].ny - SPB.qnums[qk].ny + SPB.qnums[ql].ny);
  kZ1 = (M_PI/L) * (SPB.qnums[qi].nz - SPB.qnums[qj].nz - SPB.qnums[qk].nz + SPB.qnums[ql].nz);

  kX2 = (M_PI/L) * (SPB.qnums[qi].nx - SPB.qnums[qj].nx - SPB.qnums[ql].nx + SPB.qnums[qk].nx);
  kY2 = (M_PI/L) * (SPB.qnums[qi].ny - SPB.qnums[qj].ny - SPB.qnums[ql].ny + SPB.qnums[qk].ny);
  kZ2 = (M_PI/L) * (SPB.qnums[qi].nz - SPB.qnums[qj].nz - SPB.qnums[ql].nz + SPB.qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeMtxEle(SPB.qnums[qi].m, SPB.qnums[qj].m, SPB.qnums[qk].m, SPB.qnums[ql].m);
  isoSpinEx1 = spinExchangeMtxEle(SPB.qnums[qi].t, SPB.qnums[qj].t, SPB.qnums[qk].t, SPB.qnums[ql].t);

  spinEx2 = spinExchangeMtxEle(SPB.qnums[qi].m, SPB.qnums[qj].m, SPB.qnums[ql].m, SPB.qnums[qk].m);
  isoSpinEx2 = spinExchangeMtxEle(SPB.qnums[qi].t, SPB.qnums[qj].t, SPB.qnums[ql].t, SPB.qnums[qk].t);
  
  IsIt1 = kron_del(SPB.qnums[qi].m, SPB.qnums[qk].m) * kron_del(SPB.qnums[qj].m, SPB.qnums[ql].m) * 
    kron_del(SPB.qnums[qi].t, SPB.qnums[qk].t) * kron_del(SPB.qnums[qj].t, SPB.qnums[ql].t);
  PsIt1 = spinEx1 * kron_del(SPB.qnums[qi].t, SPB.qnums[qk].t) * kron_del(SPB.qnums[qj].t, SPB.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(SPB.qnums[qi].m, SPB.qnums[qk].m)*kron_del(SPB.qnums[qj].m, SPB.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kron_del(SPB.qnums[qi].m, SPB.qnums[ql].m) * kron_del(SPB.qnums[qj].m, SPB.qnums[qk].m) * 
    kron_del(SPB.qnums[qi].t, SPB.qnums[ql].t) * kron_del(SPB.qnums[qj].t, SPB.qnums[qk].t);
  PsIt2 = spinEx2 * kron_del(SPB.qnums[qi].t, SPB.qnums[ql].t) * kron_del(SPB.qnums[qj].t, SPB.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(SPB.qnums[qi].m, SPB.qnums[ql].m) * kron_del(SPB.qnums[qj].m, SPB.qnums[qk].m) * isoSpinEx2;

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
double Coulomb_HO(int &qi, int &qj, int &qk, int &ql)
{
  int g1, g2, g3, g4, G, L;
  double dir = 0.0;
  double exch = 0.0;
  int n1, m1, n2, m2, n3, m3, n4, m4;
  double LogRatio1, LogProd2, LogRatio2, temp;
  n1 = SPB.qnums[qi].n; //1
  m1 = SPB.qnums[qi].ml;
  n2 = SPB.qnums[qj].n; //2
  m2 = SPB.qnums[qj].ml;
  n3 = SPB.qnums[ql].n; //4
  m3 = SPB.qnums[ql].ml;
  n4 = SPB.qnums[qk].n; //3
  m4 = SPB.qnums[qk].ml;
  if((m1 + m2 == m3 + m4) && SPB.qnums[qi].m == SPB.qnums[qk].m && SPB.qnums[qj].m == SPB.qnums[ql].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = j1 + j4 + (std::abs(m1) + m1 + std::abs(m4) - m4)/2;
	    g2 = j2 + j3 + (std::abs(m2) + m2 + std::abs(m3) - m3)/2;
	    g3 = j3 + j2 + (std::abs(m3) + m3 + std::abs(m2) - m2)/2;
	    g4 = j4 + j1 + (std::abs(m4) + m4 + std::abs(m1) - m1)/2;
	    G = g1 + g2 + g3 + g4;
	    LogRatio1 = logratio1(j1, j2, j3, j4);
	    LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    LogRatio2 = logratio2(G);
	    temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * std::exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    dir += (-2*((j1 + j2 + j3 + j4)%2) + 1) * std::exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    dir *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }

  n1 = SPB.qnums[qi].n; //1
  m1 = SPB.qnums[qi].ml;
  n2 = SPB.qnums[qj].n; //2
  m2 = SPB.qnums[qj].ml;
  n3 = SPB.qnums[qk].n; //3
  m3 = SPB.qnums[qk].ml;
  n4 = SPB.qnums[ql].n; //4
  m4 = SPB.qnums[ql].ml;
  if((m1 + m2 == m3 + m4) && SPB.qnums[qi].m == SPB.qnums[ql].m && SPB.qnums[qj].m == SPB.qnums[qk].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = j1 + j4 + (std::abs(m1) + m1 + std::abs(m4) - m4)/2;
	    g2 = j2 + j3 + (std::abs(m2) + m2 + std::abs(m3) - m3)/2;
	    g3 = j3 + j2 + (std::abs(m3) + m3 + std::abs(m2) - m2)/2;
	    g4 = j4 + j1 + (std::abs(m4) + m4 + std::abs(m1) - m1)/2;
	    G = g1 + g2 + g3 + g4;
	    LogRatio1 = logratio1(j1, j2, j3, j4);
	    LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    LogRatio2 = logratio2(G);
	    temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * std::exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    exch += (-2*((j1 + j2 + j3 + j4)%2) + 1) * std::exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    exch *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }
  return std::sqrt(PAR.density)*(dir - exch);
}

void Get_Fock_Matrix(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Single_Particle_States &HF, Channels &Chan, Interactions &Ints)
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
	  fock0 += SPB.qnums[t_ind].energy * HF.vectors[vector_ind + (p*nob0 + t)] * HF.vectors[vector_ind + (q*nob0 + t)];
	}
	fock += fock0;
	for(int chan = 0; chan < HF_Chan.size3; ++chan){
	  nob = HF_Chan.nob[chan];
	  for(int t = 0; t < nob; ++t){ // Sum over occupied levels
	    t_ind = HF_Chan.ob_state(chan, t).v1;
	    if(t_ind >= SPB.num_hol){ continue; } 
	    plus(tb, SPB.qnums[p_ind], SPB.qnums[t_ind]);
	    minj = abs(HF_Chan.qnums3[chan0].j - HF_Chan.qnums3[chan].j);
	    while(tb.j >= minj){
	      if( (p_ind == t_ind || q_ind == t_ind) && tb.j%4 != 0 ){ tb.j -= 2; continue; }
	      chan1 = Ind_Chan1(tb);
	      key1 = HF_Chan.tb_map[chan1][Hash(p_ind, t_ind, tb.j)];
	      key2 = HF_Chan.tb_map[chan1][Hash(q_ind, t_ind, tb.j)];
	      ind1 = HF_ME.Index[chan1] + (key1 * HF_Chan.ntb[chan1] + key2);
	      fock += (tb.j + 1.0)/(HF_Chan.qnums3[chan0].j + 1.0) * HF_ME.V[ind1];
	      tb.j -= 2;
	    }
	  }
	}
	if(p_ind < SPB.num_hol){
	  key1 = Chan.h_map[chan0][p_ind];
	  if(q_ind < SPB.num_hol){ // hh
	    key2 = Chan.h_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.hh_3_index[chan0];
	    ind3 = key1 * Chan.nh[chan0] + key2;
	    minus(tb, SPB.qnums[p_ind], SPB.qnums[q_ind]);
	    ind2 = Chan.hh1_map[Chan.ind0][Hash(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.hh_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.hh_2[ind2] = fock;
	  }
	  else{ // hp
	    key2 = Chan.p_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.hp_3_index[chan0];
	    ind3 = key1 * Chan.np[chan0] + key2;
 	    minus(tb, SPB.qnums[p_ind], SPB.qnums[q_ind]);
	    ind2 = Chan.hp1_map[Chan.ind0][Hash(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.hp_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.hp_2[ind2] = fock;
	  }
	}
	else{
	  key1 = Chan.p_map[chan0][p_ind];
	  if(q_ind < SPB.num_hol){ // ph
	    key2 = Chan.h_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.ph_3_index[chan0];
	    ind3 = key1 * Chan.nh[chan0] + key2;
 	    minus(tb, SPB.qnums[p_ind], SPB.qnums[q_ind]);
	    ind2 = Chan.ph1_map[Chan.ind0][Hash(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.ph_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.ph_2[ind2] = fock;
	  }
	  else{ // pp
	    key2 = Chan.p_map[chan0][q_ind];
	    chan_ind = Ints.Fmatrix.pp_3_index[chan0];
	    ind3 = key1 * Chan.np[chan0] + key2;
 	    minus(tb, SPB.qnums[p_ind], SPB.qnums[q_ind]);
	    ind2 = Chan.pp1_map[Chan.ind0][Hash(p_ind, q_ind, tb.j)];
	    Ints.Fmatrix.pp_3[chan_ind + ind3] = fock;
	    Ints.Fmatrix.pp_2[ind2] = fock;
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints)
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
      int ptype, qtype, rtype, stype;
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
	ptype = SPB.qnums[p].type;
	qtype = SPB.qnums[q].type;
	if(p == q){ continue; }
	r = tb2.v1;
	s = tb2.v2;
	rtype = SPB.qnums[r].type;
	stype = SPB.qnums[s].type;
	if(r == s){ continue; }
	TBME = HF_ME.V[HF_ME.Index[chan1] + (tb1_ind*ntb + tb2_ind)];
	/*if(p < q && r < s && p <= r){
	  std::cout << std::setprecision(12) << "V_hf: " << p << " " << q << " " << r << " " << s << " = " << TBME << std::endl;
	  }*/
	if(rtype == 1 && stype == 1){
	  if(ptype == 1 && qtype == 1){
	    key1 = Chan.pp_map[chan1][Hash(p, q, 0)];
	    key2 = Chan.pp_map[chan1][Hash(r, s, 0)];
	    chan_ind = Ints.Vpppp.V_1_index[chan1];
	    ind = key1 * Chan.npp[chan1] + key2;
	    Ints.Vpppp.V_1[chan_ind + ind] = TBME;
	    
	    /*chan = Ind_Chan3(SPB.qnums[r]);
	    key1 = Chan.ppp_map[chan][Hash(p, q, s, 0)];
	    key2 = Chan.p_map[chan][r];
	    chan_ind = Ints.Vpppp.V_3_3_index[chan];
	    ind = key1 * Chan.np[chan] + key2;
	    Ints.Vpppp.V_3_3[chan_ind + ind] = TBME;*/
	  }
	  else if(ptype == 0 && qtype == 1){
	    chan = Ind_Chan3(SPB.qnums[q]);
	    key1 = Chan.p_map[chan][q];
	    key2 = Chan.pph_map[chan][Hash(r, s, p, 0)];
	    chan_ind = Ints.Vhppp.V_3_2_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhppp.V_3_2[chan_ind + ind] = TBME;
	    
	    minus(tb, SPB.qnums[q], SPB.qnums[s]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.pp1_map[chan][Hash(q, s, 0)];
	    key2 = Chan.ph1_map[chan][Hash(r, p, 0)];
	    chan_ind = Ints.Vhppp.V_2_4_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhppp.V_2_4[chan_ind + ind] = TBME;
	    
	    // Vpphp -> rspq
	    chan = Ind_Chan3(SPB.qnums[p]);
	    key1 = Chan.ppp_map[chan][Hash(r, s, q, 0)];
	    key2 = Chan.h_map[chan][p];
	    chan_ind = Ints.Vpphp.V_3_3_index[chan];
	    ind = key1 * Chan.nh[chan] + key2;
	    Ints.Vpphp.V_3_3[chan_ind + ind] = TBME;
	  }
	  else if(ptype == 0 && qtype == 0){
	    key1 = Chan.hh_map[chan1][Hash(p, q, 0)];
	    key2 = Chan.pp_map[chan1][Hash(r, s, 0)];
	    chan_ind = Ints.Vhhpp.V_1_index[chan1];
	    ind = key1 * Chan.npp[chan1] + key2;
	    Ints.Vhhpp.V_1[chan_ind + ind] = TBME;
	    /*chan_ind = Ints.Vpphh.V_1_index[chan1];
	    ind = key2 * Chan.nhh[chan1] + key1;
	    Ints.Vpphh.V_1[chan_ind + ind] = TBME;*/
	    
	    chan = Ind_Chan3(SPB.qnums[p]);
	    key1 = Chan.h_map[chan][p];
	    key2 = Chan.pph_map[chan][Hash(r, s, q, 0)];
	    chan_ind = Ints.Vhhpp.V_3_1_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhhpp.V_3_1[chan_ind + ind] = TBME;
	    
	    chan = Ind_Chan3(SPB.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.pph_map[chan][Hash(r, s, p, 0)];
	    chan_ind = Ints.Vhhpp.V_3_2_index[chan];
	    ind = key1 * Chan.npph[chan] + key2;
	    Ints.Vhhpp.V_3_2[chan_ind + ind] = TBME;
	    
	    chan = Ind_Chan3(SPB.qnums[r]);
	    key1 = Chan.hhp_map[chan][Hash(p, q, s, 0)];
	    key2 = Chan.p_map[chan][r];
	    chan_ind = Ints.Vhhpp.V_3_3_index[chan];
	    ind = key1 * Chan.np[chan] + key2;
	    Ints.Vhhpp.V_3_3[chan_ind + ind] = TBME;
	    
	    minus(tb, SPB.qnums[p], SPB.qnums[s]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, s, 0)];
	    key2 = Chan.ph1_map[chan][Hash(r, q, 0)];
	    chan_ind = Ints.Vhhpp.V_2_1_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhpp.V_2_1[chan_ind + ind] = TBME;
	    
	    minus(tb, SPB.qnums[p], SPB.qnums[r]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, r, 0)];
	    key2 = Chan.ph1_map[chan][Hash(s, q, 0)];
	    chan_ind = Ints.Vhhpp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhpp.V_2_3[chan_ind + ind] = TBME;
	  } 
	}
	else if(rtype == 0 && stype == 1){
	  if(ptype == 0 && qtype == 1){
	    minus(tb, SPB.qnums[p], SPB.qnums[s]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, s, 0)];
	    key2 = Chan.hp1_map[chan][Hash(r, q, 0)];
	    chan_ind = Ints.Vhphp.V_2_1_index[chan];
	    ind = key1 * Chan.nhp1[chan] + key2;
	    Ints.Vhphp.V_2_1[chan_ind + ind] = TBME;
	    
	    minus(tb, SPB.qnums[q], SPB.qnums[r]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.ph1_map[chan][Hash(q, r, 0)];
	    key2 = Chan.ph1_map[chan][Hash(s, p, 0)];
	    chan_ind = Ints.Vhphp.V_2_2_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhphp.V_2_2[chan_ind + ind] = TBME;
	  }
	  else if(ptype == 0 && qtype == 0){
	    chan = Ind_Chan3(SPB.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.hph_map[chan][Hash(r, s, p, 0)];
	    chan_ind = Ints.Vhhhp.V_3_2_index[chan];
	    ind = key1 * Chan.nhph[chan] + key2;
	    Ints.Vhhhp.V_3_2[chan_ind + ind] = TBME;
	    
	    chan = Ind_Chan3(SPB.qnums[r]);
	    key1 = Chan.hhp_map[chan][Hash(p, q, s, 0)];
	    key2 = Chan.h_map[chan][r];
	    chan_ind = Ints.Vhhhp.V_3_3_index[chan];
	    ind = key1 * Chan.nh[chan] + key2;
	    Ints.Vhhhp.V_3_3[chan_ind + ind] = TBME;
	    
	    minus(tb, SPB.qnums[p], SPB.qnums[r]);
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hh1_map[chan][Hash(p, r, 0)];
	    key2 = Chan.ph1_map[chan][Hash(s, q, 0)];
	    chan_ind = Ints.Vhhhp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    Ints.Vhhhp.V_2_3[chan_ind + ind] = TBME;
	    
	    // Vhphh -> rspq
	    chan = Ind_Chan3(SPB.qnums[s]);
	    key1 = Chan.p_map[chan][s];
	    key2 = Chan.hhh_map[chan][Hash(p, q, r, 0)];
	    chan_ind = Ints.Vhphh.V_3_2_index[chan];
	    ind = key1 * Chan.nhhh[chan] + key2;
	    Ints.Vhphh.V_3_2[chan_ind + ind] = TBME;
	  }
	}
	else if(rtype == 0 && stype == 0){
	  if(ptype == 0 && qtype == 0){
	    key1 = Chan.hh_map[chan1][Hash(p, q, 0)];
	    key2 = Chan.hh_map[chan1][Hash(r, s, 0)];
	    chan_ind = Ints.Vhhhh.V_1_index[chan1];
	    ind = key1 * Chan.nhh[chan1] + key2;
	    Ints.Vhhhh.V_1[chan_ind + ind] = TBME;
	    
	    chan = Ind_Chan3(SPB.qnums[q]);
	    key1 = Chan.h_map[chan][q];
	    key2 = Chan.hhh_map[chan][Hash(r, s, p, 0)];
	    chan_ind = Ints.Vhhhh.V_3_2_index[chan];
	    ind = key1 * Chan.nhhh[chan] + key2;
	    Ints.Vhhhh.V_3_2[chan_ind + ind] = TBME;
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements_J(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints)
{
  struct timespec time1, time2;
  double elapsed0 = 0.0;
  clock_gettime(CLOCK_MONOTONIC, &time1);

  int ntb, length;
  two_body pq, rs;
  four_body pqrs;
  four_body *pqrs_vec;
  int *pqrs_chan;
  int *pqrs_j;
  length = 0;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    length += ntb * ntb;
  }
  pqrs_vec = new four_body[length];
  pqrs_chan = new int[length];
  pqrs_j = new int[length];
  length = 0;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    for(int pq_ind = 0; pq_ind < ntb; ++pq_ind){
      pq = HF_Chan.tb_state(chan1, pq_ind);
      pqrs.v1 = pq.v1;
      pqrs.v2 = pq.v2;
      for(int rs_ind = 0; rs_ind < ntb; ++rs_ind){
	rs = HF_Chan.tb_state(chan1, rs_ind);
	pqrs.v3 = rs.v1;
	pqrs.v4 = rs.v2;
	pqrs_vec[length] = pqrs;
	pqrs_chan[length] = chan1;
	pqrs_j[length] = HF_Chan.qnums1[chan1].j;
	++length;
      }
    }
  }

  #pragma omp parallel
  {
    int p, q, r, s, pj, qj, rj, sj, J, jmin;
    int ind, chan1, chan, chan_ind, key1, key2;
    int ptype, qtype, rtype, stype;
    double X, TBME;
    State tb;
    #pragma omp for schedule(static)
    for(int ind1 = 0; ind1 < length; ++ind1){
      p = pqrs_vec[ind1].v1;
      q = pqrs_vec[ind1].v2;
      r = pqrs_vec[ind1].v3;
      s = pqrs_vec[ind1].v4;
      ptype = SPB.qnums[p].type;
      qtype = SPB.qnums[q].type;
      rtype = SPB.qnums[r].type;
      stype = SPB.qnums[s].type;
      pj = SPB.qnums[p].j;
      qj = SPB.qnums[q].j;
      rj = SPB.qnums[r].j;
      sj = SPB.qnums[s].j;
      TBME = HF_ME.V[ind1];
      chan1 = pqrs_chan[ind1];
      J = pqrs_j[ind1];
      /*if((p <= q && r <= s && (p < r || (p == r && q <=s))) && std::fabs(TBME) > 1.0e-10){
	std::cout << std::setprecision(12) << "V_hf: " << p+1 << " " << q+1 << " " << r+1 << " " << s+1 << ", ";
	std::cout << 0.5*HF_Chan.qnums1[chan1].j << " = " << TBME << std::endl;
	}*/
      if(rtype == 1 && stype == 1){
	if(ptype == 1 && qtype == 1){
	  key1 = Chan.pp_map[chan1][Hash(p, q, J)];
	  key2 = Chan.pp_map[chan1][Hash(r, s, J)];
	  chan_ind = Ints.Vpppp.V_1_index[chan1];
	  ind = key1 * Chan.npp[chan1] + key2;
	  Ints.Vpppp.V_1[chan_ind + ind] = TBME;
	  
	  /*chan = Ind_Chan3(SPB.qnums[r]);
	  key1 = Chan.ppp_map[chan][Hash(p, q, s, J)];
	  key2 = Chan.p_map[chan][r];
	  chan_ind = Ints.Vpppp.V_3_3_index[chan];
	  ind = key1 * Chan.np[chan] + key2;
	  Ints.Vpppp.V_3_3[chan_ind + ind] = phase2(rj + sj - J) * std::sqrt((J + 1.0)/(rj + 1.0)) * TBME;*/
	}
	else if(ptype == 0 && qtype == 1){
	  chan = Ind_Chan3(SPB.qnums[q]);
	  key1 = Chan.p_map[chan][q];
	  key2 = Chan.pph_map[chan][Hash(r, s, p, J)];
	  chan_ind = Ints.Vhppp.V_3_2_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhppp.V_3_2[chan_ind + ind] = std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  //Ints.Vhppp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  
	  minus(tb, SPB.qnums[q], SPB.qnums[s]);
	  if(SPB.qnums[r].j + SPB.qnums[p].j < tb.j){ tb.j = rj + pj; }
	  jmin = std::abs(qj - sj);
	  if(std::abs(rj - pj) > jmin){ jmin = std::abs(rj - pj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.pp1_map[chan][Hash(q, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Hash(r, p, tb.j)];
	    chan_ind = Ints.Vhppp.V_2_4_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = phase2(pj + qj - J) * (J + 1.0) * CGC6(qj,pj,J,rj,sj,tb.j);
	    Ints.Vhppp.V_2_4[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }	  
	  
	  // Vpphp -> rspq
	  chan = Ind_Chan3(SPB.qnums[p]);
	  key1 = Chan.ppp_map[chan][Hash(r, s, q, J)];
	  key2 = Chan.h_map[chan][p];
	  chan_ind = Ints.Vpphp.V_3_3_index[chan];
	  ind = key1 * Chan.nh[chan] + key2;
	  Ints.Vpphp.V_3_3[chan_ind + ind] = phase2(pj + qj - J) * std::sqrt((J + 1.0)/(pj + 1.0)) * TBME;
	}
	else if(ptype == 0 && qtype == 0){
	  key1 = Chan.hh_map[chan1][Hash(p, q, J)];
	  key2 = Chan.pp_map[chan1][Hash(r, s, J)];
	  chan_ind = Ints.Vhhpp.V_1_index[chan1];
	  ind = key1 * Chan.npp[chan1] + key2;
	  Ints.Vhhpp.V_1[chan_ind + ind] = TBME;
	  /*chan_ind = Ints.Vpphh.V_1_index[chan1];
	  ind = key2 * Chan.nhh[chan1] + key1;
	  Ints.Vpphh.V_1[chan_ind + ind] = TBME;*/

	  chan = Ind_Chan3(SPB.qnums[p]);
	  key1 = Chan.h_map[chan][p];
	  key2 = Chan.pph_map[chan][Hash(r, s, q, J)];
	  chan_ind = Ints.Vhhpp.V_3_1_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhhpp.V_3_1[chan_ind + ind] = phase2(pj + qj - J) * std::sqrt((J + 1.0)/(pj + 1.0)) * TBME;
	  
	  chan = Ind_Chan3(SPB.qnums[r]);
	  key1 = Chan.hhp_map[chan][Hash(p, q, s, J)];
	  key2 = Chan.p_map[chan][r];
	  chan_ind = Ints.Vhhpp.V_3_3_index[chan];
	  ind = key1 * Chan.np[chan] + key2;
	  Ints.Vhhpp.V_3_3[chan_ind + ind] = phase2(rj + sj - J) * std::sqrt((J + 1.0)/(rj + 1.0)) * TBME;
	  
	  minus(tb, SPB.qnums[p], SPB.qnums[s]);
	  if(rj + qj < tb.j){ tb.j = rj + qj; }
	  jmin = std::abs(pj - sj);
	  if(std::abs(rj - qj) > jmin){ jmin = std::abs(rj - qj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, s, tb.j)];
	    key2 = Chan.ph1_map[chan][Hash(r, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_1_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = (J + 1.0) * CGC6(pj,qj,J,rj,sj,tb.j);
	    //X = -1.0 * (J + 1.0) * CGC6(pj,qj,J,rj,sj,tb.j);
	    Ints.Vhhpp.V_2_1[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	  
	  chan = Ind_Chan3(SPB.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.pph_map[chan][Hash(r, s, p, J)];
	  chan_ind = Ints.Vhhpp.V_3_2_index[chan];
	  ind = key1 * Chan.npph[chan] + key2;
	  Ints.Vhhpp.V_3_2[chan_ind + ind] = std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  //Ints.Vhhpp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  
	  minus(tb, SPB.qnums[p], SPB.qnums[r]);
	  if(sj + qj < tb.j){ tb.j = sj + qj; }
	  jmin = std::abs(pj - rj);
	  if(std::abs(sj - qj) > jmin){ jmin = std::abs(sj - qj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Hash(s, q, tb.j)];
	    chan_ind = Ints.Vhhpp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = phase2(rj + sj - J) * (J + 1.0) * CGC6(pj,qj,J,sj,rj,tb.j);
	    Ints.Vhhpp.V_2_3[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	}
      }
      else if(rtype == 0 && stype == 1){
	if(ptype == 0 && qtype == 1){
	  minus(tb, SPB.qnums[p], SPB.qnums[s]);
	  if(rj + qj < tb.j){ tb.j = rj + qj; }
	  jmin = std::abs(pj - sj);
	  if(std::abs(rj - qj) > jmin){ jmin = std::abs(rj - qj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hp1_map[chan][Hash(p, s, tb.j)];
	    key2 = Chan.hp1_map[chan][Hash(r, q, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_1_index[chan];
	    ind = key1 * Chan.nhp1[chan] + key2;
	    X = (J + 1.0) * CGC6(pj,qj,J,rj,sj,tb.j);
	    //X = -1.0 * (J + 1.0) * CGC6(pj,qj,J,rj,sj,tb.j);
	    Ints.Vhphp.V_2_1[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	  
	  minus(tb, SPB.qnums[q], SPB.qnums[r]);
	  if(sj + pj < tb.j){ tb.j = sj + pj; }
	  jmin = std::abs(qj - rj);
	  if(std::abs(sj - pj) > jmin){ jmin = std::abs(sj - pj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.ph1_map[chan][Hash(q, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Hash(s, p, tb.j)];
	    chan_ind = Ints.Vhphp.V_2_2_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = phase2(pj + qj + rj + sj) * (J + 1.0) * CGC6(qj,pj,J,sj,rj,tb.j);
	    //X = -1.0 * phase2(pj + qj + rj + sj) * (J + 1.0) * CGC6(qj,pj,J,sj,rj,tb.j);
	    Ints.Vhphp.V_2_2[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	}
	else if(ptype == 0 && qtype == 0){
	  chan = Ind_Chan3(SPB.qnums[r]);
	  key1 = Chan.hhp_map[chan][Hash(p, q, s, J)];
	  key2 = Chan.h_map[chan][r];
	  chan_ind = Ints.Vhhhp.V_3_3_index[chan];
	  ind = key1 * Chan.nh[chan] + key2;
	  Ints.Vhhhp.V_3_3[chan_ind + ind] = phase2(rj + sj - J) * std::sqrt((J + 1.0)/(rj + 1.0)) * TBME;
	  
	  chan = Ind_Chan3(SPB.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.hph_map[chan][Hash(r, s, p, J)];
	  chan_ind = Ints.Vhhhp.V_3_2_index[chan];
	  ind = key1 * Chan.nhph[chan] + key2;
	  Ints.Vhhhp.V_3_2[chan_ind + ind] = std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  //Ints.Vhhhp.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  
	  minus(tb, SPB.qnums[p], SPB.qnums[r]);
	  if(sj + qj < tb.j){ tb.j = sj + qj; }
	  jmin = std::abs(pj - rj);
	  if(std::abs(sj - qj) > jmin){ jmin = std::abs(sj - qj); }
	  while(tb.j >= jmin){
	    chan = Ind_Chan2(tb);
	    key1 = Chan.hh1_map[chan][Hash(p, r, tb.j)];
	    key2 = Chan.ph1_map[chan][Hash(s, q, tb.j)];
	    chan_ind = Ints.Vhhhp.V_2_3_index[chan];
	    ind = key1 * Chan.nph1[chan] + key2;
	    X = phase2(rj + sj - J) * (J + 1.0) * CGC6(pj,qj,J,sj,rj,tb.j);
	    Ints.Vhhhp.V_2_3[chan_ind + ind] += X * TBME;
	    tb.j -= 2;
	  }
	  
	  // Vhphh -> rspq
	  chan = Ind_Chan3(SPB.qnums[s]);
	  key1 = Chan.p_map[chan][s];
	  key2 = Chan.hhh_map[chan][Hash(p, q, r, J)];
	  chan_ind = Ints.Vhphh.V_3_2_index[chan];
	  ind = key1 * Chan.nhhh[chan] + key2;
	  Ints.Vhphh.V_3_2[chan_ind + ind] = std::sqrt((J + 1.0)/(sj + 1.0)) * TBME;
	  //Ints.Vhphh.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((J + 1.0)/(sj + 1.0)) * TBME;
	}
      }
      else if(rtype == 0 && stype == 0){
	if(ptype == 0 && qtype == 0){
	  key1 = Chan.hh_map[chan1][Hash(p, q, J)];
	  key2 = Chan.hh_map[chan1][Hash(r, s, J)];
	  chan_ind = Ints.Vhhhh.V_1_index[chan1];
	  ind = key1 * Chan.nhh[chan1] + key2;
	  Ints.Vhhhh.V_1[chan_ind + ind] = TBME;
	  
	  chan = Ind_Chan3(SPB.qnums[q]);
	  key1 = Chan.h_map[chan][q];
	  key2 = Chan.hhh_map[chan][Hash(r, s, p, J)];
	  chan_ind = Ints.Vhhhh.V_3_2_index[chan];
	  ind = key1 * Chan.nhhh[chan] + key2;
	  Ints.Vhhhh.V_3_2[chan_ind + ind] = std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	  //Ints.Vhhhh.V_3_2[chan_ind + ind] = -1.0 * std::sqrt((J + 1.0)/(qj + 1.0)) * TBME;
	}
      }
    }
  }
  delete[] pqrs_vec;
  delete[] pqrs_chan;
  delete[] pqrs_j;

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "Interaction load Runtime = " << elapsed0 << " sec. " << std::endl;
}

void Get_Matrix_Elements_JM(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, Channels &Chan, Interactions &Ints)
{
  int J, ntb, length;
  for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    ntb = HF_Chan.ntb[chan1];
    length = ntb * ntb;
    J = HF_Chan.qnums1[chan1].j;
    #pragma omp parallel
    {
      int p, q, r, s, pind, qind, rind, sind, tb1_ind, tb2_ind, chan_ind;
      two_body tb1, tb2;
      int pj, qj, rj, sj, pm, qm, rm, sm;
      int m1, m2, t1, t2;
      double CGC1, CGC2;
      int ind, chan, key1, key2;
      int ptype, qtype, rtype, stype;
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
	r = tb2.v1;
	s = tb2.v2;
	TBME0 = HF_ME.V[HF_ME.Index[chan1] + (tb1_ind*ntb + tb2_ind)];

	for(int M = -J; M <= J; M+=2){
	  for(int p1 = 0; p1 < SPB.shellsnum[p]; ++p1){
	    pind = SPB.shellsm[p][p1];
	    pj = SPB.shellsj[pind];
	    pm = SPB.qnums[pind].m;
	    ptype = SPB.qnums[pind].type;
	    for(int q1 = 0; q1 < SPB.shellsnum[q]; ++q1){
	      qind = SPB.shellsm[q][q1];
	      qj = SPB.shellsj[qind];
	      qm = SPB.qnums[qind].m;
	      qtype = SPB.qnums[qind].type;
	      m1 = pm + qm;
	      t1 = SPB.qnums[pind].t + SPB.qnums[qind].t;
	      if(pind == qind){ continue; }
	      for(int r1 = 0; r1 < SPB.shellsnum[r]; ++r1){
		rind = SPB.shellsm[r][r1];
		rj = SPB.shellsj[rind];
		rm = SPB.qnums[rind].m;
		rtype = SPB.qnums[rind].type;
		for(int s1 = 0; s1 < SPB.shellsnum[s]; ++s1){
		  sind = SPB.shellsm[s][s1];
		  sj = SPB.shellsj[sind];
		  sm = SPB.qnums[sind].m;
		  stype = SPB.qnums[sind].type;

		  m2 = rm + sm;
		  t2 = SPB.qnums[rind].t + SPB.qnums[sind].t;
		  if(rind == sind){ continue; }
		  
		  if(t1 != t2 || m1 != M || m2 != M){ continue; }
		  CGC1 = CGC(pj, pm, qj, qm, J, M);
		  CGC2 = CGC(rj, rm, sj, sm, J, M);
		  TBME = TBME0 * CGC1 * CGC2;

		  if(rtype == 1 && stype == 1){
		    if(ptype == 1 && qtype == 1){
		      plus(tb, SPB.qnums[pind], SPB.qnums[qind]);
		      chan = Ind_Chan1(tb);
		      key1 = Chan.pp_map[chan][Hash(pind, qind, 0)];
		      key2 = Chan.pp_map[chan][Hash(rind, sind, 0)];
		      chan_ind = Ints.Vpppp.V_1_index[chan];
		      ind = key1 * Chan.npp[chan] + key2;
		      Ints.Vpppp.V_1[chan_ind + ind] += TBME;
		      
		      /*chan = Ind_Chan3(SPB.qnums[rind]);
		      key1 = Chan.ppp_map[chan][Hash(pind, qind, sind, 0)];
		      key2 = Chan.p_map[chan][rind];
		      chan_ind = Ints.Vpppp.V_3_3_index[chan];
		      ind = key1 * Chan.np[chan] + key2;
		      Ints.Vpppp.V_3_3[chan_ind + ind] += TBME;*/
		    }
		    else if(ptype == 0 && qtype == 1){
		      chan = Ind_Chan3(SPB.qnums[qind]);
		      key1 = Chan.p_map[chan][qind];
		      key2 = Chan.pph_map[chan][Hash(rind, sind, pind, 0)];
		      chan_ind = Ints.Vhppp.V_3_2_index[chan];
		      ind = key1 * Chan.npph[chan] + key2;
		      Ints.Vhppp.V_3_2[chan_ind + ind] += TBME;
		      
		      minus(tb, SPB.qnums[qind], SPB.qnums[sind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.pp1_map[chan][Hash(qind, sind, 0)];
		      key2 = Chan.ph1_map[chan][Hash(rind, pind, 0)];
		      chan_ind = Ints.Vhppp.V_2_4_index[chan];
		      ind = key1 * Chan.nph1[chan] + key2;
		      Ints.Vhppp.V_2_4[chan_ind + ind] += TBME;
		      
		      // Vpphp -> rspq
		      chan = Ind_Chan3(SPB.qnums[pind]);
		      key1 = Chan.ppp_map[chan][Hash(rind, sind, qind, 0)];
		      key2 = Chan.h_map[chan][pind];
		      chan_ind = Ints.Vpphp.V_3_3_index[chan];
		      ind = key1 * Chan.nh[chan] + key2;
		      Ints.Vpphp.V_3_3[chan_ind + ind] += TBME;
		    }
		    else if(ptype == 0 && qtype == 0){
		      plus(tb, SPB.qnums[pind], SPB.qnums[qind]);
		      chan = Ind_Chan1(tb);
		      key1 = Chan.hh_map[chan][Hash(pind, qind, 0)];
		      key2 = Chan.pp_map[chan][Hash(rind, sind, 0)];
		      chan_ind = Ints.Vhhpp.V_1_index[chan];
		      ind = key1 * Chan.npp[chan] + key2;
		      Ints.Vhhpp.V_1[chan_ind + ind] += TBME;
		      /*chan_ind = Ints.Vpphh.V_1_index[chan];
		      ind = key2 * Chan.nhh[chan] + key1;
		      Ints.Vpphh.V_1[chan_ind + ind] += TBME;*/

		      chan = Ind_Chan3(SPB.qnums[pind]);
		      key1 = Chan.h_map[chan][pind];
		      key2 = Chan.pph_map[chan][Hash(rind, sind, qind, 0)];
		      chan_ind = Ints.Vhhpp.V_3_1_index[chan];
		      ind = key1 * Chan.npph[chan] + key2;
		      Ints.Vhhpp.V_3_1[chan_ind + ind] += TBME;
	    
		      chan = Ind_Chan3(SPB.qnums[qind]);
		      key1 = Chan.h_map[chan][qind];
		      key2 = Chan.pph_map[chan][Hash(rind, sind, pind, 0)];
		      chan_ind = Ints.Vhhpp.V_3_2_index[chan];
		      ind = key1 * Chan.npph[chan] + key2;
		      Ints.Vhhpp.V_3_2[chan_ind + ind] += TBME;
		      
		      chan = Ind_Chan3(SPB.qnums[rind]);
		      key1 = Chan.hhp_map[chan][Hash(pind, qind, sind, 0)];
		      key2 = Chan.p_map[chan][rind];
		      chan_ind = Ints.Vhhpp.V_3_3_index[chan];
		      ind = key1 * Chan.np[chan] + key2;
		      Ints.Vhhpp.V_3_3[chan_ind + ind] += TBME;
		      
		      minus(tb, SPB.qnums[pind], SPB.qnums[sind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.hp1_map[chan][Hash(pind, sind, 0)];
		      key2 = Chan.ph1_map[chan][Hash(rind, qind, 0)];
		      chan_ind = Ints.Vhhpp.V_2_1_index[chan];
		      ind = key1 * Chan.nph1[chan] + key2;
		      Ints.Vhhpp.V_2_1[chan_ind + ind] += TBME;
		      
		      minus(tb, SPB.qnums[pind], SPB.qnums[rind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.hp1_map[chan][Hash(pind, rind, 0)];
		      key2 = Chan.ph1_map[chan][Hash(sind, qind, 0)];
		      chan_ind = Ints.Vhhpp.V_2_3_index[chan];
		      ind = key1 * Chan.nph1[chan] + key2;
		      Ints.Vhhpp.V_2_3[chan_ind + ind] += TBME;
		    } 
		  }
		  else if(rtype == 0 && stype == 1){
		    if(ptype == 0 && qtype == 1){
		      minus(tb, SPB.qnums[pind], SPB.qnums[sind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.hp1_map[chan][Hash(pind, sind, 0)];
		      key2 = Chan.hp1_map[chan][Hash(rind, qind, 0)];
		      chan_ind = Ints.Vhphp.V_2_1_index[chan];
		      ind = key1 * Chan.nhp1[chan] + key2;
		      Ints.Vhphp.V_2_1[chan_ind + ind] += TBME;
		      
		      minus(tb, SPB.qnums[qind], SPB.qnums[rind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.ph1_map[chan][Hash(qind, rind, 0)];
		      key2 = Chan.ph1_map[chan][Hash(sind, pind, 0)];
		      chan_ind = Ints.Vhphp.V_2_2_index[chan];
		      ind = key1 * Chan.nph1[chan] + key2;
		      Ints.Vhphp.V_2_2[chan_ind + ind] += TBME;
		    }
		    else if(ptype == 0 && qtype == 0){
		      chan = Ind_Chan3(SPB.qnums[qind]);
		      key1 = Chan.h_map[chan][qind];
		      key2 = Chan.hph_map[chan][Hash(rind, sind, pind, 0)];
		      chan_ind = Ints.Vhhhp.V_3_2_index[chan];
		      ind = key1 * Chan.nhph[chan] + key2;
		      Ints.Vhhhp.V_3_2[chan_ind + ind] += TBME;

		      chan = Ind_Chan3(SPB.qnums[rind]);
		      key1 = Chan.hhp_map[chan][Hash(pind, qind, sind, 0)];
		      key2 = Chan.h_map[chan][rind];
		      chan_ind = Ints.Vhhhp.V_3_3_index[chan];
		      ind = key1 * Chan.nh[chan] + key2;
		      Ints.Vhhhp.V_3_3[chan_ind + ind] += TBME;
		      
		      minus(tb, SPB.qnums[pind], SPB.qnums[rind]);
		      chan = Ind_Chan2(tb);
		      key1 = Chan.hh1_map[chan][Hash(pind, rind, 0)];
		      key2 = Chan.ph1_map[chan][Hash(sind, qind, 0)];
		      chan_ind = Ints.Vhhhp.V_2_3_index[chan];
		      ind = key1 * Chan.nph1[chan] + key2;
		      Ints.Vhhhp.V_2_3[chan_ind + ind] += TBME;

		      // Vhphh -> rspq
		      chan = Ind_Chan3(SPB.qnums[sind]);
		      key1 = Chan.p_map[chan][sind];
		      key2 = Chan.hhh_map[chan][Hash(pind, qind, rind, 0)];
		      chan_ind = Ints.Vhphh.V_3_2_index[chan];
		      ind = key1 * Chan.nhhh[chan] + key2;
		      Ints.Vhphh.V_3_2[chan_ind + ind] += TBME;
		    }
		  }
		  else if(rtype == 0 && stype == 0){
		    if(ptype == 0 && qtype == 0){
		      plus(tb, SPB.qnums[pind], SPB.qnums[qind]);
		      chan = Ind_Chan1(tb);
		      key1 = Chan.hh_map[chan][Hash(pind, qind, 0)];
		      key2 = Chan.hh_map[chan][Hash(rind, sind, 0)];
		      chan_ind = Ints.Vhhhh.V_1_index[chan];
		      ind = key1 * Chan.nhh[chan] + key2;
		      Ints.Vhhhh.V_1[chan_ind + ind] += TBME;
		      
		      chan = Ind_Chan3(SPB.qnums[qind]);
		      key1 = Chan.h_map[chan][qind];
		      key2 = Chan.hhh_map[chan][Hash(rind, sind, pind, 0)];
		      chan_ind = Ints.Vhhhh.V_3_2_index[chan];
		      ind = key1 * Chan.nhhh[chan] + key2;
		      Ints.Vhhhh.V_3_2[chan_ind + ind] += TBME;
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
}

void Darmstadt_Setup(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME)
{
  int eMax = 1000, nMax = 1000, lMax = 1000, EMax = 1000, hwHO = 1000;
  std::string PATH0 = "/mnt/research/imsrg/nsuite/me/";
  std::string fullpath = PATH0 + PAR.MatrixElements;
  int stringlength = PAR.MatrixElements.length();
  std::string substring, fileend, eMaxstring;
  int filetype = 1; // type .gz
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  double scale;

  for(int i = 0; i < stringlength-2; ++i){
    if(i == stringlength-3){
      substring = PAR.MatrixElements.substr(i,3);
      fileend = substring; // bin or .gz
      if( fileend == "bin" ){ filetype = 0; }
      else if( fileend == ".gz" ){ filetype = 1; }
    }
    else{
      substring = PAR.MatrixElements.substr(i,4);
      if(substring == "eMax"){
	i += 4;
	eMax = std::stoi(PAR.MatrixElements.substr(i,2));
	eMaxstring = PAR.MatrixElements.substr(i,2);
	i += 1;
      }
      else if(substring == "lMax" || substring == "LMax"){
	i += 4;
	lMax = std::stoi(PAR.MatrixElements.substr(i,2));
	i += 1;
      }
      else if(substring == "EMax"){
	i += 4;
	EMax = std::stoi(PAR.MatrixElements.substr(i,2));
	i += 1;
      }
      else if(substring == "hwHO"){
	i += 4;
	hwHO = std::stoi(PAR.MatrixElements.substr(i,3));
	i += 2;
      }
    }
  }
  SPB.eMax = eMax;
  SPB.nMax = (eMax - eMax%2)/2;
  SPB.EMax = (EMax < 2*eMax ? EMax : 2*eMax);
  SPB.lMax = (lMax < eMax ? lMax : eMax);
  PAR.ho_energy = double(hwHO);
  PAR.ho_length = std::sqrt(2.0 * INVM / PAR.ho_energy); // 20.7355285386 = hbarc^2/(2 MN)  [MeV fm^2]
  //double COM_ho_length = std::sqrt(2.0 * INVM / COM_hw);
  std::cout << "eMax = " << SPB.eMax << ", EMax = " << SPB.EMax << ", nMax = " << SPB.nMax << ", lMax = " << SPB.lMax << ", hwHO = " << PAR.ho_energy << std::endl;
  std::cout << "HO_hw = " << PAR.ho_energy << ", HO_a = " << PAR.ho_length << " !! " << INVM/(PAR.ho_length*PAR.ho_length) << std::endl;

  //*************************//
  //SPB.eMax = 2;

  //*** determine nljMax, count states
  int lMin, n, twojMin, twojMax;
  SPB.nljMax = 0;
  for(int e0 = 0; e0 <= SPB.eMax; e0++){
    lMin = e0%2;
    for(int l = lMin; l <= std::min(e0,lMax); l += 2){
      n = (e0-l)/2;
      if(n>nMax) continue;
      twojMin = std::abs(2*l - 1);
      twojMax = 2*l + 1;
      for(int twoj = twojMin; twoj <= twojMax; twoj += 2){
	//*** n, l and j determine nlj
	SPB.nljMax++;
      }
    }
  }
  //*** allocate lookups
  SPB.nljDim = SPB.nljMax;
  SPB.e_nlj = new int[SPB.nljDim];
  SPB.n_nlj = new int[SPB.nljDim];
  SPB.l_nlj = new int[SPB.nljDim];
  SPB.twoj_nlj = new int[SPB.nljDim];

  SPB.Build_Darmstadt();
  Print_Parameters();
  HF_Chan.Build();
  HF_ME.Build(HF_Chan);

  // Count ME combinations
  int JMin, JMax, MEcount = 0;
  //*** bra loops
  for(int nlj1 = 0; nlj1 < SPB.nljMax; nlj1++){
    for(int nlj2 = 0; nlj2 <= nlj1; nlj2++){ // nlj1 >= nlj2
      if(SPB.e_nlj[nlj1] + SPB.e_nlj[nlj2] > EMax){ break; } // e1+e2 <= EMax
      //*** ket loops
      for(int nnlj1 = 0; nnlj1 <= nlj1; nnlj1++){ // nlj1*nljDim + nlj2 >= nnlj1*nljDim + nnlj2
        for(int nnlj2 = 0; nnlj2 <= (nnlj1 == nlj1 ? nlj2 : nnlj1); nnlj2++){ // nnlj1 >= nnlj2 and nlj1*nljDim + nlj2 >= nnlj1*nljDim + nnlj2
	  if( SPB.e_nlj[nnlj1] + SPB.e_nlj[nnlj2] > SPB.EMax){ break; } // ee1+ee2 <= EMax
	  if( (SPB.l_nlj[nlj1] + SPB.l_nlj[nlj2] - SPB.l_nlj[nnlj1] - SPB.l_nlj[nnlj2])%2 ){ continue; } // parity conservation
	  // triangular condition
	  JMin = std::max( std::abs(SPB.twoj_nlj[nlj1] - SPB.twoj_nlj[nlj2]), std::abs(SPB.twoj_nlj[nnlj1] - SPB.twoj_nlj[nnlj2]) );
          JMax = std::min( (SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2]), (SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2]) );
          if(JMin > JMax){ continue; }
          //*** diagonal loops
	  for(int J = JMin; J <= JMax; J += 2){
	    for(int T = 0; T <= 2; T += 2){
              for(int MT = -T; MT <= T; MT += 2){
		MEcount++; // counter
	      }
	    }
          }
        }
      }
    }
  }
  double *ME0 = new double[MEcount];
  scale = 1.0;
  //scale = 0.0;
  Darmstadt_Read(fullpath, ME0, MEcount, filetype);
  Darmstadt_Fill(HF_Chan, HF_ME, ME0, scale);

  std::string fullpath_com1;
  std::string fullpath_com2;
  // pi.pj CoM Correction for intrinsic Hamiltonian
  fullpath_com1 = PATH0 + "tpp_eMax" + eMaxstring + ".me2j.gz";
  scale = 2.0 / (PAR.A * PAR.ho_length * PAR.ho_length);
  //scale = 0.0;
  Darmstadt_Read(fullpath_com1, ME0, MEcount, 1);
  Darmstadt_Fill(HF_Chan, HF_ME, ME0, scale);

  /////////////////////////////////////////////
  /*HF_Matrix_Elements HF_ME2;
  HF_ME2.Build(HF_Chan);
  std::string fullpath_heiko = "sam_com/ho_reference/He4_hcm_eMax06_hwHO020.ham0.me2b.gz";
  scale = 1.0;
  Darmstadt_Read(fullpath_heiko, ME0, MEcount, 1);
  Darmstadt_Fill(HF_Chan, HF_ME2, ME0, scale);*/
  ////////////////////////////////////////////

  // pi.pj and ri.rj for CoM Hamiltonian
  if(PAR.CoM == 2){
    fullpath_com1 = PATH0 + "tpp_eMax" + eMaxstring + ".me2j.gz";
    fullpath_com2 = PATH0 + "r1r2_eMax" + eMaxstring + ".me2j.gz";
    scale = -2.0 * PAR.CoM_fac / (PAR.A * PAR.ho_length * PAR.ho_length);
    Darmstadt_Read(fullpath_com1, ME0, MEcount, 1);
    Darmstadt_Fill(HF_Chan, HF_ME, ME0, scale);
    scale = PAR.CoM_fac * PAR.ho_energy / PAR.A;
    Darmstadt_Read(fullpath_com2, ME0, MEcount, 1);
    Darmstadt_Fill(HF_Chan, HF_ME, ME0, scale);
  }

  /*for(int chan1 = 0; chan1 < HF_Chan.size1; ++chan1){
    for(int ind1 = 0; ind1 < HF_Chan.ntb[chan1]; ++ind1){
      for(int ind2 = ind1; ind2 < HF_Chan.ntb[chan1]; ++ind2){
	if( std::fabs(HF_ME.V[HF_ME.Index[chan1] + (ind1*HF_Chan.ntb[chan1] + ind2)]) > 1.0e-6 ){
	  std::cout << "!!  " << HF_ME.V[HF_ME.Index[chan1] + (ind1*HF_Chan.ntb[chan1] + ind2)] << std::endl;
	}
      }
    }
  }
  exit(1);*/

  delete[] ME0;
}

void Darmstadt_Read(std::string filepath, double *ME0, int count, int type)
{
  int ME2J_HEADERSIZE = 255;
  int BUFSIZE = ME2J_HEADERSIZE*sizeof(char);
  gzFile gzfd;
  FILE *fd;
  char buf[BUFSIZE], *d;

  if( type == 0 ){
    // open file
    fd = fopen(filepath.c_str(), "rb");
    if( fd == NULL ){ std::cerr << "Matrix Element file, " << filepath << ", does not exist" << std::endl; exit(1); }
    // read and check header
    fread(buf, sizeof(char), ME2J_HEADERSIZE, fd);
    if( std::strstr(buf,"me2j-f2-bin") == NULL ){ std::cerr << "Matrix Element file, " << filepath << " is wrong format" << std::endl; exit(1); }
    // read data
    fread(ME0, sizeof(double), count, fd);
    // close file
    fclose(fd);
  }
  else if( type == 1 ){
    // open file
    gzfd = gzopen(filepath.c_str(), "r");
    if( !gzfd ){ std::cerr << "Matrix Element file, " << filepath << ", does not exist" << std::endl; exit(1); }
    // read and check format identifier
    gzgets(gzfd, buf, BUFSIZE);
    if( std::strstr(buf,"me2j-f2") == NULL ){ std::cerr << "Matrix Element file, " << filepath << " is wrong format" << std::endl; exit(1); }
    // read matrix elements
    for(int i = 0; i < count; i += 10){
      gzgets(gzfd, buf, BUFSIZE);
      d = std::strtok(buf, " ");
      ME0[i] = std::atof(d);

      for(int i0 = i+1; i0 < std::min(i+10, count); ++i0){
	d = std::strtok(NULL, " ");
	ME0[i0] = std::atof(d);
      }
    }
    // close file
    gzclose(gzfd);
  }
  else{ std::cerr << "Wrong Darmstadt Matrix Element type" << std::endl; exit(1); }
}

void Darmstadt_Fill(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME, double *ME0, double factor)
{
  State statep, stateq, stater, states, tb;
  int p, q, r, s, key1, key2, key3, key4, key, ind, chan1, JMin, JMax;
  double TBME;
  int count = 0;
  //*** bra loops
  for(int nlj1 = 0; nlj1 < SPB.nljMax; nlj1++){
    statep.n = SPB.n_nlj[nlj1];
    statep.l = SPB.l_nlj[nlj1];
    statep.par = -2*(SPB.l_nlj[nlj1]%2) + 1;
    statep.j = SPB.twoj_nlj[nlj1];
    for(int nlj2 = 0; nlj2 <= nlj1; nlj2++){ // nlj1 >= nlj2
      stateq.n = SPB.n_nlj[nlj2];
      stateq.l = SPB.l_nlj[nlj2];
      stateq.par = -2*(SPB.l_nlj[nlj2]%2) + 1;
      stateq.j = SPB.twoj_nlj[nlj2];
      if(SPB.e_nlj[nlj1] + SPB.e_nlj[nlj2] > SPB.EMax){ break; } // e1+e2 <= EMax
      //*** ket loops
      for(int nnlj1 = 0; nnlj1 <= nlj1; nnlj1++){ // nlj1*nljDim + nlj2 >= nnlj1*nljDim + nnlj2
	stater.n = SPB.n_nlj[nnlj1];
	stater.l = SPB.l_nlj[nnlj1];
	stater.par = -2*(SPB.l_nlj[nnlj1]%2) + 1;
	stater.j = SPB.twoj_nlj[nnlj1];
	for(int nnlj2 = 0; nnlj2 <= (nnlj1 == nlj1 ? nlj2 : nnlj1); nnlj2++){ // nnlj1 >= nnlj2 and nlj1*nljDim + nlj2 >= nnlj1*nljDim + nnlj2
	  states.n = SPB.n_nlj[nnlj2];
	  states.l = SPB.l_nlj[nnlj2];
	  states.par = -2*(SPB.l_nlj[nnlj2]%2) + 1;
	  states.j = SPB.twoj_nlj[nnlj2];
	  if( SPB.e_nlj[nnlj1] + SPB.e_nlj[nnlj2] > SPB.EMax ){ break; } // ee1+ee2 <= EMax
	  if( (SPB.l_nlj[nlj1] + SPB.l_nlj[nlj2] - SPB.l_nlj[nnlj1] - SPB.l_nlj[nnlj2])%2 ){ continue; } // parity conservation
	  // triangular condition
	  JMin = std::max( std::abs(SPB.twoj_nlj[nlj1] - SPB.twoj_nlj[nlj2]), std::abs(SPB.twoj_nlj[nnlj1] - SPB.twoj_nlj[nnlj2]) );
	  JMax = std::min( SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2], SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2] );
	  if(JMin > JMax){ continue; }
	  //*** diagonal loops
	  for(int J = JMin; J <= JMax; J += 2){
	    for(int T = 0; T <= 2; T += 2){
	      for(int MT = -T; MT <= T; MT += 2){
		/*if(std::fabs(ME0[count]) > 1e-10){
		  std::cout << "T = " << T << "(" << MT << "), J = " << J << " : < " << nlj1 << " " << nlj2 << " | " << nnlj1 << " " << nnlj2 << " > = " << ME0[count] << std::endl;
 		  }*/
 		for(statep.t = -1; statep.t <= 1; statep.t += 2){
		  key = Ind_State(statep);
		  p = SPB.map_state[key];
		  stateq.t = -MT - statep.t; // Darmstadt uses t(p) = 1 and t(n) = -1
 		  if( std::abs(stateq.t) > 1 ){ continue; }
		  key = Ind_State(stateq);
		  q = SPB.map_state[key];
 		  if(p == q && J%4 != 0){ continue; }
		  for(stater.t = -1; stater.t <= 1; stater.t += 2){
		    key = Ind_State(stater);
		    r = SPB.map_state[key];
		    states.t = -MT - stater.t; // Darmstadt uses t(p) = 1 and t(n) = -1
		    if( std::abs(states.t) > 1 ){ continue; }
		    key = Ind_State(states);
		    s = SPB.map_state[key];
		    if(r == s && J%4 != 0){ continue; }
		    
		    TBME = PAR.tbstrength * factor * ME0[count] * CGC(1, statep.t, 1, stateq.t, T, -MT) * CGC(1, stater.t, 1, states.t, T, -MT);
		    /*if(std::fabs(TBME) > 1e-10){
		      std::cout << "T = " << T << "(" << MT << "), J = " << J << " : < " << p << " " << q << " | " << r << " " << s << " > = ";
		      std::cout << TBME << "   " << CGC(1, statep.t, 1, stateq.t, T, MT) << " " << CGC(1, stater.t, 1, states.t, T, MT) << std::endl;
		      }*/

		    plus(tb, SPB.qnums[p], SPB.qnums[q]);
		    tb.j = J;
		    chan1 = Ind_Chan1(tb);
		    key1 = HF_Chan.tb_map[chan1][Hash(p, q, tb.j)];
		    key2 = HF_Chan.tb_map[chan1][Hash(r, s, tb.j)];
		    key3 = HF_Chan.tb_map[chan1][Hash(q, p, tb.j)];
		    key4 = HF_Chan.tb_map[chan1][Hash(s, r, tb.j)];
		    
		    ind = HF_ME.Index[chan1] + (key1 * HF_Chan.ntb[chan1] + key2);
		    HF_ME.V[ind] += TBME;
		    if( nlj1 != nlj2 ){
		      ind = HF_ME.Index[chan1] + (key3 * HF_Chan.ntb[chan1] + key2);
		      HF_ME.V[ind] += -1.0 * phase2(SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2] - J) * TBME;
		    }
		    if( nnlj1 != nnlj2 ){
		      ind = HF_ME.Index[chan1] + (key1 * HF_Chan.ntb[chan1] + key4);
		      HF_ME.V[ind] += -1.0 * phase2(SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2] - J) * TBME;
		    }
		    if( nlj1 != nlj2 && nnlj1 != nnlj2 ){
		      ind = HF_ME.Index[chan1] + (key3 * HF_Chan.ntb[chan1] + key4);
		      HF_ME.V[ind] += phase2(SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2] + SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2]) * TBME;
		    }
 		    if( nlj1 != nnlj1 || nlj2 != nnlj2 ){
		      ind = HF_ME.Index[chan1] + (key2 * HF_Chan.ntb[chan1] + key1);
		      HF_ME.V[ind] += TBME;
		      if( nnlj1 != nnlj2 ){
 			ind = HF_ME.Index[chan1] + (key4 * HF_Chan.ntb[chan1] + key1);
			HF_ME.V[ind] += -1.0 * phase2(SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2] - J) * TBME;
		      }
		      if( nlj1 != nlj2 ){
 			ind = HF_ME.Index[chan1] + (key2 * HF_Chan.ntb[chan1] + key3);
			HF_ME.V[ind] += -1.0 * phase2(SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2] - J) * TBME;
		      }
		      if( nlj1 != nlj2 && nnlj1 != nnlj2 ){
 			ind = HF_ME.Index[chan1] + (key4 * HF_Chan.ntb[chan1] + key3);
			HF_ME.V[ind] += phase2(SPB.twoj_nlj[nlj1] + SPB.twoj_nlj[nlj2] + SPB.twoj_nlj[nnlj1] + SPB.twoj_nlj[nnlj2]) * TBME;
		      }
		    }
		  }
 		}
		// counter
		count++;
	      }
	    }
	  }	  
	}
      }
    }
  }
}

// one-body kinetic energy = p^2/2m
double Kinetic_Energy(int p, int q)
{
  double obcm;
  if( PAR.CoM > 0 ){ obcm = 1.0 - (1.0/PAR.A); }
  else{ obcm = 1.0; }

  int n1 = SPB.qnums[p].n;
  int n2 = SPB.qnums[q].n;
  //std::cout << "??  " << n1 << " " << n2 << " " << SPB.qnums[p].l << " " << SPB.qnums[q].l << " " << PAR.ho_energy << std::endl;
  if( SPB.qnums[p].l != SPB.qnums[q].l || SPB.qnums[p].j != SPB.qnums[q].j || SPB.qnums[p].t != SPB.qnums[q].t ){ return 0.0; }
  if(n1 == n2){ return obcm * 0.5 * PAR.ho_energy * (2*n1 + SPB.qnums[p].l + 1.5); }
  if(n1 == n2-1){ return obcm * 0.5 * PAR.ho_energy * std::sqrt((n1 + 1)*(n1 + SPB.qnums[p].l + 1.5)); }
  if(n1 == n2+1){ return obcm * 0.5 * PAR.ho_energy * std::sqrt((n2 + 1)*(n2 + SPB.qnums[p].l + 1.5)); }
  else{ return 0.0; }
}

// one-body part of CoM Hamiltonian = p^2/2mA + mw^2r^2/2A
double Hcom_1_Body(int p, int q, double hwBar)
{
  // we get a^2 from the r^2 m.e.:
  // 0.5 *m*wBar^2 * a^2 = 0.5 * m * wBar^2 * hbar^2/(m * hw) = 0.5 * hwBar^2/hw
  int n1 = SPB.qnums[p].n;
  int n2 = SPB.qnums[q].n;
  //std::cout << "??  " << n1 << " " << n2 << " " << SPB.qnums[p].l << " " << SPB.qnums[q].l << " " << PAR.ho_energy << " " << hwBar << " " << PAR.A << " " << PAR.CoM_fac << std::endl;
  if( SPB.qnums[p].l != SPB.qnums[q].l || SPB.qnums[p].j != SPB.qnums[q].j || SPB.qnums[p].t != SPB.qnums[q].t ){ return 0.0; }
  if(n1 == n2){ return PAR.CoM_fac * 0.5 * (PAR.ho_energy + hwBar*hwBar/PAR.ho_energy)/PAR.A * (2*n1 + SPB.qnums[p].l + 1.5); }
  if(n1 == n2-1){ return PAR.CoM_fac * 0.5 * (PAR.ho_energy - hwBar*hwBar/PAR.ho_energy)/PAR.A * std::sqrt((n1 + 1)*(n1 + SPB.qnums[p].l + 1.5)); }
  if(n1 == n2+1){ return PAR.CoM_fac * 0.5 * (PAR.ho_energy - hwBar*hwBar/PAR.ho_energy)/PAR.A * std::sqrt((n2 + 1)*(n2 + SPB.qnums[p].l + 1.5)); }
  else{ return 0.0; }
}

void Darmstadt3_Setup(HF_Channels &HF_Chan, HF_Matrix_Elements &HF_ME)
{
  int eMax = 1000, nMax = 1000, lMax = 1000, EMax = 1000, hwHO = 1000;
  std::string PATH0 = "/mnt/research/imsrg/nsuite/me/";
  std::string fullpath = PATH0 + PAR.MatrixElements3;
  int stringlength = PAR.MatrixElements3.length();
  std::string substring, fileend, eMaxstring;
  int filetype = 1; // type .gz
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  
  for(int i = 0; i < stringlength-2; ++i){
    if(i == stringlength-3){
      substring = PAR.MatrixElements3.substr(i,3);
      fileend = substring; // bin or .gz
      if( fileend == "bin" ){ filetype = 0; }
      else if( fileend == ".gz" ){ filetype = 1; }
    }
    else{
      substring = PAR.MatrixElements3.substr(i,4);
      if(substring == "eMax"){
	i += 4;
 	eMax = std::stoi(PAR.MatrixElements3.substr(i,2));
 	eMaxstring = PAR.MatrixElements3.substr(i,2);
 	i += 1;
      }
      else if(substring == "lMax" || substring == "LMax"){
	i += 4;
 	lMax = std::stoi(PAR.MatrixElements3.substr(i,2));
 	i += 1;
      }
      else if(substring == "EMax"){
	i += 4;
	EMax = std::stoi(PAR.MatrixElements3.substr(i,2));
	i += 1;
      }
      else if(substring == "nMax"){
	i += 4;
	nMax = std::stoi(PAR.MatrixElements3.substr(i,2));
	i += 1;
      }
      else if(substring == "hwHO"){
	i += 4;
	hwHO = std::stoi(PAR.MatrixElements3.substr(i,3));
	i += 2;
      }
    }
  }
  SPB.e3Max = eMax;
  SPB.n3Max = ((eMax - eMax%2)/2 < nMax ? (eMax - eMax%2)/2 : nMax);
  SPB.E3Max = (EMax < 3*eMax ? EMax : 3*eMax);
  SPB.l3Max = (lMax < eMax ? lMax : eMax);
  if( double(hwHO) != PAR.ho_energy ){ std::cerr << "ME-3N HO hw doesn't match, " << double(hwHO) << ", " << PAR.ho_energy << std::endl; exit(1); }

  std::cout << "eMax3 = " << SPB.e3Max << ", EMax3 = " << SPB.E3Max << ", nMax3 = " << SPB.n3Max << ", lMax3 = " << SPB.l3Max << ", hwHO3 = " << PAR.ho_energy << std::endl;
  std::cout << "HO_hw3 = " << PAR.ho_energy << ", HO_a3 = " << PAR.ho_length << " !! " << INVM/(PAR.ho_length*PAR.ho_length) << std::endl;

  //*** determine nljMax, count states
  int lMin, n, twojMin, twojMax;
  SPB.nljMax3 = 0;
  for(int e0 = 0; e0 <= SPB.e3Max; e0++){
    lMin = e0%2;
    for(int l = lMin; l <= std::min(e0,SPB.l3Max); l += 2){
      n = (e0-l)/2;
      if(n > SPB.n3Max){ continue; }
      twojMin = std::abs(2*l - 1);
      twojMax = 2*l + 1;
      for(int twoj = twojMin; twoj <= twojMax; twoj += 2){
	//*** n, l and j determine nlj
	SPB.nljMax3++;
      }
    }
  }
  //*** allocate lookups
  SPB.nljDim3 = SPB.nljMax3;
  SPB.e_nlj3 = new int[SPB.nljDim3];
  SPB.n_nlj3 = new int[SPB.nljDim3];
  SPB.l_nlj3 = new int[SPB.nljDim3];
  SPB.twoj_nlj3 = new int[SPB.nljDim3];

  //*** construct nlj index ******************************************
  int nlj = 0, ind = 0;
  for(int e0 = 0; e0 <= SPB.e3Max; e0++){
    lMin = e0%2;
    for(int l = lMin; l <= std::min(e0,SPB.l3Max); l += 2){
      n = (e0 - l)/2;
      if(n > SPB.n3Max){ continue; }
      twojMin = std::abs(2*l - 1);
      twojMax = 2*l + 1;
      for(int twoj = twojMin; twoj <= twojMax; twoj += 2){
        //*** n, l and j determine nlj
        SPB.e_nlj3[nlj]     = e0;
        SPB.n_nlj3[nlj]     = n;
        SPB.l_nlj3[nlj]     = l;
        SPB.twoj_nlj3[nlj]  = twoj;
	++nlj;
      }
    }
  }

  // reconcile existent basis from 2BME with 3BME
  // -- qnums.nlj3[ind] = -1 if 2B basis contained within 3B basis
  ind = 0;
  nlj = 0;
  for(int e0 = 0; e0 <= SPB.eMax; e0++){
    lMin = e0%2;
    for(int l = lMin; l <= std::min(e0,SPB.lMax); l += 2){
      n = (e0 - l)/2;
      if(n > SPB.nMax){ continue; }
      twojMin = std::abs(2*l - 1);
      twojMax = 2*l + 1;
      for(int twoj = twojMin; twoj <= twojMax; twoj += 2){
	for(int t = -1; t <= 1; t += 2){
	  if( e0 > SPB.e3Max || l > SPB.l3Max || n > SPB.n3Max ){ SPB.qnums[ind].nlj3 = -1; }
	  else{ SPB.qnums[ind].nlj3 = nlj; }
	  ++ind;
	}
	++nlj;
      }
    }
  }

  /*std::cout << "Basis Set..." << std::endl;
  for(int ind = 0; ind < SPB.num_states; ++ind){
    std::cout << ind << ": " << SPB.qnums[ind].nlj << ", " << SPB.qnums[ind].nlj3 << std::endl;
  }
  std::cout << std::endl;*/

  int 	nlj1, nlj2, nlj3, nnlj1, nnlj2, nnlj3, nnlj3Max;
  int 	nlj1count, nlj2count, nlj3count, nnlj1count, nnlj2count, nnlj3count;
  int 	Jab2, JabMin2, JabMax2, JJab2, JJabMin2, JJabMax2, twoJC, twoJCMin, twoJCMax;
  long long	MEcount, Dim;

  Dim = 0;
  nlj1count = 0;
  for(nlj1 = 0; nlj1 < SPB.nljMax3; nlj1++){
    if( SPB.e_nlj3[nlj1] > SPB.E3Max){ break; }
    nlj1count++;
  }
  Dim += nlj1count;

  //std::cout << "## nlj1count = " << nlj1count << std::endl;
  if(nlj1count != 0){ HF_ME.ME3Idx = new long long*****[nlj1count]; }
  if(HF_ME.ME3Idx == NULL){ std::cerr << "ME3Idx could not be allocated..." << std::endl; exit(1); }

  for(nlj1 = 0; nlj1 < SPB.nljMax3; nlj1++){
    if( SPB.e_nlj3[nlj1] > SPB.E3Max){ break; }
    nlj2count = 0;
    for(nlj2 = 0; nlj2 <= nlj1; nlj2++){
      if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] > SPB.E3Max){ break; }
      nlj2count++;
    }
    Dim += nlj2count;

    //std::cout << "## nlj2count = " << nlj2count << std::endl;
    if(nlj2count != 0){ HF_ME.ME3Idx[nlj1] = new long long****[nlj2count]; }
    if(HF_ME.ME3Idx[nlj1] == NULL){ std::cerr << "ME3Idx[nlj1] could not be allocated..." << std::endl; exit(1); }

    for(nlj2 = 0; nlj2 <= nlj1; nlj2++){
      if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] > SPB.E3Max){ break; }
      nlj3count = 0;
      for(nlj3 = 0; nlj3 <= nlj2; nlj3++){
        if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] + SPB.e_nlj3[nlj3] > SPB.E3Max){ break; }
        nlj3count++;
      }
      Dim += nlj3count;

      //std::cout << "### nlj3count = " << nlj3count << std::endl;
      if(nlj3count != 0){ HF_ME.ME3Idx[nlj1][nlj2] = new long long***[nlj3count]; }
      if(HF_ME.ME3Idx[nlj1][nlj2] == NULL){ std::cerr << "ME3Idx[nlj1][nlj2] could not be allocated..." << std::endl; exit(1); }

      for(nlj3 = 0; nlj3 <= nlj2; nlj3++){
        if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] + SPB.e_nlj3[nlj3] > SPB.E3Max){ break; }
        nnlj1count = 0;
        for(nnlj1 = 0; nnlj1 <= nlj1; nnlj1++){
          if( SPB.e_nlj3[nnlj1] > SPB.E3Max){ break; } // could be neglected
          nnlj1count++;
        }
        Dim += nnlj1count;

	//std::cout << "###@ nnlj1count = " << nnlj1count << std::endl;
	if(nnlj1count != 0){ HF_ME.ME3Idx[nlj1][nlj2][nlj3] = new long long**[nnlj1count]; }
	if(HF_ME.ME3Idx[nlj1][nlj2][nlj3] == NULL){ std::cerr << "ME3Idx[nlj1][nlj2][nlj3] could not be allocated..." << std::endl; exit(1); }

        for(nnlj1 = 0; nnlj1 <= nlj1; nnlj1++){
          if( SPB.e_nlj3[nnlj1] > SPB.E3Max){ break; } // could be neglected
          nnlj2count = 0;
          for(nnlj2 = 0; nnlj2 <= ((nlj1==nnlj1) ? nlj2 : nnlj1); nnlj2++){
            if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] > SPB.E3Max ){ break; }
            nnlj2count++;
          }
          Dim += nnlj2count;

	  //std::cout << "###@@ nnlj2count = " << nnlj2count << std::endl;
	  if(nnlj2count != 0){ HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1] = new long long*[nnlj2count]; }
	  if(HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1] == NULL){ std::cerr << "ME3Idx[nlj1][nlj2][nlj3][nnlj1] could not be allocated..." << std::endl; exit(1); }

          for(nnlj2 = 0; nnlj2 <= ((nlj1==nnlj1) ? nlj2 : nnlj1); nnlj2++){
            if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] > SPB.E3Max){ break; }
            nnlj3count = 0;
            if( nlj1 == nnlj1 && nlj2 == nnlj2 ){ nnlj3Max = nlj3; }
	    else{ nnlj3Max = nnlj2; }
            for(nnlj3 = 0; nnlj3 <= nnlj3Max; nnlj3++) {
              if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] + SPB.e_nlj3[nnlj3] > SPB.E3Max ){ break; }
              nnlj3count++;
            }
            Dim += nnlj3count;

	    //std::cout << "###@@@ nnlj3count = " << nnlj3count << std::endl;
	    if(nnlj3count != 0){ HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2] = new long long[nnlj3count]; }
	    if(HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2] == NULL){ std::cerr << "ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2] could not be allocated..." << std::endl; exit(1); }

          }//nnlj2
        }//nnlj1
      }//nlj3
    }//nlj2
  }//nlj1

  std::cout << "!!! Dim = " << Dim << std::endl;

  MEcount = 0;
  for(nlj1 = 0; nlj1 < SPB.nljMax3; nlj1++){
    if( SPB.e_nlj3[nlj1] > SPB.E3Max ){ break; }
    for(nlj2 = 0; nlj2 <= nlj1; nlj2++){
      if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] > SPB.E3Max){ break; }
      for(nlj3 = 0; nlj3 <= nlj2; nlj3++){
        if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] + SPB.e_nlj3[nlj3] > SPB.E3Max){ break; }
        for(nnlj1 = 0; nnlj1 <= nlj1; nnlj1++){
          if( SPB.e_nlj3[nnlj1] > SPB.E3Max){ break; }
          for(nnlj2 = 0; nnlj2 <= ((nlj1==nnlj1) ? nlj2 : nnlj1); nnlj2++){
            if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] > SPB.E3Max){ break; }
            if( nlj1 == nnlj1 && nlj2 == nnlj2 ){ nnlj3Max = nlj3; }
	    else{ nnlj3Max = nnlj2; }
            for(nnlj3 = 0; nnlj3 <= nnlj3Max; nnlj3++){
              if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] + SPB.e_nlj3[nnlj3] > SPB.E3Max){ break; }

              // set index to unphysical
              HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2][nnlj3] = -1;
              // check parity
              if( (SPB.l_nlj3[nlj1] + SPB.l_nlj3[nlj2] + SPB.l_nlj3[nlj3] - SPB.l_nlj3[nnlj1] - SPB.l_nlj3[nnlj2] - SPB.l_nlj3[nnlj3])%2){ continue; }
              // summation bounds
              JabMax2  = (SPB.twoj_nlj3[nlj1] + SPB.twoj_nlj3[nlj2]);
              JabMin2  = std::abs(SPB.twoj_nlj3[nlj1] - SPB.twoj_nlj3[nlj2]);
              JJabMax2 = (SPB.twoj_nlj3[nnlj1] + SPB.twoj_nlj3[nnlj2]);
              JJabMin2 = std::abs(SPB.twoj_nlj3[nnlj1] - SPB.twoj_nlj3[nnlj2]);
              if( JabMin2 > JabMax2 || JJabMin2 > JJabMax2 ){ continue; }

              //***set index
              HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2][nnlj3] = MEcount;
	      //std::cout << "***   " << nlj1 << ", " << nlj2 << ", " << nlj3 << " | " << nnlj1 << ", " << nnlj2 << ", " << nnlj3 << " -> " << MEcount << std::endl;

              // inner loops
              for(Jab2 = JabMin2; Jab2 <= JabMax2; Jab2 += 2){
                for(JJab2 = JJabMin2; JJab2 <= JJabMax2; JJab2 += 2){
                  //summation bounds
                  twoJCMin = std::max(std::abs(Jab2 - SPB.twoj_nlj3[nlj3]), std::abs(JJab2 - SPB.twoj_nlj3[nnlj3]));
                  twoJCMax = std::min(Jab2 + SPB.twoj_nlj3[nlj3], JJab2 + SPB.twoj_nlj3[nnlj3]);
		  
                  for(twoJC = twoJCMin; twoJC <= twoJCMax; twoJC += 2){
		    // the total isospin loop can be replaced by i+=5
		    MEcount += 5;
		    /*for(tab2 = 0; tab2 <= 2; tab2 += 2){
                      for(ttab2 = 0; ttab2 <= 2; ttab2 +=2){
                        //summation bounds
                        twoTMin = std::max( std::abs(tab2 - 1), std::abs(ttab2 - 1)); // twoTMin can just be used as 1
                        twoTMax = std::min( tab2 + 1, ttab2 + 1);
			
                        for(twoT = twoTMin; twoT <= twoTMax; twoT += 2){
                          ++MEcount;
			  //printf("--%i: j=%i, jj=%i, twoJ=%i, t=%i, tt=%i, twoT=%i\n", i, Jab, JJab, twoJC, tab, ttab, twoT);
                        }
                      }//ttab
		      }//tab*/
                  }//twoJ
                }//JJab
              }//Jab
            }//nnlj3
          }//nnlj2
        }//nnlj1
      }//nlj3
    }//nlj2
  }//nlj1

  std::cout << "!!! MEcount = " << MEcount << std::endl;

  std::ifstream interaction2;	// interaction file
  interaction2.open(fullpath.c_str(), std::ios::binary);
  interaction2.seekg(0, interaction2.end);
  long long length = interaction2.tellg();
  interaction2.close();
  std::cout << "!!! length = " << length << ", " << (length - 255*sizeof(char))/4 << std::endl;

  //ME3JH.iMax = i;
  HF_ME.V3 = new float[MEcount];
  Darmstadt3_Read(fullpath, HF_ME.V3, MEcount, filetype);
}

void Darmstadt3_Read(std::string filepath, float *ME0, long long count, int type)
{
  int ME3J_HEADERSIZE = 255;
  int BUFSIZE = ME3J_HEADERSIZE*sizeof(char);
  gzFile gzfd;
  FILE *fd;
  char buf[BUFSIZE], *d;

  if( type == 0 ){
    // open file
    fd = fopen(filepath.c_str(), "rb");
    if( fd == NULL ){ std::cerr << "Matrix Element file, " << filepath << ", does not exist" << std::endl; exit(1); }
    // read and check header
    fread(buf, sizeof(char), ME3J_HEADERSIZE, fd);
    if( std::strstr(buf,"me3j-f2-bin") == NULL ){ std::cerr << "Matrix Element file, " << filepath << " is wrong format" << std::endl; exit(1); }
    // read data
    fread(ME0, sizeof(float), count, fd);
    // close file
    fclose(fd);
  }
  else if( type == 1 ){
    // open file
    gzfd = gzopen(filepath.c_str(), "r");
    if( !gzfd ){ std::cerr << "Matrix Element file, " << filepath << ", does not exist" << std::endl; exit(1); }
    // read and check format identifier
    gzgets(gzfd, buf, BUFSIZE);
    if( std::strstr(buf,"me3j-f2") == NULL ){ std::cerr << "Matrix Element file, " << filepath << " is wrong format" << std::endl; exit(1); }
    // read matrix elements
    for(long long i = 0; i < count; i += 10){
      gzgets(gzfd, buf, BUFSIZE);
      d = std::strtok(buf, " ");
      ME0[i] = std::atof(d);

      for(long long i0 = i+1; i0 < std::min(i+10, count); ++i0){
	d = std::strtok(NULL, " ");
	ME0[i0] = std::atof(d);
      }
    }
    // close file
    gzclose(gzfd);
  }
  else{ std::cerr << "Wrong Darmstadt Matrix Element type" << std::endl; exit(1); }
}


//***************************************************************************
// GetME in pn format
//***************************************************************************
double ME3J_GetME_pn(HF_Matrix_Elements &HF_ME, int p, int q, int r, int s, int t, int u, int j, int jj, int J)
{
  int mt, mtt, MT, MTT;
  int tmin, tmax, ttmin, ttmax, Tmin, Tmax;
  double V, CGC1, CGC2, CGC3, CGC4;

  int p_nlj = SPB.qnums[p].nlj3;
  int q_nlj = SPB.qnums[q].nlj3;
  int r_nlj = SPB.qnums[r].nlj3;
  int s_nlj = SPB.qnums[s].nlj3;
  int t_nlj = SPB.qnums[t].nlj3;
  int u_nlj = SPB.qnums[u].nlj3;

  int p_t = SPB.qnums[p].t;
  int q_t = SPB.qnums[q].t;
  int r_t = SPB.qnums[r].t;
  int s_t = SPB.qnums[s].t;
  int t_t = SPB.qnums[t].t;
  int u_t = SPB.qnums[u].t;

  mt = p_t + q_t;
  tmin = std::abs(mt);
  tmax = 2;

  mtt = s_t + t_t;
  ttmin = std::abs(mtt);
  ttmax = 2;

  MT = mt + r_t;
  MTT = mtt + u_t;
  if( MT != MTT ){ return 0.0; }
  Tmin = std::abs(MT);
  Tmax = 3;
  
  V = 0.0;
  for(int t = tmax; t >= tmin; t -= 2){
    //std::cout << "get_CGC1( " << 1 << "," << p_t << "," << 1 << "," << q_t << "," << t << "," << mt << " ) = " << std::endl;
    //std::cout << JCOUP.get_CGC(1, p_t, 1, q_t, t, mt) << std::endl;
    CGC1 = JCOUP.get_CGC(1, p_t, 1, q_t, t, mt);
    if( CGC1 == 0.0 ){ continue; }
    for(int tt = ttmax; tt >= ttmin; tt -= 2){
      //std::cout << "get_CGC2( " << 1 << "," << s_t << "," << 1 << "," << t_t << "," << tt << "," << mtt << " ) = " << std::endl;
      //std::cout << JCOUP.get_CGC(1, s_t, 1, t_t, tt, mtt) << std::endl;
      CGC2 = JCOUP.get_CGC(1, s_t, 1, t_t, tt, mtt);
      if( CGC2 == 0.0 ){ continue; }
      for(int T = Tmax; T >= Tmin; T -= 2){
	//std::cout << "get_CGC3( " << t << "," << mt << "," << 1 << "," << r_t << "," << T << "," << MT << " ) = " << std::endl;
	//std::cout << JCOUP.get_CGC(t, mt, 1, r_t, T, MT) << std::endl;
	CGC3 = JCOUP.get_CGC(t, mt, 1, r_t, T, MT);
	if( CGC3 == 0.0 ){ continue; }
	//std::cout << "get_CGC4( " << tt << "," << mtt << "," << 1 << "," << u_t << "," << T << "," << MT << " ) = " << std::endl;
	//std::cout << JCOUP.get_CGC(tt, mtt, 1, u_t, T, MT) << std::endl;
	CGC4 = JCOUP.get_CGC(tt, mtt, 1, u_t, T, MT);
	if( CGC4 == 0.0 ){ continue; }
	V += CGC1 * CGC2 * CGC3 * CGC4 * ME3J_GetME(HF_ME, p_nlj, q_nlj, r_nlj, s_nlj, t_nlj, u_nlj, j, jj, J, t, tt, T);
      }
    }
  }

  return V;
}


//***************************************************************************
// GetME
//***************************************************************************
// get specific [(ja,jb)Jab,jc]J[(ta,tb)tab,tc]T-coupled matrix element **
//***************************************************************************
// NOTE: no assumptions are made for the index ordering - routine determines
//		 phases & possible 6j factors from ordering indices, and calls
//		 GetME_Ordered
// TODO: improve this!!!
//***************************************************************************                    2J????????????????
double ME3J_GetME(HF_Matrix_Elements &HF_ME, int nlj1, int nlj2, int nlj3, int nnlj1, int nnlj2, int nnlj3, int j12, int jj12, int twoJ, int t12, int tt12, int twoT)
{
  long long	idx, iidx;
  int		tmp, nlja, nljb, nljc, nnlja, nnljb, nnljc;
  int		twoja, twojb, twojja, twojjb;
  int		j2, jj2, j2Min, j2Max, jj2Min, jj2Max, t2, tt2, t2Min, t2Max, tt2Min, tt2Max;
  int		nlj1x, nlj2x, nlj3x, nnlj1x, nnlj2x, nnlj3x, j12x, jj12x, t12x, tt12x;

  //*** do the internal math with double accuracy, for safety reasons
  double	me, ovl1, ovl2;

  //*** E3Max cutoff must be caught here, else there might be illegal memory accesses in ME3JH
  if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] + SPB.e_nlj3[nlj3] > SPB.E3Max || SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] + SPB.e_nlj3[nnlj3] > SPB.E3Max )
    { return 0.0; }

  //*** sort indices in descending order
  nlja  = nlj1;  nljb  = nlj2;  nljc  = nlj3;
  nnlja = nnlj1; nnljb = nnlj2; nnljc = nnlj3;

  if( nlja < nljb ){ tmp = nlja; nlja = nljb; nljb = tmp; }
  if( nlja < nljc ){ tmp = nlja; nlja = nljc; nljc = tmp; }
  if( nljb < nljc ){ tmp = nljb; nljb = nljc; nljc = tmp; }
  if( nnlja < nnljb ){ tmp = nnlja; nnlja = nnljb; nnljb = tmp; }
  if( nnlja < nnljc ){ tmp = nnlja; nnlja = nnljc; nnljc = tmp; }
  if( nnljb < nnljc ){ tmp = nnljb; nnljb = nnljc; nnljc = tmp; }

  //*** lower triangle is stored - swap indices if necessary
  idx  = (nlja * SPB.nljDim3 + nljb) * SPB.nljDim3 + nljc;
  iidx = (nnlja * SPB.nljDim3 + nnljb) * SPB.nljDim3 + nnljc;

  if( iidx > idx ){
    tmp = nlja; nlja = nnlja; nnlja = tmp;
    tmp = nljb; nljb = nnljb; nnljb = tmp;
    tmp = nljc; nljc = nnljc; nnljc = tmp;

    nnlj1x = nlj1;
    nnlj2x = nlj2;
    nnlj3x = nlj3;

    nlj1x = nnlj1;
    nlj2x = nnlj2;
    nlj3x = nnlj3;

    j12x = jj12;
    jj12x = j12;

    t12x = tt12;
    tt12x = t12;
  }
  else{
    nlj1x = nlj1;
    nlj2x = nlj2;
    nlj3x = nlj3;

    nnlj1x = nnlj1;
    nnlj2x = nnlj2;
    nnlj3x = nnlj3;

    j12x = j12;
    jj12x = jj12;
    
    t12x = t12;
    tt12x = tt12;
  }

  twoja = SPB.twoj_nlj3[nlja];
  twojb = SPB.twoj_nlj3[nljb];

  twojja = SPB.twoj_nlj3[nnlja];
  twojjb = SPB.twoj_nlj3[nnljb];

  if( nljc == nlj3x ){
    j2Min = j12x;
    j2Max = j12x;
    
    t2Min = t12x;
    t2Max = t12x;
  }
  else{
    j2Min = std::abs(twoja - twojb);
    j2Max = twoja + twojb;
    
    t2Min = 0;
    t2Max = 2;
  }

  if( nnljc == nnlj3x ){
    jj2Min = jj12x;
    jj2Max = jj12x;
    
    tt2Min = tt12x;
    tt2Max = tt12x;
  }
  else{
    jj2Min = std::abs(twojja - twojjb);
    jj2Max = twojja + twojjb;
    
    tt2Min = 0;
    tt2Max = 2;
  }

  me = 0.0;
  for(j2 = j2Min; j2 <= j2Max; j2 += 2){
    for(t2 = t2Min; t2 <= t2Max; t2 += 2){
      //if( t2 < std::abs(twoT - 1) || t2 > twoT + 1 ){ continue; }
      ovl1 = ME3J_Overlap(HF_ME, nlj1x, nlj2x, nlj3x, j12x, t12x, nlja, nljb, nljc, j2, t2, twoJ, twoT);
      
      if( std::fabs(ovl1) > 1.0e-12 ){
	for(jj2 = jj2Min; jj2 <= jj2Max; jj2 += 2){
	  for(tt2 = tt2Min; tt2 <= tt2Max; tt2 += 2){
	    //if( tt2 < std::abs(twoT - 1) || tt2 > twoT + 1 ){ continue; }
	    ovl2 = ME3J_Overlap(HF_ME, nnlj1x, nnlj2x, nnlj3x, jj12x, tt12x, nnlja, nnljb, nnljc, jj2, tt2, twoJ, twoT);
	    
	    if( std::fabs(ovl1*ovl2) > 1.0e-12 ){
	      me += ovl1 * ovl2 * ME3J_GetME_Ordered(HF_ME, nlja, nljb, nljc, nnlja, nnljb, nnljc, j2, jj2, twoJ, t2, tt2, twoT);
	    }
	  }
	}
      }
    }
  }

  return me;
}


//********************************************************************
// Overlap
//********************************************************************
// returns the overlap between states with different index ordering
//********************************************************************
// NOTE: we assume nlja >= nljb >= nljc, but no particular order for
//       nlj1, nlj2, nlj3
//********************************************************************               2J????????????????????
double ME3J_Overlap(HF_Matrix_Elements &ME, int nlj1, int nlj2, int nlj3, int j12, int t12, int nlja, int nljb, int nljc, int jab, int tab, int twoJ, int twoT)
{
  int		twoj1, twoj2, twoj3;
  double	phase;

  twoj1 = SPB.twoj_nlj3[nlj1];
  twoj2 = SPB.twoj_nlj3[nlj2];
  twoj3 = SPB.twoj_nlj3[nlj3];

  phase = 1.0;

  if( nljc == nlj3 ){
    if( tab != t12 || jab != j12 ){ return 0.0; }
    if( nlja == nlj1 ){ return 1.0; }    //*** (abc) = (123)
    else{                                //*** (abc) = (213)
      if( ((j12 - (twoj1 + twoj2) + t12)/2)%2 ){ phase = -1.0; }
      return phase;
    }
  }

  if( tab < std::abs(twoT - 1) || tab > twoT + 1 ){ return 0.0; }
  
  if( nljc == nlj1 ){
    if( nlja == nlj3 ){ phase = -1.0; }  //*** (abc) = (321)
    else{                                //*** (abc) = (231)
      if( (((twoj2 + twoj3) - jab - tab + 2)/2)%2 ){ phase = -1.0; }
    }
    /*std::cout << "( " << twoj1 << "," << twoj2 << "," << j12 << "," << twoj3 << "," << twoJ << "," << jab << " ) = " << std::endl;
    std::cout << JCOUP.get_SixJ(twoj1, twoj2, j12, twoj3, twoJ, jab) << std::endl;
    std::cout << "( " << 1 << "," << 1 << "," << t12 << "," << 1 << "," << twoT << "," << tab << " ) = " << std::endl;
    std::cout << JCOUP.get_SixJ(1, 1, t12, 1, twoT, tab) << std::endl;*/

    return phase * std::sqrt((j12 + 1.0)*(jab + 1.0)*(t12 + 1.0)*(tab + 1.0)) *
      JCOUP.get_SixJ(twoj1, twoj2, j12, twoj3, twoJ, jab) * JCOUP.get_SixJ(1, 1, t12, 1, twoT, tab);
    //CGC6(twoj1, twoj2, j12, twoj3, twoJ, jab) * CGC6(1, 1, t12, 1, twoT, tab);
  }

  if( nljc == nlj2 ){
    if( nlja == nlj1 ){	                 //*** (abc) = (132)
      if( ((-(twoj2 - twoj3) + j12 - jab + t12 - tab + 2)/2)%2 ){ phase = -1.0; }
    }
    else{	                         //*** (abc) = (312)
      if( ((j12 - (twoj1 + twoj2) + t12 - 2)/2)%2 ){ phase = -1.0; }
    }
    /*std::cout << "( " << twoj2 << "," << twoj1 << "," << j12 << "," << twoj3 << "," << twoJ << "," << jab << " ) = " << std::endl;
    std::cout << JCOUP.get_SixJ(twoj2, twoj1, j12, twoj3, twoJ, jab) << std::endl;
    std::cout << "( " << 1 << "," << 1 << "," << t12 << "," << 1 << "," << twoT << "," << tab << " ) = " << std::endl;
    std::cout << JCOUP.get_SixJ(1, 1, t12, 1, twoT, tab) << std::endl;*/

    return  phase * std::sqrt((j12 + 1.0)*(jab + 1.0)*(t12 + 1.0)*(tab + 1.0)) *
      JCOUP.get_SixJ(twoj2, twoj1, j12, twoj3, twoJ, jab) * JCOUP.get_SixJ(1, 1, t12, 1, twoT, tab);
    //CGC6(twoj2, twoj1, j12, twoj3, twoJ, jab) * CGC6(1, 1, t12, 1, twoT, tab);
  }

  return 0.0;
}


//***************************************************************************
// GetME_Ordered
//***************************************************************************
// get specific [(ja,jb)Jab,jc]J[(ta,tb)tab,tc]T-coupled matrix element **
//***************************************************************************
// NOTE: this requires indices to be ordered - calls to matrix elements should
//       be made through GetME, which will ensure proper ordering
//
// TODO: further optimizations?
//***************************************************************************                     2J?????????????????????
float ME3J_GetME_Ordered(HF_Matrix_Elements &HF_ME, int nlj1, int nlj2, int nlj3, int nnlj1, int nnlj2, int nnlj3, int jab0, int jjab0, int twoJ0, int tab0, int ttab0, int twoT0)
{
  int   jab2, jabMin2, jabMax2, jjab2, jjabMin2, jjabMax2, twoJMin, twoJMax;
  long long  i;

  //*** initial index
  i = HF_ME.ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2][nnlj3];

  //*** parity or triangular condition violated
  if( i < 0 ){ return(0.0); }

  jabMax2 = (SPB.twoj_nlj3[nlj1] + SPB.twoj_nlj3[nlj2]);
  jabMin2 = std::abs(SPB.twoj_nlj3[nlj1] - SPB.twoj_nlj3[nlj2]);

  jjabMax2 = (SPB.twoj_nlj3[nnlj1] + SPB.twoj_nlj3[nnlj2]);
  jjabMin2 = std::abs(SPB.twoj_nlj3[nnlj1] - SPB.twoj_nlj3[nnlj2]);

  //*** possible isospin combinations:
  //***
  //***         tab             ttab            twoT            twoT/2
  //***         0               0               1               0
  //***         0               1               1               0
  //***         1               0               1               0
  //***         1               1               1               0
  //***         1               1               3               1

  for(jab2 = jabMin2; jab2 <= jabMax2; jab2 += 2){
    for(jjab2 = jjabMin2; jjab2 <= jjabMax2; jjab2 += 2){
      //*** diagonal loops
      twoJMin = std::max( std::abs(jab2 - SPB.twoj_nlj3[nlj3]), std::abs(jjab2 - SPB.twoj_nlj3[nnlj3]) );
      twoJMax = std::min( jab2 + SPB.twoj_nlj3[nlj3], jjab2 + SPB.twoj_nlj3[nnlj3] );

      if( twoJMin > twoJMax ){ continue; }

      if( jab2 == jab0 && jjab2 == jjab0 ){
    	//i += ((twoJ0 - twoJMin)/2)*5 + tab0*2 + ttab0 + (twoT0/2);
	i += ((twoJ0 - twoJMin)/2)*5 + tab0 + (ttab0/2) + (twoT0/2);
    	return  HF_ME.V3[i];
      }
      i += ((twoJMax - twoJMin)/2 + 1)*5;
    }
  }

  return 0.0;
}
