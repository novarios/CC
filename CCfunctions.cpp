#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"

//   Function to return Hash index for 2 indices
int Hash2(const int &p, const int &q, const int &size){ return size*p + q; }
//   Function to return Hash index for 3 indices
int Hash3(const int &p, const int &q, const int &r, const int &size){ return size*size*p + q*size + r; }

void plus(State &S, const State &S1, const State &S2){
  S.t = S1.t + S2.t;
  S.m = S1.m + S2.m;
  S.nx = S1.nx + S2.nx;
  S.ny = S1.ny + S2.ny;
  S.nz = S1.nz + S2.nz;
  S.ml = S1.ml + S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

void minus(State &S, const State &S1, const State &S2){
  S.t = S1.t - S2.t;
  S.m = S1.m - S2.m;
  S.nx = S1.nx - S2.nx;
  S.ny = S1.ny - S2.ny;
  S.nz = S1.nz - S2.nz;
  S.ml = S1.ml - S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

bool equal(const State &S1, const State &S2){
  return (S1.t == S2.t &&
	  S1.m == S2.m &&
	  S1.nx == S2.nx &&
	  S1.ny == S2.ny &&
	  S1.nz == S2.nz &&
	  S1.ml == S2.ml &&
	  S1.par == S2.par &&
	  S1.j == S2.j);
}

double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  int nhh, ind;
  State tb;
  for(int i = 0; i < Space.indhol; ++i){
    if(Parameters.basis == "finite_J"){ energy += (Space.qnums[i].j + 1) * Space.qnums[i].energy; }
    else{ energy += Space.qnums[i].energy; }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    if(nhh == 0){ continue; }
    for(int hh = 0; hh < nhh; ++hh){
      ind = hh*nhh + hh;
      if(Parameters.basis == "finite_J"){ energy -= 0.5 * (Chan.qnums1[chan].j + 1) * Ints.D_ME1.V2[chan][ind]; }
      else{ energy -= 0.5 * Ints.D_ME1.V2[chan][ind]; }
    }
  }
  return energy;
}

double Amplitudes::get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  int nhh, npp;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }    
    for(int hh = 0; hh < Chan.nhh[chan]; ++hh){
      for(int pp = 0; pp < Chan.npp[chan]; ++pp){
	if(Parameters.basis == "finite_J"){
	  energy += 0.25 * (Chan.qnums1[chan].j + 1) * D1.T1[chan][hh * Chan.npp[chan] + pp] * Ints.D_ME1.V4[chan][pp * Chan.nhh[chan] + hh];
	}
	else{
	  energy += 0.25 * D1.T1[chan][hh * Chan.npp[chan] + pp] * Ints.D_ME1.V4[chan][pp * Chan.nhh[chan] + hh];
	}
      }
    }
  }
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp1[Chan.ind0]; ++hp2){
	energy += 0.5 * S1.T1[hp1] * S1.T1[hp2] * Ints.D_ME1.V9[Chan.ind0][hp1 * Chan.nhp1[Chan.ind0] + hp2];
      }
    }
  }
  return energy;
}

void Amplitudes::copy_Amplitudes(const Input_Parameters &Parameters, const Channels &Chan, const Amplitudes &Amps)
{
  if(Parameters.approx == "doubles"){
    D1.copy_Doubles_1(Chan, Amps.D1);
  }
  else if(Parameters.approx == "singles"){
    D1.copy_Doubles_1(Chan, Amps.D1);
    S1.copy_Singles_1(Chan, Amps.S1);
  }
}

Amplitudes::Amplitudes(const Input_Parameters &Parameters, const Channels &Chan, const Amplitudes &Amps)
{
  if(Parameters.approx == "doubles"){
    D1 = Doubles_1(Chan, Amps.D1);
  }
  else if(Parameters.approx == "singles"){
    D1 = Doubles_1(Chan, Amps.D1);
    S1 = Singles_1(Chan, Amps.S1);
  }
}

Amplitudes::Amplitudes(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1 = Doubles_1(Parameters, Space, Chan);
  }
  else if(Parameters.approx == "singles"){
    D1 = Doubles_1(Parameters, Space, Chan);
    S1 = Singles_1(Parameters, Space, Chan);
  }
}

void Amplitudes::delete_struct(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.delete_struct(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.delete_struct(Chan);
    S1.delete_struct(Chan);
  }
}

void Amplitudes::zero(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.zero(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.zero(Chan);
    S1.zero(Chan);
  }
}

void Amplitudes::zero1(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D1.zero1(Chan);
  }
  else if(Parameters.approx == "singles"){
    D1.zero1(Chan);
    S1.zero(Chan);
  }
}

void Doubles_1::copy_Doubles_1(const Channels &Chan, const Doubles_1 &D)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	for(int k = 0; k < 16; ++k){ Tmap[i][16*j + k] = D.Tmap[i][16*j + k]; }
	for(int k = 0; k < 4; ++k){
	  TJnum[i][4*j + k] = D.TJnum[i][4*j + k];
	  for(int l = 0; l < 2*TJnum[i][4*j + k]; ++l){
	    TJmap[i][j][k][l] = D.TJmap[i][j][k][l];
	  }
	}
	Evec[i][j] = D.Evec[i][j];
	T1[i][j] = D.T1[i][j];
      }
    }
    length = nhh * nhh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S1[i][j] = D.S1[i][j];
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = D.Q11[i][j];
	Qmap1[i][2*j] = D.Qmap1[i][2*j];
	Qmap1[i][2*j + 1] = D.Qmap1[i][2*j + 1];
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = D.Q21[i][j];
	Qmap2[i][2*j] = D.Qmap2[i][2*j];
	Qmap2[i][2*j + 1] = D.Qmap2[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T6[i][j] = D.T6[i][j];
	T7[i][j] = D.T7[i][j];
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T8[i][j] = D.T8[i][j];
	T9[i][j] = D.T9[i][j];
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = D.Q12[i][j];
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      for(int j = 0; j < length; ++j){
	Q22[i][j] = D.Q22[i][j];
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S2[i][j] = D.S2[i][j];
	S3[i][j] = D.S3[i][j];
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S4[i][j] = D.S4[i][j];
	S5[i][j] = D.S5[i][j];
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = D.T2[i][j];
	T3[i][j] = D.T3[i][j];
	T4[i][j] = D.T4[i][j];
	T5[i][j] = D.T5[i][j];
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S6[i][j] = D.S6[i][j];
	S7[i][j] = D.S7[i][j];
      }
    }
  }
}

Doubles_1::Doubles_1(const Channels &Chan, const Doubles_1 &D)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  Tmap = new int*[Chan.size1];
  TJnum = new int*[Chan.size1];
  TJmap = new int***[Chan.size1];
  Evec = new double*[Chan.size1];
  T1 = new double*[Chan.size1];
  T2 = new double*[Chan.size2];
  T3 = new double*[Chan.size2];
  T4 = new double*[Chan.size2]; 
  T5 = new double*[Chan.size2];
  T6 = new double*[Chan.size3];
  T7 = new double*[Chan.size3];
  T8 = new double*[Chan.size3];
  T9 = new double*[Chan.size3];
  S1 = new double*[Chan.size1];
  S2 = new double*[Chan.size3];
  S3 = new double*[Chan.size3];
  S4 = new double*[Chan.size3];
  S5 = new double*[Chan.size3];
  S6 = new double*[Chan.size2];
  S7 = new double*[Chan.size2];

  Q11 = new double*[Chan.size1];
  Q21 = new double*[Chan.size1];
  Qmap1 = new int*[Chan.size1];
  Qmap2 = new int*[Chan.size1];
  Q12 = new double*[Chan.size3];
  Q22 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      Tmap[i] = new int[16 * length];
      TJnum[i] = new int[4 * length];
      TJmap[i] = new int**[length];
      Evec[i] = new double[length];
      T1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	for(int k = 0; k < 16; ++k){ Tmap[i][16*j + k] = D.Tmap[i][16*j + k]; }
	TJmap[i][j] = new int*[4];
	for(int k = 0; k < 4; ++k){
	  TJnum[i][4*j + k] = D.TJnum[i][4*j + k];
	  TJmap[i][j][k] = new int[2*TJnum[i][4*j + k]];
	  for(int l = 0; l < 2*TJnum[i][4*j + k]; ++l){
	    TJmap[i][j][k][l] = D.TJmap[i][j][k][l];
	  }
	}
	Evec[i][j] = D.Evec[i][j];
	T1[i][j] = D.T1[i][j];
      }
    }
    length = nhh * nhh;
    if(length != 0){
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S1[i][j] = D.S1[i][j];
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q11[i] = new double[length];
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = D.Q11[i][j];
	Qmap1[i][2*j] = D.Qmap1[i][2*j];
	Qmap1[i][2*j + 1] = D.Qmap1[i][2*j + 1];
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q21[i] = new double[length];
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = D.Q21[i][j];
	Qmap2[i][2*j] = D.Qmap2[i][2*j];
	Qmap2[i][2*j + 1] = D.Qmap2[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      T6[i] = new double[length];
      T7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T6[i][j] = D.T6[i][j];
	T7[i][j] = D.T7[i][j];
      }
    }
    length = nhhp * np;
    if(length != 0){
      T8[i] = new double[length];
      T9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T8[i][j] = D.T8[i][j];
	T9[i][j] = D.T9[i][j];
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q12[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = D.Q12[i][j];
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q22[i][j] = D.Q22[i][j];
      }
    }
    length = nh * nh;
    if(length != 0){
      S2[i] = new double[length];
      S3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S2[i][j] = D.S2[i][j];
	S3[i][j] = D.S3[i][j];
      }
    }
    length = np * np;
    if(length != 0){
      S4[i] = new double[length];
      S5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S4[i][j] = D.S4[i][j];
	S5[i][j] = D.S5[i][j];
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      T4[i] = new double[length];
      T5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = D.T2[i][j];
	T3[i][j] = D.T3[i][j];
	T4[i][j] = D.T4[i][j];
	T5[i][j] = D.T5[i][j];
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      S6[i] = new double[length];
      S7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S6[i][j] = D.S6[i][j];
	S7[i][j] = D.S7[i][j];
      }
    }
  }
}

Doubles_1::Doubles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  int length, length0, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  State tb;
  int ind, ind1, ind2, key1, key2, jmin, count;
  int hind1, hind2, hind3, pind1, pind2, pind3;
  int hh1, pp1, hp1;
  
  Tmap = new int*[Chan.size1];
  TJnum = new int*[Chan.size1];
  TJmap = new int***[Chan.size1];
  Evec = new double*[Chan.size1];
  T1 = new double*[Chan.size1];
  T2 = new double*[Chan.size2];
  T3 = new double*[Chan.size2];
  T4 = new double*[Chan.size2]; 
  T5 = new double*[Chan.size2];
  T6 = new double*[Chan.size3];
  T7 = new double*[Chan.size3];
  T8 = new double*[Chan.size3];
  T9 = new double*[Chan.size3];
  S1 = new double*[Chan.size1];
  S2 = new double*[Chan.size3];
  S3 = new double*[Chan.size3];
  S4 = new double*[Chan.size3];
  S5 = new double*[Chan.size3];
  S6 = new double*[Chan.size2];
  S7 = new double*[Chan.size2];

  Q11 = new double*[Chan.size1];
  Q21 = new double*[Chan.size1];
  Qmap1 = new int*[Chan.size1];
  Qmap2 = new int*[Chan.size1];
  Q12 = new double*[Chan.size3];
  Q22 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      Tmap[i] = new int[16 * length];
      TJnum[i] = new int[4 * length];
      TJmap[i] = new int**[length];
      Evec[i] = new double[length];
      T1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	for(int k = 0; k < 16; ++k){ Tmap[i][16*j + k] = -1; }
	TJmap[i][j] = new int*[4];
	for(int k = 0; k < 4; ++k){ TJnum[i][4*j + k] = 0; }
	Evec[i][j] = 0.0;
	T1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S1[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q11[i] = new double[length];
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
	Qmap1[i][2*j] = -1;
	Qmap1[i][2*j + 1] = -1;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q21[i] = new double[length];
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
	Qmap2[i][2*j] = -1;
	Qmap2[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      T6[i] = new double[length];
      T7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T6[i][j] = 0.0;
	T7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      T8[i] = new double[length];
      T9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T8[i][j] = 0.0;
	T9[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q12[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q22[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      S2[i] = new double[length];
      S3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      S4[i] = new double[length];
      S5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S4[i][j] = 0.0;
	S5[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      T4[i] = new double[length];
      T5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
	T4[i][j] = 0.0;
	T5[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      S6[i] = new double[length];
      S7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	S6[i][j] = 0.0;
	S7[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    length = Chan.nhh[i] * Chan.npp[i];
    length0 = Chan.npp[i];
    for(int hhpp = 0; hhpp < length; ++hhpp){
      pp1 = int(hhpp%length0);
      hh1 = int((hhpp - pp1)/length0);
      hind1 = Chan.hhvec[i][2*hh1];
      hind2 = Chan.hhvec[i][2*hh1 + 1];
      pind1 = Chan.ppvec[i][2*pp1];
      pind2 = Chan.ppvec[i][2*pp1 + 1];
      if(Parameters.basis == "finite_J"){
	//T2 = ((ia)(jb)')
	jmin = abs(Space.qnums[hind1].j - Space.qnums[pind1].j);
	if(abs(Space.qnums[pind2].j - Space.qnums[hind2].j) > jmin){ jmin = abs(Space.qnums[pind2].j - Space.qnums[hind2].j); }
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	if(Space.qnums[pind2].j + Space.qnums[hind2].j < tb.j){ tb.j = Space.qnums[pind2].j + Space.qnums[hind2].j; }
	TJnum[i][4*hhpp] = int(0.5*(tb.j - jmin) + 1);
	TJmap[i][hhpp][0] = new int[2 * int(0.5*(tb.j - jmin) + 1)];
	count = 0;
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind2, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  TJmap[i][hhpp][0][2*count] = ind;
	  TJmap[i][hhpp][0][2*count + 1] = ind1;
	  tb.j -= 2;
	  ++count;
	}
	//T3 = ((jb)(ia)')
	jmin = abs(Space.qnums[hind2].j - Space.qnums[pind2].j);
	if(abs(Space.qnums[pind1].j - Space.qnums[hind1].j) > jmin){ jmin = abs(Space.qnums[pind1].j - Space.qnums[hind1].j); }
	minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
	if(Space.qnums[pind1].j + Space.qnums[hind1].j < tb.j){ tb.j = Space.qnums[pind1].j + Space.qnums[hind1].j; }
	TJnum[i][4*hhpp + 1] = int(0.5*(tb.j - jmin) + 1);
	TJmap[i][hhpp][1] = new int[2 * int(0.5*(tb.j - jmin) + 1)];
	count = 0;
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  TJmap[i][hhpp][1][2*count] = ind;
	  TJmap[i][hhpp][1][2*count + 1] = ind1;
	  tb.j -= 2;
	  ++count;
	}
	//T4 = ((ib)(ja)')
	jmin = abs(Space.qnums[hind1].j - Space.qnums[pind2].j);
	if(abs(Space.qnums[pind1].j - Space.qnums[hind2].j) > jmin){ jmin = abs(Space.qnums[pind1].j - Space.qnums[hind2].j); }
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
	if(Space.qnums[pind1].j + Space.qnums[hind2].j < tb.j){ tb.j = Space.qnums[pind1].j + Space.qnums[hind2].j; }
	TJnum[i][4*hhpp + 2] = int(0.5*(tb.j - jmin) + 1);
	TJmap[i][hhpp][2] = new int[2 * int(0.5*(tb.j - jmin) + 1)];
	count = 0;
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  TJmap[i][hhpp][2][2*count] = ind;
	  TJmap[i][hhpp][2][2*count + 1] = ind1;
	  tb.j -= 2;
	  ++count;
	}
	//T5 = ((ja)(ib)')
	jmin = abs(Space.qnums[hind2].j - Space.qnums[pind1].j);
	if(abs(Space.qnums[pind2].j - Space.qnums[hind1].j) > jmin){ jmin = abs(Space.qnums[pind2].j - Space.qnums[hind1].j); }
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
	if(Space.qnums[pind2].j + Space.qnums[hind1].j < tb.j){ tb.j = Space.qnums[pind2].j + Space.qnums[hind1].j; }
	TJnum[i][4*hhpp + 3] = int(0.5*(tb.j - jmin) + 1);
	TJmap[i][hhpp][3] = new int[2 * int(0.5*(tb.j - jmin) + 1)];
	count = 0;
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  TJmap[i][hhpp][3][2*count] = ind;
	  TJmap[i][hhpp][3][2*count + 1] = ind1;
	  tb.j -= 2;
	  ++count;
	}
	//T6 = ((jab)(i))
	ind = Chan.indvec[hind1];
	key1 = Chan.hpp_map[ind][int(0.5 * Chan.qnums1[i].j * std::pow(Space.indtot, 3)) + Hash3(hind2, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	Tmap[i][16 * hhpp + 8] = ind;
	Tmap[i][16 * hhpp + 9] = ind1;
	//T7 = ((iab)(j))
	ind = Chan.indvec[hind2];
	key1 = Chan.hpp_map[ind][int(0.5 * Chan.qnums1[i].j * std::pow(Space.indtot, 3)) + Hash3(hind1, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	Tmap[i][16 * hhpp + 10] = ind;
	Tmap[i][16 * hhpp + 11] = ind1;
	//T8 = ((ijb)(a))
	ind = Chan.indvec[pind1];
	key1 = Chan.hhp_map[ind][int(0.5 * Chan.qnums1[i].j * std::pow(Space.indtot, 3)) + Hash3(hind1, hind2, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	Tmap[i][16 * hhpp + 12] = ind;
	Tmap[i][16 * hhpp + 13] = ind1;
	//T9 = ((ija)(b))
	ind = Chan.indvec[pind2];
	key1 = Chan.hhp_map[ind][int(0.5 * Chan.qnums1[i].j * std::pow(Space.indtot, 3)) + Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	Tmap[i][16 * hhpp + 14] = ind;
	Tmap[i][16 * hhpp + 15] = ind1;
      }
      else{
	TJmap[i][hhpp][0] = new int[0];
	TJmap[i][hhpp][1] = new int[0];
	TJmap[i][hhpp][2] = new int[0];
	TJmap[i][hhpp][3] = new int[0];
	//T2 = ((ia)(jb)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Tmap[i][16 * hhpp] = ind;
	Tmap[i][16 * hhpp + 1] = ind1;
	//T3 = ((jb)(ia)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Tmap[i][16 * hhpp + 2] = ind;
	Tmap[i][16 * hhpp + 3] = ind1;
	//T4 = ((ib)(ja)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Tmap[i][16 * hhpp + 4] = ind;
	Tmap[i][16 * hhpp + 5] = ind1;
	//T5 = ((ja)(ib)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Tmap[i][16 * hhpp + 6] = ind;
	Tmap[i][16 * hhpp + 7] = ind1;
	//T6 = ((jab)(i))
	ind = Chan.indvec[hind1];
	key1 = Chan.hpp_map[ind][Hash3(hind2, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind1];
	ind1 = key1 * Chan.nh[ind] + key2;
	Tmap[i][16 * hhpp + 8] = ind;
	Tmap[i][16 * hhpp + 9] = ind1;
	//T7 = ((iab)(j))
	ind = Chan.indvec[hind2];
	key1 = Chan.hpp_map[ind][Hash3(hind1, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[ind][hind2];
	ind1 = key1 * Chan.nh[ind] + key2;
	Tmap[i][16 * hhpp + 10] = ind;
	Tmap[i][16 * hhpp + 11] = ind1;
	//T8 = ((ijb)(a))
	ind = Chan.indvec[pind1];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind2, Space.indtot)];
	key2 = Chan.p_map[ind][pind1];
	ind1 = key1 * Chan.np[ind] + key2;
	Tmap[i][16 * hhpp + 12] = ind;
	Tmap[i][16 * hhpp + 13] = ind1;
	//T9 = ((ija)(b))
	ind = Chan.indvec[pind2];
	key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.p_map[ind][pind2];
	ind1 = key1 * Chan.np[ind] + key2;
	Tmap[i][16 * hhpp + 14] = ind;
	Tmap[i][16 * hhpp + 15] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    length = Chan.nhp[i] * Chan.npp[i];
    length0 = Chan.npp[i];
    for(int hppp = 0; hppp < length; ++hppp){
      pp1 = int(hppp%length0);
      hp1 = int((hppp - pp1)/length0);
      hind1 = Chan.hpvec[i][2*hp1];
      pind1 = Chan.hpvec[i][2*hp1 + 1];
      pind2 = Chan.ppvec[i][2*pp1];
      pind3 = Chan.ppvec[i][2*pp1 + 1];
      ind = Chan.indvec[pind1];
      key1 = Chan.hpp_map[ind][Hash3(hind1, pind2, pind3, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Qmap1[i][2 * hppp] = ind;
      Qmap1[i][2 * hppp + 1] = ind1;
    }
    length = Chan.nhh[i] * Chan.nhp[i];
    length0 = Chan.nhp[i];
    for(int hhhp = 0; hhhp < length; ++hhhp){
      hp1 = int(hhhp%length0);
      hh1 = int((hhhp - hp1)/length0);
      hind1 = Chan.hhvec[i][2*hh1];
      hind2 = Chan.hhvec[i][2*hh1 + 1];
      hind3 = Chan.hpvec[i][2*hp1];
      pind1 = Chan.hpvec[i][2*hp1 + 1];
      ind = Chan.indvec[hind3];
      key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
      key2 = Chan.h_map[ind][hind3];
      ind1 = key1 * Chan.nh[ind] + key2;
      Qmap2[i][2 * hhhp] = ind;
      Qmap2[i][2 * hhhp + 1] = ind1;
    }
  }
}

void Doubles_1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(nhh * npp != 0){
      for(int j = 0; j < nhh * npp; ++j){ 
	for(int k = 0; k < 4; ++k){ delete[] TJmap[i][j][k]; }
	delete[] TJmap[i][j];
      }
      delete[] Tmap[i];
      delete[] TJnum[i];
      delete[] TJmap[i];
      delete[] Evec[i];
      delete[] T1[i];
    }
    if(nhh != 0){
      delete[] S1[i];
    }
    if(nhp * npp != 0){
      delete[] Q11[i];
      delete[] Qmap1[i];
    }
    if(nhh * nhp != 0){
      delete[] Q21[i];
      delete[] Qmap2[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nhpp * nh != 0){
      delete[] T6[i];
      delete[] T7[i];
    }
    if(nhhp * np != 0){
      delete[] T8[i];
      delete[] T9[i];
    }
    if(nhpp * np != 0){
      delete[] Q12[i];
    }
    if(nhhp * nh != 0){
      delete[] Q22[i];
    }
    if(nh != 0){
      delete[] S2[i];
      delete[] S3[i];
    }
    if(np != 0){
      delete[] S4[i];
      delete[] S5[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp1 * nhp2 != 0){
      delete[] T2[i];
      delete[] T3[i];
      delete[] T4[i];
      delete[] T5[i];
    }
    if(nhp2 != 0){
      delete[] S6[i];
      delete[] S7[i];
    }
  }
  delete[] Tmap;
  delete[] TJnum;
  delete[] TJmap;
  delete[] Evec;
  delete[] T1;
  delete[] T2;
  delete[] T3;
  delete[] T4;
  delete[] T5;
  delete[] T6;
  delete[] T7;
  delete[] T8;
  delete[] T9;
  delete[] S1;
  delete[] S2;
  delete[] S3;
  delete[] S4;
  delete[] S5;
  delete[] S6;
  delete[] S7;

  delete[] Q11;
  delete[] Q21;
  delete[] Qmap1;
  delete[] Qmap2;
  delete[] Q12;
  delete[] Q22;
}

void Doubles_1::zero(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S1[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T6[i][j] = 0.0;
	T7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T8[i][j] = 0.0;
	T9[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      for(int j = 0; j < length; ++j){
	Q22[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S4[i][j] = 0.0;
	S5[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
	T4[i][j] = 0.0;
	T5[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S6[i][j] = 0.0;
	S7[i][j] = 0.0;
      }
    }
  }
}

void Doubles_1::zero1(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * nhh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S1[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T6[i][j] = 0.0;
	T7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T8[i][j] = 0.0;
	T9[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(nhhp * nh != 0){
      for(int j = 0; j < length; ++j){
	Q22[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S2[i][j] = 0.0;
	S3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S4[i][j] = 0.0;
	S5[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
	T4[i][j] = 0.0;
	T5[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	S6[i][j] = 0.0;
	S7[i][j] = 0.0;
      }
    }
  }
}

void Doubles_1::set_T(int i, int j, double T)
{
  T1[i][j] = T;
  T2[Tmap[i][16*j]][Tmap[i][16*j + 1]] = T;
  T3[Tmap[i][16*j + 2]][Tmap[i][16*j + 3]] = T;
  T4[Tmap[i][16*j + 4]][Tmap[i][16*j + 5]] = T;
  T5[Tmap[i][16*j + 6]][Tmap[i][16*j + 7]] = T;
  T6[Tmap[i][16*j + 8]][Tmap[i][16*j + 9]] = T;
  T7[Tmap[i][16*j + 10]][Tmap[i][16*j + 11]] = T;
  T8[Tmap[i][16*j + 12]][Tmap[i][16*j + 13]] = T;
  T9[Tmap[i][16*j + 14]][Tmap[i][16*j + 15]] = T;
}

double Doubles_1::get_T(int i, int j) const
{
  double tempt;
  tempt = T1[i][j];
  tempt += T2[Tmap[i][16*j]][Tmap[i][16*j + 1]];    
  tempt += T3[Tmap[i][16*j + 2]][Tmap[i][16*j + 3]];
  tempt += T4[Tmap[i][16*j + 4]][Tmap[i][16*j + 5]];
  tempt += T5[Tmap[i][16*j + 6]][Tmap[i][16*j + 7]];
  tempt += T6[Tmap[i][16*j + 8]][Tmap[i][16*j + 9]];
  tempt += T7[Tmap[i][16*j + 10]][Tmap[i][16*j + 11]];
  tempt += T8[Tmap[i][16*j + 12]][Tmap[i][16*j + 13]];
  tempt += T9[Tmap[i][16*j + 14]][Tmap[i][16*j + 15]];
  return tempt;
}

void Doubles_1::set_TJ(const Model_Space &Space, const Channels &Chan, int &chan, int &hhpp, int &i, int &j, int &a, int &b, double T)
{
  int chan1, ind1;
  double J, ij, jj, aj, bj, J1;
  J = 0.5 * Chan.qnums1[chan].j;
  ij = 0.5 * Space.qnums[i].j;
  jj = 0.5 * Space.qnums[j].j;
  aj = 0.5 * Space.qnums[a].j;
  bj = 0.5 * Space.qnums[b].j;
  T1[chan][hhpp] = T;
  for(int k = 0; k < TJnum[chan][4*hhpp]; ++k){
    chan1 = TJmap[chan][hhpp][0][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][0][2*k + 1];
    T2[chan1][ind1] += std::pow(-1.0, int(aj + bj - J)) * (2.0 * J + 1) * CGC6(bj, aj, J, ij, jj, J1) * T;
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 1]; ++k){
    chan1 = TJmap[chan][hhpp][0][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][0][2*k + 1];
    T3[chan1][ind1] += std::pow(-1.0, int(ij + jj - J)) * (2.0 * J + 1) * CGC6(aj, bj, J, jj, ij, J1) * T;
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 2]; ++k){
    chan1 = TJmap[chan][hhpp][0][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][0][2*k + 1];
    T4[chan1][ind1] += -1.0 * (2.0 * J + 1) * CGC6(aj, bj, J, ij, jj, J1) * T;
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 3]; ++k){
    chan1 = TJmap[chan][hhpp][0][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][0][2*k + 1];
    T5[chan1][ind1] += -1.0 * std::pow(-1.0, int(ij + jj + aj + bj)) * (2.0 * J + 1) * CGC6(bj, aj, J, jj, ij, J1) * T;
  }
  chan1 = Tmap[chan][16*hhpp + 8];
  ind1 = Tmap[chan][16*hhpp + 9];
  T6[chan1][ind1] = std::pow(-1.0, int(ij + jj - J)) * std::sqrt((2*J + 1)/(2*ij + 1)) * T;
  chan1 = Tmap[chan][16*hhpp + 10];
  ind1 = Tmap[chan][16*hhpp + 11];
  T7[chan1][ind1] = -1.0 * std::sqrt((2*J + 1)/(2*jj + 1)) * T;
  chan1 = Tmap[chan][16*hhpp + 12];
  ind1 = Tmap[chan][16*hhpp + 13];
  T8[chan1][ind1] = std::pow(-1.0, int(aj + bj - J)) * std::sqrt((2*J + 1)/(2*aj + 1)) * T;
  chan1 = Tmap[chan][16*hhpp + 14];
  ind1 = Tmap[chan][16*hhpp + 15];
  T9[chan1][ind1] = -1.0 * std::sqrt((2*J + 1)/(2*bj + 1)) * T;
}

double Doubles_1::get_TJ(const Model_Space &Space, const Channels &Chan, int &chan, int &hhpp, int &i, int &j, int &a, int &b) const
{
  double tempt;
  double J, ij, jj, aj, bj, J1;
  int chan1, ind1;
  J = 0.5 * Chan.qnums1[chan].j;
  ij = 0.5 * Space.qnums[i].j;
  jj = 0.5 * Space.qnums[j].j;
  aj = 0.5 * Space.qnums[a].j;
  bj = 0.5 * Space.qnums[b].j;
  tempt = T1[chan][hhpp];
  for(int k = 0; k < TJnum[chan][4*hhpp + 0]; ++k){
    chan1 = TJmap[chan][hhpp][0][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][0][2*k + 1];
    tempt += std::pow(-1.0, int(aj + bj - J)) * (2.0 * J1 + 1) * CGC6(bj, aj, J, ij, jj, J1) * T2[chan1][ind1];
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 1]; ++k){
    chan1 = TJmap[chan][hhpp][1][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][1][2*k + 1];
    tempt += std::pow(-1.0, int(ij + jj - J)) * (2.0 * J1 + 1) * CGC6(aj, bj, J, jj, ij, J1) * T3[chan1][ind1];
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 2]; ++k){
    chan1 = TJmap[chan][hhpp][2][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][2][2*k + 1];
    tempt += -1.0 * (2.0 * J1 + 1) * CGC6(aj, bj, J, ij, jj, J1) * T4[chan1][ind1];
  }
  for(int k = 0; k < TJnum[chan][4*hhpp + 3]; ++k){
    chan1 = TJmap[chan][hhpp][3][2*k];
    J1 = 0.5 * Chan.qnums2[chan1].j;
    ind1 = TJmap[chan][hhpp][3][2*k + 1];
    tempt += -1.0 * std::pow(-1.0, int(ij + jj + aj + bj)) * (2.0 * J1 + 1) * CGC6(bj, aj, J, jj, ij, J1) * T5[chan1][ind1];
  }
  chan1 = Tmap[chan][16*hhpp + 8];
  ind1 = Tmap[chan][16*hhpp + 9];
  tempt += std::pow(-1.0, int(ij + jj - J)) * std::sqrt((2*ij + 1)/(2*J + 1)) * T6[chan1][ind1];
  chan1 = Tmap[chan][16*hhpp + 10];
  ind1 = Tmap[chan][16*hhpp + 11];
  tempt += -1.0 * std::sqrt((2*jj + 1)/(2*J + 1)) * T7[chan1][ind1];
  chan1 = Tmap[chan][16*hhpp + 12];
  ind1 = Tmap[chan][16*hhpp + 13];
  tempt += std::pow(-1.0, int(aj + bj - J)) * std::sqrt((2*aj + 1)/(2*J + 1)) * T8[chan1][ind1];
  chan1 = Tmap[chan][16*hhpp + 14];
  ind1 = Tmap[chan][16*hhpp + 15];
  tempt += -1.0 * std::sqrt((2*bj + 1)/(2*J + 1)) * T9[chan1][ind1];
  return tempt;
}

void Doubles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  double fac1 = 1.0, fac2 = 0.0;
  char N = 'N';
  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], T1[i], Q11[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(T1[i], Ints.S_ME1.V20[i], Q21[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }
}

void Doubles_1::set_T_2J(const Channels &Chan, Interactions &Ints)
{
  double fac1 = 1.0, fac2 = 0.0;
  char N = 'N';
  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], T1[i], Q11[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(T1[i], Ints.S_ME1.V20[i], Q21[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }
}


void Singles_1::copy_Singles_1(const Channels &Chan, const Singles_1 &S)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    for(int i = 0; i < nhp1; ++i){
      for(int j = 0; j < 3; ++j){ Tmap[3*i + j] = S.Tmap[3*i + j]; }
      Evec[i] = S.Evec[i];
      T1[i] = S.T1[i];
      S4[i] = S.S4[i];
      for(int j = 0; j < nhp1; ++j){
	for(int k = 0; k < 18; ++k){ Tmap2[18 * (nhp1 * i + j) + k] = S.Tmap2[18 * (nhp1 * i + j) + k]; }
	for(int k = 0; k < 9; ++k){
	  TJnum[9 * (nhp1 * i + j) + k] = S.TJnum[9 * (nhp1 * i + j) + k];
	  for(int l = 0; l < 2 * TJnum[9 * (nhp1 * i + j) + k]; ++l){
	    TJmap[9 * (nhp1 * i + j) + k][l] = S.TJmap[9 * (nhp1 * i + j) + k][l];
	  }
	}
	S3[nhp1 * i +j] = S.S3[nhp1 * i +j];
      }
    }
  }
  if(nhh1 != 0){
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = S.Q31[i];
      Qmap3[2*i] = S.Qmap3[2*i];
      Qmap3[2*i + 1] = S.Qmap3[2*i + 1];
    }
  }
  if(npp1 != 0){
    for(int i = 0; i < npp1; ++i){
      Q41[i] = S.Q41[i];
      Qmap4[2*i] = S.Qmap4[2*i];
      Qmap4[2*i + 1] = S.Qmap4[2*i + 1];
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E1[i][j] = S.E1[i][j];
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q61[i][j] = S.Q61[i][j];
	Qmap6[i][2*j] = S.Qmap6[i][2*j];
	Qmap6[i][2*j + 1] = S.Qmap6[i][2*j + 1];
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q51[i][j] = S.Q51[i][j];
	Qmap5[i][2*j] = S.Qmap5[i][2*j];
	Qmap5[i][2*j + 1] = S.Qmap5[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = S.T2[i][j];
	T3[i][j] = S.T3[i][j];
      }
    }
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E6[i][j] = S.E6[i][j];
	E7[i][j] = S.E7[i][j];
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E8[i][j] = S.E8[i][j];
	E9[i][j] = S.E9[i][j];
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = S.Q11[i][j];
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = S.Q21[i][j];
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q32[i][j] = S.Q32[i][j];
	S2[i][j] = S.S2[i][j];
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q42[i][j] = S.Q42[i][j];
	S1[i][j] = S.S1[i][j];
      }
    }
    length = nhhp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q62[i][j] = S.Q62[i][j];
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q52[i][j] = S.Q52[i][j];
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Qmap1[i][2*j] = S.Qmap1[i][2*j];
	Qmap1[i][2*j + 1] = S.Qmap1[i][2*j + 1];
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Qmap2[i][2*j] = S.Qmap2[i][2*j];
	Qmap2[i][2*j + 1] = S.Qmap2[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E2[i][j] = S.E2[i][j];
	E3[i][j] = S.E3[i][j];
	E4[i][j] = S.E4[i][j];
	E5[i][j] = S.E5[i][j];
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = S.Q12[i][j];
	Q22[i][j] = S.Q22[i][j];
      }
    }
  }
}

Singles_1::Singles_1(const Channels &Chan, const Singles_1 &S)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    TJnum = new int[9 * nhp1 * nhp1];
    TJmap = new int*[9 * nhp1 * nhp1];
    Tmap = new int[3 * nhp1];
    Evec = new double[nhp1];
    T1 = new double[nhp1];
    Tmap2 = new int[18 * nhp1 * nhp1];
    S3 = new double[nhp1 * nhp1];
    S4 = new double[nhp1];
    for(int i = 0; i < nhp1; ++i){
      for(int j = 0; j < 3; ++j){ Tmap[3*i + j] = S.Tmap[3*i + j]; }
      Evec[i] = S.Evec[i];
      T1[i] = S.T1[i];
      S4[i] = S.S4[i];
      for(int j = 0; j < nhp1; ++j){
	for(int k = 0; k < 18; ++k){ Tmap2[18 * (nhp1 * i + j) + k] = S.Tmap2[18 * (nhp1 * i + j) + k]; }
	for(int k = 0; k < 9; ++k){
	  TJnum[9 * (nhp1 * i + j) + k] = S.TJnum[9 * (nhp1 * i + j) + k];
	  TJmap[9 * (nhp1 * i + j) + k] = new int[2 * TJnum[9 * (nhp1 * i + j) + k]];
	  for(int l = 0; l < 2 * TJnum[9 * (nhp1 * i + j) + k]; ++l){
	    TJmap[9 * (nhp1 * i + j) + k][l] = S.TJmap[9 * (nhp1 * i + j) + k][l];
	  }
	}
	S3[nhp1 * i +j] = S.S3[nhp1 * i +j];
      }
    }
  }
  if(nhh1 != 0){
    Q31 = new double[nhh1];
    Qmap3 = new int[2 * nhh1];
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = S.Q31[i];
      Qmap3[2*i] = S.Qmap3[2*i];
      Qmap3[2*i + 1] = S.Qmap3[2*i + 1];
    }
  }
  if(npp1 != 0){
    Q41 = new double[npp1];
    Qmap4 = new int[2 * npp1];
    for(int i = 0; i < npp1; ++i){
      Q41[i] = S.Q41[i];
      Qmap4[2*i] = S.Qmap4[2*i];
      Qmap4[2*i + 1] = S.Qmap4[2*i + 1];
    }
  }

  T2 = new double*[Chan.size3];
  T3 = new double*[Chan.size3];

  E1 = new double*[Chan.size1];
  E2 = new double*[Chan.size2];
  E3 = new double*[Chan.size2];
  E4 = new double*[Chan.size2];
  E5 = new double*[Chan.size2];
  E6 = new double*[Chan.size3];
  E7 = new double*[Chan.size3];
  E8 = new double*[Chan.size3];
  E9 = new double*[Chan.size3];

  S1 = new double*[Chan.size3];
  S2 = new double*[Chan.size3];

  Q11 = new double*[Chan.size3];
  Q21 = new double*[Chan.size3];
  Qmap1 = new int*[Chan.size3];
  Qmap2 = new int*[Chan.size3];
  Q12 = new double*[Chan.size2];
  Q22 = new double*[Chan.size2];

  Q32 = new double*[Chan.size3];
  Q42 = new double*[Chan.size3];

  Q51 = new double*[Chan.size1];
  Q61 = new double*[Chan.size1];
  Qmap5 = new int*[Chan.size1];
  Qmap6 = new int*[Chan.size1];
  Q52 = new double*[Chan.size3];
  Q62 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      E1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E1[i][j] = S.E1[i][j];
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q61[i] = new double[length];
      Qmap6[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q61[i][j] = S.Q61[i][j];
	Qmap6[i][2*j] = S.Qmap6[i][2*j];
	Qmap6[i][2*j + 1] = S.Qmap6[i][2*j + 1];
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q51[i] = new double[length];
      Qmap5[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q51[i][j] = S.Q51[i][j];
	Qmap5[i][2*j] = S.Qmap5[i][2*j];
	Qmap5[i][2*j + 1] = S.Qmap5[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = S.T2[i][j];
	T3[i][j] = S.T3[i][j];
      }
    }
    length = nhpp * nh;
    if(length != 0){
      E6[i] = new double[length];
      E7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E6[i][j] = S.E6[i][j];
	E7[i][j] = S.E7[i][j];
      }
    }
    length = nhhp * np;
    if(length != 0){
      E8[i] = new double[length];
      E9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E8[i][j] = S.E8[i][j];
	E9[i][j] = S.E9[i][j];
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Q11[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = S.Q11[i][j];
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Q21[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = S.Q21[i][j];
      }
    }
    length = nh * nh;
    if(length != 0){
      Q32[i] = new double[length];
      S2[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q32[i][j] = S.Q32[i][j];
	S2[i][j] = S.S2[i][j];
      }
    }
    length = np * np;
    if(length != 0){
      Q42[i] = new double[length];
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q42[i][j] = S.Q42[i][j];
	S1[i][j] = S.S1[i][j];
      }
    }
    length = nhhp * nh;
    if(length != 0){
      Q62[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q62[i][j] = S.Q62[i][j];
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q52[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q52[i][j] = S.Q52[i][j];
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap1[i][2*j] = S.Qmap1[i][2*j];
	Qmap1[i][2*j + 1] = S.Qmap1[i][2*j + 1];
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap2[i][2*j] = S.Qmap2[i][2*j];
	Qmap2[i][2*j + 1] = S.Qmap2[i][2*j + 1];
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      E2[i] = new double[length];
      E3[i] = new double[length];
      E4[i] = new double[length];
      E5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E2[i][j] = S.E2[i][j];
	E3[i][j] = S.E3[i][j];
	E4[i][j] = S.E4[i][j];
	E5[i][j] = S.E5[i][j];
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      Q12[i] = new double[length];
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = S.Q12[i][j];
	Q22[i][j] = S.Q22[i][j];
      }
    }
  }
}

Singles_1::Singles_1(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  int length, length0, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  State tb0, tb;
  int chan0, chan, ind, ind1, key1, key2, jmin0, jmin, count0, count;
  int hind1, hind2, hind3, pind1, pind2, pind3;
  int h1, hpp1, p1, hhp1, hp1, hp2, pp1, hh1;

  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    TJnum = new int[9 * nhp1 * nhp1];
    TJmap = new int*[9 * nhp1 * nhp1];
    Tmap = new int[3 * nhp1];
    Evec = new double[nhp1];
    T1 = new double[nhp1];
    Tmap2 = new int[18 * nhp1 * nhp1];
    S3 = new double[nhp1 * nhp1];
    S4 = new double[nhp1];
    for(int i = 0; i < nhp1; ++i){
      for(int j = 0; j < 3; ++j){ Tmap[3*i + j] = -1; }
      Evec[i] = 0.0;
      T1[i] = 0.0;
      S4[i] = 0.0;
      for(int j = 0; j < nhp1; ++j){
	for(int k = 0; k < 18; ++k){ Tmap2[18 * (nhp1 * i + j) + k] = -1; }
	for(int k = 0; k < 9; ++k){ TJnum[9 * (nhp1 * i + j) + k] = 0; }
	S3[nhp1 * i +j] = 0.0;
      }
    }
  }
  if(nhh1 != 0){
    Q31 = new double[nhh1];
    Qmap3 = new int[2 * nhh1];
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = 0.0;
      Qmap3[2*i] = -1;
      Qmap3[2*i + 1] = -1;
    }
  }
  if(npp1 != 0){
    Q41 = new double[npp1];
    Qmap4 = new int[2 * npp1];
    for(int i = 0; i < npp1; ++i){
      Q41[i] = 0.0;
      Qmap4[2*i] = -1;
      Qmap4[2*i + 1] = -1;
    }
  }

  T2 = new double*[Chan.size3];
  T3 = new double*[Chan.size3];

  E1 = new double*[Chan.size1];
  E2 = new double*[Chan.size2];
  E3 = new double*[Chan.size2];
  E4 = new double*[Chan.size2];
  E5 = new double*[Chan.size2];
  E6 = new double*[Chan.size3];
  E7 = new double*[Chan.size3];
  E8 = new double*[Chan.size3];
  E9 = new double*[Chan.size3];

  S1 = new double*[Chan.size3];
  S2 = new double*[Chan.size3];

  Q11 = new double*[Chan.size3];
  Q21 = new double*[Chan.size3];
  Qmap1 = new int*[Chan.size3];
  Qmap2 = new int*[Chan.size3];
  Q12 = new double*[Chan.size2];
  Q22 = new double*[Chan.size2];

  Q32 = new double*[Chan.size3];
  Q42 = new double*[Chan.size3];

  Q51 = new double*[Chan.size1];
  Q61 = new double*[Chan.size1];
  Qmap5 = new int*[Chan.size1];
  Qmap6 = new int*[Chan.size1];
  Q52 = new double*[Chan.size3];
  Q62 = new double*[Chan.size3];

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      E1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E1[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      Q61[i] = new double[length];
      Qmap6[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q61[i][j] = 0.0;
	Qmap6[i][2*j] = -1;
	Qmap6[i][2*j + 1] = -1;
      }
    }
    length = nhp * npp;
    if(length != 0){
      Q51[i] = new double[length];
      Qmap5[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Q51[i][j] = 0.0;
	Qmap5[i][2*j] = -1;
	Qmap5[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      T2[i] = new double[length];
      T3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
      }
    }
    length = nhpp * nh;
    if(length != 0){
      E6[i] = new double[length];
      E7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E6[i][j] = 0.0;
	E7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      E8[i] = new double[length];
      E9[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E8[i][j] = 0.0;
	E9[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Q11[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Q21[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      Q32[i] = new double[length];
      S2[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q32[i][j] = 0.0;
	S2[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      Q42[i] = new double[length];
      S1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q42[i][j] = 0.0;
	S1[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(length != 0){
      Q62[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q62[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      Q52[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q52[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      Qmap1[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap1[i][2*j] = -1;
	Qmap1[i][2*j + 1] = -1;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      Qmap2[i] = new int[2 * length];
      for(int j = 0; j < length; ++j){
	Qmap2[i][2*j] = -1;
	Qmap2[i][2*j + 1] = -1;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      E2[i] = new double[length];
      E3[i] = new double[length];
      E4[i] = new double[length];
      E5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	E2[i][j] = 0.0;
	E3[i][j] = 0.0;
	E4[i][j] = 0.0;
	E5[i][j] = 0.0;
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      Q12[i] = new double[length];
      Q22[i] = new double[length];
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
	Q22[i][j] = 0.0;
      }
    }
  }

  for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
    hind1 = Chan.hp1vec[Chan.ind0][2*hp1];
    pind1 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
    ind = Chan.indvec[hind1];
    Tmap[3 * hp1] = ind;
    key1 = Chan.p_map[ind][pind1];
    key2 = Chan.h_map[ind][hind1];
    ind1 = key1 * Chan.nh[ind] + key2;
    Tmap[3 * hp1 + 1] = ind1;
    ind1 = key2 * Chan.np[ind] + key1;
    Tmap[3 * hp1 + 2] = ind1;
  }

  length = Chan.nhp1[Chan.ind0] * Chan.nhp1[Chan.ind0];
  length0 = Chan.nhp1[Chan.ind0];
  for(int hphp = 0; hphp < length; ++hphp){
    hp2 = int(hphp%length0);
    hp1 = int((hphp - hp2)/length0);
    hind1 = Chan.hp1vec[Chan.ind0][2*hp1];
    pind2 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
    hind2 = Chan.hp1vec[Chan.ind0][2*hp2];
    pind1 = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
    if(Parameters.basis == "finite_J"){
      //E1 = ((ij)(ab))
      jmin0 = abs(Space.qnums[hind1].j - Space.qnums[hind2].j);
      if(abs(Space.qnums[pind1].j - Space.qnums[pind2].j) > jmin0){ jmin0 = abs(Space.qnums[pind1].j - Space.qnums[pind2].j); }
      plus(tb0, Space.qnums[hind1], Space.qnums[hind2]);
      if(Space.qnums[pind1].j + Space.qnums[pind2].j < tb0.j){ tb0.j = Space.qnums[pind1].j + Space.qnums[pind2].j; }
      TJnum[9 * hphp + 0] = int(0.5 * (tb0.j - jmin0) + 1);
      TJmap[9 * hphp + 0] = new int[2 * int(0.5 * (tb0.j - jmin0) + 1)];
      TJnum[9 * hphp + 5] = int(0.5 * (tb0.j - jmin0) + 1);
      TJmap[9 * hphp + 5] = new int[2 * int(0.5 * (tb0.j - jmin0) + 1)];
      TJnum[9 * hphp + 6] = int(0.5 * (tb0.j - jmin0) + 1);
      TJmap[9 * hphp + 6] = new int[2 * int(0.5 * (tb0.j - jmin0) + 1)];
      TJnum[9 * hphp + 7] = int(0.5 * (tb0.j - jmin0) + 1);
      TJmap[9 * hphp + 7] = new int[2 * int(0.5 * (tb0.j - jmin0) + 1)];
      TJnum[9 * hphp + 8] = int(0.5 * (tb0.j - jmin0) + 1);
      TJmap[9 * hphp + 8] = new int[2 * int(0.5 * (tb0.j - jmin0) + 1)];
      count0 = 0;
      while(tb0.j >= jmin0){
	chan0 = ChanInd_2b_dir(Parameters.basis, Space, tb0);
	key1 = Chan.hh_map[chan0][Hash2(hind1, hind2, Space.indtot)];
	key2 = Chan.pp_map[chan0][Hash2(pind1, pind2, Space.indtot)];
	ind = key1 * Chan.npp[chan0] + key2;
	TJmap[9 * hphp + 0][2*count0] = chan0;
	TJmap[9 * hphp + 0][2*count0 + 1] = ind;
	//E6 = ((jab)(i))
	chan = Chan.indvec[hind1];
	key1 = Chan.hpp_map[chan][int(0.5 * tb0.j * std::pow(Space.indtot, 3)) + Hash3(hind2, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[chan][hind1];
	ind = key1 * Chan.nh[chan] + key2;
	TJmap[9 * hphp + 5][2*count0] = chan;
	TJmap[9 * hphp + 5][2*count0 + 1] = ind;
	//E7 = ((iab)(j))
	chan = Chan.indvec[hind2];
	key1 = Chan.hpp_map[chan][int(0.5 * tb0.j * std::pow(Space.indtot, 3)) + Hash3(hind1, pind1, pind2, Space.indtot)];
	key2 = Chan.h_map[chan][hind2];
	ind = key1 * Chan.nh[chan] + key2;
	TJmap[9 * hphp + 6][2*count0] = chan;
	TJmap[9 * hphp + 6][2*count0 + 1] = ind;
	//E8 = ((ijb)(a))
	chan = Chan.indvec[pind1];
	key1 = Chan.hhp_map[chan][int(0.5 * tb0.j * std::pow(Space.indtot, 3)) + Hash3(hind1, hind2, pind2, Space.indtot)];
	key2 = Chan.p_map[chan][pind1];
	ind = key1 * Chan.np[chan] + key2;
	TJmap[9 * hphp + 7][2*count0] = chan;
	TJmap[9 * hphp + 7][2*count0 + 1] = ind;
	//E9 = ((ija)(b))
	chan = Chan.indvec[pind2];
	key1 = Chan.hhp_map[chan][int(0.5 * tb0.j * std::pow(Space.indtot, 3)) + Hash3(hind1, hind2, pind1, Space.indtot)];
	key2 = Chan.p_map[chan][pind2];
	ind = key1 * Chan.np[chan] + key2;
	TJmap[9 * hphp + 8][2*count0] = chan;
	TJmap[9 * hphp + 8][2*count0 + 1] = ind;

	tb0.j -= 2;
	++count0;
      }
      //E2 = ((ia)(jb)')
      jmin = abs(Space.qnums[hind1].j - Space.qnums[pind1].j);
      if(abs(Space.qnums[pind2].j - Space.qnums[hind2].j) > jmin){ jmin = abs(Space.qnums[pind2].j - Space.qnums[hind2].j); }
      minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
      if(Space.qnums[pind2].j + Space.qnums[hind2].j < tb.j){ tb.j = Space.qnums[pind2].j + Space.qnums[hind2].j; }
      TJnum[9 * hphp + 1] = int(0.5 * (tb.j - jmin) + 1);
      TJmap[9 * hphp + 1] = new int[2 * int(0.5 * (tb.j - jmin) + 1)];
      count = 0;
      while(tb.j >= jmin){
	chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind = key1 * Chan.nhp2[chan] + key2;
	TJmap[9 * hphp + 1][2*count] = chan;
	TJmap[9 * hphp + 1][2*count + 1] = ind;
	tb.j -= 2;
	++count;
      }
      //E3 = ((jb)(ia)')
      jmin = abs(Space.qnums[hind2].j - Space.qnums[pind2].j);
      if(abs(Space.qnums[pind1].j - Space.qnums[hind1].j) > jmin){ jmin = abs(Space.qnums[pind1].j - Space.qnums[hind1].j); }
      minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
      if(Space.qnums[pind1].j + Space.qnums[hind1].j < tb.j){ tb.j = Space.qnums[pind1].j + Space.qnums[hind1].j; }
      TJnum[9 * hphp + 2] = int(0.5 * (tb.j - jmin) + 1);
      TJmap[9 * hphp + 2] = new int[2 * int(0.5 * (tb.j - jmin) + 1)];
      count = 0;
      while(tb.j >= jmin){
	chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind = key1 * Chan.nhp2[ind] + key2;
	TJmap[9 * hphp + 2][2*count] = chan;
	TJmap[9 * hphp + 2][2*count + 1] = ind;
	tb.j -= 2;
	++count;
      }
      //E4 = ((ib)(ja)')
      jmin = abs(Space.qnums[hind1].j - Space.qnums[pind2].j);
      if(abs(Space.qnums[pind1].j - Space.qnums[hind2].j) > jmin){ jmin = abs(Space.qnums[pind1].j - Space.qnums[hind2].j); }
      minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
      if(Space.qnums[pind1].j + Space.qnums[hind2].j < tb.j){ tb.j = Space.qnums[pind1].j + Space.qnums[hind2].j; }
      TJnum[9 * hphp + 3] = int(0.5 * (tb.j - jmin) + 1);
      TJmap[9 * hphp + 3] = new int[2 * int(0.5 * (tb.j - jmin) + 1)];
      count = 0;
      while(tb.j >= jmin){
	chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind = key1 * Chan.nhp2[ind] + key2;
	TJmap[9 * hphp + 3][2*count] = chan;
	TJmap[9 * hphp + 3][2*count + 1] = ind;
	tb.j -= 2;
	++count;
      }
      //E5 = ((ja)(ib)')
      jmin = abs(Space.qnums[hind2].j - Space.qnums[pind1].j);
      if(abs(Space.qnums[pind2].j - Space.qnums[hind1].j) > jmin){ jmin = abs(Space.qnums[pind2].j - Space.qnums[hind1].j); }
      minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      if(Space.qnums[pind2].j + Space.qnums[hind1].j < tb.j){ tb.j = Space.qnums[pind2].j + Space.qnums[hind1].j; }
      TJnum[9 * hphp + 4] = int(0.5 * (tb.j - jmin) + 1);
      TJmap[9 * hphp + 4] = new int[2 * int(0.5 * (tb.j - jmin) + 1)];
      count = 0;
      while(tb.j >= jmin){
	chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind = key1 * Chan.nhp2[ind] + key2;
	TJmap[9 * hphp + 4][2*count] = chan;
	TJmap[9 * hphp + 4][2*count + 1] = ind;
	tb.j -= 2;
	++count;
      }
    }
    else{
      for(int k = 0; k < 9; ++k){ TJmap[9 * hphp + k] = new int[0]; }
      //E1 = ((ij)(ab))
      plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
      chan = ChanInd_2b_dir(Parameters.basis, Space, tb);
      key1 = Chan.hh_map[chan][Hash2(hind1, hind2, Space.indtot)];
      key2 = Chan.pp_map[chan][Hash2(pind1, pind2, Space.indtot)];
      ind = key1 * Chan.npp[chan] + key2;
      Tmap2[18 * hphp] = chan;
      Tmap2[18 * hphp + 1] = ind;
      //E2 = ((ia)(jb)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
      chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[chan][Hash2(hind1, pind1, Space.indtot)];
      key2 = Chan.hp2_map[chan][Hash2(hind2, pind2, Space.indtot)];
      ind = key1 * Chan.nhp2[chan] + key2;
      Tmap2[18 * hphp + 2] = chan;
      Tmap2[18 * hphp + 3] = ind;
      //E3 = ((jb)(ia)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
      chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[chan][Hash2(hind2, pind2, Space.indtot)];
      key2 = Chan.hp2_map[chan][Hash2(hind1, pind1, Space.indtot)];
      ind = key1 * Chan.nhp2[chan] + key2;
      Tmap2[18 * hphp + 4] = chan;
      Tmap2[18 * hphp + 5] = ind;
      //E4 = ((ib)(ja)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
      chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[chan][Hash2(hind1, pind2, Space.indtot)];
      key2 = Chan.hp2_map[chan][Hash2(hind2, pind1, Space.indtot)];
      ind = key1 * Chan.nhp2[chan] + key2;
      Tmap2[18 * hphp + 6] = chan;
      Tmap2[18 * hphp + 7] = ind;
      //E5 = ((ja)(ib)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      chan = ChanInd_2b_cross(Parameters.basis, Space, tb);
      key1 = Chan.hp1_map[chan][Hash2(hind2, pind1, Space.indtot)];
      key2 = Chan.hp2_map[chan][Hash2(hind1, pind2, Space.indtot)];
      ind = key1 * Chan.nhp2[chan] + key2;
      Tmap2[18 * hphp + 8] = chan;
      Tmap2[18 * hphp + 9] = ind;
      //E6 = ((jab)(i))
      chan = Chan.indvec[hind1];
      key1 = Chan.hpp_map[chan][Hash3(hind2, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[chan][hind1];
      ind = key1 * Chan.nh[chan] + key2;
      Tmap2[18 * hphp + 10] = chan;
      Tmap2[18 * hphp + 11] = ind;
      //E7 = ((iab)(j))
      chan = Chan.indvec[hind2];
      key1 = Chan.hpp_map[chan][Hash3(hind1, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[chan][hind2];
      ind = key1 * Chan.nh[chan] + key2;
      Tmap2[18 * hphp + 12] = chan;
      Tmap2[18 * hphp + 13] = ind;
      //E8 = ((ijb)(a))
      chan = Chan.indvec[pind1];
      key1 = Chan.hhp_map[chan][Hash3(hind1, hind2, pind2, Space.indtot)];
      key2 = Chan.p_map[chan][pind1];
      ind = key1 * Chan.np[chan] + key2;
      Tmap2[18 * hphp + 14] = chan;
      Tmap2[18 * hphp + 15] = ind;
      //E9 = ((ija)(b))
      chan = Chan.indvec[pind2];
      key1 = Chan.hhp_map[chan][Hash3(hind1, hind2, pind1, Space.indtot)];
      key2 = Chan.p_map[chan][pind2];
      ind = key1 * Chan.np[chan] + key2;
      Tmap2[18 * hphp + 16] = chan;
      Tmap2[18 * hphp + 17] = ind;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    length = Chan.nhpp1[i] * Chan.nh[i];
    length0 = Chan.nh[i];
    for(int hpph = 0; hpph < length; ++hpph){
      h1 = int(hpph%length0);
      hpp1 = int((hpph - h1)/length0);
      hind2 = Chan.hpp1vec[i][3*hpp1];
      pind1 = Chan.hpp1vec[i][3*hpp1 + 1];
      pind2 = Chan.hpp1vec[i][3*hpp1 + 2];
      hind1 = Chan.hvec[i][h1];
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind1].j - Space.qnums[pind1].j);
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp1[ind] + key2;
	  Qmap1[i][2 * hpph] = ind;
	  Qmap1[i][2 * hpph + 1] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	Qmap1[i][2 * hpph] = ind;
	Qmap1[i][2 * hpph + 1] = ind1;
      }
    }
    length = Chan.nhhp1[i] * Chan.np[i];
    length0 = Chan.np[i];
    for(int hhpp = 0; hhpp < length; ++hhpp){
      p1 = int(hhpp%length0);
      hhp1 = int((hhpp - p1)/length0);
      hind2 = Chan.hhp1vec[i][3*hhp1];
      hind1 = Chan.hhp1vec[i][3*hhp1 + 1];
      pind2 = Chan.hhp1vec[i][3*hhp1 + 2];
      pind1 = Chan.pvec[i][p1];
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind1].j - Space.qnums[pind1].j);
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp1[ind] + key2;
	  Qmap2[i][2 * hhpp] = ind;
	  Qmap2[i][2 * hhpp + 1] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(hind1, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	Qmap2[i][2 * hhpp] = ind;
	Qmap2[i][2 * hhpp + 1] = ind1;
      }
    }
  }

  for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){
    hind1 = Chan.hh1vec[Chan.ind0][2*hh1];
    hind2 = Chan.hh1vec[Chan.ind0][2*hh1 + 1];
    ind = Chan.indvec[hind1];
    key1 = Chan.h_map[ind][hind2];
    key2 = Chan.h_map[ind][hind1];
    ind1 = key1 * Chan.nh[ind] + key2;
    Qmap3[2 * hh1] = ind;
    Qmap3[2 * hh1 + 1] = ind1;
  }
  for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){
    pind1 = Chan.pp1vec[Chan.ind0][2*pp1];
    pind2 = Chan.pp1vec[Chan.ind0][2*pp1 + 1];
    ind = Chan.indvec[pind1];
    key1 = Chan.p_map[ind][pind1];
    key2 = Chan.p_map[ind][pind2];
    ind1 = key1 * Chan.np[ind] + key2;
    Qmap4[2 * pp1] = ind;
    Qmap4[2 * pp1 + 1] = ind1;
  }

  for(int i = 0; i < Chan.size1; ++i){
    length = Chan.nhp[i] * Chan.npp[i];
    length0 = Chan.npp[i];
    for(int hppp = 0; hppp < length; ++hppp){
      pp1 = int(hppp%length0);
      hp1 = int((hppp - pp1)/length0);
      hind1 = Chan.hpvec[i][2*hp1];
      pind1 = Chan.hpvec[i][2*hp1 + 1];
      pind2 = Chan.ppvec[i][2*pp1];
      pind3 = Chan.ppvec[i][2*pp1 + 1];
      ind = Chan.indvec[pind1];
      key1 = Chan.hpp_map[ind][Hash3(hind1, pind2, pind3, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Qmap5[i][2 * hppp] = ind;
      Qmap5[i][2 * hppp + 1] = ind1;
    }
    length = Chan.nhh[i] * Chan.nhp[i];
    length0 = Chan.nhp[i];
    for(int hhhp = 0; hhhp < length; ++hhhp){
      hp1 = int(hhhp%length0);
      hh1 = int((hhhp - hp1)/length0);
      hind1 = Chan.hhvec[i][2*hh1];
      hind2 = Chan.hhvec[i][2*hh1 + 1];
      hind3 = Chan.hpvec[i][2*hp1];
      pind1 = Chan.hpvec[i][2*hp1 + 1];
      ind = Chan.indvec[hind3];
      key1 = Chan.hhp_map[ind][Hash3(hind1, hind2, pind1, Space.indtot)];
      key2 = Chan.h_map[ind][hind3];
      ind1 = key1 * Chan.nh[ind] + key2;
      Qmap6[i][2 * hhhp] = ind;
      Qmap6[i][2 * hhhp + 1] = ind1;
    }
  }
}

void Singles_1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    for(int i = 0; i < 9 * nhp1 * nhp1; ++i){
      delete[] TJmap[i];
    }
    delete[] TJnum;
    delete[] TJmap;
    delete[] Tmap;
    delete[] Evec;
    delete[] T1;
    delete[] Tmap2;
    delete[] S3;
    delete[] S4;
  }
  if(nhh1 != 0){
    delete[] Q31;
    delete[] Qmap3;
  }
  if(npp1 != 0){
    delete[] Q41;
    delete[] Qmap4;
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(nhh * npp != 0){
      delete[] E1[i];
    }
    if(nhh * nhp != 0){
      delete[] Q61[i];
      delete[] Qmap6[i];
    }
    if(nhp * npp != 0){
      delete[] Q51[i];
      delete[] Qmap5[i];
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(nh * np != 0){
      delete[] T2[i];
      delete[] T3[i];
    }
    if(nhpp * nh != 0){
      delete[] E6[i];
      delete[] E7[i];
    }
    if(nhhp * np != 0){
      delete[] E8[i];
      delete[] E9[i];
    }
    if(nhpp1 * nh != 0){
      delete[] Q11[i];
    }
    if(nhhp1 * np != 0){
      delete[] Q21[i];
    }
    if(nh != 0){
      delete[] Q32[i];
      delete[] S2[i];
    }
    if(np != 0){
      delete[] Q42[i];
      delete[] S1[i];
    }
    if(nhpp * np != 0){
      delete[] Q52[i];
    }
    if(nhhp * nh != 0){
      delete[] Q62[i];
    }
    if(nhpp1 * nh != 0){
      delete[] Qmap1[i];
    }
    if(nhhp1 * np != 0){
      delete[] Qmap2[i];
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp1 * nhp2 != 0){
      delete[] E2[i];
      delete[] E3[i];
      delete[] E4[i];
      delete[] E5[i];
    }
    if(nhp1 != 0){
      delete[] Q12[i];
      delete[] Q22[i];
    }
  }
  delete[] T2;
  delete[] T3;
  delete[] E1;
  delete[] E2;
  delete[] E3;
  delete[] E4;
  delete[] E5;
  delete[] E6;
  delete[] E7;
  delete[] E8;
  delete[] E9;
  delete[] S1;
  delete[] S2;
  delete[] Q11;
  delete[] Q21;
  delete[] Qmap1;
  delete[] Qmap2;
  delete[] Q12;
  delete[] Q22;
  delete[] Q32;
  delete[] Q42;
  delete[] Q51;
  delete[] Q61;
  delete[] Qmap5;
  delete[] Qmap6;
  delete[] Q52;
  delete[] Q62;
}

void Singles_1::zero(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    for(int i = 0; i < nhp1; ++i){
      T1[i] = 0.0;
      S4[i] = 0.0;
      for(int j = 0; j < nhp1; ++j){
	S3[nhp1 * i +j] = 0.0;
      }
    }
  }
  if(nhh1 != 0){
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = 0.0;
    }
  }
  if(npp1 != 0){
    for(int i = 0; i < npp1; ++i){
      Q41[i] = 0.0;
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E1[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q61[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q51[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
      }
    }
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E6[i][j] = 0.0;
	E7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E8[i][j] = 0.0;
	E9[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q32[i][j] = 0.0;
	S2[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q42[i][j] = 0.0;
	S1[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q52[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q62[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E2[i][j] = 0.0;
	E3[i][j] = 0.0;
	E4[i][j] = 0.0;
	E5[i][j] = 0.0;
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
	Q22[i][j] = 0.0;
      }
    }
  }
}

void Singles_1::zero1(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1;
  nhp1 = Chan.nhp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  if(nhp1 != 0){
    for(int i = 0; i < nhp1; ++i){
      S4[i] = 0.0;
      for(int j = 0; j < nhp1; ++j){
	S3[nhp1 * i +j] = 0.0;
      }
    }
  }
  if(nhh1 != 0){
    for(int i = 0; i < nhh1; ++i){
      Q31[i] = 0.0;
    }
  }
  if(npp1 != 0){
    for(int i = 0; i < npp1; ++i){
      Q41[i] = 0.0;
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E1[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q61[i][j] = 0.0;
      }
    }
    length = nhp * npp;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q51[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    length = nh * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	T2[i][j] = 0.0;
	T3[i][j] = 0.0;
      }
    }
    length = nhpp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E6[i][j] = 0.0;
	E7[i][j] = 0.0;
      }
    }
    length = nhhp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E8[i][j] = 0.0;
	E9[i][j] = 0.0;
      }
    }
    length = nhpp1 * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q11[i][j] = 0.0;
      }
    }
    length = nhhp1 * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q21[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q32[i][j] = 0.0;
	S2[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q42[i][j] = 0.0;
	S1[i][j] = 0.0;
      }
    }
    length = nhpp * np;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q52[i][j] = 0.0;
      }
    }
    length = nhhp * nh;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q62[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    length = nhp1 * nhp2;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	E2[i][j] = 0.0;
	E3[i][j] = 0.0;
	E4[i][j] = 0.0;
	E5[i][j] = 0.0;
      }
    }
    length = nhp1 * nhp1;
    if(length != 0){
      for(int j = 0; j < length; ++j){
	Q12[i][j] = 0.0;
	Q22[i][j] = 0.0;
      }
    }
  }
}

void Singles_1::set_T(int i, double T)
{
  T1[i] = T;
  T2[Tmap[3*i]][Tmap[3*i + 1]] = T;
  T3[Tmap[3*i]][Tmap[3*i + 2]] = T;
}

double Singles_1::get_T(int i) const
{
  double tempt = T1[i];
  tempt += T2[Tmap[3*i]][Tmap[3*i + 1]];
  tempt += T3[Tmap[3*i]][Tmap[3*i + 2]];
  return tempt;
}

void Singles_1::set_TJ(const Model_Space &Space, int &hp, int &i, int &a, double T)
{
  T1[hp] = T;
  T2[Tmap[3*hp]][Tmap[3*hp + 1]] = T / std::sqrt(Space.qnums[a].j + 1);
  T3[Tmap[3*hp]][Tmap[3*hp + 1]] = T / std::sqrt(Space.qnums[i].j + 1);
}

double Singles_1::get_TJ(const Model_Space &Space, int &hp, int &i, int &a) const
{
  double tempt = T1[hp];
  tempt += T2[Tmap[3*hp]][Tmap[3*hp + 1]] * std::sqrt(Space.qnums[a].j + 1);
  tempt += T3[Tmap[3*hp]][Tmap[3*hp + 2]] * std::sqrt(Space.qnums[i].j + 1);
  return tempt;
}

void Singles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
    for(int hp2 = 0; hp2 < Chan.nhp1[Chan.ind0]; ++hp2){
      int ind0 = hp1 * Chan.nhp1[Chan.ind0] + hp2;
      double T = T1[hp1] * T1[hp2];
      E1[Tmap2[18 * ind0 + 0]][Tmap2[18 * ind0 + 1]] = T;
      E2[Tmap2[18 * ind0 + 2]][Tmap2[18 * ind0 + 3]] = T;
      E3[Tmap2[18 * ind0 + 4]][Tmap2[18 * ind0 + 5]] = T;
      E4[Tmap2[18 * ind0 + 6]][Tmap2[18 * ind0 + 7]] = T;
      E5[Tmap2[18 * ind0 + 8]][Tmap2[18 * ind0 + 9]] = T;
      E6[Tmap2[18 * ind0 + 10]][Tmap2[18 * ind0 + 11]] = T;
      E7[Tmap2[18 * ind0 + 12]][Tmap2[18 * ind0 + 13]] = T;
      E8[Tmap2[18 * ind0 + 14]][Tmap2[18 * ind0 + 15]] = T;
      E9[Tmap2[18 * ind0 + 16]][Tmap2[18 * ind0 + 17]] = T;
    }
  }

  double fac1 = 1.0, fac2 = 0.0;
  int one = 1;
  char N = 'N';
  for(int i = 0; i < Chan.size3; ++i){
    int h = Chan.nh[i];
    int p = Chan.np[i];
    int hpp1 = Chan.nhpp1[i];
    int hhp1 = Chan.nhhp1[i];
    if(h == 0 || p == 0){ continue; }
    if(hpp1 != 0){ dgemm_NN(Ints.S_ME1.V13[i], T2[i], Q11[i], &hpp1, &h, &p, &fac1, &fac2, &N, &N); }
    if(hhp1 != 0){ dgemm_NN(Ints.S_ME1.V14[i], T3[i], Q21[i], &hhp1, &p, &h, &fac1, &fac2, &N, &N); }
    for(int hpp2 = 0; hpp2 < hpp1; ++hpp2){
      for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	int ind0 = hpp2 * Chan.nh[i] + h1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hhp2 = 0; hhp2 < hhp1; ++hhp2){
      for(int p1 = 0; p1 < Chan.np[i]; ++p1){
	int ind0 = hhp2 * Chan.np[i] + p1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }

  int hp = Chan.nhp1[Chan.ind0];
  int hh = Chan.nhh1[Chan.ind0];
  int pp = Chan.npp1[Chan.ind0];
  if(hp != 0 && hh != 0){ dgemm_NN(Ints.S_ME1.V15[Chan.ind0], T1, Q31, &hh, &one, &hp, &fac1, &fac2, &N, &N); }
  if(hp != 0 && pp != 0){ dgemm_NN(Ints.S_ME1.V16[Chan.ind0], T1, Q41, &pp, &one, &hp, &fac1, &fac2, &N, &N); }
  for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){ Q32[Qmap3[2 * hh1]][Qmap3[2 * hh1 + 1]] = Q31[hh1]; }
  for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){ Q42[Qmap4[2 * pp1]][Qmap4[2 * pp1 + 1]] = Q41[pp1]; }

  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], E1[i], Q51[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(E1[i], Ints.S_ME1.V20[i], Q61[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q52[Qmap5[i][2 * ind0]][Qmap5[i][2 * ind0 + 1]] = Q51[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q62[Qmap6[i][2 * ind0]][Qmap6[i][2 * ind0 + 1]] = Q61[i][ind0];
      }
    }
  }
}

void Singles_1::set_T_2J(const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double T, ij, jj, aj, bj;
  int ind, chan1, ind1, i, j, a, b, J, J1;
  for(int hp1 = 0; hp1 < Chan.nhp1[Chan.ind0]; ++hp1){
    for(int hp2 = 0; hp2 < Chan.nhp1[Chan.ind0]; ++hp2){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      b = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      j = Chan.hp1vec[Chan.ind0][2*hp2];
      a = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
      ij = 0.5 * Space.qnums[i].j;
      bj = 0.5 * Space.qnums[b].j;
      jj = 0.5 * Space.qnums[j].j;
      aj = 0.5 * Space.qnums[a].j;
      ind = hp1 * Chan.nhp1[Chan.ind0] + hp2;
      T = T1[hp1] * T1[hp2] / std::sqrt((2.0 * aj + 1) * (2.0 * bj + 1));
      for(int k = 0; k < TJnum[9 * ind + 0]; ++k){
	chan1 = TJmap[9 * ind + 0][2*k];
	J = 0.5 * Chan.qnums1[chan1].j;
	ind1 = TJmap[9 * ind + 0][2*k + 1];
	E1[chan1][ind1] = T * std::pow(-1.0, int(aj + bj - J));
	chan1 = TJmap[9 * ind + 5][2*k];
	ind1 = TJmap[9 * ind + 5][2*k + 1];
	E6[chan1][ind1] = T * std::sqrt((2.0 * J + 1)/(2.0 * ij + 1)) * std::pow(-1.0, int(ij + jj + aj + bj));
	chan1 = TJmap[9 * ind + 6][2*k];
	ind1 = TJmap[9 * ind + 6][2*k + 1];
	E7[chan1][ind1] = T * std::sqrt((2.0 * J + 1)/(2.0 * jj + 1)) * std::pow(-1.0, int(aj + bj - J));
	chan1 = TJmap[9 * ind + 7][2*k];
	ind1 = TJmap[9 * ind + 7][2*k + 1];
	E8[chan1][ind1] = T * std::sqrt((2.0 * J + 1)/(2.0 * aj + 1));
	chan1 = TJmap[9 * ind + 8][2*k];
	ind1 = TJmap[9 * ind + 8][2*k + 1];
	E9[chan1][ind1] = T * std::sqrt((2.0 * J + 1)/(2.0 * bj + 1)) * std::pow(-1.0, int(ij + jj - J));
	for(int l = 0; l < TJnum[9 * ind + 1]; ++l){
	  chan1 = TJmap[9 * ind + 1][2*l];
	  J1 = 0.5 * Chan.qnums2[chan1].j;
	  ind1 = TJmap[9 * ind + 1][2*l + 1];
	  E2[chan1][ind1] += T * (2.0 * J + 1) * CGC6(bj, aj, J, ij, jj, J1);
	}
	for(int l = 0; l < TJnum[9 * ind + 2]; ++l){
	  chan1 = TJmap[9 * ind + 2][2*l];
	  J1 = 0.5 * Chan.qnums2[chan1].j;
	  ind1 = TJmap[9 * ind + 2][2*l + 1];
	  E3[chan1][ind1] += T * (2.0 * J + 1) * CGC6(aj, bj, J, jj, ij, J1) * std::pow(-1.0, int(ij + jj + aj + bj));
	}
	for(int l = 0; l < TJnum[9 * ind + 3]; ++l){
	  chan1 = TJmap[9 * ind + 3][2*l];
	  J1 = 0.5 * Chan.qnums2[chan1].j;
	  ind1 = TJmap[9 * ind + 3][2*l + 1];
	  E4[chan1][ind1] += T * (2.0 * J + 1) * CGC6(aj, bj, J, ij, jj, J1) * std::pow(-1.0, int(aj + bj - J));
	}
	for(int l = 0; l < TJnum[9 * ind + 4]; ++l){
	  chan1 = TJmap[9 * ind + 4][2*l];
	  J1 = 0.5 * Chan.qnums2[chan1].j;
	  ind1 = TJmap[9 * ind + 4][2*l + 1];
	  E5[chan1][ind1] += T * (2.0 * J + 1) * CGC6(bj, aj, J, jj, ij, J1) * std::pow(-1.0, int(ij + jj - J));
	}
      }
    }
  }

  double fac1 = 1.0, fac2 = 0.0;
  int one = 1;
  char N = 'N';
  for(int i = 0; i < Chan.size3; ++i){
    int h = Chan.nh[i];
    int p = Chan.np[i];
    int hpp1 = Chan.nhpp1[i];
    int hhp1 = Chan.nhhp1[i];
    if(h == 0 || p == 0){ continue; }
    if(hpp1 != 0){ dgemm_NN(Ints.S_ME1.V13[i], T2[i], Q11[i], &hpp1, &h, &p, &fac1, &fac2, &N, &N); }
    if(hhp1 != 0){ dgemm_NN(Ints.S_ME1.V14[i], T3[i], Q21[i], &hhp1, &p, &h, &fac1, &fac2, &N, &N); }
    for(int hpp2 = 0; hpp2 < hpp1; ++hpp2){
      for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	int ind0 = hpp2 * Chan.nh[i] + h1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hhp2 = 0; hhp2 < hhp1; ++hhp2){
      for(int p1 = 0; p1 < Chan.np[i]; ++p1){
	int ind0 = hhp2 * Chan.np[i] + p1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }

  int hp = Chan.nhp1[Chan.ind0];
  int hh = Chan.nhh1[Chan.ind0];
  int pp = Chan.npp1[Chan.ind0];
  if(hp != 0 && hh != 0){ dgemm_NN(Ints.S_ME1.V15[Chan.ind0], T1, Q31, &hh, &one, &hp, &fac1, &fac2, &N, &N); }
  if(hp != 0 && pp != 0){ dgemm_NN(Ints.S_ME1.V16[Chan.ind0], T1, Q41, &pp, &one, &hp, &fac1, &fac2, &N, &N); }
  for(int hh1 = 0; hh1 < Chan.nhh1[Chan.ind0]; ++hh1){ Q32[Qmap3[2 * hh1]][Qmap3[2 * hh1 + 1]] = Q31[hh1]; }
  for(int pp1 = 0; pp1 < Chan.npp1[Chan.ind0]; ++pp1){ Q42[Qmap4[2 * pp1]][Qmap4[2 * pp1 + 1]] = Q41[pp1]; }

  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.nhh[i];
    int pp = Chan.npp[i];
    int hp = Chan.nhp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    dgemm_NN(Ints.S_ME1.V19[i], E1[i], Q51[i], &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    dgemm_NN(E1[i], Ints.S_ME1.V20[i], Q61[i], &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
	int ind0 = hp1 * Chan.npp[i] + pp1;
	Q52[Qmap5[i][2 * ind0]][Qmap5[i][2 * ind0 + 1]] = Q51[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.nhp[i]; ++hp1){
	int ind0 = hh1 * Chan.nhp[i] + hp1;
	Q62[Qmap6[i][2 * ind0]][Qmap6[i][2 * ind0 + 1]] = Q61[i][ind0];
      }
    }
  }
}

Interactions::Interactions(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1 = Doubles_ME1(Chan);
  }
  else if(Parameters.approx == "singles"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
}

void Interactions::delete_struct(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1.delete_struct(Chan);;
  }
  else if(Parameters.approx == "singles"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
}

Doubles_ME1::Doubles_ME1(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  V1 = new double*[Chan.size1];
  V2 = new double*[Chan.size1];
  V3 = new double*[Chan.size2];
  V4 = new double*[Chan.size1];
  V5 = new double*[Chan.size3];
  V6 = new double*[Chan.size3];
  V7 = new double*[Chan.size3];
  V8 = new double*[Chan.size3];
  V9 = new double*[Chan.size2];
  V10 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      V1[i] = new double[npp * npp];
      for(int j = 0; j < npp * npp; ++j){ V1[i][j] = 0.0; }
    }
    if(nhh != 0){
      V2[i] = new double[nhh * nhh];
      for(int j = 0; j < nhh * nhh; ++j){ V2[i][j] = 0.0; }
    }
    if(npp * nhh != 0){
      V4[i] = new double[npp * nhh];
      for(int j = 0; j < npp * nhh; ++j){ V4[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      V5[i] = new double[nh * nhpp];
      V6[i] = new double[nh * nhpp];
      for(int j = 0; j < nh * nhpp; ++j){ V5[i][j] = 0.0; V6[i][j] = 0.0; }
    }
    if(np * nhhp != 0){
      V7[i] = new double[np * nhhp];
      V8[i] = new double[np * nhhp];
      for(int j = 0; j < np * nhhp; ++j){ V7[i][j] = 0.0; V8[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      V3[i] = new double[nhp2 * nhp2];
      for(int j = 0; j < nhp2 * nhp2; ++j){ V3[i][j] = 0.0; }
    }
    if(nhp2 * nhp1){
      V9[i] = new double[nhp2 * nhp1];
      V10[i] = new double[nhp2 * nhp1];
      for(int j = 0; j < nhp2 * nhp1; ++j){ V9[i][j] = 0.0; V10[i][j] = 0.0; }
    }
  }
}

void Doubles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      delete[] V1[i];
    }
    if(nhh != 0){
      delete[] V2[i];
    }
    if(npp * nhh != 0){
      delete[] V4[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      delete[] V5[i];
      delete[] V6[i];
    }
    if(np * nhhp != 0){
      delete[] V7[i];
      delete[] V8[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      delete[] V3[i];
    }
    if(nhp2 * nhp1){
      delete[] V9[i];
      delete[] V10[i];
    }
  }
  delete[] V1;
  delete[] V2;
  delete[] V3;
  delete[] V4;
  delete[] V5;
  delete[] V6;
  delete[] V7;
  delete[] V8;
  delete[] V9;
  delete[] V10;
}

Singles_ME1::Singles_ME1(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  V11 = new double*[Chan.size3];
  V12 = new double*[Chan.size3];
  V17 = new double*[Chan.size3];
  V18 = new double*[Chan.size3];
  V13 = new double*[Chan.size3];
  V14 = new double*[Chan.size3];
  V20 = new double*[Chan.size1];
  V19 = new double*[Chan.size1];
  V16 = new double*[Chan.size2];
  V15 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      V20[i] = new double[npp * nhp];
      for(int j = 0; j < npp * nhp; ++j){ V20[i][j] = 0.0; }
    }
    if(nhp * nhh != 0){
      V19[i] = new double[nhp * nhh];
      for(int j = 0; j < nhp * nhh; ++j){ V19[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      V11[i] = new double[np * nhpp];
      V17[i] = new double[nhpp * np];
      for(int j = 0; j < np * nhpp; ++j){ V11[i][j] = 0.0; V17[i][j] = 0.0; }
    }
    if(nh * nhhp != 0){
      V12[i] = new double[nh * nhhp];
      V18[i] = new double[nhhp * nh];
      for(int j = 0; j < nh * nhhp; ++j){ V12[i][j] = 0.0; V18[i][j] = 0.0; }
    }
    if(nhpp1 * np != 0){
      V13[i] = new double[nhpp1 * np];
      for(int j = 0; j < nhpp1 * np; ++j){ V13[i][j] = 0.0; }
    }
    if(nhhp1 * nh != 0){
      V14[i] = new double[nhhp1 * nh];
      for(int j = 0; j < nhhp1 * nh; ++j){ V14[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      V16[i] = new double[npp1 * nhp1];
      for(int j = 0; j < npp1 * nhp1; ++j){ V16[i][j] = 0.0; }
    }
    if(nhh1 * nhp1 != 0){
      V15[i] = new double[nhh1 * nhp1];
      for(int j = 0; j < nhh1 * nhp1; ++j){ V15[i][j] = 0.0; }
    }
  }
}

void Singles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      delete[] V20[i];
    }
    if(nhp * nhh != 0){
      delete[] V19[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      delete[] V11[i];
      delete[] V17[i];
    }
    if(nh * nhhp != 0){
      delete[] V12[i];
      delete[] V18[i];
    }
    if(nhpp1 * np != 0){
      delete[] V13[i];
    }
    if(nhhp1 * nh != 0){
      delete[] V14[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      delete[] V16[i];
    }
    if(nhh1 * nhp1 != 0){
      delete[] V15[i];
    }
  }
  delete[] V11;
  delete[] V12;
  delete[] V17;
  delete[] V18;
  delete[] V13;
  delete[] V14;
  delete[] V20;
  delete[] V19;
  delete[] V16;
  delete[] V15;
}

//Function to setup Channels
Channels::Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  State state;
  int key;
  int count0, ind1, ind2, h, p, hh, pp, hp, hp1, hp2, hh1, pp1, hpp, hhp, hpp1, hhp1, hhh, ppp;
  State *qnumstemp = new State[Space.indtot]; // max number of qnums groups
  int *hnum = new int[Space.indtot]; // count holes in each qnums group
  int *pnum = new int[Space.indtot]; // count particles in each qnums group
  indvec = new int[Space.indtot]; // index of qnums group for each state
  for(int i = 0; i < Space.indtot; ++i){ hnum[i] = 0; pnum[i] = 0; }

  qnumstemp[0] = Space.qnums[0];
  if(Parameters.basis != "finite_J"){ qnumstemp[0].j = 0; }
  indvec[0] = 0;
  if( Space.qnums[0].type == "hole" ){ ++hnum[0]; }
  else{ ++pnum[0]; }

  // count # of qnums groups and # of hs and ps in each qnums group, and fill indvec
  count0 = 1;
  for(int i = 1; i < Space.indtot; ++i){
    state = Space.qnums[i];
    if(Parameters.basis != "finite_J"){ state.j = 0; }
    for(int k = 0; k < count0; ++k){
      if( equal(state, qnumstemp[k]) ){
	indvec[i] = k;
	if(Space.qnums[i].type == "hole"){ ++hnum[k]; }
	else{ ++pnum[k]; }
	goto stop;
      }
      if(k == count0 - 1){
	qnumstemp[count0] = Space.qnums[i];
	if(Parameters.basis != "finite_J"){ qnumstemp[count0].j = 0; }
	indvec[i] = count0;
	if(Space.qnums[i].type == "hole"){ ++hnum[count0]; }
	else{ ++pnum[count0]; }
	++count0;
	break;
      }
    }
  stop:;
  }
  size3 = count0;

  // allocate memory for Hvec and Pvec, reset hnum and pnum
  qnums3 = new State[size3];
  hvec = new int*[size3];
  pvec = new int*[size3];
  h_map = new std::unordered_map<int,int>[size3];
  p_map = new std::unordered_map<int,int>[size3];
  nh = new int[size3];
  np = new int[size3];
  for(int i = 0; i < size3; ++i){
    h = hnum[i];
    p = pnum[i];
    qnums3[i] = qnumstemp[i];
    if(h != 0){ hvec[i] = new int[h]; }
    if(p != 0){ pvec[i] = new int[p]; }
    nh[i] = 0;
    np[i] = 0;
  }
  delete[] qnumstemp;
  delete[] hnum;
  delete[] pnum;

  // place states in appropriate Hvec or Pvec position
  for(int i = 0; i < Space.indtot; ++i){
    for(int k = 0; k < size3; ++k){
      state = Space.qnums[i];
      if(Parameters.basis != "finite_J"){ state.j = 0; }
      if( equal(state, qnums3[k]) ){
	if(Space.qnums[i].type == "hole"){
	  hvec[k][nh[k]] = i;
	  h_map[k][i] = nh[k];
	  ++nh[k];
	}
	else{
	  pvec[k][np[k]] = i;
	  p_map[k][i] = np[k];
	  ++np[k];
	}
	break;
      }
    }
  }

  /*for(int i = 0; i < size3; ++i){
    std::cout << "Chan3: " << i << ", " << qnums3[i].par << " " << qnums3[i].ml << " " << qnums3[i].m << std::endl;
    for(int j = 0; j < nh[i]; ++j){ std::cout << hvec[i][j] << " "; }
    std::cout << std::endl;
    for(int j = 0; j < np[i]; ++j){ std::cout << pvec[i][j] << " "; }
    std::cout << std::endl;
    }*/

  size1 = Space.size_2b;
  size2 = Space.size_2b;
  //std::cout << " Size1 = " << size1 << ", Size2 = " << size2 << ", Size3 = " << size3 << std::endl;

  qnums1 = new State[size1];
  qnums2 = new State[size2];

  hhvec = new int*[size1];
  ppvec = new int*[size1];
  hpvec = new int*[size1];
  hp1vec = new int*[size2];
  hp2vec = new int*[size2];
  hh1vec = new int*[size2];
  pp1vec = new int*[size2];

  hh_map = new std::unordered_map<int,int>[size1];
  pp_map = new std::unordered_map<int,int>[size1];
  hp_map = new std::unordered_map<int,int>[size1];
  hp1_map = new std::unordered_map<int,int>[size2];
  hp2_map = new std::unordered_map<int,int>[size2];
  hh1_map = new std::unordered_map<int,int>[size2];
  pp1_map = new std::unordered_map<int,int>[size2];

  nhh = new int[size1];
  npp = new int[size1];
  nhp = new int[size1];
  nhp1 = new int[size2];
  nhp2 = new int[size2];
  nhh1 = new int[size2];
  npp1 = new int[size2];

  for(int chan1 = 0; chan1 < size1; ++chan1){
    nhh[chan1] = 0;
    npp[chan1] = 0;
    nhp[chan1] = 0;
  }
  for(int chan2 = 0; chan2 < size2; ++chan2){
    nhp1[chan2] = 0;
    nhp2[chan2] = 0;
    nhh1[chan2] = 0;
    npp1[chan2] = 0;
  }

  int jmin;
  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[i].j - Space.qnums[j].j);
	plus(state, Space.qnums[i], Space.qnums[j]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  qnums1[ind1] = state;
	  ++nhh[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[i], Space.qnums[j]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  qnums2[ind2] = state;
	  ++nhh1[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[i], Space.qnums[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	qnums1[ind1] = state;
	++nhh[ind1];
	minus(state, Space.qnums[i], Space.qnums[j]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	qnums2[ind2] = state;
	++nhh1[ind2];
      }
    }
  }

  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[a].j - Space.qnums[b].j);
	plus(state, Space.qnums[a], Space.qnums[b]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  qnums1[ind1] = state;
	  ++npp[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[a], Space.qnums[b]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  qnums2[ind2] = state;
	  ++npp1[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[a], Space.qnums[b]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	qnums1[ind1] = state;
	++npp[ind1];
	minus(state, Space.qnums[a], Space.qnums[b]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	qnums2[ind2] = state;
	++npp1[ind2];
      }
    }
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int a = Space.indhol; a < Space.indtot; ++a){
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[i].j - Space.qnums[a].j);
	plus(state, Space.qnums[i], Space.qnums[a]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  qnums1[ind1] = state;
	  ++nhp[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[i], Space.qnums[a]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  qnums2[ind2] = state;
	  ++nhp1[ind2];
	  state.j -= 2;
	}
	minus(state, Space.qnums[a], Space.qnums[i]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  qnums2[ind2] = state;
	  ++nhp2[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[i], Space.qnums[a]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	qnums1[ind1] = state;
	++nhp[ind1];      
	minus(state, Space.qnums[i], Space.qnums[a]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	qnums2[ind2] = state;
	++nhp1[ind2];
	minus(state, Space.qnums[a], Space.qnums[i]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	qnums2[ind2] = state;
	++nhp2[ind2];
      }
    }
  }

  for(int i = 0; i < size1; ++i){
    hh = nhh[i];
    pp = npp[i];
    hp = nhp[i];
    if(hh != 0){ hhvec[i] = new int[2 * hh]; }
    if(pp != 0){ ppvec[i] = new int[2 * pp]; }
    if(hp != 0){ hpvec[i] = new int[2 * hp]; }
    nhh[i] = 0;
    npp[i] = 0;
    nhp[i] = 0;
  }
  for(int i = 0; i < size2; ++i){
    hp1 = nhp1[i];
    hp2 = nhp2[i];
    hh1 = nhh1[i];
    pp1 = npp1[i];
    if(hp1 != 0){ hp1vec[i] = new int[2 * hp1]; }
    if(hp2 != 0){ hp2vec[i] = new int[2 * hp2]; }
    if(hh1 != 0){ hh1vec[i] = new int[2 * hh1]; }
    if(pp1 != 0){ pp1vec[i] = new int[2 * pp1]; }
    nhp1[i] = 0;
    nhp2[i] = 0;
    nhh1[i] = 0;
    npp1[i] = 0;
  }

  for(int i = 0; i < Space.indhol; ++i){
    for(int j = 0; j < Space.indhol; ++j){
      key = Hash2(i, j, Space.indtot);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[i].j - Space.qnums[j].j);
	plus(state, Space.qnums[i], Space.qnums[j]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  hhvec[ind1][2 * nhh[ind1]] = i;
	  hhvec[ind1][2 * nhh[ind1] + 1] = j;
	  hh_map[ind1][key] = nhh[ind1];
	  ++nhh[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[i], Space.qnums[j]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  hh1vec[ind2][2 * nhh1[ind2]] = i;
	  hh1vec[ind2][2 * nhh1[ind2] + 1] = j;
	  hh1_map[ind2][key] = nhh1[ind2];
	  ++nhh1[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[i], Space.qnums[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	hhvec[ind1][2 * nhh[ind1]] = i;
	hhvec[ind1][2 * nhh[ind1] + 1] = j;
	hh_map[ind1][key] = nhh[ind1];
	++nhh[ind1];
	minus(state, Space.qnums[i], Space.qnums[j]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	hh1vec[ind2][2 * nhh1[ind2]] = i;
	hh1vec[ind2][2 * nhh1[ind2] + 1] = j;
	hh1_map[ind2][key] = nhh1[ind2];
	++nhh1[ind2];
      }
    }
  }

  for(int a = Space.indhol; a < Space.indtot; ++a){
    for(int b = Space.indhol; b < Space.indtot; ++b){
      key = Hash2(a, b, Space.indtot);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[a].j - Space.qnums[b].j);
	plus(state, Space.qnums[a], Space.qnums[b]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  ppvec[ind1][2 * npp[ind1]] = a;
	  ppvec[ind1][2 * npp[ind1] + 1] = b;
	  pp_map[ind1][key] = npp[ind1];
	  ++npp[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[a], Space.qnums[b]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  pp1vec[ind2][2 * npp1[ind2]] = a;
	  pp1vec[ind2][2 * npp1[ind2] + 1] = b;
	  pp1_map[ind2][key] = npp1[ind2];
	  ++npp1[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[a], Space.qnums[b]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	ppvec[ind1][2 * npp[ind1]] = a;
	ppvec[ind1][2 * npp[ind1] + 1] = b;
	pp_map[ind1][key] = npp[ind1];
	++npp[ind1];      
	minus(state, Space.qnums[a], Space.qnums[b]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	pp1vec[ind2][2 * npp1[ind2]] = a;
	pp1vec[ind2][2 * npp1[ind2] + 1] = b;
	pp1_map[ind2][key] = npp1[ind2];
	++npp1[ind2];
      }
    }
  }

  /*for(int i = 0; i < size1; ++i){
    std::cout << "Chan1:  " << i << " " << qnums1[i].par << " " << qnums1[i].t << " " << qnums1[i].j << std::endl;
    for(int j = 0; j < nhh[i]; ++j){
      std::cout << hhvec[i][2*j] << "," << hhvec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    for(int j = 0; j < npp[i]; ++j){
      std::cout << ppvec[i][2*j] << "," << ppvec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    }*/

  for(int i = 0; i < Space.indhol; ++i){
    for(int a = Space.indhol; a < Space.indtot; ++a){
      key = Hash2(i, a, Space.indtot);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[i].j - Space.qnums[a].j);
	plus(state, Space.qnums[i], Space.qnums[a]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  hpvec[ind1][2 * nhp[ind1]] = i;
	  hpvec[ind1][2 * nhp[ind1] + 1] = a;
	  hp_map[ind1][key] = nhp[ind1];
	  ++nhp[ind1];
	  state.j -= 2;
	}
	minus(state, Space.qnums[i], Space.qnums[a]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  hp1vec[ind2][2 * nhp1[ind2]] = i;
	  hp1vec[ind2][2 * nhp1[ind2] + 1] = a;
	  hp1_map[ind2][key] = nhp1[ind2];
	  ++nhp1[ind2];
	  state.j -= 2;
	}
	minus(state, Space.qnums[a], Space.qnums[i]);
	while(state.j >= jmin){
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  hp2vec[ind2][2 * nhp2[ind2]] = i;
	  hp2vec[ind2][2 * nhp2[ind2] + 1] = a;
	  hp2_map[ind2][key] = nhp2[ind2];
	  ++nhp2[ind2];
	  state.j -= 2;
	}
      }
      else{
	plus(state, Space.qnums[i], Space.qnums[a]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	hpvec[ind1][2 * nhp[ind1]] = i;
	hpvec[ind1][2 * nhp[ind1] + 1] = a;
	hp_map[ind1][key] = nhp[ind1];
	++nhp[ind1];
	minus(state, Space.qnums[i], Space.qnums[a]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	hp1vec[ind2][2 * nhp1[ind2]] = i;
	hp1vec[ind2][2 * nhp1[ind2] + 1] = a;
	hp1_map[ind2][key] = nhp1[ind2];
	++nhp1[ind2];
	minus(state, Space.qnums[a], Space.qnums[i]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	hp2vec[ind2][2 * nhp2[ind2]] = i;
	hp2vec[ind2][2 * nhp2[ind2] + 1] = a;
	hp2_map[ind2][key] = nhp2[ind2];
	++nhp2[ind2];
      }
    }
  }

  /*for(int i = 0; i < size2; ++i){
    std::cout << "Chan2:  " << i << " " << qnums2[i].par << " " << qnums2[i].t << " " << qnums2[i].j << std::endl;
    std::cout << nhp1[i] << " " << nhp2[i] << std::endl;
    for(int j = 0; j < nhp1[i]; ++j){
      std::cout << hp1vec[i][2*j] << "," << hp1vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    for(int j = 0; j < nhp2[i]; ++j){
      std::cout << hp2vec[i][2*j] << "," << hp2vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    }*/

  // for singles case
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    state.t = 0; //tz
    state.m = 0; //jz
    state.par = 1; //par
    state.ml = 0; //lz
    state.j = 0; //j
    ind0 = ChanInd_2b_cross(Parameters.basis, Space, state);
  }

  hppvec = new int*[size3];
  hhpvec = new int*[size3];
  hpp1vec = new int*[size3];
  hhp1vec = new int*[size3];
  hhhvec = new int*[size3];
  pppvec = new int*[size3];

  hpp_map = new std::unordered_map<int,int>[size3];
  hhp_map = new std::unordered_map<int,int>[size3];
  hpp1_map = new std::unordered_map<int,int>[size3];
  hhp1_map = new std::unordered_map<int,int>[size3];
  hhh_map = new std::unordered_map<int,int>[size3];
  ppp_map = new std::unordered_map<int,int>[size3];

  nhpp = new int[size3];
  nhhp = new int[size3];
  nhpp1 = new int[size3];
  nhhp1 = new int[size3];
  nhhh = new int[size3];
  nppp = new int[size3];

  for(int i = 0; i < size3; ++i){
    nhpp[i] = 0;
    nhhp[i] = 0;
    nhpp1[i] = 0;
    nhhp1[i] = 0;
    nhhh[i] = 0;
    nppp[i] = 0;
  }

  for(int i = 0; i < size3; ++i){
    for(int j = 0; j < size3; ++j){
      h = nh[j];
      p = np[j];
      if(Parameters.basis == "finite_J"){
	jmin = abs(qnums3[i].j - qnums3[j].j);
	plus(state, qnums3[i], qnums3[j]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  hh = nhh[ind1];
	  pp = npp[ind1];
	  hp = nhp[ind1];
	  nhhp[i] += p * hh; // <p1p2|h1h2> -> <p1|h1h2p2>
	  nhpp[i] += h * pp; // <h1h2|p1p2> -> <h1|h2p1p2>
	  nhhp1[i] += h * hp; // <h1h2|h3p1> -> <h1|h2h3p1>
	  nhpp1[i] += p * hp; // <p1p2|h1p3> -> <p1|h1p3p2>
	  nhhh[i] += h * hh; // <p1h1|h2h3> -> <p1|h1h2h3>
	  nppp[i] += p * pp; // <h1p1|p2p3> -> <h1|p1p2p3>
	  state.j -= 2;
	}
      }
      else{
	plus(state, qnums3[i], qnums3[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	hh = nhh[ind1];
	pp = npp[ind1];
	hp = nhp[ind1];
	nhhp[i] += p * hh; // <p1p2|h1h2> -> <p1|h1h2p2>
	nhpp[i] += h * pp; // <h1h2|p1p2> -> <h1|h2p1p2>
	nhhp1[i] += h * hp; // <h1h2|h3p1> -> <h1|h2h3p1>
	nhpp1[i] += p * hp; // <p1p2|h1p3> -> <p1|h1p3p2>
	nhhh[i] += h * hh; // <p1h1|h2h3> -> <p1|h1h2h3>
	nppp[i] += p * pp; // <h1p1|p2p3> -> <h1|p1p2p3>
      }
    }
  }

  for(int i = 0; i < size3; ++i){
    hpp = nhpp[i];
    hhp = nhhp[i];
    hpp1 = nhpp1[i];
    hhp1 = nhhp1[i];
    hhh = nhhh[i];
    ppp = nppp[i];
    if(hpp != 0){ hppvec[i] = new int[3 * hpp]; }
    if(hhp != 0){ hhpvec[i] = new int[3 * hhp]; }
    if(hpp1 != 0){ hpp1vec[i] = new int[3 * hpp1]; }
    if(hhp1 != 0){ hhp1vec[i] = new int[3 * hhp1]; }
    if(hhh != 0){ hhhvec[i] = new int[3 * hhh]; }
    if(ppp != 0){ pppvec[i] = new int[3 * ppp]; }
    nhpp[i] = 0;
    nhhp[i] = 0;
    nhpp1[i] = 0;
    nhhp1[i] = 0;
    nhhh[i] = 0;
    nppp[i] = 0;
  }

  for(int i = 0; i < size3; ++i){
    for(int j = 0; j < size3; ++j){
      if(Parameters.basis == "finite_J"){
	jmin = abs(qnums3[i].j - qnums3[j].j);
	plus(state, qnums3[i], qnums3[j]);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1h2> -> <p1|h1h2p2>
	    for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], pvec[j][p1], Space.indtot);
	      hhpvec[i][3 * nhhp[i]] = hhvec[ind1][2 * hh1];
	      hhpvec[i][3 * nhhp[i] + 1] = hhvec[ind1][2 * hh1 + 1];
	      hhpvec[i][3 * nhhp[i] + 2] = pvec[j][p1];
	      hhp_map[i][key] = nhhp[i];
	      ++nhhp[i];
	    }
	  }
	  for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|p1p2> -> <h1|h2p1p2>
	    for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(hvec[j][h1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	      hppvec[i][3 * nhpp[i]] = hvec[j][h1];
	      hppvec[i][3 * nhpp[i] + 1] = ppvec[ind1][2 * pp1];
	      hppvec[i][3 * nhpp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	      hpp_map[i][key] = nhpp[i];
	      ++nhpp[i];
	    }
	  }
	  for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|h3p1> -> <h1|h2h3p1>
	    for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(hvec[j][h1], hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], Space.indtot);
	      hhp1vec[i][3 * nhhp1[i]] = hvec[j][h1];
	      hhp1vec[i][3 * nhhp1[i] + 1] = hpvec[ind1][2 * hp1];
	      hhp1vec[i][3 * nhhp1[i] + 2] = hpvec[ind1][2 * hp1 + 1];
	      hhp1_map[i][key] = nhhp1[i];
	      ++nhhp1[i];
	    }
	  }
	  for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1p3> -> <p1|h1p3p2>
	    for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], pvec[j][p1], Space.indtot);
	      hpp1vec[i][3 * nhpp1[i]] = hpvec[ind1][2 * hp1];
	      hpp1vec[i][3 * nhpp1[i] + 1] = hpvec[ind1][2 * hp1 + 1];
	      hpp1vec[i][3 * nhpp1[i] + 2] = pvec[j][p1];
	      hpp1_map[i][key] = nhpp1[i];
	      ++nhpp1[i];
	    }
	  }
	  for(int h1 = 0; h1 < nh[j]; ++h1){ // <p1h1|h2h3> -> <p1|h1h2h3>
	    for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(hvec[j][h1], hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], Space.indtot);
	      hhhvec[i][3 * nhhh[i]] = hvec[j][h1];
	      hhhvec[i][3 * nhhh[i] + 1] = hhvec[ind1][2 * hh1];
	      hhhvec[i][3 * nhhh[i] + 2] = hhvec[ind1][2 * hh1 + 1];
	      hhh_map[i][key] = nhhh[i];
	      ++nhhh[i];
	    }
	  }
	  for(int p1 = 0; p1 < np[j]; ++p1){ // <h1p1|p2p3> -> <h1|p1p2p3>
	    for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	      key = int(0.5 * state.j * std::pow(Space.indtot, 3)) + Hash3(pvec[j][p1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	      pppvec[i][3 * nppp[i]] = pvec[j][p1];
	      pppvec[i][3 * nppp[i] + 1] = ppvec[ind1][2 * pp1];
	      pppvec[i][3 * nppp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	      ppp_map[i][key] = nppp[i];
	      ++nppp[i];
	    }
	  }
	  state.j -= 2;
	}
      }
      else{
	plus(state, qnums3[i], qnums3[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1h2> -> <p1|h1h2p2>
	  for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	    key = Hash3(hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], pvec[j][p1], Space.indtot);
	    hhpvec[i][3 * nhhp[i]] = hhvec[ind1][2 * hh1];
	    hhpvec[i][3 * nhhp[i] + 1] = hhvec[ind1][2 * hh1 + 1];
	    hhpvec[i][3 * nhhp[i] + 2] = pvec[j][p1];
	    hhp_map[i][key] = nhhp[i];
	    ++nhhp[i];
	  }
	}
	for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|p1p2> -> <h1|h2p1p2>
	  for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	    key = Hash3(hvec[j][h1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	    hppvec[i][3 * nhpp[i]] = hvec[j][h1];
	    hppvec[i][3 * nhpp[i] + 1] = ppvec[ind1][2 * pp1];
	    hppvec[i][3 * nhpp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	    hpp_map[i][key] = nhpp[i];
	    ++nhpp[i];
	  }
	}
	for(int h1 = 0; h1 < nh[j]; ++h1){ // <h1h2|h3p1> -> <h1|h2h3p1>
	  for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	    key = Hash3(hvec[j][h1], hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], Space.indtot);
	    hhp1vec[i][3 * nhhp1[i]] = hvec[j][h1];
	    hhp1vec[i][3 * nhhp1[i] + 1] = hpvec[ind1][2 * hp1];
	    hhp1vec[i][3 * nhhp1[i] + 2] = hpvec[ind1][2 * hp1 + 1];
	    hhp1_map[i][key] = nhhp1[i];
	    ++nhhp1[i];
	  }
	}
	for(int p1 = 0; p1 < np[j]; ++p1){ // <p1p2|h1p3> -> <p1|h1p3p2>
	  for(int hp1 = 0; hp1 < nhp[ind1]; ++hp1){
	    key = Hash3(hpvec[ind1][2 * hp1], hpvec[ind1][2 * hp1 + 1], pvec[j][p1], Space.indtot);
	    hpp1vec[i][3 * nhpp1[i]] = hpvec[ind1][2 * hp1];
	    hpp1vec[i][3 * nhpp1[i] + 1] = hpvec[ind1][2 * hp1 + 1];
	    hpp1vec[i][3 * nhpp1[i] + 2] = pvec[j][p1];
	    hpp1_map[i][key] = nhpp1[i];
	    ++nhpp1[i];
	  }
	}
	for(int h1 = 0; h1 < nh[j]; ++h1){ // <p1h1|h2h3> -> <p1|h1h2h3>
	  for(int hh1 = 0; hh1 < nhh[ind1]; ++hh1){
	    key = Hash3(hvec[j][h1], hhvec[ind1][2 * hh1], hhvec[ind1][2 * hh1 + 1], Space.indtot);
	    hhhvec[i][3 * nhhh[i]] = hvec[j][h1];
	    hhhvec[i][3 * nhhh[i] + 1] = hhvec[ind1][2 * hh1];
	    hhhvec[i][3 * nhhh[i] + 2] = hhvec[ind1][2 * hh1 + 1];
	    hhh_map[i][key] = nhhh[i];
	    ++nhhh[i];
	  }
	}
	for(int p1 = 0; p1 < np[j]; ++p1){ // <h1p1|p2p3> -> <h1|p1p2p3>
	  for(int pp1 = 0; pp1 < npp[ind1]; ++pp1){
	    key = Hash3(pvec[j][p1], ppvec[ind1][2 * pp1], ppvec[ind1][2 * pp1 + 1], Space.indtot);
	    pppvec[i][3 * nppp[i]] = pvec[j][p1];
	    pppvec[i][3 * nppp[i] + 1] = ppvec[ind1][2 * pp1];
	    pppvec[i][3 * nppp[i] + 2] = ppvec[ind1][2 * pp1 + 1];
	    ppp_map[i][key] = nppp[i];
	    ++nppp[i];
	  }
	}
      }
    }
  }

  double memory = 0.0;
  int intsize = sizeof(int);
  int doubsize = sizeof(double);
  for(int i = 0; i < size1; ++i){
    memory += (7*16 * intsize + (2*7 + 1) * doubsize) * nhh[i] * npp[i]; // Tmap, Evec, T1, V4
    memory += (7 + 1) * doubsize * nhh[i] * nhh[i]; // S1, V2
    memory += doubsize * npp[i] * npp[i]; // V1
  }
  for(int i = 0; i < size3; ++i){
    memory += (2*7 + 2) * doubsize * nhpp[i] * nh[i]; // T6, T7, V5, V6
    memory += (2*7 + 2) * doubsize * nhhp[i] * np[i]; // T8, T9, V7, V8
    memory += 2*7 * doubsize * nh[i] * nh[i]; // S2, S3
    memory += 2*7 * doubsize * np[i] * np[i]; // S4, S5
  }
  for(int i = 0; i < size2; ++i){
    memory += (4*7 + 2) * doubsize * nhp1[i] * nhp2[i]; // T2, T3, T4, T5, V9, V10
    memory += (2*7 + 1) * doubsize * nhp2[i] * nhp2[i]; // S6, S7, V3
  }

  if(Parameters.approx == "singles"){
    memory += (7*3 * intsize + 7*3 * doubsize) * nhp1[ind0]; // Tmap, Evec, T1, S4
    memory += (7*18 * intsize + 7*doubsize) * nhp1[ind0] * nhp1[ind0]; // Tmap2, S3
    memory += (7*doubsize + 7*2 * intsize) * nhh1[ind0]; // Q31, Qmap3
    memory += (7*doubsize + 7*2 * intsize) * npp1[ind0]; // Q41, Qmap4
    for(int i = 0; i < size1; ++i){
      memory += (7*doubsize + 7*2 * intsize) * nhp[i] * npp[i]; // Q11, Qmap1
      memory += (7*doubsize + 7*2 * intsize) * nhh[i] * nhp[i]; // Q21, Qmap2
      memory += 7*doubsize * nhh[i] * npp[i]; // E1
      memory += ((7 + 1) * doubsize + 7*2 * intsize) * nhh[i] * nhp[i]; // Q61, Qmap6, V19
      memory += ((7 + 1) * doubsize + 7*2 * intsize) * nhp[i] * npp[i]; // Q51, Qmap5, V20
    }
    for(int i = 0; i < size3; ++i){
      memory += (7 + 2) * doubsize * nhpp[i] * np[i]; // Q12, V11, V17
      memory += (7 + 2) * doubsize * nhhp[i] * nh[i]; // Q22, V12, V18
      memory += 7*2 * doubsize * nh[i] * np[i]; // T2, T3
      memory += 7*2 * doubsize * nhpp[i] * nh[i]; // E6, E7
      memory += 7*2 * doubsize * nhhp[i] * np[i]; // E8, E9
      memory += 7*doubsize * nhpp1[i] * nh[i]; // Q11
      memory += 7*doubsize * nhhp1[i] * np[i]; // Q21
      memory += 7*2 * doubsize * nh[i] * nh[i]; // Q32, S2
      memory += 7*2 * doubsize * np[i] * np[i]; // Q42, S1
      memory += 7*doubsize * nhhp[i] * nh[i]; // Q62
      memory += 7*doubsize * nhpp[i] * np[i]; // Q52
      memory += 7*2 * intsize * nhpp1[i] * nh[i]; // Qmap1
      memory += 7*2 * intsize * nhhp1[i] * np[i]; // Qmap2
      memory += doubsize * nhpp1[i] * np[i]; // V13
      memory += doubsize * nhhp1[i] * nh[i]; // V14
    }
    for(int i = 0; i < size2; ++i){
      memory += 7*4 * doubsize * nhp1[i] * nhp2[i]; // E2, E3, E4, E5
      memory += 7*2 * doubsize * nhp1[i] * nhp1[i]; // Q12, Q22
      memory += doubsize * npp1[i] * nhp1[i]; // V16
      memory += doubsize * nhh1[i] * nhp1[i]; // V15
    }
  }

  if(Parameters.extra == -1 || Parameters.extra == 0 || Parameters.extra == 1){
    memory += (2 * doubsize + 6 * intsize) * nhp1[ind0]; // X_ia1, Map_ia, X_ai1, Map_ai
    memory += (doubsize + 3 * intsize) * npp1[ind0]; // X_ab1, Map_ab
    memory += (2 * doubsize + 3 * intsize) * nhh1[ind0]; // X_ij1, X1_ij1, Map_ij
    memory += (doubsize + 3 * intsize) * nhp1[ind0]; // X_ai1, Map_ai
    for(int i = 0; i < size3; ++i){
      memory += 4 * doubsize * nh[i] * np[i]; // X_ia2, X_ia3, X_ai2, X_ai3
      memory += 2 * doubsize * np[i] * np[i]; // X_ab2, X_ab3
      memory += 4 * doubsize * nh[i] * nh[i]; // X_ij2, X_ij3, X1_ij2, X1_ij3
      memory += (3 * doubsize + 20 * intsize) * np[i] * nhpp[i]; // X1_iabc1, X_iabc1, Map_iabc, X_abic1, Map_abic
      memory += (4 * doubsize + 18 * intsize) * nh[i] * nhhp[i]; // X1_ijka1, X_ijka1, Map_ijka, X2_iajk1, X_iajk1, Map_iajk
      memory += 2 * doubsize * nh[i] * nppp[i]; // X1_iabc2, X_abic2
      memory += 2 * doubsize * np[i] * nhhh[i]; // X1_ijka2, X2_iajk2
      memory += 4 * doubsize * np[i] * nhpp1[i]; // X1_iabc3, X_iabc3, X_abic3, X_abic4
      memory += 2 * doubsize * nh[i] * nhhp1[i]; // X2_iajk3, X2_iajk4
      memory += 3 * doubsize * np[i] * nppp[i]; // X1_abcd2, X1_abcd3, V_abcd
      memory += 3 * doubsize * nh[i] * nhhh[i]; // X_ijkl2, X_ijkl3, V_ijkl
      memory += 4 * doubsize * nh[i] * nhpp1[i]; // X1_iajb3, X1_iajb4, X3_iajb3, X_iajb3
      memory += 3 * doubsize * np[i] * nhhp1[i]; // X1_iajb2, X3_iajb2, X3_iajb5
    }
    for(int i = 0; i < size1; ++i){
      memory += doubsize * nhh[i] * npp[i]; // X_ijab1
      memory += (2 * doubsize + 8 * intsize) * nhh[i] * nhh[i]; // X_ijkl1, X_ijkl4, Map_ijkl
      memory += (2 * doubsize + 6 * intsize) * npp[i] * npp[i]; // X1_abcd1, X_abcd1, Map_abcd
      memory += 2 * doubsize * npp[i] * nhp[i]; // X_iabc5, X_abic7
      memory += 2 * doubsize * nhh[i] * nhp[i]; // X_ijka5, X2_iajk7
    }    
    for(int i = 0; i < size2; ++i){
      memory += doubsize * npp1[i] * nhp1[i]; // X_iabc4
      memory += doubsize * nhh1[i] * nhp1[i]; // X_ijka4
      memory += (3 * doubsize + 8 * intsize) * nhp2[i] * nhp2[i]; // X1_iajb1, X3_iajb1, X_iajb1, Map_iajb
      memory += 2 * doubsize * npp1[i] * nhp2[i]; // X_abic5, X_abic6
      memory += 2 * doubsize * nhh1[i] * nhp2[i]; // X2_iajk5, X2_iajk6
    }
  }
  std::cout << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl << std::endl;
}


void Channels::delete_struct()
{
  int h, p, hh, pp, hp, hp1, hp2, hh1, pp1, hpp, hhp, hpp1, hhp1, hhh, ppp;
  for(int i = 0; i < size3; ++i){
    h = nh[i];
    p = np[i];
    if(h != 0){ delete[] hvec[i]; }
    if(p != 0){ delete[] pvec[i]; }
  }
  delete[] indvec;
  delete[] qnums3;
  delete[] hvec;
  delete[] pvec;
  delete[] h_map;
  delete[] p_map;
  delete[] nh;
  delete[] np;
  delete[] qnums1;
  delete[] qnums2;

  for(int i = 0; i < size1; ++i){
    hh = nhh[i];
    pp = npp[i];
    hp = nhp[i];
    if(hh != 0){ delete[] hhvec[i]; }
    if(pp != 0){ delete[] ppvec[i]; }
    if(hp != 0){ delete[] hpvec[i]; }
  }
  delete[] hhvec;
  delete[] ppvec;
  delete[] hpvec;
  delete[] hh_map;
  delete[] pp_map;
  delete[] hp_map;
  delete[] nhh;
  delete[] npp;
  delete[] nhp;

  for(int i = 0; i < size2; ++i){
    hp1 = nhp1[i];
    hp2 = nhp2[i];
    hh1 = nhh1[i];
    pp1 = npp1[i];
    if(hp1 != 0){ delete[] hp1vec[i]; }
    if(hp2 != 0){ delete[] hp2vec[i]; }
    if(hh1 != 0){ delete[] hh1vec[i]; }
    if(pp1 != 0){ delete[] pp1vec[i]; }
  }
  delete[] hp1vec;
  delete[] hp2vec;
  delete[] hh1vec;
  delete[] pp1vec;
  delete[] hp1_map;
  delete[] hp2_map;
  delete[] hh1_map;
  delete[] pp1_map;
  delete[] nhp1;
  delete[] nhp2;
  delete[] nhh1;
  delete[] npp1;

  for(int i = 0; i < size3; ++i){
    hpp = nhpp[i];
    hhp = nhhp[i];
    hpp1 = nhpp1[i];
    hhp1 = nhhp1[i];
    hhh = nhhh[i];
    ppp = nppp[i];
    if(hpp != 0){ delete[] hppvec[i]; }
    if(hhp != 0){ delete[] hhpvec[i]; }
    if(hpp1 != 0){ delete[] hpp1vec[i]; }
    if(hhp1 != 0){ delete[] hhp1vec[i]; }
    if(hhh != 0){ delete[] hhhvec[i]; }
    if(ppp != 0){ delete[] pppvec[i]; }
  }
  delete[] hppvec;
  delete[] hhpvec;
  delete[] hpp1vec;
  delete[] hhp1vec;
  delete[] hhhvec;
  delete[] pppvec;
  delete[] hpp_map;
  delete[] hhp_map;
  delete[] hpp1_map;
  delete[] hhp1_map;
  delete[] hhh_map;
  delete[] ppp_map;
  delete[] nhpp;
  delete[] nhhp;
  delete[] nhpp1;
  delete[] nhhp1;
  delete[] nhhh;
  delete[] nppp;
}

CC_Eff::CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  int length0, length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1, nhhh, nppp;
  State tb;
  int ind, ind1, key1, key2, jmin;
  int hind1, hind2, hind3, hind4, pind1, pind2, pind3, pind4;
  int hhp1, h1, p1, hpp1, hh1, hh2, pp1, pp2, hp1, hp2;
  nhp1 = Chan.nhp1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  nhp1 = Chan.nhp1[Chan.ind0];
  length = nhp1;
  if(length != 0){
    X_ia1 = new double[length]; // (ia)
    Map_ia = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ia1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ia[3*j + k] = -1; }
    }
  }
  length = npp1;
  if(length != 0){
    X_ab1 = new double[length]; // (ba)
    Map_ab = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ab1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ab[3*j + k] = -1; }
    }
  }
  length = nhh1;
  if(nhh1 != 0){
    X_ij1 = new double[length]; // (ji)
    X1_ij1 = new double[length]; // (ji)
    Map_ij = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ij1[j] = 0.0;
      X1_ij1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ij[3*j + k] = -1; }
    }
  }
  length = nhp1;
  if(nhp1 != 0){
    X_ai1 = new double[length]; // (ia)
    Map_ai = new int[3 * length];
    for(int j = 0; j < length; ++j){
      X_ai1[j] = 0.0;
      for(int k = 0; k < 3; ++k){ Map_ai[3*j + k] = -1; }
    }
  }

  X_ia2 = new double*[Chan.size3];
  X_ia3 = new double*[Chan.size3];

  X_ab2 = new double*[Chan.size3];
  X_ab3 = new double*[Chan.size3];

  X_ij2 = new double*[Chan.size3];
  X_ij3 = new double*[Chan.size3];
  X1_ij2 = new double*[Chan.size3];
  X1_ij3 = new double*[Chan.size3];

  X_ai2 = new double*[Chan.size3];
  X_ai3 = new double*[Chan.size3];

  X_ijab1 = new double*[Chan.size1];

  X1_iabc1 = new double*[Chan.size3];
  X1_iabc2 = new double*[Chan.size3];
  X1_iabc3 = new double*[Chan.size3];
  X_iabc1 = new double*[Chan.size3];
  X_iabc3 = new double*[Chan.size3];
  X_iabc4 = new double*[Chan.size2];
  X_iabc5 = new double*[Chan.size1];
  Map_iabc = new int*[Chan.size3];

  X1_ijka1 = new double*[Chan.size3];
  X1_ijka2 = new double*[Chan.size3];
  X_ijka1 = new double*[Chan.size3];
  X_ijka4 = new double*[Chan.size2];
  X_ijka5 = new double*[Chan.size1];
  Map_ijka = new int*[Chan.size3];

  X1_abcd1 = new double*[Chan.size1];
  X1_abcd2 = new double*[Chan.size3];
  X1_abcd3 = new double*[Chan.size3];
  X_abcd1 = new double*[Chan.size1];
  V_abcd = new double*[Chan.size3];
  Map_abcd = new int*[Chan.size1];

  X_ijkl1 = new double*[Chan.size1];
  X_ijkl2 = new double*[Chan.size3];
  X_ijkl3 = new double*[Chan.size3];
  X_ijkl4 = new double*[Chan.size1];
  V_ijkl = new double*[Chan.size3];
  Map_ijkl = new int*[Chan.size1];

  X1_iajb1 = new double*[Chan.size2];
  X1_iajb2 = new double*[Chan.size3];
  X1_iajb3 = new double*[Chan.size3];
  X1_iajb4 = new double*[Chan.size3];
  X3_iajb1 = new double*[Chan.size2];
  X3_iajb2 = new double*[Chan.size3];
  X3_iajb3 = new double*[Chan.size3];
  X3_iajb5 = new double*[Chan.size3];
  X_iajb1 = new double*[Chan.size2];
  X_iajb3 = new double*[Chan.size3];
  Map_iajb = new int*[Chan.size2];

  X_abic1 = new double*[Chan.size3];
  X_abic2 = new double*[Chan.size3];
  X_abic3 = new double*[Chan.size3];
  X_abic4 = new double*[Chan.size3];
  X_abic5 = new double*[Chan.size2];
  X_abic6 = new double*[Chan.size2];
  X_abic7 = new double*[Chan.size1];
  Map_abic = new int*[Chan.size3];

  X2_iajk1 = new double*[Chan.size3];
  X2_iajk2 = new double*[Chan.size3];
  X2_iajk3 = new double*[Chan.size3];
  X2_iajk4 = new double*[Chan.size3];
  X2_iajk5 = new double*[Chan.size2];
  X2_iajk6 = new double*[Chan.size2];
  X2_iajk7 = new double*[Chan.size1];
  X_iajk1 = new double*[Chan.size3];
  Map_iajk = new int*[Chan.size3];

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    nhhh = Chan.nhhh[i];
    nppp = Chan.nppp[i];
    length = nh * np;
    if(length != 0){
      X_ia2[i] = new double[length]; // (i,a)
      X_ia3[i] = new double[length]; // (a,i)
      X_ai2[i] = new double[length]; // (a,i)
      X_ai3[i] = new double[length]; // (i,a)
      for(int j = 0; j < length; ++j){
	X_ia2[i][j] = 0.0;
	X_ia3[i][j] = 0.0;
	X_ai2[i][j] = 0.0;
	X_ai3[i][j] = 0.0;
      }
    }
    length = np * np;
    if(length != 0){
      X_ab2[i] = new double[length]; // (a,b)
      X_ab3[i] = new double[length]; // (b,a)
      for(int j = 0; j < length; ++j){
	X_ab2[i][j] = 0.0;
	X_ab3[i][j] = 0.0;
      }
    }
    length = nh * nh;
    if(length != 0){
      X_ij2[i] = new double[length]; // (i,j)
      X_ij3[i] = new double[length]; // (j,i)
      X1_ij2[i] = new double[length]; // (i,j)
      X1_ij3[i] = new double[length]; // (j,i)
      for(int j = 0; j < length; ++j){
	X_ij2[i][j] = 0.0;
	X_ij3[i][j] = 0.0;
	X1_ij2[i][j] = 0.0;
	X1_ij3[i][j] = 0.0;
      }
    }
    length = np * nhpp;
    if(length != 0){
      X1_iabc1[i] = new double[length];
      X_iabc1[i] = new double[length];
      Map_iabc[i] = new int[8 * length];
      X_abic1[i] = new double[length];
      Map_abic[i] = new int[12 * length];
      for(int j = 0; j < length; ++j){
	X1_iabc1[i][j] = 0.0;
	X_iabc1[i][j] = 0.0;
	X_abic1[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_iabc[i][8*j + k] = -1; }
	for(int k = 0; k < 12; ++k){ Map_abic[i][12*j + k] = -1; }
      }
    }
    length = nh * nhhp;
    if(length != 0){
      X1_ijka1[i] = new double[length];
      X_ijka1[i] = new double[length];
      Map_ijka[i] = new int[6 * length];
      X2_iajk1[i] = new double[length];
      X_iajk1[i] = new double[length];
      Map_iajk[i] = new int[12 * length];
      for(int j = 0; j < length; ++j){
	X1_ijka1[i][j] = 0.0;
	X_ijka1[i][j] = 0.0;
	X2_iajk1[i][j] = 0.0;
	X_iajk1[i][j] = 0.0;
	for(int k = 0; k < 6; ++k){ Map_ijka[i][6*j + k] = -1; }
	for(int k = 0; k < 12; ++k){ Map_iajk[i][12*j + k] = -1; }
      }
    }
    length = nh * nppp;
    if(length != 0){
      X1_iabc2[i] = new double[length];
      X_abic2[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_iabc2[i][j] = 0.0;
	X_abic2[i][j] = 0.0;
      }
    }
    length = np * nhhh;
    if(length != 0){
      X1_ijka2[i] = new double[length];
      X2_iajk2[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_ijka2[i][j] = 0.0;
	X2_iajk2[i][j] = 0.0;
      }
    }
    length = np * nhpp1;
    if(length != 0){
      X1_iabc3[i] = new double[length];
      X_iabc3[i] = new double[length];
      X_abic3[i] = new double[length];
      X_abic4[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_iabc3[i][j] = 0.0;      
	X_iabc3[i][j] = 0.0;
	X_abic3[i][j] = 0.0;
	X_abic4[i][j] = 0.0;
      }
    }
    length = nh * nhhp1;
    if(length != 0){
      X2_iajk3[i] = new double[length];
      X2_iajk4[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X2_iajk3[i][j] = 0.0;
	X2_iajk4[i][j] = 0.0;
      }
    }
    length = np * nppp;
    if(length != 0){
      X1_abcd2[i] = new double[length];
      X1_abcd3[i] = new double[length];
      V_abcd[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_abcd2[i][j] = 0.0;
	X1_abcd3[i][j] = 0.0;
	V_abcd[i][j] = 0.0;
      }
    }
    length = nh * nhhh;
    if(length != 0){
      X_ijkl2[i] = new double[length];
      X_ijkl3[i] = new double[length];
      V_ijkl[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_ijkl2[i][j] = 0.0;
	X_ijkl3[i][j] = 0.0;
	V_ijkl[i][j] = 0.0;
      }
    }
    length = nh * nhpp1;
    if(length != 0){
      X1_iajb3[i] = new double[length];
      X1_iajb4[i] = new double[length];
      X3_iajb3[i] = new double[length];
      X_iajb3[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_iajb3[i][j] = 0.0;
	X1_iajb4[i][j] = 0.0;
	X3_iajb3[i][j] = 0.0;
	X_iajb3[i][j] = 0.0;
      }
    }
    length = np * nhhp1;
    if(length != 0){
      X1_iajb2[i] = new double[length];
      X3_iajb2[i] = new double[length];
      X3_iajb5[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X1_iajb2[i][j] = 0.0;
	X3_iajb2[i][j] = 0.0;
	X3_iajb5[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      X_ijab1[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_ijab1[i][j] = 0.0;
      }
    }
    length = nhh * nhh;
    if(length != 0){
      X_ijkl1[i] = new double[length];
      X_ijkl4[i] = new double[length];
      Map_ijkl[i] = new int[8 * length];
      for(int j = 0; j < length; ++j){
	X_ijkl1[i][j] = 0.0;
	X_ijkl4[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_ijkl[i][8*j + k] = -1; }
      }
    }
    length = npp * npp;
    if(length != 0){
      X1_abcd1[i] = new double[length];
      X_abcd1[i] = new double[length];
      Map_abcd[i] = new int[6 * length];
      for(int j = 0; j < length; ++j){
	X1_abcd1[i][j] = 0.0;
	X_abcd1[i][j] = 0.0;
	for(int k = 0; k < 6; ++k){ Map_abcd[i][6*j + k] = -1; }
      }
    }
    length = npp * nhp;
    if(length != 0){
      X_iabc5[i] = new double[length];
      X_abic7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_iabc5[i][j] = 0.0;
	X_abic7[i][j] = 0.0;
      }
    }
    length = nhh * nhp;
    if(length != 0){
      X_ijka5[i] = new double[length];
      X2_iajk7[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_ijka5[i][j] = 0.0;
	X2_iajk7[i][j] = 0.0;
      }
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    length = npp1 * nhp1;
    if(length != 0){
      X_iabc4[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_iabc4[i][j] = 0.0;
      }
    }
    length = nhh1 * nhp1;
    if(length != 0){
      X_ijka4[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_ijka4[i][j] = 0.0;
      }
    }
    length = nhp2 * nhp2;
    if(length != 0){
      X1_iajb1[i] = new double[length];
      X3_iajb1[i] = new double[length];
      X_iajb1[i] = new double[length];
      Map_iajb[i] = new int[8 * length];
      for(int j = 0; j < length; ++j){
	X1_iajb1[i][j] = 0.0;
	X3_iajb1[i][j] = 0.0;
	X_iajb1[i][j] = 0.0;
	for(int k = 0; k < 8; ++k){ Map_iajb[i][8*j + k] = -1; }
      }
    }
    length = npp1 * nhp2;
    if(length != 0){
      X_abic5[i] = new double[length];
      X_abic6[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X_abic5[i][j] = 0.0;
	X_abic6[i][j] = 0.0;
      }
    }
    length = nhh1 * nhp2;
    if(length != 0){
      X2_iajk5[i] = new double[length];
      X2_iajk6[i] = new double[length];
      for(int j = 0; j < length; ++j){
	X2_iajk5[i][j] = 0.0;
	X2_iajk6[i][j] = 0.0;
      }
    }
  }
 
  length = Chan.nhp1[Chan.ind0];
  for(int i = 0; i < length; ++i){
    hind1 = Chan.hp1vec[Chan.ind0][2*i];
    pind1 = Chan.hp1vec[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    key1 = Chan.h_map[ind][hind1];
    key2 = Chan.p_map[ind][pind1];
    Map_ia[3*i] = ind;
    ind1 = key1 * Chan.np[ind] + key2;
    Map_ia[3*i + 1] = ind1;
    ind1 = key2 * Chan.nh[ind] + key1;
    Map_ia[3*i + 2] = ind1;
  }
  length = Chan.npp1[Chan.ind0];
  for(int i = 0; i < length; ++i){
    pind2 = Chan.pp1vec[Chan.ind0][2*i];
    pind1 = Chan.pp1vec[Chan.ind0][2*i + 1];
    ind = Chan.indvec[pind1];
    key1 = Chan.p_map[ind][pind1];
    key2 = Chan.p_map[ind][pind2];
    Map_ab[3*i] = ind;
    ind1 = key1 * Chan.np[ind] + key2;
    Map_ab[3*i + 1] = ind1;
    ind1 = key2 * Chan.np[ind] + key1;
    Map_ab[3*i + 2] = ind1;
  }
  length = Chan.nhh1[Chan.ind0];
  for(int i = 0; i < length; ++i){
    hind2 = Chan.hh1vec[Chan.ind0][2*i];
    hind1 = Chan.hh1vec[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    key1 = Chan.h_map[ind][hind1];
    key2 = Chan.h_map[ind][hind2];
    Map_ij[3*i] = ind;
    ind1 = key1 * Chan.nh[ind] + key2;
    Map_ij[3*i + 1] = ind1;
    ind1 = key2 * Chan.nh[ind] + key1;
    Map_ij[3*i + 2] = ind1;
  }
  length = Chan.nhp1[Chan.ind0];
  for(int i = 0; i < length; ++i){
    hind1 = Chan.hp1vec[Chan.ind0][2*i];
    pind1 = Chan.hp1vec[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    key1 = Chan.p_map[ind][pind1];
    key2 = Chan.h_map[ind][hind1];
    Map_ai[3*i] = ind;
    ind1 = key1 * Chan.nh[ind] + key2;
    Map_ai[3*i + 1] = ind1;
    ind1 = key2 * Chan.np[ind] + key1;
    Map_ai[3*i + 2] = ind1;
  }

  for(int i = 0; i < Chan.size3; ++i){
    length = Chan.np[i] * Chan.nhpp[i];
    length0 = Chan.nhpp[i];
    for(int phpp = 0; phpp < length; ++phpp){
      hpp1 = int(phpp%length0);
      p1 = int((phpp - hpp1)/length0);
      pind1 = Chan.pvec[i][p1];
      hind1 = Chan.hppvec[i][3*hpp1];
      pind2 = Chan.hppvec[i][3*hpp1 + 1];
      pind3 = Chan.hppvec[i][3*hpp1 + 2];
      ind = Chan.indvec[hind1];
      Map_iabc[i][8*phpp] = ind;
      key1 = Chan.ppp_map[ind][Hash3(pind1, pind2, pind3, Space.indtot)];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_iabc[i][8*phpp + 1] = ind1;
      ind = Chan.indvec[pind2];
      Map_iabc[i][8*phpp + 2] = ind;
      key1 = Chan.hpp1_map[ind][Hash3(hind1, pind1, pind3, Space.indtot)];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_iabc[i][8*phpp + 3] = ind1;
      minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind1].j - Space.qnums[pind2].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_iabc[i][8*phpp + 4] = ind;
	  key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	  key2 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp1[ind] + key2;
	  Map_iabc[i][8*phpp + 5] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iabc[i][8*phpp + 4] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	Map_iabc[i][8*phpp + 5] = ind1;
      }
      plus(tb, Space.qnums[pind2], Space.qnums[pind3]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[pind2].j - Space.qnums[pind3].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  Map_iabc[i][8*phpp + 6] = ind;
	  key1 = Chan.pp_map[ind][Hash2(pind2, pind3, Space.indtot)];
	  key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp[ind] + key2;
	  Map_iabc[i][8*phpp + 7] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iabc[i][8*phpp + 6] = ind;
	key1 = Chan.pp_map[ind][Hash2(pind2, pind3, Space.indtot)];
	key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp[ind] + key2;
	Map_iabc[i][8*phpp + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    length = Chan.nh[i] * Chan.nhhp[i];
    length0 = Chan.nhhp[i];
    for(int hhhp = 0; hhhp < length; ++hhhp){
      hhp1 = int(hhhp%length0);
      h1 = int((hhhp - hhp1)/length0);
      hind3 = Chan.hvec[i][h1];
      hind1 = Chan.hhpvec[i][3*hhp1];
      hind2 = Chan.hhpvec[i][3*hhp1 + 1];
      pind1 = Chan.hhpvec[i][3*hhp1 + 2];
      
      ind = Chan.indvec[pind1];
      Map_ijka[i][6*hhhp] = ind;
      key1 = Chan.hhh_map[ind][Hash3(hind3, hind1, hind2, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_ijka[i][6*hhhp + 1] = ind1;
      minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind2].j - Space.qnums[pind1].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_ijka[i][6*hhhp + 2] = ind;
	  key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	  key2 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp1[ind] + key2;
	  Map_ijka[i][6*hhhp + 3] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_ijka[i][6*hhhp + 2] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	key2 = Chan.hp1_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp1[ind] + key2;
	Map_ijka[i][6*hhhp + 3] = ind1;
      }
      plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind1].j - Space.qnums[hind2].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  Map_ijka[i][6*hhhp + 4] = ind;
	  key1 = Chan.hp_map[ind][Hash2(hind3, pind1, Space.indtot)];
	  key2 = Chan.hh_map[ind][Hash2(hind1, hind2, Space.indtot)];
	  ind1 = key1 * Chan.nhh[ind] + key2;
	  Map_ijka[i][6*hhhp + 5] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_ijka[i][6*hhhp + 4] = ind;
	key1 = Chan.hp_map[ind][Hash2(hind3, pind1, Space.indtot)];
	key2 = Chan.hh_map[ind][Hash2(hind1, hind2, Space.indtot)];
	ind1 = key1 * Chan.nhh[ind] + key2;
	Map_ijka[i][6*hhhp + 5] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    length = Chan.npp[i] * Chan.npp[i];
    length0 = Chan.npp[i];
    for(int pppp = 0; pppp < length; ++pppp){
      pp2 = int(pppp%length0);
      pp1 = int((pppp - pp2)/length0);
      pind3 = Chan.ppvec[i][2*pp1];
      pind4 = Chan.ppvec[i][2*pp1 + 1];
      pind1 = Chan.ppvec[i][2*pp2];
      pind2 = Chan.ppvec[i][2*pp2 + 1];
      
      ind = Chan.indvec[pind2];
      Map_abcd[i][6*pppp] = ind;
      key1 = Chan.ppp_map[ind][Hash3(pind1, pind3, pind4, Space.indtot)];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_abcd[i][6*pppp + 1] = ind1;
      ind = Chan.indvec[pind1];
      Map_abcd[i][6*pppp + 2] = ind;
      key1 = Chan.ppp_map[ind][Hash3(pind2, pind3, pind4, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_abcd[i][6*pppp + 3] = ind1;
      ind = Chan.indvec[pind3];
      Map_abcd[i][6*pppp + 4] = ind;
      key1 = Chan.ppp_map[ind][Hash3(pind4, pind1, pind2, Space.indtot)];
      key2 = Chan.p_map[ind][pind3];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_abcd[i][6*pppp + 5] = ind1;
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    length = Chan.nhh[i] * Chan.nhh[i];
    length0 = Chan.nhh[i];
    for(int hhhh = 0; hhhh < length; ++hhhh){
      hh2 = int(hhhh%length0);
      hh1 = int((hhhh - hh2)/length0);
      hind1 = Chan.hhvec[i][2*hh1];
      hind2 = Chan.hhvec[i][2*hh1 + 1];
      hind3 = Chan.hhvec[i][2*hh2];
      hind4 = Chan.hhvec[i][2*hh2 + 1];
      
      ind = Chan.indvec[hind4];
      Map_ijkl[i][8*hhhh] = ind;
      key1 = Chan.hhh_map[ind][Hash3(hind3, hind1, hind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind4];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_ijkl[i][8*hhhh + 1] = ind1;
      ind = Chan.indvec[hind3];
      Map_ijkl[i][8*hhhh + 2] = ind;
      key1 = Chan.hhh_map[ind][Hash3(hind4, hind1, hind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind3];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_ijkl[i][8*hhhh + 3] = ind1;
      ind = i;
      Map_ijkl[i][8*hhhh + 4] = ind;
      Map_ijkl[i][8*hhhh + 5] = Chan.nhh[i]*hh2 + hh1;
      ind = Chan.indvec[hind2];
      Map_ijkl[i][8*hhhh + 6] = ind;
      key1 = Chan.hhh_map[ind][Hash3(hind1, hind3, hind4, Space.indtot)];
      key2 = Chan.h_map[ind][hind2];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_ijkl[i][8*hhhh + 7] = ind1;
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    length = Chan.nhp2[i] * Chan.nhp2[i];
    length0 = Chan.nhp2[i];
    for(int hphp = 0; hphp < length; ++hphp){
      hp2 = int(hphp%length0);
      hp1 = int((hphp - hp2)/length0);
      hind1 = Chan.hp2vec[i][2*hp1];
      pind2 = Chan.hp2vec[i][2*hp1 + 1];
      hind2 = Chan.hp2vec[i][2*hp2];
      pind1 = Chan.hp2vec[i][2*hp2 + 1];
      
      ind = Chan.indvec[pind1];
      Map_iajb[i][8*hphp] = ind;
      key1 = Chan.hhp1_map[ind][Hash3(hind1, hind2, pind2, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_iajb[i][8*hphp + 1] = ind1;
      ind = Chan.indvec[hind2];
      Map_iajb[i][8*hphp + 2] = ind;
      key1 = Chan.hpp1_map[ind][Hash3(hind1, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind2];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_iajb[i][8*hphp + 3] = ind1;
      ind = Chan.indvec[hind1];
      Map_iajb[i][8*hphp + 4] = ind;
      key1 = Chan.hpp1_map[ind][Hash3(hind2, pind2, pind1, Space.indtot)];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_iajb[i][8*hphp + 5] = ind1;
      ind = Chan.indvec[pind2];
      Map_iajb[i][8*hphp + 6] = ind;
      key1 = Chan.hhp1_map[ind][Hash3(hind2, hind1, pind1, Space.indtot)];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_iajb[i][8*hphp + 7] = ind1;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    length = Chan.nhpp[i] * Chan.np[i];
    length0 = Chan.np[i];
    for(int hppp = 0; hppp < length; ++hppp){
      p1 = int(hppp%length0);
      hpp1 = int((hppp - p1)/length0);
      hind1 = Chan.hppvec[i][3*hpp1];
      pind1 = Chan.hppvec[i][3*hpp1 + 1];
      pind2 = Chan.hppvec[i][3*hpp1 + 2];
      pind3 = Chan.pvec[i][p1];
      
      ind = Chan.indvec[hind1];
      Map_abic[i][12*hppp] = ind;
      key1 = Chan.ppp_map[ind][Hash3(pind3, pind1, pind2, Space.indtot)];
      key2 = Chan.h_map[ind][hind1];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_abic[i][12*hppp + 1] = ind1;
      ind = Chan.indvec[pind1];
      Map_abic[i][12*hppp + 2] = ind;
      key1 = Chan.hpp1_map[ind][Hash3(hind1, pind3, pind2, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_abic[i][12*hppp + 3] = ind1;
      ind = Chan.indvec[pind2];
      Map_abic[i][12*hppp + 4] = ind;
      key1 = Chan.hpp1_map[ind][Hash3(hind1, pind3, pind1, Space.indtot)];
      key2 = Chan.p_map[ind][pind2];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_abic[i][12*hppp + 5] = ind1;
      minus(tb, Space.qnums[pind3], Space.qnums[pind2]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[pind3].j - Space.qnums[pind2].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_abic[i][12*hppp + 6] = ind;
	  key1 = Chan.pp1_map[ind][Hash2(pind3, pind2, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  Map_abic[i][12*hppp + 7] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 6] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind2, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Map_abic[i][12*hppp + 7] = ind1;
      }
      minus(tb, Space.qnums[pind3], Space.qnums[pind1]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[pind3].j - Space.qnums[pind1].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_abic[i][12*hppp + 8] = ind;
	  key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  Map_abic[i][12*hppp + 9] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 8] = ind;
	key1 = Chan.pp1_map[ind][Hash2(pind3, pind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Map_abic[i][12*hppp + 9] = ind1;
      }
      plus(tb, Space.qnums[pind1], Space.qnums[pind2]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[pind1].j - Space.qnums[pind2].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  Map_abic[i][12*hppp + 10] = ind;
	  key1 = Chan.hp_map[ind][Hash2(hind1, pind3, Space.indtot)];
	  key2 = Chan.pp_map[ind][Hash2(pind1, pind2, Space.indtot)];
	  ind1 = key1 * Chan.npp[ind] + key2;
	  Map_abic[i][12*hppp + 11] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_abic[i][12*hppp + 10] = ind;
	key1 = Chan.hp_map[ind][Hash2(hind1, pind3, Space.indtot)];
	key2 = Chan.pp_map[ind][Hash2(pind1, pind2, Space.indtot)];
	ind1 = key1 * Chan.npp[ind] + key2;
	Map_abic[i][12*hppp + 11] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    length = Chan.nhhp[i] * Chan.nh[i];
    length0 = Chan.nh[i];
    for(int hhph = 0; hhph < length; ++hhph){
      h1 = int(hhph%length0);
      hhp1 = int((hhph - h1)/length0);
      hind2 = Chan.hhpvec[i][3*hhp1];
      hind3 = Chan.hhpvec[i][3*hhp1 + 1];
      pind1 = Chan.hhpvec[i][3*hhp1 + 2];
      hind1 = Chan.hvec[i][h1];
      
      ind = Chan.indvec[pind1];
      Map_iajk[i][12*hhph] = ind;
      key1 = Chan.hhh_map[ind][Hash3(hind1, hind2, hind3, Space.indtot)];
      key2 = Chan.p_map[ind][pind1];
      ind1 = key1 * Chan.np[ind] + key2;
      Map_iajk[i][12*hhph + 1] = ind1;
      ind = Chan.indvec[hind3];
      Map_iajk[i][12*hhph + 2] = ind;
      key1 = Chan.hhp1_map[ind][Hash3(hind2, hind1, pind1, Space.indtot)];
      key2 = Chan.h_map[ind][hind3];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_iajk[i][12*hhph + 3] = ind1;
      ind = Chan.indvec[hind2];
      Map_iajk[i][12*hhph + 4] = ind;
      key1 = Chan.hhp1_map[ind][Hash3(hind3, hind1, pind1, Space.indtot)];
      key2 = Chan.h_map[ind][hind2];
      ind1 = key1 * Chan.nh[ind] + key2;
      Map_iajk[i][12*hhph + 5] = ind1;
      minus(tb, Space.qnums[hind2], Space.qnums[hind1]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind2].j - Space.qnums[hind1].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_iajk[i][12*hhph + 6] = ind;
	  key1 = Chan.hh1_map[ind][Hash2(hind2, hind1, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind3, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  Map_iajk[i][12*hhph + 7] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 6] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind2, hind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind3, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Map_iajk[i][12*hhph + 7] = ind1;
      }
      minus(tb, Space.qnums[hind3], Space.qnums[hind1]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind3].j - Space.qnums[hind1].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  Map_iajk[i][12*hhph + 8] = ind;
	  key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp2[ind] + key2;
	  Map_iajk[i][12*hhph + 9] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 8] = ind;
	key1 = Chan.hh1_map[ind][Hash2(hind3, hind1, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(hind2, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp2[ind] + key2;
	Map_iajk[i][12*hhph + 9] = ind1;
      }
      plus(tb, Space.qnums[hind2], Space.qnums[hind3]);
      if(Parameters.basis == "finite_J"){
	jmin = abs(Space.qnums[hind2].j - Space.qnums[hind3].j);
	while(tb.j >= jmin){
	  ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  Map_iajk[i][12*hhph + 10] = ind;
	  key1 = Chan.hh_map[ind][Hash2(hind2, hind3, Space.indtot)];
	  key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	  ind1 = key1 * Chan.nhp[ind] + key2;
	  Map_iajk[i][12*hhph + 11] = ind1;
	  tb.j -= 2;
	}
      }
      else{
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iajk[i][12*hhph + 10] = ind;
	key1 = Chan.hh_map[ind][Hash2(hind2, hind3, Space.indtot)];
	key2 = Chan.hp_map[ind][Hash2(hind1, pind1, Space.indtot)];
	ind1 = key1 * Chan.nhp[ind] + key2;
	Map_iajk[i][12*hhph + 11] = ind1;
      }
    }
  }
}

void CC_Eff::delete_struct(const Channels &Chan)
{
  int length, nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhp2, nhh1, npp1, nhpp1, nhhp1, nhhh, nppp;

  nhp1 = Chan.nhp1[Chan.ind0];
  npp1 = Chan.npp1[Chan.ind0];
  nhh1 = Chan.nhh1[Chan.ind0];
  nhp1 = Chan.nhp1[Chan.ind0];
  if(nhp1 != 0){
    delete[] X_ia1;
    delete[] Map_ia;
  }
  if(npp1 != 0){
    delete[] X_ab1;
    delete[] Map_ab;
  }
  if(nhh1 != 0){
    delete[] X_ij1;
    delete[] X1_ij1;
    delete[] Map_ij;
  }
  if(nhp1 != 0){
    delete[] X_ai1;
    delete[] Map_ai;
  }

  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    nhhh = Chan.nhhh[i];
    nppp = Chan.nppp[i];
    length = nh * np;
    if(length != 0){
      delete[] X_ia2[i]; // (i,a)
      delete[] X_ia3[i]; // (a,i)
      delete[] X_ai2[i]; // (a,i)
      delete[] X_ai3[i]; // (i,a)
    }
    length = np * np;
    if(length != 0){
      delete[] X_ab2[i]; // (a,b)
      delete[] X_ab3[i]; // (b,a)
    }
    length = nh * nh;
    if(length != 0){
      delete[] X_ij2[i]; // (i,j)
      delete[] X_ij3[i]; // (j,i)
      delete[] X1_ij2[i]; // (i,j)
      delete[] X1_ij3[i]; // (j,i)
    }
    length = np * nhpp;
    if(length != 0){
      delete[] X1_iabc1[i];
      delete[] X1_iabc3[i];
      delete[] X_iabc1[i];
      delete[] Map_iabc[i];
      delete[] X_abic1[i];
      delete[] Map_abic[i];
    }
    length = nh * nhhp;
    if(length != 0){
      delete[] X1_ijka1[i];
      delete[] X_ijka1[i];
      delete[] Map_ijka[i];
      delete[] X2_iajk1[i];
      delete[] X_iajk1[i];
      delete[] Map_iajk[i];
    }
    length = nh * nppp;
    if(length != 0){
      delete[] X1_iabc2[i];
      delete[] X_abic2[i];
    }
    length = np * nhhh;
    if(length != 0){
      delete[] X1_ijka2[i];
      delete[] X2_iajk2[i];
    }
    length = np * nhpp1;
    if(length != 0){
      delete[] X_iabc3[i];      
      delete[] X1_iajb2[i];
      delete[] X_abic3[i];
      delete[] X_abic4[i];
    }
    length = nh * nhhp1;
    if(length != 0){
      delete[] X2_iajk3[i];
      delete[] X2_iajk4[i];
    }
    length = np * nppp;
    if(length != 0){
      delete[] X1_abcd2[i];
      delete[] X1_abcd3[i];
      delete[] V_abcd[i];
    }
    length = nh * nhhh;
    if(length != 0){
      delete[] X_ijkl2[i];
      delete[] X_ijkl3[i];
      delete[] V_ijkl[i];
    }
    length = nh * nhpp1;
    if(length != 0){
      delete[] X1_iajb3[i];
      delete[] X1_iajb4[i];
      delete[] X3_iajb3[i];
      delete[] X_iajb3[i];
    }
    length = np * nhhp1;
    if(length != 0){
      delete[] X3_iajb2[i];
      delete[] X3_iajb5[i];
    }
  }
  delete[] X_ia2;
  delete[] X_ia3;
  delete[] X_ab2;
  delete[] X_ab3;
  delete[] X_ij2;
  delete[] X_ij3;
  delete[] X1_ij2;
  delete[] X1_ij3;
  delete[] X_ai2;
  delete[] X_ai3;
  delete[] X1_iabc1;
  delete[] X1_iabc2;
  delete[] X1_iabc3;
  delete[] X_iabc1;
  delete[] X_iabc3;
  delete[] Map_iabc;
  delete[] X1_ijka1;
  delete[] X1_ijka2;
  delete[] X_ijka1;
  delete[] Map_ijka;
  delete[] X1_abcd2;
  delete[] X1_abcd3;
  delete[] V_abcd;
  delete[] X_ijkl2;
  delete[] X_ijkl3;
  delete[] V_ijkl;
  delete[] X1_iajb2;
  delete[] X1_iajb3;
  delete[] X1_iajb4;
  delete[] X3_iajb2;
  delete[] X3_iajb3;
  delete[] X3_iajb5;
  delete[] X_iajb3;
  delete[] X_abic1;
  delete[] X_abic2;
  delete[] X_abic3;
  delete[] X_abic4;
  delete[] Map_abic;
  delete[] X2_iajk1;
  delete[] X2_iajk2;
  delete[] X2_iajk3;
  delete[] X2_iajk4;
  delete[] X_iajk1;
  delete[] Map_iajk;

  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    length = nhh * npp;
    if(length != 0){
      delete[] X_ijab1[i];
    }
    length = nhh * nhh;
    if(length != 0){
      delete[] X_ijkl1[i];
      delete[] X_ijkl4[i];
      delete[] Map_ijkl[i];
    }
    length = npp * npp;
    if(length != 0){
      delete[] X1_abcd1[i];
      delete[] X_abcd1[i];
      delete[] Map_abcd[i];
    }
    length = npp * nhp;
    if(length != 0){
      delete[] X_iabc5[i];
      delete[] X_abic7[i];
    }
    length = nhh * nhp;
    if(length != 0){
      delete[] X_ijka5[i];
      delete[] X2_iajk7[i];
    }
  }
  delete[] X_ijab1;
  delete[] X_iabc5;
  delete[] X_ijka5;
  delete[] X1_abcd1;
  delete[] X_abcd1;
  delete[] Map_abcd;
  delete[] X_ijkl1;
  delete[] X_ijkl4;
  delete[] Map_ijkl;
  delete[] X_abic7;
  delete[] X2_iajk7;

  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    length = npp1 * nhp1;
    if(length != 0){
      delete[] X_iabc4[i];
    }
    length = nhh1 * nhp1;
    if(length != 0){
      delete[] X_ijka4[i];
    }
    length = nhp2 * nhp2;
    if(length != 0){
      delete[] X1_iajb1[i];
      delete[] X3_iajb1[i];
      delete[] X_iajb1[i];
      delete[] Map_iajb[i];
    }
    length = npp1 * nhp2;
    if(length != 0){
      delete[] X_abic5[i];
      delete[] X_abic6[i];
    }
    length = nhh1 * nhp2;
    if(length != 0){
      delete[] X2_iajk5[i];
      delete[] X2_iajk6[i];
    }
  }
  delete[] X_iabc4;
  delete[] X_ijka4;
  delete[] X1_iajb1;
  delete[] X3_iajb1;
  delete[] X_iajb1;
  delete[] Map_iajb;
  delete[] X_abic5;
  delete[] X_abic6;
  delete[] X2_iajk5;
  delete[] X2_iajk6;
}

void CC_Eff::set_X_ia(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhp1[Chan.ind0]; ++i){
    x = X_ia1[i];
    X_ia2[Map_ia[3*i]][Map_ia[3*i + 1]] = x;
    X_ia3[Map_ia[3*i]][Map_ia[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ab(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.npp1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ab1[i];
    x += X_ab2[Map_ab[3*i]][Map_ab[3*i + 1]];
    x += X_ab3[Map_ab[3*i]][Map_ab[3*i + 2]];
    X_ab1[i] = x;
    X_ab2[Map_ab[3*i]][Map_ab[3*i + 1]] = x;
    X_ab3[Map_ab[3*i]][Map_ab[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ij(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ij1[i];
    x += X_ij2[Map_ij[3*i]][Map_ij[3*i + 1]];
    x += X_ij3[Map_ij[3*i]][Map_ij[3*i + 2]];
    X_ij1[i] = x;
    X_ij2[Map_ij[3*i]][Map_ij[3*i + 1]] = x;
    X_ij3[Map_ij[3*i]][Map_ij[3*i + 2]] = x;
  }
}

void CC_Eff::set_X1_ij(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhh1[Chan.ind0]; ++i){
    x = 0.0;
    x += X1_ij1[i];
    x += X1_ij2[Map_ij[3*i]][Map_ij[3*i + 1]];
    x += X1_ij3[Map_ij[3*i]][Map_ij[3*i + 2]];
    X1_ij1[i] = x;
    X1_ij2[Map_ij[3*i]][Map_ij[3*i + 1]] = x;
    X1_ij3[Map_ij[3*i]][Map_ij[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ai(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.nhp1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ai1[i];
    x += X_ai2[Map_ai[3*i]][Map_ai[3*i + 1]];
    x += X_ai3[Map_ai[3*i]][Map_ai[3*i + 2]];
    X_ai1[i] = x;
    X_ai2[Map_ai[3*i]][Map_ai[3*i + 1]] = x;
    X_ai3[Map_ai[3*i]][Map_ai[3*i + 2]] = x;
  }
}

void CC_Eff::set_X1_iabc(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int p1 = 0; p1 < Chan.np[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	ind0 = p1*Chan.nhpp[i] + hpp1;
	X1_iabc2[Map_iabc[i][8*ind0]][Map_iabc[i][8*ind0 + 1]] = X1_iabc1[i][ind0];
	X1_iabc3[Map_iabc[i][8*ind0 + 2]][Map_iabc[i][8*ind0 + 3]] = X1_iabc1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X_iabc(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int p1 = 0; p1 < Chan.np[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
	ind0 = p1*Chan.nhpp[i] + hpp1;
	X_iabc3[Map_iabc[i][8*ind0 + 2]][Map_iabc[i][8*ind0 + 3]] = X_iabc1[i][ind0];
	X_iabc4[Map_iabc[i][8*ind0 + 4]][Map_iabc[i][8*ind0 + 5]] = X_iabc1[i][ind0];
	X_iabc5[Map_iabc[i][8*ind0 + 6]][Map_iabc[i][8*ind0 + 7]] = X_iabc1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X1_ijka(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	ind0 = h1*Chan.nhhp[i] + hhp1;
	X1_ijka2[Map_ijka[i][6*ind0]][Map_ijka[i][6*ind0 + 1]] = X1_ijka1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X_ijka(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
	ind0 = h1*Chan.nhhp[i] + hhp1;
	X_ijka4[Map_ijka[i][6*ind0 + 2]][Map_ijka[i][6*ind0 + 3]] = X_ijka1[i][ind0];
	X_ijka5[Map_ijka[i][6*ind0 + 4]][Map_ijka[i][6*ind0 + 5]] = X_ijka1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X1_abcd(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size1; ++i){
    for(int pp1 = 0; pp1 < Chan.npp[i]; ++pp1){
      for(int pp2 = 0; pp2 < Chan.npp[i]; ++pp2){
	x = 0.0;
	ind0 = Chan.npp[i]*pp1 + pp2;
	V_abcd[Map_abcd[i][6*ind0 + 4]][Map_abcd[i][6*ind0 + 5]] = X1_abcd1[i][ind0];
	x += X1_abcd1[i][ind0];
	x += X1_abcd2[Map_abcd[i][6*ind0]][Map_abcd[i][6*ind0 + 1]];
	x += X1_abcd3[Map_abcd[i][6*ind0 + 2]][Map_abcd[i][6*ind0 + 3]];
	X1_abcd1[i][ind0] = x;
	X1_abcd2[Map_abcd[i][6*ind0]][Map_abcd[i][6*ind0 + 1]] = x;
	X1_abcd3[Map_abcd[i][6*ind0 + 2]][Map_abcd[i][6*ind0 + 3]] = x;
      }
    }
  }
}

void CC_Eff::set_X_ijkl(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size1; ++i){
    for(int hh1 = 0; hh1 < Chan.nhh[i]; ++hh1){
      for(int hh2 = 0; hh2 < Chan.nhh[i]; ++hh2){
	x = 0.0;
	ind0 = Chan.nhh[i]*hh1 + hh2;
	V_ijkl[Map_ijkl[i][8*ind0 + 6]][Map_ijkl[i][8*ind0 + 7]] = X_ijkl1[i][ind0];
	x += X_ijkl1[i][ind0];
	x += X_ijkl2[Map_ijkl[i][8*ind0]][Map_ijkl[i][8*ind0 + 1]];
	x += X_ijkl3[Map_ijkl[i][8*ind0 + 2]][Map_ijkl[i][8*ind0 + 3]];
	x += X_ijkl4[Map_ijkl[i][8*ind0 + 4]][Map_ijkl[i][8*ind0 + 5]];
	X_ijkl1[i][ind0] = x;
	X_ijkl2[Map_ijkl[i][8*ind0]][Map_ijkl[i][8*ind0 + 1]] = x;
	X_ijkl3[Map_ijkl[i][8*ind0 + 2]][Map_ijkl[i][8*ind0 + 3]] = x;
	X_ijkl4[Map_ijkl[i][8*ind0 + 4]][Map_ijkl[i][8*ind0 + 5]] = x;
      }
    }
  }
}

void CC_Eff::set_X1_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X1_iajb1[i][ind0];
	x += X1_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]];
	x += X1_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X1_iajb1[i][ind0] = x;
	X1_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]] = x;
	X1_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
	X1_iajb4[Map_iajb[i][8*ind0 + 4]][Map_iajb[i][8*ind0 + 5]] = x;
      }
    }
  }
}

void CC_Eff::set_X3_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X3_iajb1[i][ind0];
	x += X3_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]];
	x += X3_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X3_iajb1[i][ind0] = x;
	X3_iajb2[Map_iajb[i][8*ind0]][Map_iajb[i][8*ind0 + 1]] = x;
	X3_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
	X3_iajb5[Map_iajb[i][8*ind0 + 6]][Map_iajb[i][8*ind0 + 7]] = x;
      }
    }
  }
}

void CC_Eff::set_X_iajb(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.nhp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.nhp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.nhp2[i]*hp1 + hp2;
	x += X_iajb1[i][ind0];
	x += X_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]];
	X_iajb1[i][ind0] = x;
	X_iajb3[Map_iajb[i][8*ind0 + 2]][Map_iajb[i][8*ind0 + 3]] = x;
      }
    }
  }
}

void CC_Eff::set_X_abic(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size3; ++i){
    for(int hpp1 = 0; hpp1 < Chan.nhpp[i]; ++hpp1){
      for(int p1 = 0; p1 < Chan.np[i]; ++p1){
	x = 0.0;
	ind0 = Chan.np[i]*hpp1 + p1;
	x += X_abic1[i][ind0];
	x += X_abic2[Map_abic[i][12*ind0]][Map_abic[i][12*ind0 + 1]];
	x += X_abic3[Map_abic[i][12*ind0 + 2]][Map_abic[i][12*ind0 + 3]];
	x += X_abic4[Map_abic[i][12*ind0 + 4]][Map_abic[i][12*ind0 + 5]];
	x += X_abic5[Map_abic[i][12*ind0 + 6]][Map_abic[i][12*ind0 + 7]];
	x += X_abic6[Map_abic[i][12*ind0 + 8]][Map_abic[i][12*ind0 + 9]];
	x += X_abic7[Map_abic[i][12*ind0 + 10]][Map_abic[i][12*ind0 + 11]];
	X_abic1[i][ind0] = x;
	X_abic2[Map_abic[i][12*ind0]][Map_abic[i][12*ind0 + 1]] = x;
	X_abic3[Map_abic[i][12*ind0 + 2]][Map_abic[i][12*ind0 + 3]] = x;
	X_abic4[Map_abic[i][12*ind0 + 4]][Map_abic[i][12*ind0 + 5]] = x;
	X_abic5[Map_abic[i][12*ind0 + 6]][Map_abic[i][12*ind0 + 7]] = x;
	X_abic6[Map_abic[i][12*ind0 + 8]][Map_abic[i][12*ind0 + 9]] = x;
	X_abic7[Map_abic[i][12*ind0 + 10]][Map_abic[i][12*ind0 + 11]] = x;
      }
    }
  }
}

void CC_Eff::set_X2_iajk(const Channels &Chan)
{
  int ind0;
  double x;
  for(int i = 0; i < Chan.size3; ++i){
    for(int hhp1 = 0; hhp1 < Chan.nhhp[i]; ++hhp1){
      for(int h1 = 0; h1 < Chan.nh[i]; ++h1){
	x = 0.0;
	ind0 = Chan.nh[i]*hhp1 + h1;
	x += X2_iajk1[i][ind0];
	x += X2_iajk2[Map_iajk[i][12*ind0]][Map_iajk[i][12*ind0 + 1]];
	x += X2_iajk3[Map_iajk[i][12*ind0 + 2]][Map_iajk[i][12*ind0 + 3]];
	x += X2_iajk4[Map_iajk[i][12*ind0 + 4]][Map_iajk[i][12*ind0 + 5]];
	x += X2_iajk5[Map_iajk[i][12*ind0 + 6]][Map_iajk[i][12*ind0 + 7]];
	x += X2_iajk6[Map_iajk[i][12*ind0 + 8]][Map_iajk[i][12*ind0 + 9]];
	x += X2_iajk7[Map_iajk[i][12*ind0 + 10]][Map_iajk[i][12*ind0 + 11]];
	X2_iajk1[i][ind0] = x;
	X2_iajk2[Map_iajk[i][12*ind0]][Map_iajk[i][12*ind0 + 1]] = x;
	X2_iajk3[Map_iajk[i][12*ind0 + 2]][Map_iajk[i][12*ind0 + 3]] = x;
	X2_iajk4[Map_iajk[i][12*ind0 + 4]][Map_iajk[i][12*ind0 + 5]] = x;
	X2_iajk5[Map_iajk[i][12*ind0 + 6]][Map_iajk[i][12*ind0 + 7]] = x;
	X2_iajk6[Map_iajk[i][12*ind0 + 8]][Map_iajk[i][12*ind0 + 9]] = x;
	X2_iajk7[Map_iajk[i][12*ind0 + 10]][Map_iajk[i][12*ind0 + 11]] = x;
      }
    }
  }
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

void Print_Parameters(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << std::endl;
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << "Case = " << Parameters.calc_case << ", Basis = " << Parameters.basis << ", Approximation = " << Parameters.approx << std::endl;
  if(Parameters.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << Parameters.LevelScheme << std::endl;
    if(Parameters.MatrixElements.size() > 0){ std::cout << "Interaction = " << Parameters.MatrixElements << std::endl; }
  }
  if(Parameters.calc_case == "nuclear"){
    if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){ std::cout << "Total States = " << Space.indtot << std::endl; }
    else{ std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl; }
    std::cout << "Proton Shells = " << Parameters.Pshells << ", Neutron Shells = " << Parameters.Nshells << std::endl;
    std::cout << "Protons = " << Parameters.P << ", Neutrons = " << Parameters.N << std::endl;
    if(Parameters.calc_case == "infinite"){ std::cout << "Density = " << Parameters.density << std::endl; }
  }
  else if(Parameters.calc_case == "electronic"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Density = " << Parameters.density << std::endl;
  }
  else if(Parameters.calc_case == "quantum_dot"){
    std::cout << "Number of Shells = " << Parameters.Shells << ", Total States = " << Space.indtot << std::endl;
    std::cout << "Electron Shells = " << Parameters.Pshells << ", Electrons = " << Parameters.P << std::endl;
    std::cout << "Oscillator Energy = " << Parameters.density << std::endl;
  }
  std::cout << "---------------------------------------------------------------------------" << std::endl;
  std::cout << std::endl;
}

void Model_Space::delete_struct(Input_Parameters &Parameters)
{
  delete[] qnums;
  if(Parameters.basis == "infinite"){ delete[] map_2b; }
  if(Parameters.basis == "finite_JM" && indtotj != 0){
    for(int i = 0; i < indtotj; ++i){ delete[] shellsm[i]; }
    delete[] shellsm;
  }
}

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  std::string fullpath; // Model Space file path
  std::string phline; // Std::String for each file line
  std::ifstream splevels; // Model space file

  int ind; // state index
  int l; // orbital angular momentum quantum number
  double energy; // state energy

  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count
  
  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  getline(splevels, phline);
  std::stringstream(phline) >> Space.indtot;

  // allocate memory for quntum numbers for each state
  Space.qnums = new State[Space.indtot];

  // initialize mins and maxs for each quantum number
  Space.qmins.t = 1000;
  Space.qmins.m = 1000;
  Space.qmins.nx = 1000;
  Space.qmins.ny = 1000;
  Space.qmins.nz = 1000;
  Space.qmins.ml = 1000;
  Space.qmins.j = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.t = -1000;
  Space.qmaxs.m = -1000;
  Space.qmaxs.nx = -1000;
  Space.qmaxs.ny = -1000;
  Space.qmaxs.nz = -1000;
  Space.qmaxs.ml = -1000;
  Space.qmaxs.j = -1000;
  Space.qmaxs.par = -1000;

  // States must be ordered by energy from low to high
  for(int i = 0; i < Space.indtot; ++i){
    getline(splevels, phline);
    if(Parameters.basis == "infinite"){ // ind, nx, ny, nz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].nx >> Space.qnums[i].ny >> Space.qnums[i].nz >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].ml = 0;
      Space.qnums[i].par = 1;
      Space.qnums[i].j = 0;
      if(Space.qnums[i].nx < Space.qmins.nx){ Space.qmins.nx = Space.qnums[i].nx; }
      if(Space.qnums[i].nx > Space.qmaxs.nx){ Space.qmaxs.nx = Space.qnums[i].nx; }
      if(Space.qnums[i].ny < Space.qmins.ny){ Space.qmins.ny = Space.qnums[i].ny; }
      if(Space.qnums[i].ny > Space.qmaxs.ny){ Space.qmaxs.ny = Space.qnums[i].ny; }
      if(Space.qnums[i].nz < Space.qmins.nz){ Space.qmins.nz = Space.qnums[i].nz; }
      if(Space.qnums[i].nz > Space.qmaxs.nz){ Space.qmaxs.nz = Space.qnums[i].nz; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_M"){ // ind, n, l, j, jz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].j >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      Space.qnums[i].ml = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_HO"){ // ind, n, l, lz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> Space.qnums[i].ml >> Space.qnums[i].m >> energy;
      std::cout << "! " << ind << " " << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << "  " << energy << std::endl;
      Space.qnums[i].par = 1;
      Space.qnums[i].t = -1;
      Space.qnums[i].j = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].ml < Space.qmins.ml){ Space.qmins.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].ml > Space.qmaxs.ml){ Space.qmaxs.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){ // ind, n, l, j, tz, l2n
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].j >> Space.qnums[i].t >> ind >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      Space.qnums[i].m = 0;
      Space.qnums[i].ml = 0;
      Space.qnums[i].nx = 0;
      Space.qnums[i].ny = 0;
      Space.qnums[i].nz = 0;
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].j < Space.qmins.j){ Space.qmins.j = Space.qnums[i].j; }
      if(Space.qnums[i].j > Space.qmaxs.j){ Space.qmaxs.j = Space.qnums[i].j; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    Space.qnums[i].energy = energy * Parameters.obstrength;
  }
  splevels.close();

  // Determine shell structure and assign hole/particle labels
  int *pshell;
  int *nshell;
  int pshell_num = 0;
  int nshell_num = 0;
  double p_en = -1000.0;
  double n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }
  pshell = new int[pshell_num];
  nshell = new int[nshell_num];
  pshell_num = 0;
  nshell_num = 0;
  p_en = -1000.0;
  n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      nshell[nshell_num] = i;
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }

  if(Parameters.Pshells > pshell_num){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > nshell_num){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      ++pcount;
      if(i < pshell[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      ++ncount;
      if(i < nshell[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
  }
  delete[] pshell;
  delete[] nshell;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // Find range of 2body quantum numbers from mins and maxs
  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.j = Space.qmaxs.j + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;

  // Get total number of 2body states (and setup index array for infinite case)
  if(Parameters.basis == "infinite"){
    int N2max = 0;
    int N2maxtemp;
    for(int i = 0; i < Space.indtot; ++i){
      for(int j = i+1; j < Space.indtot; ++j){
	N2maxtemp = (Space.qnums[i].nx + Space.qnums[j].nx)*(Space.qnums[i].nx + Space.qnums[j].nx) + 
	  (Space.qnums[i].ny + Space.qnums[j].ny)*(Space.qnums[i].ny + Space.qnums[j].ny) + 
	  (Space.qnums[i].nz + Space.qnums[j].nz)*(Space.qnums[i].nz + Space.qnums[j].nz);
	if(N2maxtemp > N2max){ N2max = N2maxtemp; }
      }
    }
    Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
    int count1 = 0;
    for(int nx = 2*Space.qmins.nx; nx <= 2*Space.qmaxs.nx; ++nx){
      for(int ny = 2*Space.qmins.ny; ny <= 2*Space.qmaxs.ny; ++ny){
	for(int nz = 2*Space.qmins.nz; nz <= 2*Space.qmaxs.nz; ++nz){
	  if(nx*nx + ny*ny + nz*nz <= N2max){
	    Space.map_2b[(nx - 2*Space.qmins.nx)*(Space.qsizes.ny*Space.qsizes.nz) +
			 (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = count1;
	    ++count1;
	  }
	  else{
	    Space.map_2b[(nx - 2*Space.qmins.nx)*(Space.qsizes.ny*Space.qsizes.nz) +
			 (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = -1;
	  }
	}
      }
    }
    Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;
  }
  else if(Parameters.basis == "finite_M"){
    Space.size_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_HO"){
    Space.size_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
    Space.size_2b = Space.qsizes.par * Space.qsizes.j * Space.qsizes.t;
  }
}

void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space)
{
  int ind; // state index
  int m; // angular momenum projection
  int shelllength; // number of single particle states for each shell

  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count

  State *states = new State[Space.indtot];
  int indtotj = Space.indtot;
  Space.indtotj = indtotj;
  Space.shellsm = new int*[Space.indtot];
  // reset Space.indtot with degeneracies 2j + 1
  Space.indtot = 0;
  for(int i = 0; i < indtotj; ++i){
    Space.indtot += Space.qnums[i].j + 1;
    states[i] = Space.qnums[i];
  }

  // allocate memory for quntum numbers for each state
  delete[] Space.qnums;
  Space.qnums = new State[Space.indtot];
  Space.qmins.t = 1000;
  Space.qmins.m = 1000;
  Space.qmins.j = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.t = -1000;
  Space.qmaxs.m = -1000;
  Space.qmaxs.j = -1000;
  Space.qmaxs.par = -1000;

  ind = 0;
  for(int i = 0; i < indtotj; ++i){
    shelllength = states[i].j + 1;
    Space.shellsm[i] = new int[shelllength];
    for(int j = 0; j < shelllength; ++j){
      m = -states[i].j + 2*j;
      Space.qnums[ind].par = states[i].par;
      Space.qnums[ind].j = states[i].j;
      Space.qnums[ind].m = m;
      Space.qnums[ind].t = states[i].t;
      Space.qnums[ind].energy = states[i].energy;
      Space.qnums[ind].type = states[i].type;
      Space.shellsm[i][j] = ind;
      if(Space.qnums[ind].par < Space.qmins.par){ Space.qmins.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].m < Space.qmins.m){ Space.qmins.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].t < Space.qmins.t){ Space.qmins.t = Space.qnums[ind].t; }
      if(Space.qnums[ind].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[ind].t; }
      ++ind;
    }
  }
  delete[] states;
  
  // Determine shell structure and assign hole/particle labels
  int *pshell;
  int *nshell;
  int pshell_num = 0;
  int nshell_num = 0;
  double p_en = -1000.0;
  double n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }
  pshell = new int[pshell_num];
  nshell = new int[nshell_num];
  pshell_num = 0;
  nshell_num = 0;
  p_en = -1000.0;
  n_en = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > p_en && Space.qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = Space.qnums[i].energy;
      ++pshell_num;
    }
    if(Space.qnums[i].energy > n_en && Space.qnums[i].t == 1){
      nshell[nshell_num] = i;
      n_en = Space.qnums[i].energy;
      ++nshell_num;
    }
  }

  if(Parameters.Pshells > pshell_num){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > nshell_num){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      ++pcount;
      if(i < pshell[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      ++ncount;
      if(i < nshell[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++parcount; }
    }
  }
  delete[] pshell;
  delete[] nshell;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // Find range of 2body quantum numbers from mins and maxs
  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.j = Space.qmaxs.j + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;
  Space.size_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
}

void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  double electron_prefac = hbarc_HartA * hbarc_HartA / (2.0 * m_electronc2_Hart);
  double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
  double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);

  int count = 0; // total state count
  int pcount = 0; // proton state count
  int ncount = 0; // neutron state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  int nhcount = 0; // neutron hole state count

  int shellnums [] = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23,
  		      25, 26, 27, 28, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46};

  // Find appropriate Nmax for the number of shells using shellnums
  if(Parameters.Shells > 40){ std::cerr << "Nmax too big!" << std::endl; exit(1); }
  Space.Nmax = shellnums[Parameters.Shells - 1];

  if(Parameters.calc_case == "electronic"){
    Parameters.Nshells = 0;
    double r_b = hbarc_HartA / (m_electronc2_Hart * fine_struct);
    Parameters.density = 3.0/(4.0 * M_PI * std::pow(Parameters.density * r_b, 3.0));
  }

  // Find maximum number of states (indtot) depending on Nmax and whether or not there are protons/neutrons
  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = 1, Space.qsizes.t = 3; }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1; }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){ Space.qmins.t = 1, Space.qmaxs.t = 1, Space.qsizes.t = 1; }
  else if(Parameters.calc_case == "nuclear"){ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }
  else if(Parameters.calc_case == "electronic"){ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  // Find total number of states
  Space.nmax = std::floor(std::sqrt(Space.Nmax));
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){
	for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){
	  if(shell != nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1 && Parameters.Pshells == 0){ continue; }
	      if(tz == 1 && Parameters.Nshells == 0){ continue; }
	      count++;
	    }
	  }
	}
      }
    }
  }
  Space.indtot = count;
  count = 0;

  // allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){    
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){	
	for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1){
		if(Parameters.Pshells == 0){ continue; }
		++pcount;
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
		if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++holcount; ++phcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		++ncount;
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
		if(shell < Parameters.Nshells){ Space.qnums[count].type = "hole"; ++holcount; ++nhcount; }
		else{ Space.qnums[count].type = "particle"; ++parcount; }
	      }
	      Space.qnums[count].nx = nx;
	      Space.qnums[count].ny = ny;
	      Space.qnums[count].nz = nz;
	      Space.qnums[count].m = sz;
	      Space.qnums[count].t = tz;
	      Space.qnums[count].ml = 0;
	      Space.qnums[count].par = 1;
	      Space.qnums[count].j = 0;
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.qsizes.m = 3; // -2, 0, 2
  Space.qsizes.nx = 4*Space.nmax + 1;
  Space.qsizes.ny = 4*Space.nmax + 1;
  Space.qsizes.nz = 4*Space.nmax + 1;
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  // With the number of hole states, the length scale and state energies can be calculated
  double L = pow(Space.indhol/Parameters.density, 1.0/3.0);
  if(Parameters.calc_case == "nuclear"){
    for(int i = 0; i < Space.indtot; ++i){
      if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac * M_PI*M_PI / (L*L); }
      else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac * M_PI*M_PI / (L*L); }
    }
    // Change energies to Hartree-Fock energies, E_p = E_p + V_pipi
    for(int p = 0; p < Space.indtot; ++p){
      for(int i = 0; i < Space.indhol; ++i){
	if(p == i){ continue; }
	Space.qnums[p].energy += vint_Minnesota_Momentum(Space, p, i, p, i, L);
      }
    }
  }
  else if(Parameters.calc_case == "electronic"){
    for(int i = 0; i < Space.indtot; ++i){ Space.qnums[i].energy *= electron_prefac * M_PI*M_PI / (L*L); }
    // Change energies to Hartree-Fock energies, E_p = E_p + V_pipi
    for(int p = 0; p < Space.indtot; ++p){
      for(int i = 0; i < Space.indhol; ++i){
	if(p == i){ continue; }
	Space.qnums[p].energy += Coulomb_Inf(Space, p, i, p, i, L);;
      }
    }
  }

  // Order two-body states with respect to Nx, Ny, Nz for channel index function
  int count1 = 0;
  int N2max = 4*Space.Nmax;
  Space.map_2b = new int[Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz];
  for(int nx = -2 * Space.nmax; nx <= 2 * Space.nmax; ++nx){
    for(int ny = -2 * Space.nmax; ny <= 2 * Space.nmax; ++ny){
      for(int nz = -2 * Space.nmax; nz <= 2 * Space.nmax; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max){
	  Space.map_2b[(nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (ny + 2*Space.nmax) * Space.qsizes.nz + (nz + 2*Space.nmax)] = count1;
	  ++count1;
	}
	else{
	  Space.map_2b[(nx + 2*Space.nmax) * Space.qsizes.ny*Space.qsizes.nz + (ny + 2*Space.nmax) * Space.qsizes.nz + (nz + 2*Space.nmax)] = -1;
	}
      }
    }
  }
  Space.size_2b = count1 * Space.qsizes.t * Space.qsizes.m;
}

void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  int count = 0; // total state count
  int pcount = 0; // proton state count
  int holcount = 0; // hole state count
  int parcount = 0; // particle state count
  int phcount = 0; // proton hole state count
  Parameters.Nshells = 0;
  Space.qmins.ml = 1000;
  Space.qmins.par = 1000;
  Space.qmaxs.ml = -1000;
  Space.qmaxs.par = -1000;
  Space.nmax = -1000;

  if(Parameters.Pshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1, Space.qsizes0.t = 1; }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  // Find total number of states
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int ml = -1 * shell; ml <= shell; ++ml){
      for(int n = 0; n < Parameters.Shells; ++n){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  if(n > Space.nmax){ Space.nmax = n; }
	  if(ml < Space.qmins.ml){ Space.qmins.ml = ml; }
	  if(ml > Space.qmaxs.ml){ Space.qmaxs.ml = ml; }
	  count++;
	}
      }
    }
  }
  Space.indtot = count;
  count = 0;

  // Allocate memory for quantum numbers for each state
  Space.qnums = new State[Space.indtot];
  Space.nmax = 0;
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int ml = -1 * shell; ml <= shell; ++ml){
      for(int n = 0; n < Parameters.Shells; ++n){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  E = Parameters.density*(shell + 1.0);
	  Space.qnums[count].energy = E;
	  Space.qnums[count].n = n;
	  Space.qnums[count].par = -2*(abs(ml)%2) + 1;
	  Space.qnums[count].ml = ml;
	  Space.qnums[count].m = sz;
	  Space.qnums[count].t = -1;
	  Space.qnums[count].nx = 0;
	  Space.qnums[count].ny = 0;
	  Space.qnums[count].nz = 0;
	  Space.qnums[count].j = 0;
	  if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++holcount; ++phcount; }
	  else{ Space.qnums[count].type = "particle"; ++parcount; }
	  count++;
	}
      }
    }
  }
  Space.qsizes.m = 3; // -2, 0, +2
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.par = 2; // -1, +1
  Space.size_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;
  Space.qsizes0.m = 2; // -1, +1
  Space.qsizes0.ml = Space.qmaxs.ml - Space.qmins.ml + 1;

  int key;
  for(int i = 0; i < Space.indtot; ++i){
    key = ChanInd_1b(Parameters.basis, Space, Space.qnums[i]);
    Space.map_1b[key] = i;
  }

  Space.indp = pcount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = 0;
}

 // !!!!!!!!!!!!!!!!!!!!!          FIX THIS!!!!!!!!!!!!!!!!!!!!!
int ChanInd_1b(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return (State.nx - Space.qmins.nx)*Space.qsizes0.ny*Space.qsizes0.nz*2*Space.qsizes0.t + (State.ny - Space.qmins.ny)*Space.qsizes0.nz*2*Space.qsizes0.t + (State.nz - Space.qmins.nz)*2*Space.qsizes0.t + int(0.5 * (State.m - Space.qmins.m))*Space.qsizes0.t + int(0.5 * (State.t - Space.qmins.t));
  }
  else if(basis == "finite_HO"){ // for tz = -1 only;
    return State.n*Space.qsizes.ml*2 + (State.ml - Space.qmins.ml)*2 + int(0.5*(State.m + 1));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx + 2*Space.nmax)*Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax)*Space.qsizes.nz +
			(State.nz + 2*Space.nmax)]*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t + 
                         int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - 2*Space.qmins.m)/2)*Space.qsizes.t + 
           int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - 2*Space.qmins.ml)*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t +
           int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - 2*Space.qmins.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx + 2*Space.nmax)*Space.qsizes.ny*Space.qsizes.nz + (State.ny + 2*Space.nmax)*Space.qsizes.nz +
			(State.nz + 2*Space.nmax)]*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t + 
                         int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - Space.qmins.m+Space.qmaxs.m)/2)*Space.qsizes.t + int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - Space.qmins.ml+Space.qmaxs.ml)*Space.qsizes.m*Space.qsizes.t + int((State.m + 2)/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double L = pow(Space.indhol/Parameters.density, 1.0/3.0);
  int i, j, k, l, a, b, c, d;
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(npp == 0){ goto stop1; }
    #pragma omp parallel private(a, b, c, d, TBME)
    {
      #pragma omp for schedule(static)
      for(int cd = 0; cd < npp; ++cd){
	c = Chan.ppvec[chan][2*cd];
	d = Chan.ppvec[chan][2*cd + 1];
	if(c == d){ continue; }
	for(int ab = cd; ab < npp; ++ab){
	  a = Chan.ppvec[chan][2*ab];
	  b = Chan.ppvec[chan][2*ab + 1];
	  if(a == b){ continue; }
	  TBME = Coulomb_Inf(Space, a, b, c, d, L);
	  Ints.D_ME1.V1[chan][npp*cd + ab] = TBME;
	  Ints.D_ME1.V1[chan][npp*ab + cd] = TBME;
	}
      }
    }
  stop1:;
    if(nhh == 0){ goto stop2; }
    #pragma omp parallel private(i, j, k, l, TBME)
    {
      #pragma omp for schedule(static)
      for(int ij = 0; ij < nhh; ++ij){
	i = Chan.hhvec[chan][2*ij];
	j = Chan.hhvec[chan][2*ij + 1];
	if(i == j){ continue; }
	for(int kl = ij; kl < nhh; ++kl){
	  k = Chan.hhvec[chan][2*kl];
	  l = Chan.hhvec[chan][2*kl + 1];
	  if(k == l){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, k, l, L);
	  Ints.D_ME1.V2[chan][nhh*ij + kl] = TBME;
	  Ints.D_ME1.V2[chan][nhh*kl + ij] = TBME;
	}
      }
    }
  stop2:;
    if(nhh * npp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int ab = 0; ab < npp; ++ab){
	a = Chan.ppvec[chan][2*ab];
	b = Chan.ppvec[chan][2*ab + 1];
	if(a == b){ continue; }
	for(int ij = 0; ij < nhh; ++ij){
	  i = Chan.hhvec[chan][2*ij];
	  j = Chan.hhvec[chan][2*ij + 1];
	  if(i == j){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V4[chan][nhh*ab + ij] = TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size3; ++chan){
    nh = Chan.nh[chan];
    np = Chan.np[chan];
    nhpp = Chan.nhpp[chan];
    nhhp = Chan.nhhp[chan];
    if(nh * nhpp == 0){ goto stop3; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int h = 0; h < nh; ++h){
	j = Chan.hvec[chan][h];
	for(int hpp = 0; hpp < nhpp; ++hpp){
	  i = Chan.hppvec[chan][3*hpp];
	  a = Chan.hppvec[chan][3*hpp + 1];
	  b = Chan.hppvec[chan][3*hpp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V5[chan][nhpp*h + hpp] = TBME;
	  Ints.D_ME1.V6[chan][nhpp*h + hpp] = -1.0 * TBME;
	}
      }
    }
  stop3:;
    if(np * nhhp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int p = 0; p < np; ++p){
	b = Chan.pvec[chan][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  i = Chan.hhpvec[chan][3*hhp];
	  j = Chan.hhpvec[chan][3*hhp + 1];
	  a = Chan.hhpvec[chan][3*hhp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V7[chan][nhhp*p + hhp] = TBME;
	  Ints.D_ME1.V8[chan][nhhp*p + hhp] = -1.0 * TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size2; ++chan){
    nhp1 = Chan.nhp1[chan];
    nhp2 = Chan.nhp2[chan];
    if(nhp2 == 0){ goto stop4; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp21 = 0; hp21 < nhp2; ++hp21){
	i = Chan.hp2vec[chan][2*hp21];
	b = Chan.hp2vec[chan][2*hp21 + 1];
	for(int hp22 = hp21; hp22 < nhp2; ++hp22){
	  j = Chan.hp2vec[chan][2*hp22];
	  a = Chan.hp2vec[chan][2*hp22 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, a, j, b, L);
	  Ints.D_ME1.V3[chan][nhp2*hp21 + hp22] = TBME;
	  Ints.D_ME1.V3[chan][nhp2*hp22 + hp21] = TBME;
	}
      }
    }
  stop4:;
    if(nhp1 * nhp2 == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp2 = 0; hp2 < nhp2; ++hp2){
	i = Chan.hp2vec[chan][2*hp2];
	a = Chan.hp2vec[chan][2*hp2 + 1];
	for(int hp1 = 0; hp1 < Chan.nhp1[chan]; ++hp1){
	  j = Chan.hp1vec[chan][2*hp1];
	  b = Chan.hp1vec[chan][2*hp1 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V9[chan][nhp1*hp2 + hp1] = TBME;
	  Ints.D_ME1.V10[chan][nhp1*hp2 + hp1] = -1.0 * TBME;
	}
      }
    }
  }
}

double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
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
    if(qSquared1 < 1.0e-9){ term += 0.0; }
    else{ term += 4.0 * prefactor * M_PI / qSquared1; }
  }
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-9){ term -= 0.0; }
    else{ term -= 4.0 * prefactor * M_PI / qSquared1; }
  }
  return term;
}


// Minnesota Potential for momentum basis
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql)
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

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  for(int i = 0; i < Chan.size1; ++i){
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.npp[i]; ++rs){
	shell3 = Chan.ppvec[i][2*rs];
	shell4 = Chan.ppvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V1[i][Chan.npp[i]*pq + rs] = TBME;
	Ints.D_ME1.V1[i][Chan.npp[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.nhh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec[i][2*pq];
      shell2 = Chan.hhvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V2[i][Chan.nhh[i]*pq + rs] = TBME;
	Ints.D_ME1.V2[i][Chan.nhh[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = 0; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V4[i][Chan.nhh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    for(int q = 0; q < Chan.nh[i]; ++q){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell2 = Chan.hvec[i][q];
      for(int prs = 0; prs < Chan.nhpp[i]; ++prs){
	shell1 = Chan.hppvec[i][3*prs];
	shell3 = Chan.hppvec[i][3*prs + 1];
	shell4 = Chan.hppvec[i][3*prs + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V5[i][Chan.nhpp[i]*q + prs] = TBME;
	Ints.D_ME1.V6[i][Chan.nhpp[i]*q + prs] = -1.0 * TBME;
      }
    }
    for(int s = 0; s < Chan.np[i]; ++s){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell4 = Chan.pvec[i][s];
      for(int pqr = 0; pqr < Chan.nhhp[i]; ++pqr){
	shell1 = Chan.hhpvec[i][3*pqr];
	shell2 = Chan.hhpvec[i][3*pqr + 1];
	shell3 = Chan.hhpvec[i][3*pqr + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V7[i][Chan.nhhp[i]*s + pqr] = TBME;
	Ints.D_ME1.V8[i][Chan.nhhp[i]*s + pqr] = -1.0 * TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size2; ++i){
    for(int ps = 0; ps < Chan.nhp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*ps];
      shell4 = Chan.hp2vec[i][2*ps + 1];
      for(int qr = ps; qr < Chan.nhp2[i]; ++qr){
	shell3 = Chan.hp2vec[i][2*qr];
	shell2 = Chan.hp2vec[i][2*qr + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V3[i][Chan.nhp2[i]*ps + qr] = TBME;
	Ints.D_ME1.V3[i][Chan.nhp2[i]*qr + ps] = TBME;
      }
    }
    for(int pr = 0; pr < Chan.nhp2[i]; ++pr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*pr];
      shell3 = Chan.hp2vec[i][2*pr + 1];
      for(int qs = 0; qs < Chan.nhp1[i]; ++qs){
	shell2 = Chan.hp1vec[i][2*qs];
	shell4 = Chan.hp1vec[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell4, shell3, L);
	Ints.D_ME1.V9[i][Chan.nhp1[i]*pr + qs] = TBME;
	Ints.D_ME1.V10[i][Chan.nhp1[i]*pr + qs] = -1.0 * TBME;
      }
    }
  }
}

int kron_del(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
}
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}

// Minnesota Potential for momentum basis
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
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

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  double mixmin = 0.1;
  double mixmax = 1.0;
  double mix = 1.0;
  int nhp1, nhhpp;
  int goodcount = 0;
  int goodint = 5;
  int badcount = 0;
  double width = 1.0;
  int iterations = 10;
  int Denomend = 0;
  bool denom = true;
  std::string search = "search";
  std::string stay = "stay";
  double error = 1e20, error2 = 1e20;
  int ind = 0;
  double CCoutE;
  double tempt;
  Amplitudes Amps0 = Amplitudes(Parameters, Space, Chan);
  Amplitudes Amps2 = Amplitudes(Parameters, Space, Chan);
  Amplitudes tempAmps = Amplitudes(Parameters, Space, Chan);
  //// for DIIS //////
  /*double norm, norm1, maxnorm;
  int maxind, ortho;
  double checkdot;
  double checknorm1, checknorm2;
  int nhh, npp;
  int DIISstart = 5;
  int maxl = 10;
  int N = 0;
  int P = N + 1;
  int size = Chan.size1;
  int i1, j1, a1, b1;
  if(Parameters.approx == "singles"){ size += 1; }
  double ***p = new double**[maxl];
  double ***delp = new double**[maxl];
  for(int l = 0; l < maxl; ++l){
    p[l] = new double*[size];
    delp[l] = new double*[size];
    for(int chan = 0; chan < Chan.size1; ++chan){
      p[l][chan] = new double[Chan.nhh[chan] * Chan.npp[chan]];
      delp[l][chan] = new double[Chan.nhh[chan] * Chan.npp[chan]];
      for(int hhpp = 0; hhpp < Chan.nhh[chan] * Chan.npp[chan]; ++hhpp){
	p[l][chan][hhpp] = 0.0;
	delp[l][chan][hhpp] = 0.0;
      }
    }
    if(Parameters.approx == "singles"){
      p[l][Chan.size1] = new double[Chan.nhp1[Chan.ind0]];
      delp[l][Chan.size1] = new double[Chan.nhp1[Chan.ind0]];
      for(int hp = 0; hp < Chan.nhp1[Chan.ind0]; ++hp){
	p[l][Chan.size1][hp] = 0.0;
	delp[l][Chan.size1][hp] = 0.0;
      }	
    }
  }
  double **tempdelp = new double*[size];
  for(int chan = 0; chan < Chan.size1; ++chan){
    tempdelp[chan] = new double[Chan.nhh[chan] * Chan.npp[chan]];
    for(int hhpp = 0; hhpp < Chan.nhh[chan] * Chan.npp[chan]; ++hhpp){
      tempdelp[chan][hhpp] = 0.0;
    }
  }
  if(Parameters.approx == "singles"){
    tempdelp[Chan.size1] = new double[Chan.nhp1[Chan.ind0]];
    for(int hp = 0; hp < Chan.nhp1[Chan.ind0]; ++hp){
      tempdelp[Chan.size1][hp] = 0.0;
    }	
  }
  double *B = new double[P * P];
  double *B2 = new double[P * P];
  int lwork = sizeof(double) * P;
  int *ipiv = new int[P];
  double *work = new double[sizeof(double) * P];
  int info = 0;
  B[0] = 0.0;
  B2[0] = 0.0;*/
  ////////////////////

  // Initialize Amplitudes //
  Gather_Amps0(Parameters, Space, Chan, Ints, Amps, tempAmps, mix);
  Amps0.copy_Amplitudes(Parameters, Chan, Amps);
  Amps2.copy_Amplitudes(Parameters, Chan, Amps);
  CCoutE = Amps.get_energy(Parameters, Chan, Ints);
  if(Parameters.approx == "singles"){ std::cout << "Iteration Number = " << ind << ", CCSD Energy = " << CCoutE << std::endl; }
  else{ std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCoutE << std::endl; }

  mix = mixmin;
  while((error > 1e-12 && ind < 40000) || ind < 50){
    if( ind < Denomend ){ denom = true; }
    else{ denom = false; }

    tempAmps.zero(Parameters, Chan);
    if(Parameters.basis != "finite_J"){ Doubles_Step(Space, Chan, Ints, Amps, tempAmps); }
    else{ Doubles_Step_J(Space, Chan, Ints, Amps, tempAmps); }
    if(Parameters.approx == "singles"){
      if(Parameters.basis != "finite_J"){
	Doubles_Step_2(Space, Chan, Ints, Amps, tempAmps);
	Singles_Step(Space, Chan, Ints, Amps, tempAmps);
      }
      else{
	Doubles_Step_2_J(Space, Chan, Ints, Amps, tempAmps);
	Singles_Step_J(Space, Chan, Ints, Amps, tempAmps);
      }
    }
    Amps2.zero1(Parameters, Chan);
    Gather_Amps(Parameters, Space, Chan, Ints, Amps2, tempAmps, mix, denom);
    CC_Error(Parameters, Chan, Ints, Amps, Amps2, error);
    error /= mix;

    if(error < error2 || ind < 25){
      --badcount;
      Amps0.copy_Amplitudes(Parameters, Chan, Amps);
      Amps.copy_Amplitudes(Parameters, Chan, Amps2);
      error2 = error;
    }
    else{
      ++badcount;
      mix = mixmin;
      if(error2 > 1.0){ width = 0.005; }
      else{ width = 0.005 * error2; }
      if(width < 1e-6){ width = 1e-6; }
      Random_Step(Parameters, Space, Chan, Ints, Amps0, Amps, Amps2, tempAmps, mix, denom, width, error2);
    }

    //// for DIIS //////
    /*norm1 = 0.0;
    norm = 0.0;
    for(int i = 0; i < Chan.size1; ++i){
      nhhpp = Chan.nhh[i] * Chan.npp[i];
      for(int j = 0; j < nhhpp; ++j){
	tempdelp[i][j] = Amps.D1.T1[i][j] - Amps0.D1.T1[i][j];
	norm1 += Amps.D1.T1[i][j] * Amps.D1.T1[i][j];
      }
    }
    if(Parameters.approx == "singles"){
      nhp1 = Chan.nhp1[Chan.ind0];
      for(int j = 0; j < nhp1; ++j){
	tempdelp[Chan.size1][j] = Amps.S1.T1[j] - Amps0.S1.T1[j];
	norm1 += Amps.S1.T1[j] * Amps.S1.T1[j];
      }	
    }
    norm1 = std::sqrt(norm1);
    for(int i = 0; i < Chan.size1; ++i){
      nhhpp = Chan.nhh[i] * Chan.npp[i];
      for(int j = 0; j < nhhpp; ++j){
	tempdelp[i][j] /= norm1;
	norm += tempdelp[i][j] * tempdelp[i][j];
      }
    }
    if(Parameters.approx == "singles"){
      nhp1 = Chan.nhp1[Chan.ind0];
      for(int j = 0; j < nhp1; ++j){
	tempdelp[Chan.size1][j] /= norm1;
	norm += tempdelp[Chan.size1][j] * tempdelp[Chan.size1][j];
      }	
    }

    // check orthogonality of tempdelp
    ortho = 1;
    for(int l = 0; l < N; ++l){
      checkdot = 0.0;
      checknorm1 = 0.0;
      checknorm2 = 0.0;
      for(int i = 0; i < Chan.size1; ++i){
	nhhpp = Chan.nhh[i] * Chan.npp[i];
	for(int j = 0; j < nhhpp; ++j){
	  checkdot += (tempdelp[i][j] * delp[l][i][j]);
	  checknorm1 += (tempdelp[i][j] * tempdelp[i][j]);
	  checknorm2 += (delp[l][i][j] * delp[l][i][j]);
	}
      }
      if(Parameters.approx == "singles"){
	nhp1 = Chan.nhp1[Chan.ind0];
	for(int j = 0; j < nhp1; ++j){
	  checkdot += (tempdelp[Chan.size1][j] * delp[l][Chan.size1][j]);
	  checknorm1 += (tempdelp[Chan.size1][j] * tempdelp[Chan.size1][j]);
	  checknorm2 += (delp[l][Chan.size1][j] * delp[l][Chan.size1][j]);
	}
      }
      checkdot = fabs(checkdot/(std::sqrt(checknorm1) * std::sqrt(checknorm2)));
      if(checkdot > 0.95 || norm > B[P * l + l]){ ortho = 0; break; }
    }
    if(ind < DIISstart){ ortho = 0; }

    if(ortho == 1){
      if(N < maxl){
	for(int i = 0; i < Chan.size1; ++i){
	  nhhpp = Chan.nhh[i] * Chan.npp[i];
	  for(int j = 0; j < nhhpp; ++j){
	    p[N][i][j] = Amps.D1.T1[i][j];
	    delp[N][i][j] = tempdelp[i][j];
	  }
	}
	if(Parameters.approx == "singles"){
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int j = 0; j < nhp1; ++j){
	    p[N][Chan.size1][j] = Amps.S1.T1[j];
	    delp[N][Chan.size1][j] = tempdelp[Chan.size1][j];
	  }
	}
	delete[] B2;
	B2 = new double[(P+1) * (P+1)];
	for(int j = 0; j < N; ++j){
	  for(int k = 0; k < N; ++k){
	    B2[(P+1) * j + k] = B[P * j + k];
	  }
	}
	for(int l = 0; l < P; ++l){
	  B2[(P+1) * N + l] = 0.0;
	  if(l != N){ B2[(P+1) * l + N] = 0.0; }
	  for(int i = 0; i < Chan.size1; ++i){
	    nhhpp = Chan.nhh[i] * Chan.npp[i];
	    for(int j = 0; j < nhhpp; ++j){
	      B2[(P+1) * N + l] += delp[N][i][j] * delp[l][i][j];
	      if(l != N){ B2[(P+1) * l + N] += delp[l][i][j] * delp[N][i][j]; }
	    }
	  }
	  if(Parameters.approx == "singles"){
	    nhp1 = Chan.nhp1[Chan.ind0];
	    for(int j = 0; j < nhp1; ++j){
	      B2[(P+1) * N + l] += delp[N][Chan.size1][j] * delp[l][Chan.size1][j];
	      if(l != N){ B2[(P+1) * l + N] += delp[l][Chan.size1][j] * delp[N][Chan.size1][j]; }
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
	for(int j = 0; j < P; ++j){
	  for(int k = 0; k < P; ++k){
	    B[P * j + k] = B2[P * j + k];
	  }
	}
      }
      else{
	for(int j = 0; j < P; ++j){
	  for(int k = 0; k < P; ++k){
	    B2[P * j + k] = B[P * j + k];
	  }
	}
	maxind = -1;
	maxnorm = 0.0;
	for(int j = 0; j < N; ++j){
	  norm = 0.0;
	  for(int k = 0; k < N; ++k){
	    norm += fabs(B2[P * j + k]);
	  }
	  if(norm > maxnorm){
	    maxind = j;
	    maxnorm = norm;
	  }
	}
	for(int i = 0; i < Chan.size1; ++i){
	  nhhpp = Chan.nhh[i] * Chan.npp[i];
	  for(int j = 0; j < nhhpp; ++j){
	    p[maxind][i][j] = Amps.D1.T1[i][j];
	    delp[maxind][i][j] = tempdelp[i][j];
	  }
	}
	if(Parameters.approx == "singles"){
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int j = 0; j < nhp1; ++j){
	    p[maxind][Chan.size1][j] = Amps.S1.T1[j];
	    delp[maxind][Chan.size1][j] = tempdelp[Chan.size1][j];
	  }
	}
	for(int l = 0; l < N; ++l){
	  B2[P * maxind + l] = 0.0;
	  if(l != maxind){ B2[P * l + maxind] = 0.0; }

	  for(int i = 0; i < Chan.size1; ++i){
	    nhhpp = Chan.nhh[i] * Chan.npp[i];
	    for(int j = 0; j < nhhpp; ++j){
	      B2[P * maxind + l] += delp[maxind][i][j] * delp[l][i][j];
	      if(l != maxind){ B2[P * l + maxind] += delp[l][i][j] * delp[maxind][i][j]; }
	    }
	  }
	  if(Parameters.approx == "singles"){
	    nhp1 = Chan.nhp1[Chan.ind0];
	    for(int j = 0; j < nhp1; ++j){
	      B2[P * maxind + l] += delp[maxind][Chan.size1][j] * delp[l][Chan.size1][j];
	      if(l != maxind){ B2[P * l + maxind] += delp[l][Chan.size1][j] * delp[maxind][Chan.size1][j]; }
	    }
	  }
	}
	for(int l = 0; l < N; ++l){
	  B2[P * N + l] = -1.0;
	  B2[P * l + N] = -1.0;
	}
	B2[(P+1) * P + P] = 0.0;
      }
      for(int j = 0; j < P; ++j){
	for(int k = 0; k < P; ++k){
	  B[P * j + k] = B2[P * j + k];
	}
      }
      delete[] ipiv;
      delete[] work;
      ipiv = new int[P];
      work = new double[sizeof(double) * P];
      lwork = sizeof(double) * P;
      info = 0;
      dgetrf_(&P, &P, B2, &P, ipiv, &info);
      dgetri_(&P, B2, &P, ipiv, work, &lwork, &info);
      for(int chan = 0; chan < Chan.size1; ++chan){
	nhh = Chan.npp[chan];
	npp = Chan.nhh[chan];
	if(nhh * npp == 0){ continue; }
	for(int hh = 0; hh < nhh; ++hh){
	  i1 = Chan.hhvec[chan][2*hh];
	  j1 = Chan.hhvec[chan][2*hh + 1];
	  for(int pp = 0; pp < npp; ++pp){
	    a1 = Chan.ppvec[chan][2*pp];
	    b1 = Chan.ppvec[chan][2*pp + 1];
	    tempt = 0.0;
	    ind = hh * npp + pp;
	    for(int l = 0; l < N; ++l){ tempt += -1.0 * B2[P * l + N] * p[l][chan][ind]; }
	    if(Parameters.basis != "finite_J"){ Amps.D1.set_T(chan, ind, tempt); }
	    else{ Amps.D1.set_TJ(Space, Chan, chan, ind, i1, j1, a1, b1, tempt); }
	  }
	}
      }
      if(Parameters.approx == "singles"){
	nhp1 = Chan.nhp1[Chan.ind0];
	for(int hp1 = 0; hp1 < nhp1; ++hp1){
	  i1 = Chan.hp1vec[Chan.ind0][2*hp1];
	  a1 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	  tempt = 0.0;
	  for(int l = 0; l < N; ++l){ tempt += -1.0 * B2[P * l + N] * p[l][Chan.size1][hp1]; }
	  if(Parameters.basis != "finite_J"){ Amps.S1.set_T(hp1, tempt); }
	  else{ Amps.S1.set_TJ(Space, hp1, i1, a1, tempt); }
	}
	if(Parameters.basis != "finite_J"){
	  Amps.S1.set_T_2(Chan, Ints);
	  Amps.D1.set_T_2(Chan, Ints);
	}
	else{
	  Amps.S1.set_T_2J(Space, Chan, Ints);
	  Amps.D1.set_T_2J(Chan, Ints);
	}
      }
      }*/
    ////////////////////

    CCoutE = Amps.get_energy(Parameters, Chan, Ints);
    if( !std::isfinite(CCoutE) ){ std::cerr << std::endl << ind << " : CC Solution Diverged!!" << std::endl; exit(1); }
    if(Parameters.approx == "singles"){ std::cout << "Iteration Number = " << ind+1 << ", CCSD Energy = " << CCoutE << ", error = " << error << "\r"; }
    else{ std::cout << "Iteration Number = " << ind+1 << ", CCD Energy = " << CCoutE << ", error = " << error << "\r"; }
    ++ind;
  }
  std::cout << std::endl << std::endl;
  if( error > 1e-12 ){
    std::cout << Parameters.Shells << ", " << Parameters.Pshells << ", " << Parameters.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". CC Solution Not Converged!!" << std::endl;
  }

  ///////For DIIS///////////////////////////////////
  /*delete[] ipiv;
  delete[] work;
  for(int l = 0; l < maxl; ++l){
    for(int chan = 0; chan < Chan.size1; ++chan){
      delete[] p[l][chan];
      delete[] delp[l][chan];
    }
    if(Parameters.approx == "singles"){
      delete[] p[l][Chan.size1];
      delete[] delp[l][Chan.size1];
    }
    delete[] p[l];
    delete[] delp[l];
  }
  delete[] p;
  delete[] delp;
  delete[] B;
  delete[] B2;*/
  /////////////////////////////////////////////////

  Amps0.delete_struct(Parameters, Chan);
  Amps2.delete_struct(Parameters, Chan);
  tempAmps.delete_struct(Parameters, Chan);
}

void Random_Step(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, Amplitudes &Amps2, Amplitudes &tempAmps, double &mix, bool &denom, double &width, double &error2)
{
  double error;
  int go = 1;
  int ind = 0;
  while(go == 1){
    ++ind;
    Amps.zero(Parameters, Chan);
    tempAmps.zero(Parameters, Chan);
    Randomize_Amps(Parameters, Space, Chan, Ints, Amps0, Amps, width);
    if(Parameters.basis == "finite_J"){ Doubles_Step_J(Space, Chan, Ints, Amps, tempAmps); }
    else{ Doubles_Step(Space, Chan, Ints, Amps, tempAmps); }
    if(Parameters.approx == "singles"){
      if(Parameters.basis == "finite_J"){
	Doubles_Step_2_J(Space, Chan, Ints, Amps, tempAmps);
	Singles_Step_J(Space, Chan, Ints, Amps, tempAmps);
      }
      else{
	Doubles_Step_2(Space, Chan, Ints, Amps, tempAmps);
	Singles_Step(Space, Chan, Ints, Amps, tempAmps);
      }
    }
    Amps2.zero1(Parameters, Chan);
    Gather_Amps(Parameters, Space, Chan, Ints, Amps2, tempAmps, mix, denom);
    CC_Error(Parameters, Chan, Ints, Amps, Amps2, error);
    error /= mix;

    if(error < error2){
      go = 0;
      Amps0.copy_Amplitudes(Parameters, Chan, Amps);
      Amps.copy_Amplitudes(Parameters, Chan, Amps2);
      error2 = error;
    }
    width *= 0.975;
    if(width < 10-10){ width = 10e-10; }
    if(ind >= 100){
      std::cout << Parameters.Shells << ", " << Parameters.Pshells << ", " << Parameters.density << std::endl;
      std::cerr << std::endl << "Random Step Unsuccessful!!" << std::endl; exit(1);
    }
  }
}

void Gather_Amps0(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix)
{
  int nhh, npp, a, b, i, j, ind;
  double tempt, tempen;
  double lowen = 1000;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	tempen = Space.qnums[i].energy + Space.qnums[j].energy - Space.qnums[a].energy - Space.qnums[b].energy;
	if(fabs(tempen) < lowen){ lowen = fabs(tempen); }
	Amps.D1.Evec[chan][ind] = tempen;
	Amps2.D1.Evec[chan][ind] = tempen;
	tempt = mix * Ints.D_ME1.V4[chan][pp * nhh + hh] / tempen;
	//std::cout << std::setprecision(12) << "T: " << i << " " << j << " | " << a << " " << b << " = " << tempt << " : " << tempen << std::endl;
	if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i, j, a, b, tempt); }
	else{ Amps.D1.set_T(chan, ind, tempt); }
      }
    }
  }
  if(Parameters.approx == "singles"){
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempen = Space.qnums[i].energy - Space.qnums[a].energy;
      if(fabs(tempen) < lowen){ lowen = fabs(tempen); }
      tempt = 0.0;
      Amps.S1.Evec[hp1] = tempen;
      Amps2.S1.Evec[hp1] = tempen;
      //std::cout << std::setprecision(12) << "t: " << i << " | " << a << " = " << tempt << " : " << tempen << std::endl;
      if(Parameters.basis == "finite_J"){ Amps.S1.set_TJ(Space, hp1, i, a, tempt); }
      else{ Amps.S1.set_T(hp1, tempt); }
    }
    if(Parameters.basis == "finite_J"){
      Amps.S1.set_T_2J(Space, Chan, Ints);
      Amps.D1.set_T_2J(Chan, Ints);
    }
    else{
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
  }

  if(lowen < 1.0){
    for(int chan = 0; chan < Chan.size1; ++chan){
      nhh = Chan.nhh[chan];
      npp = Chan.npp[chan];
      if(nhh * npp == 0){ continue; }
      for(int pp = 0; pp < npp; ++pp){
	a = Chan.ppvec[chan][2*pp];
	b = Chan.ppvec[chan][2*pp + 1];
	for(int hh = 0; hh < nhh; ++hh){
	  i = Chan.hhvec[chan][2*hh];
	  j = Chan.hhvec[chan][2*hh + 1];
	  ind = hh * npp + pp;
	  Amps2.D1.Evec[chan][ind] /= lowen;
	  tempt = Amps.D1.T1[chan][ind] * lowen;
	  if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i, j, a, b, tempt); }
	  else{ Amps.D1.set_T(chan, ind, tempt); }
	}
      }
    }
    if(Parameters.approx == "singles"){
      int nhp1 = Chan.nhp1[Chan.ind0];
      for(int hp1 = 0; hp1 < nhp1; ++hp1){
	Amps2.S1.Evec[hp1] /= lowen;
      }
    }
  }
}

void Gather_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix, bool &denom)
{
  int nhh, npp, a, b, i, j, ind;
  double tempt;
  double t1, t2;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	if(Parameters.basis == "finite_J"){ tempt = Amps2.D1.get_TJ(Space, Chan, chan, ind, i, j, a, b); }
	else{ tempt = Amps2.D1.get_T(chan, ind); }
	tempt += Ints.D_ME1.V4[chan][pp * nhh + hh];
	if( denom ){ tempt /= Amps2.D1.Evec[chan][ind]; }
	else{ tempt /= Amps.D1.Evec[chan][ind]; }
	tempt = mix*tempt + (1.0-mix)*Amps.D1.T1[chan][ind];
	//std::cout << std::setprecision(12) << "T: " << i << " " << j << " | " << a << " " << b << " = " << tempt << std::endl;
	if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i, j, a, b, tempt); }
	else{ Amps.D1.set_T(chan, ind, tempt); }
      }
    }
  }
  if(Parameters.approx == "singles"){  
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      if(Parameters.basis == "finite_J"){ tempt = Amps2.S1.get_TJ(Space, hp1, i, a); }
      else{ tempt = Amps2.S1.get_T(hp1); }
      if( denom ){ tempt /= Amps2.S1.Evec[hp1]; }
      else{ tempt /= Amps.S1.Evec[hp1]; }
      tempt = mix*tempt + (1.0-mix)*Amps.S1.T1[hp1];
      //std::cout << std::setprecision(12) << "t: " << i << " | " << a << " = " << tempt << std::endl;
      if(Parameters.basis == "finite_J"){ Amps.S1.set_TJ(Space, hp1, i, a, tempt); }
      else{ Amps.S1.set_T(hp1, tempt); }
    }
    if(Parameters.basis == "finite_J"){
      Amps.S1.set_T_2J(Space, Chan, Ints);
      Amps.D1.set_T_2J(Chan, Ints);
    }
    else{
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
  }
}

void Gather_Amps2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix, bool &denom, double &checkdot, int &N, double ***delp)
{
  int nhh, npp, a, b, i, j, ind;
  double tempt;
  double totalnorm1 = 0.0;
  double *totalnorm2 = new double[N];
  double *checkdot2 = new double[N];
  for(int i = 0; i < N; ++i){
    totalnorm2[i] = 0.0;
    checkdot2[i] = 0.0;
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	if(Parameters.basis == "finite_J"){ tempt = Amps2.D1.get_TJ(Space, Chan, chan, ind, i, j, a, b); }
	else{ tempt = Amps2.D1.get_T(chan, ind); }
	tempt += Ints.D_ME1.V4[chan][pp * nhh + hh];
	if( denom ){ tempt /= Amps2.D1.Evec[chan][ind]; }
	else{ tempt /= Amps.D1.Evec[chan][ind]; }
	tempt = mix*tempt + (1.0-mix)*Amps.D1.T1[chan][ind];
	totalnorm1 += (tempt - Amps.D1.T1[chan][ind]) * (tempt - Amps.D1.T1[chan][ind]);
	for(int k = 0; k < N; ++k){
	  totalnorm2[k] += delp[k][chan][ind] * delp[k][chan][ind];
	  checkdot2[k] += (tempt - Amps.D1.T1[chan][ind]) * delp[k][chan][ind];
	}
	//std::cout << std::setprecision(12) << "T: " << i << " " << j << " | " << a << " " << b << " = " << tempt << std::endl;
	if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i, j, a, b, tempt); }
	else{ Amps.D1.set_T(chan, ind, tempt); }
      }
    }
  }
  if(Parameters.approx == "singles"){
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      if(Parameters.basis == "finite_J"){ tempt = Amps2.S1.get_TJ(Space, hp1, i, a); }
      else{ tempt = Amps2.S1.get_T(hp1); }
      if( denom ){ tempt /= Amps2.S1.Evec[hp1]; }
      else{ tempt /= Amps.S1.Evec[hp1]; }
      tempt = mix*tempt + (1.0-mix)*Amps.S1.T1[hp1];
      totalnorm1 += (tempt - Amps.S1.T1[hp1]) * (tempt - Amps.S1.T1[hp1]);
      for(int k = 0; k < N; ++k){
	totalnorm2[k] += delp[k][Chan.size1][hp1] * delp[k][Chan.size1][hp1];
	checkdot2[k] += (tempt - Amps.S1.T1[hp1]) * delp[k][Chan.size1][hp1];
      }
      //std::cout << std::setprecision(12) << "t: " << i << " | " << a << " = " << tempt << std::endl;
      if(Parameters.basis == "finite_J"){ Amps.S1.set_TJ(Space, hp1, i, a, tempt); }
      else{ Amps.S1.set_T(hp1, tempt); }
    }
    if(Parameters.basis == "finite_J"){
      Amps.S1.set_T_2J(Space, Chan, Ints);
      Amps.D1.set_T_2J(Chan, Ints);
    }
    else{
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
  }
  checkdot = fabs(checkdot2[0]/std::sqrt(totalnorm1 * totalnorm2[0]));
  for(int i = 1; i < N; ++i){ checkdot = max(checkdot, fabs(checkdot2[i]/std::sqrt(totalnorm1 * totalnorm2[i]))); }
  if( !std::isfinite(checkdot) ){ checkdot = 1000; }
  delete[] totalnorm2;
  delete[] checkdot2;
}

void Gather_Amps3(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &mix, bool &denom, int &N, double ***p, double ***delp, double *B)
{
  int nhh, npp, a1, b1, i1, j1, ind;
  double tempt;
  for(int i = 0; i < N; ++i){
    B[(N + 1) * i + N - 1] = 0.0;
    B[(N + 1) * (N - 1) + i] = 0.0;
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a1 = Chan.ppvec[chan][2*pp];
      b1 = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i1 = Chan.hhvec[chan][2*hh];
	j1 = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	if(Parameters.basis == "finite_J"){ tempt = Amps2.D1.get_TJ(Space, Chan, chan, ind, i1, j1, a1, b1); }
	else{ tempt = Amps2.D1.get_T(chan, ind); }
	tempt += Ints.D_ME1.V4[chan][pp * nhh + hh];
	if( denom ){ tempt /= Amps2.D1.Evec[chan][ind]; }
	else{ tempt /= Amps.D1.Evec[chan][ind]; }
	tempt = mix*tempt + (1.0-mix)*Amps.D1.T1[chan][ind];

	p[N - 1][chan][ind] = Amps.D1.T1[chan][ind];
	delp[N - 1][chan][ind] = tempt - Amps.D1.T1[chan][ind];
	for(int i = 0; i < N; ++i){
	  B[(N + 1) * i + N - 1] += delp[i][chan][ind] * delp[N - 1][chan][ind];
	  if(i == N - 1){ break; } // don't double count
	  B[(N + 1) * (N - 1) + i] += delp[N - 1][chan][ind] * delp[i][chan][ind];
	}
	//std::cout << std::setprecision(12) << "T: " << i << " " << j << " | " << a << " " << b << " = " << tempt << std::endl;
	if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i1, j1, a1, b1, tempt); }
	else{ Amps.D1.set_T(chan, ind, tempt); }
      }
    }
  }
  if(Parameters.approx == "singles"){
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i1 = Chan.hp1vec[Chan.ind0][2*hp1];
      a1 = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      if(Parameters.basis == "finite_J"){ tempt = Amps2.S1.get_TJ(Space, hp1, i1, a1); }
      else{ tempt = Amps2.S1.get_T(hp1); }
      if( denom ){ tempt /= Amps2.S1.Evec[hp1]; }
      else{ tempt /= Amps.S1.Evec[hp1]; }
      tempt = mix*tempt + (1.0-mix)*Amps.S1.T1[hp1];

      p[N - 1][Chan.size1][hp1] = Amps.S1.T1[hp1];
      delp[N - 1][Chan.size1][hp1] = tempt - Amps.S1.T1[hp1];
      for(int i = 0; i < N; ++i){
	B[(N + 1) * i + N - 1] += delp[i][Chan.size1][hp1] * delp[N - 1][Chan.size1][hp1];
	if(i == N - 1){ break; } // don't double count
	B[(N + 1) * (N - 1) + i] += delp[N - 1][Chan.size1][hp1] * delp[i][Chan.size1][hp1];
      }
      if(Parameters.basis == "finite_J"){ Amps.S1.set_TJ(Space, hp1, i1, a1, tempt); }
      else{ Amps.S1.set_T(hp1, tempt); }
    }
    if(Parameters.basis == "finite_J"){
      Amps.S1.set_T_2J(Space, Chan, Ints);
      Amps.D1.set_T_2J(Chan, Ints);
    }
    else{
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
  }
  for(int i = 0; i < N; ++i){ B[(N + 1) * i + N] = -1.0; B[(N + 1) * N + i] = -1.0; }
  B[(N + 1) * N + N] = 0.0;
}

void CC_Error(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, Amplitudes &Amps2, double &error)
{
  int nhh, npp, a, b, i, j, ind;
  error = 0.0;
  double norm = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	//std::cout << "!error " << Amps.D1.T1[chan][ind] << " " << Amps2.D1.T1[chan][ind] << std::endl;;
	error += (Amps2.D1.T1[chan][ind] - Amps.D1.T1[chan][ind])*(Amps2.D1.T1[chan][ind] - Amps.D1.T1[chan][ind]);
	norm += Amps2.D1.T1[chan][ind] * Amps2.D1.T1[chan][ind];
      }
    }
  }
  if(Parameters.approx == "singles"){  
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      error += (Amps2.S1.T1[hp1] - Amps.S1.T1[hp1])*(Amps2.S1.T1[hp1] - Amps.S1.T1[hp1]);
      norm += Amps2.S1.T1[hp1] * Amps2.S1.T1[hp1];
    }
  }
  error = std::sqrt(error/norm);
}


void Randomize_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps0, Amplitudes &Amps, double &width)
{
  int nhh, npp, a, b, i, j, ind;
  double tempt;
  double rand;
  size_t key;
  std::unordered_map<size_t,double> t_map;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a > b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i > j){ continue; }
	tempt = Amps0.D1.T1[chan][hh * npp + pp];
	if(fabs(tempt) > 1.0e-12){
	  rand = rand_normal(0.0, width * fabs(tempt));
	  key = std::hash<float>{}(float(fabs(tempt)));
	  t_map[key] = rand; 	  
	}
      }
    }
  }
  if(Parameters.approx == "singles"){  
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      tempt = Amps0.S1.T1[hp1];
      if(fabs(tempt) > 1.0e-12){
	rand = rand_normal(0.0, width * fabs(tempt));
	key = std::hash<float>{}(float(fabs(tempt)));
	t_map[key] = rand;
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	ind = hh * npp + pp;
	tempt = Amps0.D1.T1[chan][ind];
	key = std::hash<float>{}(float(fabs(tempt)));
	if(tempt > 1.0e-12){
	  tempt += t_map[key];
	}
	else if(tempt < -1.0e-12){
	  tempt -= t_map[key];
	}
	if(Parameters.basis == "finite_J"){ Amps.D1.set_TJ(Space, Chan, chan, ind, i, j, a, b, tempt); }
	else{ Amps.D1.set_T(chan, ind, tempt); }
      }
    }
  }
  if(Parameters.approx == "singles"){  
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps0.S1.T1[hp1];
      key = std::hash<float>{}(float(fabs(tempt)));
      if(tempt > 1.0e-12){
	tempt += t_map[key];
      }
      else if(tempt < -1.0e-12){
	tempt -= t_map[key];
      }
      if(Parameters.basis == "finite_J"){ Amps.S1.set_TJ(Space, hp1, i, a, tempt); }
      else{ Amps.S1.set_T(hp1, tempt); }
    }
    if(Parameters.basis == "finite_J"){
      Amps.S1.set_T_2J(Space, Chan, Ints);
      Amps.D1.set_T_2J(Chan, Ints);
    }
    else{
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
  }
}

void Print_Amps(const Input_Parameters &Parameters, const Channels &Chan, Amplitudes &Amps)
{
  std::cout << std::endl;
  int nhh, npp, a, b, i, j;
  double tempt;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if((a >= b && Parameters.basis != "finite_J") || a > b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if((i >= j && Parameters.basis != "finite_J") || i > j){ continue; }
	tempt = Amps.D1.T1[chan][hh * npp + pp];
	if(Parameters.basis == "finite_J"){
	  std::cout << std::setprecision(12) << "T: < " << i << " " << j << " | " << a << " " << b << " >^" 
		    << 0.5*Chan.qnums1[chan].j << " = " << tempt << std::endl;
	}
	else{ std::cout << std::setprecision(12) << "T: < " << i << " " << j << " | " << a << " " << b << " > = " << tempt << std::endl; }
      }
    }
  }
  if(Parameters.approx == "singles"){  
    int nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps.S1.T1[hp1];
      std::cout << std::setprecision(12) << "t: " << i << " | " << a << " = " << tempt << std::endl;
    }
  }
  std::cout << std::endl;
}

void Doubles_Step_explicit(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2, double &error)
{
  State tb;
  int ind, ind1, ind2, key1, key2, key3, key4;
  int nh, np, nhh, npp, nhp, nhp1, nhp2, nhpp, nhhp;
  int a, b, i, j, k, l, c, d;
  double tempt, term0, term;
  double norm2 = 0.0;
  double error2 = 0.0;
  int print = 0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    nhp = Chan.nhp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }

	//if(i == 2 && j == 3 && a == 6 && b == 11){ print = 1; }
	//else{ print = 0; }
	//if(i < j && a < b){ print = 1; }
	//else{ print = 0; }

	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	  std::cout << "Term0 = < " << a << "," << b << " |V| " << i << "," << j << " > = " << Ints.D_ME1.V4[chan][pp * nhh + hh] << std::endl;
	}
	tempt = Ints.D_ME1.V4[chan][pp * nhh + hh];

	//T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  if(c == d){ continue; }
	  term0 = 0.5 * Ints.D_ME1.V1[chan][pp1 * npp + pp] * Amps1.D1.T1[chan][hh * npp + pp1];
	  term += term0;
	  if(print == 2){ std::cout << "Term1 += 0.5 * < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.D_ME1.V1[chan][pp1 * npp + pp] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	tempt += term;

	//T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int hh1 = 0; hh1 < nhh; ++hh1){
	  k = Chan.hhvec[chan][2*hh1];
	  l = Chan.hhvec[chan][2*hh1 + 1];
	  if(k == l){ continue; }
	  term0 = 0.5 * Ints.D_ME1.V2[chan][hh * nhh + hh1] * Amps1.D1.T1[chan][hh1 * npp + pp];
	  term += term0;
	  if(print == 2){ std::cout << "Term2 += 0.5 * < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.D_ME1.V2[chan][hh * nhh + hh1] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	tempt += term;

	//T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  if(c == d){ continue; }
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    term0 = 0.25 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.D1.T1[chan][hh * npp + pp1] * Amps1.D1.T1[chan][hh1 * npp + pp];
	    term += term0;
	    if(print == 2){ std::cout << "Term7 += 0.25 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.25 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	tempt += term;

	//T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || i == k){ continue; }
	  term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term3 += -1.0 * < " << k << "," << b << " |V| " << j << "," << c << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	tempt += term;

	//T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[j], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(j, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(i, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == c || j == k){ continue; }
	  term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){  std::cout << "Term4 += -1.0 * < " << k << "," << a << " |V| " << i << "," << c << " > * < " << b << "," << c << " |t| " << j << "," << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	tempt += term;

	//T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == c || i == k){ continue; }
	  term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term5 += < " << k << "," << a << " |V| " << j << "," << c << " > * < " << b << "," << c << " |t| " << i << "," << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	tempt += term;

	//T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[j], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(j, a, Space.indtot)];
 	key2 = Chan.hp2_map[ind][Hash2(i, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || j == k){ continue; }
	  term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term6 += < " << k << "," << b << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << j << "," << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	tempt += term;

	//T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || i == k){ continue; }
	  for(int hp1 = 0; hp1 < nhp1; ++hp1){
	    l = Chan.hp1vec[ind][2*hp1];
	    d = Chan.hp1vec[ind][2*hp1 + 1];
	    if(d == b || l == j || k == l || c == d){ continue; }
	    term0 = Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.D1.T2[ind][key1 * nhp2 + hp2] *  Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	tempt += term;

	//T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  d = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == d || i == k){ continue; }
	  for(int hp1 = 0; hp1 < nhp1; ++hp1){
	    l = Chan.hp1vec[ind][2*hp1];
	    c = Chan.hp1vec[ind][2*hp1 + 1];
	    if(c == a || l == j || k == l || c == d){ continue; }
	    term0 = Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.D1.T2[ind][key1 * nhp2 + hp2] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term13 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << d << " |t| " << i << "," << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	tempt += term;

	//T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[i];
	key1 = Chan.h_map[ind][i];
	key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	nh = Chan.nh[ind];
	nhpp = Chan.nhpp[ind];
	for(int h = 0; h < nh; ++h){
	  l = Chan.hvec[ind][h];
	  if(j == l){ continue; }
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(i == k || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V5[ind][h * nhpp + hpp] * Amps1.D1.T7[ind][key2 * nh + h] * Amps1.D1.T6[ind][hpp * nh + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << b << " |t| " << j << "," << l << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = -0.5 * " << Ints.D_ME1.V5[ind][h * nhpp + hpp] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " * " << Amps1.D1.T6[ind][hpp * nh + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	tempt += term;

	//T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[j];
	key1 = Chan.h_map[ind][j];
	key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	nh = Chan.nh[ind];
	nhpp = Chan.nhpp[ind];
	for(int h = 0; h < nh; ++h){
	  k = Chan.hvec[ind][h];
	  if(i == k){ continue; }
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    l = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(j == l || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V6[ind][h * nhpp + hpp] * Amps1.D1.T7[ind][key2 * nh + h] * Amps1.D1.T6[ind][hpp * nh + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term9 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << b << " |t| " << i << "," << k << " > * < " << c << "," << d << " |t| " << j << "," << l << " > = -0.5 * " << Ints.D_ME1.V6[ind][h * nhpp + hpp] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " * " << Amps1.D1.T6[ind][hpp * nh + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	tempt += term;

	//T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[a];
	key1 = Chan.p_map[ind][a];
	key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	np = Chan.np[ind];
	nhhp = Chan.nhhp[ind];
	for(int p = 0; p < np; ++p){
	  d = Chan.pvec[ind][p];
	  if(b == d){ continue; }
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(a == c || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V7[ind][p * nhhp + hhp] * Amps1.D1.T9[ind][key2 * np + p] * Amps1.D1.T8[ind][hhp * np + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term10 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V7[ind][p * nhhp + hhp] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " * " << Amps1.D1.T8[ind][hhp * np + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	tempt += term;

	//T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[b];
	key1 = Chan.p_map[ind][b];
	key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	np = Chan.np[ind];
	nhhp = Chan.nhhp[ind];
	for(int p = 0; p < np; ++p){
	  c = Chan.pvec[ind][p];
	  if(a == c){ continue; }
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    d = Chan.hhpvec[ind][3*hhp + 2];
	    if(b == d || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.D1.T9[ind][key2 * np + p] * Amps1.D1.T8[ind][hhp * np + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term11 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << j << " > * < " << b << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " * " << Amps1.D1.T8[ind][hhp * np + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	tempt += term;
	
	if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
	  if(print != 0){ std::cout << std::endl << "For Singles" << std::endl << std::endl; }

	  //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} = t1(c|i){ic}.t1(d|j){jd}.V1(ab|cd){cd,ab} (1)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    term0 = Ints.D_ME1.V1[chan][pp1 * npp + pp] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term1 += < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > = " << Ints.D_ME1.V1[chan][pp1 * npp + pp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	  tempt += term;
	  
	  //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} = V2(ij|kl){ij,kl}.t1(a|k){ka}.t1(b|l){lb} (2)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	    term0 = Ints.D_ME1.V2[chan][hh * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << Ints.D_ME1.V2[chan][hh * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	  tempt += term;

	  //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} = 0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.t1(a|k){ka}.t1(b|l){lb} (3)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      term0 = 0.5 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term3 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	  tempt += term;
	
	  //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} = 0.5 * t1(c|i){ic}.t1(d|j){jd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l){ continue; }
	      term0 = 0.5 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term4 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	  tempt += term;
	
	  //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} = t1(c|i){ic}.t1(d|j){jd}.V4(kl|cd){cd,kl}.t1(a|k){ka}.t1(b|l){lb} (5)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      term0 = Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	      term += term0;
	      if(print == 2){ std::cout << "Term5 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} = t1(c|i){ic}.t1(a|k){ka}.V3(kb|jc){kc,jb} (6)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term6 += < " << k << "," << b << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	  tempt += term;
	
	  //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} = t1(c|j){jc}.t1(b|k){kb}.V3(ka|ic){kc,ia}  (7)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(i, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term7 += < " << k << "," << a << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} = -t1(c|i){ic}.t1(b|k){kb}.V3(ka|jc){kc,ja}  (8)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term8 += -1.0 * < " << k << "," << a << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	  tempt += term;
	
	  //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} = -t1(c|j){jc}.t1(a|k){ka}.V3(kb|ic){kc,ib}  (9)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(i, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << b << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} = -t1(c|i){ic}.t1(a|k){ka}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      d = Chan.hp1vec[ind][2*hp1 + 1];
	      if(d == b || l == j || k == l || c == d){ continue; }
	      term0 = -1.0 * Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term10 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.t1(b|l){lb}.t1(d|j){jd} (11)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(a == c || i == k){ continue; }
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      d = Chan.hp1vec[ind][2*hp1 + 1];
	      if(k == l || c == d || Chan.indvec[l] != Chan.indvec[b] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      term0 = -1.0 * Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	      term += term0;
	      if(print == 2){ std::cout << "Term11 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << " |t| " << l << " > * < " << d << " |t| " << j << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} = -t1(d|i){id}.t1(b|k){kb}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    d = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      c = Chan.hp1vec[ind][2*hp1 + 1];
	      if(c == a || l == j || k == l || c == d){ continue; }
	      term0 = -1.0 * Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term12 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = -1.0 * " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.t1(a|l){la}.t1(c|j){jc} (13)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    d = Chan.hp2vec[ind][2*hp2 + 1];
	    if(b == d || i == k){ continue; }
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      c = Chan.hp1vec[ind][2*hp1 + 1];
	      if(k == l || c == d || Chan.indvec[l] != Chan.indvec[a] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      term0 = -1.0 * Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	      term += term0;
	      if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << l << " > * < " << c << " |t| " << j << " > * < " << b << "," << d << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} = -V11(ka|cd){a,kcd}.t1(c|i){ic}.T2(db|kj){kd,jb} (14)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == b || k == j || Chan.indvec[i] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term14 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << b << " |t| " << k << "," << j << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term14 = " << term << std::endl; }
	  tempt += term;

	  //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} = -V11(kb|cd){b,kcd}.t1(c|j){jc}.T2(da|ki){kd,ia} (15)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == a || k == i || Chan.indvec[j] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term15 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << a << " |t| " << k << "," << i << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term15 = " << term << std::endl; }
	  tempt += term;

	  //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} = V11(kb|cd){b,kcd}.t1(c|i){ic}.T2(da|kj){kd,ja} (16)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == a || k == j || Chan.indvec[i] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term16 += < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << a << " |t| " << k << "," << j << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term16 = " << term << std::endl; }
	  tempt += term;

	  //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} = V11(ka|cd){a,kcd}.t1(c|j){jc}.T2(db|ki){kd,ib} (17)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == b || k == i || Chan.indvec[j] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term17 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << b << " |t| " << k << "," << i << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term17 = " << term << std::endl; }
	  tempt += term;

	  //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} = -V12(kl|ic){i,klc}.t1(a|k){ka}.T2(cb|lj){lc,jb} (18)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == b || l == j || Chan.indvec[k] != Chan.indvec[a] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term18 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term18 = " << term << std::endl; }
	  tempt += term;

	  //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} = -V12(kl|jc){j,klc}.t1(b|k){kb}.T2(ca|li){lc,ia} (19)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == a || l == i || Chan.indvec[k] != Chan.indvec[b] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term19 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term19 = " << term << std::endl; }
	  tempt += term;
	  
	  //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} = V12(kl|ic){i,klc}.t1(b|k){kb}.T2(ca|lj){lc,ja} (20)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == a || l == j || Chan.indvec[k] != Chan.indvec[b] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term20 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term20 = " << term << std::endl; }
	  tempt += term;

	  //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} = V12(kl|jc){j,klc}.t1(a|k){ka}.T2(cb|li){lc,ib} (21)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == b || l == i || Chan.indvec[k] != Chan.indvec[a] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term21 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << i << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term21 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.t1(d|i){id}.t1(c|k){kc} (22)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == j){ continue; }
	    for(int p = 0; p < np; ++p){
	      d = Chan.pvec[ind][p];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hpp_map[ind][Hash3(k, c, d, Space.indtot)];
		term0 = Ints.D_ME1.V5[ind][h * nhpp + key4] * Amps1.S1.T1[key3] * Amps1.S1.T1[hp1] * Amps1.D1.T7[ind][key2 * nh + h];
		term += term0;
		if(print == 2){ std::cout << "Term22 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << j << "," << l << " > = " << Ints.D_ME1.V5[ind][h * nhpp + key4] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term22 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.t1(d|j){jd}.t1(c|l){lc} (23)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    if(k == i){ continue; }
	    for(int p = 0; p < np; ++p){
	      d = Chan.pvec[ind][p];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		l = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hpp_map[ind][Hash3(l, c, d, Space.indtot)];
		term0 = Ints.D_ME1.V6[ind][h * nhpp + key4] * Amps1.S1.T1[key3] * Amps1.S1.T1[hp1] * Amps1.D1.T7[ind][key2 * nh + h];
		term += term0;
		if(print == 2){ std::cout << "Term23 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << j << " > * < " << c << " |t| " << l << " > * < " << a << "," << b << " |t| " << i << "," << k << " > = " << Ints.D_ME1.V6[ind][h * nhpp + key4] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term23 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == j){ continue; }
	    key3 = Chan.hh1_map[Chan.ind0][Hash2(i, l, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(k == l){ continue; }
	      term0 = Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T6[ind][key2 * nh + h];
	      term += term0;
	      if(print == 2){ std::cout << "Term26 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = " << Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T6[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term26 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == i){ continue; }
	    key3 = Chan.hh1_map[Chan.ind0][Hash2(j, l, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(k == l){ continue; }
	      term0 = -1.0 * Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T6[ind][key2 * nh + h];
	      term += term0;
	      if(print == 2){ std::cout << "Term27 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = -1.0 * " << Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T6[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term27 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V17[ind][key2 * np + p] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term30 += -1.0 * < " << j << "," << c << " |V| " << a << "," << b << " > * < " << c << " |t| " << i << " > = -1.0 * " << Ints.S_ME1.V17[ind][key2 * np + p] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term30 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = Ints.S_ME1.V17[ind][key2 * np + p] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term31 += < " << i << "," << c << " |V| " << a << "," << b << " > * < " << c << " |t| " << j << " > = " << Ints.S_ME1.V17[ind][key2 * np + p] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term31 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  np = Chan.np[ind];
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    for(int p = 0; p < np; ++p){
	      c = Chan.pvec[ind][p];
	      plus(tb, Space.qnums[j], Space.qnums[c]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(j, c, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      term0 = -0.5 * Ints.S_ME1.V19[chan][key1 * nhh + hh1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term34 += -0.5 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = -0.5 * " << Ints.S_ME1.V19[chan][key1 * nhh + hh1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term34 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  np = Chan.np[ind];
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    for(int p = 0; p < np; ++p){
	      c = Chan.pvec[ind][p];
	      plus(tb, Space.qnums[i], Space.qnums[c]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(i, c, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      term0 = 0.5 * Ints.S_ME1.V19[chan][key1 * nhh + hh1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term35 += 0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.S_ME1.V19[chan][key1 * nhh + hh1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term35 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(k == l || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[l] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term38 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << a << " |t| " << l << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term38 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(k == l || Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[l] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term39 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > * < " << a << " |t| " << l << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term39 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == b){ continue; }
	    for(int h = 0; h < nh; ++h){
	      l = Chan.hvec[ind][h];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hhp_map[ind][Hash3(k, l, c, Space.indtot)];
		term0 = Ints.D_ME1.V7[ind][p * nhhp + key4] * Amps1.S1.T1[hp1] * Amps1.S1.T1[key3] * Amps1.D1.T9[ind][key2 * np + p];
		term += term0;
		if(print == 2){ std::cout << "Term24 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << a << " |t| " << l << " > * < " << b << "," << d << " |t| " << i << "," << j << " > = " << Ints.D_ME1.V7[ind][p * nhhp + key4] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term24 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    if(c == a){ continue; }
	    for(int h = 0; h < nh; ++h){
	      l = Chan.hvec[ind][h];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		d = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hhp_map[ind][Hash3(k, l, d, Space.indtot)];
		term0 = Ints.D_ME1.V8[ind][p * nhhp + key4] * Amps1.S1.T1[hp1] * Amps1.S1.T1[key3] * Amps1.D1.T9[ind][key2 * np + p];
		term += term0;
		if(print == 2){ std::cout << "Term25 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << a << "," << c << " |t| " << i << "," << j << " > = " << Ints.D_ME1.V8[ind][p * nhhp + key4] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term25 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == b){ continue; }
	    key3 = Chan.pp1_map[Chan.ind0][Hash2(d, a, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(c == d){ continue; }
	      term0 = Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T8[ind][key2 * np + p];
	      term += term0;
	      if(print == 2){ std::cout << "Term28 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = " << Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T8[ind][key2 * np + p] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term28 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  np = Chan.np[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == a){ continue; }
	    key3 = Chan.pp1_map[Chan.ind0][Hash2(d, b, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(c == d){ continue; }
	      term0 = -1.0 * Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T8[ind][key2 * np + p];
	      term += term0;
	      if(print == 2){ std::cout << "Term29 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = -1.0 * " << Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T8[ind][key2 * np + p] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term29 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V18[ind][key2 * nh + h] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term32 += -1.0 * < " << i << "," << j << " |V| " << k << "," << b << " > * < " << a << " |t| " << k << " > = -1.0 * " << Ints.S_ME1.V18[ind][key2 * nh + h] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term32 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = Ints.S_ME1.V18[ind][key2 * nh + h] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term33 += < " << i << "," << j << " |V| " << k << "," << a << " > * < " << b << " |t| " << k << " > = " << Ints.S_ME1.V18[ind][key2 * nh + h] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term33 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  nh = Chan.nh[ind];
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int h = 0; h < nh; ++h){
	      k = Chan.hvec[ind][h];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(k, b, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      term0 = -0.5 * Ints.S_ME1.V20[chan][pp1 * nhp + key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term36 += -0.5 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = -0.5 * " << Ints.S_ME1.V20[chan][pp1 * nhp + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term36 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  nh = Chan.nh[ind];
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int h = 0; h < nh; ++h){
	      k = Chan.hvec[ind][h];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(k, a, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	      term0 = 0.5 * Ints.S_ME1.V20[chan][pp1 * nhp + key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term37 += 0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.S_ME1.V20[chan][pp1 * nhp + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term37 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(c == d || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term40 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << j << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term40 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(c == d || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term41 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << j << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term41 = " << term << std::endl; }
	  tempt += term;
	}

	if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.D1.Evec[chan][hh * npp + pp] << ", " << tempt/Amps1.D1.Evec[chan][hh * npp + pp] << std::endl; }

	tempt /= Amps1.D1.Evec[chan][hh * npp + pp];
	Amps2.D1.T1[chan][hh * npp + pp] = tempt;
      }
    }
  }

  if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
    if(print != 0){ std::cout << std::endl << "For Singles" << std::endl; }
    nhp1 = Chan.nhp1[Chan.ind0];
    
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];

      //if(i == 0 && a == 10){ print = 1; }
      //else{ print = 0; }

      if(print != 0){
	std::cout << std::endl;
	std::cout << "!! < " << a << " |t| " << i << " > !!" << std::endl;
	std::cout << "Term0 += 0.0" << std::endl;
      }
      tempt = 0.0;
      
      //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hp2 = 0; hp2 < nhp1; ++hp2){
	k = Chan.hp1vec[Chan.ind0][2*hp2];
	c = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
	term0 = -1.0 * Ints.D_ME1.V3[Chan.ind0][hp1 * nhp1 + hp2] * Amps1.S1.T1[hp2];
	term += term0;
	if(print == 2){ std::cout << "Term1 += -1.0 * < " << k << "," << a << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[Chan.ind0][hp1 * nhp1 + hp2] << " * " << Amps1.S1.T1[hp2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
      tempt += term;

      //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hp2 = 0; hp2 < nhp1; ++hp2){
	k = Chan.hp1vec[Chan.ind0][2*hp2];
	c = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
	for(int hp3 = 0; hp3 < nhp1; ++hp3){
	  l = Chan.hp1vec[Chan.ind0][2*hp3];
	  d = Chan.hp1vec[Chan.ind0][2*hp3 + 1];
	  if(k == l || c == d || d == a || l == i){ continue; }
	  term0 = Ints.D_ME1.V9[Chan.ind0][hp2 * nhp1 + hp3] * Amps1.S1.T1[hp2] * Amps1.D1.T2[Chan.ind0][hp3 * nhp1 + hp1];
	  term += term0;
	  if(print == 2){ std::cout << "Term6 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > * " << Ints.D_ME1.V9[Chan.ind0][hp2 * nhp1 + hp3] << " * " << Amps1.S1.T1[hp2] << " * " << Amps1.D1.T2[Chan.ind0][hp3 * nhp1 + hp1] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = -0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      ind = Chan.indvec[a];
      key1 = Chan.p_map[ind][a];
      key2 = Chan.h_map[ind][i];
      nh = Chan.nh[ind];
      nhpp = Chan.nhpp[ind];
      np = Chan.np[ind];
      nhhp = Chan.nhhp[ind];
      for(int hpp = 0; hpp < nhpp; ++hpp){
	k = Chan.hppvec[ind][3*hpp];
	c = Chan.hppvec[ind][3*hpp + 1];
	d = Chan.hppvec[ind][3*hpp + 2];
	if(c == d || i == k){ continue; }
	term0 = -0.5 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.D1.T6[ind][hpp * nh + key2];
	term += term0;
	if(print == 2){ std::cout << "Term2 += -0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = -0.5 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.D1.T6[ind][hpp * nh + key2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} = -V11(ka|cd){a,kcd}.t1(c|i){ic}.t1(d|k){kd} (3)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hpp = 0; hpp < nhpp; ++hpp){
	k = Chan.hppvec[ind][3*hpp];
	c = Chan.hppvec[ind][3*hpp + 1];
	d = Chan.hppvec[ind][3*hpp + 2];
	if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[d]){ continue; }
	key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	key4 = Chan.hp1_map[Chan.ind0][Hash2(k, d, Space.indtot)];
	term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	term += term0;
	if(print == 2){ std::cout << "Term3 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << k << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int h = 0; h < nh; ++h){
	k = Chan.hvec[ind][h];
	for(int hpp = 0; hpp < nhpp; ++hpp){
	  l = Chan.hppvec[ind][3*hpp];
	  c = Chan.hppvec[ind][3*hpp + 1];
	  d = Chan.hppvec[ind][3*hpp + 2];
	  if(k == l || c == d || i == l || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	  term0 = -0.5 * Ints.D_ME1.V6[ind][h * nhpp + hpp] * Amps1.S1.T1[key1] * Amps1.D1.T6[ind][hpp * nh + key2];
	  term += term0;
	  if(print == 2){ std::cout << "Term7 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << l << " > = -0.5 * " << Ints.D_ME1.V6[ind][h * nhpp + hpp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T6[ind][hpp * nh + key2] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      ind = Chan.indvec[i];
      key1 = Chan.h_map[ind][i];
      key2 = Chan.p_map[ind][a];
      np = Chan.np[ind];
      nhhp = Chan.nhhp[ind];
      nh = Chan.nh[ind];
      nhpp = Chan.nhpp[ind];
      for(int hhp = 0; hhp < nhhp; ++hhp){
	k = Chan.hhpvec[ind][3*hhp];
	l = Chan.hhpvec[ind][3*hhp + 1];
	c = Chan.hhpvec[ind][3*hhp + 2];
	if(k == l || a == c){ continue; }
	term0 = -0.5 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.D1.T8[ind][hhp * np + key2];
	term += term0;
	if(print == 2){ std::cout << "Term4 += -0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.D1.T8[ind][hhp * np + key2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} = -V12(kl|ic){i,klc}.t1(a|k){ka}.t1(c|l){lc} (5)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hhp = 0; hhp < nhhp; ++hhp){
	k = Chan.hhpvec[ind][3*hhp];
	l = Chan.hhpvec[ind][3*hhp + 1];
	c = Chan.hhpvec[ind][3*hhp + 2];
	if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[c]){ continue; }
	key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	key4 = Chan.hp1_map[Chan.ind0][Hash2(l, c, Space.indtot)];
	term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	term += term0;
	if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << " |t| " << l << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int p = 0; p < np; ++p){
	c = Chan.pvec[ind][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  k = Chan.hhpvec[ind][3*hhp];
	  l = Chan.hhpvec[ind][3*hhp + 1];
	  d = Chan.hhpvec[ind][3*hhp + 2];
	  if(k == l || c == d || Chan.indvec[i] != Chan.indvec[c]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  term0 = -0.5 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.S1.T1[key1] * Amps1.D1.T8[ind][hhp * np + key2];
	  term += term0;
	  if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T8[ind][hhp * np + key2] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} = -t3(i|c){i,c}.V8(kl|cd){c,kld}.t1(a|k){ka}.t1(d|l){ld} (9)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int p = 0; p < np; ++p){
	c = Chan.pvec[ind][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  k = Chan.hhpvec[ind][3*hhp];
	  l = Chan.hhpvec[ind][3*hhp + 1];
	  d = Chan.hhpvec[ind][3*hhp + 2];
	  if(k == l || c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[d]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	  key3 = Chan.hp1_map[Chan.ind0][Hash2(l, d, Space.indtot)];
	  term0 = -1.0 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  term += term0;
	  if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << l << " > = -1.0 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
      tempt += term;

      if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.S1.Evec[hp1] << ", " << tempt/Amps1.S1.Evec[hp1] << std::endl; }
      tempt /= Amps1.S1.Evec[hp1];
      Amps2.S1.T1[hp1] = tempt;
    }
  }

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }
	tempt = Amps2.D1.T1[chan][hh * npp + pp];
	norm2 += tempt * tempt;
	error2 += (tempt - Amps1.D1.T1[chan][hh * npp + pp]) * (tempt - Amps1.D1.T1[chan][hh * npp + pp]);
	Amps1.D1.set_T(chan, hh * npp + pp, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){  
    nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps2.S1.T1[hp1];
      norm2 += tempt * tempt;
      error2 += (tempt - Amps1.S1.T1[hp1]) * (tempt - Amps1.S1.T1[hp1]);
      Amps1.S1.set_T(hp1, tempt);
    }
  }
  error = std::sqrt(error2/norm2);
}

void Doubles_Step_explicit2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps1, Amplitudes &Amps2, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, double &error)
{
  State tb;
  int ind1, ind2, ind3, ind4, Vind1, Vind2, key1, key2, key3, key4, Vkey1, Vkey2;
  int nhh, npp, nhp1, i, j, a, b;
  double tempt, term0, term, TBME;
  double norm2 = 0.0;
  double error2 = 0.0;
  int print = 0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }

	//if(i == 2 && j == 3 && a == 6 && b == 11){ print = 1; }
	//else{ print = 0; }
	//if(i < j && a < b){ print = true; }
	//else{ print = false; }

	plus(tb, Space.qnums[a], Space.qnums[b]);
	Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	plus(tb, Space.qnums[i], Space.qnums[j]);
	Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	if(Vind1 != Vind2){ continue; }
	Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	tempt = TBME;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	  std::cout << "Term0 = < " << a << "," << b << " |V| " << i << "," << j << " > = " << TBME << std::endl;
	}

	term = 0.0;
 	if(print == 2){ std::cout << std::endl; }
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(c == d){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, j, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term1 += 0.5 * < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(k, l, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, b, Space.indtot)];
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(i, j, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(k, l, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += 0.5 * < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind3 != ind4){ continue; }
	    key3 = Chan.hh_map[ind3][Hash2(k, l, Space.indtot)];
	    key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[i], Space.qnums[j]);
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2){ continue; }
		key1 = Chan.hh_map[ind1][Hash2(i, j, Space.indtot)];
		key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.25 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term7 += 0.25 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.25 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term3 += < " << k << "," << b << " |V| " << c << "," << j << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term4 += -1.0 * < " << k << "," << b << " |V| " << c << "," << i << " > * < " << a << "," << c << " |t| " << j << "," << k << " > = -1.0 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == b){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[b], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << a << " |V| " << c << "," << j << " > * < " << b << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == b){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[b], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term6 += < " << k << "," << a << " |V| " << c << "," << i << " > * < " << b << "," << c << " |t| " << j << "," << k << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  plus(tb, Space.qnums[i], Space.qnums[k]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == j){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	      key2 = Chan.pp_map[ind1][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[b], Space.qnums[d]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key3 = Chan.hh_map[ind3][Hash2(j, l, Space.indtot)];
		key4 = Chan.pp_map[ind3][Hash2(b, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << k << " > * < " << b << "," << d << " |t| " << j << "," << l << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  plus(tb, Space.qnums[j], Space.qnums[k]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == i){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	      key2 = Chan.pp_map[ind1][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[b], Space.qnums[d]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key3 = Chan.hh_map[ind3][Hash2(i, l, Space.indtot)];
		key4 = Chan.pp_map[ind3][Hash2(b, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << j << "," << k << " > * < " << b << "," << d << " |t| " << i << "," << l << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[l], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(l, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << l << "," << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[l], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(l, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == b){ continue; }
	      plus(tb, Space.qnums[b], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term9 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << c << " |t| " << l << "," << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == j){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[l], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[c]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[a], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2 || ind3 != ind4){ continue; }
		key2 = Chan.pp_map[ind2][Hash2(d, c, Space.indtot)];
		key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term10 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << "," << c << " |t| " << i << "," << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == i){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[l], Space.qnums[i]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[c]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[a], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2 || ind3 != ind4){ continue; }
		key2 = Chan.pp_map[ind2][Hash2(d, c, Space.indtot)];
		key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term11 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << "," << c << " |t| " << j << "," << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	tempt += term;

	
	if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
	  if(print != 0){ std::cout << std::endl << "For Singles" << std::endl << std::endl; }
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term1 += < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	  tempt += term;
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[i], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term2 += < " << k << "," << l << " |V| " << i << "," << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(c == d){ continue; }
		  plus(tb, Space.qnums[i], Space.qnums[j]);
		  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(c, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term3 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(k, l, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term4 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term5 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term6 += -1.0 * < " << k << "," << b << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[i]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term7 += < " << k << "," << b << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term8 += < " << k << "," << a << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[i]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << a << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == b){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[b]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term10 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == b){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[b]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term11 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == a){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[a]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == a){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[a]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == j){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[j]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[a], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term14 += < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << b << " |t| " << k << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term14 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == i){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[i]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[a], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term15 += -1.0 * < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << b << " |t| " << k << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term15 = " << term << std::endl; }
	  tempt += term;
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == j){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[j]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[b], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(b, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term16 += -1.0 * < " << b << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << a << " |t| " << k << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term16 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == i){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[i]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[b], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(b, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term17 += < " << b << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << a << " |t| " << k << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term17 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == b){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[i], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term18 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term18 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == b){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[j], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(j, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term19 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term19 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == a){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[i], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term20 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term20 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == a){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[j], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(j, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term21 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term21 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	   	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	   	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	   	for(int d = Space.indhol; d < Space.indtot; ++d){
	   	  if(Chan.indvec[d] != Chan.indvec[i]){ continue; }
	   	  key2 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	   	  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	   	  term += term0;
	   	  if(print == 2){ std::cout << "Term22 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	   	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term22 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	   	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	   	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	   	for(int d = Space.indhol; d < Space.indtot; ++d){
	   	  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	   	  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	   	  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	   	  term += term0;
	   	  if(print == 2){ std::cout << "Term23 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	   	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term23 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[a]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	for(int d = Space.indhol; d < Space.indtot; ++d){
	  	  if(d == b){ continue; }
	  	  plus(tb, Space.qnums[i], Space.qnums[j]);
	  	  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  plus(tb, Space.qnums[d], Space.qnums[b]);
	  	  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  if(ind3 != ind4){ continue; }
	  	  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	  	  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	  	  term += term0;
	  	  if(print == 2){ std::cout << "Term24 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << a << " |t| " << l << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	  	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term24 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	for(int d = Space.indhol; d < Space.indtot; ++d){
	  	  if(d == a){ continue; }
	  	  plus(tb, Space.qnums[i], Space.qnums[j]);
	  	  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  plus(tb, Space.qnums[d], Space.qnums[a]);
	  	  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  if(ind3 != ind4){ continue; }
	  	  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	  	  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	  	  term += term0;
	  	  if(print == 2){ std::cout << "Term25 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	  	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term25 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term26 += -1.0 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term26 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term27 += < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term27 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(d == b){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[d], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term28 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term28 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(d == a){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[d], Space.qnums[a]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term29 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term29 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term30 += < " << a << "," << b << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term30 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term31 += -1.0 * < " << a << "," << b << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term31 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term32 += -1.0 * < " << k << "," << b << " |V| " << i << "," << j << " > * < " << a << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term32 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term33 += < " << k << "," << a << " |V| " << i << "," << j << " > * < " << b << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term33 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term34 += 0.5 * < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term34 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term35 += -0.5 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term35 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(c == d){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[c], Space.qnums[d]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term36 += -0.5 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term36 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(c == d){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[c], Space.qnums[d]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term37 += 0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term37 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term38 += < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term38 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term39 += -1.0 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term39 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	  	key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term40 += < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term40 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
		term += term0;
		if(print == 2){ std::cout << "Term41 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << d << " |t| " << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term41 = " << term << std::endl; }
	  tempt += term;
	}

	if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.D1.Evec[chan][hh * npp + pp] << ", " << tempt/Amps1.D1.Evec[chan][hh * npp + pp] << std::endl; }
	
	tempt /= Amps1.D1.Evec[chan][hh * npp + pp];
	Amps2.D1.T1[chan][hh * npp + pp] = tempt;
      }
    }
  }

  if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
    if(print != 0){ std::cout << std::endl << "For Singles" << std::endl; }
    nhp1 = Chan.nhp1[Chan.ind0];
    
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];

      //if(i == 0 && a == 10){ print = 1; }
      //else{ print = 0; }

      if(print != 0){ 
	std::cout << std::endl;
	std::cout << "!! < " << a << " |t| " << i << " > !!" << std::endl;
	std::cout << "Term0 += 0.0" << std::endl;
      }
      tempt = 0.0;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  plus(tb, Space.qnums[a], Space.qnums[k]);
	  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[i], Space.qnums[c]);
	  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  if(Vind1 != Vind2){ continue; }
	  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  term0 = TBME * Amps1.S1.T1[key1];
	  term += term0;
	  if(print == 2){ std::cout << "Term1 += < " << a << "," << k << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
      tempt += term;      

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(k == i){ continue; }
	plus(tb, Space.qnums[i], Space.qnums[k]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(c == d){ continue; }
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[k]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += 0.5 * < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];	  
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(Chan.indvec[d] != Chan.indvec[k]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[k]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term3 += < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(k == l){ continue; }
	  plus(tb, Space.qnums[k], Space.qnums[l]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(k, l, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[c]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term4 += -0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind2][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[l]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(l, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[c]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
      tempt += term;
      
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(l == i){ continue; }
	  plus(tb, Space.qnums[l], Space.qnums[i]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(d == a){ continue; }
	      plus(tb, Space.qnums[d], Space.qnums[a]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term6 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  if(l == i){ continue; }
	  plus(tb, Space.qnums[i], Space.qnums[l]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(c == d){ continue; }
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(i, l, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term7 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(k == l){ continue; }
	  plus(tb, Space.qnums[k], Space.qnums[l]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(d == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[d]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(a, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(Chan.indvec[d] != Chan.indvec[l]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
      tempt += term;

      if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.S1.Evec[hp1] << ", " << tempt/Amps1.S1.Evec[hp1] << std::endl; }
      tempt /= Amps1.S1.Evec[hp1];
      Amps2.S1.T1[hp1] = tempt;
    }
  }

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }
	tempt = Amps2.D1.T1[chan][hh * npp + pp];
	norm2 += tempt * tempt;
	error2 += (tempt - Amps1.D1.T1[chan][hh * npp + pp]) * (tempt - Amps1.D1.T1[chan][hh * npp + pp]);
	Amps1.D1.set_T(chan, hh * npp + pp, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){  
    nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps2.S1.T1[hp1];
      norm2 += tempt * tempt;
      error2 += (tempt - Amps1.S1.T1[hp1]) * (tempt - Amps1.S1.T1[hp1]);
      Amps1.S1.set_T(hp1, tempt);
    }
  }
  error = std::sqrt(error2/norm2);
}

void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double p1 = 1.0, p12 = 0.5, p14 = 0.25, m1 = -1.0, m12 = -0.5, zero = 0.0;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh == 0 || pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &p12, &p1, &N, &N);
    //T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &p12, &p1, &N, &N);
    //T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &p1, &zero, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &p14, &p1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
    dgemm_NN(Amps1.D1.T2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &m1, &p1, &N, &N);
    //T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
    //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
    for(int i = 0; i < hp1 * hp2; ++i){
      Amps2.D1.T3[chan][i] = Amps2.D1.T2[chan][i];
      Amps2.D1.T4[chan][i] = -1.0 * Amps2.D1.T2[chan][i];
      Amps2.D1.T5[chan][i] = -1.0 * Amps2.D1.T2[chan][i];
    }
    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &p1, &p1, &N, &N);
    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &p1, &p1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.D1.T6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &m12, &p1, &N, &N);
    //T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &m12, &p1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.D1.T8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &m12, &p1, &N, &N);
    //T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &p1, &zero, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &m12, &p1, &N, &N);
  }
}

void Doubles_Step_J(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = 0.25, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh == 0 || pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac3, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac3, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac4, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
    dgemm_NN(Amps1.D1.T2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
    //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
    for(int i = 0; i < hp1 * hp2; ++i){
      Amps2.D1.T3[chan][i] = Amps2.D1.T2[chan][i];
      Amps2.D1.T4[chan][i] = Amps2.D1.T2[chan][i];
      Amps2.D1.T5[chan][i] = Amps2.D1.T2[chan][i];
    }
    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.D1.T6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac6, &fac1, &N, &N); //fac1
    //T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac6, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.D1.T8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac6, &fac1, &N, &N); //fac1
    //T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac6, &fac1, &N, &N); //fac1
  }
}

void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = -1.0, fac4 = -0.5;
  char N = 'N';
  int hp = Chan.nhp1[Chan.ind0];
  int one = 1;
  if(hp != 0){
    //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
    dgemm_NN(Ints.D_ME1.V3[Chan.ind0], Amps1.S1.T1, Amps2.S1.T1, &hp, &one, &hp, &fac3, &fac1, &N, &N); //fac1
    //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
    dgemm_NN(Ints.D_ME1.V9[Chan.ind0], Amps1.D1.T2[Chan.ind0], Amps2.S1.S3, &hp, &hp, &hp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.T1, Amps2.S1.S3, Amps2.S1.T1, &one, &hp, &hp, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p != 0 && h != 0){
      if(hpp != 0){
	//t2(a|i){a,i} = -0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.D1.T6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac4, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.S1.E6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac1, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
	dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.S1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T2[chan], Amps2.S1.S2[chan], Amps2.S1.T2[chan], &p, &h, &h, &fac4, &fac1, &N, &N); //fac1
      }
      if(hhp != 0){
	//t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.D1.T8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.S1.E8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac1, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac1, &fac1, &N, &N); //fac1
      }
    }
  }
}

void Singles_Step_J(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = -1.0, fac4 = -0.5, fac5 = 0.5;
  char N = 'N';
  int hp = Chan.nhp1[Chan.ind0];
  int one = 1;
  if(hp != 0){
    //t1(ia){ia} = V3(ka|ic){ia,kc}.t1(kc){kc} (1)
    dgemm_NN(Ints.D_ME1.V3[Chan.ind0], Amps1.S1.T1, Amps2.S1.T1, &hp, &one, &hp, &fac1, &fac1, &N, &N); //fac1
    //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
    dgemm_NN(Ints.D_ME1.V9[Chan.ind0], Amps1.D1.T2[Chan.ind0], Amps2.S1.S3, &hp, &hp, &hp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.T1, Amps2.S1.S3, Amps2.S1.T1, &one, &hp, &hp, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p != 0 && h != 0){
      if(hpp != 0){
	//t2(a|i){a,i} = 0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.D1.T6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac5, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = -V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3)
	dgemm_NN(Ints.S_ME1.V11[chan], Amps1.S1.E6[chan], Amps2.S1.T2[chan], &p, &h, &hpp, &fac3, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
	dgemm_NN(Ints.D_ME1.V6[chan], Amps1.D1.T6[chan], Amps2.S1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T2[chan], Amps2.S1.S2[chan], Amps2.S1.T2[chan], &p, &h, &h, &fac4, &fac1, &N, &N); //fac1
      }
      if(hhp != 0){
	//t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.D1.T8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5)
	dgemm_NN(Ints.S_ME1.V12[chan], Amps1.S1.E8[chan], Amps2.S1.T3[chan], &h, &p, &hhp, &fac1, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.D1.T8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9)
	dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.S1.S1[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
	dgemm_NN(Amps1.S1.T3[chan], Amps2.S1.S1[chan], Amps2.S1.T3[chan], &h, &p, &p, &fac1, &fac1, &N, &N); //fac1
      }
    }
  }
}
 
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh * pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac5, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac5, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (3)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (5)
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (6)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (7)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T3[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (8)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (9)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T5[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} (11)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.S1.E2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.S1.E2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} (14)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} (15)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} (16)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} (17)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} (18)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} (19)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} (20)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} (21)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} (22)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.S1.E6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} (23)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.S1.E6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac5, &fac1, &N, &N); //fac1
    if(p != 0){
      //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac6, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac3, &fac1, &N, &N); //fac1
      //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.S1.E8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac5, &fac1, &N, &N); //fac1
    if(h != 0){
      //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
      //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac6, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac3, &fac1, &N, &N); //fac1
      //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
    }
  }
}

void Doubles_Step_2_J(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.nhh[chan];
    int pp = Chan.npp[chan];
    if(hh * pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac5, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} (2)
    dgemm_NN(Ints.D_ME1.V2[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac5, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (3)
    dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
    dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (5)
    dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan];
    int hp2 = Chan.nhp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (6)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (7)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T3[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (8)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (9)
    dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T5[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} (11)
    dgemm_NN(Ints.D_ME1.V9[chan], Amps1.S1.E2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} (13)
    dgemm_NN(Ints.D_ME1.V10[chan], Amps1.S1.E2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} (14)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} (15)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} (16)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} (17)
    dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} (18)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} (19)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
    //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} (20)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} (21)
    dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.np[chan];
    int h = Chan.nh[chan];
    int hpp = Chan.nhpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} (22)
    dgemm_NN(Ints.D_ME1.V5[chan], Amps1.S1.E6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} (23)
    dgemm_NN(Ints.D_ME1.V6[chan], Amps1.S1.E6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
    //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
    dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac5, &fac1, &N, &N); //fac1
    if(p != 0){
      //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
      dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac6, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
      dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac3, &fac1, &N, &N); //fac1
      //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
      dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.nh[chan];
    int p = Chan.np[chan];
    int hhp = Chan.nhhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
    dgemm_NN(Ints.D_ME1.V7[chan], Amps1.S1.E8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
    dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
    dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
    //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
    dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac5, &fac1, &N, &N); //fac1
    if(h != 0){
      //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
      dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
      //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac6, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
      dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac3, &fac1, &N, &N); //fac1
      //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
      dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
    }
  }
}

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff)
{
  double fac1 = 1.0, fac5 = -1.0, fac4 = -0.5, fac3 = 0.5;
  int one = 1, hp0 = Chan.nhp1[Chan.ind0], pp0 = Chan.npp1[Chan.ind0], hh0 = Chan.nhh1[Chan.ind0];
  char N = 'N';

  if(hp0 != 0){
    // X_ia1(i|a){ia} = V9(ik|ac){ia,kc}.t1(c|k){kc}
    dgemm_NN(Ints.D_ME1.V9[Chan.ind0], Amps.S1.T1, V_Eff.X_ia1, &hp0, &one, &hp0, &fac1, &fac1, &N, &N);
    V_Eff.set_X_ia(Chan);
  }

  if(pp0 != 0){
    // X_ab1(a|b){ba} = f_ab.delta(a,b)
    for(int pp = 0; pp < pp0; ++pp){
      if(Chan.pp1vec[Chan.ind0][2*pp] == Chan.pp1vec[Chan.ind0][2*pp + 1]){
	V_Eff.X_ab1[pp] += Space.qnums[Chan.pp1vec[Chan.ind0][2*pp]].energy;
      }
    }
    if(hp0 != 0){
      // X_ab1(a|b){ba} = V16(ka|cb){ba,kc}.t1(c|k){kc}
      dgemm_NN(Ints.S_ME1.V16[Chan.ind0], Amps.S1.T1, V_Eff.X_ab1, &pp0, &one, &hp0, &fac1, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ab3(a|b){b,a} = -X_ia3(k|b){b,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X_ia3[chan], Amps.S1.T3[chan], V_Eff.X_ab3[chan], &p1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && hhp1 != 0){
      // X_ab3(a|b){b,a} = -(1/2).V8(kl|bc){b,klc}.T8(ac|kl){klc,a}
      dgemm_NN(Ints.D_ME1.V8[chan], Amps.D1.T8[chan], V_Eff.X_ab3[chan], &p1, &p1, &hhp1, &fac4, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ab(Chan);

  if(hh0 != 0){
    // X1_ij1(i|j){ji} = f_ij.delta(i,j)
    for(int hh = 0; hh < hh0; ++hh){
      if(Chan.hh1vec[Chan.ind0][2*hh] == Chan.hh1vec[Chan.ind0][2*hh + 1]){
	V_Eff.X1_ij1[hh] += Space.qnums[Chan.hh1vec[Chan.ind0][2*hh]].energy;
      }
    }
    if(hp0 != 0){
      // X1_ij1(i|j){ji} = -V15(ki|jc){ji,kc}.t1(c|k){kc}
      dgemm_NN(Ints.S_ME1.V15[Chan.ind0], Amps.S1.T1, V_Eff.X1_ij1, &hh0, &one, &hp0, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], hpp1 = Chan.nhpp[chan];
    if(h1 != 0 && hpp1 != 0){
      // X1_ij2(i|j){i,j} = (1/2).V6(ik|cd){i,kcd}.T6(cd|jk){kcd,j}
      dgemm_NN(Ints.D_ME1.V6[chan], Amps.D1.T6[chan], V_Eff.X1_ij2[chan], &h1, &h1, &hpp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_ij(Chan);

  // X_ij1(i|j){ji} = X1_ij1(i|j){ji}
  if(hh0 != 0){
    for(int hh = 0; hh < hh0; ++hh){
      V_Eff.X_ij1[hh] = V_Eff.X1_ij1[hh];
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan];
    if(h1 != 0 && p1 != 0){
      // X_ij2(i|j){i,j} = X_ia2(i|d){i,d}.t2(d|j){d,j}
      dgemm_NN(V_Eff.X_ia2[chan], Amps.S1.T2[chan], V_Eff.X_ij2[chan], &h1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ij(Chan);

  if(hp0 != 0){
    // X_ai1(a|i){ia} = -V3(ka|ic){ia,kc}.t1(c|k){kc}
    dgemm_NN(Ints.D_ME1.V3[Chan.ind0], Amps.S1.T1, V_Eff.X_ai1, &hp0, &one, &hp0, &fac5, &fac1, &N, &N);
    // X_ai1(a|i){ia} = T2(ac|ik){ia,kc}.X_ia1(k|c){kc}
    dgemm_NN(Amps.D1.T2[Chan.ind0], V_Eff.X_ia1, V_Eff.X_ai1, &hp0, &one, &hp0, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan], hpp1 = Chan.nhpp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ai2(a|i){a,i} = X_ab2(a|c){a,c}.t2(c|i){c,i}
      dgemm_NN(V_Eff.X_ab2[chan], Amps.S1.T2[chan], V_Eff.X_ai2[chan], &p1, &h1, &p1, &fac1, &fac1, &N, &N);
      if(hpp1 != 0){
	// X_ai2(a|i){a,i} = -(1/2).V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i}
	dgemm_NN(Ints.S_ME1.V11[chan], Amps.D1.T6[chan], V_Eff.X_ai2[chan], &p1, &h1, &hpp1, &fac4, &fac1, &N, &N);
      }
      // X_ai3(a|i){i,a} = -X1_ij3(k|i){i,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_ij3[chan], Amps.S1.T3[chan], V_Eff.X_ai3[chan], &h1, &p1, &h1, &fac5, &fac1, &N, &N);
      if(hhp1 != 0){
	// X_ai3(a|i){i,a} = -(1/2).V12(kl|ic){i,klc}.T8(ac|kl){klc,a}
	dgemm_NN(Ints.S_ME1.V12[chan], Amps.D1.T8[chan], V_Eff.X_ai3[chan], &h1, &p1, &hhp1, &fac4, &fac1, &N, &N);
      }
    }
  }
  V_Eff.set_X_ai(Chan);

  // X_ijab1(ij|ab){ab,ij} = V4(ij|ab){ab,ij}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.nhh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_ijab1[chan][i] = Ints.D_ME1.V4[chan][i];
      }
    }
  }

  // X(1)_iabc(ia|bc){a,ibc} = V11(ia|bc){a,ibc}
  // X(1)_ijka(ij|ka){k,ija} = V12(ij|ka){k,ija}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.np[chan] * Chan.nhpp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_iabc1[chan][i] = Ints.S_ME1.V11[chan][i];
	V_Eff.X_iabc1[chan][i] = Ints.S_ME1.V11[chan][i];
      }
    }
    length = Chan.nh[chan] * Chan.nhhp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_ijka1[chan][i] = Ints.S_ME1.V12[chan][i];
	V_Eff.X_ijka1[chan][i] = Ints.S_ME1.V12[chan][i];
      }
    }
  }

  for(int chan = 0; chan < Chan.size3; ++chan){
    int p1 = Chan.np[chan], h1 = Chan.nh[chan], hpp1 = Chan.nhpp[chan], hhp1 = Chan.nhhp[chan];
    if(p1 != 0 && h1 != 0 && hpp1 != 0){
      // X(1)_iabc(ia|bc){a,ibc} = -((1/2)).t2(a|k){a,k}.V5(ik|bc){k,ibc}
      dgemm_NN(Amps.S1.T2[chan], Ints.D_ME1.V5[chan], V_Eff.X1_iabc1[chan], &p1, &hpp1, &h1, &fac4, &fac1, &N, &N);
      dgemm_NN(Amps.S1.T2[chan], Ints.D_ME1.V5[chan], V_Eff.X_iabc1[chan], &p1, &hpp1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && h1 != 0 && hhp1 != 0){
      // X(1)_ijka(ij|ka){k,ija} = -((1/2)).t3(c|k){k,c}.V7(ij|ac){c,ija}
      dgemm_NN(Amps.S1.T3[chan], Ints.D_ME1.V7[chan], V_Eff.X1_ijka1[chan], &h1, &hhp1, &p1, &fac4, &fac1, &N, &N);
      dgemm_NN(Amps.S1.T3[chan], Ints.D_ME1.V7[chan], V_Eff.X_ijka1[chan], &h1, &hhp1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iabc(Chan);
  V_Eff.set_X1_iabc(Chan);
  V_Eff.set_X_ijka(Chan);
  V_Eff.set_X1_ijka(Chan);

  // X1_abcd(ab|cd){cd,ab} = V1(ab|cd){cd,ab}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.npp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_abcd1[chan][i] = Ints.D_ME1.V1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int ppp1 = Chan.nppp[chan], p1 = Chan.np[chan], h1 = Chan.nh[chan];
    if(ppp1 != 0 && p1 != 0 && h1 != 0){
      // X1_abcd(ab|cd){acd,b} = X1_iabc(ka|cd}{acd,k}.t3(b|k){k,b}
      dgemm_NN(V_Eff.X1_iabc2[chan], Amps.S1.T3[chan], V_Eff.X1_abcd2[chan], &ppp1, &p1, &h1, &fac1, &fac1, &N, &N);
      // X1_abcd(ab|cd){bcd,a} = -X1_iabc(kb|cd}{bcd,k}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_iabc2[chan], Amps.S1.T3[chan], V_Eff.X1_abcd3[chan], &ppp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_abcd(Chan);

  // X_abcd(ab|cd){cd,ab} = X1_abcd(ab|cd){cd,ab}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.npp[chan] * Chan.npp[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_abcd1[chan][i] = V_Eff.X1_abcd1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.npp[chan], hh1 = Chan.nhh[chan];
    if(pp1 != 0 && hh1 != 0){
      // X_abcd(ab|cd){cd,ab} = (1/2).V4(kl|cd}{cd,kl}.T1(ab|kl){kl,ab}
      dgemm_NN(Ints.D_ME1.V4[chan], Amps.D1.T1[chan], V_Eff.X_abcd1[chan], &pp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }

  // X_ijkl(ij|kl){ij,kl} = V2(ij|kl){ij,kl}
  for(int chan = 0; chan < Chan.size1; ++chan){
    int length = Chan.nhh[chan] * Chan.nhh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_ijkl1[chan][i] = Ints.D_ME1.V2[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int hhh1 = Chan.nhhh[chan], h1 = Chan.nh[chan], p1 = Chan.np[chan];
    if(hhh1 != 0 && h1 != 0 && p1 != 0){
      // X_ijkl(ij|kl){ijk,l} = X1_ijkl(ij|kc}{ijk,c}.t3(c|l){c,l}
      dgemm_NN(V_Eff.X1_ijka2[chan], Amps.S1.T2[chan], V_Eff.X_ijkl2[chan], &hhh1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X_ijkl(ij|kl){ijl,k} = -X1_ijkl(ij|lc}{ijl,c}.t3(c|k){c,k}
      dgemm_NN(V_Eff.X1_ijka2[chan], Amps.S1.T2[chan], V_Eff.X_ijkl3[chan], &hhh1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh1 = Chan.nhh[chan], pp1 = Chan.npp[chan];
    if(hh1 != 0 && pp1 != 0){
      // X_ijkl(ij|kl){kl,ij} = (1/2).T1(kl|cd}{kl,cd}.V4(ij|cd){cd,ij}
      dgemm_NN(Amps.D1.T1[chan], Ints.D_ME1.V4[chan], V_Eff.X_ijkl4[chan], &hh1, &hh1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ijkl(Chan);

  // X1_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  // X3_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  for(int chan = 0; chan < Chan.size2; ++chan){
    int length = Chan.nhp2[chan] * Chan.nhp2[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X1_iajb1[chan][i] = Ints.D_ME1.V3[chan][i];
	V_Eff.X3_iajb1[chan][i] = Ints.D_ME1.V3[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp1[chan], hpp1 = Chan.nhpp1[chan];
    if(h1 * p1 * hhp1 != 0){
      // X1_iajb(ia|jb){ijb,a} = (1/2).V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      dgemm_NN(Ints.S_ME1.V14[chan], Amps.S1.T3[chan], V_Eff.X1_iajb2[chan], &hhp1, &p1, &h1, &fac3, &fac1, &N, &N);
      // X3_iajb(ia|jb){ijb,a} = V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      dgemm_NN(Ints.S_ME1.V14[chan], Amps.S1.T3[chan], V_Eff.X3_iajb2[chan], &hhp1, &p1, &h1, &fac1, &fac1, &N, &N);
    }
    if(h1 * p1 * hpp1 != 0){
      // X1_iajb(ia|jb){iab,j} = X1_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X1_iabc3[chan], Amps.S1.T2[chan], V_Eff.X1_iajb3[chan], &hpp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X3_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X_iabc3[chan], Amps.S1.T2[chan], V_Eff.X3_iajb3[chan], &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_iajb(Chan);
  V_Eff.set_X3_iajb(Chan);

  // X_iajb(ia|jb){ib,ja} = X3_iajb(ia|jb){ib,ja}
  for(int chan = 0; chan < Chan.size2; ++chan){
    int length = Chan.nhp2[chan] * Chan.nhp2[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_iajb1[chan][i] = V_Eff.X3_iajb1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(hp1 != 0 && hp2 != 0){
      // X_iajb(ia|jb){ib,ja} = -V10(ik|cb}{ib,kc}.T5(ca|jk){kc,ja}
      dgemm_NN(Ints.D_ME1.V10[chan], Amps.D1.T5[chan], V_Eff.X_iajb1[chan], &hp2, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hpp1 = Chan.nhpp1[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      dgemm_NN(V_Eff.X_iabc3[chan], Amps.S1.T2[chan], V_Eff.X_iajb3[chan], &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iajb(Chan);

  // X_abic(ab|ic){iab,c} = V17(ab|ic){iab,c}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhpp[chan] * Chan.np[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_abic1[chan][i] = Ints.S_ME1.V17[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], ppp1 = Chan.nppp[chan], hpp1 = Chan.nhpp[chan], hpp2 = Chan.nhpp1[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_abic(ab|ic){iab,c} = -T7(ab|ik}{iab,k}.X_ia(k|c){k,c}  
      dgemm_NN(Amps.D1.T7[chan], V_Eff.X_ia2[chan], V_Eff.X_abic1[chan], &hpp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && ppp1 != 0){
      // X_abic(ab|ic){abc,i} = V*(ab|dc){abc,d}.t2(d|i){d,i}
      dgemm_NN(V_Eff.V_abcd[chan], Amps.S1.T2[chan], V_Eff.X_abic2[chan], &ppp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hpp2 != 0){
      // X_abic(ab|ic){icb',a'} = -X1_iajb(kb|ic){icb',k'}.t3(a|k){k,a}
      dgemm_NN(V_Eff.X1_iajb4[chan], Amps.S1.T3[chan], V_Eff.X_abic3[chan], &hpp2, &p1, &h1, &fac5, &fac1, &N, &N);
      // X_abic(ab|ic){ica',b'} = X1_iajb(ka|ic){ica',k'}.t3(b|k){k,b}
      dgemm_NN(V_Eff.X1_iajb4[chan], Amps.S1.T3[chan], V_Eff.X_abic4[chan], &hpp2, &p1, &h1, &fac1, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int pp1 = Chan.npp1[chan], hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(pp1 != 0 && hp1 != 0 && hp2 != 0){
      // X_abic(ab|ic){bc',ia'} = X_iabc(kb|dc){bc',kd}.T3(ad|ik){kd,ia'}
      dgemm_NN(V_Eff.X_iabc4[chan], Amps.D1.T3[chan], V_Eff.X_abic5[chan], &pp1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X_abic(ab|ic){ac',ib'} = -X_iabc(ka|dc){ac',kd}.T3(bd|ik){kd,ib'}
      dgemm_NN(V_Eff.X_iabc4[chan], Amps.D1.T3[chan], V_Eff.X_abic6[chan], &pp1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.npp[chan], hh1 = Chan.nhh[chan], hp1 = Chan.nhp[chan];
    if(pp1 != 0 && hh1 != 0 && hp1 != 0){
      // X_abic(ab|ic){ic,ab} = (1/2).X_ijka(kl|ic){ic,kl}.T1(ab|kl){kl,ab}
      dgemm_NN(V_Eff.X_ijka5[chan], Amps.D1.T1[chan], V_Eff.X_abic7[chan], &hp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_abic(Chan);

  // X2_iajk(ia|jk){jka,i} = V18(ia|jk){jka,i}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhhp[chan] * Chan.nh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X2_iajk1[chan][i] = Ints.S_ME1.V18[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhh1 = Chan.nhhh[chan], hhp1 = Chan.nhhp1[chan];
    if(h1 != 0 && p1 != 0 && hhh1 != 0){
      // X2_iajk(ia|jk){ijk,a} = -V*(il|jk){ijk,l}.t3(a|l){l,a}
      dgemm_NN(V_Eff.V_ijkl[chan], Amps.S1.T3[chan], V_Eff.X2_iajk2[chan], &hhh1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X2_iajk(ia|jk){jia',k'} = X3_iajb(ia|jd){jia',d'}.t2(d|k){d,k}
      dgemm_NN(V_Eff.X3_iajb5[chan], Amps.S1.T2[chan], V_Eff.X2_iajk3[chan], &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){kia',j'} = -X3_iajb(ia|kd){kia',d'}.t2(d|j){d,j}
      dgemm_NN(V_Eff.X3_iajb5[chan], Amps.S1.T2[chan], V_Eff.X2_iajk4[chan], &hhp1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hh1 = Chan.nhh1[chan], hp1 = Chan.nhp1[chan], hp2 = Chan.nhp2[chan];
    if(hh1 != 0 && hp1 != 0 && hp2 != 0){
      // X2_iajk(ia|jk){ij',ka'} = X_ijka(il|jc){ij',lc}.T3(ac|kl){lc,ka'}
      dgemm_NN(V_Eff.X_ijka4[chan], Amps.D1.T3[chan], V_Eff.X2_iajk5[chan], &hh1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){ik',ja'} = -X_ijka(il|kc){ik',lc}.T3(ac|jl){lc,ja'}
      dgemm_NN(V_Eff.X_ijka4[chan], Amps.D1.T3[chan], V_Eff.X2_iajk6[chan], &hh1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hp1 = Chan.nhp[chan], hh1 = Chan.nhh[chan], pp1 = Chan.npp[chan];
    if(hh1 != 0 && pp1 != 0 && hp1 != 0){
      // X2_iajk(ia|jk){jk,ia} = (1/2).T1(cd|jk){jk,cd}.X_iabc(ia|cd){cd,ia}
      dgemm_NN(Amps.D1.T1[chan], V_Eff.X_iabc5[chan], V_Eff.X2_iajk7[chan], &hh1, &hp1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X2_iajk(Chan);

  // X_iajk(ia|jk){jka,i} = X2_iajk(ia|jk){jka,i}
  for(int chan = 0; chan < Chan.size3; ++chan){
    int length = Chan.nhhp[chan] * Chan.nh[chan];
    if(length != 0){
      for(int i = 0; i < length; ++i){
	V_Eff.X_iajk1[chan][i] = V_Eff.X2_iajk1[chan][i];
      }
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.nh[chan], p1 = Chan.np[chan], hhp1 = Chan.nhhp[chan];
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X_iajk(ia|jk){jka,i} = T8(ca|jk){jka,c}.X_ia(i|c){c,i}
      dgemm_NN(Amps.D1.T8[chan], V_Eff.X_ia3[chan], V_Eff.X_iajk1[chan], &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }
}

void PA_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums)
{
  double *Ham;
  double tempen;
  State state;
  int ind;
  int *pvec2;
  int *hppvec2;
  int length;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }

    int count01 = 0;
    int count02 = 0;
    count01 = Chan.np[chan];
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // h
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // pp
	if(Chan.nh[chan1] * Chan.npp[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int pp = 0; pp < Chan.npp[chan2]; ++pp){
	      if(Chan.ppvec[chan2][2*pp] < Chan.ppvec[chan2][2*pp + 1]){
		count02 += Chan.nh[chan1];
	      }
	    }
	  }
	}
      }
    }

    pvec2 = new int[count01];
    hppvec2 = new int[3 * count02];
    count01 = 0;
    count02 = 0;
    for(int i = 0; i < Chan.np[chan]; ++i){
      pvec2[count01] = Chan.pvec[chan][i];
      ++count01;
    }
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // h
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // pp
	if(Chan.nh[chan1] * Chan.npp[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int pp = 0; pp < Chan.npp[chan2]; ++pp){
	      if(Chan.ppvec[chan2][2*pp] < Chan.ppvec[chan2][2*pp + 1]){
		for(int h = 0; h < Chan.nh[chan1]; ++h){
		  hppvec2[3*count02] = Chan.hvec[chan1][h];
		  hppvec2[3*count02 + 1] = Chan.ppvec[chan2][2*pp];
		  hppvec2[3*count02 + 2] = Chan.ppvec[chan2][2*pp + 1];
		  ++count02;
		}
	      }
	    }
	  }
	}
      }
    }

    int N = count01 + count02;
    if(N == 0){ continue; }
    Ham = new double[N*N];

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, ind0, ind;
      int bra, ket, key1, key2;
      State tb;
      length = count01 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + bra] = 0.0;
	b1 = pvec2[bra];
	k1 = pvec2[ket];
	ME = 0.0;
	ind = Chan.pp1_map[Chan.ind0][Hash2(k1, b1, Space.indtot)];
	ME += V_Eff.X_ab1[ind];
	Ham[N*ket + bra] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + (bra + count01)] = 0.0;
	b1 = hppvec2[3*bra];
	b2 = hppvec2[3*bra + 1];
	b3 = hppvec2[3*bra + 2];
	k1 = pvec2[ket];
	ME = 0.0;
	ind0 = Chan.indvec[k1];
	key1 = Chan.hpp_map[ind0][Hash3(b1, b2, b3, Space.indtot)];
	key2 = Chan.p_map[ind0][k1];
	ind = key1 * Chan.np[ind0] + key2;
	ME -= V_Eff.X_abic1[ind0][ind];
	Ham[N*ket + (bra + count01)] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%count01);
	ket = int((braket - bra)/count01);
	Ham[N*(ket + count01) + bra] = 0.0;
	k1 = hppvec2[3*ket];
	k2 = hppvec2[3*ket + 1];
	k3 = hppvec2[3*ket + 2];
	b1 = pvec2[bra];
	ME = 0.0;
	if(k2 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k3, Space.indtot)];
	  ME += V_Eff.X_ia1[ind];
	}
	if(k3 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k2, Space.indtot)];
	  ME -= V_Eff.X_ia1[ind];
	}
	ind0 = Chan.indvec[b1];
	key1 = Chan.p_map[ind0][b1];
	key2 = Chan.hpp_map[ind0][Hash3(k1, k2, k3, Space.indtot)];
	ind = key1 * Chan.nhpp[ind0] + key2;
	ME -= V_Eff.X_iabc1[ind0][ind];
	Ham[N*(ket + count01) + bra] += ME;
      }
      length = count02 * count02;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count02);
	bra = int((braket - ket)/count02);
	Ham[N*(ket + count01) + (bra + count01)] = 0.0;
	b1 = hppvec2[3*bra];
	b2 = hppvec2[3*bra + 1];
	b3 = hppvec2[3*bra + 2];
	k1 = hppvec2[3*ket];
	k2 = hppvec2[3*ket + 1];
	k3 = hppvec2[3*ket + 2];
	ME = 0.0;
	if(b1 == k1){
	  if(b3 == k3){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k2, b2, Space.indtot)];
	    ME += V_Eff.X_ab1[ind];
	  }
	  if(b2 == k2){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b3, Space.indtot)];
	    ME += V_Eff.X_ab1[ind];
	  }
	  if(b3 == k2){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b2, Space.indtot)];
	    ME -= V_Eff.X_ab1[ind];
	  }
	  if(b2 == k3){
	    ind = Chan.pp1_map[Chan.ind0][Hash2(k2, b3, Space.indtot)];
	    ME -= V_Eff.X_ab1[ind];
	  }
	}
	if(b2 == k2 && b3 == k3){
	  ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	  ME -= V_Eff.X_ij1[ind];
	}
	if(b1 == k1){
	  plus(tb, Space.qnums[b2], Space.qnums[b3]);
	  ind0 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[k2], Space.qnums[k3]);
	  if(ind0 == ChanInd_2b_dir(Parameters.basis, Space, tb)){
	    key1 = Chan.pp_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.pp_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.npp[ind0] + key2;
	    ME += V_Eff.X_abcd1[ind0][ind];
	  }
	}
	if(b3 == k3){
	  minus(tb, Space.qnums[k2], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b3 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b2], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k3){
	  minus(tb, Space.qnums[k2], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	Ham[N*(ket + count01) + (bra + count01)] += ME;
      }
    }

    double norm1p;
    int info = 0;
    if(N <= 3){
      char job1 = 'N';
      char job2 = 'V';
      int lwork = 10*N;
      double *Vl = new double[N * N];
      double *Vr = new double[N * N];
      double *wl = new double[N];
      double *wr = new double[N];
      double *work = new double[lwork];
      info = 0;
      dgeev_(&job1, &job2, &N, Ham, &N, wr, wl, Vl, &N, Vr, &N, work, &lwork, &info);

      ind = 0;
      tempen = wr[0];
      for(int i = 1; i < N; ++i){
	if(tempen - wr[i] >= 1e-8){
	  tempen = wr[i];
	  ind = i;
	}
      }
      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += Vr[N*ind + i] * Vr[N*ind + i]; }
      norm1p = std::sqrt(norm1p);

      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = wr[ind];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;
      delete[] Vl;
      delete[] Vr;
      delete[] wl;
      delete[] wr;
      delete[] work;
    }
    else{
      int ido = 0; // status integer is zero at start
      char bmat[] = "I"; // standard eigenvalue problem
      char which[] = "SR"; // smallest magnitude eigenvalues
      int nev = 1; // number of eigenvalues to calculate
      double tol = 1.0e-10; // error tolerance
      double *resid = new double[N];
      int ncv = 3*nev + 2; //5*nev; // number of lanczos vectors
      if(ncv > N){ ncv = N; }
      double *v = new double[N*ncv];
      for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
      int ldv = N;
      double *workd = new double[3*N];
      int lworkl = 4*ncv*(ncv + 2);
      double *workl = new double[lworkl];
      info = 0;
      int iparam[11];
      int ipntr[14];
      
      int ishift = 1;
      int mxiter = 5000;
      int mode = 1;
      iparam[0] = ishift;
      iparam[2] = mxiter;
      iparam[6] = mode;
      
      do{
	dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
	if(ido != -1 && ido != 1){ break; }

        #pragma omp parallel for
	for(int j = 0; j < N; ++j){
	  workd[ipntr[1]-1 + j] = 0.0;
	  for(int k = 0; k < N; ++k){
	    workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k];
	  }
	}
      } while(true);
      
      bool rvec = true;
      char howmny = 'A';
      int *select = new int[ncv];
      double *dr = new double[nev + 1];
      double *di = new double[nev + 1];
      double *z = new double[N * (nev + 1)];
      for(int i = 0; i < nev + 1; ++i){
	dr[i] = 0.0;
	di[i] = 0.0;
	for(int j = 0; j < N; ++j){
	  z[N*i + j] = 0.0;
	}
      }
      int ldz = N;
      double sigmar;
      double sigmai;
      double *workev = new double[3 * ncv];
      dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += z[i] * z[i]; }
      norm1p = std::sqrt(norm1p);
      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = dr[0];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;
   
      delete[] resid;
      delete[] v;
      delete[] workd;
      delete[] workl;
      delete[] select;
      delete[] dr;
      delete[] di;
      delete[] z;
      delete[] workev;
    }
    delete[] pvec2;
    delete[] hppvec2;
    delete[] Ham;
  }
}

void PR_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff, State *states, double *nums)
{
  double *Ham;
  double tempen;
  State state;
  int ind;
  int *hvec2;
  int *hhpvec2;
  int length;

  for(int chan = 0; chan < Chan.size3; ++chan){
    if(Chan.qnums3[chan].m != 1){ continue; }
    if(Chan.qnums3[chan].ml != 0 && Chan.qnums3[chan].ml != 1 && Chan.qnums3[chan].ml != 2){ continue; }

    int count01 = 0;
    int count02 = 0;
    count01 = Chan.nh[chan];
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
	if(Chan.np[chan1] * Chan.nhh[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int hh = 0; hh < Chan.nhh[chan2]; ++hh){
	      if(Chan.hhvec[chan2][2*hh] < Chan.hhvec[chan2][2*hh + 1]){
		count02 += Chan.np[chan1];
	      }
	    }
	  }
	}
      }
    }

    hvec2 = new int[count01];
    hhpvec2 = new int[3 * count02];
    count01 = 0;
    count02 = 0;
    for(int i = 0; i < Chan.nh[chan]; ++i){
      hvec2[count01] = Chan.hvec[chan][i];
      ++count01;
    }
    for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
      for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
	if(Chan.np[chan1] * Chan.nhh[chan2] != 0){
	  minus(state, Chan.qnums1[chan2], Chan.qnums3[chan1]);
	  if( equal(state, Chan.qnums3[chan]) ){
	    for(int hh = 0; hh < Chan.nhh[chan2]; ++hh){
	      if(Chan.hhvec[chan2][2*hh] < Chan.hhvec[chan2][2*hh + 1]){
		for(int p = 0; p < Chan.np[chan1]; ++p){
		  hhpvec2[3*count02] = Chan.hhvec[chan2][2*hh];
		  hhpvec2[3*count02 + 1] = Chan.hhvec[chan2][2*hh + 1];
		  hhpvec2[3*count02 + 2] = Chan.pvec[chan1][p];
		  ++count02;
		}
	      }
	    }
	  }
	}
      }
    }

    int N = count01 + count02;
    if(N == 0){ continue; }
    Ham = new double[N*N];

    #pragma omp parallel
    {
      double ME;
      int b1, b2, b3, k1, k2, k3, ind0, ind;
      int bra, ket, key1, key2;
      State tb;
      length = count01 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + bra] = 0.0;
	b1 = hvec2[bra];
	k1 = hvec2[ket];
	ME = 0.0;
	ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	ME -= V_Eff.X_ij1[ind];
	Ham[N*ket + bra] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count01);
	bra = int((braket - ket)/count01);
	Ham[N*ket + (bra + count01)] = 0.0;
	b1 = hhpvec2[3*bra];
	b2 = hhpvec2[3*bra + 1];
	b3 = hhpvec2[3*bra + 2];
	k1 = hvec2[ket];
	ME = 0.0;
	ind0 = Chan.indvec[k1];
	key1 = Chan.hhp_map[ind0][Hash3(b1, b2, b3, Space.indtot)];
	key2 = Chan.h_map[ind0][k1];
	ind = key1 * Chan.nh[ind0] + key2;
	ME += V_Eff.X_iajk1[ind0][ind];
	Ham[N*ket + (bra + count01)] += ME;
      }
      length = count02 * count01;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	bra = int(braket%count01);
	ket = int((braket - bra)/count01);
	Ham[N*(ket + count01) + bra] = 0.0;
	k1 = hhpvec2[3*ket];
	k2 = hhpvec2[3*ket + 1];
	k3 = hhpvec2[3*ket + 2];
	b1 = hvec2[bra];
	ME = 0.0;
	if(k1 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k2, k3, Space.indtot)];
	  ME -= V_Eff.X_ia1[ind];
	}
	if(k2 == b1){
	  ind = Chan.hp2_map[Chan.ind0][Hash2(k1, k3, Space.indtot)];
	  ME += V_Eff.X_ia1[ind];
	}
	ind0 = Chan.indvec[b1];
	key1 = Chan.h_map[ind0][b1];
	key2 = Chan.hhp_map[ind0][Hash3(k1, k2, k3, Space.indtot)];
	ind = key1 * Chan.nhhp[ind0] + key2;
	ME += V_Eff.X_ijka1[ind0][ind];
	Ham[N*(ket + count01) + bra] += ME;
      }
      length = count02 * count02;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%count02);
	bra = int((braket - ket)/count02);
	Ham[N*(ket + count01) + (bra + count01)] = 0.0;
	b1 = hhpvec2[3*bra];
	b2 = hhpvec2[3*bra + 1];
	b3 = hhpvec2[3*bra + 2];
	k1 = hhpvec2[3*ket];
	k2 = hhpvec2[3*ket + 1];
	k3 = hhpvec2[3*ket + 2];
	ME = 0.0;
	if(b3 == k3){
	  if(b2 == k2){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k1, Space.indtot)];
	    ME -= V_Eff.X_ij1[ind];
	  }
	  if(b1 == k1){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b2, k2, Space.indtot)];
	    ME -= V_Eff.X_ij1[ind];
	  }
	  if(b2 == k1){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b1, k2, Space.indtot)];
	    ME += V_Eff.X_ij1[ind];
	  }
	  if(b1 == k2){
	    ind = Chan.hh1_map[Chan.ind0][Hash2(b2, k1, Space.indtot)];
	    ME += V_Eff.X_ij1[ind];
	  }
	}
	if(b1 == k1 && b2 == k2){
	  ind = Chan.pp1_map[Chan.ind0][Hash2(k3, b3, Space.indtot)];
	  ME += V_Eff.X_ab1[ind];
	}
	if(b3 == k3){
	  plus(tb, Space.qnums[b1], Space.qnums[b2]);
	  ind0 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[k1], Space.qnums[k2]);
	  if(ind0 == ChanInd_2b_dir(Parameters.basis, Space, tb)){
	    key1 = Chan.hh_map[ind0][Hash2(k1, k2, Space.indtot)];
	    key2 = Chan.hh_map[ind0][Hash2(b1, b2, Space.indtot)];
	    ind = key1 * Chan.nhh[ind0] + key2;
	    ME += V_Eff.X_ijkl1[ind0][ind];
	  }
	}
	if(b2 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b1 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME -= V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b2 == k1){
	  minus(tb, Space.qnums[k3], Space.qnums[k2]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b1]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k2, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b1, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	if(b1 == k2){
	  minus(tb, Space.qnums[k3], Space.qnums[k1]);
	  ind0 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  minus(tb, Space.qnums[b3], Space.qnums[b2]);
	  if(ind0 == ChanInd_2b_cross(Parameters.basis, Space, tb)){
	    key1 = Chan.hp2_map[ind0][Hash2(k1, k3, Space.indtot)];
	    key2 = Chan.hp2_map[ind0][Hash2(b2, b3, Space.indtot)];
	    ind = key1 * Chan.nhp2[ind0] + key2;
	    ME += V_Eff.X_iajb1[ind0][ind];
	  }
	}
	Ham[N*(ket + count01) + (bra + count01)] += ME;
      }
    }

    double norm1p;
    int info = 0;
    if(N <= 3){
      char job1 = 'N';
      char job2 = 'V';
      int lwork = 10*N;
      double *Vl = new double[N * N];
      double *Vr = new double[N * N];
      double *wl = new double[N];
      double *wr = new double[N];
      double *work = new double[lwork];
      info = 0;
      dgeev_(&job1, &job2, &N, Ham, &N, wr, wl, Vl, &N, Vr, &N, work, &lwork, &info);

      ind = 0;
      tempen = wr[0];
      for(int i = 1; i < N; ++i){
	if(tempen - wr[i] <= 1e-8){
	  tempen = wr[i];
	  ind = i;
	}
      }
      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += Vr[N*ind + i] * Vr[N*ind + i]; }
      norm1p = std::sqrt(norm1p);      
      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = wr[ind];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;      
      delete[] Vl;
      delete[] Vr;
      delete[] wl;
      delete[] wr;
      delete[] work;
    }
    else{
      int ido = 0; // status integer is zero at start
      char bmat[] = "I"; // standard eigenvalue problem
      char which[] = "SR"; // smallest real eigenvalues
      int nev = 1; // number of eigenvalues to calculate
      double tol = 1.0e-10; // error tolerance
      double *resid = new double[N];
      int ncv = 3*nev + 2; //5*nev; // number of lanczos vectors
      if(ncv > N){ ncv = N; }
      double *v = new double[N*ncv];
      for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
      int ldv = N;
      double *workd = new double[3*N];
      int lworkl = 4*ncv*(ncv + 2);
      double *workl = new double[lworkl];
      info = 0;
      int iparam[11];
      int ipntr[14];
      
      int ishift = 1;
      int mxiter = 5000;
      int mode = 1;
      iparam[0] = ishift;
      iparam[2] = mxiter;
      iparam[6] = mode;
      
      do{
	dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	if(ido != -1 && ido != 1){ break; }

        #pragma omp parallel for
	for(int j = 0; j < N; ++j){
	  workd[ipntr[1]-1 + j] = 0.0;
	  for(int k = 0; k < N; ++k){
	    workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k];
	  }
	}
      } while(true);
      
      bool rvec = true;
      char howmny = 'A';
      int *select = new int[ncv];
      double *dr = new double[nev + 1];
      double *di = new double[nev + 1];
      double *z = new double[N * (nev + 1)];
      for(int i = 0; i < nev + 1; ++i){
	dr[i] = 0.0;
	di[i] = 0.0;
	for(int j = 0; j < N; ++j){
	  z[N*i + j] = 0.0;
	}
      }
      int ldz = N;
      double sigmar;
      double sigmai;
      double *workev = new double[3 * ncv];

      dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
      norm1p = 0.0;
      for(int i = 0; i < count01; ++i){ norm1p += z[i] * z[i]; }
      norm1p = std::sqrt(norm1p);
      states[Chan.qnums3[chan].ml] = Chan.qnums3[chan];
      nums[2*Chan.qnums3[chan].ml] = dr[0];
      nums[2*Chan.qnums3[chan].ml + 1] = norm1p;
      delete[] resid;
      delete[] v;
      delete[] workd;
      delete[] workl;
      delete[] select;
      delete[] dr;
      delete[] di;
      delete[] z;
      delete[] workev;
    }
    delete[] hvec2;
    delete[] hhpvec2;
    delete[] Ham;
  }
}

void CC_compare_JM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, std::string &inputfile)
{
  Input_Parameters Parameters2;
  Model_Space Space2;
  Channels Chan2;
  Amplitudes Amps2;
  Interactions Ints2;
  HF_Channels HF_Chan2;
  HF_Matrix_Elements HF_ME2;
  Single_Particle_States States2;
  int JM2 = 0;
  State tb1, tb2, tb3, tb4;
  int nh_j, np_j, nhh_j, npp_j, nhpp_j, nhhp_j, nhp1_j, nhp2_j, nhp0_j, chan_j, chan0_j;
  int nh_m, np_m, nhh_m, npp_m, nhpp_m, nhhp_m, nhp1_m, nhp2_m, nhp0_m, chan_m, chan0_m;
  int h_j, p_j, hh_j, pp_j, hpp_j, hhp_j, hp1_j, hp2_j;
  int h_m, p_m, hh_m, pp_m, hpp_m, hhp_m, hp1_m, hp2_m;
  double TJ;
  int ind, minj, jsize;

  int a2, b2, i2, j2;
  int k, l, c, d, k2, l2, c2, d2;

  double ja2, ma2, jb2, mb2, ji2, mi2, jj2, mj2, m12, m22;
  double tempt, tempen, term0;
  double *jvec, *tvec;
  double energy1 = 0.0;
  double energy2 = 0.0;
  double J, ij, jj, aj, bj, J1;
  int print = 0;
  int count = 0;
  int iterations = 10;

  Get_Input_Parameters(inputfile, Parameters2);
  Parameters2.basis = "finite_JM";

  Build_Model_Space(Parameters2, Space2);
  Print_Parameters(Parameters2, Space2);
  HF_Chan2 = HF_Channels(Parameters2, Space2);
  States2 = Single_Particle_States(Parameters2, Space2, HF_Chan2);
  if(Parameters2.basis == "finite_J" || Parameters2.basis == "finite_JM"){ Read_Matrix_Elements_J(Parameters2, Space2, HF_Chan2, HF_ME2); }
  else{ Read_Matrix_Elements_M(Parameters2, Space2, HF_Chan2, HF_ME2); }
  Hartree_Fock_States(Parameters2, Space2, HF_Chan2, States2, HF_ME2);
  Convert_To_HF_Matrix_Elements(HF_Chan2, States2, HF_ME2);
  States2.delete_struct(HF_Chan2);
  
  if(Parameters2.basis == "finite_JM"){
    Build_Model_Space_J2(Parameters2, Space2);
    Parameters2.basis = "finite_M";
    JM2 = 1;
  }
   
  Chan2 = Channels(Parameters2, Space2);
  Amps2 = Amplitudes(Parameters2, Space2, Chan2);
  Ints2 = Interactions(Parameters2, Chan2);
  Get_Matrix_Elements_JM(Parameters2, HF_Chan2, HF_ME2, Space2, Chan2, Ints2);
  HF_ME2.delete_struct(HF_Chan2);
  HF_Chan2.delete_struct();

  if(Parameters.approx == "singles"){
    chan0_j = Chan.ind0;
    nhp0_j = Chan.nhp1[chan0_j];
  }
  if(Parameters2.approx == "singles"){
    chan0_m = Chan2.ind0;
    nhp0_m = Chan2.nhp1[chan0_m];
  }
  
  Amplitudes tempAmps = Amplitudes(Parameters, Space, Chan);
  //Amplitudes tempAmpsJ = Amplitudes(Parameters, Space, Chan);
  Amplitudes tempAmps2 = Amplitudes(Parameters2, Space2, Chan2);
  //Amplitudes tempAmps2J = Amplitudes(Parameters2, Space2, Chan2);
  //Interactions tempInts2J = Interactions(Parameters2, Chan2);

  Amps.zero(Parameters, Chan);
  Amps2.zero(Parameters2, Chan2);

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh_j = Chan.nhh[chan];
    npp_j = Chan.npp[chan];
    if(nhh_j * npp_j == 0){ continue; }
    for(int pp = 0; pp < npp_j; ++pp){
      a2 = Chan.ppvec[chan][2*pp];
      b2 = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh_j; ++hh){
	i2 = Chan.hhvec[chan][2*hh];
	j2 = Chan.hhvec[chan][2*hh + 1];
	tempen = Space.qnums[i2].energy + Space.qnums[j2].energy - Space.qnums[a2].energy - Space.qnums[b2].energy;
	tempt = Ints.D_ME1.V4[chan][pp * nhh_j + hh] / tempen;
	//std::cout << "TJ: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > ^ " << 0.5*Chan.qnums1[chan].j << " = " << tempt << std::endl;
	ind = hh * npp_j + pp;
	Amps.D1.set_TJ(Space, Chan, chan, ind, i2, j2, a2, b2, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){
    for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
      i2 = Chan.hp1vec[chan0_j][2*hp1];
      a2 = Chan.hp1vec[chan0_j][2*hp1 + 1];
      tempen = Space.qnums[i2].energy - Space.qnums[a2].energy;
      tempt = 0.0;
      //std::cout << "tj: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
      Amps.S1.set_TJ(Space, hp1, i2, a2, tempt);
    }
    Amps.S1.set_T_2J(Space, Chan, Ints);
    //Amps.D1.set_T_2J(Chan, Ints);
  }
  for(int chan = 0; chan < Chan2.size1; ++chan){
    nhh_m = Chan2.nhh[chan];
    npp_m = Chan2.npp[chan];
    if(nhh_m * npp_m == 0){ continue; }
    for(int pp = 0; pp < npp_m; ++pp){
      a2 = Chan2.ppvec[chan][2*pp];
      b2 = Chan2.ppvec[chan][2*pp + 1];
      if(a2 == b2){ continue; }
      for(int hh = 0; hh < nhh_m; ++hh){
	i2 = Chan2.hhvec[chan][2*hh];
	j2 = Chan2.hhvec[chan][2*hh + 1];
	if(i2 == j2){ continue; }
	tempen = Space2.qnums[i2].energy + Space2.qnums[j2].energy - Space2.qnums[a2].energy - Space2.qnums[b2].energy;
	tempt = Ints2.D_ME1.V4[chan][pp * nhh_m + hh] / tempen;
	//std::cout << "TM: < " << a << "," << b << " |t| " << i << "," << j << " > = " << tempt << std::endl;
	ind = hh * npp_m + pp;
	Amps2.D1.set_T(chan, ind, tempt);
      }
    }
  }
  if(Parameters2.approx == "singles"){
    for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
      i2 = Chan2.hp1vec[chan0_m][2*hp1];
      a2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
      tempen = Space2.qnums[i2].energy - Space2.qnums[a2].energy;
      tempt = 0.0;
      //std::cout << "tm: " << i << " |t| " << a << " = " << tempt << std::endl;
      Amps2.S1.set_T(hp1, tempt);
    }
    Amps2.S1.set_T_2(Chan2, Ints2);
    //Amps2.D1.set_T_2(Chan2, Ints2);
  }
  energy1 = Amps.get_energy(Parameters, Chan, Ints);
  energy2 = Amps2.get_energy(Parameters2, Chan2, Ints2);
  std::cout << "!?!? " << energy1 << " " << energy2 << std::endl;
      
  while(count < iterations){
    ++count;
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	      tempt += Ints.D_ME1.V4[chan_j][pp_j * nhh_j + hh_j];
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term0 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		    tempt += Ints2.D_ME1.V4[chan_m][pp_m * nhh_m + hh_m];
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term0 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "?Term0 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/2 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
	      for(int pp = 0; pp < npp_j; ++pp){
		c = Chan.ppvec[chan_j][2*pp];
		d = Chan.ppvec[chan_j][2*pp + 1];
		term0 = 0.5 * Ints.D_ME1.V1[chan_j][pp * npp_j + pp_j] * Amps.D1.T1[chan_j][hh_j * npp_j + pp];
		if(print == 2){ std::cout << "Term1 += 1/2 < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = 1/2 * " << Ints.D_ME1.V1[chan_j][pp * npp_j + pp_j] << " * " << Amps.D1.T1[chan_j][hh_j * npp_j + pp] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term1 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/2 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
		    for(int pp = 0; pp < npp_m; ++pp){
		      c2 = Chan2.ppvec[chan_m][2*pp];
		      d2 = Chan2.ppvec[chan_m][2*pp + 1];
		      if(c2 == d2){ continue; }
		      term0 = 0.5 * Ints2.D_ME1.V1[chan_m][pp * npp_m + pp_m] * Amps2.D1.T1[chan_m][hh_m * npp_m + pp];
		      if(print == 2){ std::cout << "!Term1 += 1/2 < " << a2 << "," << b2 << " |V| " << c2 << "," << d2 << " > * < " << c2 << "," << d2 << " |t| " << i2 << "," << j2 << " > = 1/2 * " << Ints2.D_ME1.V1[chan_m][pp * npp_m + pp_m] << " * " << Amps2.D1.T1[chan_m][hh_m * npp_m + pp] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term1 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term1 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
    
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/2 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
	      for(int hh = 0; hh < nhh_j; ++hh){
		k = Chan.hhvec[chan_j][2*hh];
		l = Chan.hhvec[chan_j][2*hh + 1];
		term0 = 0.5 * Ints.D_ME1.V2[chan_j][hh_j * nhh_j + hh] * Amps.D1.T1[chan_j][hh * npp_j + pp_j];
		if(print == 2){ std::cout << "Term2 += 1/2 < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > ^ " << 0.5*tb1.j << " = 1/2 * " << Ints.D_ME1.V2[chan_j][hh_j * nhh_j + hh] << " * " << Amps.D1.T1[chan_j][hh * npp_j + pp_j] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term2 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/2 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
		    for(int hh = 0; hh < nhh_m; ++hh){
		      k2 = Chan2.hhvec[chan_m][2*hh];
		      l2 = Chan2.hhvec[chan_m][2*hh + 1];
		      if(k2 == l2){ continue; }
		      term0 = 0.5 * Ints2.D_ME1.V2[chan_m][hh_m * nhh_m + hh] * Amps2.D1.T1[chan_m][hh * npp_m + pp_m];
		      if(print == 2){ std::cout << "!Term2 += 1/2 < " << i2 << "," << j2 << " |V| " << k2 << "," << l2 << " > * < " << a2 << "," << b2 << " |t| " << k2 << "," << l2 << " > = 1/2 * " << Ints2.D_ME1.V2[chan_m][hh_m * nhh_m + hh] << " * " << Amps2.D1.T1[chan_m][hh * npp_m + pp_m] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term2 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term2 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
    
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/4 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
	      for(int pp = 0; pp < npp_j; ++pp){
		c = Chan.ppvec[chan_j][2*pp];
		d = Chan.ppvec[chan_j][2*pp + 1];
		for(int hh = 0; hh < nhh_j; ++hh){
		  k = Chan.hhvec[chan_j][2*hh];
		  l = Chan.hhvec[chan_j][2*hh + 1];
		  term0 = 0.25 * Ints.D_ME1.V4[chan_j][pp * nhh_j + hh] * Amps.D1.T1[chan_j][hh_j * npp_j + pp] * Amps.D1.T1[chan_j][hh * npp_j + pp_j];
		  if(print == 2){ std::cout << "Term7 += 1/4 < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 1/4 * " << Ints.D_ME1.V4[chan_j][pp * nhh_j + hh] << " * " << Amps.D1.T1[chan_j][hh_j * npp_j + pp] << " * " << Amps.D1.T1[chan_j][hh * npp_j + pp_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term7 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/4 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
		    for(int pp = 0; pp < npp_m; ++pp){
		      c2 = Chan2.ppvec[chan_m][2*pp];
		      d2 = Chan2.ppvec[chan_m][2*pp + 1];
		      if(c2 == d2){ continue; }
		      for(int hh = 0; hh < nhh_m; ++hh){
			k2 = Chan2.hhvec[chan_m][2*hh];
			l2 = Chan2.hhvec[chan_m][2*hh + 1];
			if(k2 == l2){ continue; }
			term0 = 0.25 * Ints2.D_ME1.V4[chan_m][pp * nhh_m + hh] * Amps2.D1.T1[chan_m][hh_m * npp_m + pp] * Amps2.D1.T1[chan_m][hh * npp_m + pp_m];
			if(print == 2){ std::cout << "!Term7 += 1/4 < " << k2 << "," << l2 << " |V| " << c2 << "," << d2 << " > * < " << c2 << "," << d2 << " |t| " << i2 << "," << j2 << " > * < " << a2 << "," << b2 << " |t| " << k2 << "," << l2 << " > = 1/4 * " << Ints2.D_ME1.V4[chan_m][pp * nhh_m + hh] << " * " << Amps2.D1.T1[chan_m][hh_m * npp_m + pp] << " * " << Amps2.D1.T1[chan_m][hh * npp_m + pp_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term7 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term7 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[a]);
	    minus(tb2, Space.qnums[b], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[a].j);
	    if(std::abs(Space.qnums[b].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[b].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, a, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, b, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3) J
	      //T3(ab|ij){jb,ia} = T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		term0 = Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2];
		if(print == 2){ std::cout << "Term3 += < " << c << "," << k << "^-1 |V| " << b << "," << j << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << a << "^-1 > = " << Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T2[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      tempAmps.D1.T3[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term3 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[a2]);
		    minus(tb4, Space2.qnums[b2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, a2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, b2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
		    //T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k2 = Chan2.hp2vec[chan_m][2*hp2];
		      c2 = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(a2 == c2 || i2 == k2){ continue; }
		      term0 = -1.0 * Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2];
		      if(print == 2){ std::cout << "!Term3 += -< " << c2 << "," << k2 << "^-1 |V| " << b2 << "," << j2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = -1 * " << Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    tempAmps2.D1.T3[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term3 = < " << b2 << "," << j2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + ma2;
		      m22 = mb2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, ja2 + ma2 + jj2 + mj2) * CGC(ji2, mi2, ja2, ma2, jvec[jint], m12) * CGC(jb2, mb2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term3 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
  
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[b]);
	    minus(tb2, Space.qnums[a], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[b].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, b, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, a, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5) J
	      //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		term0 = Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2];
		if(print == 2){ std::cout << "Term5 += < " << c << "," << k << "^-1 |V| " << a << "," << j << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << b << "^-1 > = " << Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T4[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      tempAmps.D1.T5[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term5 = < " << a << "," << j << "^-1 |t| " << i << "," << b << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	    
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[b2]);
		    minus(tb4, Space2.qnums[a2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, b2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, a2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
		    //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k2 = Chan2.hp2vec[chan_m][2*hp2];
		      c2 = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(b2 == c2 || i2 == k2){ continue; }
		      term0 = Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2];
		      if(print == 2){ std::cout << "!Term5 += < " << c2 << "," << k2 << "^-1 |V| " << a2 << "," << j2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| " << i2 << "," << b2 << "^-1 > = " << Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T4[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    tempAmps2.D1.T5[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term5 = < " << a2 << "," << j2 << "^-1 |t| " << i2 << "," << b2 << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mb2;
		      m22 = ma2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, jb2 + mb2 + jj2 + mj2) * CGC(ji2, mi2, jb2, mb2, jvec[jint], m12) * CGC(ja2, ma2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term5 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[a]);
	    minus(tb2, Space.qnums[b], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[a].j);
	    if(std::abs(Space.qnums[b].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[b].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, a, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, b, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp1_j; ++hp1){
		  l = Chan.hp1vec[chan_j][2*hp1];
		  d = Chan.hp1vec[chan_j][2*hp1 + 1];
		  term0 = Ints.D_ME1.V9[chan_j][hp2 * nhp1_j + hp1] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] * Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j];
		  if(print == 2){ std::cout << "Term12 += < " << l << "," << d << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << a << "^-1 > * < " << b << "," << j << "^-1 |t| " << l << "," << d << "^-1 > = " << Ints.D_ME1.V9[chan_j][hp2 * nhp1_j + hp1] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " * " << Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T2[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term12 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[a2]);
		    minus(tb4, Space2.qnums[b2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, a2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, b2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k = Chan2.hp2vec[chan_m][2*hp2];
		      c = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(k == i2 || c == a2){ continue; }
		      for(int hp1 = 0; hp1 < nhp1_m; ++hp1){
			l = Chan2.hp1vec[chan_m][2*hp1];
			d = Chan2.hp1vec[chan_m][2*hp1 + 1];
			if(l == j2 || d == b2){ continue; }
			if(k == l || c == d){ continue; }
			term0 = Ints2.D_ME1.V9[chan_m][hp2 * nhp1_m + hp1] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] * Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m];
			if(print == 2){ std::cout << "!Term12 += < " << l << "," << d << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| " << i2 << "," << a2 << "^-1 > * < " << b2 << "," << j2 << "^-1 |t| " << l << "," << d << "^-1 > = " << Ints2.D_ME1.V9[chan_m][hp2 * nhp1_m + hp1] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " * " << Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term12 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + ma2;
		      m22 = mb2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, ja2 + ma2 + jj2 + mj2) * CGC(ji2, mi2, ja2, ma2, jvec[jint], m12) * CGC(jb2, mb2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term12 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
  
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[b]);
	    minus(tb2, Space.qnums[a], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[b].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, b, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, a, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		d = Chan.hp2vec[chan_j][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp1_j; ++hp1){
		  l = Chan.hp1vec[chan_j][2*hp1];
		  c = Chan.hp1vec[chan_j][2*hp1 + 1];
		  term0 = Ints.D_ME1.V10[chan_j][hp2 * nhp1_j + hp1] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] * Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j];
		  if(print == 2){ std::cout << "Term13 += < " << l << "," << c << "^-1 |V| " << d << "," << k << "^-1 > * < " << d << "," << k << "^-1 |t| " << i << "," << b << "^-1 > * < " << a << "," << j << "^-1 |t| " << l << "," << c << "^-1 > = " << Ints.D_ME1.V10[chan_j][hp2 * nhp1_j + hp1] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " * " << Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T4[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[b2]);
		    minus(tb4, Space2.qnums[a2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, b2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, a2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k = Chan2.hp2vec[chan_m][2*hp2];
		      d = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(k == i2 || d == b2){ continue; }
		      for(int hp1 = 0; hp1 < nhp1_m; ++hp1){
			l = Chan2.hp1vec[chan_m][2*hp1];
			c = Chan2.hp1vec[chan_m][2*hp1 + 1];
			if(l == j2 || c == a2){ continue; }
			if(k == l || d == c){ continue; }
			term0 = Ints2.D_ME1.V10[chan_m][hp2 * nhp1_m + hp1] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] * Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m];
			if(print == 2){ std::cout << "!Term13 += < " << l << "," << c << "^-1 |V| " << d << "," << k << "^-1 > * < " << d << "," << k << "^-1 |t| " << i2 << "," << b2 << "^-1 > * < " << a2 << "," << j2 << "^-1 |t| " << l << "," << c << "^-1 > = " << Ints2.D_ME1.V10[chan_m][hp2 * nhp1_m + hp1] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " * " << Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T4[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term13 = < " << a << "," << j << "^-1 |t| " << i << "," << b << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mb2;
		      m22 = ma2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, jb2 + mb2 + jj2 + mj2) * CGC(ji2, mi2, jb2, mb2, jvec[jint], m12) * CGC(ja2, ma2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term13 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[i];
	    h_j = Chan.h_map[chan_j][i];
	    nh_j = Chan.nh[chan_j];
	    nhpp_j = Chan.nhpp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hpp_j = Chan.hpp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(j, a, b, Space.indtot)];
	      //T6(ab|ij){jab,i} = -1/2 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8) J
	      for(int h = 0; h < nh_j; ++h){
		l = Chan.hvec[chan_j][h];
		for(int hpp = 0; hpp < nhpp_j; ++hpp){
		  k = Chan.hppvec[chan_j][3*hpp];
		  c = Chan.hppvec[chan_j][3*hpp + 1];
		  d = Chan.hppvec[chan_j][3*hpp + 2];
		  term0 = -0.5 * Ints.D_ME1.V5[chan_j][h * nhpp_j + hpp] * Amps.D1.T7[chan_j][hpp_j * nh_j + h] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
		  if(print == 2){ std::cout << "Term8 += -1/2 < " << l << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i << " > * < " << a << "," << b << "," << j << "^-1 |t| " << l << "^-1 > = -1/2 * " << Ints.D_ME1.V5[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T7[chan_j][hpp_j * nh_j + h] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T6[chan_j][hpp_j * nh_j + h_j] += tempt;
	      if(print != 0){
		std::cout << "Term8 = < " << a << "," << b << "," << j << "^-1 |t| " << i << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[i2];
		    h_m = Chan2.h_map[chan_m][i2];
		    hpp_m = Chan2.hpp_map[chan_m][Hash3(j2, a2, b2, Space2.indtot)];
		    nh_m = Chan2.nh[chan_m];
		    nhpp_m = Chan2.nhpp[chan_m];
		    //T6(ab|ij){jab,i} = -1/2 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
		    for(int h = 0; h < nh_m; ++h){
		      l = Chan2.hvec[chan_m][h];
		      if(l == j2){ continue; }  
		      for(int hpp = 0; hpp < nhpp_m; ++hpp){
			k = Chan2.hppvec[chan_m][3*hpp];
			c = Chan2.hppvec[chan_m][3*hpp + 1];
			d = Chan2.hppvec[chan_m][3*hpp + 2];
			if(k == l || c == d || k == i2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V5[chan_m][h * nhpp_m + hpp] * Amps2.D1.T7[chan_m][hpp_m * nh_m + h] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
			if(print == 2){ std::cout << "!Term8 += -1/2 < " << l << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i2 << " > * < " << a2 << "," << b2 << "," << j2 << "^-1 |t| " << l << "^-1 > = -1/2 * " << Ints2.D_ME1.V5[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T7[chan_m][hpp_m * nh_m + h] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T6[chan_m][hpp_m * nh_m + h_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term8 = < " << a2 << "," << b2 << "," << j2 << "^-1 |t| " << i2 << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = ma2 + mb2;
		      m22 = m12 + mj2;
		      if(m22 != mi2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, jj2 + mj2) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m12) * CGC(jvec[jint], m12, jj2, mj2, ji2, mi2) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term8 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[j];
	    h_j = Chan.h_map[chan_j][j];
	    nh_j = Chan.nh[chan_j];
	    nhpp_j = Chan.nhpp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hpp_j = Chan.hpp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, a, b, Space.indtot)];
	      //T7(ab|ij){iab,j} = -1/2 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9) J
	      for(int h = 0; h < nh_j; ++h){
		k = Chan.hvec[chan_j][h];
		for(int hpp = 0; hpp < nhpp_j; ++hpp){
		  l = Chan.hppvec[chan_j][3*hpp];
		  c = Chan.hppvec[chan_j][3*hpp + 1];
		  d = Chan.hppvec[chan_j][3*hpp + 2];
		  term0 = -0.5 * Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] * Amps.D1.T7[chan_j][hpp_j * nh_j + h] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
		  if(print == 2){ std::cout << "Term9 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << j << " > * < " << a << "," << b << "," << i << "^-1 |t| " << k << "^-1 > = -1/2 * " << Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T7[chan_j][hpp_j * nh_j + h] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T7[chan_j][hpp_j * nh_j + h_j] += tempt;
	      if(print != 0){
		std::cout << "Term9 = < " << a << "," << b << "," << i << "^-1 |t| " << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = -0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[j2];
		    h_m = Chan2.h_map[chan_m][j2];
		    hpp_m = Chan2.hpp_map[chan_m][Hash3(i2, a2, b2, Space2.indtot)];
		    nh_m = Chan2.nh[chan_m];
		    nhpp_m = Chan2.nhpp[chan_m];
		    //T7(ab|ij){iab,j} = -1/2 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
		    for(int h = 0; h < nh_m; ++h){
		      k = Chan2.hvec[chan_m][h];
		      if(k == i2){ continue; }  
		      for(int hpp = 0; hpp < nhpp_m; ++hpp){
			l = Chan2.hppvec[chan_m][3*hpp];
			c = Chan2.hppvec[chan_m][3*hpp + 1];
			d = Chan2.hppvec[chan_m][3*hpp + 2];
			if(l == k || c == d || l == j2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] * Amps2.D1.T7[chan_m][hpp_m * nh_m + h] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
			if(print == 2){ std::cout << "!Term9 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << j2 << " > * < " << a2 << "," << b2 << "," << i2 << "^-1 |t| " << k << "^-1 > = -1/2 * " << Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T7[chan_m][hpp_m * nh_m + h] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T7[chan_m][hpp_m * nh_m + h_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term9 = < " << a2 << "," << b2 << "," << i2 << "^-1 |t| " << j2 << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = ma2 + mb2;
		      m22 = m12 + mi2;
		      if(m22 != mj2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, ji2 + mi2) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m12) * CGC(jvec[jint], m12, ji2, mi2, jj2, mj2) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term9 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[a];
	    p_j = Chan.p_map[chan_j][a];
	    np_j = Chan.np[chan_j];
	    nhhp_j = Chan.nhhp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hhp_j = Chan.hhp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, j, b, Space.indtot)];
	      //T8(ab|ij){ijb,a} = -1/2 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10) J
	      for(int p = 0; p < np_j; ++p){
		d = Chan.pvec[chan_j][p];
		for(int hhp = 0; hhp < nhhp_j; ++hhp){
		  k = Chan.hhpvec[chan_j][3*hhp];
		  l = Chan.hhpvec[chan_j][3*hhp + 1];
		  c = Chan.hhpvec[chan_j][3*hhp + 2];
		  term0 = -0.5 * Ints.D_ME1.V7[chan_j][p * nhhp_j + hhp] * Amps.D1.T9[chan_j][hhp_j * np_j + p] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
		  if(print == 2){ std::cout << "Term10 += -1/2 < " << k << "," << l << "," << c << "^-1 |V| " << d << " > * < " << d << " |t| " << i << "," << j << "," << b << "^-1 > * < " << a << " |t| " << k << "," << l << "," << c << "^-1 > = -1/2 * " << Ints.D_ME1.V7[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T9[chan_j][hhp_j * np_j + p] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T8[chan_j][hhp_j * np_j + p_j] += tempt;
	      if(print != 0){
		std::cout << "Term10 = < " << a << " |t| " << i << "," << j << "," << b << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[a2];
		    p_m = Chan2.p_map[chan_m][a2];
		    hhp_m = Chan2.hhp_map[chan_m][Hash3(i2, j2, b2, Space2.indtot)];
		    np_m = Chan2.np[chan_m];
		    nhhp_m = Chan2.nhhp[chan_m];
		    //T8(ab|ij){ijb,a} = -1/2 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
		    for(int p = 0; p < np_m; ++p){
		      d = Chan2.pvec[chan_m][p];
		      if(d == b2){ continue; }  
		      for(int hhp = 0; hhp < nhhp_m; ++hhp){
			k = Chan2.hhpvec[chan_m][3*hhp];
			l = Chan2.hhpvec[chan_m][3*hhp + 1];
			c = Chan2.hhpvec[chan_m][3*hhp + 2];
			if(k == l || c == d || c == a2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V7[chan_m][p * nhhp_m + hhp] * Amps2.D1.T9[chan_m][hhp_m * np_m + p] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
			if(print == 2){ std::cout << "!Term10 += -1/2 < " << k << "," << l << "," << c << "^-1 |V| " << d << " > * < " << d << " |t| " << i2 << "," << j2 << "," << b2 << "^-1 > * < " << a2 << " |t| " << k << "," << l << "," << c << "^-1 > = -1/2 * " << Ints2.D_ME1.V7[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T9[chan_m][hhp_m * np_m + p] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T8[chan_m][hhp_m * np_m + p_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term10 = < " << a2 << " |t| " << i2 << "," << j2 << "," << b2 << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = m12 + mb2;
		      if(m22 != ma2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, jb2 + mb2) * CGC(jvec[jint], m12, jb2, mb2, ja2, ma2) * CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term10 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[b];
	    p_j = Chan.p_map[chan_j][b];
	    np_j = Chan.np[chan_j];
	    nhhp_j = Chan.nhhp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hhp_j = Chan.hhp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, j, a, Space.indtot)];
	      //T9(ab|ij){ija,b} = -1/2 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11) J
	      for(int p = 0; p < np_j; ++p){
		c = Chan.pvec[chan_j][p];
		for(int hhp = 0; hhp < nhhp_j; ++hhp){
		  k = Chan.hhpvec[chan_j][3*hhp];
		  l = Chan.hhpvec[chan_j][3*hhp + 1];
		  d = Chan.hhpvec[chan_j][3*hhp + 2];
		  term0 = -0.5 * Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.D1.T9[chan_j][hhp_j * np_j + p] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
		  if(print == 2){ std::cout << "Term11 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << c << " |t| " << i << "," << j << "," << a << "^-1 > * < " << b << " |t| " << k << "," << l << "," << d << "^-1 > = -1/2 * " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T9[chan_j][hhp_j * np_j + p] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T9[chan_j][hhp_j * np_j + p_j] += tempt;
	      if(print != 0){
		std::cout << "Term11 = < " << b << " |t| " << i << "," << j << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[b2];
		    p_m = Chan2.p_map[chan_m][b2];
		    hhp_m = Chan2.hhp_map[chan_m][Hash3(i2, j2, a2, Space2.indtot)];
		    np_m = Chan2.np[chan_m];
		    nhhp_m = Chan2.nhhp[chan_m];
		    //T9(ab|ij){ija,b} = -1/2 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
		    for(int p = 0; p < np_m; ++p){
		      c = Chan2.pvec[chan_m][p];
		      if(c == a2){ continue; }  
		      for(int hhp = 0; hhp < nhhp_m; ++hhp){
			k = Chan2.hhpvec[chan_m][3*hhp];
			l = Chan2.hhpvec[chan_m][3*hhp + 1];
			d = Chan2.hhpvec[chan_m][3*hhp + 2];
			if(k == l || c == d || d == b2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.D1.T9[chan_m][hhp_m * np_m + p] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
			if(print == 2){ std::cout << "!Term11 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << c << " |t| " << i2 << "," << j2 << "," << a2 << "^-1 > * < " << b2 << " |t| " << k << "," << l << "," << d << "^-1 > = -1/2 * " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T9[chan_m][hhp_m * np_m + p] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T9[chan_m][hhp_m * np_m + p_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term11 = < " << b2 << " |t| " << i2 << "," << j2 << "," << a2 << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = m12 + ma2;
		      if(m22 != mb2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, ja2 + ma2) * CGC(jvec[jint], m12, ja2, ma2, jb2, mb2) * CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term11 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    if(Parameters.approx == "singles" && Parameters2.approx == "singles"){
      //print = 1;
      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  hp1_j = Chan.hp1_map[chan0_j][Hash2(i, a, Space.indtot)];
	  //t1(ia){ia} = V3(ka|ic){ia,kc}.t1(kc){kc} (1) J
	  for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	    k = Chan.hp1vec[chan0_j][2*hp1];
	    c = Chan.hp1vec[chan0_j][2*hp1 + 1];
	    term0 = Ints.D_ME1.V3[chan0_j][hp1_j * nhp0_j + hp1] * Amps.S1.T1[hp1];
	    if(print == 2){ std::cout << "Term_1 += < " << a << "," << i << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| > = " << Ints.D_ME1.V3[chan0_j][hp1_j * nhp0_j + hp1] << " * " << Amps.S1.T1[hp1] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T1[hp1_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_1 = < " << a << "," << i << "^-1 |t| > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
	      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = -0.5*Space2.qnums[a2].m;
	      hp1_m = Chan2.hp1_map[chan0_m][Hash2(i2, a2, Space2.indtot)];
	      //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
	      for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
		k2 = Chan2.hp1vec[chan0_m][2*hp1];
		c2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
		term0 = -1.0 * Ints2.D_ME1.V3[chan0_m][hp1_m * nhp0_m + hp1] * Amps2.S1.T1[hp1];
		if(print == 2){ std::cout << "!Term_1 += -< " << a2 << "," << i2 << "^-1 |V| " << c2 << "," << k2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| > = -1 * " << Ints2.D_ME1.V3[chan0_m][hp1_m * nhp0_m + hp1] << " * " << Amps2.S1.T1[hp1] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T1[hp1_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_1 = < " << a2 << "," << i2 << "^-1 |t| > = " << tempt << std::endl;
	      }
	      double tempt2 = std::pow(-1.0, ja2 + ma2) * CGC(ji2, mi2, ja2, ma2, 0, 0) * TJ;
	      if(print != 0){ std::cout << "? Term_1 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  hp1_j = Chan.hp1_map[chan0_j][Hash2(i, a, Space.indtot)];
	  //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6) J
	  for(int hp2 = 0; hp2 < nhp0_j; ++hp2){
	    k = Chan.hp2vec[chan0_j][2*hp2];
	    c = Chan.hp2vec[chan0_j][2*hp2 + 1];
	    for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	      l = Chan.hp1vec[chan0_j][2*hp1];
	      d = Chan.hp1vec[chan0_j][2*hp1 + 1];
	      term0 = Ints.D_ME1.V9[chan0_j][hp2 * nhp0_j + hp1] * Amps.D1.T2[chan0_j][hp1 * nhp0_j + hp1_j] * Amps.S1.T1[hp2];
	      if(print == 2){ std::cout << "Term_6 += < " << k << "," << c << "^-1 |V| " << d << "," << l << "^-1 > * < " << c << "," << k << "^-1 |t| > * < " << d << "," << l << "^-1 |t| " << i << "," << a << "^-1 > = " << Ints.D_ME1.V9[chan0_j][hp2 * nhp0_j + hp1] << " * " << Amps.D1.T2[chan0_j][hp1 * nhp0_j + hp1_j] << " * " << Amps.S1.T1[hp2] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T1[hp1_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_6 = < " << a << "," << i << "^-1 |t| > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = -0.5*Space2.qnums[a2].m;
	      hp1_m = Chan2.hp1_map[chan0_m][Hash2(i2, a2, Space2.indtot)];
	      //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
	      for(int hp2 = 0; hp2 < nhp0_m; ++hp2){
		k2 = Chan2.hp2vec[chan0_m][2*hp2];
		c2 = Chan2.hp2vec[chan0_m][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
		  l2 = Chan2.hp1vec[chan0_m][2*hp1];
		  d2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
		  term0 = Ints2.D_ME1.V9[chan0_m][hp2 * nhp0_m + hp1] * Amps2.D1.T2[chan0_m][hp1 * nhp0_m + hp1_m] * Amps2.S1.T1[hp2];
		  if(print == 2){ std::cout << "!Term_6 += < " << k2 << "," << c2 << "^-1 |V| " << d2 << "," << l2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| > * < " << d2 << "," << l2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = " << Ints2.D_ME1.V9[chan0_m][hp2 * nhp0_m + hp1] << " * " << Amps2.D1.T2[chan0_m][hp1 * nhp0_m + hp1_m] << " * " << Amps2.S1.T1[hp2] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T1[hp1_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_6 = < " << a2 << "," << i2 << "^-1 |t| > = " << tempt << std::endl;
	      }
	      double tempt2 = std::pow(-1.0, ja2 + ma2) * CGC(ji2, mi2, ja2, ma2, 0, 0) * TJ;
	      if(print != 0){ std::cout << "? Term_6 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = 1/2 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2) J
	  for(int hpp = 0; hpp < nhpp_j; ++hpp){
	    k = Chan.hppvec[chan_j][3*hpp];
	    c = Chan.hppvec[chan_j][3*hpp + 1];
	    d = Chan.hppvec[chan_j][3*hpp + 2];
	    term0 = 0.5 * Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
	    if(print == 2){ std::cout << "Term_2 += 1/2 < " << a << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i << " > = 1/2 * " << Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_2 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = -1/2 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
	      for(int hpp = 0; hpp < nhpp_m; ++hpp){
		k2 = Chan2.hppvec[chan_m][2*hpp];
		c2 = Chan2.hppvec[chan_m][2*hpp + 1];
		d2 = Chan2.hppvec[chan_m][2*hpp + 2];
		term0 = -0.5 * Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
		if(print == 2){ std::cout << "!Term_2 += -1/2 < " << a2 << " |V| " << c2 << "," << d2 << "," << k2 << "^-1 > * < " << c2 << "," << d2 << "," << k2 << "^-1 |t| " << i2 << " > = -1/2 * " << Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_2 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_2 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = -V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3) J
	  for(int hpp = 0; hpp < nhpp_j; ++hpp){
	    k = Chan.hppvec[chan_j][3*hpp];
	    c = Chan.hppvec[chan_j][3*hpp + 1];
	    d = Chan.hppvec[chan_j][3*hpp + 2];
	    term0 = -1.0 * Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] * Amps.S1.E6[chan_j][hpp * nh_j + h_j];
	    if(print == 2){ std::cout << "Term_3 += -< " << a << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |E| " << i << " > = -1 * " << Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] << " * " << Amps.S1.E6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_3 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3)
	      for(int hpp = 0; hpp < nhpp_m; ++hpp){
		k2 = Chan2.hppvec[chan_m][3*hpp];
		c2 = Chan2.hppvec[chan_m][3*hpp + 1];
		d2 = Chan2.hppvec[chan_m][3*hpp + 2];
		term0 = Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] * Amps2.S1.E6[chan_m][hpp * nh_m + h_m];
		if(print == 2){ std::cout << "!Term_3 += < " << a2 << " |V| " << c2 << "," << d2 << "," << k2 << "^-1 > * < " << c2 << "," << d2 << "," << k2 << "^-1 |E| " << i2 << " > = -1/2 * " << Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] << " * " << Amps2.S1.E6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_3 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_3 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = -1/2 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7) J
	  for(int h = 0; h < nh_j; ++h){
	    k = Chan.hvec[chan_j][h];
	    for(int hpp = 0; hpp < nhpp_j; ++hpp){
	      l = Chan.hppvec[chan_j][3*hpp];
	      c = Chan.hppvec[chan_j][3*hpp + 1];
	      d = Chan.hppvec[chan_j][3*hpp + 2];
	      term0 = -0.5 * Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] * Amps.D1.T6[chan_j][hpp * nh_j + h_j] * Amps.S1.T2[chan_j][p_j * nh_j + h];
	      if(print == 2){ std::cout << "Term_7 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << i << " > * < " << a << " |t| " << k << " > = -1/2 * " << Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " * " << Amps.S1.T2[chan_j][p_j * nh_j + h] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_7 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = -1/2 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
	      for(int h = 0; h < nh_m; ++h){
		k2 = Chan2.hvec[chan_m][h];
		for(int hpp = 0; hpp < nhpp_m; ++hpp){
		  l2 = Chan2.hppvec[chan_m][3*hpp];
		  c2 = Chan2.hppvec[chan_m][3*hpp + 1];
		  d2 = Chan2.hppvec[chan_m][3*hpp + 2];
		  term0 = -0.5 * Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m] * Amps2.S1.T2[chan_m][p_m * nh_m + h];
		  if(print == 2){ std::cout << "!Term_7 += -1/2 < " << k2 << " |V| " << c2 << "," << d2 << "," << l2 << "^-1 > * < " << c2 << "," << d2 << "," << l2 << "^-1 |t| " << i2 << " > * < " << a2 << " |t| " << k2 << " > = -1/2 * " << Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " * " << Amps2.S1.T2[chan_m][p_m * nh_m + h] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_7 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_7 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = -1/2 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4) J
	  for(int hhp = 0; hhp < nhhp_j; ++hhp){
	    k = Chan.hhpvec[chan_j][3*hhp];
	    l = Chan.hhpvec[chan_j][3*hhp + 1];
	    c = Chan.hhpvec[chan_j][3*hhp + 2];
	    term0 = -0.5 * Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
	    if(print == 2){ std::cout << "Term_4 += -1/2 < " << i << " |V| " << k << "," << l << "," << c << "^-1 > * < " << k << "," << l << "," << c << "^-1 |t| " << a << " > = -1/2 * " << Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_4 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = -1/2 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
	      for(int hhp = 0; hhp < nhhp_m; ++hhp){
		k2 = Chan2.hhpvec[chan_m][3*hhp];
		l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		c2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		term0 = -0.5 * Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
		if(print == 2){ std::cout << "!Term_4 += -1/2 < " << i2 << " |V| " << k2 << "," << l2 << "," << c2 << "^-1 > * < " << k2 << "," << l2 << "," << c2 << "^-1 |t| " << a2 << " > = -1/2 * " << Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_4 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_4 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5) J
	  for(int hhp = 0; hhp < nhhp_j; ++hhp){
	    k = Chan.hhpvec[chan_j][3*hhp];
	    l = Chan.hhpvec[chan_j][3*hhp + 1];
	    c = Chan.hhpvec[chan_j][3*hhp + 2];
	    term0 = Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] * Amps.S1.E8[chan_j][hhp * np_j + p_j];
	    if(print == 2){ std::cout << "Term_5 += < " << i << " |V| " << k << "," << l << "," << c << "^-1 > * < " << k << "," << l << "," << c << "^-1 |E| " << a << " > = " << Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] << " * " << Amps.S1.E8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_5 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5)
	      for(int hhp = 0; hhp < nhhp_m; ++hhp){
		k2 = Chan2.hhpvec[chan_m][3*hhp];
		l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		c2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		term0 = Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] * Amps2.S1.E8[chan_m][hhp * np_m + p_m];
		if(print == 2){ std::cout << "!Term_5 += < " << i2 << " |V| " << k2 << "," << l2 << "," << c2 << "^-1 > * < " << k2 << "," << l2 << "," << c2 << "^-1 |E| " << a2 << " > = " << Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] << " * " << Amps2.S1.E8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_5 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_5 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = -1/2 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8) J
	  for(int p = 0; p < np_j; ++p){
	    c = Chan.pvec[chan_j][p];
	    for(int hhp = 0; hhp < nhhp_j; ++hhp){
	      k = Chan.hhpvec[chan_j][3*hhp];
	      l = Chan.hhpvec[chan_j][3*hhp + 1];
	      d = Chan.hhpvec[chan_j][3*hhp + 2];
	      term0 = -0.5 * Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.D1.T8[chan_j][hhp * np_j + p_j] * Amps.S1.T3[chan_j][h_j * np_j + p];
	      if(print == 2){ std::cout << "Term_8 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << a << " |t| " << k << "," << l << "," << d << "^-1 > * < " << c << " |t| " << i << " > = -1/2 * " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " * " << Amps.S1.T3[chan_j][h_j * np_j + p] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_8 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = -1/2 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
	      for(int p = 0; p < np_m; ++p){
		c2 = Chan2.pvec[chan_m][p];
		for(int hhp = 0; hhp < nhhp_m; ++hhp){
		  k2 = Chan2.hhpvec[chan_m][3*hhp];
		  l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		  d2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		  term0 = -0.5 * Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.D1.T8[chan_m][hhp * np_m + p_m] * Amps2.S1.T3[chan_m][h_m * np_m + p];
		  if(print == 2){ std::cout << "!Term_8 += -1/2 < " << k2 << "," << l2 << "," << d2 << "^-1 |V| " << c2 << " > * < " << a2 << " |t| " << k2 << "," << l2 << "," << d2 << "^-1 > * < " << c2 << " |t| " << i2 << " > = -1/2 * " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " * " << Amps2.S1.T3[chan_m][h_m * np_m + p] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_8 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_8 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9) J
	  for(int p = 0; p < np_j; ++p){
	    c = Chan.pvec[chan_j][p];
	    for(int hhp = 0; hhp < nhhp_j; ++hhp){
	      k = Chan.hhpvec[chan_j][3*hhp];
	      l = Chan.hhpvec[chan_j][3*hhp + 1];
	      d = Chan.hhpvec[chan_j][3*hhp + 2];
	      term0 = Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.S1.E8[chan_j][hhp * np_j+ p_j] * Amps.S1.T3[chan_j][h_j * np_j + p];
	      if(print == 2){ std::cout << "Term_9 += < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << a << " |E| " << k << "," << l << "," << d << "^-1 > * < " << c << " |t| " << i << " > = " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.S1.E8[chan_j][hhp * np_j + p_j] << " * " << Amps.S1.T3[chan_j][h_j * np_j + p] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_9 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
	   	      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9)
	      for(int p = 0; p < np_m; ++p){
		c2 = Chan2.pvec[chan_m][p];
		for(int hhp = 0; hhp < nhhp_m; ++hhp){
		  k2 = Chan2.hhpvec[chan_m][3*hhp];
		  l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		  d2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		  term0 = Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.S1.E8[chan_m][hhp * np_m + p_m] * Amps2.S1.T3[chan_m][h_m * np_m + p];
		  if(print == 2){ std::cout << "!Term_9 += < " << k2 << "," << l2 << "," << d2 << "^-1 |V| " << c2 << " > * < " << a2 << " |E| " << k2 << "," << l2 << "," << d2 << "^-1 > * < " << c2 << " |t| " << i2 << " > = " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.S1.E8[chan_m][hhp * np_m + p_m] << " * " << Amps2.S1.T3[chan_m][h_m * np_m + p] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_9 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_9 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }
    }
    //print = 0;

    Amps.zero(Parameters, Chan);
    Amps2.zero(Parameters2, Chan2);
	
    for(int chan = 0; chan < Chan.size1; ++chan){
      nhh_j = Chan.nhh[chan];
      npp_j = Chan.npp[chan];
      if(nhh_j * npp_j == 0){ continue; }
      for(int pp = 0; pp < npp_j; ++pp){
	a2 = Chan.ppvec[chan][2*pp];
	b2 = Chan.ppvec[chan][2*pp + 1];
	for(int hh = 0; hh < nhh_j; ++hh){
	  i2 = Chan.hhvec[chan][2*hh];
	  j2 = Chan.hhvec[chan][2*hh + 1];
	  ind = hh * npp_j + pp;
	  //tempt = tempAmps.D1.T1[chan][ind];
	  tempt = tempAmps.D1.get_TJ(Space, Chan, chan, ind, i2, j2, a2, b2);
	  tempen = Space.qnums[i2].energy + Space.qnums[j2].energy - Space.qnums[a2].energy - Space.qnums[b2].energy;
	  tempt /= tempen;
	  //std::cout << "TJ: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > ^ " << 0.5*Chan.qnums1[chan].j << " = " << tempt << std::endl;
	  Amps.D1.set_TJ(Space, Chan, chan, ind, i2, j2, a2, b2, tempt);
	}
      }
    }
    if(Parameters.approx == "singles"){
      for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	i2 = Chan.hp1vec[chan0_j][2*hp1];
	a2 = Chan.hp1vec[chan0_j][2*hp1 + 1];
	//tempt = tempAmps.S1.T1[hp1];
	tempt = tempAmps.S1.get_TJ(Space, hp1, i2, a2);
	tempen = Space.qnums[i2].energy - Space.qnums[a2].energy;
	tempt /= tempen;
	//std::cout << "tj: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
	Amps.S1.set_TJ(Space, hp1, i2, a2, tempt);
      }
      Amps.S1.set_T_2J(Space, Chan, Ints);
      //Amps.D1.set_T_2J(Chan, Ints);
    }
    for(int chan = 0; chan < Chan2.size1; ++chan){
      nhh_m = Chan2.nhh[chan];
      npp_m = Chan2.npp[chan];
      if(nhh_m * npp_m == 0){ continue; }
      for(int pp = 0; pp < npp_m; ++pp){
	a2 = Chan2.ppvec[chan][2*pp];
	b2 = Chan2.ppvec[chan][2*pp + 1];
	if(a2 == b2){ continue; }
	for(int hh = 0; hh < nhh_m; ++hh){
	  i2 = Chan2.hhvec[chan][2*hh];
	  j2 = Chan2.hhvec[chan][2*hh + 1];
	  if(i2 == j2){ continue; }
	  ind = hh * npp_m + pp;
	  //tempt = tempAmps2.D1.T1[chan][ind];
	  tempt = tempAmps2.D1.get_T(chan, ind);
	  tempen = Space2.qnums[i2].energy + Space2.qnums[j2].energy - Space2.qnums[a2].energy - Space2.qnums[b2].energy;
	  tempt /= tempen;
	  //std::cout << "TM: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
	  Amps2.D1.set_T(chan, ind, tempt);
	}
      }
    }
    if(Parameters2.approx == "singles"){
      for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
	i2 = Chan2.hp1vec[chan0_m][2*hp1];
	a2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
	//tempt = tempAmps2.S1.T1[hp1];
	tempt = tempAmps2.S1.get_T(hp1);
	tempen = Space2.qnums[i2].energy - Space2.qnums[a2].energy;
	tempt /= tempen;
	//std::cout << "tm: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
	Amps2.S1.set_T(hp1, tempt);
      }
      Amps2.S1.set_T_2(Chan2, Ints2);
      //Amps2.D1.set_T_2(Chan2, Ints2);
    }
    energy1 = Amps.get_energy(Parameters, Chan, Ints);
    energy2 = Amps2.get_energy(Parameters2, Chan2, Ints2);
    std::cout << "!?!? " << energy1 << " " << energy2 << std::endl;
	
    tempAmps.zero(Parameters, Chan);
    tempAmps2.zero(Parameters2, Chan2);
  }
  tempAmps.delete_struct(Parameters, Chan);
  tempAmps2.delete_struct(Parameters2, Chan2);
  Ints2.delete_struct(Parameters2, Chan2);
  Amps2.delete_struct(Parameters2, Chan2);
  Chan2.delete_struct();
  Space2.delete_struct(Parameters2);
}
