#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

//   Function to search for index1 of p in vec1 (size num1) and index2 of q in vec2 (size num2)
//   Returns index of pq in matrix <p|q>. Index = index1*num2 + index2.
int Index11(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index11 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[i] == q){ ind2 = i; break; }
    else if(i == num2 - 1){
      for(int j = 0; j < num2; ++j){ std::cout << vec2[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index11 for " << q << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;
}

//   Function to search for index of pq in vec1 (size num1)
//   Returns index of pq in matrix <p|q>.
int Index2(const std::vector<int> &vec1, const int &p, const int &q)
{
  int ind1;
  for(int i = 0; i < int(0.5 * vec1.size()); ++i){
    if(vec1[2*i] == p && vec1[2*i + 1] == q){ ind1 = i; break; }
    else if(i == int(0.5 * vec1.size()) - 1){ return -1; }
  }
  return ind1;
}


//   Returns index of p in vector <p>.
int Index1(const std::vector<int> &vec1, const int &p)
{
  int ind1;
  for(int i = 0; i < int(vec1.size()); ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == int(vec1.size() - 1)){
      for(int j = 0; j < int(vec1.size()); ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index1 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1;
}

//   Function to search for index1 of pq in vec1 (size num1) and index2 of rs in vec2 (size num2)
//   Returns index of pqrs in matrix <pq|rs>. Index = index1*num2 + index2.
int Index22(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[2*i] == p && vec1[2*i + 1] == q){
      ind1 = i;
      break;
    }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[2*j] << "," << vec1[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index22 for " << p << " " << q << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[2*i] == r && vec2[2*i + 1] == s){
      ind2 = i;
      break;
    }
    else if(i == num2 - 1){
      for(int j = 0; j < num2; ++j){ std::cout << vec2[2*j] << "," << vec2[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index22 for " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

//   Function to search for index1 of p in vec1 (size num1) and index2 of qrs in vec2 (size num2)
//   Returns index of pqrs in matrix <p|qrs>. Index = index1*num2 + index2.
int Index13(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){ 
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index13 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[3*i] == q && vec2[3*i + 1] == r && vec2[3*i + 2] == s){ ind2 = i; break; }
    else if(i == num2 - 1){ 
      for(int j = 0; j < num2; ++j){ std::cout << vec2[3*j] << "," << vec2[3*j + 1] << "," << vec2[3*j + 2] << " "; }
      std::cout << std::endl;
      std::cerr << "Index13 for " << q << " " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

//   Function to search for index1 of pqr in vec1 (size num1) and index2 of s in vec2 (size num2)
//   Returns index of pqrs in matrix <pqr|s>. Index = index1*num2 + index2.
int Index31(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[3*i] == p && vec1[3*i + 1] == q && vec1[3*i + 2] == r){ ind1 = i; break; }
    else if(i == num1 - 1){ 
      for(int j = 0; j < num1; ++j){ std::cout << vec1[3*j] << "," << vec1[3*j + 1] << "," << vec1[3*j + 2] << " "; }
      std::cout << std::endl;
      std::cerr << "Index31 for " << p << " " << q << " " << r << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[i] == s){ ind2 = i; break; }
    else if(i == num2 - 1){ 
      for(int j = 0; j < num2; ++j){ std::cout << vec2[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index31 for " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

void plus(State &S, const State &S1, const State &S2){
  S.t = S1.t + S2.t;
  S.m = S1.m + S2.m;
  S.nx = S1.nx + S2.nx;
  S.ny = S1.ny + S2.ny;
  S.nz = S1.nz + S2.nz;
  S.ml = S1.ml + S2.ml;
  S.par = S1.par * S2.par;
}

void minus(State &S, const State &S1, const State &S2){
  S.t = S1.t - S2.t;
  S.m = S1.m - S2.m;
  S.nx = S1.nx - S2.nx;
  S.ny = S1.ny - S2.ny;
  S.nz = S1.nz - S2.nz;
  S.ml = S1.ml - S2.ml;
  S.par = S1.par * S2.par;
}

double Amplitudes::get_energy(const Input_Parameters &Parameters, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	energy += 0.25 * D1.T1[chan][hind * Chan.pp[chan] + pind] * Ints.D_ME1.V4[chan][pind * Chan.hh[chan] + hind];
      }
    }
  }
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
      for(int j = 0; j < Chan.hp1[Chan.ind0]; ++j){
	energy += 0.5 * S1.T1[i] * S1.T1[j] * Ints.D_ME1.V9[Chan.ind0][i*Chan.hp1[Chan.ind0] + j];
      }
    }
  }
  return energy;
}

void Setup_Amps(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps)
{
  if(Parameters.approx == "doubles"){
    Amps.D1 = Doubles_1(Chan, Parameters, Space);
  }
  else if(Parameters.approx == "singles"){
    Amps.D1 = Doubles_1(Chan, Parameters, Space);
    Amps.S1 = Singles_1(Chan, Parameters, Space);
  }
  else if(Parameters.approx == "triples"){
    Amps.D1 = Doubles_1(Chan, Parameters, Space);
    Amps.S1 = Singles_1(Chan, Parameters, Space);
  }
}

void Amplitudes::zero(const Channels &Chan, const Input_Parameters &Parameters)
{
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    D1.T1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    D1.S1[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    D1.Q11[i].assign(Chan.hp[i] * Chan.pp[i], 0.0);
    D1.Q21[i].assign(Chan.hh[i] * Chan.hp[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    D1.T6[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    D1.T7[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    D1.T8[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
    D1.T9[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
    D1.S2[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    D1.S3[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    D1.S4[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    D1.S5[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    D1.Q12[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    D1.Q22[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    D1.T2[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    D1.T3[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    D1.T4[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    D1.T5[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    D1.S6[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    D1.S7[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
  }
  if(Parameters.approx == "singles"){
    S1.T1.assign(Chan.hp1[Chan.ind0], 0.0);
    S1.S3.assign(Chan.hp1[Chan.ind0] * Chan.hp1[Chan.ind0], 0.0);
    S1.S4.assign(Chan.hp1[Chan.ind0], 0.0);

    S1.Q31.assign(Chan.hh1[Chan.ind0], 0.0);
    S1.Q41.assign(Chan.pp1[Chan.ind0], 0.0);
    for(int i = 0; i < Chan.size1; ++i){
      S1.E1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
      S1.Q51[i].assign(Chan.hp[i] * Chan.pp[i], 0.0);
      S1.Q61[i].assign(Chan.hh[i] * Chan.hp[i], 0.0);
    }
    for(int i = 0; i < Chan.size3; ++i){
      S1.T2[i].assign(Chan.p[i] * Chan.h[i], 0.0);
      S1.T3[i].assign(Chan.h[i] * Chan.p[i], 0.0);
      S1.S1[i].assign(Chan.p[i] * Chan.p[i], 0.0);
      S1.S2[i].assign(Chan.h[i] * Chan.h[i], 0.0);
      S1.E6[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
      S1.E7[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
      S1.E8[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
      S1.E9[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);

      S1.Q11[i].assign(Chan.hpp2[i] * Chan.h[i], 0.0);
      S1.Q21[i].assign(Chan.hhp2[i] * Chan.p[i], 0.0);

      S1.Q32[i].assign(Chan.h[i] * Chan.h[i], 0.0);
      S1.Q42[i].assign(Chan.p[i] * Chan.p[i], 0.0);

      S1.Q52[i].assign(Chan.hpp[i] * Chan.p[i], 0.0);
      S1.Q62[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
    }
    for(int i = 0; i < Chan.size2; ++i){
      S1.E2[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      S1.E3[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      S1.E4[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      S1.E5[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);

      S1.Q12[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
      S1.Q22[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
    }
  }
}

Doubles_1::Doubles_1()
{
  Tmap.resize(0);
  Evec.resize(0);
  T1.resize(0);
  T2.resize(0);
  T3.resize(0);
  T4.resize(0);
  T5.resize(0);
  T6.resize(0);
  T7.resize(0);
  T8.resize(0);
  T9.resize(0);
  S1.resize(0);
  S2.resize(0);
  S3.resize(0);
  S4.resize(0);
  S5.resize(0);
  S6.resize(0);
  S7.resize(0);
  Q11.resize(0);
  Q21.resize(0);
  Q12.resize(0);
  Q22.resize(0);
  Qmap1.resize(0);
  Qmap2.resize(0);
}

Doubles_1::Doubles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space)
{
  Tmap.resize(Chan.size1);
  Evec.resize(Chan.size1);
  T1.resize(Chan.size1);
  T2.resize(Chan.size2);
  T3.resize(Chan.size2);
  T4.resize(Chan.size2);
  T5.resize(Chan.size2);
  T6.resize(Chan.size3);
  T7.resize(Chan.size3);
  T8.resize(Chan.size3);
  T9.resize(Chan.size3);
  S1.resize(Chan.size1);
  S2.resize(Chan.size3);
  S3.resize(Chan.size3);
  S4.resize(Chan.size3);
  S5.resize(Chan.size3);
  S6.resize(Chan.size2);
  S7.resize(Chan.size2);

  Q11.resize(Chan.size1);
  Q21.resize(Chan.size1);
  Qmap1.resize(Chan.size1);
  Qmap2.resize(Chan.size1);
  Q12.resize(Chan.size3);
  Q22.resize(Chan.size3);

  for(int i = 0; i < Chan.size1; ++i){
    Tmap[i].assign(Chan.hh[i] * Chan.pp[i] * 16, 0);
    Evec[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    T1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    S1[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    Q11[i].assign(Chan.hp[i] * Chan.pp[i], 0.0);
    Q21[i].assign(Chan.hh[i] * Chan.hp[i], 0.0);
    Qmap1[i].assign(2 * Chan.hp[i] * Chan.pp[i], 0);
    Qmap2[i].assign(2 * Chan.hh[i] * Chan.hp[i], 0);
  }
  for(int i = 0; i < Chan.size3; ++i){
    T6[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    T7[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    T8[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
    T9[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
    S2[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    S3[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    S4[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    S5[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    Q12[i].assign(Chan.hpp[i] * Chan.p[i], 0.0);
    Q22[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
  }
  for(int i = 0; i < Chan.size2; ++i){
    T2[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    T3[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    T4[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    T5[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    S6[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    S7[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
  }
  for(int i = 0; i < Chan.size1; ++i){
    State tb;
    for(int hh = 0; hh < Chan.hh[i]; ++hh){
      int hind1 = Chan.hhvec1[i][2*hh];
      int hind2 = Chan.hhvec1[i][2*hh + 1];
      for(int pp = 0; pp < Chan.pp[i]; ++pp){
	int ind, ind1;
	int ind0 = hh * Chan.pp[i] + pp;
	int pind1 = Chan.ppvec1[i][2*pp];
	int pind2 = Chan.ppvec1[i][2*pp + 1];
	//T2 = ((ia)(jb)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind1, pind1, hind2, pind2);
	Tmap[i][16 * ind0] = ind;
	Tmap[i][16 * ind0 + 1] = ind1;
	//T3 = ((jb)(ia)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
    	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind2, pind2, hind1, pind1);
	Tmap[i][16 * ind0 + 2] = ind;
	Tmap[i][16 * ind0 + 3] = ind1;
	//T4 = ((ib)(ja)')
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
     	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind1, pind2, hind2, pind1);
	Tmap[i][16 * ind0 + 4] = ind;
	Tmap[i][16 * ind0 + 5] = ind1;
	//T5 = ((ja)(ib)')
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind2, pind1, hind1, pind2);
	Tmap[i][16 * ind0 + 6] = ind;
	Tmap[i][16 * ind0 + 7] = ind1;
	//T6 = ((jab)(i))
	ind = Chan.indvec[hind1];
	ind1 = Index31(Chan.hppvec1[ind], Chan.hvec1[ind], Chan.hpp[ind], Chan.h[ind], hind2, pind1, pind2, hind1);
	Tmap[i][16 * ind0 + 8] = ind;
	Tmap[i][16 * ind0 + 9] = ind1;
	//T7 = ((iab)(j))
	ind = Chan.indvec[hind2];
	ind1 = Index31(Chan.hppvec1[ind], Chan.hvec1[ind], Chan.hpp[ind], Chan.h[ind], hind1, pind1, pind2, hind2);
	Tmap[i][16 * ind0 + 10] = ind;
	Tmap[i][16 * ind0 + 11] = ind1;
	//T8 = ((ijb)(a))
	ind = Chan.indvec[pind1];
	ind1 = Index31(Chan.hhpvec1[ind], Chan.pvec1[ind], Chan.hhp[ind], Chan.p[ind], hind1, hind2, pind2, pind1);
	Tmap[i][16 * ind0 + 12] = ind;
	Tmap[i][16 * ind0 + 13] = ind1;
	//T9 = ((ija)(b))
	ind = Chan.indvec[pind2];
	ind1 = Index31(Chan.hhpvec1[ind], Chan.pvec1[ind], Chan.hhp[ind], Chan.p[ind], hind1, hind2, pind1, pind2);
	Tmap[i][16 * ind0 + 14] = ind;
	Tmap[i][16 * ind0 + 15] = ind1;
      }
    }	
  }

  for(int i = 0; i < Chan.size3; ++i){
    int ind, ind1;
    for(int hp = 0; hp < Chan.hp[i]; ++hp){
      int hind1 = Chan.hpvec1[i][2*hp];
      int pind1 = Chan.hpvec1[i][2*hp + 1];
      for(int pp = 0; pp < Chan.pp[i]; ++pp){
	int ind0 = hp * Chan.pp[i] + pp;
	int pind2 = Chan.ppvec1[i][2*pp];
	int pind3 = Chan.ppvec1[i][2*pp + 1];
	ind = Chan.indvec[pind1];
	ind1 = Index31(Chan.hppvec1[ind], Chan.pvec1[ind], Chan.hpp[ind], Chan.p[ind], hind1, pind2, pind3, pind1);
	Qmap1[i][2 * ind0] = ind;
	Qmap1[i][2 * ind0 + 1] = ind1;
      }
    }
    for(int hh = 0; hh < Chan.hh[i]; ++hh){
      int hind1 = Chan.hhvec1[i][2*hh];
      int hind2 = Chan.hhvec1[i][2*hh + 1];
      for(int hp = 0; hp < Chan.hp[i]; ++hp){
	int ind0 = hh * Chan.hp[i] + hp;
	int hind3 = Chan.hpvec1[i][2*hp];
	int pind1 = Chan.hpvec1[i][2*hp + 1];
	ind = Chan.indvec[hind3];
	ind1 = Index31(Chan.hhpvec1[ind], Chan.hvec1[ind], Chan.hhp[ind], Chan.h[ind], hind1, hind2, pind1, hind3);
	Qmap2[i][2 * ind0] = ind;
	Qmap2[i][2 * ind0 + 1] = ind1;
      }
    }
  }

  /*std::ofstream mapfile;
  mapfile.open ("mapfileN4.txt");
  for(int i = 0; i < Chan.size1; ++i){ mapfile << Chan.hh[i] * Chan.pp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size2; ++i){ mapfile << Chan.hp1[i] * Chan.hp2[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i){ mapfile << Chan.h[i] * Chan.hpp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i){ mapfile << Chan.p[i] * Chan.hhp[i] << " "; }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.hh[i]*Chan.pp[i]; ++j){
      mapfile << i << " ";
      for(int k = 0; k < 15; ++k){ mapfile << Tmap[i][15*j + k] << " "; }
      mapfile << std::endl;
    }
  }
  mapfile.close();*/
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

void Doubles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  double fac1 = 1.0, fac2 = 0.0;
  char N = 'N';
  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.hh[i];
    int pp = Chan.pp[i];
    int hp = Chan.hp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    RM_dgemm(& *Ints.S_ME1.V19[i].begin(), & *T1[i].begin(), & *Q11[i].begin(), &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    RM_dgemm(& *T1[i].begin(), & *Ints.S_ME1.V20[i].begin(), & *Q21[i].begin(), &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.hp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
	int ind0 = hp1 * Chan.pp[i] + pp1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.hp[i]; ++hp1){
	int ind0 = hh1 * Chan.hp[i] + hp1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }
}

Singles_1::Singles_1()
{
  Tmap.resize(0);
  Evec.resize(0);
  T1.resize(0);
  T2.resize(0);
  T3.resize(0);
  Tmap2.resize(0);
  E1.resize(0);
  E2.resize(0);
  E3.resize(0);
  E4.resize(0);
  E5.resize(0);
  E6.resize(0);
  E7.resize(0);
  E8.resize(0);
  E9.resize(0);
  S1.resize(0);
  S2.resize(0);
  S3.resize(0);
  S4.resize(0);
  Q11.resize(0);
  Q12.resize(0);
  Q21.resize(0);
  Q22.resize(0);
  Q31.resize(0);
  Q32.resize(0);
  Q41.resize(0);
  Q42.resize(0);
  Q51.resize(0);
  Q52.resize(0);
  Q61.resize(0);
  Q62.resize(0);
  Qmap1.resize(0);
  Qmap2.resize(0);
  Qmap3.resize(0);
  Qmap4.resize(0);
  Qmap5.resize(0);
  Qmap6.resize(0);
}


Singles_1::Singles_1(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space)
{
  Tmap.assign(3 * Chan.hp1[Chan.ind0], 0);
  Evec.assign(Chan.hp1[Chan.ind0], 0.0);
  T1.assign(Chan.hp1[Chan.ind0], 0.0);
  T2.resize(Chan.size3);
  T3.resize(Chan.size3);

  Tmap2.assign(18* Chan.hp1[Chan.ind0] * Chan.hp1[Chan.ind0], 0);

  E1.resize(Chan.size1);
  E2.resize(Chan.size2);
  E3.resize(Chan.size2);
  E4.resize(Chan.size2);
  E5.resize(Chan.size2);
  E6.resize(Chan.size3);
  E7.resize(Chan.size3);
  E8.resize(Chan.size3);
  E9.resize(Chan.size3);

  S1.resize(Chan.size3);
  S2.resize(Chan.size3);
  S3.assign(Chan.hp1[Chan.ind0] * Chan.hp1[Chan.ind0], 0.0);
  S4.assign(Chan.hp1[Chan.ind0], 0.0);

  Q11.resize(Chan.size3);
  Q21.resize(Chan.size3);
  Qmap1.resize(Chan.size3);
  Qmap2.resize(Chan.size3);
  Q12.resize(Chan.size2);
  Q22.resize(Chan.size2);

  Q31.assign(Chan.hh1[Chan.ind0], 0.0);
  Q41.assign(Chan.pp1[Chan.ind0], 0.0);
  Qmap3.assign(2 * Chan.hh1[Chan.ind0], 0);
  Qmap4.assign(2 * Chan.pp1[Chan.ind0], 0);
  Q32.resize(Chan.size3);
  Q42.resize(Chan.size3);

  Q51.resize(Chan.size1);
  Q61.resize(Chan.size1);
  Qmap5.resize(Chan.size1);
  Qmap6.resize(Chan.size1);
  Q52.resize(Chan.size3);
  Q62.resize(Chan.size3);

  for(int i = 0; i < Chan.size1; ++i){
    E1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    Q51[i].assign(Chan.hp[i] * Chan.pp[i], 0.0);
    Q61[i].assign(Chan.hh[i] * Chan.hp[i], 0.0);
    Qmap5[i].assign(2 * Chan.hp[i] * Chan.pp[i], 0);
    Qmap6[i].assign(2 * Chan.hh[i] * Chan.hp[i], 0);
  }
  for(int i = 0; i < Chan.size3; ++i){
    T2[i].assign(Chan.p[i] * Chan.h[i], 0.0);
    T3[i].assign(Chan.h[i] * Chan.p[i], 0.0);

    S1[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    S2[i].assign(Chan.h[i] * Chan.h[i], 0.0);

    E6[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    E7[i].assign(Chan.hpp[i] * Chan.h[i], 0.0);
    E8[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);
    E9[i].assign(Chan.hhp[i] * Chan.p[i], 0.0);

    Q11[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    Q21[i].assign(Chan.hhp2[i] * Chan.h[i], 0.0);
    Qmap1[i].assign(2 * Chan.hpp2[i] * Chan.h[i], 0);
    Qmap2[i].assign(2 * Chan.hhp2[i] * Chan.p[i], 0);

    Q32[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    Q42[i].assign(Chan.p[i] * Chan.p[i], 0.0);

    Q52[i].assign(Chan.hpp[i] * Chan.p[i], 0.0);
    Q62[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
  }
  for(int i = 0; i < Chan.size2; ++i){
    E2[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    E3[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    E4[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    E5[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);

    Q12[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
    Q22[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
  }

  for(int hp1 = 0; hp1 < Chan.hp1[Chan.ind0]; ++hp1){
    State tb;
    int ind, ind1;
    int hind1 = Chan.hp1vec1[Chan.ind0][2*hp1];
    int pind2 = Chan.hp1vec1[Chan.ind0][2*hp1 + 1];
    ind = Chan.indvec[hind1];
    Tmap[3 * hp1] = ind;
    ind1 = Index11(Chan.pvec1[ind], Chan.hvec1[ind], Chan.p[ind], Chan.h[ind], pind2, hind1);
    Tmap[3 * hp1 + 1] = ind1;
    ind1 = Index11(Chan.hvec1[ind], Chan.pvec1[ind], Chan.h[ind], Chan.p[ind], hind1, pind2);
    Tmap[3 * hp1 + 2] = ind1;    
    for(int hp2 = 0; hp2 < Chan.hp1[Chan.ind0]; ++hp2){
      int ind0 = hp1 * Chan.hp1[Chan.ind0] + hp2;
      int hind2 = Chan.hp1vec1[Chan.ind0][2*hp2];
      int pind1 = Chan.hp1vec1[Chan.ind0][2*hp2 + 1];
      //E1 = ((ij)(ab))
      plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
      ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
      ind1 = Index22(Chan.hhvec1[ind], Chan.ppvec1[ind], Chan.hh[ind], Chan.pp[ind], hind1, hind2, pind1, pind2);
      Tmap2[18 * ind0] = ind;
      Tmap2[18 * ind0 + 1] = ind1;
      //E2 = ((ia)(jb)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind1, pind1, hind2, pind2);
      Tmap2[18 * ind0 + 2] = ind;
      Tmap2[18 * ind0 + 3] = ind1;
      //E3 = ((jb)(ia)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind2, pind2, hind1, pind1);
      Tmap2[18 * ind0 + 4] = ind;
      Tmap2[18 * ind0 + 5] = ind1;
      //E4 = ((ib)(ja)')
      minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind1, pind2, hind2, pind1);
      Tmap2[18 * ind0 + 6] = ind;
      Tmap2[18 * ind0 + 7] = ind1;
      //E5 = ((ja)(ib)')
      minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
      ind1 = Index22(Chan.hp1vec1[ind], Chan.hp2vec1[ind], Chan.hp1[ind], Chan.hp2[ind], hind2, pind1, hind1, pind2);
      Tmap2[18 * ind0 + 8] = ind;
      Tmap2[18 * ind0 + 9] = ind1;
      //E6 = ((jab)(i))
      ind = Chan.indvec[hind1];
      ind1 = Index31(Chan.hppvec1[ind], Chan.hvec1[ind], Chan.hpp[ind], Chan.h[ind], hind2, pind1, pind2, hind1);
      Tmap2[18 * ind0 + 10] = ind;
      Tmap2[18 * ind0 + 11] = ind1;
      //E7 = ((iab)(j))
      ind = Chan.indvec[hind2];
      ind1 = Index31(Chan.hppvec1[ind], Chan.hvec1[ind], Chan.hpp[ind], Chan.h[ind], hind1, pind1, pind2, hind2);
      Tmap2[18 * ind0 + 12] = ind;
      Tmap2[18 * ind0 + 13] = ind1;
      //E8 = ((ijb)(a))
      ind = Chan.indvec[pind1];
      ind1 = Index31(Chan.hhpvec1[ind], Chan.pvec1[ind], Chan.hhp[ind], Chan.p[ind], hind1, hind2, pind2, pind1);
      Tmap2[18 * ind0 + 14] = ind;
      Tmap2[18 * ind0 + 15] = ind1;
      //E9 = ((ija)(b))
      ind = Chan.indvec[pind2];
      ind1 = Index31(Chan.hhpvec1[ind], Chan.pvec1[ind], Chan.hhp[ind], Chan.p[ind], hind1, hind2, pind1, pind2);
      Tmap2[18 * ind0 + 16] = ind;
      Tmap2[18 * ind0 + 17] = ind1;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    State tb;
    int ind, ind1;
    for(int hpp = 0; hpp < Chan.hpp2[i]; ++hpp){
      int hind2 = Chan.hpp2vec1[i][3*hpp];
      int pind1 = Chan.hpp2vec1[i][3*hpp + 1];
      int pind2 = Chan.hpp2vec1[i][3*hpp + 2];
      for(int h = 0; h < Chan.h[i]; ++h){
	int ind0 = hpp * Chan.h[i] + h;
	int hind1 = Chan.hvec1[i][h];
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp1vec1[ind], Chan.hp1[ind], Chan.hp1[ind], hind1, pind1, hind2, pind2);
	Qmap1[i][2 * ind0] = ind;
	Qmap1[i][2 * ind0 + 1] = ind1;
      }
    }
    for(int hhp = 0; hhp < Chan.hhp2[i]; ++hhp){
      int hind2 = Chan.hhp2vec1[i][3*hhp];
      int hind1 = Chan.hhp2vec1[i][3*hhp + 1];
      int pind2 = Chan.hhp2vec1[i][3*hhp + 2];
      for(int p = 0; p < Chan.p[i]; ++p){
	int ind0 = hhp * Chan.p[i] + p;
	int pind1 = Chan.pvec1[i][p];
	minus(tb, Space.qnums[hind1], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	ind1 = Index22(Chan.hp1vec1[ind], Chan.hp1vec1[ind], Chan.hp1[ind], Chan.hp1[ind], hind1, pind1, hind2, pind2);
	Qmap2[i][2 * ind0] = ind;
	Qmap2[i][2 * ind0 + 1] = ind1;
      }
    }
  }

  int ind, ind1;
  for(int hh1 = 0; hh1 < Chan.hh1[Chan.ind0]; ++hh1){
    int hind1 = Chan.hh1vec1[Chan.ind0][2*hh1];
    int hind2 = Chan.hh1vec1[Chan.ind0][2*hh1 + 1];
    ind = Chan.indvec[hind1];
    ind1 = Index11(Chan.hvec1[ind], Chan.hvec1[ind], Chan.h[ind], Chan.h[ind], hind2, hind1);
    Qmap3[2 * hh1] = ind;
    Qmap3[2 * hh1 + 1] = ind1;
  }
  for(int pp1 = 0; pp1 < Chan.pp1[Chan.ind0]; ++pp1){
    int pind1 = Chan.pp1vec1[Chan.ind0][2*pp1];
    int pind2 = Chan.pp1vec1[Chan.ind0][2*pp1 + 1];
    ind = Chan.indvec[pind1];
    ind1 = Index11(Chan.pvec1[ind], Chan.pvec1[ind], Chan.p[ind], Chan.p[ind], pind1, pind2);
    Qmap4[2 * pp1] = ind;
    Qmap4[2 * pp1 + 1] = ind1;
  }

  for(int i = 0; i < Chan.size3; ++i){
    int ind, ind1;
    for(int hp = 0; hp < Chan.hp[i]; ++hp){
      int hind1 = Chan.hpvec1[i][2*hp];
      int pind1 = Chan.hpvec1[i][2*hp + 1];
      for(int pp = 0; pp < Chan.pp[i]; ++pp){
	int ind0 = hp * Chan.pp[i] + pp;
	int pind2 = Chan.ppvec1[i][2*pp];
	int pind3 = Chan.ppvec1[i][2*pp + 1];
	ind = Chan.indvec[pind1];
	ind1 = Index31(Chan.hppvec1[ind], Chan.pvec1[ind], Chan.hpp[ind], Chan.p[ind], hind1, pind2, pind3, pind1);
	Qmap5[i][2 * ind0] = ind;
	Qmap5[i][2 * ind0 + 1] = ind1;
      }
    }
    for(int hh = 0; hh < Chan.hh[i]; ++hh){
      int hind1 = Chan.hhvec1[i][2*hh];
      int hind2 = Chan.hhvec1[i][2*hh + 1];
      for(int hp = 0; hp < Chan.hp[i]; ++hp){
	int ind0 = hh * Chan.hp[i] + hp;
	int hind3 = Chan.hpvec1[i][2*hp];
	int pind1 = Chan.hpvec1[i][2*hp + 1];
	ind = Chan.indvec[hind3];
	ind1 = Index31(Chan.hhpvec1[ind], Chan.hvec1[ind], Chan.hhp[ind], Chan.h[ind], hind1, hind2, pind1, hind3);
	Qmap6[i][2 * ind0] = ind;
	Qmap6[i][2 * ind0 + 1] = ind1;
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
  double tempt  = T1[i];
  tempt += T2[Tmap[3*i]][Tmap[3*i + 1]];
  tempt += T3[Tmap[3*i]][Tmap[3*i + 2]];
  return tempt;
}

void Singles_1::set_T_2(const Channels &Chan, Interactions &Ints)
{
  #pragma omp parallel for
  for(int hp1 = 0; hp1 < Chan.hp1[Chan.ind0]; ++hp1){
    for(int hp2 = 0; hp2 < Chan.hp1[Chan.ind0]; ++hp2){
      int ind0 = hp1 * Chan.hp1[Chan.ind0] + hp2;
      double T = T1[hp1] * T1[hp2];
      E1[Tmap2[18 * ind0]][Tmap2[18 * ind0 + 1]] = T;
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
    int h = Chan.h[i];
    int p = Chan.p[i];
    int hpp2 = Chan.hpp2[i];
    int hhp2 = Chan.hhp2[i];
    if(h == 0 || p == 0){ continue; }
    if(hhp2 != 0){ RM_dgemm(& *Ints.S_ME1.V13[i].begin(), & *T2[i].begin(), & *Q11[i].begin(), &hpp2, &h, &p, &fac1, &fac2, &N, &N); }
    if(hpp2 != 0){ RM_dgemm(& *Ints.S_ME1.V14[i].begin(), & *T3[i].begin(), & *Q21[i].begin(), &hhp2, &p, &h, &fac1, &fac2, &N, &N); }
    for(int hpp1 = 0; hpp1 < hpp2; ++hpp1){
      for(int h1 = 0; h1 < Chan.h[i]; ++h1){
	int ind0 = hpp1 * Chan.h[i] + h1;
	Q12[Qmap1[i][2 * ind0]][Qmap1[i][2 * ind0 + 1]] = Q11[i][ind0];
      }
    }
    for(int hhp1 = 0; hhp1 < hhp2; ++hhp1){
      for(int p1 = 0; p1 < Chan.p[i]; ++p1){
	int ind0 = hhp1 * Chan.p[i] + p1;
	Q22[Qmap2[i][2 * ind0]][Qmap2[i][2 * ind0 + 1]] = Q21[i][ind0];
      }
    }
  }

  int hp = Chan.hp1[Chan.ind0];
  int hh = Chan.hh1[Chan.ind0];
  int pp = Chan.pp1[Chan.ind0];
  if(hp != 0 && hh != 0){ RM_dgemm(& *Ints.S_ME1.V15[Chan.ind0].begin(), & *T1.begin(), & *Q31.begin(), &hh, &one, &hp, &fac1, &fac2, &N, &N); }
  if(hp != 0 && pp != 0){ RM_dgemm(& *Ints.S_ME1.V16[Chan.ind0].begin(), & *T1.begin(), & *Q41.begin(), &pp, &one, &hp, &fac1, &fac2, &N, &N); }
  for(int hh1 = 0; hh1 < Chan.hh1[Chan.ind0]; ++hh1){ Q32[Qmap3[2 * hh1]][Qmap3[2 * hh1 + 1]] = Q31[hh1]; }
  for(int pp1 = 0; pp1 < Chan.pp1[Chan.ind0]; ++pp1){ Q42[Qmap4[2 * pp1]][Qmap4[2 * pp1 + 1]] = Q41[pp1]; }

  for(int i = 0; i < Chan.size1; ++i){
    int hh = Chan.hh[i];
    int pp = Chan.pp[i];
    int hp = Chan.hp[i];
    if(hh == 0 || pp == 0 || hp == 0){ continue; }
    RM_dgemm(& *Ints.S_ME1.V19[i].begin(), & *E1[i].begin(), & *Q51[i].begin(), &hp, &pp, &hh, &fac1, &fac2, &N, &N);
    RM_dgemm(& *E1[i].begin(), & *Ints.S_ME1.V20[i].begin(), & *Q61[i].begin(), &hh, &hp, &pp, &fac1, &fac2, &N, &N);
    for(int hp1 = 0; hp1 < Chan.hp[i]; ++hp1){
      for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
	int ind0 = hp1 * Chan.pp[i] + pp1;
	Q52[Qmap5[i][2 * ind0]][Qmap5[i][2 * ind0 + 1]] = Q51[i][ind0];
      }
    }
    for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
      for(int hp1 = 0; hp1 < Chan.hp[i]; ++hp1){
	int ind0 = hh1 * Chan.hp[i] + hp1;
	Q62[Qmap6[i][2 * ind0]][Qmap6[i][2 * ind0 + 1]] = Q61[i][ind0];
      }
    }
  }
}


void Setup_Ints(const Input_Parameters &Parameters, const Channels &Chan, Interactions &Ints)
{
  if(Parameters.approx == "doubles"){
    Ints.D_ME1 = Doubles_ME1(Chan);
  }
  else if(Parameters.approx == "singles"){
    Ints.D_ME1 = Doubles_ME1(Chan);
    Ints.S_ME1 = Singles_ME1(Chan);
  }
  else if(Parameters.approx == "triples"){
    Ints.D_ME1 = Doubles_ME1(Chan);
    Ints.S_ME1 = Singles_ME1(Chan);
  }
}

Doubles_ME1::Doubles_ME1()
{
  V1.resize(0);
  V2.resize(0);
  V3.resize(0);
  V4.resize(0);
  V5.resize(0);
  V6.resize(0);
  V7.resize(0);
  V8.resize(0);
  V9.resize(0);
  V10.resize(0);
}

Doubles_ME1::Doubles_ME1(const Channels &Chan)
{
  V1.resize(Chan.size1);
  V2.resize(Chan.size1);
  V3.resize(Chan.size2);
  V4.resize(Chan.size1);
  V5.resize(Chan.size3);
  V6.resize(Chan.size3);
  V7.resize(Chan.size3);
  V8.resize(Chan.size3);
  V9.resize(Chan.size2);
  V10.resize(Chan.size2);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    V1[i].assign(Chan.pp[i] * Chan.pp[i], 0.0);
    V2[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    V4[i].assign(Chan.pp[i] * Chan.hh[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    V5[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
    V6[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
    V7[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
    V8[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    V3[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    V9[i].assign(Chan.hp2[i] * Chan.hp1[i], 0.0);
    V10[i].assign(Chan.hp2[i] * Chan.hp1[i], 0.0);
  }
}

Singles_ME1::Singles_ME1()
{
  V11.resize(0);
  V12.resize(0);
  V17.resize(0);
  V18.resize(0);
  V13.resize(0);
  V14.resize(0);
  V20.resize(0);
  V19.resize(0);
  V16.resize(0);
  V15.resize(0);
}

Singles_ME1::Singles_ME1(const Channels &Chan)
{
  V11.resize(Chan.size3);
  V12.resize(Chan.size3);
  V17.resize(Chan.size3);
  V18.resize(Chan.size3);
  V13.resize(Chan.size3);
  V14.resize(Chan.size3);
  V20.resize(Chan.size1);
  V19.resize(Chan.size1);
  V16.resize(Chan.size2);
  V15.resize(Chan.size2);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    V20[i].assign(Chan.pp[i] * Chan.hp[i], 0.0);
    V19[i].assign(Chan.hp[i] * Chan.hh[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    V11[i].assign(Chan.p[i] * Chan.hpp[i], 0.0);
    V12[i].assign(Chan.h[i] * Chan.hhp[i], 0.0);
    V17[i].assign(Chan.hpp[i] * Chan.p[i], 0.0);
    V18[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
    V13[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    V14[i].assign(Chan.hhp2[i] * Chan.h[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    V16[i].assign(Chan.pp1[i] * Chan.hp1[i], 0.0);
    V15[i].assign(Chan.hh1[i] * Chan.hp1[i], 0.0);
  }
}


//Function to setup Channels
void Setup_Channels(const Input_Parameters &Parameters, const Model_Space &Space, Channels &Chan)
{
  std::cout << "Building Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  //std::vector<std::vector<int> > OB_nums; // Ms,Ts,Ps for finite; Ms,Ts,Nxs,Nys,Nzs for infinite
  std::vector<std::vector<int> > Hvec, Pvec;
  State state;

  // place first state by quantum numbers
  Chan.size3 = 1;
  Hvec.resize(Chan.size3);
  Pvec.resize(Chan.size3);
  if(Parameters.basis == "infinite"){
    state.t = Space.qnums[0].t; //tz
    state.m = Space.qnums[0].m; //sz
    state.nx = Space.qnums[0].nx; //nx
    state.ny = Space.qnums[0].ny; //ny
    state.nz = Space.qnums[0].nz; //nz
  }
  else if(Parameters.basis == "finite_M"){
    state.t = Space.qnums[0].t; //tz
    state.m = Space.qnums[0].m; //jz
    state.par = Space.qnums[0].par; //par
  }
  else if(Parameters.basis == "finite_HO"){
    state.t = Space.qnums[0].t; //tz
    state.m = Space.qnums[0].m; //sz
    state.ml = Space.qnums[0].ml; //lz
  }
  else if(Parameters.basis == "finite_J"){
    state.t = Space.qnums[0].t; //tz
    state.j = Space.qnums[0].j; //j
    state.par = Space.qnums[0].par; //par
  }
  Chan.qnums3.push_back(state);
  Chan.indvec.push_back(Chan.size3 - 1);
  if(Space.qnums[0].type == "hole"){ Hvec[Chan.size3 - 1].push_back(0); }
  else{ Pvec[Chan.size3 - 1].push_back(0); }

  // place the rest of the states by their quantum numbers
  for(int i = 1; i < Space.indtot; ++i){
    if(Parameters.basis == "infinite"){
      state.t = Space.qnums[i].t; //tz
      state.m = Space.qnums[i].m; //sz
      state.nx = Space.qnums[i].nx; //nx
      state.ny = Space.qnums[i].ny; //ny
      state.nz = Space.qnums[i].nz; //nz
    }
    else if(Parameters.basis == "finite_M"){
      state.t = Space.qnums[i].t; //tz
      state.m = Space.qnums[i].m; //jz
      state.par = Space.qnums[i].par; //par
    }
    else if(Parameters.basis == "finite_HO"){
      state.t = Space.qnums[i].t; //tz
      state.m = Space.qnums[i].m; //sz
      state.ml = Space.qnums[i].ml; //lz
    }
    else if(Parameters.basis == "finite_J"){
      state.t = Space.qnums[i].t; //tz
      state.j = Space.qnums[i].j; //j
      state.par = Space.qnums[i].par; //par
    }
    for(int k = 0; k < Chan.size3; ++k){
      if(Parameters.basis == "infinite"){
	if(state.t == Chan.qnums3[k].t && state.m == Chan.qnums3[k].m && 
	   state.nx == Chan.qnums3[k].nx && state.ny == Chan.qnums3[k].ny && state.nz == Chan.qnums3[k].nz){
	  Chan.indvec.push_back(k);
	  if(Space.qnums[i].type == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  goto stop;
	}
      }
      else if(Parameters.basis == "finite_M"){
	if(state.t == Chan.qnums3[k].t && state.m == Chan.qnums3[k].m && state.par == Chan.qnums3[k].par){
	  Chan.indvec.push_back(k);
	  if(Space.qnums[i].type == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  goto stop;
	}
      }
      else if(Parameters.basis == "finite_HO"){
	if(state.t == Chan.qnums3[k].t && state.m == Chan.qnums3[k].m && state.ml == Chan.qnums3[k].ml){
	  Chan.indvec.push_back(k);
	  if(Space.qnums[i].type == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  goto stop;
	}
      }
      else if(Parameters.basis == "finite_J"){
	if(state.t == Chan.qnums3[k].t && state.j == Chan.qnums3[k].j && state.par == Chan.qnums3[k].par){
	  Chan.indvec.push_back(k);
	  if(Space.qnums[i].type == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  goto stop;
	}
      }
      if(k == Chan.size3 - 1){
	++Chan.size3;
	Hvec.resize(Chan.size3);
	Pvec.resize(Chan.size3);
	Chan.qnums3.push_back(state);
	Chan.indvec.push_back(Chan.size3 - 1);
	if(Space.qnums[i].type == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
	else{ Pvec[Chan.size3 - 1].push_back(i); }
	break;
      }
    }
  stop:;
  }

  Chan.size1 = Space.Chansize_2b;
  Chan.size2 = Space.Chansize_2b;
  std::cout << "chan.size1 = " << Chan.size1 << ", chan.size2 = " << Chan.size2 << ", chan.size3 = " << Chan.size3 << std::endl;

  Chan.qnums1.resize(Chan.size1);
  Chan.qnums2.resize(Chan.size2);

  // Make Vector of Two-Body States for Each Channel
  std::vector<std::vector<int> > tempvec1(Chan.size3);
  std::vector<std::vector<int> > tempvec2(Chan.size3);
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2;
    tempvec1[i].assign(Chan.size1, -1);
    tempvec2[i].assign(Chan.size2, -1);
    for(int j = 0; j < Chan.size3; ++j){
      if(Parameters.basis != "finite_J"){
	plus(state, Space.qnums[i], Space.qnums[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	tempvec1[i][ind1] = j;
	Chan.qnums1[ind1] = state;
	minus(state, Space.qnums[i], Space.qnums[j]);
	ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	tempvec2[i][ind2] = j;
	Chan.qnums2[ind2] = state;
      }
      else if(Parameters.basis == "finite_J"){
	plus(state, Space.qnums[i], Space.qnums[j]);
	for(int J = abs(Space.qnums[i].j - Space.qnums[j].j); J < Space.qnums[i].j + Space.qnums[j].j; J+=2){
	  state.j = J;
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  tempvec1[i][ind1] = j;
	  Chan.qnums1[ind1] = state;
	}
	minus(state, Space.qnums[i], Space.qnums[j]);
	for(int J = abs(Space.qnums[i].j - Space.qnums[j].j); J < Space.qnums[i].j + Space.qnums[j].j; J+=2){
	  state.j = J;
	  ind2 = ChanInd_2b_cross(Parameters.basis, Space, state);
	  tempvec2[i][ind2] = j;
	  Chan.qnums2[ind2] = state;
	}
      }
    }
  }

  /*std::cout << "tempvec2.size() = " << tempvec2.size() << std::endl;
  for(int p = 0; p < int(tempvec2.size()); ++p){
    std::cout << "i tempvec2[i].size() = " << p << " " << tempvec2[p].size() << std::endl;
    }*/
  
  // for singles case
  if(Parameters.approx == "singles" || Parameters.approx == "triples"){
    if(Parameters.basis == "finite_M"){
      state.t = 0; //tz
      state.m = 0; //jz
      state.par = 1; //par
    }
    else if(Parameters.basis == "finite_HO"){
      state.t = 0; //tz
      state.m = 0; //sz
      state.ml = 0; //lz
    }
    else if(Parameters.basis == "finite_J"){
      state.t = 0; //tz
      state.j = 0; //j
      state.par = 1; //par
    }
    Chan.ind0 = ChanInd_2b_cross(Parameters.basis, Space, state);
  }
  
  Chan.hhvec1.resize(Chan.size1);
  Chan.ppvec1.resize(Chan.size1);
  Chan.hpvec1.resize(Chan.size1);
  Chan.hp1vec1.resize(Chan.size2);
  Chan.hp2vec1.resize(Chan.size2);
  Chan.pp1vec1.resize(Chan.size2);
  Chan.hh1vec1.resize(Chan.size2);
  for(int i = 0; i < Chan.size3; ++i){
    //#pragma omp parallel for
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      int ind1, ind2;
      int hsize1, hsize2, psize1, psize2;
      ind1 = i;
      ind2 = tempvec1[i][j];
      hsize1 = int(Hvec[ind1].size());
      hsize2 = int(Hvec[ind2].size());
      for(int h1 = 0; h1 < hsize1; ++h1){
	for(int h2 = 0; h2 < hsize2; ++h2){
	  Chan.hhvec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hhvec1[j].push_back(Hvec[ind2][h2]);
	}
      }
      psize1 = int(Pvec[ind1].size());
      psize2 = int(Pvec[ind2].size());
      for(int p1 = 0; p1 < psize1; ++p1){
	for(int p2 = 0; p2 < psize2; ++p2){
	  Chan.ppvec1[j].push_back(Pvec[ind1][p1]);
	  Chan.ppvec1[j].push_back(Pvec[ind2][p2]);
	}
      }
      for(int h1 = 0; h1 < hsize1; ++h1){
	for(int p2 = 0; p2 < psize2; ++p2){
	  Chan.hpvec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hpvec1[j].push_back(Pvec[ind2][p2]);
	}
      }
    }
    //#pragma omp parallel for
    for(int j = 0; j < Chan.size2; ++j){
      if(tempvec2[i][j] == -1){ continue; }
      int ind1, ind2;
      int hsize1, psize1, hsize2, psize2;
      ind1 = i;
      ind2 = tempvec2[i][j];
      hsize1 = int(Hvec[ind1].size());
      psize1 = int(Pvec[ind2].size());
      for(int h1 = 0; h1 < hsize1; ++h1){
	for(int p1 = 0; p1 < psize1; ++p1){
	  Chan.hp1vec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hp1vec1[j].push_back(Pvec[ind2][p1]);
	}
      }
      hsize1 = int(Hvec[ind2].size());
      psize1 = int(Pvec[ind1].size());
      for(int h1 = 0; h1 < hsize1; ++h1){
	for(int p1 = 0; p1 < psize1; ++p1){
	  Chan.hp2vec1[j].push_back(Hvec[ind2][h1]);
	  Chan.hp2vec1[j].push_back(Pvec[ind1][p1]);
	}
      }
      hsize1 = int(Hvec[ind1].size());
      hsize2 = int(Hvec[ind2].size());
      for(int h1 = 0; h1 < hsize1; ++h1){
	for(int h2 = 0; h2 < hsize2; ++h2){
	  Chan.hh1vec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hh1vec1[j].push_back(Hvec[ind2][h2]);
	}
      }
      psize1 = int(Pvec[ind1].size());
      psize2 = int(Pvec[ind2].size());
      for(int p1 = 0; p1 < psize1; ++p1){
	for(int p2 = 0; p2 < psize2; ++p2){
	  Chan.pp1vec1[j].push_back(Pvec[ind1][p1]);
	  Chan.pp1vec1[j].push_back(Pvec[ind2][p2]);
	}
      }
    }
  }

  Chan.hvec1.resize(Chan.size3);
  Chan.pvec1.resize(Chan.size3);
  Chan.hppvec1.resize(Chan.size3);
  Chan.hhpvec1.resize(Chan.size3);
  Chan.hpp2vec1.resize(Chan.size3);
  Chan.hhp2vec1.resize(Chan.size3);
  Chan.hhhvec1.resize(Chan.size3);
  Chan.pppvec1.resize(Chan.size3);
  //#pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind2, ind3, ind4;
    int hsize1, hsize2, psize1, psize2, hsize3, psize3;
    std::vector<int> inds(3);
    hsize1 = int(Hvec[i].size());
    psize1 = int(Pvec[i].size());
    for(int h1 = 0; h1 < hsize1; ++h1){ Chan.hvec1[i].push_back(Hvec[i][h1]); }
    for(int p1 = 0; p1 < psize1; ++p1){ Chan.pvec1[i].push_back(Pvec[i][p1]); }
    for(int j = 0; j < Chan.size1; ++j){ // i+tempvec[i][j] -> j
      if(tempvec1[i][j] == -1){ continue; }
      for(int k = 0; k < Chan.size3; ++k){ // k+tempvec[k][j] -> j
	if(tempvec1[k][j] == -1){ continue; }
	ind2 = tempvec1[i][j];
	ind3 = k;
	ind4 = tempvec1[k][j];
	hsize2 = int(Hvec[ind2].size()); // <h1h2|p1p2> -> <h1|h2p1p2>
	psize1 = int(Pvec[ind3].size());
	psize2 = int(Pvec[ind4].size());
	for(int h2 = 0; h2 < hsize2; ++h2){
	  for(int p1 = 0; p1 < psize1; ++p1){
	    for(int p2 = 0; p2 < psize2; ++p2){
	      inds[0] = Hvec[ind2][h2];
	      inds[1] = Pvec[ind3][p1];
	      inds[2] = Pvec[ind4][p2];
	      Chan.hppvec1[i].insert(Chan.hppvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
	psize2 = int(Pvec[ind2].size()); // <p1p2|h1h2> -> <p1|h1h2p2>
	hsize1 = int(Hvec[ind3].size());
	hsize2 = int(Hvec[ind4].size());
	for(int h1 = 0; h1 < hsize1; ++h1){
	  for(int h2 = 0; h2 < hsize2; ++h2){
	    for(int p2 = 0; p2 < psize2; ++p2){
	      inds[0] = Hvec[ind3][h1];
	      inds[1] = Hvec[ind4][h2];
	      inds[2] = Pvec[ind2][p2];
	      Chan.hhpvec1[i].insert(Chan.hhpvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
	psize3 = int(Pvec[ind2].size()); // <p2p3|h1p1> -> <p2|h1p1p3>
	hsize1 = int(Hvec[ind3].size());
	psize1 = int(Pvec[ind4].size());
	for(int h1 = 0; h1 < hsize1; ++h1){
	  for(int p1 = 0; p1 < psize1; ++p1){
	    for(int p3 = 0; p3 < psize3; ++p3){
	      inds[0] = Hvec[ind3][h1];
	      inds[1] = Pvec[ind4][p1];
	      inds[2] = Pvec[ind2][p3];
	      Chan.hpp2vec1[i].insert(Chan.hpp2vec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
	hsize3 = int(Hvec[ind2].size()); // <h2h3|h1p1> -> <h2|h3h1p1>
	hsize1 = int(Hvec[ind3].size());
	psize1 = int(Pvec[ind4].size());
	for(int h3 = 0; h3 < hsize3; ++h3){
	  for(int h1 = 0; h1 < hsize1; ++h1){
	    for(int p1 = 0; p1 < psize1; ++p1){
	      inds[0] = Hvec[ind2][h3];
	      inds[1] = Hvec[ind3][h1];
	      inds[2] = Pvec[ind4][p1];
	      Chan.hhp2vec1[i].insert(Chan.hhp2vec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
	hsize1 = int(Hvec[ind2].size()); // <p1h1|h2h3> -> <p1|h1h2h3>
	hsize2 = int(Hvec[ind3].size());
	hsize3 = int(Hvec[ind4].size());
	for(int h1 = 0; h1 < hsize1; ++h1){
	  for(int h2 = 0; h2 < hsize2; ++h2){
	    for(int h3 = 0; h3 < hsize3; ++h3){
	      inds[0] = Hvec[ind2][h1];
	      inds[1] = Hvec[ind3][h2];
	      inds[2] = Hvec[ind4][h3];
	      Chan.hhhvec1[i].insert(Chan.hhhvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
	psize1 = int(Pvec[ind2].size()); // <h1p1|p2p3> -> <h1|p1p2p3>
	psize2 = int(Pvec[ind3].size());
	psize3 = int(Pvec[ind4].size());
	for(int p1 = 0; p1 < psize1; ++p1){
	  for(int p2 = 0; p2 < psize2; ++p2){
	    for(int p3 = 0; p3 < psize3; ++p3){
	      inds[0] = Pvec[ind2][p1];
	      inds[1] = Pvec[ind3][p2];
	      inds[2] = Pvec[ind4][p3];
	      Chan.pppvec1[i].insert(Chan.pppvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
      }
    }
  }

  double memory = 0.0;
  std::vector<int> mem1(Chan.size1, 0.0);
  std::vector<int> mem2(Chan.size2, 0.0);
  std::vector<int> mem3(Chan.size3, 0.0);

  Chan.hh.resize(Chan.size1);
  Chan.pp.resize(Chan.size1);
  Chan.hp.resize(Chan.size1);
  Chan.hp1.resize(Chan.size2);
  Chan.hp2.resize(Chan.size2);
  Chan.h.resize(Chan.size3);
  Chan.p.resize(Chan.size3);
  Chan.hpp.resize(Chan.size3);
  Chan.hhp.resize(Chan.size3);
  Chan.hpp2.resize(Chan.size3);
  Chan.hhp2.resize(Chan.size3);
  Chan.hh1.resize(Chan.size2);
  Chan.pp1.resize(Chan.size2);
  Chan.hhh.resize(Chan.size3);
  Chan.ppp.resize(Chan.size3);
  memory += (7.0 * Chan.size1 + 1) * 24.0;
  //#pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    Chan.hh[i] = int(Chan.hhvec1[i].size()/2);
    Chan.pp[i] = int(Chan.ppvec1[i].size()/2);
    Chan.hp[i] = int(Chan.hpvec1[i].size()/2);
    mem1[i] = 8.0 * (2.0 * Chan.hh[i]*Chan.hh[i] + Chan.pp[i]*Chan.pp[i] + 3.0 * Chan.hh[i]*Chan.pp[i]) + 4.0 * (15.0 * Chan.hh[i]*Chan.pp[i]);
  }
  memory += (9.0 * Chan.size2 + 1) * 24.0;
  //#pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    Chan.hp1[i] = int(Chan.hp1vec1[i].size()/2);
    Chan.hp2[i] = int(Chan.hp2vec1[i].size()/2);
    Chan.hh1[i] = int(Chan.hh1vec1[i].size()/2);
    Chan.pp1[i] = int(Chan.pp1vec1[i].size()/2);
    mem2[i] = 8.0 * (5.0 * Chan.hp1[i]*Chan.hp2[i] + 3.0 * Chan.hp1[i]*Chan.hp1[i] + Chan.hp2[i]*Chan.hp2[i]);
  }
  memory += (8.0 * Chan.size3 + 1) * 24.0;
  //#pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    Chan.h[i] = int(Chan.hvec1[i].size());
    Chan.p[i] = int(Chan.pvec1[i].size());
    Chan.hpp[i] = int(Chan.hppvec1[i].size()/3);
    Chan.hhp[i] = int(Chan.hhpvec1[i].size()/3);
    Chan.hpp2[i] = int(Chan.hpp2vec1[i].size()/3);
    Chan.hhp2[i] = int(Chan.hhp2vec1[i].size()/3);
    Chan.hhh[i] = int(Chan.hhhvec1[i].size()/3);
    Chan.ppp[i] = int(Chan.pppvec1[i].size()/3);
    mem3[i] = 8.0 * (Chan.h[i]*Chan.h[i] + 3.0 * Chan.h[i]*Chan.hpp[i] + Chan.p[i]*Chan.p[i] + 3.0 * Chan.p[i]*Chan.hhp[i]);
  }

  for(int i = 0; i < Chan.size1; ++i){ memory += mem1[i]; }
  for(int i = 0; i < Chan.size2; ++i){ memory += mem2[i]; }
  for(int i = 0; i < Chan.size3; ++i){ memory += mem3[i]; }
  std::cout << "2B1 Channels = " << Chan.size1 << ", 2B2 Channels = " << Chan.size2 << ", OB Channels = " << Chan.size3 << std::endl; 
  std::cout << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
}

CC_Eff::CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  X_ia1.assign(Chan.hp1[Chan.ind0], 0.0); // (ia)
  X_ia2.resize(Chan.size3);
  X_ia3.resize(Chan.size3);
  Map_ia.assign(3 * Chan.hp1[Chan.ind0], 0);

  X_ab1.assign(Chan.pp1[Chan.ind0], 0.0); // (ba)
  X_ab2.resize(Chan.size3);
  X_ab3.resize(Chan.size3);
  Map_ab.assign(3 * Chan.pp1[Chan.ind0], 0);

  X_ij1.assign(Chan.hh1[Chan.ind0], 0.0); // (ji)
  X_ij2.resize(Chan.size3);
  X_ij3.resize(Chan.size3);
  X1_ij1.assign(Chan.hh1[Chan.ind0], 0.0); // (ji)
  X1_ij2.resize(Chan.size3);
  X1_ij3.resize(Chan.size3);
  Map_ij.assign(3 * Chan.hh1[Chan.ind0], 0);

  X_ai1.assign(Chan.hp1[Chan.ind0], 0.0); // (ia)
  X_ai2.resize(Chan.size3);
  X_ai3.resize(Chan.size3);
  Map_ai.assign(3 * Chan.hp1[Chan.ind0], 0);

  X_ijab1.resize(Chan.size1);

  X1_iabc1.resize(Chan.size3);
  X1_iabc2.resize(Chan.size3);
  X1_iabc3.resize(Chan.size3);
  X_iabc1.resize(Chan.size3);
  X_iabc3.resize(Chan.size3);
  X_iabc4.resize(Chan.size2);
  X_iabc5.resize(Chan.size1);
  Map_iabc.resize(Chan.size3);

  X1_ijka1.resize(Chan.size3);
  X1_ijka2.resize(Chan.size3);
  X_ijka1.resize(Chan.size3);
  X_ijka4.resize(Chan.size2);
  X_ijka5.resize(Chan.size1);
  Map_ijka.resize(Chan.size3);

  X1_abcd1.resize(Chan.size1);
  X1_abcd2.resize(Chan.size3);
  X1_abcd3.resize(Chan.size3);
  X_abcd1.resize(Chan.size1);
  V_abcd.resize(Chan.size3);
  Map_abcd.resize(Chan.size1);

  X_ijkl1.resize(Chan.size1);
  X_ijkl2.resize(Chan.size3);
  X_ijkl3.resize(Chan.size3);
  X_ijkl4.resize(Chan.size1);
  V_ijkl.resize(Chan.size3);
  Map_ijkl.resize(Chan.size1);

  X1_iajb1.resize(Chan.size2);
  X1_iajb2.resize(Chan.size3);
  X1_iajb3.resize(Chan.size3);
  X1_iajb4.resize(Chan.size3);
  X3_iajb1.resize(Chan.size2);
  X3_iajb2.resize(Chan.size3);
  X3_iajb3.resize(Chan.size3);
  X3_iajb5.resize(Chan.size3);
  X_iajb1.resize(Chan.size2);
  X_iajb3.resize(Chan.size3);
  Map_iajb.resize(Chan.size2);

  X_abic1.resize(Chan.size3);
  X_abic2.resize(Chan.size3);
  X_abic3.resize(Chan.size3);
  X_abic4.resize(Chan.size3);
  X_abic5.resize(Chan.size2);
  X_abic6.resize(Chan.size2);
  X_abic7.resize(Chan.size1);
  Map_abic.resize(Chan.size3);

  X2_iajk1.resize(Chan.size3);
  X2_iajk2.resize(Chan.size3);
  X2_iajk3.resize(Chan.size3);
  X2_iajk4.resize(Chan.size3);
  X2_iajk5.resize(Chan.size2);
  X2_iajk6.resize(Chan.size2);
  X2_iajk7.resize(Chan.size1);
  X_iajk1.resize(Chan.size3);
  Map_iajk.resize(Chan.size3);

  for(int i = 0; i < Chan.size3; ++i){
    X_ia2[i].assign(Chan.h[i] * Chan.p[i], 0.0); // (i,a)
    X_ia3[i].assign(Chan.p[i] * Chan.h[i], 0.0); // (a,i)

    X_ab2[i].assign(Chan.p[i] * Chan.p[i], 0.0); // (a,b)
    X_ab3[i].assign(Chan.p[i] * Chan.p[i], 0.0); // (b,a)

    X_ij2[i].assign(Chan.h[i] * Chan.h[i], 0.0); // (i,j)
    X_ij3[i].assign(Chan.h[i] * Chan.h[i], 0.0); // (j,i)

    X_ai2[i].assign(Chan.p[i] * Chan.h[i], 0.0); // (a,i)
    X_ai3[i].assign(Chan.h[i] * Chan.p[i], 0.0); // (i,a)

    X1_ij2[i].assign(Chan.h[i] * Chan.h[i], 0.0); // (i,j)
    X1_ij3[i].assign(Chan.h[i] * Chan.h[i], 0.0); // (j,i)

    X1_iabc1[i].assign(Chan.p[i] * Chan.hpp[i], 0.0);
    X1_iabc2[i].assign(Chan.ppp[i] * Chan.h[i], 0.0);
    X1_iabc3[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    X_iabc1[i].assign(Chan.p[i] * Chan.hpp[i], 0.0);
    X_iabc3[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    Map_iabc[i].assign(8 * Chan.p[i] * Chan.hpp[i], 0);
    
    X1_ijka1[i].assign(Chan.h[i] * Chan.hhp[i], 0.0);
    X1_ijka2[i].assign(Chan.hhh[i] * Chan.p[i], 0.0);
    X_ijka1[i].assign(Chan.h[i] * Chan.hhp[i], 0.0);
    Map_ijka[i].assign(6 * Chan.h[i] * Chan.hhp[i], 0);
    
    X1_abcd2[i].assign(Chan.ppp[i] * Chan.p[i], 0.0);
    X1_abcd3[i].assign(Chan.ppp[i] * Chan.p[i], 0.0);
    V_abcd[i].assign(Chan.ppp[i] * Chan.p[i], 0.0);
    
    X_ijkl2[i].assign(Chan.hhh[i] * Chan.h[i], 0.0);
    X_ijkl3[i].assign(Chan.hhh[i] * Chan.h[i], 0.0);
    V_ijkl[i].assign(Chan.hhh[i] * Chan.h[i], 0.0);
    
    X1_iajb2[i].assign(Chan.hhp2[i] * Chan.p[i], 0.0);
    X1_iajb3[i].assign(Chan.hpp2[i] * Chan.h[i], 0.0);
    X1_iajb4[i].assign(Chan.hpp2[i] * Chan.h[i], 0.0);
    X3_iajb2[i].assign(Chan.hhp2[i] * Chan.p[i], 0.0);
    X3_iajb3[i].assign(Chan.hpp2[i] * Chan.h[i], 0.0);
    X3_iajb5[i].assign(Chan.hhp2[i] * Chan.p[i], 0.0);
    X_iajb3[i].assign(Chan.hpp2[i] * Chan.h[i], 0.0);
    
    X_abic1[i].assign(Chan.hpp[i] * Chan.p[i], 0.0);
    X_abic2[i].assign(Chan.ppp[i] * Chan.h[i], 0.0);
    X_abic3[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    X_abic4[i].assign(Chan.hpp2[i] * Chan.p[i], 0.0);
    Map_abic[i].assign(12 * Chan.hpp[i] * Chan.p[i], 0);
    
    X2_iajk1[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
    X2_iajk2[i].assign(Chan.hhh[i] * Chan.p[i], 0.0);
    X2_iajk3[i].assign(Chan.hhp2[i] * Chan.h[i], 0.0);
    X2_iajk4[i].assign(Chan.hhp2[i] * Chan.h[i], 0.0);
    X_iajk1[i].assign(Chan.hhp[i] * Chan.h[i], 0.0);
    Map_iajk[i].assign(12 * Chan.hhp[i] * Chan.h[i], 0.0);
  }

  for(int i = 0; i < Chan.size1; ++i){
    X_ijab1[i].assign(Chan.pp[i] * Chan.hh[i], 0.0);

    X_iabc5[i].assign(Chan.pp[i] * Chan.hp[i], 0.0);

    X_ijka5[i].assign(Chan.hp[i] * Chan.hh[i], 0.0);
    
    X1_abcd1[i].assign(Chan.pp[i] * Chan.pp[i], 0.0);
    X_abcd1[i].assign(Chan.pp[i] * Chan.pp[i], 0.0);
    Map_abcd[i].assign(6 * Chan.pp[i] * Chan.pp[i], 0);
    
    X_ijkl1[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    X_ijkl4[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    Map_ijkl[i].assign(8 * Chan.hh[i] * Chan.hh[i], 0);

    X_abic7[i].resize(Chan.hp[i] * Chan.pp[i], 0.0);

    X2_iajk7[i].resize(Chan.hh[i] * Chan.hp[i], 0.0);
  }

  for(int i = 0; i < Chan.size2; ++i){
    X_iabc4[i].assign(Chan.pp1[i] * Chan.hp1[i], 0.0);
    
    X_ijka4[i].assign(Chan.hh1[i] * Chan.hp1[i], 0.0);
    
    X1_iajb1[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    X3_iajb1[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    X_iajb1[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    Map_iajb[i].assign(8 * Chan.hp2[i] * Chan.hp2[i], 0);
    
    X_abic5[i].assign(Chan.pp1[i] * Chan.hp2[i], 0.0);
    X_abic6[i].assign(Chan.pp1[i] * Chan.hp2[i], 0.0);
    
    X2_iajk5[i].assign(Chan.hh1[i] * Chan.hp2[i], 0.0);
    X2_iajk6[i].assign(Chan.hh1[i] * Chan.hp2[i], 0.0);
  }
 
  State tb;
  int ind, ind1, ind0, hind1, hind2, hind3, hind4, pind1, pind2, pind3, pind4;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    hind1 = Chan.hp1vec1[Chan.ind0][2*i];
    pind1 = Chan.hp1vec1[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    Map_ia[3*i] = ind;
    ind1 = Index11(Chan.hvec1[ind], Chan.pvec1[ind], Chan.h[ind], Chan.p[ind], hind1, pind1);
    Map_ia[3*i + 1] = ind1;
    ind1 = Index11(Chan.pvec1[ind], Chan.hvec1[ind], Chan.p[ind], Chan.h[ind], pind1, hind1);
    Map_ia[3*i + 2] = ind1;
  }
  for(int i = 0; i < Chan.pp1[Chan.ind0]; ++i){
    pind2 = Chan.pp1vec1[Chan.ind0][2*i];
    pind1 = Chan.pp1vec1[Chan.ind0][2*i + 1];
    ind = Chan.indvec[pind1];
    Map_ab[3*i] = ind;
    ind1 = Index11(Chan.pvec1[ind], Chan.pvec1[ind], Chan.p[ind], Chan.p[ind], pind1, pind2);
    Map_ab[3*i + 1] = ind1;
    ind1 = Index11(Chan.pvec1[ind], Chan.pvec1[ind], Chan.p[ind], Chan.p[ind], pind2, pind1);
    Map_ab[3*i + 2] = ind1;
  }
  for(int i = 0; i < Chan.hh1[Chan.ind0]; ++i){
    hind2 = Chan.hh1vec1[Chan.ind0][2*i];
    hind1 = Chan.hh1vec1[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    Map_ij[3*i] = ind;
    ind1 = Index11(Chan.hvec1[ind], Chan.hvec1[ind], Chan.h[ind], Chan.h[ind], hind1, hind2);
    Map_ij[3*i + 1] = ind1;
    ind1 = Index11(Chan.hvec1[ind], Chan.hvec1[ind], Chan.h[ind], Chan.h[ind], hind2, hind1);
    Map_ij[3*i + 2] = ind1;
  }
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    hind1 = Chan.hp1vec1[Chan.ind0][2*i];
    pind1 = Chan.hp1vec1[Chan.ind0][2*i + 1];
    ind = Chan.indvec[hind1];
    Map_ai[3*i] = ind;
    ind1 = Index11(Chan.pvec1[ind], Chan.hvec1[ind], Chan.p[ind], Chan.h[ind], pind1, hind1);
    Map_ai[3*i + 1] = ind1;
    ind1 = Index11(Chan.hvec1[ind], Chan.pvec1[ind], Chan.h[ind], Chan.p[ind], hind1, pind1);
    Map_ai[3*i + 2] = ind1;
  }

  for(int i = 0; i < Chan.size3; ++i){
    for(int p1 = 0; p1 < Chan.p[i]; ++p1){
      pind1 = Chan.pvec1[i][p1];
      for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
	hind1 = Chan.hppvec1[i][3*hpp1];
	pind2 = Chan.hppvec1[i][3*hpp1 + 1];
	pind3 = Chan.hppvec1[i][3*hpp1 + 2];
	ind0 = Chan.hpp[i]*p1 + hpp1;

	ind = Chan.indvec[hind1];
	Map_iabc[i][8*ind0] = ind;
	ind1 = Index31(Chan.pppvec1[ind], Chan.hvec1[ind], Chan.ppp[ind], Chan.h[ind], pind1, pind2, pind3, hind1);
	Map_iabc[i][8*ind0 + 1] = ind1;
	ind = Chan.indvec[pind2];
	Map_iabc[i][8*ind0 + 2] = ind;
	ind1 = Index31(Chan.hpp2vec1[ind], Chan.pvec1[ind], Chan.hpp2[ind], Chan.p[ind], hind1, pind1, pind3, pind2);
	Map_iabc[i][8*ind0 + 3] = ind1;
	minus(tb, Space.qnums[hind1], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iabc[i][8*ind0 + 4] = ind;
	ind1 = Index22(Chan.pp1vec1[ind], Chan.hp1vec1[ind], Chan.pp1[ind], Chan.hp1[ind], pind3, pind1, hind1, pind2);
	Map_iabc[i][8*ind0 + 5] = ind1;
	plus(tb, Space.qnums[pind2], Space.qnums[pind3]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iabc[i][8*ind0 + 6] = ind;
	ind1 = Index22(Chan.ppvec1[ind], Chan.hpvec1[ind], Chan.pp[ind], Chan.hp[ind], pind2, pind3, hind1, pind1);
	Map_iabc[i][8*ind0 + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.h[i]; ++h1){
      hind3 = Chan.hvec1[i][h1];
      for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
	hind1 = Chan.hhpvec1[i][3*hhp1];
	hind2 = Chan.hhpvec1[i][3*hhp1 + 1];
	pind1 = Chan.hhpvec1[i][3*hhp1 + 2];
	ind0 = Chan.hhp[i]*h1 + hhp1;
	
	ind = Chan.indvec[pind1];
	Map_ijka[i][6*ind0] = ind;
	ind1 = Index31(Chan.hhhvec1[ind], Chan.pvec1[ind], Chan.hhh[ind], Chan.p[ind], hind3, hind1, hind2, pind1);
	Map_ijka[i][6*ind0 + 1] = ind1;
	minus(tb, Space.qnums[hind2], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_ijka[i][6*ind0 + 2] = ind;
	ind1 = Index22(Chan.hh1vec1[ind], Chan.hp1vec1[ind], Chan.hh1[ind], Chan.hp1[ind], hind3, hind1, hind2, pind1);
	Map_ijka[i][6*ind0 + 3] = ind1;
	plus(tb, Space.qnums[hind1], Space.qnums[hind2]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_ijka[i][6*ind0 + 4] = ind;
	ind1 = Index22(Chan.hpvec1[ind], Chan.hhvec1[ind], Chan.hp[ind], Chan.hh[ind], hind3, pind1, hind1, hind2);
	Map_ijka[i][6*ind0 + 5] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
      pind3 = Chan.ppvec1[i][2*pp1];
      pind4 = Chan.ppvec1[i][2*pp1 + 1];
      for(int pp2 = 0; pp2 < Chan.pp[i]; ++pp2){
	pind1 = Chan.ppvec1[i][2*pp2];
	pind2 = Chan.ppvec1[i][2*pp2 + 1];
	ind0 = Chan.pp[i]*pp1 + pp2;

	ind = Chan.indvec[pind2];
	Map_abcd[i][6*ind0] = ind;
	ind1 = Index31(Chan.pppvec1[ind], Chan.pvec1[ind], Chan.ppp[ind], Chan.p[ind], pind1, pind3, pind4, pind2);
	Map_abcd[i][6*ind0 + 1] = ind1;
	ind = Chan.indvec[pind1];
	Map_abcd[i][6*ind0 + 2] = ind;
	ind1 = Index31(Chan.pppvec1[ind], Chan.pvec1[ind], Chan.ppp[ind], Chan.p[ind], pind2, pind3, pind4, pind1);
	Map_abcd[i][6*ind0 + 3] = ind1;
	ind = Chan.indvec[pind3];
	Map_abcd[i][6*ind0 + 4] = ind;
	ind1 = Index31(Chan.pppvec1[ind], Chan.pvec1[ind], Chan.ppp[ind], Chan.p[ind], pind4, pind1, pind2, pind3);
	Map_abcd[i][6*ind0 + 5] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size1; ++i){
    for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
      hind1 = Chan.hhvec1[i][2*hh1];
      hind2 = Chan.hhvec1[i][2*hh1 + 1];
      for(int hh2 = 0; hh2 < Chan.hh[i]; ++hh2){
	hind3 = Chan.hhvec1[i][2*hh2];
	hind4 = Chan.hhvec1[i][2*hh2 + 1];
	ind0 = Chan.hh[i]*hh1 + hh2;

	ind = Chan.indvec[hind4];
	Map_ijkl[i][8*ind0] = ind;
	ind1 = Index31(Chan.hhhvec1[ind], Chan.hvec1[ind], Chan.hhh[ind], Chan.h[ind], hind3, hind1, hind2, hind4);
	Map_ijkl[i][8*ind0 + 1] = ind1;
	ind = Chan.indvec[hind3];
	Map_ijkl[i][8*ind0 + 2] = ind;
	ind1 = Index31(Chan.hhhvec1[ind], Chan.hvec1[ind], Chan.hhh[ind], Chan.h[ind], hind4, hind1, hind2, hind3);
	Map_ijkl[i][8*ind0 + 3] = ind1;
	ind = i;
	Map_ijkl[i][8*ind0 + 4] = ind;
	Map_ijkl[i][8*ind0 + 5] = Chan.hh[i]*hh2 + hh1;
	ind = Chan.indvec[hind2];
	Map_ijkl[i][8*ind0 + 6] = ind;
	ind1 = Index31(Chan.hhhvec1[ind], Chan.hvec1[ind], Chan.hhh[ind], Chan.h[ind], hind1, hind3, hind4, hind2);
	Map_ijkl[i][8*ind0 + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size2; ++i){
    for(int hp1 = 0; hp1 < Chan.hp2[i]; ++hp1){
      hind1 = Chan.hp2vec1[i][2*hp1];
      pind2 = Chan.hp2vec1[i][2*hp1 + 1];
      for(int hp2 = 0; hp2 < Chan.hp2[i]; ++hp2){
	hind2 = Chan.hp2vec1[i][2*hp2];
	pind1 = Chan.hp2vec1[i][2*hp2 + 1];
	ind0 = Chan.hp2[i]*hp1 + hp2;

	ind = Chan.indvec[pind1];
	Map_iajb[i][8*ind0] = ind;
	ind1 = Index31(Chan.hhp2vec1[ind], Chan.pvec1[ind], Chan.hhp2[ind], Chan.p[ind], hind1, hind2, pind2, pind1);
	Map_iajb[i][8*ind0 + 1] = ind1;
	ind = Chan.indvec[hind2];
	Map_iajb[i][8*ind0 + 2] = ind;
	ind1 = Index31(Chan.hpp2vec1[ind], Chan.hvec1[ind], Chan.hpp2[ind], Chan.h[ind], hind1, pind1, pind2, hind2);
	Map_iajb[i][8*ind0 + 3] = ind1;
	ind = Chan.indvec[hind1];
	Map_iajb[i][8*ind0 + 4] = ind;
	ind1 = Index31(Chan.hpp2vec1[ind], Chan.hvec1[ind], Chan.hpp2[ind], Chan.h[ind], hind2, pind2, pind1, hind1);
	Map_iajb[i][8*ind0 + 5] = ind1;
	ind = Chan.indvec[pind2];
	Map_iajb[i][8*ind0 + 6] = ind;
	ind1 = Index31(Chan.hhp2vec1[ind], Chan.pvec1[ind], Chan.hhp2[ind], Chan.p[ind], hind2, hind1, pind1, pind2);
	Map_iajb[i][8*ind0 + 7] = ind1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
      hind1 = Chan.hppvec1[i][3*hpp1];
      pind1 = Chan.hppvec1[i][3*hpp1 + 1];
      pind2 = Chan.hppvec1[i][3*hpp1 + 2];
      for(int p1 = 0; p1 < Chan.p[i]; ++p1){
	pind3 = Chan.pvec1[i][p1];
	ind0 = Chan.p[i]*hpp1 + p1;
	
	ind = Chan.indvec[hind1];
	Map_abic[i][12*ind0] = ind;
	ind1 = Index31(Chan.pppvec1[ind], Chan.hvec1[ind], Chan.ppp[ind], Chan.h[ind], pind3, pind1, pind2, hind1);
	Map_abic[i][12*ind0 + 1] = ind1;
	ind = Chan.indvec[pind1];
	Map_abic[i][12*ind0 + 2] = ind;
	ind1 = Index31(Chan.hpp2vec1[ind], Chan.pvec1[ind], Chan.hpp2[ind], Chan.p[ind], hind1, pind3, pind2, pind1);
	Map_abic[i][12*ind0 + 3] = ind1;
	ind = Chan.indvec[pind2];
	Map_abic[i][12*ind0 + 4] = ind;
	ind1 = Index31(Chan.hpp2vec1[ind], Chan.pvec1[ind], Chan.hpp2[ind], Chan.p[ind], hind1, pind3, pind1, pind2);
	Map_abic[i][12*ind0 + 5] = ind1;
	minus(tb, Space.qnums[pind3], Space.qnums[pind2]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*ind0 + 6] = ind;
	ind1 = Index22(Chan.pp1vec1[ind], Chan.hp2vec1[ind], Chan.pp1[ind], Chan.hp2[ind], pind3, pind2, hind1, pind1);
	Map_abic[i][12*ind0 + 7] = ind1;
	minus(tb, Space.qnums[pind3], Space.qnums[pind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_abic[i][12*ind0 + 8] = ind;
	ind1 = Index22(Chan.pp1vec1[ind], Chan.hp2vec1[ind], Chan.pp1[ind], Chan.hp2[ind], pind3, pind1, hind1, pind2);
	Map_abic[i][12*ind0 + 9] = ind1;
	plus(tb, Space.qnums[pind1], Space.qnums[pind2]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_abic[i][12*ind0 + 10] = ind;
	ind1 = Index22(Chan.hpvec1[ind], Chan.ppvec1[ind], Chan.hp[ind], Chan.pp[ind], hind1, pind3, pind1, pind2);
	Map_abic[i][12*ind0 + 11] = ind1;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
      hind2 = Chan.hhpvec1[i][3*hhp1];
      hind3 = Chan.hhpvec1[i][3*hhp1 + 1];
      pind1 = Chan.hhpvec1[i][3*hhp1 + 2];
      for(int h1 = 0; h1 < Chan.h[i]; ++h1){
	hind1 = Chan.hvec1[i][h1];
	ind0 = Chan.h[i]*hhp1 + h1;
	
	ind = Chan.indvec[pind1];
	Map_iajk[i][12*ind0] = ind;
	ind1 = Index31(Chan.hhhvec1[ind], Chan.pvec1[ind], Chan.hhh[ind], Chan.p[ind], hind1, hind2, hind3, pind1);
	Map_iajk[i][12*ind0 + 1] = ind1;
	ind = Chan.indvec[hind3];
	Map_iajk[i][12*ind0 + 2] = ind;
	ind1 = Index31(Chan.hhp2vec1[ind], Chan.hvec1[ind], Chan.hhp2[ind], Chan.h[ind], hind2, hind1, pind1, hind3);
	Map_iajk[i][12*ind0 + 3] = ind1;
	ind = Chan.indvec[hind2];
	Map_iajk[i][12*ind0 + 4] = ind;
	ind1 = Index31(Chan.hhp2vec1[ind], Chan.hvec1[ind], Chan.hhp2[ind], Chan.h[ind], hind3, hind1, pind1, hind2);
	Map_iajk[i][12*ind0 + 5] = ind1;
	minus(tb, Space.qnums[hind2], Space.qnums[hind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*ind0 + 6] = ind;
	ind1 = Index22(Chan.hh1vec1[ind], Chan.hp2vec1[ind], Chan.hh1[ind], Chan.hp2[ind], hind2, hind1, hind3, pind1);
	Map_iajk[i][12*ind0 + 7] = ind1;
	minus(tb, Space.qnums[hind3], Space.qnums[hind1]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	Map_iajk[i][12*ind0 + 8] = ind;
	ind1 = Index22(Chan.hh1vec1[ind], Chan.hp2vec1[ind], Chan.hh1[ind], Chan.hp2[ind], hind3, hind1, hind2, pind1);
	Map_iajk[i][12*ind0 + 9] = ind1;
	plus(tb, Space.qnums[hind2], Space.qnums[hind3]);
	ind = ChanInd_2b_dir(Parameters.basis, Space, tb);
	Map_iajk[i][12*ind0 + 10] = ind;
	ind1 = Index22(Chan.hhvec1[ind], Chan.hpvec1[ind], Chan.hh[ind], Chan.hp[ind], hind2, hind3, hind1, pind1);
	Map_iajk[i][12*ind0 + 11] = ind1;
      }
    }
  }
}

void CC_Eff::set_X_ia(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    x = 0.0;
    x += X_ia1[i];
    x += X_ia2[Map_ia[3*i]][Map_ia[3*i + 1]];
    x += X_ia3[Map_ia[3*i]][Map_ia[3*i + 2]];
    X_ia1[i] = x;
    X_ia2[Map_ia[3*i]][Map_ia[3*i + 1]] = x;
    X_ia3[Map_ia[3*i]][Map_ia[3*i + 2]] = x;
  }
}

void CC_Eff::set_X_ab(const Channels &Chan)
{
  double x;
  for(int i = 0; i < Chan.pp1[Chan.ind0]; ++i){
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
  for(int i = 0; i < Chan.hh1[Chan.ind0]; ++i){
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
  for(int i = 0; i < Chan.hh1[Chan.ind0]; ++i){
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
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
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
    for(int p1 = 0; p1 < Chan.p[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
	ind0 = p1*Chan.hpp[i] + hpp1;
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
    for(int p1 = 0; p1 < Chan.p[i]; ++p1){
      for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
	ind0 = p1*Chan.hpp[i] + hpp1;
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
    for(int h1 = 0; h1 < Chan.h[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
	ind0 = h1*Chan.hhp[i] + hhp1;
	X1_ijka2[Map_ijka[i][6*ind0]][Map_ijka[i][6*ind0 + 1]] = X1_ijka1[i][ind0];
      }
    }
  }
}

void CC_Eff::set_X_ijka(const Channels &Chan)
{
  int ind0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int h1 = 0; h1 < Chan.h[i]; ++h1){
      for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
	ind0 = h1*Chan.hhp[i] + hhp1;
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
    for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
      for(int pp2 = 0; pp2 < Chan.pp[i]; ++pp2){
	x = 0.0;
	ind0 = Chan.pp[i]*pp1 + pp2;
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
    for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
      for(int hh2 = 0; hh2 < Chan.hh[i]; ++hh2){
	x = 0.0;
	ind0 = Chan.hh[i]*hh1 + hh2;
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
    for(int hp1 = 0; hp1 < Chan.hp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.hp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.hp2[i]*hp1 + hp2;
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
    for(int hp1 = 0; hp1 < Chan.hp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.hp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.hp2[i]*hp1 + hp2;
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
    for(int hp1 = 0; hp1 < Chan.hp2[i]; ++hp1){
      for(int hp2 = 0; hp2 < Chan.hp2[i]; ++hp2){
	x = 0.0;
	ind0 = Chan.hp2[i]*hp1 + hp2;
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
    for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
      for(int p1 = 0; p1 < Chan.p[i]; ++p1){
	x = 0.0;
	ind0 = Chan.p[i]*hpp1 + p1;
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
    for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
      for(int h1 = 0; h1 < Chan.h[i]; ++h1){
	x = 0.0;
	ind0 = Chan.h[i]*hhp1 + h1;
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

/*V_Conv::V_Conv(Channels Chan, Input_Parameters Parameters, Model_Space Space, int maxJ)
{
  JME.assign(3*2*(maxJ + 1), 0.0);
  JL.resize(3*maxJ + 2);
  SzTz.resize(9);
  Skhat0.resize(2*Space.CART_tb2Indsize);
  Skhat1.resize(Space.CART_tb2Indsize);

  Y1_JLSk.resize(3);
  Y2_JLSk.resize(3);
  V_JL.resize(3);
  
  for(int i = 0; i < 3; ++i){
    if(i%2 == 1){
      Y1_JLSk[i].assign(2*Space.CART_tb2Indsize * (3*maxJ + 2), 0.0);
      Y2_JLSk[i].assign(2*Space.CART_tb2Indsize * (3*maxJ + 2), 0.0);
    }
    else{
      Y1_JLSk[i].assign(Space.CART_tb2Indsize * (3*maxJ + 2), 0.0);
      Y2_JLSk[i].assign(Space.CART_tb2Indsize * (3*maxJ + 2), 0.0);
    }
  }

  for(int i = 0; i < 3; ++i){ V_JL[i].assign((3*maxJ + 2) * (3*maxJ + 2), 0.0); }
  
  int count1 = 2;
  std::vector<int> tempvec1(2, 0);
  JL[0] = tempvec1;
  tempvec1[1] = 1;
  JL[1] = tempvec1;
  for(int J = 1; J <= maxJ; ++J){
    tempvec1[0] = J;
    for(int L = -1; L <= 1; ++L){
      tempvec1[1] = J + L;
      JL[count1] = tempvec1;
      ++count1;
    }
  }

  count1 = 0;
  for(int Sz = -1; Sz <= 1; ++Sz){
    tempvec1[0] = Sz;
    for(int Tz = -1; Tz <= 1; ++Tz){
      tempvec1[1] = Tz;
      SzTz[count1] = tempvec1;
      ++count1;
    }
  }

  count1 = 0;
  int count2 = 0;
  tempvec1.resize(4);
  int Nx2size = 4*Space.nmax + 1;
  int Ny2size, Nz2size;
  for(int S = 0; S < 2; ++S){
    tempvec1[0] = S;
    for(int i = 0; i < Nx2size; ++i){
      tempvec1[1] = Space.Nx2min + i;
      Ny2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i)))) + 1;
      for(int j = 0; j < Ny2size; ++j){
	tempvec1[2] = Space.Ny2min[i] + j;
	Nz2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j)))) + 1;
	for(int k = 0; k < Nz2size; ++k){
	  tempvec1[3] = Space.Nz2min[i][j] + k;
	  Skhat0[count1] = tempvec1;
	  ++count1;
	  if(S == 1){
	    tempvec1[3] = Space.Nz2min[i][j] + k;
	    Skhat1[count2] = tempvec1;
	    ++count2;
	  }
	}
      }
    }
  }

  }*/


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


void Print_Parameters(const Input_Parameters &Parameters)
{
  std::cout << "----------------------------------------------------------" << std::endl;
  std::cout << "Case = " << Parameters.calc_case << ", Basis = " << Parameters.basis << ", Approximation = " << Parameters.approx << std::endl;
  if(Parameters.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << Parameters.LevelScheme << std::endl;
    if(Parameters.MatrixElements.size() > 0){ std::cout << "Interaction = " << Parameters.MatrixElements << std::endl; }
  }
  else{ std::cout << "Nmax = " << Parameters.Nmax << ", Density = " << Parameters.density << std::endl; }
  std::cout << "OB strength = " << Parameters.obstrength << ", TB strength = " << Parameters.tbstrength << std::endl;
  if(Parameters.calc_case == "nuclear"){
    std::cout << "Proton Shells = " << Parameters.Pshells << ", Neutron Shells = " << Parameters.Nshells << std::endl;
    std::cout << "Protons = " << Parameters.P << ", Neutrons = " << Parameters.N << std::endl;
  }
  else if(Parameters.calc_case == "electronic"){
    std::cout << "Electron Shells = " << Parameters.Pshells << " Electrons = " << Parameters.P << std::endl;
  }
  std::cout << "----------------------------------------------------------" << std::endl;
}


void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  std::string fullpath; // Model Space file path
  std::string phline; // Std::String for each file line
  std::ifstream splevels; // Model space file
  int ind, l;
  double energy; // initialize level index, n, l, nx, ny, nz, m, tz depending on basis
  int pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;
  
  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  getline(splevels, phline);
  std::stringstream(phline) >> Space.indtot;

  Space.qnums.resize(Space.indtot);

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

  //tz,ms,nx,ny,nz,ml,n,j,par
  for(int i = 0; i < Space.indtot; ++i){
    getline(splevels, phline);
    if(Parameters.basis == "infinite"){ // ind, nx, ny, nz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].nx >> Space.qnums[i].ny >> Space.qnums[i].nz >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
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
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
      std::cout << Space.qnums[i].n << " " << Space.qnums[i].par << " " << Space.qnums[i].j << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " " << energy << std::endl;
    }
    else if(Parameters.basis == "finite_HO"){ // ind, n, l, lz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].ml >> Space.qnums[i].m >> Space.qnums[i].t >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
      if(Space.qnums[i].ml < Space.qmins.ml){ Space.qmins.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].ml > Space.qmaxs.ml){ Space.qmaxs.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_J"){ // ind, n, l, j, tz, l2n
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> l >> Space.qnums[i].j >> Space.qnums[i].t >> ind >> energy;
      Space.qnums[i].par = -2*(l%2) + 1;
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

  std::vector<int> ShellInd_P;
  std::vector<int> ShellInd_N;
  double Etemp_P = -1000.0;
  double Etemp_N = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > Etemp_P && Space.qnums[i].t == -1){
      ShellInd_P.push_back(i);
      Etemp_P = Space.qnums[i].energy;
    }
    if(Space.qnums[i].energy > Etemp_N && Space.qnums[i].t == 1){
      ShellInd_N.push_back(i);
      Etemp_N = Space.qnums[i].energy;
    }
  }
  std::cout << Parameters.Pshells << " " << int(ShellInd_P.size()) << std::endl;
  std::cout << Parameters.Nshells << " " << int(ShellInd_N.size()) << std::endl;
  if(Parameters.Pshells > int(ShellInd_P.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > int(ShellInd_N.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < ShellInd_P[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      if(i < ShellInd_N[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++ncount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].par << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
  }

  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.j = Space.qmaxs.j + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;

  //std::cout << "$ " << Space.qsizes.t << " " << Space.qsizes.m << " " << Space.qsizes.par << std::endl;

  if(Parameters.basis == "infinite"){
    int N2max = 0;
    for(int i = 0; i < Space.indtot; ++i){
      for(int j = i+1; j < Space.indtot; ++j){
	if((Space.qnums[i].nx + Space.qnums[j].nx)*(Space.qnums[i].nx + Space.qnums[j].nx) + 
	   (Space.qnums[i].ny + Space.qnums[j].ny)*(Space.qnums[i].ny + Space.qnums[j].ny) + 
	   (Space.qnums[i].nz + Space.qnums[j].nz)*(Space.qnums[i].nz + Space.qnums[j].nz) > N2max){
	  N2max = (Space.qnums[i].nx + Space.qnums[j].nx)*(Space.qnums[i].nx + Space.qnums[j].nx) + 
	    (Space.qnums[i].ny + Space.qnums[j].ny)*(Space.qnums[i].ny + Space.qnums[j].ny) + 
	    (Space.qnums[i].nz + Space.qnums[j].nz)*(Space.qnums[i].nz + Space.qnums[j].nz); }
      }
    }
    Space.map_2b.assign(Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz, -1);
    int count1 = 0;
    for(int nx = 2*Space.qmins.nx; nx <= 2*Space.qmaxs.nx; ++nx){
      for(int ny = 2*Space.qmins.ny; ny <= 2*Space.qmaxs.ny; ++ny){
	for(int nz = 2*Space.qmins.nz; nz <= 2*Space.qmaxs.nz; ++nz){
	  if(nx*nx + ny*ny + nz*nz <= N2max){
	    Space.map_2b[(nx - 2*Space.qmins.nx)*Space.qsizes.ny*Space.qsizes.nz + (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = count1;
	    ++count1;
	  }
	}
      }
    }
    Space.Chansize_2b = count1 * Space.qsizes.t * Space.qsizes.m;
  }
  else if(Parameters.basis == "finite_M"){
    Space.Chansize_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_HO"){
    Space.Chansize_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;
  }
  else if(Parameters.basis == "finite_J"){
    Space.Chansize_2b = Space.qsizes.par * Space.qsizes.j * Space.qsizes.t;
  }
}

void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space)
{
  int shelllength; // number of single particle states for each shell
  int pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;
  int ind;
  double m;

  std::vector<struct State> states = Space.qnums;
  Space.shellsm.resize(states.size());
  Space.indtot = 0;
  for(int i = 0; i < int(Space.qnums.size()); ++i){ Space.indtot += Space.qnums[i].j + 1; }

  Space.qnums.resize(Space.indtot);

  Space.qmins.t = 1000;
  Space.qmins.m = 1000;
  Space.qmins.j = 1000;
  Space.qmins.par = 1000;

  Space.qmaxs.t = -1000;
  Space.qmaxs.m = -1000;
  Space.qmaxs.j = -1000;
  Space.qmaxs.par = -1000;

  ind = 0;
  for(int i = 0; i < int(states.size()); ++i){
    shelllength = states[i].j + 1;
    Space.shellsm[i].resize(shelllength);
    for(int j = 0; j < shelllength; ++j){
      m = -states[i].j + 2*j;
      Space.qnums[ind].par = states[i].par;
      Space.qnums[ind].j = states[i].j;
      Space.qnums[ind].m = m;
      Space.qnums[ind].t = states[i].t;
      Space.qnums[ind].energy = states[i].energy;
      Space.qnums[ind].type = states[i].type;
      Space.shellsm[i][j] = ind;
      if(Space.qnums[i].par < Space.qmins.par){ Space.qmins.par = Space.qnums[i].par; }
      if(Space.qnums[i].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[i].par; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
      ++ind;
    }
  }

  std::vector<int> ShellInd_P;
  std::vector<int> ShellInd_N;
  double Etemp_P = -1000.0;
  double Etemp_N = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > Etemp_P && Space.qnums[i].t == -1){
      ShellInd_P.push_back(i);
      Etemp_P = Space.qnums[i].energy;
    }
    if(Space.qnums[i].energy > Etemp_N && Space.qnums[i].t == 1){
      ShellInd_N.push_back(i);
      Etemp_N = Space.qnums[i].energy;
    }
  }
  if(Parameters.Pshells > int(ShellInd_P.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > int(ShellInd_N.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < ShellInd_P[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      if(i < ShellInd_N[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++ncount; ++parcount; }
    }
  }

  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.j = Space.qmaxs.j + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;
  Space.Chansize_2b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
}


void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  int count = 0, pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;

  double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
  double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);

  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ 
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2*2;
    Space.qmins.t = -1, Space.qmaxs.t = 1, Space.qsizes.t = 2;
  }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1;
  }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.qmins.t = 1, Space.qmaxs.t = 1, Space.qsizes.t = 1;
  }
  else{ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }

  Space.qnums.resize(Space.indtot);
  int nmax = std::floor(std::sqrt(Parameters.Nmax));
  for(int shell = 0; shell <= Parameters.Nmax; ++shell){
    for(int nx = -nmax; nx <= nmax; ++nx){    
      for(int ny = -nmax; ny <= nmax; ++ny){	
	for(int nz = -nmax; nz <= nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz || Parameters.Nmax <= nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1){
		if(Parameters.Pshells == 0){ continue; }
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		E = 4.0*(nx*nx + ny*ny + nz*nz);
		Space.qnums[count].energy = E;
	      }
	      Space.qnums[count].nx = nx;
	      Space.qnums[count].ny = ny;
	      Space.qnums[count].nz = nz;
	      Space.qnums[count].m = sz;
	      Space.qnums[count].t = tz;
	      if(Space.qnums[count].nx < Space.qmins.nx){ Space.qmins.nx = Space.qnums[count].nx; }
	      if(Space.qnums[count].nx > Space.qmaxs.nx){ Space.qmaxs.nx = Space.qnums[count].nx; }
	      if(Space.qnums[count].ny < Space.qmins.ny){ Space.qmins.ny = Space.qnums[count].ny; }
	      if(Space.qnums[count].ny > Space.qmaxs.ny){ Space.qmaxs.ny = Space.qnums[count].ny; }
	      if(Space.qnums[count].nz < Space.qmins.nz){ Space.qmins.nz = Space.qnums[count].nz; }
	      if(Space.qnums[count].nz > Space.qmaxs.nz){ Space.qmaxs.nz = Space.qnums[count].nz; }
	      if(Space.qnums[count].m < Space.qmins.m){ Space.qmins.m = Space.qnums[count].m; }
	      if(Space.qnums[count].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[count].m; }
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.indtot = count;
  Space.qnums.resize(count);
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;

  std::vector<int> ShellInd_P;
  std::vector<int> ShellInd_N;
  double Etemp_P = -1000.0;
  double Etemp_N = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > Etemp_P && Space.qnums[i].t == -1){
      ShellInd_P.push_back(i);
      Etemp_P = Space.qnums[i].energy;
    }
    if(Space.qnums[i].energy > Etemp_N && Space.qnums[i].t == 1){
      ShellInd_N.push_back(i);
      Etemp_N = Space.qnums[i].energy;
    }
  }
  if(Parameters.Pshells >= int(ShellInd_P.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells >= int(ShellInd_N.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < ShellInd_P[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.qnums[i].t == 1){
      if(i < ShellInd_N[Parameters.Nshells]){ Space.qnums[i].type = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.qnums[i].type = "particle"; ++ncount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  double L = pow(holcount/Parameters.density, 1.0/3.0);
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac*M_PI*M_PI/(L*L); }
    else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac*M_PI*M_PI/(L*L); }
  }

  int count1 = 0;
  int N2max = 4*Parameters.Nmax;
  Space.map_2b.assign(Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz, -1);
  for(int nx = 2*Space.qmins.nx; nx <= 2*Space.qmaxs.nx; ++nx){
    for(int ny = 2*Space.qmins.ny; ny <= 2*Space.qmaxs.ny; ++ny){
      for(int nz = 2*Space.qmins.nz; nz <= 2*Space.qmaxs.nz; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max){
	  Space.map_2b[(nx - 2*Space.qmins.nx)*Space.qsizes.ny*Space.qsizes.nz + (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = count1;
	  ++count1;
	}
      }
    }
  }
  Space.Chansize_2b = count1 * Space.qsizes.t * Space.qsizes.m;
}


void EG_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  int count = 0, pcount = 0, holcount = 0, parcount = 0, phcount = 0;
  Parameters.Nshells = 0;

  double prefac = hbarc_HartA * hbarc_HartA / (2.0 * m_electronc2_Hart);
  double r_b = hbarc_HartA / (m_electronc2_Hart * fine_struct);
  Parameters.density = 3.0/(4.0 * PI * pow(Parameters.density * r_b, 3));

  if(Parameters.Pshells != 0){ 
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1;
  }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  Space.qnums.resize(Space.indtot);

  int nmax = std::floor(std::sqrt(Parameters.Nmax));
  for(int shell = 0; shell <= Parameters.Nmax; ++shell){
    for(int nx = -nmax; nx <= nmax; ++nx){    
      for(int ny = -nmax; ny <= nmax; ++ny){	
	for(int nz = -nmax; nz <= nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz || Parameters.Nmax <= nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    E = 4.0*(nx*nx + ny*ny + nz*nz);
	    Space.qnums[count].energy = E;
	    Space.qnums[count].nx = nx;
	    Space.qnums[count].ny = ny;
	    Space.qnums[count].nz = nz;
	    Space.qnums[count].m = sz;
	    Space.qnums[count].t = -1;
	    if(Space.qnums[count].nx < Space.qmins.nx){ Space.qmins.nx = Space.qnums[count].nx; }
	    if(Space.qnums[count].nx > Space.qmaxs.nx){ Space.qmaxs.nx = Space.qnums[count].nx; }
	    if(Space.qnums[count].ny < Space.qmins.ny){ Space.qmins.ny = Space.qnums[count].ny; }
	    if(Space.qnums[count].ny > Space.qmaxs.ny){ Space.qmaxs.ny = Space.qnums[count].ny; }
	    if(Space.qnums[count].nz < Space.qmins.nz){ Space.qmins.nz = Space.qnums[count].nz; }
	    if(Space.qnums[count].nz > Space.qmaxs.nz){ Space.qmaxs.nz = Space.qnums[count].nz; }
	    if(Space.qnums[count].m < Space.qmins.m){ Space.qmins.m = Space.qnums[count].m; }
	    if(Space.qnums[count].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[count].m; }
	    count++;
	  }
	}   
      } 
    }
  }
  Space.indtot = count;
  Space.qnums.resize(count);
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.nx = 2*(Space.qmaxs.nx - Space.qmins.nx) + 1;
  Space.qsizes.ny = 2*(Space.qmaxs.ny - Space.qmins.ny) + 1;
  Space.qsizes.nz = 2*(Space.qmaxs.nz - Space.qmins.nz) + 1;

  std::vector<int> ShellInd_P;
  double Etemp_P = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > Etemp_P && Space.qnums[i].t == -1){
      ShellInd_P.push_back(i);
      Etemp_P = Space.qnums[i].energy;
    }
  }
  if(Parameters.Pshells >= int(ShellInd_P.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < ShellInd_P[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indn = 0;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = 0;

  double L = pow(holcount/Parameters.density, 1.0/3.0);
  for(int i = 0; i < Space.indtot; ++i){ Space.qnums[i].energy *= prefac * PI * PI / (L * L); }

  int count1 = 0;
  int N2max = 4*Parameters.Nmax;
  Space.map_2b.assign(Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz, -1);
  for(int nx = 2*Space.qmins.nx; nx <= 2*Space.qmaxs.nx; ++nx){
    for(int ny = 2*Space.qmins.ny; ny <= 2*Space.qmaxs.ny; ++ny){
      for(int nz = 2*Space.qmins.nz; nz <= 2*Space.qmaxs.nz; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max){
	  Space.map_2b[(nx - 2*Space.qmins.nx)*Space.qsizes.ny*Space.qsizes.nz + (ny - 2*Space.qmins.ny)*Space.qsizes.nz + (nz - 2*Space.qmins.nz)] = count1;
	  ++count1;
	}
      }
    }
  }
  Space.Chansize_2b = count1 * Space.qsizes.t * Space.qsizes.m;
}


void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double E;
  int count = 0, pcount = 0, holcount = 0, parcount = 0, phcount = 0;
  Parameters.Nshells = 0;

  Space.qmins.m = 1000;
  Space.qmins.ml = 1000;
  Space.qmins.par = 1000;

  Space.qmaxs.m = -1000;
  Space.qmaxs.ml = -1000;
  Space.qmaxs.par = -1000;

  if(Parameters.Pshells != 0){ 
    Space.indtot = int(Parameters.Nmax * (Parameters.Nmax + 1));
    Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1;
  }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  Space.qnums.resize(Space.indtot);

  for(int shell = 0; shell < Parameters.Nmax; ++shell){
    for(int n = 0; n < Parameters.Nmax; ++n){
      for(int ml = -1 * int(shell - 2*n); ml <= int(shell - 2*n); ++ml){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  E = Parameters.density*(shell + 1.0);
	  Space.qnums[count].energy = E;
	  Space.qnums[count].n = n;
	  Space.qnums[count].par = -2*(abs(ml)%2) + 1;
	  Space.qnums[count].ml = ml;
	  Space.qnums[count].m = sz;
	  Space.qnums[count].t = -1;
	  if(Space.qnums[count].par < Space.qmins.par){ Space.qmins.par = Space.qnums[count].par; }
	  if(Space.qnums[count].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[count].par; }
	  if(Space.qnums[count].m < Space.qmins.m){ Space.qmins.m = Space.qnums[count].m; }
	  if(Space.qnums[count].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[count].m; }
	  if(Space.qnums[count].ml < Space.qmins.ml){ Space.qmins.ml = Space.qnums[count].ml; }
	  if(Space.qnums[count].ml > Space.qmaxs.ml){ Space.qmaxs.ml = Space.qnums[count].ml; }
	  if(Space.qnums[count].t < Space.qmins.t){ Space.qmins.t = Space.qnums[count].t; }
	  if(Space.qnums[count].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[count].t; }
	  count++;
	}
      }
    }
  }
  std::vector<int> ShellInd_P;
  double Etemp_P = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].energy > Etemp_P && Space.qnums[i].t == -1){
      ShellInd_P.push_back(i);
      Etemp_P = Space.qnums[i].energy;
    }
  }
  if(Parameters.Pshells > int(ShellInd_P.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].t == -1){
      if(i < ShellInd_P[Parameters.Pshells]){ Space.qnums[i].type = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.qnums[i].type = "particle"; ++pcount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = 0;

  for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].par << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
  }

  Space.qsizes.t = Space.qmaxs.t - Space.qmins.t + 1;
  Space.qsizes.m = Space.qmaxs.m - Space.qmins.m + 1;
  Space.qsizes.ml = 2*(Space.qmaxs.ml - Space.qmins.ml) + 1;
  Space.qsizes.par = int((Space.qmaxs.par - Space.qmins.par)/2) + 1;
  Space.Chansize_2b = Space.qsizes.ml * Space.qsizes.m * Space.qsizes.t;
}


/*double Partial_Wave_HF(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const double &maxL)
{

  for(int L = 0; L < maxL; ++L){
    for(int J1 = abs(L-1); J1 <= L+1; ++J1){
      for(int J2 = abs(L-1); J2 <= L+1; ++J2){
	for(int M1 = -J1; M1 <= J1; ++M1){
	  for(int M2 = -J2; M2 <= J2; ++M2){
	    for(int S1 = abs(J1-L); S1 >= 0; --S1){
	      for(int S2 = abs(J2-L); S2 >= 0; --S2){
		
	      }
	    }
	  }
	}
      }
    }
  }
	  
  return maxL;
  }*/

/*Model_Space CART_Build_Model_Space_Twist(Input_Parameters &Parameters, const double &thx, const double &thy, const double &thz)
{
  Model_Space Space; // Model Space information
  double E;
  int count, pcount, ncount, holcount, parcount, phcount, nhcount;
  
  double hbarc = 197.3269788; // MeVfm
  double m_neutronc2 = 939.565378; // MeV
  //double m_protonc2 = 938.272046; // MeV
  double m_protonc2 = 939.565378; // MeV
  double neutron_prefac = hbarc*hbarc/(2.0*m_neutronc2);
  double proton_prefac = hbarc*hbarc/(2.0*m_protonc2);
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1.0/3.0);

  if(Parameters.P != 0 && Parameters.N != 0){ Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2*2; }
  else if(Parameters.P != 0 || Parameters.N != 0){ Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2; }
  else{ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);
  Space.levelsnx.resize(Space.indtot);
  Space.levelsny.resize(Space.indtot);
  Space.levelsnz.resize(Space.indtot);

  count = 0;
  pcount = 0;
  ncount = 0;
  parcount = 0;
  holcount = 0;
  phcount = 0;
  nhcount = 0;

  int nmax = std::ceil(std::sqrt(Parameters.Nmax));
  for(int nx = -nmax; nx <= nmax; ++nx){    
    for(int ny = -nmax; ny <= nmax; ++ny){	
      for(int nz = -nmax; nz <= nmax; ++nz){	  
	if(Parameters.Nmax < nx*nx + ny*ny + nz*nz){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  for( int tz = -1; tz <= 1; tz = tz+2){
	    if(tz == -1){
	      if(Parameters.P == 0){ continue; }
	      E = 4.0*(nx*nx + ny*ny + nz*nz) + (4.0/M_PI)*(nx*thx + ny*thy + nz*thz) + (1.0/(M_PI*M_PI))*(thx*thx + thy*thy + thz*thz);
	      Space.levelsen[count] = E;
	    }
	    if(tz == 1){
	      if(Parameters.N == 0){ continue; }
	      E = 4.0*(nx*nx + ny*ny + nz*nz) + (4.0/M_PI)*(nx*thx + ny*thy + nz*thz) + (1.0/(M_PI*M_PI))*(thx*thx + thy*thy + thz*thz);
	      Space.levelsen[count] = E;
	    }
	    Space.levelsind[count] = count + 1;
	    Space.levelsnx[count] = nx;
	    Space.levelsny[count] = ny;
	    Space.levelsnz[count] = nz;
	    Space.levelsm[count] = sz;
	    Space.levelst[count] = tz;
	    count++;
	  }
	}   
      } 
    }
  }

  //order states by energy
  double tempE;
  int tempind;
  int nx1, ny1, nz1;
  double sz1, tz1, en1;
  for(int i = 0; i < count; ++i){
    tempE = Space.levelsen[i];
    tempind = i;
    for(int j = i+1; j < count; ++j){
      if(Space.levelsen[j] < tempE){ tempE = Space.levelsen[j]; tempind = j; }
    }
    if(tempind != i){
      nx1 = Space.levelsnx[tempind];
      ny1 = Space.levelsny[tempind];
      nz1 = Space.levelsnz[tempind];
      sz1 = Space.levelsm[tempind];
      tz1 = Space.levelst[tempind];
      en1 = Space.levelsen[tempind];
      for(int j = tempind; j > i; --j){
	Space.levelsnx[j] = Space.levelsnx[j-1];
	Space.levelsny[j] = Space.levelsny[j-1];
	Space.levelsnz[j] = Space.levelsnz[j-1];
	Space.levelsm[j] = Space.levelsm[j-1];
	Space.levelst[j] = Space.levelst[j-1];
	Space.levelsen[j] = Space.levelsen[j-1];
      }
      Space.levelsnx[i] = nx1;
      Space.levelsny[i] = ny1;
      Space.levelsnz[i] = nz1;
      Space.levelsm[i] = sz1;
      Space.levelst[i] = tz1;
      Space.levelsen[i] = en1;
    }
  }

  for(int i = 0; i < count; ++i){
    if(Space.levelst[i] == -1){
      if(phcount < Parameters.P){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.levelst[i] == 1){
      if(nhcount < Parameters.N){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
    }
  }

  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;

  for(int i = 0; i < count; ++i){
    if(Space.levelst[i] == -1){ Space.levelsen[i] *= proton_prefac*M_PI*M_PI/(L*L); }
    else if(Space.levelst[i] == 1){ Space.levelsen[i] *= neutron_prefac*M_PI*M_PI/(L*L); }
  }

  Space.levelsind.resize(count);
  Space.levelsm.resize(count);
  Space.levelst.resize(count);
  Space.levelstype.resize(count);
  Space.levelsen.resize(count);
  Space.levelsnx.resize(count);
  Space.levelsny.resize(count);
  Space.levelsnz.resize(count);
  Space.indtot = count;

  Space.nmax = int(std::floor(std::sqrt(Parameters.Nmax)));

  int count1 = 0;
  int Nxsize = 4*Space.nmax + 1;
  int Nysize, Nzsize;
  Space.Nxmin = -2*Space.nmax;
  Space.Nymin.resize(Nxsize);
  Space.Nzmin.resize(Nxsize);
  Space.CART_tb1Indvec.resize(Nxsize);
  for(int i = 0; i < Nxsize; ++i){
    Space.Nymin[i] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i))));
    Nysize = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i)))) + 1;
    Space.Nzmin[i].resize(Nysize);
    Space.CART_tb1Indvec[i].resize(Nysize);
    for(int j = 0; j < Nysize; ++j){
      Space.Nzmin[i][j] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j))));
      Nzsize = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j)))) + 1;
      Space.CART_tb1Indvec[i][j].resize(Nzsize);
      for(int k = 0; k < Nzsize; ++k){
	Space.CART_tb1Indvec[i][j][k] = count1;
	++count1;
      }
    }
  }
  Space.CART_tb1Indsize = count1;

  int count2 = 0;
  int Nx2size = 4*Space.nmax + 1;
  int Ny2size, Nz2size;
  Space.Nx2min = -2*Space.nmax;
  Space.Ny2min.resize(Nx2size);
  Space.Nz2min.resize(Nx2size);
  Space.CART_tb2Indvec.resize(Nx2size);
  for(int i = 0; i < Nx2size; ++i){
    Space.Ny2min[i] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i))));
    Ny2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i)))) + 1;
    Space.Nz2min[i].resize(Ny2size);
    Space.CART_tb2Indvec[i].resize(Ny2size);
    for(int j = 0; j < Ny2size; ++j){
      Space.Nz2min[i][j] = -int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j))));
      Nz2size = 2*int(std::floor(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j)))) + 1;
      Space.CART_tb2Indvec[i][j].resize(Nz2size);
      for(int k = 0; k < Nz2size; ++k){
	Space.CART_tb2Indvec[i][j][k] = count2;
	++count2;
      }
    }
  }
  Space.CART_tb2Indsize = count2;
  
  Space.Mmin = -2;
  Space.Msize = 3;
  Space.M2min = -2;
  Space.M2size = 3;
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.Tmin = -2;
    Space.Tsize = 3;
    Space.T2min = -2;
    Space.T2size = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.Tmin = -2;
    Space.Tsize = 1;
    Space.T2min = 0;
    Space.T2size = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.Tmin = 2;
    Space.Tsize = 1;
    Space.T2min = 0;
    Space.T2size = 1;
  }

  return Space;
}*/


int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx - 2*Space.qmins.nx)*Space.qsizes.ny*Space.qsizes.nz + (State.ny - 2*Space.qmins.ny)*Space.qsizes.nz +
			(State.nz - 2*Space.qmins.nz)]*Space.qsizes.m*Space.qsizes.t + int((State.m - 2*Space.qmins.m)/2)*Space.qsizes.t + 
                         int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - 2*Space.qmins.m)/2)*Space.qsizes.t + 
           int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - 2*Space.qmins.ml)*Space.qsizes.m*Space.qsizes.t +
      int((State.m - 2*Space.qmins.m)/2)*Space.qsizes.t + int((State.t - 2*Space.qmins.t)/2);
  }
  else if(basis == "finite_J"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - 2*Space.qmins.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State)
{
  if(basis == "infinite"){
    return Space.map_2b[(State.nx - Space.qmins.nx+Space.qmaxs.nx)*Space.qsizes.ny*Space.qsizes.nz + (State.ny - Space.qmins.ny+Space.qmaxs.ny)*Space.qsizes.nz +
			(State.nz - Space.qmins.nz+Space.qmaxs.nz)]*Space.qsizes.m*Space.qsizes.t + int((State.m - Space.qmins.m+Space.qmaxs.m)/2)*Space.qsizes.t + 
                         int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_M"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.m*Space.qsizes.t + int((State.m - Space.qmins.m+Space.qmaxs.m)/2)*Space.qsizes.t + 
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_HO"){
    return int(State.ml - Space.qmins.ml+Space.qmaxs.ml)*Space.qsizes.m*Space.qsizes.t + int((State.m - Space.qmins.m+Space.qmaxs.m)/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else if(basis == "finite_J"){
    return int((State.par - Space.qmins.par*Space.qmaxs.par)/2)*Space.qsizes.j*Space.qsizes.t + int(State.j/2)*Space.qsizes.t +
           int((State.t - Space.qmins.t+Space.qmaxs.t)/2);
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}


void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  double L = pow(Parameters.P/Parameters.density, 1.0/3.0);

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.pp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec1[i][2*pq];
      shell2 = Chan.ppvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.pp[i]; ++rs){
	shell3 = Chan.ppvec1[i][2*rs];
	shell4 = Chan.ppvec1[i][2*rs + 1];
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V1[i][Chan.pp[i]*pq + rs] = TBME;
	Ints.D_ME1.V1[i][Chan.pp[i]*rs + pq] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.hh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec1[i][2*pq];
      shell2 = Chan.hhvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.hh[i]; ++rs){
	shell3 = Chan.hhvec1[i][2*rs];
	shell4 = Chan.hhvec1[i][2*rs + 1];
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V2[i][Chan.hh[i]*pq + rs] = TBME;
	Ints.D_ME1.V2[i][Chan.hh[i]*rs + pq] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.pp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec1[i][2*pq];
      shell2 = Chan.ppvec1[i][2*pq + 1];
      for(int rs = 0; rs < Chan.hh[i]; ++rs){
	shell3 = Chan.hhvec1[i][2*rs];
	shell4 = Chan.hhvec1[i][2*rs + 1];
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V4[i][Chan.hh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int q = 0; q < Chan.h[i]; ++q){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell2 = Chan.hvec1[i][q];
      for(int prs = 0; prs < Chan.hpp[i]; ++prs){
	shell1 = Chan.hppvec1[i][3*prs];
	shell3 = Chan.hppvec1[i][3*prs + 1];
	shell4 = Chan.hppvec1[i][3*prs + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V5[i][Chan.hpp[i]*q + prs] = TBME;
	Ints.D_ME1.V6[i][Chan.hpp[i]*q + prs] = -1.0 * TBME;
      }
    }
    #pragma omp parallel for
    for(int s = 0; s < Chan.p[i]; ++s){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell4 = Chan.pvec1[i][s];
      for(int pqr = 0; pqr < Chan.hhp[i]; ++pqr){
	shell1 = Chan.hhpvec1[i][3*pqr];
	shell2 = Chan.hhpvec1[i][3*pqr + 1];
	shell3 = Chan.hhpvec1[i][3*pqr + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V7[i][Chan.hhp[i]*s + pqr] = TBME;
	Ints.D_ME1.V8[i][Chan.hhp[i]*s + pqr] = -1.0 * TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size2; ++i){
    #pragma omp parallel for
    for(int ps = 0; ps < Chan.hp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*ps];
      shell4 = Chan.hp2vec1[i][2*ps + 1];
      for(int qr = ps; qr < Chan.hp2[i]; ++qr){
	shell3 = Chan.hp2vec1[i][2*qr];
	shell2 = Chan.hp2vec1[i][2*qr + 1];
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V3[i][Chan.hp2[i]*ps + qr] = TBME;
	Ints.D_ME1.V3[i][Chan.hp2[i]*qr + ps] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pr = 0; pr < Chan.hp2[i]; ++pr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*pr];
      shell3 = Chan.hp2vec1[i][2*pr + 1];
      for(int qs = 0; qs < Chan.hp1[i]; ++qs){
	shell2 = Chan.hp1vec1[i][2*qs];
	shell4 = Chan.hp1vec1[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = Coulomb_Inf(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V9[i][Chan.hp1[i]*pr + qs] = TBME;
	Ints.D_ME1.V10[i][Chan.hp1[i]*pr + qs] = -1.0 * TBME;
      }
    }
  }
}


// Minnesota Potential for momentum basis
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
    if(qSquared1 < 1.0e-8){ term += 0.0; }
    else{ term += prefactor*4.0*M_PI/qSquared1; }
  }
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-8){ term -= 0.0; }
    else{ term -= prefactor*4.0*M_PI/qSquared1; }
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
	    //std::cout << "n:m:j:g = " << n1 << " " << n2 << " " << n3 << " " << n4 << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << " : " << j1 << " " << j2 << " " << j3 << " " << j4 << " : " << g1 << " " << g2 << " " << g3 << " " << g4 << std::endl;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    //std::cout << "logratio1(j1,j2,j3,j4) = " << LogRatio1 << std::endl;
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    //std::cout << "logprod2(n1,m1,n2,m2,n3,m3,n4,m4,j1,j2,j3,j4) = " << LogProd2 << std::endl;
	    double LogRatio2 = logratio2(G);
	    //std::cout << "G, logratio2(G) = " << G << " " << LogRatio2 << std::endl;
	    //std::cout << "1: " << j1 << " " << j2 << " " << j3 << " " << j4 << ", " << LogRatio1 << " " << LogProd2 << " " << LogRatio2 << std::endl;
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    //if(l1 != l2 || l3 != l4){ continue; }
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    //std::cout << "l = " << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1) * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		    //std::cout << "logproduct3(l1,l2,l3,l4,g1,g2,g3,g4) = " << logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) << std::endl;
		    //std::cout << "lgamma(1 + L/2) = " << lgamma(1.0 + 0.5*L) << std::endl;
		    //std::cout << "lgamma((G - L + 1)/2) = " << lgamma(0.5*(G - L + 1.0)) << std::endl;
		    //std::cout << "temp = " << temp << std::endl;
		    //std::cout << "g2,g3,l2,l3 = " << g2 << " " << g3 << " " << l2 << " " << l3 << ", " << g2+g3-l2-l3 << " " << (-2*((g2 + g3 - l2 - l3)%2) + 1) << std::endl;
		    //std::cout << "2: " << l1 << " " << l2 << " " << l3 << " " << l4 << ", " << logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) << " " << lgamma(1.0 + 0.5*L) << " " << lgamma(0.5*(G - L + 1.0)) << " " << (-2*((g2 + g3 - l2 - l3)%2) + 1) * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0))) << std::endl;
		  }
		}
	      }
	    }
	    dir += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	    //std::cout << "3: " << (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp << std::endl;
	  }
	}
      }
    }
    dir *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
    //std::cout << "! " << product1(n1, m1, n2, m2, n3, m3, n4, m4) << ", " << dir << std::endl;
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
	    //std::cout << n1 << " " << m1 << ", " << n2 << " " << m2 << ", " << n3 << " " << m3 << " : " << j1 << " " << j2 << " " << j3 << " " << j4 << " : " << g1 << " " << g2 << " " << g3 << " " << g4 << std::endl;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    double LogRatio2 = logratio2(G);
	    //std::cout << "4: " << j1 << " " << j2 << " " << j3 << " " << j4 << ", " << LogRatio1 << " " << LogProd2 << " " << LogRatio2 << std::endl;
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    //if(l1 != l2 || l3 != l4){ continue; }
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1) * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		    //std::cout << "5: " << l1 << " " << l2 << " " << l3 << " " << l4 << ", " << logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) << " " << lgamma(1.0 + 0.5*L) << " " << lgamma(0.5*(G - L + 1.0)) << " " << (-2*((g2 + g3 - l2 - l3)%2) + 1) * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0))) << std::endl;
		  }
		}
	      }
	    }
	    exch += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	    //std::cout << "6: " << (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp << std::endl;
	  }
	}
      }
    }
    exch *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
    //std::cout << "! " << product1(n1, m1, n2, m2, n3, m3, n4, m4) << ", " << dir << std::endl;
  }
  //std::cout << "dir,exch = " << dir << " " << exch << std::endl;
  return dir - exch;
}


void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.pp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec1[i][2*pq];
      shell2 = Chan.ppvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.pp[i]; ++rs){
	shell3 = Chan.ppvec1[i][2*rs];
	shell4 = Chan.ppvec1[i][2*rs + 1];
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V1[i][Chan.pp[i]*pq + rs] = TBME;
	Ints.D_ME1.V1[i][Chan.pp[i]*rs + pq] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.hh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec1[i][2*pq];
      shell2 = Chan.hhvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.hh[i]; ++rs){
	shell3 = Chan.hhvec1[i][2*rs];
	shell4 = Chan.hhvec1[i][2*rs + 1];
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V2[i][Chan.hh[i]*pq + rs] = TBME;
	Ints.D_ME1.V2[i][Chan.hh[i]*rs + pq] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.pp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec1[i][2*pq];
      shell2 = Chan.ppvec1[i][2*pq + 1];
      for(int rs = 0; rs < Chan.hh[i]; ++rs){
	shell3 = Chan.hhvec1[i][2*rs];
	shell4 = Chan.hhvec1[i][2*rs + 1];
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V4[i][Chan.hh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int q = 0; q < Chan.h[i]; ++q){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell2 = Chan.hvec1[i][q];
      for(int prs = 0; prs < Chan.hpp[i]; ++prs){
	shell1 = Chan.hppvec1[i][3*prs];
	shell3 = Chan.hppvec1[i][3*prs + 1];
	shell4 = Chan.hppvec1[i][3*prs + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V5[i][Chan.hpp[i]*q + prs] = TBME;
	Ints.D_ME1.V6[i][Chan.hpp[i]*q + prs] = -1.0 * TBME;
      }
    }
    #pragma omp parallel for
    for(int s = 0; s < Chan.p[i]; ++s){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell4 = Chan.pvec1[i][s];
      for(int pqr = 0; pqr < Chan.hhp[i]; ++pqr){
	shell1 = Chan.hhpvec1[i][3*pqr];
	shell2 = Chan.hhpvec1[i][3*pqr + 1];
	shell3 = Chan.hhpvec1[i][3*pqr + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V7[i][Chan.hhp[i]*s + pqr] = TBME;
	Ints.D_ME1.V8[i][Chan.hhp[i]*s + pqr] = -1.0 * TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size2; ++i){
    #pragma omp parallel for
    for(int ps = 0; ps < Chan.hp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*ps];
      shell4 = Chan.hp2vec1[i][2*ps + 1];
      for(int qr = ps; qr < Chan.hp2[i]; ++qr){
	shell3 = Chan.hp2vec1[i][2*qr];
	shell2 = Chan.hp2vec1[i][2*qr + 1];
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V3[i][Chan.hp2[i]*ps + qr] = TBME;
	Ints.D_ME1.V3[i][Chan.hp2[i]*qr + ps] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pr = 0; pr < Chan.hp2[i]; ++pr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*pr];
      shell3 = Chan.hp2vec1[i][2*pr + 1];
      for(int qs = 0; qs < Chan.hp1[i]; ++qs){
	shell2 = Chan.hp1vec1[i][2*qs];
	shell4 = Chan.hp1vec1[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell4, shell3, L);
	Ints.D_ME1.V9[i][Chan.hp1[i]*pr + qs] = TBME;
	Ints.D_ME1.V10[i][Chan.hp1[i]*pr + qs] = -1.0 * TBME;
      }
    }
  }
}

int kron_del(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
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


int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}


/*double vint_Conversion(const Model_Space &Space, const std::vector<std::vector<double> > &ME, const int &qi, const int &qj, const int &qk, const int &ql)
{
  int kx1, ky1, kz1, kx2, ky2, kz2, Tz, Sz, ind1, ind2, ind3, MEind;
  double TBME = 0.0, CGC10, CGC11, CGC20, CGC21;
  kx1 = Space.levelsnx[qi] - Space.levelsnx[qj];
  ky1 = Space.levelsny[qi] - Space.levelsny[qj];
  kz1 = Space.levelsnz[qi] - Space.levelsnz[qj];
  kx2 = Space.levelsnx[qk] - Space.levelsnx[ql];
  ky2 = Space.levelsny[qk] - Space.levelsny[ql];
  kz2 = Space.levelsnz[qk] - Space.levelsnz[ql];
  Tz = int(0.5*Space.levelst[qi] + 0.5*Space.levelst[qj]);
  Sz = int(0.5*Space.levelsm[qi] + 0.5*Space.levelsm[qj]);
  MEind = (Sz+1)*3 + (Tz+1);
  ind1 = CART_tbInd2_rel(Space, kx1, ky1, kz1);
  ind2 = CART_tbInd2_rel(Space, kx2, ky2, kz2);
  CGC10 = CGC(0.5,0.5*Space.levelsm[qi],0.5,0.5*Space.levelsm[qj],0,Sz);
  CGC11 = CGC(0.5,0.5*Space.levelsm[qi],0.5,0.5*Space.levelsm[qj],1,Sz);
  CGC20 = CGC(0.5,0.5*Space.levelsm[qk],0.5,0.5*Space.levelsm[ql],0,Sz);
  CGC21 = CGC(0.5,0.5*Space.levelsm[qk],0.5,0.5*Space.levelsm[ql],1,Sz);
  if(Sz == 0){
    ind3 = (2*Space.CART_tb2Indsize)*ind1 + ind2;
    TBME += ME[MEind][ind3]*CGC10*CGC20; //S=00
    ind3 = (2*Space.CART_tb2Indsize)*ind1 + (Space.CART_tb2Indsize+ind2);
    TBME += ME[MEind][ind3]*CGC10*CGC21; //S=01
    ind3 = (2*Space.CART_tb2Indsize)*(Space.CART_tb2Indsize+ind1) + ind2;
    TBME += ME[MEind][ind3]*CGC11*CGC20; //S=10
    ind3 = (2*Space.CART_tb2Indsize)*(Space.CART_tb2Indsize+ind1) + (Space.CART_tb2Indsize+ind2);
    TBME += ME[MEind][ind3]*CGC11*CGC21; //S=11
  }
  else{
    ind3 = Space.CART_tb2Indsize*ind1 + ind2;
    TBME += ME[MEind][ind3]*CGC11*CGC21; //S=11
  }
  std::cout << "** " << kx1 << " " << ky1 << " " << kz1 << " " << kx2 << " " << ky2 << " " << kz2 << " " << Sz << " " << Tz << " " << TBME << std::endl;
  return TBME;
  }*/

void Perform_CC(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps)
{
  std::cout << "Performing CCD ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  int ind = 0;
  int ind1, ind2, ind3, ind4;
  double tempen, tempt, error;
  Amplitudes Amps2 = Amps;
  double CCinE, CCoutE;
  int Stot = 0;
  int Dtot = 0;

  ////////////////////////////////////////////////////////////////

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	++Dtot;
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.qnums[ind1].energy + Space.qnums[ind2].energy - Space.qnums[ind3].energy - Space.qnums[ind4].energy;
	Amps.D1.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	tempt = Ints.D_ME1.V4[chan][pind * Chan.hh[chan] + hind] / tempen;
	Amps.D1.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	//std::cout << "T: " << ind1 << ind2 << ind3 << ind4 << " = " << tempt << " , " << tempen << std::endl;
	/*if(ind1 == ind2 || ind3 == ind4){ tempt = 0.0; }
	else if(ind1 < ind2 && ind3 < ind4){ tempt = 1000*ind3 + 100*ind4 + 10*ind1 + ind2; }
	else if(ind1 > ind2 && ind3 > ind4){ tempt = 1000*ind4 + 100*ind3 + 10*ind2 + ind1; }
	else if(ind1 < ind2 && ind3 > ind4){ tempt = -1.0*(1000*ind4 + 100*ind3 + 10*ind1 + ind2); }
	else if(ind1 > ind2 && ind3 < ind4){ tempt = -1.0*(1000*ind3 + 100*ind4 + 10*ind2 + ind1); }
	Amps.D1.set_T(chan, hind * Chan.pp[chan] + pind, tempt);*/
      }
    }
  }
  Amps2.D1.Evec = Amps.D1.Evec;

  //std::cout << std::endl;
  if(Parameters.approx == "singles"){
    for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
      ++Stot;
      ind1 = Chan.hp1vec1[Chan.ind0][2*i];
      ind2 = Chan.hp1vec1[Chan.ind0][2*i + 1];
      tempen = Space.qnums[ind1].energy - Space.qnums[ind2].energy;
      tempt = 0.0;
      //tempt = 10*ind2 + ind1;
      //tempt = 0.01;
      Amps.S1.Evec[i] = tempen;
      Amps.S1.set_T(i, tempt);
      //std::cout << "t: " << ind1 << ind2 << " = " << tempt << " " << tempen << std::endl;
    }
    Amps2.S1.Evec = Amps.S1.Evec;
    Amps.S1.set_T_2(Chan, Ints);
    Amps.D1.set_T_2(Chan, Ints);
  }
  Amps2.zero(Chan, Parameters);

  //std::cout << std::endl;
  //while(ind < 1){
  while((error > 10e-16 && ind < 1000) || ind < 10){
    Doubles_Step(Space, Chan, Ints, Amps, Amps2);
    if(Parameters.approx == "singles"){
      Doubles_Step_2(Space, Chan, Ints, Amps, Amps2);
      Singles_Step(Space, Chan, Ints, Amps, Amps2);
    }

    error = 0.0;
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  ind1 = Chan.hhvec1[chan][2*hind];
	  ind2 = Chan.hhvec1[chan][2*hind + 1];
	  ind3 = Chan.ppvec1[chan][2*pind];
	  ind4 = Chan.ppvec1[chan][2*pind + 1];
	  if(ind1 == ind2 || ind3 == ind4){
	    tempt = 0.0;
	    Amps.D1.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  }
	  else{
	    tempt = Amps2.D1.get_T(chan, hind * Chan.pp[chan] + pind);
	    tempt += Ints.D_ME1.V4[chan][pind * Chan.hh[chan] + hind];
	    tempt /= Amps2.D1.Evec[chan][hind * Chan.pp[chan] + pind];
	    error += fabs((tempt - Amps.D1.T1[chan][hind * Chan.pp[chan] + pind])/Amps.D1.T1[chan][hind * Chan.pp[chan] + pind]);
	    Amps.D1.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  }
	  //std::cout << "T: " << ind1 << ind2 << ind3 << ind4 << " = " << tempt << std::endl;
	}
      }
    }
    error /= Dtot;
    //std::cout << std::endl;
    if(Parameters.approx == "singles"){
      int error2 = 0.0;
      for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
	ind1 = Chan.hp1vec1[Chan.ind0][2*i];
	ind2 = Chan.hp1vec1[Chan.ind0][2*i + 1];
	tempt = Amps2.S1.get_T(i);
	tempt /= Amps2.S1.Evec[i];
	error2 += fabs((tempt - Amps.S1.T1[i])/Amps.S1.T1[i]);
	Amps.S1.set_T(i, tempt);
	//std::cout << "t: " << ind1 << ind2 << " = " << tempt << std::endl;
      }
      error2 /= Stot;
      error += error2;
      Amps.S1.set_T_2(Chan, Ints);
      Amps.D1.set_T_2(Chan, Ints);
    }
    CCoutE = Amps.get_energy(Parameters, Chan, Ints);
    CCinE = CCoutE;
    Amps2.zero(Chan, Parameters);

    //std::cout << std::endl;
    std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCoutE << ", error = " << error << std::endl;
    ++ind;
  }
  std::cout << std::endl << std::endl;

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	std::cout << "T: " << ind1 << ind2 << ind3 << ind4 << " = " << Amps.D1.T1[chan][hind * Chan.pp[chan] + pind] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  if(Parameters.approx == "singles"){
    for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
      ind1 = Chan.hp1vec1[Chan.ind0][2*i];
      ind2 = Chan.hp1vec1[Chan.ind0][2*i + 1];
      std::cout << "t: " << ind1 << ind2 << " = " << Amps.S1.T1[i] << std::endl;
    }
  }
}


void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  if(Parameters.basis == "infinite"){
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hh = 0; hh < Chan.hh[chan]; ++hh){
	Space.qnums[Chan.hhvec1[chan][2*hh]].energy += Ints.D_ME1.V2[chan][hh*Chan.hh[chan] + hh];
      }
    }
    for(int chan = 0; chan < Chan.size2; ++chan){
      for(int hp2 = 0; hp2 < Chan.hp2[chan]; ++hp2){
	Space.qnums[Chan.hp2vec1[chan][2*hp2+1]].energy += Ints.D_ME1.V3[chan][hp2*Chan.hp2[chan] + hp2];
      }
    }
  }
}


double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const Interactions &Ints)
{
  double energy = 0.0;
  int ind, ind1;
  State tb;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].type != "hole"){ continue; }
    energy += Space.qnums[i].energy;
    for(int j = 0; j < Space.indtot; ++j){
      if(Space.qnums[j].type != "hole" || i == j){ continue; }
      plus(tb, Space.qnums[i], Space.qnums[j]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
      ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
      energy -= 0.5 * Ints.D_ME1.V2[ind1][ind];
    }
  }
  return energy;
}

void Doubles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = 0.25, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.hh[chan];
    int pp = Chan.pp[chan];
    if(hh == 0 || pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab}
    RM_dgemm(& *Amps1.D1.T1[chan].begin(), & *Ints.D_ME1.V1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &pp, &fac3, &fac1, &N, &N);
    //T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab}
    RM_dgemm(& *Ints.D_ME1.V2[chan].begin(), & *Amps1.D1.T1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac3, &fac1, &N, &N); //fac1
    //T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab}
    RM_dgemm(& *Amps1.D1.T1[chan].begin(), & *Ints.D_ME1.V4[chan].begin(), & *Amps2.D1.S1[chan].begin(), &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps2.D1.S1[chan].begin(), & *Amps1.D1.T1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac4, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.hp1[chan];
    int hp2 = Chan.hp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb}, T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia}
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    Amps2.D1.T3[chan] = Amps2.D1.T2[chan];
    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja}, T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib}
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N);
    Amps2.D1.T5[chan] = Amps2.D1.T4[chan];
    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb}
    RM_dgemm(& *Ints.D_ME1.V9[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S6[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S6[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja}
    RM_dgemm(& *Ints.D_ME1.V10[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S7[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S7[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i}
    RM_dgemm(& *Ints.D_ME1.V5[chan].begin(), & *Amps1.D1.T6[chan].begin(), & *Amps2.D1.S2[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T7[chan].begin(), & *Amps2.D1.S2[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &h, &fac6, &fac1, &N, &N);
    //T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j}
    RM_dgemm(& *Ints.D_ME1.V6[chan].begin(), & *Amps1.D1.T6[chan].begin(), & *Amps2.D1.S3[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T7[chan].begin(), & *Amps2.D1.S3[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &h, &fac6, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a}
    RM_dgemm(& *Ints.D_ME1.V7[chan].begin(), & *Amps1.D1.T8[chan].begin(), & *Amps2.D1.S4[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T9[chan].begin(), & *Amps2.D1.S4[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &p, &fac6, &fac1, &N, &N);
    //T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b}
    RM_dgemm(& *Ints.D_ME1.V8[chan].begin(), & *Amps1.D1.T8[chan].begin(), & *Amps2.D1.S5[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T9[chan].begin(), & *Amps2.D1.S5[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &p, &fac6, &fac1, &N, &N);
  }
}

void Singles_Step(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = -1.0, fac4 = -0.5;
  char N = 'N';
  int hp = Chan.hp1[Chan.ind0];
  int one = 1;
  if(hp != 0){
    //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc}
    RM_dgemm(& *Ints.D_ME1.V3[Chan.ind0].begin(), & *Amps1.S1.T1.begin(), & *Amps2.S1.T1.begin(), &hp, &one, &hp, &fac3, &fac1, &N, &N);
    //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia}
    RM_dgemm(& *Ints.D_ME1.V9[Chan.ind0].begin(), & *Amps1.D1.T2[Chan.ind0].begin(), & *Amps2.S1.S3.begin(), &hp, &hp, &hp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.S1.T1.begin(), & *Amps2.S1.S3.begin(), & *Amps2.S1.T1.begin(), &one, &hp, &hp, &fac1, &fac1, &N, &N); //fac1
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p != 0 && h != 0){
      if(hpp != 0){
	//t2(a|i){a,i} = -0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i}
	RM_dgemm(& *Ints.S_ME1.V11[chan].begin(), & *Amps1.D1.T6[chan].begin(), & *Amps2.S1.T2[chan].begin(), &p, &h, &hpp, &fac4, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(ic|kd){kcd,i}
	RM_dgemm(& *Ints.S_ME1.V11[chan].begin(), & *Amps1.S1.E6[chan].begin(), & *Amps2.S1.T2[chan].begin(), &p, &h, &hpp, &fac1, &fac1, &N, &N); //fac1
	//t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i}
	RM_dgemm(& *Ints.D_ME1.V6[chan].begin(), & *Amps1.D1.T6[chan].begin(), & *Amps2.S1.S2[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
	RM_dgemm(& *Amps1.S1.T2[chan].begin(), & *Amps2.S1.S2[chan].begin(), & *Amps2.S1.T2[chan].begin(), &p, &h, &h, &fac4, &fac1, &N, &N);
      }
      if(hhp != 0){
	//t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a}
	RM_dgemm(& *Ints.S_ME1.V12[chan].begin(), & *Amps1.D1.T8[chan].begin(), & *Amps2.S1.T3[chan].begin(), &h, &p, &hhp, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ka|lc){klc,a}
	RM_dgemm(& *Ints.S_ME1.V12[chan].begin(), & *Amps1.S1.E8[chan].begin(), & *Amps2.S1.T3[chan].begin(), &h, &p, &hhp, &fac1, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a}
	RM_dgemm(& *Ints.D_ME1.V8[chan].begin(), & *Amps1.D1.T8[chan].begin(), & *Amps2.S1.S1[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
	RM_dgemm(& *Amps1.S1.T3[chan].begin(), & *Amps2.S1.S1[chan].begin(), & *Amps2.S1.T3[chan].begin(), &h, &p, &p, &fac4, &fac1, &N, &N); //fac1
	//t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a}
	RM_dgemm(& *Ints.D_ME1.V8[chan].begin(), & *Amps1.S1.E8[chan].begin(), & *Amps2.S1.S1[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
	RM_dgemm(& *Amps1.S1.T3[chan].begin(), & *Amps2.S1.S1[chan].begin(), & *Amps2.S1.T3[chan].begin(), &h, &p, &p, &fac1, &fac1, &N, &N); //fac1
      }
    }
  }
}
 
void Doubles_Step_2(const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.hh[chan];
    int pp = Chan.pp[chan];
    //int hp = Chan.hp[chan];
    if(hh == 0 || pp == 0){ continue; }
    //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
    RM_dgemm(& *Amps1.S1.E1[chan].begin(), & *Ints.D_ME1.V1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &pp, &fac5, &fac1, &N, &N);
    //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} (2)
    RM_dgemm(& *Ints.D_ME1.V2[chan].begin(), & *Amps1.S1.E1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac5, &fac1, &N, &N);
    //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (3)
    RM_dgemm(& *Amps1.D1.T1[chan].begin(), & *Ints.D_ME1.V4[chan].begin(), & *Amps2.D1.S1[chan].begin(), &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps2.D1.S1[chan].begin(), & *Amps1.S1.E1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac6, &fac1, &N, &N);
    //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
    RM_dgemm(& *Amps1.S1.E1[chan].begin(), & *Ints.D_ME1.V4[chan].begin(), & *Amps2.D1.S1[chan].begin(), &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps2.D1.S1[chan].begin(), & *Amps1.D1.T1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac6, &fac1, &N, &N);
    //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (5)
    RM_dgemm(& *Amps2.D1.S1[chan].begin(), & *Amps1.S1.E1[chan].begin(), & *Amps2.D1.T1[chan].begin(), &hh, &pp, &hh, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.hp1[chan];
    int hp2 = Chan.hp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (6)
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N);
    //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (7)
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T3[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (8)
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (9)
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Ints.D_ME1.V3[chan].begin(), & *Amps2.D1.T5[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
    RM_dgemm(& *Ints.D_ME1.V9[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S6[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Amps2.D1.S6[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} (11)
    RM_dgemm(& *Ints.D_ME1.V9[chan].begin(), & *Amps1.S1.E2[chan].begin(), & *Amps2.D1.S6[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S6[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
    RM_dgemm(& *Ints.D_ME1.V10[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S7[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.S1.E2[chan].begin(), & *Amps2.D1.S7[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} (13)
    RM_dgemm(& *Ints.D_ME1.V10[chan].begin(), & *Amps1.S1.E2[chan].begin(), & *Amps2.D1.S7[chan].begin(), &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T2[chan].begin(), & *Amps2.D1.S7[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N);
    //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} (14)
    RM_dgemm(& *Amps1.S1.Q12[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N);
    //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} (15)
    RM_dgemm(& *Amps1.S1.Q12[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T3[chan].begin(), &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} (16)
    RM_dgemm(& *Amps1.S1.Q12[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N);
    //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} (17)
    RM_dgemm(& *Amps1.S1.Q12[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T5[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N);
    //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} (18)
    RM_dgemm(& *Amps1.S1.Q22[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T2[chan].begin(), &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N);
    //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} (19)
    RM_dgemm(& *Amps1.S1.Q22[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T3[chan].begin(), &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N);
    //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} (20)
    RM_dgemm(& *Amps1.S1.Q22[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T4[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N);
    //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} (21)
    RM_dgemm(& *Amps1.S1.Q22[chan].begin(), & *Amps1.D1.T2[chan].begin(), & *Amps2.D1.T5[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.p[chan];
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} (22)
    RM_dgemm(& *Ints.D_ME1.V5[chan].begin(), & *Amps1.S1.E6[chan].begin(), & *Amps2.D1.S2[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T7[chan].begin(), & *Amps2.D1.S2[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &h, &fac1, &fac1, &N, &N);
    //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} (23)
    RM_dgemm(& *Ints.D_ME1.V6[chan].begin(), & *Amps1.S1.E6[chan].begin(), & *Amps2.D1.S3[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T7[chan].begin(), & *Amps2.D1.S3[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &h, &fac1, &fac1, &N, &N);
    //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
    RM_dgemm(& *Amps1.D1.T6[chan].begin(), & *Amps1.S1.Q32[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &h, &fac1, &fac1, &N, &N);
    //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
    RM_dgemm(& *Amps1.D1.T6[chan].begin(), & *Amps1.S1.Q32[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &h, &fac5, &fac1, &N, &N);
    if(p != 0){
      //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
      RM_dgemm(& *Ints.S_ME1.V17[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &p, &fac5, &fac1, &N, &N);
      //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
      RM_dgemm(& *Ints.S_ME1.V17[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &p, &fac1, &fac1, &N, &N);
      //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
      RM_dgemm(& *Amps1.D1.Q12[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &p, &fac6, &fac1, &N, &N);
      //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
      RM_dgemm(& *Amps1.D1.Q12[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &p, &fac3, &fac1, &N, &N);
      //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
      RM_dgemm(& *Amps1.S1.Q52[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T6[chan].begin(), &hpp, &h, &p, &fac1, &fac1, &N, &N);
      //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
      RM_dgemm(& *Amps1.S1.Q52[chan].begin(), & *Amps1.S1.T2[chan].begin(), & *Amps2.D1.T7[chan].begin(), &hpp, &h, &p, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
    RM_dgemm(& *Ints.D_ME1.V7[chan].begin(), & *Amps1.S1.E8[chan].begin(), & *Amps2.D1.S4[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T9[chan].begin(), & *Amps2.D1.S4[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &p, &fac1, &fac1, &N, &N);
    //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
    RM_dgemm(& *Ints.D_ME1.V8[chan].begin(), & *Amps1.S1.E8[chan].begin(), & *Amps2.D1.S5[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
    RM_dgemm(& *Amps1.D1.T9[chan].begin(), & *Amps2.D1.S5[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &p, &fac1, &fac1, &N, &N);
    //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
    RM_dgemm(& *Amps1.D1.T8[chan].begin(), & *Amps1.S1.Q42[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &p, &fac1, &fac1, &N, &N);
    //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
    RM_dgemm(& *Amps1.D1.T8[chan].begin(), & *Amps1.S1.Q42[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &p, &fac5, &fac1, &N, &N);
    if(h != 0){
      //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
      RM_dgemm(& *Ints.S_ME1.V18[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &h, &fac5, &fac1, &N, &N);
      //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
      RM_dgemm(& *Ints.S_ME1.V18[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &h, &fac1, &fac1, &N, &N);
      //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
      RM_dgemm(& *Amps1.D1.Q22[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &h, &fac6, &fac1, &N, &N);
      //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
      RM_dgemm(& *Amps1.D1.Q22[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &h, &fac3, &fac1, &N, &N);
      //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
      RM_dgemm(& *Amps1.S1.Q62[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T8[chan].begin(), &hhp, &p, &h, &fac5, &fac1, &N, &N);
      //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
      RM_dgemm(& *Amps1.S1.Q62[chan].begin(), & *Amps1.S1.T3[chan].begin(), & *Amps2.D1.T9[chan].begin(), &hhp, &p, &h, &fac1, &fac1, &N, &N);
    }
  }
}

void Build_CC_Eff(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, CC_Eff &V_Eff)
{
  double fac1 = 1.0, fac5 = -1.0, fac4 = -0.5, fac3 = 0.5; //fac2 = 0.0, fac3 = 0.5, 
  int one = 1, hp1 = Chan.hp1[Chan.ind0], pp1 = Chan.pp1[Chan.ind0], hh1 = Chan.hh1[Chan.ind0];
  char N = 'N';

  if(hp1 != 0){
    // X_ia1(i|a){ia} = V9(ik|ac){ia,kc}.t1(c|k){kc}
    RM_dgemm(& *Ints.D_ME1.V9[Chan.ind0].begin(), & *Amps.S1.T1.begin(), & *V_Eff.X_ia1.begin(), &hp1, &one, &hp1, &fac1, &fac1, &N, &N);
    V_Eff.set_X_ia(Chan);
  }

  if(pp1 != 0){
    // X_ab1(a|b){ba} = f_ab.delta(a,b)
    for(int pp = 0; pp < pp1; ++pp){
      if(Chan.pp1vec1[Chan.ind0][2*pp] == Chan.pp1vec1[Chan.ind0][2*pp + 1]){
	V_Eff.X_ab1[pp] += Space.qnums[Chan.pp1vec1[Chan.ind0][2*pp]].energy;
      }
    }
    // X_ab1(a|b){ba} = V16(ka|cb){ba,kc}.t1(c|k){kc}
    RM_dgemm(& *Ints.S_ME1.V16[Chan.ind0].begin(), & *Amps.S1.T1.begin(), & *V_Eff.X_ab1.begin(), &pp1, &one, &hp1, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hhp1 = Chan.hhp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ab3(a|b){b,a} = -X_ia3(k|b){b,k}.t3(a|k){k,a}
      RM_dgemm(& *V_Eff.X_ia3[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X_ab3[chan].begin(), &p1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && hhp1 != 0){
      // X_ab3(a|b){b,a} = -(1/2).V8(kl|bc){b,klc}.T8(ac|kl){klc,a}
      RM_dgemm(& *Ints.D_ME1.V8[chan].begin(), & *Amps.D1.T8[chan].begin(), & *V_Eff.X_ab3[chan].begin(), &p1, &p1, &hhp1, &fac4, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ab(Chan);

  if(hh1 != 0){
    // X1_ij1(i|j){ji} = f_ij.delta(i,j)
    for(int hh = 0; hh < hh1; ++hh){
      if(Chan.hh1vec1[Chan.ind0][2*hh] == Chan.hh1vec1[Chan.ind0][2*hh + 1]){
	V_Eff.X1_ij1[hh] += Space.qnums[Chan.hh1vec1[Chan.ind0][2*hh]].energy;
      }
    }
    // X1_ij1(i|j){ji} = -V15(ki|jc){ji,kc}.t1(c|k){kc}
    RM_dgemm(& *Ints.S_ME1.V15[Chan.ind0].begin(), & *Amps.S1.T1.begin(), & *V_Eff.X1_ij1.begin(), &hh1, &one, &hp1, &fac5, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], hpp1 = Chan.hpp[chan];
    if(h1 != 0 && hpp1 != 0){
      // X1_ij2(i|j){i,j} = (1/2).V6(ik|cd){i,kcd}.T6(cd|jk){kcd,j}
      RM_dgemm(& *Ints.D_ME1.V6[chan].begin(), & *Amps.D1.T6[chan].begin(), & *V_Eff.X1_ij2[chan].begin(), &h1, &h1, &hpp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_ij(Chan);
  V_Eff.X1_ij2.clear();

  // X_ij1(i|j){ji} = X1_ij1(i|j){ji}
  V_Eff.X_ij1 = V_Eff.X1_ij1;
  V_Eff.X1_ij1.clear();
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan];
    if(h1 != 0 && p1 != 0){
      // X_ij2(i|j){i,j} = X_ia2(i|d){i,d}.t2(d|j){d,j}
      RM_dgemm(& *V_Eff.X_ia2[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_ij2[chan].begin(), &h1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ij(Chan);

  if(hp1 != 0){
    // X_ai1(a|i){ia} = -V3(ka|ic){ia,kc}.t1(c|k){kc}
    RM_dgemm(& *Ints.D_ME1.V3[Chan.ind0].begin(), & *Amps.S1.T1.begin(), & *V_Eff.X_ai1.begin(), &hp1, &one, &hp1, &fac5, &fac1, &N, &N);
    // X_ai1(a|i){ia} = T2(ac|ik){ia,kc}.X_ia1(k|c){kc}
    RM_dgemm(& *Amps.D1.T2[Chan.ind0].begin(), & *V_Eff.X_ia1.begin(), & *V_Eff.X_ai1.begin(), &hp1, &one, &hp1, &fac1, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hhp1 = Chan.hhp[chan], hpp1 = Chan.hpp[chan];
    if(h1 != 0 && p1 != 0){
      // X_ai2(a|i){a,i} = X_ab2(a|c){a,c}.t2(c|i){c,i}
      RM_dgemm(& *V_Eff.X_ab2[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_ai2[chan].begin(), &p1, &h1, &p1, &fac1, &fac1, &N, &N);
      if(hpp1 != 0){
	// X_ai2(a|i){a,i} = -V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i}
	RM_dgemm(& *Ints.S_ME1.V11[chan].begin(), & *Amps.D1.T6[chan].begin(), & *V_Eff.X_ai2[chan].begin(), &p1, &h1, &hpp1, &fac4, &fac1, &N, &N);
      }
      // X_ai3(a|i){i,a} = -X1_ij3(k|i){i,k}.t3(a|k){k,a}
      RM_dgemm(& *V_Eff.X1_ij3[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X_ai3[chan].begin(), &h1, &p1, &h1, &fac5, &fac1, &N, &N);
      if(hhp1 != 0){
	// X_ai3(a|i){i,a} = -V12(kl|ic){i,klc}.T8(ac|kl){klc,a}
	RM_dgemm(& *Ints.S_ME1.V12[chan].begin(), & *Amps.D1.T8[chan].begin(), & *V_Eff.X_ai3[chan].begin(), &h1, &p1, &hhp1, &fac4, &fac1, &N, &N);
      }
    }
  }
  V_Eff.set_X_ai(Chan);
  V_Eff.X1_ij3.clear();


  // X_ijab1(ij|ab){ab,ij} = V4(ij|ab){ab,ij}
  V_Eff.X_ijab1 = Ints.D_ME1.V4;


  // X(1)_iabc(ia|bc){a,ibc} = V11(ia|bc){a,ibc}
  V_Eff.X1_iabc1 = Ints.S_ME1.V11;
  V_Eff.X_iabc1 = Ints.S_ME1.V11;
  // X(1)_ijka(ij|ka){k,ija} = V12(ij|ka){k,ija}
  V_Eff.X1_ijka1 = Ints.S_ME1.V12;
  V_Eff.X_ijka1 = Ints.S_ME1.V12;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p1 = Chan.p[chan], h1 = Chan.h[chan], hpp1 = Chan.hpp[chan], hhp1 = Chan.hhp[chan];
    if(p1 != 0 && h1 != 0 && hpp1 != 0){
      // X(1)_iabc(ia|bc){a,ibc} = -((1/2)).t2.V5
      RM_dgemm(& *Amps.S1.T2[chan].begin(), & *Ints.D_ME1.V5[chan].begin(), & *V_Eff.X1_iabc1[chan].begin(), &p1, &hpp1, &h1, &fac4, &fac1, &N, &N);
      RM_dgemm(& *Amps.S1.T2[chan].begin(), & *Ints.D_ME1.V5[chan].begin(), & *V_Eff.X_iabc1[chan].begin(), &p1, &hpp1, &h1, &fac5, &fac1, &N, &N);
    }
    if(p1 != 0 && h1 != 0 && hhp1 != 0){
      // X(1)_ijka(ij|ka){k,ija} = -((1/2)).t3.V7
      RM_dgemm(& *Amps.S1.T3[chan].begin(), & *Ints.D_ME1.V7[chan].begin(), & *V_Eff.X1_ijka1[chan].begin(), &h1, &hhp1, &p1, &fac4, &fac1, &N, &N);
      RM_dgemm(& *Amps.S1.T3[chan].begin(), & *Ints.D_ME1.V7[chan].begin(), & *V_Eff.X_ijka1[chan].begin(), &h1, &hhp1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iabc(Chan);
  V_Eff.set_X1_iabc(Chan);
  V_Eff.X1_iabc1.clear();
  V_Eff.set_X_ijka(Chan);
  V_Eff.set_X1_ijka(Chan);
  V_Eff.X1_ijka1.clear();

  // X1_abcd(ab|cd){cd,ab} = V1(ab|cd){cd,ab}
  V_Eff.X1_abcd1 = Ints.D_ME1.V1;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int ppp1 = Chan.ppp[chan], p1 = Chan.p[chan], h1 = Chan.h[chan];
    if(ppp1 != 0 && p1 != 0 && h1 != 0){
      // X1_abcd(ab|cd){acd,b} = X1_iabc(ka|cd}{acd,k}.t3(b|k){k,b}
      RM_dgemm(& *V_Eff.X1_iabc2[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X1_abcd2[chan].begin(), &ppp1, &p1, &h1, &fac1, &fac1, &N, &N);
      // X1_abcd(ab|cd){bcd,a} = -X1_iabc(kb|cd}{bcd,k}.t3(a|k){k,a}
      RM_dgemm(& *V_Eff.X1_iabc2[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X1_abcd3[chan].begin(), &ppp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_abcd(Chan);
  V_Eff.X1_abcd2.clear();
  V_Eff.X1_abcd3.clear();
  V_Eff.X1_iabc2.clear();

  // X_abcd(ab|cd){cd,ab} = X1_abcd(ab|cd){cd,ab}
  V_Eff.X_abcd1 = V_Eff.X1_abcd1;
  V_Eff.X1_abcd1.clear();
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.pp[chan], hh1 = Chan.hh[chan];
    if(pp1 != 0 && hh1 != 0){
      // X_abcd(ab|cd){cd,ab} = (1/2).V4(kl|cd}{cd,kl}.T1(ab|kl){kl,ab}
      RM_dgemm(& *Ints.D_ME1.V4[chan].begin(), & *Amps.D1.T1[chan].begin(), & *V_Eff.X_abcd1[chan].begin(), &pp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }


  // X_ijkl(ij|kl){ij,kl} = V2(ij|kl){ij,kl}
  V_Eff.X_ijkl1 = Ints.D_ME1.V2;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int hhh1 = Chan.hhh[chan], h1 = Chan.h[chan], p1 = Chan.p[chan];
    if(hhh1 != 0 && h1 != 0 && p1 != 0){
      // X1_ijkl(ij|kl){ijk,l} = X1_ijkl(ij|kc}{ijk,c}.t3(c|l){c,l}
      RM_dgemm(& *V_Eff.X1_ijka2[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_ijkl2[chan].begin(), &hhh1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X1_ijkl(ij|kl){ijl,k} = -X1_ijkl(ij|lc}{ijl,c}.t3(c|k){c,k}
      RM_dgemm(& *V_Eff.X1_ijka2[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_ijkl3[chan].begin(), &hhh1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.X1_ijka2.clear();
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh1 = Chan.hh[chan], pp1 = Chan.pp[chan];
    if(hh1 != 0 && pp1 != 0){
      // X_ijkl(ij|kl){kl,ij} = (1/2).T1(kl|cd}{kl,cd}.V4(ij|cd){cd,ij}
      RM_dgemm(& *Amps.D1.T1[chan].begin(), & *Ints.D_ME1.V4[chan].begin(), & *V_Eff.X_ijkl4[chan].begin(), &hh1, &hh1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_ijkl(Chan);


  // X1_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  V_Eff.X1_iajb1 = Ints.D_ME1.V3;
  // X3_iajb(ia|jb){ib,ja} = V(ia|jb){ib,ja}
  V_Eff.X3_iajb1 = Ints.D_ME1.V3;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hhp1 = Chan.hhp2[chan], hpp1 = Chan.hpp2[chan];
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X1_iajb(ia|jb){ijb,a} = (1/2).V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      RM_dgemm(& *Ints.S_ME1.V14[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X1_iajb2[chan].begin(), &hhp1, &p1, &h1, &fac3, &fac1, &N, &N);
      // X3_iajb(ia|jb){ijb,a} = V14(ki|jb}{ijb,k}.t3(a|k){k,a}
      RM_dgemm(& *Ints.S_ME1.V14[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X3_iajb2[chan].begin(), &hhp1, &p1, &h1, &fac1, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X1_iajb(ia|jb){iab,j} = X1_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      RM_dgemm(& *V_Eff.X1_iabc3[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X1_iajb3[chan].begin(), &hpp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X3_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      RM_dgemm(& *V_Eff.X_iabc3[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X3_iajb3[chan].begin(), &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X1_iajb(Chan);
  V_Eff.X1_iajb1.clear();
  V_Eff.X1_iajb2.clear();
  V_Eff.X1_iajb3.clear();
  V_Eff.set_X3_iajb(Chan);
  V_Eff.X3_iajb2.clear();
  V_Eff.X3_iajb3.clear();
  V_Eff.X1_iabc3.clear();


  // X_iajb(ia|jb){ib,ja} = X3_iajb(ia|jb){ib,ja}
  V_Eff.X_iajb1 = V_Eff.X3_iajb1;
  V_Eff.X3_iajb1.clear();
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.hp1[chan], hp2 = Chan.hp2[chan];
    if(hp1 != 0 && hp2 != 0){
      // X_iajb(ia|jb){ib,ja} = -V10(ik|cb}{ib,kc}.T5(ca|jk){kc,ja}
      RM_dgemm(& *Ints.D_ME1.V10[chan].begin(), & *Amps.D1.T5[chan].begin(), & *V_Eff.X_iajb1[chan].begin(), &hp2, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hpp1 = Chan.hpp2[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_iajb(ia|jb){iab,j} = (1/2).X_iabc(ia|cb}{iab,c}.t2(c|j){c,j}
      RM_dgemm(& *V_Eff.X_iabc3[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_iajb3[chan].begin(), &hpp1, &h1, &p1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_iajb(Chan);


  // X_abic(ab|ic){iab,c} = V17(ab|ic){iab,c}
  V_Eff.X_abic1 = Ints.S_ME1.V17;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], ppp1 = Chan.ppp[chan], hpp1 = Chan.hpp[chan], hpp2 = Chan.hpp2[chan];
    if(h1 != 0 && p1 != 0 && hpp1 != 0){
      // X_abic(ab|ic){iab,c} = -T7(ab|ik}{iab,k}.X_ia(k|c){k,c}  
      RM_dgemm(& *Amps.D1.T7[chan].begin(), & *V_Eff.X_ia2[chan].begin(), & *V_Eff.X_abic1[chan].begin(), &hpp1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && ppp1 != 0){
      // X_abic(ab|ic){abc,i} = V*(ab|dc){abc,d}.t2(d|i){d,i}
      RM_dgemm(& *V_Eff.V_abcd[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X_abic2[chan].begin(), &ppp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hpp2 != 0){
      // X_abic(ab|ic){icb',a'} = -X1_iajb(kb|ic){icb',k'}.t3(a|k){k,a}
      RM_dgemm(& *V_Eff.X1_iajb4[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X_abic3[chan].begin(), &hpp2, &p1, &h1, &fac5, &fac1, &N, &N);
      // X_abic(ab|ic){ica',b'} = X1_iajb(ka|ic){ica',k'}.t3(b|k){k,b}
      RM_dgemm(& *V_Eff.X1_iajb4[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X_abic4[chan].begin(), &hpp2, &p1, &h1, &fac1, &fac1, &N, &N);
    }
  }
  V_Eff.V_abcd.clear();
  V_Eff.X1_iajb4.clear();
  for(int chan = 0; chan < Chan.size2; ++chan){
    int pp1 = Chan.pp1[chan], hp1 = Chan.hp1[chan], hp2 = Chan.hp2[chan];
    if(pp1 != 0 && hp1 != 0 && hp2 != 0){
      // X_abic(ab|ic){bc',ia'} = X_iabc(kb|dc){bc',kd}.T3(ad|ik){kd,ia'}
      RM_dgemm(& *V_Eff.X_iabc4[chan].begin(), & *Amps.D1.T3[chan].begin(), & *V_Eff.X_abic5[chan].begin(), &pp1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X_abic(ab|ic){ac',ib'} = -X_iabc(ka|dc){ac',kd}.T3(bd|ik){kd,ib'}
      RM_dgemm(& *V_Eff.X_iabc4[chan].begin(), & *Amps.D1.T3[chan].begin(), & *V_Eff.X_abic6[chan].begin(), &pp1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int pp1 = Chan.pp[chan], hh1 = Chan.hh[chan], hp1 = Chan.hp[chan];
    if(pp1 != 0 && hh1 != 0 && hp1 != 0){
      // X_abic(ab|ic){ic,ab} = (1/2).X_ijka(kl|ic){ic,kl}.T1(ab|kl){kl,ab}
      RM_dgemm(& *V_Eff.X_ijka5[chan].begin(), & *Amps.D1.T1[chan].begin(), & *V_Eff.X_abic7[chan].begin(), &hp1, &pp1, &hh1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X_abic(Chan);


  // X2_iajk(ia|jk){jka,i} = V18(ia|jk){jka,i}
  V_Eff.X2_iajk1 = Ints.S_ME1.V18;
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hhh1 = Chan.hhh[chan], hhp1 = Chan.hhp2[chan];
    if(h1 != 0 && p1 != 0 && hhh1 != 0){
      // X2_iajk(ia|jk){ijk,a} = -V*(il|jk){ijk,l}.t3(a|l){l,a}
      RM_dgemm(& *V_Eff.V_ijkl[chan].begin(), & *Amps.S1.T3[chan].begin(), & *V_Eff.X2_iajk2[chan].begin(), &hhh1, &p1, &h1, &fac5, &fac1, &N, &N);
    }
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X2_iajk(ia|jk){jia',k'} = X3_iajb(ia|jd){jia',d'}.t2(d|k){d,k}
      RM_dgemm(& *V_Eff.X3_iajb5[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X2_iajk3[chan].begin(), &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){kia',j'} = -X3_iajb(ia|kd){kia',d'}.t2(d|j){d,j}
      RM_dgemm(& *V_Eff.X3_iajb5[chan].begin(), & *Amps.S1.T2[chan].begin(), & *V_Eff.X2_iajk4[chan].begin(), &hhp1, &h1, &p1, &fac5, &fac1, &N, &N);
    }
  }
  V_Eff.V_ijkl.clear();
  V_Eff.X3_iajb5.clear();
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hh1 = Chan.hh1[chan], hp1 = Chan.hp1[chan], hp2 = Chan.hp2[chan];
    if(hh1 != 0 && hp1 != 0 && hp2 != 0){
      // X2_iajk(ia|jk){ij',ka'} = X_ijka(il|jc){ij',lc}.T3(ac|kl){lc,ka'}
      RM_dgemm(& *V_Eff.X_ijka4[chan].begin(), & *Amps.D1.T3[chan].begin(), & *V_Eff.X2_iajk5[chan].begin(), &hh1, &hp2, &hp1, &fac1, &fac1, &N, &N);
      // X2_iajk(ia|jk){ik',ja'} = -X_ijka(il|kc){ik',lc}.T3(ac|jk){lc,ja'}
      RM_dgemm(& *V_Eff.X_ijka4[chan].begin(), & *Amps.D1.T3[chan].begin(), & *V_Eff.X2_iajk6[chan].begin(), &hh1, &hp2, &hp1, &fac5, &fac1, &N, &N);
    }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hp1 = Chan.hp[chan], hh1 = Chan.hh[chan], pp1 = Chan.pp[chan];
    if(hh1 != 0 && pp1 != 0 && hp1 != 0){
      // X2_iajk(ia|jk){jk,ia} = (1/2).T1(cd|jk){jk,cd}.X_iabc(ia|cd){cd,ia}
      RM_dgemm(& *Amps.D1.T1[chan].begin(), & *V_Eff.X_iabc5[chan].begin(), & *V_Eff.X2_iajk7[chan].begin(), &hh1, &hp1, &pp1, &fac3, &fac1, &N, &N);
    }
  }
  V_Eff.set_X2_iajk(Chan);
  V_Eff.X2_iajk2.clear();
  V_Eff.X2_iajk3.clear();
  V_Eff.X2_iajk4.clear();
  V_Eff.X2_iajk5.clear();
  V_Eff.X2_iajk6.clear();
  V_Eff.X2_iajk7.clear();


  // X_iajk(ia|jk){jka,i} = X2_iajk(ia|jk){jka,i}
  V_Eff.X_iajk1 = V_Eff.X2_iajk1;
  V_Eff.X2_iajk1.clear();
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h1 = Chan.h[chan], p1 = Chan.p[chan], hhp1 = Chan.hhp[chan];
    if(h1 != 0 && p1 != 0 && hhp1 != 0){
      // X_iajk(ia|jk){jka,i} = T8(ca|jk){jka,c}.X_ia(i|c){c,i}
      RM_dgemm(& *Amps.D1.T8[chan].begin(), & *V_Eff.X_ia3[chan].begin(), & *V_Eff.X_iajk1[chan].begin(), &hhp1, &h1, &p1, &fac1, &fac1, &N, &N);
    }
  }

  std::cout << "X_ia" << std::endl;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hp1vec1[Chan.ind0][2*i] << "|" << Chan.hp1vec1[Chan.ind0][2*i + 1] << "> = " << V_Eff.X_ia1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ab" << std::endl;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.pp1vec1[Chan.ind0][2*i + 1] << "|" << Chan.pp1vec1[Chan.ind0][2*i] << "> = " << V_Eff.X_ab1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ij" << std::endl;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hh1vec1[Chan.ind0][2*i + 1] << "|" << Chan.hh1vec1[Chan.ind0][2*i] << "> = " << V_Eff.X_ij1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_ai" << std::endl;
  for(int i = 0; i < Chan.hp1[Chan.ind0]; ++i){
    std::cout << "<" << Chan.hp1vec1[Chan.ind0][2*i + 1] << "|" << Chan.hp1vec1[Chan.ind0][2*i] << "> = " << V_Eff.X_ai1[i] << std::endl;
  }
  std::cout << std::endl;
  std::cout << "X_iabc" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.p[i]; ++j){
      int pind1 = Chan.pvec1[i][j];
      for(int k = 0; k < Chan.hpp[i]; ++k){
	int hind1 = Chan.hppvec1[i][3*k];
	int pind2 = Chan.hppvec1[i][3*k + 1];
	int pind3 = Chan.hppvec1[i][3*k + 2];
	std::cout << "<" << hind1 << " " << pind1 << "|" << pind2 << " " << pind3 << "> = " << V_Eff.X_iabc1[i][j*Chan.hpp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_ijka" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.h[i]; ++j){
      int hind3 = Chan.hvec1[i][j];
      for(int k = 0; k < Chan.hhp[i]; ++k){
	int hind1 = Chan.hhpvec1[i][3*k];
	int hind2 = Chan.hhpvec1[i][3*k + 1];
	int pind1 = Chan.hhpvec1[i][3*k + 2];
	std::cout << "<" << hind1 << " " << hind2 << "|" << hind3 << " " << pind1 << "> = " << V_Eff.X_ijka1[i][j*Chan.hhp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_abic" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.hpp[i]; ++j){
      int hind1 = Chan.hppvec1[i][3*j];
      int pind1 = Chan.hppvec1[i][3*j + 1];
      int pind2 = Chan.hppvec1[i][3*j + 2];
      for(int k = 0; k < Chan.p[i]; ++k){
	int pind3 = Chan.pvec1[i][k];
	std::cout << "<" << pind1 << " " << pind2 << "|" << hind1 << " " << pind3 << "> = " << V_Eff.X_abic1[i][j*Chan.p[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_iajk" << std::endl;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.hhp[i]; ++j){
      int hind2 = Chan.hhpvec1[i][3*j];
      int hind3 = Chan.hhpvec1[i][3*j + 1];
      int pind1 = Chan.hhpvec1[i][3*j + 2];
      for(int k = 0; k < Chan.h[i]; ++k){
	int hind1 = Chan.hvec1[i][k];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << hind3 << "> = " << V_Eff.X_iajk1[i][j*Chan.h[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_ijab" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.pp[i]; ++j){
      int pind1 = Chan.ppvec1[i][2*j];
      int pind2 = Chan.ppvec1[i][2*j + 1];
      for(int k = 0; k < Chan.hh[i]; ++k){
	int hind1 = Chan.hhvec1[i][2*k];
	int hind2 = Chan.hhvec1[i][2*k + 1];
	std::cout << "<" << hind1 << " " << hind2 << "|" << pind1 << " " << pind2 << "> = " << V_Eff.X_ijab1[i][j*Chan.hh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_abcd" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.pp[i]; ++j){
      int pind3 = Chan.ppvec1[i][2*j];
      int pind4 = Chan.ppvec1[i][2*j + 1];
      for(int k = 0; k < Chan.pp[i]; ++k){
	int pind1 = Chan.ppvec1[i][2*k];
	int pind2 = Chan.ppvec1[i][2*k + 1];
	std::cout << "<" << pind1 << " " << pind2 << "|" << pind3 << " " << pind4 << "> = " << V_Eff.X_abcd1[i][j*Chan.pp[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_ijkl" << std::endl;
  for(int i = 0; i < Chan.size1; ++i){
    for(int j = 0; j < Chan.hh[i]; ++j){
      int hind1 = Chan.hhvec1[i][2*j];
      int hind2 = Chan.hhvec1[i][2*j + 1];
      for(int k = 0; k < Chan.hh[i]; ++k){
	int hind3 = Chan.hhvec1[i][2*k];
	int hind4 = Chan.hhvec1[i][2*k + 1];
	std::cout << "<" << hind1 << " " << hind2 << "|" << hind3 << " " << hind4 << "> = " << V_Eff.X_ijkl1[i][j*Chan.hh[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "X_iajb" << std::endl;
  for(int i = 0; i < Chan.size2; ++i){
    for(int j = 0; j < Chan.hp2[i]; ++j){
      int hind1 = Chan.hp2vec1[i][2*j];
      int pind2 = Chan.hp2vec1[i][2*j + 1];
      for(int k = 0; k < Chan.hp2[i]; ++k){
	int hind2 = Chan.hp2vec1[i][2*k];
	int pind1 = Chan.hp2vec1[i][2*k + 1];
	std::cout << "<" << hind1 << " " << pind1 << "|" << hind2 << " " << pind2 << "> = " << V_Eff.X_iajb1[i][j*Chan.hp2[i] + k] << std::endl;
      }
    }
  }
  std::cout << std::endl;
}



void EE_EOM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const CC_Eff &V_Eff)
{
  std::cout << "EE_EOM!!" << std::endl;
  int Nbit = std::ceil(Space.indtot/64.0);
  std::vector<std::vector<std::vector<int> > > chanvec1(Chan.size2);
  std::vector<std::vector<std::vector<int> > > chanvec2(Chan.size2);
  std::vector<int> vec2(2);

  /*std::vector<unsigned long long> gstate(Nbit, 0);
  unsigned long long one = 1;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.qnums[i].type == "hole"){ gstate[std::floor(i/64.0)] += (one << (i%64)); }
    }*/

  State state;
  int ind;
  for(int chan1 = 0; chan1 < Chan.size3; ++chan1){ // p
    for(int chan2 = 0; chan2 < Chan.size3; ++chan2){ // h
      if(Chan.p[chan1] == 0 || Chan.h[chan2] == 0){ continue; }
      minus(state, Chan.qnums3[chan1], Chan.qnums3[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      vec2[0] = chan1;
      vec2[1] = chan2;
      chanvec1[ind].push_back(vec2);
    }
  }
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){ // pp
    for(int chan2 = 0; chan2 < Chan.size1; ++chan2){ // hh
      if(Chan.pp[chan1] == 0 || Chan.hh[chan2] == 0){ continue; }
      minus(state, Chan.qnums1[chan1], Chan.qnums1[chan2]);
      ind = ChanInd_2b_cross(Parameters.basis, Space, state);
      if(ind < 0 || ind >= Chan.size2){ continue; }
      vec2[0] = chan1;
      vec2[1] = chan2;
      chanvec2[ind].push_back(vec2);
    }
  }

  for(int chan = 0; chan < Chan.size2; ++chan){
    int count1 = 0;
    int count2 = 0;
    int chansize1 = int(chanvec1[chan].size());
    int chansize2 = int(chanvec2[chan].size());
    std::vector<std::vector<int> > hpvec2;
    std::vector<std::vector<int> > hhppvec2;
    std::vector<int> vec(2);

    for(int i = 0; i < chansize1; ++i){
      for(int j = 0; j < Chan.h[chanvec1[chan][i][1]]; ++j){
	for(int k = 0; k < Chan.p[chanvec1[chan][i][0]]; ++k){
	  ++count1;
	  vec[0] = Chan.hvec1[chanvec1[chan][i][1]][j];
	  vec[1] = Chan.pvec1[chanvec1[chan][i][0]][k];
	  hpvec2.push_back(vec);
	}
      }
    }
    vec.resize(4);
    for(int i = 0; i < chansize2; ++i){
      for(int j = 0; j < Chan.hh[chanvec2[chan][i][1]]; ++j){
	if(Chan.hhvec1[chanvec2[chan][i][1]][2*j] >= Chan.hhvec1[chanvec2[chan][i][1]][2*j+1]){ continue; }
	for(int k = 0; k < Chan.pp[chanvec2[chan][i][0]]; ++k){
	  if(Chan.ppvec1[chanvec2[chan][i][0]][2*k] >= Chan.ppvec1[chanvec2[chan][i][0]][2*k+1]){ continue; }
	  ++count2;
	  vec[0] = Chan.hhvec1[chanvec2[chan][i][1]][2*j];
	  vec[1] = Chan.hhvec1[chanvec2[chan][i][1]][2*j+1];
	  vec[2] = Chan.ppvec1[chanvec2[chan][i][0]][2*k];
	  vec[3] = Chan.ppvec1[chanvec2[chan][i][0]][2*k+1];
	  hhppvec2.push_back(vec);
	}
      }
    }

    std::cout << "Chan = " << chan << ": " << Chan.qnums2[chan].ml << " " << Chan.qnums2[chan].m << " " << Chan.qnums2[chan].t << std::endl;
    for(int i = 0; i < count1; ++i){
      std::cout << hpvec2[i][0] << hpvec2[i][1] << " ";
    }
    for(int i = 0; i < count2; ++i){
      std::cout << hhppvec2[i][0] << hhppvec2[i][1] << hhppvec2[i][2] << hhppvec2[i][3] << " ";
    }
    std::cout << std::endl;

    int N = count1 + count2;
    if(chan == Chan.ind0){ N += 1; }
    std::vector<double> Ham(N*N, 0.0);
    std::vector<unsigned long long> ket, bra;
    double ME, ME1, ME0;
    std::vector<int> p, q;

    for(int col = 0; col < N; ++col){
      ket.assign(Nbit, 0);
      if(col < count1){ bitsetup(hpvec2[col], ket); }
      else if(col < count1 + count2){ bitsetup(hhppvec2[(col - count1)], ket); }
      for(int row = 0; row < N; ++row){
	bra.assign(Nbit, 0);
	if(row < count1){ bitsetup(hpvec2[row], bra); }
	else if(row < count1 + count2){ bitsetup(hhppvec2[(row - count1)], bra); }
	ME = 0.0;
	std::cout << std::bitset<4>(bra[0]) << " : " << std::bitset<4>(ket[0]) << std::endl << std::endl;

	p.resize(0);
	q.resize(2);
	for(int hp1 = 0; hp1 < Chan.hp1[Chan.ind0]; ++hp1){ // {i+a} -> ia
	  ME0 = V_Eff.X_ia1[hp1];
	  if(ME0 == 0.0){ continue; }
	  q[0] = Chan.hp1vec1[Chan.ind0][2*hp1];
	  q[1] = Chan.hp1vec1[Chan.ind0][2*hp1 + 1];
	  ME1 = matrixe(bra, ket, p, q, ME0);
	  if(ME1 != 0){ std::cout << "ia = " << q[0] << q[1] << " = " << ME1 << std::endl; }
	  ME += ME1;
	}

	p.resize(1);
	q.resize(1);
	for(int pp1 = 0; pp1 < Chan.pp1[Chan.ind0]; ++pp1){ // {a+b}
	  ME0 = V_Eff.X_ab1[pp1];
	  if(ME0 == 0.0){ continue; }
	  p[0] = Chan.pp1vec1[Chan.ind0][2*pp1 + 1];
	  q[0] = Chan.pp1vec1[Chan.ind0][2*pp1];
	  ME1 = matrixe(bra, ket, p, q, ME0);
	  if(ME1 != 0){ std::cout << "ab = " << p[0] << q[0] << " = " << ME1 << std::endl; }
	  ME += ME1;
	}

	for(int hh1 = 0; hh1 < Chan.hh1[Chan.ind0]; ++hh1){ // {i+j} -> ij+ = -1 * j+i
	  ME0 = V_Eff.X_ij1[hh1];
	  if(ME0 == 0.0){ continue; }
	  p[0] = Chan.hh1vec1[Chan.ind0][2*hh1];
	  q[0] = Chan.hh1vec1[Chan.ind0][2*hh1 + 1];
	  ME1 = matrixe(bra, ket, p, q, ME0);
	  if(ME1 != 0){ std::cout << "ij = " << q[0] << p[0] << " = " << -ME1 << std::endl; }
	  ME -= ME1;
	}

	p.resize(2);
	q.resize(0);
	for(int hp1 = 0; hp1 < Chan.hp1[Chan.ind0]; ++hp1){ // {a+i} -> a+i+
	  ME0 = V_Eff.X_ai1[hp1];
	  if(ME0 == 0.0){ continue; }
	  p[0] = Chan.hp1vec1[Chan.ind0][2*hp1 + 1];
	  p[1] = Chan.hp1vec1[Chan.ind0][2*hp1];
	  ME1 = matrixe(bra, ket, p, q, ME0);
	  if(ME1 != 0){ std::cout << "ai = " << p[0] << p[1] << " = " << ME1 << std::endl; }
	  ME += ME1;
	}

	p.resize(1);
	q.resize(3);
	for(int i = 0; i < Chan.size3; ++i){ // {i+a+cb} -> ia+cb = a+ibc
	  for(int p1 = 0; p1 < Chan.p[i]; ++p1){
	    for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
	      ME0 = V_Eff.X_iabc1[i][p1*Chan.hpp[i] + hpp1];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.pvec1[i][p1];
	      q[0] = Chan.hppvec1[i][3*hpp1];
	      q[1] = Chan.hppvec1[i][3*hpp1 + 1];
	      q[2] = Chan.hppvec1[i][3*hpp1 + 2];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "iabc = " << q[0] << p[0] << q[1] << q[2] << " = " << ME1 << std::endl; }
	      ME += 0.5 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size3; ++i){ // {i+j+ak} -> ijak+ = -1 * k+ija
	  for(int h1 = 0; h1 < Chan.h[i]; ++h1){
	    for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
	      ME0 = V_Eff.X_ijka1[i][h1*Chan.hhp[i] + hhp1];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.hvec1[i][h1];
	      q[0] = Chan.hhpvec1[i][3*hhp1];
	      q[1] = Chan.hhpvec1[i][3*hhp1 + 1];
	      q[2] = Chan.hhpvec1[i][3*hhp1 + 2];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "ijka = " << q[0] << q[1] << p[0] << q[2] << " = " << -ME1 << std::endl; }
	      ME -= 0.5 * ME1;
	    }
	  }
	}

	p.resize(3);
	q.resize(1);
	for(int i = 0; i < Chan.size3; ++i){ // {a+b+ci} -> a+b+ci+ = -1 * i+a+b+c
	  for(int hpp1 = 0; hpp1 < Chan.hpp[i]; ++hpp1){
	    for(int p1 = 0; p1 < Chan.p[i]; ++p1){
	      ME0 = V_Eff.X_abic1[i][hpp1*Chan.p[i] + p1];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.hppvec1[i][3*hpp1];
	      p[1] = Chan.hppvec1[i][3*hpp1 + 1];
	      p[2] = Chan.hppvec1[i][3*hpp1 + 2];
	      q[0] = Chan.pvec1[i][p1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "abic = " << p[1] << p[2] << p[0] << q[0] << " = " << -ME1 << std::endl; }
	      ME -= 0.5 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size3; ++i){ // {i+a+kj} -> ia+k+j+ = j+k+a+i
	  for(int hhp1 = 0; hhp1 < Chan.hhp[i]; ++hhp1){
	    for(int h1 = 0; h1 < Chan.h[i]; ++h1){
	      ME0 = V_Eff.X_iajk1[i][hhp1*Chan.h[i] + h1];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.hhpvec1[i][3*hhp1];
	      p[1] = Chan.hhpvec1[i][3*hhp1 + 1];
	      p[2] = Chan.hhpvec1[i][3*hhp1 + 2];
	      q[0] = Chan.hvec1[i][h1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "iajk = " << q[0] << p[2] << p[0] << p[1] << " = " <<  ME1 << std::endl; }
	      ME += 0.5 * ME1;
	    }
	  }
	}
	
	p.resize(2);
	q.resize(2);
	for(int i = 0; i < Chan.size1; ++i){ // {a+b+dc}
	  for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
	    for(int pp2 = 0; pp2 < Chan.pp[i]; ++pp2){
	      ME0 = V_Eff.X_abcd1[i][pp1*Chan.pp[i] + pp2];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.ppvec1[i][2*pp2];
	      p[1] = Chan.ppvec1[i][2*pp2 + 1];
	      q[0] = Chan.ppvec1[i][2*pp1 + 1];
	      q[1] = Chan.ppvec1[i][2*pp1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "abcd = " << p[0] << p[1] << q[1] << q[0] << " = " <<  ME1 << std::endl; }
	      ME += 0.25 * ME1;
	    }
	  }
	}

	for(int i = 0; i < Chan.size1; ++i){ // {i+j+lk} -> ijl+k+ = l+k+ij
	  for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
	    for(int hh2 = 0; hh2 < Chan.hh[i]; ++hh2){
	      ME0 = V_Eff.X_ijkl1[i][hh1*Chan.hh[i] + hh2];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.hhvec1[i][2*hh2 + 1];
	      p[1] = Chan.hhvec1[i][2*hh2];
	      q[0] = Chan.hhvec1[i][2*hh1];
	      q[1] = Chan.hhvec1[i][2*hh1 + 1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "ijkl = " << q[0] << q[1] << p[1] << p[0] << " = " <<  ME1 << std::endl; }
	      ME += 0.25 * ME1;
	    }
	  }
	}
	
	for(int i = 0; i < Chan.size2; ++i){ // {i+a+bj} -> ia+bj+ = j+a+ib
	  for(int hp1 = 0; hp1 < Chan.hp2[i]; ++hp1){
	    for(int hp2 = 0; hp2 < Chan.hp2[i]; ++hp2){
	      ME0 = V_Eff.X_iajb1[i][hp1*Chan.hp2[i] + hp2];
	      if(ME0 == 0.0){ continue; }
	      p[0] = Chan.hp2vec1[i][2*hp2];
	      p[1] = Chan.hp2vec1[i][2*hp2 + 1];
	      q[0] = Chan.hp2vec1[i][2*hp1];
	      q[1] = Chan.hp2vec1[i][2*hp1 + 1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "iajb = " << q[0] << p[1] << p[0] << q[1] << " = " <<  ME1 << std::endl; }
	      ME += ME1;
	    }
	  }
	}

	p.resize(0);
	q.resize(4);
	for(int i = 0; i < Chan.size1; ++i){ // {i+j+ba} -> ijba
	  for(int pp1 = 0; pp1 < Chan.pp[i]; ++pp1){
	    for(int hh1 = 0; hh1 < Chan.hh[i]; ++hh1){
	      ME0 = V_Eff.X_ijab1[i][pp1*Chan.hh[i] + hh1];
	      if(ME0 == 0.0){ continue; }
	      q[0] = Chan.hhvec1[i][2*hh1];
	      q[1] = Chan.hhvec1[i][2*hh1 + 1];
	      q[2] = Chan.ppvec1[i][2*pp1 + 1];
	      q[3] = Chan.ppvec1[i][2*pp1];
	      ME1 = matrixe(bra, ket, p, q, ME0);
	      if(ME1 != 0){ std::cout << "ijab = " << q[0] << q[1] << q[2] << q[3] << " = " <<  ME1 << std::endl; }
	      ME += 0.25 * ME1;
	    }
	  }
	}

	std::cout << "ME = " << ME << std::endl << std::endl;;
	// fill column major
	Ham[N*col + row] = ME;
      }
    }
    
    for(int i = 0; i < N; ++i){
      for(int j = 0; j < N; ++j){
	std::cout << Ham[N*j + i] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    char job = 'V';
    std::vector<double> Vl(N * N);
    std::vector<double> Vr(N * N);
    std::vector<double> wr(N);
    std::vector<double> wi(N);
    std::vector<double> work(5*N);
    int lwork = 5*N, info;
    dgeev_(&job, &job, &N, & *Ham.begin(), &N, & *wr.begin(), & *wi.begin(), & *Vl.begin(), &N, & *Vr.begin(), &N, & *work.begin(), &lwork, &info);

    for(int i = 0; i < N; ++i){
      std::cout << i+1 << " :  " << wr[i] << ", " << wi[i] << std::endl;
      for(int j = 0; j < N; ++j){ std::cout << Vl[N*i + j] << " "; }
      std::cout << std::endl;
      for(int j = 0; j < N; ++j){ std::cout << Vr[N*i + j] << " "; }
      std::cout << std::endl;
    }

  }
}

void bitsetup(const std::vector<int> &vec, std::vector<unsigned long long> &state)
{
  unsigned long long one = 1;
  for(int i = 0; i < int(vec.size()); ++i){ state[std::floor(vec[i]/64.0)] += (one << (vec[i]%64)); }
}

/*void bitsetup(const std::vector<int> &ivec, const std::vector<int> &avec, std::vector<unsigned long long> &gstate)
{
  unsigned long long one = 1;
  for(int i = 0; i < int(ivec.size()); ++i){ gstate[std::floor(ivec[i]/64.0)] -= (one << (ivec[i]%64)); }
  for(int i = 0; i < int(avec.size()); ++i){ gstate[std::floor(avec[i]/64.0)] += (one << (avec[i]%64)); }
  }*/
 
double matrixe(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const std::vector<int> &p, const std::vector<int> &q, const double &ME)
{
  unsigned long long one = 1, zero = 0;
  unsigned long long comp;
  std::vector<unsigned long long> pbit(p.size());
  std::vector<unsigned long long> qbit(q.size());
  std::vector<unsigned long long> tempket = ket;
  std::vector<int> p1(p.size());
  std::vector<int> q1(q.size());
  double phase = 1.0;
  int bcount;
  int flag1 = 1, flag2 = 1;

  for(int i = 0; i < int(p.size()); ++i){
    pbit[i] = one << (p[i]%64);
    p1[i] = std::floor(p[i]/64);
    if((pbit[i] & bra[p1[i]]) == 0){ flag1 = 0; }
  }
  for(int i = 0; i < int(q.size()); ++i){
    qbit[i] = one << (q[i]%64);
    q1[i] = std::floor(q[i]/64);
    if((qbit[i] & tempket[q1[i]]) == 0){ flag2 = 0; }
  }

  if(flag1 == 1 && flag2 == 1){
    //std::cout << "!! " << ME << std::endl;
    for(int i = int(q.size())-1; i >= 0; --i){
      comp = tempket[q1[i]] & ~(~zero << (q[i]%64));
      for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
      for(int j = q1[i]-1; j >= 0; --j){ std::bitset<64> bket (tempket[j]); bcount += bket.count(); }
      phase = phase*pow(-1.0, bcount); tempket[q1[i]] = tempket[q1[i]]^qbit[i];
      //std::cout << "?? " << phase << " " << std::bitset<4>(tempket[0]) << std::endl;
    }
    for(int i = int(p.size())-1; i >= 0; --i){
      comp = tempket[p1[i]] & ~(~zero << (p[i]%64));
      for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
      for(int j = p1[i]-1; j >= 0; --j){ std::bitset<64> bket (tempket[j]); bcount += bket.count(); }
      phase = phase*pow(-1.0, bcount); tempket[p1[i]] = tempket[p1[i]]^pbit[i];
      //std::cout << "?? " << phase << " " << std::bitset<4>(tempket[0]) << std::endl;
    }
    //std::cout << "$$ " << std::bitset<4>(bra[0]) << " ?= " << std::bitset<4>(tempket[0]) << std::endl << std::endl;
    for(int i = 0; i < int(bra.size()); ++i){ if(bra[i] != tempket[i]){ return 0.0; } }
    return phase * ME;
  }
  else{ return 0.0; }
}


/*std::vector<std::vector<std::vector<double> > > p(1);
  std::vector<std::vector<std::vector<double> > > delp(1);
  std::vector<double> B;
  std::vector<double> B2;
  std::vector<int> ipiv;
  std::vector<double> work;
  int lwork;
  int info = 0;
  int maxl = 10;
  int N = 1;

  p[0] = CCin.T1;
  delp[0] = CCin.T1;
  CCout.CCDE = 0.0;
  error = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	error += 1.0e10;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	p[0][chan][hind * Chan.pp[chan] + pind] = tempt;
	delp[0][chan][hind * Chan.pp[chan] + pind] = 1.0e10;
      }
    }
  }
  CCout = CCin;
  error = std::sqrt(error)/(Parameters.P + Parameters.N);

  while((error > 10e-10 && ind < 1000) || ind < 10){
    if(N != maxl){
      p.push_back(CCin.T1);
      delp.push_back(CCin.T1);
      ++N;
      B.assign((N + 1) * (N + 1), 0.0);
      for(int i = N - 1; i >= 0; --i){
	for(int j = N - 1; j >= 0; --j){
	  B[(N + 1) * i + j] = B[N * i + j];
	}
      }
    }
    else{
      p[0] = CCin.T1;
      delp[0] = CCin.T1;
      p.resize(1);
      delp.resize(1);
      N = 1;
    }
    Doubles_Step(Space, Chan, CCME, CCin, CCout);
    CCout.CCDE = 0.0;
    error = 0.0;
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	  tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	  CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  error += (tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind]) * (tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind]);
	  p[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt;
	  delp[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt - p[N - 2][chan][hind * Chan.pp[chan] + pind];
	  for(int i = 0; i < N; ++i){
	    B[(N + 1) * i + N - 1] += delp[i][chan][hind * Chan.pp[chan] + pind] * delp[N - 1][chan][hind * Chan.pp[chan] + pind];
	    if(i == N - 1){ break; } // don't double count
	    B[(N + 1) * (N - 1) + i] += delp[N - 1][chan][hind * Chan.pp[chan] + pind] * delp[i][chan][hind * Chan.pp[chan] + pind];
	  }
	}
      }
    }
    for(int i = 0; i < N; ++i){ B[(N + 1) * i + N] = -1.0; B[(N + 1) * N + i] = -1.0; }
    B[(N + 1) * N + N] = 0.0;

    int P = N + 1;
    lwork = sizeof(double) * P;
    ipiv.resize(P);
    work.resize(sizeof(double) * P);
    B2 = B;
    dgetrf_(&P, &P, & *B2.begin(), &P, & *ipiv.begin(), &info);
    dgetri_(&P, & *B2.begin(), &P, & *ipiv.begin(), & *work.begin(), &lwork, &info);
    error = std::sqrt(error)/(Parameters.P + Parameters.N);
    
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  double tempt = 0.0;
	  for(int i = 0; i < N; ++i){ tempt += B2[(N + 1) * i + N] * p[i][chan][hind * Chan.pp[chan] + pind]; }
	  CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	}
      }
    }
    std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
    ++ind;
  }
  std::cout << std::endl << std::endl;
  
  return CCout;*/

  ////////////////////////////////////////

  /*std::vector<std::vector<std::vector<double> > > Vin(3); // 0, V/E
  std::vector<std::vector<std::vector<double> > > F(2); // Vout - Vin
  std::vector<std::vector<std::vector<double> > > delV(1);
  std::vector<std::vector<std::vector<double> > > delF(1);
  std::vector<std::vector<std::vector<double> > > u(1);
  std::vector<double> a(1, 0.0);
  std::vector<double> B(1, 0.0);
  std::vector<double> w(2, 1.0);
  std::vector<int> ipiv;
  std::vector<double> work;
  double tempt, tempen, norm;
  int lwork;
  int info = 0;
  int maxl = 10;
  int P = 1;
  int N = 2;

  double alpha = 0.9;
  double w0 = 0.01;

  Vin[0] = CCin.T1;
  Vin[1] = CCin.T1;
  F[0] = CCin.T1;
  F[1] = CCin.T1;
  delV[0] = CCin.T1;
  delF[0] = CCin.T1;
  u[0] = CCin.T1;
  norm = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    int ind1, ind2, ind3, ind4;
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = alpha * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	Vin[1][chan][hind * Chan.pp[chan] + pind] = tempt;
	F[0][chan][hind * Chan.pp[chan] + pind] = tempt;
	norm += tempt*tempt;
      }
    }
  }
  w[0] = 0.0;
  norm = std::sqrt(norm);
  CCout.Evec = CCin.Evec;

  Doubles_Step(Space, Chan, CCME, CCin, CCout);

  norm = 0.0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	F[1][chan][hind * Chan.pp[chan] + pind] = tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind];
	norm += (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind]) *
	  (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind]);
      }
    }
  }
  w[1] = 1.0;
  norm = std::sqrt(norm);

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	delV[0][chan][hind * Chan.pp[chan] + pind] = (Vin[1][chan][hind * Chan.pp[chan] + pind] - Vin[0][chan][hind * Chan.pp[chan] + pind])/norm;
	delF[0][chan][hind * Chan.pp[chan] + pind] = (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind])/norm;
	u[0][chan][hind * Chan.pp[chan] + pind] = alpha * delF[0][chan][hind * Chan.pp[chan] + pind] + delV[0][chan][hind * Chan.pp[chan] + pind];
	a[0] += w[0] * delF[0][chan][hind * Chan.pp[chan] + pind] * w[0] * delF[0][chan][hind * Chan.pp[chan] + pind];
      }
    }
  }
  a[0] += w0*w0;

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = Vin[1][chan][hind * Chan.pp[chan] + pind] + alpha * F[1][chan][hind * Chan.pp[chan] + pind];
	tempt -= w[0] * delF[0][chan][hind * Chan.pp[chan] + pind] * F[1][chan][hind * Chan.pp[chan] + pind] * u[0][chan][hind * Chan.pp[chan] + pind] / a[0];
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
      }
    }
  }
  Vin[2] = CCin.T1;

  while((error > 10e-8 && ind < 1000) || ind < 10)
    {
      Doubles_Step(Space, Chan, CCME, CCin, CCout);
      P = int(delF.size());
      N = int(delF.size() + 1);
      if(P != maxl){
	F.resize(F.size() + 1);
	Vin.resize(Vin.size() + 1);
	delV.resize(delV.size() + 1);
	delF.resize(delF.size() + 1);
	u.resize(u.size() + 1);
	w.resize(w.size() + 1);
	Vin[Vin.size() - 1] = CCout.T1;
	F[F.size() - 1] = CCout.T1;
	delV[delV.size() - 1] = CCout.T1;
	delF[delF.size() - 1] = CCout.T1;
	u[u.size() - 1] = CCout.T1;
	w[w.size() - 1] = 1;
	P = int(delF.size());
	N = int(delF.size() + 1);
	a.resize(P * P);
	for(int i = P - 2; i >= 0; --i){
	  for(int j = P - 2; j >= 0; --j){
	    a[P * i + j] = a[(P - 1) * i + j];
	  }
	}
      }
      else{
	for(int i = 0; i < (P - 1); ++i){ delV[i] = delV[i+1]; delF[i] = delF[i+1]; u[i] = u[i+1]; }
	for(int i = 0; i < (N - 1); ++i){ F[i] = F[i+1]; w[i] = w[i+1]; }
	for(int i = 0; i < N; ++i){ Vin[i] = Vin[i+1]; }
	for(int i = 0; i < (P - 1); ++i){ 
	  for(int j = 0; j < (P - 1); ++j){
	    a[P * i + j] = a[P * (i + 1) + (j + 1)];
	  }
	}
      }

      for(int i = 0; i < P; ++i){
	a[P * i + (P - 1)] = 0.0;
	a[P * (P - 1) + i] = 0.0;
      }

      CCout.CCDE = 0.0;
      error = 0.0;
      norm = 0.0;
      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	    tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	    CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    F[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind];
	    error += (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]) * (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]);
	    norm += (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind]) *
	      (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind]);
	  }
	}
      }
      w[N - 1] = 1.0/error;
      norm = std::sqrt(norm);
      error = std::sqrt(error);

      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    delV[P - 1][chan][hind * Chan.pp[chan] + pind] = (Vin[N - 1][chan][hind * Chan.pp[chan] + pind] - Vin[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    delF[P - 1][chan][hind * Chan.pp[chan] + pind] = (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    u[P - 1][chan][hind * Chan.pp[chan] + pind] = alpha * delF[P - 1][chan][hind * Chan.pp[chan] + pind] + delV[P - 1][chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      a[P * i + (P - 1)] += w[i] * delF[i][chan][hind * Chan.pp[chan] + pind] * w[P - 1] * delF[P - 1][chan][hind * Chan.pp[chan] + pind];
	      if(i == P - 1){ break; } // don't double count
	      a[P * (P - 1) + i] += w[P - 1] * delF[P - 1][chan][hind * Chan.pp[chan] + pind] * w[i] * delF[i][chan][hind * Chan.pp[chan] + pind];
	    }
	  }
	}
      }
      a[P * (P - 1) + (P - 1)] += w0 * w0;

      lwork = P;
      ipiv.resize(P);
      work.resize(sizeof(double) * P);
      B = a;
      dgetrf_(&P, &P, & *B.begin(), &P, & *ipiv.begin(), &info);
      dgetri_(&P, & *B.begin(), &P, & *ipiv.begin(), & *work.begin(), &lwork, &info);

      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = Vin[N - 1][chan][hind * Chan.pp[chan] + pind] + alpha * F[N - 1][chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      for(int j = 0; j < P; ++j){
		tempt += w[i] * delF[i][chan][hind * Chan.pp[chan] + pind] * F[N - 1][chan][hind * Chan.pp[chan] + pind] * 
		  B[P * i + j] * u[j][chan][hind * Chan.pp[chan] + pind];
	      }
	    }
	    CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  }
	}
      }
      Vin[N] = CCin.T1;
      
      std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
      ++ind;
    }
  std::cout << std::endl << std::endl;
    
return CCout;*/
