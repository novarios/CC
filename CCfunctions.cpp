#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

//Function to Search for index
int Index(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[2*i] == p && vec1[2*i + 1] == q){ ind1 = i; break; }
    else if(i == num1 - 1){
      for(int j = 0; j < num1; ++j){ std::cout << vec1[2*j] << "," << vec1[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index for " << p << " " << q << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[2*i] == r && vec2[2*i + 1] == s){ ind2 = i; break; }
    else if(i == num2 - 1){
      for(int j = 0; j < num2; ++j){ std::cout << vec2[2*j] << "," << vec2[2*j + 1] << " "; }
      std::cout << std::endl;
      std::cerr << "Index for " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

int Index2(const std::vector<int> &vec1, const std::vector<int> &vec2, const int &num1, const int &num2, const int &p, const int &q, const int &r, const int &s)
{
  int ind1, ind2;
  for(int i = 0; i < num1; ++i){
    if(vec1[i] == p){ ind1 = i; break; }
    else if(i == num1 - 1){ 
      for(int j = 0; j < num1; ++j){ std::cout << vec1[j] << " "; }
      std::cout << std::endl;
      std::cerr << "Index2 for " << p << ", not Found" << std::endl; exit(1);
    }
  }
  for(int i = 0; i < num2; ++i){
    if(vec2[3*i] == q && vec2[3*i + 1] == r && vec2[3*i + 2] == s){ ind2 = i; break; }
    else if(i == num2 - 1){ 
      for(int j = 0; j < num2; ++j){ std::cout << vec2[3*j] << "," << vec2[3*j + 1] << "," << vec2[3*j + 2] << " "; }
      std::cout << std::endl;
      std::cerr << "Index2 for " << q << " " << r << " " << s << ", not Found" << std::endl; exit(1);
    }
  }
  return ind1*num2 + ind2;  
}

CCD::CCD(const Channels &Chan, const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << "Building Amplitude Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  CCDE = 0;

  Tmap.resize(Chan.size1);
  Evec.resize(Chan.size1);
  T1.resize(Chan.size1);

  T2.resize(Chan.size2);
  T2T.resize(Chan.size2);
  T3.resize(Chan.size2);

  T4.resize(Chan.size3);
  T5.resize(Chan.size3);
  T6.resize(Chan.size3);
  T7.resize(Chan.size3);

  S1.resize(Chan.size1);
  S2.resize(Chan.size3);
  S3.resize(Chan.size3);
  S4.resize(Chan.size2);
  S4T.resize(Chan.size2);

  CCDE = 0.0;

  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    Tmap[i].resize(Chan.hh[i] * Chan.pp[i] * 15);
    Evec[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    T1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    S1[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    for(int j = 0; j < (Chan.hh[i] * Chan.pp[i]); ++j){
      Tmap[i][15*j] = j;
    }
  }

  if(Parameters.basis == "HO"){
    #pragma omp parallel for
    for(int i = 0; i < Chan.size3; ++i){
      T4[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
      T5[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
      S2[i].assign(Chan.p[i] * Chan.p[i], 0.0);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int hhp = 0; hhp < Chan.hhp[i]; ++hhp){
	hind1 = Chan.hhpvec1[i][3*hhp];
	hind2 = Chan.hhpvec1[i][3*hhp + 1];
	pind1 = Chan.hhpvec1[i][3*hhp + 2];
	ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[hind1]+Space.levelsl[hind2]), Space.levelsm[hind1]+Space.levelsm[hind2], 
			 Space.levelst[hind1]+Space.levelst[hind2]);
	for(int p2 = 0; p2 < Chan.p[i]; ++p2){
	  ++count;
	  pind2 = Chan.pvec1[i][p2];
	  if(pind1 == pind2){ continue; }
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 7] = i;
	  Tmap[ind1][15*ind + 8] = count;
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind2, pind1);
	  Tmap[ind1][15*ind + 9] = i;
	  Tmap[ind1][15*ind + 10] = count;
	}
      }
    }
      
    #pragma omp parallel for
    for(int i = 0; i < Chan.size3; ++i){
      T6[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
      T7[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
      S3[i].assign(Chan.h[i] * Chan.h[i], 0.0);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int hpp = 0; hpp < Chan.hpp[i]; ++hpp){
	hind1 = Chan.hppvec1[i][3*hpp];
	pind1 = Chan.hppvec1[i][3*hpp + 1];
	pind2 = Chan.hppvec1[i][3*hpp + 2];
	ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[pind1]+Space.levelsl[pind2]), Space.levelsm[pind1]+Space.levelsm[pind2], 
			 Space.levelst[pind1]+Space.levelst[pind2]);
	for(int h2 = 0; h2 < Chan.h[i]; ++h2){
	  ++count;
	  hind2 = Chan.hvec1[i][h2];
	  if(hind1 == hind2){ continue; }
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 11] = i;
	  Tmap[ind1][15*ind + 12] = count;
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind2, hind1, pind1, pind2);
	  Tmap[ind1][15*ind + 13] = i;
	  Tmap[ind1][15*ind + 14] = count;
	}
      }
    }
      
    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      T2[i].resize(Chan.hp1[i] * Chan.hp2[i]);
      T2T[i].resize(Chan.hp1[i] * Chan.hp2[i]);
      T3[i].resize(Chan.hp1[i] * Chan.hp2[i]);
      S4[i].resize(Chan.hp1[i] * Chan.hp1[i]);
      S4T[i].resize(Chan.hp1[i] * Chan.hp1[i]);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind1 = Chan.hp1vec1[i][2*j];
	pind1 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind2 = Chan.hp2vec1[i][2*k];
	  pind2 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[pind1]+Space.levelsl[pind2]), Space.levelsm[pind1]+Space.levelsm[pind2], 
			   Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 1] = i;
	  Tmap[ind1][15*ind + 2] = count;
	}
      }
    }

    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind2 = Chan.hp1vec1[i][2*j];
	pind2 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind1 = Chan.hp2vec1[i][2*k];
	  pind1 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[pind1]+Space.levelsl[pind2]), Space.levelsm[pind1]+Space.levelsm[pind2], 
			   Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 3] = i;
	  Tmap[ind1][15*ind + 4] = count;
	}
      }
    }
    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind2 = Chan.hp1vec1[i][2*j];
	pind1 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind1 = Chan.hp2vec1[i][2*k];
	  pind2 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[pind1]+Space.levelsl[pind2]), Space.levelsm[pind1]+Space.levelsm[pind2], 
			   Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 5] = i;
	  Tmap[ind1][15*ind + 6] = count;
	}
      }	
    }
  }
  else if(Parameters.basis == "CART"){
    #pragma omp parallel for
    for(int i = 0; i < Chan.size3; ++i){
      T4[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
      T5[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
      S2[i].assign(Chan.p[i] * Chan.p[i], 0.0);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int hhp = 0; hhp < Chan.hhp[i]; ++hhp){
	hind1 = Chan.hhpvec1[i][3*hhp];
	hind2 = Chan.hhpvec1[i][3*hhp + 1];
	pind1 = Chan.hhpvec1[i][3*hhp + 2];
	ind1 = CART_tbInd1(Space, Space.levelsnx[hind1]+Space.levelsnx[hind2], Space.levelsny[hind1]+Space.levelsny[hind2], 
			   Space.levelsnz[hind1]+Space.levelsnz[hind2], Space.levelsm[hind1]+Space.levelsm[hind2], Space.levelst[hind1]+Space.levelst[hind2]);
	for(int p2 = 0; p2 < Chan.p[i]; ++p2){
	  ++count;
	  pind2 = Chan.pvec1[i][p2];
	  if(pind1 == pind2){ continue; }
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 7] = i;
	  Tmap[ind1][15*ind + 8] = count;
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind2, pind1);
	  Tmap[ind1][15*ind + 9] = i;
	  Tmap[ind1][15*ind + 10] = count;
	}
      }
    }
    #pragma omp parallel for
    for(int i = 0; i < Chan.size3; ++i){
      T6[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
      T7[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
      S3[i].assign(Chan.h[i] * Chan.h[i], 0.0);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int hpp = 0; hpp < Chan.hpp[i]; ++hpp){
	hind1 = Chan.hppvec1[i][3*hpp];
	pind1 = Chan.hppvec1[i][3*hpp + 1];
	pind2 = Chan.hppvec1[i][3*hpp + 2];
	ind1 = CART_tbInd1(Space, Space.levelsnx[pind1]+Space.levelsnx[pind2], Space.levelsny[pind1]+Space.levelsny[pind2], 
			   Space.levelsnz[pind1]+Space.levelsnz[pind2], Space.levelsm[pind1]+Space.levelsm[pind2], Space.levelst[pind1]+Space.levelst[pind2]);
	for(int h2 = 0; h2 < Chan.h[i]; ++h2){
	  ++count;
	  hind2 = Chan.hvec1[i][h2];
	  if(hind1 == hind2){ continue; }
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 11] = i;
	  Tmap[ind1][15*ind + 12] = count;
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind2, hind1, pind1, pind2);
	  Tmap[ind1][15*ind + 13] = i;
	  Tmap[ind1][15*ind + 14] = count;
	}
      }
    }
    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      T2[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      T2T[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      T3[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
      S4[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
      S4T[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind1 = Chan.hp1vec1[i][2*j];
	pind1 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind2 = Chan.hp2vec1[i][2*k];
	  pind2 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = CART_tbInd1(Space, Space.levelsnx[pind1]+Space.levelsnx[pind2], Space.levelsny[pind1]+Space.levelsny[pind2], 
			     Space.levelsnz[pind1]+Space.levelsnz[pind2], Space.levelsm[pind1]+Space.levelsm[pind2], Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 1] = i;
	  Tmap[ind1][15*ind + 2] = count;
	}
      }
    }
    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind2 = Chan.hp1vec1[i][2*j];
	pind2 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind1 = Chan.hp2vec1[i][2*k];
	  pind1 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = CART_tbInd1(Space, Space.levelsnx[pind1]+Space.levelsnx[pind2], Space.levelsny[pind1]+Space.levelsny[pind2], 
			     Space.levelsnz[pind1]+Space.levelsnz[pind2], Space.levelsm[pind1]+Space.levelsm[pind2], Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 3] = i;
	  Tmap[ind1][15*ind + 4] = count;
	}
      }
    }
    #pragma omp parallel for
    for(int i = 0; i < Chan.size2; ++i){
      int count = -1;
      int ind, ind1;
      int hind1, hind2, pind1, pind2;
      for(int j = 0; j < Chan.hp1[i]; ++j){
	hind2 = Chan.hp1vec1[i][2*j];
	pind1 = Chan.hp1vec1[i][2*j + 1];
	for(int k = 0; k < Chan.hp2[i]; ++k){
	  ++count;
	  hind1 = Chan.hp2vec1[i][2*k];
	  pind2 = Chan.hp2vec1[i][2*k + 1];
	  if(hind1 == hind2 || pind1 == pind2){ continue; }
	  ind1 = CART_tbInd1(Space, Space.levelsnx[pind1]+Space.levelsnx[pind2], Space.levelsny[pind1]+Space.levelsny[pind2], 
			     Space.levelsnz[pind1]+Space.levelsnz[pind2], Space.levelsm[pind1]+Space.levelsm[pind2], Space.levelst[pind1]+Space.levelst[pind2]);
	  ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	  Tmap[ind1][15*ind + 5] = i;
	  Tmap[ind1][15*ind + 6] = count;
	}
      }	
    }
  }

      /*std::ofstream mapfile;
  mapfile.open ("mapfileN4.txt");
  for(int i = 0; i < Chan.size1; ++i)
    {
      mapfile << Chan.hh[i] * Chan.pp[i] << " ";
    }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size2; ++i)
    {
      mapfile << Chan.hp1[i] * Chan.hp2[i] << " ";
    }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i)
    {
      mapfile << Chan.h[i] * Chan.hpp[i] << " ";
    }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size3; ++i)
    {
      mapfile << Chan.p[i] * Chan.hhp[i] << " ";
    }
  mapfile << std::endl;
  for(int i = 0; i < Chan.size1; ++i)
    {
      for(int j = 0; j < Chan.hh[i]*Chan.pp[i]; ++j)
	{
	  mapfile << i << " ";
	  for(int k = 0; k < 15; ++k)
	    {
	      mapfile << Tmap[i][15*j + k] << " ";
	    }
	  mapfile << std::endl;
	}
    }
    mapfile.close();*/


}

void CCD::set_T(int i, int j, double T)
{
  T1[i][Tmap[i][15*j]] = T;
  T2[Tmap[i][15*j + 1]][Tmap[i][15*j + 2]] = T;
  T2T[Tmap[i][15*j + 3]][Tmap[i][15*j + 4]] = T;
  T3[Tmap[i][15*j + 5]][Tmap[i][15*j + 6]] = T;
  T4[Tmap[i][15*j + 7]][Tmap[i][15*j + 8]] = T;
  T5[Tmap[i][15*j + 9]][Tmap[i][15*j + 10]] = T;
  T6[Tmap[i][15*j + 11]][Tmap[i][15*j + 12]] = T;
  T7[Tmap[i][15*j + 13]][Tmap[i][15*j + 14]] = T;
}

double CCD::get_T(int i, int j) const
{
  double tempt;
  tempt = T1[i][Tmap[i][15*j]];
  tempt += T2[Tmap[i][15*j + 1]][Tmap[i][15*j + 2]];    
  tempt += T3[Tmap[i][15*j + 5]][Tmap[i][15*j + 6]];
  tempt += T4[Tmap[i][15*j + 7]][Tmap[i][15*j + 8]];
  tempt -= T4[Tmap[i][15*j + 9]][Tmap[i][15*j + 10]];
  tempt += T6[Tmap[i][15*j + 11]][Tmap[i][15*j + 12]];
  tempt -= T6[Tmap[i][15*j + 13]][Tmap[i][15*j + 14]];
  return tempt;
}


//Function to setup Channels
Channels HO_Setup_Channels(const Model_Space &Space)
{
  std::cout << "Building Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  Channels Chan;
  double M, T, P;

  std::vector<double> HPMs;
  std::vector<double> HPTs;
  std::vector<double> HPPs;
  std::vector<std::vector<int> > Hvec, Pvec;

  Chan.size3 = 0;
  for(int i = 0; i < Space.indtot; ++i){
    M = Space.levelsm[i];
    T = Space.levelst[i];
    P = std::pow(-1.0, Space.levelsl[i]);
    if(Chan.size3 == 0){
      ++Chan.size3;
      Hvec.resize(Chan.size3);
      Pvec.resize(Chan.size3);
      HPMs.push_back(M);
      HPTs.push_back(T);
      HPPs.push_back(P);
      Chan.indvec.push_back(Chan.size3 - 1);
      if(Space.levelstype[i] == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
      else{ Pvec[Chan.size3 - 1].push_back(i); }
      continue;
    }
    else{
      for(int k = 0; k < Chan.size3; ++k){
	if(M == HPMs[k] && T == HPTs[k] && P == HPPs[k]){
	  Chan.indvec.push_back(k);
	  if(Space.levelstype[i] == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  break;
	}
	else if(k == Chan.size3 - 1){
	  ++Chan.size3;
	  Hvec.resize(Chan.size3);
	  Pvec.resize(Chan.size3);
	  HPMs.push_back(M);
	  HPTs.push_back(T);
	  HPPs.push_back(P);
	  Chan.indvec.push_back(Chan.size3 - 1);
	  if(Space.levelstype[i] == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
	  else{ Pvec[Chan.size3 - 1].push_back(i); }
	  break;
	}
      }
    }
  }

  Chan.size1 = Space.HO_tb1Indsize * Space.Msize * Space.Tsize;
  Chan.size2 = Space.HO_tb2Indsize * Space.M2size * Space.T2size;

  std::vector<std::vector<int> > tempvec1(Chan.size3);
  std::vector<std::vector<int> > tempvec2(Chan.size3);

  //Make Vector of Two-Body States for Each Channel
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2;
    tempvec1[i].assign(Chan.size1, -1);
    tempvec2[i].assign(Chan.size2, -1);
    double M, T, M2, T2, P, P2;
    for(int j = 0; j < Chan.size3; ++j){
      M = HPMs[i] + HPMs[j];
      T = HPTs[i] + HPTs[j];
      P = HPPs[i] * HPPs[j];
      M2 = HPMs[i] - HPMs[j];
      T2 = HPTs[i] - HPTs[j];
      P2 = HPPs[i] / HPPs[j];
      ind1 = HO_tbInd1(Space, P, M, T);
      tempvec1[i][ind1] = j;
      ind2 = HO_tbInd2(Space, P2, M2, T2);
      tempvec2[i][ind2] = j;
    }
  }

  Chan.hhvec1.resize(Chan.size1);
  Chan.ppvec1.resize(Chan.size1);
  Chan.hp1vec1.resize(Chan.size2);
  Chan.hp2vec1.resize(Chan.size2);
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
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
	  if(ind1 == ind2 && h1 == h2){ continue; }
	  Chan.hhvec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hhvec1[j].push_back(Hvec[ind2][h2]);
	}
      }
      psize1 = int(Pvec[ind1].size());
      psize2 = int(Pvec[ind2].size());
      for(int p1 = 0; p1 < psize1; ++p1){
	for(int p2 = 0; p2 < psize2; ++p2){
	  if(ind1 == ind2 && p1 == p2){ continue; }
	  Chan.ppvec1[j].push_back(Pvec[ind1][p1]);
	  Chan.ppvec1[j].push_back(Pvec[ind2][p2]);
	}
      }
    }
    #pragma omp parallel for
    for(int j = 0; j < Chan.size2; ++j){
      if(tempvec2[i][j] == -1){ continue; }
      int ind1, ind2;
      int hsize1, psize1;
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
    }
  }

  Chan.hvec1.resize(Chan.size3);
  Chan.hppvec1.resize(Chan.size3);
  Chan.pvec1.resize(Chan.size3);
  Chan.hhpvec1.resize(Chan.size3);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2, ind3, ind4;
    int hsize1, hsize2, psize1, psize2;
    std::vector<int> inds(3);
    hsize1 = int(Hvec[i].size());
    if(hsize1 == 0){ goto break1; }
    for(int h1 = 0; h1 < hsize1; ++h1){ Chan.hvec1[i].push_back(Hvec[i][h1]); }
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      for(int k = 0; k < Chan.size3; ++k){
	if(tempvec1[k][j] == -1){ continue; }
	ind2 = tempvec1[i][j];
	ind3 = k;
	ind4 = tempvec1[k][j];
	hsize2 = int(Hvec[ind2].size());
	psize1 = int(Pvec[ind3].size());
	psize2 = int(Pvec[ind4].size());
	for(int h2 = 0; h2 < hsize2; ++h2){
	  for(int p1 = 0; p1 < psize1; ++p1){
	    for(int p2 = 0; p2 < psize2; ++p2){
	      if(ind3 == ind4 && p1 == p2){ continue; }
	      inds[0] = Hvec[ind2][h2];
	      inds[1] = Pvec[ind3][p1];
	      inds[2] = Pvec[ind4][p2];
	      Chan.hppvec1[i].insert(Chan.hppvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
      }
    }
  break1:
    psize1 = int(Pvec[i].size());
    if(psize1 == 0){ continue; }
    for(int p1 = 0; p1 < psize1; ++p1){ Chan.pvec1[i].push_back(Pvec[i][p1]); }
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      for(int k = 0; k < Chan.size3; ++k){
	if(tempvec1[k][j] == -1){ continue; }
	ind1 = k;
	ind2 = tempvec1[k][j];
	ind4 = tempvec1[i][j];
	hsize1 = int(Hvec[ind1].size());
	hsize2 = int(Hvec[ind2].size());
	psize2 = int(Pvec[ind4].size());
	for(int h1 = 0; h1 < hsize1; ++h1){
	  for(int h2 = 0; h2 < hsize2; ++h2){
	    if(ind1 == ind2 && h1 == h2){ continue; }
	    for(int p2 = 0; p2 < psize2; ++p2){
	      inds[0] = Hvec[ind1][h1];
	      inds[1] = Hvec[ind2][h2];
	      inds[2] = Pvec[ind4][p2];
	      Chan.hhpvec1[i].insert(Chan.hhpvec1[i].end(), inds.begin(), inds.end());
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
  memory += (7.0 * Chan.size1 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    Chan.hh[i] = int(Chan.hhvec1[i].size()/2);
    Chan.pp[i] = int(Chan.ppvec1[i].size()/2);
    mem1[i] = 8.0 * (2.0 * Chan.hh[i]*Chan.hh[i] + Chan.pp[i]*Chan.pp[i] + 3.0 * Chan.hh[i]*Chan.pp[i]) + 4.0 * (15.0 * Chan.hh[i]*Chan.pp[i]);
  }
  Chan.hp1.resize(Chan.size2);
  Chan.hp2.resize(Chan.size2);
  memory += (9.0 * Chan.size2 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    Chan.hp1[i] = int(Chan.hp1vec1[i].size()/2);
    Chan.hp2[i] = int(Chan.hp2vec1[i].size()/2);
    mem2[i] = 8.0 * (5.0 * Chan.hp1[i]*Chan.hp2[i] + 3.0 * Chan.hp1[i]*Chan.hp1[i] + Chan.hp2[i]*Chan.hp2[i]);
  }
  Chan.h.resize(Chan.size3);
  Chan.hpp.resize(Chan.size3);
  Chan.p.resize(Chan.size3);
  Chan.hhp.resize(Chan.size3);
  memory += (8.0 * Chan.size3 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    Chan.h[i] = int(Chan.hvec1[i].size());
    Chan.hpp[i] = int(Chan.hppvec1[i].size()/3);
    Chan.p[i] = int(Chan.pvec1[i].size());
    Chan.hhp[i] = int(Chan.hhpvec1[i].size()/3);
    mem3[i] = 8.0 * (Chan.h[i]*Chan.h[i] + 3.0 * Chan.h[i]*Chan.hpp[i] + Chan.p[i]*Chan.p[i] + 3.0 * Chan.p[i]*Chan.hhp[i]);
  }
  
  for(int i = 0; i < Chan.size1; ++i){ memory += mem1[i]; }
  for(int i = 0; i < Chan.size2; ++i){ memory += mem2[i]; }
  for(int i = 0; i < Chan.size3; ++i){ memory += mem3[i]; }
  std::cout << "2B1 Channels = " << Chan.size1 << ", 2B2 Channels = " << Chan.size2 << ", OB Channels = " << Chan.size3 << std::endl; 
  std::cout << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  return Chan;

}


//Function to setup Channels
Channels CART_Setup_Channels(const Model_Space &Space)
{  
  std::cout << "Building Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  Channels Chan;
  double M, T;
  int Nx, Ny, Nz;

  std::vector<double> HPMs;
  std::vector<double> HPTs;
  std::vector<int> HPNxs;
  std::vector<int> HPNys;
  std::vector<int> HPNzs;
  std::vector<std::vector<int> > Hvec, Pvec;

  Chan.size3 = 0;
  for(int i = 0; i < Space.indtot; ++i){
    M = Space.levelsm[i];
    T = Space.levelst[i];
    Nx = Space.levelsnx[i];
    Ny = Space.levelsny[i];
    Nz = Space.levelsnz[i];
    if(Chan.size3 == 0){
      ++Chan.size3;
      Hvec.resize(Chan.size3);
      Pvec.resize(Chan.size3);
      HPMs.push_back(M);
      HPTs.push_back(T);
      HPNxs.push_back(Nx);
      HPNys.push_back(Ny);
      HPNzs.push_back(Nz);
      Chan.indvec.push_back(Chan.size3 - 1);
      if(Space.levelstype[i] == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
      else{ Pvec[Chan.size3 - 1].push_back(i); }
      continue;
    }
    else{
      for(int k = 0; k < Chan.size3; ++k){
	if(M == HPMs[k] && T == HPTs[k] && Nx == HPNxs[k] && Ny == HPNys[k] && Nz == HPNzs[k]){ 
	  Chan.indvec.push_back(k);
	  if(Space.levelstype[i] == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  break;
	}
	else if(k == Chan.size3 - 1){
	  ++Chan.size3;
	  Hvec.resize(Chan.size3);
	  Pvec.resize(Chan.size3);
	  HPMs.push_back(M);
	  HPTs.push_back(T);
	  HPNxs.push_back(Nx);
	  HPNys.push_back(Ny);
	  HPNzs.push_back(Nz);
	  Chan.indvec.push_back(Chan.size3 - 1);
	  if(Space.levelstype[i] == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
	  else{ Pvec[Chan.size3 - 1].push_back(i); }
	  break;
	}
      }
    }
  }

  Chan.size1 = Space.CART_tb1Indsize * Space.Msize * Space.Tsize;
  Chan.size2 = Space.CART_tb2Indsize * Space.M2size * Space.T2size;

  std::vector<std::vector<int> > tempvec1(Chan.size3);
  std::vector<std::vector<int> > tempvec2(Chan.size3);

  //Make Vector of Two-Body States for Each Channel
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2;
    tempvec1[i].assign(Chan.size1, -1);
    tempvec2[i].assign(Chan.size2, -1);
    double M, T, M2, T2;
    int Nx, Ny, Nz, Nx2, Ny2, Nz2;
    for(int j = 0; j < Chan.size3; ++j){
      M = HPMs[i] + HPMs[j];
      T = HPTs[i] + HPTs[j];
      M2 = HPMs[i] - HPMs[j];
      T2 = HPTs[i] - HPTs[j];
      Nx = HPNxs[i] + HPNxs[j];
      Ny = HPNys[i] + HPNys[j];
      Nz = HPNzs[i] + HPNzs[j];
      Nx2 = HPNxs[i] - HPNxs[j];
      Ny2 = HPNys[i] - HPNys[j];
      Nz2 = HPNzs[i] - HPNzs[j];
      ind1 = CART_tbInd1(Space, Nx, Ny, Nz, M, T);
      tempvec1[i][ind1] = j;
      ind2 = CART_tbInd2(Space, Nx2, Ny2, Nz2, M2, T2);
      tempvec2[i][ind2] = j;
    }
  }

  Chan.hhvec1.resize(Chan.size1);
  Chan.ppvec1.resize(Chan.size1);
  Chan.hp1vec1.resize(Chan.size2);
  Chan.hp2vec1.resize(Chan.size2);
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
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
	  if(ind1 == ind2 && h1 == h2){ continue; }
	  Chan.hhvec1[j].push_back(Hvec[ind1][h1]);
	  Chan.hhvec1[j].push_back(Hvec[ind2][h2]);
	}
      }
      psize1 = int(Pvec[ind1].size());
      psize2 = int(Pvec[ind2].size());
      for(int p1 = 0; p1 < psize1; ++p1){
	for(int p2 = 0; p2 < psize2; ++p2){
	  if(ind1 == ind2 && p1 == p2){ continue; }
	  Chan.ppvec1[j].push_back(Pvec[ind1][p1]);
	  Chan.ppvec1[j].push_back(Pvec[ind2][p2]);
	}
      }
    }
    #pragma omp parallel for
    for(int j = 0; j < Chan.size2; ++j){
      if(tempvec2[i][j] == -1){ continue; }
      int ind1, ind2;
      int hsize1, psize1;
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
    }
  }
  Chan.hvec1.resize(Chan.size3);
  Chan.hppvec1.resize(Chan.size3);
  Chan.pvec1.resize(Chan.size3);
  Chan.hhpvec1.resize(Chan.size3);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2, ind3, ind4;
    int hsize1, hsize2, psize1, psize2;
    std::vector<int> inds(3);
    hsize1 = int(Hvec[i].size());
    if(hsize1 == 0){ goto break1; }
    for(int h1 = 0; h1 < hsize1; ++h1){ Chan.hvec1[i].push_back(Hvec[i][h1]); }
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      for(int k = 0; k < Chan.size3; ++k){
	if(tempvec1[k][j] == -1){ continue; }
	ind2 = tempvec1[i][j];
	ind3 = k;
	ind4 = tempvec1[k][j];
	hsize2 = int(Hvec[ind2].size());
	psize1 = int(Pvec[ind3].size());
	psize2 = int(Pvec[ind4].size());
	for(int h2 = 0; h2 < hsize2; ++h2){
	  for(int p1 = 0; p1 < psize1; ++p1){
	    for(int p2 = 0; p2 < psize2; ++p2){
	      if(ind3 == ind4 && p1 == p2){ continue; }
	      inds[0] = Hvec[ind2][h2];
	      inds[1] = Pvec[ind3][p1];
	      inds[2] = Pvec[ind4][p2];
	      Chan.hppvec1[i].insert(Chan.hppvec1[i].end(), inds.begin(), inds.end());
	    }
	  }
	}
      }
    }
  break1:
    psize1 = int(Pvec[i].size());
    if(psize1 == 0){ continue; }
    for(int p1 = 0; p1 < psize1; ++p1){ Chan.pvec1[i].push_back(Pvec[i][p1]); }
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      for(int k = 0; k < Chan.size3; ++k){
	if(tempvec1[k][j] == -1){ continue; }
	ind1 = k;
	ind2 = tempvec1[k][j];
	ind4 = tempvec1[i][j];
	hsize1 = int(Hvec[ind1].size());
	hsize2 = int(Hvec[ind2].size());
	psize2 = int(Pvec[ind4].size());
	for(int h1 = 0; h1 < hsize1; ++h1){
	  for(int h2 = 0; h2 < hsize2; ++h2){
	    if(ind1 == ind2 && h1 == h2){ continue; }
	    for(int p2 = 0; p2 < psize2; ++p2){
	      inds[0] = Hvec[ind1][h1];
	      inds[1] = Hvec[ind2][h2];
	      inds[2] = Pvec[ind4][p2];
	      Chan.hhpvec1[i].insert(Chan.hhpvec1[i].end(), inds.begin(), inds.end());
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
  memory += (7.0 * Chan.size1 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    Chan.hh[i] = int(Chan.hhvec1[i].size()/2);
    Chan.pp[i] = int(Chan.ppvec1[i].size()/2);
    mem1[i] = 8.0 * (2.0 * Chan.hh[i]*Chan.hh[i] + Chan.pp[i]*Chan.pp[i] + 3.0 * Chan.hh[i]*Chan.pp[i]) + 4.0 * (15.0 * Chan.hh[i]*Chan.pp[i]);
  }
  Chan.hp1.resize(Chan.size2);
  Chan.hp2.resize(Chan.size2);
  memory += (9.0 * Chan.size2 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    Chan.hp1[i] = int(Chan.hp1vec1[i].size()/2);
    Chan.hp2[i] = int(Chan.hp2vec1[i].size()/2);
    mem2[i] = 8.0 * (5.0 * Chan.hp1[i]*Chan.hp2[i] + 3.0 * Chan.hp1[i]*Chan.hp1[i] + Chan.hp2[i]*Chan.hp2[i]);
  }
  Chan.h.resize(Chan.size3);
  Chan.hpp.resize(Chan.size3);
  Chan.p.resize(Chan.size3);
  Chan.hhp.resize(Chan.size3);
  memory += (8.0 * Chan.size3 + 1) * 24.0;
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    Chan.h[i] = int(Chan.hvec1[i].size());
    Chan.hpp[i] = int(Chan.hppvec1[i].size()/3);
    Chan.p[i] = int(Chan.pvec1[i].size());
    Chan.hhp[i] = int(Chan.hhpvec1[i].size()/3);
    mem3[i] = 8.0 * (Chan.h[i]*Chan.h[i] + 3.0 * Chan.h[i]*Chan.hpp[i] + Chan.p[i]*Chan.p[i] + 3.0 * Chan.p[i]*Chan.hhp[i]);
  }
  
  for(int i = 0; i < Chan.size1; ++i){ memory += mem1[i]; }
  for(int i = 0; i < Chan.size2; ++i){ memory += mem2[i]; }
  for(int i = 0; i < Chan.size3; ++i){ memory += mem3[i]; }
  std::cout << "2B1 Channels = " << Chan.size1 << ", 2B2 Channels = " << Chan.size2 << ", OB Channels = " << Chan.size3 << std::endl; 
  std::cout << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  return Chan;

}


CC_Matrix_Elements::CC_Matrix_Elements(Channels Chan)
{
  HHHH.resize(Chan.size1);
  PPPP.resize(Chan.size1);
  HPHP1.resize(Chan.size2);
  HPHP2.resize(Chan.size2);
  HHPP1.resize(Chan.size1);
  HHPP2.resize(Chan.size3);
  HHPP3.resize(Chan.size3);
  HHPP4.resize(Chan.size2);
  HHPP4T.resize(Chan.size2);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    HHHH[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    PPPP[i].assign(Chan.pp[i] * Chan.pp[i], 0.0);
    HHPP1[i].assign(Chan.pp[i] * Chan.hh[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    HHPP2[i].assign(Chan.p[i] * Chan.hhp[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    HHPP3[i].assign(Chan.h[i] * Chan.hpp[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    HPHP1[i].assign(Chan.hp2[i] * Chan.hp2[i], 0.0);
    HPHP2[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
    HHPP4[i].assign(Chan.hp1[i] * Chan.hp2[i], 0.0);
    HHPP4T[i].assign(Chan.hp2[i] * Chan.hp1[i], 0.0);
  }
}


CC_Eff::CC_Eff(Channels Chan)
{
  V1.resize(Chan.size1);
  V2.resize(Chan.size1);
  V3.resize(Chan.size2);
  V4.resize(Chan.size2);
  V5.resize(Chan.size1);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    V1[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    V2[i].assign(Chan.h[i] * Chan.h[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    V3[i].assign(Chan.pp[i] * Chan.pp[i], 0.0);
    V4[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
  }
  #pragma omp parallel for
  for(int i = 0; i < Chan.size2; ++i){
    V5[i].assign(Chan.hp1[i] * Chan.hp1[i], 0.0);
  }
}


//Initialize program from input file
Input_Parameters Get_Input_Parameters(std::string &infile)
{ 
  Input_Parameters Input; // Input Parameters
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
	Input.basis = substr;
	break;
      case 2:
	Input.obstrength = atof(substr.c_str());
	break;
      case 3:
	Input.tbstrength = atof(substr.c_str());
	break;
      case 4:
	Input.P = atoi(substr.c_str());
	break;
      case 5:
	Input.N = atoi(substr.c_str());
	break;
      case 6:
	Input.LevelScheme = substr;
	break;
      case 7:
	Input.MatrixElements = substr;
	break;
      } 
    }
    else{ continue; };
  }
  return Input;
}


void Print_Parameters(const Input_Parameters &Parameters)
{
  std::cout << "----------------------------------------------------------" << std::endl;
  std::cout << "Basis = " << Parameters.basis;
  if(Parameters.LevelScheme.size() > 0){ std::cout << ", Levels Scheme = " << Parameters.LevelScheme << ", Interaction = " << Parameters.MatrixElements << std::endl; }
  else{ std::cout << ", Nmax = " << Parameters.Nmax << ", Density = " << Parameters.density << std::endl; }
  std::cout << "OB strength = " << Parameters.obstrength << ", TB strength = " << Parameters.tbstrength << std::endl;
  std::cout << "Protons = " << Parameters.P << ", Neutron = " << Parameters.N << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
}


Model_Space Build_Model_Space(const Input_Parameters &Parameters)
{
  Model_Space Space; // Model Space information
  std::string fullpath; // Model Space file path
  std::string phline; // Std::String for each file line
  std::ifstream splevels; // Model space file
  int ind, n, l, nx, ny, nz;
  double energy, m, tz; // initialize level index, n, l, nx, ny, nz, m, tz depending on basis
  int pcount, ncount, holcount, parcount;
  
  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  getline(splevels, phline);
  std::stringstream(phline) >> Space.indtot;

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);

  if(Parameters.basis == "HO"){
    Space.levelsn.resize(Space.indtot);
    Space.levelsl.resize(Space.indtot);
  }
  else if(Parameters.basis == "CART"){
    Space.levelsnx.resize(Space.indtot);
    Space.levelsny.resize(Space.indtot);
    Space.levelsnz.resize(Space.indtot);
  }

  pcount = 0;
  ncount = 0;
  parcount = 0;
  holcount = 0;

  double Mmin = 1000, Tmin = 1000, Pmin = 1000;
  double Mmax = -1000, Tmax = -1000, Pmax = -1000;
  int Nxmin = 1000, Nymin = 1000, Nzmin = 1000;
  int Nxmax = -1000, Nymax = -1000, Nzmax = -1000;

  for(int i = 0; i < Space.indtot; ++i){
    if(Parameters.basis == "HO"){
      getline(splevels, phline);
      std::stringstream(phline) >> ind >> n >> l >> m >> tz >> energy;
      Space.levelsind[i] = ind;
      Space.levelsn[i] = n;
      Space.levelsl[i] = l;
      Space.levelsm[i] = m;
      Space.levelst[i] = tz;
      Space.levelsen[i] = energy * Parameters.obstrength;
      if(m < Mmin){ Mmin = m; }
      if(tz < Tmin){ Tmin = tz; }
      if(std::pow(-1.0, l) < Pmin){ Pmin = std::pow(-1.0, l); }
      if(m > Mmax){ Mmax = m; }
      if(tz > Tmax){ Tmax = tz; }
      if(std::pow(-1.0, l) > Pmax){ Pmax = std::pow(-1.0, l); }
    }
    else if(Parameters.basis == "CART"){
      getline(splevels, phline);
      std::stringstream(phline) >> ind >> nx >> ny >> nz >> m >> tz >> energy;
      Space.levelsind[i] = ind;
      Space.levelsnx[i] = nx;
      Space.levelsny[i] = ny;
      Space.levelsnz[i] = nz;
      Space.levelsm[i] = m;
      Space.levelst[i] = tz;
      Space.levelsen[i] = energy * Parameters.obstrength;
      if(m < Mmin){ Mmin = m; }
      if(tz < Tmin){ Tmin = tz; }
      if(nx < Nxmin){ Nxmin = nx; }
      if(ny < Nymin){ Nymin = ny; }
      if(nz < Nzmin){ Nzmin = nz; }
      if(m > Mmax){ Mmax = m; }
      if(tz > Tmax){ Tmax = tz; }
      if(nx > Nxmax){ Nxmax = nx; }
      if(ny > Nymax){ Nymax = ny; }
      if(nz > Nzmax){ Nzmax = nz; }
    }
    if(tz == -1 && pcount < Parameters.P){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; }
    else if(tz == -1 && pcount >= Parameters.P){ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    if(tz == 1 && ncount < Parameters.N){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; }
    else if(tz == 1 && ncount >= Parameters.N){ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;

  if(Parameters.basis == "HO"){
    int count1 = 0;
    int Psize = int((Pmax - Pmin)/2) + 1;
    Space.Pmin = Pmin;
    Space.HO_tb1Indvec.resize(Psize);
    for(int i = 0; i < Psize; ++i){
      Space.HO_tb1Indvec[i] = count1;
      ++count1;
    }
    Space.HO_tb1Indsize = count1;
    int count2 = 0;
    int P2size = int((Pmax - Pmin)/2) + 1;
    Space.P2min = Pmin;
    Space.HO_tb2Indvec.resize(P2size);
    for(int i = 0; i < P2size; ++i){
      Space.HO_tb2Indvec[i] = count2;
      ++count2;
    }
    Space.HO_tb2Indsize = count2;
  }
  else if(Parameters.basis == "CART"){
    int count1 = 0;
    int Nxsize = int(2.0 * (Nxmax - Nxmin) + 1);
    int Nysize, Nzsize;
    Space.Nxmin = Nxmin + Nxmin;
    Space.Nymin.resize(Nxsize);
    Space.Nzmin.resize(Nxsize);
    Space.CART_tb1Indvec.resize(Nxsize);
    for(int i = 0; i < Nxsize; ++i){
      Space.Nymin[i] = Nymin + Nymin;
      Nysize = int(2.0 * (Nymax - Nymin) + 1);
      Space.Nzmin[i].resize(Nysize);
      Space.CART_tb1Indvec[i].resize(Nysize);
      for(int j = 0; j < Nysize; ++j){
	Space.Nzmin[i][j] = Nzmin + Nzmin;
	Nzsize = int(2.0 * (Nzmax - Nzmin) + 1);
	Space.CART_tb1Indvec[i][j].resize(Nzsize);
	for(int k = 0; k < Nzsize; ++k){
	  Space.CART_tb1Indvec[i][j][k] = count1;
	  ++count1;
	}
      }
    }
    Space.CART_tb1Indsize = count1;
    int count2 = 0;
    int Nx2size = int(2.0 * (Nxmax - Nxmin) + 1);
    int Ny2size, Nz2size;
    Space.Nx2min = Nxmin - Nxmax;
    Space.Ny2min.resize(Nx2size);
    Space.Nz2min.resize(Nx2size);
    Space.CART_tb2Indvec.resize(Nx2size);
    for(int i = 0; i < Nx2size; ++i){
      Space.Ny2min[i] = Nymin - Nymax;
      Ny2size = int(2.0 * (Nymax - Nymin) + 1);
      Space.Nz2min[i].resize(Ny2size);
      Space.CART_tb2Indvec[i].resize(Ny2size);
      for(int j = 0; j < Ny2size; ++j){
	Space.Nz2min[i][j] = Nzmin - Nzmax;
	Nz2size = int(2.0 * (Nzmax - Nzmin) + 1);
	Space.CART_tb2Indvec[i][j].resize(Nz2size);
	for(int k = 0; k < Nz2size; ++k){
	  Space.CART_tb2Indvec[i][j][k] = count2;
	  ++count2;
	}
      }
    }
    Space.CART_tb2Indsize = count2;
  }

  Space.Mmin = Mmin + Mmin;
  Space.Msize = int(Mmax - Mmin + 1);
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.Tmin = -2;
    Space.Tsize = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.Tmin = -2;
    Space.Tsize = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.Tmin = 2;
    Space.Tsize = 1;
  }

  Space.M2min = Mmin - Mmax;
  Space.M2size = int(Mmax - Mmin + 1);
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.T2min = -2;
    Space.T2size = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.T2min = 0;
    Space.T2size = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.T2min = 0;
    Space.T2size = 1;
  }

  splevels.close();
  
  return Space;
}


Model_Space CART_Build_Model_Space(const Input_Parameters &Parameters)
{
  Model_Space Space; // Model Space information
  double E;
  int count, pcount, ncount, holcount, parcount;
  
  double hbarc = 197.3269788; // MeVfm
  double m_neutronc2 = 939.565378; // MeV
  //double m_protonc2 = 938.272046; // MeV
  double m_protonc2 = 939.565378; // MeV
  double neutron_prefac = hbarc*hbarc/(2.*m_neutronc2);
  double proton_prefac = hbarc*hbarc/(2.*m_protonc2);
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);

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

  int shellMax = 3*Parameters.Nmax*Parameters.Nmax;
  for(int shell = 0; shell <= shellMax; ++shell){
    for(int nx = -Parameters.Nmax; nx <= Parameters.Nmax; ++nx){    
      for(int ny = -Parameters.Nmax; ny <= Parameters.Nmax; ++ny){	
	for(int nz = -Parameters.Nmax; nz <= Parameters.Nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz || Parameters.Nmax < nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = -1; tz <= 1; tz = tz+2){
	      if(tz == -1){
		if(Parameters.P == 0){ continue; }
		E = (proton_prefac*4*M_PI*M_PI/(L*L)) * (nx*nx + ny*ny + nz*nz);
		Space.levelsen[count] = E;
		if(pcount < Parameters.P){ Space.levelstype[count] = "hole"; ++pcount; ++holcount; }
		else{ Space.levelstype[count] = "particle"; ++pcount; ++parcount; }
	      }
	      if(tz == 1){
		if(Parameters.N == 0){ continue; }
		E = (neutron_prefac*4*M_PI*M_PI/(L*L)) * (nx*nx + ny*ny + nz*nz);
		Space.levelsen[count] = E;
		if(ncount < Parameters.N){ Space.levelstype[count] = "hole"; ++ncount; ++holcount; }
		else{ Space.levelstype[count] = "particle"; ++ncount; ++parcount; }
	      }
	      //std::cout << count+1 << " " << nx << " " << ny << " " << nz << " " << double(sz) << " " << double(tz) << " " << Space.levelstype[count] << std::endl;
	      Space.levelsind[count] = count + 1;
	      Space.levelsnx[count] = nx;
	      Space.levelsny[count] = ny;
	      Space.levelsnz[count] = nz;
	      Space.levelsm[count] = double(sz);
	      Space.levelst[count] = double(tz);
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;

  Space.levelsind.resize(count);
  Space.levelsm.resize(count);
  Space.levelst.resize(count);
  Space.levelstype.resize(count);
  Space.levelsen.resize(count);
  Space.levelsnx.resize(count);
  Space.levelsny.resize(count);
  Space.levelsnz.resize(count);
  Space.indtot = count;

  Space.nmax = int(std::ceil(std::sqrt(Parameters.Nmax)));

  int count1 = 0;
  int Nxsize = 4*Space.nmax + 1;
  int Nysize, Nzsize;
  Space.Nxmin = -2*Space.nmax;
  Space.Nymin.resize(Nxsize);
  Space.Nzmin.resize(Nxsize);
  Space.CART_tb1Indvec.resize(Nxsize);
  for(int i = 0; i < Nxsize; ++i){
    Space.Nymin[i] = -int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i))));
    Nysize = 2*int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i)))) + 1;
    Space.Nzmin[i].resize(Nysize);
    Space.CART_tb1Indvec[i].resize(Nysize);
    for(int j = 0; j < Nysize; ++j){
      Space.Nzmin[i][j] = -int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j))));
      Nzsize = 2*int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nxmin+i)*(Space.Nxmin+i) - (Space.Nymin[i]+j)*(Space.Nymin[i]+j)))) + 1;
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
    Space.Ny2min[i] = -int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i))));
    Ny2size = 2*int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i)))) + 1;
    Space.Nz2min[i].resize(Ny2size);
    Space.CART_tb2Indvec[i].resize(Ny2size);
    for(int j = 0; j < Ny2size; ++j){
      Space.Nz2min[i][j] = -int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j))));
      Nz2size = 2*int(std::ceil(std::sqrt(4*Parameters.Nmax - (Space.Nx2min+i)*(Space.Nx2min+i) - (Space.Ny2min[i]+j)*(Space.Ny2min[i]+j)))) + 1;
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
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.Tmin = -2;
    Space.Tsize = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.Tmin = -2;
    Space.Tsize = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.Tmin = 2;
    Space.Tsize = 1;
  }

  Space.M2min = -2;
  Space.M2size = 3;
  if(Parameters.P != 0 && Parameters.N != 0){
    Space.T2min = -2;
    Space.T2size = 3;
  }
  else if(Parameters.P != 0 && Parameters.N == 0){
    Space.T2min = 0;
    Space.T2size = 1;
  }
  else if(Parameters.P == 0 && Parameters.N != 0){
    Space.T2min = 0;
    Space.T2size = 1;
  }

  return Space;
}



int HO_tbInd1(const Model_Space &Space, const double &P, const double &M, const double &T)
{
  return Space.HO_tb1Indvec[int((P - Space.Pmin)/2)]
    + Space.HO_tb1Indsize*(int((M - Space.Mmin)/2) * (Space.Tsize) + int((T - Space.Tmin)/2));
}

int HO_tbInd2(const Model_Space &Space, const double &P2, const double &M2, const double &T2)
{
  return Space.HO_tb2Indvec[int((P2 - Space.P2min)/2)]
    + Space.HO_tb2Indsize*(int((M2 - Space.M2min)/2) * (Space.T2size) + int((T2 - Space.T2min)/2));
}

int CART_tbInd1(const Model_Space &Space, const int &Nx, const int &Ny, const int &Nz, const double &M, const double &T)
{
  return Space.CART_tb1Indvec[Nx - Space.Nxmin][Ny - Space.Nymin[Nx - Space.Nxmin]][Nz - Space.Nzmin[Nx - Space.Nxmin][Ny - Space.Nymin[Nx - Space.Nxmin]]]
    + Space.CART_tb1Indsize*(int((M - Space.Mmin)/2) * (Space.Tsize) + int((T - Space.Tmin)/2));
}

int CART_tbInd2(const Model_Space &Space, const int &Nx2, const int &Ny2, const int &Nz2, const double &M2, const double &T2)
{
  return Space.CART_tb2Indvec[Nx2 - Space.Nx2min][Ny2 - Space.Ny2min[Nx2 - Space.Nx2min]][Nz2 - Space.Nz2min[Nx2 - Space.Nx2min][Ny2 - Space.Ny2min[Nx2 - Space.Nx2min]]]
    + Space.CART_tb2Indsize*(int((M2 - Space.M2min)/2) * (Space.T2size) + int((T2 - Space.T2min)/2));
}


CC_Matrix_Elements Read_Matrix_Elements(const std::string &MEfile, const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  int NumElements; // number of ME read in from file
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::stringstream interactionstream; // stream of file line std::string
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  std::string fullpath1;
  CC_Matrix_Elements CCME(Chan);
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;

  fullpath1 = PATH + MEfile + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << MEfile << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;
 
  if(Parameters.basis == "HO"){
    for(int i = 0; i < NumElements; ++i){
      getline(interaction, interactionline);
      std::stringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
      if(std::abs(TBME) < 1.0E-10){ continue; }
      shell1 -= 1; shell2 -= 1; shell3 -= 1; shell4 -= 1;
      TBME *= Parameters.tbstrength;
      ptype = Space.levelstype[shell1];
      qtype = Space.levelstype[shell2];
      rtype = Space.levelstype[shell3];
      stype = Space.levelstype[shell4];
      
      if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[shell1]+Space.levelsl[shell2]), Space.levelsm[shell1]+Space.levelsm[shell2], 
			 Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell3, shell4);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell3, shell4);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell4, shell3);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell4, shell3);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	CCME.HHHH[ind1][ind] = TBME;
      }
      else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){ 
	ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[shell1]+Space.levelsl[shell2]), Space.levelsm[shell1]+Space.levelsm[shell2], 
			 Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell3, shell4);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell3, shell4);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell4, shell3);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell4, shell3);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell1, shell2);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell1, shell2);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell2, shell1);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell2, shell1);
	CCME.PPPP[ind1][ind] = TBME;
      }
      else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell4]-Space.levelsl[shell1]), Space.levelsm[shell4]-Space.levelsm[shell1], 
			 Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
	CCME.HPHP1[ind1][ind] = TBME;
	ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
	CCME.HPHP1[ind1][ind] = TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell1]-Space.levelsl[shell4]), Space.levelsm[shell1]-Space.levelsm[shell4], 
			 Space.levelst[shell1]-Space.levelst[shell4]);
	ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell1, shell4, shell3, shell2);
	CCME.HPHP2[ind1][ind] = TBME;
	ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell3, shell2, shell1, shell4);
	CCME.HPHP2[ind1][ind] = TBME;
      }
      else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	ind1 = HO_tbInd1(Space, std::pow(-1.0, Space.levelsl[shell1]+Space.levelsl[shell2]), Space.levelsm[shell1]+Space.levelsm[shell2], 
			 Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	CCME.HHPP1[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	CCME.HHPP1[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	CCME.HHPP1[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	CCME.HHPP1[ind1][ind] = TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell3] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
	CCME.HHPP2[ind2][ind] = TBME;
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
	CCME.HHPP2[ind2][ind] = -1.0 * TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell4] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
	CCME.HHPP2[ind2][ind] = -1.0 * TBME;
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
	CCME.HHPP2[ind2][ind] = TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell1] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
	CCME.HHPP3[ind2][ind] = TBME;
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
	CCME.HHPP3[ind2][ind] = -1.0 * TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell2] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
	CCME.HHPP3[ind2][ind] = -1.0 * TBME;
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
	CCME.HHPP3[ind2][ind] = TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell3]-Space.levelsl[shell1]), Space.levelsm[shell3]-Space.levelsm[shell1], 
			 Space.levelst[shell3]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	CCME.HHPP4[ind1][ind] = TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell3]-Space.levelsl[shell2]), Space.levelsm[shell3]-Space.levelsm[shell2], 
			 Space.levelst[shell3]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	CCME.HHPP4[ind1][ind] = -1.0 * TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell4]-Space.levelsl[shell1]), Space.levelsm[shell4]-Space.levelsm[shell1], 
			 Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	CCME.HHPP4[ind1][ind] = -1.0 * TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell4]-Space.levelsl[shell2]), Space.levelsm[shell4]-Space.levelsm[shell2], 
			 Space.levelst[shell4]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	CCME.HHPP4[ind1][ind] = TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell4]-Space.levelsl[shell2]), Space.levelsm[shell4]-Space.levelsm[shell2], 
			 Space.levelst[shell4]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	CCME.HHPP4T[ind1][ind] = TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell4]-Space.levelsl[shell1]), Space.levelsm[shell4]-Space.levelsm[shell1], 
			 Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell3]-Space.levelsl[shell2]), Space.levelsm[shell3]-Space.levelsm[shell2], 
			 Space.levelst[shell3]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
	ind1 = HO_tbInd2(Space, std::pow(-1.0, Space.levelsl[shell3]-Space.levelsl[shell1]), Space.levelsm[shell3]-Space.levelsm[shell1], 
			 Space.levelst[shell3]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	CCME.HHPP4T[ind1][ind] = TBME;
      }  
    }
  }
  else if(Parameters.basis == "CART"){
    for(int i = 0; i < NumElements; ++i){
      getline(interaction, interactionline);
      std::stringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
      if(std::abs(TBME) < 1.0E-10){ continue; }
      shell1 -= 1; shell2 -= 1; shell3 -= 1; shell4 -= 1;
      TBME *= Parameters.tbstrength;
      ptype = Space.levelstype[shell1];
      qtype = Space.levelstype[shell2];
      rtype = Space.levelstype[shell3];
      stype = Space.levelstype[shell4];
      
      if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	ind1 = CART_tbInd1(Space, Space.levelsnx[shell1]+Space.levelsnx[shell2], Space.levelsny[shell1]+Space.levelsny[shell2], 
			   Space.levelsnz[shell1]+Space.levelsnz[shell2], Space.levelsm[shell1]+Space.levelsm[shell2], 
			   Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell3, shell4);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell3, shell4);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell4, shell3);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell4, shell3);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	CCME.HHHH[ind1][ind] = TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	CCME.HHHH[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	CCME.HHHH[ind1][ind] = TBME;
      }
      else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	ind1 = CART_tbInd1(Space, Space.levelsnx[shell1]+Space.levelsnx[shell2], Space.levelsny[shell1]+Space.levelsny[shell2], 
			   Space.levelsnz[shell1]+Space.levelsnz[shell2], Space.levelsm[shell1]+Space.levelsm[shell2], 
			   Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell3, shell4);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell3, shell4);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell4, shell3);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell4, shell3);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell1, shell2);
	CCME.PPPP[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell1, shell2);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell2, shell1);
	CCME.PPPP[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell2, shell1);
	CCME.PPPP[ind1][ind] = TBME;
      }
      else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell4]-Space.levelsnx[shell1], Space.levelsny[shell4]-Space.levelsny[shell1], 
			   Space.levelsnz[shell4]-Space.levelsnz[shell1], Space.levelsm[shell4]-Space.levelsm[shell1], 
			   Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
	CCME.HPHP1[ind1][ind] = TBME;
	ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
	CCME.HPHP1[ind1][ind] = TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell1]-Space.levelsnx[shell4], Space.levelsny[shell1]-Space.levelsny[shell4], 
			   Space.levelsnz[shell1]-Space.levelsnz[shell4], Space.levelsm[shell1]-Space.levelsm[shell4], 
			   Space.levelst[shell1]-Space.levelst[shell4]);
	ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell1, shell4, shell3, shell2);
	CCME.HPHP2[ind1][ind] = TBME;
	ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell3, shell2, shell1, shell4);
	CCME.HPHP2[ind1][ind] = TBME;
      }
      else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	ind1 = CART_tbInd1(Space, Space.levelsnx[shell1]+Space.levelsnx[shell2], Space.levelsny[shell1]+Space.levelsny[shell2], 
			   Space.levelsnz[shell1]+Space.levelsnz[shell2], Space.levelsm[shell1]+Space.levelsm[shell2], 
			   Space.levelst[shell1]+Space.levelst[shell2]);
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	CCME.HHPP1[ind1][ind] = TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	CCME.HHPP1[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	CCME.HHPP1[ind1][ind] = -1.0 * TBME;
	ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	CCME.HHPP1[ind1][ind] = TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell3] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
	CCME.HHPP2[ind2][ind] = TBME;
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
	CCME.HHPP2[ind2][ind] = -1.0 * TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell4] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
	CCME.HHPP2[ind2][ind] = -1.0 * TBME;
	ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
	CCME.HHPP2[ind2][ind] = TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell1] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
	CCME.HHPP3[ind2][ind] = TBME;
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
	CCME.HHPP3[ind2][ind] = -1.0 * TBME;
	for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[shell2] == j){ ind2 = j; break; }; }
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
	CCME.HHPP3[ind2][ind] = -1.0 * TBME;
	ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
	CCME.HHPP3[ind2][ind] = TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell3]-Space.levelsnx[shell1], Space.levelsny[shell3]-Space.levelsny[shell1], 
			   Space.levelsnz[shell3]-Space.levelsnz[shell1], Space.levelsm[shell3]-Space.levelsm[shell1], 
			   Space.levelst[shell3]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	CCME.HHPP4[ind1][ind] = TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell3]-Space.levelsnx[shell2], Space.levelsny[shell3]-Space.levelsny[shell2], 
			   Space.levelsnz[shell3]-Space.levelsnz[shell2], Space.levelsm[shell3]-Space.levelsm[shell2], 
			   Space.levelst[shell3]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	CCME.HHPP4[ind1][ind] = -1.0 * TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell4]-Space.levelsnx[shell1], Space.levelsny[shell4]-Space.levelsny[shell1], 
			   Space.levelsnz[shell4]-Space.levelsnz[shell1], Space.levelsm[shell4]-Space.levelsm[shell1], 
			   Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	CCME.HHPP4[ind1][ind] = -1.0 * TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell4]-Space.levelsnx[shell2], Space.levelsny[shell4]-Space.levelsny[shell2], 
			   Space.levelsnz[shell4]-Space.levelsnz[shell2], Space.levelsm[shell4]-Space.levelsm[shell2], 
			   Space.levelst[shell4]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	CCME.HHPP4[ind1][ind] = TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell4]-Space.levelsnx[shell2], Space.levelsny[shell4]-Space.levelsny[shell2], 
			   Space.levelsnz[shell4]-Space.levelsnz[shell2], Space.levelsm[shell4]-Space.levelsm[shell2], 
			   Space.levelst[shell4]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	CCME.HHPP4T[ind1][ind] = TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell4]-Space.levelsnx[shell1], Space.levelsny[shell4]-Space.levelsny[shell1], 
			   Space.levelsnz[shell4]-Space.levelsnz[shell1], Space.levelsm[shell4]-Space.levelsm[shell1], 
			   Space.levelst[shell4]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell3]-Space.levelsnx[shell2], Space.levelsny[shell3]-Space.levelsny[shell2], 
			   Space.levelsnz[shell3]-Space.levelsnz[shell2], Space.levelsm[shell3]-Space.levelsm[shell2], 
			   Space.levelst[shell3]-Space.levelst[shell2]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
	ind1 = CART_tbInd2(Space, Space.levelsnx[shell3]-Space.levelsnx[shell1], Space.levelsny[shell3]-Space.levelsny[shell1], 
			   Space.levelsnz[shell3]-Space.levelsnz[shell1], Space.levelsm[shell3]-Space.levelsm[shell1], 
			   Space.levelst[shell3]-Space.levelst[shell1]);
	ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	CCME.HHPP4T[ind1][ind] = TBME;
      }
    }
  }
  
  interaction.close();
  
  return CCME;

}



CC_Matrix_Elements Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{

  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  CC_Matrix_Elements CCME(Chan);
  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);

  for(int i = 0; i < Chan.size1; ++i)
    {
      #pragma omp parallel for
      for(int pq = 0; pq < Chan.hh[i]; ++pq)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.hhvec1[i][2*pq];
	  shell2 = Chan.hhvec1[i][2*pq + 1];
	  for(int rs = pq; rs < Chan.hh[i]; ++rs)
	    {
	      shell3 = Chan.hhvec1[i][2*rs];
	      shell4 = Chan.hhvec1[i][2*rs + 1];
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HHHH[i][Chan.hh[i]*pq + rs] = TBME;
	      CCME.HHHH[i][Chan.hh[i]*rs + pq] = TBME;
	    }
	}
      #pragma omp parallel for
      for(int pq = 0; pq < Chan.pp[i]; ++pq)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.ppvec1[i][2*pq];
	  shell2 = Chan.ppvec1[i][2*pq + 1];
	  for(int rs = pq; rs < Chan.pp[i]; ++rs)
	    {
	      shell3 = Chan.ppvec1[i][2*rs];
	      shell4 = Chan.ppvec1[i][2*rs + 1];
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.PPPP[i][Chan.pp[i]*pq + rs] = TBME;
	      CCME.PPPP[i][Chan.pp[i]*rs + pq] = TBME;
	    }
	}
      #pragma omp parallel for
      for(int pq = 0; pq < Chan.pp[i]; ++pq)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.ppvec1[i][2*pq];
	  shell2 = Chan.ppvec1[i][2*pq + 1];
	  for(int rs = 0; rs < Chan.hh[i]; ++rs)
	    {
	      shell3 = Chan.hhvec1[i][2*rs];
	      shell4 = Chan.hhvec1[i][2*rs + 1];
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HHPP1[i][Chan.hh[i]*pq + rs] = TBME;
	    }
	}
    }

  for(int i = 0; i < Chan.size3; ++i)
    {
      #pragma omp parallel for
      for(int r = 0; r < Chan.p[i]; ++r)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell3 = Chan.pvec1[i][r];
	  for(int pqs = 0; pqs < Chan.hhp[i]; ++pqs)
	    {
	      shell1 = Chan.hhpvec1[i][3*pqs];
	      shell2 = Chan.hhpvec1[i][3*pqs + 1];
	      shell4 = Chan.hhpvec1[i][3*pqs + 2];
	      if(shell3 == shell4){ continue; }
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HHPP2[i][Chan.hhp[i]*r + pqs] = TBME;
	    }
	}
    }

  for(int i = 0; i < Chan.size3; ++i)
    {
      #pragma omp parallel for
      for(int p = 0; p < Chan.h[i]; ++p)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.hvec1[i][p];
	  for(int qrs = 0; qrs < Chan.hpp[i]; ++qrs)
	    {
	      shell2 = Chan.hppvec1[i][3*qrs];
	      shell3 = Chan.hppvec1[i][3*qrs + 1];
	      shell4 = Chan.hppvec1[i][3*qrs + 2];
	      if(shell1 == shell2){ continue; }
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HHPP3[i][Chan.hpp[i]*p + qrs] = TBME;
	    }
	}
    }

  for(int i = 0; i < Chan.size2; ++i)
    {
      #pragma omp parallel for
      for(int ps = 0; ps < Chan.hp2[i]; ++ps)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.hp2vec1[i][2*ps];
	  shell4 = Chan.hp2vec1[i][2*ps + 1];
	  for(int qr = ps; qr < Chan.hp2[i]; ++qr)
	    {
	      shell3 = Chan.hp2vec1[i][2*qr];
	      shell2 = Chan.hp2vec1[i][2*qr + 1];
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HPHP1[i][Chan.hp2[i]*ps + qr] = TBME;
	      CCME.HPHP1[i][Chan.hp2[i]*qr + ps] = TBME;
	    }
	}
      #pragma omp parallel for
      for(int ps = 0; ps < Chan.hp2[i]; ++ps)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell1 = Chan.hp2vec1[i][2*ps];
	  shell4 = Chan.hp2vec1[i][2*ps + 1];
	  for(int qs = 0; qs < Chan.hp1[i]; ++qs)
	    {
	      shell2 = Chan.hp1vec1[i][2*qs];
	      shell3 = Chan.hp1vec1[i][2*qs + 1];
	      if(shell1 == shell2 || shell3 == shell4){ continue; }
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell4, shell3, L);
	      CCME.HHPP4[i][Chan.hp1[i]*ps + qs] = TBME;
	      TBME = vint_Minnesota_Momentum(Space, shell2, shell1, shell3, shell4, L);
	      CCME.HHPP4T[i][Chan.hp1[i]*ps + qs] = TBME;
	    }
	}
      #pragma omp parallel for
      for(int qr = 0; qr < Chan.hp1[i]; ++qr)
	{
	  double TBME;
	  int shell1, shell2, shell3, shell4;
	  shell3 = Chan.hp1vec1[i][2*qr];
	  shell2 = Chan.hp1vec1[i][2*qr + 1];
	  for(int ps = qr; ps < Chan.hp1[i]; ++ps)
	    {
	      shell1 = Chan.hp1vec1[i][2*ps];
	      shell4 = Chan.hp1vec1[i][2*ps + 1];
	      TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	      CCME.HPHP2[i][Chan.hp1[i]*qr + ps] = TBME;
	      CCME.HPHP2[i][Chan.hp1[i]*ps + qr] = TBME;
	    }
	}
    }

  return CCME;

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
  double relMomTransfX1, relMomTransfY1, relMomTransfZ1, relMomTransfX2, relMomTransfY2, relMomTransfZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  relMomTransfX1 = (M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[qj] - Space.levelsnx[qk] + Space.levelsnx[ql]);
  relMomTransfY1 = (M_PI/L) * (Space.levelsny[qi] - Space.levelsny[qj] - Space.levelsny[qk] + Space.levelsny[ql]);
  relMomTransfZ1 = (M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[qj] - Space.levelsnz[qk] + Space.levelsnz[ql]);

  relMomTransfX2 = (M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[qj] - Space.levelsnx[ql] + Space.levelsnx[qk]);
  relMomTransfY2 = (M_PI/L) * (Space.levelsny[qi] - Space.levelsny[qj] - Space.levelsny[ql] + Space.levelsny[qk]);
  relMomTransfZ2 = (M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[qj] - Space.levelsnz[ql] + Space.levelsnz[qk]);

  qSquared1 = relMomTransfX1 * relMomTransfX1 + relMomTransfY1 * relMomTransfY1 + relMomTransfZ1 * relMomTransfZ1;
  qSquared2 = relMomTransfX2 * relMomTransfX2 + relMomTransfY2 * relMomTransfY2 + relMomTransfZ2 * relMomTransfZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeMtxEle(Space.levelsm[qi], Space.levelsm[qj], Space.levelsm[qk], Space.levelsm[ql]);
  isoSpinEx1 = spinExchangeMtxEle(Space.levelst[qi], Space.levelst[qj], Space.levelst[qk], Space.levelst[ql]);

  spinEx2 = spinExchangeMtxEle(Space.levelsm[qi], Space.levelsm[qj], Space.levelsm[ql], Space.levelsm[qk]);
  isoSpinEx2 = spinExchangeMtxEle(Space.levelst[qi], Space.levelst[qj], Space.levelst[ql], Space.levelst[qk]);
  
  IsIt1 = kron_del(Space.levelsm[qi], Space.levelsm[qk]) * kron_del(Space.levelsm[qj], Space.levelsm[ql]) * 
    kron_del(Space.levelst[qi], Space.levelst[qk]) * kron_del(Space.levelst[qj], Space.levelst[ql]);
  PsIt1 = spinEx1 * kron_del(Space.levelst[qi], Space.levelst[qk]) * kron_del(Space.levelst[qj], Space.levelst[ql]);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(Space.levelsm[qi], Space.levelsm[qk])*kron_del(Space.levelsm[qj], Space.levelsm[ql]) * isoSpinEx1;

  IsIt2 = kron_del(Space.levelsm[qi], Space.levelsm[ql]) * kron_del(Space.levelsm[qj], Space.levelsm[qk]) * 
    kron_del(Space.levelst[qi], Space.levelst[ql]) * kron_del(Space.levelst[qj], Space.levelst[qk]);
  PsIt2 = spinEx2 * kron_del(Space.levelst[qi], Space.levelst[ql]) * kron_del(Space.levelst[qj], Space.levelst[qk]);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(Space.levelsm[qi], Space.levelsm[ql]) * kron_del(Space.levelsm[qj], Space.levelsm[qk]) * isoSpinEx2;

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


CCD Perform_CCD(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, const Channels &Chan)
{
  std::cout << "Performing CCD ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  double error = 1000.0;
  int ind = 0;

  CCD CCin(Chan, Parameters, Space);
  CCD CCout = CCin;

  /*std::vector<std::vector<std::vector<double> > > p(2);
  std::vector<std::vector<std::vector<double> > > delp(1);
  std::vector<double> B;
  std::vector<double> B2;
  std::vector<int> ipiv;
  std::vector<double> work;
  int lwork;
  int info = 0;
  int rem = 0;
  int maxl = 8;
  int P = 1;
  int N = 2;
  double max;

  p[0] = CCin.T1;
  p[1] = CCin.T1;
  delp[0] = CCin.T1;
  for(int chan = 0; chan < Chan.size1; ++chan){
    int ind1, ind2, ind3, ind4;
    double tempen, tempt;
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	p[1][chan][hind * Chan.pp[chan] + pind] = tempt;
	delp[0][chan][hind * Chan.pp[chan] + pind] = tempt;
      }
    }
  }
  CCD CCout = CCin;

  Doubles_Step(Space, Chan, CCME, CCin, CCout);

  p.push_back(CCin.T1);
  delp.push_back(CCin.T1);
  P = int(delp.size());
  N = int(delp.size() + 1);
  B.assign(N * N, 0.0);

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	double tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	p[2][chan][hind * Chan.pp[chan] + pind] = tempt;
	delp[1][chan][hind * Chan.pp[chan] + pind] = tempt - p[1][chan][hind * Chan.pp[chan] + pind]; //CCin.T1[chan][hind * Chan.pp[chan] + pind];
	for(int i = 0; i < P; ++i){
	  for(int j = 0; j < P; ++j){
	    B[N * i + j] += delp[i][chan][hind * Chan.pp[chan] + pind] * delp[j][chan][hind * Chan.pp[chan] + pind];
	  }
	}
      }
    }
  }
  for(int i = 0; i < P; ++i){ B[N * i + P] = -1.0; B[N * P + i] = -1.0; }
  B[N * P + P] = 0.0;

  lwork = N;
  ipiv.resize(N);
  work.resize(sizeof(double) * N);
  B2 = B;
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < N; ++j){
      std::cout << B2[N*i + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  dgetrf_(&N, &N, & *B2.begin(), &N, & *ipiv.begin(), &info);
  dgetri_(&N, & *B2.begin(), &N, & *ipiv.begin(), & *work.begin(), &lwork, &info);
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < N; ++j){
      std::cout << B2[N*i + j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  error = B2[N * P + P];

  double tempt;
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = 0.0;
	for(int i = 0; i < P; ++i){
	  tempt -= B2[N * i + P] * p[i][chan][hind * Chan.pp[chan] + pind];
	}
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
      }
    }
  }

  while((error > 10e-6 && ind < 1000) || ind < 10)
    {
      Doubles_Step(Space, Chan, CCME, CCin, CCout);
      CCout.CCDE = 0.0;
      error = 0.0;
      P = int(delp.size());
      N = int(delp.size() + 1);
      rem = 0;
      max = B[0];
      for(int i = 1; i < P; ++i){ if(B[N * i + i] >= max){ max = B[N * i + i]; rem = i; }; }
      if(P != maxl){
	p.push_back(CCin.T1);
	delp.push_back(CCin.T1);
	P = int(delp.size());
	N = int(delp.size() + 1);
	B.resize(N * N);
	for(int i = P - 1; i >= 0; --i){
	  for(int j = P - 1; j >= 0; --j){
	    B[N * i + j] = B[P * i + j];
	  }
	}
      }
      else if(rem == P - 1){ 
	++maxl; 
	p.push_back(CCin.T1);
	delp.push_back(CCin.T1);
	P = int(delp.size());
	N = int(delp.size() + 1);
	B.resize(N * N);
	for(int i = P - 1; i >= 0; --i){
	  for(int j = P - 1; j >= 0; --j){
	    B[N * i + j] = B[P * i + j];
	  }
	}
      }
      else{
	for(int i = 0; i < P - 1; ++i){ 
	  if(i >= rem){ p[i] = p[i+1]; delp[i] = delp[i+1]; }
	  for(int j = 0; j < P - 1; ++j){
	    if(i >= rem || j >= rem){ B[N * i + j] = B[N * (i + 1) + (j + 1)]; }
	  }
	}
      }

      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    double tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	    tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	    p[N - 1][chan][hind * Chan.pp[chan] + pind] = tempt;
	    delp[P - 1][chan][hind * Chan.pp[chan] + pind] = tempt - p[P - 1][chan][hind * Chan.pp[chan] + pind]; //CCin.T1[chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      B[N * i + (P - 1)] += delp[i][chan][hind * Chan.pp[chan] + pind] * delp[P - 1][chan][hind * Chan.pp[chan] + pind];
	      if(i == P - 1){ break; } // don't double count
	      B[N * (P - 1) + i] += delp[P - 1][chan][hind * Chan.pp[chan] + pind] * delp[i][chan][hind * Chan.pp[chan] + pind];
	    }
	  }
	}
      }
      for(int i = 0; i < P; ++i){ B[N * i + P] = -1.0; B[N * P + i] = -1.0; }
      B[N * P + P] = 0.0;

      lwork = N;
      ipiv.resize(N);
      work.resize(sizeof(double) * N);
      B2 = B;
      for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){
	  std::cout << B2[N*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      dgetrf_(&N, &N, & *B2.begin(), &N, & *ipiv.begin(), &info);
      dgetri_(&N, & *B2.begin(), &N, & *ipiv.begin(), & *work.begin(), &lwork, &info);
      for(int i = 0; i < N; ++i){
	for(int j = 0; j < N; ++j){
	  std::cout << B2[N*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;

      error = B[N * (P - 1) + (P - 1)];
      
      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = 0.0;
	    for(int i = 0; i < P; ++i){
	      tempt += B2[N * i + P] * p[i][chan][hind * Chan.pp[chan] + pind];
	    }
	    CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  }
	}
      }
      std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
      ++ind;
    }
    std::cout << std::endl << std::endl;*/

  ////////////////////////////////////////

  /*std::vector<std::vector<std::vector<double> > > Vin(3); // 0 and V/E
  std::vector<std::vector<std::vector<double> > > F(2); // Vout - Vin
  std::vector<std::vector<std::vector<double> > > delV(1);
  std::vector<std::vector<std::vector<double> > > delF(1);
  std::vector<std::vector<std::vector<double> > > u(1);
  std::vector<double> a(1, 0.0);
  std::vector<double> B(1, 0.0);
  std::vector<int> ipiv;
  std::vector<double> work;
  double tempt, tempen, norm;
  int lwork;
  int info = 0;
  int maxl = 12;
  int P = 1;
  int N = 2;

  double alpha = 0.05;
  double w0 = 0.01;

  Vin[0] = CCin.T1;
  Vin[1] = CCin.T1;
  F[0] = CCin.T1;
  F[1] = CCin.T1;
  delV[0] = CCin.T1;
  delF[0] = CCin.T1;
  u[0] = CCin.T1;
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
      }
    }
  }
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
  norm = std::sqrt(norm);
  
  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	delV[0][chan][hind * Chan.pp[chan] + pind] = (Vin[1][chan][hind * Chan.pp[chan] + pind] - Vin[0][chan][hind * Chan.pp[chan] + pind])/norm;
	delF[0][chan][hind * Chan.pp[chan] + pind] = (F[1][chan][hind * Chan.pp[chan] + pind] - F[0][chan][hind * Chan.pp[chan] + pind])/norm;
	u[0][chan][hind * Chan.pp[chan] + pind] = alpha * delF[0][chan][hind * Chan.pp[chan] + pind] + delV[0][chan][hind * Chan.pp[chan] + pind];
	a[0] += delF[0][chan][hind * Chan.pp[chan] + pind] * delF[0][chan][hind * Chan.pp[chan] + pind];
      }
    }
  }
  a[0] += w0*w0;

  for(int chan = 0; chan < Chan.size1; ++chan){
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	tempt = Vin[1][chan][hind * Chan.pp[chan] + pind] + alpha * F[1][chan][hind * Chan.pp[chan] + pind];
	tempt -= delF[0][chan][hind * Chan.pp[chan] + pind] * F[1][chan][hind * Chan.pp[chan] + pind] * u[0][chan][hind * Chan.pp[chan] + pind] / a[0];
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
      }
    }
  }
  Vin[2] = CCin.T1;

  while((error > 10e-6 && ind < 1000) || ind < 10)
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
	Vin[Vin.size() - 1] = CCout.T1;
	F[F.size() - 1] = CCout.T1;
	delV[delV.size() - 1] = CCout.T1;
	delF[delF.size() - 1] = CCout.T1;
	u[u.size() - 1] = CCout.T1;
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
	for(int i = 0; i < (N - 1); ++i){ F[i] = F[i+1]; }
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
      norm = std::sqrt(norm);
      error = std::sqrt(error);

      std::cout << " %% " << norm << " " << error << std::endl;
      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    delV[P - 1][chan][hind * Chan.pp[chan] + pind] = (Vin[N - 1][chan][hind * Chan.pp[chan] + pind] - Vin[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    delF[P - 1][chan][hind * Chan.pp[chan] + pind] = (F[N - 1][chan][hind * Chan.pp[chan] + pind] - F[N - 2][chan][hind * Chan.pp[chan] + pind])/norm;
	    u[P - 1][chan][hind * Chan.pp[chan] + pind] = alpha * delF[P - 1][chan][hind * Chan.pp[chan] + pind] + delV[P - 1][chan][hind * Chan.pp[chan] + pind];
	    for(int i = 0; i < P; ++i){
	      a[P * i + (P - 1)] += delF[i][chan][hind * Chan.pp[chan] + pind] * delF[P - 1][chan][hind * Chan.pp[chan] + pind];
	      if(i == P - 1){ break; } // don't double count
	      a[P * (P - 1) + i] += delF[P - 1][chan][hind * Chan.pp[chan] + pind] * delF[i][chan][hind * Chan.pp[chan] + pind];
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
		tempt -= delF[i][chan][hind * Chan.pp[chan] + pind] * F[N - 1][chan][hind * Chan.pp[chan] + pind] * 
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

  double tempen, tempt;
  for(int chan = 0; chan < Chan.size1; ++chan){
    int ind1, ind2, ind3, ind4;
    for(int hind = 0; hind < Chan.hh[chan]; ++hind){
      for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	ind1 = Chan.hhvec1[chan][2*hind];
	ind2 = Chan.hhvec1[chan][2*hind + 1];
	ind3 = Chan.ppvec1[chan][2*pind];
	ind4 = Chan.ppvec1[chan][2*pind + 1];
	tempen = Space.levelsen[ind1] + Space.levelsen[ind2] - Space.levelsen[ind3] - Space.levelsen[ind4];
	tempt = CCME.HHPP1[chan][pind * Chan.hh[chan] + hind] / tempen;
	CCin.Evec[chan][hind * Chan.pp[chan] + pind] = tempen;
	CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
      }
    }
  }
  CCout.Evec = CCin.Evec;
  
  while((error > 10e-6 && ind < 1000) || ind < 10)
    {
      Doubles_Step(Space, Chan, CCME, CCin, CCout);

      CCout.CCDE = 0.0;
      for(int chan = 0; chan < Chan.size1; ++chan){
	for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	  for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	    tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	    tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	    tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	    CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	    CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  }
	}
      }
      std::cout << CCin.CCDE << " " << CCout.CCDE << std::endl;
      error = std::abs((CCout.CCDE - CCin.CCDE)/CCout.CCDE);
      CCin.CCDE = CCout.CCDE;
      
      std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
      ++ind;
    }
  std::cout << std::endl << std::endl;
    
  return CCout;

}


void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME)
{
  int ind, ind1;
  double M, T, P;
  int Nx, Ny, Nz;
  if(Parameters.basis == "HO")
    {
      for(int i = 0; i < Space.indtot; ++i)
	{
	  if(Space.levelstype[i] == "hole")
	    {
	      for(int j = 0; j < Space.indtot; ++j)
		{
		  if(Space.levelstype[j] != "hole" || i == j){ continue; }
		  M = Space.levelsm[i] + Space.levelsm[j];
		  T = Space.levelst[i] + Space.levelst[j];
		  P = std::pow(-1.0, Space.levelsl[i] + Space.levelsl[j]);
		  ind1 = HO_tbInd1(Space, P, M, T);
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
		  Space.levelsen[i] += ME.HHHH[ind1][ind];
		}
	    }
	  else
	    {
	      for(int j = 0; j < Space.indtot; ++j)
		{
		  if(Space.levelstype[j] != "hole"){ continue; }
		  M = Space.levelsm[j] - Space.levelsm[i];
		  T = Space.levelst[j] - Space.levelst[i];
		  P = std::pow(-1.0, Space.levelsl[j] - Space.levelsl[i]);
		  ind1 = HO_tbInd2(Space, P, M, T);
		  ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], j, i, j, i);
		  Space.levelsen[i] += ME.HPHP2[ind1][ind];
		}	
	    }
	}
    }
  else if(Parameters.basis == "CART")
    {
      for(int i = 0; i < Space.indtot; ++i)
	{
	  if(Space.levelstype[i] == "hole")
	    {
	      for(int j = 0; j < Space.indtot; ++j)
		{
		  if(Space.levelstype[j] != "hole" || i == j){ continue; }
		  M = Space.levelsm[i] + Space.levelsm[j];
		  T = Space.levelst[i] + Space.levelst[j];
		  Nx = Space.levelsnx[i] + Space.levelsnx[j];
		  Ny = Space.levelsny[i] + Space.levelsny[j];
		  Nz = Space.levelsnz[i] + Space.levelsnz[j];
		  ind1 = CART_tbInd1(Space, Nx, Ny, Nz, M, T);
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
		  Space.levelsen[i] += ME.HHHH[ind1][ind];
		}
	    }
	  else
	    {
	      for(int j = 0; j < Space.indtot; ++j)
		{
		  if(Space.levelstype[j] != "hole"){ continue; }
		  M = Space.levelsm[j] - Space.levelsm[i];
		  T = Space.levelst[j] - Space.levelst[i];
		  Nx = Space.levelsnx[j] - Space.levelsnx[i];
		  Ny = Space.levelsny[j] - Space.levelsny[i];
		  Nz = Space.levelsnz[j] - Space.levelsnz[i];
		  ind1 = CART_tbInd2(Space, Nx, Ny, Nz, M, T);
		  ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], j, i, j, i);
		  Space.levelsen[i] += ME.HPHP2[ind1][ind];
		}	
	    }
	}
    }
}

double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME)
{
  double energy = 0;
  int ind, ind1;
  double M, T, P;
  int Nx, Ny, Nz;
  if(Parameters.basis == "HO")
    {
      for(int i = 0; i < Space.indtot; ++i)
	{
	  if(Space.levelstype[i] != "hole"){ continue; }
	  energy += Space.levelsen[i];
	  for(int j = 0; j < Space.indtot; ++j)
	    {
	      if(Space.levelstype[j] != "hole" || i == j){ continue; }
	      M = Space.levelsm[i] + Space.levelsm[j];
	      T = Space.levelst[i] + Space.levelst[j];
	      P = std::pow(-1.0, Space.levelsl[i] + Space.levelsl[j]);
	      ind1 = HO_tbInd1(Space, P, M, T);
	      ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
	      energy -= 0.5 * ME.HHHH[ind1][ind];
	    }
	}
    }
  else if(Parameters.basis == "CART")
    {
      for(int i = 0; i < Space.indtot; ++i)
	{
	  if(Space.levelstype[i] != "hole"){ continue; }
	  energy += Space.levelsen[i];
	  for(int j = 0; j < Space.indtot; ++j)
	    {
	      if(Space.levelstype[j] != "hole" || i == j){ continue; }
	      M = Space.levelsm[i] + Space.levelsm[j];
	      T = Space.levelst[i] + Space.levelst[j];
	      Nx = Space.levelsnx[i] + Space.levelsnx[j];
	      Ny = Space.levelsny[i] + Space.levelsny[j];
	      Nz = Space.levelsnz[i] + Space.levelsnz[j];
	      ind1 = CART_tbInd1(Space, Nx, Ny, Nz, M, T);
	      ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
	      energy -= 0.5 * ME.HHHH[ind1][ind];
	    }
	}
    }
  return energy;
}

void Doubles_Step(const Model_Space &Space, const Channels &Chan, CC_Matrix_Elements &ME, CCD &CC, CCD &CC2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = 0.25, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';

  for(int chan = 0; chan < Chan.size2; ++chan)
    {
      int hp1 = Chan.hp1[chan];
      int hp2 = Chan.hp2[chan];
      if(hp1 == 0 || hp2 == 0){ continue; }
      //S4 = T2*HHPP4
      RM_dgemm(&N,&N,&hp1,&hp1,&hp2,&fac1,& *CC.T2[chan].begin(),&hp1,& *ME.HHPP4[chan].begin(),&hp2,&fac2,& *CC2.S4[chan].begin(),&hp1);
      //S4T = T2*HHPP4T
      RM_dgemm(&N,&N,&hp1,&hp1,&hp2,&fac1,& *CC.T2[chan].begin(),&hp1,& *ME.HHPP4T[chan].begin(),&hp2,&fac2,& *CC2.S4T[chan].begin(),&hp1);
      
      //T2 & T3
      RM_dgemm(&N,&N,&hp1,&hp2,&hp2,&fac1,& *CC.T2[chan].begin(),&hp1,& *ME.HPHP1[chan].begin(),&hp2,&fac2,& *CC2.T2[chan].begin(),&hp1); 
      RM_dgemm(&N,&N,&hp1,&hp2,&hp1,&fac1,& *ME.HPHP2[chan].begin(),&hp1,& *CC.T2T[chan].begin(),&hp1,&fac1,& *CC2.T2[chan].begin(),&hp1);
      CC2.T3[chan] = CC2.T2[chan]; // T3 = T2
      RM_dgemm(&N,&N,&hp1,&hp2,&hp1,&fac1,& *CC2.S4[chan].begin(),&hp1,& *CC.T2T[chan].begin(),&hp1,&fac5,& *CC2.T2[chan].begin(),&hp1);
      RM_dgemm(&N,&N,&hp1,&hp2,&hp1,&fac5,& *CC2.S4T[chan].begin(),&hp1,& *CC.T2T[chan].begin(),&hp1,&fac1,& *CC2.T3[chan].begin(),&hp1);
    }

  for(int chan = 0; chan < Chan.size1; ++chan)
    {
      int hh = Chan.hh[chan];
      int pp = Chan.pp[chan];
      if(hh == 0 || pp == 0){ continue; }
      //S1 = T1*HHPP1
      RM_dgemm(&N,&N,&hh,&hh,&pp,&fac1,& *CC.T1[chan].begin(),&hh,& *ME.HHPP1[chan].begin(),&pp,&fac2,& *CC2.S1[chan].begin(),&hh);
      
      //T1 = 0.5*T1*PPPP + 0.5*HHHH*T1 + 0.25*S1*T1
      RM_dgemm(&N,&N,&hh,&pp,&pp,&fac3,& *CC.T1[chan].begin(),&hh,& *ME.PPPP[chan].begin(),&pp,&fac2,& *CC2.T1[chan].begin(),&hh);
      RM_dgemm(&N,&N,&hh,&pp,&hh,&fac3,& *ME.HHHH[chan].begin(),&hh,& *CC.T1[chan].begin(),&hh,&fac1,& *CC2.T1[chan].begin(),&hh);
      RM_dgemm(&N,&N,&hh,&pp,&hh,&fac4,& *CC2.S1[chan].begin(),&hh,& *CC.T1[chan].begin(),&hh,&fac1,& *CC2.T1[chan].begin(),&hh);
    }
  
  for(int chan = 0; chan < Chan.size3; ++chan)
    {
      int p = Chan.p[chan];
      int hhp = Chan.hhp[chan];
      if(p == 0 || hhp == 0){ continue; }
      //S2 = HHPP2*T5
      RM_dgemm(&N,&N,&p,&p,&hhp,&fac1,& *ME.HHPP2[chan].begin(),&p,& *CC.T5[chan].begin(),&hhp,&fac2,& *CC2.S2[chan].begin(),&p);
      //T4 = -0.5*T4*S2  &  T5 
      RM_dgemm(&N,&N,&hhp,&p,&p,&fac6,& *CC.T4[chan].begin(),&hhp,& *CC2.S2[chan].begin(),&p,&fac2,& *CC2.T4[chan].begin(),&hhp);
    }
  
  for(int chan = 0; chan < Chan.size3; ++chan)
    {
      int h = Chan.h[chan];
      int hpp = Chan.hpp[chan];
      if(h == 0 || hpp == 0){ continue; }
      //S3 = HHPP3*T7
      RM_dgemm(&N,&N,&h,&h,&hpp,&fac1,& *ME.HHPP3[chan].begin(),&h,& *CC.T7[chan].begin(),&hpp,&fac2,& *CC2.S3[chan].begin(),&h);
      //T6 = -0.5*T6*S3  &  T7
      RM_dgemm(&N,&N,&hpp,&h,&h,&fac6,& *CC.T6[chan].begin(),&hpp,& *CC2.S3[chan].begin(),&h,&fac2,& *CC2.T6[chan].begin(),&hpp);
    }
}


CC_Eff Build_CC_Eff(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, CCD &CC, const Channels &Chan)
{
  CC_Eff V_Eff(Chan);
  double fac1 = 1.0, fac2 = 0.0;
  char N = 'N';

  for(int chan = 0; chan < Chan.size3; ++chan)
    {
      int p = Chan.p[chan];
      int hhp = Chan.hhp[chan];
      if(p == 0 || hhp == 0){ continue; }
      RM_dgemm(&N,&N,&p,&p,&hhp,&fac1,& *CCME.HHPP2[chan].begin(),&p,& *CC.T5[chan].begin(),&hhp,&fac2,& *V_Eff.V1[chan].begin(),&p);
      for(int i = 0; i < p; ++i){
	V_Eff.V1[chan][p * i + i] += Space.levelsen[Chan.pvec1[chan][i]];
      }
    }

  for(int chan = 0; chan < Chan.size3; ++chan)
    {
      int h = Chan.h[chan];
      int hpp = Chan.hpp[chan];
      if(h == 0 || hpp == 0){ continue; }
      RM_dgemm(&N,&N,&h,&h,&hpp,&fac1,& *CCME.HHPP3[chan].begin(),&h,& *CC.T7[chan].begin(),&hpp,&fac2,& *V_Eff.V2[chan].begin(),&h);
      for(int i = 0; i < h; ++i){
	V_Eff.V2[chan][h * i + i] += Space.levelsen[Chan.hvec1[chan][i]];
      }
    }

  for(int chan = 0; chan < Chan.size1; ++chan)
    {
      int hh = Chan.hh[chan];
      int pp = Chan.pp[chan];
      if(hh == 0 || pp == 0){ continue; }
      V_Eff.V3[chan] = CCME.PPPP[chan];
      V_Eff.V4[chan] = CCME.HHHH[chan];
      RM_dgemm(&N,&N,&pp,&pp,&hh,&fac1,& *CCME.HHPP1[chan].begin(),&pp,& *CC.T1[chan].begin(),&hh,&fac1,& *V_Eff.V3[chan].begin(),&pp);
      RM_dgemm(&N,&N,&hh,&hh,&pp,&fac1,& *CC.T1[chan].begin(),&hh,& *CCME.HHPP1[chan].begin(),&pp,&fac1,& *V_Eff.V4[chan].begin(),&hh);
    }
  
  for(int chan = 0; chan < Chan.size2; ++chan)
    {
      int hp1 = Chan.hp1[chan];
      int hp2 = Chan.hp2[chan];
      V_Eff.V5[chan] = CCME.HPHP2[chan];
      if(hp1 == 0 || hp2 == 0){ continue; }
      RM_dgemm(&N,&N,&hp1,&hp1,&hp2,&fac1,& *CC.T3[chan].begin(),&hp1,& *CCME.HHPP4[chan].begin(),&hp2,&fac1,& *V_Eff.V5[chan].begin(),&hp1);
    }

  return V_Eff;

}



void EE_EOM(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan)
{
  int Nx = Parameters.Nx;
  int Ny = Parameters.Ny;
  int Nz = Parameters.Nz;
  double M = Parameters.M;
  double T = Parameters.T;
  int ind1 = CART_tbInd1(Space, Nx, Ny, Nz, M, T);
  int count = Chan.hh[ind1] * Chan.pp[ind1];
  double memory = 24.0 + 8.0 * count * count;
  std::cout <<"Matrix for Excited State: "<< Nx <<" "<< Ny <<" "<< Nz <<" "<< M <<" "<< T <<": "<< count <<" = "<< memory/1000000.0 <<" MB"<< std::endl;

  //std::cout << std::endl << factorial(5) << " " << factorial(6.0) << " " << CGC(1.0,0.0,1.5,-0.5,0.5,-0.5) << " " << CGC(0.5,0.5,0.5,-0.5,1.0,0.0) << std::endl;
  //std::cout << std::endl << factorial2(9) << " " << factorial2(10.0) << " " << Legendre(0.25, 4, 3) << " " << Legendre(0.75, 3, -2) << std::endl;

  std::cout << SphericalY(0.5, 0.5, 4, 3) << std::endl;

}


std::vector<int> bitconfigsetup(const std::vector<int> &newconfigs, const int &N)
{
  int blength1 = int(newconfigs.size()) / N;
  std::vector<int> bconfigs;
  int btemp, bocc, shift, bsd;
  
  for(int i = 0; i < blength1; i++){
    bsd = 0;
    for(int j = 0; j < N; j++){
      bocc = newconfigs[N*i + j]; shift = bocc - 1; btemp = 1 << shift; bsd = bsd + btemp;
    }
    bconfigs.push_back(bsd);
  }
  return bconfigs;
}


double matrixe(const int &testflag, const double &strength1, const double &strength2, const std::vector<int> &indvec, const unsigned long long &bra, const unsigned long long &ket, const std::vector<double> &onebody, const std::vector<std::vector<std::vector<std::vector<int> > > > &twobodybraket, const std::vector<std::vector<std::vector<double> > > &twobody)
{
  int size = int(indvec.size());
  double temp = 0.0;
  unsigned long long one = 1, zero = 0;
  std::vector<double> tempvec(twobodybraket.size());
  
  if(bra == ket){
    for(int i = 0; i < size; ++i){
      if(((one << i) & ket) != 0){ temp = temp + onebody[indvec[i] - 1]; }
    }
  }
  
  #pragma omp parallel for
  for(int i = 0; i < int(twobodybraket.size()); ++i){
    int m, n, l, k;
    unsigned long long mbit, nbit, lbit, kbit;
    double temp2, phase;
    unsigned long long tempket;
    int flag;
    int bcount;
    unsigned long long comp;
    temp2 = 0.0;
    for(int j = 0; j < int(twobodybraket[i].size()); ++j){
      for(int q = 0; q < int(twobodybraket[i][j].size()); ++q){
	m = i + 1;
	n = i + j + 2;
	l = twobodybraket[i][j][q][0];
	k = twobodybraket[i][j][q][1];
	mbit = one << (m - 1);
	nbit = one << (n - 1);
	lbit = one << (l - 1);
	kbit = one << (k - 1);
	phase = 1.0;
	tempket = ket;
	flag = 0;
	if( (lbit & ket) != 0 && (kbit & ket) != 0 && (nbit & bra) != 0 && (mbit & bra) != 0 &&
	    ((mbit & ket) == 0 || (mbit == lbit || mbit == kbit)) && ((nbit & ket) == 0 || (nbit == lbit || nbit == kbit)) )
	  { flag = 1; }
	else if( (mbit & ket) != 0 && (nbit & ket) != 0 && (kbit & bra) != 0 && (lbit & bra) != 0 &&
		 ((lbit & ket) == 0 || (lbit == mbit || lbit == nbit)) && ((kbit & ket) == 0 || (kbit == mbit || kbit == nbit)) )
	  { flag = 1; std::swap(mbit,lbit); std::swap(m,l); std::swap(nbit,kbit); std::swap(n,k); }
	if(flag == 1){
	  comp = tempket & ~(~zero << (l - 1));
	  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
	  phase = phase*pow(-1, bcount); tempket = tempket^lbit;
	  comp = tempket & ~(~zero << (k - 1));
	  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
	  phase = phase*pow(-1, bcount); tempket = tempket^kbit;
	  comp = tempket & ~(~zero << (n - 1));
	  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
	  phase = phase*pow(-1, bcount); tempket = tempket^nbit;
	  comp = tempket & ~(~zero << (m - 1));
	  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
	  phase = phase*pow(-1, bcount); tempket = tempket^mbit;
	  if(tempket == bra){ temp2 += phase*strength2*twobody[i][j][q]; }
	}
      }
    }
    tempvec[i] = temp2;
  }
  temp = std::accumulate(tempvec.begin(),tempvec.end(),temp);
  return temp;
}
