#include "CCfunctions.hpp"
#include "MATHfunctions.hpp"

//   Function to search for index1 of pq in vec1 (size num1) and index2 of rs in vec2 (size num2)
//   Returns index of pqrs in matrix <pq|rs>. Index = index1*num2 + index2.
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

//   Function to search for index1 of p in vec1 (size num1) and index2 of qrs in vec2 (size num2)
//   Returns index of pqrs in matrix <p|qrs>. Index = index1*num2 + index2.
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
  CCDE = 0.0;
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

  #pragma omp parallel for
  for(int i = 0; i < Chan.size1; ++i){
    Tmap[i].resize(Chan.hh[i] * Chan.pp[i] * 15);
    Evec[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    T1[i].assign(Chan.hh[i] * Chan.pp[i], 0.0);
    S1[i].assign(Chan.hh[i] * Chan.hh[i], 0.0);
    for(int j = 0; j < (Chan.hh[i] * Chan.pp[i]); ++j){ Tmap[i][15*j] = j; }
  }

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
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, hind1, hind2);
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
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind1, pind2);
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
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind1, pind2);
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
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind1, pind2);
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
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind1, pind2);
	ind = Index(Chan.hhvec1[ind1], Chan.ppvec1[ind1], Chan.hh[ind1], Chan.pp[ind1], hind1, hind2, pind1, pind2);
	Tmap[ind1][15*ind + 5] = i;
	Tmap[ind1][15*ind + 6] = count;
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
Channels Setup_Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << "Building Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  Channels Chan;
  std::vector<std::vector<int> > OB_nums; // Ms,Ts,Ps for finite; Ms,Ts,Nxs,Nys,Nzs for infinite
  std::vector<int> nums; // subvector for OB_nums. size 3 for finite, size 5 for infinite.
  std::vector<std::vector<int> > Hvec, Pvec;

  if(Parameters.basis == "infinite"){ nums.resize(5); }
  else{ nums.resize(3); }

  // place first state by quantum numbers
  Chan.size3 = 1;
  Hvec.resize(Chan.size3);
  Pvec.resize(Chan.size3);
  nums[0] = Space.levelsm[0];
  nums[1] = Space.levelst[0];
  if(Parameters.basis == "infinite"){
    nums[2] = Space.levelsnx[0];
    nums[3] = Space.levelsny[0];
    nums[4] = Space.levelsnz[0];
  }
  else{
    nums[2] = -2*(Space.levelsl[0]%2) + 1;
  }
  OB_nums.push_back(nums);
  Chan.indvec.push_back(Chan.size3 - 1);
  if(Space.levelstype[0] == "hole"){ Hvec[Chan.size3 - 1].push_back(0); }
  else{ Pvec[Chan.size3 - 1].push_back(0); }

  // place the rest of the states by their quantum numbers
  for(int i = 1; i < Space.indtot; ++i){
    nums[0] = Space.levelsm[i];
    nums[1] = Space.levelst[i];
    if(Parameters.basis == "infinite"){
      nums[2] = Space.levelsnx[i];
      nums[3] = Space.levelsny[i];
      nums[4] = Space.levelsnz[i];
    }
    else{
      nums[2] = -2*(Space.levelsl[i]%2) + 1;
    }
    for(int k = 0; k < Chan.size3; ++k){
      for(int q = 0; q < int(nums.size()); ++q){
	if(nums[q] != OB_nums[k][q]){ break; }
	if(q == int(nums.size())-1){
	  Chan.indvec.push_back(k);
	  if(Space.levelstype[i] == "hole"){ Hvec[k].push_back(i); }
	  else{ Pvec[k].push_back(i); }
	  break;
	}
      }
      if(k == Chan.size3 - 1){
	++Chan.size3;
	Hvec.resize(Chan.size3);
	Pvec.resize(Chan.size3);
	OB_nums.push_back(nums);
	Chan.indvec.push_back(Chan.size3 - 1);
	if(Space.levelstype[i] == "hole"){ Hvec[Chan.size3 - 1].push_back(i); }
	else{ Pvec[Chan.size3 - 1].push_back(i); }
	break;
      }
    }
  }

  Chan.size1 = Space.Chansize_2b_dir;
  Chan.size2 = Space.Chansize_2b_cross;

  std::vector<std::vector<int> > tempvec1(Chan.size3);
  std::vector<std::vector<int> > tempvec2(Chan.size3);

  //Make Vector of Two-Body States for Each Channel
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    int ind1, ind2;
    tempvec1[i].assign(Chan.size1, -1);
    tempvec2[i].assign(Chan.size2, -1);
    for(int j = 0; j < Chan.size3; ++j){
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, i, j);
      tempvec1[i][ind1] = j;
      ind2 = ChanInd_2b_cross(Parameters.basis, Space, i, j);
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
  V1.resize(Chan.size3);
  V2.resize(Chan.size3);
  V3.resize(Chan.size1);
  V4.resize(Chan.size1);
  V5.resize(Chan.size2);
  V6.resize(Chan.size3);
  V7.resize(Chan.size3);
  #pragma omp parallel for
  for(int i = 0; i < Chan.size3; ++i){
    V1[i].assign(Chan.p[i] * Chan.p[i], 0.0);
    V2[i].assign(Chan.h[i] * Chan.h[i], 0.0);
    V6[i].assign(Chan.hpp[i] * Chan.hpp[i], 0.0);
    V7[i].assign(Chan.hhp[i] * Chan.hhp[i], 0.0);
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
	Input.calc_case = substr;
	break;
      case 2:
	Input.basis = substr;
	break;
      case 3:
	Input.obstrength = atof(substr.c_str());
	break;
      case 4:
	Input.tbstrength = atof(substr.c_str());
	break;
      case 5:
	Input.Pshells = atoi(substr.c_str());
	break;
      case 6:
	Input.Nshells = atoi(substr.c_str());
	break;
      case 7:
	Input.LevelScheme = substr;
	break;
      case 8:
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
  std::cout << "Case = " << Parameters.calc_case << ", Basis = " << Parameters.basis << std::endl;
  if(Parameters.LevelScheme.size() > 0){ 
    std::cout << "Levels Scheme = " << Parameters.LevelScheme << std::endl;
    if(Parameters.MatrixElements.size() > 0){ std::cout << "Interaction = " << Parameters.MatrixElements << std::endl; }
  }
  else{ std::cout << "Nmax = " << Parameters.Nmax << ", Density = " << Parameters.density << std::endl; }
  std::cout << "OB strength = " << Parameters.obstrength << ", TB strength = " << Parameters.tbstrength << std::endl;
  std::cout << "Proton Shells = " << Parameters.Pshells << ", Neutron Shells = " << Parameters.Nshells << std::endl;
  std::cout << "Protons = " << Parameters.P << ", Neutrons = " << Parameters.N << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
}


Model_Space Build_Model_Space(Input_Parameters &Parameters)
{
  Model_Space Space; // Model Space information
  std::string fullpath; // Model Space file path
  std::string phline; // Std::String for each file line
  std::ifstream splevels; // Model space file
  int ind, n, l, nx, ny, nz;
  double energy, j, m, tz; // initialize level index, n, l, nx, ny, nz, m, tz depending on basis
  int pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;
  
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

  if(Parameters.basis == "finite"){
    Space.levelsn.resize(Space.indtot);
    Space.levelsl.resize(Space.indtot);
    Space.levelsj.resize(Space.indtot);
    Space.Pmin = 1000, Space.Pmax = -1000;
  }
  else if(Parameters.basis == "infinite"){
    Space.levelsnx.resize(Space.indtot);
    Space.levelsny.resize(Space.indtot);
    Space.levelsnz.resize(Space.indtot);
    Space.Nxmin = 1000, Space.Nxmax = -1000;
    Space.Nymin = 1000, Space.Nymax = -1000;
    Space.Nzmin = 1000, Space.Nzmax = -1000;
  }
  Space.Mmin = 1000, Space.Mmax = -1000;
  Space.Tmin = 1000, Space.Tmax = -1000;

  for(int i = 0; i < Space.indtot; ++i){
    if(Parameters.basis == "finite"){
      getline(splevels, phline);
      std::stringstream(phline) >> ind >> n >> l >> j >> m >> tz >> energy;
      Space.levelsind[i] = ind;
      Space.levelsn[i] = n;
      Space.levelsl[i] = l;
      Space.levelsj[i] = j;
      Space.levelsm[i] = m;
      Space.levelst[i] = tz;
      Space.levelsen[i] = energy * Parameters.obstrength;
      if(m < Space.Mmin){ Space.Mmin = m; }
      if(tz < Space.Tmin){ Space.Tmin = tz; }
      if(-2*(l%2) + 1 < Space.Pmin){ Space.Pmin = -2*(l%2) + 1; }
      if(m > Space.Mmax){ Space.Mmax = m; }
      if(tz > Space.Tmax){ Space.Tmax = tz; }
      if(-2*(l%2) + 1 > Space.Pmax){ Space.Pmax = -2*(l%2) + 1; }
    }
    else if(Parameters.basis == "infinite"){
      getline(splevels, phline);
      std::stringstream(phline) >> ind >> nx >> ny >> nz >> m >> tz >> energy;
      Space.levelsind[i] = ind;
      Space.levelsnx[i] = nx;
      Space.levelsny[i] = ny;
      Space.levelsnz[i] = nz;
      Space.levelsm[i] = m;
      Space.levelst[i] = tz;
      Space.levelsen[i] = energy * Parameters.obstrength;
      if(m < Space.Mmin){ Space.Mmin = m; }
      if(tz < Space.Tmin){ Space.Tmin = tz; }
      if(nx < Space.Nxmin){ Space.Nxmin = nx; }
      if(ny < Space.Nymin){ Space.Nymin = ny; }
      if(nz < Space.Nzmin){ Space.Nzmin = nz; }
      if(m > Space.Mmax){ Space.Mmax = m; }
      if(tz > Space.Tmax){ Space.Tmax = tz; }
      if(nx > Space.Nxmax){ Space.Nxmax = nx; }
      if(ny > Space.Nymax){ Space.Nymax = ny; }
      if(nz > Space.Nzmax){ Space.Nzmax = nz; }
    }
  }

  // find N2max_dir, N2max_cross
  int N2max_dir = 0;
  int N2max_cross = 0; // find maximum nx^2 + ny^2 + nz^2
  if(Parameters.basis == "infinite"){
    for(int i = 0; i < Space.indtot; ++i){
      for(int j = i+1; j < Space.indtot; ++j){
	if((Space.levelsnx[i] + Space.levelsnx[j])*(Space.levelsnx[i] + Space.levelsnx[j]) + 
	   (Space.levelsny[i] + Space.levelsny[j])*(Space.levelsny[i] + Space.levelsny[j]) + 
	   (Space.levelsnz[i] + Space.levelsnz[j])*(Space.levelsnz[i] + Space.levelsnz[j]) > N2max_cross){
	  N2max_cross = (Space.levelsnx[i] + Space.levelsnx[j])*(Space.levelsnx[i] + Space.levelsnx[j]) + 
	    (Space.levelsny[i] + Space.levelsny[j])*(Space.levelsny[i] + Space.levelsny[j]) + 
	    (Space.levelsnz[i] + Space.levelsnz[j])*(Space.levelsnz[i] + Space.levelsnz[j]); }
	if((Space.levelsnx[i] - Space.levelsnx[j])*(Space.levelsnx[i] - Space.levelsnx[j]) + 
	   (Space.levelsny[i] - Space.levelsny[j])*(Space.levelsny[i] - Space.levelsny[j]) + 
	   (Space.levelsnz[i] - Space.levelsnz[j])*(Space.levelsnz[i] - Space.levelsnz[j]) > N2max_cross){
	  N2max_cross = (Space.levelsnx[i] - Space.levelsnx[j])*(Space.levelsnx[i] - Space.levelsnx[j]) + 
	    (Space.levelsny[i] - Space.levelsny[j])*(Space.levelsny[i] - Space.levelsny[j]) + 
	    (Space.levelsnz[i] - Space.levelsnz[j])*(Space.levelsnz[i] - Space.levelsnz[j]); }
      }
    }
  }

  std::vector<int> ShellInd;
  double Etemp = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelsen[i] > Etemp){ ShellInd.push_back(i); Etemp = Space.levelsen[i]; }
  }
  if(Parameters.Pshells >= int(ShellInd.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells >= int(ShellInd.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelst[i] == -1){
      if(i < ShellInd[Parameters.Pshells]){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.levelst[i] == 1){
      if(i < ShellInd[Parameters.Nshells]){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
    }
  }
  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  Space.Msize = int((Space.Mmax - Space.Mmin)/2) + 1;
  Space.M2size = Space.Msize;
  Space.Tsize = int((Space.Tmax - Space.Tmin)/2) + 1;
  Space.T2size = Space.Tsize;
  if(Parameters.basis == "finite"){
    Space.Psize = int((Space.Pmax - Space.Pmin)/2) + 1;
    Space.P2size = Space.Psize;
    Space.Chansize_2b_dir = Space.Psize * Space.Msize * Space.Tsize;
    Space.Chansize_2b_cross = Space.P2size * Space.M2size * Space.T2size;
  }
  else if(Parameters.basis == "infinite"){
    int count1 = 0;
    int count2 = 0;
    Space.Nxsize = 2*(Space.Nxmax - Space.Nxmin) + 1; //sizes for tb combinations in each direction tb_cross has same size
    Space.Nx2size = Space.Nxsize;
    Space.Nysize = 2*(Space.Nymax - Space.Nymin) + 1;
    Space.Ny2size = Space.Nysize;
    Space.Nzsize = 2*(Space.Nzmax - Space.Nzmin) + 1;
    Space.Nz2size = Space.Nzsize;
    Space.map_2b_dir.assign(Space.Nxsize * Space.Nysize * Space.Nzsize, -1);
    Space.map_2b_cross.assign(Space.Nx2size * Space.Ny2size * Space.Nz2size, -1);
    for(int nx = 2*Space.Nxmin; nx <= 2*Space.Nxmax; ++nx){
      for(int ny = 2*Space.Nymin; ny <= 2*Space.Nymax; ++ny){
	for(int nz = 2*Space.Nzmin; nz <= 2*Space.Nzmax; ++nz){
	  if(nx*nx + ny*ny + nz*nz <= N2max_dir){
	    Space.map_2b_dir[(nx - 2*Space.Nxmin)*(Space.Nysize*Space.Nzsize) + (ny - 2*Space.Nymin)*(Space.Nzsize) + (nz - 2*Space.Nzmin)] = count1;
	    ++count1;
	  }
	}
      }
    }
    Space.Chansize_2b_dir = count1 * Space.Msize * Space.Tsize;
    for(int nx = (Space.Nxmin - Space.Nxmax); nx <= (Space.Nxmax - Space.Nxmin); ++nx){
      for(int ny = (Space.Nymin - Space.Nymax); ny <= (Space.Nymax - Space.Nymin); ++ny){
	for(int nz = (Space.Nzmin - Space.Nzmax); nz <= (Space.Nzmax - Space.Nzmin); ++nz){
	  if(nx*nx + ny*ny + nz*nz <= N2max_cross){
	    Space.map_2b_cross[(nx - Space.Nxmin + Space.Nxmax)*(Space.Ny2size*Space.Nz2size) + (ny - Space.Nymin +Space.Nymax)*(Space.Nz2size) + 
			       (nz - Space.Nzmin + Space.Nzmax)] = count2;
	    ++count2;
	  }
	}
      }
    }
    Space.Chansize_2b_cross = count2 * Space.M2size * Space.T2size;
  }

  splevels.close();
  
  return Space;
}


Model_Space_J Build_Model_Space_J1(Input_Parameters &Parameters)
{
  Model_Space_J Space_J; // Model Space information
  std::string fullpath; // Model Space file path
  std::string phline, number; // Std::String for each file line
  std::ifstream splevels; // Model space file
  std::istringstream phstream; // Stream of file line string
  int ind, n, l, l2n;
  double energy, j, tz; // initialize level index, n, l, nx, ny, nz, m, tz depending on basis
  size_t index1, index2; // Indicies for finding parameters among file lines
  int TotOrbs; // # total orbits

  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  //Read Model Space file, get Mass, Oscillator Energy, Total number of shells
  getline(splevels, phline);
  phstream.str(phline);
  phstream >> number;
  while (number != "Total"){ getline(splevels, phline); phstream.str(phline); phstream >> number; };
  index1 = phline.find_last_of(" \t");
  index2 = phline.find_last_of("0123456789");
  TotOrbs = std::atoi( phline.substr(index1 + 1, index2 - index1).c_str() );
  getline(splevels, phline);
  phstream.str(phline);
  phstream >> number;
  while(number != "Number:"){ getline(splevels, phline); phstream.str(phline); phstream >> number; };
  
  //read rest of level scheme parameters
  Space_J.indvec.resize(TotOrbs);
  Space_J.nvec.resize(TotOrbs);
  Space_J.lvec.resize(TotOrbs);
  Space_J.jvec.resize(TotOrbs);
  Space_J.tzvec.resize(TotOrbs);
  Space_J.envec.resize(TotOrbs);
  Space_J.envec.resize(TotOrbs);
  Space_J.shellsm.resize(TotOrbs);
  
  for(int i = 0; i < TotOrbs; ++i){
    phstream.str(phline);
    phstream >> number >> ind >> n >> l >> j >> tz >> l2n >> energy;
    //energy = energy * (1.0 - (1.0/Space.A));
    Space_J.indvec[i] = ind - 1;
    Space_J.nvec[i] = n;
    Space_J.lvec[i] = l;
    Space_J.jvec[i] = j;
    Space_J.tzvec[i] = tz;
    Space_J.envec[i] = energy;
    if( i < TotOrbs - 1 ){ getline(splevels, phline); };
  }
  splevels.close();

  return Space_J;
}


Model_Space Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space_J &Space_J)
{
  Model_Space Space; // Model Space information
  int shelllength; // number of single particle states for each shell
  int pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;
  int ind;
  double m;

  Space.indtot = 0;
  for(int i = 0; i < int(Space_J.indvec.size()); ++i){ Space.indtot += Space_J.jvec[i] + 1; }

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);
  Space.levelsn.resize(Space.indtot);
  Space.levelsl.resize(Space.indtot);
  Space.levelsj.resize(Space.indtot);

  Space.Mmin = 1000, Space.Mmax = -1000;
  Space.Tmin = 1000, Space.Tmax = -1000;
  Space.Pmin = 1000, Space.Pmax = -1000;

  ind = 0;
  for(int i = 0; i < int(Space_J.jvec.size()); ++i){
    shelllength = Space_J.jvec[i] + 1;
    Space_J.shellsm[i].resize(shelllength);
    for(int j = 0; j < shelllength; ++j){
      m = -Space_J.jvec[i] + 2*j;
      Space.levelsind[ind] = ind;
      Space.levelsn[ind] = Space_J.nvec[i];
      Space.levelsl[ind] = Space_J.lvec[i];
      Space.levelsj[ind] = Space_J.jvec[i];
      Space.levelsm[ind] = m;
      Space.levelst[ind] = Space_J.tzvec[i];
      Space.levelsen[ind] = Space_J.envec[i] * Parameters.obstrength;
      Space_J.shellsm[i][j] = ind;
      if(m < Space.Mmin){ Space.Mmin = m; }
      if(Space_J.tzvec[i] < Space.Tmin){ Space.Tmin = Space_J.tzvec[i]; }
      if(-2*(Space_J.lvec[i]%2) + 1 < Space.Pmin){ Space.Pmin = -2*(Space_J.lvec[i]%2) + 1; }
      if(m > Space.Mmax){ Space.Mmax = m; }
      if(Space_J.tzvec[i] > Space.Tmax){ Space.Tmax = Space_J.tzvec[i]; }
      if(-2*(Space_J.lvec[i]%2) + 1 > Space.Pmax){ Space.Pmax = -2*(Space_J.lvec[i]%2) + 1; }
      ++ind;
    }
  }

  std::vector<int> ShellInd;
  double Etemp = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelst[i] == -1 && Space.levelsen[i] > Etemp){ ShellInd.push_back(i); Etemp = Space.levelsen[i]; }
  }
  if(Parameters.Pshells >= int(ShellInd.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelst[i] == -1){
      if(i < ShellInd[Parameters.Pshells]){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    }
  }
  ShellInd.resize(0);
  Etemp = -1000.0;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelst[i] == 1 && Space.levelsen[i] > Etemp){ ShellInd.push_back(i); Etemp = Space.levelsen[i]; }
  }
  if(Parameters.Nshells >= int(ShellInd.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelst[i] == 1){
      if(i < ShellInd[Parameters.Nshells]){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
    }
  }

  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  Space.Msize = int((Space.Mmax - Space.Mmin)/2) + 1;
  Space.M2size = Space.Msize;
  Space.Tsize = int((Space.Tmax - Space.Tmin)/2) + 1;
  Space.T2size = Space.Tsize;
  Space.Psize = int((Space.Pmax - Space.Pmin)/2) + 1;
  Space.P2size = Space.Psize;
  Space.Chansize_2b_dir = Space.Psize * Space.Msize * Space.Tsize;
  Space.Chansize_2b_cross = Space.P2size * Space.M2size * Space.T2size;
  
  return Space;
}


Model_Space CART_Build_Model_Space(Input_Parameters &Parameters)
{
  Model_Space Space; // Model Space information
  double E;
  int count = 0, pcount = 0, ncount = 0, holcount = 0, parcount = 0, phcount = 0, nhcount = 0;
  
  double hbarc = 197.3269788; // MeVfm
  double m_neutronc2 = 939.565378; // MeV
  //double m_protonc2 = 938.272046; // MeV
  double m_protonc2 = 939.565378; // MeV
  double neutron_prefac = hbarc*hbarc/(2.0*m_neutronc2);
  double proton_prefac = hbarc*hbarc/(2.0*m_protonc2);

  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ 
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2*2;
    Space.Tmin = -1, Space.Tmax = 1, Space.Tsize = 2, Space.T2size = 2;
  }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.Tmin = -1, Space.Tmax = -1, Space.Tsize = 1, Space.T2size = 1;
  }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.Tmin = 1, Space.Tmax = 1, Space.Tsize = 1, Space.T2size = 1;
  }
  else{ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);
  Space.levelsnx.resize(Space.indtot);
  Space.levelsny.resize(Space.indtot);
  Space.levelsnz.resize(Space.indtot);

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
		Space.levelsen[count] = E;
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		E = 4.0*(nx*nx + ny*ny + nz*nz);
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
  }

  std::vector<int> ShellInd;
  double Etemp = -1000.0;
  for(int i = 0; i < count; ++i){
    if(Space.levelsen[i] > Etemp){ ShellInd.push_back(i); Etemp = Space.levelsen[i]; }
  }
  if(Parameters.Pshells >= int(ShellInd.size())){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells >= int(ShellInd.size())){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < count; ++i){
    if(Space.levelst[i] == -1){
      if(i < ShellInd[Parameters.Pshells]){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
      else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
    }
    else if(Space.levelst[i] == 1){
      if(i < ShellInd[Parameters.Nshells]){ Space.levelstype[i] = "hole"; ++ncount; ++holcount; ++nhcount; }
      else{ Space.levelstype[i] = "particle"; ++ncount; ++parcount; }
    }
  }

  Space.indp = pcount;
  Space.indn = ncount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;
  Parameters.N = nhcount;

  double L = pow(holcount/Parameters.density, 1.0/3.0);
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

  Space.Mmin = -1;
  Space.Mmax = 1;
  Space.Msize = 2;
  Space.M2size = 2;
  int count1 = 0;
  int count2 = 0;
  int N2max_dir = 4*Parameters.Nmax;
  int N2max_cross = 4*Parameters.Nmax;
  Space.Nxmin = -nmax;
  Space.Nxmax = nmax;
  Space.Nymin = -nmax;
  Space.Nymax = nmax;
  Space.Nzmin = -nmax;
  Space.Nzmax = nmax;
  Space.Nxsize = 4*nmax + 1; //sizes for tb combinations in each direction tb_cross has same size
  Space.Nysize = 4*nmax + 1;
  Space.Nzsize = 4*nmax + 1;
  Space.Nx2size = 4*nmax + 1;
  Space.Ny2size = 4*nmax + 1;
  Space.Nz2size = 4*nmax + 1;
  Space.map_2b_dir.assign(Space.Nxsize * Space.Nysize * Space.Nzsize, -1);
  Space.map_2b_cross.assign(Space.Nx2size * Space.Ny2size * Space.Nz2size, -1);
  for(int nx = 2*Space.Nxmin; nx <= 2*Space.Nxmax; ++nx){
    for(int ny = 2*Space.Nymin; ny <= 2*Space.Nymax; ++ny){
      for(int nz = 2*Space.Nzmin; nz <= 2*Space.Nzmax; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max_dir){
	  Space.map_2b_dir[(nx - 2*Space.Nxmin)*(Space.Nysize*Space.Nzsize) + (ny - 2*Space.Nymin)*(Space.Nzsize) + (nz - 2*Space.Nzmin)] = count1;
	  ++count1;
	}
      }
    }
  }
  Space.Chansize_2b_dir = count1 * Space.Msize * Space.Tsize;
  for(int nx = (Space.Nxmin - Space.Nxmax); nx <= (Space.Nxmax - Space.Nxmin); ++nx){
    for(int ny = (Space.Nymin - Space.Nymax); ny <= (Space.Nymax - Space.Nymin); ++ny){
      for(int nz = (Space.Nzmin - Space.Nzmax); nz <= (Space.Nzmax - Space.Nzmin); ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max_cross){
	  Space.map_2b_cross[(nx - Space.Nxmin + Space.Nxmax)*(Space.Ny2size*Space.Nz2size) + (ny - Space.Nymin +Space.Nymax)*(Space.Nz2size) + 
			     (nz - Space.Nzmin + Space.Nzmax)] = count2;
	  ++count2;
	}
      }
    }
  }
  Space.Chansize_2b_cross = count2 * Space.M2size * Space.T2size;

  return Space;
}


Model_Space EG_Build_Model_Space(Input_Parameters &Parameters)
{
  Model_Space Space; // Model Space information
  double E;
  int count = 0, pcount = 0, holcount = 0, parcount = 0, phcount = 0;
  Parameters.Nshells = 0;

  double hbarc = 0.1973269788; // eV um
  double eVs_in_Hartree = 27.21138505; // eV
  hbarc *= 10000/eVs_in_Hartree; // Hartree Angstrom
  double massc2 = 0.5109989461; // MeV
  massc2 *= 1000000/eVs_in_Hartree; // Hartree
  double prefac = hbarc*hbarc/(2.0*massc2);

  if(Parameters.Pshells != 0){ 
    Space.indtot = (2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*(2*Parameters.Nmax+1)*2;
    Space.Tmin = -1, Space.Tmax = -1, Space.Tsize = 1, Space.T2size = 1;
  }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  Space.levelsind.resize(Space.indtot);
  Space.levelsm.resize(Space.indtot);
  Space.levelst.resize(Space.indtot);
  Space.levelstype.resize(Space.indtot);
  Space.levelsen.resize(Space.indtot);
  Space.levelsnx.resize(Space.indtot);
  Space.levelsny.resize(Space.indtot);
  Space.levelsnz.resize(Space.indtot);

  int nmax = std::floor(std::sqrt(Parameters.Nmax));
  for(int shell = 0; shell <= Parameters.Nmax; ++shell){
    for(int nx = -nmax; nx <= nmax; ++nx){    
      for(int ny = -nmax; ny <= nmax; ++ny){	
	for(int nz = -nmax; nz <= nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz || Parameters.Nmax <= nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    E = 4.0*(nx*nx + ny*ny + nz*nz);
	    Space.levelsen[count] = E;
	    Space.levelsind[count] = count + 1;
	    Space.levelsnx[count] = nx;
	    Space.levelsny[count] = ny;
	    Space.levelsnz[count] = nz;
	    Space.levelsm[count] = sz;
	    Space.levelst[count] = -1;
	    count++;
	  }   
	} 
      }
    }
  }

  std::vector<int> ShellInd;
  double Etemp = -1000.0;
  for(int i = 0; i < count; ++i){
    if(Space.levelsen[i] > Etemp){ ShellInd.push_back(i); Etemp = Space.levelsen[i]; }
  }
  if(Parameters.Pshells >= int(ShellInd.size())){ std::cerr << "Number of shells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < count; ++i){
    if(i < ShellInd[Parameters.Pshells]){ Space.levelstype[i] = "hole"; ++pcount; ++holcount; ++phcount; }
    else{ Space.levelstype[i] = "particle"; ++pcount; ++parcount; }
  }

  Space.indp = pcount;
  Space.indpar = parcount;
  Space.indhol = holcount;
  Parameters.P = phcount;

  double L = pow(holcount/Parameters.density, 1.0/3.0);
  for(int i = 0; i < count; ++i){ Space.levelsen[i] *= prefac*M_PI*M_PI/(L*L); }

  Space.levelsind.resize(count);
  Space.levelsm.resize(count);
  Space.levelst.resize(count);
  Space.levelstype.resize(count);
  Space.levelsen.resize(count);
  Space.levelsnx.resize(count);
  Space.levelsny.resize(count);
  Space.levelsnz.resize(count);
  Space.indtot = count;

  Space.Mmin = -1;
  Space.Mmax = 1;
  Space.Tmin = -1;
  Space.Tmax = -1;
  Space.Msize = 2;
  Space.M2size = 2;
  Space.Tsize = 1;
  Space.T2size = 1;
  int count1 = 0;
  int count2 = 0;
  int N2max_dir = 4*Parameters.Nmax;
  int N2max_cross = 4*Parameters.Nmax;
  Space.Nxmin = -nmax;
  Space.Nxmax = nmax;
  Space.Nymin = -nmax;
  Space.Nymax = nmax;
  Space.Nzmin = -nmax;
  Space.Nzmax = nmax;
  Space.Nxsize = 4*nmax + 1; //sizes for tb combinations in each direction tb_cross has same size
  Space.Nysize = 4*nmax + 1;
  Space.Nzsize = 4*nmax + 1;
  Space.Nx2size = 4*nmax + 1;
  Space.Ny2size = 4*nmax + 1;
  Space.Nz2size = 4*nmax + 1;
  Space.map_2b_dir.assign(Space.Nxsize * Space.Nysize * Space.Nzsize, -1);
  Space.map_2b_cross.assign(Space.Nx2size * Space.Ny2size * Space.Nz2size, -1);
  for(int nx = 2*Space.Nxmin; nx <= 2*Space.Nxmax; ++nx){
    for(int ny = 2*Space.Nymin; ny <= 2*Space.Nymax; ++ny){
      for(int nz = 2*Space.Nzmin; nz <= 2*Space.Nzmax; ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max_dir){
	  Space.map_2b_dir[(nx - 2*Space.Nxmin)*(Space.Nysize*Space.Nzsize) + (ny - 2*Space.Nymin)*(Space.Nzsize) + (nz - 2*Space.Nzmin)] = count1;
	  ++count1;
	}
      }
    }
  }
  Space.Chansize_2b_dir = count1 * Space.Msize * Space.Tsize;
  for(int nx = (Space.Nxmin - Space.Nxmax); nx <= (Space.Nxmax - Space.Nxmin); ++nx){
    for(int ny = (Space.Nymin - Space.Nymax); ny <= (Space.Nymax - Space.Nymin); ++ny){
      for(int nz = (Space.Nzmin - Space.Nzmax); nz <= (Space.Nzmax - Space.Nzmin); ++nz){
	if(nx*nx + ny*ny + nz*nz <= N2max_cross){
	  Space.map_2b_cross[(nx - Space.Nxmin + Space.Nxmax)*(Space.Ny2size*Space.Nz2size) + (ny - Space.Nymin +Space.Nymax)*(Space.Nz2size) + 
			     (nz - Space.Nzmin + Space.Nzmax)] = count2;
	  ++count2;
	}
      }
    }
  }
  Space.Chansize_2b_cross = count2 * Space.M2size * Space.T2size;
  
  return Space;
}


double Partial_Wave_HF(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, const double &maxL)
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
}



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

int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const int &ind1, const int &ind2)
{
  int M = Space.levelsm[ind1] + Space.levelsm[ind2];
  int T = Space.levelst[ind1] + Space.levelst[ind2];
  if(basis == "infinite"){
    int Nx = Space.levelsnx[ind1] + Space.levelsnx[ind2];
    int Ny = Space.levelsny[ind1] + Space.levelsny[ind2];
    int Nz = Space.levelsnz[ind1] + Space.levelsnz[ind2];
    return Space.map_2b_dir[(Nx - 2*Space.Nxmin)*(Space.Nysize*Space.Nzsize) + (Ny - 2*Space.Nymin)*(Space.Nzsize) + (Nz - 2*Space.Nzmin)]*(Space.Msize*Space.Tsize)
      + int((M - 2*Space.Mmin)/2)*Space.Tsize + int((T - 2*Space.Tmin)/2);
  }
  else{
    double P = -2*((Space.levelsl[ind1] + Space.levelsl[ind2])%2) + 1;
    return int((P - Space.Pmin)/2)*(Space.Msize*Space.Tsize) + int((M - 2*Space.Mmin)/2)*Space.Tsize + int((T - 2*Space.Tmin)/2);
  }
}

int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const int &ind1, const int &ind2)
{
  int M = Space.levelsm[ind1] - Space.levelsm[ind2];
  int T = Space.levelst[ind1] - Space.levelst[ind2];
  if(basis == "infinite"){
    int Nx = Space.levelsnx[ind1] - Space.levelsnx[ind2];
    int Ny = Space.levelsny[ind1] - Space.levelsny[ind2];
    int Nz = Space.levelsnz[ind1] - Space.levelsnz[ind2];
    return Space.map_2b_cross[(Nx - Space.Nxmin + Space.Nxmax)*(Space.Ny2size*Space.Nz2size) + (Ny - Space.Nymin + Space.Nymax)*(Space.Nzsize)
			      + (Nz - Space.Nzmin + Space.Nzmax)]*(Space.M2size*Space.T2size) + int((M - Space.Mmin + Space.Mmax)/2) * Space.T2size 
      + int((T - Space.Tmin + Space.Tmax)/2);
  }
  else{
    double P = -2*((Space.levelsl[ind1] - Space.levelsl[ind2])%2) + 1;
    return int((P - Space.Pmin)/2)*(Space.M2size*Space.T2size) + int((M - Space.Mmin + Space.Mmax)/2)*Space.T2size + int((T - Space.Tmin + Space.Tmax)/2);
  }
}

/*int CART_tbInd2_rel(const Model_Space &Space, const int &Nx2, const int &Ny2, const int &Nz2)
{
  int ind = Space.CART_tb2Indvec[Nx2-Space.Nx2min][Ny2-Space.Ny2min[Nx2-Space.Nx2min]][Nz2-Space.Nz2min[Nx2-Space.Nx2min][Ny2-Space.Ny2min[Nx2-Space.Nx2min]]];
  return ind;
  }*/


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
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
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
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2); 
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
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
      CCME.HPHP1[ind1][ind] = TBME;
      ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
      CCME.HPHP1[ind1][ind] = TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell4);
      ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell1, shell4, shell3, shell2);
      CCME.HPHP2[ind1][ind] = TBME;
      ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], shell3, shell2, shell1, shell4);
      CCME.HPHP2[ind1][ind] = TBME;
    }
    else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
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
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
      CCME.HHPP4[ind1][ind] = TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
      CCME.HHPP4[ind1][ind] = -1.0 * TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
      CCME.HHPP4[ind1][ind] = -1.0 * TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
      CCME.HHPP4[ind1][ind] = TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
      CCME.HHPP4T[ind1][ind] = TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
      CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
      CCME.HHPP4T[ind1][ind] = -1.0 * TBME;
      ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
      ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
      CCME.HHPP4T[ind1][ind] = TBME;
    }  
  }
  
  interaction.close();
  
  return CCME;

}


CC_Matrix_Elements Read_Matrix_Elements_J(const std::string &MEfile, const Input_Parameters &Parameters, const Model_Space &Space, const Model_Space_J &Space_J, const Channels &Chan)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  std::string fullpath1, fullpath2; // file path string and string with "_M" added (for J,M-Scheme)
  size_t intsize; // size of MEfile string
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME0, TBME, CMME, m1, m2, t1, t2, CGC1, CGC2, coupJ, coupT; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, par; // interaction file contents
  int pind, qind, rind, sind;
  CC_Matrix_Elements CCME(Chan);
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;

  fullpath1 = PATH + MEfile + ".int";
  intsize = MEfile.size();

  //open interaction file
  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << MEfile << ", does not exist" << std::endl; exit(1); }

  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while(number != "Total"){ 
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
    std::istringstream(interactionline) >> coupT >> par >> coupJ >> shell1 >> shell2 >> shell3 >> shell4 >> TBME0 >> CMME;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    //TBME = TBME - (Space.HOEnergy/Space.A)*CMME;
    if(shell2 < shell1){
      std::swap(shell1, shell2);
      TBME0 = TBME0 * -1.0 * pow(-1.0, int(0.5*(Space_J.jvec[shell1] + Space_J.jvec[shell2] - coupJ)));
    }
    if(shell4 < shell3){
      std::swap(shell3, shell4);
      TBME0 = TBME0 * -1.0 * pow(-1.0, int(0.5*(Space_J.jvec[shell3] + Space_J.jvec[shell4] - coupJ)));
    }
    if((shell3 < shell1) || (shell3 == shell1 && shell4 < shell2)){
      std::swap(shell1, shell3);
      std::swap(shell2, shell4);
    }
    if(shell1 == shell2){ TBME0 *= sqrt(2.0); }
    if(shell3 == shell4){ TBME0 *= sqrt(2.0); }
    for(int jz = -coupJ; jz <= coupJ; jz+=2){
      for(int p = 0; p < int(Space_J.shellsm[shell1].size()); ++p){
	for(int q = 0; q < int(Space_J.shellsm[shell2].size()); ++q){
	  for(int r = 0; r < int(Space_J.shellsm[shell3].size()); ++r){
	    for(int s = 0; s < int(Space_J.shellsm[shell4].size()); ++s){
	      pind = Space_J.shellsm[shell1][p];
	      qind = Space_J.shellsm[shell2][q];
	      rind = Space_J.shellsm[shell3][r];
	      sind = Space_J.shellsm[shell4][s];
	      m1 = Space.levelsm[pind] + Space.levelsm[qind];
	      m2 = Space.levelsm[rind] + Space.levelsm[sind];
	      t1 = Space.levelst[pind] + Space.levelst[qind];
	      t2 = Space.levelst[rind] + Space.levelst[sind];
	      if(t1 != t2 || m1 != jz || m2 != jz){ continue; }
	      if(pind >= qind && shell1 == shell2){ continue; }
	      if(rind >= sind && shell3 == shell4){ continue; }
	      if(pind > rind && shell1 == shell3 && shell2 == shell4){ continue; }
	      if(pind == rind && qind > sind && shell1 == shell3 && shell2 == shell4){ continue; }

	      CGC1 = CGC(0.5*Space.levelsj[pind], 0.5*Space.levelsm[pind], 0.5*Space.levelsj[qind], 0.5*Space.levelsm[qind], 0.5*coupJ, double(0.5*jz));
	      CGC2 = CGC(0.5*Space.levelsj[rind], 0.5*Space.levelsm[rind], 0.5*Space.levelsj[sind], 0.5*Space.levelsm[sind], 0.5*coupJ, double(0.5*jz));
	      TBME = TBME0 * (Parameters.tbstrength * CGC1 * CGC2);
	      if(fabs(TBME) < 1e-12){ continue; }
	      ptype = Space.levelstype[pind];
	      qtype = Space.levelstype[qind];
	      rtype = Space.levelstype[rind];
	      stype = Space.levelstype[sind];
	      if(ptype == "particle" && qtype == "hole"){ std::swap(pind, qind); std::swap(ptype, qtype); TBME *= -1.0; }
	      if(rtype == "particle" && stype == "hole"){ std::swap(rind, sind); std::swap(rtype, stype); TBME *= -1.0; }
	      if(ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "hole"){
		std::swap(pind, rind);
		std::swap(qind, sind);
		std::swap(ptype, rtype);
		std::swap(qtype, stype);
	      }
	      if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind, qind);
		ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], pind, qind, rind, sind);
		CCME.HHHH[ind1][ind] += TBME;
		ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], qind, pind, rind, sind);
		CCME.HHHH[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], pind, qind, sind, rind);
		CCME.HHHH[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], qind, pind, sind, rind);
		CCME.HHHH[ind1][ind] += TBME;
		if(rind != pind && sind != qind){
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], rind, sind, pind, qind);
		  CCME.HHHH[ind1][ind] += TBME;
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], sind, rind, pind, qind);
		  CCME.HHHH[ind1][ind] += -1.0 * TBME;
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], rind, sind, qind, pind);
		  CCME.HHHH[ind1][ind] += -1.0 * TBME;
		  ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], sind, rind, qind, pind);
		  CCME.HHHH[ind1][ind] += TBME;
		}
	      }
	      else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind, qind); 
		ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], pind, qind, rind, sind);
		CCME.PPPP[ind1][ind] += TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], qind, pind, rind, sind);
		CCME.PPPP[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], pind, qind, sind, rind);
		CCME.PPPP[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], qind, pind, sind, rind);
		CCME.PPPP[ind1][ind] += TBME;
		if(rind != pind && sind != qind){
		  ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], rind, sind, pind, qind);
		  CCME.PPPP[ind1][ind] += TBME;
		  ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], sind, rind, pind, qind);
		  CCME.PPPP[ind1][ind] += -1.0 * TBME;
		  ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], rind, sind, qind, pind);
		  CCME.PPPP[ind1][ind] += -1.0 * TBME;
		  ind = Index(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], sind, rind, qind, pind);
		  CCME.PPPP[ind1][ind] += TBME;
		}
	      }
	      else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, sind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], pind, sind, rind, qind);
		CCME.HPHP1[ind1][ind] += TBME;
		if(rind != pind && sind != qind){
		  ind = Index(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], rind, qind, pind, sind);
		  CCME.HPHP1[ind1][ind] += TBME;
		}
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, pind, sind);
		ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], pind, sind, rind, qind);
		CCME.HPHP2[ind1][ind] += TBME;
		if(rind != pind && sind != qind){
		  ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], rind, qind, pind, sind);
		  CCME.HPHP2[ind1][ind] += TBME;
		}
	      }
	      else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, pind, qind);
		ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], rind, sind, pind, qind);
		CCME.HHPP1[ind1][ind] += TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], sind, rind, pind, qind);
		CCME.HHPP1[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], rind, sind, qind, pind);
		CCME.HHPP1[ind1][ind] += -1.0 * TBME;
		ind = Index(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], sind, rind, qind, pind);
		CCME.HHPP1[ind1][ind] += TBME;
		for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[rind] == j){ ind2 = j; break; }; }
		ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], rind, pind, qind, sind);
		CCME.HHPP2[ind2][ind] += TBME;
		ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], rind, qind, pind, sind);
		CCME.HHPP2[ind2][ind] += -1.0 * TBME;
		for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[sind] == j){ ind2 = j; break; }; }
		ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], sind, pind, qind, rind);
		CCME.HHPP2[ind2][ind] += -1.0 * TBME;
		ind = Index2(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], sind, qind, pind, rind);
		CCME.HHPP2[ind2][ind] += TBME;
		for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[pind] == j){ ind2 = j; break; }; }
		ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], pind, qind, rind, sind);
		CCME.HHPP3[ind2][ind] += TBME;
		ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], pind, qind, sind, rind);
		CCME.HHPP3[ind2][ind] += -1.0 * TBME;
		for(int j = 0; j < Chan.size3; ++j){ if(Chan.indvec[qind] == j){ ind2 = j; break; }; }
		ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], qind, pind, rind, sind);
		CCME.HHPP3[ind2][ind] += -1.0 * TBME;
		ind = Index2(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], qind, pind, sind, rind);
		CCME.HHPP3[ind2][ind] += TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, rind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], pind, rind, qind, sind);
		CCME.HHPP4[ind1][ind] += TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, rind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], qind, rind, pind, sind);
		CCME.HHPP4[ind1][ind] += -1.0 * TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, sind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], pind, sind, qind, rind);
		CCME.HHPP4[ind1][ind] += -1.0 * TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, sind, qind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], qind, sind, pind, rind);
		CCME.HHPP4[ind1][ind] += TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, sind, qind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], qind, sind, pind, rind);
		CCME.HHPP4T[ind1][ind] += TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, sind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], pind, sind, qind, rind);
		CCME.HHPP4T[ind1][ind] += -1.0 * TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, rind, qind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], qind, rind, pind, sind);
		CCME.HHPP4T[ind1][ind] += -1.0 * TBME;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, rind, pind);
		ind = Index(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], pind, rind, qind, sind);
		CCME.HHPP4T[ind1][ind] += TBME;
	      }
	    }
	  }
	}
      }
    }
  }

  interaction.close();
  
  return CCME;

}


CC_Matrix_Elements Coulomb_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  CC_Matrix_Elements CCME(Chan);
  
  double L = pow(Parameters.P/Parameters.density, 1.0/3.0);

  for(int i = 0; i < Chan.size1; ++i){
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.hh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec1[i][2*pq];
      shell2 = Chan.hhvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.hh[i]; ++rs){
	shell3 = Chan.hhvec1[i][2*rs];
	shell4 = Chan.hhvec1[i][2*rs + 1];
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHHH[i][Chan.hh[i]*pq + rs] = TBME;
	CCME.HHHH[i][Chan.hh[i]*rs + pq] = TBME;
      }
    }
    #pragma omp parallel for
    for(int pq = 0; pq < Chan.pp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec1[i][2*pq];
      shell2 = Chan.ppvec1[i][2*pq + 1];
      for(int rs = pq; rs < Chan.pp[i]; ++rs){
	shell3 = Chan.ppvec1[i][2*rs];
	shell4 = Chan.ppvec1[i][2*rs + 1];
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.PPPP[i][Chan.pp[i]*pq + rs] = TBME;
	CCME.PPPP[i][Chan.pp[i]*rs + pq] = TBME;
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
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHPP1[i][Chan.hh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int r = 0; r < Chan.p[i]; ++r){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell3 = Chan.pvec1[i][r];
      for(int pqs = 0; pqs < Chan.hhp[i]; ++pqs){
	shell1 = Chan.hhpvec1[i][3*pqs];
	shell2 = Chan.hhpvec1[i][3*pqs + 1];
	shell4 = Chan.hhpvec1[i][3*pqs + 2];
	if(shell3 == shell4){ continue; }
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHPP2[i][Chan.hhp[i]*r + pqs] = TBME;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int p = 0; p < Chan.h[i]; ++p){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hvec1[i][p];
      for(int qrs = 0; qrs < Chan.hpp[i]; ++qrs){
	shell2 = Chan.hppvec1[i][3*qrs];
	shell3 = Chan.hppvec1[i][3*qrs + 1];
	shell4 = Chan.hppvec1[i][3*qrs + 2];
	if(shell1 == shell2){ continue; }
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHPP3[i][Chan.hpp[i]*p + qrs] = TBME;
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
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HPHP1[i][Chan.hp2[i]*ps + qr] = TBME;
	CCME.HPHP1[i][Chan.hp2[i]*qr + ps] = TBME;
      }
    }
    #pragma omp parallel for
    for(int ps = 0; ps < Chan.hp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*ps];
      shell4 = Chan.hp2vec1[i][2*ps + 1];
      for(int qs = 0; qs < Chan.hp1[i]; ++qs){
	shell2 = Chan.hp1vec1[i][2*qs];
	shell3 = Chan.hp1vec1[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = Coulomb(Space, shell1, shell2, shell4, shell3, L);
	CCME.HHPP4[i][Chan.hp1[i]*ps + qs] = TBME;
	TBME = Coulomb(Space, shell2, shell1, shell3, shell4, L);
	CCME.HHPP4T[i][Chan.hp1[i]*ps + qs] = TBME;
      }
    }
    #pragma omp parallel for
    for(int qr = 0; qr < Chan.hp1[i]; ++qr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell3 = Chan.hp1vec1[i][2*qr];
      shell2 = Chan.hp1vec1[i][2*qr + 1];
      for(int ps = qr; ps < Chan.hp1[i]; ++ps){
	shell1 = Chan.hp1vec1[i][2*ps];
	shell4 = Chan.hp1vec1[i][2*ps + 1];
	TBME = Coulomb(Space, shell1, shell2, shell3, shell4, L);
	CCME.HPHP2[i][Chan.hp1[i]*qr + ps] = TBME;
	CCME.HPHP2[i][Chan.hp1[i]*ps + qr] = TBME;
      }
    }
  }
  
  return CCME;

}


// Minnesota Potential for momentum basis
double Coulomb(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double hbarc = 0.1973269788; // eV um
  double eVs_in_Hartree = 27.21138505; // eV
  hbarc *= 10000/eVs_in_Hartree; // Hartree Angstrom
  double term = 0.0;
  double fine_struct = 1.0/137.035999139;
  double e_sq = hbarc*fine_struct;
  double prefactor = e_sq/(L*L*L);
  double qSquared1;
  double kX1, kY1, kZ1;

  if(Space.levelsnx[qi] + Space.levelsnx[qj] != Space.levelsnx[qk] + Space.levelsnx[ql] || 
     Space.levelsny[qi] + Space.levelsny[qj] != Space.levelsny[qk] + Space.levelsny[ql] || 
     Space.levelsnz[qi] + Space.levelsnz[qj] != Space.levelsnz[qk] + Space.levelsnz[ql]){ return 0.0; }
  if(Space.levelsm[qi] == Space.levelsm[qk] && Space.levelsm[qj] == Space.levelsm[ql]){
    kX1 = (2.0*M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[qk]);
    kY1 = (2.0*M_PI/L) * (Space.levelsny[qi] - Space.levelsny[qk]);
    kZ1 = (2.0*M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[qk]);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-8){ term += 0.0; }
    else{ term += prefactor*4.0*M_PI/qSquared1; }
  }
  if(Space.levelsm[qi] == Space.levelsm[ql] && Space.levelsm[qj] == Space.levelsm[qk]){
    kX1 = (2.0*M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[ql]);
    kY1 = (2.0*M_PI/L) * (Space.levelsny[qi] - Space.levelsny[ql]);
    kZ1 = (2.0*M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[ql]);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-8){ term -= 0.0; }
    else{ term -= prefactor*4.0*M_PI/qSquared1; }
  }
  return term;
}


CC_Matrix_Elements Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  CC_Matrix_Elements CCME(Chan);
  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);

  for(int i = 0; i < Chan.size1; ++i){
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
	CCME.HHHH[i][Chan.hh[i]*pq + rs] = TBME;
	CCME.HHHH[i][Chan.hh[i]*rs + pq] = TBME;
      }
    }
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
	CCME.PPPP[i][Chan.pp[i]*pq + rs] = TBME;
	CCME.PPPP[i][Chan.pp[i]*rs + pq] = TBME;
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
	CCME.HHPP1[i][Chan.hh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int r = 0; r < Chan.p[i]; ++r){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell3 = Chan.pvec1[i][r];
      for(int pqs = 0; pqs < Chan.hhp[i]; ++pqs){
	shell1 = Chan.hhpvec1[i][3*pqs];
	shell2 = Chan.hhpvec1[i][3*pqs + 1];
	shell4 = Chan.hhpvec1[i][3*pqs + 2];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHPP2[i][Chan.hhp[i]*r + pqs] = TBME;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    #pragma omp parallel for
    for(int p = 0; p < Chan.h[i]; ++p){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hvec1[i][p];
      for(int qrs = 0; qrs < Chan.hpp[i]; ++qrs){
	shell2 = Chan.hppvec1[i][3*qrs];
	shell3 = Chan.hppvec1[i][3*qrs + 1];
	shell4 = Chan.hppvec1[i][3*qrs + 2];
	if(shell1 == shell2){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	CCME.HHPP3[i][Chan.hpp[i]*p + qrs] = TBME;
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
	CCME.HPHP1[i][Chan.hp2[i]*ps + qr] = TBME;
	CCME.HPHP1[i][Chan.hp2[i]*qr + ps] = TBME;
      }
    }
    #pragma omp parallel for
    for(int ps = 0; ps < Chan.hp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec1[i][2*ps];
      shell4 = Chan.hp2vec1[i][2*ps + 1];
      for(int qs = 0; qs < Chan.hp1[i]; ++qs){
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
    for(int qr = 0; qr < Chan.hp1[i]; ++qr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell3 = Chan.hp1vec1[i][2*qr];
      shell2 = Chan.hp1vec1[i][2*qr + 1];
      for(int ps = qr; ps < Chan.hp1[i]; ++ps){
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
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  kX1 = (M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[qj] - Space.levelsnx[qk] + Space.levelsnx[ql]);
  kY1 = (M_PI/L) * (Space.levelsny[qi] - Space.levelsny[qj] - Space.levelsny[qk] + Space.levelsny[ql]);
  kZ1 = (M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[qj] - Space.levelsnz[qk] + Space.levelsnz[ql]);

  kX2 = (M_PI/L) * (Space.levelsnx[qi] - Space.levelsnx[qj] - Space.levelsnx[ql] + Space.levelsnx[qk]);
  kY2 = (M_PI/L) * (Space.levelsny[qi] - Space.levelsny[qj] - Space.levelsny[ql] + Space.levelsny[qk]);
  kZ2 = (M_PI/L) * (Space.levelsnz[qi] - Space.levelsnz[qj] - Space.levelsnz[ql] + Space.levelsnz[qk]);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
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

CCD Perform_CCD(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, const Channels &Chan)
{
  std::cout << "Performing CCD ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  
  int ind = 0;
  int ind1, ind2, ind3, ind4;
  double tempen, tempt, error;
  CCD CCin(Chan, Parameters, Space);
  CCD CCout = CCin;

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

  ////////////////////////////////////////////////////////////////

  for(int chan = 0; chan < Chan.size1; ++chan){
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
  
  while((error > 10e-16 && ind < 1000) || ind < 10){
    Doubles_Step(Space, Chan, CCME, CCin, CCout);
    
    CCout.CCDE = 0.0;
    error = 0.0;
    for(int chan = 0; chan < Chan.size1; ++chan){
      for(int hind = 0; hind < Chan.hh[chan]; ++hind){
	for(int pind = 0; pind < Chan.pp[chan]; ++pind){
	  tempt = CCout.get_T(chan, hind * Chan.pp[chan] + pind);
	  tempt += CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	  tempt /= CCout.Evec[chan][hind * Chan.pp[chan] + pind];
	  //tempt = 0.25*tempt + 0.75*CCin.T1[chan][hind * Chan.pp[chan] + pind];
	  error += (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]) * (tempt - CCin.T1[chan][hind * Chan.pp[chan] + pind]);
	  CCin.set_T(chan, hind * Chan.pp[chan] + pind, tempt);
	  CCout.CCDE += 0.25 * tempt * CCME.HHPP1[chan][pind * Chan.hh[chan] + hind];
	}
      }
    }
    error = std::sqrt(error);
    CCin.CCDE = CCout.CCDE;
    
    std::cout << "Iteration Number = " << ind << ", CCD Energy = " << CCout.CCDE << ", error = " << error << std::endl;
    ++ind;
  }
  std::cout << std::endl << std::endl;

  return CCin;

}


void HF(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME)
{
  int ind, ind1;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelstype[i] == "hole"){
      for(int j = 0; j < Space.indtot; ++j){
	if(Space.levelstype[j] != "hole" || i == j){ continue; }
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, i, j);
	ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
	Space.levelsen[i] += ME.HHHH[ind1][ind];
      }
    }
    else{
      for(int j = 0; j < Space.indtot; ++j){
	if(Space.levelstype[j] != "hole"){ continue; }
	ind1 = ChanInd_2b_cross(Parameters.basis, Space, j, i);
	ind = Index(Chan.hp1vec1[ind1], Chan.hp1vec1[ind1], Chan.hp1[ind1], Chan.hp1[ind1], j, i, j, i);
	Space.levelsen[i] += ME.HPHP2[ind1][ind];
      }	
    }
  }
}


double E_Ref(const Input_Parameters &Parameters, Model_Space &Space, const Channels &Chan, const CC_Matrix_Elements &ME)
{
  double energy = 0;
  int ind, ind1;
  for(int i = 0; i < Space.indtot; ++i){
    if(Space.levelstype[i] != "hole"){ continue; }
    energy += Space.levelsen[i];
    for(int j = 0; j < Space.indtot; ++j){
      if(Space.levelstype[j] != "hole" || i == j){ continue; }
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, i, j);
      ind = Index(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], i, j, i, j);
      energy -= 0.5 * ME.HHHH[ind1][ind];
    }
  }
  return energy;
}

void Doubles_Step(const Model_Space &Space, const Channels &Chan, CC_Matrix_Elements &ME, CCD &CC, CCD &CC2)
{
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = 0.25, fac5 = -1.0, fac6 = -0.5;
  char N = 'N';
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.hp1[chan];
    int hp2 = Chan.hp2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    //S4 = T2*HHPP4
    RM_dgemm(& *CC.T2[chan].begin(), & *ME.HHPP4[chan].begin(), & *CC2.S4[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac2, &N, &N);
    //S4T = T2*HHPP4T
    RM_dgemm(& *CC.T2[chan].begin(), & *ME.HHPP4T[chan].begin(), & *CC2.S4T[chan].begin(), &hp1, &hp1, &hp2, &fac1, &fac2, &N, &N);
    //T2 & T3
    RM_dgemm(& *CC.T2[chan].begin(), & *ME.HPHP1[chan].begin(), & *CC2.T2[chan].begin(), &hp1, &hp2, &hp2, &fac1, &fac2, &N, &N);
    RM_dgemm(& *ME.HPHP2[chan].begin(), & *CC.T2T[chan].begin(), & *CC2.T2[chan].begin(), &hp1, &hp2, &hp1, &fac1, &fac1, &N, &N);
    CC2.T3[chan] = CC2.T2[chan]; // T3 = T2
    RM_dgemm(& *CC2.S4[chan].begin(), & *CC.T2T[chan].begin(), & *CC2.T2[chan].begin(), &hp1, &hp2, &hp1, &fac1, &fac5, &N, &N);
    RM_dgemm(& *CC2.S4T[chan].begin(), & *CC.T2T[chan].begin(), & *CC2.T3[chan].begin(), &hp1, &hp2, &hp1, &fac5, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.hh[chan];
    int pp = Chan.pp[chan];
    if(hh == 0 || pp == 0){ continue; }
    //S1 = T1*HHPP1
    RM_dgemm(& *CC.T1[chan].begin(), & *ME.HHPP1[chan].begin(), & *CC2.S1[chan].begin(), &hh, &hh, &pp, &fac1, &fac2, &N, &N);
    //T1 = 0.5*T1*PPPP + 0.5*HHHH*T1 + 0.25*S1*T1
    RM_dgemm(& *CC.T1[chan].begin(), & *ME.PPPP[chan].begin(), & *CC2.T1[chan].begin(), &hh, &pp, &pp, &fac3, &fac2, &N, &N);
    RM_dgemm(& *ME.HHHH[chan].begin(), & *CC.T1[chan].begin(), & *CC2.T1[chan].begin(), &hh, &pp, &hh, &fac3, &fac1, &N, &N);
    RM_dgemm(& *CC2.S1[chan].begin(), & *CC.T1[chan].begin(), & *CC2.T1[chan].begin(), &hh, &pp, &hh, &fac4, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    //S2 = HHPP2*T5
    RM_dgemm(& *ME.HHPP2[chan].begin(), & *CC.T5[chan].begin(), & *CC2.S2[chan].begin(), &p, &p, &hhp, &fac1, &fac2, &N, &N);
    //T4 = -0.5*T4*S2  &  T5 
    RM_dgemm(& *CC.T4[chan].begin(), & *CC2.S2[chan].begin(), & *CC2.T4[chan].begin(), &hhp, &p, &p, &fac6, &fac2, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    //S3 = HHPP3*T7
    RM_dgemm(& *ME.HHPP3[chan].begin(), & *CC.T7[chan].begin(), & *CC2.S3[chan].begin(), &h, &h, &hpp, &fac1, &fac2, &N, &N);
    //T6 = -0.5*T6*S3  &  T7
    RM_dgemm(& *CC.T6[chan].begin(), & *CC2.S3[chan].begin(), & *CC2.T6[chan].begin(), &hpp, &h, &h, &fac6, &fac2, &N, &N);
  }
}


CC_Eff Build_CC_Eff(const Model_Space &Space, const Input_Parameters &Parameters, CC_Matrix_Elements &CCME, CCD &CC, const Channels &Chan)
{
  CC_Eff V_Eff(Chan);
  double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac4 = -0.5, fac5 = -1.0;
  char N = 'N';//, T = 'T';
  double ME;

  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    RM_dgemm(& *CCME.HHPP2[chan].begin(), & *CC.T5[chan].begin(), & *V_Eff.V1[chan].begin(), &p, &p, &hhp, &fac4, &fac2, &N, &N);
    for(int i = p-1; i >= 0; --i){
      for(int j = i-1; j >= 0; --j){
	ME = V_Eff.V1[chan][p*j + i];
	V_Eff.V1[chan][p*j + i] = V_Eff.V1[chan][p*i + j];
	V_Eff.V1[chan][p*i + j] = ME;
      }
    }
    for(int i = 0; i < p; ++i){ V_Eff.V1[chan][p * i + i] += Space.levelsen[Chan.pvec1[chan][i]]; }
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    RM_dgemm(& *CCME.HHPP3[chan].begin(), & *CC.T7[chan].begin(), & *V_Eff.V2[chan].begin(), &h, &h, &hpp, &fac3, &fac2, &N, &N);
    for(int i = 0; i < h; ++i){ V_Eff.V2[chan][h * i + i] += Space.levelsen[Chan.hvec1[chan][i]]; }
  }
  for(int chan = 0; chan < Chan.size1; ++chan){
    int hh = Chan.hh[chan];
    int pp = Chan.pp[chan];
    if(hh == 0 || pp == 0){ continue; }
    V_Eff.V3[chan] = CCME.PPPP[chan];
    V_Eff.V4[chan] = CCME.HHHH[chan];
    RM_dgemm(& *CCME.HHPP1[chan].begin(), & *CC.T1[chan].begin(), & *V_Eff.V3[chan].begin(), &pp, &pp, &hh, &fac3, &fac1, &N, &N);
    RM_dgemm(& *CC.T1[chan].begin(), & *CCME.HHPP1[chan].begin(), & *V_Eff.V4[chan].begin(), &hh, &hh, &pp, &fac3, &fac1, &N, &N);
    for(int i = hh-1; i >= 0; --i){
      for(int j = i-1; j >= 0; --j){
	ME = V_Eff.V4[chan][hh*j + i];
	V_Eff.V4[chan][hh*j + i] = V_Eff.V4[chan][hh*i + j];
	V_Eff.V4[chan][hh*i + j] = ME;
      }
    }
  }
  for(int chan = 0; chan < Chan.size2; ++chan){
    int hp1 = Chan.hp1[chan];
    int hp2 = Chan.hp2[chan];
    V_Eff.V5[chan] = CCME.HPHP2[chan];
    if(hp1 == 0 || hp2 == 0){ continue; }
    RM_dgemm(& *CC.T2[chan].begin(), & *CCME.HHPP4T[chan].begin(), & *V_Eff.V5[chan].begin(), &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int h = Chan.h[chan];
    int hpp = Chan.hpp[chan];
    if(h == 0 || hpp == 0){ continue; }
    RM_dgemm(& *CC.T7[chan].begin(), & *CCME.HHPP3[chan].begin(), & *V_Eff.V6[chan].begin(), &hpp, &hpp, &h, &fac1, &fac2, &N, &N);
  }
  for(int chan = 0; chan < Chan.size3; ++chan){
    int p = Chan.p[chan];
    int hhp = Chan.hhp[chan];
    if(p == 0 || hhp == 0){ continue; }
    RM_dgemm(& *CC.T5[chan].begin(), & *CCME.HHPP2[chan].begin(), & *V_Eff.V7[chan].begin(), &hhp, &hhp, &p, &fac5, &fac2, &N, &N);
  }
  return V_Eff;
}



/*void EE_EOM(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan)
{
  int Nx = Parameters.Nx;
  int Ny = Parameters.Ny;
  int Nz = Parameters.Nz;
  double M = Parameters.M;
  double T = Parameters.T;

  // Nx = Nx_pp - Nx_hh, etc...
  // (Nx_hh,Nx_pp) = (Nx_hh,Nx + Nx_hh),...
  int Nx_hh, Ny_hh, Nz_hh;
  std::vector<std::vector<int> > chanvec;
  std::vector<int> vec2(2);
  for(int i = 0; i < int(Space.CART_tb1Indvec.size()); ++i){
    Nx_hh = Space.Nxmin + i;
    for(int j = 0; j < int(Space.CART_tb1Indvec[i].size()); ++j){
      Ny_hh = Space.Nymin[i] + j;
      for(int k = 0; k < int(Space.CART_tb1Indvec[i][j].size()); ++k){
	Nz_hh = Space.Nzmin[i][j] + k;
	for(int m = Space.Mmin; m <= Space.Mmin+2*(Space.Msize-1); m=m+2){
	  for(int t = Space.Tmin; t <= Space.Tmin+2*(Space.Tsize-1); t=t+2){
	    if((Nx+Nx_hh)*(Nx+Nx_hh) + (Ny+Ny_hh)*(Ny+Ny_hh) + (Nz+Nz_hh)*(Nz+Nz_hh) <= 4*Parameters.Nmax && 
	       M+m >= Space.Mmin && M+m <= Space.Mmin+2*(Space.Msize-1) && T+t >= Space.Tmin && T+t <= Space.Tmin+2*(Space.Tsize-1)){
	      vec2[0] = CART_tbInd1(Space, Nx_hh, Ny_hh, Nz_hh, m, t);
	      vec2[1] = CART_tbInd1(Space, Nx+Nx_hh, Ny+Ny_hh, Nz+Nz_hh, M+m, T+t);
	      chanvec.push_back(vec2);
	    }
	  }
	}
      }
    }
  }

  int count = 0;
  for(int i = 0; i < int(chanvec.size()); ++i){
    count += Chan.hh[chanvec[i][0]] * Chan.pp[chanvec[i][1]];
  }
  double memory = 24.0 + 8.0 * count * count / 16.0;
  std::cout <<"Matrix for Excited State: "<< Nx <<" "<< Ny <<" "<< Nz <<" "<< M <<" "<< T <<": "<< count/4.0 <<" = "<< memory/1000000.0 <<" MB"<< std::endl;

  int chansize = int(chanvec.size());
  std::vector<int> hhppvec2;
  for(int i = 0; i < chansize; ++i){
    for(int j = 0; j < Chan.hh[chanvec[i][0]]; ++j){
      if(Chan.hhvec1[chanvec[i][0]][2*j] > Chan.hhvec1[chanvec[i][0]][2*j+1]){ continue; }
      for(int k = 0; k < Chan.pp[chanvec[i][1]]; ++k){
	if(Chan.ppvec1[chanvec[i][1]][2*k] > Chan.ppvec1[chanvec[i][1]][2*k+1]){ continue; }
	hhppvec2.push_back(Chan.hhvec1[chanvec[i][0]][2*j]);
	hhppvec2.push_back(Chan.hhvec1[chanvec[i][0]][2*j+1]);
	hhppvec2.push_back(Chan.ppvec1[chanvec[i][1]][2*k]);
	hhppvec2.push_back(Chan.ppvec1[chanvec[i][1]][2*k+1]);
      }
    }
  }
   
  int Nbit = std::ceil(Chan.size3/64.0);
  int N = int(count/4.0);
  std::vector<double> Ham(int(count*count/16.0), 0.0);
  std::vector<unsigned long long> ket;

  for(int col = 0; col < N; ++col){
    ket = bitsetup(hhppvec2[4*col], hhppvec2[4*col+1], hhppvec2[4*col+2], hhppvec2[4*col+3], Nbit);
    std::cout << "! " << col << std::endl;
    #pragma omp parallel for
    for(int row = 0; row < N; ++row){
      std::vector<unsigned long long> bra = bitsetup(hhppvec2[4*row], hhppvec2[4*row+1], hhppvec2[4*row+2], hhppvec2[4*row+3], Nbit);
      std::vector<unsigned long long> tempket = ket;
      double ME, ME1, ME0;
      ME = 0.0;
      //std::cout << "bra = " << hhppvec2[4*row] << " " << hhppvec2[4*row+1] << " " << hhppvec2[4*row+2] << " " << hhppvec2[4*row+3] << std::endl;
      //std::cout << "ket = " << hhppvec2[4*col] << " " << hhppvec2[4*col+1] << " " << hhppvec2[4*col+2] << " " << hhppvec2[4*col+3] << std::endl;
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int a = 0; a < Chan.p[chan]; ++a){
	  for(int b = 0; b < Chan.p[chan]; ++b){
	    ME0 = V_Eff.V1[chan][Chan.p[chan]*a + b];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe1(bra, ket, Chan.pvec1[chan][a], Chan.pvec1[chan][b], ME0, 0);
	    ME += ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.h[chan]; ++i){ // Normal ordering gives -1
	  for(int j = 0; j < Chan.h[chan]; ++j){
	    ME0 = V_Eff.V2[chan][Chan.h[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe1(bra, ket, Chan.hvec1[chan][j], Chan.hvec1[chan][i], ME0, 0); // {i+,j}=-j,i+
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < chansize; ++chan){
	for(int pp3 = 0; pp3 < Chan.pp[chanvec[chan][1]]; ++pp3){
	  if(Chan.ppvec1[chanvec[chan][1]][2*pp3] > Chan.ppvec1[chanvec[chan][1]][2*pp3+1]){ continue; }
	  for(int pp4 = 0; pp4 < Chan.pp[chanvec[chan][1]]; ++pp4){
	    if(Chan.ppvec1[chanvec[chan][1]][2*pp4] > Chan.ppvec1[chanvec[chan][1]][2*pp4+1]){ continue; }
	    ME0 = V_Eff.V3[chanvec[chan][1]][Chan.pp[chanvec[chan][1]]*pp4 + pp3];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.ppvec1[chanvec[chan][1]][2*pp4], Chan.ppvec1[chanvec[chan][1]][2*pp4+1], // <cd|ab> reverse
			   Chan.ppvec1[chanvec[chan][1]][2*pp3+1], Chan.ppvec1[chanvec[chan][1]][2*pp3], ME0); // <pq|rs> p+,q+,s,r
	    ME += ME1;
	  }
	}
	for(int hh3 = 0; hh3 < Chan.hh[chanvec[chan][0]]; ++hh3){
	  if(Chan.hhvec1[chanvec[chan][0]][2*hh3] > Chan.hhvec1[chanvec[chan][0]][2*hh3+1]){ continue; }
	  for(int hh4 = 0; hh4 < Chan.hh[chanvec[chan][0]]; ++hh4){ // <ij|kl>{i+,j+,l,k} = <ij|kl> l,k,i+,j+
	    if(Chan.hhvec1[chanvec[chan][0]][2*hh4] > Chan.hhvec1[chanvec[chan][0]][2*hh4+1]){ continue; }
	    ME0 = V_Eff.V4[chanvec[chan][0]][Chan.hh[chanvec[chan][0]]*hh3 + hh4];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.hhvec1[chanvec[chan][0]][2*hh4+1], Chan.hhvec1[chanvec[chan][0]][2*hh4], 
			   Chan.hhvec1[chanvec[chan][0]][2*hh3], Chan.hhvec1[chanvec[chan][0]][2*hh3+1], ME0);
	    ME += ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size2; ++chan){
	for(int hp1 = 0; hp1 < Chan.hp1[chan]; ++hp1){
	  for(int hp2 = 0; hp2 < Chan.hp1[chan]; ++hp2){ // <ia|jb>{i+,a+,b,j} = - <ia|jb> a+,j,i+,b
	    ME0 = V_Eff.V5[chan][Chan.hp1[chan]*hp1 + hp2];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.hp1vec1[chan][2*hp1+1], Chan.hp1vec1[chan][2*hp1], 
			   Chan.hp1vec1[chan][2*hp2], Chan.hp1vec1[chan][2*hp2+1], ME0);
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.hpp[chan]; ++i){
	  if(Chan.hppvec1[chan][3*i+1] > Chan.hppvec1[chan][3*i+2]){ continue; }
	  for(int j = 0; j < Chan.hpp[chan]; ++j){ // <iab|cdj>{i+,a+,b+,j,d,c} = -<iab|cdj> j,a+,b+,d,c,i+
	    if(Chan.hppvec1[chan][3*j+1] > Chan.hppvec1[chan][3*j+2]){ continue; }
	    ME0 = V_Eff.V6[chan][Chan.hpp[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe3(bra, ket, Chan.hppvec1[chan][3*i], Chan.hppvec1[chan][3*i+1], Chan.hppvec1[chan][3*i+2], Chan.hppvec1[chan][3*j+2], Chan.hppvec1[chan][3*j+1], Chan.hppvec1[chan][3*j], ME0);
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.hhp[chan]; ++i){
	  if(Chan.hhpvec1[chan][3*i] > Chan.hhpvec1[chan][3*i+1]){ continue; }
	  for(int j = 0; j < Chan.hhp[chan]; ++j){ // <ija|bkl>{i+,j+,a+,l,k,b} = <ija|bkl> a+,l,k,i+,j+,b
	    if(Chan.hhpvec1[chan][3*j] > Chan.hhpvec1[chan][3*j+1]){ continue; }
	    ME0 = V_Eff.V7[chan][Chan.h[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe3(bra, ket, Chan.hhpvec1[chan][3*i+2], Chan.hhpvec1[chan][3*i+1], Chan.hhpvec1[chan][3*i], Chan.hhpvec1[chan][3*j], Chan.hhpvec1[chan][3*j+1], Chan.hhpvec1[chan][3*j+2], ME0);
	    ME += ME1;
	  }
	}
      }
      // fill column major
      Ham[N*col + row] = ME;
    }
  }

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
  }
}

void EE_EOM_HO(const Model_Space &Space, const Input_Parameters &Parameters, const CC_Eff &V_Eff, const CCD &CC, const Channels &Chan)
{
  double Par = Parameters.Par;
  double M = Parameters.M;
  double T = Parameters.T;

  std::vector<std::vector<int> > chanvec;
  std::vector<int> vec2(2);
  // Nx = Nx_pp - Nx_hh, etc...
  // (Nx_hh,Nx_pp) = (Nx_hh,Nx + Nx_hh),...
  for(int p = Space.Pmin; p < int(Space.Pmin+2*Space.HO_tb1Indvec.size()); p+=2){ //doubles
    for(int m = Space.Mmin; m <= Space.Mmin+2*(Space.Msize-1); m+=2){
      for(int t = Space.Tmin; t <= Space.Tmin+2*(Space.Tsize-1); t+=2){
	if(M+m >= Space.Mmin && M+m <= Space.Mmin+2*(Space.Msize-1) && T+t >= Space.Tmin 
	   && T+t <= Space.Tmin+2*(Space.Tsize-1) && Par*p >= Space.Pmin && Par*p <= Space.Pmin+2*Space.HO_tb1Indvec.size()){
	  vec2[0] = HO_tbInd1(Space, p, m, t);
	  vec2[1] = HO_tbInd1(Space, Par*p, M+m, T+t);
	  chanvec.push_back(vec2);
	}
      }
    }
  }

  int count = 0;
  for(int i = 0; i < int(chanvec.size()); ++i){
    count += Chan.hh[chanvec[i][0]] * Chan.pp[chanvec[i][1]];
  }
  double memory = 24.0 + 8.0 * count * count / 16.0;
  std::cout <<"Matrix for Excited State: "<< Par <<" "<< M <<" "<< T <<": "<< count/4.0 <<" = "<< memory/1000000.0 <<" MB"<< std::endl;

  int chansize = int(chanvec.size());
  std::vector<int> hpvec2;
  std::vector<int> hhppvec2;

  std::cout << "2h2p states" << std::endl;
  for(int i = 0; i < chansize; ++i){
    for(int j = 0; j < Chan.hh[chanvec[i][0]]; ++j){
      if(Chan.hhvec1[chanvec[i][0]][2*j] > Chan.hhvec1[chanvec[i][0]][2*j+1]){ continue; }
      for(int k = 0; k < Chan.pp[chanvec[i][1]]; ++k){
	if(Chan.ppvec1[chanvec[i][1]][2*k] > Chan.ppvec1[chanvec[i][1]][2*k+1]){ continue; }
	hhppvec2.push_back(Chan.hhvec1[chanvec[i][0]][2*j]);
	hhppvec2.push_back(Chan.hhvec1[chanvec[i][0]][2*j+1]);
	hhppvec2.push_back(Chan.ppvec1[chanvec[i][1]][2*k]);
	hhppvec2.push_back(Chan.ppvec1[chanvec[i][1]][2*k+1]);
      }
    }
  }
    
  int Nbit = std::ceil(Chan.size3/64.0);
  int N = int(count/4.0);
  std::vector<double> Ham(int(count*count/16.0), 0.0);
  std::vector<unsigned long long> ket;

  for(int col = 0; col < N; ++col){
    ket = bitsetup(hhppvec2[4*col], hhppvec2[4*col+1], hhppvec2[4*col+2], hhppvec2[4*col+3], Nbit);
    #pragma omp parallel for
    for(int row = 0; row < N; ++row){
      std::vector<unsigned long long> bra = bitsetup(hhppvec2[4*row], hhppvec2[4*row+1], hhppvec2[4*row+2], hhppvec2[4*row+3], Nbit);
      std::vector<unsigned long long> tempket = ket;
      double ME, ME1, ME0;
      ME = 0.0;
      //std::cout << "bra = " << hhppvec2[4*row] << " " << hhppvec2[4*row+1] << " " << hhppvec2[4*row+2] << " " << hhppvec2[4*row+3] << std::endl;
      //std::cout << "ket = " << hhppvec2[4*col] << " " << hhppvec2[4*col+1] << " " << hhppvec2[4*col+2] << " " << hhppvec2[4*col+3] << std::endl;
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int a = 0; a < Chan.p[chan]; ++a){
	  for(int b = 0; b < Chan.p[chan]; ++b){
	    ME0 = V_Eff.V1[chan][Chan.p[chan]*a + b];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe1(bra, ket, Chan.pvec1[chan][a], Chan.pvec1[chan][b], ME0, 0);
	    ME += ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.h[chan]; ++i){ // Normal ordering gives -1
	  for(int j = 0; j < Chan.h[chan]; ++j){
	    ME0 = V_Eff.V2[chan][Chan.h[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe1(bra, ket, Chan.hvec1[chan][j], Chan.hvec1[chan][i], ME0, 0); // {i+,j}=-j,i+
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < chansize; ++chan){
	for(int pp3 = 0; pp3 < Chan.pp[chanvec[chan][1]]; ++pp3){
	  if(Chan.ppvec1[chanvec[chan][1]][2*pp3] > Chan.ppvec1[chanvec[chan][1]][2*pp3+1]){ continue; }
	  for(int pp4 = 0; pp4 < Chan.pp[chanvec[chan][1]]; ++pp4){
	    if(Chan.ppvec1[chanvec[chan][1]][2*pp4] > Chan.ppvec1[chanvec[chan][1]][2*pp4+1]){ continue; }
	    ME0 = V_Eff.V3[chanvec[chan][1]][Chan.pp[chanvec[chan][1]]*pp4 + pp3];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.ppvec1[chanvec[chan][1]][2*pp4], Chan.ppvec1[chanvec[chan][1]][2*pp4+1], // <cd|ab> reverse
			   Chan.ppvec1[chanvec[chan][1]][2*pp3+1], Chan.ppvec1[chanvec[chan][1]][2*pp3], ME0); // <pq|rs> p+,q+,s,r
	    ME += ME1;
	  }
	}
	for(int hh3 = 0; hh3 < Chan.hh[chanvec[chan][0]]; ++hh3){
	  if(Chan.hhvec1[chanvec[chan][0]][2*hh3] > Chan.hhvec1[chanvec[chan][0]][2*hh3+1]){ continue; }
	  for(int hh4 = 0; hh4 < Chan.hh[chanvec[chan][0]]; ++hh4){ // <ij|kl>{i+,j+,l,k} = <ij|kl> l,k,i+,j+
	    if(Chan.hhvec1[chanvec[chan][0]][2*hh4] > Chan.hhvec1[chanvec[chan][0]][2*hh4+1]){ continue; }
	    ME0 = V_Eff.V4[chanvec[chan][0]][Chan.hh[chanvec[chan][0]]*hh3 + hh4];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.hhvec1[chanvec[chan][0]][2*hh4+1], Chan.hhvec1[chanvec[chan][0]][2*hh4], 
			   Chan.hhvec1[chanvec[chan][0]][2*hh3], Chan.hhvec1[chanvec[chan][0]][2*hh3+1], ME0);
	    ME += ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size2; ++chan){
	for(int hp1 = 0; hp1 < Chan.hp1[chan]; ++hp1){
	  for(int hp2 = 0; hp2 < Chan.hp1[chan]; ++hp2){ // <ia|jb>{i+,a+,b,j} = - <ia|jb> a+,j,i+,b
	    ME0 = V_Eff.V5[chan][Chan.hp1[chan]*hp1 + hp2];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe2(bra, ket, Chan.hp1vec1[chan][2*hp1+1], Chan.hp1vec1[chan][2*hp1], 
			   Chan.hp1vec1[chan][2*hp2], Chan.hp1vec1[chan][2*hp2+1], ME0);
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.hpp[chan]; ++i){
	  if(Chan.hppvec1[chan][3*i+1] > Chan.hppvec1[chan][3*i+2]){ continue; }
	  for(int j = 0; j < Chan.hpp[chan]; ++j){ // <iab|cdj>{i+,a+,b+,j,d,c} = -<iab|cdj> j,a+,b+,d,c,i+
	    if(Chan.hppvec1[chan][3*j+1] > Chan.hppvec1[chan][3*j+2]){ continue; }
	    ME0 = V_Eff.V6[chan][Chan.hpp[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe3(bra, ket, Chan.hppvec1[chan][3*i], Chan.hppvec1[chan][3*i+1], Chan.hppvec1[chan][3*i+2], Chan.hppvec1[chan][3*j+2], Chan.hppvec1[chan][3*j+1], Chan.hppvec1[chan][3*j], ME0);
	    ME -= ME1;
	  }
	}
      }
      for(int chan = 0; chan < Chan.size3; ++chan){
	for(int i = 0; i < Chan.hhp[chan]; ++i){
	  if(Chan.hhpvec1[chan][3*i] > Chan.hhpvec1[chan][3*i+1]){ continue; }
	  for(int j = 0; j < Chan.hhp[chan]; ++j){ // <ija|bkl>{i+,j+,a+,l,k,b} = <ija|bkl> a+,l,k,i+,j+,b
	    if(Chan.hhpvec1[chan][3*j] > Chan.hhpvec1[chan][3*j+1]){ continue; }
	    ME0 = V_Eff.V7[chan][Chan.h[chan]*i + j];
	    if(ME0 == 0){ continue; }
	    ME1 = matrixe3(bra, ket, Chan.hhpvec1[chan][3*i+2], Chan.hhpvec1[chan][3*i+1], Chan.hhpvec1[chan][3*i], Chan.hhpvec1[chan][3*j], Chan.hhpvec1[chan][3*j+1], Chan.hhpvec1[chan][3*j+2], ME0);
	    ME += ME1;
	  }
	}
      }
      // fill column major
      Ham[N*col + row] = ME;
    }
  }

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
  }
}

std::vector<unsigned long long> bitsetup(const int &i, const int &j, const int &a, const int &b, const int &Nbit)
{
  std::vector<unsigned long long> bconfig(Nbit, 0);
  unsigned long long one = 1;
  bconfig[std::floor(i/64.0)] += (one << (i%64));
  bconfig[std::floor(j/64.0)] += (one << (j%64));
  bconfig[std::floor(a/64.0)] += (one << (a%64));
  bconfig[std::floor(b/64.0)] += (one << (b%64));
  return bconfig;
}


double matrixe1(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const double &ME, const int &test)
{
  unsigned long long one = 1, zero = 0;
  unsigned long long pbit = one << (p%64);
  unsigned long long qbit = one << (q%64);
  unsigned long long comp;
  std::vector<unsigned long long> tempket = ket;
  int p1 = std::floor(p/64);
  int q1 = std::floor(q/64);
  double phase = 1.0;
  int bcount;

  if((bra[p1] & pbit) == 0 || (tempket[q1] & qbit) == 0){ return 0.0; }
  comp = tempket[q1] & ~(~zero << (q%64));
  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
  for(int i = q1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
  phase = phase*pow(-1.0, bcount); tempket[q1] = tempket[q1]^qbit;
  comp = tempket[p1] & ~(~zero << (p%64));
  for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
  for(int i = p1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
  phase = phase*pow(-1.0, bcount); tempket[p1] = tempket[p1]^pbit;
  for(int i = 0; i < int(bra.size()); ++i){ if(bra[i] != tempket[i]){ return 0.0; } }
  return phase * ME;
}


double matrixe2(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &s, const int &r, const double &ME)
{
  unsigned long long one = 1, zero = 0;
  unsigned long long pbit = one << (p%64);
  unsigned long long qbit = one << (q%64);
  unsigned long long sbit = one << (s%64);
  unsigned long long rbit = one << (r%64);
  unsigned long long comp;
  std::vector<unsigned long long> tempket = ket;
  int p1 = std::floor(p/64);
  int q1 = std::floor(q/64);
  int s1 = std::floor(s/64);
  int r1 = std::floor(r/64);
  double phase = 1.0;
  int bcount;

  if((sbit & tempket[s1]) != 0 && (rbit & tempket[r1]) != 0 && (qbit & bra[q1]) != 0 && (pbit & bra[p1]) != 0 && 
     ((pbit & tempket[p1]) == 0 || (p == s || p == r)) && ((qbit & tempket[q1]) == 0 || (q == s || q == r))){
    comp = tempket[r1] & ~(~zero << (r%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = r1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[r1] = tempket[r1]^rbit;
    comp = tempket[s1] & ~(~zero << (s%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = s1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[s1] = tempket[s1]^sbit;
    comp = tempket[q1] & ~(~zero << (q%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = q1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[q1] = tempket[q1]^qbit;
    comp = tempket[p1] & ~(~zero << (p%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = p1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[p1] = tempket[p1]^pbit;
    for(int i = 0; i < int(bra.size()); ++i){ if(bra[i] != tempket[i]){ return 0.0; } }
    return phase * ME;
  }
  else{ return 0.0; }
}

double matrixe3(const std::vector<unsigned long long> &bra, const std::vector<unsigned long long> &ket, const int &p, const int &q, const int &r, const int &t, const int &u, const int &v, const double &ME)
{
  unsigned long long one = 1, zero = 0;
  unsigned long long pbit = one << (p%64);
  unsigned long long qbit = one << (q%64);
  unsigned long long rbit = one << (r%64);
  unsigned long long tbit = one << (t%64);
  unsigned long long ubit = one << (u%64);
  unsigned long long vbit = one << (v%64);
  unsigned long long comp;
  std::vector<unsigned long long> tempket = ket;
  int p1 = std::floor(p/64);
  int q1 = std::floor(q/64);
  int r1 = std::floor(r/64);
  int t1 = std::floor(t/64);
  int u1 = std::floor(u/64);
  int v1 = std::floor(v/64);
  double phase = 1.0;
  int bcount;

  if((tbit & tempket[t1]) != 0 && (ubit & tempket[u1]) != 0 && (vbit & tempket[v1]) != 0 && 
     (pbit & bra[p1]) != 0 && (qbit & bra[q1]) != 0 && (rbit & bra[r1]) != 0 && 
     ((pbit & tempket[p1]) == 0 || (p == t || p == u || p == v)) && 
     ((qbit & tempket[q1]) == 0 || (q == t || q == u || q == v)) &&
     ((rbit & tempket[r1]) == 0 || (r == t || r == u || r == v))){
    comp = tempket[v1] & ~(~zero << (v%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = v1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[v1] = tempket[v1]^vbit;
    comp = tempket[u1] & ~(~zero << (u%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = u1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[u1] = tempket[u1]^ubit;
    comp = tempket[t1] & ~(~zero << (t%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = t1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[t1] = tempket[t1]^tbit;
    comp = tempket[r1] & ~(~zero << (r%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = r1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[r1] = tempket[r1]^rbit;
    comp = tempket[q1] & ~(~zero << (q%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = q1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[q1] = tempket[q1]^qbit;
    comp = tempket[p1] & ~(~zero << (p%64));
    for(bcount = 0; comp; ++bcount){ comp ^= comp & -comp; }
    for(int i = p1-1; i >= 0; --i){ std::bitset<64> bket (tempket[i]); bcount += bket.count(); }
    phase = phase*pow(-1.0, bcount); tempket[p1] = tempket[p1]^pbit;
    for(int i = 0; i < int(bra.size()); ++i){ if(bra[i] != tempket[i]){ return 0.0; } }
    return phase * ME;
  }
  else{ return 0.0; }
}*/
