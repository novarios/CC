#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

struct Input_Parameters;

struct Model_Space;
struct State;
struct Channels;

void plus(State &S, const State &S1, const State &S2);
void minus(State &S, const State &S1, const State &S2);
bool equal(const State &S1, const State &S2);

int ChanInd_1b(const std::string &basis, const Model_Space &Space, const State &State);
int ChanInd_2b_dir(const std::string &basis, const Model_Space &Space, const State &State);
int ChanInd_2b_cross(const std::string &basis, const Model_Space &Space, const State &State);

void Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space);

void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);
void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space);

struct State{
  int t; // x2
  int m; // x2
  int nx;
  int ny;
  int nz;
  int ml;
  int n;
  int j; // x2
  int par; // -1,+1
  double energy;
  std::string type;
 
  State(){
    t = 0;
    m = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    ml = 0;
    n = 0;
    j = 0;
    par = 1;
    energy = -100;
    type = "none";
  };
};

//Structure for holding all model space info
struct Model_Space{
  int indp; //number of proton orbits
  int indn; //number of neutron orbits
  int indpar; //number of particle orbits
  int indhol; //number of hole orbits
  int indtot; //number of total orbits
  int indtotj;

  State *qnums;
  State qmins;
  State qmaxs;
  State qsizes;
  State qsizes0;
  int **shellsm; // for j

  int Nmax;
  int nmax;
  std::unordered_map<int,int> map_1b;
  std::unordered_map<int,int> map_2b_dir;
  std::unordered_map<int,int> map_2b_cross;
  int *map_2b;
  int size_2b;
  
  Model_Space(){};
  void delete_struct(Input_Parameters &Parameters);
};

//Structure for holding channel information
struct Channels{
  int size1;
  int size2;
  int size3;

  State *qnums1;
  State *qnums2;
  State *qnums3;
  
  int *indvec;

  int *nhh;
  int *npp;
  int *nhp;
  int *nhp1;
  int *nhp2;
  int *nh;
  int *np;
  int *nhhp;
  int *nhpp;
  int *nhhp1;
  int *nhpp1;
  int *nhh1;
  int *npp1;
  int *nhhh;
  int *nppp;

  int **hhvec;
  int **ppvec;
  int **hpvec;
  int **hp1vec;
  int **hp2vec;
  int **pvec;
  int **hvec;
  int **hhpvec;
  int **hppvec;
  int **hhp1vec;
  int **hpp1vec;
  int **hh1vec;
  int **pp1vec;
  int **hhhvec;
  int **pppvec;

  std::unordered_map<int,int> *hh_map;
  std::unordered_map<int,int> *pp_map;
  std::unordered_map<int,int> *hp_map;
  std::unordered_map<int,int> *hp1_map;
  std::unordered_map<int,int> *hp2_map;
  std::unordered_map<int,int> *p_map;
  std::unordered_map<int,int> *h_map;
  std::unordered_map<int,int> *hhp_map;
  std::unordered_map<int,int> *hpp_map;
  std::unordered_map<int,int> *hhp1_map;
  std::unordered_map<int,int> *hpp1_map;
  std::unordered_map<int,int> *hh1_map;
  std::unordered_map<int,int> *pp1_map;
  std::unordered_map<int,int> *hhh_map;
  std::unordered_map<int,int> *ppp_map;

  int ind0; // index of i-i cross channel for singles
  Channels(){};
  Channels(const Input_Parameters &Parameters, const Model_Space &Space);
  void delete_struct();
};

#endif
