#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

struct Input_Parameters;

struct Model_Space;
struct State;
struct Channels;

struct two_body;
struct three_body;

void plus(State &S, State &S1, State &S2);
void minus(State &S, State &S1, State &S2);
bool equal(State &S1, State &S2);

int ChanInd_1b(std::string &basis, Model_Space &Space, State &State);
int ChanInd_2b_dir(std::string &basis, Model_Space &Space, State &State);
int ChanInd_2b_cross(std::string &basis, Model_Space &Space, State &State);

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
    energy = -1000.0;
    type = "none";
  };

  void minimize(){
    t = -1000;
    m = -1000;
    nx = -1000;
    ny = -1000;
    nz = -1000;
    ml = -1000;
    n = -1000;
    j = -1000;
    par = -1000;
  };

  void maximize(){
    t = 1000;
    m = 1000;
    nx = 1000;
    ny = 1000;
    nz = 1000;
    ml = 1000;
    n = 1000;
    j = 1000;
    par = 1000;
  };

  void divide_spins(){
    t /= 2;
    m /= 2;
    j /= 2;
  };

  void add_one(){
    ++t;
    ++m;
    ++nx;
    ++ny;
    ++nz;
    ++ml;
    ++n;
    ++j;
  };
};

struct one_body{
  int v1;
};

struct two_body{
  int v1;
  int v2;
};

struct three_body{
  int v1;
  int v2;
  int v3;
};

//Structure for holding all model space info
struct Model_Space{
  int num_p; //number of proton orbits
  int num_n; //number of neutron orbits
  int num_par; //number of particle orbits
  int num_hol; //number of hole orbits
  int num_states; //number of total orbits
  int num_jstates; //number of orbits in j-scheme

  State *qnums;

  State qmins;
  State qmaxs;
  State qsizes;
  State qmins1;
  State qmaxs1;
  State qsizes1;
  State qmins2;
  State qmaxs2;
  State qsizes2;
  State qsizes0;

  int **shellsm; // list of m-states for each j-orbit

  int Nmax;
  int nmax;

  std::unordered_map<int,int> map_state;

  int size_1b;
  int size_2b_dir;
  int size_2b_cross;
  
  Model_Space(){};
  void Determine_Shells(Input_Parameters &Parameters);
  void Setup_Maps(Input_Parameters &Parameters);
  void delete_struct(Input_Parameters &Parameters);
  int hash2(int &p, int &q, int &j);
  int hash3(int &p, int &q, int &r, int &j);
  int ind_state(std::string &basis, State &State);
  int ind_1b(std::string &basis, State &State);
  int ind_2b_dir(std::string &basis, State &State);
  int ind_2b_cross(std::string &basis, State &State);
};

//Structure for holding channel information
struct Channels{
  int size1;
  int size2;
  int size3;

  State *qnums1;
  State *qnums2;
  State *qnums3;

  int *nh;
  one_body *h_vec;
  int *h_index;
  std::unordered_map<int,int> *h_map;
  int *np;
  one_body *p_vec;
  int *p_index;
  std::unordered_map<int,int> *p_map;
  
  int *nhh;
  two_body *hh_vec;
  int *hh_index;
  std::unordered_map<int,int> *hh_map;
  int *npp;
  two_body *pp_vec;
  int *pp_index;
  std::unordered_map<int,int> *pp_map;
  int *nhp1;
  two_body *hp1_vec;
  int *hp1_index;
  std::unordered_map<int,int> *hp1_map;
  int *nph1;
  two_body *ph1_vec;
  int *ph1_index;
  std::unordered_map<int,int> *ph1_map;

  int *nhhp;
  three_body *hhp_vec;
  int *hhp_index;
  std::unordered_map<int,int> *hhp_map;
  int *npph;
  three_body *pph_vec;
  int *pph_index;
  std::unordered_map<int,int> *pph_map;

  // For CCSD
  int *nhp;
  two_body *hp_vec;
  int *hp_index;
  std::unordered_map<int,int> *hp_map;
  int *nhh1;
  two_body *hh1_vec;
  int *hh1_index;
  std::unordered_map<int,int> *hh1_map;
  int *npp1;
  two_body *pp1_vec;
  int *pp1_index;
  std::unordered_map<int,int> *pp1_map;
  int *nhph;
  three_body *hph_vec;
  int *hph_index;
  std::unordered_map<int,int> *hph_map;
  int *nhpp;
  three_body *hpp_vec;
  int *hpp_index;
  std::unordered_map<int,int> *hpp_map;
  int *nhhh;
  three_body *hhh_vec;
  int *hhh_index;
  std::unordered_map<int,int> *hhh_map;
  int *nppp;
  three_body *ppp_vec;
  int *ppp_index;
  std::unordered_map<int,int> *ppp_map;
  int ind0; // index of i-i cross channel for singles
  //////////////////////////

  Channels(){};
  Channels(Input_Parameters &Parameters, Model_Space &Space);
  void delete_struct(Input_Parameters &Parameters);
  one_body h_state(int &chan3, int &ind3);
  one_body p_state(int &chan3, int &ind3);
  two_body hh_state(int &chan1, int &ind1);
  two_body pp_state(int &chan1, int &ind1);
  two_body hp1_state(int &chan2, int &ind2);
  two_body ph1_state(int &chan2, int &ind2);
  three_body pph_state(int &chan3, int &ind3);
  three_body hhp_state(int &chan3, int &ind3);
  two_body hp_state(int &chan1, int &ind1);
  two_body hh1_state(int &chan2, int &ind2);
  two_body pp1_state(int &chan2, int &ind2);
  three_body hph_state(int &chan3, int &ind3);
  three_body hpp_state(int &chan3, int &ind3);
  three_body hhh_state(int &chan3, int &ind3);
  three_body ppp_state(int &chan3, int &ind3);
};

#endif
