#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

struct Input_Parameters;

struct Model_Space;
struct State;
struct Channels;

struct one_body;
struct two_body;
struct three_body;
struct four_body;

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

void p_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &p_total, int *&np, int *&p_index, std::unordered_map<int,int> *&p_map, one_body *&p_vec, int &chan3size, State *qnums3, int type);
void pq_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pq_total, int *&npq, int *&pq_index, std::unordered_map<int,int> *&pq_map, two_body *&pq_vec, int &chansize, State *qnums, int &chan3size, State *qnums3, int *np, one_body *p_vec, int *p_index, int *nq, int *q_index, one_body *q_vec);
void pq1_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pq_total, int *&npq, int *&pq_index, std::unordered_map<int,int> *&pq_map, two_body *&pq_vec, int &chan2size, State *qnums2, int &chan3size, State *qnums3, int *np, one_body *p_vec, int *p_index, int *nq, int *q_index, one_body *q_vec);
void pqr_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pqr_total, int *&npqr, int *&pqr_index, std::unordered_map<int,int> *&pqr_map, three_body *&pqr_vec, State *&pqr_j, int &chan3size, State *qnums3, int *npq, two_body *pq_vec, int *pq_index, int *nr, int *r_index, one_body *r_vec);

void Map_4_count_1(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb);
void Map_4_count_2(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb);
void Map_4_count_3(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb);
void Map_4_count_4(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb);
void Map_4_count_5678(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count);
void Map_4_1(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_3(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_4(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_5(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_6(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_7(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_4_8(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J);
void Map_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, two_body &tb, int &J);
void direct_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &q, int &r, int &s, int &jmin1, State &tb1);
void cross_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &s, int &r, int &q, int &jmin2, State &tb2);

struct State{
  int t; // x2
  int m; // x2
  int nx;
  int ny;
  int nz;
  int ml;
  int n;
  int j; // x2
  int l;
  int par; // -1,+1
  double energy;
  int type;
 
  State(){
    t = -1;
    m = 0;
    nx = 0;
    ny = 0;
    nz = 0;
    ml = 0;
    n = 0;
    j = 0;
    l = 0;
    par = 1;
    energy = -1000.0;
    type = 100;
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
    l = -1000;
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
    l = 1000;
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
    ++l;
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

struct four_body{
  int v1;
  int v2;
  int v3;
  int v4;
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

  // for jm scheme
  int *shellsj; // list of j for each m-states
  int *shellsnum; // number of m-states for each j-shell
  int **shellsm; // list of m-states for each j-shell

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
  int hash2(int &p, int &q, int j);
  int hash3(int &p, int &q, int &r, int j);
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
  State *hhp_j;
  std::unordered_map<int,int> *hhp_map;
  int *npph;
  three_body *pph_vec;
  int *pph_index;
  State *pph_j;
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
  State *hph_j;
  std::unordered_map<int,int> *hph_map;
  int *nhpp;
  three_body *hpp_vec;
  int *hpp_index;
  State *hpp_j;
  std::unordered_map<int,int> *hpp_map;
  int *nhhh;
  three_body *hhh_vec;
  int *hhh_index;
  State *hhh_j;
  std::unordered_map<int,int> *hhh_map;
  int *nppp;
  three_body *ppp_vec;
  int *ppp_index;
  State *ppp_j;
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
