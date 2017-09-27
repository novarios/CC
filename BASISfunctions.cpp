#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

void plus(State &S, State &S1, State &S2){
  S.t = S1.t + S2.t;
  S.m = S1.m + S2.m;
  S.nx = S1.nx + S2.nx;
  S.ny = S1.ny + S2.ny;
  S.nz = S1.nz + S2.nz;
  S.ml = S1.ml + S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

void minus(State &S, State &S1, State &S2){
  S.t = S1.t - S2.t;
  S.m = S1.m - S2.m;
  S.nx = S1.nx - S2.nx;
  S.ny = S1.ny - S2.ny;
  S.nz = S1.nz - S2.nz;
  S.ml = S1.ml - S2.ml;
  S.par = S1.par * S2.par;
  S.j = S1.j + S2.j;
}

bool equal(State &S1, State &S2){
  return (S1.t == S2.t &&
	  S1.m == S2.m &&
	  S1.nx == S2.nx &&
	  S1.ny == S2.ny &&
	  S1.nz == S2.nz &&
	  S1.ml == S2.ml &&
	  S1.par == S2.par &&
	  S1.j == S2.j);
}

one_body Channels::h_state(int &chan3, int &ind3)
{
  return h_vec[h_index[chan3] + ind3];
}

one_body Channels::p_state(int &chan3, int &ind3)
{
  return p_vec[p_index[chan3] + ind3];
}

two_body Channels::hh_state(int &chan1, int &ind1)
{
  return hh_vec[hh_index[chan1] + ind1];
}

two_body Channels::pp_state(int &chan1, int &ind1)
{
  return pp_vec[pp_index[chan1] + ind1];
}

two_body Channels::hp1_state(int &chan2, int &ind2)
{
  return hp1_vec[hp1_index[chan2] + ind2];
}

two_body Channels::ph1_state(int &chan2, int &ind2)
{
  return ph1_vec[ph1_index[chan2] + ind2];
}

three_body Channels::pph_state(int &chan3, int &ind3)
{
  return pph_vec[pph_index[chan3] + ind3];
}

three_body Channels::hhp_state(int &chan3, int &ind3)
{
  return hhp_vec[hhp_index[chan3] + ind3];
}

two_body Channels::hp_state(int &chan1, int &ind1)
{
  return hp_vec[hp_index[chan1] + ind1];
}

two_body Channels::hh1_state(int &chan2, int &ind2)
{
  return hh1_vec[hh1_index[chan2] + ind2];
}

two_body Channels::pp1_state(int &chan2, int &ind2)
{
  return pp1_vec[pp1_index[chan2] + ind2];
}

three_body Channels::hph_state(int &chan3, int &ind3)
{
  return hph_vec[hph_index[chan3] + ind3];
}

three_body Channels::hpp_state(int &chan3, int &ind3)
{
  return hpp_vec[hpp_index[chan3] + ind3];
}

three_body Channels::hhh_state(int &chan3, int &ind3)
{
  return hhh_vec[hhh_index[chan3] + ind3];
}

three_body Channels::ppp_state(int &chan3, int &ind3)
{
  return ppp_vec[ppp_index[chan3] + ind3];
}

//   Function to return Hash index for 2 indices
int Model_Space::hash2(int &p, int &q, int j){
  return (j * std::pow(num_states, 2) / 2) + (num_states * p + q);
}

//   Function to return Hash index for 3 indices
int Model_Space::hash3(int &p, int &q, int &r, int j){
  return (j * std::pow(num_states, 3) / 2) + (num_states * num_states * p + num_states * q + r);
}

int Model_Space::ind_state(std::string &basis, State &state)
{
  if(basis == "infinite"){
    return (state.nx - qmins.nx) * qsizes.ny * qsizes.nz * qsizes.m * qsizes.t
      + (state.ny - qmins.ny) * qsizes.nz * qsizes.m * qsizes.t
      + (state.nz - qmins.nz) * qsizes.m * qsizes.t
      + int(0.5 * (state.m - qmins.m)) * qsizes.t
      + int(0.5 * (state.t - qmins.t));
  }
  else if(basis == "finite_M"){
    return (state.n - qmins.n) * qsizes.par * qsizes.m * qsizes.t
      + int(0.5 * (state.par - qmins.par)) * qsizes.m * qsizes.t
      + int(0.5 * (state.m - qmins.m)) * qsizes.t
      + int(0.5 * (state.t - qmins.t));
  }
  else if(basis == "finite_HO"){ // for tz = -1 only;
    return (state.n - qmins.n) * qsizes.ml * qsizes.m
      + (state.ml - qmins.ml) * qsizes.m
      + int(0.5 * (state.m - qmins.m));
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return (state.n - qmins.n) * qsizes.par * qsizes.j * qsizes.t
      + int(0.5 * (state.par - qmins.par)) * qsizes.j * qsizes.t
      + int(0.5 * (state.j - qmins.j)) * qsizes.t
      + int(0.5 * (state.t - qmins.t));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int Model_Space::ind_1b(std::string &basis, State &State)
{
  if(basis == "infinite"){
    return (State.nx - qmins.nx) * qsizes.ny * qsizes.nz * qsizes.m * qsizes.t
      + (State.ny - qmins.ny) * qsizes.nz * qsizes.m * qsizes.t
      + (State.nz - qmins.nz) * qsizes.m * qsizes.t
      + int(0.5 * (State.m - qmins.m)) * qsizes.t
      + int(0.5 * (State.t - qmins.t));
  }
  else if(basis == "finite_M"){
    return int(0.5 * (State.par - qmins.par)) * qsizes.m * qsizes.t
      + int(0.5 * (State.m - qmins.m)) * qsizes.t
      + int(0.5 * (State.t - qmins.t));
  }
  else if(basis == "finite_HO"){ // for tz = -1 only;
    return (State.ml - qmins.ml) * qsizes.m
      + int(0.5 * (State.m - qmins.m));
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return int(0.5 * (State.par - qmins.par)) * qsizes.j * qsizes.t
      + int(0.5 * (State.j - qmins.j)) * qsizes.t
      + int(0.5 * (State.t - qmins.t));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int Model_Space::ind_2b_dir(std::string &basis, State &State)
{
  if(basis == "infinite"){
    return (State.nx - qmins1.nx) * qsizes1.ny * qsizes1.nz * qsizes1.m * qsizes1.t
      + (State.ny - qmins1.ny) * qsizes1.nz * qsizes1.m * qsizes1.t
      + (State.nz - qmins1.nz) * qsizes1.m * qsizes1.t
      + int(0.5 * (State.m - qmins1.m)) * qsizes1.t
      + int(0.5 * (State.t - qmins1.t));
  }
  else if(basis == "finite_M"){
    return int(0.5 * (State.par - qmins1.par)) * qsizes1.m * qsizes1.t
      + int(0.5 * (State.m - qmins1.m)) * qsizes1.t
      + int(0.5 * (State.t - qmins1.t));
  }
  else if(basis == "finite_HO"){ // for tz = -1 only;
    return (State.ml - qmins1.ml) * qsizes1.m
      + int(0.5 * (State.m - qmins1.m));
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return int(0.5 * (State.par - qmins1.par)) * qsizes1.j * qsizes1.t
      + int(0.5 * (State.j - qmins1.j)) * qsizes1.t
      + int(0.5 * (State.t - qmins1.t));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

int Model_Space::ind_2b_cross(std::string &basis, State &State)
{
  if(basis == "infinite"){
    return (State.nx - qmins2.nx) * qsizes2.ny * qsizes2.nz * qsizes2.m * qsizes2.t
      + (State.ny - qmins2.ny) * qsizes2.nz * qsizes2.m * qsizes2.t
      + (State.nz - qmins2.nz) * qsizes2.m * qsizes2.t
      + int(0.5 * (State.m - qmins2.m)) * qsizes2.t
      + int(0.5 * (State.t - qmins2.t));
  }
  else if(basis == "finite_M"){
    return int(0.5 * (State.par - qmins2.par)) * qsizes2.m * qsizes2.t
      + int(0.5 * (State.m - qmins2.m)) * qsizes2.t
      + int(0.5 * (State.t - qmins2.t));
  }
  else if(basis == "finite_HO"){ // for tz = -1 only;
    return (State.ml - qmins2.ml) * qsizes2.m
      + int(0.5 * (State.m - qmins2.m));
  }
  else if(basis == "finite_J" || basis == "finite_JM"){
    return int(0.5 * (State.par - qmins2.par)) * qsizes2.j * qsizes2.t
      + int(0.5 * (State.j - qmins2.j)) * qsizes2.t
      + int(0.5 * (State.t - qmins2.t));
  }
  else{ std::cerr << "Invalid Basis Set!" << std::endl; exit(1); }
}

void Model_Space::delete_struct(Input_Parameters &Parameters)
{
  delete[] qnums;
  if(num_jstates != 0){
    delete[] shellsnum;
    for(int i = 0; i < num_jstates; ++i){ delete[] shellsm[i]; }
    delete[] shellsm;
    delete[] shellsj;
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
  int key;
  int A;

  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

  getline(splevels, phline);
  std::stringstream(phline) >> Parameters.ho_length >> Parameters.ho_energy;
  getline(splevels, phline);
  std::stringstream(phline) >> Space.num_states;

  // allocate memory for quntum numbers for each state
  Space.qnums = new State[Space.num_states];

  // initialize mins and maxs for each quantum number
  Space.qmins.maximize();
  Space.qmins1.maximize();
  Space.qmins2.maximize();
  Space.qmaxs.minimize();
  Space.qmaxs1.minimize();
  Space.qmaxs2.minimize();
  Space.num_jstates = 0;

  // States must be ordered by energy from low to high
  for(int i = 0; i < Space.num_states; ++i){
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
    }
    else if(Parameters.basis == "finite_HO"){ // ind, n, l, lz, sz, tz
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> Space.qnums[i].ml >> Space.qnums[i].m >> energy;
      if(Space.qnums[i].ml < Space.qmins.ml){ Space.qmins.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].ml > Space.qmaxs.ml){ Space.qmaxs.ml = Space.qnums[i].ml; }
      if(Space.qnums[i].m < Space.qmins.m){ Space.qmins.m = Space.qnums[i].m; }
      if(Space.qnums[i].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[i].m; }
      if(Space.qnums[i].t < Space.qmins.t){ Space.qmins.t = Space.qnums[i].t; }
      if(Space.qnums[i].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[i].t; }
    }
    else if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){ // ind, n, l, j, tz, l2n
      std::stringstream(phline) >> ind >> Space.qnums[i].n >> Space.qnums[i].l >> Space.qnums[i].j >> Space.qnums[i].t >> ind >> energy;
      Space.qnums[i].par = -2*(Space.qnums[i].l%2) + 1;
      if(Space.qnums[i].n < Space.qmins.n){ Space.qmins.n = Space.qnums[i].n; }
      if(Space.qnums[i].n > Space.qmaxs.n){ Space.qmaxs.n = Space.qnums[i].n; }
      if(Space.qnums[i].l < Space.qmins.l){ Space.qmins.l = Space.qnums[i].l; }
      if(Space.qnums[i].l > Space.qmaxs.l){ Space.qmaxs.l = Space.qnums[i].l; }
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

  // Deterime Shell Structure
  Space.Determine_Shells(Parameters);
  A = Parameters.P + Parameters.N;
  //for(int i = 0; i < Space.num_states; ++i){ Space.qnums[i].energy *= (1.0 - (1.0/A)); } // COM Correction

  // Find range of 2body quantum numbers from mins and maxs
  plus(Space.qmins1, Space.qmins, Space.qmins);
  plus(Space.qmaxs1, Space.qmaxs, Space.qmaxs);
  minus(Space.qmins2, Space.qmins, Space.qmaxs);
  minus(Space.qmaxs2, Space.qmaxs, Space.qmins);
  if(Space.qmins.par == -1 && Space.qmaxs.par == 1){
    Space.qmins1.par = -1;
    Space.qmins2.par = -1;
    Space.qmaxs1.par = 1;
    Space.qmaxs2.par = 1;
  }
  else if(Space.qmins.par == -1 && Space.qmaxs.par == -1){
    Space.qmins1.par = -1;
    Space.qmins2.par = -1;
    Space.qmaxs1.par = -1;
    Space.qmaxs2.par = -1;
  }
  else if(Space.qmins.par == 1 && Space.qmaxs.par == 1){
    Space.qmins1.par = 1;
    Space.qmins2.par = 1;
    Space.qmaxs1.par = 1;
    Space.qmaxs2.par = 1;
  }
  Space.qmins1.j = 0;
  Space.qmins2.j = 0;
  Space.qmaxs2.j = Space.qmaxs1.j;

  minus(Space.qsizes, Space.qmaxs, Space.qmins);
  minus(Space.qsizes1, Space.qmaxs1, Space.qmins1);
  minus(Space.qsizes2, Space.qmaxs2, Space.qmins2);
  if(Space.qmins.par == -1 && Space.qmaxs.par == 1){
    Space.qsizes.par = 2;
    Space.qsizes1.par = 2;
    Space.qsizes2.par = 2;
  }
  else{
    Space.qsizes.par = 1;
    Space.qsizes1.par = 1;
    Space.qsizes2.par = 1;
  }
  Space.qsizes.j = Space.qmaxs.j - Space.qmins.j;
  Space.qsizes1.j = Space.qmaxs1.j - Space.qmins1.j;
  Space.qsizes2.j = Space.qmaxs2.j - Space.qmins2.j;

  Space.qsizes.divide_spins();
  Space.qsizes1.divide_spins();
  Space.qsizes2.divide_spins();
  Space.qsizes.add_one();
  Space.qsizes1.add_one();
  Space.qsizes2.add_one();

  // Find chan3, chan1, and chan2 sizes
  if(Parameters.basis == "infinite"){
    Space.size_1b = Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz * Space.qsizes.m * Space.qsizes.t;
    Space.size_2b_dir = Space.qsizes1.nx * Space.qsizes1.ny * Space.qsizes1.nz * Space.qsizes1.m * Space.qsizes1.t;
    Space.size_2b_cross = Space.qsizes2.nx * Space.qsizes2.ny * Space.qsizes2.nz * Space.qsizes2.m * Space.qsizes2.t;
  }
  else if(Parameters.basis == "finite_M"){
    Space.size_1b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
    Space.size_2b_dir = Space.qsizes1.par * Space.qsizes1.m * Space.qsizes1.t;
    Space.size_2b_cross = Space.qsizes2.par * Space.qsizes2.m * Space.qsizes2.t;
  }
  else if(Parameters.basis == "finite_HO"){
    Space.size_1b = Space.qsizes.ml * Space.qsizes.m;
    Space.size_2b_dir = Space.qsizes1.ml * Space.qsizes1.m;
    Space.size_2b_cross = Space.qsizes2.ml * Space.qsizes2.m;
  }
  else if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
    Space.size_1b = Space.qsizes.par * Space.qsizes.j * Space.qsizes.t;
    Space.size_2b_dir = Space.qsizes1.par * Space.qsizes1.j * Space.qsizes1.t;
    Space.size_2b_cross = Space.qsizes2.par * Space.qsizes2.j * Space.qsizes2.t;
  }
 
  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }

  /*for(int i = 0; i < Space.num_states; ++i){
    std::cout << "###   " << i << ", " << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << " " << Space.qnums[i].t << ", " << Space.qnums[i].energy << std::endl;
    }*/
  for(int i = 0; i < Space.num_states; ++i){
    std::cout << "###   " << i << ", " << Space.qnums[i].n << " " << Space.qnums[i].l << " " << Space.qnums[i].j << " " << Space.qnums[i].t << ", " << Space.qnums[i].energy << std::endl;
  }
}

void Build_Model_Space_J2(Input_Parameters &Parameters, Model_Space &Space)
{
  int ind; // state index
  int m; // angular momenum projection
  int shelllength; // number of single particle states for each shell
  int key;

  State *states = new State[Space.num_states];
  int indtotj = Space.num_states;
  Space.num_jstates = indtotj;
  Space.shellsm = new int*[Space.num_states];
  Space.shellsnum = new int[Space.num_states];
  // reset Space.indtot with degeneracies 2j + 1
  Space.num_states = 0;
  for(int i = 0; i < indtotj; ++i){
    Space.num_states += Space.qnums[i].j + 1;
    Space.shellsnum[i] = Space.qnums[i].j + 1;
    states[i] = Space.qnums[i];
  }
  delete[] Space.qnums;
  Space.qnums = new State[Space.num_states];
  Space.shellsj = new int[Space.num_states];

  // initialize mins and maxs for each quantum number
  Space.qmins.maximize();
  Space.qmins1.maximize();
  Space.qmins2.maximize();
  Space.qmaxs.minimize();
  Space.qmaxs1.minimize();
  Space.qmaxs2.minimize();

  Space.num_p = 0;
  Space.num_n = 0;
  Space.num_hol = 0;
  Space.num_par = 0;

  ind = 0;
  for(int i = 0; i < indtotj; ++i){
    shelllength = states[i].j + 1;
    Space.shellsm[i] = new int[shelllength];
    for(int j = 0; j < shelllength; ++j){
      m = -states[i].j + 2*j;
      Space.qnums[ind].n = states[i].n;
      Space.qnums[ind].par = states[i].par;
      Space.qnums[ind].m = m;
      Space.qnums[ind].t = states[i].t;
      Space.qnums[ind].energy = states[i].energy;
      Space.qnums[ind].type = states[i].type;
      Space.shellsm[i][j] = ind;
      //Space.qnums[ind].j = states[i].j; // TEST!!
      Space.shellsj[ind] = states[i].j;
      if(Space.qnums[ind].par < Space.qmins.par){ Space.qmins.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].par > Space.qmaxs.par){ Space.qmaxs.par = Space.qnums[ind].par; }
      if(Space.qnums[ind].m < Space.qmins.m){ Space.qmins.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].m > Space.qmaxs.m){ Space.qmaxs.m = Space.qnums[ind].m; }
      if(Space.qnums[ind].t < Space.qmins.t){ Space.qmins.t = Space.qnums[ind].t; }
      if(Space.qnums[ind].t > Space.qmaxs.t){ Space.qmaxs.t = Space.qnums[ind].t; }
      if(Space.qnums[ind].t == -1){ ++Space.num_p; }
      else if(Space.qnums[ind].t == 1){ ++Space.num_n; }
      if(Space.qnums[ind].type == 0){ ++Space.num_hol; }
      else if(Space.qnums[ind].type == 1){ ++Space.num_par; }
      ++ind;
    }
  }
  delete[] states;

  // Find range of 2body quantum numbers from mins and maxs
  plus(Space.qmins1, Space.qmins, Space.qmins);
  plus(Space.qmaxs1, Space.qmaxs, Space.qmaxs);
  minus(Space.qmins2, Space.qmins, Space.qmaxs);
  minus(Space.qmaxs2, Space.qmaxs, Space.qmins);
  if(Space.qmins.par == -1 && Space.qmaxs.par == 1){
    Space.qmins1.par = -1;
    Space.qmins2.par = -1;
    Space.qmaxs1.par = 1;
    Space.qmaxs2.par = 1;
  }
  else if(Space.qmins.par == -1 && Space.qmaxs.par == -1){
    Space.qmins1.par = -1;
    Space.qmins2.par = -1;
    Space.qmaxs1.par = -1;
    Space.qmaxs2.par = -1;
  }
  else if(Space.qmins.par == 1 && Space.qmaxs.par == 1){
    Space.qmins1.par = 1;
    Space.qmins2.par = 1;
    Space.qmaxs1.par = 1;
    Space.qmaxs2.par = 1;
  }
  Space.qmins1.j = 0;
  Space.qmins2.j = 0;
  Space.qmaxs2.j = Space.qmaxs1.j;

  minus(Space.qsizes, Space.qmaxs, Space.qmins);
  minus(Space.qsizes1, Space.qmaxs1, Space.qmins1);
  minus(Space.qsizes2, Space.qmaxs2, Space.qmins2);
  if(Space.qmins.par == -1 && Space.qmaxs.par == 1){
    Space.qsizes.par = 2;
    Space.qsizes1.par = 2;
    Space.qsizes2.par = 2;
  }
  else{
    Space.qsizes.par = 1;
    Space.qsizes1.par = 1;
    Space.qsizes2.par = 1;
  }
  Space.qsizes.j = Space.qmaxs.j - Space.qmins.j;
  Space.qsizes1.j = Space.qmaxs1.j - Space.qmins1.j;
  Space.qsizes2.j = Space.qmaxs2.j - Space.qmins2.j;

  Space.qsizes.divide_spins();
  Space.qsizes1.divide_spins();
  Space.qsizes2.divide_spins();
  Space.qsizes.add_one();
  Space.qsizes1.add_one();
  Space.qsizes2.add_one();

  // Find chan3, chan1, and chan2 sizes
  Space.size_1b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
  Space.size_2b_dir = Space.qsizes1.par * Space.qsizes1.m * Space.qsizes1.t;
  Space.size_2b_cross = Space.qsizes2.par * Space.qsizes2.m * Space.qsizes2.t;

  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }

  /*for(int i = 0; i < Space.num_states; ++i){
    std::cout << "###   " << i << ", " << Space.qnums[i].n << " " << Space.qnums[i].j << " " << Space.qnums[i].m << " " << Space.qnums[i].par << " " << Space.qnums[i].t << ", " << Space.qnums[i].energy << " : " << Space.qnums[i].type << std::endl;
    }*/
}

void CART_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  double electron_prefac = hbarc_HartA * hbarc_HartA / (2.0 * m_electronc2_Hart);
  double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
  double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);
  int key;

  // Find appropriate Nmax for the number of shells using shellnums
  if(Parameters.Shells > 40){ std::cerr << "Nmax too big!" << std::endl; exit(1); }
  int shellnums [] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23,
  		      25, 26, 27, 28, 30, 31, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46};
  Space.Nmax = shellnums[Parameters.Shells - 1];
  Space.num_jstates = 0;

  // Set characteristic length (in Bohr radii) for electron gas
  if(Parameters.calc_case == "electronic"){
    Parameters.Nshells = 0;
    double r_b = hbarc_HartA / (m_electronc2_Hart * fine_struct);
    Parameters.density = 3.0/(4.0 * M_PI * std::pow(Parameters.density * r_b, 3.0));
  }

  // Find maximum number of states (indtot) depending on Nmax and whether or not there are protons/neutrons
  if(Parameters.Pshells != 0 && Parameters.Nshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = 1, Space.qsizes.t = 3; }
  else if(Parameters.Pshells != 0 && Parameters.Nshells == 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1; }
  else if(Parameters.Pshells == 0 && Parameters.Nshells != 0){ Space.qmins.t = 1, Space.qmaxs.t = 1, Space.qsizes.t = 1; }
  else{
    if(Parameters.calc_case == "nuclear"){ std::cerr << "No Protons or Neutrons Entered!!!" << std::endl; exit(1); }
    else if(Parameters.calc_case == "electronic"){ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }
  }

  // Find total number of states
  Space.num_states = 0;
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
	      ++Space.num_states;
	    }
	  }
	}
      }
    }
  }

  // allocate memory for quantum numbers for each state
  Space.num_p = 0;
  Space.num_n = 0;
  Space.num_hol = 0;
  Space.num_par = 0;
  Parameters.P = 0;
  Parameters.N = 0;
  Space.qnums = new State[Space.num_states];
  int count = 0;
  for(int shell = 0; shell <= Space.Nmax; ++shell){
    for(int nx = -Space.nmax; nx <= Space.nmax; ++nx){    
      for(int ny = -Space.nmax; ny <= Space.nmax; ++ny){	
	for(int nz = -Space.nmax; nz <= Space.nmax; ++nz){	  
	  if(shell != nx*nx + ny*ny + nz*nz){ continue; }
	  for(int sz = -1; sz <= 1; sz = sz+2){
	    for( int tz = Space.qmins.t; tz <= Space.qmaxs.t; tz = tz+2){
	      Space.qnums[count].energy = 4.0*(nx*nx + ny*ny + nz*nz);	
	      if(tz == -1){
		if(Parameters.Pshells == 0){ continue; }
		++Space.num_p;
		if(shell < Parameters.Pshells){ Space.qnums[count].type = 0; ++Space.num_hol; ++Parameters.P; }
		else{ Space.qnums[count].type = 1; ++Space.num_par; }
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		++Space.num_n;
		if(shell < Parameters.Nshells){ Space.qnums[count].type = 0; ++Space.num_hol; ++Parameters.N; }
		else{ Space.qnums[count].type = 1; ++Space.num_par; }
	      }
	      Space.qnums[count].nx = nx;
	      Space.qnums[count].ny = ny;
	      Space.qnums[count].nz = nz;
	      Space.qnums[count].m = sz;
	      Space.qnums[count].t = tz;
	      count++;
	    }
	  }   
	} 
      }
    }
  }
  Space.qmins.nx = -Space.nmax;
  Space.qmaxs.nx = Space.nmax;
  Space.qmins.ny = -Space.nmax;
  Space.qmaxs.ny = Space.nmax;
  Space.qmins.nz = -Space.nmax;
  Space.qmaxs.nz = Space.nmax;
  Space.qmins.m = -1;
  Space.qmaxs.m = 1;
  
  // Find range of 2body quantum numbers from mins and maxs
  plus(Space.qmins1, Space.qmins, Space.qmins);
  plus(Space.qmaxs1, Space.qmaxs, Space.qmaxs);
  minus(Space.qmins2, Space.qmins, Space.qmaxs);
  minus(Space.qmaxs2, Space.qmaxs, Space.qmins);

  minus(Space.qsizes, Space.qmaxs, Space.qmins);
  Space.qsizes.divide_spins();
  Space.qsizes.add_one();

  minus(Space.qsizes1, Space.qmaxs1, Space.qmins1);
  Space.qsizes1.divide_spins();
  Space.qsizes1.add_one();

  minus(Space.qsizes2, Space.qmaxs2, Space.qmins2);
  Space.qsizes2.divide_spins();
  Space.qsizes2.add_one();

  Space.size_1b = Space.qsizes.nx * Space.qsizes.ny * Space.qsizes.nz * Space.qsizes.m * Space.qsizes.t;
  Space.size_2b_dir = Space.qsizes1.nx * Space.qsizes1.ny * Space.qsizes1.nz * Space.qsizes1.m * Space.qsizes1.t;
  Space.size_2b_cross = Space.qsizes2.nx * Space.qsizes2.ny * Space.qsizes2.nz * Space.qsizes2.m * Space.qsizes2.t;

  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }

  // Calculate energies
  double L = pow(Space.num_hol/Parameters.density, 1.0/3.0);
  if(Parameters.calc_case == "nuclear"){
    for(int i = 0; i < Space.num_states; ++i){
      if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac * M_PI*M_PI / (L*L); }
      else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac * M_PI*M_PI / (L*L); }
    }
    // Change energies to Hartree-Fock energies, E_p = E_p + V_pipi
    for(int p = 0; p < Space.num_states; ++p){
      for(int i = 0; i < Space.num_states; ++i){
	if(Space.qnums[i].type != 0){ continue; }
	if(p == i){ continue; }
	Space.qnums[p].energy += vint_Minnesota_Momentum(Space, p, i, p, i, L);
      }
    }
  }
  else if(Parameters.calc_case == "electronic"){
    for(int i = 0; i < Space.num_states; ++i){
      Space.qnums[i].energy *= electron_prefac * M_PI*M_PI / (L*L);
    }
    // Change energies to Hartree-Fock energies, E_p = E_p + V_pipi
    for(int p = 0; p < Space.num_states; ++p){
      for(int i = 0; i < Space.num_states; ++i){
	if(Space.qnums[i].type != 0){ continue; }
	if(p == i){ continue; }
	Space.qnums[p].energy += Coulomb_Inf(Space, p, i, p, i, L);;
      }
    }
  }
}

void QD_Build_Model_Space(Input_Parameters &Parameters, Model_Space &Space)
{
  int key;
  // initialize mins and maxs for each quantum number
  Space.qmins.maximize();
  Space.qmins1.maximize();
  Space.qmins2.maximize();
  Space.qmaxs.minimize();
  Space.qmaxs1.minimize();
  Space.qmaxs2.minimize();
  Space.num_jstates = 0;

  Parameters.Nshells = 0;
  if(Parameters.Pshells != 0){ Space.qmins.t = -1, Space.qmaxs.t = -1, Space.qsizes.t = 1; }
  else{ std::cerr << "No Electrons Entered!!!" << std::endl; exit(1); }

  // Find total number of states
  Space.num_states = 0;
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int ml = -1 * shell; ml <= shell; ++ml){
      for(int n = 0; n < Parameters.Shells; ++n){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  if(n < Space.qmins.n){ Space.qmins.n = n; }
	  if(n > Space.qmaxs.n){ Space.qmaxs.n = n; }
	  if(ml < Space.qmins.ml){ Space.qmins.ml = ml; }
	  if(ml > Space.qmaxs.ml){ Space.qmaxs.ml = ml; }
	  ++Space.num_states;
	}
      }
    }
  }

  // Allocate memory for quantum numbers for each state
  Space.num_p = 0;
  Space.num_n = 0;
  Space.num_hol = 0;
  Space.num_par = 0;
  Parameters.P = 0;
  Parameters.N = 0;
  Space.nmax = 0;
  Space.qnums = new State[Space.num_states];
  int count = 0;
  for(int shell = 0; shell < Parameters.Shells; ++shell){
    for(int ml = -1 * shell; ml <= shell; ++ml){
      for(int n = 0; n < Parameters.Shells; ++n){
	if(2*n + abs(ml) != shell){ continue; }
	for(int sz = -1; sz <= 1; sz = sz+2){
	  Space.qnums[count].energy = Parameters.density*(shell + 1.0);
	  Space.qnums[count].n = n;
	  Space.qnums[count].par = -2*(abs(ml)%2) + 1;
	  Space.qnums[count].ml = ml;
	  Space.qnums[count].m = sz;
	  if(shell < Parameters.Pshells){ Space.qnums[count].type = 0; ++Space.num_hol; ++Parameters.P; }
	  else{ Space.qnums[count].type = 1; ++Space.num_par; }
	  ++Space.num_p;
	  count++;
	}
      }
    }
  }
  Space.qmins.m = -1;
  Space.qmaxs.m = 1;

  // Find range of 2body quantum numbers from mins and maxs
  plus(Space.qmins1, Space.qmins, Space.qmins);
  plus(Space.qmaxs1, Space.qmaxs, Space.qmaxs);
  minus(Space.qmins2, Space.qmins, Space.qmaxs);
  minus(Space.qmaxs2, Space.qmaxs, Space.qmins);

  minus(Space.qsizes, Space.qmaxs, Space.qmins);
  Space.qsizes.divide_spins();
  Space.qsizes.add_one();

  minus(Space.qsizes1, Space.qmaxs1, Space.qmins1);
  Space.qsizes1.divide_spins();
  Space.qsizes1.add_one();

  minus(Space.qsizes2, Space.qmaxs2, Space.qmins2);
  Space.qsizes2.divide_spins();
  Space.qsizes2.add_one();

  Space.size_1b = Space.qsizes.ml * Space.qsizes.m;
  Space.size_2b_dir = Space.qsizes1.ml * Space.qsizes1.m;
  Space.size_2b_cross = Space.qsizes2.ml * Space.qsizes2.m;

  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);    
    Space.map_state[key] = i;
    //std::cout << "###   " << i << ", " << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << ", " << Space.qnums[i].energy << ", " << key << std::endl;
  }
}

//Function to setup Channels
Channels::Channels(Input_Parameters &Parameters, Model_Space &Space)
{
  State state1, state2, state3;
  int chan1, chan2, chan3;
  int h_total, p_total;
  int hh_total, pp_total, hp1_total, ph1_total, pph_total, hhp_total;
  int hp_total, hh1_total, pp1_total, hph_total, hpp_total, hhh_total, ppp_total;

  size1 = Space.size_2b_dir;
  size2 = Space.size_2b_cross;
  size3 = Space.size_1b;

  qnums1 = new State[size1];
  qnums2 = new State[size2];
  qnums3 = new State[size3];

  p_chan_setup(Parameters,Space, h_total,nh,h_index, h_map,h_vec, size3,qnums3, 0); // h
  p_chan_setup(Parameters,Space, p_total,np,p_index, p_map,p_vec, size3,qnums3, 1); // p

  pq_chan_setup(Parameters,Space, hh_total,nhh,hh_index,hh_map,hh_vec, size1,qnums1, size3,qnums3, nh,h_vec,h_index,nh,h_index,h_vec); // hh
  pq_chan_setup(Parameters,Space, pp_total,npp,pp_index,pp_map,pp_vec, size1,qnums1, size3,qnums3, np,p_vec,p_index,np,p_index,p_vec); // pp
  pq1_chan_setup(Parameters,Space, hp1_total,nhp1,hp1_index,hp1_map,hp1_vec, size2,qnums2, size3,qnums3, nh,h_vec,h_index,np,p_index,p_vec); // hp1
  pq1_chan_setup(Parameters,Space, ph1_total,nph1,ph1_index,ph1_map,ph1_vec, size2,qnums2, size3,qnums3, np,p_vec,p_index,nh,h_index,h_vec); // ph1

  /*for(int chan = 0; chan < size1; ++chan){
    std::cout << "chan = " << chan << ", j = " << qnums1[chan].j << ", npp = " << npp[chan] << ", nhh = " << nhh[chan] << std::endl;
    for(int pp = 0; pp < npp[chan]; ++pp){
      std::cout << pp_vec[pp_index[chan] + pp].v1 << "," << pp_vec[pp_index[chan] + pp].v2 << "  ";
    }
    std::cout << std::endl;
    for(int hh = 0; hh < nhh[chan]; ++hh){
      std::cout << hh_vec[hh_index[chan] + hh].v1 << "," << hh_vec[hh_index[chan] + hh].v2 << "  ";
    }
    std::cout << std::endl;
    }*/

  pqr_chan_setup(Parameters,Space, hhp_total,nhhp,hhp_index,hhp_map,hhp_vec,hhp_j, size3,qnums3, nhh,hh_vec,hh_index, np,p_index,p_vec); // hhp
  pqr_chan_setup(Parameters,Space, pph_total,npph,pph_index,pph_map,pph_vec,pph_j, size3,qnums3, npp,pp_vec,pp_index, nh,h_index,h_vec); // pph

  // For CCSD
  if(Parameters.approx == "singles"){
    pq_chan_setup(Parameters,Space, hp_total,nhp,hp_index,hp_map,hp_vec, size1,qnums1, size3,qnums3, nh,h_vec,h_index, np,p_index,p_vec); // hp
    pq1_chan_setup(Parameters,Space, hh1_total,nhh1,hh1_index,hh1_map,hh1_vec, size2,qnums2, size3,qnums3, nh,h_vec,h_index, nh,h_index,h_vec); // hh1
    pq1_chan_setup(Parameters,Space, pp1_total,npp1,pp1_index,pp1_map,pp1_vec, size2,qnums2, size3,qnums3, np,p_vec,p_index, np,p_index,p_vec); // pp1

    pqr_chan_setup(Parameters,Space, hph_total,nhph,hph_index,hph_map,hph_vec,hph_j, size3,qnums3, nhp,hp_vec,hp_index, nh,h_index,h_vec); // hph
    pqr_chan_setup(Parameters,Space, hpp_total,nhpp,hpp_index,hpp_map,hpp_vec,hpp_j, size3,qnums3, nhp,hp_vec,hp_index, np,p_index,p_vec); // hpp
    pqr_chan_setup(Parameters,Space, hhh_total,nhhh,hhh_index,hhh_map,hhh_vec,hhh_j, size3,qnums3, nhh,hh_vec,hh_index, nh,h_index,h_vec); // hhh
    pqr_chan_setup(Parameters,Space, ppp_total,nppp,ppp_index,ppp_map,ppp_vec,ppp_j, size3,qnums3, npp,pp_vec,pp_index, np,p_index,p_vec); // ppp

    // find Chan.ind0
    state2.nx = 0;
    state2.ny = 0;
    state2.nz = 0;
    state2.t = 0;
    state2.m = 0;
    state2.par = 1;
    state2.ml = 0;
    state2.j = 0;
    ind0 = Space.ind_2b_cross(Parameters.basis, state2);
  }
  
  double memory = 0.0;
  int doubsize = sizeof(double);
  int intsize = sizeof(int);
  int mapsize = sizeof(size_t) + sizeof(void*);
  for(chan3 = 0; chan3 < size3; ++chan3){
    memory += doubsize * nh[chan3] * nh[chan3]; // F_matrix hh
    memory += doubsize * nh[chan3] * np[chan3]; // F_matrix hp
    memory += doubsize * np[chan3] * nh[chan3]; // F_matrix ph
    memory += doubsize * np[chan3] * np[chan3]; // F_matrix pp
    memory += doubsize * nh[chan3] * nhhh[chan3]; // Vhhhh3_2
    memory += doubsize * nppp[chan3] * np[chan3]; // Vpppp3_3
    memory += doubsize * nh[chan3] * npph[chan3]; // Vhhpp3_1
    memory += doubsize * nh[chan3] * npph[chan3]; // Vhhpp3_2
    memory += doubsize * nhhp[chan3] * np[chan3]; // Vhhpp3_3
    memory += doubsize * nh[chan3] * nhph[chan3]; // Vhhhp3_2
    memory += doubsize * nhhp[chan3] * nh[chan3]; // Vhhhp3_3
    memory += doubsize * np[chan3] * npph[chan3]; // Vhppp3_2
    memory += doubsize * np[chan3] * nhhh[chan3]; // Vhphh3_2
    memory += doubsize * nppp[chan3] * nh[chan3]; // Vpphp3_3
    memory += 4 * doubsize * nh[chan3] * nh[chan3]; // Xhh_3, Xhh_3od, X1hh_3, X1hh_3od
    memory += 2 * doubsize * np[chan3] * np[chan3]; // Xpp_3, Xpp_3od
    memory += doubsize * nh[chan3] * np[chan3]; // Xhp_3
    memory += doubsize * nhhh[chan3] * nh[chan3]; // Xhhhh3_3
    memory += doubsize * nhhh[chan3] * nh[chan3]; // Xhhhh3_4
    memory += doubsize * np[chan3] * nppp[chan3]; // X1pppp3_1
    memory += doubsize * np[chan3] * nppp[chan3]; // X1pppp3_2
    memory += 2 * doubsize * nh[chan3] * nhpp[chan3]; // X1hphp3_1, X2hphp3_1
    memory += 4 * doubsize * np[chan3] * nhph[chan3]; // Xhphp3_2, X1hphp3_2, X2hphp3_2, X3hphp3_2
    memory += 4 * doubsize * nhpp[chan3] * nh[chan3]; // Xhphp3_3, X1hphp3_3, X2hphp3_3, X3hphp3_3
    memory += doubsize * nhph[chan3] * np[chan3]; // X3hphp3_4
    memory += 2 * doubsize * nhhp[chan3] * nh[chan3]; // Xhhhp3_3, X1hhhp3_3
    memory += 2 * doubsize * nhhh[chan3] * np[chan3]; // X1hhhp3_4
    memory += doubsize * nh[chan3] * nppp[chan3]; // X1hppp3_1
    memory += 2 * doubsize * np[chan3] * npph[chan3]; // Xhppp3_2, X1hppp3_2
    memory += 2 * doubsize * nhpp[chan3] * np[chan3]; // Xhppp3_3, X1hppp3_3
    memory += 2 * doubsize * nh[chan3] * nhhp[chan3]; // Xhphh3_1, X1hphh3_1
    memory += 2 * doubsize * np[chan3] * nhhh[chan3]; // Xhphh3_2, X1hphh3_2
    memory += doubsize * nhph[chan3] * nh[chan3]; // Xhphh3_3
    memory += doubsize * nhph[chan3] * nh[chan3]; // Xhphh3_4
    memory += 2 * doubsize * np[chan3] * nhpp[chan3]; // Xpphp3_1, X1pphp3_1
    memory += 2 * doubsize * np[chan3] * nhpp[chan3]; // Xpphp3_2, X1pphp3_2
    memory += 2 * doubsize * nppp[chan3] * nh[chan3]; // Xpphp3_3, X1pphp3_3
    memory += 2 * doubsize * npph[chan3] * np[chan3]; // Xpphp3_4, X1pphp3_4
    memory += 6 * doubsize * np[chan3] * nhhp[chan3]; // 3 copies of T3_1, T3_2
    memory += 6 * doubsize * npph[chan3] * nh[chan3]; // 3 copies of T3_3, T3_4
    memory += (intsize * h_map[chan3].size()) + (h_map[chan3].bucket_count() * mapsize); // h_vec
    memory += (intsize * p_map[chan3].size()) + (p_map[chan3].bucket_count() * mapsize); // p_vec
    memory += (3 * intsize * hhp_map[chan3].size()) + (hhp_map[chan3].bucket_count() * mapsize); // hhp_vec
    memory += (3 * intsize * pph_map[chan3].size()) + (pph_map[chan3].bucket_count() * mapsize); // pph_vec
    memory += (3 * intsize * hpp_map[chan3].size()) + (hpp_map[chan3].bucket_count() * mapsize); // hpp_vec
    memory += (3 * intsize * hph_map[chan3].size()) + (hph_map[chan3].bucket_count() * mapsize); // hph_vec
    memory += (3 * intsize * hhh_map[chan3].size()) + (hhh_map[chan3].bucket_count() * mapsize); // hhh_vec
    memory += (3 * intsize * ppp_map[chan3].size()) + (ppp_map[chan3].bucket_count() * mapsize); // ppp_vec
  }
  for(chan1 = 0; chan1 < size1; ++chan1){
    memory += doubsize * nhh[chan1] * nhh[chan1]; // Vhhhh1
    memory += doubsize * npp[chan1] * npp[chan1]; // Vpppp1
    memory += doubsize * nhh[chan1] * npp[chan1]; // Vhhpp1
    memory += doubsize * npp[chan1] * nhh[chan1]; // Vpphh1
    memory += 3 * doubsize * nhh[chan1] * nhh[chan1]; // Xhhhh1, hhhh_fac(2)
    memory += 4 * intsize * nhh[chan1] * nhh[chan1]; // hhhh_chan(2), hhhh_ind(2)
    memory += 4 * doubsize * npp[chan1] * npp[chan1]; // Xpppp1, X1pppp1, pppp_fac(2)
    memory += 4 * intsize * npp[chan1] * npp[chan1]; // pppp_chan(2), pppp_ind(2)
    memory += 5 * doubsize * nhp[chan1] * nhp[chan1]; // hphp_fac(5)
    memory += 10 * intsize * nhp[chan1] * nhp[chan1]; // hphp_chan(5), hphp_ind(5)
    memory += 4 * doubsize * nhh[chan1] * nhp[chan1]; // Xhhhp1, hhhp_fac(3)
    memory += 6 * intsize * nhh[chan1] * nhp[chan1]; // hhhp_chan(3), hhhp_ind(3)
    memory += 5 * doubsize * nhp[chan1] * npp[chan1]; // Xhppp1, hppp_fac(4)
    memory += 8 * intsize * nhp[chan1] * npp[chan1]; // hppp_chan(4), hppp_ind(4)
    memory += 7 * doubsize * nhp[chan1] * nhh[chan1]; // Xhphh1, hphh_fac(6)
    memory += 12 * intsize * nhp[chan1] * nhh[chan1]; // hphh_chan(6), hphh_ind(6)
    memory += 7 * doubsize * npp[chan1] * nhp[chan1]; // Xpphp1, pphp_fac(6)
    memory += 12 * intsize * npp[chan1] * nhp[chan1]; // pphp_chan(6), pphp_ind(6)
    memory += 27 * doubsize * npp[chan1] * nhh[chan1]; // 3 copies of T1, T_fac(8)
    memory += 48 * intsize * npp[chan1] * nhh[chan1]; // 3 copies of T_chan(8), T_ind(8)
    memory += (2 * intsize * hh_map[chan1].size()) + (hh_map[chan1].bucket_count() * mapsize); // hh_vec
    memory += (2 * intsize * pp_map[chan1].size()) + (pp_map[chan1].bucket_count() * mapsize); // pp_vec
    memory += (2 * intsize * hp_map[chan1].size()) + (hp_map[chan1].bucket_count() * mapsize); // hp_vec
  }
  for(chan2 = 0; chan2 < size2; ++chan2){
    memory += doubsize * nhp1[chan2] * nph1[chan2]; // Vhhpp2_1
    memory += doubsize * nhp1[chan2] * nph1[chan2]; // Vhhpp2_3
    memory += doubsize * nhp1[chan2] * nhp1[chan2]; // Vhphp2_1
    memory += doubsize * nph1[chan2] * nph1[chan2]; // Vhphp2_2
    memory += doubsize * nhh1[chan2] * nph1[chan2]; // Vhhhp2_3
    memory += doubsize * npp1[chan2] * nph1[chan2]; // Vhppp2_4
    memory += 4 * doubsize * nhp1[chan2] * nph1[chan2]; // Xhphp2_1, X1hphp2_1, X2hphp2_1, X3hphp2_1
    memory += doubsize * nhh1[chan2] * nph1[chan2]; // Xhhhp2_3
    memory += doubsize * nhp1[chan2] * npp1[chan2]; // Xhppp2_3
    memory += doubsize * nhh1[chan2] * nhp1[chan2]; // Xhphh2_1
    memory += doubsize * nhh1[chan2] * nhp1[chan2]; // Xhphh2_3
    memory += doubsize * nph1[chan2] * npp1[chan2]; // Xpphp2_2
    memory += doubsize * nph1[chan2] * npp1[chan2]; // Xpphp2_3
    memory += 12 * doubsize * nph1[chan2] * nhp1[chan2]; // 3 copies of T2_1, T2_2, T2_3, T2_4
    memory += (2 * intsize * hp1_map[chan2].size()) + (hp1_map[chan2].bucket_count() * mapsize); // hp1_vec
    memory += (2 * intsize * ph1_map[chan2].size()) + (ph1_map[chan2].bucket_count() * mapsize); // ph1_vec
    memory += (2 * intsize * hh1_map[chan2].size()) + (hh1_map[chan2].bucket_count() * mapsize); // hh1_vec
    memory += (2 * intsize * pp1_map[chan2].size()) + (pp1_map[chan2].bucket_count() * mapsize); // pp1_vec
  }
  std::cout << std::endl << "Estimated Memory = " << memory/1000000.0 << " MB" << std::endl << std::endl;
}

void Channels::delete_struct(Input_Parameters &Parameters)
{
  delete[] qnums1;
  delete[] qnums2;
  delete[] qnums3;

  delete[] h_vec;
  delete[] p_vec;
  delete[] h_index;
  delete[] p_index;
  delete[] nh;
  delete[] np;
  delete[] h_map;
  delete[] p_map;

  delete[] hh_vec;
  delete[] pp_vec;
  delete[] hh_index;
  delete[] pp_index;
  delete[] nhh;
  delete[] npp;
  delete[] hh_map;
  delete[] pp_map;
  delete[] hp1_vec;
  delete[] ph1_vec;
  delete[] hp1_index;
  delete[] ph1_index;
  delete[] nhp1;
  delete[] nph1;
  delete[] hp1_map;
  delete[] ph1_map;
  delete[] pph_vec;
  delete[] hhp_vec;
  delete[] pph_j;
  delete[] hhp_j;
  delete[] pph_index;
  delete[] hhp_index;
  delete[] npph;
  delete[] nhhp;
  delete[] pph_map;
  delete[] hhp_map;
  if(Parameters.approx == "singles"){
    delete[] hp_vec;
    delete[] hp_index;
    delete[] nhp;
    delete[] hp_map;
    delete[] hh1_vec;
    delete[] pp1_vec;
    delete[] hh1_index;
    delete[] pp1_index;
    delete[] nhh1;
    delete[] npp1;
    delete[] hh1_map;
    delete[] pp1_map;
    delete[] hph_vec;
    delete[] hpp_vec;
    delete[] hhh_vec;
    delete[] ppp_vec;
    delete[] hph_j;
    delete[] hpp_j;
    delete[] hhh_j;
    delete[] ppp_j;
    delete[] hph_index;
    delete[] hpp_index;
    delete[] hhh_index;
    delete[] ppp_index;
    delete[] nhph;
    delete[] nhpp;
    delete[] nhhh;
    delete[] nppp;
    delete[] hph_map;
    delete[] hpp_map;
    delete[] hhh_map;
    delete[] ppp_map;
  }
}

void Model_Space::Determine_Shells(Input_Parameters &Parameters)
{
  int *pshell;
  int *nshell;
  int pshell_num = 0;
  int nshell_num = 0;
  double p_en = -1000.0;
  double n_en = -1000.0;
  for(int i = 0; i < num_states; ++i){
    if(qnums[i].energy != p_en && qnums[i].t == -1){
      p_en = qnums[i].energy;
      ++pshell_num;
    }
    if(qnums[i].energy != n_en && qnums[i].t == 1){
      n_en = qnums[i].energy;
      ++nshell_num;
    }
  }
  pshell = new int[pshell_num];
  nshell = new int[nshell_num];
  pshell_num = 0;
  nshell_num = 0;
  p_en = -1000.0;
  n_en = -1000.0;
  for(int i = 0; i < num_states; ++i){
    if(qnums[i].energy != p_en && qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = qnums[i].energy;
      ++pshell_num;
    }
    if(qnums[i].energy != n_en && qnums[i].t == 1){
      nshell[nshell_num] = i;
      n_en = qnums[i].energy;
      ++nshell_num;
    }
  }
  Parameters.Shells = std::max(pshell_num, nshell_num);

  num_p = 0;
  num_n = 0;
  num_hol = 0;
  num_par = 0;
  Parameters.P = 0;
  Parameters.N = 0;
  if(Parameters.Pshells > pshell_num){ std::cerr << "Pshells too big for space!" << std::endl; exit(1); }
  if(Parameters.Nshells > nshell_num){ std::cerr << "Nshells too big for space!" << std::endl; exit(1); }
  for(int i = 0; i < num_states; ++i){
    if(qnums[i].t == -1){
      ++num_p;
      if(i < pshell[Parameters.Pshells]){ qnums[i].type = 0; ++num_hol; ++Parameters.P; }
      else{ qnums[i].type = 1; ++num_par; }
    }
    if(qnums[i].t == 1){
      ++num_n;
      if(i < nshell[Parameters.Nshells]){ qnums[i].type = 0; ++num_hol; ++Parameters.N; }
      else{ qnums[i].type = 1; ++num_par; }
    }
  }
  delete[] pshell;
  delete[] nshell;
  if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
    for(int i = 0; i < num_hol; ++i){
      if(qnums[i].t == -1){ Parameters.P += qnums[i].j; }
      else if(qnums[i].t == 1){ Parameters.N += qnums[i].j; }
    }
  }
}

void p_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &p_total, int *&np, int *&p_index, std::unordered_map<int,int> *&p_map, one_body *&p_vec, int &chan3size, State *qnums3, int type){
  int np0;
  int chan3;
  one_body p_state;
  // p_total, np, p_index, p_map, p_vec
  p_total = 0;
  np = new int[chan3size];
  p_index = new int[chan3size];
  p_map = new std::unordered_map<int,int>[chan3size];
  for(chan3 = 0; chan3 < chan3size; ++chan3){ np[chan3] = 0; }
  for(int p = 0; p < Space.num_states; ++p){
    if(Space.qnums[p].type == type || type == -1){
      chan3 = Space.ind_1b(Parameters.basis, Space.qnums[p]);
      ++np[chan3];
      ++p_total;
    }
  }
  p_vec = new one_body[p_total];
  p_total = 0;
  for(chan3 = 0; chan3 < chan3size; ++chan3){
    p_index[chan3] = p_total;
    p_total += np[chan3];
    np[chan3] = 0;
  }
  for(int p = 0; p < Space.num_states; ++p){
    if(Space.qnums[p].type == type || type == -1){
      chan3 = Space.ind_1b(Parameters.basis, Space.qnums[p]);
      qnums3[chan3] = Space.qnums[p];
      p_state.v1 = p;
      np0 = np[chan3];
      p_vec[p_index[chan3] + np0] = p_state;
      p_map[chan3][p] = np0;
      ++np[chan3];
    }
  }
}

void pq_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pq_total, int *&npq, int *&pq_index, std::unordered_map<int,int> *&pq_map, two_body *&pq_vec, int &chan1size, State *qnums1, int &chan3size, State *qnums3, int *np, one_body *p_vec, int *p_index, int *nq, int *q_index, one_body *q_vec){
  int np0, nq0, npq0;
  int p_ind, q_ind;
  int p, q;
  State tb;
  int chan1;
  int jmin, key;
  two_body pq_state;
  // pq_total, npq, pq_index, pq_map, pq_vec
  pq_total = 0;
  npq = new int[chan1size];
  pq_index = new int[chan1size];
  pq_map = new std::unordered_map<int,int>[chan1size];
  for(chan1 = 0; chan1 < chan1size; ++chan1){ npq[chan1] = 0; }
  for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
    np0 = np[chan3_1];
    p_ind = p_index[chan3_1];
    for(int chan3_2 = 0; chan3_2 < chan3size; ++chan3_2){
      nq0 = nq[chan3_2];
      q_ind = q_index[chan3_2];
      plus(tb, qnums3[chan3_1], qnums3[chan3_2]);
      jmin = abs(qnums3[chan3_1].j - qnums3[chan3_2].j);
      while(tb.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	for(int p0 = 0; p0 < np0; ++p0){
	  p = p_vec[p_ind + p0].v1;
	  for(int q0 = 0; q0 < nq0; ++q0){
	    q = q_vec[q_ind + q0].v1;
	    if(p == q && tb.j%4 != 0){ continue; }
	    ++npq[chan1];
	    ++pq_total;
	  }
	}
	tb.j -= 2;
      }
    }
  }
  pq_vec = new two_body[pq_total];
  pq_total = 0;
  for(chan1 = 0; chan1 < chan1size; ++chan1){
    pq_index[chan1] = pq_total;
    pq_total += npq[chan1];
    npq[chan1] = 0;
  }
  for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
    np0 = np[chan3_1];
    p_ind = p_index[chan3_1];
    for(int chan3_2 = 0; chan3_2 < chan3size; ++chan3_2){
      nq0 = nq[chan3_2];
      q_ind = q_index[chan3_2];
      plus(tb, qnums3[chan3_1], qnums3[chan3_2]);
      jmin = abs(qnums3[chan3_1].j - qnums3[chan3_2].j);
      while(tb.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	for(int p0 = 0; p0 < np0; ++p0){
	  p = p_vec[p_ind + p0].v1;
	  for(int q0 = 0; q0 < nq0; ++q0){
	    q = q_vec[q_ind + q0].v1;
	    pq_state.v1 = p;
	    pq_state.v2 = q;
	    key = Space.hash2(p, q, tb.j);
	    chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	    qnums1[chan1] = tb;
	    if(p == q && tb.j%4 != 0){ continue; }
	    npq0 = npq[chan1];
	    pq_vec[pq_index[chan1] + npq0] = pq_state;
	    pq_map[chan1][key] = npq0;
	    ++npq[chan1];
	  }
	}
	tb.j -= 2;
      }
    }
  }
}

void pq1_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pq_total, int *&npq, int *&pq_index, std::unordered_map<int,int> *&pq_map, two_body *&pq_vec, int &chan2size, State *qnums2, int &chan3size, State *qnums3, int *np, one_body *p_vec, int *p_index, int *nq, int *q_index, one_body *q_vec){
  int np0, nq0, npq0;
  int p_ind, q_ind;
  int p, q;
  State tb;
  int chan2;
  int jmin, key;
  two_body pq_state;
  // pq_total, npq, pq_index, pq_map, pq_vec
  pq_total = 0;
  npq = new int[chan2size];
  pq_index = new int[chan2size];
  pq_map = new std::unordered_map<int,int>[chan2size];
  for(chan2 = 0; chan2 < chan2size; ++chan2){ npq[chan2] = 0; }
  for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
    np0 = np[chan3_1];
    p_ind = p_index[chan3_1];
    for(int chan3_2 = 0; chan3_2 < chan3size; ++chan3_2){
      nq0 = nq[chan3_2];
      q_ind = q_index[chan3_2];
      minus(tb, qnums3[chan3_1], qnums3[chan3_2]);
      jmin = abs(qnums3[chan3_1].j - qnums3[chan3_2].j);
      while(tb.j >= jmin){
	chan2 = Space.ind_2b_cross(Parameters.basis, tb);
	for(int p0 = 0; p0 < np0; ++p0){
	  p = p_vec[p_ind + p0].v1;
	  for(int q0 = 0; q0 < nq0; ++q0){
	    q = q_vec[q_ind + q0].v1;
	    ++npq[chan2];
	    ++pq_total;
	  }
	}
	tb.j -= 2;
      }
    }
  }
  pq_vec = new two_body[pq_total];
  pq_total = 0;
  for(chan2 = 0; chan2 < chan2size; ++chan2){
    pq_index[chan2] = pq_total;
    pq_total += npq[chan2];
    npq[chan2] = 0;
  }
  for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
    np0 = np[chan3_1];
    p_ind = p_index[chan3_1];
    for(int chan3_2 = 0; chan3_2 < chan3size; ++chan3_2){
      nq0 = nq[chan3_2];
      q_ind = q_index[chan3_2];
      minus(tb, qnums3[chan3_1], qnums3[chan3_2]);
      jmin = abs(qnums3[chan3_1].j - qnums3[chan3_2].j);
      while(tb.j >= jmin){
	chan2 = Space.ind_2b_cross(Parameters.basis, tb);
	for(int p0 = 0; p0 < np0; ++p0){
	  p = p_vec[p_ind + p0].v1;
	  for(int q0 = 0; q0 < nq0; ++q0){
	    q = q_vec[q_ind + q0].v1;
	    pq_state.v1 = p;
	    pq_state.v2 = q;
	    key = Space.hash2(p, q, tb.j);
	    chan2 = Space.ind_2b_cross(Parameters.basis, tb);
	    qnums2[chan2] = tb;
	    npq0 = npq[chan2];
	    pq_vec[pq_index[chan2] + npq0] = pq_state;
	    pq_map[chan2][key] = npq0;
	    ++npq[chan2];
	  }
	}
	tb.j -= 2;
      }
    }
  }
}

void pqr_chan_setup(Input_Parameters &Parameters, Model_Space &Space, int &pqr_total, int *&npqr, int *&pqr_index, std::unordered_map<int,int> *&pqr_map, three_body *&pqr_vec, State *&pqr_j, int &chan3size, State *qnums3, int *npq, two_body *pq_vec, int *pq_index, int *nr, int *r_index, one_body *r_vec){
  int npq0, nr0, npqr0;
  int pq_ind, r_ind;
  int p, q, r;
  State tb;
  int chan1, chan3;
  int jmin, key;
  three_body pqr_state;
  // pqr_total, npqr, pqr_index, pqr_map, pqr_vec, pqr_j
  pqr_total = 0;
  npqr = new int[chan3size];
  pqr_index = new int[chan3size];
  pqr_map = new std::unordered_map<int,int>[chan3size];
  for(chan3 = 0; chan3 < chan3size; ++chan3){ npqr[chan3] = 0; }
  for(chan3 = 0; chan3 < chan3size; ++chan3){
    for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
      nr0 = nr[chan3_1];
      plus(tb, qnums3[chan3], qnums3[chan3_1]);
      jmin = abs(qnums3[chan3].j - qnums3[chan3_1].j);
      while(tb.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	for(int r0 = 0; r0 < nr0; ++r0){
	  npqr[chan3] += npq[chan1]; // <pq|rs> -> <pqr'|s>
	  pqr_total += npq[chan1];
	}
	tb.j -= 2;
      }
    }
  }
  pqr_vec = new three_body[pqr_total];
  pqr_j = new State[pqr_total];
  pqr_total = 0;
  for(chan3 = 0; chan3 < chan3size; ++chan3){
    pqr_index[chan3] = pqr_total;
    pqr_total += npqr[chan3];
    npqr[chan3] = 0;
  }
  for(chan3 = 0; chan3 < chan3size; ++chan3){
    for(int chan3_1 = 0; chan3_1 < chan3size; ++chan3_1){
      nr0 = nr[chan3_1];
      r_ind = r_index[chan3_1];
      plus(tb, qnums3[chan3], qnums3[chan3_1]);
      jmin = abs(qnums3[chan3].j - qnums3[chan3_1].j);
      while(tb.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, tb);
	npq0 = npq[chan1];
	pq_ind = pq_index[chan1];
	for(int r0 = 0; r0 < nr0; ++r0){
	  r = r_vec[r_ind + r0].v1;
	  for(int pq0 = 0; pq0 < npq0; ++pq0){ // <pq|rs> -> <pqr'|s>
	    p = pq_vec[pq_ind + pq0].v1;
	    q = pq_vec[pq_ind + pq0].v2;
	    npqr0 = npqr[chan3];
	    key = Space.hash3(p, q, r, tb.j);
	    pqr_state.v1 = p;
	    pqr_state.v2 = q;
	    pqr_state.v3 = r;
	    pqr_vec[pqr_index[chan3] + npqr0] = pqr_state;
	    pqr_j[pqr_index[chan3] + npqr0] = tb;
	    pqr_map[chan3][key] = npqr0;
	    ++npqr[chan3];
	  }
	}
	tb.j -= 2;
      }
    }
  }
}

void Map_4_count_1(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb)
{
  State tb;
  int jmin, num;
  cross_state(Parameters, Space, fb.v1, fb.v4, fb.v3, fb.v2, jmin, tb); // 2_1
  num = (tb.j - jmin)/2 + 1;
  map_num[size * count + offset] = num;
  map_index[size * count + offset] = length;
  length += num;
}

void Map_4_count_2(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb)
{
  State tb;
  int jmin, num;
  cross_state(Parameters, Space, fb.v2, fb.v3, fb.v4, fb.v1, jmin, tb); // 2_2
  num = (tb.j - jmin)/2 + 1;
  map_num[size * count + offset] = num;
  map_index[size * count + offset] = length;
  length += num;
}

void Map_4_count_3(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb)
{
  State tb;
  int jmin, num;
  cross_state(Parameters, Space, fb.v1, fb.v3, fb.v4, fb.v2, jmin, tb); // 2_3
  num = (tb.j - jmin)/2 + 1;
  map_num[size * count + offset] = num;
  map_index[size * count + offset] = length;
  length += num;
}

void Map_4_count_4(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count, four_body &fb)
{
  State tb;
  int jmin, num;
  cross_state(Parameters, Space, fb.v2, fb.v4, fb.v3, fb.v1, jmin, tb); // 2_4
  num = (tb.j - jmin)/2 + 1;
  map_num[size * count + offset] = num;
  map_index[size * count + offset] = length;
  length += num;
}

void Map_4_count_5678(Input_Parameters &Parameters, Model_Space &Space, int *map_index, int *map_num, int size, int offset, int &length, int &count)
{
  map_num[size * count + offset] = 1;	                  // 3_1, 3_2, 3_3, 3_4
  map_index[size * count + offset] = length;
  length += 1;
}

void Map_4_1(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  State tb;
  int chan, ind, jmin, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  int num = map_num[size * count + offset];
  for(int n = 0; n < num; ++n){
    cross_state(Parameters, Space, fb[count].v1, fb[count].v4, fb[count].v3, fb[count].v2, jmin, tb);
    tb.j -= 2*n;
    chan = Space.ind_2b_cross(Parameters.basis, tb);
    key1 = map1[chan][Space.hash2(fb[count].v1, fb[count].v4, tb.j)];
    key2 = map2[chan][Space.hash2(fb[count].v3, fb[count].v2, tb.j)];
    if( fb_j[count].v1 == 0 ){ X = 1.0; }
    else{ X = -1.0 * CGC6(fb_j[count].v1, fb_j[count].v2, J[count], fb_j[count].v3, fb_j[count].v4, tb.j); }
    map_fac1[index + n] = (J[count] + 1.0) * X;
    map_fac2[index + n] = (tb.j + 1.0) * X;
    ind = key1 * num2[chan] + key2;
    map_chan[index + n] = chan;
    map_ind[index + n] = ind;
  }
}

void Map_4_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  State tb;
  int chan, ind, jmin, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  int num = map_num[size * count + offset];
  for(int n = 0; n < num; ++n){
    cross_state(Parameters, Space, fb[count].v2, fb[count].v3, fb[count].v4, fb[count].v1, jmin, tb);
    tb.j -= 2*n;
    chan = Space.ind_2b_cross(Parameters.basis, tb);
    key1 = map1[chan][Space.hash2(fb[count].v2, fb[count].v3, tb.j)];
    key2 = map2[chan][Space.hash2(fb[count].v4, fb[count].v1, tb.j)];
    if( fb_j[count].v2 == 0 ){ X = 1.0; }
    else{ X = -1.0 * phase2(fb_j[count].v1 + fb_j[count].v2 + fb_j[count].v3 + fb_j[count].v4) * CGC6(fb_j[count].v2, fb_j[count].v1, J[count], fb_j[count].v4, fb_j[count].v3, tb.j); }
    map_fac1[index + n] = (J[count] + 1.0) * X;
    map_fac2[index + n] = (tb.j + 1.0) * X;
    ind = key1 * num2[chan] + key2;
    map_chan[index + n] = chan;
    map_ind[index + n] = ind;
  }
}

void Map_4_3(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  State tb;
  int chan, ind, jmin, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  int num = map_num[size * count + offset];
  for(int n = 0; n < num; ++n){
    cross_state(Parameters, Space, fb[count].v1, fb[count].v3, fb[count].v4, fb[count].v2, jmin, tb);
    tb.j -= 2*n;
    chan = Space.ind_2b_cross(Parameters.basis, tb);
    key1 = map1[chan][Space.hash2(fb[count].v1, fb[count].v3, tb.j)];
    key2 = map2[chan][Space.hash2(fb[count].v4, fb[count].v2, tb.j)];
    if( fb_j[count].v3 == 0 ){ X = 1.0; }
    else{ X = phase2(fb_j[count].v3 + fb_j[count].v4 - J[count]) * CGC6(fb_j[count].v1, fb_j[count].v2, J[count], fb_j[count].v4, fb_j[count].v3, tb.j); }
    map_fac1[index + n] = (J[count] + 1.0) * X;
    map_fac2[index + n] = (tb.j + 1.0) * X;
    ind = key1 * num2[chan] + key2;
    map_chan[index + n] = chan;
    map_ind[index + n] = ind;
  }
}

void Map_4_4(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  State tb;
  int chan, ind, jmin, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  int num = map_num[size * count + offset];
  for(int n = 0; n < num; ++n){
    cross_state(Parameters, Space, fb[count].v2, fb[count].v4, fb[count].v3, fb[count].v1, jmin, tb);
    tb.j -= 2*n;
    chan = Space.ind_2b_cross(Parameters.basis, tb);
    key1 = map1[chan][Space.hash2(fb[count].v2, fb[count].v4, tb.j)];
    key2 = map2[chan][Space.hash2(fb[count].v3, fb[count].v1, tb.j)];
    if( fb_j[count].v4 == 0 ){ X = 1.0; }
    else{ X = phase2(fb_j[count].v1 + fb_j[count].v2 - J[count]) * CGC6(fb_j[count].v2, fb_j[count].v1, J[count], fb_j[count].v3, fb_j[count].v4, tb.j); }
    map_fac1[index + n] = (J[count] + 1.0) * X;
    map_fac2[index + n] = (tb.j + 1.0) * X;
    ind = key1 * num2[chan] + key2;
    map_chan[index + n] = chan;
    map_ind[index + n] = ind;
  }
}

void Map_4_5(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  int chan, ind, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  chan = Space.ind_1b(Parameters.basis, Space.qnums[fb[count].v1]);
  key1 = map1[chan][fb[count].v1];
  key2 = map2[chan][Space.hash3(fb[count].v3, fb[count].v4, fb[count].v2, J[count])];
  if( fb_j[count].v1 == 0.0 ){ X = 1.0; }
  else{ X = phase2(fb_j[count].v1 + fb_j[count].v2 - J[count]) * std::sqrt((J[count] + 1.0)/(fb_j[count].v1 + 1.0)); }
  map_fac1[index] = X;
  map_fac2[index] = 1.0 / X;
  ind = key1 * num2[chan] + key2;
  map_chan[index] = chan;
  map_ind[index] = ind;
}

void Map_4_6(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  int chan, ind, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  chan = Space.ind_1b(Parameters.basis, Space.qnums[fb[count].v2]);
  key1 = map1[chan][fb[count].v2];
  key2 = map2[chan][Space.hash3(fb[count].v3, fb[count].v4, fb[count].v1, J[count])];
  if( fb_j[count].v2 == 0 ){ X = 1.0; }
  else{ X = -1.0 * std::sqrt((J[count] + 1.0)/(fb_j[count].v2 + 1.0)); }
  map_fac1[index] = X;
  map_fac2[index] = 1.0 / X;
  ind = key1 * num2[chan] + key2;
  map_chan[index] = chan;
  map_ind[index] = ind;
}

void Map_4_7(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  int chan, ind, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  chan = Space.ind_1b(Parameters.basis, Space.qnums[fb[count].v3]);
  key2 = map2[chan][fb[count].v3];
  key1 = map1[chan][Space.hash3(fb[count].v1, fb[count].v2, fb[count].v4, J[count])];
  if( fb_j[count].v3 == 0 ){ X = 1.0; }
  else{ X = phase2(fb_j[count].v3 + fb_j[count].v4 - J[count]) * std::sqrt((J[count] + 1.0)/(fb_j[count].v3 + 1.0)); }
  map_fac1[index] = X;
  map_fac2[index] = 1.0 / X;
  ind = key1 * num2[chan] + key2;
  map_chan[index] = chan;
  map_ind[index] = ind;
}

void Map_4_8(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int *map_index, int *map_num, int size, int offset, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, four_body *fb, four_body *fb_j, int *J)
{
  int chan, ind, key1, key2;
  double X;
  int index = map_index[size * count + offset];
  chan = Space.ind_1b(Parameters.basis, Space.qnums[fb[count].v4]);
  key2 = map2[chan][fb[count].v4];
  key1 = map1[chan][Space.hash3(fb[count].v1, fb[count].v2, fb[count].v3, J[count])];
  if( fb_j[count].v4 == 0 ){ X = 1.0; }
  else{ X = -1.0 * std::sqrt((J[count] + 1.0)/(fb_j[count].v4 + 1.0)); }
  map_fac1[index] = X;
  map_fac2[index] = 1.0 / X;
  ind = key1 * num2[chan] + key2;
  map_chan[index] = chan;
  map_ind[index] = ind;
}

void Map_2(Input_Parameters &Parameters, Model_Space &Space, int *map_chan, int *map_ind, double *map_fac1, double *map_fac2, int &count, std::unordered_map<int,int> *map1, std::unordered_map<int,int> *map2, int *num2, two_body &tb, int &J)
{
  int chan, ind, key1, key2;
  double X;
  chan = Space.ind_1b(Parameters.basis, Space.qnums[tb.v1]);
  key1 = map1[chan][tb.v1];
  key2 = map2[chan][tb.v2];
  ind = key1 * num2[chan] + key2;
  X = std::sqrt(J + 1.0);
  map_chan[count] = chan;
  map_ind[count] = ind;
  map_fac1[count] = 1.0/X;
  map_fac2[count] = X;
}

void direct_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &q, int &r, int &s, int &jmin1, State &tb1)
{
  plus(tb1, Space.qnums[p], Space.qnums[q]);
  jmin1 = abs(Space.qnums[p].j - Space.qnums[q].j);
  if(abs(Space.qnums[r].j - Space.qnums[s].j) > jmin1){ jmin1 = abs(Space.qnums[r].j - Space.qnums[s].j); }
  if(Space.qnums[r].j + Space.qnums[s].j < tb1.j){ tb1.j = Space.qnums[r].j + Space.qnums[s].j; }
}

void cross_state(Input_Parameters &Parameters, Model_Space &Space, int &p, int &s, int &r, int &q, int &jmin2, State &tb2)
{
  minus(tb2, Space.qnums[p], Space.qnums[s]);
  jmin2 = abs(Space.qnums[p].j - Space.qnums[s].j);
  if(abs(Space.qnums[r].j - Space.qnums[q].j) > jmin2){ jmin2 = abs(Space.qnums[r].j - Space.qnums[q].j); }
  if(Space.qnums[r].j + Space.qnums[q].j < tb2.j){ tb2.j = Space.qnums[r].j + Space.qnums[q].j; }
}
