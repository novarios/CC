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
int Model_Space::hash2(int &p, int &q, int &j){
  return (j * std::pow(num_states, 2)) + (num_states * p + q);
}

//   Function to return Hash index for 3 indices
int Model_Space::hash3(int &p, int &q, int &r, int &j){
  return (j * std::pow(num_states, 3)) + (num_states * num_states * p + num_states * q + r);
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
    for(int i = 0; i < num_jstates; ++i){ delete[] shellsm[i]; }
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
  int key;

  fullpath = PATH + Parameters.LevelScheme + ".sp";
  splevels.open(fullpath.c_str());
  if (!splevels.is_open()){ std::cerr << "Level Scheme file does not exist" << std::endl; exit(1); };

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
      //std::cout << "! " << ind << " " << Space.qnums[i].n << " " << Space.qnums[i].j << " " << Space.qnums[i].m << "  " << energy << std::endl;
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
      //std::cout << "! " << ind << " " << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << "  " << energy << std::endl;
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
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }
  splevels.close();

  // Deterime Shell Structure
  Space.Determine_Shells(Parameters);

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
    std::cout << "###   " << i << ", " << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << " " << Space.qnums[i].t << ", " << Space.qnums[i].energy << std::endl;
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
  // reset Space.indtot with degeneracies 2j + 1
  Space.num_states = 0;
  for(int i = 0; i < indtotj; ++i){
    Space.num_states += Space.qnums[i].j + 1;
    states[i] = Space.qnums[i];
  }
  delete[] Space.qnums;
  Space.qnums = new State[Space.num_states];

  // initialize mins and maxs for each quantum number
  Space.qmins.maximize();
  Space.qmins1.maximize();
  Space.qmins2.maximize();
  Space.qmaxs.minimize();
  Space.qmaxs1.minimize();
  Space.qmaxs2.minimize();

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

  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }

  // Deterime Shell Structure
  Space.Determine_Shells(Parameters);

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

  // Find chan3, chan1, and chan2 sizes
  Space.size_1b = Space.qsizes.par * Space.qsizes.m * Space.qsizes.t;
  Space.size_2b_dir = Space.qsizes1.par * Space.qsizes1.m * Space.qsizes1.t;
  Space.size_2b_cross = Space.qsizes2.par * Space.qsizes2.m * Space.qsizes2.t;
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
		if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++Space.num_hol; ++Parameters.P; }
		else{ Space.qnums[count].type = "particle"; ++Space.num_par; }
	      }
	      if(tz == 1){
		if(Parameters.Nshells == 0){ continue; }
		++Space.num_n;
		if(shell < Parameters.Nshells){ Space.qnums[count].type = "hole"; ++Space.num_hol; ++Parameters.N; }
		else{ Space.qnums[count].type = "particle"; ++Space.num_par; }
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
  Space.qmins.nx = -Space.nmax;
  Space.qmaxs.nx = Space.nmax;
  Space.qmins.ny = -Space.nmax;
  Space.qmaxs.ny = Space.nmax;
  Space.qmins.nz = -Space.nmax;
  Space.qmaxs.nz = Space.nmax;
  Space.qmins.m = -1;
  Space.qmaxs.m = 1;

  for(int i = 0; i < Space.num_states; ++i){
    key = Space.ind_state(Parameters.basis, Space.qnums[i]);
    Space.map_state[key] = i;
  }
  
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

  // Calculate energies
  double L = pow(Space.num_hol/Parameters.density, 1.0/3.0);
  if(Parameters.calc_case == "nuclear"){
    for(int i = 0; i < Space.num_states; ++i){
      if(Space.qnums[i].t == -1){ Space.qnums[i].energy *= proton_prefac * M_PI*M_PI / (L*L); }
      else if(Space.qnums[i].t == 1){ Space.qnums[i].energy *= neutron_prefac * M_PI*M_PI / (L*L); }
    }
    // Change energies to Hartree-Fock energies, E_p = E_p + V_pipi
    for(int p = 0; p < Space.num_states; ++p){
      for(int i = 0; i < Space.num_hol; ++i){
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
      for(int i = 0; i < Space.num_hol; ++i){
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
	  Space.qnums[count].t = -1;
	  Space.qnums[count].nx = 0;
	  Space.qnums[count].ny = 0;
	  Space.qnums[count].nz = 0;
	  Space.qnums[count].j = 0;
	  if(shell < Parameters.Pshells){ Space.qnums[count].type = "hole"; ++Space.num_hol; ++Parameters.P; }
	  else{ Space.qnums[count].type = "particle"; ++Space.num_par; }
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
  int key, jmin;
  int chan1, chan2, chan3;
  int nh0, np0, nhh0, npp0, nhp0, nhp10, nph10, nhh10, npp10, npph0, nhhp0, nhpp0, nhph0, nhhh0, nppp0;
  int h_total, p_total;
  int hh_total, pp_total, hp1_total, ph1_total, pph_total, hhp_total;
  int hp_total, hh1_total, pp1_total, hph_total, hpp_total, hhh_total, ppp_total;
  int h1, h2, h3, p1, p2, p3;
  one_body ob;
  two_body tb;
  three_body thb;

  size1 = Space.size_2b_dir;
  size2 = Space.size_2b_cross;
  size3 = Space.size_1b;

  // allocate memory for Hvec and Pvec, reset hnum and pnum
  h_total = 0;
  p_total = 0;
  qnums3 = new State[size3];
  nh = new int[size3];
  h_index = new int[size3];
  h_map = new std::unordered_map<int,int>[size3];
  np = new int[size3];
  p_index = new int[size3];
  p_map = new std::unordered_map<int,int>[size3];
  for(chan3 = 0; chan3 < size3; ++chan3){
    nh[chan3] = 0;
    np[chan3] = 0;
  }
  for(int p = 0; p < Space.num_states; ++p){
    chan3 = Space.ind_1b(Parameters.basis, Space.qnums[p]);
    qnums3[chan3] = Space.qnums[p];
    if(Space.qnums[p].type == "hole"){
      ++nh[chan3];
      ++h_total;
    }
    else{
      ++np[chan3];
      ++p_total;
    }
  }

  h_vec = new one_body[h_total];
  p_vec = new one_body[p_total];
  h_total = 0;
  p_total = 0;
  for(chan3 = 0; chan3 < size3; ++chan3){
    h_index[chan3] = h_total;
    h_total += nh[chan3];
    nh[chan3] = 0;
    p_index[chan3] = p_total;
    p_total += np[chan3];
    np[chan3] = 0;
  }
  for(int p = 0; p < Space.num_states; ++p){
    chan3 = Space.ind_1b(Parameters.basis, Space.qnums[p]);
    qnums3[chan3] = Space.qnums[p];
    ob.v1 = p;
    if(Space.qnums[p].type == "hole"){
      nh0 = nh[chan3];
      h_vec[h_index[chan3] + nh0] = ob;
      h_map[chan3][p] = nh0;
      ++nh[chan3];
    }
    else{
      np0 = np[chan3];
      p_vec[p_index[chan3] + np0] = ob;
      p_map[chan3][p] = np0;
      ++np[chan3];
    }
  }

  /*for(int i = 0; i < size3; ++i){
    std::cout << "Chan3: " << i << ", " << qnums3[i].par << " " << qnums3[i].ml << " " << qnums3[i].m << std::endl;
    std::cout << "nh = " << nh[i] << ", np = " << np[i] << std::endl;
    for(int j = 0; j < nh[i]; ++j){ std::cout << h_vec[i][j] << " "; }
    std::cout << std::endl;
    for(int j = 0; j < np[i]; ++j){ std::cout << p_vec[i][j] << " "; }
    std::cout << std::endl;
    }*/

  size1 = Space.size_2b_dir;
  size2 = Space.size_2b_cross;
  qnums1 = new State[size1];
  qnums2 = new State[size2];
  //std::cout << " Size1 = " << size1 << ", Size2 = " << size2 << ", Size3 = " << size3 << std::endl;

  // For CCD
  hh_total = 0;
  nhh = new int[size1];
  hh_index = new int[size1];
  hh_map = new std::unordered_map<int,int>[size1];
  pp_total = 0;
  npp = new int[size1];
  pp_index = new int[size1];
  pp_map = new std::unordered_map<int,int>[size1];
  hp1_total = 0;
  nhp1 = new int[size2];
  hp1_index = new int[size2];
  hp1_map = new std::unordered_map<int,int>[size2];
  ph1_total = 0;
  nph1 = new int[size2];
  ph1_index = new int[size2];
  ph1_map = new std::unordered_map<int,int>[size2];
  for(chan1 = 0; chan1 < size1; ++chan1){
    nhh[chan1] = 0;
    npp[chan1] = 0;
  }
  for(chan2 = 0; chan2 < size2; ++chan2){
    nhp1[chan2] = 0;
    nph1[chan2] = 0;
  }
  // count nhh, npp, nhp1, nph1
  for(int p = 0; p < Space.num_states; ++p){
    for(int q = 0; q < Space.num_states; ++q){
      plus(state1, Space.qnums[p], Space.qnums[q]);
      minus(state2, Space.qnums[p], Space.qnums[q]);
      jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
      while(state1.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	qnums1[chan1] = state1;
	chan2 = Space.ind_2b_cross(Parameters.basis, state2);
	qnums2[chan2] = state2;
	if(p >= Space.num_hol && q >= Space.num_hol){ // pp
	  ++npp[chan1];
	  ++pp_total;
	}
	else if(p < Space.num_hol && q >= Space.num_hol){ // hp1
	  ++nhp1[chan2];
	  ++hp1_total;
	}
	else if(p >= Space.num_hol && q < Space.num_hol){ // ph1
	  ++nph1[chan2];
	  ++ph1_total;
	}
	else{ // hh
	  ++nhh[chan1];
	  ++hh_total;
	}
	state1.j -= 2;
	state2.j -= 2;
      }
    }
  }

  hh_vec = new two_body[hh_total];
  hh_total = 0;
  pp_vec = new two_body[pp_total];
  pp_total = 0;
  hp1_vec = new two_body[hp1_total];
  hp1_total = 0;
  ph1_vec = new two_body[ph1_total];
  ph1_total = 0;
  for(chan1 = 0; chan1 < size1; ++chan1){
    hh_index[chan1] = hh_total;
    hh_total += nhh[chan1];
    nhh[chan1] = 0;
    pp_index[chan1] = pp_total;
    pp_total += npp[chan1];
    npp[chan1] = 0;
  }
  for(chan2 = 0; chan2 < size2; ++chan2){
    hp1_index[chan2] = hp1_total;
    hp1_total += nhp1[chan2];
    nhp1[chan2] = 0;
    ph1_index[chan2] = ph1_total;
    ph1_total += nph1[chan2];
    nph1[chan2] = 0;
  }

  for(int p = 0; p < Space.num_states; ++p){
    for(int q = 0; q < Space.num_states; ++q){
      plus(state1, Space.qnums[p], Space.qnums[q]);
      minus(state2, Space.qnums[p], Space.qnums[q]);
      jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
      tb.v1 = p;
      tb.v2 = q;
      while(state1.j >= jmin){
	key = Space.hash2(p, q, state1.j);
	chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	chan2 = Space.ind_2b_cross(Parameters.basis, state2);
	if(p >= Space.num_hol && q >= Space.num_hol){ // pp
	  npp0 = npp[chan1];
	  pp_vec[pp_index[chan1] + npp0] = tb;
	  pp_map[chan1][key] = npp0;
	  ++npp[chan1];
	}
	else if(p < Space.num_hol && q >= Space.num_hol){ // hp1
	  nhp10 = nhp1[chan2];
	  hp1_vec[hp1_index[chan2] + nhp10] = tb;
	  hp1_map[chan2][key] = nhp10;
	  ++nhp1[chan2];
	}
	else if(p >= Space.num_hol && q < Space.num_hol){ // ph1
	  nph10 = nph1[chan2];
	  ph1_vec[ph1_index[chan2] + nph10] = tb;
	  ph1_map[chan2][key] = nph10;
	  ++nph1[chan2];
	}
	else{ // hh
	  nhh0 = nhh[chan1];
	  hh_vec[hh_index[chan1] + nhh0] = tb;
	  hh_map[chan1][key] = nhh0;
	  ++nhh[chan1];
	}
	state1.j -= 2;
	state2.j -= 2;
      }
    }
  }

  /*for(int i = 0; i < size1; ++i){
    std::cout << "Chan1:  " << i << " " << qnums1[i].ml << " " << qnums1[i].m << std::endl;
    std::cout << "nhh = " << nhh[i] << ", npp = " << npp[i] << std::endl;
    for(int j = 0; j < nhh[i]; ++j){
      std::cout << hh_vec[i][2*j] << "," << hh_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    for(int j = 0; j < npp[i]; ++j){
      std::cout << pp_vec[i][2*j] << "," << pp_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
  }

  for(int i = 0; i < size2; ++i){
    std::cout << "Chan2:  " << i << " " << qnums2[i].ml << " " << qnums2[i].m << std::endl;
    std::cout << "nhp1 = " << nhp1[i] << ", nph1 = " << nph1[i] << std::endl;
    for(int j = 0; j < nhp1[i]; ++j){
      std::cout << hp1_vec[i][2*j] << "," << hp1_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    for(int j = 0; j < nph1[i]; ++j){
      std::cout << ph1_vec[i][2*j] << "," << ph1_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    }*/

  hhp_total = 0;
  nhhp = new int[size3];
  hhp_index = new int[size3];
  hhp_map = new std::unordered_map<int,int>[size3];
  pph_total = 0;
  npph = new int[size3];
  pph_index = new int[size3];
  pph_map = new std::unordered_map<int,int>[size3];
  for(chan3 = 0; chan3 < size3; ++chan3){
    nhhp[chan3] = 0;
    npph[chan3] = 0;
  }
  for(chan3 = 0; chan3 < size3; ++chan3){
    for(int p = 0; p < size3; ++p){
      nh0 = nh[p];
      np0 = np[p];
      plus(state1, qnums3[chan3], qnums3[p]);
      jmin = abs(qnums3[chan3].j - qnums3[p].j);
      while(state1.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	nhhp[chan3] += nhh[chan1] * np0; // <pp|hh> -> <p|hhp'>
	hhp_total += nhh[chan1] * np0;
	npph[chan3] += npp[chan1] * nh0; // <pp|hh> -> <pph'|h>
	pph_total += npp[chan1] * nh0;
	state1.j -= 2;
      }
    }
  }

  hhp_vec = new three_body[hhp_total];
  hhp_total = 0;
  pph_vec = new three_body[pph_total];
  pph_total = 0;
  for(chan3 = 0; chan3 < size3; ++chan3){
    hhp_index[chan3] = hhp_total;
    hhp_total += nhhp[chan3];
    nhhp[chan3] = 0;
    pph_index[chan3] = pph_total;
    pph_total += npph[chan3];
    npph[chan3] = 0;
  }
  for(chan3 = 0; chan3 < size3; ++chan3){
    for(int p = 0; p < size3; ++p){
      nh0 = nh[p];
      np0 = np[p];
      plus(state1, qnums3[chan3], qnums3[p]);
      jmin = abs(qnums3[chan3].j - qnums3[p].j);
      while(state1.j >= jmin){
	chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	nhh0 = nhh[chan1];
	npp0 = npp[chan1];
	for(int hh = 0; hh < nhh0; ++hh){ // <pp|hh> -> <p|hhp'>
	  h1 = hh_vec[hh_index[chan1] + hh].v1;
	  h2 = hh_vec[hh_index[chan1] + hh].v2;
	  for(int q = 0; q < np0; ++q){
	    p1 = p_vec[p_index[p] + q].v1;
	    nhhp0 = nhhp[chan3];
	    key = Space.hash3(h1, h2, p1, state1.j);
	    thb.v1 = h1;
	    thb.v2 = h2;
	    thb.v3 = p1;
	    hhp_vec[hhp_index[chan3] + nhhp0] = thb;
	    hhp_map[chan3][key] = nhhp0;
	    ++nhhp[chan3];
	  }
	}
	for(int pp = 0; pp < npp0; ++pp){ // <pp|hh> -> <pph'|h>  
	  p1 = pp_vec[pp_index[chan1] + pp].v1;
	  p2 = pp_vec[pp_index[chan1] + pp].v2;
	  for(int q = 0; q < nh0; ++q){
	    h1 = h_vec[h_index[p] + q].v1;
	    npph0 = npph[chan3];
	    key = Space.hash3(p1, p2, h1, state1.j);
	    thb.v1 = p1;
	    thb.v2 = p2;
	    thb.v3 = h1;
	    pph_vec[pph_index[chan3] + npph0] = thb;
	    pph_map[chan3][key] = npph0;
	    ++npph[chan3];
	  }
	}
	state1.j -= 2;
      }
    }
  }

  // For CCSD
  if(Parameters.approx == "singles"){
    hp_total = 0;
    nhp = new int[size1];
    hp_index = new int[size1];
    hp_map = new std::unordered_map<int,int>[size1];
    hh1_total = 0;
    nhh1 = new int[size2];
    hh1_index = new int[size2];
    hh1_map = new std::unordered_map<int,int>[size2];
    pp1_total = 0;
    npp1 = new int[size2];
    pp1_index = new int[size2];
    pp1_map = new std::unordered_map<int,int>[size2];
    for(chan1 = 0; chan1 < size1; ++chan1){
      nhp[chan1] = 0;
    }
    for(chan2 = 0; chan2 < size2; ++chan2){
      nhh1[chan2] = 0;
      npp1[chan2] = 0;
    }

    for(int p = 0; p < Space.num_states; ++p){
      for(int q = 0; q < Space.num_states; ++q){
	plus(state1, Space.qnums[p], Space.qnums[q]);
	minus(state2, Space.qnums[p], Space.qnums[q]);
	jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
	while(state1.j >= jmin){
	  chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	  qnums1[chan1] = state1;
	  chan2 = Space.ind_2b_cross(Parameters.basis, state2);
	  qnums2[chan2] = state2;
	  if(p >= Space.num_hol && q >= Space.num_hol){ // pp1
	    ++npp1[chan2];
	    ++pp1_total;
	  }
	  else if(p < Space.num_hol && q >= Space.num_hol){ // hp
	    ++nhp[chan1];
	    ++hp_total;
	  }
	  else if(p < Space.num_hol && q < Space.num_hol){ // hh1
	    ++nhh1[chan2];
	    ++hh1_total;
	  }
	  state1.j -= 2;
	  state2.j -= 2;
	}
      }
    }

    hp_vec = new two_body[hp_total];
    hp_total = 0;
    hh1_vec = new two_body[hh1_total];
    hh1_total = 0;
    pp1_vec = new two_body[pp1_total];
    pp1_total = 0;
    for(chan1 = 0; chan1 < size1; ++chan1){
      hp_index[chan1] = hp_total;
      hp_total += nhp[chan1];
      nhp[chan1] = 0;
    }
    for(chan2 = 0; chan2 < size2; ++chan2){
      hh1_index[chan2] = hh1_total;
      hh1_total += nhh1[chan2];
      nhh1[chan2] = 0;
      pp1_index[chan2] = pp1_total;
      pp1_total += npp1[chan2];
      npp1[chan2] = 0;
    }

    for(int p = 0; p < Space.num_states; ++p){
      for(int q = 0; q < Space.num_states; ++q){
	plus(state1, Space.qnums[p], Space.qnums[q]);
	minus(state2, Space.qnums[p], Space.qnums[q]);
	jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
	tb.v1 = p;
	tb.v2 = q;
	while(state1.j >= jmin){
	  key = Space.hash2(p, q, state1.j);
	  chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	  chan2 = Space.ind_2b_cross(Parameters.basis, state2);
	  if(p >= Space.num_hol && q >= Space.num_hol){ // pp1
	    npp10 = npp1[chan2];
	    pp1_vec[pp1_index[chan2] + npp10] = tb;
	    pp1_map[chan2][key] = npp10;
	    ++npp1[chan2];
	  }
	  else if(p < Space.num_hol && q >= Space.num_hol){ // hp
	    nhp0 = nhp[chan1];
	    hp_vec[hp_index[chan1] + nhp0] = tb;
	    hp_map[chan1][key] = nhp0;
	    ++nhp[chan1];
	  }
	  else if(p < Space.num_hol && q < Space.num_hol){ // hh1
	    nhh10 = nhh1[chan2];
	    hh1_vec[hh1_index[chan2] + nhh10] = tb;
	    hh1_map[chan2][key] = nhh10;
	    ++nhh1[chan2];
	  }
	  state1.j -= 2;
	  state2.j -= 2;
	}
      }
    }

    /*for(int i = 0; i < size1; ++i){
      std::cout << "Chan1:  " << i << " " << qnums1[i].ml << " " << qnums1[i].m << std::endl;
      std::cout << "nhp = " << nhp[i] << std::endl;
      for(int j = 0; j < nhp[i]; ++j){
	std::cout << hp_vec[i][2*j] << "," << hp_vec[i][2*j + 1] << " ";
      }
      std::cout << std::endl;
    }
    
    for(int i = 0; i < size2; ++i){
      std::cout << "Chan2:  " << i << " " << qnums2[i].ml << " " << qnums2[i].m << std::endl;
      std::cout << "nhh1 = " << nhh1[i] << ", npp1 = " << npp1[i] << std::endl;
      for(int j = 0; j < nhh1[i]; ++j){
	std::cout << hh1_vec[i][2*j] << "," << hh1_vec[i][2*j + 1] << " ";
      }
      std::cout << std::endl;
      for(int j = 0; j < npp1[i]; ++j){
	std::cout << pp1_vec[i][2*j] << "," << pp1_vec[i][2*j + 1] << " ";
      }
      std::cout << std::endl;
      }*/

    hph_total = 0;
    nhph = new int[size3];
    hph_index = new int[size3];
    hph_map = new std::unordered_map<int,int>[size3];
    hpp_total = 0;
    nhpp = new int[size3];
    hpp_index = new int[size3];
    hpp_map = new std::unordered_map<int,int>[size3];
    hhh_total = 0;
    nhhh = new int[size3];
    hhh_index = new int[size3];
    hhh_map = new std::unordered_map<int,int>[size3];
    ppp_total = 0;
    nppp = new int[size3];
    ppp_index = new int[size3];
    ppp_map = new std::unordered_map<int,int>[size3];
    for(chan3 = 0; chan3 < size3; ++chan3){
      nhph[chan3] = 0;
      nhpp[chan3] = 0;
      nhhh[chan3] = 0;
      nppp[chan3] = 0;
    }
    for(chan3 = 0; chan3 < size3; ++chan3){
      for(int p = 0; p < size3; ++p){
	nh0 = nh[p];
	np0 = np[p];
	plus(state1, qnums3[chan3], qnums3[p]);
	jmin = abs(qnums3[chan3].j - qnums3[p].j);
	while(state1.j >= jmin){
	  chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	  nhph[chan3] += nhp[chan1] * nh0; // <hh|hp> -> <h|hph'>
	  hph_total += nhp[chan1] * nh0;
	  nhpp[chan3] += nhp[chan1] * np0; // <hp|pp> -> <hpp'|p>
	  hpp_total += nhp[chan1] * np0;
	  nhhh[chan3] += nhh[chan1] * nh0; // <hh|hp> -> <hhh'|p>
	  hhh_total += nhh[chan1] * nh0;
	  nppp[chan3] += npp[chan1] * np0; // <hp|pp> -> <h|ppp'>
	  ppp_total += npp[chan1] * np0;
	  state1.j -= 2;
	}
      }
    }
    hph_vec = new three_body[hph_total];
    hph_total = 0;
    hpp_vec = new three_body[hpp_total];
    hpp_total = 0;
    hhh_vec = new three_body[hhh_total];
    hhh_total = 0;
    ppp_vec = new three_body[ppp_total];
    ppp_total = 0;
    for(chan3 = 0; chan3 < size3; ++chan3){
      hph_index[chan3] = hph_total;
      hph_total += nhph[chan3];
      nhph[chan3] = 0;
      hpp_index[chan3] = hpp_total;
      hpp_total += nhpp[chan3];
      nhpp[chan3] = 0;
      hhh_index[chan3] = hhh_total;
      hhh_total += nhhh[chan3];
      nhhh[chan3] = 0;
      ppp_index[chan3] = ppp_total;
      ppp_total += nppp[chan3];
      nppp[chan3] = 0;
    }

    for(chan3 = 0; chan3 < size3; ++chan3){
      for(int p = 0; p < size3; ++p){
	nh0 = nh[p];
	np0 = np[p];
	plus(state1, qnums3[chan3], qnums3[p]);
	jmin = abs(qnums3[chan3].j - qnums3[p].j);
	while(state1.j >= jmin){
	  chan1 = Space.ind_2b_dir(Parameters.basis, state1);
	  nhh0 = nhh[chan1];
	  npp0 = npp[chan1];
	  nhp0 = nhp[chan1];
	  for(int hp = 0; hp < nhp0; ++hp){ // <hh|hp> -> <h|hph'>
	    h1 = hp_vec[hp_index[chan1] + hp].v1;
	    p1 = hp_vec[hp_index[chan1] + hp].v2;
	    for(int q = 0; q < nh0; ++q){
	      h2 = h_vec[h_index[p] + q].v1;
	      nhph0 = nhph[chan3];
	      key = Space.hash3(h1, p1, h2, state1.j);
	      thb.v1 = h1;
	      thb.v2 = p1;
	      thb.v3 = h2;
	      hph_vec[hph_index[chan3] + nhph0] = thb;
	      hph_map[chan3][key] = nhph0;
	      ++nhph[chan3];
	    }
	  }
	  for(int hp = 0; hp < nhp0; ++hp){ // <hp|pp> -> <hpp'|p>  
	    h1 = hp_vec[hp_index[chan1] + hp].v1;
	    p1 = hp_vec[hp_index[chan1] + hp].v2;
	    for(int q = 0; q < np0; ++q){
	      p2 = p_vec[p_index[p] + q].v1;
	      nhpp0 = nhpp[chan3];
	      key = Space.hash3(h1, p1, p2, state1.j);
	      thb.v1 = h1;
	      thb.v2 = p1;
	      thb.v3 = p2;
	      hpp_vec[hpp_index[chan3] + nhpp0] = thb;
	      hpp_map[chan3][key] = nhpp0;
	      ++nhpp[chan3];
	    }
	  }
	  for(int hh = 0; hh < nhh0; ++hh){ // <hh|hp> -> <hhh'|p>
	    h1 = hh_vec[hh_index[chan1] + hh].v1;
	    h2 = hh_vec[hh_index[chan1] + hh].v2;
	    for(int q = 0; q < nh0; ++q){
	      h3 = h_vec[h_index[p] + q].v1;
	      nhhh0 = nhhh[chan3];
	      key = Space.hash3(h1, h2, h3, state1.j);
	      thb.v1 = h1;
	      thb.v2 = h2;
	      thb.v3 = h3;
	      hhh_vec[hhh_index[chan3] + nhhh0] = thb;
	      hhh_map[chan3][key] = nhhh0;
	      ++nhhh[chan3];
	    }
	  }
	  for(int pp = 0; pp < npp0; ++pp){ // <hp|pp> -> <h|ppp'>
	    p1 = pp_vec[pp_index[chan1] + pp].v1;
	    p2 = pp_vec[pp_index[chan1] + pp].v2;
	    for(int q = 0; q < np0; ++q){
	      p3 = p_vec[p_index[p] + q].v1;
	      nppp0 = nppp[chan3];
	      key = Space.hash3(p1, p2, p3, state1.j);
	      thb.v1 = p1;
	      thb.v2 = p2;
	      thb.v3 = p3;
	      ppp_vec[ppp_index[chan3] + nppp0] = thb;
	      ppp_map[chan3][key] = nppp0;
	      ++nppp[chan3];
	    }
	  }
	  state1.j -= 2;
	}
      }
    }

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

  /*for(int i = 0; i < size1; ++i){
    std::cout << "Chan1:  " << i << " " << qnums1[i].par << " " << qnums1[i].t << " " << qnums1[i].j << std::endl;
    for(int j = 0; j < nhh[i]; ++j){
      std::cout << hh_vec[i][2*j] << "," << hh_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    for(int j = 0; j < npp[i]; ++j){
      std::cout << pp_vec[i][2*j] << "," << pp_vec[i][2*j + 1] << " ";
    }
    std::cout << std::endl;
    }*/
  
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
    if(qnums[i].energy > p_en && qnums[i].t == -1){
      p_en = qnums[i].energy;
      ++pshell_num;
    }
    else if(qnums[i].energy > n_en && qnums[i].t == 1){
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
    if(qnums[i].energy > p_en && qnums[i].t == -1){
      pshell[pshell_num] = i;
      p_en = qnums[i].energy;
      ++pshell_num;
    }
    else if(qnums[i].energy > n_en && qnums[i].t == 1){
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
      if(i < pshell[Parameters.Pshells]){ qnums[i].type = "hole"; ++num_hol; ++Parameters.P; }
      else{ qnums[i].type = "particle"; ++num_par; }
    }
    else if(qnums[i].t == 1){
      ++num_n;
      if(i < nshell[Parameters.Nshells]){ qnums[i].type = "hole"; ++num_hol; ++Parameters.N; }
      else{ qnums[i].type = "particle"; ++num_par; }
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
