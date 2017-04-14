#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

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


void Model_Space::delete_struct(Input_Parameters &Parameters)
{
  delete[] qnums;
  std::cout << "delete space: " << Parameters.basis << " " << indtotj << std::endl;
  if(Parameters.basis == "infinite"){ delete[] map_2b; }
  if(indtotj != 0){
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
  Space.indtotj = 0;

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
  Space.indtotj = 0;

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
  Space.indtotj = 0;

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
