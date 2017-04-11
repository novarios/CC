#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"

HF_Matrix_Elements::HF_Matrix_Elements(const HF_Channels &Chan)
{
  int ntb1;
  V = new double*[Chan.size1];
  for(int i = 0; i < Chan.size1; ++i){
    ntb1 = Chan.ntb[i];
    if(ntb1 != 0){
      ntb1 *= ntb1;
      V[i] = new double[ntb1];
      for(int j = 0; j < ntb1; ++j){
	V[i][j] = 0.0;
      }
    }
  }
}

void HF_Matrix_Elements::delete_struct(const HF_Channels &Chan)
{
  int ntb1;
  for(int i = 0; i < Chan.size1; ++i){
    ntb1 = Chan.ntb[i];
    if(ntb1 != 0){
      delete[] V[i];
    }
  }
  delete[] V;
}

HF_Channels::HF_Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  State state;
  int count0, ind1, nob1, ntb1, key;
  State *qnumstemp = new State[Space.indtot]; // max number of qnums groups
  int *obnum = new int[Space.indtot]; // count states in each qnums group
  indvec = new int[Space.indtot]; // index of qnums group for each state
  for(int i = 0; i < Space.indtot; ++i){ obnum[i] = 0; }

  qnumstemp[0] = Space.qnums[0];
  indvec[0] = 0;
  obnum[0] = 1;

  // count # of qnums groups and # of hs and ps in each qnums group, and fill indvec
  count0 = 1;
  for(int i = 1; i < Space.indtot; ++i){
    state = Space.qnums[i];
    for(int k = 0; k < count0; ++k){
      if( equal(state, qnumstemp[k]) ){
	indvec[i] = k;
	++obnum[k];
	goto stop;
      }
      else if(k == count0 - 1){
	qnumstemp[count0] = Space.qnums[i];
	indvec[i] = count0;
	++obnum[count0];
	++count0;
	break;
      }
    }
  stop:;
  }
  size3 = count0;

  // allocate memory for Hvec and Pvec, reset hnum and pnum
  qnums3 = new State[size3];
  obvec = new int*[size3];
  nob = new int[size3];
  ob_map = new std::unordered_map<int,int>[size3];
  for(int i = 0; i < size3; ++i){
    nob1 = obnum[i];
    qnums3[i] = qnumstemp[i];
    if(nob1 != 0){ obvec[i] = new int[nob1]; }
    nob[i] = 0;
  }
  delete[] qnumstemp;
  delete[] obnum;

  // place states in appropriate Hvec or Pvec position
  for(int i = 0; i < Space.indtot; ++i){
    for(int k = 0; k < size3; ++k){
      if( equal(Space.qnums[i], qnums3[k]) ){
	obvec[k][nob[k]] = i;
	ob_map[k][i] = nob[k];
	++nob[k];
	break;
      }
    }
  }

  size1 = Space.size_2b;
  qnums1 = new State[size1];
  tbvec = new int*[size1];
  tb_map = new std::unordered_map<int,int>[size1];
  ntb = new int[size1];
  for(int i = 0; i < size1; ++i){
    ntb[i] = 0;
  }

  int jmin;
  for(int p = 0; p < Space.indtot; ++p){
    for(int q = 0; q < Space.indtot; ++q){
      plus(state, Space.qnums[p], Space.qnums[q]);
      if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
	jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  qnums1[ind1] = state;
	  ++ntb[ind1];
	  state.j -= 2;
	}
      }
      else{
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	qnums1[ind1] = state;
	++ntb[ind1];
      }
    }
  }

  for(int i = 0; i < size1; ++i){
    ntb1 = ntb[i];
    if(ntb1 != 0){
      tbvec[i] = new int[2 * ntb[i]];
    }
    ntb[i] = 0;
  }

  for(int p = 0; p < Space.indtot; ++p){
    for(int q = 0; q < Space.indtot; ++q){
      key = Hash2(p, q, Space.indtot);
      plus(state, Space.qnums[p], Space.qnums[q]);
      if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
	jmin = abs(Space.qnums[p].j - Space.qnums[q].j);
	while(state.j >= jmin){
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  tbvec[ind1][2 * ntb[ind1]] = p;
	  tbvec[ind1][2 * ntb[ind1] + 1] = q;
	  tb_map[ind1][key] = ntb[ind1];
	  ++ntb[ind1];
	  state.j -= 2;
	}
      }
      else{
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	tbvec[ind1][2 * ntb[ind1]] = p;
	tbvec[ind1][2 * ntb[ind1] + 1] = q;
	tb_map[ind1][key] = ntb[ind1];
	++ntb[ind1];
      }
    }
  }

  double memory = 0.0;
  int doubsize = sizeof(double);
  for(int i = 0; i < size1; ++i){
    memory += doubsize * ntb[i] * ntb[i];
  }
}

void HF_Channels::delete_struct()
{
  int ntb1;
  for(int i = 0; i < size3; ++i){
    delete[] obvec[i];
  }
  delete[] indvec;
  delete[] qnums3;
  delete[] obvec;
  delete[] ob_map;
  delete[] nob;

  for(int i = 0; i < size1; ++i){
    ntb1 = ntb[i];
    if(ntb1 != 0){
      delete[] tbvec[i];
    }
  }
  delete[] qnums1;
  delete[] tbvec;
  delete[] tb_map;
  delete[] ntb;
}

Single_Particle_States::Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan)
{
  int ind; // count for filling states
  int nob1, h1, p1;
  double tempen, tempen1, vec1; // temp energy for ordering states
  if(Parameters.Pshells != 0){ hp = Parameters.P; }
  else{ hp = 0; }
  if(Parameters.Nshells != 0){ hn = Parameters.N; }
  else{ hn = 0; }

  vectors = new double**[Chan.size3];
  energies = new double*[Chan.size3];
  h = new int[Chan.size3];
  p = new int[Chan.size3];
  holes = new double**[Chan.size3];
  particles = new double**[Chan.size3];
  h_energies = new double*[Chan.size3];
  pt_energies = new double*[Chan.size3];

  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h[i] = 0;
    p[i] = 0;
    if(nob1 != 0){
      vectors[i] = new double*[nob1];
      energies[i] = new double[nob1];
      for(int j = 0; j < nob1; ++j){
	energies[i][j] = Space.qnums[Chan.obvec[i][j]].energy;
	vectors[i][j] = new double[nob1];
	if(Space.qnums[Chan.obvec[i][j]].type == "hole"){ ++h[i]; }
	else{ ++p[i]; }
	for(int k = 0; k < nob1; ++k){
	  if(j == k){ vectors[i][j][k] = 1.0; }
	  else{ vectors[i][j][k] = 0.0; }
	}
      }
    }
  }

  // Order states by energy
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nob[i] - 1; ++j){
      ind = j;
      tempen = energies[i][j];
      for(int k = j + 1; k < Chan.nob[i]; ++k){
	if(energies[i][k] < tempen){ tempen = energies[i][k]; ind = k; }
      }
      tempen1 = energies[i][j];
      energies[i][j] = energies[i][ind];
      energies[i][ind] = tempen1;
      for(int k = 0; k < Chan.nob[i]; ++k){
	vec1 = vectors[i][j][k];
	vectors[i][j][k] = vectors[i][ind][k];
	vectors[i][ind][k] = vec1;
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h1 = h[i];
    p1 = p[i];
    if(h1 != 0){
      holes[i] = new double*[h1];
      h_energies[i] = new double[h1];
      for(int j = 0; j < h1; ++j){
	h_energies[i][j] = energies[i][j];
	holes[i][j] = new double[nob1];
	for(int k = 0; k < nob1; ++k){
	  holes[i][j][k] = vectors[i][j][k];
	}
      }
    }
    if(p1 != 0){
      particles[i] = new double*[p1];
      pt_energies[i] = new double[p1];
      for(int j = 0; j < p1; ++j){
	pt_energies[i][j] = energies[i][j + h1];
	particles[i][j] = new double[nob1];
	for(int k = 0; k < nob1; ++k){
	  particles[i][j][k] = vectors[i][j + h1][k];
	}
      }
    }
  }
  //Separate States
  Separate(Chan);
}

void Single_Particle_States::delete_struct(const HF_Channels &Chan)
{
  int h1, p1, nob1;
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h1 = h[i];
    p1 = p[i];
    if(nob1 != 0){
      for(int j = 0; j < nob1; ++j){
	delete[] vectors[i][j];
      }
      delete[] vectors[i];
      delete[] energies[i];
    }
    if(h1 != 0){
      for(int j = 0; j < h1; ++j){
	delete[] holes[i][j];
      }
      delete[] holes[i];
      delete[] h_energies[i];
    }
    if(p1 != 0){
      for(int j = 0; j < p1; ++j){
	delete[] particles[i][j];
      }
      delete[] particles[i];
      delete[] pt_energies[i];
    }
  }
  delete[] vectors;
  delete[] energies;
  delete[] h;
  delete[] holes;
  delete[] h_energies;
  delete[] p;
  delete[] particles;
  delete[] pt_energies;
}


//ignores existing holes and particles and fills them from protons and neutrons
void Single_Particle_States::Separate(const HF_Channels &Chan)
{
  //define holes/particles
  int h1, p1, nob1;
  for(int i = 0; i < Chan.size3; ++i){
    h1 = h[i];
    p1 = p[i];
    nob1 = Chan.nob[i];
    if(h1 != 0){
      for(int j = 0; j < h1; ++j){
	delete[] holes[i][j];
      }
      delete[] holes[i];
      delete[] h_energies[i];
      holes[i] = new double*[h1];
      h_energies[i] = new double[h1];
      for(int j = 0; j < h1; ++j){
	holes[i][j] = new double[nob1];
	h_energies[i][j] = energies[i][j];
	for(int k = 0; k < nob1; ++k){
	  holes[i][j][k] = vectors[i][j][k];
	}
      }
    }
    if(p1 != 0){
      for(int j = 0; j < p1; ++j){
	delete[] particles[i][j];
      }
      delete[] particles[i];
      delete[] pt_energies[i];
      particles[i] = new double*[p1];
      pt_energies[i] = new double[p1];
      for(int j = 0; j < p1; ++j){
	particles[i][j] = new double[nob1];
	pt_energies[i][j] = energies[i][j + h1];
	for(int k = 0; k < nob1; ++k){
	  particles[i][j][k] = vectors[i][j + h1][k];
	}
      }
    }
  }
}

void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1 = PATH + Parameters.MatrixElements + ".int"; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME, hom, r2, p2; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, coupJ, coupT, par; // interaction file contents
  int ind1, ind, key1, key2, key3, key4;
  State tb;

  ME = HF_Matrix_Elements(Chan);

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> number;
  while (number != "Total"){ 
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
    std::istringstream(interactionline) >> coupT >> par >> coupJ >> shell1 >> shell2 >> shell3 >> shell4 >> TBME >> hom >> r2 >> p2;
    //TBME *= Parameters.tbstrength;
    //std::cout << coupT << " " << par << " " << coupJ << " " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if(shell1 == shell2){ TBME *= std::sqrt(2.0); } // !! check
    if(shell3 == shell4){ TBME *= std::sqrt(2.0); } // !! check
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    tb.j = coupJ;
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
    key1 = Chan.tb_map[ind1][Hash2(shell1, shell2, Space.indtot)];
    key2 = Chan.tb_map[ind1][Hash2(shell3, shell4, Space.indtot)];
    key3 = Chan.tb_map[ind1][Hash2(shell2, shell1, Space.indtot)];
    key4 = Chan.tb_map[ind1][Hash2(shell4, shell3, Space.indtot)];
    ind = key1 * Chan.ntb[ind1] + key2;
    ME.V[ind1][ind] = TBME;
    ind = key3 * Chan.ntb[ind1] + key2;
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
    ind = key1 * Chan.ntb[ind1] + key4;
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
    ind = key3 * Chan.ntb[ind1] + key4;
    ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = key2 * Chan.ntb[ind1] + key1;
      ME.V[ind1][ind] = TBME;
      ind = key4 * Chan.ntb[ind1] + key1;
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ) + 1)) * TBME;
      ind = key2 * Chan.ntb[ind1] + key3;
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ) + 1)) * TBME;
      ind = key4 * Chan.ntb[ind1] + key3;
      ME.V[ind1][ind] = pow(-1.0, int(0.5*(Space.qnums[shell1].j + Space.qnums[shell2].j + Space.qnums[shell3].j + Space.qnums[shell4].j))) * TBME;
    }
  }
  interaction.close();
}

void Read_Matrix_Elements_M(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind1, ind, key1, key2, key3, key4;
  State tb;
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;
  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
    TBME *= Parameters.tbstrength;
    //std::cout << "? " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << "  " << TBME << std::endl;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if(shell1 == shell2 || shell3 == shell4){ continue; }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
    key1 = Chan.tb_map[ind1][Hash2(shell1, shell2, Space.indtot)];
    key2 = Chan.tb_map[ind1][Hash2(shell3, shell4, Space.indtot)];
    key3 = Chan.tb_map[ind1][Hash2(shell2, shell1, Space.indtot)];
    key4 = Chan.tb_map[ind1][Hash2(shell4, shell3, Space.indtot)];
    ind = key1 * Chan.ntb[ind1] + key2;
    ME.V[ind1][ind] = TBME;
    ind = key3 * Chan.ntb[ind1] + key2;
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = key1 * Chan.ntb[ind1] + key4;
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = key3 * Chan.ntb[ind1] + key4;
    ME.V[ind1][ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = key2 * Chan.ntb[ind1] + key1;
      ME.V[ind1][ind] = TBME;
      ind = key4 * Chan.ntb[ind1] + key1;
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = key2 * Chan.ntb[ind1] + key3;
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = key4 * Chan.ntb[ind1] + key3;
      ME.V[ind1][ind] = TBME;
    }
  }
  interaction.close();
}

void Read_Matrix_Elements_QD(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::cout << "Computing Coulomb Matrix Elements for QD..." << std::endl;
  struct timespec time1, time2;
  double elapsed0 = 0.0;
  ME = HF_Matrix_Elements(Chan);
  clock_gettime(CLOCK_MONOTONIC, &time1);
  
  int pq, rs, p, q, r, s;
  int ntb, ntb0;
  int ind1, ind2, ind3, ind4;
  int length, length0;
  int *tbvec0;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    tbvec0 = new int[ntb];
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tbvec[chan][2*i] < Chan.tbvec[chan][2*i + 1]){
	tbvec0[2*ntb0] = Chan.tbvec[chan][2*i];
	tbvec0[2*ntb0 + 1] = Chan.tbvec[chan][2*i + 1];
	++ntb0;
      }
    }
    
    length = int(0.5 * ntb0 * (ntb0 + 1));
    #pragma omp parallel private(pq, rs, p, q, r, s, length0, ind1, ind2, ind3, ind4, TBME)
    {
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < length; ++pqrs){
	pq = std::floor((2*ntb0 - 1 - std::sqrt(1 + 4*ntb0 + 4*ntb0*ntb0 - 8*pqrs))/2) + 1;
	length0 = int(0.5 * pq * (2*ntb0 - pq + 1));
	rs = int(pq + pqrs - length0);
	p = tbvec0[2*pq];
	q = tbvec0[2*pq + 1];
	r = tbvec0[2*rs];
	s = tbvec0[2*rs + 1];
	ind1 = Chan.tb_map[chan][Hash2(p, q, Space.indtot)];
	ind2 = Chan.tb_map[chan][Hash2(r, s, Space.indtot)];
	ind3 = Chan.tb_map[chan][Hash2(q, p, Space.indtot)];
	ind4 = Chan.tb_map[chan][Hash2(s, r, Space.indtot)];
	TBME = Coulomb_HO(Parameters, Space, p, q, r, s);
	if(fabs(TBME) < 1.0e-14){ TBME = 0.0; }
	ME.V[chan][ind1 * ntb + ind2] = TBME;
	ME.V[chan][ind3 * ntb + ind2] = -1.0 * TBME;
	ME.V[chan][ind1 * ntb + ind4] = -1.0 * TBME;
	ME.V[chan][ind3 * ntb + ind4] = TBME;
	if(ind1 != ind2){
	  ME.V[chan][ind2 * ntb + ind1] = TBME;
	  ME.V[chan][ind4 * ntb + ind1] = -1.0 * TBME;
	  ME.V[chan][ind2 * ntb + ind3] = -1.0 * TBME;
	  ME.V[chan][ind4 * ntb + ind3] = TBME;
	}
      }
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "!! Runtime = " << elapsed0 << " sec. " << std::endl;
}

void Read_QD_ME_From_File(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  std::ifstream interaction;	// interaction file
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + "coulomb-ho2d-elements-20-shells.dat";

  interaction.open(fullpath1.c_str(), std::ios::binary);
  if(interaction.is_open()){
    //get length of file
    interaction.seekg(0, interaction.end);
    int length = interaction.tellg();
    interaction.seekg(0, interaction.beg);

    char *buffer = new char[length];
    interaction.read(buffer, length);
    interaction.close();

    double TBME;
    unsigned int neg1 = 128;
    unsigned int neg2 = 4294967040;
    unsigned int n1, n2, n3, n4;
    int ml1, ml2, ml3, ml4;
    int nmax = Parameters.Shells;
    int key, key1, key2;
    State tb;
    int ind, ind1;
    State statep, stateq, stater, states;
    int p, q, r, s;
    statep.t = -1, stateq.t = -1, stater.t = -1, states.t = -1;
    length /= 16;
    #pragma omp parallel private(n1, n2, n3, n4, ml1, ml2, ml3, ml4, key, p, q, r, s, statep, stateq, stater, states, TBME, ind, ind1, tb, key1, key2)
    {
      #pragma omp for schedule(static)
      for(int i = 0; i < length; ++i){
	n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	ml1 = 0, ml2 = 0, ml3 = 0, ml4 = 0;
	
	n1 = *(buffer + 16*i);
	ml1 = *(buffer + 16*i + 1);
	if((neg1 & ml1) != 0){ ml1 = (ml1 | neg2); }// ml1 < 0
	n2 = *(buffer + 16*i + 2);
	ml2 = *(buffer + 16*i + 3);
	if((neg1 & ml2) != 0){ ml2 = (ml2 | neg2); }// ml2 < 0
	n3 = *(buffer + 16*i + 4);
	ml3 = *(buffer + 16*i + 5);
	if((neg1 & ml3) != 0){ ml3 = (ml3 | neg2); }// ml3 < 0
	n4 = *(buffer + 16*i + 6);
	ml4 = *(buffer + 16*i + 7);
	if((neg1 & ml4) != 0){ ml4 = (ml4 | neg2); }// ml4 < 0
	TBME = *(double*)(buffer + 16*i + 8);
	if(int(2*n1 + abs(ml1)) >= nmax || int(2*n2 + abs(ml2)) >= nmax || int(2*n3 + abs(ml3)) >= nmax || int(2*n4 + abs(ml4)) >= nmax){ continue; }
	//std::cout << "< " << n1 << "," << ml1 << " ; " << n2 << "," << ml2 << " |V| " << n3 << "," << ml3 << " ; " << n4 << "," << ml4 << " > = " << TBME << std::endl;
	TBME *= std::sqrt(Parameters.density);
	statep.n = n1;
	statep.ml = ml1;
	stateq.n = n2;
	stateq.ml = ml2;
	stater.n = n3;
	stater.ml = ml3;
	states.n = n4;
	states.ml = ml4;
	for(int s1 = -1; s1 <= 1; s1 += 2){
	  statep.m = s1;
	  key = ChanInd_1b(Parameters.basis, Space, statep);
	  p = Space.map_1b[key];
	  for(int s2 = -1; s2 <= 1; s2 += 2){
	    stateq.m = s2;
	    key = ChanInd_1b(Parameters.basis, Space, stateq);
	    q = Space.map_1b[key];
	    if(p == q){ continue; }
	    for(int s3 = -1; s3 <= 1; s3 += 2){
	      stater.m = s3;
	      key = ChanInd_1b(Parameters.basis, Space, stater);
	      r = Space.map_1b[key];
	      if(s3 != s1){ continue; }
	      for(int s4 = -1; s4 <= 1; s4 += 2){
		states.m = s4;
		key = ChanInd_1b(Parameters.basis, Space, states);
		s = Space.map_1b[key];
		if(r == s || s4 != s2){ continue; }
		
		plus(tb, Space.qnums[p], Space.qnums[q]);
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);

		// C(p1q2r3s4) -> <p1q2 || r3s4>
		key1 = Chan.tb_map[ind1][Hash2(p, q, Space.indtot)];
		key2 = Chan.tb_map[ind1][Hash2(r, s, Space.indtot)];
		ind = key1 * Chan.ntb[ind1] + key2;
		ME.V[ind1][ind] += TBME;
		// C(p1q2r3s4) -> -<p1q2 || s4r3>
		key2 = Chan.tb_map[ind1][Hash2(s, r, Space.indtot)];
		ind = key1 * Chan.ntb[ind1] + key2;
		ME.V[ind1][ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(q2p1s4r3) -> <q2p1 || s4r3>
		  key1 = Chan.tb_map[ind1][Hash2(q, p, Space.indtot)];
		  key2 = Chan.tb_map[ind1][Hash2(s, r, Space.indtot)];
		  ind = key1 * Chan.ntb[ind1] + key2;
		  ME.V[ind1][ind] += TBME;
		  // C(p1q2r3s4) = C(q2p1s4r3) -> -<q2p1 || r3s4>
		  key2 = Chan.tb_map[ind1][Hash2(r, s, Space.indtot)];
		  ind = key1 * Chan.ntb[ind1] + key2;
		  ME.V[ind1][ind] -= TBME;
		}
		if(((n1 == n3 && ml1 == ml3) && (n2 == n4 && ml2 == ml4)) || ((n1 == n4 && ml1 == ml4) && (n2 == n3 && ml2 == ml3))){ continue; }
		// C(p1q2r3s4) = C(r3s4p1q2) -> <r3s4 || p1q2>
		key1 = Chan.tb_map[ind1][Hash2(r, s, Space.indtot)];
		key2 = Chan.tb_map[ind1][Hash2(p, q, Space.indtot)];
		ind = key1 * Chan.ntb[ind1] + key2;
		ME.V[ind1][ind] += TBME;
		// C(p1q2r3s4) = C(r3s4p1q2) -> -<r3s4 || q2p1>
		key2 = Chan.tb_map[ind1][Hash2(q, p, Space.indtot)];
		ind = key1 * Chan.ntb[ind1] + key2;
		ME.V[ind1][ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(s4r3q2p1) -> <s4r3 || q2p1>
		  key1 = Chan.tb_map[ind1][Hash2(s, r, Space.indtot)];
		  key2 = Chan.tb_map[ind1][Hash2(q, p, Space.indtot)];
		  ind = key1 * Chan.ntb[ind1] + key2;
		  ME.V[ind1][ind] += TBME;
		  // C(p1q2r3s4) = C(s4r3q2p1) -> -<s4r3 || p1q2>
		  key2 = Chan.tb_map[ind1][Hash2(p, q, Space.indtot)];
		  ind = key1 * Chan.ntb[ind1] + key2;
		  ME.V[ind1][ind] -= TBME;
		}
	      }
	    }
	  }
	}
      }
    }
    delete[] buffer;
  }
}

void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, const HF_Matrix_Elements &ME)
{
  //std::cout << "Diagonalizing Hartree-Fock Matrix ..." << std::endl;
  std::cout << std::setprecision(12);
  double error1; // Energy error between iterations
  double error2; // Energy error between iterations
  double Bshift; // level shift parameter
  //double width;
  int ind; // Index to keep track of iteration number
  int nob1;
  int **tempn;
  int badcount;
  Single_Particle_States HF0, HF2;

  //// for DIIS //////
  /*double norm, maxnorm;
  int maxind, ortho;
  double checkdot;
  double checknorm1, checknorm2;
  int DIISstart = 50;
  int maxl = 10;
  int N = 0;
  int P = N + 1;
  double ***p = new double**[maxl];
  double ***delp = new double**[maxl];
  for(int l = 0; l < maxl; ++l){
    p[l] = new double*[Chan.size3];
    delp[l] = new double*[Chan.size3];
    for(int chan = 0; chan < Chan.size3; ++chan){
      p[l][chan] = new double[Chan.nob[chan] * Chan.nob[chan]];
      delp[l][chan] = new double[Chan.nob[chan] * Chan.nob[chan]];
      for(int ob = 0; ob < Chan.nob[chan] * Chan.nob[chan]; ++ob){
	p[l][chan][ob] = 0.0;
	delp[l][chan][ob] = 0.0;
      }
    }
  }
  double **tempdelp = new double*[Chan.size3];
  for(int chan = 0; chan < Chan.size3; ++chan){
    tempdelp[chan] = new double[Chan.nob[chan] * Chan.nob[chan]];
    for(int ob = 0; ob < Chan.nob[chan] * Chan.nob[chan]; ++ob){
      tempdelp[chan][ob] = 0.0;
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
  
  HF0 = Single_Particle_States(Parameters, Space, Chan);
  HF2 = Single_Particle_States(Parameters, Space, Chan);
  tempn = new int*[Chan.size3];
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    HF0.h[i] = HF.h[i];
    HF2.h[i] = HF.h[i];
    tempn[i] = new int[nob1];
    for(int j = 0; j < nob1; ++j){
      tempn[i][j] = Space.qnums[Chan.obvec[i][j]].n;
      HF0.energies[i][j] = HF.energies[i][j];
      HF2.energies[i][j] = HF.energies[i][j];
      for(int k = 0; k < nob1; ++k){
	HF0.vectors[i][j][k] = HF.vectors[i][j][k];
	HF2.vectors[i][j][k] = HF.vectors[i][j][k];
      }
    }
  }

  ind = 0;
  error1 = 1000.0;
  error2 = 2000.0;
  badcount = 0;
  Bshift = 100.0;
  std::cout << std::setprecision(6);
  while((error2 > 1e-12 && ind < 20000) || ind < 20){
    ++ind;
    Hartree_Fock_Step(Parameters, Space, Chan, HF, HF2, ME, Bshift, error1);
    //// for DIIS //////
    /*norm = 0.0;
    for(int i = 0; i < Chan.size3; ++i){
      nob1 = Chan.nob[i];
      for(int j = 0; j < nob1; ++j){
	for(int k = 0; k < nob1; ++k){
	  tempdelp[i][nob1 * j + k] = HF2.vectors[i][j][k] - HF.vectors[i][j][k];
	  norm += tempdelp[i][nob1 * j + k] * tempdelp[i][nob1 * j + k];
	}
      }
    }

    // check orthogonality of tempdelp
    ortho = 1;
    for(int l = 0; l < N; ++l){
      checkdot = 0.0;
      checknorm1 = 0.0;
      checknorm2 = 0.0;
      for(int i = 0; i < Chan.size3; ++i){
	nob1 = Chan.nob[i];
	for(int j = 0; j < nob1; ++j){
	  for(int k = 0; k < nob1; ++k){
	    checkdot += (tempdelp[i][nob1 * j + k] * delp[l][i][nob1 * j + k]);
	    checknorm1 += (tempdelp[i][nob1 * j + k] * tempdelp[i][nob1 * j + k]);
	    checknorm2 += (delp[l][i][nob1 * j + k] * delp[l][i][nob1 * j + k]);
	  }
	}
      }
      checkdot = fabs(checkdot/(std::sqrt(checknorm1) * std::sqrt(checknorm2)));
      if(checkdot > 0.9 || norm > B[P * l + l]){ ortho = 0; break; }
    }
    if(ind < DIISstart){ ortho = 0; }

    if(ortho == 1){
      if(N < maxl){
	for(int i = 0; i < Chan.size3; ++i){
	  nob1 = Chan.nob[i];
	  for(int j = 0; j < nob1; ++j){
	    for(int k = 0; k < nob1; ++k){
	      p[N][i][nob1 * j + k] = HF2.vectors[i][j][k];
	      delp[N][i][nob1 * j + k] = tempdelp[i][nob1 * j + k];
	    }
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
	  for(int i = 0; i < Chan.size3; ++i){
	    nob1 = Chan.nob[i];
	    for(int j = 0; j < nob1; ++j){
	      for(int k = 0; k < nob1; ++k){
		B2[(P+1) * N + l] += delp[N][i][nob1 * j + k] * delp[l][i][nob1 * j + k];
		if(l != N){ B2[(P+1) * l + N] += delp[l][i][nob1 * j + k] * delp[N][i][nob1 * j + k]; }
	      }
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
	for(int i = 0; i < Chan.size3; ++i){
	  nob1 = Chan.nob[i];
	  for(int j = 0; j < nob1; ++j){
	    for(int k = 0; k < nob1; ++k){
	      p[maxind][i][nob1 * j + k] = HF2.vectors[i][j][k];
	      delp[maxind][i][nob1 * j + k] = tempdelp[i][nob1 * j + k];
	    }
	  }
	}
	for(int l = 0; l < N; ++l){
	  B2[P * maxind + l] = 0.0;
	  if(l != maxind){ B2[P * l + maxind] = 0.0; }
	  for(int i = 0; i < Chan.size3; ++i){
	    nob1 = Chan.nob[i];
	    for(int j = 0; j < nob1; ++j){
	      for(int k = 0; k < nob1; ++k){
		B2[P * maxind + l] += delp[maxind][i][nob1 * j + k] * delp[l][i][nob1 * j + k];
		if(l != maxind){ B2[P * l + maxind] += delp[l][i][nob1 * j + k] * delp[maxind][i][nob1 * j + k]; }
	      }
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
      for(int i = 0; i < Chan.size3; ++i){
	nob1 = Chan.nob[i];
	for(int j = 0; j < nob1; ++j){
	  for(int k = 0; k < nob1; ++k){
	    HF2.vectors[i][j][k] = 0.0;
	  }
	}
      }
      for(int l = 0; l < N; ++l){
	for(int i = 0; i < Chan.size3; ++i){
	  nob1 = Chan.nob[i];
	  for(int j = 0; j < nob1; ++j){
	    for(int k = 0; k < nob1; ++k){
	      HF2.vectors[i][j][k] += -1.0 * B2[P * l + N] * p[l][i][nob1 * j + k];
	    }
	  }
	}
      }
      for(int i = 0; i < Chan.size3; ++i){
	nob1 = Chan.nob[i];
	GramSchmidt(HF2.vectors[i], nob1);
      }
      }*/
    ////////////////////

    if(error1 > error2){
      if(error2 < 1.0e-7 && ind > 50){
	++badcount;
	if(badcount == 5){
	  for(int i = 0; i < Chan.size3; ++i){
	    nob1 = Chan.nob[i];
	    for(int j = 0; j < nob1; ++j){
	      HF.energies[i][j] = HF0.energies[i][j];
	      for(int k = 0; k < nob1; ++k){
		HF.vectors[i][j][k] = HF0.vectors[i][j][k];
	      }
	    }
	  }
	  break;
	}
      }
      for(int i = 0; i < Chan.size3; ++i){
	nob1 = Chan.nob[i];
	for(int j = 0; j < nob1; ++j){
	  HF.energies[i][j] = HF2.energies[i][j];
	  for(int k = 0; k < nob1; ++k){
	    HF.vectors[i][j][k] = HF2.vectors[i][j][k];
	  }
	}
      }
    }
    else{
      --badcount;
      for(int i = 0; i < Chan.size3; ++i){
	nob1 = Chan.nob[i];
	for(int j = 0; j < nob1; ++j){
	  HF0.energies[i][j] = HF.energies[i][j];
	  HF.energies[i][j] = HF2.energies[i][j];
	  for(int k = 0; k < nob1; ++k){
	    HF0.vectors[i][j][k] = HF.vectors[i][j][k];
	    HF.vectors[i][j][k] = HF2.vectors[i][j][k];
	  }
	}
      }
    }
    error2 = error1;
  }
  if( error2 > 1e-12 ){ std::cout << "ind = " << ind << ", error = " << error2 << ". HF Solution Not Converged!!" << std::endl; }

  ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < HF.h[i]; ++j){
      Space.qnums[Chan.obvec[i][j]] = Chan.qnums3[i];
      Space.qnums[Chan.obvec[i][j]].n = tempn[i][j];
      Space.qnums[Chan.obvec[i][j]].energy = HF.energies[i][j];
      Space.qnums[Chan.obvec[i][j]].type = "hole";
      ++ind;
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = HF.h[i]; j < HF.h[i] + HF.p[i]; ++j){
      Space.qnums[Chan.obvec[i][j]] = Chan.qnums3[i];
      Space.qnums[Chan.obvec[i][j]].n = tempn[i][j];
      Space.qnums[Chan.obvec[i][j]].energy = HF.energies[i][j];
      Space.qnums[Chan.obvec[i][j]].type = "particle";
      ++ind;
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    delete[] tempn[i];
  }
  delete[] tempn;
  HF0.delete_struct(Chan);
  HF2.delete_struct(Chan);

  /*for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].n << " " << Space.qnums[i].m << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }*/

  /*for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].n << " " << Space.qnums[i].par << " " << Space.qnums[i].j << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }*/
}

void Hartree_Fock_Step(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, const HF_Matrix_Elements &ME, const double &Bshift, double &error)
{
  double temperror;
  double tempnorm1, tempnorm2;
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  State tb;
  double *w;
  double *work;
  double *fock;
  int size1, size2;
  int ind1, ind2;
  int key1, key2;
  int mind, nind, lind, kind;
  int minj;
  int m, l;
  int length, length0;
  jobz = 'V';
  uplo = 'U';

  temperror = 0.0;
  tempnorm1 = 0.0;
  tempnorm2 = 0.0;
  //Make Fock Matrix
  for(int i = 0; i < Chan.size3; ++i){
    size1 = Chan.nob[i];
    fock = new double[size1 * size1];
    length = int(0.5 * size1 * (size1 + 1));
    #pragma omp parallel private(m, l, mind, nind, lind, kind, tb, ind1, ind2, term, length0, key1, key2, size2, minj)
    {
      #pragma omp for schedule(static)
      for(int ml = 0; ml < length; ++ml){
	m = std::floor((2*size1 - 1 - std::sqrt(1 + 4*size1 + 4*size1*size1 - 8*ml))/2) + 1;
	length0 = int(0.5 * m * (2*size1 - m + 1));
	l = int(m + ml - length0);
	mind = Chan.obvec[i][m];
	lind = Chan.obvec[i][l];
	if(m == l){ fock[size1 * m + m] = Space.qnums[mind].energy; } // Add diagonal elements to fock matrices
	else{ fock[size1 * m + l] = 0.0; }
	for(int j = 0; j < Chan.size3; ++j){
	  size2 = Chan.nob[j];
	  for(int beta = 0; beta < HF.h[j]; ++beta){ // Sum over occupied levels
	    for(int n = 0; n < size2; ++n){
	      nind = Chan.obvec[j][n];
	      for(int k = 0; k < size2; ++k){
		kind = Chan.obvec[j][k];
		if(Parameters.basis != "finite_J" && Parameters.basis != "finite_JM"){
		  plus(tb, Space.qnums[mind], Space.qnums[nind]);
		  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  key1 = Chan.tb_map[ind1][Hash2(mind, nind, Space.indtot)];
		  key2 = Chan.tb_map[ind1][Hash2(lind, kind, Space.indtot)];
		  ind2 = key1 * Chan.ntb[ind1] + key2;
		  term = HF.vectors[j][beta][n] * HF.vectors[j][beta][k] * ME.V[ind1][ind2];
		  fock[size1 * m + l] += term;
		}
		else{
		  plus(tb, Space.qnums[mind], Space.qnums[nind]);
		  minj = abs(Chan.qnums3[i].j - Chan.qnums3[j].j);
		  while(tb.j >= minj){
		    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    key1 = Chan.tb_map[ind1][Hash2(mind, nind, Space.indtot)];
		    key2 = Chan.tb_map[ind1][Hash2(lind, kind, Space.indtot)];
		    ind2 = key1 * Chan.ntb[ind1] + key2;
		    term = HF.vectors[j][beta][n] * HF.vectors[j][beta][k] * ME.V[ind1][ind2];
		    term *= (tb.j + 1.0)/(Chan.qnums3[i].j + 1.0);
		    fock[size1 * m + l] += term;
		    tb.j -= 2;
		  }
		}
	      }
	    }
	  }
	}
	for(int beta = HF.h[i]; beta < Chan.nob[i]; ++beta){ // Sum over occupied levels
	  fock[size1 * m + l] += HF.vectors[i][beta][m] * HF.vectors[i][beta][l] * Bshift;
	}
	if(m != l){ fock[size1 * l + m] = fock[size1 * m + l]; }
      }
    }

    lda = size1;
    lwork = (3+2)*size1;
    w = new double[lda];
    work = new double[lwork];
    for(int j = 0; j < size1; ++j){ w[j] = 0.0; }
    for(int j = 0; j < (3+2)*size1; ++j){ work[j] = 0.0; }
    if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, fock, &lda, w, work, &lwork, &info); }
    for(int j = HF.h[i]; j < Chan.nob[i]; ++j){ w[j] -= Bshift; } // Add back Level-shift parameter
    GramSchmidt(fock, size1);

    for(int j = 0; j < size1; ++j){
      HF2.energies[i][j] = w[j];
      if(fock[size1 * j + j] < 0.0){
	for(int k = 0; k < size1; ++k){ fock[size1 * j + k] *= -1.0; }
      }
      for(int k = 0; k < size1; ++k){
	temperror += (fabs(HF.vectors[i][j][k]) - fabs(fock[size1 * j + k])) * (fabs(HF.vectors[i][j][k]) - fabs(fock[size1 * j + k]));
	tempnorm1 += HF.vectors[i][j][k] * HF.vectors[i][j][k];
	tempnorm2 += fock[size1 * j + k] * fock[size1 * j + k];
	HF2.vectors[i][j][k] = fock[size1 * j + k];
      }
    }
    delete[] fock;
    delete[] w;
    delete[] work;
  }
  error = std::sqrt(temperror/(tempnorm1*tempnorm2));
}


void Randomize_HF(const HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, double &width)
{
  int nob1;
  double tempc;
  double rand;
  double norm;
  size_t key;
  std::unordered_map<size_t,double> c_map;
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    for(int j = 0; j < nob1; ++j){
      for(int k = 0; k < nob1; ++k){
	tempc = HF.vectors[i][j][k];
	if(fabs(tempc) > 1.0e-16){
	  rand = rand_normal(0.0, width * fabs(tempc));
	  key = std::hash<float>{}(float(fabs(tempc)));
	  c_map[key] = rand;
	}
      }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    for(int j = 0; j < nob1; ++j){
      for(int k = 0; k < nob1; ++k){
	tempc = HF.vectors[i][j][k];
	key = std::hash<float>{}(float(fabs(tempc)));
	if(tempc > 1.0e-16){ tempc += c_map[key]; }
	else if(tempc < -1.0e-16){ tempc -= c_map[key]; }
	HF2.vectors[i][j][k] = tempc;
      }
    }
    for(int j = 0; j < nob1; ++j){
      HF2.energies[i][j] = HF.energies[i][j];
      norm = 0.0;
      if(HF2.vectors[i][j][j] < 0.0){
	for(int k = 0; k < nob1; ++k){ HF2.vectors[i][j][k] *= -1.0; }
      }
    }
  }
}

void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  int length, matlength; // max length of M-Scheme indicies, length of J_ME
  double tempel;
  double *M1; // Matrices of coefficients
  double *C;
  char transa, transb;
  double alpha1, beta1;
  std::ofstream jschemefile; // file to print M-Scheme matrix elements
  int p, q, a, g;
  int ind1, ind2;
  State tb;

  for(int chan = 0; chan < Chan.size1; ++chan){
    length = Chan.ntb[chan];
    if(length == 0){ continue; }
    matlength = pow(length, 2.0);
    M1 = new double[matlength];
    C = new double[matlength];
    for(int i = 0; i < matlength; ++i){
      M1[i] = 0.0;
      C[i] = 0.0;
    }
    #pragma omp parallel private(p, q, a, g, ind1, ind2, tempel)
    {
      #pragma omp for schedule(static)
      for(int pq = 0; pq < length; ++pq){
	for(int ag = 0; ag < length; ++ag){
	  p = Chan.tbvec[chan][2*pq];
	  q = Chan.tbvec[chan][2*pq + 1];
	  a = Chan.tbvec[chan][2*ag];
	  g = Chan.tbvec[chan][2*ag + 1];
	  if(Chan.indvec[p] != Chan.indvec[a] || Chan.indvec[q] != Chan.indvec[g]){ continue; }
	  ind1 = Chan.indvec[p];
	  ind2 = Chan.indvec[q];
	  p = Chan.ob_map[ind1][p];
	  q = Chan.ob_map[ind2][q];
	  a = Chan.ob_map[ind1][a];
	  g = Chan.ob_map[ind2][g];
	  tempel = States.vectors[ind1][a][p] * States.vectors[ind2][g][q];
	  M1[length * pq + ag] = tempel;
	}
      }
    }

    transa = 'N';
    transb = 'T';
    alpha1 = 1.0;
    beta1 = 0.0;
    dgemm_NN(ME.V[chan], M1, C, &length, &length, &length, &alpha1, &beta1, &transa, &transa);
    dgemm_TN(M1, C, ME.V[chan], &length, &length, &length, &alpha1, &beta1, &transb, &transa);

    double tempV;
    for(int i = 0; i < length; ++i){
      for(int j = i; j < length; ++j){
	tempV = 0.5 * (ME.V[chan][length*i + j] + ME.V[chan][length*j + i]);
	if(fabs(tempV) < 1.0e-14){
	  ME.V[chan][length*i + j] = 0.0;
	  ME.V[chan][length*j + i] = 0.0;
	}
	else{
	  ME.V[chan][length*i + j] = tempV;
	  ME.V[chan][length*j + i] = tempV;
	}
      }
    }      
    delete[] M1;
    delete[] C;
  }
}

void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, key1, key2, key3;
  std::string ptype, qtype, rtype, stype;
  State tb;

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
      shell1 = HF_Chan.tbvec[chan][2*tb1];
      shell2 = HF_Chan.tbvec[chan][2*tb1 + 1];
      ptype = Space.qnums[shell1].type;
      qtype = Space.qnums[shell2].type;
      if(shell1 == shell2){ continue; }
      for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	shell3 = HF_Chan.tbvec[chan][2*tb2];
	shell4 = HF_Chan.tbvec[chan][2*tb2 + 1];
	rtype = Space.qnums[shell3].type;
	stype = Space.qnums[shell4].type;
	if(shell3 == shell4){ continue; }
	if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "hole"){ continue; }
	if(ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "particle"){ continue; }
	TBME = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];
	//std::cout << "V_hf: " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << " = " << TBME << std::endl;

	if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  key1 = Chan.pp_map[ind1][Hash2(shell3, shell4, Space.indtot)];
	  key2 = Chan.pp_map[ind1][Hash2(shell1, shell2, Space.indtot)];
	  ind = key1 * Chan.npp[ind1] + key2;
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  key1 = Chan.hh_map[ind1][Hash2(shell1, shell2, Space.indtot)];
	  key2 = Chan.hh_map[ind1][Hash2(shell3, shell4, Space.indtot)];
	  ind = key1 * Chan.nhh[ind1] + key2;
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind1][Hash2(shell1, shell4, Space.indtot)];
	  key2 = Chan.hp2_map[ind1][Hash2(shell3, shell2, Space.indtot)];
	  ind = key1 * Chan.nhp2[ind1] + key2;
	  Ints.D_ME1.V3[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  key1 = Chan.pp_map[ind1][Hash2(shell3, shell4, Space.indtot)];
	  key2 = Chan.hh_map[ind1][Hash2(shell1, shell2, Space.indtot)];
	  ind = key1 * Chan.nhh[ind1] + key2;
	  Ints.D_ME1.V4[ind1][ind] = TBME;
	    
	  ind1 = Chan.indvec[shell2];
	  key1 = Chan.h_map[ind1][shell2];
	  key2 = Chan.hpp_map[ind1][Hash3(shell1, shell3, shell4, Space.indtot)];
	  ind = key1 * Chan.nhpp[ind1] + key2;
	  Ints.D_ME1.V5[ind1][ind] = TBME;
	    
	  ind1 = Chan.indvec[shell1];
	  key1 = Chan.h_map[ind1][shell1];
	  key2 = Chan.hpp_map[ind1][Hash3(shell2, shell3, shell4, Space.indtot)];
	  ind = key1 * Chan.nhpp[ind1] + key2;
	  Ints.D_ME1.V6[ind1][ind] = TBME;
	    
	  ind1 = Chan.indvec[shell4];
	  key1 = Chan.p_map[ind1][shell4];
	  key2 = Chan.hhp_map[ind1][Hash3(shell1, shell2, shell3, Space.indtot)];
	  ind = key1 * Chan.nhhp[ind1] + key2;
	  Ints.D_ME1.V7[ind1][ind] = TBME;
	    
	  ind1 = Chan.indvec[shell3];
	  key1 = Chan.p_map[ind1][shell3];
	  key2 = Chan.hhp_map[ind1][Hash3(shell1, shell2, shell4, Space.indtot)];
	  ind = key1 * Chan.nhhp[ind1] + key2;
	  Ints.D_ME1.V8[ind1][ind] = TBME;

	  minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind1][Hash2(shell1, shell3, Space.indtot)];
	  key2 = Chan.hp1_map[ind1][Hash2(shell2, shell4, Space.indtot)];
	  ind = key1 * Chan.nhp1[ind1] + key2;
	  Ints.D_ME1.V9[ind1][ind] = TBME;
	    
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind1][Hash2(shell1, shell4, Space.indtot)];
	  key2 = Chan.hp1_map[ind1][Hash2(shell2, shell3, Space.indtot)];
	  ind = key1 * Chan.nhp1[ind1] + key2;
	  Ints.D_ME1.V10[ind1][ind] = TBME;
	}
	if(Parameters.approx == "singles"){
	  if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    ind1 = Chan.indvec[shell2];
	    key1 = Chan.p_map[ind1][shell2];
	    key2 = Chan.hpp_map[ind1][Hash3(shell1, shell3, shell4, Space.indtot)];
	    ind = key1 * Chan.nhpp[ind1] + key2;
	    Ints.S_ME1.V11[ind1][ind] = TBME;
	      
	    ind1 = Chan.indvec[shell3];
	    key1 = Chan.p_map[ind1][shell3];
	    key2 = Chan.hpp1_map[ind1][Hash3(shell1, shell2, shell4, Space.indtot)];
	    ind = key2 * Chan.np[ind1] + key1;
	    Ints.S_ME1.V13[ind1][ind] = TBME;
	      
	    minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    key1 = Chan.pp1_map[ind1][Hash2(shell4, shell2, Space.indtot)];
	    key2 = Chan.hp1_map[ind1][Hash2(shell1, shell3, Space.indtot)];
	    ind = key1 * Chan.nhp1[ind1] + key2;
	    Ints.S_ME1.V16[ind1][ind] = TBME;
	      
	    ind1 = Chan.indvec[shell2];
	    key1 = Chan.p_map[ind1][shell2];
	    key2 = Chan.hpp_map[ind1][Hash3(shell1, shell3, shell4, Space.indtot)];
	    ind = key2 * Chan.np[ind1] + key1;
	    Ints.S_ME1.V17[ind1][ind] = TBME;
	      
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.pp_map[ind1][Hash2(shell3, shell4, Space.indtot)];
	    key3 = Chan.hp_map[ind1][Hash2(shell1, shell2, Space.indtot)];
	    ind = key1 * Chan.nhp[ind1] + key3;
	    Ints.S_ME1.V20[ind1][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	    ind1 = Chan.indvec[shell3];
	    key1 = Chan.h_map[ind1][shell3];
	    key2 = Chan.hhp_map[ind1][Hash3(shell1, shell2, shell4, Space.indtot)];
	    ind = key1 * Chan.nhhp[ind1] + key2;
	    Ints.S_ME1.V12[ind1][ind] = TBME;
	      
	    ind1 = Chan.indvec[shell1];
	    key1 = Chan.h_map[ind1][shell1];
	    key2 = Chan.hhp1_map[ind1][Hash3(shell2, shell3, shell4, Space.indtot)];
	    ind = key2 * Chan.nh[ind1] + key1;
	    Ints.S_ME1.V14[ind1][ind] = TBME;
	      
	    minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    key1 = Chan.hh1_map[ind1][Hash2(shell3, shell2, Space.indtot)];
	    key2 = Chan.hp1_map[ind1][Hash2(shell1, shell4, Space.indtot)];
	    ind = key1 * Chan.nhp1[ind1] + key2;
	    Ints.S_ME1.V15[ind1][ind] = TBME;
	      
	    ind1 = Chan.indvec[shell3];
	    key1 = Chan.h_map[ind1][shell3];
	    key2 = Chan.hhp_map[ind1][Hash3(shell1, shell2, shell4, Space.indtot)];
	    ind = key2 * Chan.nh[ind1] + key1;
	    Ints.S_ME1.V18[ind1][ind] = TBME;
	      
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(shell1, shell2, Space.indtot)];
	    key3 = Chan.hp_map[ind1][Hash2(shell3, shell4, Space.indtot)];
	    ind = key3 * Chan.nhh[ind1] + key1;
	    Ints.S_ME1.V19[ind1][ind] = TBME;
	  }
	}
      }
    }
  }
}

void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int p, q, r, s; // interaction file contents
  int ind, ind1, key1, key2, key3, jmin;
  double pj, qj, rj, sj, tbj, J, X;
  std::string ptype, qtype, rtype, stype;
  State tb;

  //std::cout << "V_hf: " << "p" << " " << "q" << " " << "r" << " " << "s" << ", " << "ind1" << " " << "key1" << " " << "key2" << ", " << "J" << " " << "tbj" << " = " << "X * TBME" << std::endl;
  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    J = 0.5 * HF_Chan.qnums1[chan].j;
    for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
	p = HF_Chan.tbvec[chan][2*tb1];
	q = HF_Chan.tbvec[chan][2*tb1 + 1];
	ptype = Space.qnums[p].type;
	qtype = Space.qnums[q].type;
	pj = 0.5 * Space.qnums[p].j;
	qj = 0.5 * Space.qnums[q].j;
	for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	  r = HF_Chan.tbvec[chan][2*tb2];
	  s = HF_Chan.tbvec[chan][2*tb2 + 1];
	  rtype = Space.qnums[r].type;
	  stype = Space.qnums[s].type;
	  rj = 0.5 * Space.qnums[r].j;
	  sj = 0.5 * Space.qnums[s].j;
	  TBME = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];
	  //std::cout << "V_hf: " << p << " " << q << " " << r << " " << s << ", " << 0.5*HF_Chan.qnums1[chan].j << " = " << TBME << std::endl;

	  if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    key1 = Chan.pp_map[chan][Hash2(r, s, Space.indtot)];
	    key2 = Chan.pp_map[chan][Hash2(p, q, Space.indtot)];
	    ind = key1 * Chan.npp[chan] + key2;
	    Ints.D_ME1.V1[chan][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	    key1 = Chan.hh_map[chan][Hash2(p, q, Space.indtot)];
	    key2 = Chan.hh_map[chan][Hash2(r, s, Space.indtot)];
	    ind = key1 * Chan.nhh[chan] + key2;
	    Ints.D_ME1.V2[chan][ind] = TBME;
	  }
	  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	    minus(tb, Space.qnums[s], Space.qnums[p]);
	    if(Space.qnums[q].j + Space.qnums[r].j < tb.j){ tb.j = Space.qnums[q].j + Space.qnums[r].j; }
	    jmin = abs(Space.qnums[s].j - Space.qnums[p].j);
	    if(abs(Space.qnums[q].j - Space.qnums[r].j) > jmin){ jmin = abs(Space.qnums[q].j - Space.qnums[r].j); }
	    while(tb.j >= jmin){
	      tbj = 0.5 * tb.j;
	      ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      key1 = Chan.hp2_map[ind1][Hash2(p, s, Space.indtot)];
	      key2 = Chan.hp2_map[ind1][Hash2(r, q, Space.indtot)];
	      ind = key1 * Chan.nhp2[ind1] + key2;
	      X = -1.0 * std::pow(-1.0, pj + qj + rj + sj) * (2.0 * J + 1) * CGC6(qj,pj,J,sj,rj,tbj);
	      Ints.D_ME1.V3[ind1][ind] += X * TBME;
	      tb.j -= 2;
	    }
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	    key1 = Chan.pp_map[chan][Hash2(r, s, Space.indtot)];
	    key2 = Chan.hh_map[chan][Hash2(p, q, Space.indtot)];
	    ind = key1 * Chan.nhh[chan] + key2;
	    Ints.D_ME1.V4[chan][ind] = TBME;

	    ind1 = Chan.indvec[q];
	    key1 = Chan.h_map[ind1][q];
	    key2 = Chan.hpp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, r, s, Space.indtot)];
	    ind = key1 * Chan.nhpp[ind1] + key2;
	    Ints.D_ME1.V5[ind1][ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;

	    ind1 = Chan.indvec[p];
	    key1 = Chan.h_map[ind1][p];
	    key2 = Chan.hpp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(q, r, s, Space.indtot)];
	    ind = key1 * Chan.nhpp[ind1] + key2;
	    Ints.D_ME1.V6[ind1][ind] = std::pow(-1.0, pj + qj - J) * std::sqrt((2.0*J + 1)/(2.0*pj + 1)) * TBME;

	    ind1 = Chan.indvec[s];
	    key1 = Chan.p_map[ind1][s];
	    key2 = Chan.hhp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, q, r, Space.indtot)];
	    ind = key1 * Chan.nhhp[ind1] + key2;
	    Ints.D_ME1.V7[ind1][ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*sj + 1)) * TBME;

	    ind1 = Chan.indvec[r];
	    key1 = Chan.p_map[ind1][r];
	    key2 = Chan.hhp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, q, s, Space.indtot)];
	    ind = key1 * Chan.nhhp[ind1] + key2;
	    Ints.D_ME1.V8[ind1][ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;

	    minus(tb, Space.qnums[r], Space.qnums[p]);
	    if(Space.qnums[q].j + Space.qnums[s].j < tb.j){ tb.j = Space.qnums[q].j + Space.qnums[s].j; }
	    jmin = abs(Space.qnums[r].j - Space.qnums[p].j);
	    if(abs(Space.qnums[q].j - Space.qnums[s].j) > jmin){ jmin = abs(Space.qnums[q].j - Space.qnums[s].j); }
	    while(tb.j >= jmin){
	      tbj = 0.5 * tb.j;
	      ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      key1 = Chan.hp2_map[ind1][Hash2(p, r, Space.indtot)];
	      key2 = Chan.hp1_map[ind1][Hash2(q, s, Space.indtot)];
	      ind = key1 * Chan.nhp1[ind1] + key2;
	      X = std::pow(-1.0, pj + qj - J) * (2.0 * J + 1) * CGC6(qj,pj,J,rj,sj,tbj);
	      Ints.D_ME1.V9[ind1][ind] += X * TBME;
	      tb.j -= 2;
	    }

	    minus(tb, Space.qnums[s], Space.qnums[p]);
	    if(Space.qnums[q].j + Space.qnums[r].j < tb.j){ tb.j = Space.qnums[q].j + Space.qnums[r].j; }
	    jmin = abs(Space.qnums[s].j - Space.qnums[p].j);
	    if(abs(Space.qnums[q].j - Space.qnums[r].j) > jmin){ jmin = abs(Space.qnums[q].j - Space.qnums[r].j); }
	    while(tb.j >= jmin){
	      tbj = 0.5 * tb.j;
	      ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	      key1 = Chan.hp2_map[ind1][Hash2(p, s, Space.indtot)];
	      key2 = Chan.hp1_map[ind1][Hash2(q, r, Space.indtot)];
	      ind = key1 * Chan.nhp1[ind1] + key2;
	      X = -1.0 * std::pow(-1.0, pj + qj + rj + sj) * (2.0 * J + 1) * CGC6(qj,pj,J,sj,rj,tbj);
	      Ints.D_ME1.V10[ind1][ind] += X * TBME;
	      tb.j -= 2;
	    }
	  }
	  if(Parameters.approx == "singles"){
	    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	      ind1 = Chan.indvec[q];
	      key1 = Chan.p_map[ind1][q];
	      key2 = Chan.hpp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, r, s, Space.indtot)];
	      ind = key1 * Chan.nhpp[ind1] + key2;
	      Ints.S_ME1.V11[ind1][ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	      ind = key2 * Chan.np[ind1] + key1;
	      Ints.S_ME1.V17[ind1][ind] = -1.0 * std::sqrt((2.0*J + 1)/(2.0*qj + 1)) * TBME;
	      
	      ind1 = Chan.indvec[r];
	      key1 = Chan.p_map[ind1][r];
	      key2 = Chan.hpp1_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, q, s, Space.indtot)];
	      ind = key2 * Chan.np[ind1] + key1;
	      Ints.S_ME1.V13[ind1][ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	      
	      minus(tb, Space.qnums[s], Space.qnums[q]);
	      if(Space.qnums[p].j + Space.qnums[r].j < tb.j){ tb.j = Space.qnums[p].j + Space.qnums[r].j; }
	      jmin = abs(Space.qnums[s].j - Space.qnums[q].j);
	      if(abs(Space.qnums[p].j - Space.qnums[r].j) > jmin){ jmin = abs(Space.qnums[p].j - Space.qnums[r].j); }
	      while(tb.j >= jmin){
		tbj = 0.5 * tb.j;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		key1 = Chan.pp1_map[ind1][Hash2(s, q, Space.indtot)];
		key2 = Chan.hp1_map[ind1][Hash2(p, r, Space.indtot)];
		ind = key1 * Chan.nhp1[ind1] + key2;
		X = std::pow(-1.0, rj + sj - J) * (2.0 * J + 1) * CGC6(pj,qj,J,sj,rj,tbj);
		Ints.S_ME1.V16[ind1][ind] += X * TBME;
		tb.j -= 2;
	      }
	      
	      key1 = Chan.pp_map[chan][Hash2(r, s, Space.indtot)];
	      key3 = Chan.hp_map[chan][Hash2(p, q, Space.indtot)];
	      ind = key1 * Chan.nhp[chan] + key3;
	      Ints.S_ME1.V20[chan][ind] = TBME;
	    }
	    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	      ind1 = Chan.indvec[r];
	      key1 = Chan.h_map[ind1][r];
	      key2 = Chan.hhp_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(p, q, s, Space.indtot)];
	      ind = key1 * Chan.nhhp[ind1] + key2;
	      Ints.S_ME1.V12[ind1][ind] = std::pow(-1.0, rj + sj - J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	      ind = key2 * Chan.nh[ind1] + key1;
	      Ints.S_ME1.V18[ind1][ind] = std::pow(-1.0, rj + sj + J) * std::sqrt((2.0*J + 1)/(2.0*rj + 1)) * TBME;
	      
	      ind1 = Chan.indvec[p];
	      key1 = Chan.h_map[ind1][p];
	      key2 = Chan.hhp1_map[ind1][int(J * std::pow(Space.indtot, 3)) + Hash3(q, r, s, Space.indtot)];
	      ind = key2 * Chan.nh[ind1] + key1;
	      Ints.S_ME1.V14[ind1][ind] = std::pow(-1.0, pj + qj + J) * std::sqrt((2.0*J + 1)/(2.0*pj + 1)) * TBME;
	      
	      minus(tb, Space.qnums[r], Space.qnums[q]);
	      if(Space.qnums[p].j + Space.qnums[s].j < tb.j){ tb.j = Space.qnums[p].j + Space.qnums[s].j; }
	      jmin = abs(Space.qnums[r].j - Space.qnums[q].j);
	      if(abs(Space.qnums[p].j - Space.qnums[s].j) > jmin){ jmin = abs(Space.qnums[p].j - Space.qnums[s].j); }
	      while(tb.j >= jmin){
		tbj = 0.5 * tb.j;
		ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		key1 = Chan.hh1_map[ind1][Hash2(r, q, Space.indtot)];
		key2 = Chan.hp1_map[ind1][Hash2(p, s, Space.indtot)];
		ind = key1 * Chan.nhp1[ind1] + key2;
		X = -1.0 * (2.0 * J + 1) * CGC6(pj,qj,J,rj,sj,tbj);
		Ints.S_ME1.V15[ind1][ind] += X * TBME;
		tb.j -= 2;
	      }
	      
	      key1 = Chan.hh_map[chan][Hash2(p, q, Space.indtot)];
	      key3 = Chan.hp_map[chan][Hash2(r, s, Space.indtot)];
	      ind = key3 * Chan.nhh[chan] + key1;
	      Ints.S_ME1.V19[chan][ind] = TBME;
	    }
	  }
	}
      }
      //}
  }
}

void Get_Matrix_Elements_JM(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME0, TBME, m1, m2, CGC1, CGC2; // interaction two-body interaction ME and two-body COM ME
  int t1, t2;
  int pind, qind, rind, sind;
  int p, q, r, s; // interaction file contents
  double pj, qj, rj, sj, coupj, pm, qm, rm, sm, coupm;
  int ind, ind1, key1, key2, key3;
  std::string ptype, qtype, rtype, stype;
  State tb;
  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    coupj = 0.5 * HF_Chan.qnums1[chan].j;
    for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
      p = HF_Chan.tbvec[chan][2*tb1];
      q = HF_Chan.tbvec[chan][2*tb1 + 1];
      pj = 0.5 * Space.qnums[Space.shellsm[p][0]].j;
      qj = 0.5 * Space.qnums[Space.shellsm[q][0]].j;
      ptype = Space.qnums[Space.shellsm[p][0]].type;
      qtype = Space.qnums[Space.shellsm[q][0]].type;
      for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	r = HF_Chan.tbvec[chan][2*tb2];
	s = HF_Chan.tbvec[chan][2*tb2 + 1];
	rj = 0.5 * Space.qnums[Space.shellsm[r][0]].j;
	sj = 0.5 * Space.qnums[Space.shellsm[s][0]].j;
	rtype = Space.qnums[Space.shellsm[r][0]].type;
	stype = Space.qnums[Space.shellsm[s][0]].type;
	TBME0 = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];

	for(int jz = -HF_Chan.qnums1[chan].j; jz <= HF_Chan.qnums1[chan].j; jz+=2){
	  coupm = 0.5 * jz;
	  for(int p1 = 0; p1 < -1*Space.qnums[Space.shellsm[p][0]].m + 1; ++p1){
	    pind = Space.shellsm[p][p1];
	    pm = 0.5 * Space.qnums[pind].m;
	    for(int q1 = 0; q1 < -1*Space.qnums[Space.shellsm[q][0]].m + 1; ++q1){
	      qind = Space.shellsm[q][q1];
	      qm = 0.5 * Space.qnums[qind].m;
	      m1 = pm + qm;
	      t1 = Space.qnums[pind].t + Space.qnums[qind].t;
	      if(pind == qind){ continue; }
	      for(int r1 = 0; r1 < -1*Space.qnums[Space.shellsm[r][0]].m + 1; ++r1){
		rind = Space.shellsm[r][r1];
		rm = 0.5 * Space.qnums[rind].m;
		for(int s1 = 0; s1 < -1*Space.qnums[Space.shellsm[s][0]].m + 1; ++s1){
		  sind = Space.shellsm[s][s1];
		  sm = 0.5 * Space.qnums[sind].m;
		  qm = 0.5 * Space.qnums[qind].m;
		  m2 = rm + sm;
		  t2 = Space.qnums[rind].t + Space.qnums[sind].t;
		  if(rind == sind){ continue; }
		  
		  if(t1 != t2 || m1 != coupm || m2 != coupm){ continue; }
		  CGC1 = CGC(pj, pm, qj, qm, coupj, coupm);
		  CGC2 = CGC(rj, rm, sj, sm, coupj, coupm);
		  TBME = TBME0 * CGC1 * CGC2;

		  if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[pind], Space.qnums[qind]);
		    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    key1 = Chan.pp_map[ind1][Hash2(pind, qind, Space.indtot)];
		    key2 = Chan.pp_map[ind1][Hash2(rind, sind, Space.indtot)];
		    ind = key1 * Chan.npp[ind1] + key2;
		    Ints.D_ME1.V1[ind1][ind] += TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
		    plus(tb, Space.qnums[pind], Space.qnums[qind]);
		    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    key1 = Chan.hh_map[ind1][Hash2(pind, qind, Space.indtot)];
		    key2 = Chan.hh_map[ind1][Hash2(rind, sind, Space.indtot)];
		    ind = key1 * Chan.nhh[ind1] + key2;
		    Ints.D_ME1.V2[ind1][ind] += TBME;
		  }
		  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
		    minus(tb, Space.qnums[sind], Space.qnums[pind]);
		    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		    key1 = Chan.hp2_map[ind1][Hash2(pind, sind, Space.indtot)];
		    key2 = Chan.hp2_map[ind1][Hash2(rind, qind, Space.indtot)];
		    ind = key1 * Chan.nhp2[ind1] + key2;
		    Ints.D_ME1.V3[ind1][ind] += TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[pind], Space.qnums[qind]);
		    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    key1 = Chan.hh_map[ind1][Hash2(pind, qind, Space.indtot)];
		    key2 = Chan.pp_map[ind1][Hash2(rind, sind, Space.indtot)];
		    ind = key2 * Chan.nhh[ind1] + key1;
		    Ints.D_ME1.V4[ind1][ind] += TBME;
		    
		    ind1 = Chan.indvec[qind];
		    key1 = Chan.h_map[ind1][qind];
		    key2 = Chan.hpp_map[ind1][Hash3(pind, rind, sind, Space.indtot)];
		    ind = key1 * Chan.nhpp[ind1] + key2;
		    Ints.D_ME1.V5[ind1][ind] += TBME;
		    
		    ind1 = Chan.indvec[pind];
		    key1 = Chan.h_map[ind1][pind];
		    key2 = Chan.hpp_map[ind1][Hash3(qind, rind, sind, Space.indtot)];
		    ind = key1 * Chan.nhpp[ind1] + key2;
		    Ints.D_ME1.V6[ind1][ind] += TBME;
		    
		    ind1 = Chan.indvec[sind];
		    key1 = Chan.p_map[ind1][sind];
		    key2 = Chan.hhp_map[ind1][Hash3(pind, qind, rind, Space.indtot)];
		    ind = key1 * Chan.nhhp[ind1] + key2;
		    Ints.D_ME1.V7[ind1][ind] += TBME;
		    
		    ind1 = Chan.indvec[rind];
		    key1 = Chan.p_map[ind1][rind];
		    key2 = Chan.hhp_map[ind1][Hash3(pind, qind, sind, Space.indtot)];
		    ind = key1 * Chan.nhhp[ind1] + key2;
		    Ints.D_ME1.V8[ind1][ind] += TBME;
		    
		    minus(tb, Space.qnums[rind], Space.qnums[pind]);
		    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		    key1 = Chan.hp2_map[ind1][Hash2(pind, rind, Space.indtot)];
		    key2 = Chan.hp1_map[ind1][Hash2(qind, sind, Space.indtot)];
		    ind = key1 * Chan.nhp1[ind1] + key2;
		    Ints.D_ME1.V9[ind1][ind] += TBME;
		    
		    minus(tb, Space.qnums[sind], Space.qnums[pind]);
		    ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		    key1 = Chan.hp2_map[ind1][Hash2(pind, sind, Space.indtot)];
		    key2 = Chan.hp1_map[ind1][Hash2(qind, rind, Space.indtot)];
		    ind = key1 * Chan.nhp1[ind1] + key2;
		    Ints.D_ME1.V10[ind1][ind] += TBME;
		  }
		  if(Parameters.approx == "singles"){
		    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		      ind1 = Chan.indvec[qind];
		      key1 = Chan.p_map[ind1][qind];
		      key2 = Chan.hpp_map[ind1][Hash3(pind, rind, sind, Space.indtot)];
		      ind = key1 * Chan.nhpp[ind1] + key2;
		      Ints.S_ME1.V11[ind1][ind] += TBME;
		      ind = key2 * Chan.np[ind1] + key1;
		      Ints.S_ME1.V17[ind1][ind] += TBME;
		      
		      ind1 = Chan.indvec[rind];
		      key1 = Chan.p_map[ind1][rind];
		      key2 = Chan.hpp1_map[ind1][Hash3(pind, qind, sind, Space.indtot)];
		      ind = key2 * Chan.np[ind1] + key1;
		      Ints.S_ME1.V13[ind1][ind] += TBME;
		      
		      minus(tb, Space.qnums[pind], Space.qnums[rind]);
		      ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		      key1 = Chan.pp1_map[ind1][Hash2(sind, qind, Space.indtot)];
		      key2 = Chan.hp1_map[ind1][Hash2(pind, rind, Space.indtot)];
		      ind = key1 * Chan.nhp1[ind1] + key2;
		      Ints.S_ME1.V16[ind1][ind] += TBME;
		      
		      plus(tb, Space.qnums[pind], Space.qnums[qind]);
		      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		      key1 = Chan.pp_map[ind1][Hash2(rind, sind, Space.indtot)];
		      key3 = Chan.hp_map[ind1][Hash2(pind, qind, Space.indtot)];
		      ind = key1 * Chan.nhp[ind1] + key3;
		      Ints.S_ME1.V20[ind1][ind] += TBME;
		    }
		    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
		      ind1 = Chan.indvec[rind];
		      key1 = Chan.h_map[ind1][rind];
		      key2 = Chan.hhp_map[ind1][Hash3(pind, qind, sind, Space.indtot)];
		      ind = key1 * Chan.nhhp[ind1] + key2;
		      Ints.S_ME1.V12[ind1][ind] += TBME;
		      ind = key2 * Chan.nh[ind1] + key1;
		      Ints.S_ME1.V18[ind1][ind] += TBME;

		      ind1 = Chan.indvec[pind];
		      key1 = Chan.h_map[ind1][pind];
		      key2 = Chan.hhp1_map[ind1][Hash3(qind, rind, sind, Space.indtot)];
		      ind = key2 * Chan.nh[ind1] + key1;
		      Ints.S_ME1.V14[ind1][ind] += TBME;

		      minus(tb, Space.qnums[pind], Space.qnums[sind]);
		      ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
		      key1 = Chan.hh1_map[ind1][Hash2(rind, qind, Space.indtot)];
		      key2 = Chan.hp1_map[ind1][Hash2(pind, sind, Space.indtot)];
		      ind = key1 * Chan.nhp1[ind1] + key2;
		      Ints.S_ME1.V15[ind1][ind] += TBME;
		      
		      plus(tb, Space.qnums[pind], Space.qnums[qind]);
		      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		      key1 = Chan.hh_map[ind1][Hash2(pind, qind, Space.indtot)];
		      key3 = Chan.hp_map[ind1][Hash2(rind, sind, Space.indtot)];
		      ind = key3 * Chan.nhh[ind1] + key1;
		      Ints.S_ME1.V19[ind1][ind] += TBME;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
