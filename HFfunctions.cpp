#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

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
