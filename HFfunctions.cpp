#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

HF_Matrix_Elements::HF_Matrix_Elements(HF_Channels &Chan)
{
  int ntb;
  int V_total = 0;
  Index = new int[Chan.size1];
  for(int i = 0; i < Chan.size1; ++i){
    ntb = Chan.ntb[i];
    Index[i] = V_total;
    V_total += ntb * ntb;
  }
  V = new double[V_total];
  for(int i = 0; i < V_total; ++i){
    V[i] = 0.0;
  }
}

void HF_Matrix_Elements::delete_struct(HF_Channels &Chan)
{
  delete[] V;
  delete[] Index;
}

one_body HF_Channels::ob_state(int &chan3, int &ind3)
{
  return ob_vec[ob_index[chan3] + ind3];
}

two_body HF_Channels::tb_state(int &chan1, int &ind1)
{
  return tb_vec[tb_index[chan1] + ind1];
}

HF_Channels::HF_Channels(Input_Parameters &Parameters, Model_Space &Space)
{
  int chan1, chan3;
  int ob_total, tb_total;
  size1 = Space.size_2b_dir;
  size3 = Space.size_1b;
  
  qnums1 = new State[size1];
  qnums3 = new State[size3];

  p_chan_setup(Parameters,Space, ob_total,nob,ob_index,ob_map,ob_vec, size3,qnums3, -1);  // ob
  pq_chan_setup(Parameters,Space, tb_total,ntb,tb_index,tb_map,tb_vec, size1,qnums1, size3,qnums3, nob,ob_vec,ob_index,nob,ob_index,ob_vec);  // tb

  double memory = 0.0;
  int doubsize = sizeof(double);
  for(chan3 = 0; chan3 < size3; ++chan3){
    memory += 3 * doubsize * nob[chan3] * nob[chan3]; // 3 HF vectors
    memory += 3 * doubsize * nob[chan3]; // 3 HF energies
  }
  for(chan1 = 0; chan1 < size1; ++chan1){
    memory += doubsize * ntb[chan1] * ntb[chan1]; // ME.V
  }
  std::cout << std::endl << "Estimated HF-Memory = " << memory/1000000.0 << " MB" << std::endl << std::endl;
}

void HF_Channels::delete_struct()
{
  delete[] ob_vec;
  delete[] ob_index;
  delete[] qnums3;
  delete[] nob;
  delete[] ob_map;

  delete[] tb_vec;
  delete[] tb_index;
  delete[] qnums1;
  delete[] ntb;
  delete[] tb_map;
}

Single_Particle_States::Single_Particle_States(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan)
{
  int nh, nob, ob_ind, ind, chan;
  int vector_ind, energy_ind;
  int length1, length2;

  hole_size = new int[Chan.size3];
  vector_size = new int[Chan.size3];
  vector_index = new int[Chan.size3];
  energy_index = new int[Chan.size3];

  length1 = 0;
  length2 = 0;
  for(chan = 0; chan < Chan.size3; ++chan){
    nh = 0;
    nob = Chan.nob[chan];
    for(ind = 0; ind < nob; ++ind){
      ob_ind = Chan.ob_state(chan, ind).v1;
      if( Space.qnums[ob_ind].type == 0 ){ ++nh; }
    }
    hole_size[chan] = nh;
    vector_size[chan] = nob;
    vector_index[chan] = length1;
    energy_index[chan] = length2;
    length1 += nob * nob;
    length2 += nob;
  }
  vector_length = length1;
  energy_length = length2;

  vectors = new double[vector_length];
  energies = new double[energy_length];
  for(ind = 0; ind < vector_length; ++ind){ vectors[ind] = 0.0; }

  for(chan = 0; chan < Chan.size3; ++chan){
    vector_ind = vector_index[chan];
    energy_ind = energy_index[chan];
    for(ind = 0; ind < vector_size[chan]; ++ind){
      ob_ind = Chan.ob_state(chan, ind).v1;
      vectors[vector_ind + (ind * vector_size[chan] + ind)] = 1.0;
      energies[energy_ind + ind] = Space.qnums[ob_ind].energy;
    }
  }
}

void Single_Particle_States::delete_struct(HF_Channels &Chan)
{
  delete[] hole_size;
  delete[] vector_size;
  delete[] vector_index;
  delete[] energy_index;
  delete[] vectors; 
  delete[] energies;
}

void Hartree_Fock_States(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &HF, HF_Matrix_Elements &ME)
{
  double Bshift0 = 0.0;
  for(int i = 0; i < Space.num_states; ++i){ Bshift0 += Space.qnums[i].energy; }
  Bshift0 /= (Space.num_states);

  double width = 1.0;
  int DIIS_count = 0, Rand_count = 0;

  double error = 1e20, error2 = 1e20;
  int size3 = Chan.size3;
  int ind = 0;
  int *tempn;
  int ob_ind, energy_ind;
  Single_Particle_States HF0, HF2;
  int badcount = 0;
  int badlimit = 3;

  //// for DIIS /////////
  double DIISstart = 0.1;
  int maxl = 10;
  int N = 0;
  double *p = NULL;
  double *delp = NULL;
  double *tempdelp = NULL;
  double *B = NULL;
  Initialize_DIIS_HF(Chan, HF, p, delp, tempdelp, B, maxl);
  ///////////////////////

  //   Initialize Mesh for KE   //
  if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
    setup_ho_cutoff(Parameters, Space, ME);
    ME.ra = new double[ME.n_rel];
    ME.wra = new double[ME.n_rel];
    gauss_legendre(0.0, ME.cutoff, ME.ra, ME.wra, ME.n_rel);
    ME.hol = new double[ME.n_rel * (Space.qmaxs.n+1) * (Space.qmaxs.l+1)];
    ho_wfunction(Parameters, Space, ME);
  }
  ///////////////////////////////

  //   Initialize HF0, HF2, and tempn   //
  HF0 = Single_Particle_States(Parameters, Space, Chan);
  HF2 = Single_Particle_States(Parameters, Space, Chan);
  tempn = new int[HF.energy_length];
  for(int chan = 0; chan < Chan.size3; ++chan){
    for(int i = 0; i < HF.vector_size[chan]; ++i){
      tempn[HF.energy_index[chan] + i] = Space.qnums[Chan.ob_state(chan, i).v1].n;
    }
  }

  while(error > 1e-12 && ind < 10000){
    Hartree_Fock_Step(Parameters, Space, Chan, HF, HF2, ME, Bshift0, error);
    if( ind > 1000 && double(Rand_count)/double(ind) > 0.9 ){
      std::cerr << std::endl << ind << " : HF Solution Not Converged!!" << std::endl; exit(1);
    }
    std::cout << "HF-ind = " << ind << ", error = " << error << ", DIIS_count = " << DIIS_count << ", Rand_count = " << Rand_count << std::endl;

    if( error < error2 || error > 1e-4 ){
      badcount = 0;
      for(int i = 0; i < HF.energy_length; ++i){
	HF0.energies[i] = HF.energies[i];
	HF.energies[i] = HF2.energies[i];
      }
      for(int i = 0; i < HF.vector_length; ++i){
	HF0.vectors[i] = HF.vectors[i];
	HF.vectors[i] = HF2.vectors[i];
      }
      error2 = error;
      if( error < DIISstart && error > 1e-12 ){ Perform_DIIS_HF(Chan, HF, HF0, p, delp, tempdelp, B, N, maxl, DIIS_count); }
    }
    else if(badcount < badlimit){
      ++badcount;
      for(int i = 0; i < HF.energy_length; ++i){ HF.energies[i] = HF2.energies[i]; }
      for(int i = 0; i < HF.vector_length; ++i){ HF.vectors[i] = HF2.vectors[i]; }
      error2 = error;
    }
    else{
      badcount = 0;
      ++Rand_count;
      if(error2 > 1.0){ width = 0.0001; }
      else{ width = 0.0001 * error2; }
      Random_Step_HF(Parameters, Space, Chan, HF0, HF, HF2, ME, width, error2, Bshift0);
    }
    ++ind;
  }

  if( error > 1e-12 ){
    std::cout << Parameters.Shells << ", " << Parameters.Pshells << ", " << Parameters.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". HF Solution Not Converged!!" << std::endl;
  }

  for(int chan = 0; chan < size3; ++chan){
    energy_ind = HF.energy_index[chan];
    for(ind = 0; ind < HF.vector_size[chan]; ++ind){
      ob_ind = Chan.ob_state(chan, ind).v1;
      Space.qnums[ob_ind] = Chan.qnums3[chan];
      Space.qnums[ob_ind].n = tempn[energy_ind + ind];
      Space.qnums[ob_ind].energy = HF.energies[energy_ind + ind];
      if( ind < HF.hole_size[chan] ){ Space.qnums[ob_ind].type = 0; }
      else{ Space.qnums[ob_ind].type = 1; }
    }
  }
  delete[] tempn;

  HF0.delete_struct(Chan);
  HF2.delete_struct(Chan);
  Delete_DIIS_HF(Chan, p, delp, tempdelp, B, maxl);

  if(Parameters.basis == "finite_HO"){
    for(int i = 0; i < Space.num_states; ++i){
      std::cout << Space.qnums[i].n << " " << Space.qnums[i].ml << " " << Space.qnums[i].m << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }
  }
  else if(Parameters.basis == "finite_J"){
    for(int i = 0; i < Space.num_states; ++i){
      std::cout << Space.qnums[i].n << " " << Space.qnums[i].par << " " << Space.qnums[i].j << " " << Space.qnums[i].t << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }
  }
  else if(Parameters.basis == "finite_M"){
    for(int i = 0; i < Space.num_states; ++i){
      std::cout << Space.qnums[i].n << " " << Space.qnums[i].par << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " : " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }
  }

  //   Delete Mesh for KE   //
  if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){
    delete[] ME.ra;
    delete[] ME.wra;
    delete[] ME.hol;
  }
  ///////////////////////////
}

void Hartree_Fock_Step(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &Bshift0, double &error)
{
  double temperror = 0.0;
  double tempnorm = 0.0;
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  int lwork, info; // Parameters for Diagonaliztion
  double *w;
  double *work;
  double *fock;
  int size, length1;
  int vector_ind, energy_ind;
  int maxind;
  double maxc;
  jobz = 'V';
  uplo = 'U';

  double *Bshift = new double[HF.energy_length];
  for(int i = 0; i < HF.energy_length; ++i){ Bshift[i] = Bshift0; }

  //Make Fock Matrix
  for(int chan = 0; chan < Chan.size3; ++chan){
    vector_ind = HF.vector_index[chan];
    energy_ind = HF.energy_index[chan];
    size = HF.vector_size[chan];
    fock = new double[size * size];
    length1 = int(0.5 * size * (size + 1));
    #pragma omp parallel
    {
      int m, l, mind, nind, lind, kind;
      int key1, key2, chan1, ind2, vector_ind2;
      int size2, length2, minj;
      double term, KE, HF_0;
      State tb;
      #pragma omp for schedule(static)
      for(int ml = 0; ml < length1; ++ml){
	m = std::floor((2*size - 1 - std::sqrt(1 + 4*size + 4*size*size - 8*ml))/2) + 1;
	mind = Chan.ob_state(chan, m).v1;
	length2 = int(0.5 * m * (2*size - m + 1));
	l = int(m + ml - length2);
	lind = Chan.ob_state(chan, l).v1;

	HF_0 = 0.0;
	if(m == l){
	  //if(Parameters.basis == "finite_J" || Parameters.basis == "finite_JM"){ KE = kinetic_energy(Parameters, Space, ME, mind, lind); }
	  //else{ KE = Space.qnums[mind].energy; } // Add diagonal elements to fock matrices
	  KE = Space.qnums[mind].energy;
	}
	else{ KE = 0.0; }
	fock[size * m + l] = KE;
	
	for(int chan2 = 0; chan2 < Chan.size3; ++chan2){
	  vector_ind2 = HF.vector_index[chan2];
	  size2 = HF.vector_size[chan2];
	  for(int beta = 0; beta < HF.hole_size[chan2]; ++beta){ // Sum over occupied levels
	    for(int n = 0; n < size2; ++n){
	      nind = Chan.ob_state(chan2, n).v1;
	      for(int k = 0; k < size2; ++k){
		kind = Chan.ob_state(chan2, k).v1;
		plus(tb, Space.qnums[mind], Space.qnums[nind]);
		minj = abs(Chan.qnums3[chan].j - Chan.qnums3[chan2].j);
		while(tb.j >= minj){
		  if( (mind == nind || lind == kind) && tb.j%4 != 0 ){ tb.j -= 2; continue; }
		  chan1 = Space.ind_2b_dir(Parameters.basis, tb);
		  key1 = Chan.tb_map[chan1][Space.hash2(mind, nind, tb.j)];
		  key2 = Chan.tb_map[chan1][Space.hash2(lind, kind, tb.j)];
		  ind2 =  ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  term = HF.vectors[vector_ind2 + (beta * size2 + n)] * HF.vectors[vector_ind2 + (beta * size2 + k)] * ME.V[ind2];
		  term *= (tb.j + 1.0)/(Chan.qnums3[chan].j + 1.0);
		  //std::cout << "< " << mind << ", " << nind << " || " << lind << ", " << kind << " >^ " << std::setprecision(14) << tb.j << " : " << (tb.j + 1.0)/(Chan.qnums3[chan].j + 1.0) << ", " << HF.vectors[vector_ind2 + (beta * size2 + n)] << ", " << HF.vectors[vector_ind2 + (beta * size2 + k)] << " : " << ME.V[ind2] << std::endl;
		  HF_0 += term;
		  fock[size * m + l] += term;
		  tb.j -= 2;
		}
	      }
	    }
	  }
	}
	//std::cout << "a = " << mind+1 << ", c = " << lind+1 << std::setprecision(14) << " : e_kin = " << KE << ", hf = " << HF_0 << std::endl;
	for(int beta = HF.hole_size[chan]; beta < HF.vector_size[chan]; ++beta){ // Sum over unoccupied levels
	  fock[size * m + l] += HF.vectors[vector_ind + (beta * size + m)] * HF.vectors[vector_ind + (beta * size + l)] * Bshift[energy_ind + beta];
	}
	if(m != l){ fock[size * l + m] = fock[size * m + l]; }
      }
    }

    lda = size;
    lwork = (3+2)*size;
    w = new double[lda];
    work = new double[lwork];
    for(int j = 0; j < lda; ++j){ w[j] = 0.0; }
    for(int j = 0; j < lwork; ++j){ work[j] = 0.0; }
    if(size != 0){ dsyev_(&jobz, &uplo, &size, fock, &lda, w, work, &lwork, &info); }

    double deltaE;
    for(int j = HF.hole_size[chan]; j < HF.vector_size[chan]; ++j){ w[j] -= Bshift[energy_ind + j]; } // Add back Level-shift parameter
    for(int j = HF.hole_size[chan]; j < HF.vector_size[chan]; ++j){
      deltaE = fabs(w[j] - HF.energies[energy_ind + j]);
      if( deltaE > Bshift0 || HF.hole_size[chan] == 0 || HF.vector_size[chan] == HF.hole_size[chan]){ Bshift[energy_ind + j] = Bshift0 + deltaE; }
      else{ Bshift[energy_ind + j] = Bshift0 + w[HF.hole_size[chan]] - w[HF.hole_size[chan] - 1]; }
    }

    for(int j = 0; j < size; ++j){
      HF2.energies[energy_ind + j] = w[j];
      maxind = -1;
      maxc = 0.0;
      for(int k = 0; k < size; ++k){
	if( fabs(fock[size*j + k]) > maxc ){
	  maxc = fabs(fock[size*j + k]);
	  maxind = k;
	}
      }
      if( HF.vectors[vector_ind + (size*j + maxind)]/fock[size * j + maxind] < 0.0 ){
	for(int k = 0; k < size; ++k){ fock[size*j + k] *= -1.0; }
      }
    }

    for(int j = 0; j < size; ++j){
      for(int k = 0; k < size; ++k){
	temperror += std::pow(fabs(HF.vectors[vector_ind + (size*j + k)]) - fabs(fock[size * j + k]), 2.0);
	tempnorm += std::pow(fock[size * j + k], 2.0);
	HF2.vectors[vector_ind + (size*j + k)] = fock[size * j + k];
      }
    }
    delete[] fock;
    delete[] w;
    delete[] work;
  }
  delete[] Bshift;
  error = std::sqrt(temperror/tempnorm);
}

void Initialize_DIIS_HF(HF_Channels &Chan, Single_Particle_States &HF, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl)
{
  int size = HF.vector_length;
  B = new double[1];
  B[0] = 0.0;

  p = new double[maxl * size];
  delp = new double[maxl* size];
  tempdelp = new double[size];
  for(int l = 0; l < maxl; ++l){
    for(int ind = 0; ind < size; ++ind){
      p[l * size + ind] = 0.0;
      delp[l * size + ind] = 0.0;
    }
  }
  for(int ind = 0; ind < size; ++ind){
    tempdelp[ind] = 0.0;
  }
}

void Delete_DIIS_HF(HF_Channels &Chan, double *&p, double *&delp, double *&tempdelp, double *&B, int &maxl)
{
  delete[] B;
  delete[] p;
  delete[] delp;
  delete[] tempdelp;
}

void Perform_DIIS_HF(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF0, double *&p, double *&delp, double *&tempdelp, double *&B, int &N, int &maxl, int &DIIS_count)
{
  double norm = 0.0;
  double norm0 = 0.0;
  double temp;
  int size, ind, chan;
  int ortho;
  int P = N + 1;
  int lwork;
  int *ipiv;
  double *work;
  int info;
  double *B2;
  double dot_prod, norm1, norm2;
  int maxind;
  double maxc;
  size = HF.vector_length;

  // Fill temp_delp
  for(ind = 0; ind < size; ++ind){
    tempdelp[ind] = (HF.vectors[ind] - HF0.vectors[ind]);
    norm0 += HF0.vectors[ind] * HF0.vectors[ind];
  }
  norm0 = std::sqrt(norm0);
  for(ind = 0; ind < size; ++ind){
    tempdelp[ind] /= norm0;
    norm += tempdelp[ind] * tempdelp[ind];
  }
  
  // check orthogonality of tempdelp
  ortho = 1;
  for(int l = 0; l < N; ++l){
    dot_prod = 0.0;
    norm1 = 0.0;
    norm2 = 0.0;
    for(ind = 0; ind < size; ++ind){
      dot_prod += HF.vectors[ind] * p[l * size + ind];
      norm1 += HF.vectors[ind] * HF.vectors[ind];
      norm2 += p[l * size + ind] * p[l * size + ind];
    }
    dot_prod /= std::sqrt(norm1 * norm2);
    //std::cout << "!! norm <? B[P*l+l] = " << norm << " <? " << B[P*l + l] << ", dot_prod <? (1-tol) = " << dot_prod << " <? " << (1.0 - norm) << std::endl;
    if(norm > B[P * l + l] || dot_prod > (1.0 - norm)){ ortho = 0; break; }
  }
  if(ortho != 1){ return; }

  ++DIIS_count;
  if(N < maxl){ Update_B1_HF(Chan, HF, N, p, delp, tempdelp, B); }
  else{ Update_B2_HF(Chan, HF, N, p, delp, tempdelp, B); }
  P = N + 1;

  // copy B into B2 for inversion
  B2 = new double[P * P];
  for(int j = 0; j < P*P; ++j){ B2[j] = B[j]; }

  ipiv = new int[P];
  work = new double[sizeof(double) * P];
  lwork = sizeof(double) * P;
  info = 0;
  dgetrf_(&P, &P, B2, &P, ipiv, &info);
  dgetri_(&P, B2, &P, ipiv, work, &lwork, &info);

  for(ind = 0; ind < size; ++ind){
    temp = 0.0;
    for(int l = 0; l < N; ++l){ temp += -1.0 * B2[P * l + N] * p[l * size + ind]; }
    HF.vectors[ind] = temp;
  }

  for(chan = 0; chan < Chan.size3; ++chan){
    GramSchmidt(HF.vectors + HF.vector_index[chan], HF.vector_size[chan]);
    for(int j = 0; j < HF.vector_size[chan]; ++j){
      maxind = -1;
      maxc = 0.0;
      for(int k = 0; k < HF.vector_size[chan]; ++k){
	if( fabs(HF.vectors[HF.vector_index[chan] + (HF.vector_size[chan]*j + k)]) > maxc ){
	  maxc = fabs(HF.vectors[HF.vector_index[chan] + (HF.vector_size[chan]*j + k)]);
	  maxind = k;
	}
      }
      if( HF0.vectors[HF.vector_index[chan] + (HF.vector_size[chan]*j + maxind)]/HF.vectors[HF.vector_index[chan] + (HF.vector_size[chan]*j + maxind)]
	  < 0.0 ){
	for(int k = 0; k < HF.vector_size[chan]; ++k){ HF.vectors[HF.vector_index[chan] + (HF.vector_size[chan]*j + k)] *= -1.0; }
      }
    }
  }

  delete[] B2;
  delete[] ipiv;
  delete[] work;
}

void Update_B1_HF(HF_Channels &Chan, Single_Particle_States &HF, int &N, double *&p, double *&delp, double *&tempdelp, double *&B)
{
  double *B2; // used to copy B into updated B
  int P = N+1;
  int ind;
  int size = HF.vector_length;
  for(ind = 0; ind < size; ++ind){
    p[N * size + ind] = HF.vectors[ind];
    delp[N * size + ind] = tempdelp[ind];
  }
  B2 = new double[(P+1) * (P+1)];
  for(int j = 0; j < N; ++j){
    for(int k = 0; k < N; ++k){
      B2[(P+1) * j + k] = B[P * j + k];
    }
  }
  for(int l = 0; l < P; ++l){
    B2[(P+1) * N + l] = 0.0;
    if(l != N){ B2[(P+1) * l + N] = 0.0; }
    for(ind = 0; ind < size; ++ind){
      B2[(P+1) * N + l] += delp[N * size + ind] * delp[l * size + ind];
      if(l != N){ B2[(P+1) * l + N] += delp[l * size + ind] * delp[N * size + ind]; }
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
  for(int ind = 0; ind < P * P; ++ind){ B[ind] = B2[ind]; }
  delete[] B2;
}

void Update_B2_HF(HF_Channels &Chan, Single_Particle_States &HF, int &N, double *&p, double *&delp, double *&tempdelp, double *&B)
{
  int maxind, ind, size;
  double maxnorm;
  int P = N + 1;
  size = HF.vector_length;

  // Find largest norm of B to remove that vector
  maxind = -1;
  maxnorm = 0.0;
  for(int j = 0; j < N; ++j){
    if(B[P * j + j] > maxnorm){
      maxind = j;
      maxnorm = B[P * j + j];
    }
  }
  // Replace maxnorm vector with new vector
  for(ind = 0; ind < size; ++ind){
    p[maxind * size + ind] = HF.vectors[ind];
    delp[maxind * size + ind] = tempdelp[ind];
  }

  for(int l = 0; l < N; ++l){
    B[P * maxind + l] = 0.0;
    if(l != maxind){ B[P * l + maxind] = 0.0; }	
    for(ind = 0; ind < size; ++ind){
      B[P * maxind + l] += delp[maxind * size + ind] * delp[l * size + ind];
      if(l != maxind){ B[P * l + maxind] += delp[l * size + ind] * delp[maxind * size + ind]; }
    }
  }
}

void Random_Step_HF(Input_Parameters &Parameters, Model_Space &Space, HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &width, double &error2, double &Bshift)
{
  double error;
  int go = 0;
  int ind = 0;
  while(go == 0){
    ++ind;
    Randomize_HF(Chan, HF0, HF, width);
    Hartree_Fock_Step(Parameters, Space, Chan, HF, HF2, ME, Bshift, error);
    std::cout << "!!  " << ind << ", " << error << " " << error2 << ", " << width << std::endl;
    if(error < error2){
      ++go;
      for(int i = 0; i < HF.energy_length; ++i){
	HF0.energies[i] = HF.energies[i];
	HF.energies[i] = HF2.energies[i];
      }
      for(int i = 0; i < HF.vector_length; ++i){
	HF0.vectors[i] = HF.vectors[i];
	HF.vectors[i] = HF2.vectors[i];
      }
      error2 = error;
    }
    width *= 0.975;
    if(width < 10-12){ width = 10e-12; }
    if(ind >= 1000){ std::cerr << std::endl << "HF Random Step Unsuccessful!!" << std::endl; exit(1); }
  }
}

void Randomize_HF(HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, double &width)
{
  int chan, ind, vector_ind, size, maxind;
  double tempc, maxc;
  double rand;
  size_t key;
  std::unordered_map<size_t,double> c_map;
  for(ind = 0; ind < HF.vector_length; ++ind){
    tempc = HF0.vectors[ind];
    if(fabs(tempc) > 1.0e-16){
      rand = rand_normal(0.0, width * fabs(tempc));
      key = std::hash<float>{}(float(fabs(tempc)));
      c_map[key] = rand;
    }
  }

  for(ind = 0; ind < HF.vector_length; ++ind){
    tempc = HF0.vectors[ind];
    key = std::hash<float>{}(float(fabs(tempc)));
    if(tempc > 1.0e-16){ tempc += c_map[key]; }
    else if(tempc < -1.0e-16){ tempc -= c_map[key]; }
    HF.vectors[ind] = tempc;
  }

  for(chan = 0; chan < Chan.size3; ++chan){
    vector_ind = HF.vector_index[chan];
    size = HF.vector_size[chan];
    GramSchmidt(HF.vectors + vector_ind, size);
    for(int j = 0; j < size; ++j){
      maxind = -1;
      maxc = 0.0;
      for(int k = 0; k < size; ++k){
	if( fabs(HF.vectors[vector_ind + (size*j + maxind)]) > maxc ){
	  maxc = fabs(HF.vectors[vector_ind + (size*j + maxind)]);
	  maxind = k;
	}
      }
      if( HF0.vectors[vector_ind + (size*j + maxind)]/HF.vectors[vector_ind + (size*j + maxind)] < 0.0 ){
	for(int k = 0; k < size; ++k){ HF.vectors[vector_ind + (size*j + maxind)] *= -1.0; }
      }
    }
  }
}

void Convert_To_HF_Matrix_Elements(Input_Parameters &Parameters, HF_Channels &Chan, Model_Space &Space, Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  int size, length1, matlength; // max length of M-Scheme indicies, length of J_ME
  double *M1; // Matrices of coefficients
  double *C;
  char transa, transb;
  double alpha1, beta1;

  for(int chan = 0; chan < Chan.size1; ++chan){
    size = Chan.ntb[chan];
    if(size == 0){ continue; }
    matlength = size * size;
    length1 = int(0.5 * size * (size + 1));
    M1 = new double[matlength];
    C = new double[matlength];
    for(int i = 0; i < matlength; ++i){
      M1[i] = 0.0;
      C[i] = 0.0;
    }

    #pragma omp parallel
    {
      int pq, rs, length2;
      two_body pq_tb, rs_tb;
      int p, q, r, s;
      int p_chan, q_chan, r_chan, s_chan;
      int p_ind, q_ind, r_ind, s_ind;
      int vec_ind_r, vec_ind_s, p_size, q_size;
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < length1; ++pqrs){
	pq = std::floor((2*size - 1 - std::sqrt(1 + 4*size + 4*matlength - 8*pqrs))/2) + 1;
	length2 = int(0.5 * pq * (2*size - pq + 1));
	rs = int(pq + pqrs - length2);
	pq_tb = Chan.tb_state(chan, pq);
	rs_tb = Chan.tb_state(chan, rs);
	p = pq_tb.v1;
	q = pq_tb.v2;
	r = rs_tb.v1;
	s = rs_tb.v2;
	p_chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	q_chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	r_chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	s_chan = Space.ind_1b(Parameters.basis, Space.qnums[s]);
	p_ind = Chan.ob_map[p_chan][p];
	q_ind = Chan.ob_map[q_chan][q];
	r_ind = Chan.ob_map[r_chan][r];
	s_ind = Chan.ob_map[s_chan][s];
	if(p_chan != r_chan || q_chan != s_chan){ continue; }
	vec_ind_r = States.vector_index[r_chan];
	vec_ind_s = States.vector_index[s_chan];
	p_size = States.vector_size[p_chan];
	q_size = States.vector_size[q_chan];
	M1[size * pq + rs] = States.vectors[vec_ind_r + (r_ind * p_size + p_ind)] * States.vectors[vec_ind_s + (s_ind * q_size + q_ind)];
	if( pq != rs ){ M1[size * rs + pq] = States.vectors[vec_ind_r + (p_ind * p_size + r_ind)] * States.vectors[vec_ind_s + (q_ind * q_size + s_ind)]; }
      }
    }

    /*#pragma omp parallel
    {
      int pq, rs, length2;
      two_body pq_tb, rs_tb;
      int p, q, r, s;
      int p_chan, q_chan, r_chan, s_chan;
      int p_ind, q_ind, r_ind, s_ind;
      int vec_ind_r, vec_ind_s, p_size, q_size;
      double dir, exch, term, norm;
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < matlength; ++pqrs){
	rs = pqrs%size;
	pq = int((pqrs - rs)/size);
	pq_tb = Chan.tb_state(chan, pq);
	rs_tb = Chan.tb_state(chan, rs);
	p = pq_tb.v1;
	q = pq_tb.v2;
	r = rs_tb.v1;
	s = rs_tb.v2;
	norm = 1.0;
	if(p == q){ norm /= std::sqrt(2.0); }
	if(r == s){ norm /= std::sqrt(2.0); }
	ME.V[ME.Index[chan] + pq*size + rs] *= norm; // renormalize ME
	p_chan = Space.ind_1b(Parameters.basis, Space.qnums[p]);
	q_chan = Space.ind_1b(Parameters.basis, Space.qnums[q]);
	r_chan = Space.ind_1b(Parameters.basis, Space.qnums[r]);
	s_chan = Space.ind_1b(Parameters.basis, Space.qnums[s]);
	p_ind = Chan.ob_map[p_chan][p];
	q_ind = Chan.ob_map[q_chan][q];
	r_ind = Chan.ob_map[r_chan][r];
	s_ind = Chan.ob_map[s_chan][s];
	vec_ind_r = States.vector_index[r_chan];
	vec_ind_s = States.vector_index[s_chan];
	p_size = States.vector_size[p_chan];
	q_size = States.vector_size[q_chan];
	dir = 0.0;
	exch = 0.0;
	if(p_chan == r_chan && q_chan == s_chan){
	  dir = States.vectors[vec_ind_r + (r_ind * p_size + p_ind)] * States.vectors[vec_ind_s + (s_ind * q_size + q_ind)];
	}
	if(p_chan == s_chan && q_chan == r_chan && Chan.qnums1[chan].t != 0){
	  exch = States.vectors[vec_ind_r + (r_ind * q_size + q_ind)] * States.vectors[vec_ind_s + (s_ind * p_size + p_ind)];
	}
	M1[size * pq + rs] = (dir - exch * std::pow(-1.0, int(0.5*(Chan.qnums3[p_chan].j + Chan.qnums3[q_chan].j - Chan.qnums1[chan].j)))) * norm;
      }
      }*/

    /*if(chan == 20){
      std::cout << "tb_vec[" << chan << "]" << std::endl;
      for(int tb1 = 0; tb1 < size; ++tb1){ std::cout << Chan.tb_vec[Chan.tb_index[chan] + tb1].v1 << "," << Chan.tb_vec[Chan.tb_index[chan] + tb1].v2 << "  "; }
      std::cout << std::endl << std::endl;
      std::cout << std::setprecision(8);
      std::cout << "ME.V[" << chan << "]" << std::endl;
      for(int tb1 = 0; tb1 < size; ++tb1){
	for(int tb2 = 0; tb2 < size; ++tb2){
	  std::cout << ME.V[ME.Index[chan] + tb1*size + tb2] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      
      std::cout << "M1[" << chan << "]" << std::endl;
      for(int tb1 = 0; tb1 < size; ++tb1){
	for(int tb2 = 0; tb2 < size; ++tb2){
	  std::cout << M1[tb1*size + tb2] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      }*/

    transa = 'N';
    transb = 'T';
    alpha1 = 1.0;
    beta1 = 0.0;
    dgemm_NN(ME.V + ME.Index[chan], M1, C, &size, &size, &size, &alpha1, &beta1, &transa, &transa);
    dgemm_TN(M1, C, ME.V + ME.Index[chan], &size, &size, &size, &alpha1, &beta1, &transb, &transa);

    /*if(chan == 20){
      std::cout << "HF.ME.V[" << chan << "]" << std::endl;
      for(int tb1 = 0; tb1 < size; ++tb1){
	for(int tb2 = 0; tb2 < size; ++tb2){
	  std::cout << ME.V[ME.Index[chan] + tb1*size + tb2] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      }*/

    delete[] M1;
    delete[] C;
  }
    
  // print ME to file
  /*std::ofstream interactionfile; // print matrix elements to file
  int p, q, r, s;
  double tbme;
  interactionfile.open("interaction.dat");
  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    length = Chan.ntb[chan1];
    for(int pq = 0; pq < length; ++pq){
      p = Chan.tb_vec[chan1][2*pq];
      q = Chan.tb_vec[chan1][2*pq + 1];
      if(p >= q){ continue; }
      for(int rs = 0; rs < length; ++rs){
	r = Chan.tb_vec[chan1][2*rs];
	s = Chan.tb_vec[chan1][2*rs + 1];
	if(r >= s){ continue; }
	if(p > r){ continue; }
	if(p == r && q > s){ continue; }
	tbme = ME.V[chan1][pq * length + rs];
	if(std::fabs(tbme) < 1.0e-18){ continue; }
	interactionfile << std::setw(5) << p << std::setw(5) << q << std::setw(5) << r << std::setw(5) << s;
	interactionfile << std::setw(20) << std::setprecision(12) << tbme << "\n";
	//interactionfile << std::setw(4) << Space.qnums[p].n << std::setw(4) << Space.qnums[p].ml << std::setw(4) << Space.qnums[p].m;
	//interactionfile << std::setw(4) << Space.qnums[q].n << std::setw(4) << Space.qnums[q].ml << std::setw(4) << Space.qnums[q].m;
	//interactionfile << std::setw(4) << Space.qnums[r].n << std::setw(4) << Space.qnums[r].ml << std::setw(4) << Space.qnums[r].m;
	//interactionfile << std::setw(4) << Space.qnums[s].n << std::setw(4) << Space.qnums[s].ml << std::setw(4) << Space.qnums[s].m;
	//interactionfile << std::setw(24) << std::setprecision(16) << tbme << "\n";
      }
    }
  }
  interactionfile.close();*/
}
