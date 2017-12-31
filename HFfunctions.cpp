#include "HFfunctions.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "INTfunctions.hpp"
#include "DIISfunctions.hpp"
#include "MATHfunctions.hpp"
#include "AngMom.hpp"

void HF_Matrix_Elements::Build(HF_Channels &HF_Chan)
{
  int ntb;
  int V_total = 0;
  this->Index = new int[HF_Chan.size1];
  for(int i = 0; i < HF_Chan.size1; ++i){
    ntb = HF_Chan.ntb[i];
    this->Index[i] = V_total;
    V_total += ntb * ntb;
  }
  this->V = new double[V_total];
  for(int i = 0; i < V_total; ++i){
    this->V[i] = 0.0;
  }
}

void HF_Matrix_Elements::Delete()
{
  delete[] this->V;
  delete[] this->Index;

  if( PAR.ME3 == 1 ){
    delete[] this->V3;
    for(int nlj1 = 0; nlj1 < SPB.nljMax3; nlj1++){
      if( SPB.e_nlj3[nlj1] > SPB.E3Max ){ break; }
      
      for(int nlj2 = 0; nlj2 < nlj1; nlj2++){
	if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] > SPB.E3Max ){ break; }
	
	for(int nlj3 = 0; nlj3 < nlj2; nlj3++){
	  if( SPB.e_nlj3[nlj1] + SPB.e_nlj3[nlj2] + SPB.e_nlj3[nlj3] > SPB.E3Max ){ break; }
	  
	  for(int nnlj1 = 0; nnlj1 < nlj1; nnlj1++){
	    if( SPB.e_nlj3[nnlj1] > SPB.E3Max ){ break; }
	    
	    for(int nnlj2 = 0; nnlj2 < ((nlj1==nnlj1) ? nlj2 : nnlj1); nnlj2++){
	      if( SPB.e_nlj3[nnlj1] + SPB.e_nlj3[nnlj2] > SPB.E3Max ){ break; }
	      
	      if( this->ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2] != NULL ){ delete[] this->ME3Idx[nlj1][nlj2][nlj3][nnlj1][nnlj2]; }
	    }
	    if( this->ME3Idx[nlj1][nlj2][nlj3][nnlj1] != NULL ){ delete[] this->ME3Idx[nlj1][nlj2][nlj3][nnlj1]; }
	  }
	  if( this->ME3Idx[nlj1][nlj2][nlj3] != NULL ){ delete[] this->ME3Idx[nlj1][nlj2][nlj3]; }
	}
	if( this->ME3Idx[nlj1][nlj2] != NULL ){ delete[] this->ME3Idx[nlj1][nlj2]; }
      }
      if( this->ME3Idx[nlj1] != NULL ){ delete[] this->ME3Idx[nlj1]; }
    }
    if( this->ME3Idx != NULL ){ delete[] this->ME3Idx; }
    delete[] this->V3mon;
    delete[] this->V3mon_index;
  }
}

one_body HF_Channels::ob_state(int &chan3, int &ind3)
{
  return ob_vec[ob_index[chan3] + ind3];
}

two_body HF_Channels::tb_state(int &chan1, int &ind1)
{
  return tb_vec[tb_index[chan1] + ind1];
}

void HF_Channels::Build()
{
  int chan1, chan3;
  int ob_total, tb_total;
  this->size1 = SPB.size_2b_dir;
  this->size3 = SPB.size_1b;
  
  this->qnums1 = new State[this->size1];
  this->qnums3 = new State[this->size3];

  p_chan_setup(ob_total,this->nob,this->ob_index,this->ob_map,this->ob_vec, this->size3,this->qnums3, -1);  // ob
  pq_chan_setup(tb_total,this->ntb,this->tb_index,this->tb_map,this->tb_vec, this->size1,this->qnums1, this->size3,this->qnums3, this->nob,this->ob_vec,this->ob_index,this->nob,this->ob_index,this->ob_vec);  // tb

  double memory = 0.0;
  int doubsize = sizeof(double);
  for(chan3 = 0; chan3 < this->size3; ++chan3){
    memory += 3 * doubsize * this->nob[chan3] * this->nob[chan3]; // 3 HF vectors
    memory += 3 * doubsize * this->nob[chan3]; // 3 HF energies
  }
  for(chan1 = 0; chan1 < this->size1; ++chan1){
    memory += doubsize * this->ntb[chan1] * this->ntb[chan1]; // ME.V
  }
  std::cout << std::endl << "Estimated HF-Memory = " << memory/1000000.0 << " MB" << std::endl << std::endl;
}

void HF_Channels::Delete()
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

void Single_Particle_States::Build(HF_Channels &HF_Chan)
{
  int nh, nob, ob_ind, ind, chan;
  int vector_ind, energy_ind;
  int length1, length2;
  this->E3 = 0.0;

  hole_size = new int[HF_Chan.size3];
  vector_size = new int[HF_Chan.size3];
  vector_index = new int[HF_Chan.size3];
  energy_index = new int[HF_Chan.size3];

  length1 = 0;
  length2 = 0;
  for(chan = 0; chan < HF_Chan.size3; ++chan){
    nh = 0;
    nob = HF_Chan.nob[chan];
    for(ind = 0; ind < nob; ++ind){
      ob_ind = HF_Chan.ob_state(chan, ind).v1;
      if( SPB.qnums[ob_ind].type == 0 ){ ++nh; }
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

  for(chan = 0; chan < HF_Chan.size3; ++chan){
    vector_ind = vector_index[chan];
    energy_ind = energy_index[chan];
    for(ind = 0; ind < vector_size[chan]; ++ind){
      ob_ind = HF_Chan.ob_state(chan, ind).v1;
      vectors[vector_ind + (ind * vector_size[chan] + ind)] = 1.0;
      energies[energy_ind + ind] = SPB.qnums[ob_ind].energy;
    }
  }
}

void Single_Particle_States::Delete()
{
  delete[] hole_size;
  delete[] vector_size;
  delete[] vector_index;
  delete[] energy_index;
  delete[] vectors; 
  delete[] energies;
}

void Hartree_Fock_States(HF_Channels &Chan, Single_Particle_States &HF, HF_Matrix_Elements &ME)
{
  double Bshift0 = 0.0;
  for(int i = 0; i < SPB.num_states; ++i){ Bshift0 += SPB.qnums[i].energy; }
  Bshift0 /= (SPB.num_states);

  double tol = 1.0e-12;

  double width = 1.0;
  int Rand_count = 0;

  double error = 1e20, error2 = 1e20;
  int size3 = Chan.size3;
  int ind = 0;
  int ob_ind, energy_ind;
  Single_Particle_States HF0, HF2;
  int badcount = 0;
  int badlimit = 3;

  //// for DIIS /////////
  DIIS DIIS;
  int maxl = 10;
  int size = HF.vector_length;
  double DIIS_start = 0.1;
  DIIS.Build(size, maxl);

  //   Initialize HF0, HF2, and tempn   //
  HF0.Build(Chan);
  HF2.Build(Chan);

  if( PAR.ME3 == 1 ){ Three_Body_Mono(Chan, ME); }

  while(error > tol && ind < 10000){
    Hartree_Fock_Step(Chan, HF, HF2, ME, Bshift0, error);
    if( ind > 1000 && double(Rand_count)/double(ind) > 0.9 ){
      std::cerr << std::endl << ind << " : HF Solution Not Converged!!" << std::endl; exit(1);
    }
    std::cout << "HF-ind = " << ind << ", error = " << error << ", DIIS_count = " << DIIS.count << ", Rand_count = " << Rand_count << std::endl;

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
      if( error < DIIS_start && error > tol ){ DIIS.Perform(HF.vectors, HF0.vectors, 1.0); }
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
      Random_Step_HF(Chan, HF0, HF, HF2, ME, width, error2, Bshift0);
    }
    ++ind;
  }

  if( error > tol ){
    std::cout << PAR.Shells << ", " << PAR.Pshells << ", " << PAR.density << std::endl;
    std::cout << "ind = " << ind << ", error = " << error << ". HF Solution Not Converged!!" << std::endl;
  }

  for(int chan = 0; chan < size3; ++chan){
    energy_ind = HF.energy_index[chan];
    for(ind = 0; ind < HF.vector_size[chan]; ++ind){
      ob_ind = Chan.ob_state(chan, ind).v1;
      SPB.qnums[ob_ind].energy = HF.energies[energy_ind + ind];
    }
  }

  HF0.Delete();
  HF2.Delete();
  DIIS.Delete();


  std::cout << std::setprecision(8);
  if(PAR.basis == "finite_HO"){
    for(int i = 0; i < SPB.num_states; ++i){
      std::cout << SPB.qnums[i].n << " " << SPB.qnums[i].ml << " " << SPB.qnums[i].m << " : " << SPB.qnums[i].energy << " " << SPB.qnums[i].type << std::endl;
    }
  }
  else if(PAR.basis == "finite_J" || PAR.basis == "finite_JM"){
    for(int i = 0; i < SPB.num_states; ++i){
      std::cout << SPB.qnums[i].n << " " << SPB.qnums[i].l << " " << SPB.qnums[i].j << " " << SPB.qnums[i].t << " : " << SPB.qnums[i].energy << " " << SPB.qnums[i].type << std::endl;
    }
  }
  else if(PAR.basis == "finite_M"){
    for(int i = 0; i < SPB.num_states; ++i){
      std::cout << SPB.qnums[i].n << " " << SPB.qnums[i].par << " " << SPB.qnums[i].m << " " << SPB.qnums[i].t << " : " << SPB.qnums[i].energy << " " << SPB.qnums[i].type << std::endl;
    }
  }
}

void Hartree_Fock_Step(HF_Channels &Chan, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &Bshift0, double &error)
{ 
  double *errorvec = new double[Chan.size3];
  double *normvec = new double[Chan.size3];
  for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
    errorvec[chan0] = 0.0;
    normvec[chan0] = 0.0;
  }
  char jobz, uplo; // PAR for Diagonalization, Multiplication
  jobz = 'V';
  uplo = 'U';

  double *Bshift = new double[HF.energy_length];
  for(int i = 0; i < HF.energy_length; ++i){ Bshift[i] = Bshift0; }

  //Make Fock Matrix
  #pragma omp parallel
  {
    int lda; // Parameter for Diagonalization
    int lwork, info; // PAR for Diagonaliztion
    double *w;
    double *work;
    double *fock;
    int size0, length0;
    int vector_ind0, energy_ind0;
    int maxind;
    double maxc;
    double error0, norm0;

    int p, r, pind, qind, mind, rind, sind, nind;
    int key1, key2, chan1, chan3_1, chan3_2;
    int vector_ind1, vector_ind2, ind;
    int length, size1, size2, minj;
    long long count3;
    double term, KE, V, rho1, rho2, HF_0;
    State tb1, tb2;
    #pragma omp for schedule(static)
    for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
      vector_ind0 = HF.vector_index[chan0];
      energy_ind0 = HF.energy_index[chan0];
      size0 = HF.vector_size[chan0];
      fock = new double[size0 * size0];
      length0 = int(0.5 * size0 * (size0 + 1));
      count3 = 0;
      for(int pr = 0; pr < length0; ++pr){
	p = std::floor((2*size0 - 1 - std::sqrt(1 + 4*size0 + 4*size0*size0 - 8*pr))/2) + 1;
	pind = Chan.ob_state(chan0, p).v1;
	length = int(0.5 * p * (2*size0 - p + 1));
	r = int(p + pr - length);
	rind = Chan.ob_state(chan0, r).v1;
	
	HF_0 = 0.0;
	// For regular cases, SPE = diagonal	
	//if(p == r){ KE = SPB.qnums[pind].energy; }
	//else{ KE = 0.0; }

	// For HO-Basis, SPE = not diagonal
	KE = Kinetic_Energy(pind, rind) + Hcom_1_Body(pind, rind, PAR.ho_energy);
	fock[size0 * p + r] = KE;
	
	for(chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
	  vector_ind1 = HF.vector_index[chan3_1];
	  size1 = HF.vector_size[chan3_1];
	  for(int q = 0; q < size1; ++q){
	    qind = Chan.ob_state(chan3_1, q).v1;
	    for(int s = 0; s < size1; ++s){
	      sind = Chan.ob_state(chan3_1, s).v1;
	      plus(tb1, SPB.qnums[pind], SPB.qnums[qind]);
	      minj = std::abs(Chan.qnums3[chan0].j - Chan.qnums3[chan3_1].j);
	      rho1 = 0.0;
	      for(int beta = 0; beta < HF.hole_size[chan3_1]; ++beta){ // Sum over occupied levels
		rho1 += HF.vectors[vector_ind1 + (beta * size1 + q)] * HF.vectors[vector_ind1 + (beta * size1 + s)];
	      }
	      V = 0.0;
	      while( tb1.j >= minj ){
		if( (pind == qind || rind == sind) && tb1.j%4 != 0 ){ tb1.j -= 2; continue; }
		chan1 = Ind_Chan1(tb1);
		key1 = Chan.tb_map[chan1][Hash(pind, qind, tb1.j)];
		key2 = Chan.tb_map[chan1][Hash(rind, sind, tb1.j)];
		ind =  ME.Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		V += ME.V[ind] * (tb1.j + 1.0)/(Chan.qnums3[chan0].j + 1.0);
		tb1.j -= 2;
	      }
	      term = rho1 * V;
	      /*if(rho1 == 1.0){
		std::cout << "< " << pind << ", " << qind << " || " << rind << ", " << sind << " > = " << V << " * " << rho1 << std::endl;
		}*/
	      HF_0 += term;
	      fock[size0 * p + r] += term;
	    }
	  }
	}
	
	///////////////  3N  //////////////////  <pqm||rsn>
	if( PAR.ME3 == 1 && SPB.e_nlj[SPB.qnums[pind].nlj] <= SPB.E3Max && SPB.e_nlj[SPB.qnums[rind].nlj] <= SPB.E3Max ){
	  for(chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
	    vector_ind1 = HF.vector_index[chan3_1];
	    size1 = HF.vector_size[chan3_1];

	    for(int q = 0; q < size1; ++q){
	      qind = Chan.ob_state(chan3_1, q).v1;
	      if( SPB.e_nlj[SPB.qnums[pind].nlj] + SPB.e_nlj[SPB.qnums[qind].nlj] > SPB.E3Max ){ continue; }
	      for(int s = 0; s < size1; ++s){
		sind = Chan.ob_state(chan3_1, s).v1;
		if( SPB.e_nlj[SPB.qnums[rind].nlj] + SPB.e_nlj[SPB.qnums[sind].nlj] > SPB.E3Max ){ continue; }

		rho1 = 0.0;
		for(int beta1 = 0; beta1 < HF.hole_size[chan3_1]; ++beta1){ // Sum over occupied levels
		  rho1 += HF.vectors[vector_ind1 + (beta1 * size1 + q)] * HF.vectors[vector_ind1 + (beta1 * size1 + s)];
		}
		
		for(chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){
		  vector_ind2 = HF.vector_index[chan3_2];
		  size2 = HF.vector_size[chan3_2];
		  
		  for(int m = 0; m < size2; ++m){
		    mind = Chan.ob_state(chan3_2, m).v1;
		    if( SPB.e_nlj[SPB.qnums[pind].nlj] + SPB.e_nlj[SPB.qnums[qind].nlj] + SPB.e_nlj[SPB.qnums[mind].nlj] > SPB.E3Max ){ continue; }
		    for(int n = 0; n < size2; ++n){
		      nind = Chan.ob_state(chan3_2, n).v1;
		      if( SPB.e_nlj[SPB.qnums[rind].nlj] + SPB.e_nlj[SPB.qnums[sind].nlj] + SPB.e_nlj[SPB.qnums[nind].nlj] > SPB.E3Max ){ continue; }

		      rho2 = 0.0;
		      for(int beta2 = 0; beta2 < HF.hole_size[chan3_2]; ++beta2){ // Sum over occupied levels
			rho2 += HF.vectors[vector_ind2 + (beta2 * size2 + m)] * HF.vectors[vector_ind2 + (beta2 * size2 + n)];
		      }
		      if( rho1 * rho2 == 0.0 ){
			++count3;
			continue;
		      }

		      V = ME.V3mon[ME.V3mon_index[chan0] + count3];
		      //std::cout << "< " << pind << "," << qind << "," << mind << " || " << rind << "," << sind << "," << nind << " > = " << V << " * " << rho1*rho2 << " : " << count3 << std::endl;
		      ++count3;

		      term = 0.5 * rho1 * rho2 * V;
		      HF_0 += term;
		      fock[size0 * p + r] += term;
		    }
		  }
		}
	      }
	    }
	  }
	}
	
	for(int beta0 = HF.hole_size[chan0]; beta0 < size0; ++beta0){ // Sum over unoccupied levels
	  fock[size0 * p + r] += HF.vectors[vector_ind0 + (beta0 * size0 + p)] * HF.vectors[vector_ind0 + (beta0 * size0 + r)] * Bshift[energy_ind0 + beta0];
	}
	if(p != r){ fock[size0 * r + p] = fock[size0 * p + r]; }
      }
      
      /*std::cout << "!!!   Fock   !!!" << std::endl;
      for(int j = 0; j < size0; ++j){
	for(int k = 0; k < size0; ++k){
	  std::cout << fock[size0 * j + k] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
	
      lda = size0;
      lwork = (3+2)*size0;
      w = new double[lda];
      work = new double[lwork];
      for(int j = 0; j < lda; ++j){ w[j] = 0.0; }
      for(int j = 0; j < lwork; ++j){ work[j] = 0.0; }
      if(size0 != 0){ dsyev_(&jobz, &uplo, &size0, fock, &lda, w, work, &lwork, &info); }
      
      double deltaE;
      for(int j = HF.hole_size[chan0]; j < size0; ++j){ w[j] -= Bshift[energy_ind0 + j]; } // Add back Level-shift parameter
      for(int j = HF.hole_size[chan0]; j < size0; ++j){
	deltaE = fabs(w[j] - HF.energies[energy_ind0 + j]);
	if( deltaE > Bshift0 || HF.hole_size[chan0] == 0 || HF.vector_size[chan0] == HF.hole_size[chan0]){ Bshift[energy_ind0 + j] = Bshift0 + deltaE; }
	else{ Bshift[energy_ind0 + j] = Bshift0 + w[HF.hole_size[chan0]] - w[HF.hole_size[chan0] - 1]; }
      }
      
      for(int j = 0; j < size0; ++j){
	HF2.energies[energy_ind0 + j] = w[j];
	maxind = -1;
	maxc = 0.0;
	for(int k = 0; k < size0; ++k){
	  if( fabs(fock[size0 * j + k]) > maxc ){
	    maxc = fabs(fock[size0 * j + k]);
	    maxind = k;
	  }
	}
	if( HF.vectors[vector_ind0 + (size0*j + maxind)]/fock[size0 * j + maxind] < 0.0 ){
	  for(int k = 0; k < size0; ++k){ fock[size0 * j + k] *= -1.0; }
	}
      }

      error0 = 0.0;
      norm0 = 0.0;
      for(int j = 0; j < size0; ++j){
	for(int k = 0; k < size0; ++k){
	  error0 += std::pow(fabs(HF.vectors[vector_ind0 + (size0 * j + k)]) - fabs(fock[size0 * j + k]), 2.0);
	  norm0 += std::pow(fock[size0 * j + k], 2.0);
	  //temperror += std::pow(fabs(HF.vectors[vector_ind0 + (size0 * j + k)]) - fabs(fock[size0 * j + k]), 2.0);
	  //tempnorm += std::pow(fock[size0 * j + k], 2.0);
	  HF2.vectors[vector_ind0 + (size0 * j + k)] = fock[size0 * j + k];
	}
      }
      errorvec[chan0] = error0;
      normvec[chan0] = norm0;
      delete[] fock;
      delete[] w;
      delete[] work;
    }
  }
  delete[] Bshift;

  double temperror = 0.0;
  double tempnorm = 0.0;
  for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
    temperror += errorvec[chan0];
    tempnorm += normvec[chan0];
  }
  error = std::sqrt(temperror/tempnorm);
}


void Three_Body_Mono(HF_Channels &Chan, HF_Matrix_Elements &ME)
{ 
  int size0, length0, length, size1, size2;
  int p, r, pind, qind, mind, rind, sind, nind;
  int p_e, q_e, m_e, r_e, s_e, n_e;
  State tb1, tb2;
  long long count0;
  ME.V3mon_count = 0;
  ME.V3mon_index = new long long[Chan.size3];

  // count Vmon
  for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
    size0 = Chan.nob[chan0];
    length0 = int(0.5 * size0 * (size0 + 1));
    count0 = 0;

    for(int pr = 0; pr < length0; ++pr){
      p = std::floor((2*size0 - 1 - std::sqrt(1 + 4*size0 + 4*size0*size0 - 8*pr))/2) + 1;
      pind = Chan.ob_state(chan0, p).v1;
      p_e = SPB.e_nlj[SPB.qnums[pind].nlj];
      length = int(0.5 * p * (2*size0 - p + 1));
      r = int(p + pr - length);
      rind = Chan.ob_state(chan0, r).v1;
      r_e = SPB.e_nlj[SPB.qnums[rind].nlj];
      if( p_e > SPB.E3Max || r_e > SPB.E3Max ){ continue; }

      for(int chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
	size1 = Chan.nob[chan3_1];	  
	for(int q = 0; q < size1; ++q){
	  qind = Chan.ob_state(chan3_1, q).v1;
	  q_e = SPB.e_nlj[SPB.qnums[qind].nlj];
	  if( p_e + q_e > SPB.E3Max ){ continue; }
	  for(int s = 0; s < size1; ++s){
	    sind = Chan.ob_state(chan3_1, s).v1;
	    s_e = SPB.e_nlj[SPB.qnums[sind].nlj];
	    if( r_e + s_e > SPB.E3Max ){ continue; }
	    
	    for(int chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){
	      size2 = Chan.nob[chan3_2];	
	      for(int m = 0; m < size2; ++m){
		mind = Chan.ob_state(chan3_2, m).v1;
		m_e = SPB.e_nlj[SPB.qnums[mind].nlj];
		if( p_e + q_e + m_e > SPB.E3Max ){ continue; }
		for(int n = 0; n < size2; ++n){
		  nind = Chan.ob_state(chan3_2, n).v1;
		  n_e = SPB.e_nlj[SPB.qnums[nind].nlj];
		  if( r_e + s_e + n_e > SPB.E3Max ){ continue; }

		  ++count0;
		}
	      }
	    }
	  }
	}
      }
    }
    ME.V3mon_index[chan0] = ME.V3mon_count;
    ME.V3mon_count += count0;
  }
  std::cout << "3-body Mono Length = " << ME.V3mon_count << std::endl;
  ME.V3mon = new float[ME.V3mon_count];
  int *inds = new int[6 * ME.V3mon_count];
  int *js = new int[4 * ME.V3mon_count];

  // get Vmon indicies and j values
  ME.V3mon_count = 0;
  for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
    size0 = Chan.nob[chan0];
    length0 = int(0.5 * size0 * (size0 + 1));

    for(int pr = 0; pr < length0; ++pr){
      p = std::floor((2*size0 - 1 - std::sqrt(1 + 4*size0 + 4*size0*size0 - 8*pr))/2) + 1;
      pind = Chan.ob_state(chan0, p).v1;
      p_e = SPB.e_nlj[SPB.qnums[pind].nlj];
      length = int(0.5 * p * (2*size0 - p + 1));
      r = int(p + pr - length);
      rind = Chan.ob_state(chan0, r).v1;
      r_e = SPB.e_nlj[SPB.qnums[rind].nlj];
      if( p_e > SPB.E3Max || r_e > SPB.E3Max ){ continue; }

      for(int chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
	size1 = Chan.nob[chan3_1];	  
	for(int q = 0; q < size1; ++q){
	  qind = Chan.ob_state(chan3_1, q).v1;
	  q_e = SPB.e_nlj[SPB.qnums[qind].nlj];
	  if( p_e + q_e > SPB.E3Max ){ continue; }
	  for(int s = 0; s < size1; ++s){
	    sind = Chan.ob_state(chan3_1, s).v1;
	    s_e = SPB.e_nlj[SPB.qnums[sind].nlj];
	    if( r_e + s_e > SPB.E3Max ){ continue; }
	    
	    for(int chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){
	      size2 = Chan.nob[chan3_2];	
	      for(int m = 0; m < size2; ++m){
		mind = Chan.ob_state(chan3_2, m).v1;
		m_e = SPB.e_nlj[SPB.qnums[mind].nlj];
		if( p_e + q_e + m_e > SPB.E3Max ){ continue; }
		for(int n = 0; n < size2; ++n){
		  nind = Chan.ob_state(chan3_2, n).v1;
		  n_e = SPB.e_nlj[SPB.qnums[nind].nlj];
		  if( r_e + s_e + n_e > SPB.E3Max ){ continue; }

		  inds[6*ME.V3mon_count] = pind;
		  inds[6*ME.V3mon_count + 1] = qind;
		  inds[6*ME.V3mon_count + 2] = mind;
		  inds[6*ME.V3mon_count + 3] = rind;
		  inds[6*ME.V3mon_count + 4] = sind;
		  inds[6*ME.V3mon_count + 5] = nind;

		  js[4*ME.V3mon_count] = std::abs(Chan.qnums3[chan0].j - Chan.qnums3[chan3_1].j);
		  js[4*ME.V3mon_count + 1] = Chan.qnums3[chan0].j + Chan.qnums3[chan3_1].j;
		  js[4*ME.V3mon_count + 2] = Chan.qnums3[chan3_2].j;
		  js[4*ME.V3mon_count + 3] = Chan.qnums3[chan0].j;

		  ++ME.V3mon_count;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  // compute Vmon
  #pragma omp parallel for schedule(static,50)
  for(int i = 0; i < ME.V3mon_count; ++i){
    ME.V3mon[i] = 0.0;
    for(int j12 = js[4*i]; j12 <= js[4*i + 1]; j12 += 2){
      for(int J = std::abs(j12 - js[4*i + 2]); J <= (j12 + js[4*i + 2]); J += 2){
	ME.V3mon[i] += ME3J_GetME_pn(ME,inds[6*i],inds[6*i+1],inds[6*i+2],inds[6*i+3],inds[6*i+4],inds[6*i+5], j12, j12, J) * (J + 1.0)/(js[4*i + 3] + 1.0);
      }
    }
  }
  delete[] inds;
  delete[] js;
}



void Random_Step_HF(HF_Channels &Chan, Single_Particle_States &HF0, Single_Particle_States &HF, Single_Particle_States &HF2, HF_Matrix_Elements &ME, double &width, double &error2, double &Bshift)
{
  double error;
  int go = 0;
  int ind = 0;
  while(go == 0){
    ++ind;
    Randomize_HF(Chan, HF0, HF, width);
    Hartree_Fock_Step(Chan, HF, HF2, ME, Bshift, error);
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

void Convert_To_HF_Matrix_Elements(HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  int size, length1, matlength; // max length of M-Scheme indicies, length of J_ME
  double *M1; // Matrices of coefficients
  double *C;
  //double *V3;
  char transa, transb;
  double alpha1, beta1;

  for(int chan = 0; chan < Chan.size1; ++chan){
    size = Chan.ntb[chan];
    if(size == 0){ continue; }
    matlength = size * size;
    length1 = int(0.5 * size * (size + 1));
    M1 = new double[matlength];
    C = new double[matlength];
    //V3 = new double[matlength];
    for(int i = 0; i < matlength; ++i){
      M1[i] = 0.0;
      C[i] = 0.0;
      //V3[i] = 0.0;
    }

    #pragma omp parallel
    {
      int pq, rs, length2;
      two_body pq_tb, rs_tb;
      int p, q, r, s;
      int p_chan, q_chan, r_chan, s_chan;
      int p_ind, q_ind, r_ind, s_ind;
      int vec_ind_r, vec_ind_s, p_size, q_size;
      int vec_ind_B, size_B, i_ind, j_ind, jmin;
      double term, rho;
      State thb;
      #pragma omp for schedule(static)
      for(int pqrs = 0; pqrs < length1; ++pqrs){
	pq = std::floor((2*size - 1 - std::sqrt(1 + 4*size + 4*matlength - 8*pqrs))/2) + 1;
	length2 = int(0.5 * pq * (2*size - pq + 1));
	rs = int(pq + pqrs - length2);
	pq_tb = Chan.tb_state(chan, pq);
	rs_tb = Chan.tb_state(chan, rs);
	p_ind = pq_tb.v1;
	q_ind = pq_tb.v2;
	r_ind = rs_tb.v1;
	s_ind = rs_tb.v2;
	p_chan = Ind_Chan3(SPB.qnums[p_ind]);
	q_chan = Ind_Chan3(SPB.qnums[q_ind]);
	r_chan = Ind_Chan3(SPB.qnums[r_ind]);
	s_chan = Ind_Chan3(SPB.qnums[s_ind]);
	p = Chan.ob_map[p_chan][p_ind];
	q = Chan.ob_map[q_chan][q_ind];
	r = Chan.ob_map[r_chan][r_ind];
	s = Chan.ob_map[s_chan][s_ind];
	vec_ind_r = States.vector_index[r_chan];
	vec_ind_s = States.vector_index[s_chan];
	p_size = States.vector_size[p_chan];
	q_size = States.vector_size[q_chan];
	if(p_chan == r_chan && q_chan == s_chan){
	  M1[size * pq + rs] = States.vectors[vec_ind_r + (r * p_size + p)] * States.vectors[vec_ind_s + (s * q_size + q)];
	  if( pq != rs ){ M1[size * rs + pq] = States.vectors[vec_ind_r + (p * p_size + r)] * States.vectors[vec_ind_s + (q * q_size + s)]; }
	}

	////////////////  For 3N-ME, V^{pq}_{rs} <- [(J3^/J^)^2] \sum_B W^{pq^J i}_{rs^J j} * C^{i}_{B} * C^{j}_{B}
	if( PAR.ME3 == 1){
	  for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
	    vec_ind_B = States.vector_index[chan3];
	    size_B = States.vector_size[chan3];
	    for(int i = 0; i < size_B; ++i){
	      i_ind = Chan.ob_state(chan3, i).v1;
	      if( SPB.e_nlj[SPB.qnums[p_ind].nlj] + SPB.e_nlj[SPB.qnums[q_ind].nlj] + SPB.e_nlj[SPB.qnums[i_ind].nlj] > SPB.E3Max ){ continue; }
	      for(int j = 0; j < size_B; ++j){
		j_ind = Chan.ob_state(chan3, j).v1;
		if( SPB.e_nlj[SPB.qnums[r_ind].nlj] + SPB.e_nlj[SPB.qnums[s_ind].nlj] + SPB.e_nlj[SPB.qnums[j_ind].nlj] > SPB.E3Max ){ continue; }
		
		rho = 0.0;
		for(int beta = 0; beta < States.hole_size[chan3]; ++beta){ // Sum over occupied levels
		  rho += States.vectors[vec_ind_B + (beta * size_B + i)] * States.vectors[vec_ind_B + (beta * size_B + j)];
		}
		
		plus(thb, Chan.qnums1[chan], Chan.qnums3[chan3]);
		jmin = std::abs(Chan.qnums1[chan].j - Chan.qnums3[chan3].j);
		while(thb.j >= jmin){
		  term = rho * ME3J_GetME_pn(ME, p_ind, q_ind, i_ind, r_ind, s_ind, j_ind, Chan.qnums1[chan].j, Chan.qnums1[chan].j, thb.j);
		  term *= (thb.j + 1.0)/(Chan.qnums1[chan].j + 1.0);
		  //std::cout << "^^ < " << p_ind << "," << q_ind << "," << i_ind << " | " << r_ind << "," << s_ind << "," << j_ind << " > = " << ME3J_GetME_pn(ME, p_ind, q_ind, i_ind, r_ind, s_ind, j_ind, Chan.qnums1[chan].j, Chan.qnums1[chan].j, thb.j) << "  " << rho << "  " << (thb.j + 1.0)/(Chan.qnums1[chan].j + 1.0) << std::endl;
		  ME.V[ME.Index[chan] + (size * pq + rs)] += term;
		  if( pq != rs ){ ME.V[ME.Index[chan] + (size * rs + pq)] += term; }
		  thb.j -= 2;
		}
	      }
	    }
	  }
	}
	///////////////////////////////////////////////////////////////////////////////////////

      }
    }

    transa = 'N';
    transb = 'T';
    alpha1 = 1.0;
    beta1 = 0.0;
    dgemm_NN(ME.V + ME.Index[chan], M1, C, &size, &size, &size, &alpha1, &beta1, &transa, &transa);
    dgemm_TN(M1, C, ME.V + ME.Index[chan], &size, &size, &size, &alpha1, &beta1, &transb, &transa);

    //delete[] V3;
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
	//interactionfile << std::setw(4) << SPB.qnums[p].n << std::setw(4) << SPB.qnums[p].ml << std::setw(4) << SPB.qnums[p].m;
	//interactionfile << std::setw(4) << SPB.qnums[q].n << std::setw(4) << SPB.qnums[q].ml << std::setw(4) << SPB.qnums[q].m;
	//interactionfile << std::setw(4) << SPB.qnums[r].n << std::setw(4) << SPB.qnums[r].ml << std::setw(4) << SPB.qnums[r].m;
	//interactionfile << std::setw(4) << SPB.qnums[s].n << std::setw(4) << SPB.qnums[s].ml << std::setw(4) << SPB.qnums[s].m;
	//interactionfile << std::setw(24) << std::setprecision(16) << tbme << "\n";
      }
    }
  }
  interactionfile.close();*/
}


double Eref3(HF_Channels &Chan, Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  int length, vector_ind0, size0, length0, vector_ind1, size1, vector_ind2, size2;
  int i, j, iind, pind, qind, jind, rind, sind;
  long long count3;
  double rho0, rho1, rho2;
  double Eref = 0.0;
  for(int chan0 = 0; chan0 < Chan.size3; ++chan0){
    vector_ind0 = States.vector_index[chan0];
    size0 = States.vector_size[chan0];
    length0 = int(0.5 * size0 * (size0 + 1));
    count3 = 0;
    
    for(int ij = 0; ij < length0; ++ij){
      i = std::floor((2*size0 - 1 - std::sqrt(1 + 4*size0 + 4*size0*size0 - 8*ij))/2) + 1;
      iind = Chan.ob_state(chan0, i).v1;
      length = int(0.5 * i * (2*size0 - i + 1));
      j = int(i + ij - length);
      jind = Chan.ob_state(chan0, j).v1;
      if( SPB.e_nlj[SPB.qnums[iind].nlj] > SPB.E3Max || SPB.e_nlj[SPB.qnums[jind].nlj] > SPB.E3Max ){ continue; }
      
      rho0 = 0.0;
      for(int beta0 = 0; beta0 < States.hole_size[chan0]; ++beta0){ // Sum over occupied levels
	rho0 += States.vectors[vector_ind0 + (beta0 * size0 + i)] * States.vectors[vector_ind0 + (beta0 * size0 + j)];
      }
      
      for(int chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
	vector_ind1 = States.vector_index[chan3_1];
	size1 = States.vector_size[chan3_1];
	for(int p = 0; p < size1; ++p){
	  pind = Chan.ob_state(chan3_1, p).v1;
	  if( SPB.e_nlj[SPB.qnums[iind].nlj] + SPB.e_nlj[SPB.qnums[pind].nlj] > SPB.E3Max ){ continue; }
	  for(int r = 0; r < size1; ++r){
	    rind = Chan.ob_state(chan3_1, r).v1;
	    if( SPB.e_nlj[SPB.qnums[jind].nlj] + SPB.e_nlj[SPB.qnums[rind].nlj] > SPB.E3Max ){ continue; }
	    
	    rho1 = 0.0;
	    for(int beta1 = 0; beta1 < States.hole_size[chan3_1]; ++beta1){ // Sum over occupied levels
	      rho1 += States.vectors[vector_ind1 + (beta1 * size1 + p)] * States.vectors[vector_ind1 + (beta1 * size1 + r)];
	    }
	      
	    for(int chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){
	      vector_ind2 = States.vector_index[chan3_2];
	      size2 = States.vector_size[chan3_2];
	      for(int q = 0; q < size2; ++q){
		qind = Chan.ob_state(chan3_2, q).v1;
		if( SPB.e_nlj[SPB.qnums[iind].nlj] + SPB.e_nlj[SPB.qnums[pind].nlj] + SPB.e_nlj[SPB.qnums[qind].nlj] > SPB.E3Max ){ continue; }
		for(int s = 0; s < size2; ++s){
		  sind = Chan.ob_state(chan3_2, s).v1;
		  if( SPB.e_nlj[SPB.qnums[jind].nlj] + SPB.e_nlj[SPB.qnums[rind].nlj] + SPB.e_nlj[SPB.qnums[sind].nlj] > SPB.E3Max ){ continue; }
		  
		  rho2 = 0.0;
		  for(int beta2 = 0; beta2 < States.hole_size[chan3_2]; ++beta2){ // Sum over occupied levels
		    rho2 += States.vectors[vector_ind2 + (beta2 * size2 + q)] * States.vectors[vector_ind2 + (beta2 * size2 + s)];
		  }
		  
		  Eref += (1.0/6.0) * (Chan.qnums3[chan0].j + 1.0) * rho0 * rho1 * rho2 * ME.V3mon[ME.V3mon_index[chan0] + count3];
		  if( i != j ){ Eref += (1.0/6.0) * (Chan.qnums3[chan0].j + 1.0) * rho0 * rho1 * rho2 * ME.V3mon[ME.V3mon_index[chan0] + count3]; }
		  ++count3;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return Eref;
}

void HF_Matrix_Elements::Read_J(HF_Channels &Chan)
{
  std::string fullpath1 = PATH + PAR.MatrixElements + ".int"; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME, hom, r2, p2; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, coupJ, coupT, par; // interaction file contents
  int chan1, ind, key1, key2, key3, key4;
  State tb;

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << PAR.MatrixElements << ", does not exist" << std::endl; exit(1); }
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
    //if(shell1 > 44 || shell2 > 44 || shell3 > 44 || shell4 > 44){ continue; }
    //TBME *= PAR.tbstrength;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    TBME *= d_ij(shell1, shell2); // unnormalize
    TBME *= d_ij(shell3, shell4);

    //std::cout << coupT << " " << par << " " << coupJ << " " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
    plus(tb, SPB.qnums[shell1], SPB.qnums[shell2]);
    tb.j = coupJ;
    chan1 = Ind_Chan1(tb);
    key1 = Chan.tb_map[chan1][Hash(shell1, shell2, tb.j)];
    key2 = Chan.tb_map[chan1][Hash(shell3, shell4, tb.j)];
    key3 = Chan.tb_map[chan1][Hash(shell2, shell1, tb.j)];
    key4 = Chan.tb_map[chan1][Hash(shell4, shell3, tb.j)];
    ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
    this->V[ind] = TBME;
    ind = this->Index[chan1] + (key3 * Chan.ntb[chan1] + key2);
    this->V[ind] = -1.0 * phase2(SPB.qnums[shell1].j + SPB.qnums[shell2].j - coupJ) * TBME;
    ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key4);
    this->V[ind] = -1.0 * phase2(SPB.qnums[shell3].j + SPB.qnums[shell4].j - coupJ) * TBME;
    ind = this->Index[chan1] + (key3 * Chan.ntb[chan1] + key4);
    this->V[ind] = phase2(SPB.qnums[shell1].j + SPB.qnums[shell2].j + SPB.qnums[shell3].j + SPB.qnums[shell4].j) * TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = this->Index[chan1] + (key2 * Chan.ntb[chan1] + key1);
      this->V[ind] = TBME;
      ind = this->Index[chan1] + (key4 * Chan.ntb[chan1] + key1);
      this->V[ind] = -1.0 * phase2(SPB.qnums[shell3].j + SPB.qnums[shell4].j - coupJ) * TBME;
      ind = this->Index[chan1] + (key2 * Chan.ntb[chan1] + key3);
      this->V[ind] = -1.0 * phase2(SPB.qnums[shell1].j + SPB.qnums[shell2].j - coupJ) * TBME;
      ind = this->Index[chan1] + (key4 * Chan.ntb[chan1] + key3);
      this->V[ind] = phase2(SPB.qnums[shell1].j + SPB.qnums[shell2].j + SPB.qnums[shell3].j + SPB.qnums[shell4].j) * TBME;
    }
  }
  interaction.close();
}

void HF_Matrix_Elements::Read_M(HF_Channels &Chan)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int chan1, ind, key1, key2, key3, key4;
  State tb;

  fullpath1 = PATH + PAR.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << PAR.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;
  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
    TBME *= PAR.tbstrength;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;

    //std::cout << "? " << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << "  " << TBME << std::endl;
    if(shell1 == shell2 || shell3 == shell4){ continue; }
    plus(tb, SPB.qnums[shell1], SPB.qnums[shell2]);
    chan1 = Ind_Chan1(tb);
    key1 = Chan.tb_map[chan1][Hash(shell1, shell2, tb.j)];
    key2 = Chan.tb_map[chan1][Hash(shell3, shell4, tb.j)];
    key3 = Chan.tb_map[chan1][Hash(shell2, shell1, tb.j)];
    key4 = Chan.tb_map[chan1][Hash(shell4, shell3, tb.j)];
    ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
    this->V[ind] = TBME;
    ind = this->Index[chan1] + (key3 * Chan.ntb[chan1] + key2);
    this->V[ind] = -1.0 * TBME;
    ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key4);
    this->V[ind] = -1.0 * TBME;
    ind = this->Index[chan1] + (key3 * Chan.ntb[chan1] + key4);
    this->V[ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = this->Index[chan1] + (key2 * Chan.ntb[chan1] + key1);
      this->V[ind] = TBME;
      ind = this->Index[chan1] + (key4 * Chan.ntb[chan1] + key1);
      this->V[ind] = -1.0 * TBME;
      ind = this->Index[chan1] + (key2 * Chan.ntb[chan1] + key3);
      this->V[ind] = -1.0 * TBME;
      ind = this->Index[chan1] + (key4 * Chan.ntb[chan1] + key3);
      this->V[ind] = TBME;
    }
  }
  interaction.close();
}

void HF_Matrix_Elements::Read_QD(HF_Channels &Chan)
{
  std::cout << "Computing Coulomb Matrix Elements for QD..." << std::endl;
  struct timespec time1, time2;
  double elapsed0 = 0.0;
  clock_gettime(CLOCK_MONOTONIC, &time1);
  
  int pq, rs, p, q, r, s;
  int ntb, ntb0;
  int ind1, ind2, ind3, ind4;
  int length, length0;
  two_body *tbvec0;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    ntb = Chan.ntb[chan];
    if(ntb == 0){ continue; }
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tb_vec[Chan.tb_index[chan] + i].v1 < Chan.tb_vec[Chan.tb_index[chan] + i].v2){
	++ntb0;
      }
    }
    tbvec0 = new two_body[ntb0];
    ntb0 = 0;
    for(int i = 0; i < ntb; ++i){
      if(Chan.tb_vec[Chan.tb_index[chan] + i].v1 < Chan.tb_vec[Chan.tb_index[chan] + i].v2){
	tbvec0[ntb0].v1 = Chan.tb_vec[Chan.tb_index[chan] + i].v1;
	tbvec0[ntb0].v2 = Chan.tb_vec[Chan.tb_index[chan] + i].v2;
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
	p = tbvec0[pq].v1;
	q = tbvec0[pq].v2;
	r = tbvec0[rs].v1;
	s = tbvec0[rs].v2;
	ind1 = Chan.tb_map[chan][Hash(p, q, Chan.qnums1[chan].j)];
	ind2 = Chan.tb_map[chan][Hash(r, s, Chan.qnums1[chan].j)];
	ind3 = Chan.tb_map[chan][Hash(q, p, Chan.qnums1[chan].j)];
	ind4 = Chan.tb_map[chan][Hash(s, r, Chan.qnums1[chan].j)];
	TBME = Coulomb_HO(p, q, r, s);
	this->V[this->Index[chan] + (ind1 * ntb + ind2)] = TBME;
	this->V[this->Index[chan] + (ind3 * ntb + ind2)] = -1.0 * TBME;
	this->V[this->Index[chan] + (ind1 * ntb + ind4)] = -1.0 * TBME;
	this->V[this->Index[chan] + (ind3 * ntb + ind4)] = TBME;
	if(ind1 != ind2){
	  this->V[this->Index[chan] + (ind2 * ntb + ind1)] = TBME;
	  this->V[this->Index[chan] + (ind4 * ntb + ind1)] = -1.0 * TBME;
	  this->V[this->Index[chan] + (ind2 * ntb + ind3)] = -1.0 * TBME;
	  this->V[this->Index[chan] + (ind4 * ntb + ind3)] = TBME;
	}
      }
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &time2);
  elapsed0 = (time2.tv_sec - time1.tv_sec);
  elapsed0 += (time2.tv_nsec - time1.tv_nsec) / 1000000000.0;
  std::cout << std::endl << "!! Runtime = " << elapsed0 << " sec. " << std::endl;
}

void HF_Matrix_Elements::Read_QDFile(HF_Channels &Chan)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  std::ifstream interaction;	// interaction file

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
    length /= 16;

    unsigned int neg1 = 128;
    unsigned int neg2 = 4294967040;
    int nmax = PAR.Shells;
    #pragma omp parallel
    {
      double TBME;
      unsigned int n1, n2, n3, n4;
      int ml1, ml2, ml3, ml4;
      int key, key1, key2;
      State tb;
      int ind, chan1;
      State statep, stateq, stater, states;
      int p, q, r, s;
      statep.t = -1, stateq.t = -1, stater.t = -1, states.t = -1;

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
	TBME *= std::sqrt(PAR.density);
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
	  key = Ind_State(statep);
	  p = SPB.map_state[key];
	  for(int s2 = -1; s2 <= 1; s2 += 2){
	    stateq.m = s2;
	    key = Ind_State(stateq);
	    q = SPB.map_state[key];
	    if(p == q){ continue; }
	    for(int s3 = -1; s3 <= 1; s3 += 2){
	      stater.m = s3;
	      key = Ind_State(stater);
	      r = SPB.map_state[key];
	      if(s3 != s1){ continue; }
	      for(int s4 = -1; s4 <= 1; s4 += 2){
		states.m = s4;
		key = Ind_State(states);
		s = SPB.map_state[key];
		if(r == s || s4 != s2){ continue; }
		
		plus(tb, SPB.qnums[p], SPB.qnums[q]);
		chan1 = Ind_Chan1(tb);

		// C(p1q2r3s4) -> <p1q2 || r3s4>
		key1 = Chan.tb_map[chan1][Hash(p, q, 0)];
		key2 = Chan.tb_map[chan1][Hash(r, s, 0)];
		ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		this->V[ind] += TBME;
		// C(p1q2r3s4) -> -<p1q2 || s4r3>
		key2 = Chan.tb_map[chan1][Hash(s, r, 0)];
		ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		this->V[ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(q2p1s4r3) -> <q2p1 || s4r3>
		  key1 = Chan.tb_map[chan1][Hash(q, p, 0)];
		  key2 = Chan.tb_map[chan1][Hash(s, r, 0)];
		  ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  this->V[ind] += TBME;
		  // C(p1q2r3s4) = C(q2p1s4r3) -> -<q2p1 || r3s4>
		  key2 = Chan.tb_map[chan1][Hash(r, s, 0)];
		  ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  this->V[ind] -= TBME;
		}
		if(((n1 == n3 && ml1 == ml3) && (n2 == n4 && ml2 == ml4)) || ((n1 == n4 && ml1 == ml4) && (n2 == n3 && ml2 == ml3))){ continue; }
		// C(p1q2r3s4) = C(r3s4p1q2) -> <r3s4 || p1q2>
		key1 = Chan.tb_map[chan1][Hash(r, s, 0)];
		key2 = Chan.tb_map[chan1][Hash(p, q, 0)];
		ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		this->V[ind] += TBME;
		// C(p1q2r3s4) = C(r3s4p1q2) -> -<r3s4 || q2p1>
		key2 = Chan.tb_map[chan1][Hash(q, p, 0)];
		ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		this->V[ind] -= TBME;
		if((n1 != n2 || ml1 != ml2) || (n3 != n4 || ml3 != ml4)){
		  // C(p1q2r3s4) = C(s4r3q2p1) -> <s4r3 || q2p1>
		  key1 = Chan.tb_map[chan1][Hash(s, r, 0)];
		  key2 = Chan.tb_map[chan1][Hash(q, p, 0)];
		  ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  this->V[ind] += TBME;
		  // C(p1q2r3s4) = C(s4r3q2p1) -> -<s4r3 || p1q2>
		  key2 = Chan.tb_map[chan1][Hash(p, q, 0)];
		  ind = this->Index[chan1] + (key1 * Chan.ntb[chan1] + key2);
		  this->V[ind] -= TBME;
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
