#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"

HF_Matrix_Elements::HF_Matrix_Elements(const HF_Channels &Chan)
{
  V.resize(Chan.size1);
  for(int i = 0; i < Chan.size1; ++i){
    V[i].assign(Chan.tb[i] * Chan.tb[i], 0.0);
  }
}

HF_Matrix_Elements::HF_Matrix_Elements()
{
  V.resize(0);
}

//Function to setup Channels
void Setup_Channels_HF(const Input_Parameters &Parameters, const Model_Space &Space, HF_Channels &Chan)
{
  std::cout << "Building HF Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  std::vector<int> nums(3); // subvector for OB_nums. size 3 for finite_J
  std::vector<std::vector<int> > OBvec;

  State state;

  // place first state by quantum numbers
  Chan.size3 = 1;
  OBvec.resize(Chan.size3);
  if(Parameters.basis == "finite_M"){
    state.t = Space.qnums[0].t; //tz
    state.m = Space.qnums[0].m; //jz
    state.par = Space.qnums[0].par; //par
  }
  else if(Parameters.basis == "finite_HO"){
    state.t = Space.qnums[0].t; //tz
    state.m = Space.qnums[0].m; //sz
    state.ml = Space.qnums[0].ml; //lz
  }
  else if(Parameters.basis == "finite_J"){
    state.t = Space.qnums[0].t; //tz
    state.j = Space.qnums[0].j; //j
    state.par = Space.qnums[0].par; //par
  }
  Chan.qnums3.push_back(state);
  Chan.indvec.push_back(Chan.size3 - 1);
  OBvec[Chan.size3 - 1].push_back(0);

  std::cout << "test0" << std::endl;

  // place the rest of the states by their quantum numbers
  for(int i = 1; i < Space.indtot; ++i){
    if(Parameters.basis == "finite_M"){
      state.t = Space.qnums[i].t; //tz
      state.m = Space.qnums[i].m; //jz
      state.par = Space.qnums[i].par; //par
    }
    else if(Parameters.basis == "finite_HO"){
      state.t = Space.qnums[i].t; //tz
      state.m = Space.qnums[i].m; //sz
      state.ml = Space.qnums[i].ml; //lz
    }
    else if(Parameters.basis == "finite_J"){
      state.t = Space.qnums[i].t; //tz
      state.j = Space.qnums[i].j; //j
      state.par = Space.qnums[i].par; //par
    }
    for(int k = 0; k < Chan.size3; ++k){
      if(Parameters.basis == "finite_M"){
	if(state.t == Chan.qnums3[k].t && state.m == Chan.qnums3[k].m && state.par == Chan.qnums3[k].par){
	  Chan.indvec.push_back(k);
	  OBvec[k].push_back(i);
	  goto stop;
	}
      }
      else if(Parameters.basis == "finite_HO"){
	if(state.t == Chan.qnums3[k].t && state.m == Chan.qnums3[k].m && state.ml == Chan.qnums3[k].ml){
	  Chan.indvec.push_back(k);
	  OBvec[k].push_back(i);
	  goto stop;
	}
      }
      else if(Parameters.basis == "finite_J"){
	if(state.t == Chan.qnums3[k].t && state.j == Chan.qnums3[k].j && state.par == Chan.qnums3[k].par){
	  Chan.indvec.push_back(k);
	  OBvec[k].push_back(i);
	  goto stop;
	}
      }
      if(k == Chan.size3 - 1){
	++Chan.size3;
	OBvec.resize(Chan.size3);
	Chan.qnums3.push_back(state);
	Chan.indvec.push_back(Chan.size3 - 1);
	OBvec[Chan.size3 - 1].push_back(i);
	break;
      }
    }
  stop:;
  }

  std::cout << "test1, " << Space.Chansize_2b << std::endl;

  //Make Vector of Two-Body States for Each Channel
  Chan.size1 = Space.Chansize_2b;
  Chan.qnums1.resize(Chan.size1);
  std::vector<std::vector<int> > tempvec1(Chan.size3);
  for(int i = 0; i < Chan.size3; ++i){
    int ind1;
    tempvec1[i].assign(Chan.size1, -1);
    for(int j = 0; j < Chan.size3; ++j){
      if(Parameters.basis != "finite_J"){
	plus(state, Space.qnums[i], Space.qnums[j]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	std::cout << i << " " << j << " : " << ind1 << std::endl;
	std::cout << i << " " << tempvec1.size() << ", " << ind1 << " " << tempvec1[i].size() << " " << Chan.qnums1.size() << std::endl;
	tempvec1[i][ind1] = j;
	Chan.qnums1[ind1] = state;
      }
      else if(Parameters.basis == "finite_J"){
	plus(state, Space.qnums[i], Space.qnums[j]);
	for(int J = abs(Space.qnums[i].j - Space.qnums[j].j); J < Space.qnums[i].j + Space.qnums[j].j; J+=2){
	  state.j = J;
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
	  tempvec1[i][ind1] = j;
	  Chan.qnums1[ind1] = state;
	}
      }
    }
  }
  
  std::cout << "test2" << std::endl;

  Chan.tbvec1.resize(Chan.size1);
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.size1; ++j){
      if(tempvec1[i][j] == -1){ continue; }
      int ind1, ind2;
      int tbsize1, tbsize2;
      ind1 = i;
      ind2 = tempvec1[i][j];
      tbsize1 = int(OBvec[ind1].size());
      tbsize2 = int(OBvec[ind2].size());
      for(int tb1 = 0; tb1 < tbsize1; ++tb1){
	for(int tb2 = 0; tb2 < tbsize2; ++tb2){
	  Chan.tbvec1[j].push_back(OBvec[ind1][tb1]);
	  Chan.tbvec1[j].push_back(OBvec[ind2][tb2]);
	}
      }
    }
  }

  std::cout << "test3" << std::endl;

  Chan.obvec1.resize(Chan.size3);
  for(int i = 0; i < Chan.size3; ++i){
    for(int ob = 0; ob < int(OBvec[i].size()); ++ob){
      Chan.obvec1[i].push_back(OBvec[i][ob]);
    }
  }

  std::cout << "test4" << std::endl;

  Chan.tb.resize(Chan.size1);
  Chan.ob.resize(Chan.size3);
  for(int i = 0; i < Chan.size1; ++i){ Chan.tb[i] = int(Chan.tbvec1[i].size()/2); }
  for(int i = 0; i < Chan.size3; ++i){ Chan.ob[i] = int(Chan.obvec1[i].size()); }
}

//ignores existing holes and particles and fills them from protons and neutrons
void Separate_Particles_Holes(Single_Particle_States &States, const HF_Channels &Chan)
{
  //define holes/particles
  double tempen;
  int ind;
  for(int i = 0; i < Chan.size3; ++i){
    //std::cout << "** test1" << std::endl;
    States.h[i] = 0;
    States.p[i] = 0;
    for(int j = 0; j < Chan.ob[i]; ++j){
      //std::cout << "** test2" << std::endl;
      tempen = States.energies[i][j];
      ind = 0;
      for(int k = 0; k < Chan.size3; ++k){
	//std::cout << "** test3" << std::endl;
	for(int l = 0; l < Chan.ob[k]; ++l){
	  //std::cout << "** test4" << std::endl;
	  if(States.energies[k][l] <= tempen && Chan.qnums3[i].t == Chan.qnums3[k].t){ ++ind; }
	}
      }
      //std::cout << "** test5" << std::endl;
      if(ind <= States.hp && Chan.qnums3[i].t == -1){ ++States.h[i]; }
      else if(ind > States.hp && Chan.qnums3[i].t == -1){ ++States.p[i]; }
      if(ind <= States.hn && Chan.qnums3[i].t == 1){ ++States.h[i]; }
      else if(ind > States.hn && Chan.qnums3[i].t == 1){ ++States.p[i]; }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    //std::cout << "** test6" << std::endl;
    States.holes[i].resize(States.h[i]);
    States.particles[i].resize(States.p[i]);
    States.h_energies[i].resize(States.h[i]);
    States.pt_energies[i].resize(States.p[i]);
    for(int j = 0; j < Chan.ob[i]; ++j){
      //std::cout << "** test7" << std::endl;
      if(j < States.h[i]){
	States.holes[i][j] = States.vectors[i][j];
	States.h_energies[i][j] = States.energies[i][j];
      }
      else{
	States.particles[i][j - States.h[i]] = States.vectors[i][j];
	States.pt_energies[i][j - States.h[i]] = States.energies[i][j];
      }
      //std::cout << "** test8" << std::endl;
    }
  }
}

void Build_Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &States)
{
  std::cout << "Building Single-Particle States" << std::endl;
  std::cout << "-------------------------------" << std::endl;

  int ind; // count for filling states
  double tempen; // temp energy for ordering states
  if(Parameters.Pshells != 0){ States.hp = Parameters.P; }
  else{ States.hp = 0; }
  if(Parameters.Nshells != 0){ States.hn = Parameters.N; }
  else{ States.hn = 0; }

  std::cout << "States.hp, States.hn = " << States.hp << ", " << States.hn << std::endl;

  States.vectors.resize(Chan.size3);
  States.energies.resize(Chan.size3);
  States.h.resize(Chan.size3);
  States.p.resize(Chan.size3);
  States.holes.resize(Chan.size3);
  States.particles.resize(Chan.size3);
  States.h_energies.resize(Chan.size3);
  States.pt_energies.resize(Chan.size3);
  for(int i = 0; i < Chan.size3; ++i){
    States.vectors[i].resize(Chan.ob[i]);
    States.energies[i].assign(Chan.ob[i], 0.0);
    States.h[i] = 0;
    States.p[i] = 0;
    for(int j = 0; j < Chan.ob[i]; ++j){
      States.energies[i][j] = Space.qnums[Chan.obvec1[i][j]].energy;
      States.vectors[i][j].assign(Chan.ob[i], 0.0);
      States.vectors[i][j][j] = 1.0;
      if(Space.qnums[Chan.obvec1[i][j]].type == "hole"){ ++States.h[i]; }
      else if(Space.qnums[Chan.obvec1[i][j]].type == "particle"){ ++States.p[i]; }
    }
  }
  
  // Order states by energy
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i] - 1; ++j){
      ind = j;
      tempen = States.energies[i][j];
      for(int k = j + 1; k < Chan.ob[i]; ++k){
	if(States.energies[i][k] < tempen){ tempen = States.energies[i][k]; ind = k; }
      }
      std::swap(States.energies[i][j], States.energies[i][ind]);
      std::swap(States.vectors[i][j], States.vectors[i][ind]);
    }
  }

  //Separate States
  Separate_Particles_Holes(States, Chan);
}

void Read_Matrix_Elements_J(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  std::string fullpath1; // file path string and string with "_M" added (for J,M-Scheme)
  int NumElements; // number of ME read in from file
  std::string number; // string for first word of each line
  std::ifstream interaction;	// interaction file
  std::string interactionline; // interaction file line
  std::istringstream interactionstream; // stream of file line string
  size_t index1, index2; // indicies for finding parameters among file lines
  double TBME, hom, r2, p2; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4, coupJ, coupT, par; // interaction file contents
  int ind1, ind;
  State tb;

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

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
  while (number != "Tz"){
    getline(interaction, interactionline);
    interactionstream.str(interactionline);
    interactionstream >> number;
  }

  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> coupT >> par >> coupJ >> shell1 >> shell2 >> shell3 >> shell4 >> TBME >> hom >> r2 >> p2;
    coupJ *= 0.5;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if((shell1 == shell2 || shell3 == shell4) && coupJ%2 != 0){ continue; }
    //if(shell1 == shell2){ TBME *= sqrt(2.0); }
    //if(shell3 == shell4){ TBME *= sqrt(2.0); }
    if(shell1 > shell2){ std::swap(shell1, shell2); TBME *= pow(-1.0, Space.qnums[shell1].j + Space.qnums[shell2].j - coupJ); }
    if(shell3 > shell4){ std::swap(shell3, shell4); TBME *= pow(-1.0, Space.qnums[shell3].j + Space.qnums[shell4].j - coupJ); }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    tb.j = coupJ;
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
    ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell1, shell2, shell3, shell4);
    ME.V[ind1][ind] = TBME;
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
  int ind1, ind;
  State tb;
  ME = HF_Matrix_Elements(Chan);

  fullpath1 = PATH + Parameters.MatrixElements + ".int";

  interaction.open(fullpath1.c_str());
  if(!interaction.is_open()){ std::cerr << "Matrix Element file, " << Parameters.MatrixElements << ", does not exist" << std::endl; exit(1); }
  getline(interaction, interactionline);
  interactionstream.str(interactionline);
  interactionstream >> NumElements;

  //std::cout << "!! " << NumElements << " !!" << std::endl;
  for(int i = 0; i < NumElements; ++i){
    getline(interaction, interactionline);
    std::istringstream(interactionline) >> shell1 >> shell2 >> shell3 >> shell4 >> TBME;
    TBME *= Parameters.tbstrength;
    shell1 -= 1;
    shell2 -= 1;
    shell3 -= 1;
    shell4 -= 1;
    if(shell1 == shell2 || shell3 == shell4){ continue; }
    if(shell1 > shell2){ std::swap(shell1, shell2); TBME *= -1.0; }
    if(shell3 > shell4){ std::swap(shell3, shell4); TBME *= -1.0; }
    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);

    ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell1, shell2, shell3, shell4);
    ME.V[ind1][ind] = TBME;
    ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell2, shell1, shell3, shell4);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell1, shell2, shell4, shell3);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell2, shell1, shell4, shell3);
    ME.V[ind1][ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell3, shell4, shell1, shell2);
      ME.V[ind1][ind] = TBME;
      ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell3, shell4, shell2, shell1);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell4, shell3, shell1, shell2);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], shell4, shell3, shell2, shell1);
      ME.V[ind1][ind] = TBME;
    }
    //std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << ind1 << " " << ind << ", " << TBME << std::endl;
  }
  interaction.close();

  /*for(int chan = 0; chan < Chan.size1; ++chan){
    for(int tb1 = 0; tb1 < Chan.tb[chan]; ++tb1){ std::cout << Chan.tbvec1[chan][2*tb1] << Chan.tbvec1[chan][2*tb1 + 1] << " "; }
    std::cout << std::endl;
    for(int tb1 = 0; tb1 < Chan.tb[chan]; ++tb1){
      for(int tb2 = 0; tb2 < Chan.tb[chan]; ++tb2){
	std::cout << ME.V[chan][Chan.tb[chan]*tb1 + tb2] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    }*/
}

void Read_Matrix_Elements_HO(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan, HF_Matrix_Elements &ME)
{
  int pind, qind, rind, sind;
  double TBME;
  ME = HF_Matrix_Elements(Chan);
  for(int i = 0; i < Chan.size1; ++i){
    for(int pq = 0; pq < Chan.tb[i]; ++pq){
      pind = Chan.tbvec1[i][2*pq];
      qind = Chan.tbvec1[i][2*pq + 1];
      if(pind == qind){ continue; }
      for(int rs = 0; rs < Chan.tb[i]; ++rs){
	rind = Chan.tbvec1[i][2*rs];
	sind = Chan.tbvec1[i][2*rs + 1];
	if(rind == sind){ continue; }
	TBME = Coulomb_HO(Parameters, Space, pind, qind, rind, sind);
	ME.V[i][pq * Chan.tb[i] + rs] = TBME;
	if(pind < qind && rind < sind && pind <= rind && pq <= rs){
	  std::cout << pind << " " << qind << " " << rind << " " << sind << " : " << TBME << std::endl;
	}
      }
    }
  }
}

void Hartree_Fock_States_J(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, const HF_Matrix_Elements &ME)
{
  double error; // Energy error between iterations
  double Bshift; // level shift parameter
  int ind; // Index to keep track of iteration number
  std::vector<double> densmat, fock; // Fock Matrix
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  std::vector<double> w, work; // EigenEnergy Vector and work vector
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  int minj, maxj;
  State tb;

  jobz = 'V';
  uplo = 'U';
  Bshift = 50.0;

  ind = 0;
  error = 1000;
  while((error > 1e-8 && ind < 100) || ind < 10){
    //std::cout << "!! " << ind << " !!" << std::endl;

    ++ind;
    error = 0.0;
    
    //Make Fock Matrix
    for(int i = 0; i < Chan.size3; ++i){
      int size1 = Chan.ob[i];
      fock.assign(size1 * size1, 0.0);
      for(int m = 0; m < size1; ++m){
	int mind = Chan.obvec1[i][m];
	fock[size1 * m + m] += Space.qnums[mind].energy; // Add diagonal elements to fock matrices
	for(int l = 0; l < size1; ++l){
	  int lind = Chan.obvec1[i][l];
	  for(int j = 0; j < Chan.size3; ++j){
	    int size2 = Chan.ob[j];
	    minj = abs(Chan.qnums3[i].j - Chan.qnums3[j].j);
	    maxj = Chan.qnums3[i].j + Chan.qnums3[j].j;
	    for(int beta = 0; beta < HF.h[j]; ++beta){ // Sum over occupied levels
	      for(int J = minj; J <= maxj; ++J){
		for(int n = 0; n < size2; ++n){
		  int nind = Chan.obvec1[j][n];
		  for(int k = 0; k < size2; ++k){
		    int kind = Chan.obvec1[j][k];
		    plus(tb, Space.qnums[mind], Space.qnums[nind]);
		    tb.j = J;
		    int ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    int ind2 = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], mind, nind, lind, kind);
		    term = HF.vectors[j][beta][n] * HF.vectors[j][beta][k] * ME.V[ind1][ind2];
		    term *= (2*J + 1)/(2*Chan.qnums3[i].j + 1);
		    fock[size1 * m + l] += term;
		  }
		}
	      }
	    }
	  }
	  for(int beta = 0; beta < HF.h[i]; ++beta){ // Sum over occupied levels
	    fock[size1 * m + l] -= HF.vectors[i][beta][m] * HF.vectors[i][beta][l] * Bshift;
	  }
	}
      }    
      for(int j = 0; j < int(fock.size()); ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      
      lda = size1;
      lwork = (3+2)*size1;
      w.assign(size1,0);
      work.assign(lwork,0);
    
      if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, & *fock.begin(), &lda, & *w.begin(), & *work.begin(), &lwork, &info); }
      for(int j = 0; j < int(fock.size()); ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      for(int j = 0; j < HF.h[i]; ++j){ w[j] += Bshift; } //Add back Level-shift parameter
      
      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  HF.vectors[i][j][k] = fock[size1 * j + k];
	}
      }
      
      int ind2;
      double tempen2;
      // Order states by energy
      for(int j = 0; j < size1 - 1; ++j){
	ind2 = j;
	tempen2 = w[j];
	for(int k = j + 1; k < size1; ++k){
	  if(w[k] < tempen2){ tempen2 = w[k]; ind2 = k; }
	}
	std::swap(w[j], w[ind2]);
	std::swap(HF.vectors[i][j], HF.vectors[i][ind2]);
      }
      
      for(int j = 0; j < size1; ++j){
	error += fabs(HF.energies[i][j] - w[j])/Chan.ob[i]; 
	HF.energies[i][j] = w[j];
      }
      GramSchmidt(HF.vectors[i]);
      Separate_Particles_Holes(HF, Chan);
      //Separate_Particles_Holes(HF, Space.Pocc, Space.Nocc, Space);
      
      //std::cout << "error = " << error << std::endl << std::endl;
    }
  }

  ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i]; ++j){
      Space.qnums[ind] = Chan.qnums3[i];
      Space.qnums[ind].energy = HF.energies[i][j];
      if(j < HF.h[i]){ Space.qnums[ind].type = "hole"; }
      else{ Space.qnums[ind].type = "particle"; }
      ++ind;
    }
  }

  /*std::ofstream HFlevelfile;
  std::string filename = PATH + Parameters.LevelScheme + "_HF.sp";
  HFlevelfile.open(filename.c_str());
  HFlevelfile << "Mass number A of chosen nucleus (important for CoM corrections): \t" << Space.A << "\n";
  HFlevelfile << "Oscillator energy: \t" << Space.HOEnergy << "\n";
  HFlevelfile << "Total number of single-particle orbits: \t" << Space.shelltot << "\n";
  HFlevelfile << "Legend:   \tn \tl \t2j \ttz \t2n+l \tHO-energy \tevalence \tparticle/hole \tinside/outside \n";
  for(int i = 0; i < int(HF.protons.size()); ++i){
    HFlevelfile << "Number:   " << i+1 << "\t" << p_n[i] << "\t" << HF.p_l[i] << "\t" << int(2*HF.p_j[i]) << "\t";
    HFlevelfile << "-1" << "\t" << 2*p_n[i]+HF.p_l[i] << "\t" << std::setprecision(8) << HF.p_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Pocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  for(int i = 0; i < int(HF.neutrons.size()); ++i){
    HFlevelfile << "Number:   " << i+int(HF.protons.size())+1 << "\t" << n_n[i] << "\t" << HF.n_l[i] << "\t" << int(2*HF.n_j[i]) << "\t";
    HFlevelfile << "1" << "\t" << 2*n_n[i]+HF.n_l[i] << "\t" << std::setprecision(8) << HF.n_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Nocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  
  HFlevelfile.close();*/
  
}

void Hartree_Fock_States(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, const HF_Matrix_Elements &ME)
{
  double error; // Energy error between iterations
  double Bshift; // level shift parameter
  int ind; // Index to keep track of iteration number
  std::vector<double> densmat, fock; // Fock Matrix
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  std::vector<double> w, work; // EigenEnergy Vector and work vector
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  Single_Particle_States States = HF;
  State tb;
  
  jobz = 'V';
  uplo = 'U';
  Bshift = 0.0;

  ind = 0;
  error = 1000;

  /*for(int i = 0; i < Chan.size3; ++i){
    int size1 = Chan.ob[i];
    std::cout << "h,p = " << HF.h[i] << " " << HF.p[i] << std::endl;
    for(int j = 0; j < size1; ++j){
      std::cout << "E = " << HF.energies[i][j] << std::endl;
      for(int k = 0; k < size1; ++k){
	std::cout << HF.vectors[i][j][k] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    }*/

  while((error > 1e-18 && ind < 1000) || ind < 10){ //10
    //std::cout << "!! " << ind << " !!" << std::endl;

    ++ind;
    error = 0.0;
    
    //Make Fock Matrix
    for(int i = 0; i < Chan.size3; ++i){
      int size1 = Chan.ob[i];
      fock.assign(size1 * size1, 0.0);
      //std::cout << "chan = " << i << std::endl;
      for(int m = 0; m < size1; ++m){
	int mind = Chan.obvec1[i][m];
	fock[size1 * m + m] += Space.qnums[mind].energy; // Add diagonal elements to fock matrices
	//std::cout << "mind = " << mind << ", <m|f|m> += " << Space.qnums[mind].energy << std::endl;
	for(int l = 0; l < size1; ++l){
	  int lind = Chan.obvec1[i][l];
	  //std::cout << "lind = " << lind << std::endl;
	  for(int j = 0; j < Chan.size3; ++j){
	    int size2 = Chan.ob[j];
	    //std::cout << "chan2 = " << j << std::endl;
	    for(int beta = 0; beta < HF.h[j]; ++beta){ // Sum over occupied levels
	      for(int n = 0; n < size2; ++n){
		int nind = Chan.obvec1[j][n];
		for(int k = 0; k < size2; ++k){
		  int kind = Chan.obvec1[j][k];
		  plus(tb, Space.qnums[mind], Space.qnums[nind]);
		  int ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  int ind2 = Index22(Chan.tbvec1[ind1], Chan.tbvec1[ind1], Chan.tb[ind1], Chan.tb[ind1], mind, nind, lind, kind);
		  term = HF.vectors[j][beta][n] * HF.vectors[j][beta][k] * ME.V[ind1][ind2];
		  //std::cout << "B, nind, kind = " << beta << " " << nind << " " << kind << std::endl;
		  //std::cout << "HF[B][n], HF[B][k], ME.V, term = " << HF.vectors[j][beta][n] << " " << HF.vectors[j][beta][k] << " " << ME.V[ind1][ind2] << " " << term << std::endl;
		  fock[size1 * m + l] += term;
		}
	      }
	    }
	  }
	  for(int beta = 0; beta < HF.h[i]; ++beta){ // Sum over occupied levels
	    fock[size1 * m + l] -= HF.vectors[i][beta][m] * HF.vectors[i][beta][l] * Bshift;
	  }
	}
      }    
      for(int j = 0; j < int(fock.size()); ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}

      /*std::cout << "chan = " << i << std::endl;
      for(int x = 0; x < size1; ++x){
	for(int y = 0; y < size1; ++y){
	  std::cout << fock[size1*x + y] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl*/
	
      lda = size1;
      lwork = (3+2)*size1;
      w.assign(size1,0);
      work.assign(lwork,0);
    
      if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, & *fock.begin(), &lda, & *w.begin(), & *work.begin(), &lwork, &info); }
      for(int j = 0; j < int(fock.size()); ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      for(int j = 0; j < HF.h[i]; ++j){ w[j] += Bshift; } //Add back Level-shift parameter
      
      /*for(int x = 0; x < size1; ++x){
	for(int y = 0; y < size1; ++y){
	  std::cout << fock[size1*x + y] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
	
      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  error += fabs((States.vectors[i][j][k] - fock[size1 * j + k])/States.vectors[i][j][k]);
	  States.vectors[i][j][k] = fock[size1 * j + k];
	}
      }

      int ind2;
      double tempen2;
      // Order states by energy
      for(int j = 0; j < size1 - 1; ++j){
	ind2 = j;
	tempen2 = w[j];
	for(int k = j + 1; k < size1; ++k){
	  if(w[k] < tempen2){ tempen2 = w[k]; ind2 = k; }
	}
	std::swap(w[j], w[ind2]);
	std::swap(States.vectors[i][j], States.vectors[i][ind2]);
      }
      
      for(int j = 0; j < size1; ++j){
	States.energies[i][j] = w[j];
      }
      GramSchmidt(States.vectors[i]);
      Separate_Particles_Holes(States, Chan);
      
      //std::cout << "error = " << error << std::endl << std::endl;
    }
    error /= (Space.indtot * Space.indtot);
    HF = States;
  }

  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i]; ++j){
      Space.qnums[Chan.obvec1[i][j]] = Chan.qnums3[i];
      Space.qnums[Chan.obvec1[i][j]].energy = HF.energies[i][j];
      if(j < HF.h[i]){ Space.qnums[Chan.obvec1[i][j]].type = "hole"; }
      else{ Space.qnums[Chan.obvec1[i][j]].type = "particle"; }
    }
  }

  for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].n << " " << Space.qnums[i].par << " " << Space.qnums[i].j << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
  }


  /*std::ofstream HFlevelfile;
  std::string filename = PATH + Parameters.LevelScheme + "_HF.sp";
  HFlevelfile.open(filename.c_str());
  HFlevelfile << "Mass number A of chosen nucleus (important for CoM corrections): \t" << Space.A << "\n";
  HFlevelfile << "Oscillator energy: \t" << Space.HOEnergy << "\n";
  HFlevelfile << "Total number of single-particle orbits: \t" << Space.shelltot << "\n";
  HFlevelfile << "Legend:   \tn \tl \t2j \ttz \t2n+l \tHO-energy \tevalence \tparticle/hole \tinside/outside \n";
  for(int i = 0; i < int(HF.protons.size()); ++i){
    HFlevelfile << "Number:   " << i+1 << "\t" << p_n[i] << "\t" << HF.p_l[i] << "\t" << int(2*HF.p_j[i]) << "\t";
    HFlevelfile << "-1" << "\t" << 2*p_n[i]+HF.p_l[i] << "\t" << std::setprecision(8) << HF.p_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Pocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  for(int i = 0; i < int(HF.neutrons.size()); ++i){
    HFlevelfile << "Number:   " << i+int(HF.protons.size())+1 << "\t" << n_n[i] << "\t" << HF.n_l[i] << "\t" << int(2*HF.n_j[i]) << "\t";
    HFlevelfile << "1" << "\t" << 2*n_n[i]+HF.n_l[i] << "\t" << std::setprecision(8) << HF.n_energies[i] << "\t" << "0.000000" << "\t";
    if(i < Space.Nocc){ HFlevelfile << "hole    " << "\t" << "inside" << "\n"; }
    else{ HFlevelfile << "particle" << "\t" << "inside" << "\n"; }
  }
  
  HFlevelfile.close();*/
  
}


void Convert_To_HF_Matrix_Elements(const HF_Channels &Chan, const Single_Particle_States &States, HF_Matrix_Elements &ME)
{
  std::cout << "Converting Matrix Elements to HF basis" << std::endl;
  std::cout << "--------------------------------------" << std::endl;

  int length, matlength; // max length of M-Scheme indicies, length of J_ME
  double tempel;
  std::vector<double> M1, M2; // Matrices of coefficients
  //std::vector<double> V;
  std::vector<double> C;
  char transa;
  double alpha1, beta1;
  std::ofstream jschemefile; // file to print M-Scheme matrix elements

  /*for(int chan = 0; chan < Chan.size3; ++chan){
    for(int i = 0; i < Chan.ob[chan]; ++i){
      std::cout << "! " << i << " : ";
      for(int j = 0; j < Chan.ob[chan]; ++j){
	std::cout << States.vectors[chan][i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    }*/
      
  for(int chan = 0; chan < Chan.size1; ++chan){
    length = Chan.tb[chan];
    matlength = pow(length, 2.0);
    M1.assign(matlength, 0.0);
    M2.assign(matlength, 0.0);
    C.assign(matlength, 0.0);

    for(int pq = 0; pq < length; ++pq){
      for(int ag = 0; ag < length; ++ag){
	int pind = Chan.tbvec1[chan][2*pq];
	int qind = Chan.tbvec1[chan][2*pq + 1];
	int aind = Chan.tbvec1[chan][2*ag];
	int gind = Chan.tbvec1[chan][2*ag + 1];
	//td::cout << "!! " << pind << " " << qind << " " << aind << " " << gind << std::endl;
	if(Chan.indvec[pind] != Chan.indvec[aind] || Chan.indvec[qind] != Chan.indvec[gind]){ continue; }
	int ind1 = Chan.indvec[pind];
	int ind2 = Chan.indvec[qind];
	//std::cout << "pass! " << ind1 << " " << ind2 << std::endl;
	pind = Index1(Chan.obvec1[ind1], pind);
	qind = Index1(Chan.obvec1[ind2], qind);
	aind = Index1(Chan.obvec1[ind1], aind);
	gind = Index1(Chan.obvec1[ind2], gind);
	//std::cout << "?? " << pind << " " << aind << " = " << States.vectors[ind1][pind][aind] << std::endl;
	//std::cout << "?? " << qind << " " << gind << " = " << States.vectors[ind2][qind][gind] << std::endl;
	tempel = States.vectors[ind1][pind][aind] * States.vectors[ind2][qind][gind];
	M1[length * pq + ag] = tempel;
	M2[length * ag + pq] = tempel;
      }
    }

    /*std::cout << "M1" << std::endl;
    for(int i = 0; i < length; ++i){ std::cout << Chan.tbvec1[chan][2*i] << Chan.tbvec1[chan][2*i+1] << " "; }
    std::cout << std::endl;
    for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
	std::cout << M1[length*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "V0" << std::endl;
    for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
	std::cout << ME.V[chan][length*i + j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;*/

    /*for(int ag = 0; ag < length; ++ag){
      int aind = Chan.tbvec1[chan][2*ag];
      int gind = Chan.tbvec1[chan][2*ag + 1];
      for(int bd = 0; bd < length; ++bd){
	int bind = Chan.tbvec1[chan][2*bd];
	int dind = Chan.tbvec1[chan][2*bd + 1];
	tempel = ME.V[chan][ag*Chan.tb[chan] + bd];
	if(aind == gind){ tempel /= sqrt(2.0); }
	if(bind == dind){ tempel /= sqrt(2.0); }
	V[length * ag + bd] = tempel;
      }
      }*/
      
    transa = 'N';
    alpha1 = 1.0;
    beta1 = 0.0;
    
    if(length != 0){
      RM_dgemm(& *M1.begin(), & *ME.V[chan].begin(), & *C.begin(), &length, &length, &length, &alpha1, &beta1, &transa, &transa);
      /*std::cout << "V1" << std::endl;
      for(int i = 0; i < length; ++i){
	for(int j = 0; j < length; ++j){
	  std::cout << C[length*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
      RM_dgemm(& *C.begin(), & *M2.begin(), & *ME.V[chan].begin(), &length, &length, &length, &alpha1, &beta1, &transa, &transa);
      /*std::cout << "V2" << std::endl;
      for(int i = 0; i < length; ++i){
	for(int j = 0; j < length; ++j){
	  std::cout << ME.V[chan][length*i + j] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/
    }
  }
    /*for(int p = 0; p < plength; ++p){
      for(int q = p; q < plength; ++q){
	for(int r = p; r < plength; ++r){
	  for(int s = r; s < plength; ++s){
	    num1 = number_diff(p, q, plength);
	    num2 = number_diff(r, s, plength);
	    if(num1 > num2){ std::swap(num1, num2); }
	    tempel = PP_V[pmatlength * num1 + num2];
	    if(p == q){ tempel /= sqrt(2.0); }
	    if(r == s){ tempel /= sqrt(2.0); }
	    if(fabs(tempel) > 1.0e-10){ 
	      HF_ME.set_ppJME(p, q, r, s, J, tempel);
	    }
	  }
	}
      }
    }
  }
  
    HF_ME.cut_ppJME();
    
    ind = 0;
    for(int J = 0; J <= Space.max2J; ++J){
      for(int m = 0; m < plength; ++m){
	for(int n = m; n < plength; ++n){
	  for(int l = m; l < plength; ++l){
	    for(int k = l; k < plength; ++k){
	      tempel = HF_ME.get_ppJME(m, n, l, k, J);
	      if(fabs(tempel) >= 1.0e-10){ ++ind; }
	    }
	  }
	}
      }
    }
    }*/

  // print J_ME to file
  /*jschemefile.open((PATH + MatrixElements + "_HF.int").c_str());
  jschemefile << "Total number of twobody matx elements:" << "\t" << ind << "\n";
  jschemefile << "----> Interaction part\n";   
  jschemefile << "Nucleon-Nucleon interaction model:n3lo\n";            
  jschemefile << "Type of calculation: nocore\n";              
  jschemefile << "Number and value of starting energies:   1	0.000000E+00\n";
  jschemefile << "Total number of twobody matx elements:\t" << ind << "\n";
  jschemefile << "Tz      Par      2J      a      b      c      d      <ab|V|cd>\n";
  for(int J = 0; J <= Space.max2J; ++J){
    for(int m = 0; m < plength; ++m){
      for(int n = m; n < plength; ++n){
	for(int l = m; l < plength; ++l){
	  for(int k = l; k < plength; ++k){
	    tempel = HF_ME.get_ppJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t -1\t" << pow(-1.0,States.p_l[m]+States.p_l[n]) << "\t" << 2*J << "\t"
			<< m + 1 << "\t" << n + 1 << "\t" << l + 1 << "\t" << k + 1
			<< "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
    for(int m = 0; m < plength; ++m){
      for(int n = 0; n < nlength; ++n){
	for(int l = m; l < plength; ++l){
	  for(int k = 0; k < nlength; ++k){
	    tempel = HF_ME.get_pnJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t  0\t" << pow(-1.0,States.p_l[m]+States.n_l[n]) << "\t" << 2*J << "\t"
			<< m + 1 << "\t" << n + plength + 1 << "\t" << l + 1 << "\t" << k + plength + 1 
			<< "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
    for(int m = 0; m < nlength; ++m){
      for(int n = m; n < nlength; ++n){
	for(int l = m; l < nlength; ++l){
	  for(int k = l; k < nlength; ++k){
	    tempel = HF_ME.get_nnJME(m, n, l, k, J);
	    if(fabs(tempel) < 1.0e-10){ continue; }
	    jschemefile << "\t  1\t" << pow(-1.0,States.n_l[m]+States.n_l[n]) << "\t" << 2*J << "\t"
			<< m + plength + 1 << "\t" << n + plength + 1 << "\t" << l + plength + 1 << "\t"
			<< k + plength + 1 << "\t" << std::setprecision(8) << tempel << "\n";
	  }
	}
      }
    }
  }
  
  jschemefile.close();*/
}

/*void Setup_HF_Space(Model_Space &Space, const Single_Particle_States &States, const HF_Channels &Chan)
{
  int ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i]; ++j){
      Space.qnums[ind] = Chan.qnums3[i];
      Space.qnums[ind].energy = HF.energies[i][j];
      if(j < Chan.h[i]){ Space.qnums[ind].type = "hole"; }
      else{ Space.qnums[ind].type = "particle"; }
      ++ind;
    }
  }

  double tempen;
  // Order states by energy
  for(int i = 0; i < Space.indtot - 1; ++i){
    ind = i;
    tempen = Space.qnums[i].energy;
    for(int j = i + 1; j < Space.indtot; ++j){
      if(Space.qnums[j].energy < tempen){ tempen = Space.qnums[j].energy; ind = j; }
    }
    std::swap(Space.qnums[i], Space.qnums[ind]);
  }

  // find j
  ind = 0;
  for(int i = 0; i < Space.indtot; ++i){
    tempen2 = en[i];
    ind2 = 0;
    for(int j = i + 1; j < Space.indtot; ++j){
      if(en[j] == tempen2){ ++ind2; }
      else{
	for(int k = i; k < j; ++k){ J[j] = ind2; }
	i = j - 1;
      }
    }
  }

  for(int i = 0; i < Space.indtot; ++i){
    for(int k = 0; k < Chan.ob[chan[i]]; ++k){
      if(abs(2*Space.levelsl[Chan.obvec1[chan[i]][k]] - J[ind]) == 1){
	l[i] = Space.levelsl[Chan.obvec1[chan[i]][k]];
	break;
      }
    }
  }

  for(int i = 0; i < Space.indtot; ++i){
    int tempn = -1;
    int templ = l[i];
    double tempj = J[i];
    for(int j = 0; j <= i; ++j){
      if(l[j] == templ && J[j] == tempj){ ++tempn; }
    }
    n[i] = tempn;

  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i]; ++j){
      int ind = Chan.obvec1[i][j];
      Space.qnums[ind].energy = States.energies[i][j];
      if(j < States.h[i]){ Space.qnums[ind].type = "hole"; }
      else{ Space.qnums[ind].type = "particle"; }
    }
  }

  // find j
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.ob[i]; ++j){
      int ind1 = Chan.obvec1[i][j];
      double tempen2 = Space.qnums[ind1].energy;
      int ind = 0;
      for(int k = 0; k < Chan.size3; ++k){
	for(int l = 0; l < Chan.ob[k]; ++l){
	  int ind2 = Chan.obvec1[k][l];
	  if(Space.qnums[ind2].energy == tempen2){ ++ind; break; }
	}
      }
      Space.qnums[ind1].j = ind - 1;
    }
    }
    }*/

void Get_Matrix_Elements(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME; // interaction two-body interaction ME and two-body COM ME
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;
  State tb;

  /*for(int i = 0; i < Space.indtot; ++i){
    std::cout << Space.qnums[i].par << " " << Space.qnums[i].m << " " << Space.qnums[i].t << " " << Space.qnums[i].energy << " " << Space.qnums[i].type << std::endl;
    }*/

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    for(int tb1 = 0; tb1 < HF_Chan.tb[chan]; ++tb1){
      shell1 = HF_Chan.tbvec1[chan][2*tb1];
      shell2 = HF_Chan.tbvec1[chan][2*tb1 + 1];
      ptype = Space.qnums[shell1].type;
      qtype = Space.qnums[shell2].type;
      //std::cout << "shell1, shell2 = " << shell1 << " " << shell2 << std::endl;
      if(ptype == qtype && shell1 >= shell2){ continue; }
      if(ptype == "particle" && qtype == "hole"){ continue; }
      for(int tb2 = 0; tb2 < HF_Chan.tb[chan]; ++tb2){
	shell3 = HF_Chan.tbvec1[chan][2*tb2];
	shell4 = HF_Chan.tbvec1[chan][2*tb2 + 1];
	rtype = Space.qnums[shell3].type;
	stype = Space.qnums[shell4].type;
	//std::cout << "shell3, shell4 = " << shell3 << " " << shell4 << std::endl;
	if(rtype == stype && shell3 >= shell4){ continue; }
	if(rtype == "particle" && stype == "hole"){ continue; }
	if(ptype == rtype && shell1 > shell3){ continue; }
	if(ptype == rtype && qtype == stype && shell1 == shell3 && shell2 > shell4){ continue; }
	if(ptype == "particle" && rtype == "hole"){ continue; }
	//std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << std::endl;
	TBME = HF_ME.V[chan][tb1*HF_Chan.tb[chan] + tb2];
	std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
	if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	}
	else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
	  Ints.D_ME1.V3[ind1][ind] = TBME;
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V3[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V4[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V4[ind1][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell2];
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V5[ind2][ind] = TBME;
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell1];
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V5[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell1];
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V6[ind2][ind] = TBME;
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell2];
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V6[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell4];
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
	  Ints.D_ME1.V7[ind2][ind] = TBME;
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
	  Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell3];
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
	  Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V7[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell3];
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
	  Ints.D_ME1.V8[ind2][ind] = TBME;
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell4];
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
	  Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
	  Ints.D_ME1.V8[ind2][ind] = TBME;
	  
	  minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	  Ints.D_ME1.V9[ind1][ind] = TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	  Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	  Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	  Ints.D_ME1.V9[ind1][ind] = TBME;
	  
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
	  Ints.D_ME1.V10[ind1][ind] = TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
	  Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
	  Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
	  Ints.D_ME1.V10[ind1][ind] = TBME;
	}
	if(Parameters.approx == "singles"){
	  if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.S_ME1.V11[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.S_ME1.V11[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell4, shell3);
	    Ints.S_ME1.V13[ind2][ind] = TBME;
	    ind2 = Chan.indvec[shell4];
	    ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell3, shell4);
	    Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
	    
	    minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell4, shell2, shell1, shell3);
	    Ints.S_ME1.V16[ind2][ind] = TBME;
	    minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V16[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell2];
	    ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell3, shell4, shell2);
	    Ints.S_ME1.V17[ind2][ind] = TBME;
	    ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell4, shell3, shell2);
	    Ints.S_ME1.V17[ind2][ind] = -1.0 * TBME;
	    
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell3, shell4, shell1, shell2);
	    Ints.S_ME1.V20[ind2][ind] = TBME;
	    ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell4, shell3, shell1, shell2);
	    Ints.S_ME1.V20[ind2][ind] = -1.0 * TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.S_ME1.V12[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V12[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell1];
	    ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell2, shell3, shell4, shell1);
	    Ints.S_ME1.V14[ind2][ind] = TBME;
	    ind2 = Chan.indvec[shell2];
	    ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell1, shell3, shell4, shell2);
	    Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;
	    
	    minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V15[ind2][ind] = TBME;
	    minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell1, shell2, shell4);
	    Ints.S_ME1.V15[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell1, shell2, shell4, shell3);
	    Ints.S_ME1.V18[ind2][ind] = TBME;
	    ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell2, shell1, shell4, shell3);
	    Ints.S_ME1.V18[ind2][ind] = -1.0 * TBME;
	    
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell1, shell2);
	    Ints.S_ME1.V19[ind2][ind] = TBME;
	    ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell2, shell1);
	    Ints.S_ME1.V19[ind2][ind] = -1.0 * TBME;
	  }
	}
      }
    }
  }
  /*for(int i = 0; i < Chan.size1; ++i){
    for(int pp = 0; pp < Chan.pp[i]; ++pp){ std::cout << Chan.ppvec1[i][2*pp] << Chan.ppvec1[i][2*pp + 1] << " "; }
    std::cout << std::endl;
    for(int hh = 0; hh < Chan.hh[i]; ++hh){ std::cout << Chan.hhvec1[i][2*hh] << Chan.hhvec1[i][2*hh + 1] << " "; }
    std::cout << std::endl;
    for(int pp = 0; pp < Chan.pp[i]; ++pp){
      for(int hh = 0; hh < Chan.hh[i]; ++hh){
	std::cout << Ints.D_ME1.V4[i][pp*Chan.hh[i] + hh] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    }*/
}

/*void Get_Matrix_Elements_J(const Input_Parameters &Parameters, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double TBME0, TBME, m1, m2, t1, t2, CGC1, CGC2; // interaction two-body interaction ME and two-body COM ME
  int pind, qind, rind, sind;
  int shell1, shell2, shell3, shell4; // interaction file contents
  int ind, ind1, ind2;
  std::string ptype, qtype, rtype, stype;
  State tb;

  for(int chan = 0; chan < HF_Chan.size1; ++chan){
    for(int tb1 = 0; tb1 < HF_Chan.tb[chan]; ++tb1){
      shell1 = HF_Chan.tbvec1[chan][2*tb1];
      shell2 = HF_Chan.tbvec1[chan][2*tb1 + 1];
      ptype = Space.qnums[shell1].type;
      qtype = Space.qnums[shell2].type;
      if(ptype == qtype && shell1 > shell2){ continue; }
      if(ptype == "particle" && qtype == "hole"){ continue; }
      for(int tb2 = tb1; tb2 < HF_Chan.tb[chan]; ++tb2){
	shell3 = HF_Chan.tbvec1[chan][2*tb2];
	shell4 = HF_Chan.tbvec1[chan][2*tb2 + 1];
	rtype = Space.qnums[shell3].type;
	stype = Space.qnums[shell4].type;
	if(rtype == stype && shell3 > shell4){ continue; }
	if(rtype == "particle" && stype == "hole"){ continue; }
	if(ptype == rtype && qtype == stype && shell1 == shell3 && shell2 > shell4){ continue; }
	TBME0 = HF_ME.V[chan][tb1*HF_Chan.tb[chan] + tb2];
	//if(shell1 == shell2){ TBME0 *= sqrt(2.0); }
	//if(shell3 == shell4){ TBME0 *= sqrt(2.0); }
	for(int jz = -HF_Chan.qnums1[chan].j; jz <= HF_Chan.qnums1[chan].j; jz+=2){
	  for(int p = 0; p < int(Space.shellsm[shell1].size()); ++p){
	    for(int q = 0; q < int(Space.shellsm[shell2].size()); ++q){
	      for(int r = 0; r < int(Space.shellsm[shell3].size()); ++r){
		for(int s = 0; s < int(Space.shellsm[shell4].size()); ++s){
		  pind = Space.shellsm[shell1][p];
		  qind = Space.shellsm[shell2][q];
		  rind = Space.shellsm[shell3][r];
		  sind = Space.shellsm[shell4][s];
		  m1 = Space.qnums[pind].m + Space.qnums[qind].m;
		  m2 = Space.qnums[rind].m + Space.qnums[sind].m;
		  t1 = Space.qnums[pind].t + Space.qnums[qind].t;
		  t2 = Space.qnums[rind].t + Space.qnums[sind].t;
		  if(t1 != t2 || m1 != jz || m2 != jz){ continue; }
		  if(pind >= qind && shell1 == shell2){ continue; }
		  if(rind >= sind && shell3 == shell4){ continue; }
		  if(pind > rind && shell1 == shell3 && shell2 == shell4){ continue; }
		  if(pind == rind && qind > sind && shell1 == shell3 && shell2 == shell4){ continue; }
		  CGC1 = CGC(0.5*Space.qnums[pind].j, 0.5*Space.qnums[pind].m, 0.5*Space.qnums[qind].j, 0.5*Space.qnums[qind].m, 0.5*HF_Chan.qnums1[chan].j, double(0.5*jz));
		  CGC2 = CGC(0.5*Space.qnums[rind].j, 0.5*Space.qnums[rind].m, 0.5*Space.qnums[sind].j, 0.5*Space.qnums[sind].m, 0.5*HF_Chan.qnums1[chan].j, double(0.5*jz));
		  TBME = TBME0 * CGC1 * CGC2;
		  if(fabs(TBME) < 1e-12){ continue; }
		  ptype = Space.qnums[pind].type;
		  qtype = Space.qnums[qind].type;
		  rtype = Space.qnums[rind].type;
		  stype = Space.qnums[sind].type;
		  if(ptype == "particle" && qtype == "hole"){ std::swap(pind, qind); std::swap(ptype, qtype); TBME *= -1.0; }
		  if(rtype == "particle" && stype == "hole"){ std::swap(rind, sind); std::swap(rtype, stype); TBME *= -1.0; }
		  if((ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "hole") || 
		     (ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "hole") ||
		     (ptype == "particle" && qtype == "particle" && rtype == "hole" && stype == "particle")){
		    std::swap(pind, rind);
		    std::swap(qind, sind);
		    std::swap(ptype, rtype);
		    std::swap(qtype, stype);
		  }
		  if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.hhvec1[ind1], Chan.hhvec1[ind1], Chan.hh[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V2[ind1][ind] = TBME;
		  }
		  else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2); 
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.ppvec1[ind1], Chan.pp[ind1], Chan.pp[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V1[ind1][ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell1, shell4, shell3, shell2);
		    Ints.D_ME1.V3[ind1][ind] = TBME;
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp2vec1[ind1], Chan.hp2[ind1], Chan.hp2[ind1], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V3[ind1][ind] = TBME;
		  }
		  else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
		    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell1, shell2);
		    Ints.D_ME1.V4[ind1][ind] = TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell1, shell2);
		    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell3, shell4, shell2, shell1);
		    Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
		    ind = Index22(Chan.ppvec1[ind1], Chan.hhvec1[ind1], Chan.pp[ind1], Chan.hh[ind1], shell4, shell3, shell2, shell1);
		    Ints.D_ME1.V4[ind1][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell3];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		    Ints.D_ME1.V5[ind2][ind] = TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell4];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
		    Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
		    Ints.D_ME1.V5[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell4];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell1, shell2, shell3);
		    Ints.D_ME1.V6[ind2][ind] = TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell4, shell2, shell1, shell3);
		    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell3];
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		    Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.pvec1[ind2], Chan.hhpvec1[ind2], Chan.p[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		    Ints.D_ME1.V6[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell1];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V7[ind2][ind] = TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell2];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V7[ind2][ind] = TBME;
		    
		    ind2 = Chan.indvec[shell2];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		    Ints.D_ME1.V8[ind2][ind] = TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
		    ind2 = Chan.indvec[shell1];
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell3, shell4);
		    Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
		    ind = Index13(Chan.hvec1[ind2], Chan.hppvec1[ind2], Chan.h[ind2], Chan.hpp[ind2], shell1, shell2, shell4, shell3);
		    Ints.D_ME1.V8[ind2][ind] = TBME;
		    
		    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
		    Ints.D_ME1.V9[ind1][ind] = TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
		    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
		    Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
		    Ints.D_ME1.V9[ind1][ind] = TBME;
		    
		    minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell4, shell2, shell3);
		    Ints.D_ME1.V10[ind1][ind] = TBME;
		    minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell4, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell4, shell1, shell3);
		    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell1);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell1, shell3, shell2, shell4);
		    Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
		    minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
		    ind1 = ChanInd_2b(Parameters.basis, Space, tb);
		    //ind1 = ChanInd_2b_cross(Parameters.basis, Space, shell3, shell2);
		    ind = Index22(Chan.hp2vec1[ind1], Chan.hp1vec1[ind1], Chan.hp2[ind1], Chan.hp1[ind1], shell2, shell3, shell1, shell4);
		    Ints.D_ME1.V10[ind1][ind] = TBME;
		  }
		  if(Parameters.approx == "singles"){
		    if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
		      ind2 = Chan.indvec[shell2];
		      ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell3, shell4);
		      Ints.S_ME1.HPPP[ind2][ind] = TBME;
		      ind = Index13(Chan.pvec1[ind2], Chan.hppvec1[ind2], Chan.p[ind2], Chan.hpp[ind2], shell2, shell1, shell4, shell3);
		      Ints.S_ME1.HPPP[ind2][ind] = -1.0 * TBME;
		      
		      ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell3, shell4, shell2);
		      Ints.S_ME1.HPPP2[ind2][ind] = TBME;
		      ind = Index31(Chan.hppvec1[ind2], Chan.pvec1[ind2], Chan.hpp[ind2], Chan.p[ind2], shell1, shell4, shell3, shell2);
		      Ints.S_ME1.HPPP2[ind2][ind] = -1.0 * TBME;
		      
		      ind2 = Chan.indvec[shell3];
		      ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell4, shell3);
		      Ints.S_ME1.V13[ind2][ind] = TBME;
		      ind2 = Chan.indvec[shell4];
		      ind = Index31(Chan.hpp2vec1[ind2], Chan.pvec1[ind2], Chan.hpp2[ind2], Chan.p[ind2], shell1, shell2, shell3, shell4);
		      Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
		      
		      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		      ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell3, shell4, shell1, shell2);
		      Ints.S_ME1.HPPP4[ind2][ind] = TBME;
		      ind = Index22(Chan.ppvec1[ind2], Chan.hpvec1[ind2], Chan.pp[ind2], Chan.hp[ind2], shell4, shell3, shell1, shell2);
		      Ints.S_ME1.HPPP4[ind2][ind] = -1.0 * TBME;

		      minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell3);
		      ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell4, shell2, shell1, shell3);
		      Ints.S_ME1.HPPP5[ind2][ind] = TBME;
		      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell4);
		      ind = Index22(Chan.pp1vec1[ind2], Chan.hp1vec1[ind2], Chan.pp1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HPPP5[ind2][ind] = -1.0 * TBME;
		    }
		    else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
		      ind2 = Chan.indvec[shell3];
		      ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell1, shell2, shell4);
		      Ints.S_ME1.HHHP[ind2][ind] = TBME;
		      ind = Index13(Chan.hvec1[ind2], Chan.hhpvec1[ind2], Chan.h[ind2], Chan.hhp[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HHHP[ind2][ind] = -1.0 * TBME;
		      
		      ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell1, shell2, shell4, shell3);
		      Ints.S_ME1.HHHP2[ind2][ind] = TBME;
		      ind = Index31(Chan.hhpvec1[ind2], Chan.hvec1[ind2], Chan.hhp[ind2], Chan.h[ind2], shell2, shell1, shell4, shell3);
		      Ints.S_ME1.HHHP2[ind2][ind] = -1.0 * TBME;
		      
		      ind2 = Chan.indvec[shell1];
		      ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell2, shell3, shell4, shell1);
		      Ints.S_ME1.V14[ind2][ind] = TBME;
		      ind2 = Chan.indvec[shell2];
		      ind = Index31(Chan.hhp2vec1[ind2], Chan.hvec1[ind2], Chan.hhp2[ind2], Chan.h[ind2], shell1, shell3, shell4, shell2);
		      Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;

		      plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);		      
		      //ind2 = ChanInd_2b_dir(Parameters.basis, Space, shell1, shell2);
		      ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell1, shell2);
		      Ints.S_ME1.HHHP4[ind2][ind] = TBME;
		      ind = Index22(Chan.hpvec1[ind2], Chan.hhvec1[ind2], Chan.hp[ind2], Chan.hh[ind2], shell3, shell4, shell2, shell1);
		      Ints.S_ME1.HHHP4[ind2][ind] = -1.0 * TBME;
		      
		      minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell1, shell4);
		      ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell2, shell1, shell4);
		      Ints.S_ME1.HHHP5[ind2][ind] = TBME;
		      minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
		      ind2 = ChanInd_2b(Parameters.basis, Space, tb);
		      //ind2 = ChanInd_2b_cross(Parameters.basis, Space, shell2, shell4);
		      ind = Index22(Chan.hh1vec1[ind2], Chan.hp1vec1[ind2], Chan.hh1[ind2], Chan.hp1[ind2], shell3, shell1, shell2, shell4);
		      Ints.S_ME1.HHHP5[ind2][ind] = -1.0 * TBME;
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
}*/
