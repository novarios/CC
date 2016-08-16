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
      V[i] = new double[ntb1 * ntb1];
      for(int j = 0; j < ntb1*ntb1; ++j){
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

//Function to setup Channels
HF_Channels::HF_Channels(const Input_Parameters &Parameters, const Model_Space &Space)
{
  std::cout << "Building HF Channels ... " << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;

  State state;
  int count0, ind1, nob1, ntb1;
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
      if(k == count0 - 1){
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
	++nob[k];
      }
    }
  }

  /*for(int i = 0; i < size3; ++i){
    nob1 = obnum[i];
    std::cout << "Chan3 = " << i << " : " << qnums3[i].ml << " " << qnums3[i].m << " " << qnums3[i].ml << ", " << nob1 << std::endl;
    if(nob1 != 0){
      std::cout << "obvec = ";
      for(int j = 0; j < nob1; ++j){ std::cout << obvec[i][j] << " "; }
      std::cout << std::endl;
    }
    }*/

  size1 = Space.size_2b;
  qnums1 = new State[size1];
  tbvec = new int*[size1];
  ntb = new int[size1];

  for(int i = 0; i < size1; ++i){
    ntb[i] = 0;
  }

  for(int p = 0; p < Space.indtot; ++p){
    for(int q = 0; q < Space.indtot; ++q){
      plus(state, Space.qnums[p], Space.qnums[q]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      qnums1[ind1] = state;
      ++ntb[ind1];
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
      plus(state, Space.qnums[p], Space.qnums[q]);
      ind1 = ChanInd_2b_dir(Parameters.basis, Space, state);
      tbvec[ind1][2 * ntb[ind1]] = p;
      tbvec[ind1][2 * ntb[ind1] + 1] = q;
      ++ntb[ind1];
    }
  }

  /*for(int i = 0; i < size1; ++i){
    std::cout << "Chan1 = " << i << " : " << qnums1[i].ml << " " << qnums1[i].m << " " << qnums1[i].ml << std::endl;
    if(ntb1 != 0){
      std::cout << "tbvec = ";
      for(int j = 0; j < ntb1; ++j){ std::cout << tbvec[i][2*j] << tbvec[i][2*j + 1] << " "; }
      std::cout << std::endl;
    }
    }*/
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
  delete[] nob;

  for(int i = 0; i < size1; ++i){
    ntb1 = ntb[i];
    if(ntb1 != 0){
      delete[] tbvec[i];
    }
  }
  delete[] qnums1;
  delete[] tbvec;
  delete[] ntb;
}

//ignores existing holes and particles and fills them from protons and neutrons
void Single_Particle_States::Separate(const HF_Channels &Chan)
{
  //define holes/particles
  double tempen;
  int ind, h1, p1, nob1;

  for(int i = 0; i < Chan.size3; ++i){
    h1 = h[i];
    p1 = p[i];
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

  for(int i = 0; i < Chan.size3; ++i){
    h[i] = 0;
    p[i] = 0;
    for(int j = 0; j < Chan.nob[i]; ++j){
      tempen = energies[i][j];
      ind = 0;
      for(int k = 0; k < Chan.size3; ++k){
	if(Chan.qnums3[i].t == Chan.qnums3[k].t){
	  for(int l = 0; l < Chan.nob[k]; ++l){
	    if(energies[k][l] <= tempen){ ++ind; }
	  }
	}
      }
      if(Chan.qnums3[i].t == -1){
	if(ind <= hp){ ++h[i]; }
	else{ ++p[i]; }
      }
      else if(Chan.qnums3[i].t == 1){
	if(ind <= hn){ ++h[i]; }
	else{ ++p[i]; }
      }
    }
  }

  for(int i = 0; i < Chan.size3; ++i){
    h1 = h[i];
    p1 = p[i];
    nob1 = Chan.nob[i];
    if(h1 != 0){
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

  /*for(int i = 0; i < Chan.size3; ++i){
    nob1 = Chan.nob[i];
    h1 = h[i];
    p1 = p[i];
    std::cout << "Chan3 = " << i << " : " << nob1 << " " << h1 << " " << p1 << std::endl;
    if(h1 != 0){
      for(int j = 0; j < h1; ++j){
	std::cout << "h: " << Chan.obvec[i][j] << " : " << h_energies[i][j] << std::endl;
      }
    }
    if(p1 != 0){
      for(int j = 0; j < p1; ++j){
	std::cout << "p: " << Chan.obvec[i][j+h1] << " : " << pt_energies[i][j] << std::endl;
      }
    }
    }*/

}

Single_Particle_States::Single_Particle_States(const Input_Parameters &Parameters, const Model_Space &Space, const HF_Channels &Chan)
{
  std::cout << "Building Single-Particle States" << std::endl;
  std::cout << "-------------------------------" << std::endl;

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

  /*for(int i = 0; i < Chan.size3; ++i){
    std::cout << "Chan3 = " << i << " : " << h[i] << " " << p[i] << std::endl;
    for(int j = 0; j < Chan.nob[i]; ++j){
      std::cout << Chan.obvec[i][j] << " : " << energies[i][j] << std::endl;
    }
  }
  std::cout << std::endl;*/

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
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell3, shell4);
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

    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell3, shell4);
    ME.V[ind1][ind] = TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell3, shell4);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell1, shell2, shell4, shell3);
    ME.V[ind1][ind] = -1.0 * TBME;
    ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell2, shell1, shell4, shell3);
    ME.V[ind1][ind] = TBME;
    if(shell1 != shell3 || shell2 != shell4){
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell1, shell2);
      ME.V[ind1][ind] = TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell3, shell4, shell2, shell1);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell1, shell2);
      ME.V[ind1][ind] = -1.0 * TBME;
      ind = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], shell4, shell3, shell2, shell1);
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
    for(int pq = 0; pq < Chan.ntb[i]; ++pq){
      pind = Chan.tbvec[i][2*pq];
      qind = Chan.tbvec[i][2*pq + 1];
      if(pind == qind){ continue; }
      for(int rs = 0; rs < Chan.ntb[i]; ++rs){
	rind = Chan.tbvec[i][2*rs];
	sind = Chan.tbvec[i][2*rs + 1];
	if(rind == sind){ continue; }
	TBME = Coulomb_HO(Parameters, Space, pind, qind, rind, sind);
	ME.V[i][pq * Chan.ntb[i] + rs] = TBME;
	/*if(pind < qind && rind < sind && pind <= rind && pq <= rs){
	  std::cout << pind << " " << qind << " " << rind << " " << sind << " : " << TBME << std::endl;
	  }*/
      }
    }
  }
}

void Hartree_Fock_States_J(const Input_Parameters &Parameters, Model_Space &Space, const HF_Channels &Chan, Single_Particle_States &HF, const HF_Matrix_Elements &ME)
{
  double error; // Energy error between iterations
  double Bshift; // level shift parameter
  int ind; // Index to keep track of iteration number
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  int minj, maxj;
  State tb;

  double *w;
  double *work;

  //double *densmat;
  double *fock;

  jobz = 'V';
  uplo = 'U';
  Bshift = 50.0;

  ind = 0;
  error = 1000;
  while((error > 1e-8 && ind < 100) || ind < 10){
    ++ind;
    error = 0.0;
    
    //Make Fock Matrix
    for(int i = 0; i < Chan.size3; ++i){
      int size1 = Chan.nob[i];
      fock = new double[size1 * size1];
      for(int m = 0; m < size1; ++m){
	int mind = Chan.obvec[i][m];
	for(int l = 0; l < size1; ++l){
	  int lind = Chan.obvec[i][l];
	  if(m == l){ fock[size1 * m + m] = Space.qnums[mind].energy; }	// Add diagonal elements to fock matrices
	  else{ fock[size1 * m + l] = 0.0; }
	  for(int j = 0; j < Chan.size3; ++j){
	    int size2 = Chan.nob[j];
	    minj = abs(Chan.qnums3[i].j - Chan.qnums3[j].j);
	    maxj = Chan.qnums3[i].j + Chan.qnums3[j].j;
	    for(int beta = 0; beta < HF.h[j]; ++beta){ // Sum over occupied levels
	      for(int J = minj; J <= maxj; ++J){
		for(int n = 0; n < size2; ++n){
		  int nind = Chan.obvec[j][n];
		  for(int k = 0; k < size2; ++k){
		    int kind = Chan.obvec[j][k];
		    plus(tb, Space.qnums[mind], Space.qnums[nind]);
		    tb.j = J;
		    int ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		    int ind2 = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], mind, nind, lind, kind);
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
      for(int j = 0; j < size1*size1; ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      
      lda = size1;
      lwork = (3+2)*size1;
      w = new double[lda];
      work = new double[lwork];
      for(int j = 0; j < size1; ++j){ w[j] = 0.0; }
      for(int j = 0; j < (3+2)*size1; ++j){ work[j] = 0.0; }
    
      if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, fock, &lda, w, work, &lwork, &info); }
      for(int j = 0; j < size1*size1; ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      for(int j = 0; j < HF.h[i]; ++j){ w[j] += Bshift; } //Add back Level-shift parameter
      
      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  HF.vectors[i][j][k] = fock[size1 * j + k];
	}
      }
      
      delete[] fock;
      delete[] w;
      delete[] work;

      int ind2;
      double tempen2;
      double tempen3;
      double vec3;
      // Order states by energy
      for(int j = 0; j < size1 - 1; ++j){
	ind2 = j;
	tempen2 = w[j];
	for(int k = j + 1; k < size1; ++k){
	  if(w[k] < tempen2){ tempen2 = w[k]; ind2 = k; }
	}
	tempen3 = w[j];
	w[j] = w[ind2];
	w[ind2] = tempen3;
	for(int k = 0; k < size1; ++k){
	  vec3 = HF.vectors[i][j][k];
	  HF.vectors[i][j][k] = HF.vectors[i][ind2][k];
	  HF.vectors[i][ind2][k] = vec3;
	}
      }
      
      for(int j = 0; j < size1; ++j){
	error += fabs(HF.energies[i][j] - w[j])/Chan.nob[i]; 
	HF.energies[i][j] = w[j];
      }
      GramSchmidt(HF.vectors[i], size1);
      HF.Separate(Chan);
    }
  }

  ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < Chan.nob[i]; ++j){
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
  Single_Particle_States HFtemp(Parameters, Space, Chan);
  double error; // Energy error between iterations
  double Bshift; // level shift parameter
  int ind; // Index to keep track of iteration number
  char jobz, uplo; // Parameters for Diagonalization, Multiplication
  int lda; // Parameter for Diagonalization
  int lwork, info; // Parameters for Diagonaliztion
  double term;
  State tb;
  int nob1;
  int ind2;
  double tempen2;
  double tempen3;
  double vec3;

  double *w;
  double *work;

  double *fock;

  jobz = 'V';
  uplo = 'U';
  Bshift = 50.0;

  ind = 0;
  error = 1000;

  /*for(int i = 0; i < Chan.size3; ++i){
    int size1 = Chan.nob[i];
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
    ++ind;
    error = 0.0;
    
    //Make Fock Matrix
    for(int i = 0; i < Chan.size3; ++i){
      int size1 = Chan.nob[i];
      //std::cout << "chan = " << i << ", size = " << size1 << ", h = " << HFtemp.h[i] << ", p = " << HFtemp.p[i] << std::endl;
      fock = new double[size1 * size1];
      for(int m = 0; m < size1; ++m){
	int mind = Chan.obvec[i][m];
	for(int l = 0; l < size1; ++l){
	  int lind = Chan.obvec[i][l];
	  if(m == l){ fock[size1 * m + m] = Space.qnums[mind].energy; }	// Add diagonal elements to fock matrices
	  else{ fock[size1 * m + l] = 0.0; }
	  for(int j = 0; j < Chan.size3; ++j){
	    int size2 = Chan.nob[j];
	    for(int beta = 0; beta < HF.h[j]; ++beta){ // Sum over occupied levels
	      for(int n = 0; n < size2; ++n){
		int nind = Chan.obvec[j][n];
		for(int k = 0; k < size2; ++k){
		  int kind = Chan.obvec[j][k];
		  plus(tb, Space.qnums[mind], Space.qnums[nind]);
		  int ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  int ind2 = Index22(Chan.tbvec[ind1], Chan.tbvec[ind1], Chan.ntb[ind1], Chan.ntb[ind1], mind, nind, lind, kind);
		  term = HFtemp.vectors[j][beta][n] * HFtemp.vectors[j][beta][k] * ME.V[ind1][ind2];
		  fock[size1 * m + l] += term;
		}
	      }
	    }
	  }
	  for(int beta = 0; beta < HFtemp.h[i]; ++beta){ // Sum over occupied levels
	    fock[size1 * m + l] -= HFtemp.vectors[i][beta][m] * HFtemp.vectors[i][beta][l] * Bshift;
	  }
	}
      }
      for(int j = 0; j < size1*size1; ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}

      /*for(int x = 0; x < size1; ++x){
	for(int y = 0; y < size1; ++y){
	  std::cout << fock[size1*x + y] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/

      lda = size1;
      lwork = (3+2)*size1;
      w = new double[lda];
      work = new double[lwork];
      for(int j = 0; j < size1; ++j){ w[j] = 0.0; }
      for(int j = 0; j < (3+2)*size1; ++j){ work[j] = 0.0; }
    
      if(size1 != 0){ dsyev_(&jobz, &uplo, &size1, fock, &lda, w, work, &lwork, &info); }
      for(int j = 0; j < size1*size1; ++j){ if(fabs(fock[j]) < 10e-13){ fock[j] = 0.0; }}
      for(int j = 0; j < HFtemp.h[i]; ++j){ w[j] += Bshift; } //Add back Level-shift parameter
      
      for(int j = 0; j < size1; ++j){
	for(int k = 0; k < size1; ++k){
	  HF.vectors[i][j][k] = fock[size1 * j + k];
	}
      }

      /*for(int x = 0; x < size1; ++x){
	std::cout << "E = " << w[x] << std::endl;
	for(int y = 0; y < size1; ++y){
	  std::cout << fock[size1*x + y] << " ";
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;*/

      // Order states by energy
      for(int j = 0; j < size1 - 1; ++j){
	ind2 = j;
	tempen2 = w[j];
	for(int k = j + 1; k < size1; ++k){
	  if(w[k] < tempen2){ tempen2 = w[k]; ind2 = k; }
	}
	tempen3 = w[j];
	w[j] = w[ind2];
	w[ind2] = tempen3;
	for(int k = 0; k < size1; ++k){
	  vec3 = HF.vectors[i][j][k];
	  HF.vectors[i][j][k] = HF.vectors[i][ind2][k];
	  HF.vectors[i][ind2][k] = vec3;
	}
      }

      for(int j = 0; j < size1; ++j){
	error += fabs(HF.energies[i][j] - w[j])/Chan.nob[i]; 
	HF.energies[i][j] = w[j];
      }
      GramSchmidt(HF.vectors[i], size1);

      delete[] fock;
      delete[] w;
      delete[] work;
    }
    HF.Separate(Chan);

    // HFtemp = HF
    for(int i = 0; i < Chan.size3; ++i){
      nob1 = Chan.nob[i];
      HFtemp.h[i] = HF.h[i];
      HFtemp.p[i] = HF.p[i];
      if(nob1 != 0){
	for(int j = 0; j < nob1; ++j){
	  HFtemp.energies[i][j] = HF.energies[i][j];
	  for(int k = 0; k < nob1; ++k){
	    HFtemp.vectors[i][j][k] = HF.vectors[i][j][k];
	  }
	}
      }
    }
  }

  HFtemp.delete_struct(Chan);

  ind = 0;
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = 0; j < HF.h[i]; ++j){
      Space.qnums[ind] = Chan.qnums3[i];
      Space.qnums[ind].energy = HF.energies[i][j];
      Space.qnums[ind].type = "hole";
      ++ind;
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    for(int j = HF.h[i]; j < HF.h[i] + HF.p[i]; ++j){
      Space.qnums[ind] = Chan.qnums3[i];
      Space.qnums[ind].energy = HF.energies[i][j];
      Space.qnums[ind].type = "particle";
      ++ind;
    }
  }

  // Order holes by energy
  for(int i = 0; i < Space.indhol; ++i){
    ind2 = i;
    tempen2 = Space.qnums[i].energy;
    for(int j = i + 1; j < Space.indhol; ++j){
      if((tempen2 - Space.qnums[j].energy) > 1e-6){
	tempen2 = Space.qnums[j].energy;
	ind2 = j;
      }
    }
    tb = Space.qnums[ind2];
    for(int j = ind2; j > i; --j){
      Space.qnums[j] = Space.qnums[j - 1];
    }
    Space.qnums[i] = tb;
  }
  // Order particles by energy
  for(int i = Space.indhol; i < Space.indtot; ++i){
    ind2 = i;
    tempen2 = Space.qnums[i].energy;
    for(int j = i + 1; j < Space.indtot; ++j){
      if((tempen2 - Space.qnums[j].energy) > 1e-6){
	tempen2 = Space.qnums[j].energy;
	ind2 = j;
      }
    }
    tb = Space.qnums[ind2];
    for(int j = ind2; j > i; --j){
      Space.qnums[j] = Space.qnums[j - 1];
    }
    Space.qnums[i] = tb;
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
  double *M1, *M2; // Matrices of coefficients
  double *C;

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
    length = Chan.ntb[chan];
    if(length == 0){ continue; }
    matlength = pow(length, 2.0);

    M1 = new double[matlength];
    M2 = new double[matlength];
    C = new double[matlength];
    for(int i = 0; i < matlength; ++i){
      M1[i] = 0.0;
      M2[i] = 0.0;
      C[i] = 0.0;
    }

    for(int pq = 0; pq < length; ++pq){
      for(int ag = 0; ag < length; ++ag){
	int p = Chan.tbvec[chan][2*pq];
	int q = Chan.tbvec[chan][2*pq + 1];
	int a = Chan.tbvec[chan][2*ag];
	int g = Chan.tbvec[chan][2*ag + 1];
	//td::cout << "!! " << p << " " << q << " " << a << " " << g << std::endl;
	if(Chan.indvec[p] != Chan.indvec[a] || Chan.indvec[q] != Chan.indvec[g]){ continue; }
	int ind1 = Chan.indvec[p];
	int ind2 = Chan.indvec[q];
	//std::cout << "pass! " << ind1 << " " << ind2 << std::endl;
	p = Index1(Chan.obvec[ind1], Chan.nob[ind1], p);
	q = Index1(Chan.obvec[ind2], Chan.nob[ind2], q);
	a = Index1(Chan.obvec[ind1], Chan.nob[ind1], a);
	g = Index1(Chan.obvec[ind2], Chan.nob[ind2], g);
	//std::cout << "?? " << p << " " << a << " = " << States.vectors[ind1][p][a] << std::endl;
	//std::cout << "?? " << q << " " << g << " = " << States.vectors[ind2][q][g] << std::endl;
	tempel = States.vectors[ind1][p][a] * States.vectors[ind2][q][g];
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
    
    RM_dgemm(M1, ME.V[chan], C, &length, &length, &length, &alpha1, &beta1, &transa, &transa);
    /*std::cout << "V1" << std::endl;
      for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
      std::cout << C[length*i + j] << " ";
      }
      std::cout << std::endl;
      }
      std::cout << std::endl;*/
    RM_dgemm(C, M2, ME.V[chan], &length, &length, &length, &alpha1, &beta1, &transa, &transa);
    /*std::cout << "V2" << std::endl;
      for(int i = 0; i < length; ++i){
      for(int j = 0; j < length; ++j){
      std::cout << ME.V[chan][length*i + j] << " ";
      }
      std::cout << std::endl;
      }
      std::cout << std::endl;*/

    delete[] M1;
    delete[] M2;
    delete[] C;
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
    for(int tb1 = 0; tb1 < HF_Chan.ntb[chan]; ++tb1){
      shell1 = HF_Chan.tbvec[chan][2*tb1];
      shell2 = HF_Chan.tbvec[chan][2*tb1 + 1];
      ptype = Space.qnums[shell1].type;
      qtype = Space.qnums[shell2].type;
      //std::cout << "shell1, shell2 = " << shell1 << " " << shell2 << std::endl;
      if(ptype == qtype && shell1 >= shell2){ continue; }
      if(ptype == "particle" && qtype == "hole"){ continue; }
      for(int tb2 = 0; tb2 < HF_Chan.ntb[chan]; ++tb2){
	shell3 = HF_Chan.tbvec[chan][2*tb2];
	shell4 = HF_Chan.tbvec[chan][2*tb2 + 1];
	rtype = Space.qnums[shell3].type;
	stype = Space.qnums[shell4].type;
	//std::cout << "shell3, shell4 = " << shell3 << " " << shell4 << std::endl;
	if(rtype == stype && shell3 >= shell4){ continue; }
	if(rtype == "particle" && stype == "hole"){ continue; }
	if(ptype == rtype && shell1 > shell3){ continue; }
	if(ptype == rtype && qtype == stype && shell1 == shell3 && shell2 > shell4){ continue; }
	if(ptype == "particle" && rtype == "hole"){ continue; }
	//std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << std::endl;
	TBME = HF_ME.V[chan][tb1*HF_Chan.ntb[chan] + tb2];
	//std::cout << shell1 << " " << shell2 << " " << shell3 << " " << shell4 << ", " << TBME << std::endl;
	if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "hole"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V2[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.hhvec[ind1], Chan.hhvec[ind1], Chan.nhh[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V2[ind1][ind] = TBME;
	}
	else if(ptype == "particle" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V1[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.ppvec[ind1], Chan.npp[ind1], Chan.npp[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V1[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "particle" && rtype == "hole" && stype == "particle"){
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell1, shell4, shell3, shell2);
	  Ints.D_ME1.V3[ind1][ind] = TBME;
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp2vec[ind1], Chan.nhp2[ind1], Chan.nhp2[ind1], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V3[ind1][ind] = TBME;
	}
	else if(ptype == "hole" && qtype == "hole" && rtype == "particle" && stype == "particle"){
	  plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell1, shell2);
	  Ints.D_ME1.V4[ind1][ind] = TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell1, shell2);
	  Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell3, shell4, shell2, shell1);
	  Ints.D_ME1.V4[ind1][ind] = -1.0 * TBME;
	  ind = Index22(Chan.ppvec[ind1], Chan.hhvec[ind1], Chan.npp[ind1], Chan.nhh[ind1], shell4, shell3, shell2, shell1);
	  Ints.D_ME1.V4[ind1][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell2];
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V5[ind2][ind] = TBME;
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell1];
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V5[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V5[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell1];
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell3, shell4);
	  Ints.D_ME1.V6[ind2][ind] = TBME;
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell1, shell2, shell4, shell3);
	  Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell2];
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	  Ints.D_ME1.V6[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.hvec[ind2], Chan.hppvec[ind2], Chan.nh[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	  Ints.D_ME1.V6[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell4];
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	  Ints.D_ME1.V7[ind2][ind] = TBME;
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	  Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell3];
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	  Ints.D_ME1.V7[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V7[ind2][ind] = TBME;
	  
	  ind2 = Chan.indvec[shell3];
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	  Ints.D_ME1.V8[ind2][ind] = TBME;
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	  Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	  ind2 = Chan.indvec[shell4];
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell1, shell2, shell3);
	  Ints.D_ME1.V8[ind2][ind] = -1.0 * TBME;
	  ind = Index13(Chan.pvec[ind2], Chan.hhpvec[ind2], Chan.np[ind2], Chan.nhhp[ind2], shell4, shell2, shell1, shell3);
	  Ints.D_ME1.V8[ind2][ind] = TBME;
	  
	  minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	  Ints.D_ME1.V9[ind1][ind] = TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	  Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	  Ints.D_ME1.V9[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	  Ints.D_ME1.V9[ind1][ind] = TBME;
	  
	  minus(tb, Space.qnums[shell4], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell4, shell2, shell3);
	  Ints.D_ME1.V10[ind1][ind] = TBME;
	  minus(tb, Space.qnums[shell4], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell4, shell1, shell3);
	  Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell1]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell1, shell3, shell2, shell4);
	  Ints.D_ME1.V10[ind1][ind] = -1.0 * TBME;
	  minus(tb, Space.qnums[shell3], Space.qnums[shell2]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  ind = Index22(Chan.hp2vec[ind1], Chan.hp1vec[ind1], Chan.nhp2[ind1], Chan.nhp1[ind1], shell2, shell3, shell1, shell4);
	  Ints.D_ME1.V10[ind1][ind] = TBME;
	}
	if(Parameters.approx == "singles"){
	  if(ptype == "hole" && qtype == "particle" && rtype == "particle" && stype == "particle"){
	    ind2 = Chan.indvec[shell2];
	    ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell3, shell4);
	    Ints.S_ME1.V11[ind2][ind] = TBME;
	    ind = Index13(Chan.pvec[ind2], Chan.hppvec[ind2], Chan.np[ind2], Chan.nhpp[ind2], shell2, shell1, shell4, shell3);
	    Ints.S_ME1.V11[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell4, shell3);
	    Ints.S_ME1.V13[ind2][ind] = TBME;
	    ind2 = Chan.indvec[shell4];
	    ind = Index31(Chan.hpp1vec[ind2], Chan.pvec[ind2], Chan.nhpp1[ind2], Chan.np[ind2], shell1, shell2, shell3, shell4);
	    Ints.S_ME1.V13[ind2][ind] = -1.0 * TBME;
	    
	    minus(tb, Space.qnums[shell1], Space.qnums[shell3]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell4, shell2, shell1, shell3);
	    Ints.S_ME1.V16[ind2][ind] = TBME;
	    minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.pp1vec[ind2], Chan.hp1vec[ind2], Chan.npp1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V16[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell2];
	    ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell3, shell4, shell2);
	    Ints.S_ME1.V17[ind2][ind] = TBME;
	    ind = Index31(Chan.hppvec[ind2], Chan.pvec[ind2], Chan.nhpp[ind2], Chan.np[ind2], shell1, shell4, shell3, shell2);
	    Ints.S_ME1.V17[ind2][ind] = -1.0 * TBME;
	    
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell3, shell4, shell1, shell2);
	    Ints.S_ME1.V20[ind2][ind] = TBME;
	    ind = Index22(Chan.ppvec[ind2], Chan.hpvec[ind2], Chan.npp[ind2], Chan.nhp[ind2], shell4, shell3, shell1, shell2);
	    Ints.S_ME1.V20[ind2][ind] = -1.0 * TBME;
	  }
	  else if(ptype == "hole" && qtype == "hole" && rtype == "hole" && stype == "particle"){
	    ind2 = Chan.indvec[shell3];
	    ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell1, shell2, shell4);
	    Ints.S_ME1.V12[ind2][ind] = TBME;
	    ind = Index13(Chan.hvec[ind2], Chan.hhpvec[ind2], Chan.nh[ind2], Chan.nhhp[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V12[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell1];
	    ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell2, shell3, shell4, shell1);
	    Ints.S_ME1.V14[ind2][ind] = TBME;
	    ind2 = Chan.indvec[shell2];
	    ind = Index31(Chan.hhp1vec[ind2], Chan.hvec[ind2], Chan.nhhp1[ind2], Chan.nh[ind2], shell1, shell3, shell4, shell2);
	    Ints.S_ME1.V14[ind2][ind] = -1.0 * TBME;
	    
	    minus(tb, Space.qnums[shell1], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell2, shell1, shell4);
	    Ints.S_ME1.V15[ind2][ind] = TBME;
	    minus(tb, Space.qnums[shell2], Space.qnums[shell4]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hh1vec[ind2], Chan.hp1vec[ind2], Chan.nhh1[ind2], Chan.nhp1[ind2], shell3, shell1, shell2, shell4);
	    Ints.S_ME1.V15[ind2][ind] = -1.0 * TBME;
	    
	    ind2 = Chan.indvec[shell3];
	    ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell1, shell2, shell4, shell3);
	    Ints.S_ME1.V18[ind2][ind] = TBME;
	    ind = Index31(Chan.hhpvec[ind2], Chan.hvec[ind2], Chan.nhhp[ind2], Chan.nh[ind2], shell2, shell1, shell4, shell3);
	    Ints.S_ME1.V18[ind2][ind] = -1.0 * TBME;
	    
	    plus(tb, Space.qnums[shell1], Space.qnums[shell2]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell1, shell2);
	    Ints.S_ME1.V19[ind2][ind] = TBME;
	    ind = Index22(Chan.hpvec[ind2], Chan.hhvec[ind2], Chan.nhp[ind2], Chan.nhh[ind2], shell3, shell4, shell2, shell1);
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
