#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

Interactions::Interactions(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1 = Doubles_ME1(Chan);
  }
  else if(Parameters.approx == "singles"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1 = Doubles_ME1(Chan);
    S_ME1 = Singles_ME1(Chan);
  }
}

void Interactions::delete_struct(const Input_Parameters &Parameters, const Channels &Chan)
{
  if(Parameters.approx == "doubles"){
    D_ME1.delete_struct(Chan);;
  }
  else if(Parameters.approx == "singles"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
  else if(Parameters.approx == "triples"){
    D_ME1.delete_struct(Chan);
    S_ME1.delete_struct(Chan);
  }
}

Doubles_ME1::Doubles_ME1(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  V1 = new double*[Chan.size1];
  V2 = new double*[Chan.size1];
  V3 = new double*[Chan.size2];
  V4 = new double*[Chan.size1];
  V5 = new double*[Chan.size3];
  V6 = new double*[Chan.size3];
  V7 = new double*[Chan.size3];
  V8 = new double*[Chan.size3];
  V9 = new double*[Chan.size2];
  V10 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      V1[i] = new double[npp * npp];
      for(int j = 0; j < npp * npp; ++j){ V1[i][j] = 0.0; }
    }
    if(nhh != 0){
      V2[i] = new double[nhh * nhh];
      for(int j = 0; j < nhh * nhh; ++j){ V2[i][j] = 0.0; }
    }
    if(npp * nhh != 0){
      V4[i] = new double[npp * nhh];
      for(int j = 0; j < npp * nhh; ++j){ V4[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      V5[i] = new double[nh * nhpp];
      V6[i] = new double[nh * nhpp];
      for(int j = 0; j < nh * nhpp; ++j){ V5[i][j] = 0.0; V6[i][j] = 0.0; }
    }
    if(np * nhhp != 0){
      V7[i] = new double[np * nhhp];
      V8[i] = new double[np * nhhp];
      for(int j = 0; j < np * nhhp; ++j){ V7[i][j] = 0.0; V8[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      V3[i] = new double[nhp2 * nhp2];
      for(int j = 0; j < nhp2 * nhp2; ++j){ V3[i][j] = 0.0; }
    }
    if(nhp2 * nhp1){
      V9[i] = new double[nhp2 * nhp1];
      V10[i] = new double[nhp2 * nhp1];
      for(int j = 0; j < nhp2 * nhp1; ++j){ V9[i][j] = 0.0; V10[i][j] = 0.0; }
    }
  }
}

void Doubles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    if(npp != 0){
      delete[] V1[i];
    }
    if(nhh != 0){
      delete[] V2[i];
    }
    if(npp * nhh != 0){
      delete[] V4[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    if(nh * nhpp != 0){
      delete[] V5[i];
      delete[] V6[i];
    }
    if(np * nhhp != 0){
      delete[] V7[i];
      delete[] V8[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhp2 = Chan.nhp2[i];
    if(nhp2 != 0){
      delete[] V3[i];
    }
    if(nhp2 * nhp1){
      delete[] V9[i];
      delete[] V10[i];
    }
  }
  delete[] V1;
  delete[] V2;
  delete[] V3;
  delete[] V4;
  delete[] V5;
  delete[] V6;
  delete[] V7;
  delete[] V8;
  delete[] V9;
  delete[] V10;
}

Singles_ME1::Singles_ME1(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  V11 = new double*[Chan.size3];
  V12 = new double*[Chan.size3];
  V17 = new double*[Chan.size3];
  V18 = new double*[Chan.size3];
  V13 = new double*[Chan.size3];
  V14 = new double*[Chan.size3];
  V20 = new double*[Chan.size1];
  V19 = new double*[Chan.size1];
  V16 = new double*[Chan.size2];
  V15 = new double*[Chan.size2];
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      V20[i] = new double[npp * nhp];
      for(int j = 0; j < npp * nhp; ++j){ V20[i][j] = 0.0; }
    }
    if(nhp * nhh != 0){
      V19[i] = new double[nhp * nhh];
      for(int j = 0; j < nhp * nhh; ++j){ V19[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      V11[i] = new double[np * nhpp];
      V17[i] = new double[nhpp * np];
      for(int j = 0; j < np * nhpp; ++j){ V11[i][j] = 0.0; V17[i][j] = 0.0; }
    }
    if(nh * nhhp != 0){
      V12[i] = new double[nh * nhhp];
      V18[i] = new double[nhhp * nh];
      for(int j = 0; j < nh * nhhp; ++j){ V12[i][j] = 0.0; V18[i][j] = 0.0; }
    }
    if(nhpp1 * np != 0){
      V13[i] = new double[nhpp1 * np];
      for(int j = 0; j < nhpp1 * np; ++j){ V13[i][j] = 0.0; }
    }
    if(nhhp1 * nh != 0){
      V14[i] = new double[nhhp1 * nh];
      for(int j = 0; j < nhhp1 * nh; ++j){ V14[i][j] = 0.0; }
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      V16[i] = new double[npp1 * nhp1];
      for(int j = 0; j < npp1 * nhp1; ++j){ V16[i][j] = 0.0; }
    }
    if(nhh1 * nhp1 != 0){
      V15[i] = new double[nhh1 * nhp1];
      for(int j = 0; j < nhh1 * nhp1; ++j){ V15[i][j] = 0.0; }
    }
  }
}

void Singles_ME1::delete_struct(const Channels &Chan)
{
  int nhh, npp, nhp, nh, np, nhpp, nhhp, nhp1, nhh1, npp1, nhpp1, nhhp1;
  for(int i = 0; i < Chan.size1; ++i){
    nhh = Chan.nhh[i];
    npp = Chan.npp[i];
    nhp = Chan.nhp[i];
    if(npp * nhp != 0){
      delete[] V20[i];
    }
    if(nhp * nhh != 0){
      delete[] V19[i];
    }
  }
  for(int i = 0; i < Chan.size3; ++i){
    nh = Chan.nh[i];
    np = Chan.np[i];
    nhpp = Chan.nhpp[i];
    nhhp = Chan.nhhp[i];
    nhpp1 = Chan.nhpp1[i];
    nhhp1 = Chan.nhhp1[i];
    if(np * nhpp != 0){
      delete[] V11[i];
      delete[] V17[i];
    }
    if(nh * nhhp != 0){
      delete[] V12[i];
      delete[] V18[i];
    }
    if(nhpp1 * np != 0){
      delete[] V13[i];
    }
    if(nhhp1 * nh != 0){
      delete[] V14[i];
    }
  }
  for(int i = 0; i < Chan.size2; ++i){
    nhp1 = Chan.nhp1[i];
    nhh1 = Chan.nhh1[i];
    npp1 = Chan.npp1[i];
    if(npp1 * nhp1 != 0){
      delete[] V16[i];
    }
    if(nhh1 * nhp1 != 0){
      delete[] V15[i];
    }
  }
  delete[] V11;
  delete[] V12;
  delete[] V17;
  delete[] V18;
  delete[] V13;
  delete[] V14;
  delete[] V20;
  delete[] V19;
  delete[] V16;
  delete[] V15;
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

void Coulomb_Inf_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  double L = pow(Space.indhol/Parameters.density, 1.0/3.0);
  int i, j, k, l, a, b, c, d;
  int nhh, npp, nh, np, nhpp, nhhp, nhp1, nhp2;
  double TBME;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(npp == 0){ goto stop1; }
    #pragma omp parallel private(a, b, c, d, TBME)
    {
      #pragma omp for schedule(static)
      for(int cd = 0; cd < npp; ++cd){
	c = Chan.ppvec[chan][2*cd];
	d = Chan.ppvec[chan][2*cd + 1];
	if(c == d){ continue; }
	for(int ab = cd; ab < npp; ++ab){
	  a = Chan.ppvec[chan][2*ab];
	  b = Chan.ppvec[chan][2*ab + 1];
	  if(a == b){ continue; }
	  TBME = Coulomb_Inf(Space, a, b, c, d, L);
	  Ints.D_ME1.V1[chan][npp*cd + ab] = TBME;
	  Ints.D_ME1.V1[chan][npp*ab + cd] = TBME;
	}
      }
    }
  stop1:;
    if(nhh == 0){ goto stop2; }
    #pragma omp parallel private(i, j, k, l, TBME)
    {
      #pragma omp for schedule(static)
      for(int ij = 0; ij < nhh; ++ij){
	i = Chan.hhvec[chan][2*ij];
	j = Chan.hhvec[chan][2*ij + 1];
	if(i == j){ continue; }
	for(int kl = ij; kl < nhh; ++kl){
	  k = Chan.hhvec[chan][2*kl];
	  l = Chan.hhvec[chan][2*kl + 1];
	  if(k == l){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, k, l, L);
	  Ints.D_ME1.V2[chan][nhh*ij + kl] = TBME;
	  Ints.D_ME1.V2[chan][nhh*kl + ij] = TBME;
	}
      }
    }
  stop2:;
    if(nhh * npp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int ab = 0; ab < npp; ++ab){
	a = Chan.ppvec[chan][2*ab];
	b = Chan.ppvec[chan][2*ab + 1];
	if(a == b){ continue; }
	for(int ij = 0; ij < nhh; ++ij){
	  i = Chan.hhvec[chan][2*ij];
	  j = Chan.hhvec[chan][2*ij + 1];
	  if(i == j){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V4[chan][nhh*ab + ij] = TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size3; ++chan){
    nh = Chan.nh[chan];
    np = Chan.np[chan];
    nhpp = Chan.nhpp[chan];
    nhhp = Chan.nhhp[chan];
    if(nh * nhpp == 0){ goto stop3; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int h = 0; h < nh; ++h){
	j = Chan.hvec[chan][h];
	for(int hpp = 0; hpp < nhpp; ++hpp){
	  i = Chan.hppvec[chan][3*hpp];
	  a = Chan.hppvec[chan][3*hpp + 1];
	  b = Chan.hppvec[chan][3*hpp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V5[chan][nhpp*h + hpp] = TBME;
	  Ints.D_ME1.V6[chan][nhpp*h + hpp] = -1.0 * TBME;
	}
      }
    }
  stop3:;
    if(np * nhhp == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int p = 0; p < np; ++p){
	b = Chan.pvec[chan][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  i = Chan.hhpvec[chan][3*hhp];
	  j = Chan.hhpvec[chan][3*hhp + 1];
	  a = Chan.hhpvec[chan][3*hhp + 2];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V7[chan][nhhp*p + hhp] = TBME;
	  Ints.D_ME1.V8[chan][nhhp*p + hhp] = -1.0 * TBME;
	}
      }
    }
  }
  
  for(int chan = 0; chan < Chan.size2; ++chan){
    nhp1 = Chan.nhp1[chan];
    nhp2 = Chan.nhp2[chan];
    if(nhp2 == 0){ goto stop4; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp21 = 0; hp21 < nhp2; ++hp21){
	i = Chan.hp2vec[chan][2*hp21];
	b = Chan.hp2vec[chan][2*hp21 + 1];
	for(int hp22 = hp21; hp22 < nhp2; ++hp22){
	  j = Chan.hp2vec[chan][2*hp22];
	  a = Chan.hp2vec[chan][2*hp22 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, a, j, b, L);
	  Ints.D_ME1.V3[chan][nhp2*hp21 + hp22] = TBME;
	  Ints.D_ME1.V3[chan][nhp2*hp22 + hp21] = TBME;
	}
      }
    }
  stop4:;
    if(nhp1 * nhp2 == 0){ continue; }
    #pragma omp parallel private(i, j, a, b, TBME)
    {
      #pragma omp for schedule(static)
      for(int hp2 = 0; hp2 < nhp2; ++hp2){
	i = Chan.hp2vec[chan][2*hp2];
	a = Chan.hp2vec[chan][2*hp2 + 1];
	for(int hp1 = 0; hp1 < Chan.nhp1[chan]; ++hp1){
	  j = Chan.hp1vec[chan][2*hp1];
	  b = Chan.hp1vec[chan][2*hp1 + 1];
	  if(i == j || a == b){ continue; }
	  TBME = Coulomb_Inf(Space, i, j, a, b, L);
	  Ints.D_ME1.V9[chan][nhp1*hp2 + hp1] = TBME;
	  Ints.D_ME1.V10[chan][nhp1*hp2 + hp1] = -1.0 * TBME;
	}
      }
    }
  }
}

double Coulomb_Inf(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double term = 0.0;
  double e_sq = hbarc_HartA * fine_struct;
  double prefactor = e_sq/(L*L*L);
  double qSquared1;
  double kX1, kY1, kZ1;
  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx || 
     Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny || 
     Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qk].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qk].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qk].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-9){ term += 0.0; }
    else{ term += 4.0 * prefactor * M_PI / qSquared1; }
  }
  if(Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    kX1 = (2.0*M_PI/L) * (Space.qnums[qi].nx - Space.qnums[ql].nx);
    kY1 = (2.0*M_PI/L) * (Space.qnums[qi].ny - Space.qnums[ql].ny);
    kZ1 = (2.0*M_PI/L) * (Space.qnums[qi].nz - Space.qnums[ql].nz);
    qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
    if(qSquared1 < 1.0e-9){ term -= 0.0; }
    else{ term -= 4.0 * prefactor * M_PI / qSquared1; }
  }
  return term;
}


// Minnesota Potential for momentum basis
double Coulomb_HO(const Input_Parameters &Parameters, const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql)
{
  int g1, g2, g3, g4, G, L;
  double dir = 0.0;
  double exch = 0.0;
  int n1, m1, n2, m2, n3, m3, n4, m4;
  n1 = Space.qnums[qi].n; //1
  m1 = Space.qnums[qi].ml;
  n2 = Space.qnums[qj].n; //2
  m2 = Space.qnums[qj].ml;
  n3 = Space.qnums[ql].n; //4
  m3 = Space.qnums[ql].ml;
  n4 = Space.qnums[qk].n; //3
  m4 = Space.qnums[qk].ml;
  if((m1 + m2 == m3 + m4) && Space.qnums[qi].m == Space.qnums[qk].m && Space.qnums[qj].m == Space.qnums[ql].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = int(j1 + j4 + 0.5*(abs(m1) + m1) + 0.5*(abs(m4) - m4));
	    g2 = int(j2 + j3 + 0.5*(abs(m2) + m2) + 0.5*(abs(m3) - m3));
	    g3 = int(j3 + j2 + 0.5*(abs(m3) + m3) + 0.5*(abs(m2) - m2));
	    g4 = int(j4 + j1 + 0.5*(abs(m4) + m4) + 0.5*(abs(m1) - m1));
	    G = g1 + g2 + g3 + g4;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    double LogRatio2 = logratio2(G);
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    dir += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    dir *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }

  n1 = Space.qnums[qi].n; //1
  m1 = Space.qnums[qi].ml;
  n2 = Space.qnums[qj].n; //2
  m2 = Space.qnums[qj].ml;
  n3 = Space.qnums[qk].n; //3
  m3 = Space.qnums[qk].ml;
  n4 = Space.qnums[ql].n; //4
  m4 = Space.qnums[ql].ml;
  if((m1 + m2 == m3 + m4) && Space.qnums[qi].m == Space.qnums[ql].m && Space.qnums[qj].m == Space.qnums[qk].m){
    for(int j1 = 0; j1 <= n1; ++j1){
      for(int j2 = 0; j2 <= n2; ++j2){
	for(int j3 = 0; j3 <= n3; ++j3){
	  for(int j4 = 0; j4 <= n4; ++j4){
	    g1 = int(j1 + j4 + 0.5*(abs(m1) + m1) + 0.5*(abs(m4) - m4));
	    g2 = int(j2 + j3 + 0.5*(abs(m2) + m2) + 0.5*(abs(m3) - m3));
	    g3 = int(j3 + j2 + 0.5*(abs(m3) + m3) + 0.5*(abs(m2) - m2));
	    g4 = int(j4 + j1 + 0.5*(abs(m4) + m4) + 0.5*(abs(m1) - m1));
	    G = g1 + g2 + g3 + g4;
	    double LogRatio1 = logratio1(j1, j2, j3, j4);
	    double LogProd2 = logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4);
	    double LogRatio2 = logratio2(G);
	    double temp = 0.0;
	    for(int l1 = 0; l1 <= g1; ++l1){
	      for(int l2 = 0; l2 <= g2; ++l2){
		for(int l3 = 0; l3 <= g3; ++l3){
		  for(int l4 = 0; l4 <= g4; ++l4){
		    if(l1 + l2 != l3 + l4){ continue; }
		    L = l1 + l2 + l3 + l4;
		    temp += (-2*((g2 + g3 - l2 - l3)%2) + 1)
		      * exp(logproduct3(l1, l2, l3, l4, g1, g2, g3, g4) + lgamma(1.0 + 0.5*L) + lgamma(0.5*(G - L + 1.0)));
		  }
		}
	      }
	    }
	    exch += (-2*((j1 + j2 + j3 + j4)%2) + 1) * exp(LogRatio1 + LogProd2 + LogRatio2) * temp;
	  }
	}
      }
    }
    exch *= product1(n1, m1, n2, m2, n3, m3, n4, m4);
  }
  return std::sqrt(Parameters.density)*(dir - exch);
}

void Minnesota_Matrix_Elements(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints)
{
  std::cout << "Building Interaction Matrices ... " << std::endl;  
  double L = pow((Parameters.P + Parameters.N)/Parameters.density, 1./3.);
  for(int i = 0; i < Chan.size1; ++i){
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.npp[i]; ++rs){
	shell3 = Chan.ppvec[i][2*rs];
	shell4 = Chan.ppvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V1[i][Chan.npp[i]*pq + rs] = TBME;
	Ints.D_ME1.V1[i][Chan.npp[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.nhh[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hhvec[i][2*pq];
      shell2 = Chan.hhvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = pq; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V2[i][Chan.nhh[i]*pq + rs] = TBME;
	Ints.D_ME1.V2[i][Chan.nhh[i]*rs + pq] = TBME;
      }
    }
    for(int pq = 0; pq < Chan.npp[i]; ++pq){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.ppvec[i][2*pq];
      shell2 = Chan.ppvec[i][2*pq + 1];
      if(shell1 == shell2){ continue; }
      for(int rs = 0; rs < Chan.nhh[i]; ++rs){
	shell3 = Chan.hhvec[i][2*rs];
	shell4 = Chan.hhvec[i][2*rs + 1];
	if(shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V4[i][Chan.nhh[i]*pq + rs] = TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size3; ++i){
    for(int q = 0; q < Chan.nh[i]; ++q){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell2 = Chan.hvec[i][q];
      for(int prs = 0; prs < Chan.nhpp[i]; ++prs){
	shell1 = Chan.hppvec[i][3*prs];
	shell3 = Chan.hppvec[i][3*prs + 1];
	shell4 = Chan.hppvec[i][3*prs + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V5[i][Chan.nhpp[i]*q + prs] = TBME;
	Ints.D_ME1.V6[i][Chan.nhpp[i]*q + prs] = -1.0 * TBME;
      }
    }
    for(int s = 0; s < Chan.np[i]; ++s){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell4 = Chan.pvec[i][s];
      for(int pqr = 0; pqr < Chan.nhhp[i]; ++pqr){
	shell1 = Chan.hhpvec[i][3*pqr];
	shell2 = Chan.hhpvec[i][3*pqr + 1];
	shell3 = Chan.hhpvec[i][3*pqr + 2];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V7[i][Chan.nhhp[i]*s + pqr] = TBME;
	Ints.D_ME1.V8[i][Chan.nhhp[i]*s + pqr] = -1.0 * TBME;
      }
    }
  }
  
  for(int i = 0; i < Chan.size2; ++i){
    for(int ps = 0; ps < Chan.nhp2[i]; ++ps){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*ps];
      shell4 = Chan.hp2vec[i][2*ps + 1];
      for(int qr = ps; qr < Chan.nhp2[i]; ++qr){
	shell3 = Chan.hp2vec[i][2*qr];
	shell2 = Chan.hp2vec[i][2*qr + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell3, shell4, L);
	Ints.D_ME1.V3[i][Chan.nhp2[i]*ps + qr] = TBME;
	Ints.D_ME1.V3[i][Chan.nhp2[i]*qr + ps] = TBME;
      }
    }
    for(int pr = 0; pr < Chan.nhp2[i]; ++pr){
      double TBME;
      int shell1, shell2, shell3, shell4;
      shell1 = Chan.hp2vec[i][2*pr];
      shell3 = Chan.hp2vec[i][2*pr + 1];
      for(int qs = 0; qs < Chan.nhp1[i]; ++qs){
	shell2 = Chan.hp1vec[i][2*qs];
	shell4 = Chan.hp1vec[i][2*qs + 1];
	if(shell1 == shell2 || shell3 == shell4){ continue; }
	TBME = vint_Minnesota_Momentum(Space, shell1, shell2, shell4, shell3, L);
	Ints.D_ME1.V9[i][Chan.nhp1[i]*pr + qs] = TBME;
	Ints.D_ME1.V10[i][Chan.nhp1[i]*pr + qs] = -1.0 * TBME;
      }
    }
  }
}

int kron_del(const int &i, const int &j)
{
  if(i != j){ return 0; }
  return 1;
}
int spinExchangeMtxEle(const int &i, const int &j, const int &k, const int &l)
{
  if(i == l && j == k){ return 1; }
  else{ return 0; }
}

// Minnesota Potential for momentum basis
double vint_Minnesota_Momentum(const Model_Space &Space, const int &qi, const int &qj, const int &qk, const int &ql, const double &L)
{
  double V_R1, V_T1, V_S1, V_R2, V_T2, V_S2;
  double V_0R, V_0T, V_0S;
  double kappa_R, kappa_T, kappa_S;
  double kX1, kY1, kZ1, kX2, kY2, kZ2;
  double qSquared1, spinEx1, isoSpinEx1, qSquared2, spinEx2, isoSpinEx2;
  double IsIt1, PsIt1, PsPt1, IsPt1, IsIt2, PsIt2, PsPt2, IsPt2;
  V_0R = 200; //MeV
  V_0T = 178; //MeV
  V_0S = 91.85; //MeV
  kappa_R = 1.487; //fm^-2
  kappa_T = 0.639; //fm^-2
  kappa_S = 0.465; //fm^-2

  if(Space.qnums[qi].nx + Space.qnums[qj].nx != Space.qnums[qk].nx + Space.qnums[ql].nx){ return 0.0; }
  if(Space.qnums[qi].ny + Space.qnums[qj].ny != Space.qnums[qk].ny + Space.qnums[ql].ny){ return 0.0; }
  if(Space.qnums[qi].nz + Space.qnums[qj].nz != Space.qnums[qk].nz + Space.qnums[ql].nz){ return 0.0; }
  if(Space.qnums[qi].m + Space.qnums[qj].m != Space.qnums[qk].m + Space.qnums[ql].m){ return 0.0; }
  if(Space.qnums[qi].t + Space.qnums[qj].t != Space.qnums[qk].t + Space.qnums[ql].t){ return 0.0; }

  kX1 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[qk].nx + Space.qnums[ql].nx);
  kY1 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[qk].ny + Space.qnums[ql].ny);
  kZ1 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[qk].nz + Space.qnums[ql].nz);

  kX2 = (M_PI/L) * (Space.qnums[qi].nx - Space.qnums[qj].nx - Space.qnums[ql].nx + Space.qnums[qk].nx);
  kY2 = (M_PI/L) * (Space.qnums[qi].ny - Space.qnums[qj].ny - Space.qnums[ql].ny + Space.qnums[qk].ny);
  kZ2 = (M_PI/L) * (Space.qnums[qi].nz - Space.qnums[qj].nz - Space.qnums[ql].nz + Space.qnums[qk].nz);

  qSquared1 = kX1 * kX1 + kY1 * kY1 + kZ1 * kZ1;
  qSquared2 = kX2 * kX2 + kY2 * kY2 + kZ2 * kZ2;
  
  V_R1 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared1/(4*kappa_R));
  V_T1 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared1/(4*kappa_T));
  V_S1 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared1/(4*kappa_S));

  V_R2 = V_0R/(L*L*L)*pow(M_PI/kappa_R,1.5) * exp(-qSquared2/(4*kappa_R));
  V_T2 = -V_0T/(L*L*L)*pow(M_PI/kappa_T,1.5) * exp(-qSquared2/(4*kappa_T));
  V_S2 = -V_0S/(L*L*L)*pow(M_PI/kappa_S,1.5) * exp(-qSquared2/(4*kappa_S));
  
  spinEx1 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[qk].m, Space.qnums[ql].m);
  isoSpinEx1 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[qk].t, Space.qnums[ql].t);

  spinEx2 = spinExchangeMtxEle(Space.qnums[qi].m, Space.qnums[qj].m, Space.qnums[ql].m, Space.qnums[qk].m);
  isoSpinEx2 = spinExchangeMtxEle(Space.qnums[qi].t, Space.qnums[qj].t, Space.qnums[ql].t, Space.qnums[qk].t);
  
  IsIt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m) * kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsIt1 = spinEx1 * kron_del(Space.qnums[qi].t, Space.qnums[qk].t) * kron_del(Space.qnums[qj].t, Space.qnums[ql].t);
  PsPt1 = spinEx1 * isoSpinEx1;
  IsPt1 = kron_del(Space.qnums[qi].m, Space.qnums[qk].m)*kron_del(Space.qnums[qj].m, Space.qnums[ql].m) * isoSpinEx1;

  IsIt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * 
    kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsIt2 = spinEx2 * kron_del(Space.qnums[qi].t, Space.qnums[ql].t) * kron_del(Space.qnums[qj].t, Space.qnums[qk].t);
  PsPt2 = spinEx2 * isoSpinEx2;
  IsPt2 = kron_del(Space.qnums[qi].m, Space.qnums[ql].m) * kron_del(Space.qnums[qj].m, Space.qnums[qk].m) * isoSpinEx2;

  return 0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * IsIt1 + 
    0.25 * (V_T1 - V_S1) * PsIt1 - 
    0.5 * (V_R1 + 0.5*V_T1 + 0.5*V_S1) * PsPt1 - 
    0.25 * (V_T1 - V_S1) * IsPt1 -
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * IsIt2 - 
    0.25 * (V_T2 - V_S2) * PsIt2 +
    0.5 * (V_R2 + 0.5*V_T2 + 0.5*V_S2) * PsPt2 + 
    0.25 * (V_T2 - V_S2) * IsPt2;
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
