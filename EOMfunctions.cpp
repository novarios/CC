#include "EOMfunctions.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "CCfunctions.hpp"
#include "INTfunctions.hpp"
#include "EffINTfunctions.hpp"
#include "MATHfunctions.hpp"

void EOM::PA1_EOM(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  double *Ham;
  State state;
  int length;
  int np0, npph0, npph, np_tot, npph_tot;
  three_body *pph_vec;
  State *J_state;
  int *J_vec;
  int *a_vec;
  int *b_vec;
  int *i_vec;
  int chan_j;
  int chan, N;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int p1, p2, h1;

  // Count N_states and setup EOM structures //
  count_states(Chan, 1);
  np_tot = 0;
  npph_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan]; 
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (PAR.basis == "finite_J" && p1 == p2) ){ ++npph0; }
      //++npph0;
    }
    nob[n] = np0;
    nthb[n] = npph0;
    nstate[n] = np0 + npph0;
    ob_index[n] = np_tot;
    thb_index[n] = npph_tot;
    state_index[n] = np_tot + npph_tot;
    np_tot += np0;
    npph_tot += npph0;
  }
  ob_vec = new one_body[np_tot];
  thb_vec = new three_body[npph_tot];
  thb_qnums = new State[npph_tot];
  state_vec_R = new double[np_tot + npph_tot];
  state_vec_L = new double[np_tot + npph_tot];

  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan];
    chan_j = Chan.qnums3[chan].j;
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (PAR.basis == "finite_J" && p1 == p2) ){
	++npph0;
      }
    }

    pph_vec = new three_body[npph0];
    J_state = new State[npph0];
    J_vec = new int[npph0];
    a_vec = new int[npph0];
    b_vec = new int[npph0];
    i_vec = new int[npph0];

    // set EOM array
    npph0 = 0;
    for(int p = 0; p < np0; ++p){ ob_vec[ob_index[n] + p] = Chan.p_state(chan, p); }
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      h1 = Chan.pph_state(chan, pph).v3;
      if( (p1 < p2) || (PAR.basis == "finite_J" && p1 == p2) ){
	pph_vec[npph0] = Chan.pph_state(chan, pph);
	J_state[npph0] = Chan.pph_j[Chan.pph_index[chan] + pph];
	J_vec[npph0] = Chan.pph_j[Chan.pph_index[chan] + pph].j;
	a_vec[npph0] = SPB.qnums[p1].j;
	b_vec[npph0] = SPB.qnums[p2].j;
	i_vec[npph0] = SPB.qnums[h1].j;
	// set EOM arrays
	thb_vec[thb_index[n] + npph0] = Chan.pph_state(chan, pph);
	thb_qnums[thb_index[n] + npph0] = Chan.pph_j[Chan.pph_index[chan] + pph];
	++npph0;
      }
    }

    N = np0 + npph0;
    Ham = new double[N*N];
    eigenvalues = new double[1];
    eigenvectors_R = new double[N];
    eigenvectors_L = new double[N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }

    #pragma omp parallel
    {
      double ME, X;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind, jmin;
      int bra, ket, key1, key2;
      State tb1, tb2;
      double fac1, fac2;

      length = N * N;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%N);
	bra = int((braket - ket)/N);
	fac2 = 1.0;
	ME = 0.0;

	if( bra >= np0 ){
	  bra -= np0;
	  b1 = pph_vec[bra].v1;
	  b2 = pph_vec[bra].v2;
	  b3 = pph_vec[bra].v3;
	  if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	  if( ket >= np0 ){  //  < pph | pph >
	    ket -= np0;
	    k1 = pph_vec[ket].v1;
	    k2 = pph_vec[ket].v2;
	    k3 = pph_vec[ket].v3;
	    if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      if(b1 == k1 && equal(SPB.qnums[b2], SPB.qnums[k2])){
		minus(tb1, SPB.qnums[b2], SPB.qnums[k2]);
		ind = Chan.pp1_map[Chan.ind0][Hash(b2, k2, 0)];
		ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(b_vec[bra] + 1.0);
		//std::cout << "1: " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xpp.X_2[ind] / std::sqrt(b_vec[bra] + 1.0) << std::endl;
	      }
	      if(b1 == k2 && equal(SPB.qnums[b2], SPB.qnums[k1])){
		minus(tb1, SPB.qnums[b2], SPB.qnums[k1]);
		ind = Chan.pp1_map[Chan.ind0][Hash(b2, k1, 0)];
		ME -= Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(b_vec[bra] + 1.0);
		//std::cout << "2: " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(b_vec[bra] + 1.0) << std::endl;
	      }
	      if(b2 == k2 && equal(SPB.qnums[b1], SPB.qnums[k1])){
		minus(tb1, SPB.qnums[b1], SPB.qnums[k1]);
		ind = Chan.pp1_map[Chan.ind0][Hash(b1, k1, 0)];
		ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[bra] + 1.0);
		//std::cout << "3: " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[bra] + 1.0) << std::endl;
	      }
	      if(b2 == k1 && equal(SPB.qnums[b1], SPB.qnums[k2])){
		minus(tb1, SPB.qnums[b1], SPB.qnums[k2]);
		ind = Chan.pp1_map[Chan.ind0][Hash(b1, k2, 0)];
		ME -= Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(a_vec[bra] + 1.0);
		//std::cout << "4: " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xpp.X_2[ind] * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]) / std::sqrt(a_vec[bra] + 1.0) << std::endl;
	      }
	    }

	    if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) && equal(SPB.qnums[b3], SPB.qnums[k3])){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[b3]);
	      ind = Chan.hh1_map[Chan.ind0][Hash(k3, b3, 0)];
	      ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0);
	      //std::cout << "5_1: " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0) << std::endl;
	      if(b1 == b2){
		ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0);
		//std::cout << "5_2: " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0) << std::endl;
	      }
	    }

	    if(b1 == k1){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k2]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k2].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b2].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "6_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b1 == k2){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k1]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k1].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b2].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "7_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k2){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k1]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k1].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b1].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(a_vec[bra] + b_vec[bra] - J_vec[bra] + a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "8_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k1){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k2]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k2].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b1].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "9_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      chan0 = Ind_Chan1(J_state[bra]);
	      key1 = Chan.pp_map[chan0][Hash(b1, b2, J_state[bra].j)];
	      key2 = Chan.pp_map[chan0][Hash(k1, k2, J_state[ket].j)];
	      chan_ind = Eff_Ints.Xpppp.X_1_index[chan0];
	      ind = key1 * Chan.npp[chan0] + key2;
	      ME += Eff_Ints.Xpppp.X_1[chan_ind + ind];
	      //std::cout << "10: " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xpppp.X_1[chan_ind + ind] << std::endl;
	    }
	    ket += np0;
	  }
	  else{  //  < pph | p >
	    k1 = Chan.p_state(chan, ket).v1;
	    
	    chan0 = Ind_Chan3(SPB.qnums[b3]);
	    key1 = Chan.ppp_map[chan0][Hash(b1, b2, k1, J_vec[bra])];
	    key2 = Chan.h_map[chan0][b3];
	    chan_ind = Eff_Ints.Xpphp.X_3_3_index[chan0];
	    ind = key1 * Chan.nh[chan0] + key2;
	    ME -= Eff_Ints.Xpphp.X_3_3[chan_ind + ind] * std::sqrt((i_vec[bra] + 1.0)/(chan_j + 1.0)) * phase2(i_vec[bra] + chan_j - J_vec[bra]);
	    //std::cout << "11: " << bra+np0 << "|" << ket << " = < " << b1 << "," << b2 << "," << b3 << "|" << k1 << "> = " << -fac2 * Eff_Ints.Xpphp.X_3_3[chan_ind + ind] * std::sqrt((i_vec[bra] + 1.0)/(chan_j + 1.0)) * phase2(i_vec[bra] + chan_j - J_vec[bra]) << std::endl;
	  }
	  bra += np0;
	}
	else if( ket >= np0 ){  //  < p | pph >
	  ket -= np0;
	  k1 = pph_vec[ket].v1;
	  k2 = pph_vec[ket].v2;
	  k3 = pph_vec[ket].v3;
	  b1 = Chan.p_state(chan, bra).v1;
	  if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	  if(k1 == b1 && equal(SPB.qnums[k3], SPB.qnums[k2])){
	    minus(tb1, SPB.qnums[k3], SPB.qnums[k2]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k3, k2, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (b_vec[ket] + 1.0))) * phase2(a_vec[ket] + b_vec[ket] - J_vec[ket]);
	    ME += Eff_Ints.Xhp.X_2[ind] * X;
	    //std::cout << "12: " << bra << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xhp.X_2[ind] * X << std::endl;
	  }
	  if(k2 == b1 && equal(SPB.qnums[k3], SPB.qnums[k1])){
	    minus(tb1, SPB.qnums[k3], SPB.qnums[k1]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k3, k1, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (a_vec[ket] + 1.0)));
	    ME -= Eff_Ints.Xhp.X_2[ind] * X;
	    //std::cout << "13: " << bra << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhp.X_2[ind] * X << std::endl;
	  }

	  plus(tb1, SPB.qnums[k1], SPB.qnums[k2]);
	  key1 = Chan.p_map[chan][b1];
	  key2 = Chan.pph_map[chan][Hash(k1, k2, k3, J_vec[ket])];
	  chan_ind = Eff_Ints.Xhppp.X_3_2_index[chan];
	  ind = key1 * npph + key2;
	  ME -= Eff_Ints.Xhppp.X_3_2[chan_ind + ind];
	  //std::cout << "14: " << bra << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhppp.X_3_2[chan_ind + ind] << std::endl;
	  ket += np0;
	}
	else{  //  < p | p >
	  b1 = Chan.p_state(chan, bra).v1;
	  k1 = Chan.p_state(chan, ket).v1;

	  minus(tb1, SPB.qnums[b1], SPB.qnums[k1]);
	  ind = Chan.pp1_map[Chan.ind0][Hash(b1, k1, 0)];
	  ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(chan_j + 1.0);
	  //std::cout << "15: " << bra << "|" << ket << " = " << fac2 * Eff_Ints.Xpp.X_2[ind] / std::sqrt(chan_j + 1.0) << std::endl;
	}
	Ham[N*ket + bra] += fac2 * ME;
      }
    }

    if(N <= 1){ Asym_Diagonalize1(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }

    ////////////////////////////////////////////
    /*PA PA_Amps, PA_Amps2;
    int count0, key0;
    double fac0;
    PA_Amps.Setup(Chan, chan);
    PA_Amps2.Setup(Chan, PA_Amps);
    
    int num = 2;
    bool rvec = true;
    char howmny = 'A';
    char which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
    int mxiter = 5000;
    int nev = num, ncv = 3*nev + 2; // number of eigenvalues, lanczos vectors to calculate
    if(ncv > N){ ncv = N; }
    int ldv = N, ldz = N, lworkl = 4*ncv*(ncv + 2);
    int mode = 1, ishift = 1, info = 0, ido = 0; // status integer is zero at start
    double sigmar, sigmai, tol = 1.0e-10; // error tolerance
    int *select = new int[ncv];
    double *resid = new double[N];
    double *workev = new double[3 * ncv];
    double *workd = new double[3*N];
    double *workd2 = new double[3*N]; ////////    for testing!!!
    double *workl = new double[lworkl];
    double *v = new double[N*ncv];
    double *dr = new double[nev + 1];
    double *di = new double[nev + 1];
    double *z = new double[N * (nev + 1)];
    int iparam[11], ipntr[14];
    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }
    
    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      for(int j = 0; j < N; ++j){ workd[ipntr[0]-1 + j] = 1.0; }
      std::cout << "0**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[0]-1 + j] << " "; }
      std::cout << std::endl;
      for(int j = 0; j < N; ++j){ // right eigenproblem
	workd[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k]; }
      }
      for(int j = 0; j < N; ++j){ // left eigenproblem
	workd2[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd2[ipntr[1]-1 + j] += Ham[N*j + k] * workd2[ipntr[0]-1 + k]; }
      }
      std::cout << "1R**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      std::cout << "1L**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd2[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      ///////////////////////////////////////////////////////////
      for(int j = 0; j < PA_Amps.np0; ++j){ PA_Amps.R[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){
	PA_Amps.R3[PA_Amps.pph1_ind1[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j];
	PA_Amps.R3[PA_Amps.pph1_ind2[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j] * PA_Amps.pph1_fac1[j];
      }
      PA_Amps2.Zero_R();
      PA_Amps.Update_R1(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      PA_Amps.Update_R2(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      std::cout << "2R**  " << std::setprecision(6);
      for(int j = 0; j < PA_Amps.np0; ++j){ std::cout << PA_Amps2.R[j] << " "; }
      for(int j = 0; j < PA_Amps.npph1; ++j){ std::cout << PA_Amps2.R3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j] << " "; }
      std::cout << std::endl;

      for(int j = 0; j < PA_Amps.np0; ++j){ PA_Amps.L[j] = workd2[ipntr[0]-1 + j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){
	PA_Amps.L3[PA_Amps.pph1_ind1[j]] = workd2[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j];
	PA_Amps.L3[PA_Amps.pph1_ind2[j]] = workd2[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j] * PA_Amps.pph1_fac1[j];
      }
      PA_Amps2.Zero_L();
      PA_Amps.Update_L1(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      PA_Amps.Update_L2(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      std::cout << "2L**  " << std::setprecision(6);
      for(int j = 0; j < PA_Amps.np0; ++j){ std::cout << PA_Amps2.L[j] << " "; }
      for(int j = 0; j < PA_Amps.npph1; ++j){ std::cout << PA_Amps2.L3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j] << " "; }
      std::cout << std::endl;
      //////////////////////////////////////////////////////////
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    
    //for(int i = 0; i < num; ++i){
    //std::cout << "!!  " << dr[i] << ", " << N << std::endl;
    //}
    //std::cout << std::endl;

    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;
    
    PA_Amps.Delete();
    PA_Amps2.Delete();*/
    ////////////////////////////////////////////////////

    norm1p = 0.0;
    for(int i = 0; i < np0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);

    del_E[n] = eigenvalues[0];
    norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
    delete[] pph_vec;
    delete[] J_state;
    delete[] J_vec;
    delete[] a_vec;
    delete[] b_vec;
    delete[] i_vec;
    delete[] Ham;
  }
}

void EOM::PR1_EOM(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  double *Ham;
  State state;
  int length;
  int nh0, nhhp0, nhhp, nh_tot, nhhp_tot;//, nhhh;
  three_body *hhp_vec;
  State *J_state;
  int *J_vec;
  int *i_vec;
  int *j_vec;
  int *a_vec;
  int chan_j;
  int chan, N;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int h1, h2, p1;

  // Count N_states and setup EOM structures //
  count_states(Chan, 0);
  nh_tot = 0;
  nhhp_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan]; 
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (PAR.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    nob[n] = nh0;
    nthb[n] = nhhp0;
    nstate[n] = nh0 + nhhp0;
    ob_index[n] = nh_tot;
    thb_index[n] = nhhp_tot;
    state_index[n] = nh_tot + nhhp_tot;
    nh_tot += nh0;
    nhhp_tot += nhhp0;
  }
  ob_vec = new one_body[nh_tot];
  thb_vec = new three_body[nhhp_tot];
  thb_qnums = new State[nhhp_tot];
  state_vec_R = new double[nh_tot + nhhp_tot];
  state_vec_L = new double[nh_tot + nhhp_tot];
  /////////////////////////////////////////////
  std::cout << "1PR-EOM" << std::endl;

  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan];
    chan_j = Chan.qnums3[chan].j;
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (PAR.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    hhp_vec = new three_body[nhhp0];
    J_state = new State[nhhp0];
    J_vec = new int[nhhp0];
    i_vec = new int[nhhp0];
    j_vec = new int[nhhp0];
    a_vec = new int[nhhp0];

    // set EOM array    
    nhhp0 = 0;
    for(int h = 0; h < nh0; ++h){ ob_vec[ob_index[n] + h] = Chan.h_state(chan, h); }
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      p1 = Chan.hhp_state(chan, hhp).v3;
      if( (h1 < h2) || (PAR.basis == "finite_J" && h1 == h2) ){
	hhp_vec[nhhp0] = Chan.hhp_state(chan, hhp);
	J_state[nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp];
	J_vec[nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp].j;
	i_vec[nhhp0] = SPB.qnums[h1].j;
	j_vec[nhhp0] = SPB.qnums[h2].j;
	a_vec[nhhp0] = SPB.qnums[p1].j;
	// set EOM arrays
	thb_vec[thb_index[n] + nhhp0] = Chan.hhp_state(chan, hhp);
	thb_qnums[thb_index[n] + nhhp0] = Chan.hhp_j[Chan.hhp_index[chan] + hhp];
	++nhhp0;
      }
    }

    N = nh0 + nhhp0;
    Ham = new double[N*N];
    eigenvalues = new double[1];
    eigenvectors_R = new double[N];
    eigenvectors_L = new double[N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }
    
    #pragma omp parallel
    {
      double ME, X;
      int b1, b2, b3, k1, k2, k3, chan0, ind, chan_ind, jmin;
      int bra, ket, key1, key2;
      State tb1, tb2;
      double fac1, fac2;

      length = N * N;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%N);
	bra = int((braket - ket)/N);
	fac2 = 1.0;
	ME = 0.0;
	  
	if( bra >= nh0 ){
	  bra -= nh0;
	  b1 = hhp_vec[bra].v1;
	  b2 = hhp_vec[bra].v2;
	  b3 = hhp_vec[bra].v3;
	  if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	  if( ket >= nh0 ){  //  < hhp | hhp >
	    ket -= nh0;
	    k1 = hhp_vec[ket].v1;
	    k2 = hhp_vec[ket].v2;
	    k3 = hhp_vec[ket].v3;
	    if(k1 == k2){ fac2 /= std::sqrt(2.0); }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      if(b1 == k1 && equal(SPB.qnums[b2], SPB.qnums[k2])){
		minus(tb1, SPB.qnums[k2], SPB.qnums[b2]);
		ind = Chan.hh1_map[Chan.ind0][Hash(k2, b2, 0)];
		ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(j_vec[bra] + 1.0);
	      }
	      if(b1 == k2 && equal(SPB.qnums[b2], SPB.qnums[k1])){
		minus(tb1, SPB.qnums[k1], SPB.qnums[k2]);
		ind = Chan.hh1_map[Chan.ind0][Hash(k1, b2, 0)];
		ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(j_vec[bra] + 1.0);
	      }
	      if(b2 == k2 && equal(SPB.qnums[b1], SPB.qnums[k1])){
		minus(tb1, SPB.qnums[k1], SPB.qnums[b1]);
		ind = Chan.hh1_map[Chan.ind0][Hash(k1, b1, 0)];
		ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[bra] + 1.0);
	      }
	      if(b2 == k1 && equal(SPB.qnums[b1], SPB.qnums[k2])){
		minus(tb1, SPB.qnums[k2], SPB.qnums[b1]);
		ind = Chan.hh1_map[Chan.ind0][Hash(k2, b1, 0)];
		ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(i_vec[bra] + 1.0);
	      }
	    }

	    if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) && equal(SPB.qnums[b3], SPB.qnums[k3])){
	      minus(tb1, SPB.qnums[b3], SPB.qnums[k3]);
	      ind = Chan.pp1_map[Chan.ind0][Hash(b3, k3, 0)];
	      ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0);
	      if(b1 == b2){ ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0); }
	    }

	    if(b1 == k1){
	      minus(tb1, SPB.qnums[k2], SPB.qnums[k3]);
	      minus(tb2, SPB.qnums[b2], SPB.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k2].j - SPB.qnums[k3].j);
		if(abs(SPB.qnums[b2].j - SPB.qnums[b3].j) > jmin){ jmin = abs(SPB.qnums[b2].j - SPB.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(j_vec[bra] + j_vec[ket] + a_vec[bra] + a_vec[ket]);
		  X *= CGC6(i_vec[bra],j_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(i_vec[ket],j_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k2, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b2, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b1 == k2){
	      minus(tb1, SPB.qnums[k1], SPB.qnums[k3]);
	      minus(tb2, SPB.qnums[b2], SPB.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k1].j - SPB.qnums[k3].j);
		if(abs(SPB.qnums[b2].j - SPB.qnums[b3].j) > jmin){ jmin = abs(SPB.qnums[b2].j - SPB.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);;
		  X *= phase2(j_vec[bra] + i_vec[ket] + a_vec[bra] + a_vec[ket] + i_vec[ket] + j_vec[ket] - J_vec[ket]);
		  X *= CGC6(i_vec[bra],j_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(j_vec[ket],i_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k1, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b2, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k2){
	      minus(tb1, SPB.qnums[k1], SPB.qnums[k3]);
	      minus(tb2, SPB.qnums[b1], SPB.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k1].j - SPB.qnums[k3].j);
		if(abs(SPB.qnums[b1].j - SPB.qnums[b3].j) > jmin){ jmin = abs(SPB.qnums[b1].j - SPB.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(i_vec[bra] + i_vec[ket] + a_vec[bra] + a_vec[ket]);
		  X *= phase2(i_vec[bra] + j_vec[bra] - J_vec[bra] + i_vec[ket] + j_vec[ket] - J_vec[ket]);
		  X *= CGC6(j_vec[bra],i_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(j_vec[ket],i_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k1, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b1, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k1){
	      minus(tb1, SPB.qnums[k2], SPB.qnums[k3]);
	      minus(tb2, SPB.qnums[b1], SPB.qnums[b3]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k2].j - SPB.qnums[k3].j);
		if(abs(SPB.qnums[b1].j - SPB.qnums[b3].j) > jmin){ jmin = abs(SPB.qnums[b1].j - SPB.qnums[b3].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(i_vec[bra] + j_vec[ket] + a_vec[bra] + a_vec[ket] + i_vec[bra] + j_vec[bra] - J_vec[bra]);
		  X *= CGC6(j_vec[bra],i_vec[bra],J_vec[bra],a_vec[bra],chan_j,tb1.j) * CGC6(i_vec[ket],j_vec[ket],J_vec[ket],a_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k2, k3, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b1, b3, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }

	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      chan0 = Ind_Chan1(J_state[bra]);
	      key1 = Chan.hh_map[chan0][Hash(k1, k2, J_state[ket].j)];
	      key2 = Chan.hh_map[chan0][Hash(b1, b2, J_state[bra].j)];
	      chan_ind = Eff_Ints.Xhhhh.X_1_index[chan0];
	      ind = key1 * Chan.nhh[chan0] + key2;
	      ME += Eff_Ints.Xhhhh.X_1[chan_ind + ind];
	    }
	    ket += nh0;
	  }
	  else{  //  < hhp | h >
	    k1 = Chan.h_state(chan, ket).v1;

	    chan0 = Ind_Chan3(SPB.qnums[b3]);
	    key1 = Chan.p_map[chan0][b3];
	    key2 = Chan.hhh_map[chan0][Hash(b1, b2, k1, J_vec[bra])];
	    chan_ind = Eff_Ints.Xhphh.X_3_2_index[chan0];
	    ind = key1 * Chan.nhhh[chan0] + key2;
	    ME -= Eff_Ints.Xhphh.X_3_2[chan_ind + ind] * std::sqrt((a_vec[bra] + 1.0) / (chan_j + 1.0)) * phase2(chan_j + a_vec[bra] - J_vec[bra]);
	    //std::cout << "11: " << bra+nh0 << "|" << ket << " = < " << b1 << "," << b2 << "," << b3 << "|" << k1 << "> = " << -Eff_Ints.Xhphh.X_3_2[chan_ind + ind] * std::sqrt((a_vec[bra] + 1.0) / (chan_j + 1.0)) * phase2(chan_j + a_vec[bra] - J_vec[bra]) << std::endl;
	  }
	  bra += nh0;
	}
	else if( ket >= nh0 ){  //  < h | hhp >
	  ket -= nh0;
	  k1 = hhp_vec[ket].v1;
	  k2 = hhp_vec[ket].v2;
	  k3 = hhp_vec[ket].v3;
	  b1 = Chan.h_state(chan, bra).v1;
	  if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	    
	  if(k1 == b1 && equal(SPB.qnums[k2], SPB.qnums[k3])){
	    minus(tb1, SPB.qnums[k2], SPB.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k2, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (j_vec[ket] + 1.0))) * phase2(i_vec[ket] + j_vec[ket] - J_vec[ket]);
	    ME += Eff_Ints.Xhp.X_2[ind] * X;
	  }
	  if(k2 == b1 && equal(SPB.qnums[k1], SPB.qnums[k3])){
	    minus(tb1, SPB.qnums[k1], SPB.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k1, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (i_vec[ket] + 1.0)));
	    ME -= Eff_Ints.Xhp.X_2[ind] * X;
	  }

	  key1 = Chan.hhp_map[chan][Hash(k1, k2, k3, J_vec[ket])];
	  key2 = Chan.h_map[chan][b1];
	  chan_ind = Eff_Ints.Xhhhp.X_3_3_index[chan];
	  ind = key1 * nh0 + key2;
	  ME += -1.0 * Eff_Ints.Xhhhp.X_3_3[chan_ind + ind];
	  ket += nh0;
	}
	else{  //  < h | h >
	  b1 = Chan.h_state(chan, bra).v1;
	  k1 = Chan.h_state(chan, ket).v1;

	  minus(tb1, SPB.qnums[k1], SPB.qnums[b1]);
	  ind = Chan.hh1_map[Chan.ind0][Hash(k1, b1, 0)];
	  ME += -1.0 * Eff_Ints.Xhh.X_2[ind] / std::sqrt(chan_j + 1.0);
	}
	Ham[N*ket + bra] += fac2 * ME;
      }
    }
    
    ////////////////////////////////////////////
    /*PR PR_Amps, PR_Amps2;
    int count0, key0;
    double fac0;
    PR_Amps.Setup(Chan, chan);
    PR_Amps2.Setup(Chan, PR_Amps);
    
    int num = 2;
    bool rvec = true;
    char howmny = 'A';
    char which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
    int mxiter = 5000;
    int nev = num, ncv = 3*nev + 2; // number of eigenvalues, lanczos vectors to calculate
    if(ncv > N){ ncv = N; }
    int ldv = N, ldz = N, lworkl = 4*ncv*(ncv + 2);
    int mode = 1, ishift = 1, info = 0, ido = 0; // status integer is zero at start
    double sigmar, sigmai, tol = 1.0e-10; // error tolerance
    int *select = new int[ncv];
    double *resid = new double[N];
    double *workev = new double[3 * ncv];
    double *workd = new double[3*N];
    double *workd2 = new double[3*N]; ////////    for testing!!!
    double *workl = new double[lworkl];
    double *v = new double[N*ncv];
    double *dr = new double[nev + 1];
    double *di = new double[nev + 1];
    double *z = new double[N * (nev + 1)];
    int iparam[11], ipntr[14];
    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }
    
    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      //for(int j = 0; j < N; ++j){ workd[ipntr[0]-1 + j] = 1.0; }
      std::cout << "0**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[0]-1 + j] << " "; }
      std::cout << std::endl;
      for(int j = 0; j < N; ++j){ // right eigenproblem
	workd[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k]; }
      }
      for(int j = 0; j < N; ++j){ // left eigenproblem
	workd2[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd2[ipntr[1]-1 + j] += Ham[N*j + k] * workd[ipntr[0]-1 + k]; }
      }
      std::cout << "1R**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      std::cout << "1L**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd2[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      ///////////////////////////////////////////////////////////
      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.R[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.R3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.R3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_R();
      PR_Amps.Update_R1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_R2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      std::cout << "2R**  " << std::setprecision(6);
      for(int j = 0; j < PR_Amps.nh0; ++j){ std::cout << PR_Amps2.R[j] << " "; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ std::cout << PR_Amps2.R3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j] << " "; }
      std::cout << std::endl;

      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.L[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.L3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.L3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_L();
      PR_Amps.Update_L1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_L2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      std::cout << "2L**  " << std::setprecision(6);
      for(int j = 0; j < PR_Amps.nh0; ++j){ std::cout << PR_Amps2.L[j] << " "; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ std::cout << PR_Amps2.L3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j] << " "; }
      std::cout << std::endl;
      //////////////////////////////////////////////////////////
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    
    //for(int i = 0; i < num; ++i){
    //std::cout << "!!  " << dr[i] << ", " << N << std::endl;
    //}
    //std::cout << std::endl;

    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;
    
    PR_Amps.Delete();
    PR_Amps2.Delete();*/
    ////////////////////////////////////////////////////

    if(N <= 1){ Asym_Diagonalize1(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }

    norm1p = 0.0;
    for(int i = 0; i < nh0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);

    del_E[n] = eigenvalues[0];
    norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
    delete[] hhp_vec;
    delete[] J_state;
    delete[] J_vec;
    delete[] i_vec;
    delete[] j_vec;
    delete[] a_vec;
    delete[] Ham;
  }
}

/*void EOM::PN_EOM(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  double *Ham;
  State state;
  int length;
  int nph0, npphh0, npphh, nph_tot, npphh_tot;
  four_body *pphh_vec;
  State *J_state1;
  State *J_state2;
  int *J_vec1;
  int *J_vec2;
  int *a1_vec;
  int *i1_vec;
  int *a_vec;
  int *b_vec;
  int *i_vec;
  int *j_vec;
  int chan_j;
  int chan, N;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int p1, p2, h1, h2;
  int *npphh_vec;
  int *npphh0_vec;

  // Count N_states and setup EOM structures //
  count_states2(Chan);
  npphh_vec = new int[N_states];
  npphh0_vec = new int[N_states];
  nph_tot = 0;
  npphh_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nph0 = Chan.nph1[chan];
    npphh0 = 0;
    npphh = 0;
    for(int chan1_1 = 0; chan1_1 < Chan.size1; ++chan1_1){
      if( Chan.npp[chan1_1] == 0 ){ continue; }
      for(int chan1_2 = 0; chan1_2 < Chan.size1; ++chan1_2){
	if( Chan.nhh[chan1_2] == 0 ){ continue; }

	minus(fb, Chan.qnums1[chan1_1], Chan.qnums1[chan1_2]);	
	jmin = std::abs(Chan.qnums1[chan1_1].j - Chan.qnums1[chan1_2].j);
	while(fb.j >= jmin){
	  if( Ind_Chan2(fb) == chan ){
	    for(int pp = 0; pp < Chan.npp[chan1_1]; ++pp){
	      p1 = Chan.pp_state(chan, pp).v1;
	      p2 = Chan.pp_state(chan, pp).v2;
	      for(int hh = 0; hh < Chan.nhh[chan1_2]; ++hh){
		h1 = Chan.hh_state(chan, hh).v1;
		h2 = Chan.hh_state(chan, hh).v2;
		+npphh;
		if( (p1 < p2 && h1 < h2) || (PAR.basis == "finite_J" && (p1 <= p2 && h1 <= h2)) ){ ++npphh0; }
	      }
	    }
	  }
	  fb.j -= 2;
	}
      }
    }
    npphh_vec[n] = npphh;
    npphh0_vec[n] = npphh0;

    ntb[n] = nph0;
    nfb[n] = npphh0;
    nstate[n] = nph0 + npphh0;
    tb_index[n] = nph_tot;
    fb_index[n] = npphh_tot;
    state_index[n] = nph_tot + npphh_tot;
    nph_tot += nph0;
    npphh_tot += npphh0;
  }
  tb_vec = new one_body[nph_tot];
  fb_vec = new three_body[npphh_tot];
  fb1_qnums = new State[npphh_tot];
  fb2_qnums = new State[npphh_tot];
  state_vec_R = new double[nph_tot + npphh_tot];
  state_vec_L = new double[nph_tot + npphh_tot];
  /////////////////////////////////////////////
  std::cout << "PN-EOM" << std::endl;

  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nph0 = Chan.nph1[chan];
    npphh0 = npphh0_vec[n];
    npphh = npphh_vec[n];
    chan_j = Chan.qnums2[chan].j;

    a1_vec = new int[nph0];
    i1_vec = new int[nph0];
    pphh_vec = new four_body[npphh0];
    J_state1 = new State[npphh0];
    J_state2 = new State[npphh0];
    J_vec1 = new int[npphh0];
    J_vec1 = new int[npphh0];
    a_vec = new int[npphh0];
    b_vec = new int[npphh0];
    i_vec = new int[npphh0];
    j_vec = new int[npphh0];

    // set EOM array    
    npphh0 = 0;
    for(int ph1 = 0; ph1 < nph0; ++ph1){ tb_vec[tb_index[n] + ph1] = Chan.ph1_state(chan, ph1); }
    for(int chan1_1 = 0; chan1_1 < Chan.size1; ++chan1_1){
      if( Chan.npp[chan1_1] == 0 ){ continue; }
      for(int chan1_2 = 0; chan1_2 < Chan.size1; ++chan1_2){
	if( Chan.nhh[chan1_2] == 0 ){ continue; }

	minus(fb, Chan.qnums1[chan1_1], Chan.qnums1[chan1_2]);	
	jmin = std::abs(Chan.qnums1[chan1_1].j - Chan.qnums1[chan1_2].j);
	while(fb.j >= jmin){
	  if( Ind_Chan2(fb) == chan ){
	    for(int pp = 0; pp < Chan.npp[chan1_1]; ++pp){
	      p1 = Chan.pp_state(chan, pp).v1;
	      p2 = Chan.pp_state(chan, pp).v2;
	      for(int hh = 0; hh < Chan.nhh[chan1_2]; ++hh){
		h1 = Chan.hh_state(chan, hh).v1;
		h2 = Chan.hh_state(chan, hh).v2;
		if( (p1 < p2 && h1 < h2) || (PAR.basis == "finite_J" && (p1 <= p2 && h1 <= h2)) ){
		  pphh_vec[npphh0].v1 = Chan.pp_state(chan1_1, pp).v1;
		  pphh_vec[npphh0].v2 = Chan.pp_state(chan1_1, pp).v2;
		  pphh_vec[npphh0].v3 = Chan.hh_state(chan1_2, hh).v1;
		  pphh_vec[npphh0].v4 = Chan.hh_state(chan1_2, hh).v2;
		  J_state1[npphh0] = Chan.qnums1[chan1_1];
		  J_state2[npphh0] = Chan.qnums1[chan1_2];
		  J_vec1[npphh0] = Chan.qnums1[chan1_1].j;
		  J_vec2[npphh0] = Chan.qnums1[chan1_2].j;
		  a_vec[npphh0] = SPB.qnums[p1].j;
		  b_vec[npphh0] = SPB.qnums[p2].j;
		  i_vec[npphh0] = SPB.qnums[h1].j;
		  j_vec[npphh0] = SPB.qnums[h2].j;
		  // set EOM arrays
		  fb_vec[fb_index[n] + npphh0].v1 = Chan.pp_state(chan1_1, pp).v1;
		  fb_vec[fb_index[n] + npphh0].v2 = Chan.pp_state(chan1_1, pp).v2;
		  fb_vec[fb_index[n] + npphh0].v3 = Chan.hh_state(chan1_2, hh).v1;
		  fb_vec[fb_index[n] + npphh0].v4 = Chan.hh_state(chan1_2, hh).v2;
		  fb1_qnums[fb_index[n] + npphh0] = Chan.qnums1[chan1_1];
		  fb2_qnums[fb_index[n] + npphh0] = Chan.qnums1[chan1_2];
		  ++npphh0;
		}
	      }
	    }
	  }
	  fb.j -= 2;
	}
      }
    }

    N = nph0 + npphh0;
    Ham = new double[N*N];
    eigenvalues = new double[1];
    eigenvectors_R = new double[N];
    eigenvectors_L = new double[N];
    for(int ind = 0; ind < N*N; ++ind){ Ham[ind] = 0.0; }
    
    #pragma omp parallel
    {
      double ME, X;
      int b1, b2, b3, b4, k1, k2, k3, k4, chan0, ind, chan_ind, jmin;
      int bra, ket, key1, key2;
      State tb1, tb2;
      double fac1, fac2;

      length = N * N;
      #pragma omp for schedule(static)
      for(int braket = 0; braket < length; ++braket){
	ket = int(braket%N);
	bra = int((braket - ket)/N);
	fac2 = 1.0;
	ME = 0.0;
	  
	if( bra >= nph0 ){
	  bra -= nph0;
	  b1 = pphh_vec[bra].v1;
	  b2 = pphh_vec[bra].v2;
	  b3 = pphh_vec[bra].v3;
	  b4 = pphh_vec[bra].v3;
	  if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	  if(b3 == b4){ fac2 /= std::sqrt(2.0); }
	  if( ket >= nph0 ){  //  < pphh | pphh >
	    ket -= nph0;
	    k1 = pphh_vec[ket].v1;
	    k2 = pphh_vec[ket].v2;
	    k3 = pphh_vec[ket].v3;
	    k4 = pphh_vec[ket].v4;
	    if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	    if(k3 == k4){ fac2 /= std::sqrt(2.0); }

	    if(b3 == k3 && b4 == k4  && b2 == k2 && equal(J_state2[bra], J_state2[ket]) && equal(SPB.qnums[b1], SPB.qnums[k1])){
	      minus(tb1, SPB.qnums[b1], SPB.qnums[k1]);
	      ind = Chan.pp1_map[Chan.ind0][Hash(b1, k1, 0)];
	      ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0);
	      if(b3 == b4){ ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0); }
	    }
	    if(b3 == k3 && b4 == k4  && b1 == k1 && equal(J_state2[bra], J_state2[ket]) && equal(SPB.qnums[b2], SPB.qnums[k2])){
	      minus(tb1, SPB.qnums[b2], SPB.qnums[k2]);
	      ind = Chan.pp1_map[Chan.ind0][Hash(b2, k2, 0)];
	      ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(b_vec[ket] + 1.0);
	      if(b3 == b4){ ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(b_vec[ket] + 1.0); }
	    }

	    if(b1 == k1 && b2 == k2  && b4 == k4 && equal(J_state1[bra], J_state1[ket]) && equal(SPB.qnums[b3], SPB.qnums[k3])){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[b3]);
	      ind = Chan.hh1_map[Chan.ind0][Hash(k3, b3, 0)];
	      ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0);
	      if(b1 == b2){ ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[ket] + 1.0); }
	    }
	    if(b1 == k1 && b2 == k2  && b3 == k3 && equal(J_state1[bra], J_state1[ket]) && equal(SPB.qnums[b4], SPB.qnums[k4])){
	      minus(tb1, SPB.qnums[k4], SPB.qnums[b4]);
	      ind = Chan.hh1_map[Chan.ind0][Hash(k4, b4, 0)];
	      ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(j_vec[ket] + 1.0);
	      if(b1 == b2){ ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(j_vec[ket] + 1.0); }
	    }

	    if(b1 == k1 && b3 == k3){}
	    if(b1 == k1 && b4 == k4){}
	    if(b1 == k1 && b3 == k4){}
	    if(b1 == k1 && b4 == k3){}

	    if(b2 == k2 && b3 == k3){}
	    if(b2 == k2 && b4 == k4){}
	    if(b2 == k2 && b3 == k4){}
	    if(b2 == k2 && b4 == k3){}

	    if(b1 == k2 && b3 == k3){}
	    if(b1 == k2 && b4 == k4){}
	    if(b1 == k2 && b3 == k4){}
	    if(b1 == k2 && b4 == k3){}

	    if(b2 == k1 && b3 == k3){}
	    if(b2 == k1 && b4 == k4){}
	    if(b2 == k1 && b3 == k4){}
	    if(b2 == k1 && b4 == k3){}


	    if(b1 == k1 && b4 == k4){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k2]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k2].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b2].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b1 == k2){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k1]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b2]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k1].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b2].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b2].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(a_vec[bra],b_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b2, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "7_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k2){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k1]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k1].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b1].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0);
		  X *= phase2(a_vec[bra] + b_vec[bra] - J_vec[bra] + a_vec[ket] + b_vec[ket] - J_vec[ket]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(b_vec[ket],a_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k1, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME -= Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "8_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << -fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }
	    if(b2 == k1){
	      minus(tb1, SPB.qnums[k3], SPB.qnums[k2]);
	      minus(tb2, SPB.qnums[b3], SPB.qnums[b1]);
	      if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	      else{ tb2.j = tb1.j; }
	      if( equal(tb1, tb2) ){
		jmin = abs(SPB.qnums[k3].j - SPB.qnums[k2].j);
		if(abs(SPB.qnums[b3].j - SPB.qnums[b1].j) > jmin){ jmin = abs(SPB.qnums[b3].j - SPB.qnums[b1].j); }
		while(tb1.j >= jmin){
		  X = std::sqrt((J_vec[bra] + 1.0) * (J_vec[ket] + 1.0)) * (tb1.j + 1.0) * phase2(a_vec[bra] + b_vec[bra] - J_vec[bra]);
		  X *= CGC6(b_vec[bra],a_vec[bra],J_vec[bra],i_vec[bra],chan_j,tb1.j) * CGC6(a_vec[ket],b_vec[ket],J_vec[ket],i_vec[ket],chan_j,tb1.j);
		  chan0 = Ind_Chan2(tb1);
		  key1 = Chan.hp1_map[chan0][Hash(k3, k2, tb1.j)];
		  key2 = Chan.hp1_map[chan0][Hash(b3, b1, tb1.j)];
		  chan_ind = Eff_Ints.Xhphp.X_2_1_index[chan0];
		  ind = key1 * Chan.nhp1[chan0] + key2;
		  ME += Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X;
		  //std::cout << "9_" << tb1.j << ": " << bra+np0 << "|" << ket+np0 << " = " << fac2 * Eff_Ints.Xhphp.X_2_1[chan_ind + ind] * X << std::endl;
		  tb1.j -= 2;
		}
	      }
	    }


	    if(b3 == k3 && equal(J_state[bra], J_state[ket])){
	      chan0 = Ind_Chan1(J_state[bra]);
	      key1 = Chan.hh_map[chan0][Hash(k1, k2, J_state[ket].j)];
	      key2 = Chan.hh_map[chan0][Hash(b1, b2, J_state[bra].j)];
	      chan_ind = Eff_Ints.Xhhhh.X_1_index[chan0];
	      ind = key1 * Chan.nhh[chan0] + key2;
	      ME += Eff_Ints.Xhhhh.X_1[chan_ind + ind];
	    }
	    ket += nh0;
	  }
	  else{  //  < hhp | h >
	    k1 = Chan.h_state(chan, ket).v1;

	    chan0 = Ind_Chan3(SPB.qnums[b3]);
	    key1 = Chan.p_map[chan0][b3];
	    key2 = Chan.hhh_map[chan0][Hash(b1, b2, k1, J_vec[bra])];
	    chan_ind = Eff_Ints.Xhphh.X_3_2_index[chan0];
	    ind = key1 * Chan.nhhh[chan0] + key2;
	    ME -= Eff_Ints.Xhphh.X_3_2[chan_ind + ind] * std::sqrt((a_vec[bra] + 1.0) / (chan_j + 1.0)) * phase2(chan_j + a_vec[bra] - J_vec[bra]);
	    //std::cout << "11: " << bra+nh0 << "|" << ket << " = < " << b1 << "," << b2 << "," << b3 << "|" << k1 << "> = " << -Eff_Ints.Xhphh.X_3_2[chan_ind + ind] * std::sqrt((a_vec[bra] + 1.0) / (chan_j + 1.0)) * phase2(chan_j + a_vec[bra] - J_vec[bra]) << std::endl;
	  }
	  bra += nh0;
	}
	else if( ket >= nh0 ){  //  < h | hhp >
	  ket -= nh0;
	  k1 = hhp_vec[ket].v1;
	  k2 = hhp_vec[ket].v2;
	  k3 = hhp_vec[ket].v3;
	  b1 = Chan.h_state(chan, bra).v1;
	  if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	    
	  if(k1 == b1 && equal(SPB.qnums[k2], SPB.qnums[k3])){
	    minus(tb1, SPB.qnums[k2], SPB.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k2, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (j_vec[ket] + 1.0))) * phase2(i_vec[ket] + j_vec[ket] - J_vec[ket]);
	    ME += Eff_Ints.Xhp.X_2[ind] * X;
	  }
	  if(k2 == b1 && equal(SPB.qnums[k1], SPB.qnums[k3])){
	    minus(tb1, SPB.qnums[k1], SPB.qnums[k3]);
	    ind = Chan.hp1_map[Chan.ind0][Hash(k1, k3, 0)];
	    X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (i_vec[ket] + 1.0)));
	    ME -= Eff_Ints.Xhp.X_2[ind] * X;
	  }

	  key1 = Chan.hhp_map[chan][Hash(k1, k2, k3, J_vec[ket])];
	  key2 = Chan.h_map[chan][b1];
	  chan_ind = Eff_Ints.Xhhhp.X_3_3_index[chan];
	  ind = key1 * nh0 + key2;
	  ME += -1.0 * Eff_Ints.Xhhhp.X_3_3[chan_ind + ind];
	  ket += nh0;
	}
	else{  //  < h | h >
	  b1 = Chan.h_state(chan, bra).v1;
	  k1 = Chan.h_state(chan, ket).v1;

	  minus(tb1, SPB.qnums[k1], SPB.qnums[b1]);
	  ind = Chan.hh1_map[Chan.ind0][Hash(k1, b1, 0)];
	  ME += -1.0 * Eff_Ints.Xhh.X_2[ind] / std::sqrt(chan_j + 1.0);
	}
	Ham[N*ket + bra] += fac2 * ME;
      }
    }
    
    ////////////////////////////////////////////
    PR PR_Amps, PR_Amps2;
    int count0, key0;
    double fac0;
    PR_Amps.Setup(Chan, chan);
    PR_Amps2.Setup(Chan, PR_Amps);
    
    int num = 2;
    bool rvec = true;
    char howmny = 'A';
    char which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
    int mxiter = 5000;
    int nev = num, ncv = 3*nev + 2; // number of eigenvalues, lanczos vectors to calculate
    if(ncv > N){ ncv = N; }
    int ldv = N, ldz = N, lworkl = 4*ncv*(ncv + 2);
    int mode = 1, ishift = 1, info = 0, ido = 0; // status integer is zero at start
    double sigmar, sigmai, tol = 1.0e-10; // error tolerance
    int *select = new int[ncv];
    double *resid = new double[N];
    double *workev = new double[3 * ncv];
    double *workd = new double[3*N];
    double *workd2 = new double[3*N]; ////////    for testing!!!
    double *workl = new double[lworkl];
    double *v = new double[N*ncv];
    double *dr = new double[nev + 1];
    double *di = new double[nev + 1];
    double *z = new double[N * (nev + 1)];
    int iparam[11], ipntr[14];
    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }
    
    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      //for(int j = 0; j < N; ++j){ workd[ipntr[0]-1 + j] = 1.0; }
      std::cout << "0**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[0]-1 + j] << " "; }
      std::cout << std::endl;
      for(int j = 0; j < N; ++j){ // right eigenproblem
	workd[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k]; }
      }
      for(int j = 0; j < N; ++j){ // left eigenproblem
	workd2[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd2[ipntr[1]-1 + j] += Ham[N*j + k] * workd[ipntr[0]-1 + k]; }
      }
      std::cout << "1R**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      std::cout << "1L**  " << std::setprecision(6);
      for(int j = 0; j < N; ++j){ std::cout << workd2[ipntr[1]-1 + j] << " "; }
      std::cout << std::endl;
      ///////////////////////////////////////////////////////////
      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.R[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.R3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.R3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_R();
      PR_Amps.Update_R1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_R2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      std::cout << "2R**  " << std::setprecision(6);
      for(int j = 0; j < PR_Amps.nh0; ++j){ std::cout << PR_Amps2.R[j] << " "; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ std::cout << PR_Amps2.R3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j] << " "; }
      std::cout << std::endl;

      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.L[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.L3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.L3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_L();
      PR_Amps.Update_L1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_L2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      std::cout << "2L**  " << std::setprecision(6);
      for(int j = 0; j < PR_Amps.nh0; ++j){ std::cout << PR_Amps2.L[j] << " "; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ std::cout << PR_Amps2.L3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j] << " "; }
      std::cout << std::endl;
      //////////////////////////////////////////////////////////
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    
    //for(int i = 0; i < num; ++i){
    //std::cout << "!!  " << dr[i] << ", " << N << std::endl;
    //}
    //std::cout << std::endl;

    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;
    
    PR_Amps.Delete();
    PR_Amps2.Delete();
    ////////////////////////////////////////////////////

    if(N <= 1){ Asym_Diagonalize1(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }
    else{ Asym_Diagonalize2(Ham, N, eigenvalues, eigenvectors_L, eigenvectors_R, 1); }

    norm1p = 0.0;
    for(int i = 0; i < nh0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);

    del_E[n] = eigenvalues[0];
    norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
    delete[] hhp_vec;
    delete[] J_state;
    delete[] J_vec;
    delete[] i_vec;
    delete[] j_vec;
    delete[] a_vec;
    delete[] Ham;
  }
}*/

void EOM::count_states(Channels &Chan, int type)
{
  int count1 = 0, count2 = 0, count_limit1 = 0, count_limit2 = 0;
  int num, ind, index;
  double en, tempen;
  double limit;
  int p_limit = 6;
  int n_limit = 0;
  if( type == 1 ){ limit = 1.0e10; }
  else{ limit = -1.0e10; } // type = 0

  if( PAR.Pshells != 0 ){
    for(int chan = 0; chan < Chan.size3; ++chan){
      if( type == 1 ){ num = Chan.np[chan]; }
      else{ num = Chan.nh[chan]; } // type = 0
      if( num == 0 || Chan.qnums3[chan].t == 1 ){ continue; } // skip if no particles/holes in this channel
      if(PAR.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
      else if(PAR.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
      ++count_limit1;
    }
    if( count_limit1 > p_limit ){ count_limit1 = p_limit; }
  }
  if( PAR.Nshells != 0 ){
    for(int chan = 0; chan < Chan.size3; ++chan){
      if( type == 1 ){ num = Chan.np[chan]; }
      else{ num = Chan.nh[chan]; } // type = 0
      if( num == 0 || Chan.qnums3[chan].t == -1 ){ continue; } // skip if no particles/holes in this channel
      if(PAR.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
      else if(PAR.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
      ++count_limit2;
    }
    if( count_limit2 > n_limit ){ count_limit2 = n_limit; }
  }
  N_states = count_limit1 + count_limit2;   // EOM variable
  chan_vec = new int[N_states];             // EOM variable
  qnums = new State[N_states];              // EOM variable
  del_E = new double[N_states];             // EOM variable
  norm_1p = new double[N_states];           // EOM variable
  nob = new int[N_states];                  // EOM variable
  ob_index = new int[N_states];             // EOM variable
  nthb = new int[N_states];                 // EOM variable
  thb_index = new int[N_states];            // EOM variable
  nstate = new int[N_states];               // EOM variable
  state_index = new int[N_states];          // EOM variable

  if( PAR.Pshells != 0 ){
    while( count1 < count_limit1 ){
      en = limit;
      ind = -1;
      for(int chan = 0; chan < Chan.size3; ++chan){
	if( type == 1 ){ num = Chan.np[chan]; index = 0; }
	else{ num = Chan.nh[chan]; index = num-1; } // type = 0
	if( num == 0 || Chan.qnums3[chan].t == 1 ){ continue; } // skip if no particles/holes in this channel
	if(PAR.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
        else if(PAR.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
	for(int n = 0; n < count1; ++n){ if( chan_vec[n] == chan ){ goto stop1; } }
	if( type == 1 ){ tempen = SPB.qnums[Chan.p_state(chan, index).v1].energy; }
	else{ tempen = SPB.qnums[Chan.h_state(chan, index).v1].energy; } // type = 0
	if( (type == 1 && tempen < en) || (type == 0 && tempen > en) ){
	  en = tempen;
	  ind = chan;
	}
      stop1:;
      }
      chan_vec[count1] = ind;
      qnums[count1] = Chan.qnums3[ind];
      ++count1;
    }
  }
  if( PAR.Nshells != 0 ){
    while( count2 < count_limit2 ){
      en = limit;
      ind = -1;
      for(int chan = 0; chan < Chan.size3; ++chan){
	if( type == 1 ){ num = Chan.np[chan]; index = 0; }
	else{ num = Chan.nh[chan]; index = num-1; } // type = 0
	if( num == 0 || Chan.qnums3[chan].t == -1 ){ continue; } // skip if no particles/holes in this channel
	if(PAR.basis == "finite_HO" && (Chan.qnums3[chan].m != 1 || Chan.qnums3[chan].ml < 0 || Chan.qnums3[chan].ml > 2)){ continue; }
        else if(PAR.basis == "finite_M" && (Chan.qnums3[chan].m != 1)){ continue; }
	for(int n = count1; n < count1+count2; ++n){ if( chan_vec[n] == chan ){ goto stop2; } }
	if( type == 1 ){ tempen = SPB.qnums[Chan.p_state(chan, index).v1].energy; }
	else{ tempen = SPB.qnums[Chan.h_state(chan, index).v1].energy; } // type = 0
	if( (type == 1 && tempen < en) || (type == 0 && tempen > en) ){
	  en = tempen;
	  ind = chan;
	}
      stop2:;
      }
      chan_vec[count1 + count2] = ind;
      qnums[count1 + count2] = Chan.qnums3[ind];
      ++count2;
    }
  }
}

void EOM::count_states2(Channels &Chan)
{
  if( PAR.Pshells == 0 || PAR.Nshells == 0 ){ std::cerr << "No ph-channels" << std::endl; exit(1); }
  int count1 = 0, count2 = 0, count_limit1 = 0, count_limit2 = 0;
  int num, p, h, ind, index;
  double en, tempen;
  int p_limit = 6;
  int n_limit = 0;
  double limit = 1.0e10;

  double *energies = new double[Chan.size2];
  for(int chan = 0; chan < Chan.size2; ++chan){
    energies[chan] = limit;
    num = Chan.nph1[chan];
    if( num == 0 ){ continue; } // skip if no ph-states in this channel
    if(PAR.basis == "finite_M" && (Chan.qnums2[chan].m != 0)){ continue; }
    if( Chan.qnums2[chan].t == -2 ){ ++count_limit1; } // +P-N
    if( Chan.qnums2[chan].t == 2 ){ ++count_limit2; }  // +N-P
    for(int ph1 = 0; ph1 < num; ++ph1){
      p = Chan.ph1_state(chan, ph1).v1;
      h = Chan.ph1_state(chan, ph1).v2;
      if( SPB.qnums[p].energy - SPB.qnums[h].energy < energies[chan] ){
	energies[chan] = SPB.qnums[p].energy - SPB.qnums[h].energy;
      }
    }
  }
  if( count_limit1 > p_limit ){ count_limit1 = p_limit; }
  if( count_limit2 > n_limit ){ count_limit2 = n_limit; }

  N_states = count_limit1 + count_limit2;   // EOM variable
  chan_vec = new int[N_states];             // EOM variable
  qnums = new State[N_states];              // EOM variable
  del_E = new double[N_states];             // EOM variable
  norm_1p = new double[N_states];           // EOM variable
  ntb = new int[N_states];                  // EOM variable
  tb_index = new int[N_states];             // EOM variable
  nfb = new int[N_states];                  // EOM variable
  fb_index = new int[N_states];             // EOM variable
  nstate = new int[N_states];               // EOM variable
  state_index = new int[N_states];          // EOM variable
  
  count1 = 0;
  while( count1 < count_limit1 ){
    en = limit;
    ind = -1;
    for(int chan = 0; chan < Chan.size2; ++chan){
      num = Chan.nph1[chan];
      if( num == 0 ){ continue; } // skip if no ph-states in this channel
      if(PAR.basis == "finite_M" && (Chan.qnums2[chan].m != 0)){ continue; }
      if( Chan.qnums2[chan].t == -2 ){ // +P-N
	if( energies[chan] < en ){
	  en = energies[chan];
	  ind = chan;
	}
      }
    }
    energies[ind] = limit;
    chan_vec[count1] = ind;
    qnums[count1] = Chan.qnums2[ind];
    ++count1;
  }
  
  count2 = 0;
  while( count2 < count_limit2 ){
    en = limit;
    ind = -1;
    for(int chan = 0; chan < Chan.size2; ++chan){
      num = Chan.nph1[chan];
      if( num == 0 ){ continue; } // skip if no ph-states in this channel
      if(PAR.basis == "finite_M" && (Chan.qnums2[chan].m != 0)){ continue; }
      if( Chan.qnums2[chan].t == 2 ){ // +N-P
	if( energies[chan] < en ){
	  en = energies[chan];
	  ind = chan;
	}
      }
    }
    energies[ind] = limit;
    chan_vec[count1 + count2] = ind;
    qnums[count1 + count2] = Chan.qnums2[ind];
    ++count2;
  }
}

void EOM::Print_EOM_1P(double Energy)
{
  std::cout << std::fixed;
  if(PAR.basis == "finite_HO"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << PAR.Shells << std::setw(5) << PAR.Pshells << std::setw(5) << qnums[i].ml;
      std::cout << std::setw(5) << qnums[i].m << std::setprecision(2) << std::setw(8) << PAR.density;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
  else if(PAR.basis == "finite_J"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << PAR.Shells << std::setw(5) << PAR.Pshells << std::setw(5) << PAR.Nshells;
      std::cout << std::setw(5) << qnums[i].j << std::setw(5) << qnums[i].par << std::setw(5) << qnums[i].t;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
  else if(PAR.basis == "finite_M"){
    for(int i = 0; i < N_states; ++i){
      std::cout << std::setw(5) << PAR.Shells << std::setw(5) << PAR.Pshells << std::setw(5) << PAR.Nshells;
      std::cout << std::setw(5) << qnums[i].m << std::setw(5) << qnums[i].par << std::setw(5) << qnums[i].t;
      std::cout << std::setprecision(9) << std::setw(17) << Energy << std::setw(17) << del_E[i];
      std::cout << std::setw(17) << Energy + del_E[i] << std::setw(17) << norm_1p[i] << std::endl;
    }
  }
}

void EOM::One_Body_Transition_1PA(Channels &Chan, OP_1b &OP)
{
  int ind, ind1, ind2, length;
  int nob1, nob2, nthb1, nthb2, nstate1, nstate2;
  int ob_ind1, ob_ind2, thb_ind1, thb_ind2, state_ind1, state_ind2;
  int chan3_1, chan3_2, chan1_j, chan2_j, lambda = 2;
  double Operator;
  for(ind1 = 0; ind1 < N_states; ++ind1){
    for(ind2 = 0; ind2 < N_states; ++ind2){
      for(ind = 0; ind < OP.chan3_length; ++ind){
	chan3_1 = chan_vec[ind1];
	chan3_2 = chan_vec[ind2];
	if(chan3_1 != OP.chan3_vec1[ind] || chan3_2 != OP.chan3_vec2[ind]){ continue; }
	chan1_j = Chan.qnums3[chan3_1].j;
	nob1 = nob[ind1];
	nthb1 = nthb[ind1];
	nstate1 = nstate[ind1];
	ob_ind1 = ob_index[ind1];
	thb_ind1 = thb_index[ind1];
	state_ind1 = state_index[ind1];
	chan2_j = Chan.qnums3[chan3_2].j;
	nob2 = nob[ind2];
	nthb2 = nthb[ind2];
	nstate2 = nstate[ind2];
	ob_ind2 = ob_index[ind2];
	thb_ind2 = thb_index[ind2];
	state_ind2 = state_index[ind2];

	Operator = 0.0;

        #pragma omp parallel
	{
	  double op;
	  int b1, b2, b3, k1, k2, k3, index;
	  int bra, ket, key1, key2;
	  State tb1, tb2;
	  double fac2;//, X;
	  
	  length = nstate1 * nstate2;
          #pragma omp for schedule(static)
	  for(int braket = 0; braket < length; ++braket){
	    ket = int(braket%nstate2);
	    bra = int((braket - ket)/nstate2);
	    fac2 = 1.0;
	    op = 0.0;
	    
	    if( bra >= nob1 ){
	      bra -= nob1;
	      b1 = thb_vec[thb_index[ind1] + bra].v1;
	      b2 = thb_vec[thb_index[ind1] + bra].v2;
	      b3 = thb_vec[thb_index[ind1] + bra].v3;
	      if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	      if( ket >= nob2 ){  //  < pph | pph >
		ket -= nob2;
		k1 = thb_vec[thb_index[ind2] + ket].v1;
		k2 = thb_vec[thb_index[ind2] + ket].v2;
		k3 = thb_vec[thb_index[ind2] + ket].v3;
		if(k1 == k2){ fac2 /= std::sqrt(2.0); }

		/*if( b3 == k3 ){
		  if( b1 == k1 ){
		    key1 = Chan.p_map[chan3_1][b2];
		    key2 = Chan.p_map[chan3_2][k2];
		    index = key1 * Chan.np[chan3_2] + key2;
		    if( PAR.basis == "finite_J" ){
		      X = std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0) * (J_vec[bra] + 1.0) * (J_vec[ket] + 1.0) * (lambda + 1.0));
		      X *= CGC6(J_vec[bra],lambda,J_vec[ket], b_vec[ket],b_vec[bra],a_vec[bra]) * CGC6(J_vec[bra],lambda,J_vec[ket], chan2_j,i_vec[bra],chan1_j);
		      X *= phase2(chan2_j - a_vec[ket] + b_vec[ket] - i_vec[ket]);
		      op += OP.pp_3[OP.pp_3_index[ind] + index] * X;
		    }
		    else if( PAR.basis == "finite_M" ){ op += OP.pp_3[OP.pp_3_index[ind] + index]; }
		  }
		  if( b1 == k2 ){
		    key1 = Chan.p_map[chan3_1][b2];
		    key2 = Chan.p_map[chan3_2][k1];
		    index = key1 * Chan.np[chan3_2] + key2;
		    if( PAR.basis == "finite_J" ){
		      X = std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0) * (J_vec[bra] + 1.0) * (J_vec[ket] + 1.0) * (lambda + 1.0));
		      X *= CGC6(J_vec[bra],lambda,J_vec[ket], a_vec[ket],b_vec[bra],a_vec[bra]) * CGC6(J_vec[bra],lambda,J_vec[ket], chan2_j,i_vec[bra],chan1_j);
		      X *= phase2(chan2_j - i_vec[ket] - J_vec[ket]);
		      op += OP.pp_3[OP.pp_3_index[ind] + index] * X;
		    }
		    else if( PAR.basis == "finite_M" ){ op -= OP.pp_3[OP.pp_3_index[ind] + index]; }
		  }
		  if( b2 == k2 ){
		    key1 = Chan.p_map[chan3_1][b1];
		    key2 = Chan.p_map[chan3_2][k1];
		    index = key1 * Chan.np[chan3_2] + key2;
		    if( PAR.basis == "finite_J" ){
		      X = std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0) * (J_vec[bra] + 1.0) * (J_vec[ket] + 1.0) * (lambda + 1.0));
		      X *= CGC6(J_vec[bra],lambda,J_vec[ket], a_vec[ket],a_vec[bra],b_vec[bra]) * CGC6(J_vec[bra],lambda,J_vec[ket], chan2_j,i_vec[bra],chan1_j);
		      X *= -1.0 * phase2(chan2_j + a_vec[bra] + b_vec[bra] - i_vec[bra] - J_vec[bra] - J_vec[ket]);
		      op += OP.pp_3[OP.pp_3_index[ind] + index] * X;
		    }
		    else if( PAR.basis == "finite_M" ){ op += OP.pp_3[OP.pp_3_index[ind] + index]; }
		  }
		  if( b2 == k1 ){
		    key1 = Chan.p_map[chan3_1][b1];
		    key2 = Chan.p_map[chan3_2][k2];
		    index = key1 * Chan.np[chan3_2] + key2;
		    if( PAR.basis == "finite_J" ){
		      X = std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0) * (J_vec[bra] + 1.0) * (J_vec[ket] + 1.0) * (lambda + 1.0));
		      X *= CGC6(J_vec[bra],lambda,J_vec[ket], b_vec[ket],a_vec[bra],b_vec[bra]) * CGC6(J_vec[bra],lambda,J_vec[ket], chan2_j,i_vec[bra],chan1_j);
		      X *= phase2(chan2_j + a_vec[bra] + b_vec[ket] - i_vec[bra] - J_vec[bra]);
		      op += OP.pp_3[OP.pp_3_index[ind] + index] * X;
		    }
		    else if( PAR.basis == "finite_M" ){ op -= OP.pp_3[OP.pp_3_index[ind] + index]; }
		  }
		}
		
		if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) ){
		  key1 = Chan.h_map[chan3_2][k3];
		  key2 = Chan.h_map[chan3_1][b3];
		  index = key1 * Chan.nh[chan3_1] + key2;
		  if( PAR.basis == "finite_J" ){
		    X = std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0) * (lambda + 1.0));
		    X *= CGC6(i_vec[bra],lambda,i_vec[ket], chan2_j,J_vec[bra],chan1_j);
		    X *= -1.0 * phase2(i_vec[bra] + chan1_j - J_vec[bra]);
		    op += OP.hh_3[OP.hh_3_index[ind] + index] * X;
		    if( b1 == b2 ){ op += OP.hh_3[OP.hh_3_index[ind] + index] * X; }
		  }
		  else if( PAR.basis == "finite_M" ){
		    op -= OP.hh_3[OP.hh_3_index[ind] + index];
		    if(b1 == b2){ op -= OP.hh_3[OP.hh_3_index[ind] + index]; }
		  }
		  }*/
		ket += nob2;
	      }
	      else{  //  < pph | p >
		k1 = Chan.p_state(chan3_2, ket).v1;

	      }
	      bra += nob1;
	    }
	    else if( ket >= nob2 ){  //  < p | pph >
	      ket -= nob2;
	      k1 = thb_vec[thb_index[ind2] + ket].v1;
	      k2 = thb_vec[thb_index[ind2] + ket].v2;
	      k3 = thb_vec[thb_index[ind2] + ket].v3;
	      b1 = Chan.p_state(chan3_1, bra).v1;
	      if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	      
	    }
	    else{  //  < p | p >
	      b1 = ob_vec[ob_index[ind1] + bra].v1;
	      k1 = ob_vec[ob_index[ind2] + ket].v1;
	      
	      key1 = Chan.p_map[chan3_1][b1];
	      key2 = Chan.p_map[chan3_2][k1];
	      index = key1 * Chan.np[chan3_2] + key2;
	      if( PAR.basis == "finite_J" ){
		op += OP.pp_3[OP.pp_3_index[ind] + index] * (lambda + 1.0);
	      }
	      else if( PAR.basis == "finite_M" ){
		op += OP.pp_3[OP.pp_3_index[ind] + index];
	      }
	    }
	    op *= fac2 * state_vec_R[state_index[ind1] + bra] * state_vec_L[state_index[ind2] + ket];
	    if(std::fabs(op) > 1e-10){
	      if( PAR.basis == "finite_J" ){ std::cout << "1: < " << b1 << " |o| " << k1 << " > = " << op << std::endl; }
	      else if( PAR.basis == "finite_M" ){ std::cout << "2: < " << b1 << " |o| " << k1 << " > = " << op << std::endl; }
	    }
	    Operator += op;
	  }
	}
	std::cout << "Operator = " << Operator << std::endl;
      }
    }
  }
}

void EOM::One_Body_Transition_1PR(Channels &Chan, OP_1b &OP)
{
  int ind, ind1, ind2, length;
  int nob1, nob2, nthb1, nthb2, nstate1, nstate2;
  int ob_ind1, ob_ind2, thb_ind1, thb_ind2, state_ind1, state_ind2;
  int chan3_1, chan3_2, chan1_j, chan2_j, lambda = 2;
  double Operator;
  for(ind1 = 0; ind1 < N_states; ++ind1){
    for(ind2 = 0; ind2 < N_states; ++ind2){
      for(ind = 0; ind < OP.chan3_length; ++ind){
	chan3_1 = chan_vec[ind1];
	chan3_2 = chan_vec[ind2];
	if(chan3_1 != OP.chan3_vec1[ind] || chan3_2 != OP.chan3_vec2[ind]){ continue; }
	chan1_j = Chan.qnums3[chan3_1].j;
	nob1 = nob[ind1];
	nthb1 = nthb[ind1];
	nstate1 = nstate[ind1];
	ob_ind1 = ob_index[ind1];
	thb_ind1 = thb_index[ind1];
	state_ind1 = state_index[ind1];
	chan2_j = Chan.qnums3[chan3_2].j;
	nob2 = nob[ind2];
	nthb2 = nthb[ind2];
	nstate2 = nstate[ind2];
	ob_ind2 = ob_index[ind2];
	thb_ind2 = thb_index[ind2];
	state_ind2 = state_index[ind2];

	Operator = 0.0;

        #pragma omp parallel
	{
	  double op;
	  int b1, b2, b3, k1, k2, k3, index;
	  int bra, ket, key1, key2;
	  State tb1, tb2;
	  double fac2;
	  
	  length = nstate1 * nstate2;
          #pragma omp for schedule(static)
	  for(int braket = 0; braket < length; ++braket){
	    ket = int(braket%nstate2);
	    bra = int((braket - ket)/nstate2);
	    fac2 = 1.0;
	    op = 0.0;
	    
	    if( bra >= nob1 ){
	      bra -= nob1;
	      b1 = thb_vec[thb_index[ind1] + bra].v1;
	      b2 = thb_vec[thb_index[ind1] + bra].v2;
	      b3 = thb_vec[thb_index[ind1] + bra].v3;
	      if(b1 == b2){ fac2 /= std::sqrt(2.0); }
	      if( ket >= nob2 ){  //  < hhp | hhp >
		ket -= nob2;
		k1 = thb_vec[thb_index[ind2] + ket].v1;
		k2 = thb_vec[thb_index[ind2] + ket].v2;
		k3 = thb_vec[thb_index[ind2] + ket].v3;
		if(k1 == k2){ fac2 /= std::sqrt(2.0); }
		
		/*if(b3 == k3 && equal(J_state[bra], J_state[ket])){
		  if(b1 == k1 && equal(SPB.qnums[b2], SPB.qnums[k2])){
		    minus(tb1, SPB.qnums[k2], SPB.qnums[b2]);
		    ind = Chan.hh1_map[Chan.ind0][Hash(k2, b2, 0)];
		    ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(j_vec[bra] + 1.0);
		  }
		  if(b1 == k2 && equal(SPB.qnums[b2], SPB.qnums[k1])){
		    minus(tb1, SPB.qnums[k1], SPB.qnums[k2]);
		    ind = Chan.hh1_map[Chan.ind0][Hash(k1, b2, 0)];
		    ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(j_vec[bra] + 1.0);
		  }
		  if(b2 == k2 && equal(SPB.qnums[b1], SPB.qnums[k1])){
		    minus(tb1, SPB.qnums[k1], SPB.qnums[b1]);
		    ind = Chan.hh1_map[Chan.ind0][Hash(k1, b1, 0)];
		    ME -= Eff_Ints.Xhh.X_2[ind] / std::sqrt(i_vec[bra] + 1.0);
		  }
		  if(b2 == k1 && equal(SPB.qnums[b1], SPB.qnums[k2])){
		    minus(tb1, SPB.qnums[k2], SPB.qnums[b1]);
		    ind = Chan.hh1_map[Chan.ind0][Hash(k2, b1, 0)];
		    ME += Eff_Ints.Xhh.X_2[ind] * phase2(i_vec[bra] + j_vec[bra] - J_vec[bra]) / std::sqrt(i_vec[bra] + 1.0);
		  }
		}
		
		if(b1 == k1 && b2 == k2 && equal(J_state[bra], J_state[ket]) && equal(SPB.qnums[b3], SPB.qnums[k3])){
		  minus(tb1, SPB.qnums[b3], SPB.qnums[k3]);
		  ind = Chan.pp1_map[Chan.ind0][Hash(b3, k3, 0)];
		  ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0);
		  if(b1 == b2){ ME += Eff_Ints.Xpp.X_2[ind] / std::sqrt(a_vec[ket] + 1.0); }
		  }*/
		ket += nob2;
	      }
	      else{  //  < hhp | h >
		k1 = ob_vec[ob_index[ind2] + ket].v1;

	      }
	      bra += nob1;
	    }
	    else if( ket >= nob2 ){  //  < h | hhp >
	      ket -= nob2;
	      k1 = thb_vec[thb_index[ind2] + ket].v1;
	      k2 = thb_vec[thb_index[ind2] + ket].v2;
	      k3 = thb_vec[thb_index[ind2] + ket].v3;
	      b1 = ob_vec[ob_index[ind1] + bra].v1;
	      if(k1 == k2){ fac2 /= std::sqrt(2.0); }
	    
	      /*if(k1 == b1 && equal(SPB.qnums[k2], SPB.qnums[k3])){
		minus(tb1, SPB.qnums[k2], SPB.qnums[k3]);
		ind = Chan.hp1_map[Chan.ind0][Hash(k2, k3, 0)];
		X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (j_vec[ket] + 1.0))) * phase2(i_vec[ket] + j_vec[ket] - J_vec[ket]);
		ME += Eff_Ints.Xhp.X_2[ind] * X;
	      }
	      if(k2 == b1 && equal(SPB.qnums[k1], SPB.qnums[k3])){
		minus(tb1, SPB.qnums[k1], SPB.qnums[k3]);
		ind = Chan.hp1_map[Chan.ind0][Hash(k1, k3, 0)];
		X = std::sqrt((J_vec[ket] + 1.0) / ((chan_j + 1.0) * (i_vec[ket] + 1.0)));
		ME -= Eff_Ints.Xhp.X_2[ind] * X;
		}*/

	      ket += nob2;
	    }
	    else{  //  < h | h >
	      b1 = ob_vec[ob_index[ind1] + bra].v1;
	      k1 = ob_vec[ob_index[ind2] + ket].v1;

	      key1 = Chan.h_map[chan3_1][b1];
	      key2 = Chan.h_map[chan3_2][k1];
	      index = key1 * Chan.nh[chan3_2] + key2;
	      op += phase2(lambda + chan2_j - chan1_j) * OP.hh_3[OP.hh_3_index[ind] + index] / std::sqrt((chan1_j + 1.0) * (chan2_j + 1.0));
	    }
	    op *= fac2 * state_vec_R[state_index[ind1] + bra] * state_vec_L[state_index[ind2] + ket];
	    Operator += op;
	  }
	}
      }
    }
  }
}

void EOM::Delete()
{
  if(N_states != 0){
    delete[] qnums;
    delete[] chan_vec;
    delete[] del_E;
    delete[] norm_1p;
    delete[] nob;
    //delete[] ob_vec;
    delete[] ob_index;
    delete[] nthb;
    //delete[] thb_vec;
    delete[] thb_index;
    //delete[] thb_qnums;
    delete[] nstate;
    //delete[] state_vec_R;
    //delete[] state_vec_L;
    delete[] state_index;
  }
}

bool Fermi_0(Channels &Chan, int chan3_1, int chan3_2)
{
  return (Chan.qnums3[chan3_1].t != Chan.qnums3[chan3_2].t && Chan.qnums3[chan3_1].par == Chan.qnums3[chan3_2].par &&
	  ((PAR.basis == "finite_J" && Chan.qnums3[chan3_1].j == Chan.qnums3[chan3_2].j) ||
	   (PAR.basis == "finite_M" && Chan.qnums3[chan3_1].m == Chan.qnums3[chan3_2].m)));
}

double Fermi_1(int v1, int v2)
{
  if(SPB.qnums[v1].n != SPB.qnums[v2].n || SPB.qnums[v1].l != SPB.qnums[v2].l){ return 0.0; }
  else if(PAR.basis == "finite_J"){ return std::sqrt(SPB.qnums[v1].j + 1.0); }
  else if(PAR.basis == "finite_M" && SPB.shellsj[v1] == SPB.shellsj[v2]){ return std::sqrt(SPB.shellsj[v1] + 1.0); }
  else{ return 0.0; }
}

bool Gamow_Teller_0(Channels &Chan, int chan3_1, int chan3_2)
{
  return (Chan.qnums3[chan3_1].t != Chan.qnums3[chan3_2].t && Chan.qnums3[chan3_1].par == Chan.qnums3[chan3_2].par &&
	  ((PAR.basis == "finite_J" && std::abs(Chan.qnums3[chan3_1].j - Chan.qnums3[chan3_2].j) <= 2) ||
	   (PAR.basis == "finite_M" && std::abs(Chan.qnums3[chan3_1].m - Chan.qnums3[chan3_2].m) <= 2)));
}

double Gamow_Teller_1(int v1, int v2)
{
  double term;
  int l1 = SPB.qnums[v1].l;
  if(SPB.qnums[v1].n != SPB.qnums[v2].n || l1 != SPB.qnums[v2].l){ return 0.0; }
  else if(PAR.basis == "finite_J"){
    term = std::sqrt(2.0*(SPB.qnums[v1].j + 1.0)*(SPB.qnums[v2].j + 1.0)) * phase2(2*l1 + SPB.qnums[v1].j + 3) * CGC6(1,1,2, SPB.qnums[v2].j,SPB.qnums[v1].j,2*l1);
    /*if( v1 == 4 && v2 == 5 ){
      std::cout << "!! " << v1 << " " << v2 << " : " << std::sqrt(2.0*(SPB.qnums[v1].j + 1.0)*(SPB.qnums[v2].j + 1.0)) << " " << phase2(2*l1 + SPB.qnums[v1].j + 3) << " " << CGC6(1,1,2, SPB.qnums[v2].j,SPB.qnums[v1].j,2*l1) << " : " << SPB.qnums[v1].j << " " << SPB.qnums[v2].j << " " << l1 << std::endl;
      }*/
    return term;
  }
  else if(PAR.basis == "finite_M" && SPB.shellsj[v1] == SPB.shellsj[v2]){
    term = std::sqrt(2.0*(SPB.shellsj[v1] + 1.0)*(SPB.shellsj[v2] + 1.0)) * phase2(2*l1 + SPB.shellsj[v1] + 3) * CGC6(1,1,2, SPB.shellsj[v2],SPB.shellsj[v1],2*l1);
    /*if( v1 == 4 && v2 == 8 ){
      std::cout << "!! " << v1 << " " << v2 << " : " << term << " : CGC( " << SPB.shellsj[v2] << "," << SPB.qnums[v2].m << "," << 2 << "," << SPB.qnums[v1].m-SPB.qnums[v2].m << "," << SPB.shellsj[v1] << "," << SPB.qnums[v1].m << ") = " << CGC(SPB.shellsj[v2],SPB.qnums[v2].m, 2,SPB.qnums[v1].m-SPB.qnums[v2].m, SPB.shellsj[v1],SPB.qnums[v1].m) << " : / " << std::sqrt(SPB.shellsj[v1] + 1.0) << std::endl;
      }*/
    return term * CGC(SPB.shellsj[v2],SPB.qnums[v2].m, 2,SPB.qnums[v1].m-SPB.qnums[v2].m, SPB.shellsj[v1],SPB.qnums[v1].m) / std::sqrt(SPB.shellsj[v1] + 1.0);
  }
  else{ return 0.0; }
}

OP_1b::OP_1b(Channels &Chan, Amplitudes &Amps, int type)
{
  int length, ind, ind0, chan3_1, chan3_2;
  int h1, h2, p1, p2;
  double term = 0.0;
  int length2, chan2, *chan2_vec0 = new int[Chan.size2];
  State tb2;
  int lambda;
  if( type == 0 ){ lambda = 0; }
  else{ lambda = 2; } // type = 1
  //////////
  int nh1, np1, nh2, np2, nph, nhp, chan_ind1, chan_ind2, chan_ind3, one = 1;
  double plus1 = 1.0, minus1 = -1.0;
  double *T2, *t3;
  char N = 'N';

  /*if(PAR.basis == "finite_M"){
    for(ind = 0; ind < SPB.num_states; ++ind){
      std::cout << "! " << ind << " : " << SPB.shellsj[ind] << std::endl;
    }
    }*/

  length = 0;
  length2 = 0;
  for(chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
    for(chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){	
      if((type == 0 && Fermi_0(Chan, chan3_1, chan3_2)) ||
	 (type == 1 && Gamow_Teller_0(Chan, chan3_1, chan3_2))){
	++length;

	minus(tb2, Chan.qnums3[chan3_1], Chan.qnums3[chan3_2]);
	if( type == 0 ){ tb2.j = lambda; }
	else if( type == 1 ){ tb2.j = lambda; }
	chan2 = Ind_Chan2(tb2);
	for(int i = 0; i < length2; ++i){ if(chan2 == chan2_vec0[i]){ goto stop; } }
	chan2_vec0[length2] = chan2;
	++length2;
      stop:;
	tb2.j -= 2;
      }
    }
  }

  // set chan3 vecs
  chan3_length = length;
  chan3_vec1 = new int[length];
  chan3_vec2 = new int[length];
  length = 0;
  for(chan3_1 = 0; chan3_1 < Chan.size3; ++chan3_1){
    for(chan3_2 = 0; chan3_2 < Chan.size3; ++chan3_2){
      if((type == 0 && Fermi_0(Chan, chan3_1, chan3_2)) ||
	 (type == 1 && Gamow_Teller_0(Chan, chan3_1, chan3_2))){
	chan3_vec1[length] = chan3_1;
	chan3_vec2[length] = chan3_2;
	++length;
      }
    }
  }
  // set chan2 vec
  chan2_length = length2;
  chan2_vec = new int[length2];
  for(ind = 0; ind < length2; ++ind){ chan2_vec[ind] = chan2_vec0[ind]; }
  delete[] chan2_vec0;

  // setup O_hh
  length = 0;
  hh_3_index = new int[chan3_length];
  for(ind = 0; ind < chan3_length; ++ind){
    hh_3_index[ind] = length;
    length += Chan.nh[chan3_vec1[ind]] * Chan.nh[chan3_vec2[ind]];
  }
  hh_3 = new double[length];
  hh_map_ind = new int[length];
  hh_3_length = length;
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int hind2 = 0; hind2 < Chan.nh[chan3_2]; ++hind2){
	h2 = Chan.h_vec[Chan.h_index[chan3_2] + hind2].v1;
	if(type == 0){ term = Fermi_1(h1, h2); }
	else if(type == 1){ term = Gamow_Teller_1(h1, h2); }
	hh_3[hh_3_index[ind] + (Chan.nh[chan3_2]*hind1 + hind2)] = term;
	//if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << h1 << " |o| " << h2 << " > = " << term << std::endl; }
      }
    }
  }
  length2 = 0;
  hh_2_index = new int[chan2_length];
  for(ind = 0; ind < chan2_length; ++ind){
    hh_2_index[ind] = length2;
    length2 += Chan.nhh1[chan2_vec[ind]];
  }
  hh_2 = new double[length2];
  hh_2_length = length2;
  for(ind = 0; ind < hh_2_length; ++ind){ hh_2[ind] = 0.0; }
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int hind2 = 0; hind2 < Chan.nh[chan3_2]; ++hind2){
	h2 = Chan.h_vec[Chan.h_index[chan3_2] + hind2].v1;
	ind0 = hh_3_index[ind] + hind1*Chan.nh[chan3_2] + hind2;
	//Map_2_1b(hh_map_ind, chan2_vec,chan2_length,hh_2_index, Chan.hh1_map, ind0,h1,h2, lambda);
      }
    }
  }

  // setup O_pp
  length = 0;
  pp_3_index = new int[chan3_length];
  for(ind = 0; ind < chan3_length; ++ind){
    pp_3_index[ind] = length;
    length += Chan.np[chan3_vec1[ind]] * Chan.np[chan3_vec2[ind]];
  }
  pp_3 = new double[length];
  pp_map_ind = new int[length];
  pp_3_length = length;
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int pind2 = 0; pind2 < Chan.np[chan3_2]; ++pind2){
	p2 = Chan.p_vec[Chan.p_index[chan3_2] + pind2].v1;
	if(type == 0){ term = Fermi_1(p1, p2); }
	else if(type == 1){ term = Gamow_Teller_1(p1, p2); }
	pp_3[pp_3_index[ind] + (Chan.np[chan3_2]*pind1 + pind2)] = term;
	//if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << p1 << " |o| " << p2 << " > = " << term << std::endl; }
      }
    }
  }
  length2 = 0;
  pp_2_index = new int[chan2_length];
  for(ind = 0; ind < chan2_length; ++ind){
    pp_2_index[ind] = length2;
    length2 += Chan.npp1[chan2_vec[ind]];
  }
  pp_2 = new double[length2];
  pp_2_length = length2;
  for(ind = 0; ind < pp_2_length; ++ind){ pp_2[ind] = 0.0; }
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int pind2 = 0; pind2 < Chan.np[chan3_2]; ++pind2){
	p2 = Chan.p_vec[Chan.p_index[chan3_2] + pind2].v1;
	ind0 = pp_3_index[ind] + pind1*Chan.np[chan3_2] + pind2;
	//Map_2_1b(pp_map_ind, chan2_vec,chan2_length,pp_2_index, Chan.pp1_map, ind0,p1,p2, lambda);
      }
    }
  }

  // setup O_ph
  length = 0;
  ph_3_index = new int[chan3_length];
  for(ind = 0; ind < chan3_length; ++ind){
    ph_3_index[ind] = length;
    length += Chan.np[chan3_vec1[ind]] * Chan.nh[chan3_vec2[ind]];
  }
  ph_3 = new double[length];
  ph_map_ind = new int[length];
  ph_3_length = length;
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int hind1 = 0; hind1 < Chan.nh[chan3_2]; ++hind1){
	h1 = Chan.h_vec[Chan.h_index[chan3_2] + hind1].v1;
	if(type == 0){ term = Fermi_1(p1, h1); }
	else if(type == 1){ term = Gamow_Teller_1(p1, h1); }
	ph_3[ph_3_index[ind] + (Chan.nh[chan3_2]*pind1 + hind1)] = term;
	//if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << p1 << " |o| " << h1 << " > = " << term << std::endl; }
      }
    }
  }
  length2 = 0;
  ph_2_index = new int[chan2_length];
  for(ind = 0; ind < chan2_length; ++ind){
    ph_2_index[ind] = length2;
    length2 += Chan.nph1[chan2_vec[ind]];
  }
  ph_2 = new double[length2];
  ph_2_length = length2;
  for(ind = 0; ind < ph_2_length; ++ind){ ph_2[ind] = 0.0; }
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int hind1 = 0; hind1 < Chan.nh[chan3_2]; ++hind1){
	h1 = Chan.h_vec[Chan.h_index[chan3_2] + hind1].v1;
	ind0 = ph_3_index[ind] + pind1*Chan.nh[chan3_2] + hind1;
	//Map_2_1b(ph_map_ind, chan2_vec,chan2_length,ph_2_index, Chan.ph1_map, ind0,p1,h1, lambda);
      }
    }
  }

  // setup O_hp
  length = 0;
  hp_3_index = new int[chan3_length];
  for(ind = 0; ind < chan3_length; ++ind){
    hp_3_index[ind] = length;
    length += Chan.nh[chan3_vec1[ind]] * Chan.np[chan3_vec2[ind]];
  }
  hp_3 = new double[length];
  hp_map_ind = new int[length];
  hp_3_length = length;
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int pind1 = 0; pind1 < Chan.np[chan3_2]; ++pind1){
	p1 = Chan.p_vec[Chan.p_index[chan3_2] + pind1].v1;
	if(type == 0){ term = Fermi_1(h1, p1); }
	else if(type == 1){ term = Gamow_Teller_1(h1, p1); }
	hp_3[hp_3_index[ind] + (Chan.np[chan3_2]*hind1 + pind1)] = term;
	//if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << h1 << " |o| " << p1 << " > = " << term << std::endl; }
      }
    }
  }
  length2 = 0;
  hp_2_index = new int[chan2_length];
  for(ind = 0; ind < chan2_length; ++ind){
    hp_2_index[ind] = length2;
    length2 += Chan.nhp1[chan2_vec[ind]];
  }
  hp_2 = new double[length2];
  hp_2_length = length2;
  for(ind = 0; ind < hp_2_length; ++ind){ hp_2[ind] = 0.0; }
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int pind1 = 0; pind1 < Chan.np[chan3_2]; ++pind1){
	p1 = Chan.p_vec[Chan.p_index[chan3_2] + pind1].v1;
	ind0 = hp_3_index[ind] + hind1*Chan.np[chan3_2] + pind1;
	//Map_2_1b(hp_map_ind, chan2_vec,chan2_length,hp_2_index, Chan.hp1_map, ind0,h1,p1, lambda);
      }
    }
  }

  /*for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int hind2 = 0; hind2 < Chan.nh[chan3_2]; ++hind2){
	h2 = Chan.h_vec[Chan.h_index[chan3_2] + hind2].v1;
	term = hh_3[hh_3_index[ind] + (Chan.nh[chan3_2]*hind1 + hind2)];
	if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << h1 << " | o | " << h2 << " > = " << term << std::endl; }
      }
    }
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int pind2 = 0; pind2 < Chan.np[chan3_2]; ++pind2){
	p2 = Chan.p_vec[Chan.p_index[chan3_2] + pind2].v1;
	term = pp_3[pp_3_index[ind] + (Chan.np[chan3_2]*pind1 + pind2)];
	if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << p1 << " | o | " << p2 << " > = " << term << std::endl; }
      }
    }
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int hind1 = 0; hind1 < Chan.nh[chan3_2]; ++hind1){
	h1 = Chan.h_vec[Chan.h_index[chan3_2] + hind1].v1;
	term = ph_3[ph_3_index[ind] + (Chan.nh[chan3_2]*pind1 + hind1)];
	if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << p1 << " | o | " << h1 << " > = " << term << std::endl; }
      }
    }
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int pind1 = 0; pind1 < Chan.np[chan3_2]; ++pind1){
	p1 = Chan.p_vec[Chan.p_index[chan3_2] + pind1].v1;
	term = hp_3[hp_3_index[ind] + (Chan.np[chan3_2]*hind1 + pind1)];
	if( std::fabs(term) > 1.0e-8 ){ std::cout << "< " << h1 << " | o | " << p1 << " > = " << term << std::endl; }
      }
    }
    }*/

  // Similarity transformed operator
  for(ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    nh1 = Chan.nh[chan3_1];
    np1 = Chan.np[chan3_1];
    nh2 = Chan.nh[chan3_2];
    np2 = Chan.np[chan3_2];

    chan_ind1 = Amps.S1.t3_index[chan3_1];
    chan_ind2 = hh_3_index[ind];
    chan_ind3 = ph_3_index[ind];
    if(np1 * nh2 * nh1 != 0){
      t3 = new double[np1 * nh1];
      Amps.S1.Set_t3(t3, np1*nh1, chan_ind1);
      //OP(a|i){a,i}  <-  -t3(a|k){a,k}.OP_hh3(k|i){k,i}
      dgemm_NN(t3, (hh_3 + chan_ind2), (ph_3 + chan_ind3), &np1, &nh2, &nh1, &minus1, &plus1, &N, &N);
      delete[] t3;
    }

    chan_ind1 = Amps.S1.t3_index[chan3_2];
    chan_ind2 = pp_3_index[ind];
    chan_ind3 = ph_3_index[ind];
    if(np1 * nh2 * np2 != 0){
      t3 = new double[np2 * nh2];
      Amps.S1.Set_t3(t3, np2*nh2, chan_ind1);
      //OP(a|i){a,i}  <-  OP_pp3(a|c){a,c}.t3(c|i){c,i}
      dgemm_NN((pp_3 + chan_ind2), t3, (ph_3 + chan_ind3), &np1, &nh2, &np2, &plus1, &plus1, &N, &N);
      delete[] t3;
    }
  }
  //gather1(Chan);
  for(ind = 0; ind < chan2_length; ++ind){
    chan2 = chan2_vec[ind];
    nph = Chan.nph1[chan2];
    nhp = Chan.nhp1[chan2];
    chan_ind1 = Amps.D1.T2_3_index[chan2];
    chan_ind2 = hp_2_index[ind];
    chan_ind3 = ph_2_index[ind];
    if(nph * nhp != 0){
      T2 = new double[nph * nhp];
      Amps.D1.Set_T2_3(T2, nph*nhp, chan_ind1);
      //OP(a|i){ai'}  <-  + T2_3(ac|ik){ai',kc'}.OP_hp2(k|c){kc'}
      dgemm_NN(T2, (hp_2 + chan_ind2), (ph_2 + chan_ind3), &nph, &one, &nhp, &plus1, &plus1, &N, &N);

      delete[] T2;
    }
  }
  //gather2(Chan);

}

void OP_1b::Gather1(Channels &Chan)
{
  int ind0, h1, p1, chan3_1, chan3_2;
  double x;
  for(int ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int pind1 = 0; pind1 < Chan.np[chan3_2]; ++pind1){
	p1 = Chan.p_vec[Chan.p_index[chan3_2] + pind1].v1;
	ind0 = hp_3_index[ind] + hind1*Chan.np[chan3_2] + pind1;
	x = hp_3[ind0];
	hp_3[ind0] = 0.0;
	hp_2[hp_map_ind[ind0]] = x;
      }
    }
  }
}

void OP_1b::Gather2(Channels &Chan)
{
  int ind0, h1, h2, p1, p2, chan3_1, chan3_2;
  double x;
  for(int ind = 0; ind < chan3_length; ++ind){
    chan3_1 = chan3_vec1[ind];
    chan3_2 = chan3_vec2[ind];
    for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int hind2 = 0; hind2 < Chan.nh[chan3_2]; ++hind2){
	h2 = Chan.h_vec[Chan.h_index[chan3_2] + hind2].v1;
	ind0 = hh_3_index[ind] + hind1*Chan.nh[chan3_2] + hind2;
	x = hh_3[ind0];
	//x += hh_2[hh_map_ind[ind0]];
	hh_3[ind0] = x;
	//hh_2[hh_map_ind[ind0]] = x;
      }
    }
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int pind2 = 0; pind2 < Chan.np[chan3_2]; ++pind2){
	p2 = Chan.p_vec[Chan.p_index[chan3_2] + pind2].v1;
	ind0 = pp_3_index[ind] + pind1*Chan.np[chan3_2] + pind2;
	x = pp_3[ind0];
	//x += pp_2[pp_map_ind[ind0]];
	pp_3[ind0] = x;
	//pp_2[pp_map_ind[ind0]] = x;
      }
    }
    for(int pind1 = 0; pind1 < Chan.np[chan3_1]; ++pind1){
      p1 = Chan.p_vec[Chan.p_index[chan3_1] + pind1].v1;
      for(int hind1 = 0; hind1 < Chan.nh[chan3_2]; ++hind1){
	h1 = Chan.h_vec[Chan.h_index[chan3_2] + hind1].v1;
	ind0 = ph_3_index[ind] + pind1*Chan.nh[chan3_2] + hind1;
	x = ph_3[ind0];
	//x += ph_2[ph_map_ind[ind0]];
	ph_3[ind0] = x;
	//ph_2[ph_map_ind[ind0]] = x;
      }
    }
    /*for(int hind1 = 0; hind1 < Chan.nh[chan3_1]; ++hind1){
      h1 = Chan.h_vec[Chan.h_index[chan3_1] + hind1].v1;
      for(int pind1 = 0; pind1 < Chan.np[chan3_2]; ++pind1){
	p1 = Chan.p_vec[Chan.p_index[chan3_2] + pind1].v1;
	ind0 = hp_3_index[ind] + hind1*Chan.np[chan3_2] + pind1;
	x = hp_3[ind0];
	x += hp_2[hp_map_ind[ind0]];
	hp_3[ind0] = x;
	hp_2[hp_map_ind[ind0]] = x;
      }
      }*/
  }
}

void OP_1b::Delete()
{
  if(chan3_length != 0){
    delete[] chan3_vec1;
    delete[] chan3_vec2;
  }
  if(hh_3_length != 0){
    delete[] hh_3_index;
    delete[] hh_3;
    delete[] hh_map_ind;
  }
  if(pp_3_length != 0){
    delete[] pp_3_index;
    delete[] pp_3;
    delete[] pp_map_ind;
  }
  if(ph_3_length != 0){
    delete[] ph_3_index;
    delete[] ph_3;
    delete[] ph_map_ind;
  }
  if(hp_3_length != 0){
    delete[] hp_3_index;
    delete[] hp_3;
    delete[] hp_map_ind;
  }
  if(chan2_length != 0){
    delete[] chan2_vec;
  }
  if(hh_2_length != 0){
    delete[] hh_2_index;
    delete[] hh_2;
  }
  if(pp_2_length != 0){
    delete[] pp_2_index;
    delete[] pp_2;
  }
  if(ph_2_length != 0){
    delete[] ph_2_index;
    delete[] ph_2;
  }
  if(hp_2_length != 0){
    delete[] hp_2_index;
    delete[] hp_2;
  }
}

void PA::Setup(Channels &Chan, int chan30)
{
  int chan1, chan2, chan3, chan31, chan32;
  int ind, jmin, length;
  int a, b, i;
  State ab, ib, ia;
  int x0 = Chan.qnums3[chan30].j;
  three_body abi, abi_j;
  three_body *thb_ind;
  three_body *thb_j;
  int *J;
  std::unordered_map<int,int> *R21_map, *R22_map, *L21_map, *L22_map;

  this->chan30 = chan30;
  this->nh0 = Chan.nh[chan30];
  this->np0 = Chan.np[chan30];
  this->npph0 = Chan.npph[chan30];

  this->R = new double[np0];
  this->L = new double[np0];
  this->R_length = np0;
  this->L_length = np0;
  for(ind = 0; ind < np0; ++ind){
  this->R[ind] = 0.0;
  this->L[ind] = 0.0;
  }
  this->R3 = new double[npph0];
  this->L3 = new double[npph0];
  this->R3_length = npph0;
  this->L3_length = npph0;
  for(ind = 0; ind < npph0; ++ind){
    this->R3[ind] = 0.0;
    this->L3[ind] = 0.0;
  }

  this->chan31_num = new int[Chan.size1 * Chan.size3];
  this->chan32_num = new int[Chan.size2 * Chan.size3];
  for(ind = 0; ind < Chan.size1 * Chan.size3; ++ind){ this->chan31_num[ind] = 0; }
  for(ind = 0; ind < Chan.size2 * Chan.size3; ++ind){ this->chan32_num[ind] = 0; }
  
  this->npph1 = 0;
  for(int pph = 0; pph < this->npph0; ++pph){
    a = Chan.pph_state(chan30, pph).v1;
    b = Chan.pph_state(chan30, pph).v2;
    i = Chan.pph_state(chan30, pph).v3;
    ab = Chan.pph_j[Chan.pph_index[chan30] + pph];
    if( (a < b) || (PAR.basis == "finite_J" && a == b) ){ ++this->npph1; }

    chan3 = Ind_Chan3(SPB.qnums[i]);
    chan1 = Ind_Chan1(ab);
    chan31 = Chan.size3 * chan1 + chan3;
    this->chan31_num[chan31] = Chan.npp[chan1] * Chan.nh[chan3];

    chan3 = Ind_Chan3(SPB.qnums[a]);
    minus(ib, SPB.qnums[i], SPB.qnums[b]);
    jmin = std::abs(SPB.qnums[i].j - SPB.qnums[b].j);
    while(ib.j >= jmin){
      chan2 = Ind_Chan2(ib);
      chan32 = Chan.size3 * chan2 + chan3;
      this->chan32_num[chan32] += Chan.np[chan3] * Chan.nhp1[chan2];
      ib.j -= 2;
    }
  }

  this->pph1_ind1 = new int[this->npph1];    // for pph1 -> pph0, map(p1,p2)
  this->pph1_ind2 = new int[this->npph1];    // for pph1 -> pph0, map(p2,p1)
  this->pph1_fac1 = new double[this->npph1]; // for pph1 -> pph0, 1.0        if p1==p2, -1*(-1)^(ja+jb-J) otherwise
  this->pph1_fac0 = new double[this->npph1]; // for pph1 -> pph0, sqrt(2)    if p1==p2, 1                 otherwise

  length = 0;
  this->R1_index = new int[Chan.size1 * Chan.size3];
  this->L1_index = new int[Chan.size1 * Chan.size3];
  for(int ind = 0; ind < Chan.size1 * Chan.size3; ++ind){
    if( this->chan31_num[ind] == 0 ){
      this->R1_index[ind] = -1;
      this->L1_index[ind] = -1;
      continue;
    }
    this->R1_index[ind] = length;
    this->L1_index[ind] = length;
    length += this->chan31_num[ind];
  }
  this->R1_length = length;
  this->L1_length = length;

  length = 0;
  this->R2_1_index = new int[Chan.size2 * Chan.size3];
  this->R2_2_index = new int[Chan.size2 * Chan.size3];
  this->L2_1_index = new int[Chan.size2 * Chan.size3];
  this->L2_2_index = new int[Chan.size2 * Chan.size3];
  for(chan32 = 0; chan32 < Chan.size2 * Chan.size3; ++chan32){
    if( this->chan32_num[chan32] == 0 ){
      this->R2_1_index[chan32] = -1;
      this->R2_2_index[chan32] = -1;
      this->L2_1_index[chan32] = -1;
      this->L2_2_index[chan32] = -1;
      continue;
    }
    this->R2_1_index[chan32] = length;
    this->R2_2_index[chan32] = length;
    this->L2_1_index[chan32] = length;
    this->L2_2_index[chan32] = length;
    length += this->chan32_num[chan32];
  }
  this->R2_1_length = length;
  this->R2_2_length = length;
  this->L2_1_length = length;
  this->L2_2_length = length;
  this->Rmap21_num = new int[length];
  this->Rmap22_num = new int[length];
  this->Rmap21_index = new int[length];
  this->Rmap22_index = new int[length];
  this->Lmap21_num = new int[length];
  this->Lmap22_num = new int[length];
  this->Lmap21_index = new int[length];
  this->Lmap22_index = new int[length];
  R21_map = new std::unordered_map<int,int>[length];
  R22_map = new std::unordered_map<int,int>[length];
  L21_map = new std::unordered_map<int,int>[length];
  L22_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap21_num[ind] = 0;
    this->Rmap22_num[ind] = 0;
    this->Lmap21_num[ind] = 0;
    this->Lmap22_num[ind] = 0;
  }

  thb_ind = new three_body[this->npph0];
  thb_j = new three_body[this->npph0];
  J = new int[this->npph0];
  this->npph1 = 0;
  for(int pph = 0; pph < this->npph0; ++pph){
    a = Chan.pph_state(chan30, pph).v1;
    b = Chan.pph_state(chan30, pph).v2;
    i = Chan.pph_state(chan30, pph).v3;
    abi.v1 = a;
    abi.v2 = b;
    abi.v3 = i;
    abi_j.v1 = SPB.qnums[a].j;
    abi_j.v2 = SPB.qnums[b].j;
    abi_j.v3 = SPB.qnums[i].j;
    thb_ind[pph] = abi;
    thb_j[pph] = abi_j;
    J[pph] = Chan.pph_j[Chan.pph_index[chan30] + pph].j;
    if( (a < b) || (PAR.basis == "finite_J" && a == b) ){
      this->pph1_ind1[this->npph1] = pph;
      this->pph1_ind2[this->npph1] = Chan.pph_map[chan30][Hash(b, a, i, J[pph])];
      if( a == b ){
	this->pph1_fac0[this->npph1] = std::sqrt(2.0);
	this->pph1_fac1[this->npph1] = 1.0;
      }
      else{
	this->pph1_fac0[this->npph1] = 1.0;
	this->pph1_fac1[this->npph1] = -1.0 * phase2(abi_j.v1 + abi_j.v2 - J[pph]);
      }
      ++this->npph1;
    }
    Map_Count_PAR_3_to_21(this->Rmap21_num, this->R2_1_index, Chan.p_map, Chan.hp1_map, Chan.nhp1, abi, J[pph], Chan.size3, R21_map);
    Map_Count_PAR_3_to_22(this->Rmap22_num, this->R2_2_index, Chan.p_map, Chan.hp1_map, Chan.nhp1, abi, J[pph], Chan.size3, R22_map);
    Map_Count_PAL_3_to_21(this->Lmap21_num, this->L2_1_index, Chan.hp1_map, Chan.p_map, Chan.np, abi, J[pph], Chan.size3, L21_map);
    Map_Count_PAL_3_to_22(this->Lmap22_num, this->L2_2_index, Chan.hp1_map, Chan.p_map, Chan.np, abi, J[pph], Chan.size3, L22_map);
  }

  this->Rmap1_ind = new int[this->R1_length];
  for(int ind = 0; ind < this->R1_length; ++ind){
    this->Rmap1_ind[ind] = 0;
  }

  length = 0;
  for(int ind = 0; ind < this->R2_1_length; ++ind){
    this->Rmap21_index[ind] = length;
    length += this->Rmap21_num[ind];
  }
  this->Rmap21_ind = new int[length];
  this->Rmap21_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap21_ind[ind] = 0;
    this->Rmap21_fac[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->R2_2_length; ++ind){
    this->Rmap22_index[ind] = length;
    length += this->Rmap22_num[ind];
  }
  this->Rmap22_ind = new int[length];
  this->Rmap22_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap22_ind[ind] = 0;
    this->Rmap22_fac[ind] = 0.0;
  }

  this->Lmap1_ind = new int[this->L1_length];
  for(int ind = 0; ind < this->L1_length; ++ind){
    this->Lmap1_ind[ind] = 0;
  }

  length = 0;
  for(int ind = 0; ind < this->L2_1_length; ++ind){
    this->Lmap21_index[ind] = length;
    length += this->Lmap21_num[ind];
  }
  this->Lmap21_ind = new int[length];
  this->Lmap21_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Lmap21_ind[ind] = 0;
    this->Lmap21_fac[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->L2_2_length; ++ind){
    this->Lmap22_index[ind] = length;
    length += this->Lmap22_num[ind];
  }
  this->Lmap22_ind = new int[length];
  this->Lmap22_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Lmap22_ind[ind] = 0;
    this->Lmap22_fac[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int pph = 0; pph < this->npph0; ++pph){
    Map_PAR_3_to_1(this->Rmap1_ind,this->R1_index,pph, Chan.pp_map,Chan.h_map,Chan.nh, thb_ind,thb_j,J,x0,Chan.size3);
    Map_PAR_3_to_21(this->Rmap21_ind,this->Rmap21_fac,this->Rmap21_index,this->R2_1_index,pph, Chan.p_map,Chan.hp1_map,Chan.nhp1, thb_ind,thb_j,J,x0,Chan.size3, R21_map);
    Map_PAR_3_to_22(this->Rmap22_ind,this->Rmap22_fac,this->Rmap22_index,this->R2_2_index,pph, Chan.p_map,Chan.hp1_map,Chan.nhp1, thb_ind,thb_j,J,x0,Chan.size3, R22_map);
    Map_PAL_3_to_1(this->Lmap1_ind,this->L1_index,pph, Chan.h_map,Chan.pp_map,Chan.npp, thb_ind,thb_j,J,x0,Chan.size3);
    Map_PAL_3_to_21(this->Lmap21_ind,this->Lmap21_fac,this->Lmap21_index,this->L2_1_index,pph, Chan.hp1_map,Chan.p_map,Chan.np, thb_ind,thb_j,J,x0,Chan.size3, L21_map);
    Map_PAL_3_to_22(this->Lmap22_ind,this->Lmap22_fac,this->Lmap22_index,this->L2_2_index,pph, Chan.hp1_map,Chan.p_map,Chan.np, thb_ind,thb_j,J,x0,Chan.size3, L22_map);
  }

  delete[] thb_ind;
  delete[] thb_j;
  delete[] J;
  delete[] R21_map;
  delete[] R22_map;
  delete[] L21_map;
  delete[] L22_map;
}

void PA::Set_R1(double *R1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Rmap1_ind[ind + offset];
    R1[ind] = this->R3[index];
  }
}

void PA::Set_R2_1(double *R21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    R21[ind] = 0.0;
    for(int n = 0; n < this->Rmap21_num[ind + offset]; ++n){
      index = this->Rmap21_index[ind + offset] + n;
      index1 = this->Rmap21_ind[index];
      R21[ind] += this->Rmap21_fac[index] * this->R3[index1];
    }
  }
}

void PA::Set_R2_2(double *R22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    R22[ind] = 0.0;
    for(int n = 0; n < this->Rmap22_num[ind + offset]; ++n){
      index = this->Rmap22_index[ind + offset] + n;
      index1 = this->Rmap22_ind[index];
      R22[ind] += this->Rmap22_fac[index] * this->R3[index1];
    }
  }
}

void PA::Gather_R1(double *R1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Rmap1_ind[ind + offset];
    this->R3[index] += R1[ind];
    R1[ind] = 0.0;
  }
}

void PA::Gather_R2_1(double *R21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Rmap21_num[ind + offset]; ++n){
      index = this->Rmap21_index[ind + offset] + n;
      index1 = this->Rmap21_ind[index];
      this->R3[index1] += this->Rmap21_fac[index] * R21[ind];
    }
    R21[ind] = 0.0;
  }
}

void PA::Gather_R2_2(double *R22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Rmap22_num[ind + offset]; ++n){
      index = this->Rmap22_index[ind + offset] + n;
      index1 = this->Rmap22_ind[index];
      this->R3[index1] += this->Rmap22_fac[index] * R22[ind];
    }
    R22[ind] = 0.0;
  }
}

void PA::Set_L1(double *L1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Lmap1_ind[ind + offset];
    L1[ind] = this->L3[index];
  }
}

void PA::Set_L2_1(double *L21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    L21[ind] = 0.0;
    for(int n = 0; n < this->Lmap21_num[ind + offset]; ++n){
      index = this->Lmap21_index[ind + offset] + n;
      index1 = this->Lmap21_ind[index];
      L21[ind] += this->Lmap21_fac[index] * this->L3[index1];
    }
  }
}

void PA::Set_L2_2(double *L22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    L22[ind] = 0.0;
    for(int n = 0; n < this->Lmap22_num[ind + offset]; ++n){
      index = this->Lmap22_index[ind + offset] + n;
      index1 = this->Lmap22_ind[index];
      L22[ind] += this->Lmap22_fac[index] * this->L3[index1];
      //std::cout << "Set22[" << ind << "] += " << this->Lmap22_fac[index] << " * " << this->L3[index1] << " = " << this->Lmap22_fac[index] * this->L3[index1] << std::endl;
    }
  }
}

void PA::Gather_L1(double *L1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Lmap1_ind[ind + offset];
    this->L3[index] += L1[ind];
    L1[ind] = 0.0;
  }
}

void PA::Gather_L2_1(double *L21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Lmap21_num[ind + offset]; ++n){
      index = this->Lmap21_index[ind + offset] + n;
      index1 = this->Lmap21_ind[index];
      this->L3[index1] += this->Lmap21_fac[index] * L21[ind];
      //std::cout << "Gather21[" << index1 << "] = " << this->Lmap21_fac[index] << " * " << L21[ind] << " = " << this->Lmap21_fac[index] * L21[ind] << std::endl;
    }
    L21[ind] = 0.0;
  }
}

void PA::Gather_L2_2(double *L22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Lmap22_num[ind + offset]; ++n){
      index = this->Lmap22_index[ind + offset] + n;
      index1 = this->Lmap22_ind[index];
      this->L3[index1] += this->Lmap22_fac[index] * L22[ind];
      //std::cout << "Gather22[" << index1 << "] = " << this->Lmap22_fac[index] << " * " << L22[ind] << " = " << this->Lmap22_fac[index] * L22[ind] << std::endl;
    }
    L22[ind] = 0.0;
  }
}

void PA::Setup(Channels &Chan, PA &PA1)
{
  int chan31, chan32, ind, length;
  this->chan30 = PA1.chan30;
  this->nh0 = PA1.nh0;
  this->np0 = PA1.np0;
  this->npph0 = PA1.npph0;
  this->npph1 = PA1.npph1;
  this->pph1_ind1 = new int[PA1.npph1];
  this->pph1_ind2 = new int[PA1.npph1];
  this->pph1_fac1 = new double[PA1.npph1];
  this->pph1_fac0 = new double[PA1.npph1];
  for(ind = 0; ind < PA1.npph1; ++ind){
    this->pph1_ind1[ind] = PA1.pph1_ind1[ind];
    this->pph1_ind2[ind] = PA1.pph1_ind2[ind];
    this->pph1_fac1[ind] = PA1.pph1_fac1[ind];
    this->pph1_fac0[ind] = PA1.pph1_fac0[ind];
  }

  this->R = new double[PA1.np0];
  this->L = new double[PA1.np0];
  this->R_length = PA1.np0;
  this->L_length = PA1.np0;
  for(ind = 0; ind < PA1.np0; ++ind){
    this->R[ind] = PA1.R[ind];
    this->L[ind] = PA1.L[ind];
  }
  this->R3 = new double[PA1.npph0];
  this->L3 = new double[PA1.npph0];
  this->R3_length = PA1.npph0;
  this->L3_length = PA1.npph0;
  for(ind = 0; ind < PA1.npph0; ++ind){
    this->R3[ind] = PA1.R3[ind];
    this->L3[ind] = PA1.L3[ind];
  }

  this->chan31_num = new int[Chan.size1 * Chan.size3];
  this->chan32_num = new int[Chan.size2 * Chan.size3];
  for(ind = 0; ind < Chan.size1 * Chan.size3; ++ind){ this->chan31_num[ind] = PA1.chan31_num[ind]; }
  for(ind = 0; ind < Chan.size2 * Chan.size3; ++ind){ this->chan32_num[ind] = PA1.chan32_num[ind]; }
  
  this->R1_index = new int[Chan.size1 * Chan.size3];
  this->L1_index = new int[Chan.size1 * Chan.size3];
  for(chan31 = 0; chan31 < Chan.size1 * Chan.size3; ++chan31){
    this->R1_index[chan31] = PA1.R1_index[chan31];
    this->L1_index[chan31] = PA1.L1_index[chan31];
  }
  this->R1_length = PA1.R1_length;
  this->L1_length = PA1.L1_length;

  this->R2_1_index = new int[Chan.size2 * Chan.size3];
  this->R2_2_index = new int[Chan.size2 * Chan.size3];
  this->L2_1_index = new int[Chan.size2 * Chan.size3];
  this->L2_2_index = new int[Chan.size2 * Chan.size3];
  for(chan32 = 0; chan32 < Chan.size2 * Chan.size3; ++chan32){
    this->R2_1_index[chan32] = PA1.R2_1_index[chan32];
    this->R2_2_index[chan32] = PA1.R2_2_index[chan32];
    this->L2_1_index[chan32] = PA1.L2_1_index[chan32];
    this->L2_2_index[chan32] = PA1.L2_2_index[chan32];
  }
  this->R2_1_length = PA1.R2_1_length;
  this->R2_2_length = PA1.R2_2_length;
  this->L2_1_length = PA1.L2_1_length;
  this->L2_2_length = PA1.L2_2_length;

  this->Rmap1_ind = new int[PA1.R1_length];
  for(ind = 0; ind < PA1.R1_length; ++ind){
    this->Rmap1_ind[ind] = PA1.Rmap1_ind[ind];
  }

  this->Rmap21_num = new int[PA1.R2_1_length];
  this->Rmap21_index = new int[PA1.R2_1_length];
  length = 0;
  for(ind = 0; ind < PA1.R2_1_length; ++ind){
    this->Rmap21_num[ind] = PA1.Rmap21_num[ind];
    this->Rmap21_index[ind] = PA1.Rmap21_index[ind];
    length += PA1.Rmap21_num[ind];
  }
  this->Rmap21_ind = new int[length];
  this->Rmap21_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Rmap21_ind[ind] = PA1.Rmap21_ind[ind];
    this->Rmap21_fac[ind] = PA1.Rmap21_fac[ind];
  }

  this->Rmap22_num = new int[PA1.R2_2_length];
  this->Rmap22_index = new int[PA1.R2_2_length];
  length = 0;
  for(ind = 0; ind < PA1.R2_2_length; ++ind){
    this->Rmap22_num[ind] = PA1.Rmap22_num[ind];
    this->Rmap22_index[ind] = PA1.Rmap22_index[ind];
    length += PA1.Rmap22_num[ind];
  }
  this->Rmap22_ind = new int[length];
  this->Rmap22_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Rmap22_ind[ind] = PA1.Rmap22_ind[ind];
    this->Rmap22_fac[ind] = PA1.Rmap22_fac[ind];
  }

  this->Lmap1_ind = new int[PA1.L1_length];
  for(ind = 0; ind < PA1.L1_length; ++ind){
    this->Lmap1_ind[ind] = PA1.Lmap1_ind[ind];
  }

  this->Lmap21_num = new int[PA1.L2_1_length];
  this->Lmap21_index = new int[PA1.L2_1_length];
  length = 0;
  for(ind = 0; ind < PA1.L2_1_length; ++ind){
    this->Lmap21_num[ind] = PA1.Lmap21_num[ind];
    this->Lmap21_index[ind] = PA1.Lmap21_index[ind];
    length += PA1.Lmap21_num[ind];
  }
  this->Lmap21_ind = new int[length];
  this->Lmap21_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Lmap21_ind[ind] = PA1.Lmap21_ind[ind];
    this->Lmap21_fac[ind] = PA1.Lmap21_fac[ind];
  }

  this->Lmap22_num = new int[PA1.L2_2_length];
  this->Lmap22_index = new int[PA1.L2_2_length];
  length = 0;
  for(ind = 0; ind < PA1.L2_2_length; ++ind){
    this->Lmap22_num[ind] = PA1.Lmap22_num[ind];
    this->Lmap22_index[ind] = PA1.Lmap22_index[ind];
    length += PA1.Lmap22_num[ind];
  }
  this->Lmap22_ind = new int[length];
  this->Lmap22_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Lmap22_ind[ind] = PA1.Lmap22_ind[ind];
    this->Lmap22_fac[ind] = PA1.Lmap22_fac[ind];
  }
}

void PA::Delete()
{
  delete[] R;
  delete[] R3;
  delete[] L;
  delete[] L3;

  delete[] pph1_ind1;
  delete[] pph1_ind2;
  delete[] pph1_fac1;
  delete[] pph1_fac0;

  delete[] chan31_num;
  delete[] chan32_num;

  delete[] R1_index;
  delete[] R2_1_index;
  delete[] R2_2_index;

  delete[] L1_index;
  delete[] L2_1_index;
  delete[] L2_2_index;

  delete[] Rmap1_ind;

  delete[] Rmap21_index;
  delete[] Rmap21_num;
  delete[] Rmap21_ind;
  delete[] Rmap21_fac;

  delete[] Rmap22_index;
  delete[] Rmap22_num;
  delete[] Rmap22_ind;
  delete[] Rmap22_fac;

  delete[] Lmap1_ind;

  delete[] Lmap21_index;
  delete[] Lmap21_num;
  delete[] Lmap21_ind;
  delete[] Lmap21_fac;

  delete[] Lmap22_index;
  delete[] Lmap22_num;
  delete[] Lmap22_ind;
  delete[] Lmap22_fac;
}

void PA::Zero_R()
{
  for(int ind = 0; ind < this->R_length; ++ind){ this->R[ind] = 0.0; }
  for(int ind = 0; ind < this->R3_length; ++ind){ this->R3[ind] = 0.0; }
}

void PA::Zero_L()
{
  for(int ind = 0; ind < this->L_length; ++ind){ this->L[ind] = 0.0; }
  for(int ind = 0; ind < this->L3_length; ++ind){ this->L3[ind] = 0.0; }
}

void PA::Update_R1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2)
{
  double p1 = 1.0, m12 = -0.5;
  int one = 1;
  int chan_ind1;
  char N = 'N';
  int chan32;
  int chan30 = this->chan30;
  int np0 = this->np0;
  int npph0 = this->npph0;
  int nhp0 = Chan.nhp1[Chan.ind0];
  double *X3, *R2;

  if( np0 != 0 ){
    chan_ind1 = Eff_Ints.Xpp.X_3_index[chan30];

    X3 = new double[np0 * np0];
    Eff_Ints.Xpp.Set_X_3(X3, np0 * np0, chan_ind1);
    //  w*r(a|){a}  =  Xpp3(a|c){a,c}.r(c|){c}
    dgemm_NN(X3, this->R, PA2.R, &np0, &one, &np0, &p1, &p1, &N, &N);

    delete[] X3;
  }

  chan32 = Chan.size3 * Chan.ind0 + chan30;
  if( this->chan32_num[chan32] != 0 ){
    chan_ind1 = this->R2_1_index[chan32];
    
    R2 = new double[np0 * nhp0];
    this->Set_R2_1(R2, np0 * nhp0, chan_ind1);
    //  w*r(a|){a}  =  r(ac|k){a,kc'}.Xhp2(k|c){kc'}
    dgemm_NN(R2, Eff_Ints.Xhp.X_2, PA2.R, &np0, &one, &nhp0, &p1, &p1, &N, &N);
    
    delete[] R2;
  }
    
  if( np0 * npph0 != 0 ){
    chan_ind1 = Eff_Ints.Xhppp.X_3_2_index[chan30];
    
    //  w*r(a|){a}  =  -(1/2).Xhppp3_2(ka|cd){a,cdk'}.r(cd|k){cdk'}
    dgemm_NN((Eff_Ints.Xhppp.X_3_2 + chan_ind1), this->R3, PA2.R, &np0, &one, &npph0, &m12, &p1, &N, &N);
  }
}

void PA::Update_R2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int nh, np, npp, nhp1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';
  int chan31, chan32;
  int chan30 = this->chan30;
  int nh0 = this->nh0;
  int np0 = this->np0;
  int npph0 = this->npph0;
  double *X3, *RR2, *R2, *RR1, *R1, *T3, *Q3;

  if( np0 * npph0 != 0 ){
    chan_ind1 = Eff_Ints.Xpphp.X_3_4_index[chan30];
    
    X3 = new double[npph0 * np0];
    Eff_Ints.Xpphp.Set_X_3_4(X3, npph0 * np0, chan_ind1);

    //  w*r(ab|i){abi'}  =  -Xpphp3_4(ab|ic){abi',c}.r(c|){c}
    dgemm_NN(X3, this->R, PA2.R3, &npph0, &one, &np0, &m1, &p1, &N, &N);
    
    delete[] X3;
  }

  /*if( nh0 * npph0 != 0 ){
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan30];
    chan_ind2 = Amps.D1.T3_3_index[chan30];
    
    T3 = new double[npph0 * nh0];
    Amps.D1.Set_T3_3(T3, npph0 * nh0, chan_ind2);
    Q3 = new double[nh0];
    for(int ind = 0; ind < nh0; ++ind){ Q3[ind] = 0.0; }
    //w*r(ab|i){abi'}  =  -(1/2).T3_3(ab|ki){abi',k}.Vhhpp3_1(kl|cd){k,cdl'}.r(cd|l){cdl'}
    dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), this->R3, Q3, &nh0, &one, &npph0, &p1, &p1, &N, &N);
    dgemm_NN(T3, Q3, PA2.R3, &npph0, &one, &nh0, &m12, &p1, &N, &N);
    
    delete[] T3;
    delete[] Q3;
    }*/
    
  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    if( nhp1 == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      np = Chan.np[chan3];
      if( np == 0 ){ continue; }
      chan32 = Chan.size3 * chan2 + chan3;
      
      if( this->chan32_num[chan32] != 0 ){
	chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
	chan_ind2 = this->R2_1_index[chan32];
	chan_ind3 = this->R2_2_index[chan32];
	chan_ind4 = Eff_Ints.Xhphp.X_2_1_index[chan2];
	
	X3 = new double[np * np];
	Eff_Ints.Xpp.Set_X_3(X3, np*np, chan_ind1);
	RR2 = new double[np * nhp1];
        this->Set_R2_2(RR2, np*nhp1, chan_ind3);
	R2 = new double[np * nhp1];
	for(int ind = 0; ind < np * nhp1; ++ind){ R2[ind] = 0.0; }
	//  w*r(ab|i){a,ib'}  =  -Xpp3(a|c){a,c}.r(bc|i){c,ib'}
	dgemm_NN(X3, RR2, R2, &np, &nhp1, &np, &m1, &p1, &N, &N);
	//  w*r(ab|i){a,ib'}  =  r(ca|k){a,kc'}.Xhphp2_1(kb|ic){kc',ib'}
	dgemm_NN(RR2, (Eff_Ints.Xhphp.X_2_1 + chan_ind4), R2, &np, &nhp1, &nhp1, &p1, &p1, &N, &N);
	PA2.Gather_R2_1(R2, np * nhp1, chan_ind2);

	//  w*r(ab|i){b,ia'}  =  Xpp3(b|c){b,c}.r(ac|i){c,ia'}
	dgemm_NN(X3, RR2, R2, &np, &nhp1, &np, &p1, &p1, &N, &N);
	//  w*r(ab|i){b,ia'}  =  -r(cb|k){b,kc'}.Xhphp2_1(ka|ic){kc',ia'}
	dgemm_NN(RR2, (Eff_Ints.Xhphp.X_2_1 + chan_ind4), R2, &np, &nhp1, &nhp1, &m1, &p1, &N, &N);
	PA2.Gather_R2_2(R2, np * nhp1, chan_ind3);

	delete[] X3;
	delete[] R2;
	delete[] RR2;
      }
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    if( npp == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      nh = Chan.nh[chan3];
      if( nh == 0 ){ continue; }
      chan31 = Chan.size3 * chan1 + chan3;

      if( this->chan31_num[chan31] != 0 ){
	chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
	chan_ind2 = this->R1_index[chan31];
	chan_ind3 = Eff_Ints.Xpppp.X_1_index[chan1];
	
	X3 = new double[nh * nh];
	Eff_Ints.Xhh.Set_X_3(X3, nh*nh, chan_ind1);
	RR1 = new double[npp * nh];
	this->Set_R1(RR1, npp*nh, chan_ind2);
	R1 = new double[npp * nh];
	for(int ind = 0; ind < npp * nh; ++ind){ R1[ind] = 0.0; }
	//  w*r(ab|i){ab,i}  =  -r(ab|k){ab,k}.Xhh3(k|i){k,i}
	dgemm_NN(RR1, X3, R1, &npp, &nh, &nh, &m1, &p1, &N, &N);
	//  w*r(ab|i){ab,i}  =  (1/2).Xpppp1(ab|cd){ab,cd}.r(cd|i){cd,i}
	dgemm_NN((Eff_Ints.Xpppp.X_1 + chan_ind3), RR1, R1, &npp, &nh, &npp, &p12, &p1, &N, &N);
	PA2.Gather_R1(R1, npp*nh, chan_ind2);

	delete[] X3;
	delete[] R1;
	delete[] RR1;
      }
    }
  }
}

void PA::Update_L1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2)
{
  double p1 = 1.0, m12 = -0.5;
  int one = 1;
  int chan_ind1;
  char N = 'N';
  int chan30 = this->chan30;
  int np0 = this->np0;
  int npph0 = this->npph0;
  double *X3;

  //  E*l(|a){a}  =  E0.l(|a){a}  
  double E0 = Amps.dE;
  for(int p = 0; p < np0; ++p){
    PA2.L[p] += E0 * this->L[p];
  }

  if( np0 != 0 ){
    chan_ind1 = Eff_Ints.Xpp.X_3_index[chan30];

    X3 = new double[np0 * np0];
    Eff_Ints.Xpp.Set_X_3(X3, np0 * np0, chan_ind1);
    //  E*l(|a){a}  =  l(|c){c}.Xpp3(c|a){c,a}
    dgemm_NN(this->L, X3, PA2.L, &one, &np0, &np0, &p1, &p1, &N, &N);

    delete[] X3;
  }

  if( np0 * npph0 != 0 ){
    chan_ind1 = Eff_Ints.Xpphp.X_3_4_index[chan30];
    
    X3 = new double[npph0 * np0];
    Eff_Ints.Xpphp.Set_X_3_4(X3, npph0 * np0, chan_ind1);
    //  E*l(|a){a}  =  -(1/2).l(k|cd){cdk'}.Xpphp3_4(cd|ka){cdk',a}
    dgemm_NN(this->L3, X3, PA2.L, &one, &np0, &npph0, &m12, &p1, &N, &N);
    
    delete[] X3;
  }
}

void PA::Update_L2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int nh, np, npp, nhp1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';
  int chan31, chan32;
  int chan30 = this->chan30;
  int nh0 = this->nh0;
  int np0 = this->np0;
  int npph0 = this->npph0;
  int nhp0 = Chan.nhp1[Chan.ind0];
  double *X3, *L2_1, *L2_2, *LL2, *L2, *LL1, *L1, *T3, *Q3;

  //  E*l(i|ab){abi'}  =  E0.l(i|ab){abi'}
  double E0 = Amps.dE;
  for(int pph = 0; pph < npph0; ++pph){
    PA2.L3[pph] += E0 * this->L3[pph];
  }

  chan32 = Chan.size3 * Chan.ind0 + chan30;
  if( this->chan32_num[chan32] != 0 ){
    chan_ind1 = this->L2_1_index[chan32];
    chan_ind2 = this->L2_2_index[chan32];
    
    L2_1 = new double[nhp0 * np0];
    L2_2 = new double[nhp0 * np0];
    for(int ind = 0; ind < nhp0 * np0; ++ind){
      L2_1[ind] = 0.0;
      L2_2[ind] = 0.0;
    }
    for(int hp1 = 0; hp1 < nhp0; ++hp1){
      for(int p = 0; p < np0; ++p){
	//  E*l(i|ab){ib',a}  =  Xhp2(i|b){ib'}.l(|a){a}
	L2_1[hp1*np0 + p] += Eff_Ints.Xhp.X_2[hp1] * this->L[p];
	//  E*l(i|ab){ia',b}  =  -Xhp2(i|a){ia'}.l(|b){b}
	L2_2[hp1*np0 + p] -= Eff_Ints.Xhp.X_2[hp1] * this->L[p];
      }
    }
    PA2.Gather_L2_1(L2_1, nhp0*np0, chan_ind1);
    PA2.Gather_L2_2(L2_2, nhp0*np0, chan_ind2);
    delete[] L2_1;
    delete[] L2_2;
  }

  if( np0 * npph0 != 0 ){
    chan_ind1 = Eff_Ints.Xhppp.X_3_2_index[chan30];
    
    //  E*l(i|ab){abi'}  =  -l(|c){,c}.Xhppp3_2(ic|ab){c,abi'}
    dgemm_NN(this->L, (Eff_Ints.Xhppp.X_3_2 + chan_ind1), PA2.L3, &one, &npph0, &np0, &m1, &p1, &N, &N);
  }

  /*if( nh0 * npph0 != 0 ){
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan30];
    chan_ind2 = Amps.D1.T3_3_index[chan30];
    
    T3 = new double[npph0 * nh0];
    Amps.D1.Set_T3_3(T3, npph0 * nh0, chan_ind2);
    Q3 = new double[nh0];
    for(int ind = 0; ind < nh0; ++ind){ Q3[ind] = 0.0; }
    //E*l(i|ab){abi'}  =  -(1/2).l(l|cd){cdl'}.T3_3(cd|kl){cdl',k}.Vhhpp3_1(ki|ab){k,abi'}
    dgemm_NN(this->L3, T3, Q3, &one, &nh0, &npph0, &p1, &p1, &N, &N);
    dgemm_NN(Q3, (Ints.Vhhpp.V_3_1 + chan_ind1), PA2.L3, &one, &npph0, &nh0, &m12, &p1, &N, &N);
    
    delete[] T3;
    delete[] Q3;
    }*/

  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nhp1 = Chan.nhp1[chan2];
    if( nhp1 == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      np = Chan.np[chan3];
      if( np == 0 ){ continue; }
      chan32 = Chan.size3 * chan2 + chan3;
      
      if( this->chan32_num[chan32] != 0 ){
	chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
	chan_ind2 = this->L2_1_index[chan32];
	chan_ind3 = this->L2_2_index[chan32];
	chan_ind4 = Eff_Ints.Xhphp.X_2_1_index[chan2];
	
	X3 = new double[np * np];
	Eff_Ints.Xpp.Set_X_3(X3, np*np, chan_ind1);
	LL2 = new double[nhp1 * np];
	this->Set_L2_2(LL2, nhp1*np, chan_ind3);
	L2 = new double[nhp1 * np];
	for(int ind = 0; ind < np * nhp1; ++ind){ L2[ind] = 0.0; }
	//  E*l(i|ab){ib',a}  =  -l(i|ac){ia',c}.Xpp3(c|b){c,b}
	dgemm_NN(LL2, X3, L2, &nhp1, &np, &np, &m1, &p1, &N, &N);
	//  E*l(i|ab){ib',a}  =  Xhphp2_1(ic|kb){ib',kc'}.l(k|ca){kc',a}
	dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind4), LL2, L2, &nhp1, &np, &nhp1, &p1, &p1, &N, &N);
	PA2.Gather_L2_1(L2, nhp1 * np, chan_ind2);

	//  E*l(i|ab){ia',b}  =  l(i|bc){ib',c}.Xpp3(c|a){c,a}
	dgemm_NN(LL2, X3, L2, &nhp1, &np, &np, &p1, &p1, &N, &N);
	//  E*l(i|ab){ia',b}  =  -Xhphp2_1(ic|ka){ia',kc'}.l(k|cb){kc',b}
	dgemm_NN((Eff_Ints.Xhphp.X_2_1 + chan_ind4), LL2, L2, &nhp1, &np, &nhp1, &m1, &p1, &N, &N);
	PA2.Gather_L2_2(L2, nhp1 * np, chan_ind3);

	delete[] X3;
	delete[] L2;
	delete[] LL2;
      }
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    npp = Chan.npp[chan1];
    if( npp == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      nh = Chan.nh[chan3];
      if( nh == 0 ){ continue; }
      chan31 = Chan.size3 * chan1 + chan3;

      if( this->chan31_num[chan31] != 0 ){
	chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
	chan_ind2 = this->L1_index[chan31];
	chan_ind3 = Eff_Ints.Xpppp.X_1_index[chan1];
	
	X3 = new double[nh * nh];
	Eff_Ints.Xhh.Set_X_3(X3, nh*nh, chan_ind1);
	LL1 = new double[nh * npp];
	this->Set_L1(LL1, nh*npp, chan_ind2);
	L1 = new double[nh * npp];
	for(int ind = 0; ind < nh * npp; ++ind){ L1[ind] = 0.0; }
	//  E*l(i|ab){i,ab}  =  -Xhh3(i|k){i,k}.l(k|ab){k,ab}
	dgemm_NN(X3, LL1, L1, &nh, &npp, &nh, &m1, &p1, &N, &N);
	//  E*l(i|ab){i,ab}  =  (1/2).l(i|cd){i,cd}.Xpppp1(cd|ab){cd,ab}
	dgemm_NN(LL1, (Eff_Ints.Xpppp.X_1 + chan_ind3), L1, &nh, &npp, &npp, &p12, &p1, &N, &N);
	PA2.Gather_L1(L1, nh*npp, chan_ind2);

	delete[] X3;
	delete[] L1;
	delete[] LL1;
      }
    }
  }
}

void PR::Setup(Channels &Chan, int chan30)
{
  int chan1, chan2, chan3, chan31, chan32;
  int ind, jmin, length;
  int i, j, a;
  State ij, aj, ai;
  int x0 = Chan.qnums3[chan30].j;
  three_body ija, ija_j;
  three_body *thb_ind;
  three_body *thb_j;
  int *J;
  std::unordered_map<int,int> *R21_map, *R22_map, *L21_map, *L22_map;

  this->chan30 = chan30;
  this->np0 = Chan.np[chan30];
  this->nh0 = Chan.nh[chan30];
  this->nhhp0 = Chan.nhhp[chan30];

  this->R = new double[nh0];
  this->L = new double[nh0];
  this->R_length = nh0;
  this->L_length = nh0;
  for(ind = 0; ind < nh0; ++ind){
  this->R[ind] = 0.0;
  this->L[ind] = 0.0;
  }
  this->R3 = new double[nhhp0];
  this->L3 = new double[nhhp0];
  this->R3_length = nhhp0;
  this->L3_length = nhhp0;
  for(ind = 0; ind < nhhp0; ++ind){
    this->R3[ind] = 0.0;
    this->L3[ind] = 0.0;
  }

  this->chan31_num = new int[Chan.size1 * Chan.size3];
  this->chan32_num = new int[Chan.size2 * Chan.size3];
  for(ind = 0; ind < Chan.size1 * Chan.size3; ++ind){ this->chan31_num[ind] = 0; }
  for(ind = 0; ind < Chan.size2 * Chan.size3; ++ind){ this->chan32_num[ind] = 0; }
  
  this->nhhp1 = 0;
  for(int hhp = 0; hhp < this->nhhp0; ++hhp){
    i = Chan.hhp_state(chan30, hhp).v1;
    j = Chan.hhp_state(chan30, hhp).v2;
    a = Chan.hhp_state(chan30, hhp).v3;
    ij = Chan.hhp_j[Chan.hhp_index[chan30] + hhp];
    if( (i < j) || (PAR.basis == "finite_J" && i == j) ){ ++this->nhhp1; }

    chan3 = Ind_Chan3(SPB.qnums[a]);
    chan1 = Ind_Chan1(ij);
    chan31 = Chan.size3 * chan1 + chan3;
    this->chan31_num[chan31] = Chan.nhh[chan1] * Chan.np[chan3];

    chan3 = Ind_Chan3(SPB.qnums[i]);
    minus(aj, SPB.qnums[a], SPB.qnums[j]);
    jmin = std::abs(SPB.qnums[a].j - SPB.qnums[j].j);
    while(aj.j >= jmin){
      chan2 = Ind_Chan2(aj);
      chan32 = Chan.size3 * chan2 + chan3;
      this->chan32_num[chan32] += Chan.nph1[chan2] * Chan.nh[chan3];
      aj.j -= 2;
    }
  }

  this->hhp1_ind1 = new int[this->nhhp1];    // for hhp1 -> hhp0, map(h1,h2)
  this->hhp1_ind2 = new int[this->nhhp1];    // for hhp1 -> hhp0, map(h2,h1)
  this->hhp1_fac1 = new double[this->nhhp1]; // for hhp1 -> hhp0, 1.0        if h1==h2, -1*(-1)^(ji+jj-J) otherwise
  this->hhp1_fac0 = new double[this->nhhp1]; // for hhp1 -> hhp0, sqrt(2)    if h1==h2, 1                 otherwise

  length = 0;
  this->R1_index = new int[Chan.size1 * Chan.size3];
  this->L1_index = new int[Chan.size1 * Chan.size3];
  for(int ind = 0; ind < Chan.size1 * Chan.size3; ++ind){
    if( this->chan31_num[ind] == 0 ){
      this->R1_index[ind] = -1;
      this->L1_index[ind] = -1;
      continue;
    }
    this->R1_index[ind] = length;
    this->L1_index[ind] = length;
    length += this->chan31_num[ind];
  }
  this->R1_length = length;
  this->L1_length = length;

  length = 0;
  this->R2_1_index = new int[Chan.size2 * Chan.size3];
  this->R2_2_index = new int[Chan.size2 * Chan.size3];
  this->L2_1_index = new int[Chan.size2 * Chan.size3];
  this->L2_2_index = new int[Chan.size2 * Chan.size3];
  for(chan32 = 0; chan32 < Chan.size2 * Chan.size3; ++chan32){
    if( this->chan32_num[chan32] == 0 ){
      this->R2_1_index[chan32] = -1;
      this->R2_2_index[chan32] = -1;
      this->L2_1_index[chan32] = -1;
      this->L2_2_index[chan32] = -1;
      continue;
    }
    this->R2_1_index[chan32] = length;
    this->R2_2_index[chan32] = length;
    this->L2_1_index[chan32] = length;
    this->L2_2_index[chan32] = length;
    length += this->chan32_num[chan32];
  }
  this->R2_1_length = length;
  this->R2_2_length = length;
  this->L2_1_length = length;
  this->L2_2_length = length;
  this->Rmap21_num = new int[length];
  this->Rmap22_num = new int[length];
  this->Rmap21_index = new int[length];
  this->Rmap22_index = new int[length];
  this->Lmap21_num = new int[length];
  this->Lmap22_num = new int[length];
  this->Lmap21_index = new int[length];
  this->Lmap22_index = new int[length];
  R21_map = new std::unordered_map<int,int>[length];
  R22_map = new std::unordered_map<int,int>[length];
  L21_map = new std::unordered_map<int,int>[length];
  L22_map = new std::unordered_map<int,int>[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap21_num[ind] = 0;
    this->Rmap22_num[ind] = 0;
    this->Lmap21_num[ind] = 0;
    this->Lmap22_num[ind] = 0;
  }

  thb_ind = new three_body[this->nhhp0];
  thb_j = new three_body[this->nhhp0];
  J = new int[this->nhhp0];
  this->nhhp1 = 0;
  for(int hhp = 0; hhp < this->nhhp0; ++hhp){
    i = Chan.hhp_state(chan30, hhp).v1;
    j = Chan.hhp_state(chan30, hhp).v2;
    a = Chan.hhp_state(chan30, hhp).v3;
    ija.v1 = i;
    ija.v2 = j;
    ija.v3 = a;
    ija_j.v1 = SPB.qnums[i].j;
    ija_j.v2 = SPB.qnums[j].j;
    ija_j.v3 = SPB.qnums[a].j;
    thb_ind[hhp] = ija;
    thb_j[hhp] = ija_j;
    J[hhp] = Chan.hhp_j[Chan.hhp_index[chan30] + hhp].j;
    if( (i < j) || (PAR.basis == "finite_J" && i == j) ){
      this->hhp1_ind1[this->nhhp1] = hhp;
      this->hhp1_ind2[this->nhhp1] = Chan.hhp_map[chan30][Hash(j, i, a, J[hhp])];
      if( i == j ){
	this->hhp1_fac0[this->nhhp1] = std::sqrt(2.0);
	this->hhp1_fac1[this->nhhp1] = 1.0;
      }
      else{
	this->hhp1_fac0[this->nhhp1] = 1.0;
	this->hhp1_fac1[this->nhhp1] = -1.0 * phase2(ija_j.v1 + ija_j.v2 - J[hhp]);
      }
      ++this->nhhp1;
    }
    Map_Count_PRR_3_to_21(this->Rmap21_num, this->R2_1_index, Chan.ph1_map, Chan.h_map, Chan.nh, ija, J[hhp], Chan.size3, R21_map);
    Map_Count_PRR_3_to_22(this->Rmap22_num, this->R2_2_index, Chan.ph1_map, Chan.h_map, Chan.nh, ija, J[hhp], Chan.size3, R22_map);
    Map_Count_PRL_3_to_21(this->Lmap21_num, this->L2_1_index, Chan.h_map, Chan.ph1_map, Chan.nph1, ija, J[hhp], Chan.size3, L21_map);
    Map_Count_PRL_3_to_22(this->Lmap22_num, this->L2_2_index, Chan.h_map, Chan.ph1_map, Chan.nph1, ija, J[hhp], Chan.size3, L22_map);
  }

  this->Rmap1_ind = new int[this->R1_length];
  for(int ind = 0; ind < this->R1_length; ++ind){
    this->Rmap1_ind[ind] = 0;
  }

  length = 0;
  for(int ind = 0; ind < this->R2_1_length; ++ind){
    this->Rmap21_index[ind] = length;
    length += this->Rmap21_num[ind];
  }
  this->Rmap21_ind = new int[length];
  this->Rmap21_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap21_ind[ind] = 0;
    this->Rmap21_fac[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->R2_2_length; ++ind){
    this->Rmap22_index[ind] = length;
    length += this->Rmap22_num[ind];
  }
  this->Rmap22_ind = new int[length];
  this->Rmap22_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Rmap22_ind[ind] = 0;
    this->Rmap22_fac[ind] = 0.0;
  }

  this->Lmap1_ind = new int[this->L1_length];
  for(int ind = 0; ind < this->L1_length; ++ind){
    this->Lmap1_ind[ind] = 0;
  }

  length = 0;
  for(int ind = 0; ind < this->L2_1_length; ++ind){
    this->Lmap21_index[ind] = length;
    length += this->Lmap21_num[ind];
  }
  this->Lmap21_ind = new int[length];
  this->Lmap21_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Lmap21_ind[ind] = 0;
    this->Lmap21_fac[ind] = 0.0;
  }

  length = 0;
  for(int ind = 0; ind < this->L2_2_length; ++ind){
    this->Lmap22_index[ind] = length;
    length += this->Lmap22_num[ind];
  }
  this->Lmap22_ind = new int[length];
  this->Lmap22_fac = new double[length];
  for(int ind = 0; ind < length; ++ind){
    this->Lmap22_ind[ind] = 0;
    this->Lmap22_fac[ind] = 0.0;
  }

  #pragma omp parallel for schedule(static)
  for(int hhp = 0; hhp < this->nhhp0; ++hhp){
    Map_PRR_3_to_1(this->Rmap1_ind,this->R1_index,hhp, Chan.p_map,Chan.hh_map,Chan.nhh, thb_ind,thb_j,J,x0,Chan.size3);
    Map_PRR_3_to_21(this->Rmap21_ind,this->Rmap21_fac,this->Rmap21_index,this->R2_1_index,hhp, Chan.ph1_map,Chan.h_map,Chan.nh, thb_ind,thb_j,J,x0,Chan.size3, R21_map);
    Map_PRR_3_to_22(this->Rmap22_ind,this->Rmap22_fac,this->Rmap22_index,this->R2_2_index,hhp, Chan.ph1_map,Chan.h_map,Chan.nh, thb_ind,thb_j,J,x0,Chan.size3, R22_map);
    Map_PRL_3_to_1(this->Lmap1_ind,this->L1_index,hhp, Chan.hh_map,Chan.p_map,Chan.np, thb_ind,thb_j,J,x0,Chan.size3);
    Map_PRL_3_to_21(this->Lmap21_ind,this->Lmap21_fac,this->Lmap21_index,this->L2_1_index,hhp, Chan.h_map,Chan.ph1_map,Chan.nph1, thb_ind,thb_j,J,x0,Chan.size3, L21_map);
    Map_PRL_3_to_22(this->Lmap22_ind,this->Lmap22_fac,this->Lmap22_index,this->L2_2_index,hhp, Chan.h_map,Chan.ph1_map,Chan.nph1, thb_ind,thb_j,J,x0,Chan.size3, L22_map);
  }

  delete[] thb_ind;
  delete[] thb_j;
  delete[] J;
  delete[] R21_map;
  delete[] R22_map;
  delete[] L21_map;
  delete[] L22_map;
}

void PR::Set_R1(double *R1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Rmap1_ind[ind + offset];
    R1[ind] = this->R3[index];
  }
}

void PR::Set_R2_1(double *R21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    R21[ind] = 0.0;
    for(int n = 0; n < this->Rmap21_num[ind + offset]; ++n){
      index = this->Rmap21_index[ind + offset] + n;
      index1 = this->Rmap21_ind[index];
      R21[ind] += this->Rmap21_fac[index] * this->R3[index1];
    }
  }
}

void PR::Set_R2_2(double *R22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    R22[ind] = 0.0;
    for(int n = 0; n < this->Rmap22_num[ind + offset]; ++n){
      index = this->Rmap22_index[ind + offset] + n;
      index1 = this->Rmap22_ind[index];
      R22[ind] += this->Rmap22_fac[index] * this->R3[index1];
    }
  }
}

void PR::Gather_R1(double *R1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Rmap1_ind[ind + offset];
    this->R3[index] += R1[ind];
    R1[ind] = 0.0;
  }
}

void PR::Gather_R2_1(double *R21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Rmap21_num[ind + offset]; ++n){
      index = this->Rmap21_index[ind + offset] + n;
      index1 = this->Rmap21_ind[index];
      this->R3[index1] += this->Rmap21_fac[index] * R21[ind];
    }
    R21[ind] = 0.0;
  }
}

void PR::Gather_R2_2(double *R22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Rmap22_num[ind + offset]; ++n){
      index = this->Rmap22_index[ind + offset] + n;
      index1 = this->Rmap22_ind[index];
      this->R3[index1] += this->Rmap22_fac[index] * R22[ind];
    }
    R22[ind] = 0.0;
  }
}

void PR::Set_L1(double *L1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Lmap1_ind[ind + offset];
    L1[ind] = this->L3[index];
  }
}

void PR::Set_L2_1(double *L21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    L21[ind] = 0.0;
    for(int n = 0; n < this->Lmap21_num[ind + offset]; ++n){
      index = this->Lmap21_index[ind + offset] + n;
      index1 = this->Lmap21_ind[index];
      L21[ind] += this->Lmap21_fac[index] * this->L3[index1];
    }
  }
}

void PR::Set_L2_2(double *L22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    L22[ind] = 0.0;
    for(int n = 0; n < this->Lmap22_num[ind + offset]; ++n){
      index = this->Lmap22_index[ind + offset] + n;
      index1 = this->Lmap22_ind[index];
      L22[ind] += this->Lmap22_fac[index] * this->L3[index1];
    }
  }
}

void PR::Gather_L1(double *L1, int length, int offset)
{
  int index;
  for(int ind = 0; ind < length; ++ind){
    index = this->Lmap1_ind[ind + offset];
    this->L3[index] += L1[ind];
    L1[ind] = 0.0;
  }
}

void PR::Gather_L2_1(double *L21, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Lmap21_num[ind + offset]; ++n){
      index = this->Lmap21_index[ind + offset] + n;
      index1 = this->Lmap21_ind[index];
      this->L3[index1] += this->Lmap21_fac[index] * L21[ind];
    }
    L21[ind] = 0.0;
  }
}

void PR::Gather_L2_2(double *L22, int length, int offset)
{
  int index, index1;
  for(int ind = 0; ind < length; ++ind){
    for(int n = 0; n < this->Lmap22_num[ind + offset]; ++n){
      index = this->Lmap22_index[ind + offset] + n;
      index1 = this->Lmap22_ind[index];
      this->L3[index1] += this->Lmap22_fac[index] * L22[ind];
    }
    L22[ind] = 0.0;
  }
}

void PR::Setup(Channels &Chan, PR &PR1)
{
  int chan31, chan32, ind, length;
  this->chan30 = PR1.chan30;
  this->np0 = PR1.np0;
  this->nh0 = PR1.nh0;
  this->nhhp0 = PR1.nhhp0;
  this->nhhp1 = PR1.nhhp1;
  this->hhp1_ind1 = new int[PR1.nhhp1];
  this->hhp1_ind2 = new int[PR1.nhhp1];
  this->hhp1_fac1 = new double[PR1.nhhp1];
  this->hhp1_fac0 = new double[PR1.nhhp1];
  for(ind = 0; ind < PR1.nhhp1; ++ind){
    this->hhp1_ind1[ind] = PR1.hhp1_ind1[ind];
    this->hhp1_ind2[ind] = PR1.hhp1_ind2[ind];
    this->hhp1_fac1[ind] = PR1.hhp1_fac1[ind];
    this->hhp1_fac0[ind] = PR1.hhp1_fac0[ind];
  }

  this->R = new double[PR1.nh0];
  this->L = new double[PR1.nh0];
  this->R_length = PR1.nh0;
  this->L_length = PR1.nh0;
  for(ind = 0; ind < PR1.nh0; ++ind){
    this->R[ind] = PR1.R[ind];
    this->L[ind] = PR1.L[ind];
  }
  this->R3 = new double[PR1.nhhp0];
  this->L3 = new double[PR1.nhhp0];
  this->R3_length = PR1.nhhp0;
  this->L3_length = PR1.nhhp0;
  for(ind = 0; ind < PR1.nhhp0; ++ind){
    this->R3[ind] = PR1.R3[ind];
    this->L3[ind] = PR1.L3[ind];
  }

  this->chan31_num = new int[Chan.size1 * Chan.size3];
  this->chan32_num = new int[Chan.size2 * Chan.size3];
  for(ind = 0; ind < Chan.size1 * Chan.size3; ++ind){ this->chan31_num[ind] = PR1.chan31_num[ind]; }
  for(ind = 0; ind < Chan.size2 * Chan.size3; ++ind){ this->chan32_num[ind] = PR1.chan32_num[ind]; }
  
  this->R1_index = new int[Chan.size1 * Chan.size3];
  this->L1_index = new int[Chan.size1 * Chan.size3];
  for(chan31 = 0; chan31 < Chan.size1 * Chan.size3; ++chan31){
    this->R1_index[chan31] = PR1.R1_index[chan31];
    this->L1_index[chan31] = PR1.L1_index[chan31];
  }
  this->R1_length = PR1.R1_length;
  this->L1_length = PR1.L1_length;

  this->R2_1_index = new int[Chan.size2 * Chan.size3];
  this->R2_2_index = new int[Chan.size2 * Chan.size3];
  this->L2_1_index = new int[Chan.size2 * Chan.size3];
  this->L2_2_index = new int[Chan.size2 * Chan.size3];
  for(chan32 = 0; chan32 < Chan.size2 * Chan.size3; ++chan32){
    this->R2_1_index[chan32] = PR1.R2_1_index[chan32];
    this->R2_2_index[chan32] = PR1.R2_2_index[chan32];
    this->L2_1_index[chan32] = PR1.L2_1_index[chan32];
    this->L2_2_index[chan32] = PR1.L2_2_index[chan32];
  }
  this->R2_1_length = PR1.R2_1_length;
  this->R2_2_length = PR1.R2_2_length;
  this->L2_1_length = PR1.L2_1_length;
  this->L2_2_length = PR1.L2_2_length;

  this->Rmap1_ind = new int[PR1.R1_length];
  for(ind = 0; ind < PR1.R1_length; ++ind){
    this->Rmap1_ind[ind] = PR1.Rmap1_ind[ind];
  }

  this->Rmap21_num = new int[PR1.R2_1_length];
  this->Rmap21_index = new int[PR1.R2_1_length];
  length = 0;
  for(ind = 0; ind < PR1.R2_1_length; ++ind){
    this->Rmap21_num[ind] = PR1.Rmap21_num[ind];
    this->Rmap21_index[ind] = PR1.Rmap21_index[ind];
    length += PR1.Rmap21_num[ind];
  }
  this->Rmap21_ind = new int[length];
  this->Rmap21_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Rmap21_ind[ind] = PR1.Rmap21_ind[ind];
    this->Rmap21_fac[ind] = PR1.Rmap21_fac[ind];
  }

  this->Rmap22_num = new int[PR1.R2_2_length];
  this->Rmap22_index = new int[PR1.R2_2_length];
  length = 0;
  for(ind = 0; ind < PR1.R2_2_length; ++ind){
    this->Rmap22_num[ind] = PR1.Rmap22_num[ind];
    this->Rmap22_index[ind] = PR1.Rmap22_index[ind];
    length += PR1.Rmap22_num[ind];
  }
  this->Rmap22_ind = new int[length];
  this->Rmap22_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Rmap22_ind[ind] = PR1.Rmap22_ind[ind];
    this->Rmap22_fac[ind] = PR1.Rmap22_fac[ind];
  }

  this->Lmap1_ind = new int[PR1.L1_length];
  for(ind = 0; ind < PR1.L1_length; ++ind){
    this->Lmap1_ind[ind] = PR1.Lmap1_ind[ind];
  }

  this->Lmap21_num = new int[PR1.L2_1_length];
  this->Lmap21_index = new int[PR1.L2_1_length];
  length = 0;
  for(ind = 0; ind < PR1.L2_1_length; ++ind){
    this->Lmap21_num[ind] = PR1.Lmap21_num[ind];
    this->Lmap21_index[ind] = PR1.Lmap21_index[ind];
    length += PR1.Lmap21_num[ind];
  }
  this->Lmap21_ind = new int[length];
  this->Lmap21_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Lmap21_ind[ind] = PR1.Lmap21_ind[ind];
    this->Lmap21_fac[ind] = PR1.Lmap21_fac[ind];
  }

  this->Lmap22_num = new int[PR1.L2_2_length];
  this->Lmap22_index = new int[PR1.L2_2_length];
  length = 0;
  for(ind = 0; ind < PR1.L2_2_length; ++ind){
    this->Lmap22_num[ind] = PR1.Lmap22_num[ind];
    this->Lmap22_index[ind] = PR1.Lmap22_index[ind];
    length += PR1.Lmap22_num[ind];
  }
  this->Lmap22_ind = new int[length];
  this->Lmap22_fac = new double[length];
  for(ind = 0; ind < length; ++ind){
    this->Lmap22_ind[ind] = PR1.Lmap22_ind[ind];
    this->Lmap22_fac[ind] = PR1.Lmap22_fac[ind];
  }
}

void PR::Delete()
{
  delete[] R;
  delete[] R3;
  delete[] L;
  delete[] L3;

  delete[] hhp1_ind1;
  delete[] hhp1_ind2;
  delete[] hhp1_fac1;
  delete[] hhp1_fac0;

  delete[] chan31_num;
  delete[] chan32_num;

  delete[] R1_index;
  delete[] R2_1_index;
  delete[] R2_2_index;

  delete[] L1_index;
  delete[] L2_1_index;
  delete[] L2_2_index;

  delete[] Rmap1_ind;

  delete[] Rmap21_index;
  delete[] Rmap21_num;
  delete[] Rmap21_ind;
  delete[] Rmap21_fac;

  delete[] Rmap22_index;
  delete[] Rmap22_num;
  delete[] Rmap22_ind;
  delete[] Rmap22_fac;

  delete[] Lmap1_ind;

  delete[] Lmap21_index;
  delete[] Lmap21_num;
  delete[] Lmap21_ind;
  delete[] Lmap21_fac;

  delete[] Lmap22_index;
  delete[] Lmap22_num;
  delete[] Lmap22_ind;
  delete[] Lmap22_fac;
}

void PR::Zero_R()
{
  for(int ind = 0; ind < this->R_length; ++ind){ this->R[ind] = 0.0; }
  for(int ind = 0; ind < this->R3_length; ++ind){ this->R3[ind] = 0.0; }
}

void PR::Zero_L()
{
  for(int ind = 0; ind < this->L_length; ++ind){ this->L[ind] = 0.0; }
  for(int ind = 0; ind < this->L3_length; ++ind){ this->L3[ind] = 0.0; }
}

void PR::Update_R1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2)
{
  double p1 = 1.0, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int chan_ind1;
  char N = 'N';
  int chan32;
  int chan30 = this->chan30;
  int nh0 = this->nh0;
  int nhhp0 = this->nhhp0;
  int nph0 = Chan.nph1[Chan.ind0];
  int pind1, hind1, index;
  double *X2, *X3, *R2;

  if( nh0 != 0 ){
    chan_ind1 = Eff_Ints.Xhh.X_3_index[chan30];

    X3 = new double[nh0 * nh0];
    Eff_Ints.Xhh.Set_X_3(X3, nh0 * nh0, chan_ind1);
    //  w*r(|i){,i}  =  -r(|k){,k}.Xhh3(k|i){k,i}
    dgemm_NN(this->R, X3, PR2.R, &one, &nh0, &nh0, &m1, &p1, &N, &N);

    delete[] X3;
  }

  chan32 = Chan.size3 * Chan.ind0 + chan30;
  if( this->chan32_num[chan32] != 0 ){
    chan_ind1 = this->R2_1_index[chan32];

    X2 = new double[nph0]; // build (Xhp_2' = X(i|a){,ai'})
    for(int ph1 = 0; ph1 < nph0; ++ph1){
      pind1 = Chan.ph1_state(Chan.ind0, ph1).v1;
      hind1 = Chan.ph1_state(Chan.ind0, ph1).v2;
      index = Chan.hp1_map[Chan.ind0][Hash(hind1, pind1, 0)];
      X2[ph1] = Eff_Ints.Xhp.X_2[index];
    }
    R2 = new double[nph0 * nh0];
    this->Set_R2_1(R2, nph0 * nh0, chan_ind1);
    //  w*r(|i){,i}  =  Xhp2'(k|c){,ck'}.r(c|ik){ck',i}
    dgemm_NN(X2, R2, PR2.R, &one, &nh0, &nph0, &p1, &p1, &N, &N);
    
    delete[] X2;
    delete[] R2;
  }
   
  if( nhhp0 * nh0 != 0 ){
    chan_ind1 = Eff_Ints.Xhhhp.X_3_3_index[chan30];

    //  w*r(|i){,i}  =  -(1/2).r(c|kl){,klc'}.Xhhhp3_3(kl|ic){klc',i}
    dgemm_NN(this->R3, (Eff_Ints.Xhhhp.X_3_3 + chan_ind1), PR2.R, &one, &nh0, &nhhp0, &m12, &p1, &N, &N);
  }
}

void PR::Update_R2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int nh, np, nhh, nph1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';
  int chan31, chan32;
  int chan30 = this->chan30;
  int np0 = this->np0;
  int nh0 = this->nh0;
  int nhhp0 = this->nhhp0;
  int pind1, pind2, hind1, hind2, pj1, pj2, hj1, hj2, key1, key2, index;
  double *X2, *X3, *RR2, *R2, *RR1, *R1, *T3, *Q3;

  if( nh0 * nhhp0 != 0 ){
    chan_ind1 = Eff_Ints.Xhphh.X_3_1_index[chan30];
    
    X3 = new double[nh0 * nhhp0];
    Eff_Ints.Xhphh.Set_X_3_1(X3, nh0 * nhhp0, chan_ind1);

    //  w*r(a|ij){,ija'}  =  -r(|k){,k}.Xhphh3_1(ka|ij){k,ija'}
    dgemm_NN(this->R, X3, PR2.R3, &one, &nhhp0, &nh0, &m1, &p1, &N, &N);

    delete[] X3;
  }

  /*if( nh0 * npph0 != 0 ){
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan30];
    chan_ind2 = Amps.D1.T3_3_index[chan30];

    T3 = new double[npph0 * nh0];
    Amps.D1.Set_T3_3(T3, npph0 * nh0, chan_ind2);
    Q3 = new double[nh0];
    for(int ind = 0; ind < nh0; ++ind){ Q3[ind] = 0.0; }
    //w*r(ab|i){abi'}  =  -(1/2).T3_3(ab|ki){abi',k}.Vhhpp3_1(kl|cd){k,cdl'}.r(cd|l){cdl'}
    dgemm_NN((Ints.Vhhpp.V_3_1 + chan_ind1), this->R3, Q3, &nh0, &one, &npph0, &p1, &p1, &N, &N);
    dgemm_NN(T3, Q3, PA2.R3, &npph0, &one, &nh0, &m12, &p1, &N, &N);

    delete[] T3;
    delete[] Q3;
    }*/

  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    if( nph1 == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      nh = Chan.nh[chan3];
      if( nh == 0 ){ continue; }
      chan32 = Chan.size3 * chan2 + chan3;

      if( this->chan32_num[chan32] != 0 ){
	chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
	chan_ind2 = this->R2_1_index[chan32];
	chan_ind3 = this->R2_2_index[chan32];
	chan_ind4 = Eff_Ints.Xhphp.X_2_2_index[chan2];
	
	X3 = new double[nh * nh];
	Eff_Ints.Xhh.Set_X_3(X3, nh * nh, chan_ind1);
	
	X2 = new double[nph1 * nph1];
	Eff_Ints.Xhphp.Set_X_2_2(X2, nph1 * nph1, chan_ind4);

	RR2 = new double[nph1 * nh];
        this->Set_R2_2(RR2, nph1 * nh, chan_ind3);
	R2 = new double[nph1 * nh];
	for(int ind = 0; ind < nph1 * nh; ++ind){ R2[ind] = 0.0; }
	//  w*r(a|ij){aj',i}  =  r(a|jk){aj',k}.Xhh3(k|i){k,i}
	dgemm_NN(RR2, X3, R2, &nph1, &nh, &nh, &p1, &p1, &N, &N);
	//  w*r(a|ij){aj',i}  =  Xhphp2_2(ka|jc){aj',ck'}.r(c|ki){ck',i}
	dgemm_NN(X2, RR2, R2, &nph1, &nh, &nph1, &p1, &p1, &N, &N);
	PR2.Gather_R2_1(R2, nph1 * nh, chan_ind2);

	//  w*r(a|ij){ai',j}  =  -r(a|ik){ai',k}.Xhh3(k|j){k,j}
	dgemm_NN(RR2, X3, R2, &nph1, &nh, &nh, &m1, &p1, &N, &N);
	//  w*r(a|ij){ai',j}  =  -Xhphp2_2(ka|ic){ai',ck'}.r(c|kj){ck',j}
	dgemm_NN(X2, RR2, R2, &nph1, &nh, &nph1, &m1, &p1, &N, &N);
	PR2.Gather_R2_2(R2, nph1 * nh, chan_ind3);	

	delete[] X2;
	delete[] X3;
	delete[] R2;
	delete[] RR2;
      }
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    if( nhh == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      np = Chan.np[chan3];
      if( np == 0 ){ continue; }
      chan31 = Chan.size3 * chan1 + chan3;

      if( this->chan31_num[chan31] != 0 ){
	chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
	chan_ind2 = this->R1_index[chan31];
	chan_ind3 = Eff_Ints.Xhhhh.X_1_index[chan1];
	
	X3 = new double[np * np];
	Eff_Ints.Xpp.Set_X_3(X3, np * np, chan_ind1);
	RR1 = new double[nhh * np];
	this->Set_R1(RR1, nhh * np, chan_ind2);
	R1 = new double[nhh * np];
	for(int ind = 0; ind < nhh * np; ++ind){ R1[ind] = 0.0; }
	//  w*r(a|ij){a,ij}  =  Xpp3(a|c){a,c}.r(c|ij){c,ij}
	dgemm_NN(X3, RR1, R1, &np, &nhh, &np, &p1, &p1, &N, &N);
	//  w*r(a|ij){a,ij}  =  (1/2).r(a|kl){a,kl}.Xhhhh1(kl|ij){kl,ij}
	dgemm_NN(RR1, (Eff_Ints.Xhhhh.X_1 + chan_ind3), R1, &np, &nhh, &nhh, &p12, &p1, &N, &N);
	PR2.Gather_R1(R1, nhh * np, chan_ind2);	

	delete[] X3;
	delete[] R1;
	delete[] RR1;
      }
    }
  }
}

void PR::Update_L1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int chan_ind1;
  char N = 'N';
  int chan30 = this->chan30;
  int nh0 = this->nh0;
  int nhhp0 = this->nhhp0;
  double *X3;

  //  E*l(i|){i,}  =  E0.l(i|){i,}  
  double E0 = Amps.dE;
  for(int h = 0; h < nh0; ++h){
    PR2.L[h] += E0 * this->L[h];
  }

  if( nh0 != 0 ){
    chan_ind1 = Eff_Ints.Xhh.X_3_index[chan30];

    X3 = new double[nh0 * nh0];
    Eff_Ints.Xhh.Set_X_3(X3, nh0 * nh0, chan_ind1);
    //  E*l(i|){i,}  =  -Xhh3(i|k){i,k}.l(k|){k,}
    dgemm_NN(X3, this->L, PR2.L, &nh0, &one, &nh0, &m1, &p1, &N, &N);

    delete[] X3;
  }

  if( nh0 * nhhp0 != 0 ){
    chan_ind1 = Eff_Ints.Xhphh.X_3_1_index[chan30];
    
    X3 = new double[nh0 * nhhp0];
    Eff_Ints.Xhphh.Set_X_3_1(X3, nh0 * nhhp0, chan_ind1);
    //  E*l(i|){i,}  =  -(1/2).Xhphh3_1(ic|kl){i,klc'}.l(kl|c){klc',}
    dgemm_NN(X3, this->L3, PR2.L, &nh0, &one, &nhhp0, &m12, &p1, &N, &N);
    
    delete[] X3;
  }
}

void PR::Update_L2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2)
{
  double p1 = 1.0, p12 = 0.5, m1 = -1.0, m12 = -0.5;
  int one = 1;
  int nh, np, nhh, nph1;
  int chan_ind1, chan_ind2, chan_ind3, chan_ind4;
  char N = 'N';
  int chan31, chan32;
  int chan30 = this->chan30;
  int np0 = this->np0;
  int nh0 = this->nh0;
  int nhhp0 = this->nhhp0;
  int nph0 = Chan.nph1[Chan.ind0];
  int pind1, pind2, hind1, hind2, key1, key2, index;
  double *X2, *X3, *L2_1, *L2_2, *LL2, *L2, *LL1, *L1, *T3, *Q3;

  //  E*l(ij|a){ija',}  =  E0.l(ij|a){ija',}
  double E0 = Amps.dE;
  for(int hhp = 0; hhp < nhhp0; ++hhp){
    PR2.L3[hhp] += E0 * this->L3[hhp];
  }

  chan32 = Chan.size3 * Chan.ind0 + chan30;
  if( this->chan32_num[chan32] != 0 ){
    chan_ind1 = this->L2_1_index[chan32];
    chan_ind2 = this->L2_2_index[chan32];
    
    L2_1 = new double[nh0 * nph0];
    L2_2 = new double[nh0 * nph0];
    for(int ind = 0; ind < nh0 * nph0; ++ind){
      L2_1[ind] = 0.0;
      L2_2[ind] = 0.0;
    }

    X2 = new double[nph0]; // build (Xhp_2' = X(i|a){,ai'})
    for(int ph1 = 0; ph1 < nph0; ++ph1){
      pind1 = Chan.ph1_state(Chan.ind0, ph1).v1;
      hind1 = Chan.ph1_state(Chan.ind0, ph1).v2;
      index = Chan.hp1_map[Chan.ind0][Hash(hind1, pind1, 0)];
      X2[ph1] = Eff_Ints.Xhp.X_2[index];
    }

    for(int h = 0; h < nh0; ++h){
      for(int ph1 = 0; ph1 < nph0; ++ph1){
	//  E*l(ij|a){i,aj'}  =  l(i|){i,}.Xhp2'(i|b){ib'}
	L2_1[h*nph0 + ph1] += this->L[h] * X2[ph1];
	//  E*l(ij|a){j,ai'}  = -l(j|){j,}.Xhp2'(i|a){ia'}
	L2_2[h*nph0 + ph1] -= this->L[h] * X2[ph1];
      }
    }
    PR2.Gather_L2_1(L2_1, nh0*nph0, chan_ind1);
    PR2.Gather_L2_2(L2_2, nh0*nph0, chan_ind2);

    delete[] X2;
    delete[] L2_1;
    delete[] L2_2;
  }

  if( nhhp0 * nh0 != 0 ){
    chan_ind1 = Eff_Ints.Xhhhp.X_3_3_index[chan30];
    
    //  E*l(ij|a){ija',}  =  -Xhhhp3_3(ij|ka){ija',k}.l(k|){k,}
    dgemm_NN((Eff_Ints.Xhhhp.X_3_3 + chan_ind1), this->L, PR2.L3, &nhhp0, &one, &nh0, &m1, &p1, &N, &N);
  }

  /*if( nh0 * npph0 != 0 ){
    chan_ind1 = Ints.Vhhpp.V_3_1_index[chan30];
    chan_ind2 = Amps.D1.T3_3_index[chan30];
    
    T3 = new double[npph0 * nh0];
    Amps.D1.Set_T3_3(T3, npph0 * nh0, chan_ind2);
    Q3 = new double[nh0];
    for(int ind = 0; ind < nh0; ++ind){ Q3[ind] = 0.0; }
    //E*l(i|ab){abi'}  =  -(1/2).l(l|cd){cdl'}.T3_3(cd|kl){cdl',k}.Vhhpp3_1(ki|ab){k,abi'}
    dgemm_NN(this->L3, T3, Q3, &one, &nh0, &npph0, &p1, &p1, &N, &N);
    dgemm_NN(Q3, (Ints.Vhhpp.V_3_1 + chan_ind1), PA2.L3, &one, &npph0, &nh0, &m12, &p1, &N, &N);

    delete[] T3;
    delete[] Q3;
    }*/

  for(int chan2 = 0; chan2 < Chan.size2; ++chan2){
    nph1 = Chan.nph1[chan2];
    if( nph1 == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      nh = Chan.nh[chan3];
      if( nh == 0 ){ continue; }
      chan32 = Chan.size3 * chan2 + chan3;

      if( this->chan32_num[chan32] != 0 ){
	chan_ind1 = Eff_Ints.Xhh.X_3_index[chan3];
	chan_ind2 = this->L2_1_index[chan32];
	chan_ind3 = this->L2_2_index[chan32];
	chan_ind4 = Eff_Ints.Xhphp.X_2_2_index[chan2];
	
	X3 = new double[nh * nh];
	Eff_Ints.Xhh.Set_X_3(X3, nh * nh, chan_ind1);

	X2 = new double[nph1 * nph1];
	Eff_Ints.Xhphp.Set_X_2_2(X2, nph1 * nph1, chan_ind4);

	LL2 = new double[nph1 * nh];
	this->Set_L2_2(LL2, nph1 * nh, chan_ind3);
	L2 = new double[nph1 * nh];
	for(int ind = 0; ind < nph1 * nh; ++ind){ L2[ind] = 0.0; }
	//  E*l(ij|a){i,aj'}  =  Xhh3(i|k){i,k}.l(jk|a){k,aj'}
	dgemm_NN(X3, LL2, L2, &nh, &nph1, &nh, &p1, &p1, &N, &N);
	//  E*l(ij|a){i,aj'}  =  l(ki|c){i,ck'}.Xhphp2_2(jc|ka){ck',aj'}
	dgemm_NN(LL2, X2, L2, &nh, &nph1, &nph1, &p1, &p1, &N, &N);
	PR2.Gather_L2_1(L2, nh * nph1, chan_ind2);

	//  E*l(ij|a){j,ai'}  =  -Xhh3(i|k){j,k}.l(ik|a){k,ai'}
	dgemm_NN(X3, LL2, L2, &nh, &nph1, &nh, &m1, &p1, &N, &N);
	//  E*l(ij|a){j,ai'}  =  -l(kj|c){j,ck'}.Xhphp2_2(ic|ka){ck',ai'}
	dgemm_NN(LL2, X2, L2, &nh, &nph1, &nph1, &m1, &p1, &N, &N);
	PR2.Gather_L2_2(L2, nh * nph1, chan_ind3);

	delete[] X2;
	delete[] X3;
	delete[] L2;
	delete[] LL2;
      }
    }
  }

  for(int chan1 = 0; chan1 < Chan.size1; ++chan1){
    nhh = Chan.nhh[chan1];
    if( nhh == 0 ){ continue; }
    for(int chan3 = 0; chan3 < Chan.size3; ++chan3){
      np = Chan.np[chan3];
      if( np == 0 ){ continue; }
      chan31 = Chan.size3 * chan1 + chan3;

      if( this->chan31_num[chan31] != 0 ){
	chan_ind1 = Eff_Ints.Xpp.X_3_index[chan3];
	chan_ind2 = this->L1_index[chan31];
	chan_ind3 = Eff_Ints.Xhhhh.X_1_index[chan1];
	
	X3 = new double[np * np];
	Eff_Ints.Xpp.Set_X_3(X3, np*np, chan_ind1);
	LL1 = new double[nhh * np];
	this->Set_L1(LL1, np * nhh, chan_ind2);
	L1 = new double[nhh * np];
	for(int ind = 0; ind < np * nhh; ++ind){ L1[ind] = 0.0; }
	//  E*l(ij|a){ij,a}  =  l(ij|c){ij,c}.Xpp3(c|a){c,a}
	dgemm_NN(LL1, X3, L1, &nhh, &np, &np, &p1, &p1, &N, &N);
	//  E*l(ij|a){ij,a}  =  (1/2).Xhhhh1(ij|kl){ij,kl}.l(kl|a){kl,a}
	dgemm_NN((Eff_Ints.Xhhhh.X_1 + chan_ind3), LL1, L1, &nhh, &np, &nhh, &p12, &p1, &N, &N);
	PR2.Gather_L1(L1, nhh * np, chan_ind2);	

	delete[] X3;
	delete[] L1;
	delete[] LL1;
      }
    }
  }
}















void EOM::EOM_1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  std::cout << "1PA-EOM" << std::endl;
  int chan, N, np_tot, npph_tot, np0, npph0, npph;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int p1, p2, h1;

  // Count N_states and setup EOM structures //
  this->count_states(Chan, 1);
  np_tot = 0;
  npph_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan]; 
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (PAR.basis == "finite_J" && p1 == p2) ){ ++npph0; }
    }
    this->nob[n] = np0;
    this->nthb[n] = npph0;
    this->nstate[n] = np0 + npph0;
    this->ob_index[n] = np_tot;
    this->thb_index[n] = npph_tot;
    this->state_index[n] = np_tot + npph_tot;
    np_tot += np0;
    npph_tot += npph0;
  }
  this->ob_vec = new one_body[np_tot];
  this->thb_vec = new three_body[npph_tot];
  this->thb_qnums = new State[npph_tot];
  this->state_vec_R = new double[np_tot + npph_tot];
  this->state_vec_L = new double[np_tot + npph_tot];

  this->EOM_PA(Chan, Ints, Amps, Eff_Ints);
}

void EOM::EOM_2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  std::cout << "1PR-EOM" << std::endl;
  int chan, N, nh_tot, nhhp_tot, nh0, nhhp0, nhhp;
  double norm1p;
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
  int h1, h2, p1;

  // Count N_states and setup EOM structures //
  this->count_states(Chan, 0);
  nh_tot = 0;
  nhhp_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan]; 
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (PAR.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    this->nob[n] = nh0;
    this->nthb[n] = nhhp0;
    this->nstate[n] = nh0 + nhhp0;
    this->ob_index[n] = nh_tot;
    this->thb_index[n] = nhhp_tot;
    this->state_index[n] = nh_tot + nhhp_tot;
    nh_tot += nh0;
    nhhp_tot += nhhp0;
  }
  this->ob_vec = new one_body[nh_tot];
  this->thb_vec = new three_body[nhhp_tot];
  this->thb_qnums = new State[nhhp_tot];
  this->state_vec_R = new double[nh_tot + nhhp_tot];
  this->state_vec_L = new double[nh_tot + nhhp_tot];

  this->EOM_PR(Chan, Ints, Amps, Eff_Ints);
}

void EOM::EOM_PA(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  int nev = 1;
  int ncv, ldv, ldz, lworkl;
  bool rvec = true;
  char howmny = 'A', which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
  int mxiter = 5000;
  int mode, ishift, info, ido;
  double sigmar, sigmai, tol = 1.0e-10; // error tolerance
  int *select;
  double *resid, *workev, *workd, *workl, *v, *dr, *di, *z;
  int iparam[11], ipntr[14];
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
 
  State state;
  int length;
  int np0, npph0, npph, np_tot, npph_tot;
  three_body *pph_vec;
  State *J_state;
  int *J_vec;
  int *a_vec;
  int *b_vec;
  int *i_vec;
  int chan_j;
  int chan, N;
  double norm, norm1p;
  int p1, p2, h1;

  PA PA_Amps, PA_Amps2;
  int count0, key0;
  double fac0;

  // Count N_states and setup EOM structures //
  count_states(Chan, 1);
  np_tot = 0;
  npph_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    //PA_Amps.Setup(Chan, chan);
    //PA_Amps2.Setup(Chan, PA_Amps);
    np0 = Chan.np[chan];
    npph0 = 0;
    npph = Chan.npph[chan]; 
    for(int pph = 0; pph < npph; ++pph){
      p1 = Chan.pph_state(chan, pph).v1;
      p2 = Chan.pph_state(chan, pph).v2;
      if( (p1 < p2) || (PAR.basis == "finite_J" && p1 == p2) ){ ++npph0; }
      //++npph0;
    }
    nob[n] = np0;
    nthb[n] = npph0;
    nstate[n] = np0 + npph0;
    ob_index[n] = np_tot;
    thb_index[n] = npph_tot;
    state_index[n] = np_tot + npph_tot;
    np_tot += np0;
    npph_tot += npph0;
  }
  ob_vec = new one_body[np_tot];
  thb_vec = new three_body[npph_tot];
  thb_qnums = new State[npph_tot];
  state_vec_R = new double[np_tot + npph_tot];
  state_vec_L = new double[np_tot + npph_tot];

  for(int n = 0; n < this->N_states; ++n){
    chan = this->chan_vec[n];
    PA_Amps.Setup(Chan, chan);
    PA_Amps2.Setup(Chan, PA_Amps);

    N = PA_Amps.np0 + PA_Amps.npph1;
    eigenvalues = new double[1];
    eigenvectors_L = new double[N];
    eigenvectors_R = new double[N];
    ncv = 3*nev + 2;
    if(ncv > N){ ncv = N; }
    ldv = N;
    ldz = N;
    lworkl = 4*ncv*(ncv + 2);

    mode = 1;
    ishift = 1;
    info = 0;
    ido = 0; // status integer is zero at start
  
    select = new int[ncv];
    resid = new double[N];
    workev = new double[3 * ncv];
    workd = new double[3*N];
    workl = new double[lworkl];
    v = new double[N*ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    z = new double[N * (nev + 1)];

    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }

    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      for(int j = 0; j < PA_Amps.np0; ++j){ PA_Amps.R[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){
	PA_Amps.R3[PA_Amps.pph1_ind1[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j];
	PA_Amps.R3[PA_Amps.pph1_ind2[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j] * PA_Amps.pph1_fac1[j];
      }
      PA_Amps2.Zero_R();
      PA_Amps.Update_R1(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      PA_Amps.Update_R2(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      for(int j = 0; j < PA_Amps.np0; ++j){ workd[ipntr[1]-1 + j] = PA_Amps2.R[j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){ workd[ipntr[1]-1 + PA_Amps.np0 + j] = PA_Amps2.R3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j]; }
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    std::cout << "%%R  " << dr[0] << std::endl;
    //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[j] << " "; }
    //std::cout << std::endl;
    eigenvalues[0] = dr[0];
    for(int j = 0; j < PA_Amps.np0; ++j){ eigenvectors_R[j] = PA_Amps2.R[j]; }
    for(int j = 0; j < PA_Amps.npph1; ++j){ eigenvectors_R[PA_Amps.np0 + j] = PA_Amps2.R3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j]; }
    //for(int i = 0; i < N; ++i){ eigenvectors_R[i] = z[i]; }
    
    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;

    ////////////////////////////////////////////////////////////////////////////////

    mode = 1;
    ishift = 1;
    info = 0;
    ido = 0; // status integer is zero at start
  
    select = new int[ncv];
    resid = new double[N];
    workev = new double[3 * ncv];
    workd = new double[3*N];
    workl = new double[lworkl];
    v = new double[N*ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    z = new double[N * (nev + 1)];

    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }

    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      for(int j = 0; j < PA_Amps.np0; ++j){ PA_Amps.L[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){
	PA_Amps.L3[PA_Amps.pph1_ind1[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j];
	PA_Amps.L3[PA_Amps.pph1_ind2[j]] = workd[ipntr[0]-1 + PA_Amps.np0 + j] * PA_Amps.pph1_fac0[j] * PA_Amps.pph1_fac1[j];
      }
      PA_Amps2.Zero_L();
      PA_Amps.Update_L1(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      PA_Amps.Update_L2(Chan, Ints, Amps, Eff_Ints, PA_Amps2);
      for(int j = 0; j < PA_Amps.np0; ++j){ workd[ipntr[1]-1 + j] = PA_Amps2.L[j]; }
      for(int j = 0; j < PA_Amps.npph1; ++j){ workd[ipntr[1]-1 + PA_Amps.np0 + j] = PA_Amps2.L3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j]; }
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    std::cout << "%%L  " << dr[0] - Amps.dE << std::endl;
    //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[j] << " "; }
    //std::cout << std::endl;
    for(int j = 0; j < PA_Amps.np0; ++j){ eigenvectors_L[j] = PA_Amps2.L[j]; }
    for(int j = 0; j < PA_Amps.npph1; ++j){ eigenvectors_L[PA_Amps.np0 + j] = PA_Amps2.L3[PA_Amps.pph1_ind1[j]] / PA_Amps.pph1_fac0[j]; }
    //for(int i = 0; i < N; ++i){ eigenvectors_L[i] = z[i]; }
    
    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;

    norm = 0.0;
    for(int i = 0; i < N; ++i){ norm += eigenvectors_L[i] * eigenvectors_R[i]; }
    if( norm < 0.0 ){
      for(int i = 0; i < N; ++i){ eigenvectors_R[i] *= -1.0; }
      norm *= -1.0;
    }
    norm = std::sqrt(norm);
    for(int i = 0; i < N; ++i){
      eigenvectors_L[i] /= norm;
      eigenvectors_R[i] /= norm;
      if( std::fabs(eigenvectors_L[i]) < 1.0e-10 ){ eigenvectors_L[i] = 0.0; }
      if( std::fabs(eigenvectors_R[i]) < 1.0e-10 ){ eigenvectors_R[i] = 0.0; }
    }

    norm1p = 0.0;
    for(int i = 0; i < PA_Amps.np0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);
    
    this->del_E[n] = eigenvalues[0];
    this->norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      this->state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      this->state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }
    
    PA_Amps.Delete();
    PA_Amps2.Delete();

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
  }
}


void EOM::EOM_PR(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints)
{
  int nev = 1;
  int ncv, ldv, ldz, lworkl;
  bool rvec = true;
  char howmny = 'A', which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
  int mxiter = 5000;
  int mode, ishift, info, ido;
  double sigmar, sigmai, tol = 1.0e-10; // error tolerance
  int *select;
  double *resid, *workev, *workd, *workl, *v, *dr, *di, *z;
  int iparam[11], ipntr[14];
  double *eigenvalues, *eigenvectors_L, *eigenvectors_R;
 
  State state;
  int length;
  int nh0, nhhp0, nhhp, nh_tot, nhhp_tot;
  three_body *hhp_vec;
  State *J_state;
  int *J_vec;
  int *i_vec;
  int *j_vec;
  int *a_vec;
  int chan_j;
  int chan, N;
  double norm, norm1p;
  int h1, h2, p1;

  PR PR_Amps, PR_Amps2;
  int count0, key0;
  double fac0;

  // Count N_states and setup EOM structures //
  count_states(Chan, 0);
  nh_tot = 0;
  nhhp_tot = 0;
  for(int n = 0; n < N_states; ++n){
    chan = chan_vec[n];
    nh0 = Chan.nh[chan];
    nhhp0 = 0;
    nhhp = Chan.nhhp[chan]; 
    for(int hhp = 0; hhp < nhhp; ++hhp){
      h1 = Chan.hhp_state(chan, hhp).v1;
      h2 = Chan.hhp_state(chan, hhp).v2;
      if( (h1 < h2) || (PAR.basis == "finite_J" && h1 == h2) ){ ++nhhp0; }
    }
    nob[n] = nh0;
    nthb[n] = nhhp0;
    nstate[n] = nh0 + nhhp0;
    ob_index[n] = nh_tot;
    thb_index[n] = nhhp_tot;
    state_index[n] = nh_tot + nhhp_tot;
    nh_tot += nh0;
    nhhp_tot += nhhp0;
  }
  ob_vec = new one_body[nh_tot];
  thb_vec = new three_body[nhhp_tot];
  thb_qnums = new State[nhhp_tot];
  state_vec_R = new double[nh_tot + nhhp_tot];
  state_vec_L = new double[nh_tot + nhhp_tot];

  for(int n = 0; n < this->N_states; ++n){
    chan = this->chan_vec[n];
    PR_Amps.Setup(Chan, chan);
    PR_Amps2.Setup(Chan, PR_Amps);

    N = PR_Amps.nh0 + PR_Amps.nhhp1;
    eigenvalues = new double[1];
    eigenvectors_L = new double[N];
    eigenvectors_R = new double[N];
    ncv = 3*nev + 2;
    if(ncv > N){ ncv = N; }
    ldv = N;
    ldz = N;
    lworkl = 4*ncv*(ncv + 2);

    mode = 1;
    ishift = 1;
    info = 0;
    ido = 0; // status integer is zero at start
  
    select = new int[ncv];
    resid = new double[N];
    workev = new double[3 * ncv];
    workd = new double[3*N];
    workl = new double[lworkl];
    v = new double[N*ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    z = new double[N * (nev + 1)];

    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }

    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.R[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.R3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.R3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_R();
      PR_Amps.Update_R1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_R2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      for(int j = 0; j < PR_Amps.nh0; ++j){ workd[ipntr[1]-1 + j] = PR_Amps2.R[j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ workd[ipntr[1]-1 + PR_Amps.nh0 + j] = PR_Amps2.R3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j]; }
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    std::cout << "%%R  " << dr[0] << std::endl;
    //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[j] << " "; }
    //std::cout << std::endl;
    eigenvalues[0] = dr[0];
    for(int j = 0; j < PR_Amps.nh0; ++j){ eigenvectors_R[j] = PR_Amps2.R[j]; }
    for(int j = 0; j < PR_Amps.nhhp1; ++j){ eigenvectors_R[PR_Amps.nh0 + j] = PR_Amps2.R3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j]; }
    //for(int i = 0; i < N; ++i){ eigenvectors_R[i] = z[i]; }
    
    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;

    ////////////////////////////////////////////////////////////////////////////////

    mode = 1;
    ishift = 1;
    info = 0;
    ido = 0; // status integer is zero at start
  
    select = new int[ncv];
    resid = new double[N];
    workev = new double[3 * ncv];
    workd = new double[3*N];
    workl = new double[lworkl];
    v = new double[N*ncv];
    dr = new double[nev + 1];
    di = new double[nev + 1];
    z = new double[N * (nev + 1)];

    iparam[0] = ishift;
    iparam[2] = mxiter;
    iparam[6] = mode;
    for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
    for(int i = 0; i < nev + 1; ++i){
      dr[i] = 0.0;
      di[i] = 0.0;
      for(int j = 0; j < N; ++j){ z[N*i + j] = 0.0; }
    }

    do{
      dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
      if(ido != -1 && ido != 1){ break; }
      for(int j = 0; j < PR_Amps.nh0; ++j){ PR_Amps.L[j] = workd[ipntr[0]-1 + j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){
	PR_Amps.L3[PR_Amps.hhp1_ind1[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j];
	PR_Amps.L3[PR_Amps.hhp1_ind2[j]] = workd[ipntr[0]-1 + PR_Amps.nh0 + j] * PR_Amps.hhp1_fac0[j] * PR_Amps.hhp1_fac1[j];
      }
      PR_Amps2.Zero_L();
      PR_Amps.Update_L1(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      PR_Amps.Update_L2(Chan, Ints, Amps, Eff_Ints, PR_Amps2);
      for(int j = 0; j < PR_Amps.nh0; ++j){ workd[ipntr[1]-1 + j] = PR_Amps2.L[j]; }
      for(int j = 0; j < PR_Amps.nhhp1; ++j){ workd[ipntr[1]-1 + PR_Amps.nh0 + j] = PR_Amps2.L3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j]; }
    } while(true);
    dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    std::cout << "%%L  " << dr[0] - Amps.dE << std::endl;
    //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[j] << " "; }
    //std::cout << std::endl;
    for(int j = 0; j < PR_Amps.nh0; ++j){ eigenvectors_L[j] = PR_Amps2.L[j]; }
    for(int j = 0; j < PR_Amps.nhhp1; ++j){ eigenvectors_L[PR_Amps.nh0 + j] = PR_Amps2.L3[PR_Amps.hhp1_ind1[j]] / PR_Amps.hhp1_fac0[j]; }
    //for(int i = 0; i < N; ++i){ eigenvectors_L[i] = z[i]; }
    
    delete[] select;
    delete[] resid;
    delete[] workev;
    delete[] workd;
    delete[] workl;
    delete[] v;
    delete[] dr;
    delete[] di;
    delete[] z;

    norm = 0.0;
    for(int i = 0; i < N; ++i){ norm += eigenvectors_L[i] * eigenvectors_R[i]; }
    if( norm < 0.0 ){
      for(int i = 0; i < N; ++i){ eigenvectors_R[i] *= -1.0; }
      norm *= -1.0;
    }
    norm = std::sqrt(norm);
    for(int i = 0; i < N; ++i){
      eigenvectors_L[i] /= norm;
      eigenvectors_R[i] /= norm;
      if( std::fabs(eigenvectors_L[i]) < 1.0e-10 ){ eigenvectors_L[i] = 0.0; }
      if( std::fabs(eigenvectors_R[i]) < 1.0e-10 ){ eigenvectors_R[i] = 0.0; }
    }

    norm1p = 0.0;
    for(int i = 0; i < PR_Amps.nh0; ++i){ norm1p += eigenvectors_L[i] * eigenvectors_R[i]; }
    norm1p = std::sqrt(norm1p);
    
    this->del_E[n] = eigenvalues[0];
    this->norm_1p[n] = norm1p;
    for(int i = 0; i < N; ++i){
      this->state_vec_R[state_index[n] + i] = eigenvectors_R[i];
      this->state_vec_L[state_index[n] + i] = eigenvectors_L[i];
    }
    
    PR_Amps.Delete();
    PR_Amps2.Delete();

    delete[] eigenvalues;
    delete[] eigenvectors_L;
    delete[] eigenvectors_R;
  }
}
