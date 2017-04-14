#include "CCfunctions.hpp"
#include "HFfunctions.hpp"
#include "MATHfunctions.hpp"
#include "INTfunctions.hpp"
#include "TESTfunctions.hpp"
#include "BASISfunctions.hpp"

void Doubles_Step_explicit(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps1, Amplitudes &Amps2, double &error)
{
  State tb;
  int ind, ind1, ind2, key1, key2, key3, key4;
  int nh, np, nhh, npp, nhp, nhp1, nhp2, nhpp, nhhp;
  int a, b, i, j, k, l, c, d;
  double tempt, term0, term;
  double norm2 = 0.0;
  double error2 = 0.0;
  int print = 0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    nhp = Chan.nhp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }

	//if(i == 2 && j == 3 && a == 6 && b == 11){ print = 1; }
	//else{ print = 0; }
	//if(i < j && a < b){ print = 1; }
	//else{ print = 0; }

	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	  std::cout << "Term0 = < " << a << "," << b << " |V| " << i << "," << j << " > = " << Ints.D_ME1.V4[chan][pp * nhh + hh] << std::endl;
	}
	tempt = Ints.D_ME1.V4[chan][pp * nhh + hh];

	//T1(ab|ij){ij,ab} = 0.5 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  if(c == d){ continue; }
	  term0 = 0.5 * Ints.D_ME1.V1[chan][pp1 * npp + pp] * Amps1.D1.T1[chan][hh * npp + pp1];
	  term += term0;
	  if(print == 2){ std::cout << "Term1 += 0.5 * < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.D_ME1.V1[chan][pp1 * npp + pp] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	tempt += term;

	//T1(ab|ij){ij,ab} = 0.5 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int hh1 = 0; hh1 < nhh; ++hh1){
	  k = Chan.hhvec[chan][2*hh1];
	  l = Chan.hhvec[chan][2*hh1 + 1];
	  if(k == l){ continue; }
	  term0 = 0.5 * Ints.D_ME1.V2[chan][hh * nhh + hh1] * Amps1.D1.T1[chan][hh1 * npp + pp];
	  term += term0;
	  if(print == 2){ std::cout << "Term2 += 0.5 * < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.D_ME1.V2[chan][hh * nhh + hh1] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	tempt += term;

	//T1(ab|ij){ij,ab} = 0.25 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int pp1 = 0; pp1 < npp; ++pp1){
	  c = Chan.ppvec[chan][2*pp1];
	  d = Chan.ppvec[chan][2*pp1 + 1];
	  if(c == d){ continue; }
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    term0 = 0.25 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.D1.T1[chan][hh * npp + pp1] * Amps1.D1.T1[chan][hh1 * npp + pp];
	    term += term0;
	    if(print == 2){ std::cout << "Term7 += 0.25 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.25 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	tempt += term;

	//T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || i == k){ continue; }
	  term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term3 += -1.0 * < " << k << "," << b << " |V| " << j << "," << c << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	tempt += term;

	//T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[j], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(j, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(i, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == c || j == k){ continue; }
	  term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){  std::cout << "Term4 += -1.0 * < " << k << "," << a << " |V| " << i << "," << c << " > * < " << b << "," << c << " |t| " << j << "," << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	tempt += term;

	//T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == c || i == k){ continue; }
	  term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term5 += < " << k << "," << a << " |V| " << j << "," << c << " > * < " << b << "," << c << " |t| " << i << "," << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	tempt += term;

	//T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[j], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(j, a, Space.indtot)];
 	key2 = Chan.hp2_map[ind][Hash2(i, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || j == k){ continue; }
	  term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	  term += term0;
	  if(print == 2){ std::cout << "Term6 += < " << k << "," << b << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << j << "," << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key2] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	}
	if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	tempt += term;

	//T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[a]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  c = Chan.hp2vec[ind][2*hp2 + 1];
	  if(a == c || i == k){ continue; }
	  for(int hp1 = 0; hp1 < nhp1; ++hp1){
	    l = Chan.hp1vec[ind][2*hp1];
	    d = Chan.hp1vec[ind][2*hp1 + 1];
	    if(d == b || l == j || k == l || c == d){ continue; }
	    term0 = Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.D1.T2[ind][key1 * nhp2 + hp2] *  Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	tempt += term;

	//T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	minus(tb, Space.qnums[i], Space.qnums[b]);
	ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	nhp1 = Chan.nhp1[ind];
	nhp2 = Chan.nhp2[ind];
	for(int hp2 = 0; hp2 < nhp2; ++hp2){
	  k = Chan.hp2vec[ind][2*hp2];
	  d = Chan.hp2vec[ind][2*hp2 + 1];
	  if(b == d || i == k){ continue; }
	  for(int hp1 = 0; hp1 < nhp1; ++hp1){
	    l = Chan.hp1vec[ind][2*hp1];
	    c = Chan.hp1vec[ind][2*hp1 + 1];
	    if(c == a || l == j || k == l || c == d){ continue; }
	    term0 = Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.D1.T2[ind][key1 * nhp2 + hp2] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term13 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << d << " |t| " << i << "," << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	tempt += term;

	//T6(ab|ij){jab,i} = -0.5 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[i];
	key1 = Chan.h_map[ind][i];
	key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	nh = Chan.nh[ind];
	nhpp = Chan.nhpp[ind];
	for(int h = 0; h < nh; ++h){
	  l = Chan.hvec[ind][h];
	  if(j == l){ continue; }
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(i == k || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V5[ind][h * nhpp + hpp] * Amps1.D1.T7[ind][key2 * nh + h] * Amps1.D1.T6[ind][hpp * nh + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << b << " |t| " << j << "," << l << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = -0.5 * " << Ints.D_ME1.V5[ind][h * nhpp + hpp] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " * " << Amps1.D1.T6[ind][hpp * nh + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	tempt += term;

	//T7(ab|ij){iab,j} = -0.5 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[j];
	key1 = Chan.h_map[ind][j];
	key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	nh = Chan.nh[ind];
	nhpp = Chan.nhpp[ind];
	for(int h = 0; h < nh; ++h){
	  k = Chan.hvec[ind][h];
	  if(i == k){ continue; }
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    l = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(j == l || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V6[ind][h * nhpp + hpp] * Amps1.D1.T7[ind][key2 * nh + h] * Amps1.D1.T6[ind][hpp * nh + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term9 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << b << " |t| " << i << "," << k << " > * < " << c << "," << d << " |t| " << j << "," << l << " > = -0.5 * " << Ints.D_ME1.V6[ind][h * nhpp + hpp] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " * " << Amps1.D1.T6[ind][hpp * nh + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	tempt += term;

	//T8(ab|ij){ijb,a} = -0.5 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[a];
	key1 = Chan.p_map[ind][a];
	key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	np = Chan.np[ind];
	nhhp = Chan.nhhp[ind];
	for(int p = 0; p < np; ++p){
	  d = Chan.pvec[ind][p];
	  if(b == d){ continue; }
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(a == c || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V7[ind][p * nhhp + hhp] * Amps1.D1.T9[ind][key2 * np + p] * Amps1.D1.T8[ind][hhp * np + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term10 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V7[ind][p * nhhp + hhp] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " * " << Amps1.D1.T8[ind][hhp * np + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	tempt += term;

	//T9(ab|ij){ija,b} = -0.5 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	ind = Chan.indvec[b];
	key1 = Chan.p_map[ind][b];
	key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	np = Chan.np[ind];
	nhhp = Chan.nhhp[ind];
	for(int p = 0; p < np; ++p){
	  c = Chan.pvec[ind][p];
	  if(a == c){ continue; }
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    d = Chan.hhpvec[ind][3*hhp + 2];
	    if(b == d || k == l || c == d){ continue; }
	    term0 = -0.5 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.D1.T9[ind][key2 * np + p] * Amps1.D1.T8[ind][hhp * np + key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term11 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << j << " > * < " << b << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " * " << Amps1.D1.T8[ind][hhp * np + key1] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	tempt += term;
	
	if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
	  if(print != 0){ std::cout << std::endl << "For Singles" << std::endl << std::endl; }

	  //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} = t1(c|i){ic}.t1(d|j){jd}.V1(ab|cd){cd,ab} (1)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    term0 = Ints.D_ME1.V1[chan][pp1 * npp + pp] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term1 += < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > = " << Ints.D_ME1.V1[chan][pp1 * npp + pp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	  tempt += term;
	  
	  //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} = V2(ij|kl){ij,kl}.t1(a|k){ka}.t1(b|l){lb} (2)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	    term0 = Ints.D_ME1.V2[chan][hh * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << Ints.D_ME1.V2[chan][hh * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	  tempt += term;

	  //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} = 0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.t1(a|k){ka}.t1(b|l){lb} (3)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      term0 = 0.5 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term3 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	  tempt += term;
	
	  //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} = 0.5 * t1(c|i){ic}.t1(d|j){jd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l){ continue; }
	      term0 = 0.5 * Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term4 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	  tempt += term;
	
	  //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} = t1(c|i){ic}.t1(d|j){jd}.V4(kl|cd){cd,kl}.t1(a|k){ka}.t1(b|l){lb} (5)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	    for(int hh1 = 0; hh1 < nhh; ++hh1){
	      k = Chan.hhvec[chan][2*hh1];
	      l = Chan.hhvec[chan][2*hh1 + 1];
	      if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      term0 = Ints.D_ME1.V4[chan][pp1 * nhh + hh1] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	      term += term0;
	      if(print == 2){ std::cout << "Term5 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << Ints.D_ME1.V4[chan][pp1 * nhh + hh1] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} = t1(c|i){ic}.t1(a|k){ka}.V3(kb|jc){kc,jb} (6)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term6 += < " << k << "," << b << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	  tempt += term;
	
	  //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} = t1(c|j){jc}.t1(b|k){kb}.V3(ka|ic){kc,ia}  (7)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(i, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term7 += < " << k << "," << a << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > = " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} = -t1(c|i){ic}.t1(b|k){kb}.V3(ka|jc){kc,ja}  (8)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term8 += -1.0 * < " << k << "," << a << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	  tempt += term;
	
	  //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} = -t1(c|j){jc}.t1(a|k){ka}.V3(kb|ic){kc,ib}  (9)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp2_map[ind][Hash2(i, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = -1.0 * Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << b << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[ind][hp2 * nhp2 + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} = -t1(c|i){ic}.t1(a|k){ka}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(i, a, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(j, b, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      d = Chan.hp1vec[ind][2*hp1 + 1];
	      if(d == b || l == j || k == l || c == d){ continue; }
	      term0 = -1.0 * Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term10 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.t1(b|l){lb}.t1(d|j){jd} (11)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    c = Chan.hp2vec[ind][2*hp2 + 1];
	    if(a == c || i == k){ continue; }
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      d = Chan.hp1vec[ind][2*hp1 + 1];
	      if(k == l || c == d || Chan.indvec[l] != Chan.indvec[b] || Chan.indvec[j] != Chan.indvec[d]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      term0 = -1.0 * Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	      term += term0;
	      if(print == 2){ std::cout << "Term11 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << " |t| " << l << " > * < " << d << " |t| " << j << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V9[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} = -t1(d|i){id}.t1(b|k){kb}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key1 = Chan.hp1_map[ind][Hash2(i, b, Space.indtot)];
	  key2 = Chan.hp2_map[ind][Hash2(j, a, Space.indtot)];
	  nhp1 = Chan.nhp1[ind];
	  nhp2 = Chan.nhp2[ind];
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    d = Chan.hp2vec[ind][2*hp2 + 1];
	    if(Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      c = Chan.hp1vec[ind][2*hp1 + 1];
	      if(c == a || l == j || k == l || c == d){ continue; }
	      term0 = -1.0 * Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][hp1 * nhp2 + key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term12 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = -1.0 * " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][hp1 * nhp2 + key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	  tempt += term;
	
	  //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.t1(a|l){la}.t1(c|j){jc} (13)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int hp2 = 0; hp2 < nhp2; ++hp2){
	    k = Chan.hp2vec[ind][2*hp2];
	    d = Chan.hp2vec[ind][2*hp2 + 1];
	    if(b == d || i == k){ continue; }
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      l = Chan.hp1vec[ind][2*hp1];
	      c = Chan.hp1vec[ind][2*hp1 + 1];
	      if(k == l || c == d || Chan.indvec[l] != Chan.indvec[a] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      term0 = -1.0 * Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4] * Amps1.D1.T2[ind][key1 * nhp2 + hp2];
	      term += term0;
	      if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << l << " > * < " << c << " |t| " << j << " > * < " << b << "," << d << " |t| " << i << "," << k << " > = -1.0 * " << Ints.D_ME1.V10[ind][hp2 * nhp1 + hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " * " << Amps1.D1.T2[ind][key1 * nhp2 + hp2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	  tempt += term;
	
	  //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} = -V11(ka|cd){a,kcd}.t1(c|i){ic}.T2(db|kj){kd,jb} (14)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == b || k == j || Chan.indvec[i] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term14 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << b << " |t| " << k << "," << j << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term14 = " << term << std::endl; }
	  tempt += term;

	  //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} = -V11(kb|cd){b,kcd}.t1(c|j){jc}.T2(da|ki){kd,ia} (15)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == a || k == i || Chan.indvec[j] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term15 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << a << " |t| " << k << "," << i << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term15 = " << term << std::endl; }
	  tempt += term;

	  //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} = V11(kb|cd){b,kcd}.t1(c|i){ic}.T2(da|kj){kd,ja} (16)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == a || k == j || Chan.indvec[i] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term16 += < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << a << " |t| " << k << "," << j << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term16 = " << term << std::endl; }
	  tempt += term;

	  //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} = V11(ka|cd){a,kcd}.t1(c|j){jc}.T2(db|ki){kd,ib} (17)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    minus(tb, Space.qnums[k], Space.qnums[d]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(c == d || d == b || k == i || Chan.indvec[j] != Chan.indvec[c] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(k, d, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term17 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << b << " |t| " << k << "," << i << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term17 = " << term << std::endl; }
	  tempt += term;

	  //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} = -V12(kl|ic){i,klc}.t1(a|k){ka}.T2(cb|lj){lc,jb} (18)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == b || l == j || Chan.indvec[k] != Chan.indvec[a] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term18 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term18 = " << term << std::endl; }
	  tempt += term;

	  //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} = -V12(kl|jc){j,klc}.t1(b|k){kb}.T2(ca|li){lc,ia} (19)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == a || l == i || Chan.indvec[k] != Chan.indvec[b] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term19 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term19 = " << term << std::endl; }
	  tempt += term;
	  
	  //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} = V12(kl|ic){i,klc}.t1(b|k){kb}.T2(ca|lj){lc,ja} (20)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[i], Space.qnums[b]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(j, a, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == a || l == j || Chan.indvec[k] != Chan.indvec[b] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term20 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term20 = " << term << std::endl; }
	  tempt += term;

	  //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} = V12(kl|jc){j,klc}.t1(a|k){ka}.T2(cb|li){lc,ib} (21)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  minus(tb, Space.qnums[j], Space.qnums[a]);
	  ind1 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	  key2 = Chan.hp2_map[ind1][Hash2(i, b, Space.indtot)];
	  nhp2 = Chan.nhp2[ind1];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    minus(tb, Space.qnums[l], Space.qnums[c]);
	    ind2 = ChanInd_2b_cross(Parameters.basis, Space, tb);
	    if(k == l || c == b || l == i || Chan.indvec[k] != Chan.indvec[a] || ind2 != ind1){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key4 = Chan.hp1_map[ind2][Hash2(l, c, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.D1.T2[ind1][key4 * nhp2 + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term21 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << i << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T2[ind1][key4 * nhp2 + key2] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term21 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.t1(d|i){id}.t1(c|k){kc} (22)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == j){ continue; }
	    for(int p = 0; p < np; ++p){
	      d = Chan.pvec[ind][p];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hpp_map[ind][Hash3(k, c, d, Space.indtot)];
		term0 = Ints.D_ME1.V5[ind][h * nhpp + key4] * Amps1.S1.T1[key3] * Amps1.S1.T1[hp1] * Amps1.D1.T7[ind][key2 * nh + h];
		term += term0;
		if(print == 2){ std::cout << "Term22 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << j << "," << l << " > = " << Ints.D_ME1.V5[ind][h * nhpp + key4] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term22 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.t1(d|j){jd}.t1(c|l){lc} (23)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    if(k == i){ continue; }
	    for(int p = 0; p < np; ++p){
	      d = Chan.pvec[ind][p];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		l = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hpp_map[ind][Hash3(l, c, d, Space.indtot)];
		term0 = Ints.D_ME1.V6[ind][h * nhpp + key4] * Amps1.S1.T1[key3] * Amps1.S1.T1[hp1] * Amps1.D1.T7[ind][key2 * nh + h];
		term += term0;
		if(print == 2){ std::cout << "Term23 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << j << " > * < " << c << " |t| " << l << " > * < " << a << "," << b << " |t| " << i << "," << k << " > = " << Ints.D_ME1.V6[ind][h * nhpp + key4] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T7[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term23 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == j){ continue; }
	    key3 = Chan.hh1_map[Chan.ind0][Hash2(i, l, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(k == l){ continue; }
	      term0 = Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T6[ind][key2 * nh + h];
	      term += term0;
	      if(print == 2){ std::cout << "Term26 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = " << Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T6[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term26 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    l = Chan.hvec[ind][h];
	    if(l == i){ continue; }
	    key3 = Chan.hh1_map[Chan.ind0][Hash2(j, l, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(k == l){ continue; }
	      term0 = -1.0 * Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T6[ind][key2 * nh + h];
	      term += term0;
	      if(print == 2){ std::cout << "Term27 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = -1.0 * " << Ints.S_ME1.V15[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T6[ind][key2 * nh + h] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term27 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  key2 = Chan.hpp_map[ind][Hash3(j, a, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V17[ind][key2 * np + p] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term30 += -1.0 * < " << j << "," << c << " |V| " << a << "," << b << " > * < " << c << " |t| " << i << " > = -1.0 * " << Ints.S_ME1.V17[ind][key2 * np + p] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term30 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  key2 = Chan.hpp_map[ind][Hash3(i, a, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhpp = Chan.nhpp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = Ints.S_ME1.V17[ind][key2 * np + p] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term31 += < " << i << "," << c << " |V| " << a << "," << b << " > * < " << c << " |t| " << j << " > = " << Ints.S_ME1.V17[ind][key2 * np + p] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term31 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  np = Chan.np[ind];
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    for(int p = 0; p < np; ++p){
	      c = Chan.pvec[ind][p];
	      plus(tb, Space.qnums[j], Space.qnums[c]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(j, c, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      term0 = -0.5 * Ints.S_ME1.V19[chan][key1 * nhh + hh1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term34 += -0.5 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = -0.5 * " << Ints.S_ME1.V19[chan][key1 * nhh + hh1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term34 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  np = Chan.np[ind];
	  for(int hh1 = 0; hh1 < nhh; ++hh1){
	    k = Chan.hhvec[chan][2*hh1];
	    l = Chan.hhvec[chan][2*hh1 + 1];
	    if(k == l){ continue; }
	    for(int p = 0; p < np; ++p){
	      c = Chan.pvec[ind][p];
	      plus(tb, Space.qnums[i], Space.qnums[c]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(i, c, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      term0 = 0.5 * Ints.S_ME1.V19[chan][key1 * nhh + hh1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh1 * npp + pp];
	      term += term0;
	      if(print == 2){ std::cout << "Term35 += 0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << Ints.S_ME1.V19[chan][key1 * nhh + hh1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh1 * npp + pp] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term35 = " << term << std::endl; }
	  tempt += term;

	  //T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[j];
	  key1 = Chan.h_map[ind][j];
	  nhhp = Chan.nhhp[ind];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(k == l || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[l] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	    term0 = Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term38 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << a << " |t| " << l << " > = " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term38 = " << term << std::endl; }
	  tempt += term;

	  //T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[i];
	  key1 = Chan.h_map[ind][i];
	  nhhp = Chan.nhhp[ind];
	  for(int hhp = 0; hhp < nhhp; ++hhp){
	    k = Chan.hhpvec[ind][3*hhp];
	    l = Chan.hhpvec[ind][3*hhp + 1];
	    c = Chan.hhpvec[ind][3*hhp + 2];
	    if(k == l || Chan.indvec[j] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[l] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term39 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > * < " << a << " |t| " << l << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term39 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == b){ continue; }
	    for(int h = 0; h < nh; ++h){
	      l = Chan.hvec[ind][h];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hhp_map[ind][Hash3(k, l, c, Space.indtot)];
		term0 = Ints.D_ME1.V7[ind][p * nhhp + key4] * Amps1.S1.T1[hp1] * Amps1.S1.T1[key3] * Amps1.D1.T9[ind][key2 * np + p];
		term += term0;
		if(print == 2){ std::cout << "Term24 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << a << " |t| " << l << " > * < " << b << "," << d << " |t| " << i << "," << j << " > = " << Ints.D_ME1.V7[ind][p * nhhp + key4] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term24 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  nh = Chan.nh[ind];
	  np = Chan.np[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    c = Chan.pvec[ind][p];
	    if(c == a){ continue; }
	    for(int h = 0; h < nh; ++h){
	      l = Chan.hvec[ind][h];
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int hp1 = 0; hp1 < nhp1; ++hp1){
		k = Chan.hp1vec[Chan.ind0][2*hp1];
		d = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
		if(k == l || c == d){ continue; }
		key4 = Chan.hhp_map[ind][Hash3(k, l, d, Space.indtot)];
		term0 = Ints.D_ME1.V8[ind][p * nhhp + key4] * Amps1.S1.T1[hp1] * Amps1.S1.T1[key3] * Amps1.D1.T9[ind][key2 * np + p];
		term += term0;
		if(print == 2){ std::cout << "Term25 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << a << "," << c << " |t| " << i << "," << j << " > = " << Ints.D_ME1.V8[ind][p * nhhp + key4] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.D1.T9[ind][key2 * np + p] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term25 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  np = Chan.np[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == b){ continue; }
	    key3 = Chan.pp1_map[Chan.ind0][Hash2(d, a, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(c == d){ continue; }
	      term0 = Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T8[ind][key2 * np + p];
	      term += term0;
	      if(print == 2){ std::cout << "Term28 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = " << Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T8[ind][key2 * np + p] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term28 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  np = Chan.np[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int p = 0; p < np; ++p){
	    d = Chan.pvec[ind][p];
	    if(d == a){ continue; }
	    key3 = Chan.pp1_map[Chan.ind0][Hash2(d, b, Space.indtot)];
	    for(int hp1 = 0; hp1 < nhp1; ++hp1){
	      k = Chan.hp1vec[Chan.ind0][2*hp1];
	      c = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
	      if(c == d){ continue; }
	      term0 = -1.0 * Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] * Amps1.S1.T1[hp1] * Amps1.D1.T8[ind][key2 * np + p];
	      term += term0;
	      if(print == 2){ std::cout << "Term29 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = -1.0 * " << Ints.S_ME1.V16[Chan.ind0][key3 * nhp1 + hp1] << " * " << Amps1.S1.T1[hp1] << " * " << Amps1.D1.T8[ind][key2 * np + p] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term29 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, b, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V18[ind][key2 * nh + h] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term32 += -1.0 * < " << i << "," << j << " |V| " << k << "," << b << " > * < " << a << " |t| " << k << " > = -1.0 * " << Ints.S_ME1.V18[ind][key2 * nh + h] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term32 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  key2 = Chan.hhp_map[ind][Hash3(i, j, a, Space.indtot)];
	  nh = Chan.nh[ind];
	  nhhp = Chan.nhhp[ind];
	  nhp1 = Chan.nhp1[Chan.ind0];
	  for(int h = 0; h < nh; ++h){
	    k = Chan.hvec[ind][h];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    term0 = Ints.S_ME1.V18[ind][key2 * nh + h] * Amps1.S1.T1[key3];
	    term += term0;
	    if(print == 2){ std::cout << "Term33 += < " << i << "," << j << " |V| " << k << "," << a << " > * < " << b << " |t| " << k << " > = " << Ints.S_ME1.V18[ind][key2 * nh + h] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term33 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  nh = Chan.nh[ind];
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int h = 0; h < nh; ++h){
	      k = Chan.hvec[ind][h];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(k, b, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	      term0 = -0.5 * Ints.S_ME1.V20[chan][pp1 * nhp + key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term36 += -0.5 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = -0.5 * " << Ints.S_ME1.V20[chan][pp1 * nhp + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term36 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  nh = Chan.nh[ind];
	  for(int pp1 = 0; pp1 < npp; ++pp1){
	    c = Chan.ppvec[chan][2*pp1];
	    d = Chan.ppvec[chan][2*pp1 + 1];
	    if(c == d){ continue; }
	    for(int h = 0; h < nh; ++h){
	      k = Chan.hvec[ind][h];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != chan){ continue; }
	      key1 = Chan.hp_map[ind1][Hash2(k, a, Space.indtot)];
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	      term0 = 0.5 * Ints.S_ME1.V20[chan][pp1 * nhp + key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[chan][hh * npp + pp1];
	      term += term0;
	      if(print == 2){ std::cout << "Term37 += 0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << Ints.S_ME1.V20[chan][pp1 * nhp + key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[chan][hh * npp + pp1] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term37 = " << term << std::endl; }
	  tempt += term;

	  //T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[b];
	  key1 = Chan.p_map[ind][b];
	  nhpp = Chan.nhpp[ind];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(c == d || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term40 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << j << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term40 = " << term << std::endl; }
	  tempt += term;

	  //T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  ind = Chan.indvec[a];
	  key1 = Chan.p_map[ind][a];
	  nhpp = Chan.nhpp[ind];
	  for(int hpp = 0; hpp < nhpp; ++hpp){
	    k = Chan.hppvec[ind][3*hpp];
	    c = Chan.hppvec[ind][3*hpp + 1];
	    d = Chan.hppvec[ind][3*hpp + 2];
	    if(c == d || Chan.indvec[k] != Chan.indvec[b] || Chan.indvec[i] != Chan.indvec[d] || Chan.indvec[j] != Chan.indvec[c]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
	    key4 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    term0 = Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	    term += term0;
	    if(print == 2){ std::cout << "Term41 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << c << " |t| " << j << " > = " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term41 = " << term << std::endl; }
	  tempt += term;
	}

	if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.D1.Evec[chan][hh * npp + pp] << ", " << tempt/Amps1.D1.Evec[chan][hh * npp + pp] << std::endl; }

	tempt /= Amps1.D1.Evec[chan][hh * npp + pp];
	Amps2.D1.T1[chan][hh * npp + pp] = tempt;
      }
    }
  }

  if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
    if(print != 0){ std::cout << std::endl << "For Singles" << std::endl; }
    nhp1 = Chan.nhp1[Chan.ind0];
    
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];

      //if(i == 0 && a == 10){ print = 1; }
      //else{ print = 0; }

      if(print != 0){
	std::cout << std::endl;
	std::cout << "!! < " << a << " |t| " << i << " > !!" << std::endl;
	std::cout << "Term0 += 0.0" << std::endl;
      }
      tempt = 0.0;
      
      //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hp2 = 0; hp2 < nhp1; ++hp2){
	k = Chan.hp1vec[Chan.ind0][2*hp2];
	c = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
	term0 = -1.0 * Ints.D_ME1.V3[Chan.ind0][hp1 * nhp1 + hp2] * Amps1.S1.T1[hp2];
	term += term0;
	if(print == 2){ std::cout << "Term1 += -1.0 * < " << k << "," << a << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > = -1.0 * " << Ints.D_ME1.V3[Chan.ind0][hp1 * nhp1 + hp2] << " * " << Amps1.S1.T1[hp2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
      tempt += term;

      //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hp2 = 0; hp2 < nhp1; ++hp2){
	k = Chan.hp1vec[Chan.ind0][2*hp2];
	c = Chan.hp1vec[Chan.ind0][2*hp2 + 1];
	for(int hp3 = 0; hp3 < nhp1; ++hp3){
	  l = Chan.hp1vec[Chan.ind0][2*hp3];
	  d = Chan.hp1vec[Chan.ind0][2*hp3 + 1];
	  if(k == l || c == d || d == a || l == i){ continue; }
	  term0 = Ints.D_ME1.V9[Chan.ind0][hp2 * nhp1 + hp3] * Amps1.S1.T1[hp2] * Amps1.D1.T2[Chan.ind0][hp3 * nhp1 + hp1];
	  term += term0;
	  if(print == 2){ std::cout << "Term6 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > * " << Ints.D_ME1.V9[Chan.ind0][hp2 * nhp1 + hp3] << " * " << Amps1.S1.T1[hp2] << " * " << Amps1.D1.T2[Chan.ind0][hp3 * nhp1 + hp1] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = -0.5 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      ind = Chan.indvec[a];
      key1 = Chan.p_map[ind][a];
      key2 = Chan.h_map[ind][i];
      nh = Chan.nh[ind];
      nhpp = Chan.nhpp[ind];
      np = Chan.np[ind];
      nhhp = Chan.nhhp[ind];
      for(int hpp = 0; hpp < nhpp; ++hpp){
	k = Chan.hppvec[ind][3*hpp];
	c = Chan.hppvec[ind][3*hpp + 1];
	d = Chan.hppvec[ind][3*hpp + 2];
	if(c == d || i == k){ continue; }
	term0 = -0.5 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.D1.T6[ind][hpp * nh + key2];
	term += term0;
	if(print == 2){ std::cout << "Term2 += -0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = -0.5 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.D1.T6[ind][hpp * nh + key2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} = -V11(ka|cd){a,kcd}.t1(c|i){ic}.t1(d|k){kd} (3)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hpp = 0; hpp < nhpp; ++hpp){
	k = Chan.hppvec[ind][3*hpp];
	c = Chan.hppvec[ind][3*hpp + 1];
	d = Chan.hppvec[ind][3*hpp + 2];
	if(c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[d]){ continue; }
	key3 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	key4 = Chan.hp1_map[Chan.ind0][Hash2(k, d, Space.indtot)];
	term0 = -1.0 * Ints.S_ME1.V11[ind][key1 * nhpp + hpp] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	term += term0;
	if(print == 2){ std::cout << "Term3 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << k << " > = -1.0 * " << Ints.S_ME1.V11[ind][key1 * nhpp + hpp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
      tempt += term;

      //t2(a|i){a,i} = -0.5 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int h = 0; h < nh; ++h){
	k = Chan.hvec[ind][h];
	for(int hpp = 0; hpp < nhpp; ++hpp){
	  l = Chan.hppvec[ind][3*hpp];
	  c = Chan.hppvec[ind][3*hpp + 1];
	  d = Chan.hppvec[ind][3*hpp + 2];
	  if(k == l || c == d || i == l || Chan.indvec[k] != Chan.indvec[a]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	  term0 = -0.5 * Ints.D_ME1.V6[ind][h * nhpp + hpp] * Amps1.S1.T1[key1] * Amps1.D1.T6[ind][hpp * nh + key2];
	  term += term0;
	  if(print == 2){ std::cout << "Term7 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << l << " > = -0.5 * " << Ints.D_ME1.V6[ind][h * nhpp + hpp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T6[ind][hpp * nh + key2] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = -0.5 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      ind = Chan.indvec[i];
      key1 = Chan.h_map[ind][i];
      key2 = Chan.p_map[ind][a];
      np = Chan.np[ind];
      nhhp = Chan.nhhp[ind];
      nh = Chan.nh[ind];
      nhpp = Chan.nhpp[ind];
      for(int hhp = 0; hhp < nhhp; ++hhp){
	k = Chan.hhpvec[ind][3*hhp];
	l = Chan.hhpvec[ind][3*hhp + 1];
	c = Chan.hhpvec[ind][3*hhp + 2];
	if(k == l || a == c){ continue; }
	term0 = -0.5 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.D1.T8[ind][hhp * np + key2];
	term += term0;
	if(print == 2){ std::cout << "Term4 += -0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.D1.T8[ind][hhp * np + key2] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} = -V12(kl|ic){i,klc}.t1(a|k){ka}.t1(c|l){lc} (5)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int hhp = 0; hhp < nhhp; ++hhp){
	k = Chan.hhpvec[ind][3*hhp];
	l = Chan.hhpvec[ind][3*hhp + 1];
	c = Chan.hhpvec[ind][3*hhp + 2];
	if(k == l || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[c]){ continue; }
	key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	key4 = Chan.hp1_map[Chan.ind0][Hash2(l, c, Space.indtot)];
	term0 = -1.0 * Ints.S_ME1.V12[ind][key1 * nhhp + hhp] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
	term += term0;
	if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << " |t| " << l << " > = -1.0 * " << Ints.S_ME1.V12[ind][key1 * nhhp + hhp] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
      }
      if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = -0.5 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int p = 0; p < np; ++p){
	c = Chan.pvec[ind][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  k = Chan.hhpvec[ind][3*hhp];
	  l = Chan.hhpvec[ind][3*hhp + 1];
	  d = Chan.hhpvec[ind][3*hhp + 2];
	  if(k == l || c == d || Chan.indvec[i] != Chan.indvec[c]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  term0 = -0.5 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.S1.T1[key1] * Amps1.D1.T8[ind][hhp * np + key2];
	  term += term0;
	  if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T8[ind][hhp * np + key2] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
      tempt += term;

      //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} = -t3(i|c){i,c}.V8(kl|cd){c,kld}.t1(a|k){ka}.t1(d|l){ld} (9)
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int p = 0; p < np; ++p){
	c = Chan.pvec[ind][p];
	for(int hhp = 0; hhp < nhhp; ++hhp){
	  k = Chan.hhpvec[ind][3*hhp];
	  l = Chan.hhpvec[ind][3*hhp + 1];
	  d = Chan.hhpvec[ind][3*hhp + 2];
	  if(k == l || c == d || Chan.indvec[i] != Chan.indvec[c] || Chan.indvec[k] != Chan.indvec[a] || Chan.indvec[l] != Chan.indvec[d]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	  key3 = Chan.hp1_map[Chan.ind0][Hash2(l, d, Space.indtot)];
	  term0 = -1.0 * Ints.D_ME1.V8[ind][p * nhhp + hhp] * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  term += term0;
	  if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << l << " > = -1.0 * " << Ints.D_ME1.V8[ind][p * nhhp + hhp] << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
      tempt += term;

      if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.S1.Evec[hp1] << ", " << tempt/Amps1.S1.Evec[hp1] << std::endl; }
      tempt /= Amps1.S1.Evec[hp1];
      Amps2.S1.T1[hp1] = tempt;
    }
  }

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }
	tempt = Amps2.D1.T1[chan][hh * npp + pp];
	norm2 += tempt * tempt;
	error2 += (tempt - Amps1.D1.T1[chan][hh * npp + pp]) * (tempt - Amps1.D1.T1[chan][hh * npp + pp]);
	Amps1.D1.set_T(chan, hh * npp + pp, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){  
    nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps2.S1.T1[hp1];
      norm2 += tempt * tempt;
      error2 += (tempt - Amps1.S1.T1[hp1]) * (tempt - Amps1.S1.T1[hp1]);
      Amps1.S1.set_T(hp1, tempt);
    }
  }
  error = std::sqrt(error2/norm2);
}

void Doubles_Step_explicit2(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Amplitudes &Amps1, Amplitudes &Amps2, const HF_Channels &HF_Chan, const HF_Matrix_Elements &HF_ME, double &error)
{
  State tb;
  int ind1, ind2, ind3, ind4, Vind1, Vind2, key1, key2, key3, key4, Vkey1, Vkey2;
  int nhh, npp, nhp1, i, j, a, b;
  double tempt, term0, term, TBME;
  double norm2 = 0.0;
  double error2 = 0.0;
  int print = 0;
  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }

	//if(i == 2 && j == 3 && a == 6 && b == 11){ print = 1; }
	//else{ print = 0; }
	//if(i < j && a < b){ print = true; }
	//else{ print = false; }

	plus(tb, Space.qnums[a], Space.qnums[b]);
	Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	plus(tb, Space.qnums[i], Space.qnums[j]);
	Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	if(Vind1 != Vind2){ continue; }
	Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	tempt = TBME;
	if(print != 0){
	  std::cout << std::endl;
	  std::cout << "!! < " << a << "," << b << " |t| " << i << "," << j << " > !!" << std::endl;
	  std::cout << "Term0 = < " << a << "," << b << " |V| " << i << "," << j << " > = " << TBME << std::endl;
	}

	term = 0.0;
 	if(print == 2){ std::cout << std::endl; }
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(c == d){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, j, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term1 += 0.5 * < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(k, l, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, b, Space.indtot)];
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(i, j, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(k, l, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += 0.5 * < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind3 != ind4){ continue; }
	    key3 = Chan.hh_map[ind3][Hash2(k, l, Space.indtot)];
	    key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[i], Space.qnums[j]);
		ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2){ continue; }
		key1 = Chan.hh_map[ind1][Hash2(i, j, Space.indtot)];
		key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.25 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term7 += 0.25 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.25 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term3 += < " << k << "," << b << " |V| " << c << "," << j << " > * < " << a << "," << c << " |t| " << i << "," << k << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term4 += -1.0 * < " << k << "," << b << " |V| " << c << "," << i << " > * < " << a << "," << c << " |t| " << j << "," << k << " > = -1.0 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == b){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[b], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << a << " |V| " << c << "," << j << " > * < " << b << "," << c << " |t| " << i << "," << k << " > = -1.0 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == b){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[b], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term6 += < " << k << "," << a << " |V| " << c << "," << i << " > * < " << b << "," << c << " |t| " << j << "," << k << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
	if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  plus(tb, Space.qnums[i], Space.qnums[k]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == j){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	      key2 = Chan.pp_map[ind1][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[b], Space.qnums[d]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key3 = Chan.hh_map[ind3][Hash2(j, l, Space.indtot)];
		key4 = Chan.pp_map[ind3][Hash2(b, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << i << "," << k << " > * < " << b << "," << d << " |t| " << j << "," << l << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  plus(tb, Space.qnums[j], Space.qnums[k]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == i){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[l]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	      key2 = Chan.pp_map[ind1][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[b], Space.qnums[d]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key3 = Chan.hh_map[ind3][Hash2(i, l, Space.indtot)];
		key4 = Chan.pp_map[ind3][Hash2(b, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << j << "," << k << " > * < " << b << "," << d << " |t| " << i << "," << l << " > = " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[l], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(l, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << "," << c << " |t| " << l << "," << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  for(int l = 0; l < Space.indhol; ++l){
	    if(k == l){ continue; }
	    plus(tb, Space.qnums[l], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(l, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(c == b){ continue; }
	      plus(tb, Space.qnums[b], Space.qnums[c]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind1 != ind2){ continue; }
	      key2 = Chan.pp_map[ind2][Hash2(b, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind3 != ind4){ continue; }
		key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term9 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << b << "," << c << " |t| " << l << "," << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == i){ continue; }
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == j){ continue; }
	    plus(tb, Space.qnums[i], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[l], Space.qnums[j]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[c]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[a], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2 || ind3 != ind4){ continue; }
		key2 = Chan.pp_map[ind2][Hash2(d, c, Space.indtot)];
		key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term10 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << "," << c << " |t| " << i << "," << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	tempt += term;

	term = 0.0;
	if(print == 2){ std::cout << std::endl; }
	for(int k = 0; k < Space.indhol; ++k){
	  if(k == j){ continue; }
	  for(int l = 0; l < Space.indhol; ++l){
	    if(l == i){ continue; }
	    plus(tb, Space.qnums[j], Space.qnums[k]);
	    ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[l], Space.qnums[i]);
	    ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    key1 = Chan.hh_map[ind1][Hash2(j, k, Space.indtot)];
	    key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(c == d){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[c]);
		ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[a], Space.qnums[b]);
		ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind1 != ind2 || ind3 != ind4){ continue; }
		key2 = Chan.pp_map[ind2][Hash2(d, c, Space.indtot)];
		key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		term += term0;
		if(print == 2){ std::cout << "Term11 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << d << "," << c << " |t| " << j << "," << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	}
	if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	tempt += term;

	
	if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
	  if(print != 0){ std::cout << std::endl << "For Singles" << std::endl << std::endl; }
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term1 += < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
	  tempt += term;
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[i], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term2 += < " << k << "," << l << " |V| " << i << "," << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(c == d){ continue; }
		  plus(tb, Space.qnums[i], Space.qnums[j]);
		  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(c, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term3 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(k, l, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term4 += 0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key3 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key4 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];	    
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3] * Amps1.S1.T1[key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term5 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " * " << Amps1.S1.T1[key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term6 += -1.0 * < " << k << "," << b << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[b]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[i]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term7 += < " << k << "," << b << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[j]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term8 += < " << k << "," << a << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[a]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[i]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];	    
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	      term += term0;
	      if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << a << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	    }
	  }
	  if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == b){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[b]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term10 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term10 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == b){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[b]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term11 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << d << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term11 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == a){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[a]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term12 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term12 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
		key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		for(int d = Space.indhol; d < Space.indtot; ++d){
		  if(d == a){ continue; }
		  plus(tb, Space.qnums[d], Space.qnums[a]);
		  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(ind3 != ind4){ continue; }
		  key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
		  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
		  term += term0;
		  if(print == 2){ std::cout << "Term13 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << b << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
		}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term13 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == j){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[j]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[a], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term14 += < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << b << " |t| " << k << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term14 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == i){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[i]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == b){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[a], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term15 += -1.0 * < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << b << " |t| " << k << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term15 = " << term << std::endl; }
	  tempt += term;
	  
	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == j){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[j]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[b], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(b, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term16 += -1.0 * < " << b << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << "," << a << " |t| " << k << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term16 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(k == i){ continue; }
	    plus(tb, Space.qnums[k], Space.qnums[i]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(d == a){ continue; }
		plus(tb, Space.qnums[d], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(k, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[b], Space.qnums[k]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(b, k, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term17 += < " << b << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << j << " > * < " << d << "," << a << " |t| " << k << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term17 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == b){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[i], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term18 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term18 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == b){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[b]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[j], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(j, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term19 += < " << k << "," << l << " |V| " << j << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term19 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == a){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[i], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term20 += < " << k << "," << l << " |V| " << i << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term20 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
		if(c == a){ continue; }
		plus(tb, Space.qnums[c], Space.qnums[a]);
		ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(ind2 != ind3){ continue; }
		key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
		key3 = Chan.pp_map[ind3][Hash2(c, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[j], Space.qnums[c]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(j, c, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
		term += term0;
		if(print == 2){ std::cout << "Term21 += -1.0 * < " << k << "," << l << " |V| " << j << "," << c << " > * < " << b << " |t| " << k << " > * < " << c << "," << a << " |t| " << l << "," << i << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term21 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(l, j, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	   	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	   	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	   	for(int d = Space.indhol; d < Space.indtot; ++d){
	   	  if(Chan.indvec[d] != Chan.indvec[i]){ continue; }
	   	  key2 = Chan.hp1_map[Chan.ind0][Hash2(i, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	   	  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	   	  term += term0;
	   	  if(print == 2){ std::cout << "Term22 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << " |t| " << i << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	   	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term22 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[a], Space.qnums[b]);
	      ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind3 != ind4){ continue; }
	      key3 = Chan.hh_map[ind3][Hash2(l, i, Space.indtot)];
	      key4 = Chan.pp_map[ind4][Hash2(a, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	   	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	   	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	   	for(int d = Space.indhol; d < Space.indtot; ++d){
	   	  if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	   	  key2 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	   	  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	   	  term += term0;
	   	  if(print == 2){ std::cout << "Term23 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << " |t| " << j << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	   	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term23 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[a]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, a, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	for(int d = Space.indhol; d < Space.indtot; ++d){
	  	  if(d == b){ continue; }
	  	  plus(tb, Space.qnums[i], Space.qnums[j]);
	  	  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  plus(tb, Space.qnums[d], Space.qnums[b]);
	  	  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  if(ind3 != ind4){ continue; }
	  	  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	  	  key4 = Chan.pp_map[ind4][Hash2(d, b, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	  term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	  	  term += term0;
	  	  if(print == 2){ std::cout << "Term24 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << a << " |t| " << l << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	  	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term24 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key2 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	for(int d = Space.indhol; d < Space.indtot; ++d){
	  	  if(d == a){ continue; }
	  	  plus(tb, Space.qnums[i], Space.qnums[j]);
	  	  ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  plus(tb, Space.qnums[d], Space.qnums[a]);
	  	  ind4 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	  if(ind3 != ind4){ continue; }
	  	  key3 = Chan.hh_map[ind3][Hash2(i, j, Space.indtot)];
	  	  key4 = Chan.pp_map[ind4][Hash2(d, a, Space.indtot)];
		  plus(tb, Space.qnums[k], Space.qnums[l]);
		  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  plus(tb, Space.qnums[c], Space.qnums[d]);
		  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		  if(Vind1 != Vind2){ continue; }
		  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	  term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4];
	  	  term += term0;
	  	  if(print == 2){ std::cout << "Term25 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << b << " |t| " << l << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.D1.T1[ind3][key3 * Chan.npp[ind3] + key4] << " = " << term0 << std::endl; }
	  	}
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term25 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == j){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[j]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(l, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term26 += -1.0 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term26 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(l == i){ continue; }
	      plus(tb, Space.qnums[l], Space.qnums[i]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term27 += < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << k << " > * < " << a << "," << b << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term27 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(d == b){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[d], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(d, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term28 += < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << b << " |t| " << i << "," << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term28 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[k]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(d == a){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[d], Space.qnums[a]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term29 += -1.0 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << i << "," << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term29 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term30 += < " << a << "," << b << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term30 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[i]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term31 += -1.0 * < " << a << "," << b << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term31 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[b]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term32 += -1.0 * < " << k << "," << b << " |V| " << i << "," << j << " > * < " << a << " |t| " << k << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term32 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[a]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[j]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, j, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1];
	    term += term0;
	    if(print == 2){ std::cout << "Term33 += < " << k << "," << a << " |V| " << i << "," << j << " > * < " << b << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " = " << term0 << std::endl; }
	  }
	  if(print != 0){ std::cout << "Term33 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term34 += 0.5 * < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term34 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    for(int l = 0; l < Space.indhol; ++l){
	      if(k == l){ continue; }
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
	  	plus(tb, Space.qnums[a], Space.qnums[b]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(a, b, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term35 += -0.5 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term35 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(c == d){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[c], Space.qnums[d]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term36 += -0.5 * < " << k << "," << b << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term36 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(c == d){ continue; }
	  	plus(tb, Space.qnums[i], Space.qnums[j]);
	  	ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	plus(tb, Space.qnums[c], Space.qnums[d]);
	  	ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  	if(ind2 != ind3){ continue; }
	  	key2 = Chan.hh_map[ind2][Hash2(i, j, Space.indtot)];
	  	key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = 0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term37 += 0.5 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << b << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << j << " > = 0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term37 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[j]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, j, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term38 += < " << k << "," << l << " |V| " << c << "," << j << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term38 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int l = 0; l < Space.indhol; ++l){
	      if(Chan.indvec[l] != Chan.indvec[b]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, b, Space.indtot)];
	      for(int c = Space.indhol; c < Space.indtot; ++c){
	  	if(Chan.indvec[c] != Chan.indvec[j]){ continue; }
	  	key1 = Chan.hp1_map[Chan.ind0][Hash2(j, c, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[l]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[i]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, i, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term39 += -1.0 * < " << k << "," << l << " |V| " << c << "," << i << " > * < " << c << " |t| " << j << " > * < " << a << " |t| " << k << " > * < " << b << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term39 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
	  	if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
	  	key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[b]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, b, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  	term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	  	term += term0;
	  	if(print == 2){ std::cout << "Term40 += < " << k << "," << b << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << j << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term40 = " << term << std::endl; }
	  tempt += term;

	  term = 0.0;
	  if(print == 2){ std::cout << std::endl; }
	  for(int k = 0; k < Space.indhol; ++k){
	    if(Chan.indvec[k] != Chan.indvec[b]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, b, Space.indtot)];
	    for(int c = Space.indhol; c < Space.indtot; ++c){
	      if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	      key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	      for(int d = Space.indhol; d < Space.indtot; ++d){
		if(Chan.indvec[d] != Chan.indvec[j]){ continue; }
		key3 = Chan.hp1_map[Chan.ind0][Hash2(j, d, Space.indtot)];
		plus(tb, Space.qnums[k], Space.qnums[a]);
		Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		plus(tb, Space.qnums[c], Space.qnums[d]);
		Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
		if(Vind1 != Vind2){ continue; }
		Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, a, Space.indtot)];
		Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
		TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
		term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
		term += term0;
		if(print == 2){ std::cout << "Term41 += -1.0 * < " << k << "," << a << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << b << " |t| " << k << " > * < " << d << " |t| " << j << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	      }
	    }
	  }
	  if(print != 0){ std::cout << "Term41 = " << term << std::endl; }
	  tempt += term;
	}

	if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.D1.Evec[chan][hh * npp + pp] << ", " << tempt/Amps1.D1.Evec[chan][hh * npp + pp] << std::endl; }
	
	tempt /= Amps1.D1.Evec[chan][hh * npp + pp];
	Amps2.D1.T1[chan][hh * npp + pp] = tempt;
      }
    }
  }

  if(Parameters.approx == "singles" && Chan.nhp1[Chan.ind0] != 0){
    if(print != 0){ std::cout << std::endl << "For Singles" << std::endl; }
    nhp1 = Chan.nhp1[Chan.ind0];
    
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];

      //if(i == 0 && a == 10){ print = 1; }
      //else{ print = 0; }

      if(print != 0){ 
	std::cout << std::endl;
	std::cout << "!! < " << a << " |t| " << i << " > !!" << std::endl;
	std::cout << "Term0 += 0.0" << std::endl;
      }
      tempt = 0.0;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	  plus(tb, Space.qnums[a], Space.qnums[k]);
	  Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  plus(tb, Space.qnums[i], Space.qnums[c]);
	  Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  if(Vind1 != Vind2){ continue; }
	  Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	  Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	  TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	  term0 = TBME * Amps1.S1.T1[key1];
	  term += term0;
	  if(print == 2){ std::cout << "Term1 += < " << a << "," << k << " |V| " << i << "," << c << " > * < " << c << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << " = " << term0 << std::endl; }
	}
      }
      if(print != 0){ std::cout << "Term1 = " << term << std::endl; }
      tempt += term;      

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(k == i){ continue; }
	plus(tb, Space.qnums[i], Space.qnums[k]);
	ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(c == d){ continue; }
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(i, k, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(c, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[k]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = 0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term2 += 0.5 * < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << k << " > = 0.5 * " << TBME << " * " << Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term2 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int c = Space.indhol; c < Space.indtot; ++c){
	  if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	  key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];	  
	  for(int d = Space.indhol; d < Space.indtot; ++d){
	    if(Chan.indvec[d] != Chan.indvec[k]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(k, d, Space.indtot)];
	    plus(tb, Space.qnums[a], Space.qnums[k]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[c], Space.qnums[d]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(a, k, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term3 += < " << a << "," << k << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << d << " |t| " << k << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term3 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(k == l){ continue; }
	  plus(tb, Space.qnums[k], Space.qnums[l]);
	  ind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(c == a){ continue; }
	    plus(tb, Space.qnums[a], Space.qnums[c]);
	    ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(ind1 != ind2){ continue; }
	    key1 = Chan.hh_map[ind1][Hash2(k, l, Space.indtot)];
	    key2 = Chan.pp_map[ind2][Hash2(a, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[c]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -0.5 * TBME * Amps1.D1.T1[ind1][key1 * Chan.npp[ind1] + key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term4 += -0.5 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << "," << c << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.D1.T1[ind2][key1 * Chan.npp[ind1] + key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term4 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[l]){ continue; }
	    key2 = Chan.hp1_map[Chan.ind0][Hash2(l, c, Space.indtot)];
	    plus(tb, Space.qnums[k], Space.qnums[l]);
	    Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    plus(tb, Space.qnums[i], Space.qnums[c]);
	    Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	    if(Vind1 != Vind2){ continue; }
	    Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	    Vkey2 = HF_Chan.tb_map[Vind2][Hash2(i, c, Space.indtot)];
	    TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	    term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2];
	    term += term0;
	    if(print == 2){ std::cout << "Term5 += -1.0 * < " << k << "," << l << " |V| " << i << "," << c << " > * < " << a << " |t| " << k << " > * < " << c << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " = " << term0 << std::endl; }
	  }
	}
      }
      if(print != 0){ std::cout << "Term5 = " << term << std::endl; }
      tempt += term;
      
      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(l == i){ continue; }
	  plus(tb, Space.qnums[l], Space.qnums[i]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[k] != Chan.indvec[c]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(k, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(d == a){ continue; }
	      plus(tb, Space.qnums[d], Space.qnums[a]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(l, i, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(d, a, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term6 += < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << k << " > * < " << d << "," << a << " |t| " << l << "," << i << " > = " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term6 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key1 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  if(l == i){ continue; }
	  plus(tb, Space.qnums[i], Space.qnums[l]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(c == d){ continue; }
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(i, l, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(c, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term7 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << a << " |t| " << k << " > * < " << c << "," << d << " |t| " << i << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term7 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	for(int l = 0; l < Space.indhol; ++l){
	  if(k == l){ continue; }
	  plus(tb, Space.qnums[k], Space.qnums[l]);
	  ind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(d == a){ continue; }
	      plus(tb, Space.qnums[a], Space.qnums[d]);
	      ind3 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(ind2 != ind3){ continue; }
	      key2 = Chan.hh_map[ind2][Hash2(k, l, Space.indtot)];
	      key3 = Chan.pp_map[ind3][Hash2(a, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -0.5 * TBME * Amps1.S1.T1[key1] * Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term8 += -0.5 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << "," << d << " |t| " << k << "," << l << " > = -0.5 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.D1.T1[ind2][key2 * Chan.npp[ind2] + key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term8 = " << term << std::endl; }
      tempt += term;

      term = 0.0;
      if(print == 2){ std::cout << std::endl; }
      for(int k = 0; k < Space.indhol; ++k){
	if(Chan.indvec[k] != Chan.indvec[a]){ continue; }
	key2 = Chan.hp1_map[Chan.ind0][Hash2(k, a, Space.indtot)];
	for(int l = 0; l < Space.indhol; ++l){
	  for(int c = Space.indhol; c < Space.indtot; ++c){
	    if(Chan.indvec[c] != Chan.indvec[i]){ continue; }
	    key1 = Chan.hp1_map[Chan.ind0][Hash2(i, c, Space.indtot)];
	    for(int d = Space.indhol; d < Space.indtot; ++d){
	      if(Chan.indvec[d] != Chan.indvec[l]){ continue; }
	      key3 = Chan.hp1_map[Chan.ind0][Hash2(l, d, Space.indtot)];
	      plus(tb, Space.qnums[k], Space.qnums[l]);
	      Vind1 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      plus(tb, Space.qnums[c], Space.qnums[d]);
	      Vind2 = ChanInd_2b_dir(Parameters.basis, Space, tb);
	      if(Vind1 != Vind2){ continue; }
	      Vkey1 = HF_Chan.tb_map[Vind1][Hash2(k, l, Space.indtot)];
	      Vkey2 = HF_Chan.tb_map[Vind2][Hash2(c, d, Space.indtot)];
	      TBME = HF_ME.V[Vind1][Vkey1 * HF_Chan.ntb[Vind1] + Vkey2];
	      term0 = -1.0 * TBME * Amps1.S1.T1[key1] * Amps1.S1.T1[key2] * Amps1.S1.T1[key3];
	      term += term0;
	      if(print == 2){ std::cout << "Term9 += -1.0 * < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << " |t| " << i << " > * < " << a << " |t| " << k << " > * < " << d << " |t| " << l << " > = -1.0 * " << TBME << " * " << Amps1.S1.T1[key1] << " * " << Amps1.S1.T1[key2] << " * " << Amps1.S1.T1[key3] << " = " << term0 << std::endl; }
	    }
	  }
	}
      }
      if(print != 0){ std::cout << "Term9 = " << term << std::endl; }
      tempt += term;

      if(print != 0){ std::cout << std::endl << "tempt, E, tempt/E = " << tempt << ", " << Amps1.S1.Evec[hp1] << ", " << tempt/Amps1.S1.Evec[hp1] << std::endl; }
      tempt /= Amps1.S1.Evec[hp1];
      Amps2.S1.T1[hp1] = tempt;
    }
  }

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh = Chan.nhh[chan];
    npp = Chan.npp[chan];
    if(nhh * npp == 0){ continue; }
    for(int pp = 0; pp < npp; ++pp){
      a = Chan.ppvec[chan][2*pp];
      b = Chan.ppvec[chan][2*pp + 1];
      if(a == b){ continue; }
      for(int hh = 0; hh < nhh; ++hh){
	i = Chan.hhvec[chan][2*hh];
	j = Chan.hhvec[chan][2*hh + 1];
	if(i == j){ continue; }
	tempt = Amps2.D1.T1[chan][hh * npp + pp];
	norm2 += tempt * tempt;
	error2 += (tempt - Amps1.D1.T1[chan][hh * npp + pp]) * (tempt - Amps1.D1.T1[chan][hh * npp + pp]);
	Amps1.D1.set_T(chan, hh * npp + pp, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){  
    nhp1 = Chan.nhp1[Chan.ind0];
    for(int hp1 = 0; hp1 < nhp1; ++hp1){
      i = Chan.hp1vec[Chan.ind0][2*hp1];
      a = Chan.hp1vec[Chan.ind0][2*hp1 + 1];
      tempt = Amps2.S1.T1[hp1];
      norm2 += tempt * tempt;
      error2 += (tempt - Amps1.S1.T1[hp1]) * (tempt - Amps1.S1.T1[hp1]);
      Amps1.S1.set_T(hp1, tempt);
    }
  }
  error = std::sqrt(error2/norm2);
}

void CC_compare_JM(const Input_Parameters &Parameters, const Model_Space &Space, const Channels &Chan, Interactions &Ints, Amplitudes &Amps, std::string &inputfile)
{
  Input_Parameters Parameters2;
  Model_Space Space2;
  Channels Chan2;
  Amplitudes Amps2;
  Interactions Ints2;
  HF_Channels HF_Chan2;
  HF_Matrix_Elements HF_ME2;
  Single_Particle_States States2;
  int JM2 = 0;
  State tb1, tb2, tb3, tb4;
  int nh_j, np_j, nhh_j, npp_j, nhpp_j, nhhp_j, nhp1_j, nhp2_j, nhp0_j, chan_j, chan0_j;
  int nh_m, np_m, nhh_m, npp_m, nhpp_m, nhhp_m, nhp1_m, nhp2_m, nhp0_m, chan_m, chan0_m;
  int h_j, p_j, hh_j, pp_j, hpp_j, hhp_j, hp1_j, hp2_j;
  int h_m, p_m, hh_m, pp_m, hpp_m, hhp_m, hp1_m, hp2_m;
  double TJ;
  int ind, minj, jsize;

  int a2, b2, i2, j2;
  int k, l, c, d, k2, l2, c2, d2;

  double ja2, ma2, jb2, mb2, ji2, mi2, jj2, mj2, m12, m22;
  double tempt, tempen, term0;
  double *jvec, *tvec;
  double energy1 = 0.0;
  double energy2 = 0.0;
  int print = 0;
  int count = 0;
  int iterations = 10;

  Get_Input_Parameters(inputfile, Parameters2);
  Parameters2.basis = "finite_JM";

  Build_Model_Space(Parameters2, Space2);
  Print_Parameters(Parameters2, Space2);
  HF_Chan2 = HF_Channels(Parameters2, Space2);
  States2 = Single_Particle_States(Parameters2, Space2, HF_Chan2);
  Read_Matrix_Elements_J(Parameters2, Space2, HF_Chan2, HF_ME2);
  Hartree_Fock_States(Parameters2, Space2, HF_Chan2, States2, HF_ME2);
  Convert_To_HF_Matrix_Elements(HF_Chan2, States2, HF_ME2);
  States2.delete_struct(HF_Chan2);
  
  Build_Model_Space_J2(Parameters2, Space2);
  Parameters2.basis = "finite_M";
  JM2 = 1;
   
  Chan2 = Channels(Parameters2, Space2);
  Amps2 = Amplitudes(Parameters2, Space2, Chan2);
  Ints2 = Interactions(Parameters2, Chan2);
  Get_Matrix_Elements_JM(Parameters2, HF_Chan2, HF_ME2, Space2, Chan2, Ints2);
  HF_ME2.delete_struct(HF_Chan2);
  HF_Chan2.delete_struct();

  if(Parameters.approx == "singles"){
    chan0_j = Chan.ind0;
    nhp0_j = Chan.nhp1[chan0_j];
  }
  if(Parameters2.approx == "singles"){
    chan0_m = Chan2.ind0;
    nhp0_m = Chan2.nhp1[chan0_m];
  }
  
  Amplitudes tempAmps = Amplitudes(Parameters, Space, Chan);
  //Amplitudes tempAmpsJ = Amplitudes(Parameters, Space, Chan);
  Amplitudes tempAmps2 = Amplitudes(Parameters2, Space2, Chan2);
  //Amplitudes tempAmps2J = Amplitudes(Parameters2, Space2, Chan2);
  //Interactions tempInts2J = Interactions(Parameters2, Chan2);

  Amps.zero(Parameters, Chan);
  Amps2.zero(Parameters2, Chan2);

  for(int chan = 0; chan < Chan.size1; ++chan){
    nhh_j = Chan.nhh[chan];
    npp_j = Chan.npp[chan];
    if(nhh_j * npp_j == 0){ continue; }
    for(int pp = 0; pp < npp_j; ++pp){
      a2 = Chan.ppvec[chan][2*pp];
      b2 = Chan.ppvec[chan][2*pp + 1];
      for(int hh = 0; hh < nhh_j; ++hh){
	i2 = Chan.hhvec[chan][2*hh];
	j2 = Chan.hhvec[chan][2*hh + 1];
	tempen = Space.qnums[i2].energy + Space.qnums[j2].energy - Space.qnums[a2].energy - Space.qnums[b2].energy;
	tempt = Ints.D_ME1.V4[chan][pp * nhh_j + hh] / tempen;
	//std::cout << "TJ: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > ^ " << 0.5*Chan.qnums1[chan].j << " = " << tempt << std::endl;
	ind = hh * npp_j + pp;
	Amps.D1.set_TJ(Space, Chan, chan, ind, i2, j2, a2, b2, tempt);
      }
    }
  }
  if(Parameters.approx == "singles"){
    for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
      i2 = Chan.hp1vec[chan0_j][2*hp1];
      a2 = Chan.hp1vec[chan0_j][2*hp1 + 1];
      tempen = Space.qnums[i2].energy - Space.qnums[a2].energy;
      tempt = 0.0;
      //std::cout << "tj: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
      Amps.S1.set_TJ(Space, hp1, i2, a2, tempt);
    }
    Amps.S1.set_T_2J(Space, Chan, Ints);
    //Amps.D1.set_T_2J(Chan, Ints);
  }
  for(int chan = 0; chan < Chan2.size1; ++chan){
    nhh_m = Chan2.nhh[chan];
    npp_m = Chan2.npp[chan];
    if(nhh_m * npp_m == 0){ continue; }
    for(int pp = 0; pp < npp_m; ++pp){
      a2 = Chan2.ppvec[chan][2*pp];
      b2 = Chan2.ppvec[chan][2*pp + 1];
      if(a2 == b2){ continue; }
      for(int hh = 0; hh < nhh_m; ++hh){
	i2 = Chan2.hhvec[chan][2*hh];
	j2 = Chan2.hhvec[chan][2*hh + 1];
	if(i2 == j2){ continue; }
	tempen = Space2.qnums[i2].energy + Space2.qnums[j2].energy - Space2.qnums[a2].energy - Space2.qnums[b2].energy;
	tempt = Ints2.D_ME1.V4[chan][pp * nhh_m + hh] / tempen;
	//std::cout << "TM: < " << a << "," << b << " |t| " << i << "," << j << " > = " << tempt << std::endl;
	ind = hh * npp_m + pp;
	Amps2.D1.set_T(chan, ind, tempt);
      }
    }
  }
  if(Parameters2.approx == "singles"){
    for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
      i2 = Chan2.hp1vec[chan0_m][2*hp1];
      a2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
      tempen = Space2.qnums[i2].energy - Space2.qnums[a2].energy;
      tempt = 0.0;
      //std::cout << "tm: " << i << " |t| " << a << " = " << tempt << std::endl;
      Amps2.S1.set_T(hp1, tempt);
    }
    Amps2.S1.set_T_2(Chan2, Ints2);
    //Amps2.D1.set_T_2(Chan2, Ints2);
  }
  energy1 = Amps.get_energy(Parameters, Chan, Ints);
  energy2 = Amps2.get_energy(Parameters2, Chan2, Ints2);
  std::cout << "!?!? " << energy1 << " " << energy2 << std::endl;
      
  while(count < iterations){
    ++count;
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	      tempt += Ints.D_ME1.V4[chan_j][pp_j * nhh_j + hh_j];
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term0 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		    tempt += Ints2.D_ME1.V4[chan_m][pp_m * nhh_m + hh_m];
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term0 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "?Term0 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/2 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
	      for(int pp = 0; pp < npp_j; ++pp){
		c = Chan.ppvec[chan_j][2*pp];
		d = Chan.ppvec[chan_j][2*pp + 1];
		term0 = 0.5 * Ints.D_ME1.V1[chan_j][pp * npp_j + pp_j] * Amps.D1.T1[chan_j][hh_j * npp_j + pp];
		if(print == 2){ std::cout << "Term1 += 1/2 < " << a << "," << b << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = 1/2 * " << Ints.D_ME1.V1[chan_j][pp * npp_j + pp_j] << " * " << Amps.D1.T1[chan_j][hh_j * npp_j + pp] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term1 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/2 * T1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
		    for(int pp = 0; pp < npp_m; ++pp){
		      c2 = Chan2.ppvec[chan_m][2*pp];
		      d2 = Chan2.ppvec[chan_m][2*pp + 1];
		      if(c2 == d2){ continue; }
		      term0 = 0.5 * Ints2.D_ME1.V1[chan_m][pp * npp_m + pp_m] * Amps2.D1.T1[chan_m][hh_m * npp_m + pp];
		      if(print == 2){ std::cout << "!Term1 += 1/2 < " << a2 << "," << b2 << " |V| " << c2 << "," << d2 << " > * < " << c2 << "," << d2 << " |t| " << i2 << "," << j2 << " > = 1/2 * " << Ints2.D_ME1.V1[chan_m][pp * npp_m + pp_m] << " * " << Amps2.D1.T1[chan_m][hh_m * npp_m + pp] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term1 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term1 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
    
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/2 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
	      for(int hh = 0; hh < nhh_j; ++hh){
		k = Chan.hhvec[chan_j][2*hh];
		l = Chan.hhvec[chan_j][2*hh + 1];
		term0 = 0.5 * Ints.D_ME1.V2[chan_j][hh_j * nhh_j + hh] * Amps.D1.T1[chan_j][hh * npp_j + pp_j];
		if(print == 2){ std::cout << "Term2 += 1/2 < " << i << "," << j << " |V| " << k << "," << l << " > * < " << a << "," << b << " |t| " << k << "," << l << " > ^ " << 0.5*tb1.j << " = 1/2 * " << Ints.D_ME1.V2[chan_j][hh_j * nhh_j + hh] << " * " << Amps.D1.T1[chan_j][hh * npp_j + pp_j] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term2 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/2 * V2(ij|kl){ij,kl}.T1(ab|kl){kl,ab} (2)
		    for(int hh = 0; hh < nhh_m; ++hh){
		      k2 = Chan2.hhvec[chan_m][2*hh];
		      l2 = Chan2.hhvec[chan_m][2*hh + 1];
		      if(k2 == l2){ continue; }
		      term0 = 0.5 * Ints2.D_ME1.V2[chan_m][hh_m * nhh_m + hh] * Amps2.D1.T1[chan_m][hh * npp_m + pp_m];
		      if(print == 2){ std::cout << "!Term2 += 1/2 < " << i2 << "," << j2 << " |V| " << k2 << "," << l2 << " > * < " << a2 << "," << b2 << " |t| " << k2 << "," << l2 << " > = 1/2 * " << Ints2.D_ME1.V2[chan_m][hh_m * nhh_m + hh] << " * " << Amps2.D1.T1[chan_m][hh * npp_m + pp_m] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term2 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term2 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
    
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_dir(Parameters.basis, Space, tb1);
	      hh_j = Chan.hh_map[chan_j][Hash2(i, j, Space.indtot)];
	      pp_j = Chan.pp_map[chan_j][Hash2(a, b, Space.indtot)];
	      nhh_j = Chan.nhh[chan_j];
	      npp_j = Chan.npp[chan_j];
	    
	      //T1(ab|ij){ij,ab} = 1/4 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
	      for(int pp = 0; pp < npp_j; ++pp){
		c = Chan.ppvec[chan_j][2*pp];
		d = Chan.ppvec[chan_j][2*pp + 1];
		for(int hh = 0; hh < nhh_j; ++hh){
		  k = Chan.hhvec[chan_j][2*hh];
		  l = Chan.hhvec[chan_j][2*hh + 1];
		  term0 = 0.25 * Ints.D_ME1.V4[chan_j][pp * nhh_j + hh] * Amps.D1.T1[chan_j][hh_j * npp_j + pp] * Amps.D1.T1[chan_j][hh * npp_j + pp_j];
		  if(print == 2){ std::cout << "Term7 += 1/4 < " << k << "," << l << " |V| " << c << "," << d << " > * < " << c << "," << d << " |t| " << i << "," << j << " > * < " << a << "," << b << " |t| " << k << "," << l << " > = 1/4 * " << Ints.D_ME1.V4[chan_j][pp * nhh_j + hh] << " * " << Amps.D1.T1[chan_j][hh_j * npp_j + pp] << " * " << Amps.D1.T1[chan_j][hh * npp_j + pp_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T1[chan_j][hh_j * npp_j + pp_j] += tempt;
	      if(print != 0){
		std::cout << "Term7 = < " << a << "," << b << " |t| " << i << "," << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j = tb1.j - 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	  
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_dir(Parameters2.basis, Space2, tb3);
		    hh_m = Chan2.hh_map[chan_m][Hash2(i2, j2, Space2.indtot)];
		    pp_m = Chan2.pp_map[chan_m][Hash2(a2, b2, Space2.indtot)];
		    nhh_m = Chan2.nhh[chan_m];
		    npp_m = Chan2.npp[chan_m];
		  
		    //T1(ab|ij){ij,ab} = 1/4 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(kl|ab){kl,ab} (7)
		    for(int pp = 0; pp < npp_m; ++pp){
		      c2 = Chan2.ppvec[chan_m][2*pp];
		      d2 = Chan2.ppvec[chan_m][2*pp + 1];
		      if(c2 == d2){ continue; }
		      for(int hh = 0; hh < nhh_m; ++hh){
			k2 = Chan2.hhvec[chan_m][2*hh];
			l2 = Chan2.hhvec[chan_m][2*hh + 1];
			if(k2 == l2){ continue; }
			term0 = 0.25 * Ints2.D_ME1.V4[chan_m][pp * nhh_m + hh] * Amps2.D1.T1[chan_m][hh_m * npp_m + pp] * Amps2.D1.T1[chan_m][hh * npp_m + pp_m];
			if(print == 2){ std::cout << "!Term7 += 1/4 < " << k2 << "," << l2 << " |V| " << c2 << "," << d2 << " > * < " << c2 << "," << d2 << " |t| " << i2 << "," << j2 << " > * < " << a2 << "," << b2 << " |t| " << k2 << "," << l2 << " > = 1/4 * " << Ints2.D_ME1.V4[chan_m][pp * nhh_m + hh] << " * " << Amps2.D1.T1[chan_m][hh_m * npp_m + pp] << " * " << Amps2.D1.T1[chan_m][hh * npp_m + pp_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T1[chan_m][hh_m * npp_m + pp_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term7 = < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = ma2 + mb2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term7 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[a]);
	    minus(tb2, Space.qnums[b], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[a].j);
	    if(std::abs(Space.qnums[b].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[b].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, a, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, b, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3) J
	      //T3(ab|ij){jb,ia} = T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		term0 = Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2];
		if(print == 2){ std::cout << "Term3 += < " << c << "," << k << "^-1 |V| " << b << "," << j << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << a << "^-1 > = " << Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T2[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      tempAmps.D1.T3[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term3 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[a2]);
		    minus(tb4, Space2.qnums[b2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, a2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, b2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (3)
		    //T3(ab|ij){jb,ia} = -T2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (4)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k2 = Chan2.hp2vec[chan_m][2*hp2];
		      c2 = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(a2 == c2 || i2 == k2){ continue; }
		      term0 = -1.0 * Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2];
		      if(print == 2){ std::cout << "!Term3 += -< " << c2 << "," << k2 << "^-1 |V| " << b2 << "," << j2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = -1 * " << Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    tempAmps2.D1.T3[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term3 = < " << b2 << "," << j2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + ma2;
		      m22 = mb2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, ja2 + ma2 + jj2 + mj2) * CGC(ji2, mi2, ja2, ma2, jvec[jint], m12) * CGC(jb2, mb2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term3 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
  
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[b]);
	    minus(tb2, Space.qnums[a], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[b].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, b, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, a, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5) J
	      //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		term0 = Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2];
		if(print == 2){ std::cout << "Term5 += < " << c << "," << k << "^-1 |V| " << a << "," << j << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << b << "^-1 > = " << Ints.D_ME1.V3[chan_j][hp2 * nhp2_j + hp2_j] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps.D1.T4[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      tempAmps.D1.T5[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term5 = < " << a << "," << j << "^-1 |t| " << i << "," << b << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	    
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[b2]);
		    minus(tb4, Space2.qnums[a2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, b2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, a2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T4(ab|ij){ib,ja} = T2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (5)
		    //T5(ab|ij){ja,ib} = T2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (6)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k2 = Chan2.hp2vec[chan_m][2*hp2];
		      c2 = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(b2 == c2 || i2 == k2){ continue; }
		      term0 = Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2];
		      if(print == 2){ std::cout << "!Term5 += < " << c2 << "," << k2 << "^-1 |V| " << a2 << "," << j2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| " << i2 << "," << b2 << "^-1 > = " << Ints2.D_ME1.V3[chan_m][hp2 * nhp2_m + hp2_m] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " = " << term0 << std::endl; }
		      tempt += term0;
		    }
		    tempAmps2.D1.T4[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    tempAmps2.D1.T5[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term5 = < " << a2 << "," << j2 << "^-1 |t| " << i2 << "," << b2 << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mb2;
		      m22 = ma2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, jb2 + mb2 + jj2 + mj2) * CGC(ji2, mi2, jb2, mb2, jvec[jint], m12) * CGC(ja2, ma2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term5 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[a]);
	    minus(tb2, Space.qnums[b], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[a].j);
	    if(std::abs(Space.qnums[b].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[b].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, a, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, b, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		c = Chan.hp2vec[chan_j][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp1_j; ++hp1){
		  l = Chan.hp1vec[chan_j][2*hp1];
		  d = Chan.hp1vec[chan_j][2*hp1 + 1];
		  term0 = Ints.D_ME1.V9[chan_j][hp2 * nhp1_j + hp1] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] * Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j];
		  if(print == 2){ std::cout << "Term12 += < " << l << "," << d << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| " << i << "," << a << "^-1 > * < " << b << "," << j << "^-1 |t| " << l << "," << d << "^-1 > = " << Ints.D_ME1.V9[chan_j][hp2 * nhp1_j + hp1] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " * " << Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T2[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      if(print != 0){
		std::cout << "Term12 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[a2]);
		    minus(tb4, Space2.qnums[b2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, a2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, b2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T2(ab|ij){ia,jb} = T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (12)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k = Chan2.hp2vec[chan_m][2*hp2];
		      c = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(k == i2 || c == a2){ continue; }
		      for(int hp1 = 0; hp1 < nhp1_m; ++hp1){
			l = Chan2.hp1vec[chan_m][2*hp1];
			d = Chan2.hp1vec[chan_m][2*hp1 + 1];
			if(l == j2 || d == b2){ continue; }
			if(k == l || c == d){ continue; }
			term0 = Ints2.D_ME1.V9[chan_m][hp2 * nhp1_m + hp1] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] * Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m];
			if(print == 2){ std::cout << "!Term12 += < " << l << "," << d << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| " << i2 << "," << a2 << "^-1 > * < " << b2 << "," << j2 << "^-1 |t| " << l << "," << d << "^-1 > = " << Ints2.D_ME1.V9[chan_m][hp2 * nhp1_m + hp1] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " * " << Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term12 = < " << b << "," << j << "^-1 |t| " << i << "," << a << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + ma2;
		      m22 = mb2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, ja2 + ma2 + jj2 + mj2) * CGC(ji2, mi2, ja2, ma2, jvec[jint], m12) * CGC(jb2, mb2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term12 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }
  
    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    minus(tb1, Space.qnums[i], Space.qnums[b]);
	    minus(tb2, Space.qnums[a], Space.qnums[j]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[b].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[j].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[j].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      chan_j = ChanInd_2b_cross(Parameters.basis, Space, tb1);
	      hp1_j = Chan.hp1_map[chan_j][Hash2(i, b, Space.indtot)];
	      hp2_j = Chan.hp2_map[chan_j][Hash2(j, a, Space.indtot)];
	      nhp1_j = Chan.nhp1[chan_j];
	      nhp2_j = Chan.nhp2[chan_j];
	      //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13) J
	      for(int hp2 = 0; hp2 < nhp2_j; ++hp2){
		k = Chan.hp2vec[chan_j][2*hp2];
		d = Chan.hp2vec[chan_j][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp1_j; ++hp1){
		  l = Chan.hp1vec[chan_j][2*hp1];
		  c = Chan.hp1vec[chan_j][2*hp1 + 1];
		  term0 = Ints.D_ME1.V10[chan_j][hp2 * nhp1_j + hp1] * Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] * Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j];
		  if(print == 2){ std::cout << "Term13 += < " << l << "," << c << "^-1 |V| " << d << "," << k << "^-1 > * < " << d << "," << k << "^-1 |t| " << i << "," << b << "^-1 > * < " << a << "," << j << "^-1 |t| " << l << "," << c << "^-1 > = " << Ints.D_ME1.V10[chan_j][hp2 * nhp1_j + hp1] << " * " << Amps.D1.T2[chan_j][hp1_j * nhp2_j + hp2] << " * " << Amps.D1.T2[chan_j][hp1 * nhp2_j + hp2_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T4[chan_j][hp1_j * nhp2_j + hp2_j] += tempt;
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    minus(tb3, Space2.qnums[i2], Space2.qnums[b2]);
		    minus(tb4, Space2.qnums[a2], Space2.qnums[j2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = ChanInd_2b_cross(Parameters2.basis, Space2, tb3);
		    hp1_m = Chan2.hp1_map[chan_m][Hash2(i2, b2, Space2.indtot)];
		    hp2_m = Chan2.hp2_map[chan_m][Hash2(j2, a2, Space2.indtot)];
		    nhp1_m = Chan2.nhp1[chan_m];
		    nhp2_m = Chan2.nhp2[chan_m];
		    //T4(ab|ij){ib,ja} = T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (13)
		    for(int hp2 = 0; hp2 < nhp2_m; ++hp2){
		      k = Chan2.hp2vec[chan_m][2*hp2];
		      d = Chan2.hp2vec[chan_m][2*hp2 + 1];
		      if(k == i2 || d == b2){ continue; }
		      for(int hp1 = 0; hp1 < nhp1_m; ++hp1){
			l = Chan2.hp1vec[chan_m][2*hp1];
			c = Chan2.hp1vec[chan_m][2*hp1 + 1];
			if(l == j2 || c == a2){ continue; }
			if(k == l || d == c){ continue; }
			term0 = Ints2.D_ME1.V10[chan_m][hp2 * nhp1_m + hp1] * Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] * Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m];
			if(print == 2){ std::cout << "!Term13 += < " << l << "," << c << "^-1 |V| " << d << "," << k << "^-1 > * < " << d << "," << k << "^-1 |t| " << i2 << "," << b2 << "^-1 > * < " << a2 << "," << j2 << "^-1 |t| " << l << "," << c << "^-1 > = " << Ints2.D_ME1.V10[chan_m][hp2 * nhp1_m + hp1] << " * " << Amps2.D1.T2[chan_m][hp1_m * nhp2_m + hp2] << " * " << Amps2.D1.T2[chan_m][hp1 * nhp2_m + hp2_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T4[chan_m][hp1_m * nhp2_m + hp2_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term13 = < " << a << "," << j << "^-1 |t| " << i << "," << b << "^-1 > = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mb2;
		      m22 = ma2 + mj2;
		      if(m12 != m22 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, jb2 + mb2 + jj2 + mj2) * CGC(ji2, mi2, jb2, mb2, jvec[jint], m12) * CGC(ja2, ma2, jj2, mj2, jvec[jint], m22) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term13 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[i];
	    h_j = Chan.h_map[chan_j][i];
	    nh_j = Chan.nh[chan_j];
	    nhpp_j = Chan.nhpp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hpp_j = Chan.hpp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(j, a, b, Space.indtot)];
	      //T6(ab|ij){jab,i} = -1/2 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8) J
	      for(int h = 0; h < nh_j; ++h){
		l = Chan.hvec[chan_j][h];
		for(int hpp = 0; hpp < nhpp_j; ++hpp){
		  k = Chan.hppvec[chan_j][3*hpp];
		  c = Chan.hppvec[chan_j][3*hpp + 1];
		  d = Chan.hppvec[chan_j][3*hpp + 2];
		  term0 = -0.5 * Ints.D_ME1.V5[chan_j][h * nhpp_j + hpp] * Amps.D1.T7[chan_j][hpp_j * nh_j + h] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
		  if(print == 2){ std::cout << "Term8 += -1/2 < " << l << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i << " > * < " << a << "," << b << "," << j << "^-1 |t| " << l << "^-1 > = -1/2 * " << Ints.D_ME1.V5[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T7[chan_j][hpp_j * nh_j + h] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T6[chan_j][hpp_j * nh_j + h_j] += tempt;
	      if(print != 0){
		std::cout << "Term8 = < " << a << "," << b << "," << j << "^-1 |t| " << i << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = -0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[i2];
		    h_m = Chan2.h_map[chan_m][i2];
		    hpp_m = Chan2.hpp_map[chan_m][Hash3(j2, a2, b2, Space2.indtot)];
		    nh_m = Chan2.nh[chan_m];
		    nhpp_m = Chan2.nhpp[chan_m];
		    //T6(ab|ij){jab,i} = -1/2 * T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.T6(cd|ik){kcd,i} (8)
		    for(int h = 0; h < nh_m; ++h){
		      l = Chan2.hvec[chan_m][h];
		      if(l == j2){ continue; }  
		      for(int hpp = 0; hpp < nhpp_m; ++hpp){
			k = Chan2.hppvec[chan_m][3*hpp];
			c = Chan2.hppvec[chan_m][3*hpp + 1];
			d = Chan2.hppvec[chan_m][3*hpp + 2];
			if(k == l || c == d || k == i2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V5[chan_m][h * nhpp_m + hpp] * Amps2.D1.T7[chan_m][hpp_m * nh_m + h] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
			if(print == 2){ std::cout << "!Term8 += -1/2 < " << l << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i2 << " > * < " << a2 << "," << b2 << "," << j2 << "^-1 |t| " << l << "^-1 > = -1/2 * " << Ints2.D_ME1.V5[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T7[chan_m][hpp_m * nh_m + h] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T6[chan_m][hpp_m * nh_m + h_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term8 = < " << a2 << "," << b2 << "," << j2 << "^-1 |t| " << i2 << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = ma2 + mb2;
		      m22 = m12 + mj2;
		      if(m22 != mi2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, jj2 + mj2) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m12) * CGC(jvec[jint], m12, jj2, mj2, ji2, mi2) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term8 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[j];
	    h_j = Chan.h_map[chan_j][j];
	    nh_j = Chan.nh[chan_j];
	    nhpp_j = Chan.nhpp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hpp_j = Chan.hpp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, a, b, Space.indtot)];
	      //T7(ab|ij){iab,j} = -1/2 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9) J
	      for(int h = 0; h < nh_j; ++h){
		k = Chan.hvec[chan_j][h];
		for(int hpp = 0; hpp < nhpp_j; ++hpp){
		  l = Chan.hppvec[chan_j][3*hpp];
		  c = Chan.hppvec[chan_j][3*hpp + 1];
		  d = Chan.hppvec[chan_j][3*hpp + 2];
		  term0 = -0.5 * Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] * Amps.D1.T7[chan_j][hpp_j * nh_j + h] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
		  if(print == 2){ std::cout << "Term9 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << j << " > * < " << a << "," << b << "," << i << "^-1 |t| " << k << "^-1 > = -1/2 * " << Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T7[chan_j][hpp_j * nh_j + h] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T7[chan_j][hpp_j * nh_j + h_j] += tempt;
	      if(print != 0){
		std::cout << "Term9 = < " << a << "," << b << "," << i << "^-1 |t| " << j << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = -0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[j2];
		    h_m = Chan2.h_map[chan_m][j2];
		    hpp_m = Chan2.hpp_map[chan_m][Hash3(i2, a2, b2, Space2.indtot)];
		    nh_m = Chan2.nh[chan_m];
		    nhpp_m = Chan2.nhpp[chan_m];
		    //T7(ab|ij){iab,j} = -1/2 * T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.T6(cd|jl){lcd,j} (9)
		    for(int h = 0; h < nh_m; ++h){
		      k = Chan2.hvec[chan_m][h];
		      if(k == i2){ continue; }  
		      for(int hpp = 0; hpp < nhpp_m; ++hpp){
			l = Chan2.hppvec[chan_m][3*hpp];
			c = Chan2.hppvec[chan_m][3*hpp + 1];
			d = Chan2.hppvec[chan_m][3*hpp + 2];
			if(l == k || c == d || l == j2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] * Amps2.D1.T7[chan_m][hpp_m * nh_m + h] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
			if(print == 2){ std::cout << "!Term9 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << j2 << " > * < " << a2 << "," << b2 << "," << i2 << "^-1 |t| " << k << "^-1 > = -1/2 * " << Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T7[chan_m][hpp_m * nh_m + h] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T7[chan_m][hpp_m * nh_m + h_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term9 = < " << a2 << "," << b2 << "," << i2 << "^-1 |t| " << j2 << " > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = ma2 + mb2;
		      m22 = m12 + mi2;
		      if(m22 != mj2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, ji2 + mi2) * CGC(ja2, ma2, jb2, mb2, jvec[jint], m12) * CGC(jvec[jint], m12, ji2, mi2, jj2, mj2) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term9 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[a];
	    p_j = Chan.p_map[chan_j][a];
	    np_j = Chan.np[chan_j];
	    nhhp_j = Chan.nhhp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hhp_j = Chan.hhp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, j, b, Space.indtot)];
	      //T8(ab|ij){ijb,a} = -1/2 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10) J
	      for(int p = 0; p < np_j; ++p){
		d = Chan.pvec[chan_j][p];
		for(int hhp = 0; hhp < nhhp_j; ++hhp){
		  k = Chan.hhpvec[chan_j][3*hhp];
		  l = Chan.hhpvec[chan_j][3*hhp + 1];
		  c = Chan.hhpvec[chan_j][3*hhp + 2];
		  term0 = -0.5 * Ints.D_ME1.V7[chan_j][p * nhhp_j + hhp] * Amps.D1.T9[chan_j][hhp_j * np_j + p] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
		  if(print == 2){ std::cout << "Term10 += -1/2 < " << k << "," << l << "," << c << "^-1 |V| " << d << " > * < " << d << " |t| " << i << "," << j << "," << b << "^-1 > * < " << a << " |t| " << k << "," << l << "," << c << "^-1 > = -1/2 * " << Ints.D_ME1.V7[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T9[chan_j][hhp_j * np_j + p] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T8[chan_j][hhp_j * np_j + p_j] += tempt;
	      if(print != 0){
		std::cout << "Term10 = < " << a << " |t| " << i << "," << j << "," << b << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = 0.5*Space2.qnums[a2].m;
		    mb2 = -0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[a2];
		    p_m = Chan2.p_map[chan_m][a2];
		    hhp_m = Chan2.hhp_map[chan_m][Hash3(i2, j2, b2, Space2.indtot)];
		    np_m = Chan2.np[chan_m];
		    nhhp_m = Chan2.nhhp[chan_m];
		    //T8(ab|ij){ijb,a} = -1/2 * T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.T8(ac|kl){klc,a} (10)
		    for(int p = 0; p < np_m; ++p){
		      d = Chan2.pvec[chan_m][p];
		      if(d == b2){ continue; }  
		      for(int hhp = 0; hhp < nhhp_m; ++hhp){
			k = Chan2.hhpvec[chan_m][3*hhp];
			l = Chan2.hhpvec[chan_m][3*hhp + 1];
			c = Chan2.hhpvec[chan_m][3*hhp + 2];
			if(k == l || c == d || c == a2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V7[chan_m][p * nhhp_m + hhp] * Amps2.D1.T9[chan_m][hhp_m * np_m + p] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
			if(print == 2){ std::cout << "!Term10 += -1/2 < " << k << "," << l << "," << c << "^-1 |V| " << d << " > * < " << d << " |t| " << i2 << "," << j2 << "," << b2 << "^-1 > * < " << a2 << " |t| " << k << "," << l << "," << c << "^-1 > = -1/2 * " << Ints2.D_ME1.V7[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T9[chan_m][hhp_m * np_m + p] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T8[chan_m][hhp_m * np_m + p_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term10 = < " << a2 << " |t| " << i2 << "," << j2 << "," << b2 << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = m12 + mb2;
		      if(m22 != ma2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += std::pow(-1.0, jb2 + mb2) * CGC(jvec[jint], m12, jb2, mb2, ja2, ma2) * CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term10 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    for(int i = 0; i < Space.indhol; ++i){
      for(int j = 0; j < Space.indhol; ++j){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  for(int b = Space.indhol; b < Space.indtot; ++b){
	    plus(tb1, Space.qnums[i], Space.qnums[j]);
	    plus(tb2, Space.qnums[a], Space.qnums[b]);
	    if(tb2.j < tb1.j){ tb1.j = tb2.j; }
	    else{ tb2.j = tb1.j; }
	    minj = std::abs(Space.qnums[i].j - Space.qnums[j].j);
	    if(std::abs(Space.qnums[a].j - Space.qnums[b].j) > minj){ minj = std::abs(Space.qnums[a].j - Space.qnums[b].j); }
	    if( !(equal(tb1, tb2)) ){ continue; }
	    jsize = int(0.5*(tb1.j - minj) + 1);
	    jvec = new double[jsize];
	    tvec = new double[jsize];
	    chan_j = Chan.indvec[b];
	    p_j = Chan.p_map[chan_j][b];
	    np_j = Chan.np[chan_j];
	    nhhp_j = Chan.nhhp[chan_j];
	    while(tb1.j >= minj){
	      tempt = 0.0;
	      hhp_j = Chan.hhp_map[chan_j][int(0.5 * tb1.j * std::pow(Space.indtot, 3)) + Hash3(i, j, a, Space.indtot)];
	      //T9(ab|ij){ija,b} = -1/2 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11) J
	      for(int p = 0; p < np_j; ++p){
		c = Chan.pvec[chan_j][p];
		for(int hhp = 0; hhp < nhhp_j; ++hhp){
		  k = Chan.hhpvec[chan_j][3*hhp];
		  l = Chan.hhpvec[chan_j][3*hhp + 1];
		  d = Chan.hhpvec[chan_j][3*hhp + 2];
		  term0 = -0.5 * Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.D1.T9[chan_j][hhp_j * np_j + p] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
		  if(print == 2){ std::cout << "Term11 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << c << " |t| " << i << "," << j << "," << a << "^-1 > * < " << b << " |t| " << k << "," << l << "," << d << "^-1 > = -1/2 * " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T9[chan_j][hhp_j * np_j + p] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps.D1.T9[chan_j][hhp_j * np_j + p_j] += tempt;
	      if(print != 0){
		std::cout << "Term11 = < " << b << " |t| " << i << "," << j << "," << a << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
	      }
	      jvec[int(0.5*(tb1.j - minj))] = 0.5 * tb1.j;
	      tvec[int(0.5*(tb1.j - minj))] = tempt;
	      tb1.j -= 2;
	    }
	    if(print != 0){ std::cout << std::endl; }
	      
	    for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	      for(int j1 = 0; j1 < (Space.qnums[j].j + 1); ++j1){
		for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
		  for(int b1 = 0; b1 < (Space.qnums[b].j + 1); ++b1){
		    tempt = 0.0;
		    i2 = Space2.shellsm[i][i1];
		    j2 = Space2.shellsm[j][j1];
		    ji2 = 0.5*Space2.qnums[i2].j;
		    jj2 = 0.5*Space2.qnums[j2].j;
		    mi2 = 0.5*Space2.qnums[i2].m;
		    mj2 = 0.5*Space2.qnums[j2].m;
		    if(i2 == j2){ continue; }
		    a2 = Space2.shellsm[a][a1];
		    b2 = Space2.shellsm[b][b1];
		    ja2 = 0.5*Space2.qnums[a2].j;
		    jb2 = 0.5*Space2.qnums[b2].j;
		    ma2 = -0.5*Space2.qnums[a2].m;
		    mb2 = 0.5*Space2.qnums[b2].m;
		    if(a2 == b2){ continue; }
		    plus(tb3, Space2.qnums[i2], Space2.qnums[j2]);
		    plus(tb4, Space2.qnums[a2], Space2.qnums[b2]);
		    tb3.j = 0;
		    tb4.j = 0;
		    if( !(equal(tb3, tb4)) ){ continue; }
		    chan_m = Chan2.indvec[b2];
		    p_m = Chan2.p_map[chan_m][b2];
		    hhp_m = Chan2.hhp_map[chan_m][Hash3(i2, j2, a2, Space2.indtot)];
		    np_m = Chan2.np[chan_m];
		    nhhp_m = Chan2.nhhp[chan_m];
		    //T9(ab|ij){ija,b} = -1/2 * T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.T8(bd|kl){kld,b} (11)
		    for(int p = 0; p < np_m; ++p){
		      c = Chan2.pvec[chan_m][p];
		      if(c == a2){ continue; }  
		      for(int hhp = 0; hhp < nhhp_m; ++hhp){
			k = Chan2.hhpvec[chan_m][3*hhp];
			l = Chan2.hhpvec[chan_m][3*hhp + 1];
			d = Chan2.hhpvec[chan_m][3*hhp + 2];
			if(k == l || c == d || d == b2){ continue; }
			term0 = -0.5 * Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.D1.T9[chan_m][hhp_m * np_m + p] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
			if(print == 2){ std::cout << "!Term11 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << c << " |t| " << i2 << "," << j2 << "," << a2 << "^-1 > * < " << b2 << " |t| " << k << "," << l << "," << d << "^-1 > = -1/2 * " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T9[chan_m][hhp_m * np_m + p] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
			tempt += term0;
		      }
		    }
		    tempAmps2.D1.T9[chan_m][hhp_m * np_m + p_m] += tempt;
		    if(print != 0){
		      std::cout << "!Term11 = < " << b2 << " |t| " << i2 << "," << j2 << "," << a2 << "^-1 > ^ " << 0.5*tb1.j << " = " << tempt << std::endl;
		    }
		    double tempt2 = 0.0;
		    for(int jint = 0; jint < jsize; ++jint){
		      m12 = mi2 + mj2;
		      m22 = m12 + ma2;
		      if(m22 != mb2 || fabs(m12) > jvec[jint]){ continue; }
		      tempt2 += -1.0 * std::pow(-1.0, ja2 + ma2) * CGC(jvec[jint], m12, ja2, ma2, jb2, mb2) * CGC(ji2, mi2, jj2, mj2, jvec[jint], m12) * tvec[jint];
		    }
		    if(print != 0){ std::cout << "? Term11 = " << tempt2 << std::endl; }
		  }
		}
	      }
	    }
	    if(print != 0){ std::cout << std::endl; }
	    delete[] jvec;
	    delete[] tvec;
	  }
	}
      }
    }

    if(Parameters.approx == "singles" && Parameters2.approx == "singles"){
      //print = 1;
      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  hp1_j = Chan.hp1_map[chan0_j][Hash2(i, a, Space.indtot)];
	  //t1(ia){ia} = V3(ka|ic){ia,kc}.t1(kc){kc} (1) J
	  for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	    k = Chan.hp1vec[chan0_j][2*hp1];
	    c = Chan.hp1vec[chan0_j][2*hp1 + 1];
	    term0 = Ints.D_ME1.V3[chan0_j][hp1_j * nhp0_j + hp1] * Amps.S1.T1[hp1];
	    if(print == 2){ std::cout << "Term_1 += < " << a << "," << i << "^-1 |V| " << c << "," << k << "^-1 > * < " << c << "," << k << "^-1 |t| > = " << Ints.D_ME1.V3[chan0_j][hp1_j * nhp0_j + hp1] << " * " << Amps.S1.T1[hp1] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T1[hp1_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_1 = < " << a << "," << i << "^-1 |t| > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
	      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = -0.5*Space2.qnums[a2].m;
	      hp1_m = Chan2.hp1_map[chan0_m][Hash2(i2, a2, Space2.indtot)];
	      //t1(ia){ia} = -V3(ka|ic){ia,kc}.t1(kc){kc} (1)
	      for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
		k2 = Chan2.hp1vec[chan0_m][2*hp1];
		c2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
		term0 = -1.0 * Ints2.D_ME1.V3[chan0_m][hp1_m * nhp0_m + hp1] * Amps2.S1.T1[hp1];
		if(print == 2){ std::cout << "!Term_1 += -< " << a2 << "," << i2 << "^-1 |V| " << c2 << "," << k2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| > = -1 * " << Ints2.D_ME1.V3[chan0_m][hp1_m * nhp0_m + hp1] << " * " << Amps2.S1.T1[hp1] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T1[hp1_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_1 = < " << a2 << "," << i2 << "^-1 |t| > = " << tempt << std::endl;
	      }
	      double tempt2 = std::pow(-1.0, ja2 + ma2) * CGC(ji2, mi2, ja2, ma2, 0, 0) * TJ;
	      if(print != 0){ std::cout << "? Term_1 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  hp1_j = Chan.hp1_map[chan0_j][Hash2(i, a, Space.indtot)];
	  //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6) J
	  for(int hp2 = 0; hp2 < nhp0_j; ++hp2){
	    k = Chan.hp2vec[chan0_j][2*hp2];
	    c = Chan.hp2vec[chan0_j][2*hp2 + 1];
	    for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	      l = Chan.hp1vec[chan0_j][2*hp1];
	      d = Chan.hp1vec[chan0_j][2*hp1 + 1];
	      term0 = Ints.D_ME1.V9[chan0_j][hp2 * nhp0_j + hp1] * Amps.D1.T2[chan0_j][hp1 * nhp0_j + hp1_j] * Amps.S1.T1[hp2];
	      if(print == 2){ std::cout << "Term_6 += < " << k << "," << c << "^-1 |V| " << d << "," << l << "^-1 > * < " << c << "," << k << "^-1 |t| > * < " << d << "," << l << "^-1 |t| " << i << "," << a << "^-1 > = " << Ints.D_ME1.V9[chan0_j][hp2 * nhp0_j + hp1] << " * " << Amps.D1.T2[chan0_j][hp1 * nhp0_j + hp1_j] << " * " << Amps.S1.T1[hp2] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T1[hp1_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_6 = < " << a << "," << i << "^-1 |t| > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = -0.5*Space2.qnums[a2].m;
	      hp1_m = Chan2.hp1_map[chan0_m][Hash2(i2, a2, Space2.indtot)];
	      //t1(ia){ia} = t1(kc){kc}.V9(kl|cd){kc,ld}.T2(da|li){ld,ia} (6)
	      for(int hp2 = 0; hp2 < nhp0_m; ++hp2){
		k2 = Chan2.hp2vec[chan0_m][2*hp2];
		c2 = Chan2.hp2vec[chan0_m][2*hp2 + 1];
		for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
		  l2 = Chan2.hp1vec[chan0_m][2*hp1];
		  d2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
		  term0 = Ints2.D_ME1.V9[chan0_m][hp2 * nhp0_m + hp1] * Amps2.D1.T2[chan0_m][hp1 * nhp0_m + hp1_m] * Amps2.S1.T1[hp2];
		  if(print == 2){ std::cout << "!Term_6 += < " << k2 << "," << c2 << "^-1 |V| " << d2 << "," << l2 << "^-1 > * < " << c2 << "," << k2 << "^-1 |t| > * < " << d2 << "," << l2 << "^-1 |t| " << i2 << "," << a2 << "^-1 > = " << Ints2.D_ME1.V9[chan0_m][hp2 * nhp0_m + hp1] << " * " << Amps2.D1.T2[chan0_m][hp1 * nhp0_m + hp1_m] << " * " << Amps2.S1.T1[hp2] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T1[hp1_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_6 = < " << a2 << "," << i2 << "^-1 |t| > = " << tempt << std::endl;
	      }
	      double tempt2 = std::pow(-1.0, ja2 + ma2) * CGC(ji2, mi2, ja2, ma2, 0, 0) * TJ;
	      if(print != 0){ std::cout << "? Term_6 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = 1/2 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2) J
	  for(int hpp = 0; hpp < nhpp_j; ++hpp){
	    k = Chan.hppvec[chan_j][3*hpp];
	    c = Chan.hppvec[chan_j][3*hpp + 1];
	    d = Chan.hppvec[chan_j][3*hpp + 2];
	    term0 = 0.5 * Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] * Amps.D1.T6[chan_j][hpp * nh_j + h_j];
	    if(print == 2){ std::cout << "Term_2 += 1/2 < " << a << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |t| " << i << " > = 1/2 * " << Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_2 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = -1/2 * V11(ka|cd){a,kcd}.T6(cd|ik){kcd,i} (2)
	      for(int hpp = 0; hpp < nhpp_m; ++hpp){
		k2 = Chan2.hppvec[chan_m][2*hpp];
		c2 = Chan2.hppvec[chan_m][2*hpp + 1];
		d2 = Chan2.hppvec[chan_m][2*hpp + 2];
		term0 = -0.5 * Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m];
		if(print == 2){ std::cout << "!Term_2 += -1/2 < " << a2 << " |V| " << c2 << "," << d2 << "," << k2 << "^-1 > * < " << c2 << "," << d2 << "," << k2 << "^-1 |t| " << i2 << " > = -1/2 * " << Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_2 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_2 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = -V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3) J
	  for(int hpp = 0; hpp < nhpp_j; ++hpp){
	    k = Chan.hppvec[chan_j][3*hpp];
	    c = Chan.hppvec[chan_j][3*hpp + 1];
	    d = Chan.hppvec[chan_j][3*hpp + 2];
	    term0 = -1.0 * Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] * Amps.S1.E6[chan_j][hpp * nh_j + h_j];
	    if(print == 2){ std::cout << "Term_3 += -< " << a << " |V| " << c << "," << d << "," << k << "^-1 > * < " << c << "," << d << "," << k << "^-1 |E| " << i << " > = -1 * " << Ints.S_ME1.V11[chan_j][p_j * nhpp_j + hpp] << " * " << Amps.S1.E6[chan_j][hpp * nh_j + h_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_3 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = V11(ka|cd){a,kcd}.E6(cd|ik){kcd,i} (3)
	      for(int hpp = 0; hpp < nhpp_m; ++hpp){
		k2 = Chan2.hppvec[chan_m][3*hpp];
		c2 = Chan2.hppvec[chan_m][3*hpp + 1];
		d2 = Chan2.hppvec[chan_m][3*hpp + 2];
		term0 = Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] * Amps2.S1.E6[chan_m][hpp * nh_m + h_m];
		if(print == 2){ std::cout << "!Term_3 += < " << a2 << " |V| " << c2 << "," << d2 << "," << k2 << "^-1 > * < " << c2 << "," << d2 << "," << k2 << "^-1 |E| " << i2 << " > = -1/2 * " << Ints2.S_ME1.V11[chan_m][p_m * nhpp_m + hpp] << " * " << Amps2.S1.E6[chan_m][hpp * nh_m + h_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_3 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_3 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhpp_j = Chan.nhpp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t2(a|i){a,i} = -1/2 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7) J
	  for(int h = 0; h < nh_j; ++h){
	    k = Chan.hvec[chan_j][h];
	    for(int hpp = 0; hpp < nhpp_j; ++hpp){
	      l = Chan.hppvec[chan_j][3*hpp];
	      c = Chan.hppvec[chan_j][3*hpp + 1];
	      d = Chan.hppvec[chan_j][3*hpp + 2];
	      term0 = -0.5 * Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] * Amps.D1.T6[chan_j][hpp * nh_j + h_j] * Amps.S1.T2[chan_j][p_j * nh_j + h];
	      if(print == 2){ std::cout << "Term_7 += -1/2 < " << k << " |V| " << c << "," << d << "," << l << "^-1 > * < " << c << "," << d << "," << l << "^-1 |t| " << i << " > * < " << a << " |t| " << k << " > = -1/2 * " << Ints.D_ME1.V6[chan_j][h * nhpp_j + hpp] << " * " << Amps.D1.T6[chan_j][hpp * nh_j + h_j] << " * " << Amps.S1.T2[chan_j][p_j * nh_j + h] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T2[chan_j][p_j * nh_j + h_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_7 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t2(a|i){a,i} = -1/2 * t2(a|k){a,k}.V6(kl|cd){k,lcd}.T6(cd|il){lcd,i} (7)
	      for(int h = 0; h < nh_m; ++h){
		k2 = Chan2.hvec[chan_m][h];
		for(int hpp = 0; hpp < nhpp_m; ++hpp){
		  l2 = Chan2.hppvec[chan_m][3*hpp];
		  c2 = Chan2.hppvec[chan_m][3*hpp + 1];
		  d2 = Chan2.hppvec[chan_m][3*hpp + 2];
		  term0 = -0.5 * Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] * Amps2.D1.T6[chan_m][hpp * nh_m + h_m] * Amps2.S1.T2[chan_m][p_m * nh_m + h];
		  if(print == 2){ std::cout << "!Term_7 += -1/2 < " << k2 << " |V| " << c2 << "," << d2 << "," << l2 << "^-1 > * < " << c2 << "," << d2 << "," << l2 << "^-1 |t| " << i2 << " > * < " << a2 << " |t| " << k2 << " > = -1/2 * " << Ints2.D_ME1.V6[chan_m][h * nhpp_m + hpp] << " * " << Amps2.D1.T6[chan_m][hpp * nh_m + h_m] << " * " << Amps2.S1.T2[chan_m][p_m * nh_m + h] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T2[chan_m][p_m * nh_m + h_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_7 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_7 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = -1/2 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4) J
	  for(int hhp = 0; hhp < nhhp_j; ++hhp){
	    k = Chan.hhpvec[chan_j][3*hhp];
	    l = Chan.hhpvec[chan_j][3*hhp + 1];
	    c = Chan.hhpvec[chan_j][3*hhp + 2];
	    term0 = -0.5 * Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] * Amps.D1.T8[chan_j][hhp * np_j + p_j];
	    if(print == 2){ std::cout << "Term_4 += -1/2 < " << i << " |V| " << k << "," << l << "," << c << "^-1 > * < " << k << "," << l << "," << c << "^-1 |t| " << a << " > = -1/2 * " << Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_4 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = -1/2 * V12(kl|ic){i,klc}.T8(ac|kl){klc,a} (4)
	      for(int hhp = 0; hhp < nhhp_m; ++hhp){
		k2 = Chan2.hhpvec[chan_m][3*hhp];
		l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		c2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		term0 = -0.5 * Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] * Amps2.D1.T8[chan_m][hhp * np_m + p_m];
		if(print == 2){ std::cout << "!Term_4 += -1/2 < " << i2 << " |V| " << k2 << "," << l2 << "," << c2 << "^-1 > * < " << k2 << "," << l2 << "," << c2 << "^-1 |t| " << a2 << " > = -1/2 * " << Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_4 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_4 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5) J
	  for(int hhp = 0; hhp < nhhp_j; ++hhp){
	    k = Chan.hhpvec[chan_j][3*hhp];
	    l = Chan.hhpvec[chan_j][3*hhp + 1];
	    c = Chan.hhpvec[chan_j][3*hhp + 2];
	    term0 = Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] * Amps.S1.E8[chan_j][hhp * np_j + p_j];
	    if(print == 2){ std::cout << "Term_5 += < " << i << " |V| " << k << "," << l << "," << c << "^-1 > * < " << k << "," << l << "," << c << "^-1 |E| " << a << " > = " << Ints.S_ME1.V12[chan_j][h_j * nhhp_j + hhp] << " * " << Amps.S1.E8[chan_j][hhp * np_j + p_j] << " = " << term0 << std::endl; }
	    tempt += term0;
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_5 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhpp_m = Chan2.nhpp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = V12(kl|ic){i,klc}.E8(ac|kl){klc,a} (5)
	      for(int hhp = 0; hhp < nhhp_m; ++hhp){
		k2 = Chan2.hhpvec[chan_m][3*hhp];
		l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		c2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		term0 = Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] * Amps2.S1.E8[chan_m][hhp * np_m + p_m];
		if(print == 2){ std::cout << "!Term_5 += < " << i2 << " |V| " << k2 << "," << l2 << "," << c2 << "^-1 > * < " << k2 << "," << l2 << "," << c2 << "^-1 |E| " << a2 << " > = " << Ints2.S_ME1.V12[chan_m][h_m * nhhp_m + hhp] << " * " << Amps2.S1.E8[chan_m][hhp * np_m + p_m] << " = " << term0 << std::endl; }
		tempt += term0;
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_5 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_5 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = -1/2 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8) J
	  for(int p = 0; p < np_j; ++p){
	    c = Chan.pvec[chan_j][p];
	    for(int hhp = 0; hhp < nhhp_j; ++hhp){
	      k = Chan.hhpvec[chan_j][3*hhp];
	      l = Chan.hhpvec[chan_j][3*hhp + 1];
	      d = Chan.hhpvec[chan_j][3*hhp + 2];
	      term0 = -0.5 * Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.D1.T8[chan_j][hhp * np_j + p_j] * Amps.S1.T3[chan_j][h_j * np_j + p];
	      if(print == 2){ std::cout << "Term_8 += -1/2 < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << a << " |t| " << k << "," << l << "," << d << "^-1 > * < " << c << " |t| " << i << " > = -1/2 * " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.D1.T8[chan_j][hhp * np_j + p_j] << " * " << Amps.S1.T3[chan_j][h_j * np_j + p] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_8 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }

	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = -1/2 * t3(i|c){i,c}.V8(kl|cd){c,kld}.T8(ad|kl){kld,a} (8)
	      for(int p = 0; p < np_m; ++p){
		c2 = Chan2.pvec[chan_m][p];
		for(int hhp = 0; hhp < nhhp_m; ++hhp){
		  k2 = Chan2.hhpvec[chan_m][3*hhp];
		  l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		  d2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		  term0 = -0.5 * Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.D1.T8[chan_m][hhp * np_m + p_m] * Amps2.S1.T3[chan_m][h_m * np_m + p];
		  if(print == 2){ std::cout << "!Term_8 += -1/2 < " << k2 << "," << l2 << "," << d2 << "^-1 |V| " << c2 << " > * < " << a2 << " |t| " << k2 << "," << l2 << "," << d2 << "^-1 > * < " << c2 << " |t| " << i2 << " > = -1/2 * " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.D1.T8[chan_m][hhp * np_m + p_m] << " * " << Amps2.S1.T3[chan_m][h_m * np_m + p] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_8 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_8 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }

      for(int i = 0; i < Space.indhol; ++i){
	for(int a = Space.indhol; a < Space.indtot; ++a){
	  if(Chan.indvec[i] != Chan.indvec[a]){ continue; }
	  tempt = 0.0;
	  chan_j = Chan.indvec[i];
	  nh_j = Chan.nh[chan_j];
	  np_j = Chan.np[chan_j];
	  nhhp_j = Chan.nhhp[chan_j];
	  h_j = Chan.h_map[chan_j][i];
	  p_j = Chan.p_map[chan_j][a];
	  //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9) J
	  for(int p = 0; p < np_j; ++p){
	    c = Chan.pvec[chan_j][p];
	    for(int hhp = 0; hhp < nhhp_j; ++hhp){
	      k = Chan.hhpvec[chan_j][3*hhp];
	      l = Chan.hhpvec[chan_j][3*hhp + 1];
	      d = Chan.hhpvec[chan_j][3*hhp + 2];
	      term0 = Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] * Amps.S1.E8[chan_j][hhp * np_j+ p_j] * Amps.S1.T3[chan_j][h_j * np_j + p];
	      if(print == 2){ std::cout << "Term_9 += < " << k << "," << l << "," << d << "^-1 |V| " << c << " > * < " << a << " |E| " << k << "," << l << "," << d << "^-1 > * < " << c << " |t| " << i << " > = " << Ints.D_ME1.V8[chan_j][p * nhhp_j + hhp] << " * " << Amps.S1.E8[chan_j][hhp * np_j + p_j] << " * " << Amps.S1.T3[chan_j][h_j * np_j + p] << " = " << term0 << std::endl; }
	      tempt += term0;
	    }
	  }
	  tempAmps.S1.T3[chan_j][h_j * np_j + p_j] += tempt;
	  if(print != 0){
	    std::cout << "Term_9 = < " << a << " |t| " << i << " > = " << tempt << std::endl;
	  }
	  TJ = tempt;
	  if(print != 0){ std::cout << std::endl; }
	   	      
	  for(int i1 = 0; i1 < (Space.qnums[i].j + 1); ++i1){
	    for(int a1 = 0; a1 < (Space.qnums[a].j + 1); ++a1){
	      tempt = 0.0;
	      i2 = Space2.shellsm[i][i1];
	      a2 = Space2.shellsm[a][a1];
	      if(Chan2.indvec[i2] != Chan2.indvec[a2]){ continue; }
	      ji2 = 0.5*Space2.qnums[i2].j;
	      ja2 = 0.5*Space2.qnums[a2].j;
	      mi2 = 0.5*Space2.qnums[i2].m;
	      ma2 = 0.5*Space2.qnums[a2].m;
	      chan_m = Chan2.indvec[i2];
	      nh_m = Chan2.nh[chan_m];
	      np_m = Chan2.np[chan_m];
	      nhhp_m = Chan2.nhhp[chan_m];
	      h_m = Chan2.h_map[chan_m][i2];
	      p_m = Chan2.p_map[chan_m][a2];
	      //t3(i|a){i,a} = t3(i|c){i,c}.V8(kl|cd){c,kld}.E8(ad|kl){kld,a} (9)
	      for(int p = 0; p < np_m; ++p){
		c2 = Chan2.pvec[chan_m][p];
		for(int hhp = 0; hhp < nhhp_m; ++hhp){
		  k2 = Chan2.hhpvec[chan_m][3*hhp];
		  l2 = Chan2.hhpvec[chan_m][3*hhp + 1];
		  d2 = Chan2.hhpvec[chan_m][3*hhp + 2];
		  term0 = Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] * Amps2.S1.E8[chan_m][hhp * np_m + p_m] * Amps2.S1.T3[chan_m][h_m * np_m + p];
		  if(print == 2){ std::cout << "!Term_9 += < " << k2 << "," << l2 << "," << d2 << "^-1 |V| " << c2 << " > * < " << a2 << " |E| " << k2 << "," << l2 << "," << d2 << "^-1 > * < " << c2 << " |t| " << i2 << " > = " << Ints2.D_ME1.V8[chan_m][p * nhhp_m + hhp] << " * " << Amps2.S1.E8[chan_m][hhp * np_m + p_m] << " * " << Amps2.S1.T3[chan_m][h_m * np_m + p] << " = " << term0 << std::endl; }
		  tempt += term0;
		}
	      }
	      tempAmps2.S1.T3[chan_m][h_m * np_m + p_m] += tempt;
	      if(print != 0){
		std::cout << "!Term_9 = < " << a2 << " |t| " << i2 << " > = " << tempt << std::endl;
	      }
	      double tempt2 = CGC(ja2, ma2, 0, 0, ji2, mi2) * TJ;
	      if(print != 0){ std::cout << "? Term_9 = " << tempt2 << std::endl; }
	    }
	  }
	}
      }
    }
    //print = 0;

    /*double fac1 = 1.0, fac2 = 0.0, fac3 = 0.5, fac5 = -1.0, fac6 = -0.5;
    char N = 'N';
    for(int chan = 0; chan < Chan.size1; ++chan){
      int hh = Chan.nhh[chan];
      int pp = Chan.npp[chan];
      if(hh * pp == 0){ continue; }
      //T1(ab|ij){ij,ab} = -E1(cd|ij){ij,cd}.V1(ab|cd){cd,ab} (1)
      dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V1[chan], Amps2.D1.T1[chan], &hh, &pp, &pp, &fac5, &fac1, &N, &N); //fac1
      //T1(ab|ij){ij,ab} = -V2(ij|kl){ij,kl}.E1(ab|kl){kl,ab} (2)
      dgemm_NN(Ints.D_ME1.V2[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac5, &fac1, &N, &N); //fac1
      //T1(ab|ij){ij,ab} = -0.5 * T1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (3)
      dgemm_NN(Amps1.D1.T1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
      //T1(ab|ij){ij,ab} = -0.5 * E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.T1(ab|kl){kl,ab} (4)
      dgemm_NN(Amps1.S1.E1[chan], Ints.D_ME1.V4[chan], Amps2.D1.S1[chan], &hh, &hh, &pp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps2.D1.S1[chan], Amps1.D1.T1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac6, &fac1, &N, &N); //fac1
      //T1(ab|ij){ij,ab} = E1(cd|ij){ij,cd}.V4(kl|cd){cd,kl}.E1(ab|kl){kl,ab} (5)
      dgemm_NN(Amps2.D1.S1[chan], Amps1.S1.E1[chan], Amps2.D1.T1[chan], &hh, &pp, &hh, &fac1, &fac1, &N, &N); //fac1
    }
    for(int chan = 0; chan < Chan.size2; ++chan){
      int hp1 = Chan.nhp1[chan];
      int hp2 = Chan.nhp2[chan];
      if(hp1 == 0 || hp2 == 0){ continue; }
      //T2(ab|ij){ia,jb} = E2(ac|ik){ia,kc}.V3(kb|jc){kc,jb} (6)
      dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
      //T3(ab|ij){jb,ia} = E2(bc|jk){jb,kc}.V3(ka|ic){kc,ia} (7)
      dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T3[chan], &hp1, &hp2, &hp2, &fac1, &fac1, &N, &N); //fac1
      //T4(ab|ij){ib,ja} = -E2(bc|ik){ib,kc}.V3(ka|jc){kc,ja} (8)
      dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T5(ab|ij){ja,ib} = -E2(ac|jk){ja,kc}.V3(kb|ic){kc,ib} (9)
      dgemm_NN(Amps1.S1.E2[chan], Ints.D_ME1.V3[chan], Amps2.D1.T5[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T2(ab|ij){ia,jb} = -E2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.T2(db|lj){ld,jb} (10)
      dgemm_NN(Ints.D_ME1.V9[chan], Amps1.D1.T2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T2(ab|ij){ia,jb} = -T2(ac|ik){ia,kc}.V9(kl|cd){kc,ld}.E2(db|lj){ld,jb} (11)
      dgemm_NN(Ints.D_ME1.V9[chan], Amps1.S1.E2[chan], Amps2.D1.S6[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S6[chan], Amps2.D1.T2[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T4(ab|ij){ib,ja} = -E2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.T2(ca|lj){lc,ja} (12)
      dgemm_NN(Ints.D_ME1.V10[chan], Amps1.D1.T2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.S1.E2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T4(ab|ij){ib,ja} = -T2(bd|ik){ib,kd}.V10(kl|cd){kd,lc}.E2(ca|lj){lc,ja} (13)
      dgemm_NN(Ints.D_ME1.V10[chan], Amps1.S1.E2[chan], Amps2.D1.S7[chan], &hp2, &hp2, &hp1, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T2[chan], Amps2.D1.S7[chan], Amps2.D1.T4[chan], &hp1, &hp2, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T2(ab|ij){ia,jb} = -Q12(ad|ik){ia,kd}.T2(db|kj){kd,jb} (14)
      dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T3(ab|ij){jb,ia} = -Q12(bd|jk){jb,kd}.T2(da|ki){kd,ia} (15)
      dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T4(ab|ij){ib,ja} = Q12(bd|ik){ib,kd}.T2(da|kj){kd,ja} (16)
      dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
      //T5(ab|ij){ja,ib} = Q12(ad|jk){ja,kd}.T2(db|ki){kd,ib} (17)
      dgemm_NN(Amps1.S1.Q12[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
      //T2(ab|ij){ia,jb} = -Q22(ac|il){ia,kc}.T2(cb|lj){kc,jb} (18)
      dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T2[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T3(ab|ij){jb,ia} = -Q22(bc|jl){jb,kc}.T2(ca|li){kc,ia} (19)
      dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T3[chan], &hp1, &hp1, &hp2, &fac5, &fac1, &N, &N); //fac1
      //T4(ab|ij){ib,ja} = Q22(bc|il){ib,kc}.T2(ca|lj){kc,ja} (20)
      dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T4[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
      //T5(ab|ij){ja,ib} = Q22(ac|jl){ja,kc}.T2(cb|li){kc,ib} (21)
      dgemm_NN(Amps1.S1.Q22[chan], Amps1.D1.T2[chan], Amps2.D1.T5[chan], &hp1, &hp1, &hp2, &fac1, &fac1, &N, &N); //fac1
    }
    for(int chan = 0; chan < Chan.size3; ++chan){
      int p = Chan.np[chan];
      int h = Chan.nh[chan];
      int hpp = Chan.nhpp[chan];
      if(h == 0 || hpp == 0){ continue; }
      //T6(ab|ij){jab,i} = T7(ab|jl){jab,l}.V5(kl|cd){l,kcd}.E6(cd|ik){kcd,i} (22)
      dgemm_NN(Ints.D_ME1.V5[chan], Amps1.S1.E6[chan], Amps2.D1.S2[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S2[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = T7(ab|ik){iab,k}.V6(kl|cd){k,lcd}.E6(cd|jl){lcd,j} (23)
      dgemm_NN(Ints.D_ME1.V6[chan], Amps1.S1.E6[chan], Amps2.D1.S3[chan], &h, &h, &hpp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T7[chan], Amps2.D1.S3[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
      //T6(ab|ij){jab,i} = T6(ab|lj){jab,l}.Q32(l|i){l,i} (26)
      dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T6[chan], &hpp, &h, &h, &fac1, &fac1, &N, &N); //fac1
      //T7(ab|ij){iab,j} = -T6(ab|li){iab,l}.Q32(l|j){l,j} (27)
      dgemm_NN(Amps1.D1.T6[chan], Amps1.S1.Q32[chan], Amps2.D1.T7[chan], &hpp, &h, &h, &fac5, &fac1, &N, &N); //fac1
      if(p != 0){
	//T6(ab|ij){jab,i} = -V17(jc|ab){jab,c}.t2(c|i){c,i} (30)
	dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
	//T7(ab|ij){iab,j} = V17(ic|ab){iab,c}.t2(c|j){c,j} (31)
	dgemm_NN(Ints.S_ME1.V17[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
	//T6(ab|ij){jab,i} = -0.5 * DQ12(ab|jc){jab,c}.t2(c|i){c,i} (34)
	dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac6, &fac1, &N, &N); //fac1
	//T7(ab|ij){iab,j} = 0.5 * DQ12(ab|ic){iab,c}.t2(c|j){c,j} (35)
	dgemm_NN(Amps1.D1.Q12[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac3, &fac1, &N, &N); //fac1
	//T6(ab|ij){jab,i} = Q52(jc|ab){jab,c}.t2(c|i){c,i} (38)
	dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T6[chan], &hpp, &h, &p, &fac1, &fac1, &N, &N); //fac1
	//T7(ab|ij){iab,j} = -Q52(ic|ab){iab,c}.t2(c|j){c,j} (39)
	dgemm_NN(Amps1.S1.Q52[chan], Amps1.S1.T2[chan], Amps2.D1.T7[chan], &hpp, &h, &p, &fac5, &fac1, &N, &N); //fac1
      }
    }
    for(int chan = 0; chan < Chan.size3; ++chan){
      int h = Chan.nh[chan];
      int p = Chan.np[chan];
      int hhp = Chan.nhhp[chan];
      if(p == 0 || hhp == 0){ continue; }
      //T8(ab|ij){ijb,a} = T9(bd|ij){ijb,d}.V7(kl|cd){d,klc}.E8(ac|kl){klc,a} (24)
      dgemm_NN(Ints.D_ME1.V7[chan], Amps1.S1.E8[chan], Amps2.D1.S4[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S4[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = T9(ac|ij){ija,c}.V8(kl|cd){c,kld}.E8(bd|kl){kld,b} (25)
      dgemm_NN(Ints.D_ME1.V8[chan], Amps1.S1.E8[chan], Amps2.D1.S5[chan], &p, &p, &hhp, &fac1, &fac2, &N, &N);
      dgemm_NN(Amps1.D1.T9[chan], Amps2.D1.S5[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
      //T8(ab|ij){ijb,a} = T8(db|ij){ijb,d}*Q42(d|a){d,a} (28)
      dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T8[chan], &hhp, &p, &p, &fac1, &fac1, &N, &N); //fac1
      //T9(ab|ij){ija,b} = -T8(da|ij){ija,d}*Q42(d|b){d,b} (29)
      dgemm_NN(Amps1.D1.T8[chan], Amps1.S1.Q42[chan], Amps2.D1.T9[chan], &hhp, &p, &p, &fac5, &fac1, &N, &N); //fac1
      if(h != 0){
	//T8(ab|ij){ijb,a} = -V18(ij|kb){ijb,k}.t3(k|a){k,a} (32)
	dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
	//T9(ab|ij){ija,b} = V18(ij|ka){ija,k}.t3(k|b){k,b} (33)
	dgemm_NN(Ints.S_ME1.V18[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
	//T8(ab|ij){ijb,a} = -0.5 * DQ22(kb|ij){ijb,k}.t3(k|a){k,a} (36)
	dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac6, &fac1, &N, &N); //fac1
	//T9(ab|ij){ija,b} = 0.5 * DQ22(ka|ij){ija,k}.t3(k|b){k,b} (37)
	dgemm_NN(Amps1.D1.Q22[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac3, &fac1, &N, &N); //fac1
	//T8(ab|ij){ijb,a} = -Q62(ij|kb){ijb,k}.t3(k|a){k,a} (40)
	dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T8[chan], &hhp, &p, &h, &fac5, &fac1, &N, &N); //fac1
	//T9(ab|ij){ija,b} = Q62(ij|ka){ija,k}.t3(k|b){k,b} (41)
	dgemm_NN(Amps1.S1.Q62[chan], Amps1.S1.T3[chan], Amps2.D1.T9[chan], &hhp, &p, &h, &fac1, &fac1, &N, &N); //fac1
      }
      }*/

    Amps.zero(Parameters, Chan);
    Amps2.zero(Parameters2, Chan2);
	
    for(int chan = 0; chan < Chan.size1; ++chan){
      nhh_j = Chan.nhh[chan];
      npp_j = Chan.npp[chan];
      if(nhh_j * npp_j == 0){ continue; }
      for(int pp = 0; pp < npp_j; ++pp){
	a2 = Chan.ppvec[chan][2*pp];
	b2 = Chan.ppvec[chan][2*pp + 1];
	for(int hh = 0; hh < nhh_j; ++hh){
	  i2 = Chan.hhvec[chan][2*hh];
	  j2 = Chan.hhvec[chan][2*hh + 1];
	  ind = hh * npp_j + pp;
	  //tempt = tempAmps.D1.T1[chan][ind];
	  tempt = tempAmps.D1.get_TJ(Space, Chan, chan, ind, i2, j2, a2, b2);
	  tempen = Space.qnums[i2].energy + Space.qnums[j2].energy - Space.qnums[a2].energy - Space.qnums[b2].energy;
	  tempt /= tempen;
	  //std::cout << "TJ: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > ^ " << 0.5*Chan.qnums1[chan].j << " = " << tempt << std::endl;
	  Amps.D1.set_TJ(Space, Chan, chan, ind, i2, j2, a2, b2, tempt);
	}
      }
    }
    if(Parameters.approx == "singles"){
      for(int hp1 = 0; hp1 < nhp0_j; ++hp1){
	i2 = Chan.hp1vec[chan0_j][2*hp1];
	a2 = Chan.hp1vec[chan0_j][2*hp1 + 1];
	//tempt = tempAmps.S1.T1[hp1];
	tempt = tempAmps.S1.get_TJ(Space, hp1, i2, a2);
	tempen = Space.qnums[i2].energy - Space.qnums[a2].energy;
	tempt /= tempen;
	//std::cout << "tj: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
	Amps.S1.set_TJ(Space, hp1, i2, a2, tempt);
      }
      Amps.S1.set_T_2J(Space, Chan, Ints);
      //Amps.D1.set_T_2J(Chan, Ints);
    }
    for(int chan = 0; chan < Chan2.size1; ++chan){
      nhh_m = Chan2.nhh[chan];
      npp_m = Chan2.npp[chan];
      if(nhh_m * npp_m == 0){ continue; }
      for(int pp = 0; pp < npp_m; ++pp){
	a2 = Chan2.ppvec[chan][2*pp];
	b2 = Chan2.ppvec[chan][2*pp + 1];
	if(a2 == b2){ continue; }
	for(int hh = 0; hh < nhh_m; ++hh){
	  i2 = Chan2.hhvec[chan][2*hh];
	  j2 = Chan2.hhvec[chan][2*hh + 1];
	  if(i2 == j2){ continue; }
	  ind = hh * npp_m + pp;
	  //tempt = tempAmps2.D1.T1[chan][ind];
	  tempt = tempAmps2.D1.get_T(chan, ind);
	  tempen = Space2.qnums[i2].energy + Space2.qnums[j2].energy - Space2.qnums[a2].energy - Space2.qnums[b2].energy;
	  tempt /= tempen;
	  //std::cout << "TM: < " << a2 << "," << b2 << " |t| " << i2 << "," << j2 << " > = " << tempt << std::endl;
	  Amps2.D1.set_T(chan, ind, tempt);
	}
      }
    }
    if(Parameters2.approx == "singles"){
      for(int hp1 = 0; hp1 < nhp0_m; ++hp1){
	i2 = Chan2.hp1vec[chan0_m][2*hp1];
	a2 = Chan2.hp1vec[chan0_m][2*hp1 + 1];
	//tempt = tempAmps2.S1.T1[hp1];
	tempt = tempAmps2.S1.get_T(hp1);
	tempen = Space2.qnums[i2].energy - Space2.qnums[a2].energy;
	tempt /= tempen;
	//std::cout << "tm: " << i2 << " |t| " << a2 << " = " << tempt << std::endl;
	Amps2.S1.set_T(hp1, tempt);
      }
      Amps2.S1.set_T_2(Chan2, Ints2);
      //Amps2.D1.set_T_2(Chan2, Ints2);
    }
    energy1 = Amps.get_energy(Parameters, Chan, Ints);
    energy2 = Amps2.get_energy(Parameters2, Chan2, Ints2);
    std::cout << "!?!? " << energy1 << " " << energy2 << std::endl;
	
    tempAmps.zero(Parameters, Chan);
    tempAmps2.zero(Parameters2, Chan2);
  }
  tempAmps.delete_struct(Parameters, Chan);
  tempAmps2.delete_struct(Parameters2, Chan2);
  Ints2.delete_struct(Parameters2, Chan2);
  Amps2.delete_struct(Parameters2, Chan2);
  Chan2.delete_struct();
  Space2.delete_struct(Parameters2);
}
