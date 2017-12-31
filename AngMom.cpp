#include "AngMom.hpp"
#include "CoupledCluster.hpp"
#include "BASISfunctions.hpp"
#include "MATHfunctions.hpp"

/*void AngMom::Build()
{
  int ind1, ind2, indL, indR, indJ;
  int mm2, maxj, maxJ, Jmin1, Jmax1, Jmin2, Jmax2;
  long long cgcsize, sixjsize;

  int maxj2 = 1;
  int maxj3 = 1;
  for(int i = 0; i < SPB.num_states; ++i){
    for(int j = i; j < SPB.num_states; ++j){
      if( SPB.qnums[i].j + SPB.qnums[j].j > maxj2 ){
	maxj2 = SPB.qnums[i].j + SPB.qnums[j].j;
      }
      for(int k = j; k < SPB.num_states; ++k){
	if( SPB.qnums[i].j + SPB.qnums[j].j + SPB.qnums[k].j > maxj3 ){
	  maxj3 = SPB.qnums[i].j + SPB.qnums[j].j + SPB.qnums[k].j;
	}
      }
    }
  }

  maxj = maxj3;
  maxJ = maxj2;
  this->jsize = (maxj + 1)/2;
  this->jsize = maxJ + 1;
  //threejsize = (this->jsize + 1)*(this->jsize + 2)*(this->jsize + 3)*(5 + 11*this->jsize + 4*this->jsize*this->jsize)/30;
  //sixjsize = ((this->jsize*this->jsize)*(this->jsize*this->jsize + 1)/4)*this->Jsize*this->Jsize;
  std::cout << "maxj = " << maxj << std::endl;
  std::cout << "maxJ = " << maxJ << std::endl;
  std::cout << "jsize = " << this->jsize << std::endl;
  std::cout << "Jsize = " << this->Jsize << std::endl;
  // std::cout << "CGC Size = " << cgcsize << std::endl;
  //std::cout << "Six_J Size = " << sixjsize << std::endl;

  //  Setup CGC
  this->CGC_vec = new double**[2 * this->jsize * (this->jsize + 1)]; //  j1, -j1 <= m <= j1
  for(int j1 = 0; j1 < this->jsize; ++j1){
    for(int m1 = 0; m1 < 2*(j1 + 1); ++m1){
      ind1 = 2*j1*this->jsize + m1;
      this->CGC_vec[ind1] = new double*[j1 + 1];

      for(int j2 = 0; j2 <= j1; ++j2){  //  j2 <= j1
	Jmax1 = std::min(j1 + j2 + 1, maxJ/2);
	Jmin1 = std::abs(j1 - j2);
	this->CGC_vec[ind1][j2] = new double[this->Jsize * this->Jsize];
	for(indJ = 0; indJ < this->Jsize*this->Jsize; ++indJ){ this->CGC_vec[ind1][j2][indJ] = 0.0; }

        #pragma omp parallel for collapse(2)
	for(int J = Jmin1; J <= Jmax1; ++J){
	  for(int M = 0; M <= Jmax1; ++M){
	    mm2 = 2*M - (2*(m1 - j1) - 1);  // m-projection, not index
	    if( std::abs(mm2) > (2*j2 + 1) || M > J ){ continue; }
	    if( logfac(2*j1 + 1) + logfac(2*j2 + 1) + logfac(2*J) < 30.0 ){
	      this->CGC_vec[ind1][j2][J*this->Jsize + M] = this->CGC_1(2*j1 + 1, 2*(m1 - j1) - 1, 2*j2 + 1, mm2, 2*J, 2*M);
	    }
	    else{
	      this->CGC_vec[ind1][j2][J*this->Jsize + M] = this->CGC_large_1(2*j1 + 1, 2*(m1 - j1) - 1, 2*j2 + 1, mm2, 2*J, 2*M);
	    }
	  }
	}
      }
    }
  }

  //  Setup SixJ
  this->SixJ_vec = new double**[this->jsize * this->jsize];
  for(int j1 = 0; j1 < this->jsize; ++j1){
    for(int j2 = 0; j2 < this->jsize; ++j2){
      ind1 = j1*this->jsize + j2;
      this->SixJ_vec[ind1] = new double*[ind1 + 1];

      for(int jj1 = 0; jj1 < this->jsize; ++jj1){
	indL = j1*this->jsize + jj1;
	for(int jj2 = 0; jj2 < this->jsize; ++jj2){
	  ind2 = jj1*this->jsize + jj2;
	  indR = j2*this->jsize + jj2;
	  if( ind2 > ind1 ){ continue; }
	  this->SixJ_vec[ind1][ind2] = new double[this->Jsize * this->Jsize];
	  for(indJ = 0; indJ < this->Jsize*this->Jsize; ++indJ){ this->SixJ_vec[ind1][ind2][indJ] = 0.0; }
	  if( indR > indL ){ continue; }

	  Jmax1 = std::min(j1 + j2 + 1, jj1 + jj2 + 1);
	  Jmax1 = std::min(Jmax1, maxJ/2);
	  Jmin1 = std::max(std::abs(j1 - j2), std::abs(jj1 - jj2));
	  if( Jmin1 > Jmax1 ){ continue; }

	  Jmax2 = std::min(j1 + jj2 + 1, jj1 + j2 + 1);
	  Jmax2 = std::min(Jmax2, maxJ/2);
	  Jmin2 = std::max(std::abs(j1 - jj2), std::abs(jj1 - j2));
	  if( Jmin2 > Jmax2 ){ continue; }

          #pragma omp parallel for collapse(2)
	  for(int J = Jmin1; J <= Jmax1; ++J){
	    for(int JJ = Jmin2; JJ <= Jmax2; ++JJ){
	      this->SixJ_vec[ind1][ind2][J*this->Jsize + JJ] = this->CGC6_1(2*j1 + 1, 2*j2 + 1, 2*J, 2*jj1 + 1, 2*jj2 + 1, 2*JJ);
	    }
	  }
	}
      }
    }
  }
}

double AngMom::get_CGC(int j1, int m1, int j2, int m2, int jtot, int mtot)
{
  if( this->CGC_vec == NULL ){ return 0.0; }
  if( j1%2 && j2%2 && !(jtot%2) ){ return this->get_CGC0(j1, m1, j2, m2, jtot, mtot); }
  else if( !(j1%2) && j2%2 && jtot%2 ){ return this->get_CGC0(jtot, -mtot, j2, m2, j1, -m1) * std::sqrt((jtot + 1.0)/(j1 + 1.0)) * phase2(j2 + m2); }
  else if( j1%2 && !(j2%2) && jtot%2 ){ return this->get_CGC0(j1, m1, jtot, -mtot, j2, -m2) * std::sqrt((jtot + 1.0)/(j2 + 1.0)) * phase2(j1 - m1); }
  else{ std::cerr << "AngMom::get_CGC Error!!!" << std::endl; exit(1); }
}

double AngMom::get_CGC0(int j1, int m1, int j2, int m2, int jtot, int mtot)
{
  if(m1 + m2 - mtot != 0 || std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(mtot) > jtot || std::abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  int ind1, indJ;
  int jj1, mm1, jj2, mm2, mmtot;
  if( j1 >= j2 ){
    jj1 = j1;
    jj2 = j2;
    if( mtot >= 0 ){
      mm1 = m1;
      mm2 = m2;
      mmtot = mtot;

      ind1 = (jj1 - 1)*this->jsize + ((mm1 + jj1)/2);
      indJ = (jtot/2)*this->Jsize + (mmtot/2);
      return this->CGC_vec[ind1][(jj2 - 1)/2][indJ];
    }
    else{
      mm1 = -m1;
      mm2 = -m2;
      mmtot = -mtot;

      ind1 = (jj1 - 1)*this->jsize + ((mm1 + jj1)/2);
      indJ = (jtot/2)*this->Jsize + (mmtot/2);
      return this->CGC_vec[ind1][(jj2 - 1)/2][indJ] * phase2(j1 + j2 - jtot);
    }
  }
  else{
    jj1 = j2;
    jj2 = j1;
    if( mtot >= 0 ){
      mm1 = m2;
      mm2 = m1;
      mmtot = mtot;

      ind1 = (jj1 - 1)*this->jsize + ((mm1 + jj1)/2);
      indJ = (jtot/2)*this->Jsize + (mmtot/2);
      return this->CGC_vec[ind1][(jj2 - 1)/2][indJ] * phase2(j1 + j2 - jtot);
    }
    else{
      mm1 = -m2;
      mm2 = -m1;
      mmtot = -mtot;

      ind1 = (jj1 - 1)*this->jsize + ((mm1 + jj1)/2);
      indJ = (jtot/2)*this->Jsize + (mmtot/2);
      return this->CGC_vec[ind1][(jj2 - 1)/2][indJ];
    }
  }
}

double AngMom::get_SixJ(int j1, int j2, int J, int jj1, int jj2, int JJ)
{
  if( j1 == 0 || this->SixJ_vec == NULL ){ return 1.0; } // no angular momentum coupling
  int ind1_1 = ((j1 - 1)/2)*this->jsize + ((j2 - 1)/2);
  int ind2_1 = ((jj1 - 1)/2)*this->jsize + ((jj2 - 1)/2);
  int ind1_2 = ((j2 - 1)/2)*this->jsize + ((j1 - 1)/2);
  int ind2_2 = ((jj2 - 1)/2)*this->jsize + ((jj1 - 1)/2);
  int indL_1 = ((j1 - 1)/2)*this->jsize + ((jj1 - 1)/2);
  int indR_1 = ((j2 - 1)/2)*this->jsize + ((jj2 - 1)/2);
  int indL_2 = ((jj1 - 1)/2)*this->jsize + ((j1 - 1)/2);
  int indR_2 = ((jj2 - 1)/2)*this->jsize + ((j2 - 1)/2);
  int indJ = (J/2)*this->Jsize + (JJ/2);
  if( ind1_1 >= ind2_1 && indL_1 >= indR_1 ){
    return this->SixJ_vec[ind1_1][ind2_1][indJ];
  }
  else if( ind1_1 < ind2_1 && indL_2 >= indR_2 ){ return this->SixJ_vec[ind2_1][ind1_1][indJ]; }
  else if( ind1_2 >= ind2_2 && indL_1 < indR_1 ){ return this->SixJ_vec[ind1_2][ind2_2][indJ]; }
  else if( ind1_2 < ind2_2 && indL_2 < indR_2 ){ return this->SixJ_vec[ind2_2][ind1_2][indJ]; }
  else{ return 0.0; }
}

double AngMom::get_NineJ(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9)
{
  double ninej0, ninej;
  int jmin = std::max(std::abs(j1 - j9), std::abs(j2 - j6));
  jmin = std::max(jmin, std::abs(j4 - j8));
  int jmax = std::min(j1 + j9, j2 + j6);
  jmax = std::min(jmax, j4 + j8);

  ninej = 0.0;
  for(int x = jmin; x <= jmax; x += 2){
    ninej0 = (x + 1.0) * phase2(x);
    ninej0 *= this->get_SixJ(j1, j4, j7, j8, j9, x);
    ninej0 *= this->get_SixJ(j2, j5, j8, j4, x, j6);
    ninej0 *= this->get_SixJ(j3, j6, j9, x, j1, j2);
    ninej += ninj0;
  }
  return ninej;
}

void AngMom::Delete()
{
  int ind1, ind2;
  for(int j1 = 0; j1 < this->jsize; ++j1){
    for(int m1 = 0; m1 <= 2*j1; ++m1){
      ind1 = 2*j1*this->jsize + m1;
      for(int j2 = 0; j2 <= j1; ++j2){
	if( this->CGC_vec[ind1][j1] != NULL ){ delete[] this->CGC_vec[ind1][j2]; }
      }
      if( this->CGC_vec[ind1] != NULL ){ delete[] this->CGC_vec[ind1]; }
    }
  }
  if( this->CGC_vec != NULL ){ delete[] this->CGC_vec; }

  for(int j1 = 0; j1 < this->jsize; ++j1){
    for(int j2 = 0; j2 < this->jsize; ++j2){
      ind1 = j1*this->jsize + j2;
      for(int jj1 = 0; jj1 < this->jsize; ++jj1){
	for(int jj2 = 0; jj2 < this->jsize; ++jj2){
	  ind2 = jj1*this->jsize + jj2;
	  if( ind2 > ind1 ){ continue; }
	  if( this->SixJ_vec[ind1][ind2] != NULL ){ delete[] this->SixJ_vec[ind1][ind2]; }
	}
      }
      if( this->SixJ_vec[ind1] != NULL ){ delete[] this->SixJ_vec[ind1]; }
    }
  }
  if( this->SixJ_vec != NULL ){ delete[] this->SixJ_vec; }
}

// for 2*j, assume j1 >= j2, abs(m1) <= j1, abs(m2) <= j2, abs(mtot) <= jtot, abs(j1 - j2) <= jtot <= (j1 + j2), mtot = m1 + m2
double AngMom::CGC_1(int j1, int m1, int j2, int m2, int jtot, int mtot)
{
  double fac1, fac2, fac3, CGC;
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  long long num1, den1, n2_1, n2_2, n2_3;

  num1 = (jtot + 1.0) * fac((jtot + j1 - j2)/2) * fac((jtot - j1 + j2)/2) * fac((j1 + j2 - jtot)/2);
  den1 = fac((j1 + j2 + jtot + 2)/2);
  fac1 = std::sqrt(num1) / std::sqrt(den1);
  n2_1 = fac((jtot + mtot)/2) * fac((jtot - mtot)/2);
  n2_2 = fac((j1 - m1)/2) * fac((j1 + m1)/2);
  n2_3 = fac((j2 - m2)/2) * fac((j2 + m2)/2);
  fac2 = std::sqrt(n2_1) * std::sqrt(n2_2) * std::sqrt(n2_3);
  maxk = std::min(j1 + j2 - jtot, j1 - m1);
  maxk = std::min(maxk, j2 + m2);

  fac3 = 0.0;
  for(int k = 0; k <= maxk; k += 2){
    d3_1 = j1 + j2 - jtot - k;
    d3_2 = j1 - m1 - k;
    d3_3 = j2 + m2 - k;
    d3_4 = jtot - j2 + m1 + k;
    d3_5 = jtot - j1 - m2 + k;
    if (d3_1 >= 0 && d3_2 >= 0 && d3_3 >= 0 && d3_4 >= 0 && d3_5 >= 0){
      fac3 += ((((((phase2(k) / fac(k/2)) / fac(d3_1/2)) / fac(d3_2/2)) / fac(d3_3/2)) / fac(d3_4/2)) / fac(d3_5/2));
    }
  }
  CGC = fac1*fac2*fac3;
  return CGC;
}

double AngMom::CGC_large_1(int j1, int m1, int j2, int m2, int jtot, int mtot)
{
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  double fac1, fac2, fac3, fac3_1, CGC, num1, den1, n2_1, n2_2, n2_3;

  num1 = std::log(jtot + 1) + logfac((jtot + j1 - j2)/2) + logfac((jtot - j1 + j2)/2) + logfac((j1 + j2 - jtot)/2);
  den1 = logfac((j1 + j2 + jtot + 2)/2);
  fac1 = std::exp(0.5 * (num1 - den1));
  n2_1 = logfac((jtot + mtot)/2) + logfac((jtot - mtot)/2);
  n2_2 = logfac((j1 - m1)/2) + logfac((j1 + m1)/2);
  n2_3 = logfac((j2 - m2)/2) + logfac((j2 + m2)/2);
  fac2 = std::exp(0.5 * (n2_1 + n2_2 + n2_3));
  maxk = std::min(j1 + j2 - jtot, j1 - m1);
  maxk = std::min(maxk, j2 + m2);

  fac3 = 0.0;
  for(int k = 0; k <= maxk; k += 2){
    d3_1 = j1 + j2 - jtot - k;
    d3_2 = j1 - m1 - k;
    d3_3 = j2 + m2 - k;
    d3_4 = jtot - j2 + m1 + k;
    d3_5 = jtot - j1 - m2 + k;
    if (d3_1 >= 0 && d3_2 >= 0 && d3_3 >= 0 && d3_4 >= 0 && d3_5 >= 0){
      fac3_1 = -logfac(k/2) - logfac(d3_1/2) - logfac(d3_2/2) - logfac(d3_3/2) - logfac(d3_4/2) - logfac(d3_5/2);
      fac3 += phase2(k) * std::exp(fac3_1);
    }
  }
  CGC = fac1*fac2*fac3;
  return CGC;
}

// for 2*j, assume j1 >= j2, abs(m1) <= j1, abs(m2) <= j2, abs(mtot) <= jtot, abs(j1 - j2) <= jtot <= (j1 + j2)
double AngMom::CGC3_1(int j1, int m1, int j2, int m2, int jtot, int mtot) // for 2*j
{
  if(m1 + m2 + mtot != 0 || abs(m1) > j1 || abs(m2) > j2 || abs(mtot) > jtot || abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  else{ return phase2(j1 - j2 - mtot) * this->get_CGC(j1,m1,j2,m2,jtot,-mtot) / std::sqrt(jtot + 1.0); }
}

double AngMom::CGC6_1(int j1, int j2, int j3, int j4, int j5, int j6) // for 2*j
{
  double sixj = 0.0;
  double phase0 = phase2(j1 + j2 + j3 + j4 + j5 + j6);
  int m1, m2, m3, m4, m5, m6;
  int min, max;
  double CGC1;
  for(m1 = -j1; m1 <= j1; m1 += 2){
    for(m2 = -j2; m2 <= j2; m2 += 2){
      m3 = m1 + m2;
      CGC1 = CGC3_1(j1,m1,j2,m2,j3,-m3);
      min = std::max(-j4,-j5-m3);
      max = std::min(j4,j5-m3);
      for(m4 = min; m4 <= max; m4 += 2){
	m5 = m3 + m4;
	m6 = m1 - m5;
	sixj += this->CGC3_1(j1,-m1,j5,m5,j6,m6) * this->CGC3_1(j4,m4,j5,-m5,j3,m3) * this->CGC3_1(j4,-m4,j2,-m2,j6,-m6)
	  * phase0 * phase2(2*m1 + m2 + m5) * CGC1;
      }
    }
  }
  return sixj;
}*/


int count1(int N)
{
  return int(2 + 0.5*(std::floor(0.5*(N - 1))*(3 + std::floor(0.5*(N - 1))) + std::floor(0.5*N)*(3 + std::floor(0.5*N))));
}

int count2(int N)
{
  return int(3 + std::floor(0.5*(N - 1))*(3 + std::floor(0.5*(N - 1))) + std::floor(0.5*N)*(2 + std::floor(0.5*N)));
}

void AngMom::Build()
{
  int ind1, ind2, ind3;
  int m3, Jmin1, Jmax1, Jmin2, Jmax2;

  this->jmax = 1;
  if( PAR.ME3 == 1 ){
    for(int i = 0; i < SPB.num_states; ++i){
      for(int j = i; j < SPB.num_states; ++j){
	if( SPB.qnums[i].j + SPB.qnums[j].j > this->jmax ){
	  this->jmax = SPB.qnums[i].j + SPB.qnums[j].j;
	}
	for(int k = j; k < SPB.num_states; ++k){
	  if( SPB.e_nlj3[SPB.qnums[i].nlj] + SPB.e_nlj3[SPB.qnums[j].nlj] + SPB.e_nlj3[SPB.qnums[k].nlj] <= SPB.E3Max){
	    if( SPB.qnums[i].j + SPB.qnums[j].j + SPB.qnums[k].j > this->jmax ){
	      this->jmax = SPB.qnums[i].j + SPB.qnums[j].j + SPB.qnums[k].j;
	    }
	  }
	}
      }
    }
  }
  else{
    for(int i = 0; i < SPB.num_states; ++i){
      if( 2 * SPB.qnums[i].j > this->jmax ){
	this->jmax = 2 * SPB.qnums[i].j;
      }
    }
  }
  
  this->jmax = 21;

  this->jsize = this->jmax + 1;
  std::cout << "maxj = " << this->jmax << std::endl;
  std::cout << "jsize = " << this->jsize << std::endl;

  long long threej_count = 0;
  long long sixj_count = 0;

  //  Setup ThreeJ
  this->ThreeJ_vec = new double**[count1(this->jmax)];  //  sum_j=0^maxj(sum_m=0^j)
  for(int j1 = 0; j1 <= this->jmax; ++j1){  //  j1 <= maxj, 0 <= m1 <= j1
    for(int m1 = j1%2; m1 <= j1; m1 += 2){
      ind1 = count1(j1 - 1) + (m1 - m1%2)/2;

      this->ThreeJ_vec[ind1] = new double*[count2(j1)];  //  sum_j=0^maxj(sum_m=-j^j) = (maxj + 1)^2
      for(int j2 = 0; j2 <= j1; ++j2){  //  j2 <= j1, -j2 <= m2 <= j2
	for(int m2 = -j2; m2 <= j2; m2 += 2){
	  ind2 = count2(j2 - 1) + (m2 + j2)/2;

	  m3 = -m1 - m2;
	  Jmax1 = std::min(j1 + j2, this->jmax);
	  if( std::abs(m3) > Jmax1 ){ continue; }
	  Jmin1 = std::max(std::abs(j1 - j2), std::abs(m3));

	  this->ThreeJ_vec[ind1][ind2] = new double[(Jmax1 - Jmin1)/2 + 1];
          #pragma omp parallel for reduction(+:threej_count)
	  for(int j3 = Jmin1; j3 <= Jmax1; j3 += 2){
	    if( logfac(j1) + logfac(j2) + logfac(j3) < 30.0 ){
	      this->ThreeJ_vec[ind1][ind2][(j3 - Jmin1)/2] = this->CGC_1(j1, m1, j2, m2, j3, -m3) * phase2(j1 - j2 - m3) / std::sqrt(j3 + 1);
	    }
	    else{
	      this->ThreeJ_vec[ind1][ind2][(j3 - Jmin1)/2] = this->CGC_large_1(j1, m1, j2, m2, j3, -m3) * phase2(j1 - j2 - m3) / std::sqrt(j3 + 1);
	    }
	    ++threej_count;
	  }
	}
      }
    }
  }
  std::cout << "ThreeJ count = " << threej_count << std::endl;
  
  //  Setup SixJ
  this->SixJ_vec = new double**[(this->jmax + 1) * (this->jmax + 2) / 2];  //  sum_j1=0^maxj(sum_j4=0^j1) = (maxj + 1)*(maxj + 2)/2
  for(int j1 = 0; j1 <= this->jmax; ++j1){  //  j1 <= maxj, j4 <= j1
    for(int j4 = 0; j4 <= j1; ++j4){
      ind1 = (j1 * (j1 + 1) / 2) + j4;

      this->SixJ_vec[ind1] = new double*[(j1 + 1) * (j1 + 2) / 2];  //  sum_j2=0^j1(sum_j5=0^j2) = (j1 + 1)*(j1 + 2)/2
      for(int j2 = 0; j2 <= j1; ++j2){  //  j2 <= j1, j5 <= j2
	for(int j5 = 0; j5 <= j2; ++j5){
	  ind2 = (j2 * (j2 + 1) / 2) + j5;

	  if( ((j1 + j2)%2 != (j4 + j5)%2) || ((j1 + j5)%2 != (j2 + j4)%2) ){ continue; }
	  Jmin1 = std::max(std::abs(j1 - j2), std::abs(j4 - j5));
	  Jmax1 = std::min(j1 + j2, j4 + j5);
	  Jmax1 = std::min(Jmax1, j2);
	  if( Jmin1 > Jmax1 ){ continue; }
	  Jmin2 = std::max(std::abs(j1 - j5), std::abs(j2 - j4));
	  Jmax2 = std::min(j1 + j5, j2 + j4);
	  Jmax2 = std::min(Jmax2, j2);
	  if( Jmin2 > Jmax2 ){ continue; }

	  this->SixJ_vec[ind1][ind2] = new double[((Jmax1 - Jmin1)/2 + 1) * ((Jmax2 - Jmin2)/2 + 1)];
          #pragma omp parallel for collapse(2) private(ind3) reduction(+:sixj_count)
	  for(int j3 = Jmin1; j3 <= Jmax1; j3 += 2){
	    for(int j6 = Jmin2; j6 <= Jmax2; j6 += 2){
	      ind3 = ((j3 - Jmin1)/2) * ((Jmax2 - Jmin2)/2 + 1) + (j6 - Jmin2)/2;
	      this->SixJ_vec[ind1][ind2][ind3] = this->SixJ_1(j1, j2, j3, j4, j5, j6);
	      std::cout << "SixJ[" << ind1 << "][" << ind2 << "][" << ind3 << "]: (" << j1 << "," << j2 << "," << j3 << "," << j4 << "," << j5 << "," << j6 << ") = " << this->SixJ_vec[ind1][ind2][ind3] << std::endl;
	      ++sixj_count;
	    }
	  }
	}
      }
    }
  }
  std::cout << "SixJ count = " << sixj_count << std::endl;

  std::cout << "!! NineJ = " << get_NineJ(9, 5, 4, 2, 1, 3, 9, 4, 5) << std::endl;
  std::cout << "!! NineJ = " << get_NineJ(9, 8, 9, 5, 5, 10, 4, 3, 1) << std::endl;
}

double AngMom::get_CGC(int j1, int m1, int j2, int m2, int j3, int m3)
{
  return phase2(-j1 + j2 - m3) * std::sqrt(j3 + 1.0) * this->get_ThreeJ(j1, m1, j2, m2, j3, -m3);
}

double AngMom::get_ThreeJ(int j1, int m1, int j2, int m2, int j3, int m3)
{
  if( m1 + m2 + m3 != 0 || std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(m3) > j3 || std::abs(j1 - j2) > j3 || j1 + j2 < j3){ return 0.0; }
  if( std::abs(m1)%2 != j1%2 || std::abs(m2)%2 != j2%2 || std::abs(m3)%2 != j3%2 ){ return 0.0; }
  int ind1, ind2, ind3, Jmin1;
  int jj1 = j1;
  int mm1 = m1;
  int jj2 = j2;
  int mm2 = m2;
  int jj3 = j3;
  int mm3 = m3;
  int J = j1 + j2 + j3;
  double phase = 1.0;

  // need j1 >= j2 >= j3
  if( jj3 > jj2 ){
    std::swap(jj2, jj3);
    std::swap(mm2, mm3);
    phase *= phase2(J);
  }
  if( jj2 > jj1 ){
    std::swap(jj1, jj2);
    std::swap(mm1, mm2);
    phase *= phase2(J);
  }
  if( jj3 > jj2 ){
    std::swap(jj2, jj3);
    std::swap(mm2, mm3);
    phase *= phase2(J);
  }
  // need m1 >= 0
  if( mm1 < 0 ){
    mm1 *= -1;
    mm2 *= -1;
    mm3 *= -1;
    phase *= phase2(J);
  }
  // get inds
  Jmin1 = std::max(std::abs(jj1 - jj2), std::abs(mm3));
  ind1 = count1(jj1 - 1) + (mm1 - mm1%2)/2;
  ind2 = count2(jj2 - 1) + (mm2 + jj2)/2;
  ind3 = (jj3 - Jmin1)/2;

  return phase * this->ThreeJ_vec[ind1][ind2][ind3];
}

double AngMom::get_SixJ(int j1, int j2, int j3, int j4, int j5, int j6)
{
  if( (j3 > j1 + j2) || (j3 > j4 + j5) || (j6 > j1 + j5) || (j6 > j4 + j2) ){ return 0.0; }
  if( (j3 < std::abs(j1 - j2)) || (j3 < std::abs(j4 - j5)) || (j6 < std::abs(j1 - j5)) || (j6 < std::abs(j4 - j2)) ){ return 0.0; }
  if( ((j1 + j2)%2 != (j4 + j5)%2) || ((j1 + j5)%2 != (j2 + j4)%2) ){ return 0.0; }

  int jj1 = j1;
  int jj2 = j2;
  int jj3 = j3;
  int jj4 = j4;
  int jj5 = j5;
  int jj6 = j6;
  int jmax1 = std::max(j1, j4);
  int jmax2 = std::max(j2, j5);
  int jmax3 = std::max(j3, j6);
  int ind1, ind2, ind3, Jmin1, Jmin2, Jmax2;

  // need jmax1 >= jmax2 >= jmax3
  if( jmax3 > jmax2 ){
    std::swap(jj2, jj3);
    std::swap(jj5, jj6);
    std::swap(jmax2, jmax3);
  }
  if( jmax2 > jmax1 ){
    std::swap(jj1, jj2);
    std::swap(jj4, jj5);
    std::swap(jmax1, jmax2);
  }
  if( jmax3 > jmax2 ){
    std::swap(jj2, jj3);
    std::swap(jj5, jj6);
    std::swap(jmax2, jmax3);
  }
  // need j1 >= j4
  if( jj4 > jj1 ){
    std::swap(jj1, jj4);
    std::swap(jj2, jj5);
  }
  // need j2 >= j5
  if( jj5 > jj2 ){
    std::swap(jj2, jj5);
    std::swap(jj3, jj6);
  }

  ind1 = (jj1 * (jj1 + 1) / 2) + jj4;
  ind2 = (jj2 * (jj2 + 1) / 2) + jj5;
  Jmin1 = std::max(std::abs(jj1 - jj2), std::abs(jj4 - jj5));
  Jmin2 = std::max(std::abs(jj1 - jj5), std::abs(jj2 - jj4));
  Jmax2 = std::min(jj1 + jj5, jj2 + jj4);
  Jmax2 = std::min(Jmax2, jj2);
  ind3 = ((jj3 - Jmin1)/2) * ((Jmax2 - Jmin2)/2 + 1) + (jj6 - Jmin2)/2;

  if( j1 == 14 && j2 == 15 && j3 == 3 && j4 == 11 && j5 == 10 && j6 == 4 ){
    std::cout << std::endl;
    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
    std::cout << ind1 << "  " << ind2 << "  " << ind3 << std::endl << std::endl;
    std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
  }

  return this->SixJ_vec[ind1][ind2][ind3];
}

double AngMom::get_NineJ(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9)
{
  double ninej0, ninej;
  int jmin = std::max(std::abs(j1 - j9), std::abs(j2 - j6));
  jmin = std::max(jmin, std::abs(j4 - j8));
  int jmax = std::min(j1 + j9, j2 + j6);
  jmax = std::min(jmax, j4 + j8);

  ninej = 0.0;
  for(int x = jmin; x <= jmax; x += 2){
    ninej0 = (x + 1.0) * phase2(2*x);
    ninej0 *= this->get_SixJ(j1, j4, j7, j8, j9, x);
    ninej0 *= this->get_SixJ(j2, j5, j8, j4, x, j6);
    ninej0 *= this->get_SixJ(j3, j6, j9, x, j1, j2);
    std::cout << "x = " << x << ", fac = " << (x + 1.0) * phase2(2*x) << ", sixj(" << j1 << "," << j4 << "," << j7 << "," << j8 << "," << j9 << "," << x << ") = " << this->get_SixJ(j1, j4, j7, j8, j9, x) << ", sixj(" << j2 << "," << j5 << "," << j8 << "," << j4 << "," << x << "," << j6 << ") = " << this->get_SixJ(j2, j5, j8, j4, x, j6) << ", sixj(" << j3 << "," << j6 << "," << j9 << "," << x << "," << j1 << "," << j2 << ") = " << this->get_SixJ(j3, j6, j9, x, j1, j2) << std::endl;
    ninej += ninej0;
  }
  return ninej;
}

// for 2*j, assume j1 >= j2, abs(m1) <= j1, abs(m2) <= j2, abs(mtot) <= jtot, abs(j1 - j2) <= jtot <= (j1 + j2), mtot = m1 + m2
double AngMom::CGC_1(int j1, int m1, int j2, int m2, int j3, int m3)
{
  double fac1, fac2, fac3, CGC;
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  long long num1, den1, n2_1, n2_2, n2_3;
  double phase = 1.0;
  if( m3 < 0 ){
    m1 *= -1;
    m2 *= -1;
    m3 *= -1;
    phase *= phase2(j1 + j2 - j3);
  }

  num1 = (j3 + 1) * fac((j3 + j1 - j2)/2) * fac((j3 - j1 + j2)/2) * fac((j1 + j2 - j3)/2);
  den1 = fac((j1 + j2 + j3 + 2)/2);
  fac1 = std::sqrt(num1) / std::sqrt(den1);
  n2_1 = fac((j3 + m3)/2) * fac((j3 - m3)/2);
  n2_2 = fac((j1 - m1)/2) * fac((j1 + m1)/2);
  n2_3 = fac((j2 - m2)/2) * fac((j2 + m2)/2);
  fac2 = std::sqrt(n2_1) * std::sqrt(n2_2) * std::sqrt(n2_3);
  maxk = std::min(j1 + j2 - j3, j1 - m1);
  maxk = std::min(maxk, j2 + m2);

  fac3 = 0.0;
  for(int k = 0; k <= maxk; k += 2){
    d3_1 = j1 + j2 - j3 - k;
    d3_2 = j1 - m1 - k;
    d3_3 = j2 + m2 - k;
    d3_4 = j3 - j2 + m1 + k;
    d3_5 = j3 - j1 - m2 + k;
    if (d3_1 >= 0 && d3_2 >= 0 && d3_3 >= 0 && d3_4 >= 0 && d3_5 >= 0){
      fac3 += ((((((phase2(k) / fac(k/2)) / fac(d3_1/2)) / fac(d3_2/2)) / fac(d3_3/2)) / fac(d3_4/2)) / fac(d3_5/2));
    }
  }
  CGC = fac1*fac2*fac3;
  return CGC * phase;
}

double AngMom::CGC_large_1(int j1, int m1, int j2, int m2, int j3, int m3)
{
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  double fac1, fac2, fac3, fac3_1, CGC, num1, den1, n2_1, n2_2, n2_3;
  double phase = 1.0;
  if( m3 < 0 ){
    m1 *= -1;
    m2 *= -1;
    m3 *= -1;
    phase *= phase2(j1 + j2 - j3);
  }

  num1 = std::log(j3 + 1) + logfac((j3 + j1 - j2)/2) + logfac((j3 - j1 + j2)/2) + logfac((j1 + j2 - j3)/2);
  den1 = logfac((j1 + j2 + j3 + 2)/2);
  fac1 = std::exp(0.5 * (num1 - den1));
  n2_1 = logfac((j3 + m3)/2) + logfac((j3 - m3)/2);
  n2_2 = logfac((j1 - m1)/2) + logfac((j1 + m1)/2);
  n2_3 = logfac((j2 - m2)/2) + logfac((j2 + m2)/2);
  fac2 = std::exp(0.5 * (n2_1 + n2_2 + n2_3));
  maxk = std::min(j1 + j2 - j3, j1 - m1);
  maxk = std::min(maxk, j2 + m2);

  fac3 = 0.0;
  for(int k = 0; k <= maxk; k += 2){
    d3_1 = j1 + j2 - j3 - k;
    d3_2 = j1 - m1 - k;
    d3_3 = j2 + m2 - k;
    d3_4 = j3 - j2 + m1 + k;
    d3_5 = j3 - j1 - m2 + k;
    if (d3_1 >= 0 && d3_2 >= 0 && d3_3 >= 0 && d3_4 >= 0 && d3_5 >= 0){
      fac3_1 = -logfac(k/2) - logfac(d3_1/2) - logfac(d3_2/2) - logfac(d3_3/2) - logfac(d3_4/2) - logfac(d3_5/2);
      fac3 += phase2(k) * std::exp(fac3_1);
    }
  }
  CGC = fac1*fac2*fac3;
  return CGC * phase;
}

double AngMom::SixJ_1(int j1, int j2, int j3, int j4, int j5, int j6) // for 2*j
{
  double SixJ = 0.0;
  int J = j1 + j2 + j3 + j4 + j5 + j6;
  int m1, m2, m3, m4, m5, m6;
  int min, max;
  double ThreeJ1;
  for(m1 = -j1; m1 <= j1; m1 += 2){
    for(m2 = -j2; m2 <= j2; m2 += 2){
      m3 = m1 + m2;
      ThreeJ1 = this->get_ThreeJ(j1,m1,j2,m2,j3,-m3);
      min = std::max(-j4,-j5-m3);
      max = std::min(j4,j5-m3);
      for(m4 = min; m4 <= max; m4 += 2){
	m5 = m3 + m4;
	m6 = m1 - m5;
	SixJ += phase2(J - m3 - 2*m5 - m6) * ThreeJ1
	  * this->get_ThreeJ(j1,-m1,j5,m5,j6,m6)
	  * this->get_ThreeJ(j4,m4,j5,-m5,j3,m3)
	  * this->get_ThreeJ(j4,-m4,j2,-m2,j6,-m6);
	if( j1 == 15 && j2 == 14 && j3 == 3 && j4 == 10 && j5 == 11 && j6 == 4){
	  std::cout << "threej1(" << j1 << "," << m1 << "," << j2 << "," << m2 << "," << j3 << "," << -m3 << ") = " << this->get_ThreeJ(j1,m1,j2,m2,j3,-m3) << std::endl;
	  std::cout << "threej2(" << j1 << "," << -m1 << "," << j5 << "," << m5 << "," << j6 << "," << m6 << ") = " << this->get_ThreeJ(j1,-m1,j5,m5,j6,m6) << std::endl;
	  std::cout << "threej3(" << j4 << "," << m4 << "," << j5 << "," << -m5 << "," << j3 << "," << m3 << ") = " << this->get_ThreeJ(j4,m4,j5,-m5,j3,m3) << std::endl;
	  std::cout << "threej4(" << j4 << "," << -m4 << "," << j2 << "," << -m2 << "," << j6 << "," << -m6 << ") = " << this->get_ThreeJ(j4,-m4,j2,-m2,j6,-m6) << std::endl;
	  std::cout << "SixJ = " << SixJ << ", phase = " << phase2(J - m3 - 2*m5 - m6) << std::endl;
	}
      }
    }
  }
  return SixJ;
}

void AngMom::Delete()
{
  int ind1, ind2, m3, Jmin1, Jmax1, Jmin2, Jmax2;
  for(int j1 = 0; j1 <= this->jsize; ++j1){  //  j1 <= maxj, 0 <= m1 <= j1
    for(int m1 = 0; m1 <= j1; ++m1){
      ind1 = (j1 * (j1 + 1) / 2) + m1;

      for(int j2 = 0; j2 <= j1; ++j2){  //  j2 <= j1, -j2 <= m2 <= j2
	for(int m2 = -j2; m1 <= j2; ++m2){
	  ind2 = (j2 * j2) + (m2 + j2);

	  m3 = -m1 - m2;
	  Jmax1 = j1 + j2;
	  if( std::abs(m3) > Jmax1 ){ continue; }
	  if( this->ThreeJ_vec[ind1][ind2] != NULL ){ delete[] this->ThreeJ_vec[ind1][ind2]; }
	}
      }
      if( this->ThreeJ_vec[ind1] != NULL ){ delete[] this->ThreeJ_vec[ind1]; }
    }
  }
  if( this->ThreeJ_vec != NULL ){ delete[] this->ThreeJ_vec; }
  
  this->SixJ_vec = new double**[(this->jmax + 1) * (this->jmax + 2) / 2];  //  sum_j1=0^maxj(sum_j4=0^j1) = (maxj + 1)*(maxj + 2)/2
  for(int j1 = 0; j1 <= this->jsize; ++j1){  //  j1 <= maxj, j4 <= j1
    for(int j4 = 0; j4 <= j1; ++j4){
      ind1 = (j1 * (j1 + 1) / 2) + j4;

      this->SixJ_vec[ind1] = new double*[(j1 + 1) * (j1 + 2) / 2];  //  sum_j2=0^j1(sum_j5=0^j2) = (j1 + 1)*(j1 + 2)/2
      for(int j2 = 0; j2 <= j1; ++j2){  //  j2 <= j1, j5 <= j2
	for(int j5 = 0; j5 <= j2; ++j5){
	  ind2 = (j2 * (j2 + 1) / 2) + j5;

	  Jmin1 = std::max(std::abs(j1 - j2), std::abs(j4 - j5));
	  Jmax1 = std::min(j1 + j2, j4 + j5);
	  if( Jmin1 > Jmax1 ){ continue; }

	  Jmin2 = std::max(std::abs(j1 - j5), std::abs(j2 - j4));
	  Jmax2 = std::min(j1 + j5, j2 + j4);
	  if( Jmin2 > Jmax2 ){ continue; }

	  if( this->SixJ_vec[ind1][ind2] != NULL ){ delete[] this->SixJ_vec[ind1][ind2]; }
	}
      }
      if( this->SixJ_vec[ind1] != NULL ){ delete[] this->SixJ_vec[ind1]; }
    }
  }
  if( this->SixJ_vec != NULL ){ delete[] this->SixJ_vec; }
}
