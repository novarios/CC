#include "MATHfunctions.hpp"
#include "CoupledCluster.hpp"

double phase(int arg)
{
  return std::pow(-1.0, arg);
}

double phase2(int arg)
{
  return phase(arg/2);
}

double d_ij(int i, int j)
{
  if(i == j){ return std::sqrt(2.0); }
  else{ return 1.0; }
}

double delta_ij(int i, int j)
{
  if(i == j){ return 1.0; }
  else{ return 0.0; }
}

long long fac(int n)
{
  if(n < 0){ std::cerr << n << " : Factorial intput should be >= 0" << std::endl; exit(1); }
  long long fac = 1;
  for(int a = 2; a <= n; a++){ fac *= a; }
  return fac;
}

long long fac2(int n)
{
  if(n < 0){ std::cerr << n << " : Factorial2 intput should be >= 0" << std::endl; exit(1); }
  else if(n == 0){ return 1; }
  long long intfactorial;
  if(n % 2 == 0){ intfactorial = 2; for(int a = 4; a <= n; a = a + 2){ intfactorial = intfactorial * a; } }
  else{ intfactorial = 1; for(int a = 3; a <= n; a = a + 2){ intfactorial = intfactorial * a; } }
  return intfactorial;
}

double logfac(int n)
{
  if(n < 0){ std::cerr << n << " : LogFactorial intput should be >= 0" << std::endl; exit(1); }
  double fac = 0.0;
  for(int a = 2; a < n+1; a++){ fac += std::log(a); }
  return fac;
}

double logfac2(int n)
{
  if(n < 0){ std::cerr << n << " : LogFactorial2 intput should be >= 0" << std::endl; exit(1); }
  else if(n == 0){ return 0.0; }
  double fac = 0.0;
  if(n % 2 == 0){ fac = log(2); for(int a = 4; a <= n; a += 2){ fac += log(a); } }
  else{ fac = 0.0; for(int a = 3; a <= n; a += 2){ fac += log(a); } }
  return fac;
}

int choose(int int1, int int2)
{
  return fac(int1) / (fac(int2) * fac(int1 - int2));
}

double logratio1(int int1, int int2, int int3, int int4)
{
  return -logfac(int1) - logfac(int2) - logfac(int3) - logfac(int4);
}

double logratio2(int G)
{
  return -0.5 * (G + 1) * std::log(2);
}

double product1(int n1, int m1, int n2, int m2, int n3, int m3, int n4, int m4)
{
  double prod = logfac(n1) + logfac(n2) + logfac(n3) + logfac(n4);
  prod -= (logfac(n1 + std::abs(m1)) + logfac(n2 + std::abs(m2)) + logfac(n3 + std::abs(m3)) + logfac(n4 + std::abs(m4)));
  prod *= 0.5;
  return std::exp(prod);
}

double logproduct2(int n1, int m1, int n2, int m2, int n3, int m3, int n4, int m4, int j1, int j2, int j3, int j4)
{
  double prod = logfac(n1 + std::abs(m1)) + logfac(n2 + std::abs(m2)) + logfac(n3 + std::abs(m3)) + logfac(n4 + std::abs(m4));
  prod -= (logfac(n1 - j1) + logfac(n2 - j2) + logfac(n3 - j3) + logfac(n4 - j4));
  prod -= (logfac(j1 + std::abs(m1)) + logfac(j2 + std::abs(m2)) + logfac(j3 + std::abs(m3)) + logfac(j4 + std::abs(m4)));
  return prod;
}

double logproduct3(int l1, int l2, int l3, int l4, int g1, int g2, int g3, int g4)
{
  double prod = logfac(g1) + logfac(g2) + logfac(g3) + logfac(g4);
  prod -= (logfac(l1) + logfac(l2) + logfac(l3) + logfac(l4));
  prod -= (logfac(g1 - l1) + logfac(g2 - l2) + logfac(g3 - l3) + logfac(g4 - l4));
  return prod;
}

double loggamma(double x){
  if(x <= 0.0){ std::cerr << x << " : LogGamma intput should be > 0" << std::endl; exit(1); }  
  else if(x == 1.0 || x == 2.0){ return 0.0; }
  double x0, x2, xp, gl0, gl;
  int n = 0;
  if(x <= 7.0){
    n = int(7 - x);
    x0 = x + n;
  }

  x0 = x;
  x2 = 1.0 / (x0 * x0);
  xp = 2.0 * PI;
  gl0 = -1.39243221690590;
  gl0 = gl0 * x2 + 1.796443723688307e-01;
  gl0 = gl0 * x2 - 2.955065359477124e-02;
  gl0 = gl0 * x2 + 6.410256410256410e-03;
  gl0 = gl0 * x2 - 1.917526917526918e-03;
  gl0 = gl0 * x2 + 8.417508417508418e-04;
  gl0 = gl0 * x2 - 5.952380952380952e-04;
  gl0 = gl0 * x2 + 7.936507936507937e-04;
  gl0 = gl0 * x2 - 2.777777777777778e-03;
  gl0 = gl0 * x2 + 8.333333333333333e-02;
  gl = gl0/x0 + 0.5*log(xp) + (x0 - 0.5)*log(x0) - x0;
  if(x <= 7.0){
    for(int i = 1; i <= n; ++i){
      gl -= log(x0 - 1.0);
      x0 -= 1.0;
    }
  }
  return gl;
}

double CGC(int j1, int m1, int j2, int m2, int jtot, int mtot) // for 2*j
{
  if(m1 + m2 - mtot != 0 || std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(mtot) > jtot || std::abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  double fac1, fac2, fac3, num1, den1, n2_1, n2_2, n2_3, CGC;
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  int change1 = 0, change2 = 0; // flags to change from general formula

  if(j1 < j2){ std::swap(j1,j2); std::swap(m1,m2); change1 = 1; };  
  if(mtot < 0){ m1 = -m1; m2 = -m2; mtot = std::abs(mtot); change2 = 1; }

  num1 = (jtot + 1.0) * fac((jtot + j1 - j2)/2) * fac((jtot - j1 + j2)/2) * fac((j1 + j2 - jtot)/2);
  den1 = fac((j1 + j2 + jtot + 2)/2);
  fac1 = std::sqrt(num1 / den1);
  n2_1 = fac((jtot + mtot)/2) * fac((jtot - mtot)/2);
  n2_2 = fac((j1 - m1)/2) * fac((j1 + m1)/2);
  n2_3 = fac((j2 - m2)/2) * fac((j2 + m2)/2);
  fac2 = std::sqrt(n2_1 * n2_2 * n2_3);
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
      fac3 += phase2(k) / (fac(k/2) * fac(d3_1/2) * fac(d3_2/2) * fac(d3_3/2) * fac(d3_4/2) * fac(d3_5/2));
    }
  }
  CGC = fac1*fac2*fac3;
  if (change1 == 1){ CGC *= phase2(jtot - j2 - j1); };
  if (change2 == 1){ CGC *= phase2(jtot - j1 - j2); };
  return CGC;
}

double CGC3(int j1, int m1, int j2, int m2, int jtot, int mtot) // for 2*j
{
  if(m1 + m2 + mtot != 0 || abs(m1) > j1 || abs(m2) > j2 || abs(mtot) > jtot || abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  else{ return phase2(j1 - j2 - mtot) * CGC(j1,m1,j2,m2,jtot,-mtot) / std::sqrt(jtot + 1.0); }
}

double CGC6(int j1, int j2, int j3, int j4, int j5, int j6) // for 2*j
{
  double sixj = 0.0;
  double phase0 = phase2(j1 + j2 + j3 + j4 + j5 + j6);
  int m1, m2, m3, m4, m5, m6;
  int min, max;
  double CGC1;
  for(m1 = -j1; m1 <= j1; m1 += 2){
    for(m2 = -j2; m2 <= j2; m2 += 2){
      m3 = m1 + m2;
      CGC1 = CGC3(j1,m1,j2,m2,j3,-m3);
      min = std::max(-j4,-j5-m3);
      max = std::min(j4,j5-m3);
      for(m4 = min; m4 <= max; m4 += 2){
	m5 = m3 + m4;
	m6 = m1 - m5;
	sixj += phase2(2*m1 + m2 + m5) * CGC1 * CGC3(j1,-m1,j5,m5,j6,m6) * CGC3(j4,m4,j5,-m5,j3,m3) * CGC3(j4,-m4,j2,-m2,j6,-m6);
	if( j1 == 3 && j2 == 2 && j3 == 1 && j4 == 0 && j5 == 1 && j6 == 2){
	  std::cout << "threej1(" << j1 << "," << m1 << "," << j2 << "," << m2 << "," << j3 << "," << -m3 << ") = " << CGC3(j1,m1,j2,m2,j3,-m3) << std::endl;
	  std::cout << "threej2(" << j1 << "," << -m1 << "," << j5 << "," << m5 << "," << j6 << "," << m6 << ") = " << CGC3(j1,-m1,j5,m5,j6,m6) << std::endl;
	  std::cout << "threej3(" << j4 << "," << m4 << "," << j5 << "," << -m5 << "," << j3 << "," << m3 << ") = " << CGC3(j4,m4,j5,-m5,j3,m3) << std::endl;
	  std::cout << "threej4(" << j4 << "," << -m4 << "," << j2 << "," << -m2 << "," << j6 << "," << -m6 << ") = " << CGC3(j4,-m4,j2,-m2,j6,-m6) << std::endl;
	  std::cout << "SixJ = " << sixj << ", phase = " << phase0 * phase2(2*m1 + m2 + m5) << std::endl;
	}
      }
    }
  }
  return phase0 * sixj;
}

/*double CGC(int j1, int m1, int j2, int m2, int jtot, int mtot) // for 2*j
{
  if(m1 + m2 - mtot != 0 || std::abs(m1) > j1 || std::abs(m2) > j2 || std::abs(mtot) > jtot || std::abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  double fac1, fac2, fac3, fac30, num1, den1, n2_1, n2_2, n2_3, CGC;
  int maxk, d3_1, d3_2, d3_3, d3_4, d3_5;
  int change1 = 0, change2 = 0; // flags to change from general formula

  if(j1 < j2){ std::swap(j1,j2); std::swap(m1,m2); change1 = 1; };  
  if(mtot < 0){ m1 = -m1; m2 = -m2; mtot = std::abs(mtot); change2 = 1; }

  num1 = std::log(jtot + 1.0) + logfac((jtot + j1 - j2)/2) + logfac((jtot - j1 + j2)/2) + logfac((j1 + j2 - jtot)/2);
  den1 = logfac((j1 + j2 + jtot + 2)/2);
  fac1 = 0.5 * (num1 - den1);
  n2_1 = logfac((jtot + mtot)/2) + logfac((jtot - mtot)/2);
  n2_2 = logfac((j1 - m1)/2) + logfac((j1 + m1)/2);
  n2_3 = logfac((j2 - m2)/2) + logfac((j2 + m2)/2);
  fac2 = 0.5 * (n2_1 + n2_2 + n2_3);
  maxk = std::min(j1 + j2 - jtot, j1 - m1);
  maxk = std::min(maxk, j2 + m2);

  fac3 = 0.0;
  for(int k = 0; k <= maxk; k += 2){
    d3_1 = j1 + j2 - jtot - k;
    d3_2 = j1 - m1 - k;
    d3_3 = j2 + m2 - k;
    d3_4 = jtot - j2 + m1 + k;
    d3_5 = jtot - j1 - m2 + k;
    if( d3_1 >= 0 && d3_2 >= 0 && d3_3 >= 0 && d3_4 >= 0 && d3_5 >= 0 ){
      fac30 = -logfac(k/2) - logfac(d3_1/2) - logfac(d3_2/2) - logfac(d3_3/2) - logfac(d3_4/2) -logfac(d3_5/2);
      fac3 += phase2(k) * std::pow(e, fac30);
    }
  }
  CGC = std::pow(e, fac1 + fac2) * fac3;
  if( change1 == 1 ){ CGC *= phase2(jtot - j2 - j1); };
  if( change2 == 1 ){ CGC *= phase2(jtot - j1 - j2); };
  return CGC;
}

double CGC3(int j1, int m1, int j2, int m2, int jtot, int mtot) // for 2*j
{
  if(m1 + m2 + mtot != 0 || abs(m1) > j1 || abs(m2) > j2 || abs(mtot) > jtot || abs(j1 - j2) > jtot || j1 + j2 < jtot){ return 0.0; }
  else{ return phase2(j1 - j2 - mtot) * CGC(j1,m1,j2,m2,jtot,-mtot) / std::sqrt(jtot + 1.0); }
}

double CGC6(int j1, int j2, int j3, int j4, int j5, int j6) // for 2*j
{
  double sixj = 0.0;
  int S, m1, m2, m3, m4, m5, m6;
  for(m1 = -j1; m1 <= j1; m1 += 2){
    for(m2 = -j2; m2 <= j2; m2 += 2){
      for(m4 = -j4; m4 <= j4; m4 += 2){
	for(m5 = -j5; m5 <= j5; m5 += 2){
	  m3 = m2 + m1;
	  m6 = m1 - m5;
	  if(m3 == m5 - m4 && m6 == -m2 - m4){
       	    S = j1 - m1 + j2 - m2 + j3 - m3 + j4 - m4 + j5 - m5 + j6 - m6;
	    sixj += phase2(S) * CGC3(j1,m1,j2,m2,j3,-m3) * CGC3(j1,-m1,j5,m5,j6,m6) * CGC3(j4,m4,j5,-m5,j3,m3) * CGC3(j4,-m4,j2,-m2,j6,-m6);
	  }
	}
      }
    }
  }
  return sixj;
  }*/

double CGC9_0(int j1, int j2, int j3, int j4, int j5, int j6)
{
  return phase2(j2 + j5 + j3 + j6) * CGC6(j1, j2, j5, j4, j3, j6);
}

double Pandya(int j1, int j2, int j3, int j4, int j13)
{
  double p = 0.0;
  int j12min = int(abs(j1 - j2));
  int j12max = j1 + j2;
  for(int j = j12min; j <= j12max; j += 2){ p += phase2(j3 + j4 + j) * (j + 1.0) * CGC6(j1, j2, j, j4, j3, j13); }
  return p;
}

double Chi_J(int a, int b, int c, int d, int J1, int J2)
{
  double fac1 = std::sqrt((J1 + 1.0) * (J2 + 1.0));
  return fac1 * phase2(b + J1 + c + J2) * CGC6(a, b, J1, d, c, J2);
}

double Legendre(double x, int l, int m) // l, m x2
{
  if(std::abs(m) > l){ return 0.0; }
  int M = abs(m);
  int L = M;
  double P;
  if(M == 0){ P = 1.0; }
  else{ P = std::pow(-1.0, M/2) * fac2(M - 1) * std::pow(1.0 - x*x, 0.25*M); } // P_M^M
  double P1 = 0.0;
  double tempP; //P_L-1^M, reset P1
  while(L != l){
    tempP = P;
    if(L == M){ P = x * (L + 1) * P; }
    else{ P = ((L + 1) * x * P - ((L + M)/2) * P1)/((L - M)/2 + 1); }
    P1 = tempP;
    L += 2;
  }
  if(M != m){ P *= std::pow(-1.0, M/2) * fac((l - M)/2) / fac((l + M)/2); }
  return P;
}

std::complex<double> SphericalY_C(double theta, double phi, int l, int m) // l, m x2
{
  if(std::abs(m) > l){ return 0.0; }
  std::complex<double> I (0.0, 1.0), M (m, 0.0), PHI (phi, 0.0);
  std::complex<double> fac1 (std::pow(-1.0, m/2), 0.0), fac2 (std::sqrt(((l + 1) * fac((l - m)/2))/(4.0 * PI * fac((l + m)/2))), 0.0);
  std::complex<double> fac3 (Legendre(std::cos(theta), l, m), 0.0);
  return fac1 * fac2 * fac3 * pow(e, I * M * PHI);
}

double SphericalY(double theta, double phi, int l, int m) // l, m x2
{
  if(std::abs(m) > l){ return 0.0; }
  if(m < 0){ return std::sqrt(2.0) * std::pow(-1.0, m/2) * SphericalY_C(theta, phi, l, std::abs(m)).imag(); }
  else if(m == 0){ return SphericalY_C(theta, phi, l, 0).real(); }
  else{ return std::sqrt(2.0) * std::pow(-1.0, m/2) * SphericalY_C(theta, phi, l, m).real(); }
}

double SphericalYTens(double theta, double phi, int j, int l, int s, int ms) // j, l, s, ms x2
{
  if(std::abs(ms) > s){ return 0.0; }
  double Y = 0.0;
  for(int ml = -l; ml <= l; ml += 2){ Y += SphericalY(theta, phi, l/2, ml/2) * CGC(l, ml, s, ms, j, ml + ms); }
  return Y;
}

double Erf(double z)
{
  double t = 1.0/(1 + 0.5*fabs(z));
  double tau = t * std::pow(e, -1.0*z*z - 1.26551223 + 1.00002368*t + 0.37409196*t*t + 0.09678418*t*t*t
		       - 0.18628806*t*t*t*t + 0.27886807*t*t*t*t*t - 1.13520398*t*t*t*t*t*t
		       + 1.48851587*t*t*t*t*t*t*t - 0.82215223*t*t*t*t*t*t*t*t + 0.17087277*t*t*t*t*t*t*t*t*t);
  if(z >= 0.0){ return 1.0 - tau; }
  else{ return tau - 1.0; }
}

void projection(double *u, double *v, double *proj, int size)
{
  double innerprod = 0.0;
  double norm = 0.0;
  for(int i = 0; i < size; ++i){ innerprod += v[i]*u[i]; norm += u[i]*u[i]; }
  for(int i = 0; i < size; ++i){ proj[i] = -1.0 * (innerprod/norm)*u[i]; }
}

void GramSchmidt(double *Vectors, int size)
{
  if(size == 0){ return; }
  double norm;
  double innerprod;
  double *V = new double[size*size];
  for(int i = 0; i < size; ++i){
    for(int j = 0; j < size; ++j){ V[size * i + j] = Vectors[size * i + j]; }
  }
  norm = 0.0;
  for(int i = 0; i < size; ++i){ norm += V[i] * V[i]; }
  for(int i = 0; i < size; ++i){ Vectors[i] = V[i]/std::sqrt(norm); }
  for(int i = 1; i < size; ++i){
    for(int j = 0; j < size; ++j){ Vectors[size * i + j] = V[size * i + j]; }
    for(int j = 0; j < i; ++j){
      innerprod = 0.0;
      for(int k = 0; k < size; ++k){ innerprod += Vectors[size * i + k] * Vectors[size * j + k]; }
      for(int k = 0; k < size; ++k){ Vectors[size * i + k] -= innerprod * Vectors[size * j + k]; }
    }
    norm = 0.0;
    for(int j = 0; j < size; ++j){ norm += Vectors[size * i + j] * Vectors[size * i + j]; }
    for(int j = 0; j < size; ++j){ Vectors[size * i + j] /= std::sqrt(norm); }
  }
  delete[] V;
}

//box muller method
double rand_normal(double mean, double stddev)
{
  static double n2 = 0.0;
  static int n2_cached = 0;
  if(!n2_cached){
    double x, y, r;
    do{
      x = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
      y = 2.0 * double(rand())/double(RAND_MAX) - 1.0;
      r = x*x + y*y;
    }
    while(r == 0.0 || r > 1.0);
    double d = std::sqrt(-2.0*log(r)/r);
    double n1 = x*d;
    n2 = y*d;
    double result = n1*stddev + mean;
    n2_cached = 1;
    return result;
  }
  else{
    n2_cached = 0;
    return n2*stddev + mean;
  }
}

void Asym_Diagonalize1(double *Ham, int &N, double *eigenvalues, double *eigenvectors_L, double *eigenvectors_R, int num)
{
  if(num > N){ std::cerr << "Asym_Diagonalize1: number of vectors > size" << std::endl; exit(1); }
  int info = 0;
  int ind = -1, ind0 = -1;
  char job1 = 'V';
  char job2 = 'V';
  int lwork = 10*N;
  double *Vl = new double[N * N];
  double *Vr = new double[N * N];
  double *wr = new double[N];
  double *wi = new double[N];
  double *work = new double[lwork];
  double tempen, templ, tempr;
  dgeev_(&job1, &job2, &N, Ham, &N, wr, wi, Vl, &N, Vr, &N, work, &lwork, &info);
  
  ind0 = 0;
  while(ind0 < num){
    tempen = 1.0e10;
    for(int i = ind0; i < N; ++i){
      if(wr[i] < tempen){
	tempen = wr[i];
	ind = i;
      }
    }
    wr[ind] = wr[ind0];
    wr[ind0] = tempen;
    eigenvalues[ind0] = tempen;
    for(int i = 0; i < N; ++i){
      tempr = Vr[N*ind + i];
      Vr[N*ind + i] = Vr[N*ind0 + i];
      Vr[N*ind0 + i] = tempr;
      eigenvectors_R[N*ind0 + i] = tempr;
      templ = Vl[N*ind + i];
      Vl[N*ind + i] = Vl[N*ind0 + i];
      Vl[N*ind0 + i] = templ;
      eigenvectors_L[N*ind0 + i] = templ;
    }
    ++ind0;
  }

  // Biorthogonalize
  double norm;
  for(int i = 0; i < num; ++i){
    norm = 0.0;
    for(int j = 0; j < N; ++j){ norm += eigenvectors_L[N*i + j] * eigenvectors_R[N*i + j]; }
    if( norm < 0.0 ){
      norm *= -1.0;
      for(int j = 0; j < N; ++j){ eigenvectors_R[N*i + j] *= -1.0; }
    }
    norm = std::sqrt(norm);
    for(int j = 0; j < N; ++j){
      eigenvectors_L[N*i + j] /= norm;
      eigenvectors_R[N*i + j] /= norm;
    }
  }

  delete[] Vl;
  delete[] Vr;
  delete[] wr;
  delete[] wi;
  delete[] work;
}

void Asym_Diagonalize2_0(double *Ham, int &N, double *eigenvalues, double *eigenvectors, int num, char type)
{
  bool rvec = true;
  char howmny = 'A';
  char which[] = "SR", bmat[] = "I"; // standard eigenvalue problem
  int mxiter = 5000;
  int nev = num, ncv = 3*nev + 2; // number of eigenvalues, lanczos vectors to calculate
  if(ncv > N){ ncv = N; }
  int ldv = N, ldz = N, lworkl = 4*ncv*(ncv + 2);
  int mode = 1, ishift = 1, info = 0, ido = 0; // status integer is zero at start
  double sigmar, sigmai, tol = 1.0e-12; // error tolerance
  int *select = new int[ncv];
  double *resid = new double[N];
  double *workev = new double[3 * ncv];
  double *workd = new double[3*N];
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
    if(type == 'R'){
      #pragma omp parallel for
      for(int j = 0; j < N; ++j){
	workd[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k]; }
      }
    }
    else if(type == 'L'){
      #pragma omp parallel for
      for(int j = 0; j < N; ++j){
	workd[ipntr[1]-1 + j] = 0.0;
	for(int k = 0; k < N; ++k){ workd[ipntr[1]-1 + j] += Ham[N*j + k] * workd[ipntr[0]-1 + k]; }
      }
    }
  } while(true);
  dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  
  for(int i = 0; i < num; ++i){
    if(type == 'R'){
      std::cout << "%%R  " << dr[i] << std::endl;
      //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[N*i + j] << " "; }
      //std::cout << std::endl;
      eigenvalues[i] = dr[i];
    }
    else{
      std::cout << "%%L  " << dr[i] << std::endl;
      //for(int j = 0; j < N; ++j){ std::cout << std::setprecision(5) << z[N*i + j] << " "; }
      //std::cout << std::endl;
    }
    for(int j = 0; j < N; ++j){ eigenvectors[N*i + j] = z[N*i + j]; }
  }
  delete[] select;
  delete[] resid;
  delete[] workev;
  delete[] workd;
  delete[] workl;
  delete[] v;
  delete[] dr;
  delete[] di;
  delete[] z;
}

void Asym_Diagonalize2(double *Ham, int &N, double *eigenvalues, double *eigenvectors_L, double *eigenvectors_R, int num)
{
  if(num > N){ std::cerr << "Asym_Diagonalize1: number of vectors > size" << std::endl; exit(1); }
  Asym_Diagonalize2_0(Ham, N, eigenvalues, eigenvectors_L, num, 'L');
  Asym_Diagonalize2_0(Ham, N, eigenvalues, eigenvectors_R, num, 'R');
 
  // Biorthogonalize
  double norm;
  for(int i = 0; i < num; ++i){
    norm = 0.0;
    for(int j = 0; j < N; ++j){ norm += eigenvectors_L[N*i + j] * eigenvectors_R[N*i + j]; }
    if( norm < 0.0 ){
      norm *= -1.0;
      for(int j = 0; j < N; ++j){ eigenvectors_R[N*i + j] *= -1.0; }
    }
    norm = std::sqrt(norm);
    for(int j = 0; j < N; ++j){
      eigenvectors_L[N*i + j] /= norm;
      eigenvectors_R[N*i + j] /= norm;
      if( std::fabs(eigenvectors_L[N*i + j]) < 1.0e-10 ){ eigenvectors_L[N*i + j] = 0.0; }
      if( std::fabs(eigenvectors_R[N*i + j]) < 1.0e-10 ){ eigenvectors_R[N*i + j] = 0.0; }
    }
  }
}

void Asym_Diagonalize2_2(double *Ham, int &N, double *eigenvalues, double *eigenvectors_L, double *eigenvectors_R, int num)
{
  if(num > N){ std::cerr << "Asym_Diagonalize1: number of vectors > size" << std::endl; exit(1); }
  Asym_Diagonalize2_0(Ham, N, eigenvalues, eigenvectors_L, num, 'L');
  Asym_Diagonalize2_0(Ham, N, eigenvalues, eigenvectors_R, num, 'R');
 
  // Biorthogonalize
  double norm;
  for(int i = 0; i < num; ++i){
    norm = 0.0;
    for(int j = 0; j < N; ++j){ norm += eigenvectors_L[N*i + j] * eigenvectors_R[N*i + j]; }
    if( norm < 0.0 ){
      norm *= -1.0;
      for(int j = 0; j < N; ++j){ eigenvectors_R[N*i + j] *= -1.0; }
    }
    norm = std::sqrt(norm);
    for(int j = 0; j < N; ++j){
      eigenvectors_L[N*i + j] /= norm;
      eigenvectors_R[N*i + j] /= norm;
      if( std::fabs(eigenvectors_L[N*i + j]) < 1.0e-10 ){ eigenvectors_L[N*i + j] = 0.0; }
      if( std::fabs(eigenvectors_R[N*i + j]) < 1.0e-10 ){ eigenvectors_R[N*i + j] = 0.0; }
    }
  }
}

/*void Arnoldi_1(double *Ham, int &N, double *eigenvalues, double *eigenvectors, int num, char type)
{
  int ido = 0; // status flag
  char which[] = "SR"; // target smallest real part
  char bmat[] = "I"; // standard eigenvalue problem
  int nev = 1; // calculate lowest eigenpair
  double tol = 1.0e-10; // error tolerance
  double *resid = new double[N];
  int ncv = 3*nev + 2; // number of eigenvalues, lanczos vectors to calculate
  double *v = new double[N*ncv];
  int ldv = N;
  int mode = 1; // standard eigen-problem A*x = lambda*x
  int ishift = 1; // exact implicit shifts with respect to the current Hessenberg matrix
  int info = 0; // status integer is zero at start
  int iparam[11], ipntr[14];
  iparam[0] = ishift;
  iparam[2] = mxiter;
  iparam[6] = mode;
  int lworkl = 4*ncv*(ncv + 2);
  double *workd = new double[3*N];
  double *workl = new double[lworkl];

  dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
  if(ido != -1 && ido != 1){ break; }
  for(int j = 0; j < N; ++j){ workd[ipntr[1]-1 + j] = vec[j]; }
  }*/

//  input: x1   : lower limit of the integration interval                      
//         x2   : upper limit ---------- "" -------------                      
//         n    : the desired number of mesh points                            
//  ouput: x    : gauss-legendre mesh points on the interval (x1,x2)          
//         w    : the corresponding weights                                   
//  From : Numerical recipes, F90 version : M. Hjorth-Jensen
void gauss_legendre(double x1, double x2, double *x, double *w, int n)
{
  int m;
  double eps = 3.0e-14;
  double pp = 0.0;
  double p1, p2, p3, xl, xm, z, z1;
  m = (n+1)/ 2 ;
  xm = 0.5 * (x2+x1);
  xl = 0.5 * (x2-x1);
  for(int i = 0; i < m; ++i){
    z1 = 0.0;
    //z = std::cos(PI*(i - 0.25)/(n + 0.5));
    z = std::cos(PI*(i + 0.75)/(n + 0.5));
    while( std::fabs(z-z1) > eps){
      p1 = 1.0;
      p2 = 0.0;
      for(int j = 1; j <= n; ++j){
	p3 = p2;
	p2 = p1;
	p1 = double(((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3) / j);
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      z1 = z;
      z = z - p1/pp;
    }
    x[i] = xm - xl*z;
    x[n - 1 - i] = xm + xl*z;
    w[i] = 2.0 * xl / ((1.0 - z*z)*pp*pp);
    w[n - 1 - i] = w[i];
  }
}

// L_n^alpha(x),  cx = new double[n + 1]
void laguerre_general(int n, double alpha, double x, double *cx)
{
  if( alpha <= -1.0 ){ std::cerr << "LAGUERRE GENERAL ERROR, ALPHA <= -1" << std::endl; exit(1); }
  if( n < 0 ){ return; }
  cx[0] = 1.0;
  if( n == 0 ){ return; }
  cx[1] = 1.0 + alpha - x;
  for(int i = 2; i <= n; ++i){ cx[i] = double( ( (2*i - 1 + alpha - x )*cx[i-1] + (1 - i - alpha)*cx[i-2] ) / i); }
}

/*void setup_ho_cutoff(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME)
{
  int number_of_iterations, int_points;
  double sigma, norm, oscl_r, sum_hf;
  double z_rel, factor, xp, contrib;
  double *q_points, *q_weights;
  int lmax = Space.qmaxs.l;
  int nmax = Space.qmaxs.n;
  double *sum_norm = new double[lmax + 1];
  double *cx = new double[nmax + 1];

  sigma = 1.0;
  norm = 1.0;
  int_points = 10;
  number_of_iterations = 0;
  oscl_r = Parameters.ho_length * std::sqrt(2.0);
  //  for momentum space if we set z= 0.5*(oscl_r*k)**2 > 60-70, then EXP(-60) < 10E-27.
  //  Making it twice as large ensures that we account for extensions due to the Laguerre polynoms which depend on z**(2n).

  ME.cutoff = 2.0 * std::sqrt(60.0)/Parameters.ho_length;
  while( number_of_iterations < 20 && std::fabs(sigma) > 1.0e-4 ){
    q_points = new double[int_points];
    q_weights = new double[int_points];
    gauss_legendre(0.0, ME.cutoff, q_points, q_weights, int_points);
    for(int h = 0; h <= lmax; ++h){
      sum_hf = 0.0;
      factor = 0.5 * ((nmax + h + 2) * std::log(2.0) + logfac(nmax) - logfac2(2*nmax + 2*h + 1) - 0.5 * std::log(PI));
      factor = std::pow(e, factor);
      for(int iq = 0; iq < int_points; ++iq){
	z_rel = q_points[iq] * q_points[iq] * oscl_r * oscl_r;
	laguerre_general( nmax, h + 0.5, z_rel, cx );
	xp = std::pow(e, -0.5*z_rel) * std::pow(q_points[iq]*oscl_r, h) * cx[nmax];
	contrib = xp * factor * std::pow(oscl_r, 1.5); 
	sum_hf += q_weights[iq] * std::pow(contrib * q_points[iq], 2);
      }
      sum_norm[h] = sum_hf;
    }

    sigma = 0.0;
    for(int h = 0; h <= lmax; ++h){ sigma += std::fabs(sum_norm[h] - norm); }
    sigma /= (lmax + 1);

    ++number_of_iterations;

    delete[] q_points;
    delete[] q_weights;
    if( std::fabs(sigma) > 1.0e-4 ){ int_points += 10; }
  }
  delete[] sum_norm;
  delete[] cx;

  ME.n_rel = int_points;
  //std::cout << "???        points = " << ME.n_rel << ", cutoff = " << ME.cutoff << std::endl;
}

void ho_wfunction(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME)
{
  double ph, factor, z_lab, xp;
  double sum_rel;
  int lmax = Space.qmaxs.l;
  int nmax = Space.qmaxs.n;
  double oscl = Parameters.ho_length;
  double *cx = new double[nmax + 1];

  for(int n = 0; n <= nmax; ++n){
    ph = std::pow(-1.0, n);
    for(int l = 0; l <= lmax; ++l){
      factor = 0.5 * ((n + l + 2) * std::log(2.0) + logfac(n) - logfac2(2*n + 2*l + 1) - 0.5 * std::log(PI));
      factor = std::pow(e, factor);
      sum_rel = 0.0;
      for(int i = 0; i < ME.n_rel; ++i){
	z_lab = ME.ra[i] * ME.ra[i] * oscl * oscl;
        laguerre_general( n, l + 0.5, z_lab, cx );
	xp = std::pow(e, -0.5 * z_lab) * std::pow(ME.ra[i] * oscl, l) * cx[n];
	ME.hol[(nmax+1)*(lmax+1)*i + (nmax+1)*l + n] = xp * ph * factor * std::pow(oscl, 1.5);
	sum_rel += ME.wra[i] * std::pow(ME.hol[(nmax+1)*(lmax+1)*i + (nmax+1)*l + n] * ME.ra[i], 2);
      }
      //std::cout << "HO:  n = " << n << ", l = " << l << ", norm = " << sum_rel << std::endl;
    }
  }
  delete[] cx;
}

double kinetic_energy(Input_Parameters &Parameters, Model_Space &Space, HF_Matrix_Elements &ME, int a, int c)
{
  double e_kin;
  int lmax = Space.qmaxs.l;
  int nmax = Space.qmaxs.n;
  double *int_factor = new double[ME.n_rel];
  double *kin_energy = new double[ME.n_rel];
  double mass;
  if( Space.qnums[a].t == -1 ){ mass = m_neutronc2; }
  else if( Space.qnums[a].t == 1 ){ mass = m_protonc2; }
  for(int i = 0; i < ME.n_rel; ++i){
    int_factor[i] = ME.wra[i] * ME.ra[i] * ME.ra[i];
    kin_energy[i] = 0.5 * ME.ra[i] * ME.ra[i] * hbarc_MeVfm * hbarc_MeVfm / mass;
  }

  e_kin = 0.0;
  for(int k = 0; k < ME.n_rel; ++k){
    e_kin += 2 * ME.hol[(nmax+1)*(lmax+1)*k + (nmax+1)*Space.qnums[a].l + Space.qnums[a].n]
      * ME.hol[(nmax+1)*(lmax+1)*k + (nmax+1)*Space.qnums[c].l + Space.qnums[c].n] * int_factor[k] * kin_energy[k];
  }
  delete[] int_factor;
  delete[] kin_energy;
  
  return e_kin;
  }*/
