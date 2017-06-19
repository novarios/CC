#include "MATHfunctions.hpp"

long long factorial(const int &n)
{
  if(n < 0){ std::cerr << n << " : Factorial intput should be >= 0" << std::endl; exit(1); }
  long long intfactorial = 1;
  for(int a = 2; a <= n; a++){ intfactorial = intfactorial * a; }
  return intfactorial;
}

long long factorial(const double &n)
{
  if(n + 0.2 < 0){ std::cerr << n << " : Factorial intput should be >= 0" << std::endl; exit(1); }
  else if(fabs(n) < 0.1){ return 1; }
  long long intfactorial = 1;
  for(int a = 2; a <= n + 0.1; a++){ intfactorial = intfactorial * a; }
  return intfactorial;
}

double logfac(const int &n)
{
  if(n < 0){ std::cerr << n << " : LogFactorial intput should be >= 0" << std::endl; exit(1); }
  double fac = 0.0;
  for(int a = 2; a < n+1; a++){ fac += log(a); }
  return fac;
}

long long factorial2(const int &n)
{
  if(n < 0){ std::cerr << n << " : Factorial2 intput should be >= 0" << std::endl; exit(1); }
  else if(n == 0){ return 1; }
  long long intfactorial;
  if(n % 2 == 0){ intfactorial = 2; for(int a = 4; a <= n; a = a + 2){ intfactorial = intfactorial * a; } }
  else{ intfactorial = 1; for(int a = 3; a <= n; a = a + 2){ intfactorial = intfactorial * a; } }
  return intfactorial;
}

long long factorial2(const double &n)
{
  if(n + 0.2 < 0){ std::cerr << n << " : Factorial2 intput should be >= 0" << std::endl; exit(1); }
  else if(fabs(n) < 0.1){ return 1; }
  long long intfactorial;
  if(int(n + 0.1) % 2 == 0){ intfactorial = 2; for(int a = 4; a <= n + 0.1; a = a + 2){ intfactorial = intfactorial * a; } }
  else{ intfactorial = 1; for(int a = 3; a <= n + 0.1; a = a + 2){ intfactorial = intfactorial * a; } }
  return intfactorial;
}

int choose(const int &int1, const int &int2)
{
  return factorial(int1) / (factorial(int2) * factorial(int1 - int2));
}

double logratio1(const int &int1, const int &int2, const int &int3, const int &int4)
{
  return -logfac(int1) - logfac(int2) - logfac(int3) - logfac(int4);
}

double logratio2(const int &G)
{
  return -0.5 * (G + 1) * log(2);
}

double product1(const int &n1, const int &m1, const int &n2, const int &m2, const int &n3, const int &m3, const int &n4, const int &m4)
{
  double prod = logfac(n1) + logfac(n2) + logfac(n3) + logfac(n4);
  prod -= (logfac(n1 + abs(m1)) + logfac(n2 + abs(m2)) + logfac(n3 + abs(m3)) + logfac(n4 + abs(m4)));
  prod *= 0.5;
  return exp(prod);
}

double logproduct2(const int &n1, const int &m1, const int &n2, const int &m2, const int &n3, const int &m3, const int &n4, const int &m4, const int &j1, const int &j2, const int &j3, const int &j4)
{
  double prod = logfac(n1 + abs(m1)) + logfac(n2 + abs(m2)) + logfac(n3 + abs(m3)) + logfac(n4 + abs(m4));
  prod -= (logfac(n1 - j1) + logfac(n2 - j2) + logfac(n3 - j3) + logfac(n4 - j4));
  prod -= (logfac(j1 + abs(m1)) + logfac(j2 + abs(m2)) + logfac(j3 + abs(m3)) + logfac(j4 + abs(m4)));
  return prod;
}

double logproduct3(const int &l1, const int &l2, const int &l3, const int &l4, const int &g1, const int &g2, const int &g3, const int &g4)
{
  double prod = logfac(g1) + logfac(g2) + logfac(g3) + logfac(g4);
  prod -= (logfac(l1) + logfac(l2) + logfac(l3) + logfac(l4));
  prod -= (logfac(g1 - l1) + logfac(g2 - l2) + logfac(g3 - l3) + logfac(g4 - l4));
  return prod;
}

double loggamma(const double &x){
  if(x <= 0.0){ std::cerr << x << " : LogGamma intput should be >= 0" << std::endl; exit(1); }  
  else if(x == 1.0 || x == 2.0){ return 0.0; }
  double x0, x2, xp, gl0, gl;
  int n;
  x0 = x;
  if(x <= 7.0){
    n = int(7 - x);
    x0 = x + n;
  }
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

double CGC(double j1, double m1, double j2, double m2, double jtot, double mtot)
{
  if(fabs(m1 + m2 - mtot) > 0.1){ std::cout << "CGC1" << std::endl; return 0.0; } //projections must add correctly
  else if((jtot < fabs(j1 - j2)) || (jtot > j1 + j2)){ std::cout << "CGC2" << std::endl; return 0.0; } //triangle rule
  else if((fabs(m1) > j1) || (fabs(m2) > j2) || (fabs(mtot) > jtot)){ std::cout << "CGC3" << std::endl; return 0.0; } //unphysical

  double num1, den1, fac1, num2_1, num2_2, num2_3, fac2, den3_2, den3_1, den3_3, den3_4, den3_5, fac3; //numerators, denominators, and factors for CGC
  int change1 = 0, change2 = 0; //flags to change from general formula
  double maxk1, maxk2; //variables to find maximum sum
  double CGC; //clebsch-gordon coefficient
  if(j1 < j2){ std::swap(j1,j2); std::swap(m1,m2); change1 = 1; };  
  if(mtot < 0){ m1 = -m1; m2 = -m2; change2 = 1; }
  mtot = fabs(mtot);
  num1 = (2 * jtot + 1) * factorial(jtot + j1 - j2) * factorial(jtot - j1 + j2) * factorial(j1 + j2 - jtot);
  den1 = double(factorial(j1 + j2 + jtot + 1));
  fac1 = std::sqrt(num1 / den1);
  num2_1 = double(factorial(jtot + mtot) * factorial(jtot - mtot));
  num2_2 = double(factorial(j1 - m1) * factorial(j1 + m1));
  num2_3 = double(factorial(j2 - m2) * factorial(j2 + m2));
  fac2 = std::sqrt(num2_1 * num2_2 * num2_3);
  maxk1 = std::min(j1 + j2 - jtot, j1 - m1);
  maxk2 = std::min(maxk1, j2 + m2);
  fac3 = 0.0;
  for(int k = 0; k <= maxk2+1; ++k){
    den3_1 = j1 + j2 - jtot - k;
    den3_2 = j1 - m1 - k;
    den3_3 = j2 + m2 - k;
    den3_4 = jtot - j2 + m1 + k;
    den3_5 = jtot - j1 - m2 + k;
    if (den3_1 > -0.1 && den3_2 > -0.1 && den3_3 > -0.1 && den3_4 > -0.1 && den3_5 > -0.1){
      fac3 = fac3 + (pow(-1.0, k) / (factorial(k)*factorial(den3_1)*factorial(den3_2)*factorial(den3_3)*factorial(den3_4)*factorial(den3_5)));
    }
  }
  CGC = fac1*fac2*fac3;
  if (change1 == 1){ CGC = CGC*pow(-1.0, int(abs(jtot - j2 - j1) + 0.1)); };
  if (change2 == 1){ CGC = CGC*pow(-1.0, int(abs(jtot - j1 - j2) + 0.1)); };
  return CGC;
}

double CGC3(const double &j1, const double &m1, const double &j2, const double &m2, const double &jtot, const double &mtot)
{
  double threej;  
  if(m1 + m2 + mtot != 0 || abs(m1) > j1 || abs(m2) > j2 || abs(mtot) > jtot || abs(j1 - j2) > jtot || j1 + j2 < jtot){ threej = 0.0; }
  else{ threej = (pow(-1.0, int(abs(j1 - j2 - mtot) + 0.1)) / sqrt(2.0 * jtot + 1.0)) * CGC(j1,m1,j2,m2,jtot,-mtot); }
  return threej;
}

double CGC6(const double &j1, const double &j2, const double &j3, const double &j4, const double &j5, const double &j6)
{
  double sixj = 0.0;
  int S;  
  for(double m1 = -j1; m1 <= j1; m1 = m1 + 1.0){
    for(double m2 = -j2; m2 <= j2; m2 = m2 + 1.0){
      for(double m4 = -j4; m4 <= j4; m4 = m4 + 1.0){
	for(double m5 = -j5; m5 <= j5; m5 = m5 + 1.0){
	  double m3 = m2 + m1, m6 = m1 - m5;
	  if(m3 != m5 - m4 || m6 != -m2 - m4){ continue; }
	  else{
       	    S = int((j1 - m1) + (j2 - m2) + (j3 - m3) + (j4 - m4) + (j5 - m5) + (j6 - m6) + 0.01);
	    sixj = sixj + pow(-1.0, S) * CGC3(j1,m1,j2,m2,j3,-m3) * CGC3(j1,-m1,j5,m5,j6,m6) * CGC3(j4,m4,j5,-m5,j3,m3) * CGC3(j4,-m4,j2,-m2,j6,-m6);
	  }
	}
      }
    }
  }
  return sixj;
}

double CGC9_0(const double &j1, const double &j2, const double &j3, const double &j4, const double &j5, const double &j6)
{
  return pow(-1.0, int(j2 + j5 + j3 + j6)) * CGC6(j1, j2, j5, j4, j3, j6);
}

double Pandya(const double &j1, const double &j2, const double &j3, const double &j4, const double &j13)
{
  double p = 0.0;
  int j12min = int(abs(j1 - j2));
  int j12max = j1 + j2;
  for(int j = j12min; j < j12max; ++j){ p += pow(-1.0, int(j3 + j4 + j)) * (2 * j + 1) * CGC6(j1, j2, j, j4, j3, j13); }
  return p;
}

double Chi_J(const int &a, const int &b, const int &c, const int &d, const int &J1, const int &J2)
{
  double a2 = 0.5*a;
  double b2 = 0.5*b;
  double c2 = 0.5*c;
  double d2 = 0.5*d;
  double J1_2 = 0.5*J1;
  double J2_2 = 0.5*J2;
  double fac1 = std::sqrt((2*J1_2 + 1) * (2*J2_2 + 1));
  return fac1 * pow(-1.0, int(b2 + J1_2 + c2 + J2_2)) * CGC6(a2, b2, J1_2, d2, c2, J2_2);
}

double Legendre(const double &x, const int &l, const int &m)
{
  if(abs(m) > l){ return 0.0; }
  int M = abs(m);
  int L = M;
  double P;
  if(M == 0){ P = 1.0; }
  else{ P = pow(-1.0, M) * factorial2(2*M - 1) * pow(1.0 - x*x, M/2.0); } //P_M^M
  double P1, tempP; //P_L-1^M, reset P1
  while(L != l){
    tempP = P;
    if(L == M){ P = x * (2*L + 1) * P; }
    else{ P = ((2*L + 1) * x * P - (L + M) * P1)/(L - M + 1); }
    P1 = tempP;
    L += 1;
  }
  if(M != m){ P *= pow(-1.0, M) * factorial(l - M) / factorial(l + M); }
  return P;
}

std::complex<double> SphericalY_C(const double &theta, const double &phi, const int &l, const int &m)
{
  if(abs(m) > l){ return 0.0; }
  std::complex<double> I (0.0, 1.0), M (m, 0.0), PHI (phi, 0.0);
  std::complex<double> fac1 (pow(-1.0, m), 0.0), fac2 (sqrt(((2*l + 1)*factorial(l - m))/(4*PI*factorial(l + m))), 0.0);
  std::complex<double> fac3 (Legendre(cos(theta), l, m), 0.0);
  return fac1 * fac2 * fac3 * pow(e, I * M * PHI);
}

double SphericalY(const double &theta, const double &phi, const int &l, const int &m)
{
  if(abs(m) > l){ return 0.0; }
  if(m < 0){ return sqrt(2)*pow(-1.0, m)*SphericalY_C(theta, phi, l, abs(m)).imag(); }
  else if(m == 0){ return SphericalY_C(theta, phi, l, 0).real(); }
  else{ return sqrt(2)*pow(-1.0, m)*SphericalY_C(theta, phi, l, m).real(); }
}

double SphericalYTens(const double &theta, const double &phi, const double &j, const int &l, const int &s, const int &ms)
{
  if(abs(ms) > s){ return 0.0; }
  double Y = 0.0;
  for(int ml = -l; ml <= l; ++ml){ Y += SphericalY(theta, phi, l, ml) * CGC(l, ml, s, ms, j, ml + ms); }
  return Y;
}

double Erf(const double &z)
{
  double t = 1/(1 + 0.5*fabs(z));
  double tau = t * pow(e, -1.0*z*z - 1.26551223 + 1.00002368*t + 0.37409196*t*t + 0.09678418*t*t*t
		       - 0.18628806*t*t*t*t + 0.27886807*t*t*t*t*t - 1.13520398*t*t*t*t*t*t
		       + 1.48851587*t*t*t*t*t*t*t - 0.82215223*t*t*t*t*t*t*t*t + 0.17087277*t*t*t*t*t*t*t*t*t);
  if(z >= 0.0){ tau = 1 - tau; }
  else{ tau = tau - 1; }
  return tau;
}

double Laguerre(const int &k, const double &alpha, const double &x)
{
  if(k == 0){ return 1.0; }
  else if(k == 1){ return 1.0 + alpha - x; }
  else{
    int k0 = 1;
    double L0 = 1.0;
    double L1 = 1.0 + alpha - x;
    double L2 = ((2*k0 + 1.0 + alpha - x)*L1 - (k0 + alpha)*L0)/(k0 + 1);
    L0 = L1;
    L1 = L2;
    ++k0;
    while(k0 < k){
      L2 = ((2*k0 + 1.0 + alpha - x)*L1 - (k0 + alpha)*L0)/(k0 + 1);
      L0 = L1;
      L1 = L2;
      ++k0;
    }
    return L2;
  }
}

double HOfunction(const double &hw, const int &k, const int &l, const int &m, const double &r, const double &theta, const double &phi){
  double nu = 1000000 * m_electronc2 * hw / (2 * hbarc_eVum * hbarc_eVum);
  double N = std::sqrt(std::sqrt(2 * nu * nu * nu / PI) * (std::pow(2, k + 2*l + 3) * factorial(k) * std::pow(nu, l)) / factorial2(2*k + 2*l + 1));
  return N * std::pow(r, l) * std::pow(e, -1.0 * nu * r * r) * Laguerre(k, l+0.5, 2 * nu * r * r) * SphericalY(theta, phi, l, m);
}

void projection(double *u, double *v, double *proj, const int &size)
{
  double innerprod = 0.0;
  double norm = 0.0;
  for(int i = 0; i < size; ++i){ innerprod += v[i]*u[i]; norm += u[i]*u[i]; }
  for(int i = 0; i < size; ++i){ proj[i] = -1.0 * (innerprod/norm)*u[i]; }
}

void GramSchmidt(double **Vectors, const int &size)
{
  double norm;
  double innerprod;
  double **V = new double*[size];
  for(int i = 0; i < size; ++i){
    V[i] = new double[size];
    for(int j = 0; j < size; ++j){ V[i][j] = Vectors[i][j]; }
  }

  norm = 0.0;
  for(int i = 0; i < size; ++i){ norm += V[0][i] * V[0][i]; }
  for(int i = 0; i < size; ++i){ Vectors[0][i] = V[0][i]/std::sqrt(norm); }
  for(int i = 1; i < size; ++i){
    for(int j = 0; j < size; ++j){ Vectors[i][j] = V[i][j]; }
    for(int j = 0; j < i; ++j){
      innerprod = 0.0;
      for(int k = 0; k < size; ++k){ innerprod += Vectors[i][k] * Vectors[j][k]; }
      for(int k = 0; k < size; ++k){ Vectors[i][k] -= innerprod * Vectors[j][k]; }
    }
    norm = 0.0;
    for(int j = 0; j < size; ++j){ norm += Vectors[i][j] * Vectors[i][j]; }
    for(int j = 0; j < size; ++j){ Vectors[i][j] /= std::sqrt(norm); }
  }
  for(int i = 0; i < size; ++i){
    delete[] V[i];
  }
  delete[] V;
}

void GramSchmidt(double *Vectors, const int &size)
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

void Asym_Diagonalize1(double *Ham, int &N, double &eigenvalue, double &norm1p, int &np0)
{
  int info = 0;
  int ind;
  char job1 = 'N';
  char job2 = 'V';
  int lwork = 10*N;
  double *Vl = new double[N * N];
  double *Vr = new double[N * N];
  double *wr = new double[N];
  double *wi = new double[N];
  double *work = new double[lwork];
  double tempen;
  info = 0;
  dgeev_(&job1, &job2, &N, Ham, &N, wr, wi, Vl, &N, Vr, &N, work, &lwork, &info);
  
  ind = -1;
  tempen = 1000;
  for(int i = 0; i < N; ++i){
    if(wr[i] < tempen){
      tempen = wr[i];
      ind = i;
    }
  }
  
  norm1p = 0.0;
  for(int i = 0; i < np0; ++i){ norm1p += Vr[N*ind + i] * Vr[N*ind + i]; }
  norm1p = std::sqrt(norm1p);
  eigenvalue = wr[ind];
  
  //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << wr[ind] << " : " << norm1p << std::endl;
  delete[] Vl;
  delete[] Vr;
  delete[] wr;
  delete[] wi;
  delete[] work;
}

void Asym_Diagonalize2(double *Ham, int &N, double &eigenvalue, double &norm1p, int &np0)
{
  int info = 0;
  int ido = 0; // status integer is zero at start
  char bmat[] = "I"; // standard eigenvalue problem
  char which[] = "SR"; // smallest magnitude eigenvalues
  int nev = 1; // number of eigenvalues to calculate
  double tol = 1.0e-10; // error tolerance
  double *resid = new double[N];
  int ncv = 3*nev + 2; // number of lanczos vectors
  if(ncv > N){ ncv = N; }
  double *v = new double[N*ncv];
  for(int i = 0; i < N*ncv; ++i){ v[i] = 0.0; }
  int ldv = N;
  double *workd = new double[3*N];
  int lworkl = 4*ncv*(ncv + 2);
  double *workl = new double[lworkl];
  info = 0;
  int iparam[11];
  int ipntr[14];
  int ishift = 1;
  int mxiter = 5000;
  int mode = 1;
  iparam[0] = ishift;
  iparam[2] = mxiter;
  iparam[6] = mode;

  do{
    dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);	
    if(ido != -1 && ido != 1){ break; }
    
    #pragma omp parallel for
    for(int j = 0; j < N; ++j){
      workd[ipntr[1]-1 + j] = 0.0;
      for(int k = 0; k < N; ++k){
	workd[ipntr[1]-1 + j] += Ham[N*k + j] * workd[ipntr[0]-1 + k];
      }
    }
  } while(true);
  
  bool rvec = true;
  char howmny = 'A';
  int *select = new int[ncv];
  double *dr = new double[nev + 1];
  double *di = new double[nev + 1];
  double *z = new double[N * (nev + 1)];
  for(int i = 0; i < nev + 1; ++i){
    dr[i] = 0.0;
    di[i] = 0.0;
    for(int j = 0; j < N; ++j){
      z[N*i + j] = 0.0;
    }
  }
  int ldz = N;
  double sigmar;
  double sigmai;
  double *workev = new double[3 * ncv];
  
  dneupd_(&rvec, &howmny, select, dr, di, z, &ldz, &sigmar, &sigmai, workev, bmat, &N, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
  
  norm1p = 0.0;
  for(int i = 0; i < np0; ++i){ norm1p += z[i] * z[i]; }
  norm1p = std::sqrt(norm1p);
  eigenvalue = dr[0];
  
  //std::cout << "! " << Chan.qnums3[chan].ml << " " << Chan.qnums3[chan].m << " : " << dr[0] << " : " << norm1p << std::endl;   
  delete[] resid;
  delete[] v;
  delete[] workd;
  delete[] workl;
  delete[] select;
  delete[] dr;
  delete[] di;
  delete[] z;
  delete[] workev;
}
