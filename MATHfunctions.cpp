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

double CGC(double j1, double m1, double j2, double m2, double jtot, double mtot)
{
  //std::cout << "! " << j1 << " " << m1 << " " << j2 << " " << m2 << " " << jtot << " " << mtot << std::endl;
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
  fac1 = sqrt(num1 / den1);
  num2_1 = double(factorial(jtot + mtot) * factorial(jtot - mtot));
  num2_2 = double(factorial(j1 - m1) * factorial(j1 + m1));
  num2_3 = double(factorial(j2 - m2) * factorial(j2 + m2));
  fac2 = sqrt(num2_1 * num2_2 * num2_3);

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
  
  //std::cout << fac1 << " " << fac2 << " " << fac3 << std::endl;

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
  //std::cout << theta << " " << phi << " " << l << " " << m << std::endl;
  if(abs(m) > l){ return 0.0; }
  std::complex<double> I (0.0, 1.0), M (m, 0.0), PHI (phi, 0.0);
  std::complex<double> fac1 (pow(-1.0, m), 0.0), fac2 (sqrt(((2*l + 1)*factorial(l - m))/(4*PI*factorial(l + m))), 0.0);
  std::complex<double> fac3 (Legendre(cos(theta), l, m), 0.0);
  return fac1 * fac2 * fac3 * pow(e, I * M * PHI);
}

double SphericalY(const double &theta, const double &phi, const int &l, const int &m)
{
  //std::cout << theta << " " << phi << " " << l << " " << m << std::endl;
  if(abs(m) > l){ return 0.0; }
  if(m < 0){ return sqrt(2)*pow(-1.0, m)*SphericalY_C(theta, phi, l, abs(m)).imag(); }
  else if(m == 0){ return SphericalY_C(theta, phi, l, 0).real(); }
  else{ return sqrt(2)*pow(-1.0, m)*SphericalY_C(theta, phi, l, m).real(); }
}

double SphericalYTens(const double &theta, const double &phi, const double &j, const int &l, const int &s, const int &ms)
{
  if(abs(ms) > s){ return 0.0; }
  double Y = 0.0;
  for(int ml = -l; ml <= l; ++ml){
    Y += SphericalY(theta, phi, l, ml) * CGC(l, ml, s, ms, j, ml + ms);
  }
  //std::cout << theta << " " << phi << " " << j << " " << l << " " << s << " " << ms << std::endl;
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

