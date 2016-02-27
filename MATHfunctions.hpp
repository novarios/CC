#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <complex>

#define PI 3.141592654
#define e 2.718281828

long long factorial(const int &n);
long long factorial(const double &n);
long long factorial2(const int &n);
long long factorial2(const double &n);
double CGC(double j1, double m1, double j2, double m2, double jtot, double mtot);
double CGC3(const double &j1, const double &m1, const double &j2, const double &m2, const double &jtot, const double &mtot);
double CGC6(const double &j1, const double &j2, const double &j3, const double &j4, const double &j5, const double &j6);
double Legendre(const double &x, const int &l, const int &m);
std::complex<double> SphericalY_C(const double &theta, const double &phi, const int &l, const int &m);
double SphericalY(const double &theta, const double &phi, const int &l, const int &m);
double SphericalYTens(const double &theta, const double &phi, const double &j, const int &l, const int &s, const int &ms);
double Erf(const double &z);

#endif
