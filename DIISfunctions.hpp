#ifndef DIISFUNCTIONS_H
#define DIISFUNCTIONS_H

#include <cmath>

struct DIIS{
  int size;
  int count;
  int maxl;
  int N;
  double *p;
  double *delp;
  double *tempdelp;
  double *B;

  DIIS(){};
  void Build(int size, int maxl);
  void Delete();
  void Perform(double *vec, double *vec0, double mix);
  void Update_B1(double *vec);
  void Update_B2(double *vec);
};

#endif
