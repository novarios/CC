#ifndef ANGMOM_H
#define ANGMOM_H

struct Input_Parameters;
struct Model_Space;
struct AngMom;

int count1(int N);
int count2(int N);

struct AngMom{
  int jmax;
  int jsize;
  //int Jsize;

  //double ***CGC_vec;
  double ***ThreeJ_vec;
  double ***SixJ_vec;

  AngMom(){};
  void Build();
  void Delete();
  //double CGC_1(int j1, int m1, int j2, int m2, int jtot, int mtot);
  //double CGC_large_1(int j1, int m1, int j2, int m2, int jtot, int mtot);
  double CGC_1(int j1, int m1, int j2, int m2, int j3, int m3);
  double CGC_large_1(int j1, int m1, int j2, int m2, int j3, int m3);
  double SixJ_1(int j1, int j2, int j3, int j4, int j5, int j6);
  //double CGC3_1(int j1, int m1, int j2, int m2, int jtot, int mtot);
  //double CGC6_1(int j1, int j2, int j3, int j4, int j5, int j6);
  double get_CGC(int j1, int m1, int j2, int m2, int j3, int m3);
  //double get_CGC0(int j1, int m1, int j2, int m2, int jtot, int mtot);
  double get_ThreeJ(int j1, int m1, int j2, int m2, int j3, int m3);
  double get_SixJ(int j1, int j2, int j3, int j4, int j5, int j6);
  double get_NineJ(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9);
};
extern struct AngMom JCOUP;

#endif
