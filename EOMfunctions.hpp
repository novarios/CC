#ifndef EOMFUNCTIONS_H
#define EOMFUNCTIONS_H

struct State;
struct Channels;
struct Amplitudes;
struct Interactions;
struct Eff_Interactions;
struct OP_1b;
struct one_body;
struct three_body;
struct EOM;
struct PA;
struct PR;

struct EOM{
  int N_states;                         //  number of total EOM states
  State *qnums;                         //  quantum numbers for N states
  int *chan_vec;                        //  1b-channels corresponding to each state
  double *del_E;                        //  energy difference
  double *norm_1p;                      //  1-body character of state

  int *nob;                             //  number of one-body components in each of N states
  one_body *ob_vec;                     //  contains nob one-body indices for each of N states
  int *ob_index;                        //  index in ob_vec for start of each sub_vec

  int *ntb;                             //  number of one-body components in each of N states
  one_body *tb_vec;                     //  contains nob one-body indices for each of N states
  int *tb_index;                        //  index in ob_vec for start of each sub_vec

  int *nthb;                            //  number of three-body components in each of N states
  three_body *thb_vec;                  //  contains nthb three-body indices for each of N states
  int *thb_index;                       //  index in thb_vec for start of each sub_vec
  State *thb_qnums;                     //  quantum numbers for 2-body state in 3-body state

  int *nfb;                            //  number of three-body components in each of N states
  three_body *fb_vec;                  //  contains nthb three-body indices for each of N states
  int *fb_index;                       //  index in thb_vec for start of each sub_vec
  State *fb1_qnums;                     //  quantum numbers for 2-body state in 3-body state
  State *fb2_qnums;                     //  quantum numbers for 2-body state in 3-body state

  int *nstate;                          //  number of components in each of N states (nob + nthb)
  double *state_vec_R;                  //  contains nstate R coefficients for each of N states
  double *state_vec_L;                  //  contains nstate L coefficients for each of N states
  int *state_index;                     //  index in state_vec for start of each sub_vec

  EOM(){};
  void count_states(Channels &Chan, int type);
  void count_states2(Channels &Chan);
  void PA1_EOM(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);
  void PR1_EOM(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);
  void Print_EOM_1P(double Energy);

  void PA1_EOM_2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);

  void EOM_1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);
  void EOM_PA(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);
  void EOM_2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);
  void EOM_PR(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints);

  void One_Body_Transition_1PA(Channels &Chan, OP_1b &OP);
  void One_Body_Transition_1PR(Channels &Chan, OP_1b &OP);
  void Delete();
};

struct OP_1b{
  int chan3_length;   //  length of chan_vec1, chan_vec2
  int *chan3_vec1;    //  chan3 of bra
  int *chan3_vec2;    //  chan3 of ket

  int chan2_length;   //  length of chan2_vec
  int *chan2_vec;    //  chan2 of bra-ket

  int *hh_3_index;
  double *hh_3;
  int hh_3_length;
  int *hh_2_index;
  double *hh_2;
  int hh_2_length;
  int *hh_map_ind;

  int *pp_3_index;
  double *pp_3;
  int pp_3_length;
  int *pp_2_index;
  double *pp_2;
  int pp_2_length;
  int *pp_map_ind;

  int *ph_3_index;
  double *ph_3;
  int ph_3_length;
  int *ph_2_index;
  double *ph_2;
  int ph_2_length;
  int *ph_map_ind;

  int *hp_3_index;
  double *hp_3;
  int hp_3_length;
  int *hp_2_index;
  double *hp_2;
  int hp_2_length;
  int *hp_map_ind;

  OP_1b(){};
  OP_1b(Channels &Chan, Amplitudes &Amps, int type); // 0 = Fermi, 1 = Gamow-Teller
  void Gather1(Channels &Chan);
  void Gather2(Channels &Chan);
  void Delete();
};

struct PA{
  double *R;     // r(a|)     ->  chan3_k
  double *R3;    // r(abi'|)  ->  chan3_k
  // R1          // r(ab|i)   ->  chan1_ab,chan3_i
  // R2_1        // r(a|ib')  ->  chan3_a,chan2_ib'
  // R2_2        // r(b|ia')  ->  chan3_b,chan2_ia'
  double *L;     // l(|a)     ->  chan3_k
  double *L3;    // l(|abi')  ->  chan3_k
  // L1          // l(i|ab)   ->  chan3_i,chan1_ab
  // L2_1        // l(ib'|a)  ->  chan2_ib',chan3_a
  // L2_2        // l(ia'|b)  ->  chan2_ia',chan3_b

  int chan30;
  int nh0;
  int np0;
  int npph0;
  int npph1;
  int *pph1_ind1;
  int *pph1_ind2;
  double *pph1_fac1;
  double *pph1_fac0;
  
  int *chan31_num;
  int *chan32_num;

  int R_length;
  int R3_length;
  int R1_length;
  int R2_1_length;
  int R2_2_length;

  int *R1_index;
  int *R2_1_index;
  int *R2_2_index;

  int *Rmap1_ind;

  int *Rmap21_index;
  int *Rmap21_num;
  int *Rmap21_ind;
  double *Rmap21_fac;

  int *Rmap22_index;
  int *Rmap22_num;
  int *Rmap22_ind;
  double *Rmap22_fac;

  int L_length;
  int L3_length;
  int L1_length;
  int L2_1_length;
  int L2_2_length;

  int *L1_index;
  int *L2_1_index;
  int *L2_2_index;

  int *Lmap1_ind;

  int *Lmap21_index;
  int *Lmap21_num;
  int *Lmap21_ind;
  double *Lmap21_fac;

  int *Lmap22_index;
  int *Lmap22_num;
  int *Lmap22_ind;
  double *Lmap22_fac;

  void Setup(Channels &Chan, int chan30);
  void Setup(Channels &Chan, PA &PA1);
  void Delete();
  void Zero_R();
  void Zero_L();
  void Set_R1(double *R1, int length, int offset);
  void Set_R2_1(double *R21, int length, int offset);
  void Set_R2_2(double *R22, int length, int offset);
  void Gather_R1(double *R1, int length, int offset);
  void Gather_R2_1(double *R21, int length, int offset);
  void Gather_R2_2(double *R22, int length, int offset);
  void Set_L1(double *L1, int length, int offset);
  void Set_L2_1(double *L21, int length, int offset);
  void Set_L2_2(double *L22, int length, int offset);
  void Gather_L1(double *L1, int length, int offset);
  void Gather_L2_1(double *L21, int length, int offset);
  void Gather_L2_2(double *L22, int length, int offset);

  void Update_R1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2);
  void Update_R2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2);
  void Update_L1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2);
  void Update_L2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PA &PA2);
};

struct PR{
  double *R;     // r(|i)     ->  chan3_k
  double *R3;    // r(|ija')  ->  chan3_k
  // R1          // r(a|ij)   ->  chan3_a,chan1_ij
  // R2_1        // r(aj'|i)  ->  chan2_aj',chan3_i
  // R2_2        // r(ai'|j)  ->  chan2_ai',chan3_j
  double *L;     // l(i|)     ->  chan3_k
  double *L3;    // l(ija'|)  ->  chan3_k
  // L1          // l(ij|a)   ->  chan1_ij,chan3_a
  // L2_1        // l(i|aj')  ->  chan3_i,chan2_aj'
  // L2_2        // l(j|ai')  ->  chan3_j,chan2_ai'

  int chan30;
  int np0;
  int nh0;
  int nhhp0;
  int nhhp1;
  int *hhp1_ind1;
  int *hhp1_ind2;
  double *hhp1_fac1;
  double *hhp1_fac0;
  
  int *chan31_num;
  int *chan32_num;

  int R_length;
  int R3_length;
  int R1_length;
  int R2_1_length;
  int R2_2_length;

  int *R1_index;
  int *R2_1_index;
  int *R2_2_index;

  int *Rmap1_ind;

  int *Rmap21_index;
  int *Rmap21_num;
  int *Rmap21_ind;
  double *Rmap21_fac;

  int *Rmap22_index;
  int *Rmap22_num;
  int *Rmap22_ind;
  double *Rmap22_fac;

  int L_length;
  int L3_length;
  int L1_length;
  int L2_1_length;
  int L2_2_length;

  int *L1_index;
  int *L2_1_index;
  int *L2_2_index;

  int *Lmap1_ind;

  int *Lmap21_index;
  int *Lmap21_num;
  int *Lmap21_ind;
  double *Lmap21_fac;

  int *Lmap22_index;
  int *Lmap22_num;
  int *Lmap22_ind;
  double *Lmap22_fac;

  void Setup(Channels &Chan, int chan30);
  void Setup(Channels &Chan, PR &PR1);
  void Delete();
  void Zero_R();
  void Zero_L();
  void Set_R1(double *R1, int length, int offset);
  void Set_R2_1(double *R21, int length, int offset);
  void Set_R2_2(double *R22, int length, int offset);
  void Gather_R1(double *R1, int length, int offset);
  void Gather_R2_1(double *R21, int length, int offset);
  void Gather_R2_2(double *R22, int length, int offset);
  void Set_L1(double *L1, int length, int offset);
  void Set_L2_1(double *L21, int length, int offset);
  void Set_L2_2(double *L22, int length, int offset);
  void Gather_L1(double *L1, int length, int offset);
  void Gather_L2_1(double *L21, int length, int offset);
  void Gather_L2_2(double *L22, int length, int offset);

  void Update_R1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2);
  void Update_R2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2);
  void Update_L1(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2);
  void Update_L2(Channels &Chan, Interactions &Ints, Amplitudes &Amps, Eff_Interactions &Eff_Ints, PR &PR2);
};

#endif
