#include "DIISfunctions.hpp"
#include "MATHfunctions.hpp"

void DIIS::Build(int size, int maxl)
{
  this->count = 0;
  this->N = 0;
  this->size = size;
  this->maxl = maxl;
  B = new double[1];
  B[0] = 0.0;

  this->p = new double[maxl * size];
  this->delp = new double[maxl * size];
  this->tempdelp = new double[size];
  for(int l = 0; l < maxl; ++l){
    for(int ind = 0; ind < size; ++ind){
      this->p[l * ind] = 0.0;
      this->delp[l * ind] = 0.0;
    }
  }
  for(int ind = 0; ind < size; ++ind){
    this->tempdelp[ind] = 0.0;
  }
}

void DIIS::Delete()
{
  delete[] B;
  delete[] p;
  delete[] delp;
  delete[] tempdelp;
}

void DIIS::Perform(double *vec, double *vec0, double mix)
{
  double norm = 0.0;
  double norm0 = 0.0;
  double tempt;
  int ind;
  bool ortho;
  int P = this->N + 1;
  int lwork;
  int *ipiv;
  double *work;
  int info;
  double *B2;
  double dot_prod, norm1, norm2;  
  
  // Fill temp_delp
  for(ind = 0; ind < this->size; ++ind){
    this->tempdelp[ind] = (vec[ind] - vec0[ind]) / mix;
    norm0 += vec0[ind] * vec0[ind];
  }
  norm0 = std::sqrt(norm0);
  for(ind = 0; ind < this->size; ++ind){
    this->tempdelp[ind] /= norm0;
    norm += this->tempdelp[ind] * this->tempdelp[ind];
  }

  // check orthogonality of tempdelp
  ortho = true;
  for(int l = 0; l < this->N; ++l){
    dot_prod = 0.0;
    norm1 = 0.0;
    norm2 = 0.0;
    for(ind = 0; ind < this->size; ++ind){
      dot_prod += vec[ind] * this->p[l * this->size + ind];
      norm1 += vec[ind] * vec[ind];
      norm2 += this->p[l * this->size + ind] * this->p[l * this->size + ind];
    }
    dot_prod /= std::sqrt(norm1 * norm2);
    //std::cout << "!! norm, dot_prod, B[P*l+l] " << norm << " " << dot_prod << " " << B[P*l + l] << std::endl;
    if(norm > B[P * l + l] || dot_prod > (1.0 - std::pow(norm, 1.5))){ ortho = false; break; }
  }
  if( !ortho ){ return; }

  ++this->count;
  if(this->N < this->maxl){ this->Update_B1(vec); }
  else{ this->Update_B2(vec); }
  P = this->N + 1;

  // copy B into B2 for inversion
  B2 = new double[P * P];
  for(int j = 0; j < P*P; ++j){ B2[j] = this->B[j]; }
  
  ipiv = new int[P];
  work = new double[sizeof(double) * P];
  lwork = sizeof(double) * P;
  info = 0;
  dgetrf_(&P, &P, B2, &P, ipiv, &info);
  dgetri_(&P, B2, &P, ipiv, work, &lwork, &info);
  for(ind = 0; ind < this->size; ++ind){
    vec[ind] = 0.0;
    tempt = 0.0;
    for(int l = 0; l < this->N; ++l){ tempt += -1.0 * B2[P * l + this->N] * this->p[l * this->size + ind]; }
    vec[ind] = tempt;
  }
  delete[] B2;
  delete[] ipiv;
  delete[] work;
}

void DIIS::Update_B1(double *vec)
{
  double *B2; // used to copy B into updated B
  int ind;
  int P = this->N+1;
  for(ind = 0; ind < this->size; ++ind){
    this->p[N * size + ind] = vec[ind];
    this->delp[N * size + ind] = this->tempdelp[ind];
  }

  B2 = new double[(P+1) * (P+1)];
  for(int j = 0; j < this->N; ++j){
    for(int k = 0; k < this->N; ++k){
      B2[(P+1) * j + k] = B[P * j + k];
    }
  }
  for(int l = 0; l < P; ++l){
    B2[(P+1) * this->N + l] = 0.0;
    if(l != this->N){ B2[(P+1) * l + this->N] = 0.0; }
    for(ind = 0; ind < size; ++ind){
      B2[(P+1) * this->N + l] += this->delp[this->N * this->size + ind] * this->delp[l * this->size + ind];
      if(l != N){ B2[(P+1) * l + this->N] += this->delp[l * this->size + ind] * this->delp[this->N * this->size + ind]; }
    }
  }
  for(int l = 0; l < P; ++l){
    B2[(P+1) * P + l] = -1.0;
    B2[(P+1) * l + P] = -1.0;
  }
  B2[(P+1) * P + P] = 0.0;
  ++N;
  ++P;
  delete[] B;
  this->B = new double[P * P];
  for(int ind = 0; ind < P * P; ++ind){ this->B[ind] = B2[ind]; }
  delete[] B2;
}

void DIIS::Update_B2(double *vec)
{
  int maxind, ind;
  double maxnorm;
  int P = this->N + 1;

  // Find largest norm of B to remove that vector
  maxind = -1;
  maxnorm = 0.0;
  for(int j = 0; j < this->N; ++j){
    if(std::fabs(this->B[P * j + j]) > maxnorm){
      maxind = j;
      maxnorm = std::fabs(this->B[P * j + j]);
    }
  }
  // Replace maxnorm vector with new vector
  for(ind = 0; ind < this->size; ++ind){
    this->p[maxind * this->size + ind] = vec[ind];
    this->delp[maxind * this->size + ind] = this->tempdelp[ind];
  }
  for(int l = 0; l < this->N; ++l){
    this->B[P * maxind + l] = 0.0;
    if(l != maxind){ this->B[P * l + maxind] = 0.0; }	
    for(ind = 0; ind < this->size; ++ind){
      this->B[P * maxind + l] += this->delp[maxind * this->size + ind] * this->delp[l * this->size + ind];
      if(l != maxind){ this->B[P * l + maxind] += this->delp[l * this->size + ind] * this->delp[maxind * this->size + ind]; }
    }
  }
}
