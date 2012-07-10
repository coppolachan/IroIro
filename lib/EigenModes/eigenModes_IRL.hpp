//----------------------------------------------------------------------
// eigenModes_IRL.h
//----------------------------------------------------------------------
#ifndef EIGENMODES_IRL_INCLUDED
#define EIGENMODES_IRL_INCLUDED
#include <vector>
#include "include/fopr.h"

class Field;
class SortEigen;

/*
  @brief Class to calculate EigenModes using Implicit restarted Lanzcos method
*/
class EigenModes_IRL{
private:
  const Fopr* opr_;
  const SortEigen* sort_;
  int Nk_;
  int Nm_;
  double enorm_;
  double vthrs_;
  int Niter_;

  void setUnit(std::vector<double>& Qt)const;

  void decomp_Lanczos(std::vector<double>& lmda,
		      std::vector<double>& lmdb, 
		      Field& f,
		      std::vector<Field>& v,int k)const;

 void qr_decomp(std::vector<double>& lmda,
		std::vector<double>& lmdb,
		std::vector<double>& Qt,
		int Nk,double sft,int k_min,int k_max)const;

  void diagonalize(std::vector<double>& lmda,
		   std::vector<double>& lmdb, 
		   std::vector<double>& Qt,int Nk)const;

  void orthogonalize(Field& w,const std::vector<Field>& evec,int k)const;
public:
  EigenModes_IRL(const Fopr* fopr,const SortEigen* sort,
		 int Nk,int Np,double enorm,double vthrs,int Niter)
    :opr_(fopr),sort_(sort),
     Nk_(Nk),Nm_(Nk+Np),enorm_(enorm),vthrs_(vthrs),
     Niter_(Niter){}

  void calc(std::vector<double>& lmd,std::vector<Field>& evec,int& Nsbt)const;
};

#endif
