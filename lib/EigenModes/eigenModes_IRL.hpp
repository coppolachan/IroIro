//----------------------------------------------------------------------
// eigenModes_IRL.h
//----------------------------------------------------------------------
#ifndef EIGENMODES_IRL_INCLUDED
#define EIGENMODES_IRL_INCLUDED
#include <vector>
#include <algorithm>

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
  int Np_;
  double enorm_;
  double vthrs_;
  int Niter_;

  void setUnit_Qt(int Nn, std::vector<double>& Qt)const;
  void step(std::vector<double>& lmda,
	    std::vector<double>& lmdb, 
	    std::vector<Field>& evec,
	    Field& f,int Nm,int k)const;

 void qr_decomp(std::vector<double>& lmda,
		std::vector<double>& lmdb,
		int Nk,
		int Nm,
		std::vector<double>& Qt,
		double Dsft, 
		int kmin,
		int kmax)const;

  void diagonalize(std::vector<double>& lmda,
		   std::vector<double>& lmdb, 
		   int Nm2,
		   int Nm,
		   std::vector<double>& Qt)const;

  void orthogonalize(Field& w,
		     const std::vector<Field>& evec,
		     int k)const;

public:
  EigenModes_IRL(const Fopr* fopr,
		 const SortEigen* sort,
		 int Nk,
		 int Np,
		 double enorm,
		 double vthrs,
		 int Niter)
    :opr_(fopr),sort_(sort),
     Nk_(Nk),
     Np_(Np),
     enorm_(enorm),
     vthrs_(vthrs),
     Niter_(Niter){}

  void calc(std::vector<double>& lmd,
	    std::vector<Field>& evec,
	    const Field& b,
	    int& Nsbt,
	    int& Nconv) const;
};

#endif
