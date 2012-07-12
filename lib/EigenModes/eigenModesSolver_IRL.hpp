/*! @file eigenModesSolver_IRL.hpp
 *  @brief Class to calculate EigenModes using Implicit 
 *  restarted Lanzcos method 
 */
#ifndef EIGENMODESSOLVER_IRL_INCLUDED
#define EIGENMODESSOLVER_IRL_INCLUDED

#include "include/fopr.h"
#include "eigenModesSolver.hpp"
#include "eigenSorter.hpp"
#include "include/pugi_interface.h"

class EigenSorter;

namespace EigenModes{
  struct LowestModes{};
  struct HighestModes{};
}

class EigenModesSolver_IRL :public EigenModesSolver{
private:
  const Fopr_Herm* opr_;
  const EigenSorter* esorter_;
  int Nk_;
  int Nm_;
  double prec_;
  double thrs_;
  int Niter_;
  
  void setUnit(std::vector<double>& Qt)const;

  void lanczos_init(std::vector<double>& ta,std::vector<double>& tb, 
		    std::vector<Field>& V)const;

  void lanczos_ext(std::vector<double>& ta,std::vector<double>& tb, 
		   std::vector<Field>& V,Field& f)const;

 void QRfact_Givens(std::vector<double>& ta,std::vector<double>& tb,
		      std::vector<double>& Qt,
		      int Nk,double sft,int k_min,int k_max)const;

  void diagonalize(std::vector<double>& ta,std::vector<double>& tb, 
		   std::vector<double>& Qt,int Nk)const;

  void orthogonalize(Field& w,const std::vector<Field>& evec,int k)const;
  
public:
  EigenModesSolver_IRL(const Fopr_Herm* fopr,XML::node node);

  EigenModesSolver_IRL(const Fopr_Herm* fopr,
		       int Nk,int Np,double prec,double thrs,int Niter,
		       EigenModes::LowestModes)
    :opr_(fopr),esorter_(new EigenSorter_low),
     Nk_(Nk),Nm_(Nk+Np),prec_(prec),thrs_(opr_->func(thrs)),
     Niter_(Niter){}

  EigenModesSolver_IRL(const Fopr_Herm* fopr,
		       int Nk,int Np,double prec,double thrs,int Niter,
		       EigenModes::HighestModes)
    :opr_(fopr),esorter_(new EigenSorter_high),
     Nk_(Nk),Nm_(Nk+Np),prec_(prec),thrs_(opr_->func(thrs)),
     Niter_(Niter){}
  
  ~EigenModesSolver_IRL(){delete esorter_;}

  void calc(std::vector<double>& lmd,std::vector<Field>& evec,int& Nsbt)const;
};

#endif
