/*! @file eigenModesSolver_IRL.hpp
 *  @brief Class to calculate EigenModes using Implicit 
 *  restarted Lanzcos method 
 */
#ifndef EIGENMODESSOLVER_IRL_INCLUDED
#define EIGENMODESSOLVER_IRL_INCLUDED

#include "include/fopr.h"
#include "include/pugi_interface.h"
#include "Communicator/comm_io.hpp"
#include "eigenModesSolver.hpp"
#include "eigenSorter.hpp"

class EigenSorter;

class EigenModesSolver_IRL :public EigenModesSolver{
private:
  const Fopr_Herm* opr_;
  const EigenSorter* esorter_;
  int Nthrs_;
  int Nk_;
  int Nm_;
  double prec_;
  double scale_;
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
  EigenModesSolver_IRL(const Fopr_Herm* fopr,const EigenSorter* esorter,
		       XML::node node)
    :opr_(fopr),esorter_(esorter){
    CCIO::cout<<"EigenModesSolver_IRL being constructed\n";

    XML::read(node,"Nthreshold",Nthrs_,MANDATORY);
    XML::read(node,"Nk",Nk_,MANDATORY);
    int Np;
    XML::read(node,"Np",Np,MANDATORY);
    Nm_= Nk_+Np;

    XML::read(node,"precision",prec_,MANDATORY);
    XML::read(node,"max_iter",Niter_,MANDATORY);
    
    double scale = fopr->func(esorter->thrs());
    CCIO::cout<<"scaling factor = "<<scale<<"\n";
    prec_*= scale;
    CCIO::cout<<"net precision = "<<prec_<<"\n";
  }
  EigenModesSolver_IRL(const Fopr_Herm* fopr,const EigenSorter* esorter,
		       int Nk,int Np,double prec,int Niter,int Nthrs=10000)
    :opr_(fopr),esorter_(esorter),
     Nthrs_(Nthrs),Nk_(Nk),Nm_(Nk+Np),Niter_(Niter),
     prec_(prec*fopr->func(esorter->thrs())){}
  
  void calc(std::vector<double>& lmd,std::vector<Field>& evec,int& Nsbt)const;
};

#endif
