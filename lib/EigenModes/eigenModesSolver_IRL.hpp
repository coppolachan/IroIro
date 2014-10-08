/*! @file eigenModesSolver_IRL.hpp
 *  @brief Class to calculate EigenModes using Implicit 
 *  restarted Lanzcos method 
 */
#ifndef EIGENMODESSOLVER_IRL_INCLUDED
#define EIGENMODESSOLVER_IRL_INCLUDED

#include "Fopr/fopr.h"
#include "include/pugi_interface.h"
#include "Communicator/comm_io.hpp"
#include "eigenModesSolver.hpp"
#include "eigenModes.hpp"
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
  double ngPrec_;
  double scale_;
  int Niter_;
  
  void setUnit(std::vector<double>& Qt)const;

  void lanczos_ext(std::vector<double>& ta,std::vector<double>& tb, 
		   std::vector<Field>& V,Field& f,int ini,int fin)const;

 void QRfact_Givens(std::vector<double>& ta,std::vector<double>& tb,
		    std::vector<double>& Qt,
		    int Nk,double sft,int k_min,int k_max)const;

  int diagonalize(std::vector<double>& ta,std::vector<double>& tb, 
		   std::vector<double>& Qt,int Nk)const;

  void orthogonalize(Field& w,const std::vector<Field>& evec,int k)const;
  //void projectOut(std::vector<Field>& evec)const;
public:
  EigenModesSolver_IRL(const Fopr_Herm* fopr,
		       const EigenSorter* esorter,
		       XML::node node)
    :opr_(fopr),esorter_(esorter),ngPrec_(1.0e+5){

    XML::read(node,"Nthreshold",Nthrs_,MANDATORY);
    XML::read(node,"Nk",Nk_,MANDATORY);
    int Np;
    XML::read(node,"Np",Np,MANDATORY);
    Nm_= Nk_+Np;
    XML::read(node,"precision",prec_,MANDATORY);
    XML::read(node,"NGprec_factor",ngPrec_);
    XML::read(node,"max_iter",Niter_,MANDATORY);

    assert(Nk_>Nthrs_);
    CCIO::cout<<"EigenModesSolver_IRL constructed\n";
  }

  EigenModesSolver_IRL(const Fopr_Herm* fopr,const EigenSorter* esorter,
		       int Nk,int Np,double prec,int Niter,
		       double ngPrec=1.0e+5,int Nthrs=10000)
    :opr_(fopr),esorter_(esorter),
     Nthrs_(Nthrs),Nk_(Nk),Nm_(Nk+Np),Niter_(Niter),
     prec_(prec),
     ngPrec_(ngPrec){
    ///
    assert(Nm_>Nk_);
    assert(Nk_>Nthrs_);
    CCIO::cout<<"EigenModesSolver_IRL constructed\n";
  }
  
  void calc(std::vector<double>& lmd,std::vector<Field>& evec,int& Neigen,
	    std::vector<Field>* exvec=NULL)const;
};

#endif
