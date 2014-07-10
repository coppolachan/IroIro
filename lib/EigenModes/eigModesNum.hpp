#ifndef FOPR_EIGMODESNUM_INCLUDED
#define FOPR_EIGMODESNUM_INCLUDED
#include "Fopr/fopr.h"
#include <vector>

class RandNum;

class EigModesNum{
  double prec_;
  double epsilon_;
  double Mratio_;
  int max_iter_;
  int Nrand_;
  int Npoly_;

  mutable std::vector<double> cutoffs_;
  const Fopr_DdagD* DdagD_;
  const RandNum* rng_;  

  std::vector<double> coeff_;
  std::vector<double> pole_;
  static double func_sgn(double,double,const std::vector<double>&);
  static double func_rh4(double,void*);
  static double func_sqrtInv(double,void*);
  void init_ChebyshevApprox();
  void hproj(Field& f,double Msq)const;
public:
  EigModesNum(double prec,int max_iter,int Nrand,int Npoly,double epsilon,
	      const std::vector<double>& cutoff,
	      const Dirac* D,const RandNum* rng)
    :prec_(prec),max_iter_(max_iter),Nrand_(Nrand),
     Npoly_(Npoly),coeff_(Npoly+1),
     epsilon_(epsilon),cutoffs_(cutoff),
     DdagD_(new Fopr_DdagD(D)),rng_(rng){ init_ChebyshevApprox();}

  EigModesNum(XML::node node,const Dirac* D,const RandNum* rng)
    :DdagD_(new Fopr_DdagD(D)),rng_(rng){
    
    XML::read(node,"precision",    prec_,    MANDATORY);
    XML::read(node,"max_iteration",max_iter_,MANDATORY);
    XML::read(node,"Nrand",        Nrand_,   MANDATORY);
    XML::read(node,"epsilon",      epsilon_, MANDATORY);
    XML::read(node,"Npoly",        Npoly_,   MANDATORY);
    XML::read_array(node,"cutoffs",cutoffs_);
    
    coeff_.resize(Npoly_+1);
    init_ChebyshevApprox();
  }

  ~EigModesNum(){ if(DdagD_) delete DdagD_;  }

  void do_count()const;
};

#endif
