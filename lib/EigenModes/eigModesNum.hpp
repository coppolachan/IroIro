#ifndef FOPR_EIGMODESNUM_INCLUDED
#define FOPR_EIGMODESNUM_INCLUDED

#include "fopr.h"
#include <vector>

class RandNum;

class EigModesNum{

  double prec_;
  double epsilon_;
  int max_iter_;
  int Nrand_;
  int Ndeg_;

  mutable std::vector<double> cutoffs_;
  const Fopr_DdagD* DdagD_;
  const RandNum* rng_;  

  double r0_;
  std::vector<double> res_;
  std::vector<double> pole_;
  void init_RationalApprox();
public:
  EigModesNum(double prec,int max_iter,int Nrand,int Ndeg,double epsilon,
	      const std::vector<double>& cutoff,
	      const Dirac* D,const RandNum* rng)
    :prec_(prec),max_iter_(max_iter),Nrand_(Nrand),Ndeg_(Ndeg),
     epsilon_(epsilon),cutoffs_(cutoff),
     DdagD_(new Fopr_DdagD(D)),rng_(rng){ init_RationalApprox();}

  EigModesNum(XML::node node,const Dirac* D,const RandNum* rng)
    :DdagD_(new Fopr_DdagD(D)),rng_(rng){
    
    XML::read(node,"precision",    prec_,    MANDATORY);
    XML::read(node,"max_iteration",max_iter_,MANDATORY);
    XML::read(node,"Nrand",        Nrand_,   MANDATORY);
    XML::read(node,"epsilon",      epsilon_, MANDATORY);
    XML::read(node,"Ndegree",      Ndeg_,    MANDATORY);
    XML::read_array(node,"cutoffs",cutoffs_);

    init_RationalApprox();
  }

  ~EigModesNum(){ if(DdagD_) delete DdagD_;  }

  void count()const;
  void count(const std::vector<double>&)const;
};

#endif
