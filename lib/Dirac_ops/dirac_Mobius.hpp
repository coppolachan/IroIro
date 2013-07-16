/*! @file dirac_Mobius_EvenOdd.hpp
 *  @brief Definition of Mobius operator
 Time-stamp: <2013-07-02 13:47:55 noaki>
*/

#ifndef DIRAC_MOBIUS_INCLUDED
#define DIRAC_MOBIUS_INCLUDED

#include "dirac_WilsonLike.hpp"
class Solver;

class Dirac_Mobius : public DiracWilsonLike {
private:
  const DiracWilsonLike* Dd_; // Dirac_op for denominator
  const Solver* slv_; // solver for denominator 
  double c1_,c2_; /*!< @brief D_Mobius = c1_ + c2_ D^-1 (massd_) 
 		              where c1_=b/c, c2_ = -2b/(c^2 * (mass_d+4)*/
  double mass_;
  Dirac_Mobius(const Dirac_Mobius&); /*!< @brief simple copy is prohibited */
public:
  Dirac_Mobius(const XML::node& DMnode,
	       const DiracWilsonLike* D,const Solver* slv)
    :Dd_(D),slv_(slv){

    double b,c;
    XML::node Dnode = DMnode;
    XML::descend(Dnode,"Dirac_denominator",MANDATORY);
    XML::read(Dnode,"mass",mass_,MANDATORY);
    XML::read(DMnode,"b",b,MANDATORY);
    XML::read(DMnode,"c",c,MANDATORY);

    if(c == 0.0){
      std::cout<<"c must be positive !\n";
      abort();
    }else{
      c1_=b/c;
      c2_=-2.0*b/(c*c*(mass_+4.0));
      // remove  *(mass_ +4.0)  after D_W is canonically normalized !!!!  
      mass_-= 2.0/c;
    }
  }
  
  Dirac_Mobius(const DiracWilsonLike* D,const Solver* slv,
	       double b,double c,double massd)
    :Dd_(D),slv_(slv),
     c1_(b/c),c2_(-2.0*b/(c*c*(massd+4.0))),mass_(massd-2.0/c){}

  size_t fsize() const{return Dd_->fsize();}
  size_t gsize() const{return Dd_->gsize();}
  double getMass()const{return mass_;}

  const Field* getGaugeField_ptr()const{ return Dd_->getGaugeField_ptr(); }

  const Field gamma5(const Field& f) const{return Dd_->gamma5(f);}
  const Field mult(const Field&) const;
  const Field mult_dag(const Field&) const;
};

#endif
