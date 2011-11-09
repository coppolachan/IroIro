/*!
 * @file fopr_chebyshev_DdagD.h
 *
 * @brief Definition of Chebyshev operator
 *
 */

#ifndef FOPR_CHEBYSHEV_INCLUDED
#define FOPR_CHEBYSHEV_INCLUDED

#include "include/pugi_interface.h"
#include "fopr.h"

#include <vector>
#include <iomanip>

struct Fopr_Chebyshev_DdagDParams{
  int Npoly;
  double Fcb1,Fcb2;

  Fopr_Chebyshev_DdagDParams(XML::node node) {
    double vthrs, vmax;
    XML::read(node, "NumPolynomials", Npoly);
    XML::read(node, "Threshold", vthrs);
    XML::read(node, "MaxEigenmode", vmax);
    
    Fcb1 = 2.0/(vmax*vmax-vthrs*vthrs);
    Fcb2 = -(vmax*vmax+vthrs*vthrs)/(vmax*vmax-vthrs*vthrs);
  }

  Fopr_Chebyshev_DdagDParams(int Npoly_, 
			     double vthrs_, 
			     double vmax_) {
    Npoly = Npoly_;
    Fcb1 = 2.0/(vmax_*vmax_-vthrs_*vthrs_);
    Fcb2 = -(vmax_*vmax_+vthrs_*vthrs_)/(vmax_*vmax_-vthrs_*vthrs_);
  }
};

class Fopr_Chebyshev_DdagD : public Fopr {
private:
  const Fopr_Chebyshev_DdagDParams Params;
  const Dirac* D_;
public:
  Fopr_Chebyshev_DdagD(int Npoly,
		       double vthrs,
		       double vmax,
		       const Dirac* D)
    :Params(Npoly, vthrs, vmax),
     D_(D){}

  Fopr_Chebyshev_DdagD(XML::node Chebyshev_node,
		       const Dirac* D)
    :Params(Chebyshev_node),
     D_(D){}

  ~Fopr_Chebyshev_DdagD(){delete D_;}

  const Field mult(const Field& f) const;
  const Field mult_dag(const Field& f) const{ return mult(f);}
  double mult(double x) const;
  size_t fsize()const {return D_->fsize();}
};
#endif
