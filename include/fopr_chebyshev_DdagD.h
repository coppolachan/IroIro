/*! @file fopr_chebyshev_DdagD.h
 * @brief Definition of Chebyshev operator
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

  Fopr_Chebyshev_DdagDParams(XML::node node){
    double vmin, vmax;
    XML::read(node, "Npoly", Npoly);
    XML::read(node, "vmin", vmin);
    XML::read(node, "vmax", vmax);
    
    Fcb1 = 2.0/(vmax*vmax-vmin*vmin);
    Fcb2 = -(vmax*vmax+vmin*vmin)/(vmax*vmax-vmin*vmin);
  }

  Fopr_Chebyshev_DdagDParams(int Npoly_,double vmin_,double vmax_){
    Npoly = Npoly_;
    Fcb1 = 2.0/(vmax_*vmax_-vmin_*vmin_);
    Fcb2 = -(vmax_*vmax_+vmin_*vmin_)/(vmax_*vmax_-vmin_*vmin_);
  }
};

class Fopr_Chebyshev_DdagD : public Fopr_Herm {
private:
  const Fopr_Chebyshev_DdagDParams Params;
  const Dirac* D_;
public:
  Fopr_Chebyshev_DdagD(int Npoly,double vmin,double vmax,const Dirac* D)
    :Params(Npoly,vmin,vmax),D_(D){}

  Fopr_Chebyshev_DdagD(XML::node Chebyshev_node,const Dirac* D)
    :Params(Chebyshev_node),D_(D){}
  
  ~Fopr_Chebyshev_DdagD(){}

  double func(double x) const;
  const Field mult(const Field& f) const;
  const Field mult_dag(const Field& f) const{ return mult(f);}

  size_t fsize()const {return D_->fsize();}
};
#endif
