/*! @file topologyGeom.cpp
 *  @brief implementation of the TopologyGeom class
 */
#include "topologyGeom.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/numerical_const.hpp"
#include "staples.hpp"
#include <numeric>

using namespace SUNmatUtils;
using namespace FieldUtils;

double TopologyGeom::get_Q(const GaugeField& G) const{
  std::vector<double> q(Nvol_);
  calc_Q(q,G);
  double Q = accumulate(q.begin(),q.end(),0.0);
  return Communicator::instance()->reduce_sum(Q);
}

double TopologyGeom::get_Q(std::vector<double>& q,const GaugeField& G) const{
  q.resize(Nvol_);
  calc_Q(q,G);
  double Q = accumulate(q.begin(),q.end(),0.0);
  return Communicator::instance()->reduce_sum(Q);
}

void TopologyGeom::calc_Q(std::vector<double>& q,const GaugeField& G) const{
  
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  GaugeField1D Bx = stpl_.fieldStrength(G,YDIR,ZDIR);
  GaugeField1D By = stpl_.fieldStrength(G,ZDIR,XDIR);
  GaugeField1D Bz = stpl_.fieldStrength(G,XDIR,YDIR);

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  GaugeField1D Ex = stpl_.fieldStrength(G,TDIR,XDIR);
  GaugeField1D Ey = stpl_.fieldStrength(G,TDIR,YDIR);
  GaugeField1D Ez = stpl_.fieldStrength(G,TDIR,ZDIR);
  
  double coeff = 8.0/(32.0*PI*PI);
  SUNmat u,v;

  // q is the local topology density 
  for(int site=0; site<Nvol_; ++site){
    u = mat(Bx,site); u*= mat(Ex,site);
    v = mat(By,site); v*= mat(Ey,site);
    u+= v;
    v = mat(Bz,site); v*= mat(Ez,site);
    u+= v;
    q[site] = coeff*ReTr(u);
  }
}
