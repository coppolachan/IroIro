/*! @file topologyGeom.cpp
 *  @brief implementation of the TopologyGeom class
 */
#include "topologyGeom.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "include/numerical_const.hpp"
#include "staples.hpp"

using namespace SUNmatUtils;
using namespace FieldUtils;

double TopologyGeom::getQ(const GaugeField& G) const{
  
  // Bx = -iF(y,z), By = -iF(z,y), Bz = -iF(x,y)
  GaugeField1D Bx = stpl_.fieldStrength(G,YDIR,ZDIR);
  GaugeField1D By = stpl_.fieldStrength(G,ZDIR,XDIR);
  GaugeField1D Bz = stpl_.fieldStrength(G,XDIR,YDIR);

  // Ex = -iF(t,x), Ey = -iF(t,y), Ez = -iF(t,z)
  GaugeField1D Ex = stpl_.fieldStrength(G,TDIR,XDIR);
  GaugeField1D Ey = stpl_.fieldStrength(G,TDIR,YDIR);
  GaugeField1D Ez = stpl_.fieldStrength(G,TDIR,ZDIR);

  double Q=0.0;   // Q = ReTr F~*F 
  SUNmat u,v;
  for(int site=0; site<Nvol_; ++site){
    u = mat(Bx,site); u*= mat(Ex,site);
    v = mat(By,site); v*= mat(Ey,site);
    u+= v;
    v = mat(Bz,site); v*= mat(Ez,site);
    u+= v;
    Q += ReTr(u);
  }
  Q = Communicator::instance()->reduce_sum(Q);
  return Q*8.0/(32.0*PI*PI);
}
