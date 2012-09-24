/*!--------------------------------------------------------------------------
 * @file boundaryCond.hpp
 * @brief Definition of the BoundaryCond class
 *-------------------------------------------------------------------------*/
#ifndef BOUNDARY_COND_INCLUDED
#define BOUNDARY_COND_INCLUDED

#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Tools/sunMat.hpp"
#include <complex> 


/*!@brief interface class*/
class BoundaryCond{
public:
  virtual void apply_bc(GaugeField& u)const = 0;
  virtual void apply_bc(GaugeField& ue,GaugeField& uo)const = 0;
  virtual ~BoundaryCond(){}
};

/*!@brief periodic boundary condition (nothing happens)*/
class BoundaryCond_periodic: public BoundaryCond{
public:
  BoundaryCond_periodic(){}
  BoundaryCond_periodic(const XML::node&){}

  void apply_bc(GaugeField& u)const{}
  void apply_bc(GaugeField& ue,GaugeField& uo)const{}
};

/*!@brief anti-periodic boundary condition */
class BoundaryCond_antiPeriodic: public BoundaryCond{
  int dir_;
public:
  BoundaryCond_antiPeriodic(int dir):dir_(dir){}
  BoundaryCond_antiPeriodic(const XML::node&);

  void apply_bc(GaugeField& u)const;
  void apply_bc(GaugeField& ue,GaugeField& uo)const;
};

/*!@brief boundary condition with U(1) phase*/
class BoundaryCond_U1phase: public BoundaryCond{
  int dir_;
  std::complex<double> bc_;
public:
  BoundaryCond_U1phase(double theta):bc_(cos(theta),sin(theta)){}
  BoundaryCond_U1phase(const XML::node&);

  void apply_bc(GaugeField& u)const;
  void apply_bc(GaugeField& ue,GaugeField& uo)const;
};

/*!@brief boundary condition with SU(N) matrix*/
class BoundaryCond_SUNmat: public BoundaryCond{
  int dir_;
  SUNmat bm_;
public:
  BoundaryCond_SUNmat(const SUNmat& bm):bm_(bm){}
  BoundaryCond_SUNmat(const XML::node&);

  void apply_bc(GaugeField& u)const;
  void apply_bc(GaugeField& ue,GaugeField& uo)const;
};

///////////////////////////////
BoundaryCond* createBC(const XML::node&);

#endif
