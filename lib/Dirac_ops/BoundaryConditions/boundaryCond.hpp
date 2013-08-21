/*!--------------------------------------------------------------------------
 * @file boundaryCond.hpp
 * @brief Definition of the BoundaryCond class
 *-------------------------------------------------------------------------*/
#ifndef BOUNDARY_COND_INCLUDED
#define BOUNDARY_COND_INCLUDED

#include "include/pugi_interface.h"
#include "include/common_fields.hpp"
#include "Tools/sunMat.hpp"
#include "antiPeriodicBC.hpp"
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
  AntiPeriodicBC<GaugeField>* apbc_;
public:
  BoundaryCond_antiPeriodic(int dir)
  :apbc_(new AntiPeriodicBC<GaugeField>(dir)){}
 
  BoundaryCond_antiPeriodic(const XML::node& bcnode):apbc_(NULL){
    int dir;
    const char* dir_name = bcnode.attribute("dir").value();
    if(     !strcmp(dir_name,"X")) dir= XDIR;
    else if(!strcmp(dir_name,"Y")) dir= YDIR;
    else if(!strcmp(dir_name,"Z")) dir= ZDIR;
    else if(!strcmp(dir_name,"T")) dir= TDIR;
    else {
      CCIO::cout<<"No valid direction available\n";
      abort();
    }
    apbc_=new AntiPeriodicBC<GaugeField>(dir);
  }
  ~BoundaryCond_antiPeriodic(){if(apbc_) delete apbc_;}

  void apply_bc(GaugeField& u)const{apbc_->apply_bc(u);}
  void apply_bc(GaugeField& ue,GaugeField& uo)const{apbc_->apply_bc(ue,uo);}
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
