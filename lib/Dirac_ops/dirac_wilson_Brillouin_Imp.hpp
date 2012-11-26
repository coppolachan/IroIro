//----------------------------------------------------------------------
// dirac_wilson_Brillouin_Imp.hpp
//----------------------------------------------------------------------
#ifndef DIRAC_WILSON_BRILLOUIN_IMP_INCLUDED
#define DIRAC_WILSON_BRILLOUIN_IMP_INCLUDED

#include "dirac.hpp"
#include "include/format_F.h"
#include "include/format_G.h"
#include "include/pugi_interface.h"

#include "Main/Geometry/siteIndex.hpp"
#include "Main/Geometry/siteMap.hpp"

typedef Format::Format_F ffmt_t;
typedef Format::Format_G gfmt_t;

class Dirac_Wilson_Brillouin_Imp: public DiracWilsonLike{
private:
  const Field* const u_;
  double kimp_,mass_;
  int Nx_,Ny_,Nz_,Nt_,Nvol_;

  const ffmt_t ff_;
  const gfmt_t gf_;

  const size_t fsize_;
  const size_t gsize_;
  const Communicator* comm_;

  int r0(int c)const{return 2*c;}
  int r1(int c)const{return 2*(NC_+c);}
  int r2(int c)const{return 2*(2*NC_+c);}
  int r3(int c)const{return 2*(3*NC_+c);} 

  int i0(int c)const{return 2*c+1;}
  int i1(int c)const{return 2*(NC_+c)+1;}
  int i2(int c)const{return 2*(2*NC_+c)+1;}
  int i3(int c)const{return 2*(3*NC_+c)+1;} 

  int sr(int s,int c)const{return 2*(s*NC_+c);}
  int si(int s,int c)const{return 2*(s*NC_+c)+1;}

  int re(int c1,int c2)const{return 2*(NC_*c1+c2);}
  int im(int c1,int c2)const{return 2*(NC_*c1+c2)+1;}
  
  // full site operation
  int xsl(int x,int n,int dir)const{return SiteMap::shiftSite.xslice(x,n,dir);}
  int slsize(int dir)const{return SiteMap::shiftSite.slice_size(0,dir);}
  //  int slsize(int x,int dir)const{return SiteMap::shiftSite.slice_size(x,dir);}

  const Field sft_p(const Field& f,int dir)const;
  const Field sft_m(const Field& f,int dir)const;
  const Field delg(const Field& f,int dir,double a,double b)const;
  const Field lap(const Field& f,int dir,double a,double b)const;

  const Field gamma_x(const Field&)const;
  const Field gamma_y(const Field&)const;
  const Field gamma_z(const Field&)const;
  const Field gamma_t(const Field&)const;

  static const Field (Dirac_Wilson_Brillouin_Imp::*gm[])(const Field&)const;
  
  int sgm(int,int,int)const;
public:
 /*! @brief constructor to create instance with normal site indexing */
  Dirac_Wilson_Brillouin_Imp(double mass,const Field* u)
    :kimp_(1.0/(mass+225.0/128.0)),u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){}

  Dirac_Wilson_Brillouin_Imp(const XML::node& node,const Field* u)
    :u_(u),
     Nx_(CommonPrms::instance()->Nx()),
     Ny_(CommonPrms::instance()->Ny()),
     Nz_(CommonPrms::instance()->Nz()),
     Nt_(CommonPrms::instance()->Nt()),
     Nvol_(CommonPrms::instance()->Nvol()),
     comm_(Communicator::instance()),
     ff_(Nvol_),fsize_(ff_.size()),
     gf_(Nvol_),gsize_(gf_.size()){
    //
    double mass;
    XML::read(node, "mass", mass);
    kimp_= 1.0/(mass+225.0/128.0);
  }
  
  size_t fsize() const{return fsize_;}
  size_t gsize() const{return gsize_;}

  const Field mult(const Field&)const;
  const Field mult_dag(const Field&)const;
  //  const Field lap_bri(const Field&)const;
  const Field lap_hop(const Field&)const;
  const Field del_iso(const Field&)const;

  ////////////////////////////////////////Preconditioned versions
  // Wilson operator has no defined preconditioner now 
  const Field mult_prec     (const Field& f)const{return f;}
  const Field mult_dag_prec (const Field& f)const{return f;}
  const Field left_prec     (const Field& f)const{return f;}
  const Field right_prec    (const Field& f)const{return f;}
  const Field left_dag_prec (const Field& f)const{return f;}
  const Field right_dag_prec(const Field& f)const{return f;}
  //////////////////////////////////////////////////////////////

  const Field gamma5(const Field&) const;
  /*
  void gamma5_mult(Field&, const Field&) const;
  void gamma5_ptr(double*, double* const) const;
  const Field proj_p(const Field&) const;
  const Field proj_m(const Field&) const;
  void proj_p(Field&, const Field&, int) const;
  void proj_m(Field&, const Field&, int) const;
  */
  const Field md_force(const Field& , const Field&)const{}
  void md_force_p(Field&,const Field&,const Field&)const{}
  void md_force_m(Field&,const Field&,const Field&)const{}
  
  double getKappa() const {return 1.0/(2.0/kimp_-17.0/4.0);}  
  const ffmt_t get_fermionFormat() const {return ff_;}
  const std::vector<int> get_gsite() const;

  void update_internal_state(){}
};

#endif
