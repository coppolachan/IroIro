/*!
 * @file source_types.hpp
 * @brief Definition of source types for propagator calculations
 * Include this file if you need to specify explicitly the source type
 */
#ifndef SOURCES_TYPES_H_
#define SOURCES_TYPES_H_

#include "include/numerical_const.hpp"
#include "Geometry/siteIndex.hpp"
#include "Tools/randNum.h"
#include "Tools/randNum_MP.h"
#include "source.hpp"
#include "Tools/EnumToString.hpp"

#include <cmath>
#include <iostream>

////// Source_local -------------------------------------------------------

template<typename FMT> class Source_local :public Source{
private:
  std::vector<int> sp_;
  FMT* ff_;
  bool has_source_;
  int loc_;
  void set_src();
public:
  Source_local(std::vector<int> sp, int Nvol)
    :sp_(sp),ff_(new FMT(Nvol)),has_source_(false),loc_(0){
    set_src();
  }
  ~Source_local(){ delete ff_;}
  const Field mksrc(int s, int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename FMT> 
void Source_local<FMT>::set_src(){

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  int Ndim = CommonPrms::instance()->Ndim();

  int nid = Communicator::instance()->nodeid(sp_[XDIR]/Nx,sp_[YDIR]/Ny,
					     sp_[ZDIR]/Nz,sp_[TDIR]/Nt);
  if(Communicator::instance()->nodeid() == nid){
    has_source_= true;
    loc_= SiteIndex::instance()->site(sp_[XDIR]%Nx,sp_[YDIR]%Ny,
				      sp_[ZDIR]%Nz,sp_[TDIR]%Nt);
  }
}

template<typename FMT> const Field Source_local<FMT>::
mksrc(int s,int c)const{
  Field src(ff_->size());
  if(has_source_) src.set(ff_->index_r(c,s,loc_),1.0);
  return src;
}

template<typename FMT> const Field Source_local<FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}

////// Source_exp   -----------------------------------------------------
template<typename FMT> class Source_exp :public Source{
private:
  std::vector<int> sp_;
  FMT* ff_;
  std::valarray<double> smr_;
  double alpha_;
  void set_src();
public:
  Source_exp(std::vector<int> sp,double al,int Nvol)
    :sp_(sp),ff_(new FMT(Nvol)),
     alpha_(al),smr_(0.0,ff_->Nvol()){
    set_src();
  }
  ~Source_exp(){ delete ff_;}
  const Field mksrc(int s,int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename FMT> void Source_exp<FMT>::set_src(){

  int tb = SiteIndex::instance()->global_t(0);
  int tt = SiteIndex::instance()->global_t(SiteIndex::instance()->Bdir(TDIR));

  if(tb <= sp_[TDIR] && sp_[TDIR]<= tt){
    int ti = sp_[TDIR]%CommonPrms::instance()->Nt();

    for(int sl= 0; sl<SiteIndex::instance()->slsize(ti,TDIR); ++sl){

      int site = SiteIndex::instance()->slice_t(ti,sl);
      int gsite = SiteIndex::instance()->get_gsite(site);

      int xi = SiteIndex::instance()->g_x(gsite) -sp_[XDIR];
      int yi = SiteIndex::instance()->g_y(gsite) -sp_[YDIR];
      int zi = SiteIndex::instance()->g_z(gsite) -sp_[ZDIR];

      int xr = CommonPrms::instance()->Lx() -xi;
      int yr = CommonPrms::instance()->Ly() -yi;
      int zr = CommonPrms::instance()->Lz() -zi;
     
      smr_[site] = exp(-alpha_*sqrt(double(xi*xi +yi*yi +zi*zi)))
	+          exp(-alpha_*sqrt(double(xr*xr +yi*yi +zi*zi)))
	+          exp(-alpha_*sqrt(double(xi*xi +yr*yr +zi*zi)))
	+          exp(-alpha_*sqrt(double(xi*xi +yi*yi +zr*zr)))
	+          exp(-alpha_*sqrt(double(xi*xi +yr*yr +zr*zr)))
	+          exp(-alpha_*sqrt(double(xr*xr +yi*yi +zr*zr)))
	+          exp(-alpha_*sqrt(double(xr*xr +yr*yr +zi*zi)))
	+          exp(-alpha_*sqrt(double(xr*xr +yr*yr +zr*zr)));
    }
  }
  double nrm = (smr_*smr_).sum();
  smr_/= Communicator::instance()->reduce_sum(nrm);
}

template<typename FMT> const Field Source_exp<FMT>::
mksrc(int s,int c)const{
  Field src(ff_->size());
  for(int site=0; site<ff_->Nvol(); ++site)
    src.set(ff_->index_r(c,s,site),smr_[site]);
  return src;
}

template<typename FMT> const Field Source_exp<FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}

////// Source_Gauss ----------------------------------------------------
template<typename FMT> class Source_Gauss :public Source{
private:
  std::vector<int> sp_;
  FMT* ff_;
  std::valarray<double> smr_;
  double alpha_;
  void set_src();
public:
  Source_Gauss(std::vector<int> source_position, double al, int Nvol)
    :sp_(source_position),ff_(new FMT(Nvol)),
     alpha_(al),smr_(0.0,ff_->Nvol()){
    set_src();
  }
  ~Source_Gauss(){ delete ff_;}
  const Field mksrc(int s,int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename FMT> void Source_Gauss<FMT>::set_src(){

  int tb = SiteIndex::instance()->global_t(0);
  int tt = SiteIndex::instance()->global_t(SiteIndex::instance()->Bdir(TDIR));

  if(tb <= sp_[TDIR] && sp_[TDIR]<= tt){
    
    int ti = sp_[TDIR]%CommonPrms::instance()->Nt();

    for(int sl= 0; sl<SiteIndex::instance()->slsize(ti,TDIR); ++sl){

      int site = SiteIndex::instance()->slice_t(ti,sl);
      int gsite = SiteIndex::instance()->get_gsite(site);

      int xi = SiteIndex::instance()->g_x(gsite) -sp_[XDIR];
      int yi = SiteIndex::instance()->g_y(gsite) -sp_[YDIR];
      int zi = SiteIndex::instance()->g_z(gsite) -sp_[ZDIR];

      int xr = CommonPrms::instance()->Lx() -xi;
      int yr = CommonPrms::instance()->Ly() -yi;
      int zr = CommonPrms::instance()->Lz() -zi;
     
      smr_[site] = exp(-alpha_*double(xi*xi +yi*yi +zi*zi))
	+          exp(-alpha_*double(xr*xr +yi*yi +zi*zi))
	+          exp(-alpha_*double(xi*xi +yr*yr +zi*zi))
	+          exp(-alpha_*double(xi*xi +yi*yi +zr*zr))
	+          exp(-alpha_*double(xi*xi +yr*yr +zr*zr))
	+          exp(-alpha_*double(xr*xr +yi*yi +zr*zr))
	+          exp(-alpha_*double(xr*xr +yr*yr +zi*zi))
	+          exp(-alpha_*double(xr*xr +yr*yr +zr*zr));
    }
  }
  double nrm = (smr_*smr_).sum();
  smr_/= Communicator::instance()->reduce_sum(nrm);
}

template<typename FMT> const Field Source_Gauss<FMT>::
mksrc(int s,int c)const{
  Field src(ff_->size());
  for(int site=0; site<ff_->Nvol(); ++site)
      src.set(ff_->index_r(c,s,site),smr_[site]);
  return src;
}

template<typename FMT> const Field Source_Gauss<FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}

////// Source_wall ----------------------------------------------------
template<typename FMT> class Source_wall :public Source{
private:
  int spt_;
  FMT* ff_;
  std::valarray<double> smr_;
  void set_src();
public:
  Source_wall(int spt, int Nvol)
    :spt_(spt),ff_(new FMT(Nvol)),smr_(0.0,ff_->Nvol()){ set_src(); }
  ~Source_wall(){ delete ff_;}
  const Field mksrc(int s,int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename FMT> void Source_wall<FMT>::set_src(){
  for(int site=0; site<ff_->Nvol(); ++site){
    int gsite = SiteIndex::instance()->get_gsite(site);
    int t = SiteIndex::instance()->g_t(gsite);
    if(t == spt_) smr_[site] = 1.0;
  }
  double nrm = (smr_*smr_).sum();
  smr_/= Communicator::instance()->reduce_sum(nrm);
}

template<typename FMT> const Field Source_wall<FMT>::
mksrc(int s,int c)const{
  Field src(ff_->size());
  for(int site=0; site<ff_->Nvol(); ++site)
    src.set(ff_->index_r(c,s,site),smr_[site]);
  return Field(src);
}

template<typename FMT> const Field Source_wall<FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}

////// Source_wnoise----------------------------------------------------
template<typename IDX,typename FMT> class Source_wnoise :public Source{
private:
  const RandNum& rand_;
  FMT* ff_;
  IDX* idx_;
  std::valarray<double> src_;
public:
  Source_wnoise(const RandNum& rand, int Nvol)
    :rand_(rand),idx_(IDX::instance()),ff_(new FMT(Nvol)),
     src_(0.0,ff_->size()){
    //
    std::vector<int> gsite = idx_->get_gsite();
    MPrand::mp_get(src_,rand_,gsite,*ff_);
  }

  ~Source_wnoise(){ delete ff_;}
  const Field mksrc(int s,int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename IDX,typename FMT> 
const Field Source_wnoise<IDX,FMT>::mksrc(int s,int c)const{
  Field wns(ff_->size());
  wns.set(ff_->cs_slice(c,s),src_[ff_->cs_slice(c,s)]);
  return wns;
}

template<typename IDX,typename FMT> 
const Field Source_wnoise<IDX,FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}

////// Source_Z2----------------------------------------------------
/*! 
 * @brief Source type: Z2 noise source
 */
enum Z2Type {StandardZ2, ComplexZ2};
Begin_Enum_String( Z2Type ) {
  Enum_String( StandardZ2 );
  Enum_String( ComplexZ2 );
} 
End_Enum_String;

template<typename IDX,typename FMT> class Source_Z2noise :public Source{
private:
  const RandNum& rand_generator_;/*!< @brief Random number generator (to be provided)*/
  Z2Type Type_;
  IDX* idx_;
  FMT* ff_;/*!< @brief %Field %Format specifier */
  std::valarray<double> src_;/*!< @brief %Source field result */
public:
  /*!  @brief Constructor */
  Source_Z2noise(const RandNum& rand,int Nvol,Z2Type Type)
    :rand_generator_(rand),
     idx_(IDX::instance()),
     ff_(new FMT(Nvol)),
     Type_(Type),src_(0.0,ff_->size()){
    CCIO::cout<<"creating Source_Z2noise\n";
    setup_source();
  }

  /*! @brief Destructor */
  ~Source_Z2noise(){ delete ff_;}

  void setup_source(){
    CCIO::cout<<"setting up Source_Z2noise\n";
    std::valarray<double> white_noise(0.0,ff_->size());
    double cosine;
    
    std::vector<int> gsite = idx_->get_gsite();
    MPrand::mp_get(white_noise,rand_generator_,gsite,*ff_);
    
    for (int idx = 0; idx <white_noise.size()/2; ++idx){
      cosine = cos(2.0*PI*white_noise[idx]);
      src_[2*idx] = copysign(1.0,cosine);
      src_[2*idx+1] = 0;	
      if (Type_ == ComplexZ2) {
       	src_[2*idx]   = src_[2*idx]*INVSQRT2;
       	src_[2*idx+1] = src_[2*idx]*INVSQRT2;
      }		   
    }
  }	   	   
  /*! @brief Constructs the source vector for specific spin and colour*/
  const Field mksrc(int s,int c)const;
  const Field mksrc(const std::vector<int>& lv,int s,int c)const;
};

template<typename IDX,typename FMT>
const Field Source_Z2noise<IDX,FMT>::mksrc(int s,int c)const{
  Field wns(ff_->size());
  wns.set(ff_->cs_slice(c,s), src_[ff_->cs_slice(c,s)]);
  return wns;
}

template<typename IDX,typename FMT> 
const Field Source_Z2noise<IDX,FMT>::
mksrc(const std::vector<int>& lv,int s,int c)const{
  return Field(mksrc(s,c)[ff_->get_sub(lv)]);
}



#endif
