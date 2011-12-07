/*!
 *
 * @file source_types.hpp
 *
 * @brief Definition of source types for propagator calculations
 *
 * Include this file if you need to specify explicitly the source type
 */

#include "Main/Geometry/siteIndex.h"
#include "Tools/randNum.h"
#include "Tools/randNum_MP.h"
#include "source.hpp"

#include <cmath>

#include <iostream>

#define INVSQRT2 0.7071067811865475
////// Source_local -------------------------------------------------------

template<typename FMT> class Source_local :public Source{
private:
  std::vector<int> sp_;
  std::valarray<double> src_;
  FMT* ff_;
  bool has_source_;
  int loc_;
  void set_src();
public:
  Source_local(std::vector<int> sp, int Nvol)
    :sp_(sp),ff_(new FMT(Nvol)),has_source_(false),loc_(0){
    src_.resize(ff_->size());
    set_src();
  }
  ~Source_local(){ delete ff_;}
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT> 
void Source_local<FMT>::set_src(){

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  int Ndim = CommonPrms::instance()->Ndim();

  int nid= Communicator::
    instance()->nodeid(sp_[0]/Nx,sp_[1]/Ny,sp_[2]/Nz,sp_[3]/Nt);

  if(Communicator::instance()->nodeid() == nid){
    has_source_= true;
    loc_= SiteIndex::instance()->site(sp_[0]%Nx,sp_[1]%Ny,
				      sp_[2]%Nz,sp_[3]%Nt);
  }
}

template<typename FMT> 
const Field Source_local<FMT>::mksrc(int s,int c){
  src_=0.0;
  if(has_source_) src_[ff_->index_r(c,s,loc_)] =1.0;
  return Field(src_);
}

template<typename FMT> 
const Field Source_local<FMT>::mksrc(const std::vector<int>& lv,int s,int c){
  Field fsrc = mksrc(s,c);
  return Field(fsrc[ff_->get_sub(lv)]);
}

////// Source_exp   -----------------------------------------------------
template<typename FMT> class Source_exp :public Source{
private:
  std::vector<int> sp_;
  FMT* ff_;
  std::valarray<double> src_;
  static std::valarray<double> smr_;
  double alpha_;
  void set_src();
public:
  Source_exp(std::vector<int> sp, double al, int Nvol)
    :sp_(sp),ff_(new FMT(Nvol)),alpha_(al){
    src_.resize(ff_->size());    
    set_src();
  }
  ~Source_exp(){ delete ff_;}
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT>
std::valarray<double> Source_exp<FMT>::smr_;

template<typename FMT>
void Source_exp<FMT>::set_src(){
  smr_=0.0;

  int NPx =(Communicator::instance()->ipe(0))*CommonPrms::instance()->NPEx();
  int NPy =(Communicator::instance()->ipe(1))*CommonPrms::instance()->NPEy();
  int NPz =(Communicator::instance()->ipe(2))*CommonPrms::instance()->NPEz();
  int NPt =(Communicator::instance()->ipe(3))*CommonPrms::instance()->NPEt();

  for(int site= 0; site< ff_->Nvol(); ++site){
    int x = NPx+SiteIndex::instance()->x(site);
    int y = NPy+SiteIndex::instance()->y(site);
    int z = NPz+SiteIndex::instance()->z(site);
    int t = NPt+SiteIndex::instance()->t(site);

    double r = sqrt((x-sp_[0])*(x-sp_[0]) +(y-sp_[1])*(y-sp_[1])
		    +(z-sp_[2])*(z-sp_[2]) +(t-sp_[3])*(t-sp_[3]));
    smr_[site]= exp(-alpha_*r);
  }
}

template<typename FMT> 
const Field Source_exp<FMT>::mksrc(int s, int c){
  src_=0.0;
  for(int site=0; site< ff_->Nvol(); ++site) 
    src_[ff_->index_r(c,s,site)] =smr_[site];
  return Field(src_);
}

template<typename FMT> 
const Field Source_exp<FMT>::mksrc(const std::vector<int>& lv,int s, int c){
  Field fsrc = mksrc(s,c);
  return Field(fsrc[ff_->get_sub(lv)]);
}

////// Source_Gauss ----------------------------------------------------
template<typename FMT> class Source_Gauss :public Source{
private:
  std::vector<int> sp_;
  FMT* ff_;
  std::valarray<double> src_;
  static std::valarray<double> smr_;
  double alpha_;
  void set_src();
public:
  Source_Gauss(std::vector<int> source_position, double al, int Nvol)
    :sp_(source_position),ff_(new FMT(Nvol)),alpha_(al){
    src_.resize(ff_->size());
    set_src();
  }
  ~Source_Gauss(){ delete ff_;}
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT>
void Source_Gauss<FMT>::set_src(){
  smr_.resize(ff_->Nvol(),0.0);
  int NPx =(Communicator::instance()->ipe(0))*CommonPrms::instance()->NPEx();
  int NPy =(Communicator::instance()->ipe(1))*CommonPrms::instance()->NPEy();
  int NPz =(Communicator::instance()->ipe(2))*CommonPrms::instance()->NPEz();
  int NPt =(Communicator::instance()->ipe(3))*CommonPrms::instance()->NPEt();

  for(int site= 0; site< ff_->Nvol(); ++site){
    int x = NPx+SiteIndex::instance()->x(site);
    int y = NPy+SiteIndex::instance()->y(site);
    int z = NPz+SiteIndex::instance()->z(site);
    int t = NPt+SiteIndex::instance()->t(site);

    double rsq = (x-sp_[0])*(x-sp_[0]) +(y-sp_[1])*(y-sp_[1])
                +(z-sp_[2])*(z-sp_[2]) +(t-sp_[3])*(t-sp_[3]);
    smr_[site]= exp(-alpha_*rsq);
  }
}

template<typename FMT>
std::valarray<double> Source_Gauss<FMT>::smr_;

template<typename FMT> 
const Field Source_Gauss<FMT>::mksrc(int s, int c){
  src_=0.0;
  src_[ff_->cs_slice(c,s)] =smr_;
  return Field(src_);
}
template<typename FMT> 
const Field Source_Gauss<FMT>::mksrc(const std::vector<int>& lv,int s,int c){
  Field fsrc = mksrc(s,c);
  return Field(fsrc[ff_->get_sub(lv)]);
}

////// Source_wall ----------------------------------------------------
template<typename FMT> class Source_wall :public Source{
private:
  int spt_;
  FMT* ff_;
  std::valarray<double> src_;
  static std::valarray<double> smr_;
  void set_src();
public:
  Source_wall(int spt, int Nvol)
    :spt_(spt),ff_(new FMT(Nvol)){
    src_.resize(ff_->size());
    set_src();
  }
  ~Source_wall(){ delete ff_;}
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT>
void Source_wall<FMT>::set_src(){
  smr_=0.0;

  int NPt =(Communicator::instance()->ipe(3))*CommonPrms::instance()->NPEt();
  for(int site= 0; site< ff_->Nvol(); ++site){
    int t = NPt+SiteIndex::instance()->t(site);
    if(t == spt_) smr_[site]= 1.0;
  }
}

template<typename FMT>
std::valarray<double> Source_wall<FMT>::smr_;

template<typename FMT>
const Field Source_wall<FMT>::mksrc(int s, int c){
  src_=0.0;
  for(int site=0; site< ff_->Nvol(); ++site)
    src_[ff_->index_r(c,s,site)] =smr_[site];

  return Field(src_);
}

template<typename FMT> 
const Field Source_wall<FMT>::mksrc(const std::vector<int>& lv,int s,int c){
  Field fsrc = mksrc(s,c);
  return Field(fsrc[ff_->get_sub(lv)]);
}

////// Source_wnoise----------------------------------------------------
template<typename FMT> class Source_wnoise :public Source{
private:
  const RandNum& rand_;
  FMT* ff_;
  std::valarray<double> src_;
public:
  Source_wnoise(const RandNum& rand, int Nvol)
    :rand_(rand),ff_(new FMT(Nvol)){
    //
    src_.resize(ff_->size());
    MPrand::mp_get(src_,rand_,*ff_);
  }

  ~Source_wnoise(){ delete ff_;}
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT>
const Field Source_wnoise<FMT>::mksrc(int s, int c){
  Field wns(ff_->size());
  wns.set(ff_->cs_slice(c,s), src_[ff_->cs_slice(c,s)]);
  return wns;
}

template<typename FMT> 
const Field Source_wnoise<FMT>::
mksrc(const std::vector<int>& lv,int s,int c){
  Field scs = mksrc(s,c);
  return Field(scs[ff_->get_sub(lv)]);
}

////// Source_Z2----------------------------------------------------
/*! 
 * @brief Source type: Z2 noise source
 */
enum Z2Type {StandardZ2, ComplexZ2};

template<typename FMT> class Source_Z2noise :public Source{
private:
  const RandNum& rand_generator_;/*!< @brief Random number generator (to be provided)*/
  Z2Type Type_;
  FMT* ff_;/*!< @brief %Field %Format specifier */
  std::valarray<double> source_;/*!< @brief %Source field result */
public:
  /*! 
   * @brief Constructor
   */
  Source_Z2noise(const RandNum& rand, int Nvol, Z2Type Type)
    :rand_generator_(rand),
     ff_(new FMT(Nvol)),
     Type_(Type){
    setup_source();
  }

  /*! 
   * @brief Destructor
   */
  ~Source_Z2noise(){ delete ff_;}

  const void setup_source(){
    std::valarray<double> white_noise(ff_->size());
    double cosine;
    source_.resize(ff_->size());
    MPrand::mp_get(white_noise,rand_generator_,*ff_);
    
    for (int idx = 0; idx <white_noise.size()/2; ++idx){
      cosine = cos(white_noise[idx]);
      source_[idx] = copysign(1.0,cosine);
      source_[2*idx+1] = 0;	
      if (Type_ == ComplexZ2) {
       	source_[2*idx]   = source_[2*idx]*INVSQRT2;
       	source_[2*idx+1] = source_[2*idx];
      }		   
    }
  }	   
  
	   
  /*!
   * @brief Constructs the source vector for specific spin and colour
   */
  const Field mksrc(int s, int c);
  const Field mksrc(const std::vector<int>& lv, int s, int c);
};

template<typename FMT>
const Field Source_Z2noise<FMT>::mksrc(int s, int c){
  Field wns(ff_->size());
  wns.set(ff_->cs_slice(c,s), source_[ff_->cs_slice(c,s)]);
  return wns;
}

template<typename FMT> 
const Field Source_Z2noise<FMT>::
mksrc(const std::vector<int>& lv,int s,int c){
  Field scs = mksrc(s,c);
  return Field(scs[ff_->get_sub(lv)]);
}

