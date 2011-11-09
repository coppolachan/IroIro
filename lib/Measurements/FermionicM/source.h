//----------------------------------------------------------------------
// source.h
//----------------------------------------------------------------------
#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED

#ifndef SITEINDEX_INCLUDED
#include "Main/Geometry/siteIndex.h"
#endif

#ifndef RANDNUM_INCLUDED
#include "Tools/randNum.h"
#endif

#ifndef FIELD_INCLUDED
#include "include/field.h"
#endif

class Source{
public:
  virtual ~Source(){}
  virtual const Field mksrc(int s, int c) = 0;
};

////// Source_local -----------------------------------------------------
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
};

template<typename FMT> 
void Source_local<FMT>::set_src(){

  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nt = CommonPrms::instance()->Nt();
  int Ndim = CommonPrms::instance()->Ndim();

#ifdef USE_MPI
  const int size[Ndim] = { Nx_,Ny_,Nz_,Nt_, };
  int coord[Ndim];
  for (int i=0; i<Ndim; ++i) coord[i] = sp_[i]/size[i];
  int nid = Communicator::instance()->nodeid(Ndim, coord);
#endif
#ifndef USE_MPI
  int nid= Communicator::
    instance()->nodeid(sp_[0]/Nx,sp_[1]/Ny,sp_[2]/Nz,sp_[3]/Nt);
#endif

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
  src_[ff_->cs_slice(c,s)] =smr_;
  return Field(src_);
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
  Source_Gauss(std::vector<int> sp, double al, int Nvol)
    :sp_(sp),ff_(new FMT(Nvol)),alpha_(al){
    src_.resize(ff_->size());
    set_src();
  }
  ~Source_Gauss(){ delete ff_;}
  const Field mksrc(int s, int c);
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
  src_[ff_->cs_slice(c,s)] =smr_;
  return Field(src_);
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
    src_.resize(ff_->size());
  }
  ~Source_wnoise(){ delete ff_;}
  const Field mksrc(int s, int c);
};

template<typename FMT>
const Field Source_wnoise<FMT>::mksrc(int s, int c){
  src_=0.0;

  std::valarray<double> wnz_(ff_->Nvol());
  rand_.get(wnz_);

  //  for(int i =0; i<wnz_.size(); ++i) 
  //    std::cout<<"wnz_["<<i<<"]="<<wnz_[i]<<std::endl;

  src_[ff_->cs_slice(c,s)] =wnz_;
  return Field(src_);
}
#endif
