//---------------------------------------------------------------------
// shiftField.h
//---------------------------------------------------------------------
#ifndef SHIFTFIELD_INCLUDED
#define SHIFTFIELD_INCLUDED

#ifndef COMMOMPRMS_INCLUDED
#include "include/commonPrms.h"
#endif

#ifndef SITEINDEX_INCLUDED
#include "siteIndex.h"
#endif

#ifndef FIELD_INCLUDED
#include "include/field.h"
#endif

class ShiftField{
public:
  virtual ~ShiftField(){}
  virtual const std::valarray<double> cv(int in,int site,int ex=0)const =0;
  virtual const std::valarray<double> iv(int site, int ex=0)const =0;
  virtual const std::valarray<double> getva()const =0;
  virtual double re(int c,int s,int site,int ex=0)const =0;
  virtual double im(int c,int s,int site,int ex=0)const =0;
  virtual bool on_bdry(int site,int ex=0)const =0;
  virtual double re_on_bdry(int c,int s,int site,int ex=0)const =0;
  virtual double im_on_bdry(int c,int s,int site,int ex=0)const =0;
  virtual double re_on_bulk(int c,int s,int site,int ex=0)const =0;
  virtual double im_on_bulk(int c,int s,int site,int ex=0)const =0;
  virtual void setf(const Field& field) =0; 
  virtual void setf(const std::valarray<double>& field) =0; 
};

// ShiftField_up  ------------------------
template<typename FMT> class ShiftField_up :public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;
  
  Communicator* com_;
  SiteIndex* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  const std::vector<int> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_up(const FMT* fmt,int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bd_(fmt_->get_sub(idx_->bdlw(dir))),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bdry_(bdfmt_->size()){}

  ShiftField_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bd_(fmt_->get_sub(idx_->bdlw(dir))),
     bdry_(bdfmt_->size()){ com_->transfer_fw(bdry_,*field_,bd_,dir);}

  ShiftField_up(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bd_(fmt_->get_sub(idx_->bdlw(dir))),
     bdry_(bdfmt_->size()){ com_->transfer_fw(bdry_,field,bd_,dir);}

  ~ShiftField_up(){delete bdfmt_;}

  void setf(const Field& field){ 
    field_= &(field.getva());
    com_->transfer_fw(bdry_,*field_,bd_,dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_fw(bdry_,field,bd_,dir_);
  }

  const std::valarray<double> cv(int in,int site,int ex=0) const;
  const std::valarray<double> iv(int site,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int site,int ex=0) const;
  double im(int c,int s,int site,int ex=0) const;
  bool on_bdry(int site,int ex=0) const;
  double re_on_bdry(int c,int s,int site,int ex=0) const;
  double im_on_bdry(int c,int s,int site,int ex=0) const;
  double re_on_bulk(int c,int s,int site,int ex=0) const;
  double im_on_bulk(int c,int s,int site,int ex=0) const;
};

template<typename FMT> 
const std::valarray<double>
ShiftField_up<FMT>::cv(int in,int site,int ex) const {
  if(idx_->cmp(site,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->cslice(in,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_up<FMT>::iv(int site,int ex) const {
  if(idx_->cmp(site,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->islice(idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_up<FMT>::on_bdry(int site,int ex) const {
  return (idx_->cmp(site,dir_) == idx_->Bdir(dir_));
}  

template<typename FMT> 
double ShiftField_up<FMT>::re(int c,int s,int site,int ex) const {
  if(idx_->cmp(site,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_r(c,s,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_up<FMT>::im(int c,int s,int site,int ex) const {
  if(idx_->cmp(site,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_i(c,s,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_up<FMT>::re_on_bdry(int c,int s,int site,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->x_b(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_up<FMT>::im_on_bdry(int c,int s,int site,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->x_b(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_up<FMT>::re_on_bulk(int c,int s,int site,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_up<FMT>::im_on_bulk(int c,int s,int site,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->x_p(site,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> ShiftField_up<FMT>::getva() const {

  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex< fmt_->Nex(); ++ex)
    for(int site=0; site< fmt_->Nvol(); ++site)
      va[fmt_->islice(site,ex)] = iv(site,ex);
  
  return va;
}

// ShiftField_dn  ------------------------
template<typename FMT> class ShiftField_dn :public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;  

  Communicator* com_;
  SiteIndex* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  const std::vector<int> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_dn(const FMT* fmt, int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bd_(fmt_->get_sub(idx_->bdup(dir))),
     bdry_(bdfmt_->size()){}

  ShiftField_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bd_(fmt_->get_sub(idx_->bdup(dir))),
     bdry_(bdfmt_->size()){ com_->transfer_bk(bdry_,*field_,bd_,dir);}

  ShiftField_dn(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     bdfmt_(new FMT(idx_->Vdir(dir),fmt_->Nex())),
     bd_(fmt_->get_sub(idx_->bdup(dir))),
     bdry_(bdfmt_->size()){ com_->transfer_bk(bdry_,field,bd_,dir);}

  ~ShiftField_dn(){delete bdfmt_;}

  void setf(const Field& field){ 
    field_= &(field.getva());
    com_->transfer_bk(bdry_,*field_,bd_,dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_bk(bdry_,field,bd_,dir_);
  }
  
  const std::valarray<double> cv(int in,int site,int ex=0) const;
  const std::valarray<double> iv(int site,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int site,int ex=0) const;
  double im(int c,int s,int site,int ex=0) const;
  bool on_bdry(int site,int ex=0) const;
  double re_on_bdry(int c,int s,int site,int ex=0) const;
  double im_on_bdry(int c,int s,int site,int ex=0) const;
  double re_on_bulk(int c,int s,int site,int ex=0) const;
  double im_on_bulk(int c,int s,int site,int ex=0) const;
};

template<typename FMT>
const std::valarray<double>
ShiftField_dn<FMT>::cv(int in,int site,int ex) const {
  if(idx_->cmp(site,dir_) == 0) 
    return bdry_[     bdfmt_->cslice(in,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->x_m(site,dir_),ex)];
}  

template<typename FMT>
const std::valarray<double> 
ShiftField_dn<FMT>::iv(int site,int ex) const {
  if(idx_->cmp(site,dir_) == 0)
    return bdry_[     bdfmt_->islice(idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->x_m(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::re(int c,int s,int site,int ex) const {
  if(idx_->cmp(site,dir_) == 0)
    return bdry_[     bdfmt_->index_r(c,s,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->x_m(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::im(int c,int s,int site,int ex) const {
  if(idx_->cmp(site,dir_) == 0)
    return bdry_[     bdfmt_->index_i(c,s,idx_->x_b(site,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->x_m(site,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_dn<FMT>::on_bdry(int site,int ex) const {
  return (idx_->cmp(site,dir_) == 0);
}  

template<typename FMT> 
double ShiftField_dn<FMT>::re_on_bdry(int c,int s,int site,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->x_b(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::im_on_bdry(int c,int s,int site,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->x_b(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::re_on_bulk(int c,int s,int site,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->x_m(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::im_on_bulk(int c,int s,int site,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->x_m(site,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> ShiftField_dn<FMT>::getva() const {

  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex< fmt_->Nex(); ++ex)
    for(int site=0; site< fmt_->Nvol(); ++site)
      va[fmt_->islice(site,ex)] = iv(site,ex);
  
  return va;
}
#endif
