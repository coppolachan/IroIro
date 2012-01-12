//---------------------------------------------------------------------
// shiftField_eo.h
//---------------------------------------------------------------------
#ifndef SHIFTFIELD_EO_INCLUDED
#define SHIFTFIELD_EO_INCLUDED

#ifndef SITEINDEX_EO_INCLUDED
#include "siteIndex_eo.h"
#endif

#ifndef SHIFTFIELD_INCLUDED
#include "shiftField.h"
#endif

#ifndef FIELD_INCLUDED
#include "field.h"
#endif

/******************** ShiftField_even_up ************************/
template<typename FMT> class ShiftField_even_up:public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;

  Communicator* com_;
  SiteIndex_eo* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  const std::valarray<size_t> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_even_up(const FMT* fmt, int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdlw(dir))),
     bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),
     bdry_(bdfmt_->size()){}

  ShiftField_even_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdlw(dir))),
     bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),
     bdry_(bdfmt_->size()){ 
    com_->transfer_fw(bdry_,(*field_)[bd_],dir);
  }

  ShiftField_even_up(const std::valarray<double>& field,
		     const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdlw(dir))),
     bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),
     bdry_(bdfmt_->size()){ 
    com_->transfer_fw(bdry_,field[bd_],dir);
  }

  ~ShiftField_even_up(){ delete bdfmt_;}
  
  void setf(const Field& field){
    field_= &(field.getva());
    com_->transfer_fw(bdry_,(*field_)[bd_],dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_fw(bdry_,field[bd_],dir_);
  }

  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  const double* get_bdry_addr(int site,int ex=0);
  const double* get_bulk_addr(int site,int ex=0);

  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
};

template<typename FMT> 
const std::valarray<double>
ShiftField_even_up<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->cslice(in,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_even_up<FMT>::iv(int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->islice(idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ox_p(   hs,dir_),ex)];
} 
 
template<typename FMT> 
double ShiftField_even_up<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_r(c,s,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_i(c,s,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_even_up<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_));
}  

template<typename FMT> 
const double* ShiftField_even_up<FMT>::get_bdry_addr(int hs,int ex){
  return &bdry_[bdfmt_->index_r(0,0,idx_->obid_up(hs,dir_),ex)];
}

template<typename FMT> 
const double* ShiftField_even_up<FMT>::get_bulk_addr(int hs,int ex){
  return &(*const_cast<std::valarray<double> *>(field_))[fmt_->index_r(0,0,idx_->ox_p(hs,dir_),ex)];
}


template<typename FMT> 
double ShiftField_even_up<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->obid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->obid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ox_p(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ox_p(hs,dir_),ex)];
}  

template<typename FMT>
const std::valarray<double> ShiftField_even_up<FMT>::getva() const {
 std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<fmt_->Nex(); ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return va;
}

/******************** ShiftField_even_dn ************************/
template<typename FMT> class ShiftField_even_dn:public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;

  Communicator* com_;
  SiteIndex_eo* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  const std::valarray<size_t> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_even_dn(const FMT* fmt, int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdup(dir))),
     bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){}

  ShiftField_even_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdup(dir))),
          bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){ 
    com_->transfer_bk(bdry_,(*field_)[bd_],dir);
  }

  ShiftField_even_dn(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->ebdup(dir))),
     bdfmt_(new FMT(idx_->Ve_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){ 
    com_->transfer_bk(bdry_,field[bd_],dir);
  }

  ~ShiftField_even_dn(){ delete bdfmt_;}

  void setf(const Field& field){
    field_= &(field.getva());
    com_->transfer_bk(bdry_,(*field_)[bd_],dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_bk(bdry_,field[bd_],dir_);
  }

  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  const double* get_bdry_addr(int site,int ex=0);
  const double* get_bulk_addr(int site,int ex=0);

  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
};

template<typename FMT> 
const std::valarray<double>
ShiftField_even_dn<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->cslice(in,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_even_dn<FMT>::iv(int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->islice(idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_r(c,s,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_i(c,s,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_even_dn<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->o_cmp(hs,dir_) == 0);
}  

template<typename FMT> 
const double* ShiftField_even_dn<FMT>::get_bdry_addr(int hs,int ex) {
  return &bdry_[bdfmt_->index_r(0,0,idx_->obid_lw(hs,dir_),ex)];
}

template<typename FMT> 
const double* ShiftField_even_dn<FMT>::get_bulk_addr(int hs,int ex) {
  return &(*const_cast<std::valarray<double> *>(field_))[fmt_->index_r(0,0,idx_->ox_m(hs,dir_),ex)];
}

template<typename FMT> 
double ShiftField_even_dn<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->obid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->obid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ox_m(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ox_m(hs,dir_),ex)];
}  

template<typename FMT>
const std::valarray<double> ShiftField_even_dn<FMT>::getva() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<fmt_->Nex(); ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return va;
}

/******************** ShiftField_odd_up ************************/
template<typename FMT> class ShiftField_odd_up:public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;

  Communicator* com_;
  SiteIndex_eo* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  //  const std::vector<int> bd_;
  const std::valarray<size_t> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_odd_up(const FMT* fmt, int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdlw(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){}

  ShiftField_odd_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdlw(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){ 
    com_->transfer_fw(bdry_,(*field_)[bd_],dir);
  }

  ShiftField_odd_up(const std::valarray<double>& field,
		    const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdlw(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){ 
    com_->transfer_fw(bdry_,field[bd_],dir);
  }

  ~ShiftField_odd_up(){ delete bdfmt_;}

  void setf(const Field& field){
    field_= &(field.getva());
    com_->transfer_fw(bdry_,(*field_)[bd_],dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_fw(bdry_,field[bd_],dir_);
  }

  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  const double* get_bdry_addr(int site,int ex=0);
  const double* get_bulk_addr(int site,int ex=0);
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
};

template<typename FMT> 
const std::valarray<double>
ShiftField_odd_up<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->cslice(in,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_odd_up<FMT>::iv(int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->islice(idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ex_p(   hs,dir_),ex)];
}  
 
template<typename FMT> 
double ShiftField_odd_up<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_odd_up<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_));
}  

template<typename FMT> 
const double* ShiftField_odd_up<FMT>::get_bdry_addr(int hs,int ex){
  return &bdry_[bdfmt_->index_r(0,0,idx_->ebid_up(hs,dir_),ex)];
}

template<typename FMT> 
const double* ShiftField_odd_up<FMT>::get_bulk_addr(int hs,int ex){
  return &(*const_cast<std::valarray<double> *>(field_))[fmt_->index_r(0,0,idx_->ex_p(hs,dir_),ex)];
}


template<typename FMT> 
double ShiftField_odd_up<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->ebid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->ebid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ex_p(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ex_p(hs,dir_),ex)];
}  

template<typename FMT>
const std::valarray<double> ShiftField_odd_up<FMT>::getva() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<fmt_->Nex(); ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return va;
}

/******************** ShiftField_odd_dn ************************/
template<typename FMT> class ShiftField_odd_dn:public ShiftField{
private:
  const FMT* fmt_;
  const int dir_;

  Communicator* com_;
  SiteIndex_eo* idx_;

  const FMT* bdfmt_;
  const std::valarray<double>* field_;
  const std::valarray<size_t> bd_;
  std::valarray<double> bdry_;

public:
  ShiftField_odd_dn(const FMT* fmt, int dir)
    :field_(0),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdup(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){}

  ShiftField_odd_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdup(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){
    com_->transfer_bk(bdry_,(*field_)[bd_],dir);
  }

  ShiftField_odd_dn(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),fmt_(fmt),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     bd_(fmt_->get_sub(idx_->obdup(dir))),
     bdfmt_(new FMT(idx_->Vo_dir(dir),fmt_->Nex())),    
     bdry_(bdfmt_->size()){ 
    com_->transfer_bk(bdry_,field[bd_],dir);
  }

  ~ShiftField_odd_dn(){ delete bdfmt_;}

  void setf(const Field& field){
    field_= &(field.getva());
    com_->transfer_bk(bdry_,(*field_)[bd_],dir_);
  }

  void setf(const std::valarray<double>& field){
    field_= &field;
    com_->transfer_bk(bdry_,field[bd_],dir_);
  }

  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  const std::valarray<double> getva() const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  const double* get_bdry_addr(int site,int ex=0);
  const double* get_bulk_addr(int site,int ex=0);
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
  const Field field() const;
};

template<typename FMT> 
const std::valarray<double>
ShiftField_odd_dn<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->cslice(in,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_odd_dn<FMT>::iv(int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->islice(idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_odd_dn<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->e_cmp(hs,dir_) == 0);
}  

template<typename FMT> 
const double* ShiftField_odd_dn<FMT>::get_bdry_addr(int hs,int ex) {
  return &bdry_[bdfmt_->index_r(0,0,idx_->ebid_lw(hs,dir_),ex)];
}

template<typename FMT> 
const double* ShiftField_odd_dn<FMT>::get_bulk_addr(int hs,int ex) {
  return &(*const_cast<std::valarray<double> *>(field_))[fmt_->index_r(0,0,idx_->ex_m(hs,dir_),ex)];
}

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_r(c,s,idx_->ebid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[bdfmt_->index_i(c,s,idx_->ebid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ex_m(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ex_m(hs,dir_),ex)];
}  

template<typename FMT>
const std::valarray<double> ShiftField_odd_dn<FMT>::getva() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<fmt_->Nex(); ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return va;
}

#endif
