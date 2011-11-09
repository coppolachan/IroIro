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
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex_eo* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;

  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_even_up(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){}

  ShiftField_even_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ShiftField_even_up(const std::valarray<double>& field,
		     const FMT* fmt,int dir)
    :field_(&field),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ~ShiftField_even_up(){ delete bdfmt_;}
  
  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;

  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;

  const Field field() const;
};

template<typename FMT> 
void ShiftField_even_up<FMT>::setup(){

  std::vector<int> b(idx_->ebdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_even_up<FMT>::setup(const Field& field,int dir){

  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->ebdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_even_up<FMT>::
setup(const std::valarray<double>& field,int dir){

  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->ebdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
const std::valarray<double>
ShiftField_even_up<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->cslice(in,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_even_up<FMT>::iv(int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->islice(idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ex_p(   hs,dir_),ex)];
} 
 
template<typename FMT> 
double ShiftField_even_up<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_even_up<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->e_cmp(hs,dir_) == idx_->Bdir(dir_));
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_up<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ex_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const Field ShiftField_even_up<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<Nex_; ++ex)
    for(int hs=0; hs<fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return Field(va);
}

/******************** ShiftField_even_dn ************************/
template<typename FMT> class ShiftField_even_dn:public ShiftField{
private:
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex_eo* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;
  
  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_even_dn(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){}

  ShiftField_even_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ShiftField_even_dn(const std::valarray<double>& field,
		     const FMT* fmt,int dir)
    :field_(&field),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ~ShiftField_even_dn(){ delete bdfmt_;}

  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
  const Field field() const;
};

template<typename FMT> 
void ShiftField_even_dn<FMT>::setup(){
  
  std::vector<int> b(idx_->ebdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_even_dn<FMT>::setup(const Field& field,int dir){

  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->ebdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_even_dn<FMT>::
setup(const std::valarray<double>& field,int dir){

  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->ebdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
const std::valarray<double>
ShiftField_even_dn<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->cslice(in,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_even_dn<FMT>::iv(int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->islice(idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->e_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_even_dn<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->e_cmp(hs,dir_) == 0);
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_r(c,s,idx_->ebid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_i(c,s,idx_->ebid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_even_dn<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ex_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const Field ShiftField_even_dn<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<Nex_; ++ex)
    for(int hs=0; hs<fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return Field(va);
}

/******************** ShiftField_odd_up ************************/
template<typename FMT> class ShiftField_odd_up:public ShiftField{
private:
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex_eo* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;

  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_odd_up(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){}
  
  ShiftField_odd_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ShiftField_odd_up(const std::valarray<double>& field,
		     const FMT* fmt,int dir)
    :field_(&field),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ~ShiftField_odd_up(){ delete bdfmt_;}

  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
  const Field field() const;
};

template<typename FMT> 
void ShiftField_odd_up<FMT>::setup(){
  
  std::vector<int> b(idx_->obdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_odd_up<FMT>::setup(const Field& field,int dir){

  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->obdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_odd_up<FMT>::
setup(const std::valarray<double>& field,int dir){

  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->obdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
const std::valarray<double>
ShiftField_odd_up<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->cslice(in,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_odd_up<FMT>::iv(int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->islice(idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ox_p(   hs,dir_),ex)];
}  
 
template<typename FMT> 
double ShiftField_odd_up<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_r(c,s,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_)) 
    return bdry_[     bdfmt_->index_i(c,s,idx_->obid_up(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_odd_up<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->o_cmp(hs,dir_) == idx_->Bdir(dir_));
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_r(c,s,idx_->obid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_i(c,s,idx_->obid_up(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_up<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ox_p(   hs,dir_),ex)];
}  

template<typename FMT> 
const Field ShiftField_odd_up<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<Nex_; ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return Field(va);
}

/******************** ShiftField_odd_dn ************************/
template<typename FMT> class ShiftField_odd_dn:public ShiftField{
private:
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex_eo* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;

  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_odd_dn(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){}
  
  ShiftField_odd_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ShiftField_odd_dn(const std::valarray<double>& field,
		    const FMT* fmt,int dir)
    :field_(&field),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex_eo::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()){setup();}

  ~ShiftField_odd_dn(){ delete bdfmt_;}

  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int hs,int ex=0) const;
  const std::valarray<double> iv(int hs,int ex=0) const;
  double re(int c,int s,int hs,int ex=0) const;
  double im(int c,int s,int hs,int ex=0) const;
  bool on_bdry(int hs,int ex=0) const;
  double re_on_bdry(int c,int s,int hs,int ex=0) const;
  double im_on_bdry(int c,int s,int hs,int ex=0) const;
  double re_on_bulk(int c,int s,int hs,int ex=0) const;
  double im_on_bulk(int c,int s,int hs,int ex=0) const;
  const Field field() const;
};

template<typename FMT> 
void ShiftField_odd_dn<FMT>::setup(){
  
  std::vector<int> b(idx_->obdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_odd_dn<FMT>::setup(const Field& field,int dir){

  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->obdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_odd_dn<FMT>::
setup(const std::valarray<double>& field,int dir){

  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->obdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  bdfmt_= new FMT(bsize,Nex_);
}

template<typename FMT> 
const std::valarray<double>
ShiftField_odd_dn<FMT>::cv(int in,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->cslice(in,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->cslice(in,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const std::valarray<double> 
ShiftField_odd_dn<FMT>::iv(int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->islice(idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->islice(idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_r(c,s,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_r(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im(int c,int s,int hs,int ex) const {
  if(idx_->o_cmp(hs,dir_) == 0)
    return bdry_[     bdfmt_->index_i(c,s,idx_->obid_lw(hs,dir_),ex)];
  else return (*field_)[fmt_->index_i(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
bool ShiftField_odd_dn<FMT>::on_bdry(int hs,int ex) const {
  return (idx_->o_cmp(hs,dir_) == 0);
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_r(c,s,idx_->obid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im_on_bdry(int c,int s,int hs,int ex) const {
  return bdry_[     bdfmt_->index_i(c,s,idx_->obid_lw(hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::re_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_r(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_odd_dn<FMT>::im_on_bulk(int c,int s,int hs,int ex) const {
  return (*field_)[fmt_->index_i(c,s,idx_->ox_m(   hs,dir_),ex)];
}  

template<typename FMT> 
const Field ShiftField_odd_dn<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex<Nex_; ++ex)
    for(int hs=0; hs< fmt_->Nvol(); ++hs)
      va[fmt_->islice(hs,ex)] = iv(hs,ex);

  return Field(va);
}

#endif
