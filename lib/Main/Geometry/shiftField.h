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
  virtual void setup(const Field& field,int dir) =0;
  virtual void setup(const std::valarray<double>& field,int dir) =0;
  virtual const std::valarray<double> cv(int in,int site,int ex=0)const =0;
  virtual const std::valarray<double> iv(int site, int ex=0)const =0;
  virtual double re(int c,int s,int site,int ex=0)const =0;
  virtual double im(int c,int s,int site,int ex=0)const =0;
  virtual bool on_bdry(int site,int ex=0)const =0;
  virtual double re_on_bdry(int c,int s,int site,int ex=0)const =0;
  virtual double im_on_bdry(int c,int s,int site,int ex=0)const =0;
  virtual double re_on_bulk(int c,int s,int site,int ex=0)const =0;
  virtual double im_on_bulk(int c,int s,int site,int ex=0)const =0;
  virtual const Field field()const =0;
};

// ShiftField_up  ------------------------
template<typename FMT> class ShiftField_up :public ShiftField{
private:
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;

  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_up(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()),
     bdfmt_(new FMT){}

  ShiftField_up(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),
     dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),
     Nin_(fmt->Nin()),
     Nex_(fmt->Nex()),
     bdfmt_(new FMT){setup();}

  ShiftField_up(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()),
     bdfmt_(new FMT){setup();}

  ~ShiftField_up(){delete bdfmt_;}

  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int site,int ex=0) const;
  const std::valarray<double> iv(int site,int ex=0) const;
  double re(int c,int s,int site,int ex=0) const;
  double im(int c,int s,int site,int ex=0) const;
  bool on_bdry(int site,int ex=0) const;
  double re_on_bdry(int c,int s,int site,int ex=0) const;
  double im_on_bdry(int c,int s,int site,int ex=0) const;
  double re_on_bulk(int c,int s,int site,int ex=0) const;
  double im_on_bulk(int c,int s,int site,int ex=0) const;
  const Field field() const;
};

template<typename FMT> 
void ShiftField_up<FMT>::setup(){
  
  std::vector<int> b(idx_->bdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_up<FMT>::setup(const Field& field,int dir){
  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->bdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_up<FMT>::setup(const std::valarray<double>& field,int dir){
  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->bdlw(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_fw(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

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
const Field ShiftField_up<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex< Nex_; ++ex)
    for(int site=0; site< fmt_->Nvol(); ++site)
      va[fmt_->islice(site,ex)] = iv(site,ex);
  
  return Field(va);
}

// ShiftField_dn  ------------------------
template<typename FMT> class ShiftField_dn :public ShiftField{
private:
  const std::valarray<double>* field_;
  Communicator* com_;
  SiteIndex* idx_;
  const FMT* fmt_;
  int Nin_;
  int Nex_;
  int dir_;  

  std::valarray<double> bdry_;
  FMT* bdfmt_;
  void setup();
public:
  explicit ShiftField_dn(const FMT* fmt)
    :field_(0),dir_(-1),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()),
     bdfmt_(new FMT){}

  ShiftField_dn(const Field& field,const FMT* fmt,int dir)
    :field_(&(field.getva())),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()),
     bdfmt_(new FMT){setup();}

  ShiftField_dn(const std::valarray<double>& field,const FMT* fmt,int dir)
    :field_(&field),dir_(dir),
     com_(Communicator::instance()),
     idx_(SiteIndex::instance()),
     fmt_(fmt),Nin_(fmt->Nin()),Nex_(fmt->Nex()),
     bdfmt_(new FMT){setup();}

  ~ShiftField_dn(){delete bdfmt_;}
  
  void setup(const Field& field,int dir);
  void setup(const std::valarray<double>& field,int dir);
  const std::valarray<double> cv(int in,int site,int ex=0) const;
  const std::valarray<double> iv(int site,int ex=0) const;
  double re(int c,int s,int site,int ex=0) const;
  double im(int c,int s,int site,int ex=0) const;
  bool on_bdry(int site,int ex=0) const;
  double re_on_bdry(int c,int s,int site,int ex=0) const;
  double im_on_bdry(int c,int s,int site,int ex=0) const;
  double re_on_bulk(int c,int s,int site,int ex=0) const;
  double im_on_bulk(int c,int s,int site,int ex=0) const;
  const Field field() const;
};

template<typename FMT>
void ShiftField_dn<FMT>::setup(){

  std::vector<int> b(idx_->bdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_dn<FMT>::setup(const Field& field,int dir){
  field_= &(field.getva());
  dir_= dir;
  std::vector<int> b(idx_->bdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

template<typename FMT> 
void ShiftField_dn<FMT>::setup(const std::valarray<double>& field,int dir){
  field_= &field;
  dir_= dir;
  std::vector<int> b(idx_->bdup(dir_));
  int bsize = b.size();
  int size = Nin_*Nex_*bsize;

  std::valarray<size_t> bd(size);
  int k=0;
  for(int ex=0; ex<Nex_; ++ex)
    for(int ss=0; ss<bsize; ++ss)
      for(int in=0; in<Nin_; ++in) bd[k++] =fmt_->index(in,b[ss],ex);
  
  bdry_.resize(size);
  com_->transfer_bk(bdry_,(*field_)[bd],size,dir_);
  //bdfmt_= new FMT(bsize,Nex_);
  bdfmt_->setup(bsize,Nex_);
}

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
  return bdry_[     bdfmt_->index_r(c,s,idx_->x_b(site,dir_),ex)];
}  

template<typename FMT> 
double ShiftField_dn<FMT>::im_on_bdry(int c,int s,int site,int ex) const {
  return bdry_[     bdfmt_->index_i(c,s,idx_->x_b(site,dir_),ex)];
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
const Field ShiftField_dn<FMT>::field() const {
  std::valarray<double> va(fmt_->size());

  for(int ex=0; ex< Nex_; ++ex)
    for(int site=0; site< fmt_->Nvol(); ++site)
      va[fmt_->islice(site,ex)] = iv(site,ex);
  
  return Field(va);
}
#endif
