//---------------------------------------------------------------------
// field.h
//---------------------------------------------------------------------
#ifndef FIELD_INCLUDED
#define FIELD_INCLUDED

#include <valarray>
#include <vector>
#include <iostream>

#include "Communicator/communicator.h"
#include <fstream>

class Field{
private:
  std::valarray<double> field_;
public:
  Field(){}
  explicit Field(const std::valarray<double>& va):field_(va){}
  Field(std::size_t size, double val=0.0):field_(val,size){}
  Field(const Field& f):field_(f.field_){}

  std::size_t size() const {return field_.size();}
  void resize(std::size_t size, double val=0.0){field_.resize(size,val);}

  double operator[](const size_t i) const { return field_[i]; }
  const double& operator[](const size_t i) { return field_[i]; }

  std::valarray<double> operator[](const std::slice& sl) const { 
    return field_[sl];}
  std::slice_array<double> operator[](const std::slice& sl) {
    return field_[sl];}
  std::valarray<double> operator[](const std::gslice& sl) const{ 
    return field_[sl];}
  std::gslice_array<double> operator[](const std::gslice& sl){
    return field_[sl];}
  std::valarray<double> operator[](const std::valarray<size_t>& va) const{
    return field_[va];}
  std::indirect_array<double> operator[](const std::valarray<size_t>& va){
    return field_[va];}

  const std::valarray<double>& getva()const{ return field_;}

  void set(const std::size_t i, double val){field_[i]= val;}
  void set(const std::slice& sl, const std::valarray<double>& va){ 
    field_[sl]= va;}
  void set(const std::gslice& sl, const std::valarray<double>& va){ 
    field_[sl]= va;}
  void set(const std::valarray<std::size_t>& idx, 
	   const std::valarray<double>& va){ 
    field_[idx]= va;
  }
  void add(const std::size_t i, double val){ field_[i]+= val;}
  void add(const std::slice& sl, const std::valarray<double>& va){
    field_[sl]+= va;}
  void add(const std::gslice& sl, const std::valarray<double>& va){ 
    field_[sl]+= va;}
  void add(const std::valarray<std::size_t>& idx, 
	   const std::valarray<double>& va){ 
    field_[idx]+= va;
  }

  double norm() const {
    double a = (field_*field_).sum();
    double b = Communicator::instance()->reduce_sum(a);
    return sqrt(b);
  }

  double average_abs() const{
    double a = 0.0;
    for(int i=0; i<field_.size()/2; ++i) 
      a += sqrt(field_[2*i]*field_[2*i] +field_[2*i+1]*field_[2*i+1]);
    
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field_.size()/2);
    return b/c;
  }

  double average_real() const{
    double a = 0.0;
    for(int i=0; i<field_.size()/2; ++i) a += field_[2*i];
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field_.size()/2);
    return b/c;
  }

  double average_imag() const{
    double a = 0.0;
    for(int i=0; i<field_.size()/2; ++i) a += field_[2*i+1];
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field_.size()/2);
    return b/c;
  }

  double max_element() const{
    std::valarray<double> fab = field_*field_;
    double max= fab[0]+fab[1];
    int idx = 0;
    for(int i=1; i< field_.size()/2; ++i){

      double tmp = fab[2*i] +fab[2*i+1];
      if(max < tmp){
	max = tmp;
	idx = i;
      }
    }	
    int mn = Communicator::instance()->reduce_max(max,idx,field_.size()/2);
    return sqrt(max);
  }

  double min_element() const{
    std::valarray<double> fab = field_*field_;
    double min= fab[0]+fab[1];
    int idx = 0;
    for(int i=1; i< field_.size()/2; ++i){

      double tmp = fab[2*i] +fab[2*i+1];
      if(min > tmp){
	min = tmp;
	idx = i;
      }
    }	
    int mn = Communicator::instance()->reduce_min(min,idx,field_.size()/2);
    return sqrt(min);
  }
  
  template<typename T>
  Field& operator=(const T& rhs){
    *this = rhs.eval();
    return *this;
  }
  template<typename T>
  Field& operator+=(const T& rhs){
    *this += rhs.eval();
    return *this;
  }
  template<typename T>
  Field& operator-=(const T& rhs){
    *this -= rhs.eval();
     return *this;
  }

  Field& operator-();
  Field& operator=(const Field&);
  Field& operator=(const double);
  Field& operator=(const std::valarray<double>&);

  Field& operator+=(const Field&);
  Field& operator+=(const std::valarray<double>&);

  Field& operator-=(const Field&);
  Field& operator-=(const std::valarray<double>&);

  Field& operator*=(const double);
  Field& operator/=(const double);

  double operator*(const Field&) const;
  double operator*(const std::valarray<double>&) const;

  double im_prod(const Field&) const;
  double im_prod(const std::valarray<double>&) const;

  void write_stream(std::ofstream& out) {
    out.write((char*)&field_[0], sizeof(double)*field_.size());
  }
  void read_stream(std::ifstream& in) {
    in.read((char*)&field_[0], sizeof(double)*field_.size());
  }
};

inline Field& Field::operator-(){
  field_= -field_;
  return *this;
}
inline Field& Field::operator=(const Field& rhs){
  field_= rhs.field_;
  return *this;
}
inline Field& Field::operator=(const double r){
  field_= r;
  return *this;
}
inline Field& Field::operator=(const std::valarray<double>& rhs){
  field_= rhs;
  return *this;
}

inline Field& Field::operator+=(const Field& rhs){
  field_+= rhs.field_;
  return *this;
}
inline Field& Field::operator+=(const std::valarray<double>& rhs){
  field_+= rhs;
  return *this;
}

inline Field& Field::operator-=(const Field& rhs){
  field_-= rhs.field_;
  return *this;
}
inline Field& Field::operator-=(const std::valarray<double>& rhs){
  field_-= rhs;
  return *this;
}

inline Field& Field::operator*=(const double rhs){
  field_*= rhs;
  return *this;
}

inline Field& Field::operator/=(const double rhs){
  field_/= rhs;
  return *this;
}

inline double Field::operator*(const Field& rhs) const{
  double a = (field_*rhs.field_).sum();
  //double b = Communicator::reduce_sum(a);
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}
inline double Field::operator*(const std::valarray<double>& rhs) const{
  double a = (field_*rhs).sum();
  //double b = Communicator::reduce_sum(a);
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}

inline double Field::im_prod(const Field& rhs) const{
  std::slice re(0,field_.size()/2,2);
  std::slice im(1,field_.size()/2,2);
  
  std::valarray<double> lhs_im = field_[im];
  std::valarray<double> lhs_re = field_[re];

  double a = (lhs_re*rhs[im]).sum()-(lhs_im*rhs[re]).sum();
  //double b = Communicator::reduce_sum(a);
  double b = Communicator::instance()->reduce_sum(a);
  return b;
}
inline double Field::im_prod(const std::valarray<double>& rhs) const{
  std::slice re(0,field_.size()/2,2);
  std::slice im(1,field_.size()/2,2);

  std::valarray<double> lhs_im = field_[im];
  std::valarray<double> lhs_re = field_[re];

  double a = (lhs_re*rhs[im]).sum()-(lhs_im*rhs[re]).sum();
  //double b =  Communicator::reduce_sum(a);
  double b =  Communicator::instance()->reduce_sum(a);
  return b;
}

namespace FieldExpression{
//////// expression templates ///////
  struct Add{
    static Field calc(const Field& lhs, const Field& rhs){
      return Field(lhs)+=rhs;
    }
  };
  
  struct Sub{
    static Field calc(const Field& lhs, const Field& rhs){
      return Field(lhs)-=rhs;
    }
  };

  struct Mul{
    static Field calc(const Field& lhs, const double& rhs){
      return Field(lhs)*=rhs;
    }
    static Field calc(const double& lhs, const Field& rhs){
      return Field(rhs)*=lhs;
    }
  };
  
  template<typename L, typename Op, typename R>
  class AdSbMc{
  private:  
    const L& lhs;
    const R& rhs;
  public:
    AdSbMc(const L& Lhs,const R& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Op::calc(lhs.eval(),rhs.eval());}
  };

  template<typename L, typename Op>
  class AdSbMc<L,Op,Field>{
  private:  
    const L& lhs;
    const Field& rhs;
  public:
    AdSbMc(const L& Lhs,const Field& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Op::calc(lhs.eval(),rhs);}
  };

  template<typename Op, typename R>
  class AdSbMc<Field,Op,R>{
  private:  
    const Field& lhs;
    const R& rhs;
  public:
    AdSbMc(const Field& Lhs,const R& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Op::calc(lhs,rhs.eval());}
  };
  
  template<typename L>
  class AdSbMc<L,Mul,double>{
  private:  
    const L& lhs;
    const double& rhs;
  public:
    AdSbMc(const L& Lhs,const double& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Mul::calc(lhs.eval(),rhs);}
  };
  
  template<>
  class AdSbMc<Field,Mul,double>{
  private:  
    const Field& lhs;
    const double& rhs;
  public:
    AdSbMc(const Field& Lhs,const double& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Mul::calc(lhs,rhs);}
  };

  template<typename R>
  class AdSbMc<double,Mul,R>{
  private:  
    const double& lhs;
    const R& rhs;
  public:
    AdSbMc(const double& Lhs,const R& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Mul::calc(lhs,rhs.eval());}
  };

  template<>
  class AdSbMc<double,Mul,Field>{
  private:  
    const double& lhs;
    const Field& rhs;
  public:
    AdSbMc(const double& Lhs,const Field& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Mul::calc(lhs,rhs);}
  };

  template<typename Op>
  class AdSbMc<Field,Op,Field>{
  private:  
    const Field& lhs;
    const Field& rhs;
  public:
    AdSbMc(const Field& Lhs,const Field& Rhs):lhs(Lhs),rhs(Rhs){}
    Field eval() const{return Op::calc(lhs,rhs);}
  };

  template<typename L,typename R>
  AdSbMc<L,Add,R> operator+(const L& lhs, const R& rhs){
    return AdSbMc<L,Add,R>(lhs, rhs);
  }

  template<typename L,typename R>
  AdSbMc<L,Sub,R> operator-(const L& lhs, const R& rhs){
    return AdSbMc<L,Sub,R>(lhs, rhs);
  }

  template<typename R>
  AdSbMc<double,Mul,R> operator*(const double& lhs, const R& rhs){
    return AdSbMc<double,Mul,R>(lhs, rhs);
  }
}




#endif
