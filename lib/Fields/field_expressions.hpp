#ifndef FIELD_EXPR_H_
#define FIELD_EXPR_H_

#include "include/field.h"

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
