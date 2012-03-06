/*!
  @file field.h
  @brief Definition of a general class storing a field
*/
#ifndef FIELD_INCLUDED
#define FIELD_INCLUDED
#include "Communicator/communicator.h"
#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>

class Field{
private:
  std::valarray<double> field;
public:

  Field(){}
  explicit Field(const std::valarray<double>& va):field(va){}
  Field(std::size_t size, double val=0.0):field(val,size){}
  Field(const Field& f):field(f.field){}

  std::size_t size() const {return field.size();}
  void resize(std::size_t size, double val=0.0){field.resize(size,val);}

  double operator[](const size_t i) const { return field[i]; }
  const double& operator[](const size_t i) { return field[i]; }

  std::valarray<double> operator[](const std::slice& sl) const { 
    return field[sl];}
  std::slice_array<double> operator[](const std::slice& sl) {
    return field[sl];}
  std::valarray<double> operator[](const std::gslice& sl) const{ 
    return field[sl];}
  std::gslice_array<double> operator[](const std::gslice& sl){
    return field[sl];}
  std::valarray<double> operator[](const std::valarray<size_t>& va) const{
    return field[va];}
  std::indirect_array<double> operator[](const std::valarray<size_t>& va){
    return field[va];}

  const std::valarray<double>& getva()const{ return field;}
  
  //direct address access
  #ifndef HITACHISR16K
  const double* getaddr(const size_t i)const{return &field[i];} 
  #else
  double* getaddr(const size_t i){return &field[i];}   
  #endif

  void set(const std::size_t i, double val){field[i]= val;}
  void set(const std::slice& sl, const std::valarray<double>& va){ 
    field[sl]= va;}
  void set(const std::gslice& sl, const std::valarray<double>& va){ 
    field[sl]= va;}
  void set(const std::valarray<std::size_t>& idx, 
	   const std::valarray<double>& va){ 
    field[idx]= va;
  }
  void add(const std::size_t i, double val){ field[i]+= val;}
  void add(const std::slice& sl, const std::valarray<double>& va){
    field[sl]+= va;}
  void add(const std::gslice& sl, const std::valarray<double>& va){ 
    field[sl]+= va;}
  void add(const std::valarray<std::size_t>& idx, 
	   const std::valarray<double>& va){ 
    field[idx]+= va;
  }

  double norm() const {
    double a = (field*field).sum();
    double b = Communicator::instance()->reduce_sum(a);
    return sqrt(b);
  }

  double average_abs() const{
    double a = 0.0;
    for(int i=0; i<field.size()/2; ++i) 
      a += sqrt(field[2*i]*field[2*i] +field[2*i+1]*field[2*i+1]);
    
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field.size()/2);
    return b/c;
  }

  double average_real() const{
    double a = 0.0;
    for(int i=0; i<field.size()/2; ++i) a += field[2*i];
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field.size()/2);
    return b/c;
  }

  double average_imag() const{
    double a = 0.0;
    for(int i=0; i<field.size()/2; ++i) a += field[2*i+1];
    double b = Communicator::instance()->reduce_sum(a);
    double c = Communicator::instance()->reduce_sum(field.size()/2);
    return b/c;
  }

  double max_element() const{
    std::valarray<double> fab = field*field;
    double max= fab[0]+fab[1];
    int idx = 0;
    for(int i=1; i< field.size()/2; ++i){

      double tmp = fab[2*i] +fab[2*i+1];
      if(max < tmp){
	max = tmp;
	idx = i;
      }
    }	
    int mn = Communicator::instance()->reduce_max(max,idx,field.size()/2);
    return sqrt(max);
  }

  double min_element() const{
    std::valarray<double> fab = field*field;
    double min= fab[0]+fab[1];
    int idx = 0;
    for(int i=1; i< field.size()/2; ++i){

      double tmp = fab[2*i] +fab[2*i+1];
      if(min > tmp){
	min = tmp;
	idx = i;
      }
    }	
    int mn = Communicator::instance()->reduce_min(min,idx,field.size()/2);
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
    out.write((char*)&field[0], sizeof(double)*field.size());
  }
  void read_stream(std::ifstream& in) {
    in.read((char*)&field[0], sizeof(double)*field.size());
  }
};

#include "Fields/field.inl"  //definition of inlined functions
#endif
