/*!@file eigenSorter.h pp
 * @brief declaration of the EigenSorter classes
 */
#ifndef EIGENSORTER_INCLUDED
#define EIGENSORTER_INCLUDED

#include "Communicator/comm_io.hpp"
#include <vector>
#include <utility>

class Field;

class EigenSorter{
public:
  virtual void push(std::vector<double>& lmd,int N) const = 0;
  virtual void push(std::vector<double>& lmd,
		    std::vector<Field>& evec,int N) const = 0;
  virtual bool beyond_thrs(double lmd) const =0;
  virtual double thrs()const = 0;
  virtual ~EigenSorter(){}
};

class EigenSorter_high :public EigenSorter{
private:
  double thrs_;
  static bool greater_lmd(double left,double right);
  static bool greater_pair(std::pair<double,Field>& left, 
			   std::pair<double,Field>& right);
public:
  EigenSorter_high(double thrs):thrs_(thrs){
    CCIO::cout<<"EigenSorter_high is constructed. Lower threshold = "<<thrs_<<"\n";
  }
  void push(std::vector<double>& lmd,int N) const;
  void push(std::vector<double>& lmd,std::vector<Field>& evec,int N) const;
  bool beyond_thrs(double lmd) const;
  double thrs()const{return thrs_;}
};

class EigenSorter_low :public EigenSorter{
private:
  double thrs_;
  static bool less_lmd(double left,double right);
  static bool less_pair(std::pair<double,Field>& left, 
			std::pair<double,Field>& right);
public:
  EigenSorter_low(double thrs):thrs_(thrs){
    CCIO::cout<<"EigenSorter_low is constructed. Higher threshold = "<<thrs_<<"\n";
  }
  void push(std::vector<double>& lmd,int N) const;
  void push(std::vector<double>& lmd,std::vector<Field>& evec,int N) const;
  bool beyond_thrs(double lmd) const;
  double thrs()const{return thrs_;}
};

#endif
