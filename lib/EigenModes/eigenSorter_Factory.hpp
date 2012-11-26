/*! @file eigenSorter_Factory.hpp
 *  @brief factory of EigenSorter
 */
#ifndef EIGENSORTER_FACTORY_INCLUDED
#define EIGENSORTER_FACTORY_INCLUDED

#include "eigenSorter.hpp"
#include "include/foprHermFactory.hpp"
#include "include/pugi_interface.h"
#include <memory>

class EigenSorterFactory{
public:
  virtual EigenSorter* getEigenSorter(const Fopr_Herm*)const =0;
  virtual EigenSorter* getEigenSorter()const =0;
  virtual ~EigenSorterFactory(){}
};

class EigenSorterFactory_low: public EigenSorterFactory{
  double thrs_;
public:
  EigenSorterFactory_low(double thrs):thrs_(thrs){}

  EigenSorter* getEigenSorter(const Fopr_Herm* opr)const{
    return new EigenSorter_low(opr->func(thrs_));
  }
  EigenSorter* getEigenSorter()const{
    return new EigenSorter_low(thrs_);
  }
};

class EigenSorterFactory_high: public EigenSorterFactory{
  double thrs_;
public:
  EigenSorterFactory_high(double thrs):thrs_(thrs){}

  EigenSorter* getEigenSorter(const Fopr_Herm* opr)const{
    return new EigenSorter_high(opr->func(thrs_));
  }
  EigenSorter* getEigenSorter()const{
    return new EigenSorter_high(thrs_);
  }
};



#endif
