/*! @file eigenSorter_Factory.hpp
 *  @brief factory of SortEigen 
 */
#ifndef EIGENSORTER_FACTORY_INCLUDED
#define EIGENSORTER_FACTORY_INCLUDED

#include "eigenSorter.hpp"
#include "include/pugi_interface.h"

class EigenSorterFactory{
public:
  virtual EigenSorter* getEigenSorter() = 0;
  virtual ~EigenSorterFactory(){}
};

class EigenSorterFactory_high: public EigenSorterFactory {
public:
  EigenSorter* getEigenSorter(){ return new EigenSorter_high;}
};

class EigenSorterFactory_low: public EigenSorterFactory {
public:
  EigenSorter* getEigenSorter(){ return new EigenSorter_low;}
};

namespace EigenModes{
  EigenSorterFactory* createEigenSorterFactory(XML::node node);
}
#endif
