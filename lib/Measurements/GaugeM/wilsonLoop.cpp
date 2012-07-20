/*! @file wilsonLoop.cpp
 *  @brief implementation of wilsonLoop.cpp
 */
#include "wilsonLoop.hpp"

using namespace std;

vector<double> WilsonLoop::calc(const GaugeField& u)const{
  GaugeField1D ut = DirSlice(u,mu_dir_);
  for(int nu=0; nu<NDIM_;++nu){
    if(nu != mu_dir_){
      GaugeField1D v = DirSlice(u,]-)
      uarm = shift(,mu,Forward());
  
}
