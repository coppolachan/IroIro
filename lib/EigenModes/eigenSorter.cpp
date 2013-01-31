/*!@file eigenSorter.cpp
 * @brief implementation of the eigenSorter class 
 */
#include "eigenSorter.hpp"

#ifndef FIELD_INCLUDED
#include "include/field.h"
#endif

#include <algorithm>
#include <cmath>

using namespace std;

//*** EigenSorter_high ***
bool EigenSorter_high::greater_pair(pair<double,Field>& left,
				  pair<double,Field>& right){
  return fabs(left.first) > fabs(right.first);}  

bool EigenSorter_high::greater_lmd(double left,double right){
  return fabs(left) > fabs(right);}  

bool EigenSorter_high::beyond_thrs(double lmd) const{
  return fabs(lmd) < fabs(thrs_);}

void EigenSorter_high::push(vector<double>& lmd,int N) const{
  std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),greater_lmd);}

void EigenSorter_high::push(vector<double>& lmd,vector<Field>& evec,int N)const{
  vector<pair<double, Field> > emod;
  vector<pair<double, Field> >::iterator it;

  for(int i=0;i<lmd.size();++i)
    emod.push_back(pair<double,Field>(lmd[i],evec[i]));
		   
  partial_sort(emod.begin(),emod.begin()+N,emod.end(),greater_pair);
  it=emod.begin();
  for(int i=0;i<N;++i){
    lmd[i]=it->first;
    evec[i]=it->second;
    ++it;
  }
}

//*** EigenSorter_low ***
bool EigenSorter_low::less_pair(pair<double,Field>& left,
			      pair<double,Field>& right){
  return fabs(left.first) < fabs(right.first);}  

bool EigenSorter_low::less_lmd(double left,double right){
  return fabs(left) < fabs(right);}  

bool EigenSorter_low::beyond_thrs(double lmd) const{
  return fabs(lmd) > fabs(thrs_);}

void EigenSorter_low::push(vector<double>& lmd,int N) const{
  std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),less_lmd);}

void EigenSorter_low::push(vector<double>& lmd,vector<Field>& evec,int N) const{
  vector<pair<double, Field> > emod;
  vector<pair<double, Field> >::iterator it;

  for(int i=0;i<lmd.size();++i)
    emod.push_back(pair<double,Field>(lmd[i],evec[i]));

  partial_sort(emod.begin(),emod.begin()+N,emod.end(),less_pair);
  it=emod.begin();
  for(int i=0;i<N;++i){
    lmd[i]=it->first;
    evec[i]=it->second;
    ++it;
  }
}


