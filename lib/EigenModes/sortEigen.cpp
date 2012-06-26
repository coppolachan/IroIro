//--------------------------------------------------------------------
// sortEigen.cpp
//--------------------------------------------------------------------
#include "sortEigen.h"

#ifndef FIELD_INCLUDED
#include "include/field.h"
#endif

#include <algorithm>
#include <cmath>

using namespace std;

//*** SortEighen_high ***
bool SortEigen_high::greater_pair(pair<double,Field>& left,
				  pair<double,Field>& right){
  return fabs(left.first) > fabs(right.first);}  

bool SortEigen_high::greater_lmd(double left,double right){
  return fabs(left) > fabs(right);}  

bool SortEigen_high::beyond_thrs(double lmd,double thrs) const{
  return fabs(lmd) < fabs(thrs);}

void SortEigen_high::push(vector<double>& lmd,int N) const{
  std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),greater_lmd);}

void SortEigen_high::push(vector<double>& lmd,vector<Field>& evec,int N)const{
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

//*** SortEighen_low ***
bool SortEigen_low::less_pair(pair<double,Field>& left,
			      pair<double,Field>& right){
  return fabs(left.first) < fabs(right.first);}  

bool SortEigen_low::less_lmd(double left,double right){
  return fabs(left) < fabs(right);}  

bool SortEigen_low::beyond_thrs(double lmd,double thrs) const{
  return fabs(lmd) > fabs(thrs);}

void SortEigen_low::push(vector<double>& lmd,int N) const{
  std::partial_sort(lmd.begin(),lmd.begin()+N,lmd.end(),less_lmd);}

void SortEigen_low::push(vector<double>& lmd,vector<Field>& evec,int N) const{
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


