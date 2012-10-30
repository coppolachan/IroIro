#ifndef BARYON_CONTRACT_INCLUDED
#define BARYON_CONTRACT_INCLUDED
#include "include/macros.hpp"
#include <stdlib.h>

namespace ColorUtils{

  class ColorEpsilon{
    std::vector<int> idx_;
  public:
    ColorEpsilon(){
      int idx[] = {0,1,2,
		   2,0,1,
		   1,2,0,
		   1,0,2,
		   2,1,0,
		   0,2,1};
      for(int i=0;i<6*3; ++i) idx_.push_back(idx[i]);
    }
    void get_elm(double& pm,std::vector<int>& ep,int id)const{
      pm =1.0-2.0*(id/NC_);
      for(int i=0; i<NC_; ++i) ep[i] = idx_[NC_*id +i];
    }      
  };
}

#endif
