#include "randNum_MP.h" 

void MPrand::mp_get(std::valarray<double>& rn,const RandNum& rand){
  rn=0.0;
  
  if (rand.parallel_safe()){
    rand.get(rn);
    Communicator::instance()->sync();	
  } else {
    
    int NP = CommonPrms::instance()->NP();
    
    std::valarray<double> Rn_source(rn.size());
    for(int node=0; node<NP; ++node) {
      if(Communicator::instance()->primaryNode()) 
	rand.get(Rn_source);	
      
      Communicator::instance()->sync();
      Communicator::instance()->send_1to1(rn,Rn_source,rn.size(),node,0,node);
    }
  }
}



void MPrand::mp_get_gauss(std::valarray<double>& rn,const RandNum& rand){
  rn=0.0;
  
  if (rand.parallel_safe()){
    rand.get_gauss(rn);
    Communicator::instance()->sync();	
  } else {
    
    int NP = CommonPrms::instance()->NP();
    
    std::valarray<double> Rn_source(rn.size());
    for(int node=0; node<NP; ++node) {
	if(Communicator::instance()->primaryNode()) 
	  rand.get_gauss(Rn_source);	
	
	Communicator::instance()->sync();
	Communicator::instance()->send_1to1(rn,Rn_source,rn.size(),node,0,node);
    }
  }
}

