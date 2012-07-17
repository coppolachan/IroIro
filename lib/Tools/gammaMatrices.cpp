#include "gammaMatrices.hpp"

namespace GammaMatrices{
  
  GammaResult Gamma::operator() (int spinor_in) const{
    GammaResult res;
    res.spinor_index = indexes[spinor_in];
    res.complex_factor[0] = factors[2*spinor_in];
    res.complex_factor[1] = factors[2*spinor_in+1];
    return res;
  }
  
  Unit::Unit(){
    indexes[0] = 0;
    indexes[1] = 1;
    indexes[2] = 2;
    indexes[3] = 3;
      
    factors[0] = 1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] = 1.0;  factors[5] = 0.0;
    factors[6] = 1.0;  factors[7] = 0.0;
  }

  Gamma1::Gamma1(){
    indexes[0] = 3;
    indexes[1] = 2;
    indexes[2] = 1;
    indexes[3] = 0;
      
    factors[0] = 0.0;  factors[1] =-1.0;
    factors[2] = 0.0;  factors[3] =-1.0;
    factors[4] = 0.0;  factors[5] = 1.0;
    factors[6] = 0.0;  factors[7] = 1.0;
  }

  Gamma2::Gamma2(){
    indexes[0] = 3;
    indexes[1] = 2;
    indexes[2] = 1;
    indexes[3] = 0;
    
    factors[0] =-1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] = 1.0;  factors[5] = 0.0;
    factors[6] =-1.0;  factors[7] = 0.0;
  }

  Gamma3::Gamma3(){
    indexes[0] = 2;
    indexes[1] = 3;
    indexes[2] = 0;
    indexes[3] = 1;
      
    factors[0] = 0.0;  factors[1] =-1.0;
    factors[2] = 0.0;  factors[3] = 1.0;
    factors[4] = 0.0;  factors[5] = 1.0;
    factors[6] = 0.0;  factors[7] =-1.0;
  }

  Gamma4::Gamma4(){
    indexes[0] = 0;
    indexes[1] = 1;
    indexes[2] = 2;
    indexes[3] = 3;
    
    factors[0] = 1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] =-1.0;  factors[5] = 0.0;
    factors[6] =-1.0;  factors[7] = 0.0;
  }

  Gamma5::Gamma5(){
    indexes[0] = 2;
    indexes[1] = 3;
    indexes[2] = 0;
    indexes[3] = 1;
      
    factors[0] = 1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] = 1.0;  factors[5] = 0.0;
    factors[6] = 1.0;  factors[7] = 0.0;
  }

  Gamma1_5::Gamma1_5(){
    indexes[0] = 1;
    indexes[1] = 0;
    indexes[2] = 3;
    indexes[3] = 2;
      
    factors[0] = 0.0;  factors[1] =-1.0;
    factors[2] = 0.0;  factors[3] =-1.0;
    factors[4] = 0.0;  factors[5] = 1.0;
    factors[6] = 0.0;  factors[7] = 1.0;
  }

  Gamma2_5::Gamma2_5(){
    indexes[0] = 1;
    indexes[1] = 0;
    indexes[2] = 3;
    indexes[3] = 2;
      
    factors[0] =-1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] = 1.0;  factors[5] = 0.0;
    factors[6] =-1.0;  factors[7] = 0.0;
  }

  Gamma3_5::Gamma3_5(){
    indexes[0] = 0;
    indexes[1] = 1;
    indexes[2] = 2;
    indexes[3] = 3;
      
    factors[0] = 0.0;  factors[1] =-1.0;
    factors[2] = 0.0;  factors[3] = 1.0;
    factors[4] = 0.0;  factors[5] = 1.0;
    factors[6] = 0.0;  factors[7] =-1.0;
  }

  Gamma4_5::Gamma4_5(){
    indexes[0] = 2;
    indexes[1] = 3;
    indexes[2] = 0;
    indexes[3] = 1;
      
    factors[0] = 1.0;  factors[1] = 0.0;
    factors[2] = 1.0;  factors[3] = 0.0;
    factors[4] =-1.0;  factors[5] = 0.0;
    factors[6] =-1.0;  factors[7] = 0.0;
  }

  ChargeConj::ChargeConj(){ //Gamma_2*Gamma_4
    indexes[0] = 3;
    indexes[1] = 2;
    indexes[2] = 1;
    indexes[3] = 0;
      
    factors[0] = 1.0;  factors[1] = 0.0;
    factors[2] =-1.0;  factors[3] = 0.0;
    factors[4] = 1.0;  factors[5] = 0.0;
    factors[6] =-1.0;  factors[7] = 0.0;
  }
};
