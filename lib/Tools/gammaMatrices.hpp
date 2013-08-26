/*!
 * @file gammaMatrices.hpp
 * @brief Declares gamma matrices classes 
 */
#ifndef GAMMA_MATR_HPP_
#define GAMMA_MATR_HPP_
#include "include/macros.hpp"
#include <vector>

namespace GammaMatrices {
  //Chiral representation

  //maybe decided at compile time

  struct GammaResult{
    int spn;
    double facr,faci;
  };
  
  class Gamma {
  protected:
    std::vector<int> idx_;
    std::vector<double> fac_;
  public:
    Gamma():fac_(2*NDIM_,0.0){}
    GammaResult operator() (int sp) const;
  };

  struct Unit: public Gamma{Unit();};
  struct Gamma1: public Gamma{Gamma1();};
  struct Gamma2: public Gamma{Gamma2();};
  struct Gamma3: public Gamma{Gamma3();};
  struct Gamma4: public Gamma{Gamma4();};
  struct Gamma5: public Gamma{Gamma5();};

  struct Gamma1_5: public Gamma {Gamma1_5();};
  struct Gamma2_5: public Gamma {Gamma2_5();};
  struct Gamma3_5: public Gamma {Gamma3_5();};
  struct Gamma4_5: public Gamma {Gamma4_5();};

  struct Sigma12_5: public Gamma {Sigma12_5();};
  struct Sigma13_5: public Gamma {Sigma13_5();};
  struct Sigma14_5: public Gamma {Sigma14_5();};
  struct Sigma23_5: public Gamma {Sigma23_5();};
  struct Sigma24_5: public Gamma {Sigma24_5();};
  struct Sigma34_5: public Gamma {Sigma34_5();};

  struct CConj:   public Gamma {CConj();};     //C= Gamma2*Gamma4
  struct CGamma1: public Gamma {CGamma1();};   //C*Gamma1
  struct CGamma2: public Gamma {CGamma2();};   //C*Gamma2
  struct CGamma3: public Gamma {CGamma3();};   //C*Gamma3

  struct CGamma5: public Gamma {CGamma5();};   //C*Gamma5
}

#endif // GAMMA_MATR_HPP_
