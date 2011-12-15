/*!
 * @file gammaMatrices.hpp
 *
 * @brief Declares gamma matrices classes 
 *
 */
#ifndef GAMMA_MATR_HPP_
#define GAMMA_MATR_HPP_

namespace GammaMatrices {
  //Chiral representation

  //maybe decided at compile time

  struct GammaResult{
    int spinor_index;
    double complex_factor[2];
  };
  
  class Gamma {
  protected:
    int indexes[4];
    double factors[8];
  public:
    GammaResult operator() (int spinor_in) const;
  };

  class Unit: public Gamma{
  public:
    Unit();
    ~Unit(){};
  };

  class Gamma1: public Gamma{
  public:
    Gamma1();
    ~Gamma1(){};
  };

  class Gamma2: public Gamma{
  public:
    Gamma2();
    ~Gamma2(){};
  };

  class Gamma3: public Gamma{
  public:
    Gamma3();
    ~Gamma3(){};
  };

  class Gamma4: public Gamma{
  public:
    Gamma4();
    ~Gamma4(){};
  };

  class Gamma5: public Gamma{
  public:
    Gamma5();
    ~Gamma5(){};
  };

  class Gamma1_5: public Gamma {
  public:
    Gamma1_5();
    ~Gamma1_5(){};

  };

  class Gamma2_5: public Gamma {
  public:
    Gamma2_5();
    ~Gamma2_5(){};

  };

  class Gamma3_5: public Gamma {
  public:
    Gamma3_5();
    ~Gamma3_5(){};

  };

  class Gamma4_5: public Gamma {
  public:
    Gamma4_5();
    ~Gamma4_5(){};

  };

  class ChargeConj: public Gamma {
    //Gamma4*Gamma2*Gamma5
  public:
    ChargeConj();
    ~ChargeConj(){};

  };

}

#endif // GAMMA_MATR_HPP_
