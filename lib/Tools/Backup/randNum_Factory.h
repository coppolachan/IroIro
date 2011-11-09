/*!
  @file randNum_Factory.h
  
  @brief Defines the Factories for Random Number Generators
 
*/

#ifndef RANDNUM_FACT_
#define RANDNUM_FACT_

#include "randNum.h"
#include "randNum_MT19937.h"
#include "include/singleton.h"
#include "include/pugi_interface.h"

/*!
 *@class RandomNumberCreator
 *
 *@brief Abstract Factory class for Random Number Generators
 *
 */
class RandomNumberCreator {
public:
  /*!
   Virtual function returning the Random Number Generator
   */
  virtual RandNum* getRandomNumGenerator() = 0;
};

class RandNum_MT19937_Creator : public RandomNumberCreator {
  unsigned long *init;
  int length;
  
public:
  RandNum* getRandomNumGenerator(){
    return createRNG();
  };
  RandNum_MT19937_Creator(unsigned long *in, int l):
    init(in),
    length(l){};

  
  RandNum_MT19937_Creator(XML::node node) {
    std::vector<unsigned long> inputs;
    inputs.resize(4);
    int it;
    
    XML::read_array(node, "init", inputs);
    init = new unsigned long[inputs.size()];

    for (it = 0; it < inputs.size(); it++) {
      init[it] = inputs[it];
    }
    length = inputs.size();
  }

  
private:  
  RandNum_MT19937* createRNG(){
    return new RandNum_MT19937(init, length); 
  };
  
};

namespace RNG_Env {
  static RandomNumberCreator* RNG;
  RandomNumberCreator* createRNGfactory(XML::node);  
}


#endif
