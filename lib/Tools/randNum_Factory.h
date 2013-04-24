/*! @file randNum_Factory.h
    @brief Defines the Factories for Random Number Generators
*/
#ifndef RANDNUM_FACT_
#define RANDNUM_FACT_

#include "include/singleton.h"
#include "include/pugi_interface.h"
#include "Communicator/communicator.h"

#include "randNum.h"
#include "RandomNumGen/randNum_MT19937.h"
#include "RandomNumGen/dcmt_wrapper.hpp"


/*!
 *@class RandomNumberCreator
 *@brief Abstract Factory class for Random Number Generators
 */
class RandomNumberCreator {
public:
  /*! Virtual function returning the Random Number Generator */
  virtual RandNum* getRandomNumGenerator() = 0;
};

//Specific factories 
////////////////////////////////////////////////////////////////////
class RandNum_MT19937_Creator : public RandomNumberCreator {
  unsigned long *init_;
  int length_;
  std::string filename_;
  std::vector<unsigned long> inputs_;
  bool fromfile_;

  RandNum_MT19937* createRNG(){
    if(fromfile_) return new RandNum_MT19937(filename_);
    else          return new RandNum_MT19937(init_,length_); 
  }
public:
  explicit RandNum_MT19937_Creator(std::string filename)
    :filename_(filename),fromfile_(true){}

  RandNum_MT19937_Creator(unsigned long *in, int l)
    :init_(in),length_(l),fromfile_(false){}

  RandNum_MT19937_Creator(XML::node node){
    if(node.child("seedFile")!= NULL){
      fromfile_= true;
      XML::read(node,"seedFile",filename_,MANDATORY); 
    }else{
      fromfile_= false;
      inputs_.resize(4);
      XML::read_array(node,"init",inputs_,MANDATORY);
      init_= &inputs_[0];
      length_= inputs_.size();
    }
  }
  RandNum* getRandomNumGenerator(){ return createRNG();}
};

class RandNum_DCMT_Creator : public RandomNumberCreator{
  unsigned long generatorSeed_;
  unsigned long seed_;
  std::string filename_;
  bool fromfile_;

  RandNum_DCMT* createRNG(){
    if(fromfile_) return new RandNum_DCMT(filename_);
    else          return new RandNum_DCMT(generatorSeed_,
					  seed_,
					  Communicator::instance()->id()); 
  }
 public:
  explicit RandNum_DCMT_Creator(std::string filename)
    :filename_(filename),fromfile_(true){}
  
  RandNum_DCMT_Creator(XML::node node){
    if(node.child("seedFile")!= NULL){
      fromfile_= true;
      XML::read(node,"seedFile",filename_,MANDATORY); 
    }else{
      fromfile_= false;
      XML::read(node,"genSeed",generatorSeed_,MANDATORY);
      XML::read(node,"seed",seed_,MANDATORY);
    }
  }
   RandNum_DCMT_Creator(unsigned long gS, int s)
     :generatorSeed_(gS), seed_(s), fromfile_(false){}

  RandNum* getRandomNumGenerator(){ return createRNG();}
};


//////////////////////////////////////////////////////////////
namespace RNG_Env {
  static RandomNumberCreator* RNG;
  RandomNumberCreator* createRNGfactory(XML::node);  
}

#endif
