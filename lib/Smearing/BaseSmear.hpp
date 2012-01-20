/*
  @file BaseSmear.hpp

  @brief Declares base smearing class Smear

 */

#ifndef BASE_SMEAR_H
#define BASE_SMEAR_H

class Field;

class Smear {
public:
  virtual ~Smear(){}
  
  virtual void smear(Field&, const Field&) = 0;
};


#endif
