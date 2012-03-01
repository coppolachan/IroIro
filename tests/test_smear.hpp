/*!
 * @file test_smear.hpp
 *
 * @brief Tests for the smearing classes APE, Stout
 *
 */
#ifndef TEST_SMEAR_H_
#define TEST_SMEAR_H_

#include "include/common_code.hpp"
#include "tests.hpp"

class Test_Smear : public TestGeneral{
private:
  const char* test_name;
  XML::node Smear_node_;
  GaugeField conf_;

  GaugeField smeared_u_;  /*!< Smeared field */
  GaugeField previous_u_; /*!< Temporary field */
  Field* gauge_pointer;   /*!< Pointer to a gauge field */

public:
  Test_Smear(const XML::node node, GaugeField& conf)
    :Smear_node_(node),
     conf_(conf){
    test_name = "TestSmear";
    XML::descend(Smear_node_, test_name, MANDATORY);        
  }

  ~Test_Smear(){}
  
  int run();
};

#endif
