/*! @file tests.hpp
 * @brief Abstract class for tests
 */
#ifndef TESTS_GEN_H_
#define TESTS_GEN_H_

#define TEST_PASSED 1
#define TEST_FAILED 0

#include "include/commandline.hpp"
#include "lib/Measurements/measGeneral.hpp"
#include "lib/Tools/jobUtils.hpp"

class TestGeneral {
public:
  virtual int run() = 0;
};

namespace TestEnv{

  template <class TestClass>
  TestClass StartUp(int argc, char* argv[]){
    CommandOptions Options = ReadCmdLine(argc, argv);

    //Reading input file
    XML::node top_node = XML::getInputXML(Options.filename);  

    //Initializing geometry using XML input
    Geometry geom(top_node);

    //Initialize GaugeField using XML input
    GaugeGlobal GaugeF(geom);

    GaugeF.initialize(top_node);
    
    // Echo of input xml
    JobUtils::echo_input(Options.filename);

    return TestClass(top_node, GaugeF);
  }

  template <class TestClass>
  int StartRun(int argc, char* argv[]){
    CommandOptions Options = ReadCmdLine(argc, argv);

    //Reading input file   
    XML::node top_node = XML::getInputXML(Options.filename);  

    MeasGeneral meas(top_node,Geometry(top_node));

    // Echo of input xml
    JobUtils::echo_input(Options.filename);

    meas.do_meas<TestClass>();
    
    return 0;
  }
}

#define CREATE_RUN_TEST(Class){ int result = TestEnv::StartRun<Class>(argc,argv)
#define CLEAR_TEST }

#define CREATE_TEST(Class) { Class TestObject = TestEnv::StartUp<Class>(argc, argv)
#define RUN_TEST TestObject.run()
#define CLEAR_TEST }

#endif //TESTS_GEN_H_
