/*!
 * JLQCD code for measurements for the Hadron correlators for 2+1 flavors
 * 
//------------------------------------------------------------------------
/*!
 * @file main_hadron_correlator_Nf2p1.cpp 
 * @brief Main source code for Hadron correlators (2+1 flavors)
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 * @author Jun-Ichi Noaki
 */
//------------------------------------------------------------------------
#include "include/iroiro_code.hpp"
#include "include/commandline.hpp"
#include "include/geometry.hpp"
#include "Main/gaugeGlobal.hpp"
#include "lib/Tools/jobUtils.hpp"
#include "HadronSpectrum_Nf2p1.hpp"

using namespace XML;

int run_hadron_corr(MeasGeneral&);

int main(int argc, char* argv[]){
  int status;
  CommandOptions Options = ReadCmdLine(argc, argv);
  
  //Reading input file   
  XML::node top_node = XML::getInputXML(Options.filename);  
  
  //Initializing Measurements using XML input
  MeasGeneral meas(top_node,Geometry(top_node));
  
  // Echo of input xml
  JobUtils::echo_input(Options.filename);
  
  status = run_hadron_corr(meas);//meas.do_meas<TestClass>();
  
  return status;
}


int run_hadron_corr(MeasGeneral& measurement) {
  
  CCIO::header(IROIRO_PACKAGE_STRING);

  measurement.do_meas<HadronSpectrum_Nf2p1>();

}
