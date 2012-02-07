/*!
 * @mainpage KEK common code for %Lattice QCD simulations
 *
 * \image html keklogo-c.jpg 
 * \image html JICFUSsymbolmark170px.png 
 *
 * JLQCD branch of the code for lattice simulations of QCD  
 *
 *
 * \section feat Features
 *
 * Current implementation:
 * - Actions (Gauge: Wilson, Rectangle, Fermion: 2 flavors, 2 flavors Ratio, Overlap)
 * - %Dirac operators (Wilson, Clover, Overlap, Even-odd preconditioned Wilson, Generalized Domain Wall (4d - 5d) )
 * - Linear Solvers (Conjugate Gradient Unpreconditioned, Conjugate Gradient Preconditioned, BiConjugate Gradient)
 * - Measurements (Quark Propagator [Wilson, Domain Wall], Gauge Quantities)
 * - Smearing (APE, Stout analytic), HMC Smeared runs
 * - Random Number Generators (Mersenne Twister)
 * - %XML control of program behavior
 *
 *
 * \section Plat Supported platforms
 * 
 * - <a href="http://gcc.gnu.org/">GNU compiler</a> (tested with g++ version 4.6.1)
 * - Multicore version tested with <a href="http://www.open-mpi.org/">openMPI</a> (version 1.4.4)  
 * - <a href="http://software.intel.com/en-us/articles/c-compilers/">INTEL compiler</a> (tested with icpc version 12.1.2)
 * - <a href="http://www-01.ibm.com/software/awdtools/xlcpp/">IBM XLC compiler</a> (tested with xlC version 11.1, cross platform compilation)
 *
 * \authors {<a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>,  Shoji Hashimoto, Jun-Ichi Noaki}
 *
 */
#include "documentation_pages.h"
//------------------------------------------------------------------------
/*!
 * @file main_hmc.cpp 
 * @brief Main source code for running HMC updates
 *
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 */
//------------------------------------------------------------------------



#include "include/common_code.hpp"
#include "include/commandline.hpp"
#include "HMC/hmcGeneral.hpp"

using namespace XML;

int run(GaugeField, node);

int main(int argc, char* argv[]){
  int status;
  CommandOptions Options = ReadCmdLine(argc, argv);
 
  CCIO::header(PACKAGE_STRING);
 
  //Reading input file
  node top_node = getInputXML(Options.filename);  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  //Initialize GaugeField using XML input
  GaugeField GaugeF(geom);
  GaugeF.initialize(top_node);

  node HMC_node = top_node;
  descend(HMC_node, "HMC");
  
  status = run(GaugeF, HMC_node);

  return status;
}


int run(GaugeField Gfield_, node HMC_node) {
  CCIO::header("Starting HMC updater");

  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_node);
  Integrators::Integr = 
    Integrators::createIntegratorFactory(HMC_node, Gfield_.Format);

  //Initialization of class
  HMCgeneral hmc_general(HMC_node);

  ////////////// HMC calculation /////////////////
  double elapsed_time;
  TIMING_START;
  try{
    CCIO::cout<< "-----  HMC starts\n"<<std::endl;
    hmc_general.evolve(Gfield_.U);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }

  TIMING_END(elapsed_time);
  ////////////// HMC calculation end /////////////////

  CCIO::cout << "Total elapsed time (s): "<< elapsed_time/1000.0 << "\n";

  CCIO::cout << "Saving configuration on disk in binary format\n";
  CCIO::SaveOnDisk< Format::Format_G >(Gfield_.U, "final_conf.bin");
  return 0;
}

