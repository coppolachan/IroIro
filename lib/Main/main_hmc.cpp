/*!
 * @mainpage KEK code for %Lattice QCD simulations
 *
 * \image html LogoIroIro170px.png 
 * \image html keklogo-c.jpg 
 * \image html JICFUSsymbolmark170px.png 
 * 
 *
 * JLQCD code for lattice simulations of QCD  
 * 
 *
 * \section feat Features
 *
 * Current implementation:
 * - Actions (Gauge: Wilson, Rectangle, Fermion: 2 flavors, 2 flavors Ratio, RHMC Nf flavors, RHMC Nf flavors ratio, Overlap)
 * - %Dirac operators (Wilson, Clover, Staggered, Adjoint Staggered, Overlap, Even-odd preconditioned Wilson, Generalized Domain Wall (4d-5d), Wilson Brillouin, Hybrid DomainWall-Overlap, Moebius Kernel)
 * - Linear %Solvers (Conjugate Gradient, BiCG Stabilized, Rational-Multishift)
 * - Measurements (Quark propagators [Wilson, Domain Wall], momentum space propagators, Meson and Baryon correlators, Eigenmodes, Low mode preconditioning, Gauge fixing - Coulomb & Landau,Gauge quantities [Plaquette, Polyakov Loop, Wilson Loop], Topological charge, Wilson Flow)
 * - Smearing (APE, Stout analytic), HMC Smeared runs
 * - Random Number Generators (Mersenne Twister, Dynamic Creation Mersenne Twister)
 * - I/O support (Plain ASCII, Plain binary, ILDG, NERSC, MILC, JLQCD-legacy)
 * - Peter Boyle's <a href="http://www2.ph.ed.ac.uk/~paboyle/bagel/">BAGEL</a>/BFM integration
 * - %XML control of program behavior
 * - Lightweight, compiles in a minute...
 *
 *
 * \section Plat Supported platforms
 * 
 * - <a href="http://gcc.gnu.org/">GNU compiler</a> (tested with g++ version 4.6.1)
 * - Multicore version tested with <a href="http://www.open-mpi.org/">openMPI</a> (version 1.4.4)  
 * - <a href="http://software.intel.com/en-us/articles/c-compilers/">INTEL compiler</a> (tested with icpc version 12.1.2)
 * - <a href="http://www-01.ibm.com/software/awdtools/xlcpp/">IBM XLC compiler</a> (tested with xlC version 12.1, cross platform compilation)
 *
 * \section Dev Development
 * 
 * - Main developers: <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>, J. Noaki
 * - Specific routine development: Y. Cho, H. Fukaya 
 * - Code mantainers: G. Cossu, J. Noaki 
 * - Performance optimization:  P. Boyle (BFM), G. Cossu, J. Doi (IBM libs)
 * - Testers:  S. Hashimoto,  T. Kaneko, J. Noaki
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
#include "include/iroiro_code.hpp"
#include "include/commandline.hpp"
#include "include/geometry.hpp"
#include "Main/gaugeGlobal.hpp"
#include "HMC/hmcGeneral.hpp"
#include "lib/Tools/jobUtils.hpp"

using namespace XML;

int run_hmc(GaugeField, TrajInfo*, node);

int main(int argc, char* argv[]){
  int status;
  CommandOptions Options = ReadCmdLine(argc, argv);
 
   //Reading input file
  node top_node = getInputXML(Options.filename);  

  //Initializing geometry using XML input
  Geometry geom(top_node);

  GaugeGlobal GaugeF(geom);
  TrajInfo Info = GaugeF.initialize(top_node);

  node HMC_node = top_node;
  descend(HMC_node, "HMC");

  // Echo of input xml
  JobUtils::echo_input(Options.filename);
  
  status = run_hmc(GaugeF, &Info, HMC_node);

  return status;
}


int run_hmc(GaugeField GaugeF_, TrajInfo* Info, node HMC_node) {

  CCIO::header(IROIRO_PACKAGE_STRING);
 
  
  CCIO::header("Starting HMC updater");

  RNG_Env::RNG = RNG_Env::createRNGfactory(HMC_node);
 
  CCIO::cout << "Creating integrator factory\n";

  Integrators::Integr = 
    Integrators::createIntegratorFactory(HMC_node);

  //Initialization of HMC classes from XML file
  CCIO::cout << "Initialization of HMC classes from XML file\n";
  HMCgeneral hmc_general(HMC_node, Info);

  ////////////// HMC calculation starts /////////////////
  double elapsed_time;
  TIMING_START;
  try{
    CCIO::cout<< "-------------------  HMC starts\n"<<std::endl;
    hmc_general.evolve(GaugeF_);
  }catch(const char* error){
    CCIO::cerr << error << std::endl;
    return EXIT_FAILURE;
  }
  TIMING_END(elapsed_time);
  ////////////// HMC calculation ends /////////////////

  CCIO::cout << "Total elapsed time (s): "<< elapsed_time/1000.0 << "\n";

  CCIO::cout << "Saving configuration on disk in binary format\n";
  CCIO::SaveOnDisk< Format::Format_G >(GaugeF_.data, "final_conf.bin");
  
  return 0;
}

