/*!
 * @file test_MesonSpectrum.cpp
 * @brief Definition of classes for calculating meson correlators
 */
#include "test_MesonSpectrum.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/GaugeM/staples.hpp"

using namespace std;

int Test_MesonSpectrum::run(){
  //// Quark Propagator ////
  XML::descend(node_,"QuarkProp");
  QuarkPropagatorFactory* 
    qpfact = QuarkPropagators::createQuarkPropagatorFactory(node_);
  QuarkPropagator* qprop = qpfact->getQuarkProp(conf_);
  
  //// source creation ////
  XML::next_sibling(node_,"Source");
  SourceFactory* SrcFactory 
    = Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_);
  Source* src = SrcFactory->getSource();

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  // meson correlators
  CCIO::cout << " ---- Making up meson correlators\n";
  MesonCorrelator pp(Pion), v1v1(Vector1);
  vector<double> Cpp   = pp.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq);  

  // output
  CCIO::cout << " ---- Output in "<< output_.c_str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(output_.c_str());

    writer<< setiosflags(  ios_base::scientific);
    writer<<"---pp meson correlator---"<<endl;
    for(int t=0; t<Cpp.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(20)<<setiosflags(ios_base::left )<< Cpp[t]<<endl;
    }
    writer<<"---v1v1 meson correlator---"<<endl;
    for(int t=0; t<Cv1v1.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(20)<<setiosflags(ios_base::left )<< Cv1v1[t]<<endl;
    }
    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  return 0;
}
