/*!
 * @file test_MesonSpectrum.cpp
 * @brief implementation of test_MesonSpectrum class
 */
#include "test_MesonSpectrum.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
//#include "Measurements/GaugeM/staples.hpp"
#include <memory>

using namespace std;

int Test_MesonSpectrum::run(){

  //// Quark Propagator ////
  XML::descend(node_,"QuarkProp");
  
  auto_ptr<QuarkPropagatorFactory> qpfact(QuarkPropagators::createQuarkPropagatorFactory(node_));
  auto_ptr<QuarkPropagator> qprop(qpfact->getQuarkProp(conf_));
  
  //// source creation ////
  XML::next_sibling(node_,"Source");
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  // meson correlators
  CCIO::cout << " ---- Making up meson correlators\n";
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  vector<double> Cpp   = pp.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq,sq);  

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
    writer<<"---v2v2 meson correlator---"<<endl;
    for(int t=0; t<Cv2v2.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(20)<<setiosflags(ios_base::left )<< Cv2v2[t]<<endl;
    }
    writer<<"---v3v3 meson correlator---"<<endl;
    for(int t=0; t<Cv3v3.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(20)<<setiosflags(ios_base::left )<< Cv3v3[t]<<endl;
    }
    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  
  return 0;
}
