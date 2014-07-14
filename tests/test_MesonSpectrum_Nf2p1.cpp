/*!
 * @file test_MesonSpectrum_Nf2p1.cpp
 * @brief implementation of test_MesonSpectrum_Nf2p1 class
 */
#include "test_MesonSpectrum_Nf2p1.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
//#include "Measurements/GaugeM/staples.hpp"
#include <memory>

using namespace std;

int Test_MesonSpectrum_Nf2p1::run(){
  XML::node node = input_.node;
  InputConfig config = input_.getConfig();

  //// Quark Propagator ////
  XML::descend(node,"QuarkPropUpDown");
  auto_ptr<QuarkPropagatorFactory> qpf_ud(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_ud(qpf_ud->getQuarkProp(config));

  //// Quark Propagator(strange) ////
  XML::next_sibling(node,"QuarkPropStrange");
  auto_ptr<QuarkPropagatorFactory> qpf_s(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_s(qpf_s->getQuarkProp(config));
  
  //// source creation ////
  XML::next_sibling(node,"Source");
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq_ud, sq_s;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating ud-propagator\n";
  qprop_ud->calc(sq_ud,*src);

  CCIO::cout << " ---- Calculating s-propagator\n";
  qprop_s->calc(sq_s,*src);

  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  // S=0 correlators
  CCIO::cout << " ---- Making meson correlators\n";
  vector<double> Cpp   = pp.calculate<Format::Format_F>(  sq_ud,sq_ud);  
  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq_ud,sq_ud);  
  vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq_ud,sq_ud);  
  vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq_ud,sq_ud);  

  vector<double> Kpp   = pp.calculate<Format::Format_F>(  sq_s,sq_ud);  
  vector<double> Kv1v1 = v1v1.calculate<Format::Format_F>(sq_s,sq_ud);  
  vector<double> Kv2v2 = v2v2.calculate<Format::Format_F>(sq_s,sq_ud);  
  vector<double> Kv3v3 = v3v3.calculate<Format::Format_F>(sq_s,sq_ud);  

  vector<double> Spp   = pp.calculate<Format::Format_F>(  sq_s,sq_s);  
  vector<double> Sv1v1 = v1v1.calculate<Format::Format_F>(sq_s,sq_s);  
  vector<double> Sv2v2 = v2v2.calculate<Format::Format_F>(sq_s,sq_s);  
  vector<double> Sv3v3 = v3v3.calculate<Format::Format_F>(sq_s,sq_s);  

  // output
  CCIO::cout << " ---- Output in "<< input_.output.c_str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(input_.output.c_str());
    
    writer<< setiosflags(  ios_base::scientific);
    writer<<"==== ud-ud content ===="<<endl;
    writer<<"---pp correlator---"<<endl;
    for(int t=0; t<Cpp.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cpp[t]<<endl;
    }
    writer<<"---v1v1 correlator---"<<endl;
    for(int t=0; t<Cv1v1.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv1v1[t]<<endl;
    }
    writer<<"---v2v2 correlator---"<<endl;
    for(int t=0; t<Cv2v2.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv2v2[t]<<endl;
    }
    writer<<"---v3v3 correlator---"<<endl;
    for(int t=0; t<Cv3v3.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv3v3[t]<<endl;
    }
    writer<<"==== ud-s content ===="<<endl;
    writer<<"---pp correlator---"<<endl;
    for(int t=0; t<Kpp.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Kpp[t]<<endl;
    }
    writer<<"---v1v1 correlator---"<<endl;
    for(int t=0; t<Kv1v1.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Kv1v1[t]<<endl;
    }
    writer<<"---v2v2 correlator---"<<endl;
    for(int t=0; t<Kv2v2.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Kv2v2[t]<<endl;
    }
    writer<<"---v3v3 correlator---"<<endl;
    for(int t=0; t<Kv3v3.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Kv3v3[t]<<endl;
    }
    writer<<"==== s-s content ===="<<endl;
    writer<<"---pp correlator---"<<endl;
    for(int t=0; t<Spp.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Spp[t]<<endl;
    }
    writer<<"---v1v1 correlator---"<<endl;
    for(int t=0; t<Sv1v1.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Sv1v1[t]<<endl;
    }
    writer<<"---v2v2 correlator---"<<endl;
    for(int t=0; t<Sv2v2.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Sv2v2[t]<<endl;
    }
    writer<<"---v3v3 correlator---"<<endl;
    for(int t=0; t<Sv3v3.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Sv3v3[t]<<endl;
    }
    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();

  return 0;
}
