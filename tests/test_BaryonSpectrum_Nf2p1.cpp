/*!
 * @file test_BaryonSpectrum_Nf2p1.cpp
 * @brief implementation of test_BaryonSpectrum_Nf2p1 class
 */
#include "include/factories.hpp"
#include "test_BaryonSpectrum_Nf2p1.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
//#include "Measurements/GaugeM/staples.hpp"

using namespace std;

int Test_BaryonSpectrum_Nf2p1::run(){
  //// Quark Propagator ////
  XML::descend(node_,"QuarkPropUpDown");
  auto_ptr<QuarkPropagatorFactory> 
    qpf_ud(QuarkPropagators::createQuarkPropagatorFactory(node_));
  auto_ptr<QuarkPropagator> qprop_ud(qpf_ud->getQuarkProp(conf_));

  //// Quark Propagator(strange) ////
  XML::next_sibling(node_,"QuarkPropStrange");
  auto_ptr<QuarkPropagatorFactory> 
    qpf_s(QuarkPropagators::createQuarkPropagatorFactory(node_));
  auto_ptr<QuarkPropagator> qprop_s(qpf_s->getQuarkProp(conf_));
  
  //// source creation ////
  XML::next_sibling(node_,"Source");
  auto_ptr<SourceFactory> 
    SrcFactory(Sources::createSourceFactory<SiteIndex,Format::Format_F>(node_));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq_ud, sq_s;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating ud-propagator\n";
  qprop_ud->calc(sq_ud,*src);

  CCIO::cout << " ---- Calculating s-propagator\n";
  qprop_s->calc(sq_s,*src);

  BaryonCorrelator baryon_correl(sq_ud,sq_s);

  CCIO::cout << " ---- Output in "<< output_.c_str()<<"\n";

  //// begin output //////////////////////
  ofstream writer(output_.c_str());

  output(writer,baryon_correl.nucleon(UP),baryon_correl.nucleon(DN),
	 "---nucleon---");
  output(writer,baryon_correl.sigma8(UP), baryon_correl.sigma8(DN),
	 "---sigma8---");
  output(writer,baryon_correl.xi8(UP),    baryon_correl.xi8(DN),
	 "---xi8---");
  output(writer,baryon_correl.lambda(UP), baryon_correl.lambda(DN),
	 "---lambda---");

  output(writer,baryon_correl.delta(XDIR,UP),baryon_correl.delta(XDIR,DN),
	 "---delta_X---");
  output(writer,baryon_correl.delta(YDIR,UP),baryon_correl.delta(YDIR,DN),
	 "---delta_Y---");
  output(writer,baryon_correl.delta(ZDIR,UP),baryon_correl.delta(ZDIR,DN),
	 "---delta_Z---");

  output(writer,baryon_correl.omega(XDIR,UP),baryon_correl.omega(XDIR,DN),
	 "---omega_X---");
  output(writer,baryon_correl.omega(YDIR,UP),baryon_correl.omega(YDIR,DN),
	 "---omega_Y---");
  output(writer,baryon_correl.omega(ZDIR,UP),baryon_correl.omega(ZDIR,DN),
	 "---omega_Z---");

  output(writer,baryon_correl.sigma10(XDIR,UP),baryon_correl.sigma10(XDIR,DN),
	 "---sigma10_X---");
  output(writer,baryon_correl.sigma10(YDIR,UP),baryon_correl.sigma10(YDIR,DN),
	 "---sigma10_Y---");
  output(writer,baryon_correl.sigma10(ZDIR,UP),baryon_correl.sigma10(ZDIR,DN),
	 "---sigma10_Z---");

  output(writer,baryon_correl.xi10(XDIR,UP),baryon_correl.xi10(XDIR,DN),
	 "---xi10_X---");
  output(writer,baryon_correl.xi10(YDIR,UP),baryon_correl.xi10(YDIR,DN),
	 "---xi10_Y---");
  output(writer,baryon_correl.xi10(ZDIR,UP),baryon_correl.xi10(ZDIR,DN),
	 "---xi10_Z---");
  return 0;
}

void Test_BaryonSpectrum_Nf2p1::
output(ofstream& writer,const correl_t& Upper,const correl_t& Lower,
       string msg){
  int Lt = CommonPrms::instance()->Lt();
  
  if(Communicator::instance()->primaryNode()){
    writer<< setiosflags(ios_base::scientific);
    writer<<msg.c_str()<<endl;
    for(int t=0; t<Lt; ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Upper[t].real()
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Lower[t].real()
	    <<endl;
    }
    writer<< resetiosflags(ios_base::scientific);
  }
}
