/*!
 * @file test_HadronSpectrum_Nf2p1.cpp
 * @brief implementation of test_HadronSpectrum_Nf2p1 class
 */
#include "include/factories.hpp"
#include "test_HadronSpectrum_Nf2p1.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
//#include "Measurements/GaugeM/staples.hpp"

using namespace std;

int Test_HadronSpectrum_Nf2p1::run(){
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

  CCIO::cout << " ---- Output in "<< output_.c_str()<<"\n";
  
  ofstream writer(output_.c_str());

  ////////// Meson Correlation Functions ///////////
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);

  CCIO::cout << " ---- Making meson correlators\n";
  if(Communicator::instance()->primaryNode()) writer<<"=== meson correlators ===\n";

  output_meson(writer,  pp.calculate<Format::Format_F>(sq_ud,sq_ud),"---pion---");
  output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_ud),"---kaon---");
  output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_s), "---eta_s---");
  
  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_ud,sq_ud),"---rho_X---");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_ud,sq_ud),"---rho_Y---");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_ud,sq_ud),"---rho_Z---");

  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_ud),"---k_star_X---");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_ud),"---k_star_Y---");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_ud),"---k_star_Z---");

  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_s),"---phi_X---");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_s),"---phi_Y---");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_s),"---phi_Z---");

  ////////// Baryon Correlation Functions ///////////
  BaryonCorrelator baryons(sq_ud,sq_s);

  CCIO::cout << " ---- Making baryon correlators\n";
  if(Communicator::instance()->primaryNode()) writer<<"=== baryon correlators ===\n";

  output_baryon(writer,baryons.nucleon(UP),     baryons.nucleon(DN),     "---nucleon---");
  output_baryon(writer,baryons.sigma8(UP),      baryons.sigma8(DN),      "---sigma8---");
  output_baryon(writer,baryons.xi8(UP),         baryons.xi8(DN),         "---xi8---");
  output_baryon(writer,baryons.lambda(UP),      baryons.lambda(DN),      "---lambda---");
  
  output_baryon(writer,baryons.delta(XDIR,UP),  baryons.delta(XDIR,DN),  "---delta_X---");
  output_baryon(writer,baryons.delta(YDIR,UP),  baryons.delta(YDIR,DN),  "---delta_Y---");
  output_baryon(writer,baryons.delta(ZDIR,UP),  baryons.delta(ZDIR,DN),  "---delta_Z---");
  
  output_baryon(writer,baryons.omega(XDIR,UP),  baryons.omega(XDIR,DN),  "---omega_X---");
  output_baryon(writer,baryons.omega(YDIR,UP),  baryons.omega(YDIR,DN),  "---omega_Y---");
  output_baryon(writer,baryons.omega(ZDIR,UP),  baryons.omega(ZDIR,DN),  "---omega_Z---");
  
  output_baryon(writer,baryons.sigma10(XDIR,UP),baryons.sigma10(XDIR,DN),"---sigma10_X---");
  output_baryon(writer,baryons.sigma10(YDIR,UP),baryons.sigma10(YDIR,DN),"---sigma10_Y---");
  output_baryon(writer,baryons.sigma10(ZDIR,UP),baryons.sigma10(ZDIR,DN),"---sigma10_Z---");
  
  output_baryon(writer,baryons.xi10(XDIR,UP),   baryons.xi10(XDIR,DN),   "---xi10_X---");
  output_baryon(writer,baryons.xi10(YDIR,UP),   baryons.xi10(YDIR,DN),   "---xi10_Y---");
  output_baryon(writer,baryons.xi10(ZDIR,UP),   baryons.xi10(ZDIR,DN),   "---xi10_Z---");

  return 0;
}

void Test_HadronSpectrum_Nf2p1::
output_meson(ofstream& writer,const vector<double>& Correl,string msg){
  int Lt = CommonPrms::instance()->Lt();
  
  if(Communicator::instance()->primaryNode()){
    writer<< setiosflags(ios_base::scientific);
    writer<<msg.c_str()<<endl;   // print of the separator

    for(int t=0; t<Lt; ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Correl[t]
	    << endl;
    }
    writer<< resetiosflags(ios_base::scientific);
  }
}

void Test_HadronSpectrum_Nf2p1::
output_baryon(ofstream& writer,const correl_t& Upper,const correl_t& Lower,
	      string msg){
  int Lt = CommonPrms::instance()->Lt();
  
  if(Communicator::instance()->primaryNode()){
    writer<< setiosflags(ios_base::scientific);
    writer<<msg.c_str()<<endl;  // print of the separator

    for(int t=0; t<Lt; ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Upper[t].real()
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Lower[t].real()
	    << endl;
    }
    writer<< resetiosflags(ios_base::scientific);
  }
}


