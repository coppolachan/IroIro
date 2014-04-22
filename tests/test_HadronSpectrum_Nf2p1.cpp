/*!
 * @file test_HadronSpectrum_Nf2p1.cpp
 * @brief Implementation of test_HadronSpectrum_Nf2p1 class
 *
 * Time-stamp: <2014-04-21 09:13:40 noaki>
 */
#include "include/factories.hpp"
#include "include/messages_macros.hpp"
#include "test_HadronSpectrum_Nf2p1.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"

using namespace std;

int Test_HadronSpectrum_Nf2p1::run(){
  XML::node node = input_.node;
  InputConfig config = input_.getConfig();

  //////////////////// Quark Propagator ////////////////////
  XML::descend(node,"QuarkPropUpDown");
  auto_ptr<QuarkPropagatorFactory> 
    qpf_ud(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_ud(qpf_ud->getQuarkProp(config));

  ////////////////// Quark Propagator(strange) ///////////////
  XML::next_sibling(node,"QuarkPropStrange");
  auto_ptr<QuarkPropagatorFactory> 
    qpf_s(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_s(qpf_s->getQuarkProp(config));
  
  ///////////////////// source creation ////////////////////////
  XML::next_sibling(node,"Source");
  auto_ptr<SourceFactory> 
    SrcFactory(Sources::createSourceFactory<SiteIndex,Format::Format_F>(node));

  _Message(DEBUG_VERB_LEVEL,"SrcFactory created\n");
  auto_ptr<Source> src(SrcFactory->getSource());
  _Message(DEBUG_VERB_LEVEL,"Source created\n");

  ////////////////////// Propagator actual calculation ///////////////////////
  prop_t sq_ud, sq_s;  //Defines a vector of fields
  CCIO::cout << "            _ \n";
  CCIO::cout << "  _   _  __| |      _ __  _ __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __ \n";
  CCIO::cout << " | | | |/ _` |_____| '_ \\| '__/ _ \\| '_ \\ / _` |/ _` |/ _` | __/ _ \\| '__|\n";
  CCIO::cout << " | |_| | (_| |_____| |_) | | | (_) | |_) | (_| | (_| | (_| | || (_) | |   \n";
  CCIO::cout << "  \\__,_|\\__,_|     | .__/|_|  \\___/| .__/ \\__,_|\\__, |\\__,_|\\__\\___/|_|   \n";
  CCIO::cout << "                   |_|             |_|          |___/                     \n\n";

  qprop_ud->calc(sq_ud,*src);

  CCIO::cout << "                                                    _              \n";
  CCIO::cout << "  ___       _ __  _ __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __  \n";
  CCIO::cout << " / __|_____| '_ \\| '__/ _ \\| '_ \\ / _` |/ _` |/ _` | __/ _ \\| '__| \n";
  CCIO::cout << " \\__ \\_____| |_) | | | (_) | |_) | (_| | (_| | (_| | || (_) | |    \n";
  CCIO::cout << " |___/     | .__/|_|  \\___/| .__/ \\__,_|\\__, |\\__,_|\\__\\___/|_|    \n";
  CCIO::cout << "           |_|             |_|          |___/                      \n\n";

  qprop_s->calc(sq_s,*src);

  CCIO::cout << "\n ::::::::::::: Output saved in "<< input_.output.c_str()<<"\n";
  ofstream writer(input_.output.c_str());

  /////////////////////////////// Meson Correlation Functions /////////////////////////////////////
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);

  CCIO::cout << " ::::::::::::: Performing contractions for meson correlators\n";
  if(Communicator::instance()->primaryNode()) writer<<"====== meson correlators ======\n";

  output_meson(writer,  pp.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Pion ------");
  output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_ud),"------ Kaon ------");
  output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_s), "------ Eta_s ------");
  
  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_X ------");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_Y ------");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_Z ------");

  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_X ------");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_Y ------");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_Z ------");

  output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_X ------");
  output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_Y ------");
  output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_Z ------");

  ////////// Baryon Correlation Functions ///////////
  BaryonCorrelator baryons(sq_ud,sq_s);

  CCIO::cout << " ::::::::::::: Performing contractions for baryon correlators\n";
  if(Communicator::instance()->primaryNode()) writer<<"====== baryon correlators ======\n";

  output_baryon(writer,baryons.nucleon(UP),     baryons.nucleon(DN),     "------ Nucleon ------");
  output_baryon(writer,baryons.sigma8(UP),      baryons.sigma8(DN),      "------ Sigma8 ------");
  output_baryon(writer,baryons.xi8(UP),         baryons.xi8(DN),         "------ Xi8 ------");
  output_baryon(writer,baryons.lambda(UP),      baryons.lambda(DN),      "------ Lambda ------");
  
  output_baryon(writer,baryons.delta(XDIR,UP),  baryons.delta(XDIR,DN),  "------ Delta_X ------");
  output_baryon(writer,baryons.delta(YDIR,UP),  baryons.delta(YDIR,DN),  "------ Delta_Y ------");
  output_baryon(writer,baryons.delta(ZDIR,UP),  baryons.delta(ZDIR,DN),  "------ Delta_Z ------");
  
  output_baryon(writer,baryons.omega(XDIR,UP),  baryons.omega(XDIR,DN),  "------ Omega_X ------");
  output_baryon(writer,baryons.omega(YDIR,UP),  baryons.omega(YDIR,DN),  "------ Omega_Y ------");
  output_baryon(writer,baryons.omega(ZDIR,UP),  baryons.omega(ZDIR,DN),  "------ Omega_Z ------");
  
  output_baryon(writer,baryons.sigma10(XDIR,UP),baryons.sigma10(XDIR,DN),"------ Sigma10_X ------");
  output_baryon(writer,baryons.sigma10(YDIR,UP),baryons.sigma10(YDIR,DN),"------ Sigma10_Y ------");
  output_baryon(writer,baryons.sigma10(ZDIR,UP),baryons.sigma10(ZDIR,DN),"------ Sigma10_Z ------");
  
  output_baryon(writer,baryons.xi10(XDIR,UP),   baryons.xi10(XDIR,DN),   "------ Xi10_X ------");
  output_baryon(writer,baryons.xi10(YDIR,UP),   baryons.xi10(YDIR,DN),   "------ Xi10_Y ------");
  output_baryon(writer,baryons.xi10(ZDIR,UP),   baryons.xi10(ZDIR,DN),   "------ Xi10_Z ------");

  /////////////////////////////// Extra Meson Channels /////////////////////////////////////
  MesonCorrelator ss(Scalar),a4a4(AVector4),a4p(AV4_PS),pa4(PS_AV4); 

  CCIO::cout << " ::::::::::::: Performing contractions for extra meson channels\n";
  if(Communicator::instance()->primaryNode()) writer<<"====== extra meson channels ======\n";

  output_meson(writer, ss.calculate<Format::Format_F>(sq_ud,sq_ud),"------ SS(ud,ud) ------");
  output_meson(writer, ss.calculate<Format::Format_F>(sq_s, sq_ud),"------ SS(s,ud) ------");
  output_meson(writer, ss.calculate<Format::Format_F>(sq_s, sq_s), "------ SS(s,s) ------");
  
  output_meson(writer,a4a4.calculate<Format::Format_F>(sq_ud,sq_ud),"------ A4A4(ud,ud) ------");
  output_meson(writer,a4a4.calculate<Format::Format_F>(sq_s, sq_ud),"------ A4A4(s,ud) ------");
  output_meson(writer,a4a4.calculate<Format::Format_F>(sq_s, sq_s), "------ A4A4(s,s) ------");
  /// <A4(t)P(0)>
  output_meson(writer,a4p.calculate<Format::Format_F>(sq_ud,sq_ud),"------ A4P(ud,ud) ------");
  output_meson(writer,a4p.calculate<Format::Format_F>(sq_s, sq_ud),"------ A4P(s,ud) ------");
  output_meson(writer,a4p.calculate<Format::Format_F>(sq_s, sq_s), "------ A4P(s,s) ------");
  /// <P(t)A4(0)>
  output_meson(writer,pa4.calculate<Format::Format_F>(sq_ud,sq_ud),"------ PA4(ud,ud) ------");
  output_meson(writer,pa4.calculate<Format::Format_F>(sq_s, sq_ud),"------ PA4(s,ud) ------");
  output_meson(writer,pa4.calculate<Format::Format_F>(sq_s, sq_s), "------ PA4(s,s) ------");
  
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


