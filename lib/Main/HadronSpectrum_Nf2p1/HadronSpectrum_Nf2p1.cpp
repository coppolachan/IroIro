/*!
 * @file HadronSpectrum_Nf2p1.cpp
 * @brief Implementation of HadronSpectrum_Nf2p1 class
 *
 * Time-stamp: <2014-10-05 11:36:16 noaki>
 */

#include "HadronSpectrum_Nf2p1.hpp"
#include "include/factories.hpp"
#include "include/messages_macros.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "Measurements/FermionicM/outputUtils.hpp"
using namespace std;

int HadronSpectrum_Nf2p1::run(){
  XML::node node = input_.node;
  InputConfig config = input_.config;

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
    SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));

  _Message(DEBUG_VERB_LEVEL,"SrcFactory created\n");
  auto_ptr<Source> src(SrcFactory->getSource());
  _Message(DEBUG_VERB_LEVEL,"Source created\n");

  ////////////////////// Propagator actual calculation ///////////////////////
  prop_t sq_ud, sq_s;  //Defines a vector of fields
  CCIO::cout<<"            _ \n";
  CCIO::cout<<"  _   _  __| |      _ __  _ __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __ \n";
  CCIO::cout<<" | | | |/ _` |_____| '_ \\| '__/ _ \\| '_ \\ / _` |/ _` |/ _` | __/ _ \\| '__|\n";
  CCIO::cout<<" | |_| | (_| |_____| |_) | | | (_) | |_) | (_| | (_| | (_| | || (_) | |   \n";
  CCIO::cout<<"  \\__,_|\\__,_|     | .__/|_|  \\___/| .__/ \\__,_|\\__, |\\__,_|\\__\\___/|_|\n";
  CCIO::cout<<"                   |_|             |_|          |___/                     \n\n";

  qprop_ud->calc(sq_ud,*src);

  CCIO::cout<<"                                                    _              \n";
  CCIO::cout<<"  ___       _ __  _ __ ___  _ __   __ _  __ _  __ _| |_ ___  _ __  \n";
  CCIO::cout<<" / __|_____| '_ \\| '__/ _ \\| '_ \\ / _` |/ _` |/ _` | __/ _ \\| '__|\n";
  CCIO::cout<<" \\__ \\_____| |_) | | | (_) | |_) | (_| | (_| | (_| | || (_) | |    \n";
  CCIO::cout<<" |___/     | .__/|_|  \\___/| .__/ \\__,_|\\__, |\\__,_|\\__\\___/|_|\n";
  CCIO::cout<<"           |_|             |_|          |___/                      \n\n";

  qprop_s->calc(sq_s,*src);

  CCIO::cout << "\n ::::::::::::: Output saved in "<< input_.output.c_str()<<"\n";
  ofstream writer(input_.output.c_str());

  Hadrons::mesonProp(sq_ud,sq_s,writer);
  Hadrons::baryonProp(sq_ud,sq_s,writer);
  Hadrons::mesonExtraProp(sq_ud,sq_s,writer);
}
