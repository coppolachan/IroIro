/*!
 * @file test_HadronSpectrum_MultiQprop.cpp
 * @brief Implementation of test_HadronSpectrum_MultiQprop class
 * This class obtains hadron correlators from combinations of 
 * Source-smeared propagators and multi-mass propagators
 * Time-stamp: <2014-07-10 16:36:19 noaki>
 */
#include "include/factories.hpp"
#include "include/messages_macros.hpp"
#include "test_HadronSpectrum_MultiQprop.hpp"
#include "Fopr/foprHermFuncFactoryCreator.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "Measurements/FermionicM/outputUtils.hpp"
#include <typeinfo>

using namespace std;
int Test_HadronSpectrum_MultiQprop::run(){

  InputConfig config = input_.getConfig();

  /////////////// Creation of factory & objects //////////////
  XML::node node = input_.node;
  XML::descend(node,"Source",MANDATORY);              // basic source 
  auto_ptr<SourceFactory> 
    SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());    
  CCIO::cout<<"base source ready\n";
  
  XML::next_sibling(node,"SourceSinkSmearedPropagator",MANDATORY);
  XML::descend(node,"QuarkPropagator",MANDATORY);     // Quark Propagator(smeared)
  auto_ptr<QuarkPropagatorFactory> 
    qps(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_smrd(qps->getQuarkProp(config));
  CCIO::cout<<"QuarkPropagator(smrd) ready\n";

  vector<FoprHermFactory*> sff;                       // Smearing operators
  for(XML::node snode=node.next_sibling("Smearing");    
      snode; 
      snode=snode.next_sibling("Smearing"))
    sff.push_back(FuncHermite::createHermOpFuncFactory(snode));
  vector<Fopr_Herm*> fs;                           
  for(int it=0; it<sff.size(); ++it) fs.push_back(sff[it]->getFoprHerm(config));    
  CCIO::cout<<"Smearing operators ready\n";

  node = input_.node;  
  XML::descend(node,"MultiMassPropagator",MANDATORY); // Quark Propagator(non-smeared)
  vector<double> mass_list;                          
  XML::read_array(node,"masses",mass_list,MANDATORY); // quark mass list
  CCIO::cout<<"mass-list ready\n";

  XML::descend(node,"QuarkPropagator",MANDATORY);

  ///////////// Actual calculation & output ///////////////
  for(int mm=0; mm<mass_list.size(); ++mm){
    CCIO::cout<<"*=============*\n";
    CCIO::cout<<"  mass = "<<mass_list[mm]<<" \n";
    CCIO::cout<<"*=============*\n";

    auto_ptr<QuarkPropagatorFactory> 
      mff(QuarkPropagators::createQuarkPropagatorFactory(node,mass_list[mm])); 
    prop_t qpBase;
    auto_ptr<QuarkPropagator> qprop_mmass(mff->getQuarkProp(config));// basic prop
    CCIO::cout<<"QuarkPropagator(non-smrd) ready\n";
    
    qprop_mmass->calc(qpBase,*src);

    stringstream outputfile;
    outputfile<<input_.output<<"_m"<<mm<<flush;
    ofstream writer(outputfile.str().c_str());           // creating the output stream 

    CCIO::cout<<"\n ::::::::::::: Output saving in "<< outputfile.str().c_str()<<"\n";

    for(int s0=0; s0<fs.size(); ++s0){
      Source_wrapper<Format::Format_F> sw(src.get(),fs[s0]);  // source smearing
      
      CCIO::cout<<" [solving at source smearing "<<s0<<" ]\n";
      prop_t qpSrc;
      qprop_smrd->calc(qpSrc,sw);
    
      for(int s1=0; s1<fs.size(); ++s1){
	CCIO::cout<<" [sink smearing "<<s1<<" ]\n";
	prop_t qpSmrd;                                        // sink smearing
	for(int cs=0; cs<qpSrc.size(); ++cs) qpSmrd.push_back(fs[s1]->mult(qpSrc[cs]));

	writer<<":::::: (sink, source) = ("<<s1<<", "<<s0<<") ::::::\n";
	Hadrons::mesonProp(qpSmrd,qpBase,writer);
	Hadrons::baryonProp(qpSmrd,qpBase,writer);
	Hadrons::mesonExtraProp(qpSmrd,qpBase,writer);
      }
    }
  }
  for(int sr=0; sr<fs.size(); ++sr){ delete fs[sr]; delete sff[sr]; }
  return 0;
}

