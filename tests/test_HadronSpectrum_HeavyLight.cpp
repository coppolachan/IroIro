/*!
 * @file test_HadronSpectrum_HeavyLight.cpp
 * @brief Implementation of test_HadronSpectrum_HeavyLight class
 * This class obtains hadron correlators from combinations of 
 * Source-smeared propagators and multi-mass propagators
 * Time-stamp: <2014-08-07 13:29:07 noaki>
 */
#include "include/factories.hpp"
#include "include/messages_macros.hpp"
#include "test_HadronSpectrum_HeavyLight.hpp"
#include "Fopr/foprHermFuncFactoryCreator.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "Measurements/FermionicM/outputUtils.hpp"
#include <typeinfo>

using namespace std;
int Test_HadronSpectrum_HeavyLight::run(){

  InputConfig config = input_.getConfig();

  XML::node node = input_.node;
  XML::descend(node,"Source",MANDATORY);             // basic source 
  auto_ptr<SourceFactory> 
    SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());    
  CCIO::cout<<"base source ready\n";

  /*=-=-=-=-=-=-=-= light hadron propagators =-=-=-=-=-=-=-=-=-=*/  
  XML::next_sibling(node,"SmearedPropagator",MANDATORY);

  XML::descend(node,"UpDownPropagator",MANDATORY);      // u/d propagator(smeared)
  auto_ptr<QuarkPropagatorFactory> 
    qpudF(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_ud(qpudF->getQuarkProp(config));
  CCIO::cout<<"u/d propagator(smrd) ready\n";

  XML::next_sibling(node,"StrangePropagator",MANDATORY);// strange propagator(smeared)
  auto_ptr<QuarkPropagatorFactory> 
    qpsF(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop_s(qpsF->getQuarkProp(config));
  CCIO::cout<<"strange propagator(smrd) ready\n";

  vector<FoprHermFactory*> sff;                        // Smearing operators
  for(XML::node snode=node.next_sibling("Smearing");    
      snode; 
      snode=snode.next_sibling("Smearing"))
    sff.push_back(FuncHermite::createHermOpFuncFactory(snode));
  vector<Fopr_Herm*> fs;                           
  for(int it=0; it<sff.size(); ++it) fs.push_back(sff[it]->getFoprHerm(config));    
  CCIO::cout<<"Smearing operators ready\n";

  stringstream outputfile;
  outputfile<<input_.output<<"_ll"<<flush;
  ofstream writerLL(outputfile.str().c_str());             // output stream 

  CCIO::cout<<":::::: light-light meson correlators ::::::\n";
  
  vector<prop_t> qpUpDown;                                  // smeared u/d props
  vector<prop_t> qpStrange;                                 // smeared strange props
  for(int s0=0; s0<fs.size(); ++s0){
    Source_wrapper<Format::Format_F> sw(src.get(),fs[s0]);  // source smearing
                                                           
    CCIO::cout<<" [up/down-propagator, smearing "<<s0<<" ]\n";
    prop_t qpud;  qprop_ud->calc(qpud,sw);                  // actual calc (u/d)
    CCIO::cout<<" [strange-propagator, smearing "<<s0<<" ]\n";
    prop_t qps;   qprop_s ->calc(qps,sw);                   // actual calc (strange)

    CCIO::cout<<"\n ::::::::::: Output saving in "<< outputfile.str().c_str()<<"\n";
    for(int s1=0; s1<fs.size(); ++s1){
      prop_t qpSud, qpSst;                                 // sink smearing
      for(int cs=0; cs<qpud.size(); ++cs){
	qpSud.push_back(fs[s1]->mult(qpud[cs]));
	qpSst.push_back(fs[s1]->mult(qps[cs]));
      }      
      CCIO::cout<<" [contraction with (sink, source) = ("<<s1<<", "<< s0<<")]\n";

      if(Communicator::instance()->primaryNode())
	writerLL<<":::::: ud-ud mesons (sink, source) = ("<<s1<<","<<s0<<") ::::::\n";
      Hadrons::mesonPropGeneral(qpSud,qpSud,writerLL);    
      
      if(Communicator::instance()->primaryNode())
	writerLL<<":::::: ud-s mesons (sink, source) = ("<<s1<<","<<s0<<") ::::::\n";
      Hadrons::mesonPropGeneral(qpSud,qpSst,writerLL);

      if(Communicator::instance()->primaryNode())
	writerLL<<":::::: s-s mesons (sink, source) = ("<<s1<<","<<s0<<") ::::::\n";
      Hadrons::mesonPropGeneral(qpSst,qpSst,writerLL);    
      
      if(Communicator::instance()->primaryNode())
	writerLL<<":::::: baryons (sink, source) = ("<<s1<<","<<s0<<") ::::::\n";
      Hadrons::baryonProp(qpSud,qpSst,writerLL);
    }
    qpUpDown.push_back(qpud);  // saving solved propagators
    qpStrange.push_back(qps);  
  }
  outputfile.str("");

  /*=-=-=-=-=-=-=-= heavy-light hadron propagators =-=-=-=-=-=-=-=-=-=*/  
  node = input_.node;  
  XML::descend(node,"MultiMassPropagator",MANDATORY);// Quark Propagator(non-smeared)
  vector<double> mass_list;                          
  XML::read_array(node,"masses",mass_list,MANDATORY);// quark mass list
  CCIO::cout<<"mass-list ready\n";

  XML::descend(node,"HeavyQuarkPropagator",MANDATORY);
  
  CCIO::cout<<":::::: heavy-light meson correlators ::::::\n";

  for(int mm=0; mm<mass_list.size(); ++mm){  

    auto_ptr<QuarkPropagatorFactory> 
      mff(QuarkPropagators::createQuarkPropagatorFactory(node,mass_list[mm])); 
    auto_ptr<QuarkPropagator> qprop_mmass(mff->getQuarkProp(config));// basic prop
    CCIO::cout<<"QuarkPropagator(non-smrd) ready\n";

    CCIO::cout<<" [solving at mass "<<mm<<" ]\n";    
    prop_t qpHeavy;
    qprop_mmass->calc(qpHeavy,*src);
    
    outputfile<<input_.output<<"_m"<<mm<<flush;
    ofstream writerHL(outputfile.str().c_str());      // output stream 
    
    CCIO::cout<<"\n ::::::::::: Output saving in "<< outputfile.str().c_str()<<"\n";

    for(int s0=0; s0<fs.size(); ++s0){                 // smeared-unsmeared
      for(int s1=0; s1<fs.size(); ++s1){
	prop_t qpSud, qpSst;                           // sink smearing
	for(int cs=0; cs<qpUpDown[s0].size(); ++cs){
	  /*
	  Field tmp = fs[s1]->mult(qpUpDown[s0][cs]);
	  tmp /= tmp.norm();
	  qpSud.push_back(tmp);

	  tmp = fs[s1]->mult(qpStrange[s0][cs]);
	  tmp /= tmp.norm();
	  qpSst.push_back(tmp);
	  */
	  qpSud.push_back(fs[s1]->mult(qpUpDown[s0][cs]));
	  qpSst.push_back(fs[s1]->mult(qpStrange[s0][cs]));
	}
	CCIO::cout<<" [contraction with (sink, source) = ("<<s1<<", "<< s0<<")]\n";

	if(Communicator::instance()->primaryNode())
	  writerHL<<":::::: heavy-ud mesons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	Hadrons::mesonPropGeneral(qpSud,qpHeavy,writerHL);    

	if(Communicator::instance()->primaryNode())
	  writerHL<<":::::: heavy-s mesons  (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	Hadrons::mesonPropGeneral(qpSst,qpHeavy,writerHL);    

	if(Communicator::instance()->primaryNode())
	  writerHL<<":::::: heavy-ud baryons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	Hadrons::baryonProp(qpSud,qpHeavy,writerHL);

	if(Communicator::instance()->primaryNode())
	  writerHL<<":::::: heavy-s baryons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	Hadrons::baryonProp(qpSst,qpHeavy,writerHL);
      }
    }
    if(Communicator::instance()->primaryNode())
      writerHL<<":::::: heaby-heaby mesons ::::::\n";	
    Hadrons::mesonPropGeneral(qpHeavy,qpHeavy,writerHL);    
    
    if(Communicator::instance()->primaryNode())
      writerHL<<":::::: heaby-heaby baryons ::::::\n";	
    Hadrons::baryonProp(qpHeavy,qpHeavy,writerHL);
    outputfile.str("");
  }

  /// gavage collection ///
  for(int sr=0; sr<fs.size(); ++sr){ delete fs[sr]; delete sff[sr]; }
  return 0;
}

