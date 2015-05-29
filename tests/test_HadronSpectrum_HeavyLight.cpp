/*!
 * @file test_HadronSpectrum_HeavyLight.cpp
 * @brief Implementation of test_HadronSpectrum_HeavyLight class
 * This class obtains hadron correlators from combinations of 
 * Source-smeared propagators and multi-mass propagators
 * Time-stamp: <2015-05-20 17:19:21 noaki>
 */
#include "factories.hpp"
#include "messages_macros.hpp"
#include "test_HadronSpectrum_HeavyLight.hpp"
#include "Fopr/foprHermFuncFactoryCreator.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "Measurements/FermionicM/baryonCorrelator.hpp"
#include "Measurements/FermionicM/outputUtils.hpp"
#include "timings.hpp"
#include <typeinfo>

using namespace std;
int Test_HadronSpectrum_HeavyLight::run(){

  InputConfig config = input_.config;

  XML::node node = input_.node;
  vector<vector<int> > lpos;
  XML::descend(node,"LocalSource",MANDATORY);           // basic source 
  for(XML::iterator it=node.begin(); it!=node.end(); ++it)
    lpos.push_back(XML::child_array(*it));

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
  
  for(int lp=0; lp<lpos.size(); ++lp){   /// loop over source locations
    CCIO::cout<<"Source location = ("
	      <<lpos[lp][XDIR]
	      <<","      <<lpos[lp][YDIR]
	      <<","      <<lpos[lp][ZDIR]
	      <<","      <<lpos[lp][TDIR]
	      <<")\n";

    Source_local<Format::Format_F> sloc(lpos[lp]);
    
    stringstream outputfile;
    outputfile<<input_.output
	      <<"_ll_src"<<lpos[lp][XDIR]
	      <<"_"      <<lpos[lp][YDIR]
	      <<"_"      <<lpos[lp][ZDIR]
	      <<"_"      <<lpos[lp][TDIR]<<flush;

    ofstream writerLL(outputfile.str().c_str());              // output stream 

    CCIO::cout<<":::::: light-light meson correlators ::::::\n";
    
    vector<vector<prop_t> > qpUpDown(fs.size());              // smeared u/d props
    vector<vector<prop_t> > qpStrange(fs.size());             // smeared strange props
    long double propagator_timer, source_timer, save_timer, smearing_timer;

    for(int s0=0; s0<fs.size(); ++s0){
      FINE_TIMING_START(source_timer);
      Source_wrapper<Format::Format_F> sw(&sloc,fs[s0]);      // source smearing
      FINE_TIMING_END(source_timer);
      CCIO::cout << " Timing - Source         :  "<< source_timer << " seconds\n";

      FINE_TIMING_START(propagator_timer);                                       
      CCIO::cout<<" [up/down-propagator, smearing "<<s0<<" ]\n";
      prop_t qpud;  qprop_ud->calc(qpud,sw);                  // actual calc (u/d)
      FINE_TIMING_END(propagator_timer);
      CCIO::cout << " Timing - Propagator ud  :  "<< propagator_timer << " seconds\n";

      FINE_TIMING_START(propagator_timer);                
      CCIO::cout<<" [strange-propagator, smearing "<<s0<<" ]\n";
      prop_t qps;   qprop_s ->calc(qps,sw);                   // actual calc (strange)
      FINE_TIMING_END(propagator_timer);
      CCIO::cout << " Timing - Propagator s   :  "<< propagator_timer << " seconds\n";
      
      CCIO::cout<<"\n ::::::::::: Output saving in "<< outputfile.str().c_str()<<"\n";
      for(int s1=0; s1<fs.size(); ++s1){
	prop_t qpSud, qpSst;                                 // sink smearing
	CCIO::cout<<"\n ::::::::::: Smearing sink... "; 
	FINE_TIMING_START(smearing_timer);
	for(int cs=0; cs<qpud.size(); ++cs){
	  qpSud.push_back(fs[s1]->mult(qpud[cs]));
	  qpSst.push_back(fs[s1]->mult(qps[cs]));
	}    
	FINE_TIMING_END(smearing_timer);
	CCIO::cout << " done in "<< smearing_timer << " seconds\n";
  
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
	
    
	FINE_TIMING_START(save_timer);  
	CCIO::cout<<"\n ::::::::::: Saving smeared propagators... "; 
	qpUpDown[s0].push_back(qpSud);
	qpStrange[s0].push_back(qpSst);     
	FINE_TIMING_END(save_timer);
	CCIO::cout<< " done in "<< save_timer <<" seconds\n";
      }
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
      qprop_mmass->calc(qpHeavy,sloc);
    
      outputfile<<input_.output<<"_m"<<mm
		<<"_src"<<lpos[lp][XDIR]
		<<"_"   <<lpos[lp][YDIR]
		<<"_"   <<lpos[lp][ZDIR]
		<<"_"   <<lpos[lp][TDIR]<<flush;

      ofstream writerHL(outputfile.str().c_str());      // output stream 
    
      CCIO::cout<<"\n ::::::::::: Output saving in "<< outputfile.str().c_str()<<"\n";
      
      for(int s0=0; s0<fs.size(); ++s0){                 // smeared-unsmeared
	for(int s1=0; s1<fs.size(); ++s1){
	
	  CCIO::cout<<" [contraction with (sink, source) = ("<<s1<<", "<< s0<<")]\n";

	  if(Communicator::instance()->primaryNode())
	    writerHL<<":::::: heavy-ud mesons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	  Hadrons::mesonPropGeneral(qpUpDown[s0][s1],qpHeavy,writerHL);    
	  
	  if(Communicator::instance()->primaryNode())
	    writerHL<<":::::: heavy-s mesons  (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	  Hadrons::mesonPropGeneral(qpStrange[s0][s1],qpHeavy,writerHL);    
	  
	  if(Communicator::instance()->primaryNode())
	    writerHL<<":::::: heavy-ud baryons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	  Hadrons::baryonProp(qpUpDown[s0][s1],qpHeavy,writerHL);

	  if(Communicator::instance()->primaryNode())
	    writerHL<<":::::: heavy-s baryons (sink,source) = ("<<s1<<","<<s0<<") ::::::\n";
	  Hadrons::baryonProp(qpStrange[s0][s1],qpHeavy,writerHL);
	}
      }
      if(Communicator::instance()->primaryNode())
	writerHL<<":::::: heavy-heavy mesons ::::::\n";	
      Hadrons::mesonPropGeneral(qpHeavy,qpHeavy,writerHL);    
    
      if(Communicator::instance()->primaryNode())
	writerHL<<":::::: heavy-heavy baryons ::::::\n";	
      Hadrons::baryonProp(qpHeavy,qpHeavy,writerHL);
      outputfile.str("");
    }
    
    /// garbage collection ///
    for(int sr=0; sr<fs.size(); ++sr){ delete fs[sr]; delete sff[sr]; }
    return 0;
  }
}

