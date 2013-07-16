/*!
 * @file test_MesonSpectrum.cpp
 * @brief implementation of test_MesonSpectrum class
 */
#include "test_MesonSpectrum_SpMom.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
//#include "Measurements/GaugeM/staples.hpp"
#include <memory>

using namespace std;

int Test_MesonSpectrum_SpMom::run(){
  XML::node node = input_.node;
  InputConfig config = input_.getConfig();  

  //// Quark Propagator ////
  XML::descend(node,"QuarkProp");
  
  auto_ptr<QuarkPropagatorFactory> qpfact(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop(qpfact->getQuarkProp(config));
  
  //// source creation ////
  XML::next_sibling(node,"Source");
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<SiteIndex,Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);


  // meson correlators
  CCIO::cout << " ---- Making up meson correlators\n";
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  
  // output file open
  ofstream writer(input_.output.c_str());  

  // Fourier Tr
  int Lx = CommonPrms::instance()->Lx();
  int Ly = CommonPrms::instance()->Ly();
  int Lz = CommonPrms::instance()->Lz();
  const int Mx = Lx/2/PI;
  const int My = Ly/2/PI;
  const int Mz = Lz/2/PI;

  for(int nx = -Mx;nx<=Mx;++nx){
    for(int ny = -My;ny<=My;++ny){
      for(int nz = -Mz;nz<=Mz;++nz){

	double px = 2*PI*nx/Lx;
	double py = 2*PI*ny/Ly;
	double pz = 2*PI*nz/Lz;

	vector<double> Cpp = pp.calculate_mom<Format::Format_F>(sq,sq,px,py,pz);
	//	vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq,px,py,pz);  
	//	vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq,sq,px,py,pz);  
	//	vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq,sq,px,py,px);

	// output
	//CCIO::cout << " ---- Output p=("<<px<<","<<py<<","<<pz<<") in "<< input_.output.c_str()<<"\n";

	if(Communicator::instance()->primaryNode()){
	  CCIO::cout << "calculation p=("<<px<<","<<py<<","<<pz<<")\n";

	  
	  writer<< setiosflags(  ios_base::scientific);
	  writer<<"---pp meson correlator in p=("<<px<<","<<py<<","<<pz<<")"<<endl;
	  writer<<px<<"  "<<py<<"  "<<pz<<endl;    
	  for(int t=0; t<Cpp.size()/2; ++t){
	    writer<< setw(2) <<setiosflags(ios_base::right)<< t;
	    writer<< setw(20)<<setiosflags(ios_base::left )<< Cpp[2*t];
	    writer<< setw(20)<<setiosflags(ios_base::left )<< Cpp[2*t+1]<<endl;
	  }
	  /*    writer<<"---v1v1 meson correlator---"<<endl;
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
		}*/
	  writer<< resetiosflags(ios_base::scientific);
	  writer<<endl;
	}
      }
    }
  }
  writer.close();
  Communicator::instance()->sync();
  
  return 0;
}
