/*!
 * @file test_MesonSpectrum.cpp
 * @brief implementation of test_MesonSpectrum class
 */
#include "test_MesonSpectrum.hpp"
#include "include/factories.hpp"
#include "Measurements/FermionicM/quark_prop_meas_factory.hpp"
#include "Measurements/FermionicM/qprop_mom.hpp"
#include "Measurements/FermionicM/meson_correlator.hpp"
#include "include/errors.hpp"
#include <memory>

using namespace std;

int get_Direction(XML::node node) {
  if (!XML::attribute_compare(node, "dir", "X"))
    return XDIR;
  if (!XML::attribute_compare(node, "dir", "Y"))
    return YDIR;
  if (!XML::attribute_compare(node, "dir", "Z"))
    return ZDIR;
  if (!XML::attribute_compare(node, "dir", "T"))
    return TDIR;

  ErrorString msg;
  msg << "Direction not valid\n" << 
    "Specify the attribute chosing in {X,Y,Z,T}";
  Errors::XMLerr(msg);
}


int Test_MesonSpectrum::run(){
  XML::node node = input_.node;
  //// Quark Propagator ////
  XML::descend(node,"QuarkProp");
  XML::node qprop_node = node;

  // Select direction
  int dir =  get_Direction(qprop_node);
  CCIO::cout << "Selected dir="<<dir << "\n";
  
  auto_ptr<QuarkPropagatorFactory> qpfact(QuarkPropagators::createQuarkPropagatorFactory(node));
  InputConfig config = input_.getConfig();
  auto_ptr<QuarkPropagator> qprop(qpfact->getQuarkProp(config));
  
  //// source creation ////
  XML::next_sibling(node,"Source");
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<SiteIndex,Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  // meson correlators
  CCIO::cout << " ---- Making up meson correlators in direction "<<dir<<"\n";
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  vector<double> Cpp   = pp.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq,sq,dir);  

  // output
  CCIO::cout << " ---- Output in "<< input_.output.c_str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(input_.output.c_str());

    writer<< setiosflags(  ios_base::scientific);
    writer<<"---pp meson correlator---"<<endl;
    for(int t=0; t<Cpp.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cpp[t]
	    << endl;

    writer<<"---v1v1 meson correlator---"<<endl;
    for(int t=0; t<Cv1v1.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv1v1[t]
	    << endl;

    writer<<"---v2v2 meson correlator---"<<endl;
    for(int t=0; t<Cv2v2.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv2v2[t]
	    << endl;

    writer<<"---v3v3 meson correlator---"<<endl;
    for(int t=0; t<Cv3v3.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cv3v3[t]
	    << endl;
    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  
  return 0;
}
