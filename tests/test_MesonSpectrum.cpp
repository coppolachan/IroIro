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

bool is_ScreeningMass(XML::node node) {
  if (!XML::attribute_compare(node, "type", "Mass"))
    return false;
  if (!XML::attribute_compare(node, "type", "ScreeningMass"))
    return true;

  ErrorString msg;
  msg << "Type not set\n" << 
    "Specify the attribute type chosing in {Mass, ScreeningMass}";
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
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  // meson correlators
  CCIO::cout << " ---- Making up meson correlators \n";
  MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  MesonCorrelator ss(Scalar), Av1Av1(AVector1),Av2Av2(AVector2),Av3Av3(AVector2);


  
  // assumes that all the spatial sizes are the same
  int space_size = CommonPrms::instance()->global_size(XDIR);
  vector<double> temp, temp2, Cpp(space_size, 0.0), Css(space_size,0.0);
  for (int dir1 = 0; dir1 < 3; dir1++){
    vector<double> temp  = pp.calculate<Format::Format_F>(sq,sq,dir1);  
    vector<double> temp2 = ss.calculate<Format::Format_F>(sq,sq,dir1);  
    for (int t=0; t< space_size; t++){
      Cpp[t] +=  temp[t]*0.333333333333333333333;
      Css[t] += temp2[t]*0.333333333333333333333;
    }
  }
  // Susceptibility
  double chi_pi, chi_delta;
  for (int t=0; t< space_size; t++){
    chi_pi += Cpp[t];
    chi_delta -= Css[t];
  }

  CCIO::cout << "Susceptibilities:\n";
  CCIO::cout << "---- chi_pi    :  "<< chi_pi << "\n";
  CCIO::cout << "---- chi_delta :  "<< chi_delta << "\n";

  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq,sq,dir);  

  vector<double> CAv1v1 = Av1Av1.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> CAv2v2 = Av2Av2.calculate<Format::Format_F>(sq,sq,dir);  
  vector<double> CAv3v3 = Av3Av3.calculate<Format::Format_F>(sq,sq,dir);  


  // output
  CCIO::cout << " ---- Output written in "<< input_.output <<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(input_.output.c_str());

    writer<< setiosflags(  ios_base::scientific);
    writer<<"---Pseudoscalar (connected) meson correlator---\n";
    for(int t=0; t<Cpp.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cpp[t]
	    << endl;

    writer<<"---Scalar (connected) meson correlator---\n";
    for(int t=0; t<Cpp.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Css[t]
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





   writer<<"---Av1v1 meson correlator---"<<endl;
    for(int t=0; t<Cv1v1.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << CAv1v1[t]
	    << endl;

    writer<<"---Av2v2 meson correlator---"<<endl;
    for(int t=0; t<Cv2v2.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << CAv2v2[t]
	    << endl;

    writer<<"---Av3v3 meson correlator---"<<endl;
    for(int t=0; t<Cv3v3.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << CAv3v3[t]
	    << endl;


    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  
  return 0;
}
