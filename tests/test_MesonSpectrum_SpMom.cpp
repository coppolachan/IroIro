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
  //// Quark Propagator ////
  XML::descend(node,"QuarkProp");
  
  auto_ptr<QuarkPropagatorFactory> qpfact(QuarkPropagators::createQuarkPropagatorFactory(node));
  auto_ptr<QuarkPropagator> qprop(qpfact->getQuarkProp(input_.config));
  
  //// source creation ////
  XML::next_sibling(node,"Source");
  auto_ptr<SourceFactory> SrcFactory(Sources::createSourceFactory<Format::Format_F>(node));
  auto_ptr<Source> src(SrcFactory->getSource());

  prop_t sq;  //Defines a vector of fields
  CCIO::cout << " ---- Calculating propagator\n";
  qprop->calc(sq,*src);

  // meson correlators
  CCIO::cout << " ---- Making up meson correlators\n";
  MesonCorrelator pp(Pion), sc(Scalar);
  MesonCorrelator v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);
  MesonCorrelator av1av1(AVector1),av2av2(AVector2),av3av3(AVector3),av4av4(AVector4);
  MesonCorrelator t12t12(Tensor12), t13t13(Tensor13),t23t23(Tensor23);
  
  // output file open
  ofstream writer(input_.output.c_str());  
  writer<< setiosflags(ios_base::scientific);
  
  // Fourier Tr
  int Lx = CommonPrms::instance()->Lx();
  int Ly = CommonPrms::instance()->Ly();
  int Lz = CommonPrms::instance()->Lz();
  const int Mx = Lx/2/PI;
  const int My = Ly/2/PI;
  const int Mz = Lz/2/PI;

  writer<<"---pseudo scalar meson correlator"<<endl;
  /*
  for(int nx=-Mx; nx<=Mx; ++nx){
    double px = 2*PI*nx/Lx;	

    for(int ny=-My; ny<=My; ++ny){
      double py = 2*PI*ny/Ly;	

      for(int nz=-Mz; nz<=Mz; ++nz){
	double pz = 2*PI*nz/Lz;
	vector<double> Cpp = pp.calculate_PSmom<Format::Format_F>(sq,sq,px,py,pz);
  */
  for(int nx=0; nx<=2; ++nx){
    double px = 2*PI*nx/Lx;	

    for(int ny=0; ny<=2; ++ny){
      double py = 2*PI*ny/Ly;	

      for(int nz=0; nz<=2; ++nz){
	double pz = 2*PI*nz/Lz;

	vector<double> Cpp = pp.calculate_PSmom<Format::Format_F>(sq,sq,px,py,pz);

	if(Communicator::instance()->primaryNode()){
	  CCIO::cout<<"momentum (nx,ny,nz) = ("<<nx<<","<<ny<<","<<nz<<")"<<endl;

	  writer<<"---momentum index---"<<endl;
	  writer<<setw(4)<<setiosflags(ios_base::right)<<nx;
	  writer<<setw(4)<<setiosflags(ios_base::right)<<ny;
	  writer<<setw(4)<<setiosflags(ios_base::right)<<nz<<endl;
	  writer<<"---momentum---"<<endl;
	  writer<<setw(20)<<setprecision(8)<<setiosflags(ios_base::right)<<px;
	  writer<<setw(20)<<setprecision(8)<<setiosflags(ios_base::right)<<py;
	  writer<<setw(20)<<setprecision(8)<<setiosflags(ios_base::right)<<pz<<endl;

	  for(int t=0; t<Cpp.size(); ++t)
	    writer<< setw(2) <<setiosflags(ios_base::right)<< t
		  << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
		  << Cpp[t]
		  << endl;
	}
      }
    }
  }

  /// vector mesons
  vector<double> Cv1v1 = v1v1.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv2v2 = v2v2.calculate<Format::Format_F>(sq,sq);  
  vector<double> Cv3v3 = v3v3.calculate<Format::Format_F>(sq,sq);

  if(Communicator::instance()->primaryNode()){
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
  }
  
  /// scalar mesons
  vector<double> Csc = sc.calculate<Format::Format_F>(sq,sq);
  if(Communicator::instance()->primaryNode()){
    writer<<"---scalar meson correlator---"<<endl;
    for(int t=0; t<Csc.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Csc[t]
	    << endl;
  }
  
  /// axial vector mesons
  vector<double> Cav1av1 = av1av1.calculate<Format::Format_F>(sq,sq);
  vector<double> Cav2av2 = av2av2.calculate<Format::Format_F>(sq,sq);
  vector<double> Cav3av3 = av3av3.calculate<Format::Format_F>(sq,sq);
  vector<double> Cav4av4 = av4av4.calculate<Format::Format_F>(sq,sq);

  if(Communicator::instance()->primaryNode()){
    writer<<"---av1av1 meson correlator---"<<endl;
    for(int t=0; t<Cav1av1.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cav1av1[t]
	    << endl;
    
    writer<<"---av2av2 meson correlator---"<<endl;
    for(int t=0; t<Cav2av2.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cav2av2[t]
	    << endl;

    writer<<"---av3av3 meson correlator---"<<endl;
    for(int t=0; t<Cav3av3.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cav3av3[t]
	    << endl;

    writer<<"---av4av4 meson correlator---"<<endl;
    for(int t=0; t<Cav4av4.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Cav4av4[t]
	    << endl;
  }

  /// tensor mesons
  vector<double> Ct12t12 = t12t12.calculate<Format::Format_F>(sq,sq);
  vector<double> Ct13t13 = t13t13.calculate<Format::Format_F>(sq,sq);
  vector<double> Ct23t23 = t23t23.calculate<Format::Format_F>(sq,sq);
  if(Communicator::instance()->primaryNode()){
    writer<<"---t12t12 meson correlator---"<<endl;
    for(int t=0; t<Ct12t12.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Ct12t12[t]
	    << endl;
    
    writer<<"---t13t13 meson correlator---"<<endl;
    for(int t=0; t<Ct13t13.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Ct13t13[t]
	    << endl;

    writer<<"---t23t23 meson correlator---"<<endl;
    for(int t=0; t<Ct23t23.size(); ++t)
      writer<< setw(2) <<setiosflags(ios_base::right)<< t
	    << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	    << Ct23t23[t]
	    << endl;
  }
  CCIO::cout<<"finished making up other meson correlators without momenta"<<endl;

  writer<< resetiosflags(ios_base::scientific);
  writer.close();
  Communicator::instance()->sync();
  
  return 0;
}
