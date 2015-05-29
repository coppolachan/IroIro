#include "test_QnumSuscept.hpp"
#include "include/field.h"
#include "Measurements/FermionicM/sources_factory.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Solver/solver_Factory.hpp"
#include "PugiXML/xmlUtilities.hpp"
#include <memory>

using namespace std;

int Test_QnumSuscept::run(){
  CCIO::cout<<"Test_QnumSuscept::run() called\n";

  XML::node qnsus_node= input_.node;
  XML::descend(qnsus_node,"QuarkNumberSusceptibility");

  //////// Dirac Operator ////////
  XML::node dirac_node= qnsus_node;
  XML::descend(dirac_node,"DiracWFD");

  auto_ptr<DiracWilsonLikeFiniteDensityFactory> 
    Dfact(Diracs::createDiracWilsonLikeFiniteDensityFactory(dirac_node));
  
  auto_ptr<DiracWilsonLikeFiniteDensity> D(Dfact->getDirac(input_.config));

  CCIO::cout<<"Dirac operator created\n";

  ////////// Solver ///////////
  XML::node slv_node= qnsus_node;
  XML::descend(slv_node,"Solver");

  auto_ptr<SolverFactory> slv_fact(Solvers::createSolverFactory(slv_node));
  Fopr_DdagD DdagD(D.get());
  auto_ptr<Solver> slv(slv_fact->getSolver(&DdagD));

  CCIO::cout<<"Solver created\n";
  
  //////// Noise Source /////////
  int Nnoise;
  XML::read(qnsus_node,"Nnoise",Nnoise,MANDATORY);

  XML::descend(qnsus_node,"NoiseSource");  
  auto_ptr<SourceFactory> 
    src_fact(Sources::createSourceFactory<Format::Format_F>(qnsus_node));
  auto_ptr<Source> src(src_fact->getSource());

  CCIO::cout<<"Noise source created\n";

  /////// Actual Calculation ///////  
  const char* dataName[]= {"Di",
			   "DsDi",
			   "DtpDi",
			   "DtmDi",
			   "DexDi",
			   "DiDtmDi",
			   "DsDiDtmDi",
			   "DtpDiDtmDi",
			   "DtmDiDtmDi",
			   "DexDiDtmDi"};
  vector<complex<double> > eXe(10);

  std::stringstream ofile;
  ofile << input_.output.c_str();
  ofstream writer(ofile.str().c_str());

  for(int ni=0; ni<Nnoise; ++ni){ // loop over noise 

    src->refresh();
    for(int e=0;e<eXe.size();++e) eXe[e] = 0.0;

    CCIO::cout<<"[ source "<<ni<<" ]\n";
    
    for(int s=0; s<ND_; ++s){  // full color-spin dilution
      for(int c=0; c<NC_; ++c){  		  
	
	CCIO::cout<<"s="<<s<<" c="<<c<<"\n";	
	
	Field eta = src->mksrc(s,c);
	Field Mi_e(eta.size());

	//// one-inverce 
	Field tmp = D->mult_dag(eta);
	SolverOutput info =  slv->solve(Mi_e,tmp); // D^{-1}*eta

	CCIO::cout<<" 1st inversion: num of iter= "<<info.Iterations
		  <<" residual= "<<info.diff
		  <<" timing (msec)= "<<info.timing<<"\n";

	eXe[0] += complex<double>(eta*Mi_e,eta.im_prod(Mi_e));
	
	tmp = D->mult_Ds(Mi_e);    // Ds*D^{-1}*eta
	eXe[1] += complex<double>(eta*tmp,eta.im_prod(tmp));

	tmp = D->mult_Dtp(Mi_e);   // Dtp*D^{-1}*eta
	eXe[2] += complex<double>(eta*tmp,eta.im_prod(tmp));
	
	tmp = D->mult_Ex(Mi_e);  // Dex*D^{-1}*eta, must precede eXe[3]
	eXe[4] += complex<double>(eta*tmp,eta.im_prod(tmp));

	tmp = D->mult_Dtm(Mi_e);   // Dtm*D^{-1}*eta
	eXe[3] += complex<double>(eta*tmp,eta.im_prod(tmp));

	//// two-inverces
	Mi_e = tmp;
	tmp = D->mult_dag(Mi_e);
	info = slv->solve(Mi_e,tmp); // D^{-1}*Dtm*D^{-1}*eta

	CCIO::cout<<" 2nd inversion: num of iter= "<<info.Iterations
		  <<" residual= "<<info.diff
		  <<" timing (msec)= "<<info.timing<<"\n";
	
	eXe[5] += complex<double>(eta*Mi_e,eta.im_prod(Mi_e));

	tmp = D->mult_Ds(Mi_e);  // Ds*D^{-1}*Dtm*D^{-1}*eta
	eXe[6] += complex<double>(eta*tmp,eta.im_prod(tmp));

	tmp = D->mult_Dtp(Mi_e);  // Dtp*D^{-1}*Dtm*D^{-1}*eta
	eXe[7] += complex<double>(eta*tmp,eta.im_prod(tmp));

	tmp = D->mult_Dtm(Mi_e);  // Dtm*D^{-1}*Dtm*D^{-1}*eta
	eXe[8] += complex<double>(eta*tmp,eta.im_prod(tmp));

	tmp = D->mult_Ex(Mi_e); // Dex*D^{-1}*Dtm*D^{-1}*eta
	eXe[9] += complex<double>(eta*tmp,eta.im_prod(tmp));
      }
    } 
    ////// output /////
    CCIO::cout << " ---- Output in "<< ofile.str()<<"\n";

    if(Communicator::instance()->primaryNode()){
      writer<<"[noise "<<ni<<" ]\n";
      writer<<setiosflags(ios_base::scientific)
	    <<setiosflags(ios_base::left);

      for(int e=0; e<eXe.size(); ++e) 
	writer<< setw(12)<< dataName[e]
	      << setw(25)<< setprecision(16)<< eXe[e].real()
	      << setw(25)<< setprecision(16)<< eXe[e].imag()<<endl;
      writer<< resetiosflags(ios_base::scientific);
    } 
  }
  return 0;
}
