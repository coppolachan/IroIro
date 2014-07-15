/*!
 * @file test_ReadEigenModes.cpp
 * @brief Reads eigenmodes
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2014-06-12 14:30:41 neo>
 */

#include "test_ReadEigenModes.hpp"
#include "include/format_F.h"
#include "Measurements/GaugeM/staples.hpp"
#include "Measurements/GaugeM/topologyGeom.hpp"
#include "Measurements/GaugeM/polyakovLoop.hpp"
#include "Tools/sunMatUtils.hpp"
#include "Tools/fieldUtils.hpp"
#include "EigenModes/eigenModes.hpp"
#include <fstream> 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


////////////////////////////////////////
// IROIRO Uses DIRAC representation
////////////////////////////////////////

using namespace FieldUtils;
using namespace SUNmatUtils;

void print_va(std::valarray<double> va) {
  for (int i = 0; i < va.size(); i++){
    std::cout << i << "  "<< va[i] << "\n";
  }
}


double plaquette_site(const GaugeField& F, Staples stp, int Nvol_) {
  std::ofstream plaqf;
  CommonPrms* cmprms = CommonPrms::instance();
  int Nt = cmprms->Nt();
  int Nx = cmprms->Nx();
  int Ny = cmprms->Ny();
  int Nz = cmprms->Nz();

  double plaq = 0.0;
  double min_plaq = 0.0;
  int min_loc[4];
  GaugeField1D stpl(Nvol_);
  double plaqv[Nvol_];
  for(int i=0;i<NDIM_-1;++i){
    int j = (i+1)%(NDIM_-1);
    stpl = stp.lower(F,i,j);
    for(int site=0; site<Nvol_; ++site){
      int gsite = SiteIndex::instance()->get_gsite(site);
      plaq = ReTr(mat(F,site,i)*mat_dag(stpl,site));  // P_ij
      int x = SiteIndex::instance()->g_x(gsite);
      int y = SiteIndex::instance()->g_y(gsite);
      int z = SiteIndex::instance()->g_z(gsite);
      int t = SiteIndex::instance()->g_t(gsite);

      plaqv[site] += plaq;
      
      if (plaq < min_plaq){
	min_plaq = plaq;
	min_loc[0] = x;	    min_loc[1] = y;
	min_loc[2] = z;	    min_loc[3] = t;
      }
    }
  }

 for (int t = 0 ; t < Nt; t++){
    std::stringstream ss;
    ss << "ScalarGball_sites_t" << t; 
    plaqf.open(ss.str().c_str()); 

    
    plaqf << "x  y  z  norm\n";
    
    for (int x = 0 ; x < Nx; x++){
      for (int y = 0 ; y < Ny; y++){
	for (int z = 0 ; z < Nz; z++){
	  int site = SiteIndex::instance()->site(x,y,z,t);
	  plaqf << x << " " << y << " "
		<< z << " " << plaqv[site] << "\n";    
	}
      }
    }
    plaqf.close();
  }
 

  std::cout << "Plaq Min at "<<min_loc[0]<< " "<<min_loc[1] << " "<<min_loc[2] << " "<<min_loc[3]<<"  "<< min_plaq<<"\n";  

}





int Test_ReadEigenModes::run(){

  CCIO::cout<<"Test_ReadEigenModes::run() called\n";

  Communicator* comm = Communicator::instance();
  CommonPrms* cmprms = CommonPrms::instance();

  XML::node ReadEig_node = input_.node;
  InputConfig config = input_.getConfig();

  XML::descend(ReadEig_node,"ReadEigenMode",MANDATORY);

  int Nev; // stores  the number of eigenmodes to read
  int EigMode; // number of eigenmode to store
  std::string Emod_filename;
  std::string Eval_filename;

  // Check eigenvector file to read
  XML::read(ReadEig_node, "ModesFilename", Emod_filename,MANDATORY);
  XML::read(ReadEig_node, "EvalFilename", Eval_filename,MANDATORY);

  // Number of eigenvectors to read
  XML::read(ReadEig_node, "NumEigenModes",Nev, MANDATORY);

  XML::read(ReadEig_node, "EigenMode",EigMode, MANDATORY);


  // Initialization and checks
  // Lattice size parameters
  int Nvol = cmprms->Nvol();
  int Nt = cmprms->Nt();
  int Nx = cmprms->Nx();
  int Ny = cmprms->Ny();
  int Nz = cmprms->Nz();
  Field* u = &((input_.gconf)->data);
  Staples stpl;
  double plq = stpl.plaquette(*input_.gconf);
  CCIO::cout<<" Plaquette ="<< plq <<std::endl;
 
  //plaquette_site(*input_.gconf, stpl, Nvol);
 
  // Polyakov Loop field
  PolyakovLoop PLoop(TDIR);

  std::ofstream plf, plif, topf;
  std::stringstream sspl,sstop; 
  sspl << Emod_filename << "_PolyakovLoop"; 
  sstop << Emod_filename << "_TopCharge"; 
  plf.open(sspl.str().c_str());
  topf.open(sstop.str().c_str());
  sspl.str("");
  sspl << Emod_filename << "_PolyakovLoopIm";
  plif.open(sspl.str().c_str());

 
  GaugeField1D PL_field = PLoop.get_PLField(*input_.gconf);

  TopologyGeom topG; 
  std::vector<double> q_loc;
  topG.get_Q(q_loc, *input_.gconf);
  
  plf << "x  y  z  PL\n";
  topf << "x  y  z  TC\n";

  for(int s=0; s<SiteIndex::instance()->slsize(0,TDIR); ++s){
    SUNmat PL = mat(PL_field,s);
    double pl_re = SUNmatUtils::ReTr(PL);
    double pl_im = SUNmatUtils::ImTr(PL);
    int gsite = SiteIndex::instance()->get_gsite(s);
    int x = SiteIndex::instance()->g_x(gsite);
    int y = SiteIndex::instance()->g_y(gsite);
    int z = SiteIndex::instance()->g_z(gsite);

    plf << x << " "
	<< y << " "
	<< z << " "
	<< pl_re << " "
	<< pl_im << "\n";

    plif << x << " "
	 << y << " "
	 << z << " "
	 << pl_im <<"\n";
    
    topf << x << " "
	 << y << " "
	 << z << " "
	 << q_loc[s] << "\n";


  }
  plf.close();
  plif.close();
  topf.close();

  Format::Format_F FF(Nvol);

  // Allocate the required eigenvector fields 
  EigenModes EigM(ReadEig_node);//must specify threshold

  // Read the required eigenvectors and store in fields
  EigM.initialize< Format::Format_F >(Eval_filename, Emod_filename);




  // Calculate the local norm of the eigenvector
  CCIO::cout << "Norm: "<< EigM.evecs_[EigMode].norm() << "\n";


  for (int eval = 0; eval < EigM.evecs_.size(); eval++){
    std::ofstream outf, outc;
    
    double ipr=0.0; //inverse participatio ratio
    
    int max_loc[4];
    double max= 0.0;
    
    int min_loc[4];
    double min = 0.0;
    
    
    for (int t = 0 ; t < Nt; t++){
      std::stringstream ss;
      ss << Emod_filename << "_"<< eval << "_norm_t" << t; 
      
      std::stringstream sc;
      sc << Emod_filename << "_"<< eval << "_chirality_t" << t;   
      
      if(comm->primaryNode()) {
	outf.open(ss.str().c_str()); 
	outc.open(sc.str().c_str());
      }
      
      outf << "x  y  z  norm\n";
      outc << "x  y  z  chi_norm\n";
      
      for (int x = 0 ; x < Nx; x++){
	for (int y = 0 ; y < Ny; y++){
	  for (int z = 0 ; z < Nz; z++){
	    
	    
	    int site = SiteIndex::instance()->site(x,y,z,t);
	    int site0 = SiteIndex::instance()->site(x,y,z,0);
	    //Local norm;
	    double site_norm = 0.0;
	    SUNmat PL = mat(PL_field,site0);
	    double pl_re = SUNmatUtils::ReTr(PL);
	    double pl_im = SUNmatUtils::ImTr(PL);
	    
	    std::valarray<double> site_Array(EigM.evecs_[eval][FF.islice(site)]);
	    
	    //std::cout << "spinor\n";
	    //print_va(site_Array);
	    
	    // DIRAC representation
	    std::valarray<double> spinor_up(12);
	    std::valarray<double> spinor_dn(12);
	    
	    std::valarray<double> chiral(24);
	    std::valarray<double> chiral_left(24);
	    std::valarray<double> chiral_right(24);
	    
	    for (int s=0; s< 6; s++){
	      spinor_up[2*s]   = site_Array[2*s];
	      spinor_up[2*s+1] = site_Array[2*s+1];
	      
	      spinor_dn[2*s]   = site_Array[2*(s+ND_*NC_/2)];
	      spinor_dn[2*s+1] = site_Array[2*(s+ND_*NC_/2)+1];
	      
	      chiral[2*s]   = spinor_dn[2*s];
	      chiral[2*s+1] = spinor_dn[2*s+1];
	      
	      chiral[(s+ND_*NC_/2)*2]   = spinor_up[2*s];
	      chiral[(s+ND_*NC_/2)*2+1] = spinor_up[2*s+1];
	      
	      chiral_left[2*s]   = 0.5*(spinor_up[2*s] + spinor_dn[2*s]);
	      chiral_left[2*s+1] = 0.5*(spinor_up[2*s+1] + spinor_dn[2*s+1]);
	      chiral_left[(s+ND_*NC_/2)*2]   = 0.5*(spinor_up[2*s] + spinor_dn[2*s]);
	      chiral_left[(s+ND_*NC_/2)*2+1] = 0.5*(spinor_up[2*s+1] + spinor_dn[2*s+1]);
	      
	      chiral_right[2*s]   = 0.5*(spinor_up[2*s] - spinor_dn[2*s]);
	      chiral_right[2*s+1] = 0.5*(spinor_up[2*s+1] - spinor_dn[2*s+1]);
	      chiral_right[(s+ND_*NC_/2)*2]   = 0.5*(-spinor_up[2*s] + spinor_dn[2*s]);
	      chiral_right[(s+ND_*NC_/2)*2+1] = 0.5*(-spinor_up[2*s+1] + spinor_dn[2*s+1]);
	      
	    }
	    
	    //std::cout << "gamma5 spinor\n";
	    //print_va(chiral);
	    
	    double chi_norm = (site_Array*chiral).sum();
	    double chil_norm = (site_Array*chiral_left).sum();
	    double chir_norm = (site_Array*chiral_right).sum();
	    
	    //std::cout << "chi_norm: "<<chi_norm<<"\n";
	    
	    site_Array *=site_Array;
	    site_norm = sqrt(site_Array.sum());
	    
	    outf << SiteIndex::instance()->global_x(x) << " "
		 << SiteIndex::instance()->global_y(y) << " "
		 << SiteIndex::instance()->global_z(z) << " "
		 << site_norm << " "<< pl_re << " "<< pl_im <<"\n";
	    
	    ipr += site_norm*site_norm*site_norm*site_norm;
	    
	  
	    
	    
	    outc << SiteIndex::instance()->global_x(x) << " "
		 << SiteIndex::instance()->global_y(y) << " "
		 << SiteIndex::instance()->global_z(z) << " "
		 << chi_norm << " " << chil_norm << " " <<chir_norm 
		 << " " << (4/M_PI*atan(sqrt(chil_norm/chir_norm))-1.0) <<"\n";
	    
	    if (chi_norm > max) {
	      max = chi_norm;
	      max_loc[0] = x;	    max_loc[1] = y;
	      max_loc[2] = z;	    max_loc[3] = t;
	    }
	    
	    if (chi_norm < min) {
	      min = chi_norm;
	      min_loc[0] = x;	    min_loc[1] = y;
	      min_loc[2] = z;	    min_loc[3] = t;
	    }

	    
	    
	  }
	}
      }
      
      
      comm->sync();    
      
      
    
      if(comm->primaryNode()) {
	outf.close();
	outc.close();
      }
  }
    
    
    
    std::cout << "IPR ["<< eval << "]: "<<ipr*Nvol << "\n";
    
    
    
    std::cout << "Max at "<<max_loc[0]<< " "<<max_loc[1] << " "<<max_loc[2] << " "<<max_loc[3]<<"  "<< max<<"\n";
    std::cout << "Min at "<<min_loc[0]<< " "<<min_loc[1] << " "<<min_loc[2] << " "<<min_loc[3]<<"  "<< min<<"\n";
    
  }



}
