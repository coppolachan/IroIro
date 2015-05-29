/*!
 * @file test_EigenModesPropagator.cpp
 * @brief Reads eigenmodes and calculated the meson correlators
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2015-02-02 18:29:55 neo>
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
#include <omp.h>

////////////////////////////////////////
// IROIRO Uses DIRAC representation
////////////////////////////////////////

using namespace FieldUtils;
using namespace SUNmatUtils;

void print_va(std::valarray<double> va){
  for(int i=0; i<va.size(); ++i) std::cout<< i<<"  "<< va[i]<<"\n";
}

double plaquette_site(const GaugeField& F,Staples stp,int Nvol_){
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

  for(int t=0; t<Nt; ++t){
    std::stringstream ss;
    ss<< "ScalarGball_sites_t" << t; 
    plaqf.open(ss.str().c_str()); 
    
    plaqf <<"x  y  z  norm\n";
    
    for(int x = 0; x < Nx; ++x){
      for(int y = 0; y < Ny; ++y){
	for(int z = 0; z < Nz; ++z){
	  int site = SiteIndex::instance()->site(x,y,z,t);
	  plaqf<< x <<" "<< y <<" "<< z <<" "<< plaqv[site]<<"\n";    
	}
      }
    }
    plaqf.close();
  }
  std::cout<<"Plaq Min at "<<min_loc[0]<<" "<<min_loc[1]<<" "<<min_loc[2]<<" "<<min_loc[3]
	   <<"  "<< min_plaq<<"\n";  
}

int Test_ReadEigenModes::run(){
  CCIO::cout<<"Test_ReadEigenModes::run() called\n";
  bool Gamma5 = 0;


  Communicator* comm = Communicator::instance();
  CommonPrms* cmprms = CommonPrms::instance();
  SiteIndex* SIdx    = SiteIndex::instance();
  XML::node ReadEig_node = input_.node;
  XML::descend(ReadEig_node,"ReadEigenMode",MANDATORY);

  // Mass parameter
  double mass; // number of eigenmode to store
  XML::read(ReadEig_node,"Mass",mass, MANDATORY);

  // Check eigenvector file to read
  std::string Emod_filename;
  XML::read(ReadEig_node, "ModesFilename", Emod_filename,MANDATORY);
  std::string Eval_filename;
  XML::read(ReadEig_node, "EvalFilename", Eval_filename,MANDATORY);
  
  int Nev;     // Number of eigenvectors to read
  XML::read(ReadEig_node,"NumEigenModes",Nev, MANDATORY);
  int EigMode; // number of eigenmode to store
  XML::read(ReadEig_node,"EigenMode",EigMode, MANDATORY);
  XML::read(ReadEig_node,"OnlyGamma5", Gamma5);


  CCIO::cout << "Calculating only gamma_5 : "<< Gamma5 << "\n";
  // Prefactor
  double prefactor = 1.0/((1.0-mass)*(1.0-mass));


  //// Initialization and checks ////

  // Lattice size parameters
  int Nvol = cmprms->Nvol();
  int Nt = cmprms->Nt();
  int Nx = cmprms->Nx();
  int Ny = cmprms->Ny();
  int Nz = cmprms->Nz();

  int Lt = cmprms->Lt();
  int Lx = cmprms->Lx();
  int Ly = cmprms->Ly();
  int Lz = cmprms->Lz();
  Staples stpl;
  double plq = stpl.plaquette(*(input_.config.gconf));
  CCIO::cout<<" Plaquette ="<< plq <<std::endl;
 
  Format::Format_F FF(Nvol);

  // Allocate the required eigenvector fields 
  EigenModes EigM;
  

  Eigen::Predic* eig_pred = Eigen::predFactory(ReadEig_node.child("ReadCondition"));
  // Read the required eigenvectors and store in fields
  Eigen::initFromFile<Format::Format_F>(EigM,eig_pred,Eval_filename, Emod_filename);
  
  int iev = EigM.evals.size();
  CCIO::cout<<"Found "<< iev <<" eigenvalues in "<< Eval_filename<<"\n";
  
  
  FILE* fhandle = fopen("correlators.data", "w"); 

  // Finite temperature meson spatial correlators
  double pion[Nx];
  double pion_im[Nx];
  double delta[Nx];
  double delta_im[Nx];

  double* block_ps = (double*)calloc(Nvol*2, sizeof(double));
  double* block_s  = (double*)calloc(Nvol*2, sizeof(double));
  
  for (int dist = 1; dist < Nx ; dist++){
    pion[dist]     = 0.0;
    pion_im[dist]  = 0.0;
    delta[dist]    = 0.0;
    delta_im[dist] = 0.0;
  }
  
  int source_num = 4;
  int skip_modes = 0;


  // DIRAC representation
  std::valarray<double> spinor_up(12);
  std::valarray<double> spinor_dn(12);

  
  // Cycle over the eigenmodes
  for(int eid1=skip_modes; eid1<std::min(Nev,iev); eid1++){
    double eval1 = (EigM.evals[eid1]);
    double chirality = 0.0;
    double delta_gamma5 = 0.0;
    for(int eid2=skip_modes; eid2<std::min(Nev,iev); eid2++){
      double eval2 = (EigM.evals[eid2]);
      
      double check_r = 0.0;
      double check_i = 0.0;
      double chir    = 0.0;
      double chii    = 0.0;

      for (int site = 0; site < Nvol ; site++){
	block_ps[2*site  ] = 0.0;
	block_ps[2*site+1] = 0.0;

	block_s[2*site  ]  = 0.0;
	block_s[2*site+1]  = 0.0;
	std::valarray<double> nX(EigM.evecs[eid1][FF.islice(site)]);
	std::valarray<double> mX(EigM.evecs[eid2][FF.islice(site)]);

	// nX^dag * mX = \psi_n(x)^dag \psi_m(x)
	for (int s=0; s< 6; s++){
	  int re = 2*s; int im = 2*s + 1;
	  
	  spinor_up[re]  = mX[2*s];
	  spinor_up[im]  = mX[2*s+1];
	  
	  spinor_dn[re]  = mX[2*s+ND_*NC_];
	  spinor_dn[im]  = mX[2*s+ND_*NC_+1];

	  //pseudoscalar ok
	  block_ps[2*site  ] +=  nX[re]*mX[re] + nX[im]*mX[im];  //real 
	  block_ps[2*site+1] += -nX[re]*mX[im] + nX[im]*mX[re];  //imag
	  block_ps[2*site  ] +=  nX[re+ND_*NC_]*mX[re+ND_*NC_] + nX[im+ND_*NC_]*mX[im+ND_*NC_];  //real 
	  block_ps[2*site+1] += -nX[re+ND_*NC_]*mX[im+ND_*NC_] + nX[im+ND_*NC_]*mX[re+ND_*NC_];  //imag

	  //scalar  -  matrix elements of gamma_5
	  block_s[2*site  ] +=  nX[re]*spinor_dn[re] + nX[im]*spinor_dn[im];  //real 
	  block_s[2*site+1] += -nX[re]*spinor_dn[im] + nX[im]*spinor_dn[re];  //imag
	  block_s[2*site  ] +=  nX[re+ND_*NC_]*spinor_up[re] + nX[im+ND_*NC_]*spinor_up[im];  //real 
	  block_s[2*site+1] += -nX[re+ND_*NC_]*spinor_up[im] + nX[im+ND_*NC_]*spinor_up[re];  //imag
	}
	
	check_r += block_ps[2*site];
	check_i += block_ps[2*site+1];
	
	chir += block_s[2*site];
	chii += block_s[2*site+1];

      }
      // check orthogonality sum(block_ps) =0 if n neq m
      printf("Gamma5 %5d\t%5d\t%25.15e\t%25.15e\t%25.15e\t%25.15e\n", eid1, eid2, check_r, check_i, chir, chii);

      delta_gamma5 += (1.0+mass)*(eval1+eval2)*(chir*chir+chii*chii); 

      if (eid1 == eid2)
	chirality = chir;

      
      if (!Gamma5){
	
	// Correlator
	
	//source
	
	int x1 = 0;
	int y1 = 0;
	int z1 = 0;
	int t1 = 0;
	
	
	for (int x1 = 0; x1  < Nx ; x1+=Nx/source_num){
	  y1= x1; z1=x1;
	  
	  for (int x  = 0; x  < Nx ; x++){
	    //sources
	    int sx = SIdx->site(x,y1,z1,t1);
	    int sy = SIdx->site(x1,x,z1,t1);
	    int sz = SIdx->site(x1,y1,x,t1);
	    
	    //sink
	    for (int dist = 1; dist < Nx ; dist++){
	      for (int sink = 0; sink < SIdx->slsize(0,0); sink++){
		double psr, psi;
		double sr, si;
		int s2 = (SIdx->*SiteIndex::slice_dir[0])((x+dist)%Nx,sink);//sink1
		
		psr  = block_ps[2*sx  ]*block_ps[2*s2] + block_ps[2*sx+1]*block_ps[2*s2+1];
		psi  = block_ps[2*sx+1]*block_ps[2*s2] - block_ps[2*sx  ]*block_ps[2*s2+1];
		
		sr  = block_s[2*sx  ]*block_s[2*s2] + block_s[2*sx+1]*block_s[2*s2+1];
		si  = block_s[2*sx+1]*block_s[2*s2] - block_s[2*sx  ]*block_s[2*s2+1];
		
		s2 = (SIdx->*SiteIndex::slice_dir[1])((x+dist)%Ny,sink);//sink2
		psr += block_ps[2*sy  ]*block_ps[2*s2] + block_ps[2*sy+1]*block_ps[2*s2+1];
		psi += block_ps[2*sy+1]*block_ps[2*s2] - block_ps[2*sy  ]*block_ps[2*s2+1];
		
		sr  += block_s[2*sy  ]*block_s[2*s2] + block_s[2*sy+1]*block_s[2*s2+1];
		si  += block_s[2*sy+1]*block_s[2*s2] - block_s[2*sy  ]*block_s[2*s2+1];
		
		s2 = (SIdx->*SiteIndex::slice_dir[2])((x+dist)%Nz,sink);//sink3
		psr += block_ps[2*sz  ]*block_ps[2*s2] + block_ps[2*sz+1]*block_ps[2*s2+1];
		psi += block_ps[2*sz+1]*block_ps[2*s2] - block_ps[2*sz  ]*block_ps[2*s2+1];
		
		sr  += block_s[2*sz  ]*block_s[2*s2] + block_s[2*sz+1]*block_s[2*s2+1];
		si  += block_s[2*sz+1]*block_s[2*s2] - block_s[2*sz  ]*block_s[2*s2+1];
		
		pion[dist]    += psr/(eval1*eval2);
		pion_im[dist] += psi/(eval1*eval2);
		
		delta[dist]    += sr/(eval1*eval2);
		delta_im[dist] += si/(eval1*eval2);
	      }//end sink
	      
	    }//end dist
	  }//end source xyz
	}//end x1
      }
    }//end eid2

    double delta_n = -2.0*(eval1*eval1+mass) + chirality*(1.0+mass)*eval1*2.0 ;
    delta_gamma5 -= 2.0*chirality*(eval1*eval1+mass);
    printf("Delta\t%5d\t%25.15e\t%25.15e\n", eid1, delta_n, delta_gamma5);

    double mres_numerator   = (1.0+eval1*eval1)/(eval1*eval1)*delta_n - 2.0/eval1*delta_gamma5;
    double mres_denominator = (1.0+eval1*eval1)/(eval1*eval1) - 2.0*chirality/eval1;

   printf("Mres\t%5d\t%25.15e\t%25.15e\n", eid1, mres_numerator, mres_denominator);


    if (!Gamma5){
      for (int s=1; s<Nx; ++s){
	printf("%5d\t%5d\t%20.15e\t%20.15e\n", eid1, s, prefactor*pion[s]/(source_num*3*Nx), prefactor*pion_im[s]/(source_num*3*Nx)); // normalize by the number of sources
      }
      for (int s=1; s<Nx; ++s){
	printf(" - %5d\t%5d\t%20.15e\t%20.15e\n", eid1, s, prefactor*delta[s]/(source_num*3*Nx), prefactor*delta_im[s]/(source_num*3*Nx)); // normalize by the number of sources
      }
      
    }
    fflush(stdout);
  }//end eid1
  

  if (!Gamma5){
    fprintf(fhandle, "#### PseudoScalar\n");
    for (int s=1; s<Nx; ++s){
      fprintf(fhandle, "%5d\t%20.15e\t%20.15e\n", s, prefactor*pion[s]/(source_num*3*Nx), prefactor*pion_im[s]/(source_num*3*Nx)); // normalize by the number of sources
    }
    fprintf(fhandle, "#### Scalar\n");
    for (int s=1; s<Nx; ++s){
      fprintf(fhandle, "%5d\t%20.15e\t%20.15e\n", s, prefactor*delta[s]/(source_num*3*Nx), prefactor*delta_im[s]/(source_num*3*Nx)); // normalize by the number of sources
    }
    
  }
  
  fclose(fhandle);
  free(block_ps);
  free(block_s);
  
  
}
