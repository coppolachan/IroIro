/*!
 * @file test_CorrelateEigenmodes.cpp
 * @brief Calculates the inner product of two eigenmodes, even from different files.
 *
 * @author Guido Cossu
 * @author <a href="http://suchix.kek.jp/guido_cossu/">Guido Cossu</a>
 *
 * Time-stamp: <2015-03-16 18:04:44 neo>
 */

#include "test_CorrelateEigenModes.hpp"
#include "include/format_F.h"
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

void print_va(std::valarray<double> va){
  for(int i=0; i<va.size(); ++i) std::cout<< i<<"  "<< va[i]<<"\n";
}



int Test_CorrelateEigenModes::run(){
  CCIO::cout<<"Test_ReadEigenModes::run() called\n";

  Communicator* comm = Communicator::instance();
  CommonPrms* cmprms = CommonPrms::instance();

  XML::node ReadEig_node = input_.node;
  XML::descend(ReadEig_node,"ReadEigenMode",MANDATORY);

  // Check eigenvector file to read
  std::string Emod_filename_1;
  XML::read(ReadEig_node, "ModesFilename_1", Emod_filename_1,MANDATORY);
  std::string Eval_filename_1;
  XML::read(ReadEig_node, "EvalFilename_1", Eval_filename_1,MANDATORY);

  std::string Emod_filename_2;
  XML::read(ReadEig_node, "ModesFilename_2", Emod_filename_2,MANDATORY);
  std::string Eval_filename_2;
  XML::read(ReadEig_node, "EvalFilename_2", Eval_filename_2,MANDATORY);

  
  int Nev_1;     // Which eigenvectors on first file
  int Nev_2;     // Which eigenvectors on second file

  XML::read(ReadEig_node,"NumEigenMode1",Nev_1, MANDATORY);
  XML::read(ReadEig_node,"NumEigenMode2",Nev_2, MANDATORY);
  int EigMode; // number of eigenmode to store
  XML::read(ReadEig_node,"EigenMode",EigMode, MANDATORY);

  //// Initialization and checks ////

  // Lattice size parameters
  int Nvol = cmprms->Nvol();
  int Nt = cmprms->Nt();
  int Nx = cmprms->Nx();
  int Ny = cmprms->Ny();
  int Nz = cmprms->Nz();
 
  // Allocate the required eigenvector fields 
  EigenModes EigM;
  Eigen::Predic* eig_pred = Eigen::predFactory(ReadEig_node.child("ReadCondition"));

  // Get number of stored eigenvalues
  std::ifstream reader1(Eval_filename_1.c_str()); 
  int idummy, iev=0;
  double ev=0.0;
  while(reader1>>idummy && reader1>>ev) iev++;
  CCIO::cout<<"Found "<< iev <<" eigenvalues in "<< Eval_filename_1<<"\n";

  std::ifstream reader2(Eval_filename_2.c_str()); 
  idummy, iev=0;
  ev=0.0;
  while(reader2>>idummy && reader2>>ev) iev++;
  CCIO::cout<<"Found "<< iev <<" eigenvalues in "<< Eval_filename_2<<"\n";

  Format::Format_F FF(Nvol);
  // Read eigenmodes from the two files
  double eval_1;
  Field evec_1(FF.size());
  double eval_2;
  Field evec_2(FF.size());

  Eigen::pickUpFromFile<Format::Format_F>(eval_1,evec_1,Nev_1,eig_pred,
					  Eval_filename_1,Emod_filename_1);
  
  CCIO::cout<< "Eigenvalue 1 : "<< eval_1<<"\n";
  CCIO::cout<< "Norm 1 : "<< evec_1.norm()<<"\n";
  
  Eigen::pickUpFromFile<Format::Format_F>(eval_2,evec_2,Nev_2,eig_pred,
					  Eval_filename_2,Emod_filename_2);

  CCIO::cout<< "Eigenvalue 2 : "<< eval_2<<"\n";
  CCIO::cout<< "Norm 2 : "<< evec_2.norm()<<"\n";
  
  
  
  std::ofstream outf, outc;
  
  
  double site_innerprod = 0.0;   

 
  
  for(int t=0; t<Nt; ++t){
    for(int x=0; x<Nx; ++x){
      for(int y=0; y<Ny; ++y){
	for(int z=0; z<Nz; ++z){
	  
	  int site = SiteIndex::instance()->site(x,y,z,t);
	  int site0 = SiteIndex::instance()->site(x,y,z,0);
	  //Local norm;

	  
	  std::valarray<double> site_Array_1(evec_1[FF.islice(site)]);
	  std::valarray<double> site_Array_2(evec_2[FF.islice(site)]);
	    
	  //std::cout << "spinor\n";
	  //print_va(site_Array);
	  
	  site_Array_2 *=site_Array_1;  // overwrites array_2
	  site_innerprod += site_Array_2.sum();
	  
	  
	  
	  
	}//Z
      }//Y
    }//X
  }//T
 

  CCIO::cout << "Inner product :  " << Nev_1 << "  "<< Nev_2 << "  " << eval_1 << "  " << eval_2 << " " << site_innerprod << "\n";
}


