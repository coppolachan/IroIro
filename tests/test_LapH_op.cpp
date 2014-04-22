/* 
 * 
 * @file test_LapH_op.cpp
 * @brief Tests for the LapH meson operator construction
 */
#include <stdio.h>

#include "test_LapH.hpp"
#include "Dirac_ops/dirac_Operator_FactoryCreator.hpp"
#include "Measurements/GaugeM/staples.hpp"

#include "include/timings.hpp"
#include "include/messages_macros.hpp"

using namespace std;
using namespace Format;
using namespace EvenOddUtils;

typedef FermionField1sp ScalarField;

int LapH_Pion() {

  int Nvol = CommonPrms::instance()->Nvol();
  int Nt = CommonPrms::instance()->Nt();
  int Nx = CommonPrms::instance()->Nx();
  int Ny = CommonPrms::instance()->Ny();
  int Nz = CommonPrms::instance()->Nz();
  int Nvol3D = Nx*Ny*Nz;


  char filename[128];
  
  int Ndilution = 6;
  int Nev = 120;
  int time_slice = 0;


  std::vector<double> source_meson_operator;

  std::vector<FermionField> srcfield_0(ND_);
  std::vector<FermionField> srcfield_1(ND_);

  std::vector<FermionField> phi_0(ND_);
  std::vector<FermionField> smeared_phi_0(ND_);
  std::vector<FermionField> phi_1(ND_);
  std::vector<FermionField> smeared_phi_1(ND_);
  std::vector<ScalarField> Sfield(Nev); 
  for (int i = 0; i < Nev; i++){
    Sfield[i].resize(Nvol3D*2);//resize 
  }


  // !!!!!!!! WARNING:  assuming chiral representation in the contraction


  // Load all the Laplacian eigenvectors to smear the sinks
  for(int eigvec_number=0; eigvec_number < Nev; eigvec_number++)
    {     
      CCIO::ReadFromDisk3D< Format::Format_S >(Sfield[eigvec_number].data, 
					       "eigen_t0_5830", 
					       eigvec_number, 0, "Binary", false);
      
    }
  
  
  
  for (int i = 0 ; i < Ndilution; i++){
    
    for (int spin = 0; spin < ND_ ; spin++){
      
      // Load [0] source 
      sprintf(filename, "smeared_rho_0_NvI6_t0_i%d", i*ND_+spin);
      CCIO::ReadFromDisk<FermionField>(srcfield_0[spin], filename);
      sprintf(filename, "phi_0_NvI6_t0_i%d", i*ND_+spin);
      CCIO::ReadFromDisk<FermionField>(phi_0[spin], filename);
      
      
      //Smear the sink fields moe to external routine
      std::vector<double> tempv(Nev*2);
      for (int ev= 0; ev < Nev; ev+=2){
	for (int site = 0 ; site  < Nvol3D; site++) {
	  tempv[ev] += Sfield[ev][site*2]*phi_0[spin][(Nvol3D*time_slice+site)*NC_*ND_*2]+
	    Sfield[ev][site*2+1]*phi_0[spin][(Nvol3D*time_slice+site)*NC_*ND_*2+1];
	  tempv[ev+1] += Sfield[ev][site*2]*phi_0[spin][(Nvol3D*time_slice+site)*NC_*ND_*2+1]-
	    Sfield[ev][site*2+1]*phi_0[spin][(Nvol3D*time_slice+site)*NC_*ND_*2];
	}
      }

      for (int site = 0 ; site  < Nvol3D; site++) {    
	for (int ev= 0; ev < Nev; ev++){
	  smeared_phi_0[spin].add((Nvol3D*time_slice+site)*NC_*ND_*2, Sfield[ev][site*2]*tempv[ev] - Sfield[ev][site*2+1]*tempv[ev+1]);
	  smeared_phi_0[spin].add((Nvol3D*time_slice+site)*NC_*ND_*2+1, Sfield[ev][site*2]*tempv[ev+1] + Sfield[ev][site*2+1]*tempv[ev]);
	}
      }
      // got the smeared sink 


    }//ndilution 
    
    for (int j = 0 ; j < Ndilution; j++){

      for (int spin = 0; spin < ND_ ; spin++){
	// Load [1] source 
	sprintf(filename, "smeared_rho_1_NvI6_t0_i%d", j*ND_+spin);
	CCIO::ReadFromDisk<FermionField>(srcfield_1[spin], filename);
	sprintf(filename, "phi_1_NvI6_t0_i%d", i*ND_+spin);
	CCIO::ReadFromDisk<FermionField>(phi_1[spin], filename);

	//smear the sink phi_1 with the same routine as before


      }

      // Contraction
      double ij_element = 0.0;
      for (int spin = 0; spin < ND_ ; spin++){
	ij_element = srcfield_0[spin]*srcfield_1[spin]; //pion
      }
      source_meson_operator.append(ij_element);//indexed 1-d vector i*(Ndilution)+j
      

    }
  }


}
