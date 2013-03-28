//---------------------------------------------------------------------
/*! @file sunRepresentations.hpp
  @brief \f$SU(N)\f$ Generators in various representations

  Class declarations
*/ 
//---------------------------------------------------------------------
#ifndef SUNREP_INCLUDED
#define SUNREP_INCLUDED

#include <iostream>
#include "include/macros.hpp"
#include "sunMat.hpp"


#define ADJCOL COLORS*COLORS-1

template <size_t COLORS>
class SUNRep {
  SUNmatrix<ADJCOL> AdjointGen[ADJCOL];
  SUNmatrix<COLORS> FundamentalGen[ADJCOL];
  
  void GenerateRep();
public:
  SUNRep() {
    GenerateRep();
  }

  SUNmatrix<ADJCOL> lambda_adj(int a){ return AdjointGen[a]; }
  SUNmatrix<COLORS> lambda_fund(int a){ return FundamentalGen[a]; }

};

template <size_t COLORS>
void SUNRep<COLORS>::GenerateRep() {
  // Generating the fundamental representation
  ////////////////////////////////////////////
  double rsq = 1.0/sqrt(2.0);
  // initialization
  _Message(0, "SUNRep Message: Generating fundamental representation for SU("<<COLORS<<")...\n");
  for (int a = 0; a < ADJCOL; a++) 
    FundamentalGen[a].zero();
  
  // setting off-diagonal generators
  _Message(DEBUG_VERB_LEVEL, "Off diagonal generators...\n");
  int a = 0;
  for (int i = 0; i < COLORS-1; i++) {
    for (int j = i+1; j < COLORS; j++) {
      FundamentalGen[a].set(i,j, rsq, 0.0);
      FundamentalGen[a].set(j,i, rsq, 0.0); 

      a++;

      FundamentalGen[a].set(i,j, 0.0, -rsq);
      FundamentalGen[a].set(j,i, 0.0,  rsq);
      a++;	
    }
  }
 _Message(DEBUG_VERB_LEVEL, "Wrote "<<a<<" off-diagonal generators\n");

  for (int j = 0; j < COLORS-1; j++) {
    for (int i = 0; i < j; i++) {
      FundamentalGen[a].set(i,i, 1.0/sqrt((double)((j+1)*(j+2))), 0.0);
    }
    FundamentalGen[a].set(j+2,j+2, -double(j+1)/sqrt((double)((j+1)*(j+2))), 0.0);
    a++;
  }

  _Message(DEBUG_VERB_LEVEL, "Wrote "<<a<<" fundamental generators\n");

  // Generating the adjoint representation
  ////////////////////////////////////////////
  _Message(0, "SUNRep Message: Generating adjoint representation for SU("<<COLORS<<")...\n");
  SUNmatrix<ADJCOL> temp;
  for(int ia = 0; ia < ADJCOL; ia++) {

    for (int ic = 0; ic < ADJCOL; ic++) {
      for (int ib = 0; ib < ADJCOL; ib++) {
	
	temp.zero();
	for (int j = 0; j < COLORS; j++){ 
	  for (int i = 0; i < COLORS; i++) {
	    for (int k = 0; k < COLORS; k++) {
	      temp.set(i,j, 
		       temp.r(i,j)+
		       (FundamentalGen[ia].r(i,k)*FundamentalGen[ib].i(k,j)+
			FundamentalGen[ia].i(i,k)*FundamentalGen[ib].r(k,j)-
			FundamentalGen[ib].r(i,k)*FundamentalGen[ia].i(k,j)-
			FundamentalGen[ib].i(i,k)*FundamentalGen[ia].r(k,j)),
			temp.i(i,j) + 
		       (FundamentalGen[ia].i(i,k)*FundamentalGen[ib].i(k,j)-
			FundamentalGen[ia].r(i,k)*FundamentalGen[ib].r(k,j)+
			FundamentalGen[ib].r(i,k)*FundamentalGen[ia].r(k,j)-
			FundamentalGen[ib].i(i,k)*FundamentalGen[ia].i(k,j)));
	    }
	  }
	}

	for (int i = 0; i < COLORS; i++) {
	  for (int k = 0; k < COLORS; k++) {
	    AdjointGen[ia].set(ib,ic,
			       AdjointGen[ia].r(ib,ic) +
			       rsq*(temp.r(i,k)*FundamentalGen[ic].r(k,i)-temp.i(i,k)*FundamentalGen[ic].i(k,i)),
			       AdjointGen[ia].i(ib,ic) +
			       rsq*(temp.r(i,k)*FundamentalGen[ic].i(k,i)+temp.i(i,k)*FundamentalGen[ic].r(k,i)));
	  }
	}


      }
    }
    _Message(DEBUG_VERB_LEVEL, "Generated adjoint generator number "<<ia<<"\n");
  }






}

#endif
