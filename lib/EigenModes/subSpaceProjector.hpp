/*!@file subSpaceProjector.hpp
 * @brief declaration of the utility to do projection of the fermion field 
 *  acorrding to given subspace.
 * Time-stamp: <2014-10-05 11:10:46 noaki>
 */
#ifndef SUBSPACEPROJECTOR_INCLUDED
#define SUBSPACEPROJECTOR_INCLUDED

#include "eigenModes.hpp"

namespace SubSpace{
  void project_real(Field& w,const Field& f,
		    const vector_Field& vsub,const std::vector<double>& q);
  void project_complex(Field& w,const Field& f,
		       const vector_Field& vsub,const std::vector<double>& q);

  void project(Field& w,const Field& f,const vector_Field& vsub);
  void projectOut(Field& w,const Field& f,const vector_Field& vsub);
}

#endif
