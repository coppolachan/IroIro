/*! 
  @file findminmax.cpp

  @brief Find the minimum and maximum eigenmode of a real matrix

*/

#include "findminmax.hpp"
#include "Tools/randNum_Factory.h"

findMinMax::findMinMax(Fopr* K, 
		       const RandNum* R,
		       size_t KernelVectSize):
  Kernel_(K),
  rand_(R),
  vect_size(KernelVectSize){}

double findMinMax::findMax() const {
  std::valarray<double> gauss_va;
  int loop_count;
  double norm, epsilon, old_norm;
  double target_prec = 1e-5;
  double max;

  gauss_va.resize(vect_size);
  CCIO::cout << "Getting gauss... ";
  rand_->get_gauss(gauss_va); 
  CCIO::cout << "ok\n";

  Field gauss_vect(gauss_va);

  loop_count = 0;
  do {
    old_norm = gauss_vect.norm();
    CCIO::cout << "Gauss vect old_norm = "<< old_norm << "\n";
    gauss_vect /= old_norm;

    gauss_vect = Kernel_->mult(gauss_vect);

    norm = gauss_vect.norm();
    CCIO::cout << "Gauss vect norm = "<< norm << "\n";
 
    
    epsilon = fabs(old_norm-norm)/norm;
    
    loop_count++;
    CCIO::cout << "Epsilon = "<< epsilon << "\n";
    CCIO::cout << "Loop #"<<loop_count << "\n";
  } while (epsilon > target_prec);

  return max = norm/old_norm;
  


}
