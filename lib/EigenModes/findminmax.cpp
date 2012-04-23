/*! 
  @file findminmax.cpp

  @brief Find the minimum and maximum eigenmode of a real matrix

*/

#include "findminmax.hpp"
#include "include/messages_macros.hpp"

findMinMax::findMinMax(Fopr* K, 
		       const RandNum* R,
		       size_t KernelVectSize):
  Kernel_(K),
  rand_(R),
  vect_size(KernelVectSize){}


MinMaxOut findMinMax::findExtrema() const {
  MinMaxOut res;
  res.max = findMax();
  CCIO::cout << "Maximum eigenvalue : " << res.max << "\n";
  res.min = findMin(res.max);
  CCIO::cout << "Minimum eigenvalue : " << res.min << "\n";
  return res;
}

double findMinMax::findMax() const {
  std::valarray<double> gauss_va;
  int loop_count;
  double norm, epsilon, old_norm;
  double target_prec = 1e-5;

  gauss_va.resize(vect_size);
  rand_->get_gauss(gauss_va); 

  Field gauss_vect(gauss_va);
  old_norm = gauss_vect.norm();
  gauss_vect /= old_norm;
  old_norm = 1.0;
  loop_count = 0;

  CCIO::cout << "Calculating maximum eigenmode...";
  do {
    _Message(DEBUG_VERB_LEVEL,"Gauss vect old_norm = "<< old_norm << "\n");

    gauss_vect = Kernel_->mult(gauss_vect);

    norm = gauss_vect.norm();
    _Message(DEBUG_VERB_LEVEL,"Gauss vect norm = "<< norm << "\n");
    gauss_vect /= norm;
    
    epsilon = fabs(old_norm-norm)/norm;
    old_norm = norm;    
    loop_count++;
    _Message(DEBUG_VERB_LEVEL,"Epsilon = "<< epsilon << "\n");
    _Message(DEBUG_VERB_LEVEL,"Loop #"<<loop_count << "\n");

  } while (epsilon > target_prec);

  CCIO::cout << "done\n";
  return norm;
}

double findMinMax::findMin(const double& max) const {
  std::valarray<double> gauss_va;
  int loop_count;
  double norm, epsilon, old_norm;
  double target_prec = 1e-5;

  gauss_va.resize(vect_size);
  rand_->get_gauss(gauss_va); 

  Field gauss_vect(gauss_va), temp_vect;
  temp_vect.resize(gauss_vect.size());
  old_norm = gauss_vect.norm();
  gauss_vect /= old_norm;
  old_norm = 1.0;
  loop_count = 0;

  CCIO::cout << "Calculating minimum eigenmode...";
  do {
    _Message(DEBUG_VERB_LEVEL,"Gauss vect old_norm = "<< old_norm << "\n");

    temp_vect = gauss_vect;
    temp_vect *= max;
    temp_vect -= Kernel_->mult(gauss_vect);
    gauss_vect = temp_vect;

    norm = gauss_vect.norm();
    _Message(DEBUG_VERB_LEVEL,"Gauss vect norm = "<< norm << "\n");
    gauss_vect /= norm;
    
    epsilon = fabs(old_norm-norm)/norm;
    old_norm = norm;    
    loop_count++;
    _Message(DEBUG_VERB_LEVEL,"Epsilon = "<< epsilon << "\n");
    _Message(DEBUG_VERB_LEVEL,"Loop #"<<loop_count << "\n");

  } while (epsilon > target_prec);

  CCIO::cout << "done\n";
  return (max-norm);
}
