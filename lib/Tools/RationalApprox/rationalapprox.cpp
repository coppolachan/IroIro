/*!

 * @file rationalapprox.cpp

 * @brief Calculates the rational approximation for a given function

 */

#include <assert.h>
#include "rationalapprox.hpp"
#include "Tools/Remez/alg_remez.hpp"
#include "Communicator/comm_io.hpp"

// Standard Constructor 
RationalApprox::RationalApprox(RationalApprox_params Par):Params(Par)
{
  fill();
}

// XML constructor
RationalApprox::RationalApprox(const XML::node Approx_node):
  Params(Approx_node)
{
  fill();
}

void RationalApprox::fill() {

  CCIO::cout << "Calculating Rational Approximation of the function:\n";
  CCIO::cout << "f(x) = x^("<<Params.exponent_num<<"/"<<Params.exponent_den<<") ";
  CCIO::cout << "in the interval ["<<Params.lambda_low<<";"<<Params.lambda_high<<"]\n";

  AlgRemez RemezApprox(Params.lambda_low, Params.lambda_high, Params.gmp_remez_precision);

  double error = RemezApprox.generateApprox(Params.numerator_deg,
					    Params.denominator_deg,
					    Params.exponent_num,
					    Params.exponent_den);
  
  CCIO::cout << "Approximation error: "<< error << "\n";

  // Find the partial fraction expansion of the approximation 
  // to the function x^{a/b} (this only works currently for 
  // the special case that num_deg = den_deg)
  
  // temporary arrays
  double *res = new double[Params.numerator_deg];
  double *den = new double[Params.denominator_deg];

  assert (Params.numerator_deg == Params.denominator_deg);
  RemezApprox.getPFE(res, den, &RA_a0);

  // Fill the vectors
  RA_res.resize(Params.numerator_deg);
  RA_pole.resize(Params.denominator_deg);


  for (int i = 0; i < Params.numerator_deg;  i++) {
    RA_res[i]  = res[i];
    RA_pole[i] = den[i];

    CCIO::cout << "Res["<<i<<"] = "<< std::setw(25)<< RA_res[i] 
	       << "   Pole["<<i<<"] = "<< std::setw(25)<< RA_pole[i] << "\n";
  }

  delete[] res;
  delete[] den;

}







































/*

                                             // RESCALED COEFFICIENTS

void RationalApprox::first_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::first_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;

 // normalized coefficients
 *this=*first_inv_approx_norm_coeff;

 findminmax(min, max);

 min*=0.95;
 max*=1.05;

 epsilon=min/max;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_inv_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = ((double) no_flavours) /8./ ((double)no_ps);
 epsilon=pow(max, exponent);
 RA_a0*=epsilon;
 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::first_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


void RationalApprox::md_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::md_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;

 // normalized coefficients
 *this=*md_inv_approx_norm_coeff;

 findminmax(min, max);

 min*=0.95;
 max*=1.05;

 epsilon=min/max;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_md_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = -((double) no_flavours)/4./((double)no_ps);
 epsilon=pow(max, exponent);
 RA_a0*=epsilon;
 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::md_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


void RationalApprox::last_inv_approx_coeff(void)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside RationalApprox::last_inv_approx_coef ..."<<endl;
 #endif
 int i;
 REAL min, max, epsilon;

 // normalized coefficients
 *this=*last_inv_approx_norm_coeff;

 findminmax(min, max);

 min*=0.95;
 max*=1.05;

 epsilon=min/max;

 if(epsilon<min_epsilon)
   {
   ofstream err_file;
   err_file.open(QUOTEME(ERROR_FILE), ios::app);   
   err_file << "WARNING: in first_last_approx_coef epsilon="<<epsilon<<" < min_epsilon="<<min_epsilon;
   err_file << " at update_iteration="<<update_iteration <<endl;
   err_file.close();
   } 

  
 // rescale coeff.
 double exponent = -((double) no_flavours)/4./((double)no_ps);
 epsilon=pow(max, exponent);
 RA_a0*=epsilon;
 for(i=0; i<approx_order; i++)
    {
    RA_a[i]*=(max*epsilon);
    RA_b[i]*=max;
    }

 #ifdef DEBUG_MODE
 cout << "\tterminated RationalApprox::last_inv_approx_coeff, epsilon="<< min/max<< " min_epsilon="<<min_epsilon<<endl;
 #endif
 }


                                     // EXPLICIT CALCULATIONS ON FERMIONS


void first_inv_approx_calc(REAL res)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside first_inv_approx_calc ..."<<endl;
 #endif

 RationalApprox approx;
 approx.first_inv_approx_coeff(); 

 #ifdef USE_GPU
 int order;
 double constant, numerators[max_approx_order];

 cu_multips_shifted_invert (res, approx);

 get_order(order, approx);
 get_const(constant, approx);
 get_numerators(numerators, approx);

 cuda_sol_sum(order, constant, numerators);

 smartunpack_multifermion(fermion_chi, chi_packed);
 #endif

 #ifndef USE_GPU
 int iter, pseudofermion;
 long int i;
 Vec3 vr_1;

 multips_shifted_invert (fermion_shiftmulti, fermion_phi, res, approx);

 for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++)
    {
    for(i=0; i<sizeh; i++)
       {
       vr_1=(approx.RA_a0)*(fermion_phi->fermion[pseudofermion][i]);
       for(iter=0; iter<(approx.approx_order); iter++)
          {
          vr_1+=(approx.RA_a[iter])*(fermion_shiftmulti->fermion[pseudofermion][iter][i]);
          }
       fermion_chi->fermion[pseudofermion][i]=vr_1;
       }
    }
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated first_inv_approx_calc"<<endl;
 #endif
 }



void last_inv_approx_calc(REAL res)
 {
 #ifdef DEBUG_MODE
 cout << "DEBUG: inside  last_inv_approx_calc  .."<<endl;
 #endif

 RationalApprox approx;
 approx.last_inv_approx_coeff(); 

 #ifdef USE_GPU
 int order;
 double constant, numerators[max_approx_order];

 cu_multips_shifted_invert (res, approx);

 get_order(order, approx);
 get_const(constant, approx);
 get_numerators(numerators, approx);

 cuda_sol_sum(order, constant, numerators);

 smartunpack_multifermion(fermion_phi, chi_packed);
 #endif

 #ifndef USE_GPU
 int iter, pseudofermion;
 long int i;
 Vec3 vr_1;

 multips_shifted_invert (fermion_shiftmulti, fermion_chi, res, approx);

 for(pseudofermion=0; pseudofermion<no_ps; pseudofermion++)
    {
    for(i=0; i<sizeh; i++)
       {
       vr_1=(approx.RA_a0)*(fermion_chi->fermion[pseudofermion][i]);
       for(iter=0; iter<(approx.approx_order); iter++)
          {
          vr_1+=(approx.RA_a[iter])*(fermion_shiftmulti->fermion[pseudofermion][iter][i]);
          }
       fermion_phi->fermion[pseudofermion][i]=vr_1;
       }
    }
 #endif

 #ifdef DEBUG_MODE
 cout << "\tterminated root_1_12_calc"<<endl;
 #endif
 }

*/
