/*!
        @file    $Id:: gaugeFixing_Coulomb.h #$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru) 
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2012-03-15 17:01:15 #$

        @version $LastChangedRevision: 137 $
*/

#ifndef GAUGEFIXING_COULOMB_INCLUDED
#define GAUGEFIXING_COULOMB_INCLUDED

#include <valarray>

#include "index_eo.h"
#include "field_G.h"

class RandomNumbers;

//! Coulomb gauge fixing.
/*
  This class fix the gauge of configuration to the Coulomb gauge.
  The implementation assumes that the dimension is 4 and the
  Coulomb gauge fixing is performed within each time slice.
  The algorithm is that developed by the Los Alamos group [see the
  implementation note].
  Overrelaxation is incorporated.
  To escape the Gribov copy, if convergence is not reached on some
  timeslices within the iterations specified by Nreset, random
  gauge transformation is performed to reset the configuration on
  that timeslice.
  This is the reason that random number generator is needed at the
  construction of this class.

  The implementation is not complete:
  - only applies to SU(3) case: because of specific implementation
      of maxTr function (Cabibbo-Marinari maximization).
  - unnecessary arithmetic operations exist for the timeslices
    on which the gauge is already fixed to good precision.
  These should be improved in the version beyond test phase.
                                        [16 Feb 2012 H.Matsufuru]
 */
class GaugeFixing_Coulomb{

 private:
  int    d_Niter;   // max iteration number
  int    d_Nnaive;  // number of naive iterations
  int    d_Nmeas;   // interval of measurements
  int    d_Nreset;  // Number of iteration to reset the config.
  double d_Enorm;   // convergence criterion
  double d_wp;      // overrelaxation parameter
  RandomNumbers* d_rnd;
  Index_eo d_index;

 public:
  GaugeFixing_Coulomb(RandomNumbers* rnd){
    d_rnd = rnd;
    d_Niter = 0;
    //set_prms();
  };

  ~GaugeFixing_Coulomb(){
  };

  void set_prms(int Niter, int Nnaive, int Nmeas,
                int Nreset, double Enorm, double wp);

  void set_prms();

  void fix(Field_G& Ufix, const Field_G& Uorg);

  void gauge_trf(Field_G& Ue, Field_G& Uo, Field_G& Geo, int Ieo);

  void set_randomGtrf(std::valarray<double>& sg, Field_G& Geo);

  //! one step of gauge fixing with overrelaxation parameter wp.
  void gfix_step(std::valarray<double>& sg,
                 Field_G& Ue, Field_G& Uo, double wp);

  void calc_SG(std::valarray<double>& sg, std::valarray<double>& Fval,
               Field_G& Ue, Field_G& Uo);
  void calc_W(Field_G& Weo, Field_G& Ue, Field_G& Uo, int Ieo);
  void calc_DLT(Field_G& Weo, Field_G& Ue, Field_G& Uo, int Ieo);

  void maxTr(Field_G&, Field_G&);
  void maxTr1(Field_G&, Field_G&);
  void maxTr2(Field_G&, Field_G&);
  void maxTr3(Field_G&, Field_G&);

  void sum_global_t(std::valarray<double>& val_global,
                    std::valarray<double>& val_local);

};

#endif  
