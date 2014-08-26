/*!
 * @file BFM_HDCG.hpp
 * @brief Declares classes for P. Boyle HDCG inverter
 * Time-stamp: <2014-08-07 15:25:22 neo>
 */
#ifndef _BFM_HDCG_EXT_H_
#define _BFM_HDCG_EXT_H_

#include <BfmHDCG.h>



template <class Float> 
class BFM_HDCG_Extend : public BfmHDCG < Float >
{
public:
  BFM_HDCG_Extend(int _N5,
		  int _Nvec,
		  int _BlockSize[5],
		  int _QuadrantSize[4],
		  int _GlobalSize[4],
		  int _SubgridSize[4],
		  int _NodeCoord[4],
		  bfm_internal<double> * _linop_d,
		  bfm_internal<float>  * _linop_f) : BfmHDCG <Float> (_N5,
								      _Nvec,
								      _BlockSize,
								      _QuadrantSize,
								      _GlobalSize,
								      _SubgridSize,
								      _NodeCoord,
								      _linop_d,
								      _linop_f){};



  int MyNodeNumber();
  int NodeFromCoord(int g[4]);
};

#endif
