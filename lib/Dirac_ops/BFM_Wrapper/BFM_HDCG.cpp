/*!
 * @file BFM_HDCG.cpp
 * @brief Declares classes for P. Boyle HDCG inverter
 * Time-stamp: <2014-10-16 14:34:31 cossu>
 */

#include "BFM_HDCG.hpp"
#include "Communicator/communicator.hpp"
#include "Geometry/siteIndex.hpp"
#include "Communicator/comm_io.hpp"

int MyNodeNumber(void){
  CCIO::cout << "My node number overloaded\n";
  return Communicator::instance()->nodeid();
}



int NodeFromCoord(int g[4]){
  //Move this into the site index class/////////////////////////

  SiteIndex* SIdx = SiteIndex::instance();
  CommonPrms* CP = CommonPrms::instance();
  Communicator* Comm = Communicator::instance();

  // Linearize global coordinate
  int gs = (g[0] + CP->Lx()*(g[1] + 
			     CP->Ly()*( g[2] +
					CP->Lz()*g[3])));

  // Return linear PE index

  /*
  return (SIdx->process_id_x(gs)+ CP->NPEx()*( SIdx->process_id_y(gs) +
						 CP->NPEy()*( SIdx->process_id_z(gs) +
								CP->NPEz()*SIdx->process_id_t(gs))));
  */

  return (SIdx->process_id_x(g[0])+ CP->NPEx()*( SIdx->process_id_y(g[1]) +
						 CP->NPEy()*( SIdx->process_id_z(g[2]) +
								CP->NPEz()*SIdx->process_id_t(g[3]))));


}





template class BFM_HDCG_Extend<double>;
template class BFM_HDCG_Extend<float>;
