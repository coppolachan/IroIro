#include "map_scalar.hpp"
#include "siteIndex.h"

Map::Map(int dir, ShiftSign sign) {
  sites = CommonPrms::instance()->Lvol(); //assuming single node
  CCIO::cout << "Mapper scalar " << dir << " constructing. Sites: "<< sites<<"\n";
  int Block[4];
  Block[0] =          CommonPrms::instance()->Nx();
  Block[1] = Block[0]*CommonPrms::instance()->Ny();
  Block[2] = Block[1]*CommonPrms::instance()->Nz();
  Block[3] = Block[2]*CommonPrms::instance()->Nt();  

  SiteIndex* Sindex = SiteIndex::instance();
  site_map.resize(sites);

  for (int site = 0; site < sites; ++site) {
    if (sign == Forward) {
      site_map[site] = Sindex->x_p(site,dir);
      if(Sindex->cmp(site,dir)==Sindex->Bdir(dir)) site_map[site] -= Block[dir];
    }
    if (sign == Backward) {
      site_map[site] = Sindex->x_m(site,dir);
      if(Sindex->cmp(site,dir)==0) site_map[site] += Block[dir];
    }   
  }
}

