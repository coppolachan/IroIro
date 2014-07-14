#include "test_Wilson_FiniteDensity.hpp"
#include "Dirac_ops/dirac_wilson_FiniteDensity.hpp"
#include "Dirac_ops/dirac_clover.hpp"
using namespace std;

void Test_Wilson_FiniteDensity::run(){
  AntiPeriodicBC<GaugeField> apbc(TDIR);
  apbc.apply_bc(conf_);


  Dirac_Wilson_FiniteDensity Dw(node_,&(conf_.data));

  CCIO::cout<<"Dirac_Wilson_FiniteDensity created\n";  

  Field b(Dw.fsize(),1.0);
  Format::Format_F ff(CommonPrms::Nvol());
    
  Field Db = Dw.mult(b);
  double Dbnrm = Db.norm();

  Db = Dw.mult_dag(b);
  double Ddbnrm = Db.norm();

  CCIO::cout<<"Dbnrm= "<<Dbnrm<<" Ddbnrm= "<<Ddbnrm<<"\n";

  Dirac_Clover Dcl(&Dw,1.0,&(conf_.data));
  Db = Dcl.mult(b);
  Dbnrm = Db.norm();

  Db = Dcl.mult_dag(b);
  Ddbnrm = Db.norm();

  CCIO::cout<<"Clover: Dbnrm= "<<Dbnrm<<" Ddbnrm= "<<Ddbnrm<<"\n";

  /*
  int col_i = 0;    int col_j = 0; 
  int spn_i = 0;    int spn_j = 0; 

  for(int i=0; i<ff.Nvol();++i){
    for(int j=0; j<ff.Nvol();++j){

      //  int i = 384; int j = 448;

      b = 0.0;
      b.set(ff.index_r(col_j,spn_j,j),1.0);
      b.set(ff.index_i(col_j,spn_j,j),0.0);
      Field Db = Dw.mult(b);
      double Db_r = Db[ff.index_r(col_i,spn_i,i)];
      double Db_i = Db[ff.index_i(col_i,spn_i,i)];

      b = 0.0;
      b.set(ff.index_r(col_i,spn_i,i),1.0);
      b.set(ff.index_i(col_i,spn_i,i),0.0);
      Db = Dw.mult_dag(b);
      double Ddb_r = Db[ff.index_r(col_j,spn_j,j)];
      double Ddb_i = Db[ff.index_i(col_j,spn_j,j)];

      if(Db_r != Ddb_r || Db_i != -Ddb_i) 
	CCIO::cout<<"("<<i<<","<<j<<")"
		  <<" Dij= ("    <<Db_r <<","<<Db_i <<")"
		  <<" Ddag_ji= ("<<Ddb_r<<","<<Ddb_i<<")\n";

    }
  }
  */
}
