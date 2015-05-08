#include "outputUtils.hpp"
#include "include/timings.hpp"
using namespace std;

namespace Hadrons{

  void output_meson(ofstream& writer,const vector<double>& Correl,string msg){
    int Lt = CommonPrms::instance()->Lt();
  
    if(Communicator::instance()->primaryNode()){
      writer<< setiosflags(ios_base::scientific);
      writer<<msg.c_str()<<endl;   // print of the separator

      for(int t=0; t<Lt; ++t){
	writer<< setw(2) <<setiosflags(ios_base::right)<< t
	      << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << Correl[t]
	      << endl;
      }
      writer<< resetiosflags(ios_base::scientific);
    }
  }

  void output_baryon(ofstream& writer,
		     const correl_t& Upper,const correl_t& Lower,
		     string msg){
    int Lt = CommonPrms::instance()->Lt();
    
    if(Communicator::instance()->primaryNode()){
      writer<< setiosflags(ios_base::scientific);
      writer<<msg.c_str()<<endl;  // print of the separator
      
      for(int t=0; t<Lt; ++t){
	writer<< setw(2) <<setiosflags(ios_base::right)<< t
	      << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << Upper[t].real()
	      << setw(25)<<setprecision(16)<<setiosflags(ios_base::left )
	      << Lower[t].real()
	      << endl;
      }
      writer<< resetiosflags(ios_base::scientific);
    }
  }
  
  /// Meson Channels ///
  void mesonProp(const prop_t& sq_ud,const prop_t& sq_s,ofstream& writer){
    long double total_timer, meson_timer;
    FINE_TIMING_START(total_timer);

    MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);

    CCIO::cout << " ::::::::::: light-light meson correlators\n";
    if(Communicator::instance()->primaryNode()) writer<<"====== meson correlators ======\n";

    FINE_TIMING_START(meson_timer);
    output_meson(writer,  pp.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Pion ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - Pion            :  "<< meson_timer << " seconds\n";
    FINE_TIMING_START(meson_timer);
    output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_ud),"------ Kaon ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - Kaon            :  "<< meson_timer << " seconds\n";
    FINE_TIMING_START(meson_timer);
    output_meson(writer,  pp.calculate<Format::Format_F>(sq_s, sq_s), "------ Eta_s ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - Eta_s           :  "<< meson_timer << " seconds\n";
 


    FINE_TIMING_START(meson_timer);
    output_meson(writer,v1v1.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_X ------");
    output_meson(writer,v2v2.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_Y ------");
    output_meson(writer,v3v3.calculate<Format::Format_F>(sq_ud,sq_ud),"------ Rho_Z ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - Rho(X,Y,Z)      :  "<< meson_timer << " seconds\n";

    FINE_TIMING_START(meson_timer);
    output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_X ------");
    output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_Y ------");
    output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_ud),"------ K_star_Z ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - K_star(X,Y,Z)   :  "<< meson_timer << " seconds\n";

    FINE_TIMING_START(meson_timer);
    output_meson(writer,v1v1.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_X ------");
    output_meson(writer,v2v2.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_Y ------");
    output_meson(writer,v3v3.calculate<Format::Format_F>(sq_s,sq_s),"------ Phi_Z ------");
    FINE_TIMING_END(meson_timer);
    CCIO::cout << " Timing - Phi(X,Y,Z)      :  "<< meson_timer << " seconds\n";


    FINE_TIMING_END(total_timer);
    CCIO::cout << " Timing - Cumulative      :  "<< total_timer << " seconds\n";

  }

  /// General Mesons ///
  void mesonPropGeneral(const prop_t& sq1,const prop_t& sq2,ofstream& writer){
    CCIO::cout << " ::::::::::: meson correlators\n";
    if(Communicator::instance()->primaryNode()) writer<<"====== meson correlators ======\n";

    long double total_timer, meson_timer;
    FINE_TIMING_START(total_timer);


    { /// basic channels

      MesonCorrelator pp(Pion), v1v1(Vector1), v2v2(Vector2), v3v3(Vector3);    

      FINE_TIMING_START(meson_timer);
      output_meson(writer,  pp.calculate<Format::Format_F>(sq1,sq2),"------PP------"  );
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - PP              :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer,v1v1.calculate<Format::Format_F>(sq1,sq2),"------V1V1------");
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - V1V1            :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer,v2v2.calculate<Format::Format_F>(sq1,sq2),"------V2V2------");
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - V2V2            :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer,v3v3.calculate<Format::Format_F>(sq1,sq2),"------V3V3------");
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - V3V3            :  "<< meson_timer << " seconds\n";
  
    }
    { /// extra channels
      MesonCorrelator ss(Scalar),a4a4(AVector4),a4p(AV4_PS),pa4(PS_AV4); 
      FINE_TIMING_START(meson_timer);
      output_meson(writer,  ss.calculate<Format::Format_F>(sq1,sq2),"------SS------"  );
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - SS              :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer,a4a4.calculate<Format::Format_F>(sq1,sq2),"------A4A4------");
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - A4A4            :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer, a4p.calculate<Format::Format_F>(sq1,sq2),"------A4P------" );
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - A4P             :  "<< meson_timer << " seconds\n";
      FINE_TIMING_START(meson_timer);
      output_meson(writer, pa4.calculate<Format::Format_F>(sq1,sq2),"------PA4------" );
      FINE_TIMING_END(meson_timer);
      CCIO::cout << " Timing - PA4              :  "<< meson_timer << " seconds\n";

    }

    FINE_TIMING_END(total_timer);
    CCIO::cout << " Timing - Cumulative      :  "<< total_timer << " seconds\n";

  }

  /// Extra Meson Channels ///
  void mesonExtraProp(const prop_t& sq_ud,const prop_t& sq_s,ofstream& writer){
    MesonCorrelator ss(Scalar),a4a4(AVector4),a4p(AV4_PS),pa4(PS_AV4); 

    CCIO::cout << " ::::::::::: Performing contractions for extra meson channels\n";
    if(Communicator::instance()->primaryNode()) writer<<"====== extra meson channels ======\n";

    output_meson(writer, ss.calculate<Format::Format_F>(sq_ud,sq_ud),"------ SS(ud,ud) ------");
    output_meson(writer, ss.calculate<Format::Format_F>(sq_s, sq_ud),"------ SS(s,ud) ------");
    output_meson(writer, ss.calculate<Format::Format_F>(sq_s, sq_s), "------ SS(s,s) ------");
    
    output_meson(writer,a4a4.calculate<Format::Format_F>(sq_ud,sq_ud),"------ A4A4(ud,ud) ------");
    output_meson(writer,a4a4.calculate<Format::Format_F>(sq_s, sq_ud),"------ A4A4(s,ud) ------");
    output_meson(writer,a4a4.calculate<Format::Format_F>(sq_s, sq_s), "------ A4A4(s,s) ------");
    /// <A4(t)P(0)>
    output_meson(writer,a4p.calculate<Format::Format_F>(sq_ud,sq_ud),"------ A4P(ud,ud) ------");
    output_meson(writer,a4p.calculate<Format::Format_F>(sq_s, sq_ud),"------ A4P(s,ud) ------");
    output_meson(writer,a4p.calculate<Format::Format_F>(sq_s, sq_s), "------ A4P(s,s) ------");
    /// <P(t)A4(0)>
    output_meson(writer,pa4.calculate<Format::Format_F>(sq_ud,sq_ud),"------ PA4(ud,ud) ------");
    output_meson(writer,pa4.calculate<Format::Format_F>(sq_s, sq_ud),"------ PA4(s,ud) ------");
    output_meson(writer,pa4.calculate<Format::Format_F>(sq_s, sq_s), "------ PA4(s,s) ------");
  }

  /// Baryon Correlation Functions ///
  void baryonProp(const prop_t& sq_ud,const prop_t& sq_s,ofstream& writer){
    long double total_timer, baryon_timer;
    FINE_TIMING_START(total_timer);

    BaryonCorrelator baryons(sq_ud,sq_s);

    CCIO::cout << " ::::::::::: Baryon correlators\n";
    if(Communicator::instance()->primaryNode()) writer<<"====== baryon correlators ======\n";


    FINE_TIMING_START(baryon_timer);
    output_baryon(writer,baryons.nucleon(UP),     baryons.nucleon(DN),     "------ Nucleon ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Nucleon              :  "<< baryon_timer << " seconds\n";
    FINE_TIMING_START(baryon_timer);
    output_baryon(writer,baryons.sigma8(UP),      baryons.sigma8(DN),      "------ Sigma8 ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Sigma8               :  "<< baryon_timer << " seconds\n";
    FINE_TIMING_START(baryon_timer);
    output_baryon(writer,baryons.xi8(UP),         baryons.xi8(DN),         "------ Xi8 ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Xi8                  :  "<< baryon_timer << " seconds\n";
    FINE_TIMING_START(baryon_timer);
    output_baryon(writer,baryons.lambda(UP),      baryons.lambda(DN),      "------ Lambda ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Lambda               :  "<< baryon_timer << " seconds\n";



    FINE_TIMING_START(baryon_timer);
    output_baryon(writer,baryons.delta(XDIR,UP),  baryons.delta(XDIR,DN),  "------ Delta_X ------");
    output_baryon(writer,baryons.delta(YDIR,UP),  baryons.delta(YDIR,DN),  "------ Delta_Y ------");
    output_baryon(writer,baryons.delta(ZDIR,UP),  baryons.delta(ZDIR,DN),  "------ Delta_Z ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Delta(X,Y,Z)         :  "<< baryon_timer << " seconds\n";


    FINE_TIMING_START(baryon_timer); 
    output_baryon(writer,baryons.omega(XDIR,UP),  baryons.omega(XDIR,DN),  "------ Omega_X ------");
    output_baryon(writer,baryons.omega(YDIR,UP),  baryons.omega(YDIR,DN),  "------ Omega_Y ------");
    output_baryon(writer,baryons.omega(ZDIR,UP),  baryons.omega(ZDIR,DN),  "------ Omega_Z ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Omega(X,Y,Z)         :  "<< baryon_timer << " seconds\n";


    FINE_TIMING_START(baryon_timer);     
    output_baryon(writer,baryons.sigma10(XDIR,UP),baryons.sigma10(XDIR,DN),"------ Sigma10_X ------");
    output_baryon(writer,baryons.sigma10(YDIR,UP),baryons.sigma10(YDIR,DN),"------ Sigma10_Y ------");
    output_baryon(writer,baryons.sigma10(ZDIR,UP),baryons.sigma10(ZDIR,DN),"------ Sigma10_Z ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Sigma10(X,Y,Z)       :  "<< baryon_timer << " seconds\n";


    FINE_TIMING_START(baryon_timer); 
    output_baryon(writer,baryons.xi10(XDIR,UP),   baryons.xi10(XDIR,DN),   "------ Xi10_X ------");
    output_baryon(writer,baryons.xi10(YDIR,UP),   baryons.xi10(YDIR,DN),   "------ Xi10_Y ------");
    output_baryon(writer,baryons.xi10(ZDIR,UP),   baryons.xi10(ZDIR,DN),   "------ Xi10_Z ------");
    FINE_TIMING_END(baryon_timer);
    CCIO::cout << " Timing - Xi10(X,Y,Z)          :  "<< baryon_timer << " seconds\n";   

    
    FINE_TIMING_END(total_timer);
    CCIO::cout << " Timing - Cumulative      :  "<< total_timer << " seconds\n";
  }


}
