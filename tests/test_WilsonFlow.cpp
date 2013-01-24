/*! @file test_WilsonFlow.cpp
    @brief implementation of test code for the wilson flow
*/
#include "test_WilsonFlow.hpp"
#include "Measurements/GaugeM/wilsonFlow.hpp"
#include "Measurements/GaugeM/topologyGeom.hpp"
#include <string>

using namespace std;

int Test_WilsonFlow::run(){
  XML::descend(node_,"WilsonFlow");               // object creation
  WilsonFlow wflow(node_,conf_);
  int Nvol=CommonPrms::instance()->Nvol();
  vector<double> q(Nvol);                            // topology density

  CCIO::cout<<" ---- Calculating Wilson flow\n"; 
  vector<double> tau;
  vector<double> ttEstd;
  vector<double> ttEsym;

  int Mstep = wflow.MonitorStep();
  
  for(int t=0; t<wflow.Nstep(); ++t){            // wilson flow 
    tau.push_back(wflow.tau(t));
    ttEstd.push_back(wflow.Edens_plaq(t));
    ttEsym.push_back(wflow.Edens_clover(t));
    wflow.evolve_step();
    if(!((t+1)%Mstep)){
      CCIO::cout<<"Monitor at t="<<t+1<<"\n";
      monitor(wflow.getU(),q);                      // topology monitor
      topologyoutput(output_,t+1,q);                // saving topology density
      /*
      std::stringstream ofile_tmp;
      ofile_tmp << output_.c_str() <<"_cooledconf_at_t" << t+1;
      wflow.save_config(ofile_tmp.str().c_str());
      */
    }
  }

  ////// File Output //////
  wflow.save_config(output_.c_str());    // saving evolved config (if required)
  std::stringstream ofile;
  ofile << output_.c_str() <<"_ttE";     // output of ttE

  CCIO::cout << " ---- Output in "<< ofile.str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(ofile.str().c_str());
    for(int t=0; t<wflow.Nstep(); ++t)
      writer<<setw(10)<<fixed<<setprecision( 6)<<setiosflags(ios_base::left)<<tau[t]
	    <<setw(20)<<fixed<<setprecision(16)<<setiosflags(ios_base::left)<<ttEstd[t]
	    <<setw(20)<<fixed<<setprecision(16)<<setiosflags(ios_base::left)<<ttEsym[t]
	    <<endl;
    writer.close();
  }
  Communicator::instance()->sync();
  return 0;
}

void Test_WilsonFlow::monitor(const GaugeField& U, vector<double>& q)const{

  TopologyGeom tg;                             // geometrical topology
  double Qt = tg.get_Q(q,U);                   
  CCIO::cout<<" GeometricalTopology = "<<Qt<<"\n";

  // check 
  //int Nvol=CommonPrms::instance()->Nvol();
  //double Qcheck = 0.0;  
  //for(int site=0; site<Nvol; ++site){
  //  Qcheck += q[site];
  //}
  //Qcheck = Communicator::instance()->reduce_sum(Qcheck);
  //CCIO::cout<<"Check ! GeometricalTopology = "<<Qcheck<<"\n";
  // check

  CCIO::cout<<" Admissibility check:\n"; 
  Staples stpl;                                // check of the smoothness 
  static double sp_adm = 0.067;                // admissible threshold
  static double pl_adm = 1.0-sp_adm/NC_;
  CCIO::cout << "   (pl_adm =" << pl_adm << ")\n";

  double sp_max = NC_*(1.0-stpl.plaq_min(U,pl_adm));
  double sp_ave = NC_*(1.0-stpl.plaquette(U));

  CCIO::cout<<"   sp_max = "        << sp_max <<"\n";
  CCIO::cout<<"   sp_ave = "        << sp_ave <<"\n"; 
  CCIO::cout<<"   (sp_admissible = "<< sp_adm <<")\n"; 
  CCIO::cout<<"   sp_admissible - sp_max = "<<sp_adm-sp_max <<"\n";
  CCIO::cout<<"   sp_admissible - sp_ave = "<<sp_adm-sp_ave <<"\n";
}

void Test_WilsonFlow::topologyoutput(string& fname, int t, vector<double>& q)const{

  int Lvol = CommonPrms::instance()->Lvol();
  int Nvol=CommonPrms::instance()->Nvol();

  // gathering topology density
  std::vector<double> q_tmp(Lvol,0.0);
  for(int site=0; site<Nvol; ++site){
    q_tmp[SiteIndex::instance()->get_gsite(site)] = q[site];
  }
  std::vector<double> q_all(Lvol);
  for(int site=0; site<Lvol; ++site)
    q_all[site] = Communicator::instance()->reduce_sum(q_tmp[site]);

  // output
  std::stringstream ofile;
  ofile << fname.c_str() <<"_topdensity_at_t"<< t;     // output of ttE

  CCIO::cout << " ---- Output in "<< ofile.str()<<"\n";
  if(Communicator::instance()->primaryNode()){

    /* out put in txt
    ofstream writer(ofile.str().c_str());
    for(int site=0; site<Lvol; ++site)
      writer<<setw(10)<<fixed<<setprecision( 6)<<setiosflags(ios_base::left)<<site
	    <<setw(20)<<fixed<<setprecision(16)<<setiosflags(ios_base::left)<<q_all[site]
	    <<endl;
    writer.close();
    out put in txt */

    FILE* outFile;
    outFile = fopen(ofile.str().c_str(),"w");
    fwrite((const char*)&q_all[0], sizeof(double), Lvol, outFile);

    /* check start
    outFile = fopen(ofile.str().c_str(),"r");
    std::vector<double> q_check(Lvol);
    fread((char*)&q_check[0], sizeof(double), Lvol, outFile);
    ofile << "_check";
    ofstream writer(ofile.str().c_str());
    for(int site=0; site<Lvol; ++site)
      writer<<setw(10)<<fixed<<setprecision( 6)<<setiosflags(ios_base::left)<<site
	    <<setw(20)<<fixed<<setprecision(16)<<setiosflags(ios_base::left)<<q_all[site]<<q_check[site]
	    <<q_all[site]-q_check[site]
	    <<endl;
    writer.close();
     check end */ 
    
  }
  Communicator::instance()->sync();
}
