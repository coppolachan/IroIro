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

  CCIO::cout<<" ---- Calculating Wilson flow\n"; 
  vector<double> tau;
  vector<double> ttEstd;
  vector<double> ttEsym;

  int Mstep = wflow.MonitorStep();
  
  for(int t=0; t<wflow.Nstep(); ++t){            // wilson flow 
    tau.push_back(wflow.tau(t));
    ttEstd.push_back(wflow.Edens_plaq(t));
    ttEsym.push_back(wflow.Edens_clover(t));

    if(!(t%Mstep)){
      CCIO::cout<<"Monitor at t="<<t<<"\n";
      monitor(wflow.getU());
    }
    wflow.evolve_step();
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

void Test_WilsonFlow::monitor(const GaugeField& U)const{
  TopologyGeom tg;                             // geometrical topology
  double Qt = tg.getQ(U);
  CCIO::cout<<" GeometricalTopology = "<<Qt<<"\n";

  Staples stpl;                                // check of the smoothness 
  double sp_max = NC_*(1.0-stpl.plaq_min(U));
  double sp_ave = NC_*(1.0-stpl.plaquette(U));
  static double sp_adm = 0.067;                // admissible threshold

  CCIO::cout<<" sp_max = "        << sp_max <<"\n";
  CCIO::cout<<" sp_ave = "        << sp_ave <<"\n"; 
  CCIO::cout<<" (sp_admissible = "<< sp_adm <<")\n"; 
  CCIO::cout<<" sp_admissible - sp_max = "<<sp_adm-sp_max <<"\n";
  CCIO::cout<<" sp_admissible - sp_ave = "<<sp_adm-sp_ave <<"\n";
}
