/*! @file test_WilsonFlow.cpp
    @brief implementation of test code for the wilson flow
*/
#include "test_WilsonFlow.hpp"
#include "Measurements/GaugeM/wilsonFlow.hpp"

using namespace std;

int Test_WilsonFlow::run(){
  // object creation
  XML::descend(node_,"WilsonFlow");
  WilsonFlow wflow(node_,conf_);

  // wilson flow 
  CCIO::cout<<" ---- Calculating Wilson flow\n";
  vector<double> tau = wflow.get_t();
  vector<double> ttE = wflow.evolve();

  // output
  CCIO::cout << " ---- Output in "<< output_.c_str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(output_.c_str());

    //writer<< setiosflags(  ios_base::scientific);
    for(int t=0; t<ttE.size(); ++t)
      writer<<setw(10)<<fixed<<setprecision( 6)<<setiosflags(ios_base::left)<<tau[t]
	    <<setw(20)<<fixed<<setprecision(16)<<setiosflags(ios_base::left)<<ttE[t]
	    <<endl;
    //writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  return 0;
}
