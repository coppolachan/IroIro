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
  vector<double> Eplaq = wflow.evolve();

  // output
  CCIO::cout << " ---- Output in "<< output_.c_str()<<"\n";
  if(Communicator::instance()->primaryNode()){
    ofstream writer(output_.c_str());

    writer<< setiosflags(  ios_base::scientific);
    for(int t=0; t<Eplaq.size(); ++t){
      writer<< setw(2) <<setiosflags(ios_base::right)<< t;
      writer<< setw(20)<<setiosflags(ios_base::left )<< Eplaq[t]<<endl;
    }
    writer<< resetiosflags(ios_base::scientific);
    writer.close();
  }
  Communicator::instance()->sync();
  return 0;
}
