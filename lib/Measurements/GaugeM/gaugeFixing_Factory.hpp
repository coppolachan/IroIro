
class GaugeFixingFactory{
  virtual ~GFixFactory(){}
  virtual GaugeFixing* getGaugeFixing() = 0;
};

class GFixFactory_Free: public GFixFactory {
  GaugeFixing* getGaugeFixing(){
    return new GaugeFixing_Free;
  }
};

class GFixFactory_Landau: public GFixFactory {
  RaiiFactoryObj<RandomNumberCreator> RNGobj_;
  RaiiFactoryObj<RandomNumber> RNG_;
  const XML::node gfix_node_;
public:
  GaugeFixingFactory(XML::node node):gfix_node_(node){
    XML::descend(node,"RandomNumberGen");
    RNGobj_.save(RNG_Env::createRNGfactory(node));
    XML::descend(gfix_node_,"Fixing");
  }  
  
  GaugeFixing* getGaugeFixing(){
    RNG_.save(RNGobj_.get()->getRandomNumberGen());
    return new GaugeFixing_Landau(*(RNG_.get()),gfix_node_);
  }
};

class GFixFactory_Coulomb: public GFixFactory {
  RaiiFactoryObj<RandomNumberCreator> RNGobj_;
  RaiiFactoryObj<RandomNumber> RNG_;
  const XML::node gfix_node_;
public:
  GaugeFixingFactory(XML::node node):gfix_node_(node){
    XML::descend(node,"RandomNumberGen");
    RNGobj_.save(RNG_Env::createRNGfactory(node));
    XML::descend(gfix_node_,"Fixing");
  }  
  
  GaugeFixing* getGaugeFixing(){
    RNG_.save(RNGobj_.get()->getRandomNumberGen());
    return new GaugeFixing_Coulomb(*(RNG_.get()),gfix_node_);
  }
};

