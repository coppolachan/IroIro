/*!@file eigenModes.cpp
 * @brief implementation of the factory function for Eigen::Predic
 */

#include "eigenModes.hpp"
#include "PugiXML/xmlUtilities.hpp"

namespace Eigen{
  Predic* predFactory(XML::node node){
    XML::nullCheck(node,"eigen-reading condition ");
    const char* condition = node.attribute("name").value();
    if(!strcmp(condition, "BelowEq"))    return new BelowEq(node);
    if(!strcmp(condition, "Beyond"))     return new Beyond(node);
    if(!strcmp(condition, "SmallerEq"))  return new SmallerEq(node);
    if(!strcmp(condition, "LargerThan")) return new LargerThan(node);
    return new Predic;  /// default behavior (read as it is)
  }
}
