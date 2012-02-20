#ifndef __PLUMED_Value_h
#define __PLUMED_Value_h

#include <vector>
#include <string>
#include <map>
#include "PlumedException.h"
#include "Tools.h"
#include "AtomNumber.h"
#include "Vector.h"

namespace PLMD{

class ActionWithValue;

/// Class containing a value which can be addressed by PLMD::ActionWithArguments.
/// It also contains the derivative of this value with respect to
/// an arbitrary number of parameters.
/// Typically, an object of type PLMD::ActionWithValue will contain one or
/// more objects of type PLUMD::Value, one per component.
class Value{
  ActionWithValue&action;
  double value;
  double inputForce;
  bool forced;
  std::vector<double> derivatives;
  std::map<AtomNumber,Vector> gradients;
  std::string name;
  bool deriv;
  enum {unset,periodic,notperiodic} periodicity;
  double min,max;
  double max_minus_min;
  double inv_max_minus_min;
public:
  void setGradients();
  Value(ActionWithValue&action,const std::string& name);
  void set(double);
  double get()const;
  void setPeriodicity(bool);
  void setDomain(double,double);
  bool isPeriodic()const;
  void getDomain(double&,double&)const;
  const std::string& getName()const;
  const std::string getFullName()const;
  void enableDerivatives();
  bool hasDerivatives()const;
  void setNumberOfParameters(int n);
  void setDerivatives(int i,double d);
  void clearInputForce();
  void clearDerivatives();
  double getForce()const;
  void  addForce(double f);
  const std::vector<double> &  getDerivatives()const;
  ActionWithValue& getAction();

  double difference(double)const;
  double difference(double,double)const;

/// check if a force has been added at this step
  bool checkForced()const;
  static double projection(const Value&,const Value&);

};

inline
void Value::set(double v){
  value=v;
}

inline
double Value::get()const{
  return value;
}

inline
const std::string& Value::getName()const{
  return name;
}

inline
ActionWithValue& Value::getAction(){
  return action;
}

inline
double Value::getForce()const{
  return inputForce;
}

inline
const std::vector<double> & Value::getDerivatives()const{
  return derivatives;
}

inline
bool Value::hasDerivatives()const{
  return deriv;
}

inline
void Value::setNumberOfParameters(int n){
  if(deriv)derivatives.resize(n);
}

inline
void Value::setDerivatives(int i,double d){
  derivatives[i]=d;
}

inline
void Value::clearInputForce(){
  forced=false;
  inputForce=0.0;
}

inline
void Value::clearDerivatives(){
  derivatives.assign(derivatives.size(),0.0);
}

inline
bool Value::checkForced()const{
  return forced;
}

inline
void Value::addForce(double f){
  plumed_massert(hasDerivatives(),"forces can only be added to values with derivatives");
  forced=true;
  inputForce+=f;
}

inline
double Value::difference(double d1,double d2)const{
  if(periodicity==notperiodic){
    return d2-d1;
  }else if(periodicity==periodic){
    double s=(d2-d1)*inv_max_minus_min;
    s=Tools::pbc(s);
    return s*max_minus_min;
  } else plumed_merror("periodicity should be set to compute differences");
}

inline
double Value::difference(double d)const{
  return difference(get(),d);
}

}

#endif

