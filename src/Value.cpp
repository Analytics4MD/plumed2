#include "Value.h"
#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"
#include "ActionWithVirtualAtom.h"
#include "PlumedException.h"
#include "Atoms.h"
#include "PlumedMain.h"

using namespace PLMD;

Value::Value(ActionWithValue&action,const std::string& name):
  action(action),
  value(0.0),
  name(name),
  deriv(false),
  periodicity(unset),
  min(0.0),
  max(0.0),
  max_minus_min(0.0),
  inv_max_minus_min(0.0)
{}

bool Value::isPeriodic()const{
  plumed_massert(periodicity!=unset,"periodicity should be set");
  return periodicity==periodic;
}

void Value::getDomain(double&min,double&max)const{
  min=this->min;
  max=this->max;
}

void Value::setPeriodicity(bool p){
  if(p) periodicity=periodic;
  else periodicity=notperiodic;
}

void Value::setDomain(double min,double max){
  this->min=min;
  this->max=max;
  max_minus_min=max-min;
  if(max_minus_min!=0.0) inv_max_minus_min=1.0/max_minus_min;
}


const std::string Value::getFullName()const{
  if(name.length()==0) return action.getLabel();
  else return action.getLabel()+"."+name;
}

void Value::enableDerivatives()
{
  deriv=true;derivatives.resize(action.getNumberOfParameters());
}

void Value::setGradients(){
  gradients.clear();
  ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(&action);
  ActionWithArguments*aw=dynamic_cast<ActionWithArguments*>(&action);
  if(aa){
    Atoms&atoms((aa->plumed).getAtoms());
    for(unsigned j=0;j<aa->getNumberOfAtoms();++j){
      AtomNumber an=aa->getAbsoluteIndex(j);
      if(atoms.isVirtualAtom(an.index())){
        const ActionWithVirtualAtom* a=atoms.getVirtualAtomsAction(an.index());
        for(std::map<AtomNumber,Tensor>::const_iterator p=a->getGradients().begin();p!=a->getGradients().end();++p){
// controllare l'ordine del matmul:
          gradients[(*p).first]+=matmul(Vector(derivatives[3*j],derivatives[3*j+1],derivatives[3*j+2]),(*p).second);
        }
      } else {
        for(unsigned i=0;i<3;i++) gradients[an][i]+=derivatives[3*j+i];
      }
    }
  } else if(aw){
    std::vector<Value*> values=aw->getArguments();
    for(unsigned j=0;j<derivatives.size();j++){
      for(std::map<AtomNumber,Vector>::const_iterator p=values[j]->gradients.begin();p!=values[j]->gradients.end();++p){
        AtomNumber iatom=(*p).first;
        gradients[iatom]+=(*p).second*derivatives[j];
      }
    }
  } else plumed_error();
}

double Value::projection(const Value& v1,const Value&v2){
  double proj=0.0;
  const std::map<AtomNumber,Vector> & grad1(v1.gradients);
  const std::map<AtomNumber,Vector> & grad2(v2.gradients);
  for(std::map<AtomNumber,Vector>::const_iterator p1=grad1.begin();p1!=grad1.end();++p1){
    AtomNumber a=(*p1).first;
    std::map<AtomNumber,Vector>::const_iterator p2=grad2.find(a);
    if(p2!=grad2.end()){
      proj+=dotProduct((*p1).second,(*p2).second);
    }
  }
  return proj;
}




