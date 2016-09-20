/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include <cmath>

using namespace std;

namespace PLMD{
namespace bias{

//+PLUMEDOC BIAS EMRESTRAINT
/*
Put the doc here

*/
//+ENDPLUMEDOC


class EMrestraint : public Bias
{
  // temperature in kbt
  double kbt_;
  // exp data points
  vector<double> ovdd_;
  
public:
  EMrestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(EMrestraint,"EMRESTRAINT")

void EMrestraint::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TEMP","temperature in energy units");
  componentsAreNotOptional(keys); 
}

EMrestraint::EMrestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
{

  parse("TEMP",kbt_);
  checkRead();

  // last half of the arguments inside ovdd_;
  for(unsigned i=getNumberOfArguments()/2; i<getNumberOfArguments();++i){
    ovdd_.push_back(getArgument(i));
  }

  log.printf("  temperature of the system %f\n",kbt_);
  log.printf("  number of experimental data points %u\n",static_cast<unsigned>(ovdd_.size()));
  
}


void EMrestraint::calculate(){
   
  double ene = 0.0;
  unsigned int ndata = getNumberOfArguments()/2;
  
  vector<double> ene_der(ndata);
  
  // cycle on arguments
  // count number of non-zero overlaps
  double ndata_zero = 0.0;
  for(unsigned i=0;i<ndata;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0){
     // individual term
     ene_der[i] = std::log(ovmd/ovdd_[i]);
     // increment energy
     ene += ene_der[i] * ene_der[i];
     // increment counter
     ndata_zero += 1.0;
    }
  };
  
  // constant factor
  double fact = kbt_ * 0.5 * ndata_zero;

  // get derivatives
  for(unsigned i=0;i<ndata;++i){
    // check for zero overlaps
    double ovmd = getArgument(i);
    if(ovmd > 0.0 && ene > 0.0){
     // calculate derivative
     double der = 2.0 * fact / ene * ene_der[i] / ovmd;
     // set forces
     setOutputForce(i, -der);
    }
  }

  // set value of the bias
  setBias(fact * std::log(ene));
}


}
}


