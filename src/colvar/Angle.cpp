/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Angle.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR ANGLE
/*
Calculate an angle.

This command can be used to compute the angle between three atoms. Alternatively
if four atoms appear in the atom
specification it calculates the angle between the vector joining atoms 1 and 2 and that joining
atoms 3 and 4.

\par Examples

This command tells plumed to calculate the angle between the vector connecting atom 1 to atom 2 and
the vector connecting atom 2 to atom 3 and to print it on file COLVAR1. At the same time,
the angle between vector connecting atom 1 to atom 2 and the vector connecting atom 3 to atom 4 is printed
on file COLVAR2.
\verbatim

a: ANGLE ATOMS=1,2,3 
# equivalently one could state:
# a: ANGLE ATOMS=1,2,2,3

b: ANGLE ATOMS=1,2,3,4

PRINT ARG=a FILE=COLVAR1
PRINT ARG=b FILE=COLVAR2
\endverbatim
(see also \ref PRINT)

*/
//+ENDPLUMEDOC
   
class Angle : public Colvar {
  bool pbc;

public:
  Angle(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Angle,"ANGLE")

void Angle::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the list of atoms involved in this collective variable (either 3 or 4 atoms)");
}

Angle::Angle(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(atoms.size()==3){
    log.printf("  between atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  }else if(atoms.size()==4){
    log.printf("  between lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  }else error("Number of specified atoms should be either 3 or 4");

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
  checkRead();
}

// calculator
void Angle::calculate(){

  Vector dij,dik;
  if(pbc){
    dij=pbcDistance(getPosition(2),getPosition(3));
    dik=pbcDistance(getPosition(1),getPosition(0));
  } else {
    dij=delta(getPosition(2),getPosition(3));
    dik=delta(getPosition(1),getPosition(0));
  }
  Vector ddij,ddik;
  PLMD::Angle a;
  double angle=a.compute(dij,dik,ddij,ddik);
  setAtomsDerivatives(0,ddik);
  setAtomsDerivatives(1,-ddik);
  setAtomsDerivatives(2,-ddij);
  setAtomsDerivatives(3,ddij);
  setValue           (angle);
  setBoxDerivatives  (-(Tensor(dij,ddij)+Tensor(dik,ddik)));
}

}
}



