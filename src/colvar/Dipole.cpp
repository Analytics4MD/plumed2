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

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR DIPOLE 
/*
Calcualte the dipole moment for a group of atoms.

\par Examples
The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\verbatim
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=5
\endverbatim
(see also \ref PRINT)

\attention 
If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
where N is the number of atoms. This implies that the dipole (which for a charged system depends
on the position) is computed on the geometric center of the group.


*/
//+ENDPLUMEDOC
   
class Dipole : public Colvar {
  vector<AtomNumber> ga_lista;
public:
  Dipole(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Dipole,"DIPOLE")

void Dipole::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.remove("NOPBC");
}

Dipole::Dipole(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  parseAtomList("GROUP",ga_lista);
  checkRead();
  addValueWithDerivatives(); setNotPeriodic();

  log.printf("  of %d atoms\n",ga_lista.size());
  for(unsigned int i=0;i<ga_lista.size();++i){
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  requestAtoms(ga_lista);
}

// calculator
void Dipole::calculate()
{
 double dipole=0.;
 vector<Vector> deriv(getNumberOfAtoms());
 Vector dipje;
 vector<double> charges(getNumberOfAtoms());

 double ctot(0.0);
 for(unsigned i=0;i<charges.size();++i){
   charges[i]=getCharge(i);
   ctot+=charges[i];
 }

 ctot/=charges.size();

 for(unsigned i=0;i<charges.size();++i) charges[i]-=ctot;

// deriv.resize(getPositions().size());
// deriv.resize(getNumberOfAtoms());
 for(unsigned int i=0;i<ga_lista.size();i++) {
   dipje += charges[i]*getPosition(i);
 }
 dipole = dipje.modulo();

 for(unsigned int i=0;i<ga_lista.size();i++) {
   double dfunc=charges[i]/dipole;
   deriv[i] = deriv[i] + (dfunc)*dipje;
 }

// for(unsigned i=0;i<getPositions().size();++i) setAtomsDerivatives(i,deriv[i]);
 for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(i,deriv[i]);
 setValue           (dipole);
 setBoxDerivativesNoPbc();
}

}
}
