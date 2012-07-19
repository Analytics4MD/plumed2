/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The PLUMED team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of PLUMED, version 2.0.

   PLUMED is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PLUMED is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with PLUMED.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_Grex_h
#define __PLUMED_Grex_h

#include "WithCmd.h"
#include <string>
#include <vector>

namespace PLMD{

class PlumedMain;
class Atoms;
class PlumedCommunicator;

class GREX:
  public WithCmd
{
  bool initialized;
  PlumedCommunicator& intracomm;
  PlumedCommunicator& intercomm;
  PlumedMain& plumedMain;
  Atoms&      atoms;
  int partner;
  double localDeltaBias;
  double foreignDeltaBias;
  std::vector<double> allDeltaBias;
  std::string buffer;
  int myreplica;
public:
  GREX(PlumedMain&);
  ~GREX();
  void cmd(const std::string&key,void*val=NULL);
  void calculate();
  void savePositions();
};

}

#endif
