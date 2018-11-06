/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
   (see the PEOPLE target at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This target is part of plumed, version 2.

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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/Units.h"
#include <cstdio>
#include <memory>
#include "core/SetupMolInfo.h"
#include "core/ActionSet.h"
#include "retrieve.h"
#include <vector>
#include <tuple>
#include <chrono>

#ifndef __PLUMED_DOES_LICHENS_DISPATCH
#define __PLUMED_DOES_LICHENS_DISPATCH
#endif
#if defined(__PLUMED_DOES_LICHENS_DISPATCH)
//#include <lichenstarget/lichens_target.h>
#endif

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::milliseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

using namespace std;

namespace PLMD
{
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPATOMS
/*

*/
//+ENDPLUMEDOC

class DispatchAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  double lenunit;
  Retrieve* retrieve_ptr;
#if defined(__PLUMED_DOES_LICHENS_DISPATCH)
  int temp=0;
#endif

public:
  explicit DispatchAtoms(const ActionOptions&);
  ~DispatchAtoms();
  static void registerKeywords( Keywords& keys );
  void calculate() {}
  void apply() {}
  void update();
};

PLUMED_REGISTER_ACTION(DispatchAtoms,"DISPATCHATOMS")

void DispatchAtoms::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the atoms should be output");
  keys.add("atoms", "ATOMS", "the atom indices whose positions you would like to print out");
  keys.add("compulsory", "TARGET", "target on which to output coordinates; extension is automatically detected");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
}

DispatchAtoms::DispatchAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao)
{
  retrieve_ptr = new Retrieve();
  vector<AtomNumber> atoms;
  string target;
  parse("TARGET",target);
  if(target.length()==0) error("name out output target was not specified");
  
  parseAtomList("ATOMS",atoms);

  std::string unitname; parse("UNITS",unitname);
  lenunit=1.0;

  checkRead();

  log.printf("  dispatching the following atoms in %s :", unitname.c_str() );
  for(unsigned i=0; i<atoms.size(); ++i) log.printf(" %d",atoms[i].serial() );
  log.printf("\n");
  requestAtoms(atoms);
}

void DispatchAtoms::update() {
  //Get the number of processes
  TimeVar t1=timeNow();
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  //printf("Number of atom positions being dispatched: %d rank %d/%d\n",getNumberOfAtoms(),world_rank,world_size);
  if (world_rank ==0)
  {
    const Tensor & t(getPbc().getBox());
    double lx, ly, lz, xy, xz, yz; //xy, xz, yz are tilt factors 
    lx = lenunit*t(0,0);
    ly = lenunit*t(1,1);
    lz = lenunit*t(2,2);
    xy = lenunit*t(0,1); // 0 for orthorhombic
    xz = lenunit*t(0,2); // 0 for orthorhombic
    yz = lenunit*t(1,2); // 0 for orthorhombic

    int atom_count = getNumberOfAtoms();

    std::vector<std::tuple<double, double, double> > positions;
    int types[atom_count];
    for(int i=0; i<atom_count; ++i) {
      auto pos_tuple = std::make_tuple(lenunit*getPosition(i)(0), lenunit*getPosition(i)(1), lenunit*getPosition(i)(2));
      positions.push_back(pos_tuple);
      types[i] = 0;
    }
    //printf("Calling retrieve call\n");
    char* module_name = "calc_voronoi_for_frame";
    char* function_name = "analyze";

    retrieve_ptr->analyze_frame(module_name,
                           function_name,
                           types,
                           positions,
                           lx,
                           ly,
                           lz,
                           xy,
                           xz,
                           yz);
  }
  int time_elapsed = duration(timeNow()-t1);
  printf("DISPATCH_UPDATE_TIME_rank_%d : %d\n",world_rank,time_elapsed);
}

DispatchAtoms::~DispatchAtoms() {
}


}
}
