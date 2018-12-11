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
#include "dataspaces_writer.h"
//#include "plumed_chunker.h"
#include "chunk.h"
#include "common.h"
#include <vector>
#include <tuple>
#include <chrono>
#include <string>

#ifndef __PLUMED_DOES_LICHENS_DISPATCH
#define __PLUMED_DOES_LICHENS_DISPATCH
#endif
#if defined(__PLUMED_DOES_LICHENS_DISPATCH)
//#include <lichenstarget/lichens_target.h>
#endif

typedef std::chrono::high_resolution_clock::time_point TimeVar;
typedef std::chrono::duration<double, std::milli> DurationMilli;
#define timeNow() std::chrono::high_resolution_clock::now()

using namespace std;

namespace PLMD
{
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPATOMS
/*
The analysis function should have the following syntax:
"analyze(types, points, box_points, step)" where 
types is an array of integers,\n\
points is an array of array of cartesian coordinates [[x1,y1,z1],[x2,y2,z3],..[xn,yn.zn]]\n\
box_points is an array of double values representing the simulation box. The values are: [Lx,Ly,Ly,xy,yz,xz]\"
*/
//+ENDPLUMEDOC

class DispatchAtoms:
  public ActionAtomistic,
  public ActionPilot
{
  double lenunit;
  Retrieve* retrieve_ptr;
  DataSpacesWriter* dataspaces_writer_ptr;
  //PlumedChunker* chunker_ptr;
  int dispatch_method = 0; // 1: a4md, 2:python 
#if defined(__PLUMED_DOES_LICHENS_DISPATCH)
  int temp=0;
#endif
  int world_size;
  int world_rank;
  double total_simulation_time_ms = 0.0;
  double simulation_time_ms = 0.0;
  double total_plumed_time_ms = 0.0;
  DurationMilli plumed_time_ms;
  //int stride = 0;
  string target;
  string python_function;
  int total_steps = 0;
  int nstride = 0;
  TimeVar t_start;
  unsigned long int current_chunk_id = 0;
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
  keys.add("compulsory","TOTAL_STEPS","1","the total number of time steps being simulated");
  keys.add("compulsory", "TARGET","a4md","target on which to output coordinates; Options are \"a4md\" or a python module (provide the name of the .py file including the .py extension. for e.g. md_analysis.py)");
  keys.add("compulsory", "PYTHON_FUNCTION", "NONE", "name of the python module to load and execute the analysis code. Applicable only if TARGET is a .py file");
  keys.add("compulsory", "UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
}

DispatchAtoms::DispatchAtoms(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionPilot(ao)
{
  //Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  vector<AtomNumber> atoms;

  parse("TOTAL_STEPS",total_steps);
  parse("TARGET",target);
  //parse("STRIDE",nst);
  parse("PYTHON_FUNCTION",python_function);
  log.printf("TOTAL_STEPS: %i\n",total_steps);
  nstride = getStride();
  log.printf("STRIDE: %i\n",nstride);
  log.printf("TARGET: %s\n",target.c_str());
  log.printf("PYTHON_FUNCTION: %s\n",python_function.c_str());
  if(target.length()==0) error("name out output target was not specified");

  
  if (target.find("py") != std::string::npos)
  {
    printf("---====== DispatchAtoms found target to be a python module %s\n",target.c_str()); 
    if(python_function.length()==0) error("TARGET was specified as \"file with .py extension\", but PYTHON_MODULE was not specified");
    if (world_rank == 0)
    {
      retrieve_ptr = new Retrieve((char*)target.c_str(), (char*)python_function.c_str());
    }
    dispatch_method = 2;
  }
  else if (target == "a4md")
  {
    //chunker_ptr = new PlumedChunker();
    char* temp_var_name = "test_var";
    unsigned long int total_chunks = total_steps/nstride + 1;// +1 for the first call before starting simulation
    dataspaces_writer_ptr = new DataSpacesWriter(temp_var_name, total_chunks);
    dispatch_method = 1;
  }
  else
    log.printf(" ERROR: Unknown TARGET specified in DispatchAtoms Action.\n");

  parseAtomList("ATOMS",atoms);

  std::string unitname; parse("UNITS",unitname);
  lenunit=1.0;

  checkRead();

  //log.printf("  dispatching the following atoms in %s :", unitname.c_str() );
  //for(unsigned i=0; i<atoms.size(); ++i) log.printf(" %d",atoms[i].serial() );
  //log.printf("\n");
  requestAtoms(atoms);
}

void DispatchAtoms::update() {
    //printf("Number of atom positions being dispatched: %d rank %d/%d\n",getNumberOfAtoms(),world_rank,world_size);
  if (world_rank == 0)
  {
    auto step = getStep(); 
    if (step > 0)
    {
        DurationMilli round_about_time_ms = timeNow() - t_start;
        auto prev_simulation_time_ms = round_about_time_ms.count() - plumed_time_ms.count();
        total_simulation_time_ms += prev_simulation_time_ms;
        total_plumed_time_ms += plumed_time_ms.count();
    }
    if (step == total_steps)
    {
        log.printf("total_dispatch_action_time_ms : %f\n",total_plumed_time_ms);
        log.printf("total_simulation_time_ms : %f\n",total_simulation_time_ms);
        log.printf("total_time_steps : %d\n",total_steps);
    }
    t_start = timeNow();

    if (dispatch_method == 1)//a4md
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

      //printf("In DISPATCH step: %d\n",step);
      bool wf3 = true;
      if (dispatch_method==1)
        {
          //printf("Preparing data to send to dataspace, step %i\n",step);
          //std::vector<std::vector<double> > positions;
          std::vector<double> x_positions;
          std::vector<double> y_positions;
          std::vector<double> z_positions;

          std::vector<int> types;
          for(int i=0; i<atom_count; ++i) 
          {
            //std::vector<double> pos_tuple = {lenunit*getPosition(i)(0), lenunit*getPosition(i)(1), lenunit*getPosition(i)(2)};
            x_positions.push_back(lenunit*getPosition(i)(0));
            y_positions.push_back(lenunit*getPosition(i)(1));
            z_positions.push_back(lenunit*getPosition(i)(2));
            //positions.push_back(pos_tuple);
            types.push_back(0);
          }
          PLMDChunk plmd_chunk(current_chunk_id++,
                               step,
                               types,
                               x_positions,
                               y_positions,
                               z_positions,
                               lx,
                               ly,
                               lz,
                               xy,
                               xz,
                               yz);
          Chunk* chunk = &plmd_chunk; 
          std::vector<Chunk*> chunks = {chunk};
          dataspaces_writer_ptr->write_chunks(chunks);

        }
        else if (dispatch_method == 2)
        {
          std::vector<double> x_positions;
          std::vector<double> y_positions;
          std::vector<double> z_positions;

          std::vector<int> types;
          for(int i=0; i<atom_count; ++i) 
          {
            //std::vector<double> pos_tuple = {lenunit*getPosition(i)(0), lenunit*getPosition(i)(1), lenunit*getPosition(i)(2)};
            x_positions.push_back(lenunit*getPosition(i)(0));
            y_positions.push_back(lenunit*getPosition(i)(1));
            z_positions.push_back(lenunit*getPosition(i)(2));
            //positions.push_back(pos_tuple);
            types.push_back(0);
          }
          retrieve_ptr->analyze_frame(types,
                                      x_positions,
                                      y_positions,
                                      z_positions,
                                      lx,
                                      ly,
                                      lz,
                                      xy,
                                      xz,
                                      yz,
                                      step);
        }
        else
          log.printf(" ERROR: Unknown TARGET specified in DispatchAtoms Action.\n");
    }
    plumed_time_ms = timeNow()-t_start;
  }
}

DispatchAtoms::~DispatchAtoms() {
  dataspaces_writer_ptr = NULL;
  delete dataspaces_writer_ptr;
  retrieve_ptr = NULL;
  delete retrieve_ptr;

}
}
}
