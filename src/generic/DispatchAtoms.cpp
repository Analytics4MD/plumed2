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
#include "dataspaces_reader.h"
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
  DataSpacesReader* dataspaces_reader_ptr;
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
  double total_retrieve_time_ms = 0.0;
  double total_read_plumed_data_time_ms = 0.0;
  TimeVar t_start;
  DurationMilli plumed_time_ms;
  //int stride = 0;
  string target;
  string python_module;
  string python_function;
  string data_stage;
  int total_steps = 0;
  int nstride = 0;
  unsigned long int current_chunk_id = 0;
  std::vector<std::string> ds_write_time_stamps;
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
  keys.add("compulsory", "TARGET","NONE","target on which to output coordinates; Options are \"a4md\" or a python module (provide the name of the .py file including the .py extension. for e.g. md_analysis.py)");
  keys.add("compulsory", "STAGE_DATA_IN", "NONE", "name of the python module to load and execute the analysis code. Applicable only if TARGET is a .py file");
  keys.add("compulsory", "PYTHON_MODULE", "NONE", "name of the python module to load and execute the analysis code. Applicable only if TARGET is a .py file");
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
  parse("STAGE_DATA_IN", data_stage);
  parse("PYTHON_MODULE",python_module);
  parse("PYTHON_FUNCTION",python_function);
  log.printf("TOTAL_STEPS: %i\n",total_steps);
  nstride = getStride();
  log.printf("STRIDE: %i\n",nstride);
  log.printf("TARGET: %s\n",target.c_str());
  log.printf("PYTHON_FUNCTION: %s\n",python_function.c_str());
  if(target == "NONE") error("name out output target was not specified");

  if (target == "py")
  {
    if(python_module=="NONE") error("TARGET was specified as \"py\", but PYTHON_MODULE was not specified");
    if(python_function=="NONE") error("TARGET was specified as \"py\", but PYTHON_FUNCTION was not specified");
    retrieve_ptr = new Retrieve((char*)python_module.c_str(), (char*)python_function.c_str());
    printf("----===== Initialized Retrieve ====----\n");
    if (data_stage == "dataspaces")
    {
      dispatch_method = 2;
      char* temp_var_name = "test_var";
      unsigned long int total_chunks = total_steps/nstride + 1;// +1 for the first call before starting simulation
      printf("----===== Initializing DataSpaces Reader and Writer ====----\n");
      dataspaces_writer_ptr = new DataSpacesWriter(temp_var_name, total_chunks);
      dataspaces_reader_ptr = new DataSpacesReader(temp_var_name, total_chunks);
      printf("----===== Initialized DataSpaces Reader and Writer ====----\n");
    }
    else
    { 
      dispatch_method = 1;
    }
  }
  else if (target == "a4md")
  {
    if (data_stage == "dataspaces")
    {
      char* temp_var_name = "test_var";
      unsigned long int total_chunks = total_steps/nstride + 1;// +1 for the first call before starting simulation
      dataspaces_writer_ptr = new DataSpacesWriter(temp_var_name, total_chunks);
      dispatch_method = 3;
    }
    else
      error("Invalid option for target = a4md. Currently USE_DATASPACES is the only valid option for a4md target.");
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
  auto step = getStep(); 
  //printf("---=== DispatchAtoms::update %i ===---\n",step);
  if (step > 0)
  {
      DurationMilli round_about_time_ms = timeNow() - t_start;
      auto prev_simulation_time_ms = round_about_time_ms.count() - plumed_time_ms.count();
      total_simulation_time_ms += prev_simulation_time_ms;
  }
  t_start = timeNow();

  if (dispatch_method > 0)
  {
    TimeVar t_start_read_frame = timeNow();
    const Tensor & t(getPbc().getBox());
    double lx, ly, lz, xy, xz, yz; //xy, xz, yz are tilt factors 
    lx = lenunit*t(0,0);
    ly = lenunit*t(1,1);
    lz = lenunit*t(2,2);
    xy = lenunit*t(0,1); // 0 for orthorhombic
    xz = lenunit*t(0,2); // 0 for orthorhombic
    yz = lenunit*t(1,2); // 0 for orthorhombic

    int atom_count = getNumberOfAtoms();
    
    std::vector<double> x_positions;
    std::vector<double> y_positions;
    std::vector<double> z_positions;

    std::vector<int> types;
    for(int i=0; i<atom_count; ++i) 
    {
      x_positions.push_back(lenunit*getPosition(i)(0));
      y_positions.push_back(lenunit*getPosition(i)(1));
      z_positions.push_back(lenunit*getPosition(i)(2));
      types.push_back(0);
    }
    DurationMilli read_plumed_data_time_ms = timeNow()-t_start_read_frame;
    total_read_plumed_data_time_ms += read_plumed_data_time_ms.count();
    //printf("In DISPATCH step: %d\n",step);
    bool wf3 = true;
    if (dispatch_method==1) // plumed
    {
      TimeVar t_start_retrieve = timeNow();
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
      
      DurationMilli retrieve_time_ms = timeNow()-t_start_retrieve;
      total_retrieve_time_ms += retrieve_time_ms.count();
    }
    else if (dispatch_method > 1) //plumed_ds
    {
      //printf("Preparing data to send to dataspace, step %i\n",step);
      PLMDChunk plmd_chunk(current_chunk_id,
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
      //printf("----===== Writing Chunk %i to DataSpaces START====----\n",current_chunk_id);
      dataspaces_writer_ptr->write_chunks(chunks);
      //printf("----===== Writing Chunk %i to DataSpaces STOP====----\n",current_chunk_id);
      if (dispatch_method == 2)
      {
        TimeVar t_start_retrieve = timeNow();
        //printf("----===== Reading Chunk %i from DataSpaces START====----\n",current_chunk_id);
        std::vector<Chunk*> in_chunks = dataspaces_reader_ptr->get_chunks(current_chunk_id, current_chunk_id);
        //printf("----===== Reading Chunk %i from DataSpaces STOP====----\n",current_chunk_id);
        for (Chunk* chunk: in_chunks)
        {
              PLMDChunk *plmdchunk = dynamic_cast<PLMDChunk *>(chunk);
              //printf("Printing typecasted chunk\n");
              //chunk->print();
              auto x_positions = plmdchunk->get_x_positions();
              auto y_positions = plmdchunk->get_y_positions();
              auto z_positions = plmdchunk->get_z_positions();
              auto types_vector = plmdchunk->get_types();
             
              double lx, ly, lz, xy, xz, yz; //xy, xz, yz are tilt factors 
              lx = plmdchunk->get_box_lx();
              ly = plmdchunk->get_box_ly();
              lz = plmdchunk->get_box_lz();
              xy = plmdchunk->get_box_xy(); // 0 for orthorhombic
              xz = plmdchunk->get_box_xz(); // 0 for orthorhombic
              yz = plmdchunk->get_box_yz(); // 0 for orthorhombic
              int step = plmdchunk->get_timestep();
              
              retrieve_ptr->analyze_frame(types_vector,
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
        DurationMilli retrieve_time_ms = timeNow()-t_start_retrieve;
        total_retrieve_time_ms += retrieve_time_ms.count();
      }
      current_chunk_id++;
    }
    else
      log.printf(" ERROR: Unknown TARGET specified in DispatchAtoms Action.\n");
  }
  else
      log.printf(" ERROR: Unknown TARGET specified in DispatchAtoms Action.\n");

  plumed_time_ms = timeNow()-t_start;
  total_plumed_time_ms += plumed_time_ms.count();

  //high_resolution_clock::time_point p = timeNow();
  //milliseconds ms = duration_cast<milliseconds>(p.time_since_epoch());
  //std::size_t fractional_seconds = ms.count() % 1000;
  //seconds s = duration_cast<seconds>(ms);
  //std::time_t t = s.count();
  //char mbstr[100];
  //if (std::strftime(mbstr, sizeof(mbstr), "%A %c", std::localtime(&t))==0)
  //{
  //    std::cout << "ERROR while trying to convert time to string in DispatchAtoms" << '\n';
  //}
  //std::string timeval = std::string(mbstr)+ ","+std::to_string(fractional_seconds);
  //ds_write_time_stamps.push_back(timeval);

  if (step == total_steps)
  {
      log.printf("total_dispatch_action_time_ms : %f\n",total_plumed_time_ms);
      log.printf("total_simulation_time_ms : %f\n",total_simulation_time_ms);
      log.printf("total_retriever_time_ms : %f\n",total_retrieve_time_ms);
      log.printf("total_read_plumed_data_time_ms : %f\n",total_read_plumed_data_time_ms);
      log.printf("total_time_steps : %d\n",total_steps);
      //std::ofstream outFile("ds_write_time_stamps.txt");
      //for (const auto &e : ds_write_time_stamps) outFile << e << "\n";
  }
}

DispatchAtoms::~DispatchAtoms() {
  dataspaces_writer_ptr = NULL;
  delete dataspaces_writer_ptr;
  dataspaces_reader_ptr = NULL;
  delete dataspaces_reader_ptr;
  retrieve_ptr = NULL;
  delete retrieve_ptr;

}
}
}
