/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#ifndef __PLUMED_mapping_Mapping_h
#define __PLUMED_mapping_Mapping_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "vesselbase/ActionWithVessel.h"
#include "FakeFrame.h"
#include <vector>

namespace PLMD {

class PDB;
class ReferenceConfiguration;

namespace mapping {

class Mapping :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue,
  public vesselbase::ActionWithVessel
  {
private:
//  The derivative wrt to the distance from the frame
  std::vector<double> dfframes;
/// These are the configurations that serve as references
  std::vector<ReferenceConfiguration*> frames; 
/// The names of the projection coordinates
  std::vector<std::string> property;
/// The forces on each of the derivatives (used in apply)
  std::vector<double> forcesToApply;
protected:
/// The (transformed) distance from each frame
  std::vector<double> fframes;
/// These are where the reference configurations should be projected
  std::vector< std::vector<double> > low_dim;  // N.B. These are protected for a reason
/// Get the number of frames in the path
  unsigned getNumberOfReferencePoints() const ;
/// Calculate the value of the distance from the ith frame
  double calculateDistanceFunction( const unsigned& ifunc, const bool& squared );
/// Store the distance function
  void storeDistanceFunction( const unsigned& ifunc );
public:
  static void registerKeywords( Keywords& keys );
  Mapping(const ActionOptions&);
  ~Mapping();
/// Overload the virtual functions that appear in both ActionAtomistic and ActionWithArguments
  void calculateNumericalDerivatives( ActionWithValue* a=NULL );
  void lockRequests();
  void unlockRequests();
/// Distance from a point is never periodic
  bool isPeriodic(){ return false; }
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives();  // N.B. This is replacing the virtual function in ActionWithValue
/// Get the value of lambda for paths and property maps 
  virtual double getLambda();
/// This does the transformation of the distance by whatever function is required
  virtual double transformHD( const double& dist, double& df )=0;
/// Get the number of properties we are projecting onto
  unsigned getNumberOfProperties() const ;
/// Get the name of the ith property we are projecting
  std::string getPropertyName( const unsigned& iprop ) const ;
/// Get the index of a particular named property 
  unsigned getPropertyIndex( const std::string& name ) const ;
/// Get the name of the ith argument
  std::string getArgumentName( unsigned& iarg );
/// Get the value of the ith property for the current frame
  double getPropertyValue( const unsigned& iprop ) const ;
/// Return the current value of the high dimensional function
  double getCurrentHighDimFunctionValue( const unsigned& ider ) const ;
/// Perform chain rule for derivatives
  void mergeDerivatives( const unsigned& , const double& );
/// Apply the forces 
  void apply();
};

inline
unsigned Mapping::getNumberOfReferencePoints() const {
  return low_dim.size();
}

inline
unsigned Mapping::getNumberOfDerivatives(){
  unsigned nat=getNumberOfAtoms();
  if(nat>0) return 3*nat + 9 + getNumberOfArguments();
  return getNumberOfArguments();
}

inline
void Mapping::lockRequests(){
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

inline
void Mapping::unlockRequests(){
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

inline
unsigned Mapping::getNumberOfProperties() const {
  return property.size();
}

inline
std::string Mapping::getPropertyName( const unsigned& iprop ) const {
  plumed_dbg_assert( iprop<property.size() );
  return property[iprop];
}

inline
double Mapping::getPropertyValue( const unsigned& iprop ) const {
  plumed_dbg_assert( iprop<property.size() );
  return low_dim[current][iprop];
}

inline 
double Mapping::getCurrentHighDimFunctionValue( const unsigned& ider ) const {
  plumed_dbg_assert( ider<2 );
  return fframes[ider*low_dim.size() + current];
}

inline
void Mapping::storeDistanceFunction( const unsigned& ifunc ){
  plumed_dbg_assert( ifunc<low_dim.size() );
  unsigned storef=low_dim.size()+ifunc;
  fframes[storef]=fframes[ifunc]; dfframes[storef]=dfframes[ifunc];
  frames[storef]->copyDerivatives( frames[ifunc] );
}

}
}
#endif
