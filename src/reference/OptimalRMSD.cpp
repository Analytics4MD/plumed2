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
#include "MetricRegister.h"
#include "RMSDBase.h"
#include "tools/Matrix.h"

namespace PLMD{

class OptimalRMSD : public RMSDBase {
public:
  OptimalRMSD(const ReferenceConfigurationOptions& ro);
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );

  template <bool safe,bool alEqDis>
  double optimalAlignment(const  std::vector<double>  & align,
                          const  std::vector<double>  & displace,
                          const std::vector<Vector> & positions,
                          bool squared=false);
};

PLUMED_REGISTER_METRIC(OptimalRMSD,"OPTIMAL")

OptimalRMSD::OptimalRMSD(const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
RMSDBase(ro)
{
}

void OptimalRMSD::read( const PDB& pdb ){
  readAtomsFromPDB( pdb ); 
}

double OptimalRMSD::calc( const std::vector<Vector>& pos, const bool& squared ){
  if( getAlign()==getDisplace() ) return optimalAlignment<false,true>(getAlign(),getDisplace(),pos,squared); 
  return optimalAlignment<false,false>(getAlign(),getDisplace(),pos,squared);
}

template <bool safe,bool alEqDis>
double OptimalRMSD::optimalAlignment(const  std::vector<double>  & align,
                                     const  std::vector<double>  & displace,
                                     const std::vector<Vector> & positions,
                                     bool squared) {
  double dist(0);
  double norm(0);
  double dnorm(0);
// This is the trace of positions*positions + reference*reference
  double sum00w(0);
  double sum11w(0);
// This is positions*reference
  Tensor sum01w;

  Vector cpositions;
  Vector creference;

  unsigned n=getNumberOfReferencePositions();
  // Number of members of align and displace is the number of reference positions
  plumed_dbg_assert( n==align.size() && n==displace.size() );
  // Positions array might contain vectors that are not particularly interesting
  plumed_dbg_assert( positions.size()==getNumberOfAtoms() );

// first expensive loop: compute centers
  for(unsigned iat=0;iat<n;iat++){
    unsigned iatom=getAtomIndex(iat);
    double w=align[iat]; norm+=w;
    if(!alEqDis) dnorm+=displace[iatom];
    cpositions+=positions[iatom]*w;
    creference+=getReferencePosition(iat)*w;
  }
  double invnorm=1.0/norm;
  double invdnorm;
  if(!alEqDis) invdnorm=1.0/dnorm;

  cpositions*=invnorm;
  creference*=invnorm;

// second expensive loop: compute second moments wrt centers
  for(unsigned iat=0;iat<n;iat++){
    unsigned iatom=getAtomIndex(iat);
    double w=align[iat]; Vector ref=getReferencePosition(iat);
    sum00w+=dotProduct(positions[iatom]-cpositions,positions[iatom]-cpositions)*w;
    sum11w+=dotProduct(ref-creference,ref-creference)*w;
    sum01w+=Tensor(positions[iatom]-cpositions,ref-creference)*w;
  }

  double rr00=sum00w*invnorm;
  Tensor rr01=sum01w*invnorm;
  double rr11=sum11w*invnorm;

  Matrix<double> m=Matrix<double>(4,4);
  m[0][0]=rr00+rr11+2.0*(-rr01[0][0]-rr01[1][1]-rr01[2][2]);
  m[1][1]=rr00+rr11+2.0*(-rr01[0][0]+rr01[1][1]+rr01[2][2]);
  m[2][2]=rr00+rr11+2.0*(+rr01[0][0]-rr01[1][1]+rr01[2][2]);
  m[3][3]=rr00+rr11+2.0*(+rr01[0][0]+rr01[1][1]-rr01[2][2]);
  m[0][1]=2.0*(-rr01[1][2]+rr01[2][1]);
  m[0][2]=2.0*(+rr01[0][2]-rr01[2][0]);
  m[0][3]=2.0*(-rr01[0][1]+rr01[1][0]);
  m[1][2]=2.0*(-rr01[0][1]-rr01[1][0]);
  m[1][3]=2.0*(-rr01[0][2]-rr01[2][0]);
  m[2][3]=2.0*(-rr01[1][2]-rr01[2][1]);
  m[1][0] = m[0][1];
  m[2][0] = m[0][2];
  m[2][1] = m[1][2];
  m[3][0] = m[0][3];
  m[3][1] = m[1][3];
  m[3][2] = m[2][3];

  std::vector<double> eigenvals; Matrix<double> eigenvecs;
  int diagerror=diagMat(m, eigenvals, eigenvecs );

  if (diagerror!=0){
    std::string sdiagerror;
    Tools::convert(diagerror,sdiagerror);
    std::string msg="DIAGONALIZATION FAILED WITH ERROR CODE "+sdiagerror;
    plumed_merror(msg);
  }

  dist=eigenvals[0];

  Matrix<double> ddist_dm(4,4);

  Vector4d q(eigenvecs[0][0],eigenvecs[0][1],eigenvecs[0][2],eigenvecs[0][3]);

// This is the rotation matrix that brings reference to positions
// i.e. matmul(rotation,reference[iat])+shift is fitted to positions[iat]

  Tensor rotation;
  rotation[0][0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotation[1][1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotation[2][2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  rotation[0][1]=2*(+q[0]*q[3]+q[1]*q[2]);
  rotation[0][2]=2*(-q[0]*q[2]+q[1]*q[3]);
  rotation[1][2]=2*(+q[0]*q[1]+q[2]*q[3]);
  rotation[1][0]=2*(-q[0]*q[3]+q[1]*q[2]);
  rotation[2][0]=2*(+q[0]*q[2]+q[1]*q[3]);
  rotation[2][1]=2*(-q[0]*q[1]+q[2]*q[3]);

  double prefactor=2.0*invnorm;
  Vector shift=cpositions-matmul(rotation,creference);

  if(!squared && alEqDis) prefactor*=0.5/sqrt(dist);

// if "safe", recompute dist here to a better accuracy
  if(safe || !alEqDis) dist=0.0;

// If safe is set to "false", MSD is taken from the eigenvalue of the M matrix
// If safe is set to "true", MSD is recomputed from the rotational matrix
// For some reason, this last approach leads to less numerical noise but adds an overhead

// third expensive loop: derivatives
  for(unsigned iat=0;iat<n;iat++){
    unsigned iatom=getAtomIndex(iat);
// there is no need for derivatives of rotation and shift here as it is by construction zero
// (similar to Hellman-Feynman forces)
    Vector d(positions[iatom]-shift - matmul(rotation,getReferencePosition(iat)));
    if(alEqDis){
       addAtomicDerivatives( iat, prefactor*align[iat]*d );
       if(safe) dist+=align[iat]*invnorm*modulo2(d);
    } else {
      dist+=displace[iat]*invdnorm*modulo2(d);
// these are the derivatives assuming the roto-translation as frozen
      addAtomicDerivatives( iat, 2*displace[iat]*invdnorm*d );
    }
  }

  if(!alEqDis){
// Here we have to recover the correct derivatives.
// Instead of computing explicitly the derivatives of the rotational matrix
// we enforce that the final result is invariant wrt translation and rotation.
// If we interpret these derivatives as forces, this amounts in
// setting to zero the total force and torque.

// This is the inertia matrix:
    Tensor I;
    for(unsigned iat=0;iat<n;iat++){
      unsigned iatom=getAtomIndex(iat);
      Vector p=positions[iatom]-cpositions;
      I+=align[iat]*(Tensor::identity()*modulo2(p)-Tensor(p,p));
    }
// Total force:
    Vector lin; for(unsigned iat=0;iat<getNumberOfAtoms();iat++) lin+=atom_ders[iat];
// Remove from each atom a force proportional to its weight
    for(unsigned iat=0;iat<n;iat++) addAtomicDerivatives( iat, -lin*align[iat]*invnorm);
// Total torque:
    Vector omega; 
    for(unsigned iat=0;iat<n;iat++){
        unsigned iatom=getAtomIndex(iat);
        omega+=crossProduct(positions[iatom]-cpositions,retrieveAtomicDerivatives(iat) );
    }
    omega=matmul(inverse(I),omega);
// Remove from each atom a torque proportional to its weight
    for(unsigned iat=0;iat<n;iat++){
       unsigned iatom=getAtomIndex(iat);
       addAtomicDerivatives( iat, -crossProduct(omega,positions[iatom]-cpositions)*align[iat]);
    }
  }

  if(!squared){
    dist=sqrt(dist);
    if(!alEqDis){
      double xx=0.5/dist;
      for(unsigned iat=0;iat<atom_ders.size();iat++) atom_ders[iat]*=xx;
    }
  }

  return dist;
}

}
