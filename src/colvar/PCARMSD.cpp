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
#include "core/Atoms.h"
#include "core/PlumedMain.h"
#include "ActionRegister.h"
#include "tools/PDB.h"
#include "tools/RMSD.h"
#include "tools/Tools.h"

using namespace std;

namespace PLMD{
namespace colvar{
   
class PCARMSD : public Colvar {
	
  PLMD::RMSD* rmsd;
  bool squared; 
  std::vector< std::vector<Vector> > eigenvectors;
  std::vector<PDB> pdbv;
  std::vector<string> pca_names;
public:
  PCARMSD(const ActionOptions&);
  ~PCARMSD();
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};


using namespace std;

//+PLUMEDOC DCOLVAR PCARMSD 
/*
Calculate the PCA components for a number of provided eigenvectors and an average structure. Performs optimal alignment at every step. 

\par Examples

\verbatim
PCARMSD AVERAGE=file.pdb EIGENVECTORS=eigenvectors.pdb 
\endverbatim

...

*/
//+ENDPLUMEDOC

PLUMED_REGISTER_ACTION(PCARMSD,"PCARMSD")

void PCARMSD::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("compulsory","AVERAGE","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  keys.add("compulsory","EIGENVECTORS","a file in pdb format containing the reference structure and the atoms involved in the CV.");
  useCustomisableComponents(keys);
  keys.addOutputComponent("err","COMPONENTS","the error component ");
  //keys.add("compulsory","TYPE","SIMPLE","the manner in which RMSD alignment is performed.  Should be OPTIMAL or SIMPLE.");
  //keys.addFlag("SQUARED",false," This should be setted if you want MSD instead of RMSD ");
}

PCARMSD::PCARMSD(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),squared(true)
{
  string f_average;
  parse("AVERAGE",f_average);
  string type;	
  type.assign("OPTIMAL");
  string f_eigenvectors;
  parse("EIGENVECTORS",f_eigenvectors);
  checkRead();

  PDB pdb;

  // read everything in ang and transform to nm if we are not in natural units
  if( !pdb.read(f_average,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength()) )
      error("missing input file " + f_average );


  rmsd = new RMSD();
  bool remove_com=true;
  bool normalize_weights=true;
  // here align and displace are a simple vector of ones
  std::vector<double> align; align=pdb.getOccupancy();for(unsigned i=0;i<align.size();i++){align[i]=1.;} ;
  std::vector<double> displace;  displace=pdb.getBeta();for(unsigned i=0;i<displace.size();i++){displace[i]=1.;} ; 
  log.printf("SIZE OF ALIGN: %d \n",align.size());
  for(unsigned i=0;i<align.size();i++){log.printf("AL %f\n",align[i]);}
  log.printf("SIZE OF DISPLACE: %d \n",align.size());
  for(unsigned i=0;i<displace.size();i++){log.printf("DIS %f\n",displace[i]);}
  // reset again to reimpose unifrom weights (safe to disable this)
  rmsd->set(align,displace,pdb.getPositions(),type,remove_com,normalize_weights);
  requestAtoms( pdb.getAtomNumbers() );

  addComponentWithDerivatives("err"); componentIsNotPeriodic("err"); 

  log.printf("  average from file %s\n",f_average.c_str());
  log.printf("  which contains %d atoms\n",getNumberOfAtoms());
  log.printf("  method for alignment : %s \n",type.c_str() );

  // now get the eigenvectors 
  // open the file
  FILE* fp=fopen(f_eigenvectors.c_str(),"r");
  std::vector<AtomNumber> aaa;
  unsigned neigenvects;
  neigenvects=0;
  if (fp!=NULL)
  {
    log<<"  Opening the eigenvectors file "<<f_eigenvectors.c_str()<<"\n";
    bool do_read=true;
    while (do_read){
         PDB mypdb; 
	 // check the units for reading this file: how can they make sense? 
         do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
         if(do_read){
            unsigned nat=0;
            neigenvects++;
            if(mypdb.getAtomNumbers().size()==0) error("number of atoms in a frame should be more than zero");
            if(nat==0) nat=mypdb.getAtomNumbers().size();
            if(nat!=mypdb.getAtomNumbers().size()) error("frames should have the same number of atoms");
            if(aaa.empty()) aaa=mypdb.getAtomNumbers();
            if(aaa!=mypdb.getAtomNumbers()) error("frames should contain same atoms in same order");
            log<<"  Found eigenvector: "<<neigenvects<<" containing  "<<mypdb.getAtomNumbers().size()<<" atoms\n"; 
	    pdbv.push_back(mypdb); 
            eigenvectors.push_back(mypdb.getPositions()); 		
         }else{break ;}
    }
    fclose (fp);
    log<<"  Found total "<<neigenvects<< " eigenvectors in the file "<<f_eigenvectors.c_str()<<" \n"; 
    if(neigenvects==0) error("at least one eigenvector is expected");
  } 
  // the components 
  for(unsigned i=0;i<neigenvects;i++){
        string xx; Tools::convert(i,xx);
        string name; name=string("pca_")+xx;
	pca_names.push_back(name);
	addComponentWithDerivatives(name.c_str()); componentIsNotPeriodic(name.c_str());	
  }  
  turnOnDerivatives();

}

PCARMSD::~PCARMSD(){
  delete rmsd;
}


// calculator
void PCARMSD::calculate(){
        Tensor rotation,invrotation;
        Matrix<std::vector<Vector> > drotdpos(3,3);
        std::vector<Vector> DDistDRef;
        std::vector<Vector> alignedpos;
        std::vector<Vector> centeredpos;
        std::vector<Vector> centeredref;
        std::vector<Vector> ddistdpos;
        std::vector<Vector> derivatives;
        double r=rmsd->calc_PCAelements( getPositions(), ddistdpos, rotation ,  drotdpos , alignedpos ,centeredpos, centeredref ,squared);
	invrotation=rotation.transpose();
	
	Value* verr=getPntrToComponent("err");
	verr->set(r);
	for(unsigned iat=0;iat<getNumberOfAtoms();iat++){
		        setAtomsDerivatives (verr,iat,ddistdpos[iat]);
	}

	std::vector< Vector > der;
	der.resize(getNumberOfAtoms());


	// (type OPTIMAL with homogeneous weights, normalized weights, squared=true)
	// 1) centeredpos
	// 2) drotdpos   
	// 3) invrotation
	// 4) ddistdpos   OK
	// 5) centeredref OK
	// 6) alignedpos OK
	// 7) err OK
	for(unsigned i=0;i<eigenvectors.size();i++){
		Value* value=getPntrToComponent(pca_names[i].c_str());
		double val;val=0.;
		for(unsigned iat=0;iat<getNumberOfAtoms();iat++){
			val+=dotProduct(alignedpos[iat]-centeredref[iat],eigenvectors[i][iat]);	
		}
		value->set(val);
		// here the loop is reversed to better suit the structure of the derivative of the rotation matrix
		double tmp1;
		der.clear();
		for(unsigned a=0;a<3;a++){
			for(unsigned b=0;b<3;b++){
				for(unsigned iat=0;iat<getNumberOfAtoms();iat++){
					tmp1=0.;
					for(unsigned n=0;n<getNumberOfAtoms();n++){
						tmp1+=centeredpos[n][b]*eigenvectors[i][n][a];
					}
					der[iat]+=drotdpos[a][b][iat]*tmp1;	
		//			log.printf("XXXXX  %d %d %d : %f %f %f\n",a,b,iat,drotdpos[a][b][iat][0],drotdpos[a][b][iat][1],drotdpos[a][b][iat][2]);
				}
			}
		}
		Vector v1;
		for(unsigned n=0;n<getNumberOfAtoms();n++){
				v1+=(1./getNumberOfAtoms())*matmul(invrotation,eigenvectors[i][n]);	
		}	
		for(unsigned iat=0;iat<getNumberOfAtoms();iat++){
			der[iat]+=matmul(invrotation,eigenvectors[i][iat])-v1;	
		        setAtomsDerivatives (value,iat,der[iat]);
		}		
 	 }

  setBoxDerivativesNoPbc();

}

}
}



