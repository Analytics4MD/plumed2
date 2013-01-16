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
#include "File.h"
#include "Exception.h"
#include "core/Action.h"
#include "core/PlumedMain.h"
#include "core/Value.h"
#include "Communicator.h"
#include "Tools.h"
#include <cstdarg>
#include <cstring>

#include <iostream>
#include <string>

namespace PLMD{

size_t IFile::llread(char*ptr,size_t s){
  plumed_assert(fp);
  size_t r;
  r=fread(ptr,1,s,fp);
  if(feof(fp))   eof=true;
  if(ferror(fp)) err=true;
  return r;
}

IFile& IFile::advanceField(){
  plumed_assert(!inMiddleOfField);
  std::string line;
  bool done=false;
  while(!done){
    getline(line);
    if(!*this){return *this;}
    std::vector<std::string> words=Tools::getWords(line);
    if(words.size()>=2 && words[0]=="#!" && words[1]=="FIELDS"){
      fields.clear();
      for(unsigned i=2;i<words.size();i++){
        Field field;
        field.name=words[i];
        fields.push_back(field);
      }
    } else if(words.size()==4 && words[0]=="#!" && words[1]=="SET"){
      Field field;
      field.name=words[2];
      field.value=words[3];
      field.constant=true;
      fields.push_back(field);
    } else {
      unsigned nf=0;
      for(unsigned i=0;i<fields.size();i++) if(!fields[i].constant) nf++;
      Tools::trimComments(line);
      words=Tools::getWords(line);
      if( words.size()==nf ){
          unsigned j=0;
          for(unsigned i=0;i<fields.size();i++){
            if(fields[i].constant) continue;
            fields[i].value=words[j];
            fields[i].read=false;
            j++;
          }
          done=true;
      } else if( words.size()!=0 ) {
          plumed_merror("mismatch between number of fields in file and expected number");
      }
    }
  }
  inMiddleOfField=true;
  return *this;
}

IFile& IFile::open(const std::string&name){
  FileBase::open(name,"r");
  return *this;
}

IFile& IFile::scanFieldList(std::vector<std::string>&s){
  if(!inMiddleOfField) advanceField();
  if(!*this) return *this;
  s.clear();
  for(unsigned i=0;i<fields.size();i++)
    s.push_back(fields[i].name);
  return *this;
}

bool IFile::FieldExist(const std::string& s){
     std::vector<std::string> slist;
     scanFieldList(slist);
     int mycount = (int) std::count(slist.begin(), slist.end(), s);
     if(mycount>0) return true;
     else return false;
}

IFile& IFile::scanField(const std::string&name,std::string&str){
  if(!inMiddleOfField) advanceField();
  if(!*this) return *this;
  unsigned i=findField(name);
  str=fields[i].value;
  fields[i].read=true;
  return *this;
}

IFile& IFile::scanField(const std::string&name,double &x){
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(const std::string&name,int &x){
  std::string str;
  scanField(name,str);
  if(*this) Tools::convert(str,x);
  return *this;
}

IFile& IFile::scanField(Value* val){
  double ff; scanField(  val->getName(), ff );
  val->set( ff );
  if( FieldExist("min_" + val->getName() ) ){ 
      std::string min, max;
      scanField("min_" + val->getName(), min );
      scanField("max_" + val->getName(), max );
      val->setDomain( min, max ); 
  } else {
      val->setNotPeriodic();
  }
  return *this;
}

IFile& IFile::scanField(){
  if(!ignoreFields){
     for(unsigned i=0;i<fields.size();i++){
       plumed_massert(fields[i].read,"field "+fields[i].name+" was not read: all the fields need to be read otherwise you could miss important infos" );
     }
  }
  inMiddleOfField=false;
  return *this;
}

IFile::IFile():
  inMiddleOfField(false),
  ignoreFields(false)
{
}

IFile::~IFile(){
  if(inMiddleOfField) std::cerr<<"WARNING: IFile closed in the middle of reading. seems strange!\n";
}

IFile& IFile::getline(std::string &str){
  char tmp;
  str="";
  fpos_t pos;
  fgetpos(fp,&pos);
  while(llread(&tmp,1)==1 && tmp && tmp!='\n' && !eof && !err){
    str+=tmp;
  }
  if(err || eof || tmp!='\n'){
    eof = true;
    str="";
    fsetpos(fp,&pos);
  }
  return *this;
}

unsigned IFile::findField(const std::string&name)const{
  unsigned i;
  for(i=0;i<fields.size();i++) if(fields[i].name==name) break;
  if(i>=fields.size()) plumed_merror(name);
  return i;
}

void IFile::reset(bool reset){
 eof = reset;
 err = reset;
 if(!reset) clearerr(fp);
 return;
} 

void IFile::allowIgnoredFields(){
  ignoreFields=true;
}

bool IFile::findCvsAndPeriodic(std::string filename, std::vector<std::string> &cvs,std::vector<std::string> &pmin,std::vector<std::string> &pmax, bool &multivariate){
       std::vector<std::string> fields;
       if(FileExist(filename)){
          cvs.clear(); pmin.clear(); pmax.clear(); 
          open(filename);
          scanFieldList(fields);
          size_t founds,foundm,foundp;
          bool before_sigma=true;
          for(int i=0;i<fields.size();i++){
              size_t pos = 0;
              //found=(fields[i].find("sigma_", pos) || fields[i].find("min_", pos) || fields[i].find("max_", pos) ) ;
              founds=fields[i].find("sigma_", pos)  ;
              foundm=fields[i].find("min_", pos)  ;
              foundp=fields[i].find("max_", pos)  ;
              if (founds!=std::string::npos || foundm!=std::string::npos ||  foundp!=std::string::npos )before_sigma=false;
              // cvs are after time and before sigmas 
              size_t  found; 
              found=fields[i].find("time", pos); 
              if( found==std::string::npos && before_sigma){
                   cvs.push_back(fields[i]);
                   std::cerr<<"found variable number  "<<cvs.size()<<" :  "<<cvs.back()<<std::endl;
                   // get periodicity
                   pmin.push_back("none");
                   pmax.push_back("none");
                   if(FieldExist("min_"+cvs.back())){
              		std::string val;
              		scanField("min_"+cvs.back(),val);
                        pmin[pmin.size()-1]=val; 
                       // std::cerr<<"found min   :  "<<pmin.back()<<std::endl;
                   }
     	           if(FieldExist("max_"+cvs.back())){
              		std::string val;
              		scanField("max_"+cvs.back(),val);
                        pmax[pmax.size()-1]=val; 
                       // std::cerr<<"found max   :  "<<pmax.back()<<std::endl;
                   }
              }
          }
          // is multivariate ???
          std::string sss;
          multivariate=false;
          if(FieldExist("multivariate")){;
         	 scanField("multivariate",sss);
         	 if(sss=="true"){ multivariate=true;}
         	 else if(sss=="false"){ multivariate=false;}
          }
          close();
	  return true;
       }else { 
	return false;
       }
}

}
