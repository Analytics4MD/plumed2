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
#ifndef __RANDOM_H
#define __RANDOM_H

#include <string>
#include <iosfwd>

namespace PLMD{

/// \ingroup TOOLBOX
class Random{
	static const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
	static const int NDIV=(1+(IM-1)/NTAB);
	static const double EPS;
	static const double AM;
	static const double RNMX;
	static const double fact;
	bool incPrec;
        bool switchGaussian;
        double saveGaussian;
	int iy;
	int iv[NTAB];
	int idum;
	std::string name;
public:
	Random(const std::string & name="");
	void setSeed(int idum);
	double RandU01();
	double U01();
	double U01d();
	void WriteStateFull(std::ostream &)const;
	void ReadStateFull (std::istream &);
	void fromString(const std::string & str);
	void toString(std::string & str)const;
	friend std::ostream & operator<<(std::ostream & out,const Random & rng){
		rng.WriteStateFull(out); return out;
	};
	friend std::istream & operator>>(std::istream & in,Random & rng){
		rng.ReadStateFull(in); return in;
	};
	double Gaussian();
	void IncreasedPrecis(bool i){incPrec=i;};
};

}

#endif

