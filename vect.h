/*
    comet - an actin-based motility simulator
    Copyright (C) 2005 Mark J Dayel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#ifndef vect_H
#define vect_H

//#include "stdafx.h"


                                                                                 

class vect
{
 public:

	double x;
    double y;
    double z;
    
    inline vect(void)
	{
	    x=y=z=0.0;
	}
    
    inline vect(const double &a, const  double &b, const  double &c) 
	{
	    x=a;
	    y=b;
	    z=c;
	}
    
    virtual ~vect(void)
	{
	}
    
    inline void zero()
	{
	    x=y=z=0;
	    return;
	}


    inline vect operator-(const vect &param) const
	{
        return vect (*this) -= param;
	}

    inline vect operator+(const vect &param) const
	{
        return vect (*this) += param;
	}

    inline vect operator*(const double &param) const
	{
        return vect (*this) *= param;
	}

    inline vect operator/(const double &param) const
	{
        return vect (*this) /= param;
	}


    inline vect& operator+=(const vect &param)
	{
	    x += param.x;
	    y += param.y;
	    z += param.z;
	    
	    return *this;
	}
    
    inline vect& operator-=(const vect &param)
	{
	    x -= param.x;
	    y -= param.y;
	    z -= param.z;
	    
	    return *this;
	}
    
    inline vect& operator*=(const double &scale)
	{
	    x *= scale;
	    y *= scale;
	    z *= scale;
	    
	    return *this;
	}

    inline vect& operator/=(const double &scale)
	{   
        double tmpscale = 1/scale;
	    x *= tmpscale;
	    y *= tmpscale;
	    z *= tmpscale;
	    
	    return *this;
	}

   
    
    inline double dot(const vect &v) const
	{
	    return ((x * v.x) + (y * v.y) + (z * v.z));
	}
    
    inline double length() const
	{
	    return calcdist(x,y,z);
	}



    inline double recip_length() const
	{
	    return recipcalcdist(x,y,z);
	}
    
    
    const inline vect cross(const vect &v) const
	{
	    vect temp;
	    
	    temp.x = y*v.z - z*v.y;
	    temp.y = z*v.x - x*v.z;
	    temp.z = x*v.y - y*v.x;
	    
	    return temp;
	}
    
    const inline vect unitvec() const
	{
	    vect temp;
	    double recip_len = recipcalcdist(x,y,z);
	    
	    temp.x = x*recip_len;
	    temp.y = y*recip_len;
	    temp.z = z*recip_len;
	    
	    return temp;
	}
    

/*    
      void operator <<(ostream& out)
      {
      
      out << x << "," << y << "," << z;
      
      }
*/  
    inline double sqrlength(void) const
	{
	    return x*x + y*y + z*z;
	}

    // FIXME: Nasty hack here, these are friends to allow placement within the class, not
    // for encapsulation (encapsulation!).  Defining these outside the class leads to linkage problems
    // because of the dependancies between stdafx/vect and the client code. (ML)
    friend ostream &operator<<(ostream &stm, vect const &v){
	stm << v.x << " " << v.y << " " << v.z;
	return stm;
    }
    
    friend istream &operator>>(istream &stm, vect &v){
	stm >> v.x >> v.y >> v.z;
	return stm;
    }
    
};
#endif
