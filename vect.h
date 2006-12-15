/*
Copyright (C) 2005 Mark J Dayel

You may not distribute, deploy, provide copies of, rent, loan, lease, 
transfer, or grant any rights in the software or derivative works thereof 
in any form to any person.  Reproduction, adaptation, or translation of 
this program without prior written permission from the author is 
prohibited.  Proofs of infringement of copyright violations include but 
not limited to similar code style and structures, similar layout and 
design, similar algorithm design, and containing parts of the original 
software source code.  Copyright notice must remain intact and cannot be 
removed without prior written permission from the author.
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

    inline vect operator+(const vect &param) const
	{
	    vect tmp;
	    
	    tmp.x = x + param.x;
	    tmp.y = y + param.y;
	    tmp.z = z + param.z;
	    
	    return tmp;
	}

    
    inline vect operator-(const vect &param) const
	{
	    vect tmp;
	    
	    tmp.x = x - param.x;
	    tmp.y = y - param.y;
	    tmp.z = z - param.z;
	    
	    return tmp;
	}
    
    inline vect operator-(void) const
	{
	    vect tmp;
	    
	    tmp.x = -x;
	    tmp.y = -y;
	    tmp.z = -z;
	    
	    return tmp;
	}
    
    
    inline vect operator*(const double &scale) const
	{
	    vect tmp;
	    
	    tmp.x = x * scale;
	    tmp.y = y * scale;
	    tmp.z = z * scale;
	    
	    return tmp;
	}

	inline vect operator/(const double &scale) const
	{
	    vect tmp;
		double recip_scale = 1/scale;
	    
	    tmp.x = x * recip_scale;
	    tmp.y = y * recip_scale;
	    tmp.z = z * recip_scale;
	    
	    return tmp;
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
    
    //inline double dot(vect &a, vect &b)
    //{
    //	return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
    //}
    
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
    
    //inline vect cross(const vect &a, const vect &b)
    //{
    //	vect temp;
    
    //	temp.x = a.y*b.z - a.z*b.y;
    //	temp.y = a.z*b.x - a.x*b.z;
    //	temp.z = a.x*b.y - a.y*b.x;
    
    //	return temp;
    //}
    
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
