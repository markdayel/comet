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

	inline vect(void)
	{
		x=y=z=0.0;
	}

	inline vect(const MYDOUBLE &a, const  MYDOUBLE &b, const  MYDOUBLE &c) 
	{
		x=a;
		y=b;
		z=c;
	}

	virtual ~vect(void)
	{
	}

    MYDOUBLE x;
	MYDOUBLE y;
	MYDOUBLE z;

	inline vect operator+(const vect &param)
	{
		vect tmp;

		tmp.x = x + param.x;
		tmp.y = y + param.y;
		tmp.z = z + param.z;

		return tmp;
	}

	inline vect operator-(const vect &param)
	{
		vect tmp;

		tmp.x = x - param.x;
		tmp.y = y - param.y;
		tmp.z = z - param.z;

		return tmp;
	}

	inline vect operator-(void)
	{
		vect tmp;

		tmp.x = -x;
		tmp.y = -y;
		tmp.z = -z;

		return tmp;
	}


	inline vect operator*(const MYDOUBLE &scale)
	{
		vect tmp;

		tmp.x = x * scale;
		tmp.y = y * scale;
		tmp.z = z * scale;

		return tmp;
	}




	inline vect operator+=(const vect &param)
	{
		x += param.x;
		y += param.y;
		z += param.z;

		return *this;
	}

	inline vect operator-=(const vect &param)
	{
		x -= param.x;
		y -= param.y;
		z -= param.z;

		return *this;
	}

	inline vect operator*=(const MYDOUBLE &scale)
	{
		x *= scale;
		y *= scale;
		z *= scale;

		return *this;
	}

	inline MYDOUBLE dot(vect &a, vect &b)
	{
		return ((a.x * b.x) + (a.y * b.y) + (a.z * b.z));
	}

	inline MYDOUBLE length()
	{
		return calcdist(x,y,z);
	}

	inline vect cross(const vect &a, const vect &b)
	{
		vect temp;

		temp.x = a.y*b.z - a.z-b.y;
		temp.y = a.z*b.x - a.x-b.z;
		temp.z = a.x*b.y - a.y-b.x;

		return temp;
	}

	inline vect zero()
	{
		x=y=z=0;
		return *this;
	}

	void operator <<(ostream& out)
	{

	out << x << "," << y << "," << z;

	}

	inline MYDOUBLE sqrlength(void)
	{
		return x*x + y*y + z*z;
	}
};

#endif
