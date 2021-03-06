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

    When using the code in academic work please cite:

      Dayel MJ, Akin O, Landeryou M, Risca V, Mogilner A, et al. (2009) 
      In Silico Reconstitution of Actin-Based Symmetry Breaking and Motility. 
      PLoS Biol 7(9):e1000201. doi:10.1371/journal.pbio.1000201

    and include any modifications to the source code with the published work.

*/


#ifndef rotationmatrix_H
#define rotationmatrix_H

#include "vect.h"

class rotationmatrix
{
    
 public:
 //   enum axis {
	//xaxis,
	//yaxis,
	//zaxis };
    
    rotationmatrix(void);
    rotationmatrix(const double &angle, const projection &a);
    
    ~rotationmatrix(void);

	void settoidentity(void);
    
    double xx,xy,xz;
    double yx,yy,yz;
    double zx,zy,zz;
    
    void rotatematrix(const rotationmatrix& rotmatrix);
    void rotatematrix(const double &angle, const projection &a);
    void rotatematrix(const double &x_angle, const double &y_angle, const double &z_angle);
    void rotatematrixrevorder(const double &x_angle, const double &y_angle, const double &z_angle);
    //void rotatematrix(const vect & v);
    rotationmatrix inverse();

    rotationmatrix& operator*=(const rotationmatrix& rotmatrix);
    
    //void rotate(double &x, double &y, double &z);
    
    //void rotate(vect& v);
    
    inline void rotate(double &x, double &y, double &z) const
	{
	    //double tmp_x, tmp_y;
	    
	    double tmp_x = x*xx + y*yx + z*zx;
	    double tmp_y = x*xy + y*yy + z*zy;
	               z = x*xz + y*yz + z*zz;
	    
	    x = tmp_x;
	    y = tmp_y;
	    
	    return;
	    
	}

	inline void rotate(float &x, float &y, float &z) const
	{
	    //double tmp_x, tmp_y;
	    
	    float tmp_x = x*(float)xx + y*(float)yx + z*(float)zx;
	    float tmp_y = x*(float)xy + y*(float)yy + z*(float)zy;
	              z = x*(float)xz + y*(float)yz + z*(float)zz;
	    
	    x = tmp_x;
	    y = tmp_y;
	    
	    return;
	    
	}
    
    inline void rotate(vect& v) const
	{
		rotate(v.x,v.y,v.z);
	    return;
	}

    inline void rotate(vect* v) const
	{
		rotate(v->x,v->y,v->z);
	    return;
	}

    // these are for rotation the vtk coordinates which are in arrays

    inline void rotate(float f[])
    {
        rotate( f[0], f[1], f[2] );
    }

    inline void rotate(double f[])
    {
        rotate( f[0], f[1], f[2] );
    }

    inline rotationmatrix operator*(const double &scale) const
	{
	    rotationmatrix tmp;
	    
	    tmp.xx = xx * scale;
	    tmp.xy = xy * scale;
	    tmp.xz = xz * scale;

        tmp.yx = yx * scale;
	    tmp.yy = yy * scale;
	    tmp.yz = yz * scale;

        tmp.zx = zx * scale;
	    tmp.zy = zy * scale;
	    tmp.zz = zz * scale;
	    
	    return tmp;
	}
    
    void getangles(double& x_angle, double& y_angle, double& z_angle) const; 
    void getaxis_and_angle(double& x, double& y, double& z, double& angle) const;
};

// hmm, we don't need to be a friend (as data is public) but bear it in mind if we 
// intoduce encapsulation later on :)
ostream& operator<<(ostream &ostm, rotationmatrix const &rm);
istream& operator>>(istream &istm, rotationmatrix &rm);

#endif
