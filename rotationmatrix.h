#ifndef rotationmatrix_H
#define rotationmatrix_H

#include "vect.h"

class rotationmatrix
{
    
 public:
    enum axis {
	xaxis,
	yaxis,
	zaxis };
    
    rotationmatrix(void);
    rotationmatrix(double angle, axis a);
    
    ~rotationmatrix(void);

	void settoidentity(void);
    
    double xx,xy,xz;
    double yx,yy,yz;
    double zx,zy,zz;
    
    
    void rotatematrix(double angle, axis a);
    void rotatematrix(const double &x_angle, const double &y_angle, const double &z_angle);
    
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
    
 //   inline rotationmatrix rotate(double &x, double &y, double &z) const
	//{
	//    static double tmp_x, tmp_y;
	//    
	//    tmp_x = x*xx + y*yx + z*zx;
	//    tmp_y = x*xy + y*yy + z*zy;
	//        z = x*xz + y*yz + z*zz;
	//    
	//    x = tmp_x;
	//    y = tmp_y;
	//    
	//    return *this;
	//    
	//}
 //   
 //   inline rotationmatrix rotate(vect& v) const
	//{
	//    return rotate(v.x,v.y,v.z);
	//}

    void getangles(double& x_angle, double& y_angle, double& z_angle);    
};

// hmm, we don't need to be a friend (as data is public) but bear it in mind if we 
// intoduce encapsulation later on :)
ostream& operator<<(ostream &ostm, rotationmatrix const &rm);
istream& operator>>(istream &istm, rotationmatrix &rm);

#endif
