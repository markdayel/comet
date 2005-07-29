#ifndef rotationmatrix_H
#define rotationmatrix_H


class rotationmatrix
{
public:
	enum axis {
		xaxis,
		yaxis,
		zaxis };

	rotationmatrix(void);
	rotationmatrix(MYDOUBLE angle, axis a);

	~rotationmatrix(void);

	MYDOUBLE xx,xy,xz;
	MYDOUBLE yx,yy,yz;
	MYDOUBLE zx,zy,zz;


	void rotatematrix(MYDOUBLE angle, axis a);
	void rotatematrix(MYDOUBLE x_angle, MYDOUBLE y_angle, MYDOUBLE z_angle);

	//void rotate(MYDOUBLE &x, MYDOUBLE &y, MYDOUBLE &z);

	//void rotate(vect& v);

	inline void rotate(MYDOUBLE &x, MYDOUBLE &y, MYDOUBLE &z)
	{
		static MYDOUBLE tmp_x, tmp_y;

		tmp_x = x*xx + y*yx + z*zx;
		tmp_y = x*xy + y*yy + z*zy;
			z = x*xz + y*yz + z*zz;

		x = tmp_x;
		y = tmp_y;

	}

	inline void rotate(vect& v)
	{
		rotate(v.x,v.y,v.z);
	}

};

#endif