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

#include "stdafx.h"
#include "nucleator.h"

nucleator::nucleator(void)
{
	
	radius = RADIUS;
	segment = SEGMENT;
	//P_NUC = (MYDOUBLE) 0.8 / (4*PI*radius*radius);
	geometry = sphere;
	position.x=position.y=position.z=0;
	surf_area = 4 * PI * radius * radius;
	movability = (MYDOUBLE) 0.05 * NODE_INCOMPRESSIBLE_RADIUS / radius;  // inverse ratio of 'viscosities' of node and bead

	radial_rep_distrib_x.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_y.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_z.reserve(RADIAL_SEGMENTS);

}

nucleator::~nucleator(void)
{
}

nucleator::nucleator(shape set_geometry, actin *actinptr)
{
	
	radius = RADIUS;
	segment = SEGMENT;
	//P_NUC = (MYDOUBLE)0.8 / (4*PI*radius*radius);
	geometry = set_geometry;

	if (geometry==sphere)
	{
		surf_area = 4 * PI * radius * radius;
		movability = (MYDOUBLE) 0.05 * NODE_INCOMPRESSIBLE_RADIUS / radius;
	}
	else
	{
		surf_area = 4 * PI * radius * radius  +  2 * PI * radius * segment;
		movability =(MYDOUBLE) 0.05 * NODE_INCOMPRESSIBLE_RADIUS / (radius + segment/4);   // what should this be??
	}
	ptheactin = actinptr;
	ptheactin->nucleation_object = this;
	position.x=position.y=position.z=0;
	definenucleatorgrid();

	radial_rep_distrib_x.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_y.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_z.reserve(RADIAL_SEGMENTS);
}

int nucleator::addnodes(void)
{
	int nodesadded = 0;
	switch (geometry)
	{
	case (sphere):
		{
		nodesadded = addnodessphere();
		break;
		}
	case (capsule):
		{
		nodesadded = addnodescapsule();
		break;
		}
	}
	return nodesadded;
}

int nucleator::addnodessphere(void)
{
	int nodesadded = 0;
	MYDOUBLE x,y,z,r,theta;
	
	MYDOUBLE floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;

	if (( floatingnodestoadd - nodestoadd ) > ( (MYDOUBLE) rand()) / (MYDOUBLE)RAND_MAX)
		nodestoadd++;

	for (int i=0; i< nodestoadd; i++)
	{
			z = (((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX/2)) - 1 ;		// random number -1 to 1
			theta = (2 * PI * ((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX));  // circle vector
			
			if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
			{
				r = radius * mysqrt(1 - z*z);		// radius of circle
			}
			else
			{
				r = radius;  
			}

			x =  r * cos(theta) * (MYDOUBLE)1.001;				// x and y of point
			y =  r * sin(theta) * (MYDOUBLE)1.001;
			z*=  radius * (MYDOUBLE) 1.001;						// z just scaled by radius

			if (ASYMMETRIC_NUCLEATION!=0)
			{
				if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation abouve z=0
					if (z<0) continue;
				if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
					if (z < (radius) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/2) - 1))
						continue;
				if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
					if (z < (radius) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/4) - 3))
						continue;
			}

			ptheactin->node[ptheactin->highestnodecount++].polymerize(x,y,z);
			nodesadded++;
	}
	return nodesadded;
}

int nucleator::addnodescapsule(void)
{
	int nodesadded = 0;
	MYDOUBLE x,y,z,r,theta;
	bool onseg;

	MYDOUBLE rad = radius * (MYDOUBLE) 1.001;

	MYDOUBLE floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;

	if (( floatingnodestoadd - nodestoadd ) > ( (MYDOUBLE) rand()) / (MYDOUBLE)RAND_MAX)
		nodestoadd++;

	for (int i=0; i< nodestoadd; i++)
	{
		//	Pick a random point on capsule:

		// first choose whether on sphere or capsule:
		// p(sphere) = Area of sphere/area of capsule
		//           = 4.PI.r^2 / (4.PI.r^2 + 2.PI.r.h)
		//           = r/(r+2h)


		//onseg = (((2 * rad)/(2 * rad + 3 * segment)) < (((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX)));

		z = (((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX/2)) - 1 ;  // random number -1 to 1
		theta = (2 * PI * ((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX));  // circle vector
		
		onseg = ( (segment /(2*radius+segment)) > ( ((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX) ) ); // on ends or on segment?
		
		if (onseg)
		{
			r = rad;
			z*= segment/2;
		}
		else
		{
			if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
			{
				r = rad * mysqrt(1 - z*z);		// radius of circle
			}
			else
			{
				r = rad;  
			}

			z*=rad;
						
			if (z>0)
				z+= (segment/2); 
			else
				z-= (segment/2);
		}

		x =  r * cos(theta); // x and y of point
		y =  r * sin(theta);
		
		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation abouve z=0
				if (z<0) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (segment/2 + rad) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/2) - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (segment/2 + rad) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/4) - 3))
					continue;
		}

		ptheactin->node[ptheactin->highestnodecount++].polymerize(x,y,z);
		nodesadded++;

	}
	return nodesadded;
}

int nucleator::save(ofstream *outputstream) 
{

MYDOUBLE x,y,z,r;

if (geometry==sphere)
{

 // sphere

	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<=1; z1+= RAD_INCOMP/10)
		{
			r = RAD_INCOMP * mysqrt(1 - z1*z1);		// radius of circle
			
			x = r * cos(theta);				// x and y of point
			y = r * sin(theta);
			z = z1 * RAD_INCOMP;						// z just scaled by radius

			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
			            << 1 << "," << 1 << "," 
						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
						<< 0 << "," << 0 << "," << 0
						<< endl;
		}
}
else
{
	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<=1; z1+= RAD_INCOMP/10)
		{
			r = RAD_INCOMP * mysqrt(1 - z1*z1);		// radius of circle
			
			x = r * cos(theta);				// x and y of point
			y = r * sin(theta);
			z = z1 * RAD_INCOMP;	// z just scaled by radius

			if (z>0)
				z+= (segment/2); 
			else
				z-= (segment/2);

			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
			            << 1 << "," << 1 << "," 
						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
						<< 0 << "," << 0 << "," << 0
						<< endl;

		}
		
	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<=1; z1+= RAD_INCOMP/10)
		{
					
			x = RAD_INCOMP * cos(theta);				// x and y of point
			y = RAD_INCOMP * sin(theta);
			z = z1 * (segment/2);						// z just scaled by segment


			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
			            << 1 << "," << 1 << "," 
						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
						<< 0 << "," << 0 << "," << 0
						<< endl;

		}

}
	return 0;
}

int nucleator::savevrml(ofstream *outputstream) 
{

MYDOUBLE x,y,z,r;

int numpoints = 0;

if (geometry==sphere)
{

 // sphere

	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<=1; z1+= RAD_INCOMP/10)
		{
			r = RAD_INCOMP * mysqrt(1 - z1*z1);		// radius of circle
			
			x = r * cos(theta);				// x and y of point
			y = r * sin(theta);
			z = z1 * RAD_INCOMP;						// z just scaled by radius

			if (numpoints==0)
				*outputstream	<< x << "," << y << "," << z;
			else
				*outputstream	<< "," << x << "," << y << "," << z;

			numpoints++;
		}
}
else
{
	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<=1; z1+= RAD_INCOMP/10)
		{
			r = RAD_INCOMP * mysqrt(1 - z1*z1);		// radius of circle
			
			x = r * cos(theta);				// x and y of point
			y = r * sin(theta);
			z = z1 * RAD_INCOMP;						// z just scaled by radius

			if (z>0)
				z+= (segment/2); 
			else
				z-= (segment/2);

			if (numpoints==0)
				*outputstream	<< x << "," << y << "," << z;
			else
				*outputstream	<< "," << x << "," << y << "," << z;

			numpoints++;
		}
		
	for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/20)
		for (MYDOUBLE z1=-1; z1<1; z1+= radius/10)
		{
					
			x = RAD_INCOMP * cos(theta);				// x and y of point
			y = RAD_INCOMP * sin(theta);
			z = z1 * (segment/2);						// z just scaled by radius

			*outputstream	<< "," << x << "," << y << "," << z;
			numpoints++;
		}


	}



	*outputstream << "] }" << endl;
	*outputstream << "        color Color { color [ ";

	for (int i=0; i<numpoints; i++)
	{
			if (i==0)
				*outputstream	<< "1 1 1";
			else
				*outputstream	<< ",1 1 1";
	}

	return 0;
}

int nucleator::definenucleatorgrid(void)
{

	int_vect gridpos, lastgridpoint;
	bool newgridpoint;

	lastgridpoint.x = lastgridpoint.y = lastgridpoint.z = -1;

	// traverse cuboid around shape, and add to nucleatorgrid

	MYDOUBLE scalefacrad = ((MYDOUBLE)1 + ((MYDOUBLE)2)*(MYDOUBLE)GRIDRES/RAD_INCOMP);

	cout << "Defining nucleator grid for ";
	if (geometry==capsule)
		cout << "capsule...";
	else
		cout << "sphere...";
	cout.flush();

	for (MYDOUBLE x=-(MYDOUBLE)1.2*(RAD_INCOMP+GRIDRES); x<(MYDOUBLE)1.2*(RAD_INCOMP+GRIDRES); x+=GRIDRES/10)
		for (MYDOUBLE y=-(MYDOUBLE)1.2*(RAD_INCOMP+GRIDRES); y<(MYDOUBLE)1.2*(RAD_INCOMP+GRIDRES); y+=GRIDRES/10)
			for (MYDOUBLE z=-(MYDOUBLE)1.2*(RAD_INCOMP+segment+GRIDRES); z<(MYDOUBLE)1.2*(RAD_INCOMP+segment+GRIDRES);z+=GRIDRES/10)
				{
					gridpos.x = (((int)(x / GRIDRES)) + (GRIDSIZE/2) );
					gridpos.y = (((int)(y / GRIDRES)) + (GRIDSIZE/2) );
					gridpos.z = (((int)(z / GRIDRES)) + (GRIDSIZE/2) );

					if ((gridpos.x==lastgridpoint.x) &&   // is it same point as last one?
						(gridpos.y==lastgridpoint.y) &&
						(gridpos.z==lastgridpoint.z))
						continue;

					lastgridpoint = gridpos;

					if (iswithinnucleator(x/scalefacrad,y/scalefacrad,z/scalefacrad))
					{
						newgridpoint = true;

						for (vector <int_vect>::iterator i=ptheactin->nucleatorgrid.begin(); i<ptheactin->nucleatorgrid.end() ; i++ )
						{	 
							if ((i->x == gridpos.x) &&
								(i->y == gridpos.y) &&
								(i->z == gridpos.z))
							{
								newgridpoint = false;
								break;
							}
						}
	
						if (newgridpoint)
						{
							ptheactin->nucleatorgrid.push_back(gridpos);
						}
					}

					
				}

	cout << "done." << endl;
	cout << (int) ptheactin->nucleatorgrid.size() << " nucleator gridpoints created" << endl << endl;

	return 0;
}


bool nucleator::iswithinnucleator(MYDOUBLE x, MYDOUBLE y, MYDOUBLE z)
{
	switch (geometry)
	{
	case (sphere):
		{
		return (((x*x)+(y*y)+(z*z))<(RAD_INCOMP*RAD_INCOMP));
		break;
		}
	case (capsule):
		{
		if (fabs(z)<(segment/2))
			return (((x*x)+(y*y))<RAD_INCOMP*RAD_INCOMP); // cylinder segment
		else
		{
			return (((x*x)+(y*y)+((fabs(z)-(segment/2))*(fabs(z)-(segment/2))))<RAD_INCOMP*RAD_INCOMP);
		}
		break;
		}
	}

	return false;
}

int nucleator::collision(MYDOUBLE &x, MYDOUBLE &y, MYDOUBLE &z)
{
	// node has entered nucleator,
	// return co-ords of node pushed to surface...

	// MYDOUBLE r2,theta,phi;
	MYDOUBLE z2, r,scale;
	MYDOUBLE oldx,oldy,oldz;

	oldx = x;
	oldy = y;
	oldz = z;

	MYDOUBLE rad = RAD_INCOMP * (MYDOUBLE) 1.001; // needed to prevent rounding errors putting back inside nuclator

	switch (geometry)
	{
	case (sphere):
		{

		r = calcdist(x,y,z);

		scale = rad / r;
		x*=scale;
		y*=scale;
		z*=scale;

		radial_rep_distrib_z[(int)(((atan2(y,x)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(z))*(rad-r);
		radial_rep_distrib_x[(int)(((atan2(y,z)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(x))*(rad-r);
		radial_rep_distrib_y[(int)(((atan2(z,x)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(y))*(rad-r);

		break;
		}

	case (capsule):
		{
		if ((fabs(z)<(segment/2)))
			{ 
				// on the cylinder

				r = calcdist(x,y);
				scale = rad / r;
				x*=scale;
				y*=scale;
			}
			else
			{	// on the ends
				
				if (z<0)  // make into a sphere again
					z2 = z + (segment/2);
				else
					z2 = z - (segment/2);

				// calculate theta, phi :

				r = calcdist(x,y,z2);
				scale = rad / r;
				x  *= scale;
				y  *= scale;
				z2 *= scale;

				if (z<0)  
					z = z2 - (segment/2);
				else
					z = z2 + (segment/2);
			
			}

			break;
		}

	}

if (((x-oldx) > NODE_INCOMPRESSIBLE_RADIUS) ||
	((y-oldy) > NODE_INCOMPRESSIBLE_RADIUS) ||
	((z-oldz) > NODE_INCOMPRESSIBLE_RADIUS))
{
	cout << "node nucleus ejection too great " <<  endl;
	return 1;
}


position.x -= (x - oldx) * movability;  // move the nucleator
position.y -= (y - oldy) * movability;
position.z -= (z - oldz) * movability;

x -= (x - oldx) * movability;	// and move node by same amount so that
y -= (y - oldy) * movability;	// still on surface
z -= (z - oldz) * movability;

//x = oldx + ((x - oldx) * (1-movability));  // move the nodes so that
//y = oldy + ((y - oldy) * (1-movability));  // exactly on nucleator surface
//z = oldz + ((z - oldz) * (1-movability));


	return 0;
}

int nucleator::saveradialsegments(ofstream *outputstream) 
{

	for (int i=0; i<RADIAL_SEGMENTS; i++)
	{
		if (i>0)
			*outputstream << ",";
		*outputstream << radial_rep_distrib_x[i];
	}

	for (int i=0; i<RADIAL_SEGMENTS; i++)
	{
		if (i>0)
			*outputstream << ",";
		*outputstream << radial_rep_distrib_y[i];
	}

	for (int i=0; i<RADIAL_SEGMENTS; i++)
	{
		if (i>0)
			*outputstream << ",";
		*outputstream << radial_rep_distrib_z[i];
	}

	return 0;
}

int nucleator::clearradialsegments()
{
	for (int i=0; i<RADIAL_SEGMENTS; i++)
	{
		radial_rep_distrib_x[i]=0;
		radial_rep_distrib_y[i]=0;
		radial_rep_distrib_z[i]=0;
	}
	return 0;
}
