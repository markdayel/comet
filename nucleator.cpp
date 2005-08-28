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

#ifdef SEED_INSIDE			
double NUCPOINT_SCALE = 0.9;
#else
double NUCPOINT_SCALE = 1.001;
#endif
			
nucleator::nucleator(void)
{  
	// this should never be called anyhow.  Needs pointer to actin

	//radius = RADIUS;
	//segment = 2*CAPSULE_HALF_LINEAR;
	////P_NUC = (double) 0.8 / (4*PI*radius*radius);
	//geometry = sphere;
	//position.zero();
	//surf_area = 4 * PI * radius * radius;
 //   // inverse ratio of 'viscosities' of node and bead
 //   movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / radius;    

	//deltanucposn.zero();

	//torque.zero();

	//centerofmass.zero();

	//momentofinertia.x = 1000 * MofI;  // mark:  todo: calculate these numbers
	//momentofinertia.y = 1000 * MofI;
	//momentofinertia.z = 500 * MofI;

 //   definecagepoints();
	//colour.r = colour.g = colour.b = 1.0;

}

nucleator::~nucleator(void)
{
}

nucleator::nucleator(shape set_geometry, actin *actinptr)
{
	radius = RADIUS;
	//segment = 2*CAPSULE_HALF_LINEAR;
	//P_NUC = (double)0.8 / (4*PI*radius*radius);
	geometry = set_geometry;

	if (geometry==sphere)
	{
		surf_area = 4 * PI * radius * radius;
		movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / radius;
	}
	else
	{
		surf_area = 4 * PI * radius * radius  +  4 * PI * radius * CAPSULE_HALF_LINEAR;
		movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / (radius + CAPSULE_HALF_LINEAR/2);   // what should this be??
	}

	ptheactin = actinptr;
	ptheactin->p_nuc = this;

	position.zero();
	definenucleatorgrid();

	deltanucposn.zero();

	centerofmass.zero();;

	torque.zero();;

	momentofinertia.x = 1000 * MofI;  // mark:  todo: calculate these numbers
	momentofinertia.y = 1000 * MofI;
	momentofinertia.z = 500 * MofI;

	//if( geometry==sphere ){
	//radial_rep_distrib_x.reserve(RADIAL_SEGMENTS);
	//radial_rep_distrib_y.reserve(RADIAL_SEGMENTS);
	//radial_rep_distrib_z.reserve(RADIAL_SEGMENTS);
	//} else {
	//    set_rep_bins();
	//}

	definecagepoints();
	colour.r = colour.g = colour.b = 1.0;

	segs.setupsegments(this, ptheactin);

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
	double x,y,z,r,theta;
	
	double floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;

	if (( floatingnodestoadd - nodestoadd ) * RAND_MAX > rand() )
		nodestoadd++;

	for (int i=0; i< nodestoadd; i++)
	{
		z = (2 * rand() * RECIP_RAND_MAX ) - 1 ;		// random number -1 to 1
		theta = (2 * PI * rand() * RECIP_RAND_MAX);  // circle vector
		
		if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
		{
			r = radius * sqrt(1 - z*z);		// radius of circle
		}
		else
		{
			r = radius;  
		}

		x =  r * cos(theta) * NUCPOINT_SCALE;				// x and y of point
		y =  r * sin(theta) * NUCPOINT_SCALE;
		z*=  radius * NUCPOINT_SCALE;					// z just scaled by radius

		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation above z=0
				if ((y<0) || (fabs(x+z)>0.5)) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (radius) *( 2 * rand() * RECIP_RAND_MAX - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (radius) *( 4 * rand() * RECIP_RAND_MAX - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4) { // fixed random location
			static double fixed_x = x;
			static double fixed_y = y;
			static double fixed_z = z;
			x = fixed_x;
			y = fixed_y;
			z = fixed_z;				
		}
			if (ASYMMETRIC_NUCLEATION==7)  // half caps one side
				if ( ((y>0)&&(z>0)) || ((y<0)&&(z<0)) || (z<0))
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
	double x,y,z,r,theta;
	bool onseg;

	double rad = radius * NUCPOINT_SCALE;

	double floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;

	if (( floatingnodestoadd - nodestoadd ) > ( rand() * RECIP_RAND_MAX ))
		nodestoadd++;

	for (int i=0; i< nodestoadd; i++)
	{
		//	Pick a random point on capsule:

		// first choose whether on sphere or capsule:
		// p(sphere) = Area of sphere/area of capsule
		//           = 4.PI.r^2 / (4.PI.r^2 + 2.PI.r.h)
		//           = r/(r+2h)


		//onseg = (((2 * rad)/(2 * rad + 3 * segment)) < (((double) rand()) / (double)(RAND_MAX)));

		z = ( 2 * rand() * RECIP_RAND_MAX ) - 1 ;  // random number -1 to 1
		theta = (2 * PI * rand() * RECIP_RAND_MAX);  // circle vector
		
		onseg = ( (CAPSULE_HALF_LINEAR /(radius+CAPSULE_HALF_LINEAR)) > ( rand() * RECIP_RAND_MAX )); // on ends or on segment?
		
		if (onseg)
		{
			r = rad;
			z*= CAPSULE_HALF_LINEAR;
		}
		else
		{
			if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
			{
				r = rad * sqrt(1 - z*z);		// radius of circle
			}
			else
			{
				r = rad;  
			}

			z*=rad;
						
			if (z>0)
				z+= (CAPSULE_HALF_LINEAR); 
			else
				z-= (CAPSULE_HALF_LINEAR);
		}

		x =  r * cos(theta); // x and y of point
		y =  r * sin(theta);
		
		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation above z=0
				if (z<0) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 2*  rand() * RECIP_RAND_MAX - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 4 * rand() * RECIP_RAND_MAX - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4)  // caps only
				if (fabs(z) < (CAPSULE_HALF_LINEAR))
					continue;
			if (ASYMMETRIC_NUCLEATION==5)  // half caps only
				if ( (fabs(z) < (CAPSULE_HALF_LINEAR)) || ((x<0)&&(z>0)) || ((x>0)&&(z<0)) )
					continue;
			if (ASYMMETRIC_NUCLEATION==6)  // caps only
				if (fabs(z) < 0.5 * (CAPSULE_HALF_LINEAR + rad) * ( 4 * rand() * RECIP_RAND_MAX - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==7)  // half caps one side
				if ( (fabs(z) < (0.7*CAPSULE_HALF_LINEAR)) || ((x>0)&&(z>0)) || ((x<0)&&(z<0)) || (z>0))
					continue;
			if (ASYMMETRIC_NUCLEATION==8)  //  cap one side
				if ( (fabs(z) < (CAPSULE_HALF_LINEAR)) || (z<0))
					continue;
		}

		//}

		
		ptheactin->node[ptheactin->highestnodecount++].polymerize(x,y,z);
		nodesadded++;

	}
	return nodesadded;
}

//
//int nucleator::save(ofstream *outputstream) 
//{
//
//double x,y,z,r;
//
//if (geometry==sphere)
//{
//
// // sphere
//
//	for (double theta=-PI; theta<PI; theta+=2*PI/20)
//		for (double z1=-1; z1<=1; z1+= RAD_INCOMP/10)
//		{
//			r = RAD_INCOMP * sqrt(1 - z1*z1);		// radius of circle
//			
//			x = r * cos(theta);				// x and y of point
//			y = r * sin(theta);
//			z = z1 * RAD_INCOMP;						// z just scaled by radius
//
//			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
//			            << 1 << "," << 1 << "," 
//						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
//						<< 0 << "," << 0 << "," << 0
//						<< endl;
//		}
//}
//else
//{
//	for (double theta=-PI; theta<PI; theta+=2*PI/20)
//		for (double z1=-1; z1<=1; z1+= RAD_INCOMP/10)
//		{
//			r = RAD_INCOMP * sqrt(1 - z1*z1);		// radius of circle
//			
//			x = r * cos(theta);				// x and y of point
//			y = r * sin(theta);
//			z = z1 * RAD_INCOMP;	// z just scaled by radius
//
//			if (z>0)
//				z+= (CAPSULE_HALF_LINEAR); 
//			else
//				z-= (CAPSULE_HALF_LINEAR);
//
//			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
//			            << 1 << "," << 1 << "," 
//						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
//						<< 0 << "," << 0 << "," << 0
//						<< endl;
//
//		}
//		
//	for (double theta=-PI; theta<PI; theta+=2*PI/20)
//		for (double z1=-1; z1<=1; z1+= RAD_INCOMP/10)
//		{
//					
//			x = RAD_INCOMP * cos(theta);				// x and y of point
//			y = RAD_INCOMP * sin(theta);
//			z = z1 * CAPSULE_HALF_LINEAR;						// z just scaled by segment
//
//
//			*outputstream	<< 0 << "," << x << "," << y << "," << z << "," 
//			            << 1 << "," << 1 << "," 
//						<< 1 << "," << 1 << "," << 1<< "," << 0 << "," << 0 << ","
//						<< 0 << "," << 0 << "," << 0
//						<< endl;
//
//		}
//
//}
//	return 0;
//}

int nucleator::savevrml(ofstream *outputstream) 
{

int numpoints = 0;

for (vector <vect>::iterator point=cagepoints.begin(); 
	      point<cagepoints.end() ; ++point )
	{
		if (numpoints==0)
			*outputstream	<< point->x << "," << point->y << "," << point->z;
		else
			*outputstream	<< "," << point->x << "," << point->y << "," << point->z;

		numpoints++;
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

	double scalefacrad = ((double)1 + ((double)2)*(double)GRIDRES/RAD_INCOMP);

	cout << "Defining nucleator grid for ";
	if (geometry==capsule)
		cout << "capsule...";
	else
		cout << "sphere...";
	cout.flush();

	for (double x=-(double)1.2*(RAD_INCOMP+GRIDRES); x<(double)1.2*(RAD_INCOMP+GRIDRES); x+=GRIDRES/10)
		for (double y=-(double)1.2*(RAD_INCOMP+GRIDRES); y<(double)1.2*(RAD_INCOMP+GRIDRES); y+=GRIDRES/10)
			for (double z=-(double)1.2*(RAD_INCOMP+CAPSULE_HALF_LINEAR*2+GRIDRES); z<(double)1.2*(RAD_INCOMP+CAPSULE_HALF_LINEAR*2+GRIDRES);z+=GRIDRES/10)
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


bool nucleator::iswithinnucleator(const double& x, const double& y, const double& z)
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
		if (fabs(z)<(CAPSULE_HALF_LINEAR))
			return (((x*x)+(y*y))<RAD_INCOMP*RAD_INCOMP); // cylinder segment
		else
		{
			return (((x*x)+(y*y)+((fabs(z)-(CAPSULE_HALF_LINEAR))*(fabs(z)-(CAPSULE_HALF_LINEAR))))<RAD_INCOMP*RAD_INCOMP);
		}
		break;
		}
	}

	return false;
}


bool nucleator::collision(nodes &node)//(double &x, double &y, double &z)
{  // returns true if succeeds, false if fails due to too great node ejection

	// node has entered nucleator,
	// return co-ords of node pushed to surface...

    double r, scale, z2;
	vect node_disp;
	vect oldpos;

	oldpos = node;

	double rad = RAD_INCOMP * (double) 1.001; // needed to prevent rounding errors putting back inside nuclator

    // FIXME: add no movement of nodes to outside the nucleator when SEED_INSIDE is set (ML)? 
	switch (geometry)
	{
	case (sphere):
		{
		r = node.length();;

		scale = rad / r;

		node *= scale;
		
        node.nucleator_impacts += (rad-r);

		break;
		}

	case (capsule):

		{
		if ((fabs(node.z)<(CAPSULE_HALF_LINEAR)))
			{ 
				// on the cylinder

				r = calcdist(node.x,node.y);
				scale = rad / r;
				node.x *= scale;
				node.y *= scale;
				
                node.nucleator_impacts += (rad-r);
			}
			else
			{	// on the ends
				
				if (node.z<0)  // make into a sphere again
					z2 = node.z + (CAPSULE_HALF_LINEAR);
				else
					z2 = node.z - (CAPSULE_HALF_LINEAR);

				// calculate theta, phi :

				r = calcdist(node.x,node.y,z2);
				scale = rad / r;
				node.x *= scale;
				node.y *= scale;
				z2 *= scale;

				if (node.z<0)  
					node.z = z2 - (CAPSULE_HALF_LINEAR);
				else
					node.z = z2 + (CAPSULE_HALF_LINEAR);
			
				node.nucleator_impacts += (rad-r);
			}

			break;
		}

	}

node_disp = node - oldpos;

if ((fabs(node_disp.x) > NODE_INCOMPRESSIBLE_RADIUS) ||
	(fabs(node_disp.y) > NODE_INCOMPRESSIBLE_RADIUS) ||
	(fabs(node_disp.z) > NODE_INCOMPRESSIBLE_RADIUS))
{
	cout << "node nucleus ejection too great " <<  endl;
	cout << "old (x,y,z): " <<  oldpos.x << ", " << oldpos.y << ", " << oldpos.z << endl;
	cout << "new (x,y,z): " <<  node.x << ", " << node.y << ", " << node.z <<  endl;
	return false;  // failed
}

#ifndef SEED_INSIDE

deltanucposn -= node_disp * movability;  // move the nucleator

node -= node_disp * movability;	

#endif

// rotate the nucleator

// calculate torque (as pos x disp)

if (ROTATION)
{
	torque.x += ((oldpos.y-centerofmass.y)*node_disp.z - (oldpos.z-centerofmass.z)*node_disp.y);
	torque.y += ((oldpos.z-centerofmass.z)*node_disp.x - (oldpos.x-centerofmass.x)*node_disp.z);
	torque.z += ((oldpos.x-centerofmass.x)*node_disp.y - (oldpos.y-centerofmass.y)*node_disp.x);
}

	return true; // sucessful node ejection
}


int nucleator::save_data(ofstream &ostr) 
{
    ostr << geometry << endl;
    ostr << radius << endl;
    ostr << CAPSULE_HALF_LINEAR << endl;
    ostr << surf_area << endl;
    ostr << movability << endl;
    ostr << position << endl;
    ostr << deltanucposn << endl;
    ostr << torque << endl;
    ostr << centerofmass << endl;
    ostr << momentofinertia << endl;
    ostr << nucleator_rotation << endl;
    // omitting
    //  colour
    //  cagepoints
    //  radial_rep_distrib_x
    //  radial_rep_distrib_y
    //  radial_rep_distrib_z
    //  nbdy_segs
    //  ncap_segs
    //  fbar_cap_x;
    //  fbar_cap_y;
    //  fbar_cap_ang;
    //  fbar_bdy_x;
    //  fbar_bdy_y;
    
    /*
    ostr << position.x << ","
	 << position.y << ","
	 << position.z << endl;
    */
    return 0;
}

int nucleator::load_data(ifstream &istr) 
{
    int geom;
    istr >> geom;
    if(geom == 0)
	geometry= sphere;
    else
	geometry = capsule;
    istr >> radius;
    istr >> CAPSULE_HALF_LINEAR;
    istr >> surf_area;
    istr >> movability;
    istr >> position;
    istr >> deltanucposn;
    istr >> torque;
    istr >> centerofmass;
    istr >> momentofinertia;
    istr >> nucleator_rotation;
    
    return 0;
}

//bool nucleator::is_sphere()
//{
//    return geometry == sphere;
//}
//
//bool nucleator::is_capsule()
//{
//    return geometry == capsule;
//}


void nucleator::definecagepoints(void)
{

	const double pointdensity = 20;
	double pointspacing;

	cagepoints.resize(0);

	double r,xx,yy,zz;

	if (geometry==sphere)
	{

	// sphere

		for (double theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (double phi=-PI; phi<PI; phi+=2*PI/pointdensity)
			{
				
				r = RAD_INCOMP * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RAD_INCOMP;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
			}


	}
	else
	{


		for (double theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (double phi=-PI; phi<PI; phi+=2*PI/pointdensity)
			{
				
				r = RAD_INCOMP * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RAD_INCOMP;						// z just scaled by radius

				if (zz>0)
					zz+= (CAPSULE_HALF_LINEAR); 
				else
					zz-= (CAPSULE_HALF_LINEAR);

				cagepoints.push_back(vect(xx,yy,zz));
			}

		// cylinder
			
		pointspacing = (RAD_INCOMP * 2 *PI) / pointdensity;

		for (double theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (double z1=(-(CAPSULE_HALF_LINEAR)); z1<(0.001+CAPSULE_HALF_LINEAR); z1+= pointspacing)
			{
						
				xx = RAD_INCOMP * cos(theta);				// x and y of point
				yy = RAD_INCOMP * sin(theta);
				zz = -z1;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
				
			}

	}
}
