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
MYDOUBLE NUCPOINT_SCALE = 0.9;
#else
MYDOUBLE NUCPOINT_SCALE = 1.001;
#endif
			
nucleator::nucleator(void)
{
	radius = RADIUS;
	segment = SEGMENT;
	//P_NUC = (MYDOUBLE) 0.8 / (4*PI*radius*radius);
	geometry = sphere;
	position.zero();
	surf_area = 4 * PI * radius * radius;
    // inverse ratio of 'viscosities' of node and bead
    movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / radius;    

	deltanucposn.zero();

	torque.zero();

	centerofmass.zero();

	momentofinertia.x = 1000 * MofI;  // mark:  todo: calculate these numbers
	momentofinertia.y = 1000 * MofI;
	momentofinertia.z = 500 * MofI;

    definecagepoints();
	colour.r = colour.g = colour.b = 1.0;

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
		movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / radius;
	}
	else
	{
		surf_area = 4 * PI * radius * radius  +  2 * PI * radius * segment;
		movability = FORCE_SCALE_FACT * NODE_INCOMPRESSIBLE_RADIUS / (radius + segment/4);   // what should this be??
	}

	ptheactin = actinptr;
	ptheactin->nucleation_object = this;
	position.zero();
	definenucleatorgrid();

	deltanucposn.zero();

	centerofmass.zero();;

	torque.zero();;

	momentofinertia.x = 1000 * MofI;  // mark:  todo: calculate these numbers
	momentofinertia.y = 1000 * MofI;
	momentofinertia.z = 500 * MofI;

	if( geometry==sphere ){
	radial_rep_distrib_x.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_y.reserve(RADIAL_SEGMENTS);
	radial_rep_distrib_z.reserve(RADIAL_SEGMENTS);
	} else {
	    set_rep_bins();
	}

	definecagepoints();
	colour.r = colour.g = colour.b = 1.0;

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
					if (z < (radius) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/2) - 1))
						continue;
				if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
					if (z < (radius) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/4) - 3))
						continue;
			    if (ASYMMETRIC_NUCLEATION==4) { // fixed random location
				static MYDOUBLE fixed_x = x;
				static MYDOUBLE fixed_y = y;
				static MYDOUBLE fixed_z = z;
				x = fixed_x;
				y = fixed_y;
				z = fixed_z;				
			}
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

	MYDOUBLE rad = radius * NUCPOINT_SCALE;

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
				r = rad * sqrt(1 - z*z);		// radius of circle
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
			if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation above z=0
				if (z<0) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (segment/2 + rad) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/2) - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (segment/2 + rad) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/4) - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4)  // caps only
				if (fabs(z) < (segment/2))
					continue;
			if (ASYMMETRIC_NUCLEATION==5)  // half caps only
				if ( (fabs(z) < (segment/2)) || ((x<0)&&(z>0)) || ((x>0)&&(z<0)) )
					continue;
			if (ASYMMETRIC_NUCLEATION==6)  // caps only
				if (fabs(z) < 0.5 * (segment/2 + rad) *( (MYDOUBLE) rand() / (MYDOUBLE)(RAND_MAX/4) - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==7)  // half caps one side
				if ( (fabs(z) < (segment/2)) || ((x<0)&&(z>0)) || ((x>0)&&(z<0)) || (z>0))
					continue;
			if (ASYMMETRIC_NUCLEATION==8)  //  cap one side
				if ( (fabs(z) < (segment/2)) || (z<0))
					continue;
		}

		//}

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
			r = RAD_INCOMP * sqrt(1 - z1*z1);		// radius of circle
			
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
			r = RAD_INCOMP * sqrt(1 - z1*z1);		// radius of circle
			
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


bool nucleator::iswithinnucleator(const MYDOUBLE& x, const MYDOUBLE& y, const MYDOUBLE& z)
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

void nucleator::set_rep_bins()
{
    // vector<MYDOUBLE> radial_rep_bin_x;
    // vector<MYDOUBLE> radial_rep_bin_y;
		
    // bit of a guess, let's partition the samples as 1/2, 1/4, 1/4
    // body segments

	// try to keep equidistant:

	nbdy_segs = int( (2*segment) / ( (2 * PI * radius) / RADIAL_SEGMENTS ) );

    //nbdy_segs = int(0.5*RADIAL_SEGMENTS);
    if(nbdy_segs%2 != 0)
		nbdy_segs++;

    // body segments
    // clear everything
    fbar_bdy_x.clear();
    fbar_bdy_y.clear();
    
    MYDOUBLE seg_length = SEGMENT / (nbdy_segs/2.0 - 1);

    MYDOUBLE x,y;
    cout << "body_segments:" << nbdy_segs << endl;
    
    // partition points on the capsule body
    for(int i=0; i<(nbdy_segs/2.0); i++){
	fbar_bdy_x.push_back( -RADIUS );
	fbar_bdy_y.push_back( -SEGMENT/2.0 + i*seg_length );
    }
    for(int i=0; i<(nbdy_segs/2.0); i++){
	fbar_bdy_x.push_back( +RADIUS );
	fbar_bdy_y.push_back( -SEGMENT/2.0 + i*seg_length );
    }
    
    // cap segments
    // clear everything
    fbar_cap_x.clear();
    fbar_cap_y.clear();
    fbar_cap_ang.clear();
    
    ncap_segs = RADIAL_SEGMENTS;// - nbdy_segs;
    if(ncap_segs%2 != 0)
		ncap_segs++;

    cout << "cap_segments:" << ncap_segs << endl;
    
    MYDOUBLE seg_angle = PI/(ncap_segs/2.0);
    MYDOUBLE angle;
    
    for(int i=0; i<(ncap_segs); i++){
	angle = i*seg_angle + seg_angle/2.0;
	x = RADIUS * cos(angle);
	y = RADIUS * sin(angle);
	fbar_cap_x.push_back( x );
	fbar_cap_y.push_back( y );
	fbar_cap_ang.push_back( atan2(y, x) );	
    }
    
    radial_rep_distrib_x.reserve(nbdy_segs + ncap_segs);
    //radial_rep_distrib_x.clear();
    fill(radial_rep_distrib_x.begin(), radial_rep_distrib_x.end(), 0);

    radial_rep_distrib_y.reserve(nbdy_segs + ncap_segs);
    //radial_rep_distrib_y.clear();
    fill(radial_rep_distrib_y.begin(), radial_rep_distrib_y.end(), 0);

    radial_rep_distrib_z.reserve(nbdy_segs + ncap_segs);
    //radial_rep_distrib_z.clear();
    fill(radial_rep_distrib_z.begin(), radial_rep_distrib_z.end(), 0);
    
}

bool nucleator::collision(nodes &node)//(MYDOUBLE &x, MYDOUBLE &y, MYDOUBLE &z)
{  // returns true if succeeds, false if fails due to too great node ejection

if (USE_THREADS)
	pthread_mutex_lock(&beadmovelock_mutex); // lock the beadmove mutex

	// node has entered nucleator,
	// return co-ords of node pushed to surface...

    MYDOUBLE r, scale, z2;
	vect node_disp;
	vect oldpos;

	oldpos = node;

	MYDOUBLE rad = RAD_INCOMP * (MYDOUBLE) 1.001; // needed to prevent rounding errors putting back inside nuclator

    // FIXME: add no movement of nodes to outside the nucleator when SEED_INSIDE is set (ML)? 
	switch (geometry)
	{
	case (sphere):
		{
		r = node.length();;

		scale = rad / r;

		node*=scale;
		
		if (NUCLEATOR_FORCES)
		{
			radial_rep_distrib_z[(int)(((atan2(node.y,node.x)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(node.z))*(rad-r);
			radial_rep_distrib_x[(int)(((atan2(node.y,node.z)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(node.x))*(rad-r);
			radial_rep_distrib_y[(int)(((atan2(node.z,node.x)/PI)+1)*RADIAL_SEGMENTS)%RADIAL_SEGMENTS]+= (1-fabs(node.y))*(rad-r);
		}

		break;
		}

	case (capsule):
		//{
	 //   if(fabs(z)<=(segment/2.0)){
		//// on the cylinder body						
		//// eject from the cylinder, and place back on surface
		//		r = calcdist(x,y);
		//		scale = rad / r;
		//x *= scale; // x/rad;
		//y *= scale; // y/rad;
		//		
		//radial_rep_distrib_x[get_zbin(y,z)]   += (1-fabs(x/rad))*(rad-r); // scaled
		//radial_rep_distrib_y[get_zbin(x,z)]   += (1-fabs(y/rad))*(rad-r);
		//radial_rep_distrib_z[get_angbin(x,y)] += (rad-r); // not scaled
		///*
		//if((1-fabs(x))*(rad-r)<0)
		//    cout << "a" << x << ":" << (rad-r)<< ":" <<(1-fabs(x))*(rad-r) << endl;
		//if((1-fabs(y))*(rad-r)<0)
		//    cout << "b" << (1-fabs(y))*(rad-r) << endl;
		//if((rad-r)<0)
		//    cout << "c" << (rad-r) << endl;
		//*/
	 //   } else {
		//// on the endcaps
		//// treat locally as a sphere
  //              // move back to the surface of the nucleator
		//MYDOUBLE lz = 0;
		//lz = (z>0)?(z-SEGMENT/2.0):(lz+SEGMENT/2.0);		

		//r = calcdist(x,y,lz);
		//scale = rad / r;
		//x  *= scale;
		//y  *= scale;
		//z*=scale;
		//z = (z>0)?lz+SEGMENT/2.0:lz-SEGMENT/2.0;

  //              // on the end caps
		//radial_rep_distrib_x[nbdy_segs + get_angbin(y,lz)] += (1-fabs(x/rad))*(rad-r); // scaled
		//radial_rep_distrib_y[nbdy_segs + get_angbin(x,lz)] += (1-fabs(y/rad))*(rad-r);
		//radial_rep_distrib_z[get_angbin(x,y)]              += (1-fabs(lz/rad*scale))*(rad-r);
		///*
		//if((1-fabs(x))*(rad-r)<0)
		//    cout << "1a"<< (1-fabs(x))*(rad-r) << endl;
		//if((1-fabs(y))*(rad-r)<0)
		//    cout << "1b"<< (1-fabs(y))*(rad-r) << endl;
		//if((1-fabs(lz*scale))*(rad-r)<0)
		//    cout << "1c" << (1-fabs(lz*scale))*(rad-r) << endl;
		//*/
		//	}
		//	break;
		//}

		{
		if ((fabs(node.z)<(segment/2)))
			{ 
				// on the cylinder

				r = calcdist(node.x,node.y);
				scale = rad / r;
				node.x*=scale;
				node.y*=scale;
				
				if (NUCLEATOR_FORCES)
				{
					radial_rep_distrib_x[get_zbin(node.y,node.z)]   += (1-fabs(node.x/rad))*(rad-r); // scaled
					radial_rep_distrib_y[get_zbin(node.x,node.z)]   += (1-fabs(node.y/rad))*(rad-r);
					radial_rep_distrib_z[get_angbin(node.x,node.y)] += (rad-r); // not scaled
				}

			}
			else
			{	// on the ends
				
				if (node.z<0)  // make into a sphere again
					z2 = node.z + (segment/2);
				else
					z2 = node.z - (segment/2);

				// calculate theta, phi :

				r = calcdist(node.x,node.y,z2);
				scale = rad / r;
				node.x  *= scale;
				node.y  *= scale;
				z2 *= scale;

				if (node.z<0)  
					node.z = z2 - (segment/2);
				else
					node.z = z2 + (segment/2);
			

				if (NUCLEATOR_FORCES)
				{
					radial_rep_distrib_x[nbdy_segs + get_angbin(node.y,z2)] += (1-fabs(node.x/rad))*(rad-r); // scaled
					radial_rep_distrib_y[nbdy_segs + get_angbin(node.x,z2)] += (1-fabs(node.y/rad))*(rad-r);
					radial_rep_distrib_z[get_angbin(node.x,node.y)]         += (1-fabs(z2/rad*scale))*(rad-r);
				}

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

//torque += delta_torque;

//if (ptheactin->iteration_num % 2000 == 0)
//{
//	cout << setprecision(20) << torque.x << " " << torque.y << " " << torque.z << endl;
//}

if (USE_THREADS)
	pthread_mutex_unlock(&beadmovelock_mutex); // unlock the beadmove mutex

	return true; // sucessful node ejection
}

int nucleator::get_zbin(const MYDOUBLE x, const MYDOUBLE y)
{
    // FIXME:
    // Confusingly this is x,y where y is moving up capsule
    // and x is moving away (ie y = z, x = x|y) this made sense at the time.
    int indx = 0;
    int np = (int)fbar_bdy_y.size();
    MYDOUBLE mindist = fabs(y - fbar_bdy_y[indx]);;
    MYDOUBLE dist;
    
    for(int i=0; i<np; i++)
	{
		if(x * fbar_bdy_x[i] >= 0)
		{ // same side
			
			dist = fabs(y - fbar_bdy_y[i]);

			if(dist < mindist)
			{
				mindist = dist;
				indx = i;
			}
		}
    }
    /*
    cout << "Zbin  " 
	 << "input x, y" << x << " " << y
	 << "  matches: "
	 << fbar_cap_x[indx] << " "
	 << fbar_cap_y[indx] 
    	 << " (" << indx << ")" << endl;
    */
    return indx;
}
 
int nucleator::get_angbin(const MYDOUBLE x, const MYDOUBLE y)
{
    // angular bin, assumed to be on a circle
    int indx = -1;
    int np = (int)fbar_cap_ang.size();
    MYDOUBLE ang;
    MYDOUBLE mindiff = 2*PI;
    MYDOUBLE diff;

    ang = atan2(y, x);
    for(int i=0; i<np; i++){
	diff =  fabs(ang - fbar_cap_ang[i]);
	if(diff<mindiff){
	    mindiff = diff;
	    indx=i;
	}
    }
    /*
    cout << "AngBin  "
	 << "input x, y, angle : (" << x << ", " << y << ") "<< ang
	 << "  matches: ("
	 << fbar_cap_x[indx] << ", " 
	 << fbar_cap_y[indx] << ") " 
	 << fbar_cap_ang[indx]
	 << " indx: " << indx
	 << endl;
    */
    return indx;
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
	if (geometry == sphere)
	{	
		for (int i=0; i<RADIAL_SEGMENTS; ++i)
		{
			radial_rep_distrib_x[i]=0;
			radial_rep_distrib_y[i]=0;
			radial_rep_distrib_z[i]=0;
		}
	}
	else
	{
		for (int i=0; i<(nbdy_segs + ncap_segs); ++i)
		{
			radial_rep_distrib_x[i]=0;
			radial_rep_distrib_y[i]=0;
			radial_rep_distrib_z[i]=0;
		}
	}

	return 0;
}

int nucleator::save_data(ofstream &ostr) 
{
    ostr << geometry << endl;
    ostr << radius << endl;
    ostr << segment << endl;
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
    istr >> segment;
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

bool nucleator::is_sphere()
{
    return geometry == sphere;
}

bool nucleator::is_capsule()
{
    return geometry == capsule;
}

int nucleator::n_force_segments()
{
    return nbdy_segs + ncap_segs;
}

void nucleator::definecagepoints(void)
{

	const MYDOUBLE pointdensity = 20;
	MYDOUBLE pointspacing;

	cagepoints.resize(0);

	MYDOUBLE r,xx,yy,zz;

	if (geometry==sphere)
	{

	// sphere

		for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (MYDOUBLE phi=-PI; phi<PI; phi+=2*PI/pointdensity)
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


		for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (MYDOUBLE phi=-PI; phi<PI; phi+=2*PI/pointdensity)
			{
				
				r = RAD_INCOMP * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RAD_INCOMP;						// z just scaled by radius

				if (zz>0)
					zz+= (segment/2); 
				else
					zz-= (segment/2);

				cagepoints.push_back(vect(xx,yy,zz));
			}

		// cylinder
			
		pointspacing = (RAD_INCOMP * 2 *PI) / pointdensity;

		for (MYDOUBLE theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (MYDOUBLE z1=(-(segment/2)); z1<(0.001+segment/2); z1+= pointspacing)
			{
						
				xx = RAD_INCOMP * cos(theta);				// x and y of point
				yy = RAD_INCOMP * sin(theta);
				zz = -z1;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
				
			}

	}
}
