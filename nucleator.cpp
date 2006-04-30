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

// NUCPOINT_SCALE determines where to put the nodes after they are ejected
// slightly above 1 to stop rounding errors

#ifdef SEED_INSIDE			
double NUCPOINT_SCALE = 0.9;				 
#else
double NUCPOINT_SCALE = 1.001;
#endif
			
nucleator::nucleator(void)
{  
}

nucleator::~nucleator(void)
{
}

nucleator::nucleator(shape set_geometry)//, actin *actinptr)
{
	//radius = RADIUS;
	//segment = 2*CAPSULE_HALF_LINEAR;
	//P_NUC = (double)0.8 / (4*PI*radius*radius);

	geometry = set_geometry;

	if (geometry==sphere)
	{
		surf_area = 4 * PI * RADIUS * RADIUS;
		movability = 0.25 / (RADIUS * NUCLEATOR_INERTIA);

		momentofinertia.x = 1000 * MofI;  // mark:  todo: calculate these numbers
		momentofinertia.y = 1000 * MofI;
		momentofinertia.z = 1000 * MofI;
	}
	else
	{
		surf_area = 4 * PI * RADIUS * RADIUS  +  4 * PI * RADIUS * CAPSULE_HALF_LINEAR;
		movability = 0.25 / (NUCLEATOR_INERTIA * (RADIUS + CAPSULE_HALF_LINEAR/2));   // what should this be??
	
	    momentofinertia.x = 1000 * CAPSULE_HALF_LINEAR * MofI;  // mark:  todo: calculate these numbers
		momentofinertia.y = 1000 * CAPSULE_HALF_LINEAR * MofI;
		momentofinertia.z = 1000 * MofI;
	}



	//ptheactin = actinptr;
	ptheactin->p_nuc = this;

	position.zero();
	//definenucleatorgrid();

	deltanucposn.zero();

	centerofmass.zero();

	torque.zero();

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

	int nodestoadd = (int) floatingnodestoadd;  // actual nodes have to be integer

	if (( floatingnodestoadd - nodestoadd ) * RAND_MAX > rand() )
		nodestoadd++;   // deal with fractional nodes using probability

	for (int i=0; i != nodestoadd; ++i)
	{
		z = (2 * (rand() * RECIP_RAND_MAX) ) -1 ;		// random number -1 to 1
		theta = (2 * PI * (rand() * RECIP_RAND_MAX));  // circle vector
		
		if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
		{
			r = RADIUS * sqrt(1 - z*z);		// radius of circle
		}
		else
		{
			r = RADIUS;  
		}

        //cout << "RECIP_RAND_MAX*RAND_MAX: " <<RAND_MAX*RECIP_RAND_MAX << " RAND_MAX:" <<RAND_MAX << endl;

		x =  r * cos(theta) * NUCPOINT_SCALE;   // x and y of point
		y =  r * sin(theta) * NUCPOINT_SCALE;

		z *=  RADIUS * NUCPOINT_SCALE;          // z just scaled by radius

		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION==1)  // no nucleation above z=0
				if ((y<0) || (fabs(x+z)>0.5)) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (RADIUS) *( 2 * (rand() * RECIP_RAND_MAX) - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (RADIUS) *( 4 * (rand() * RECIP_RAND_MAX) - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4) 
            { // fixed random location
			static double fixed_x = x;
			static double fixed_y = y;
			static double fixed_z = z;
			x = fixed_x;
			y = fixed_y;
			z = fixed_z;				
		    }
			if (ASYMMETRIC_NUCLEATION==7)  // half caps one side
				if ( ( (y>0) && (z>0) ) || 
                     ( (y<0) && (z<0) ) || 
                       (z<0) )
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

	double rad = RADIUS * NUCPOINT_SCALE;

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

		z = ( 2 * (rand() * RECIP_RAND_MAX) ) - 1 ;  // random number -1 to 1
		theta = (2 * PI * (rand() * RECIP_RAND_MAX));  // circle vector
		
		onseg = ( (CAPSULE_HALF_LINEAR /(RADIUS+CAPSULE_HALF_LINEAR)) > ( rand() * RECIP_RAND_MAX )); // on ends or on segment?
		
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

			z *= rad;
						
			if (z>0)
				z += CAPSULE_HALF_LINEAR; 
			else
				z -= CAPSULE_HALF_LINEAR;
		}

		x =  r * cos(theta); // x and y of point
		y =  r * sin(theta);
		
		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION==1)  /// no nucleation above z=0
				if (z<0) continue;
			if (ASYMMETRIC_NUCLEATION==2)  // linear degredation to zero
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 2*  (rand() * RECIP_RAND_MAX) - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 4 * (rand() * RECIP_RAND_MAX) - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4)  // caps only
				if (fabs(z) < (CAPSULE_HALF_LINEAR))
					continue;
			if (ASYMMETRIC_NUCLEATION==5)  // half caps only
				if ( (fabs(z) < (CAPSULE_HALF_LINEAR)) || ((x<0)&&(z>0)) || ((x>0)&&(z<0)) )
					continue;
			if (ASYMMETRIC_NUCLEATION==6)  // caps only
				if (fabs(z) < 0.5 * (CAPSULE_HALF_LINEAR + rad) * ( 4 * (rand() * RECIP_RAND_MAX) - 3))
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

bool nucleator::collision(nodes &node)//(double &x, double &y, double &z)
{  // returns true if succeeds, false if fails due to too great node ejection

	// node has entered nucleator,
	// return co-ords of node pushed to surface...

    double r, scale, z2;
	vect node_disp;
	vect oldpos;

	oldpos = node;

	double rad = RADIUS * NUCPOINT_SCALE; // needed to prevent rounding errors putting back inside nuclator

    // FIXME: add no movement of nodes to outside the nucleator when SEED_INSIDE is set (ML)? 
	switch (geometry)
	{
	    case (sphere):
	    {
		    r = node.length();

			if ((!node.stucktonucleator) && (STICK_TO_NUCLEATOR)) // re-stick to nucleator if come off
			{
				node.stucktonucleator = true;
				node.nucleator_stuck_position = node * (1/r); // link to point *on* the nucleator surface
			}

			scale = rad / r;

		    node *= scale;
        	
            node.nucleator_impacts += (rad-r);


		    break;
	    }

	    case (capsule):

	    {
		    if ((fabs(node.z)<CAPSULE_HALF_LINEAR))
		    { 
			    // on the cylinder

			    r = calcdist(node.x,node.y);

				if ((!node.stucktonucleator) && (STICK_TO_NUCLEATOR) && (RESTICK_TO_NUCLEATOR)) // re-stick to nucleator if come off
				{
					node.stucktonucleator = true;
					node.nucleator_stuck_position.x = node.x * (1/r);
					node.nucleator_stuck_position.y = node.y * (1/r);
					node.nucleator_stuck_position.z = node.z;// link to point *on* the nucleator surface
				}

			    scale = rad * (1/r);
			    node.x *= scale;
			    node.y *= scale;
    			
                node.nucleator_impacts += (rad-r);


		    }
		    else
		    {	// on the ends
    			
			    // calculate theta, phi :

				if (node.z<0)  // make into a sphere again
				    z2 = node.z + CAPSULE_HALF_LINEAR;
			    else
				    z2 = node.z - CAPSULE_HALF_LINEAR;

				r = calcdist(node.x,node.y,z2);

				if ((!node.stucktonucleator) && (STICK_TO_NUCLEATOR) && (RESTICK_TO_NUCLEATOR)) // re-stick to nucleator if come off
				{
					node.stucktonucleator = true;
					node.nucleator_stuck_position.x = node.x * (1/r);
					node.nucleator_stuck_position.y = node.y * (1/r);
					
					if (node.z<0)  // make into a sphere again
						node.nucleator_stuck_position.z = z2 * (1/r) - (CAPSULE_HALF_LINEAR);// link to point *on* the nucleator surface
					else
						node.nucleator_stuck_position.z = z2 * (1/r) + (CAPSULE_HALF_LINEAR);// link to point *on* the nucleator surface
					
				}

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

    if ((fabs(node_disp.x) > 0.2*RADIUS) ||
	    (fabs(node_disp.y) > 0.2*RADIUS) ||
	    (fabs(node_disp.z) > 0.2*RADIUS))
    {
		cout << endl;
		cout << "node " << node.nodenum << " nucleus ejection too great. radius:" << oldpos.length() <<  endl;
	    cout << "old (x,y,z): " <<  oldpos.x << ", " << oldpos.y << ", " << oldpos.z << endl;
	    cout << "new (x,y,z): " <<  node.x << ", " << node.y << ", " << node.z <<  endl;
		cout << "#links: " << (int) node.listoflinks.size() << endl;
	    return false;  // failed
    }


#ifndef SEED_INSIDE

    node -= node_disp * movability;  // move the node *back* scaled by movibility (since nucleator should move this much

    move_nuc(oldpos,node_disp);

#endif

	return true; // sucessful node ejection
}

void nucleator::move_nuc(vect& origin_of_movement, vect& tomove)
{
	/// updates the vectors for moving and rotating the nucleator
	/// taking tomove as the movement vector from the origin_of_movement
	/// doesn't actually move anything here---that is done in actin::applyforces()

    // displacement is easy:
    deltanucposn -= tomove * movability;

    // rotate the nucleator

    // calculate torque (as pos x disp)

    if (ROTATION)
    {
        vect lever_arm = origin_of_movement - centerofmass;
        torque += lever_arm.cross(tomove);
    }

}

int nucleator::save_data(ofstream &ostr) 
{
    ostr << geometry << endl;
    ostr << RADIUS << endl;
    ostr << CAPSULE_HALF_LINEAR << endl;
    ostr << surf_area << endl;
    ostr << movability << endl;
    ostr << position << endl;
    ostr << deltanucposn << endl;
    ostr << torque << endl;
    ostr << centerofmass << endl;
    ostr << momentofinertia << endl;
//    ostr << nucleator_rotation << endl;
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
    istr >> RADIUS;
    istr >> CAPSULE_HALF_LINEAR;
    istr >> surf_area;
    istr >> movability;
    istr >> position;
    istr >> deltanucposn;
    istr >> torque;
    istr >> centerofmass;
    istr >> momentofinertia;
//    istr >> nucleator_rotation;
    
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
				
				r = RADIUS * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RADIUS;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
			}


	}
	else
	{


		for (double theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (double phi=-PI; phi<PI; phi+=2*PI/pointdensity)
			{
				
				r = RADIUS * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RADIUS;						// z just scaled by radius

				if (zz>0)
					zz+= (CAPSULE_HALF_LINEAR); 
				else
					zz-= (CAPSULE_HALF_LINEAR);

				cagepoints.push_back(vect(xx,yy,zz));
			}

		// cylinder
			
		pointspacing = (RADIUS * 2 *PI) / pointdensity;

		for (double theta=-PI; theta<PI; theta+=2*PI/pointdensity)
			for (double z1=(-(CAPSULE_HALF_LINEAR)); z1<(0.001+CAPSULE_HALF_LINEAR); z1+= pointspacing)
			{
						
				xx = RADIUS * cos(theta);				// x and y of point
				yy = RADIUS * sin(theta);
				zz = -z1;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
				
			}

	}
}
