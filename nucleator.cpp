/*

    comet - an actin-based motility simulator
    Copyright (C) 2009 Mark J Dayel

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

    When using the code in academic work please cite Dayel et al. 2009 
    and include any modifications to the source code with the published work.

*/

#include "stdafx.h"
#include "nucleator.h"

// NUCPOINT_SCALE determines where to put the nodes after they are ejected
// slightly above 1 to stop rounding errors

#ifdef SEED_INSIDE			
double NUCPOINT_SCALE = 0.9;				 
#else
double NUCPOINT_SCALE = 1.0001;
#endif
			
//nucleator::nucleator(void)
//{  
//}

nucleator::~nucleator(void)
{
}


double arctanh(double x) 
{ 
  if (fabs(x)>0.9999999) 
    return (x/fabs(x))*1e18;

  if (x < 0.000000000001) 
      return 0.0;

  return log((1.0+x)/(1.0-x)) / 2.0 ;
}

nucleator::nucleator(void)//, actin *actinptr)
{

	if (NUCSHAPE==sphere)
	{
		surf_area = 4 * PI * RADIUS * RADIUS;
        
        if (VARY_INERT_W_RAD)
        {   // movability is a *fraction* representing the ratio of how movable the nucleator is compared to the nodes
            inertia = 4 * (RADIUS * NUCLEATOR_INERTIA);
        }
        else
        {
            inertia = 4 * (2.5 * NUCLEATOR_INERTIA);
        }

		momentofinertia.x = 1000 * RADIUS * MOFI;  // mark:  todo: calculate these numbers
		momentofinertia.y = 1000 * RADIUS * MOFI;
		momentofinertia.z = 1000 * RADIUS * MOFI;
	}
	else if (NUCSHAPE==ellipsoid)
	{

        double a = RADIUS;
        //double b = RADIUS;
        double c = RADIUS * ELLIPSOID_STRETCHFACTOR;

        if (ELLIPSOID_STRETCHFACTOR > 1)
        {   // prolate
            double OE = acos( a / c );
            surf_area = 2 * PI * ( a*a + c*c * OE / tan(OE));
        }
        else
        {   // oblate
            double OE = acos( c / a );
            surf_area = 2 * PI * ( a*a + c*c * arctanh(sin(OE)) / sin(OE));
        }


		//surf_area = 4 * PI * RADIUS * RADIUS * ELLIPSOID_STRETCHFACTOR;
        
        if (VARY_INERT_W_RAD)
        {   // movability is a *fraction* representing the ratio of how movable the nucleator is compared to the nodes
            inertia = 4 * (RADIUS * NUCLEATOR_INERTIA * ELLIPSOID_STRETCHFACTOR) ;
        }
        else
        {
            inertia = 4 * (2.5 * NUCLEATOR_INERTIA);
        }

		momentofinertia.x = 1000 * RADIUS * MOFI * ELLIPSOID_STRETCHFACTOR;  // mark:  todo: calculate these numbers
		momentofinertia.y = 1000 * RADIUS * MOFI * ELLIPSOID_STRETCHFACTOR;
		momentofinertia.z = 1000 * RADIUS * MOFI;
	}
	else
	{
		surf_area = 4 * PI * RADIUS * RADIUS  +  4 * PI * RADIUS * CAPSULE_HALF_LINEAR;
		
        if (VARY_INERT_W_RAD)
        {
            inertia = 4 * (NUCLEATOR_INERTIA * (RADIUS + CAPSULE_HALF_LINEAR/2));   // what should this be??
        }
        else
        {
            inertia = 4 * (NUCLEATOR_INERTIA * (2.5 + 1.5/2));   // what should this be??
        }
        

	    momentofinertia.x = 1000 * RADIUS * MOFI * CAPSULE_HALF_LINEAR;  // mark:  todo: calculate these numbers
		momentofinertia.y = 1000 * RADIUS * MOFI * CAPSULE_HALF_LINEAR;
		momentofinertia.z = 1000 * RADIUS * MOFI;
	}



	//ptheactin = actinptr;
	ptheactin->p_nuc = this;

	position.zero();
    last_delta_position.zero();

    deltanucposn_sum.zero();

    //position = vect(5,5,5);

	//definenucleatorgrid();

	deltanucposn.zero();

	centerofmass.zero();

	torque.zero();

    last_torque_rotate.settoidentity();

	definecagepoints();
	colour.r = colour.g = colour.b = 1.0;

	segs.setupsegments(this, ptheactin);

}

int nucleator::addnodes(void) const
{
	int nodesadded = 0;
	switch (NUCSHAPE)
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
    case (ellipsoid):
		{
		nodesadded = addnodesellipsoid();
		break;
		}
	}
	return nodesadded;
}

int nucleator::addnodessphere(void) const
{
	int nodesadded = 0;
	double x,y,z,r,theta;
	
	double floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;  // actual nodes have to be integer

	if ( prob_to_bool( floatingnodestoadd - nodestoadd ) )
		nodestoadd++;   // deal with fractional nodes using probability

	for (int i=0; i != nodestoadd; ++i)
	{
		z = 2 * rand_0to1() - 1 ;		// random number -1 to 1
		theta = 2 * PI * rand_0to1();  // circle vector
		
		if (z*z<1) // avoid floating exception due to rounding errors causing -ve sqrt
		{
			r = RADIUS * sqrt(1 - z*z);		// radius of circle
		}
		else
		{
			r = RADIUS;  
		}

		x =  r * cos(theta) * NUCPOINT_SCALE;   // x and y of point
		y =  r * sin(theta) * NUCPOINT_SCALE;

		z *=  RADIUS * NUCPOINT_SCALE;          // z just scaled by radius


        // modifier for asymmetric nucleation

		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION == 1)  // no nucleation above z=0
				if ((y < 0) || (fabs(x+z) > 0.5)) continue;
			if (ASYMMETRIC_NUCLEATION == 2)  // linear degredation to zero
				if (z < RADIUS * ( 2 * rand_0to1() - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION == 3)  // linear degredation
				if (z < RADIUS * ( 4 * rand_0to1() - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION == 4) 
            { // fixed random location
			    static double fixed_x = x;
			    static double fixed_y = y;
			    static double fixed_z = z;
			    x = fixed_x;
			    y = fixed_y;
			    z = fixed_z;				
		    }
			if (ASYMMETRIC_NUCLEATION == 7)  // half caps one side
				if ( ( (y > 0) && (z > 0) ) || 
                     ( (y < 0) && (z < 0) ) || 
                       (z < 0) )
					continue;
		}
    
        // modifier for +ve feedback in polymerisation

        if (POLY_FEEDBACK)
        {   // calculate probability of polymerization due to +ve feedback (based on nearby nodes)    
            if ( prob_to_bool( polyfeedbackprob(x,y,z) ) )
            {   
                // don't add node
                ptheactin->attemptedpolrate++;  // an aborted node should appear as an attempt
                continue;   // don't polymerise
            }
        }

        // transform co-ords to lab frame
        vect nucframepos(x,y,z);
        vect worldframepos=nucframepos;

        ptheactin->nuc_to_world_frame(worldframepos);

        double COVERSLIP_DIST_NO_POL = RADIUS / 5.0; // shut off polymerization within 1/10 diameter of coverslip

        if ((COVERSLIPGAP > 0) &&
            ((worldframepos.x * 2.0 >   (COVERSLIPGAP - COVERSLIP_DIST_NO_POL) ) ||
             (worldframepos.x * 2.0 < - (COVERSLIPGAP - COVERSLIP_DIST_NO_POL) )) )  
        {   // skip of outside coverslip
            continue;
        }




        // add the new node:
        
        ptheactin->node[ptheactin->highestnodecount].polymerize(worldframepos);             

        // set the stuck position to the nuc frame pos'n
        ptheactin->node[ptheactin->highestnodecount].nucleator_stuck_position = nucframepos;

        ptheactin->highestnodecount++;

        nodesadded++;
	}

	return nodesadded;
}

int nucleator::addnodesellipsoid(void) const
{
	int nodesadded = 0;
	double x,y,z,r,theta;
    double sz; // sphere co-ords

	double floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;  // actual nodes have to be integer

	if ( prob_to_bool( floatingnodestoadd - nodestoadd ) )
		nodestoadd++;   // deal with fractional nodes using probability

	while (nodesadded < nodestoadd)
	{
		sz = 2 * rand_0to1() - 1 ;		// random number -1 to 1
		theta = 2 * PI * rand_0to1();  // circle vector
		
		if (sz*sz<1) // avoid floating exception due to rounding errors causing -ve sqrt
		{
			r = RADIUS * sqrt(1 - sz*sz);		// radius of circle
		}
		else
		{
			r = RADIUS;  
		}

		x =  r * cos(theta) * NUCPOINT_SCALE;   // x and y of point
		y =  r * sin(theta) * NUCPOINT_SCALE;

		sz *=  RADIUS * NUCPOINT_SCALE;          // z just scaled by radius

        z = sz * ELLIPSOID_STRETCHFACTOR; // ellipse z is sphere z stretched by ELLIPSOID_STRETCHFACTOR
        
        // reject point with probability proportional to distance moved 
        if (ELLIPSOID_STRETCHFACTOR >= 1.0)
        {   // prolate
            if (prob_to_bool(fabs( (sz - z) / (RADIUS * ELLIPSOID_STRETCHFACTOR) )))
                continue;  
        } else
        {   // oblate
            if (prob_to_bool(1 - fabs( (sz - z) / (RADIUS) )))
                continue;
            //if (prob_to_bool( pow(((r - fabs(z)) / RADIUS) ,0.5) ))
            //    continue;

            //if (prob_to_bool( pow( 1- fabs( (sz - z) / (RADIUS * ELLIPSOID_STRETCHFACTOR) ) , 2) ))
            //    continue;
        }

        // modifier for asymmetric nucleation

		if (ASYMMETRIC_NUCLEATION!=0)
		{
			if (ASYMMETRIC_NUCLEATION == 1)  // no nucleation above z=0
				if ((y < 0) || (fabs(x+z) > 0.5)) continue;
			if (ASYMMETRIC_NUCLEATION == 2)  // linear degredation to zero
				if (z < RADIUS * ( 2 * rand_0to1() - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION == 3)  // linear degredation
				if (z < RADIUS * ( 4 * rand_0to1() - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION == 4) 
            { // fixed random location
			    static double fixed_x = x;
			    static double fixed_y = y;
			    static double fixed_z = z;
			    x = fixed_x;
			    y = fixed_y;
			    z = fixed_z;				
		    }
			if (ASYMMETRIC_NUCLEATION == 7)  // half caps one side
				if ( ( (y > 0) && (z > 0) ) || 
                     ( (y < 0) && (z < 0) ) || 
                       (z < 0) )
					continue;
		}

        // modifier for +ve feedback in polymerisation

        if (POLY_FEEDBACK)
        {   // calculate probability of polymerization due to +ve feedback (based on nearby nodes)    
            if ( prob_to_bool( polyfeedbackprob(x,y,z) ) )
            {   
                // don't add node
                ptheactin->attemptedpolrate++;  // an aborted node should appear as an attempt
                continue;   // don't polymerise
            }
        }

        // transform co-ords to lab frame
        vect nucframepos(x,y,z);
        vect worldframepos=nucframepos;

        ptheactin->nuc_to_world_frame(worldframepos);

        if ((COVERSLIPGAP > 0) &&
            ((worldframepos.x * 2 >  COVERSLIPGAP) ||
             (worldframepos.x * 2 < -COVERSLIPGAP)) )  
        {   // skip of outside coverslip
            continue;
        }


        // add the new node:
        
        ptheactin->node[ptheactin->highestnodecount].polymerize(worldframepos);             

        // set the stuck position to the nuc frame pos'n
        ptheactin->node[ptheactin->highestnodecount].nucleator_stuck_position = nucframepos;
        ptheactin->highestnodecount++;
        nodesadded++;
	} 

	return nodesadded;
}

int nucleator::addnodescapsule(void) const
{
	int nodesadded = 0;
	double x,y,z,r,theta;
	bool onseg;

	double rad = RADIUS * NUCPOINT_SCALE;

	double floatingnodestoadd = DELTA_T * P_NUC * surf_area;  // number of nodes to add

	int nodestoadd = (int) floatingnodestoadd;

	if ( prob_to_bool( floatingnodestoadd - nodestoadd ) )
		nodestoadd++;

	for (int i=0; i != nodestoadd; i++)
	{
		//	Pick a random point on capsule:

		// first choose whether on sphere or capsule:
		// p(sphere) = Area of sphere/area of capsule
		//           = 4.PI.r^2 / (4.PI.r^2 + 2.PI.r.h)
		//           = r/(r+2h)


		z =  2 * rand_0to1() - 1 ;  // random number -1 to 1
		theta = 2 * PI * rand_0to1() ;  // circle vector
		
		onseg = ( prob_to_bool(CAPSULE_HALF_LINEAR /(RADIUS+CAPSULE_HALF_LINEAR)) ); // on ends or on segment?
		
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
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 2 * rand_0to1() - 1))
					continue;
			if (ASYMMETRIC_NUCLEATION==3)  // linear degredation
				if (z < (CAPSULE_HALF_LINEAR + rad) *( 4 * rand_0to1() - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==4)  // caps only
				if (fabs(z) < (CAPSULE_HALF_LINEAR))
					continue;
			if (ASYMMETRIC_NUCLEATION==5)  // half caps only
				if ( (fabs(z) < (CAPSULE_HALF_LINEAR)) || ((x<0)&&(z>0)) || ((x>0)&&(z<0)) )
					continue;
			if (ASYMMETRIC_NUCLEATION==6)  // caps only
				if (fabs(z) < 0.5 * (CAPSULE_HALF_LINEAR + rad) * ( 4 * rand_0to1() - 3))
					continue;
			if (ASYMMETRIC_NUCLEATION==7)  // half caps one side
				if ( (fabs(z) < (0.7*CAPSULE_HALF_LINEAR)) || ((x>0)&&(z>0)) || ((x<0)&&(z<0)) || (z>0))
					continue;
			if (ASYMMETRIC_NUCLEATION==8)  //  cap one side
				if ( (fabs(z) < (CAPSULE_HALF_LINEAR)) || (z<0))
					continue;
            if (ASYMMETRIC_NUCLEATION==9)  /// no nucleation above z=0 , fin on one side :)
            {
                if (z<0) continue;

                //double FINPITCH =  0.1; // 0.1  ; // turns per length
                //double FINWIDTHANGLE = 30 * PI / 180; // the angular width of the 
                //double FINRATIO = 0.1; // 0 to 1, ratio of fin to non-fin

                double finangle = z * (FINPITCH * ( 2 * PI / 2 * CAPSULE_HALF_LINEAR ) );
                bool finnoskip = prob_to_bool(FINRATIO);

                if ( fabs( fmod((theta - finangle * PI / 180 ), 2*PI )) * 2 < FINWIDTHANGLE  )
                {   // on the fin
                    if (!finnoskip)
                        continue;
                }
                else
                {
                    if (finnoskip)
                       continue;
                }

            }

            if (ASYMMETRIC_NUCLEATION==10)  /// no nucleation above z=0 , fin on one side :)
            {
                if (z<-CAPSULE_HALF_LINEAR/2) continue;

                //double FINPITCH =  0.1; // 0.1  ; // turns per length
                //double FINWIDTHANGLE = 30 * PI / 180; // the angular width of the 
                //double FINRATIO = 0.1; // 0 to 1, ratio of fin to non-fin

                double finangle = z * (FINPITCH * ( 2 * PI / 2 * CAPSULE_HALF_LINEAR ) );
                bool finnoskip = prob_to_bool(FINRATIO);

                if ( fabs( fmod((theta - finangle * PI / 180 ), 2*PI )) * 2 < FINWIDTHANGLE  )
                {   // on the fin
                    if (!finnoskip)
                        continue;
                }
                else
                {
                    if (finnoskip)
                       continue;
                }

            }


		}


        // modifier for +ve feedback in polymerisation

        if (POLY_FEEDBACK)
        {   // calculate probability of polymerization due to +ve feedback (based on nearby nodes)    
            if ( prob_to_bool( polyfeedbackprob(x,y,z) ) )
            {   
                // don't add node
                ptheactin->attemptedpolrate++;  // an aborted node should appear as an attempt
                continue;   // don't polymerise
            }
        }

        // transform co-ords to lab frame
        vect nucframepos(x,y,z);
        vect worldframepos=nucframepos;

        ptheactin->nuc_to_world_frame(worldframepos);

        // add the new node:
        
        ptheactin->node[ptheactin->highestnodecount].polymerize(worldframepos);

        // set the stuck position to the nuc frame pos'n
        ptheactin->node[ptheactin->highestnodecount].nucleator_stuck_position = nucframepos;
        ptheactin->highestnodecount++;
        nodesadded++;

	}

	return nodesadded;
}


double nucleator::polyfeedbackprob(const double& x, const double& y, const double& z) const
{  ///  calculate probability of polymerization due to +ve feedback (based on nearby nodes and distances)
    
    // find nodes within radius
    // create a dummy temp node for the findnearbynodes call
    // which needs to know the gridcoords for the node

    nodes tempnode;

    tempnode.x = x;
    tempnode.y = y;
    tempnode.z = z;

    tempnode.setgridcoords();

    int POL_XLINK_GRIDSEARCH = (int) floor( POLY_FEEDBACK_DIST / GRIDRES ) + 1;

    ptheactin->findnearbynodes(tempnode,POL_XLINK_GRIDSEARCH,0);

    // found nodes, now sum weights

    double prob_weight_sum = 0.0;

    for (vector <nodes*>::iterator nearnode  = ptheactin->recti_near_nodes[0].begin();
                                   nearnode != ptheactin->recti_near_nodes[0].end(); 
                                 ++nearnode)
    {
        if (  (!(*nearnode)->harbinger) &&
              ( (*nearnode)->polymer) )  // only count real nodes
        {  
            // calculate distance
            double dist = ((*(*nearnode)) - tempnode).length();

            if (dist < POLY_FEEDBACK_DIST)
            {   // if within the range
                prob_weight_sum += 1 / dist;
            }
        }
    }

    if (prob_weight_sum < SQRT_ACCURACY_LOSS)
        return POLY_FEEDBACK_MIN_PROB;  // no real nodes were within range

    return mymax(POLY_FEEDBACK_MIN_PROB, 1 - (POLY_FEEDBACK_FACTOR / prob_weight_sum)); 

}


int nucleator::savevrml(ofstream *outputstream) 
{  /// sends the cage points for the nucleator to a vrml stream

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

	for (int i=0; i != numpoints; i++)
	{
			if (i==0)
				*outputstream	<< "1 1 1";
			else
				*outputstream	<< ",1 1 1";
	}

	return 0;
}


vect nucleator::eject_point(const vect& inputpoint) const
{   /// returns the closest point on the surface of the nucleator to the point given

    switch (NUCSHAPE)
	{
	    case (sphere):
	    {
            return eject_sphere(inputpoint);
		    break;
	    }

        case (ellipsoid):
	    {
            return eject_ellipsoid(inputpoint);
		    break;
	    }

	    case (capsule):

	    {
            return eject_capsule(inputpoint);
		    break;
	    }

    }

    cerr << "Unknown nucleator shape!" << endl;

    return inputpoint; // shouldn't get here!
}


vect nucleator::eject_sphere(const vect& inputpoint) const
{
    return inputpoint.unitvec() * (RADIUS * NUCPOINT_SCALE);
}

vect nucleator::eject_ellipsoid(const vect& inputpoint) const
{
    double nlen1 = RADIUS * NUCPOINT_SCALE;
    //double nlen2 = rad;

    vect new_posonellipse, posonsphere, iterpos, ejection_normal;

    iterpos = inputpoint;

    while  (nlen1 > ((RADIUS * NUCPOINT_SCALE) * 0.002))   // relates to how close before we accept it
    {

        // map onto sphere
        posonsphere.x = iterpos.x;
        posonsphere.y = iterpos.y;
        posonsphere.z = iterpos.z / ELLIPSOID_STRETCHFACTOR;

        // eject to sphere
        new_posonellipse = posonsphere.unitvec() * (RADIUS * NUCPOINT_SCALE); // this is sphere not really ellipsoid yet
        new_posonellipse.z *= ELLIPSOID_STRETCHFACTOR; // stretch z to map onto ellipsoid
        
        // this is not right yet---the ejection normals are not correct
        // but we can transform them below

        // find the ejection normal, shrink it along ellipse axis, and stretch along it's length

        ejection_normal = new_posonellipse - iterpos;
        nlen1 = ejection_normal.length();

        ejection_normal.z /= (ELLIPSOID_STRETCHFACTOR * ELLIPSOID_STRETCHFACTOR);
        //nlen2 = ejection_normal.length(); // squish along ellipse axis

        //ejection_normal *= 0.9 * nlen1 / nlen2; // scale length 
        

        iterpos = iterpos + ejection_normal;

    }

    return iterpos;
}

vect nucleator::eject_capsule(const vect& inputpoint) const
{
    vect ejectedpoint=inputpoint;
    double z2, scale;

    if ((fabs(ejectedpoint.z) < CAPSULE_HALF_LINEAR))
    { 
	    // on the cylinder

        // note the nucleator_stuck_position is used to produce the patterned speckle
        // tracks and should represent the last nucleator collision position

	    scale = (RADIUS * NUCPOINT_SCALE) / calcdist(ejectedpoint.x, ejectedpoint.y);

	    ejectedpoint.x *= scale;
	    ejectedpoint.y *= scale;
    }
    else
    {	// on the ends
		
	    // calculate theta, phi :

		if (ejectedpoint.z<0)  // make into a sphere again
		    z2 = ejectedpoint.z + CAPSULE_HALF_LINEAR;
	    else
		    z2 = ejectedpoint.z - CAPSULE_HALF_LINEAR;

	    scale = (RADIUS * NUCPOINT_SCALE) / calcdist(ejectedpoint.x,ejectedpoint.y,z2);
	    ejectedpoint.x *= scale;
	    ejectedpoint.y *= scale;
	    z2 *= scale;

		if (ejectedpoint.z<0)
		    ejectedpoint.z = z2 - CAPSULE_HALF_LINEAR;
	    else
		    ejectedpoint.z = z2 + CAPSULE_HALF_LINEAR;
    }

    return ejectedpoint;
}

bool nucleator::collision(nodes &node_world)//(double &x, double &y, double &z)
{   /// returns true if succeeds, false if fails due to too great node ejection

	/// node has entered nucleator,
	/// push to surface along normal
    /// and use the push vector to move and rotate nucleator

	const vect oldpos = node_world.pos_in_nuc_frame;

    vect nodepos = node_world;
   
    // FIXME: add no movement of nodes to outside the nucleator when SEED_INSIDE is set (ML)? 
	if (NUC_FRICTION_COEFF < 0.0001)  
    {   // if no friction, then just eject point
        node_world.pos_in_nuc_frame = eject_point(node_world.pos_in_nuc_frame);
    }
    else
    {   // friction
        double normaldist; 
        double len_ratio;
        vect tangdisp, frictdisp; 
        vect inffrictionpoint, zerofrictpoint;

        nodepos -= node_world.delta; // find previous node position (world frame)

        nodepos -= position; // center it on the nucleator

        // rotate the node along with the nucleator by last rotation
        // as though it were stuck to it
        last_torque_rotate.inverse().rotate(nodepos);    

        ptheactin->world_to_nuc_rot.rotate(nodepos); // and bring into nucleator frame

        //vect rotated_last_delta_position = last_delta_position;
        //ptheactin->world_to_nuc_rot.rotate(rotated_last_delta_position);

        //nodepos = eject_point(nodepos);
        nodepos += last_delta_position;  // move by the nucleator motion

        // the combination of these two rotations brings the point into the nucleator frame
        // and also means that if the nucleator rotated or moved last time
        // the point will be swept along (as though infinte friction)
        
        // (assuming that since it collided, the point was very close to the surface last time
        // ideally this should be the intersection of the vector with the sphere
        // but this is a good approximation)

        // so, point on surface if infinite friction, i.e. project last node position onto surface
            
        inffrictionpoint = eject_point(nodepos);

        // point on surface if zero friction, i.e. project current position onto surface
	    zerofrictpoint = eject_point(node_world.pos_in_nuc_frame);

        // calculate normal displacement
        normaldist = (zerofrictpoint - node_world.pos_in_nuc_frame).length(); 
        
        // calculate the tangential displacement for infinite friciton
        tangdisp = inffrictionpoint - zerofrictpoint; 
       
        // friction is the normal displacement * frict coeff in the direction of the movement
        frictdisp = tangdisp.unitvec() * (normaldist * NUC_FRICTION_COEFF);   // in direction of node movement

        // find ratio of max movment to proposed friction movement
        len_ratio = tangdisp.length() / frictdisp.length();

        if (len_ratio < 1.0)
        {   // friction can't be greater than tangential force
            frictdisp *= len_ratio;
        }

        // move node to final point                           

        node_world.pos_in_nuc_frame = zerofrictpoint + frictdisp; // tangdisp * 0.9 ; 
    }    

    //node_world.previous_pos_in_nuc_frame = oldpos;

    if ((!node_world.stucktonucleator) && STICK_TO_NUCLEATOR && RESTICK_TO_NUCLEATOR) // re-stick to nucleator if come off
	{
		node_world.stucktonucleator = true;
		node_world.nucleator_stuck_position = node_world.pos_in_nuc_frame; // link to point *on* the nucleator surface
	}

    vect node_disp = node_world.pos_in_nuc_frame - oldpos;  // check now much we have moved
    
    node_world.nucleator_impacts += node_disp;

#ifndef SEED_INSIDE

    node_world.pos_in_nuc_frame -= node_disp / inertia;  // move the node *back* scaled by movibility (since nucleator should move this much)

    move_nuc(oldpos,node_disp);  // move the nucleator

#endif

    node_world.x = node_world.pos_in_nuc_frame.x; // set the position of the node
    node_world.y = node_world.pos_in_nuc_frame.y;
    node_world.z = node_world.pos_in_nuc_frame.z;

    ptheactin->nuc_to_world_frame(node_world);  // move node back to the world frame of ref

    // check if moved too much
    if ((fabs(node_disp.x) > 0.2*RADIUS) ||
	    (fabs(node_disp.y) > 0.2*RADIUS) ||
	    (fabs(node_disp.z) > 0.2*RADIUS*ELLIPSOID_STRETCHFACTOR))
    {
		cout << setprecision(3);
        cout << "Warning: Node " << node_world.nodenum << " nucleus ejection too great: " << node_disp.length() << " uM" << endl;
	    return false;  // failed (node still moved, though)
    }

	return true; // sucessful node ejection
}


void nucleator::move_nuc(const vect& origin_of_movement, const vect& tomove)
{
	/// updates the vectors for moving and rotating the nucleator

	/// Takes tomove as the movement vector from the origin_of_movement
	/// doesn't actually move anything here---that is done in actin::applyforces()

    // displacement is easy:
    deltanucposn -= tomove / inertia;

    // rotate the nucleator

    // calculate torque (as pos x disp)

    if (ROTATION)
    {
        vect lever_arm = origin_of_movement - centerofmass;
        torque += lever_arm.cross(tomove); // todo: is this right? or should it be scaled by movilility? (yes?)
    }

}

int nucleator::save_data(ofstream &ostr) 
{
    ostr << NUCSHAPE << endl;
    ostr << RADIUS << endl;
    ostr << CAPSULE_HALF_LINEAR << endl;
    ostr << surf_area << endl;
    ostr << (1/inertia) << endl;
    ostr << position << endl;
    ostr << deltanucposn_sum << endl;
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
    //if(geom == 0)
	//geometry= sphere;
    //else
	//geometry = capsule;

    double movability;
    
    istr >> RADIUS;
    istr >> CAPSULE_HALF_LINEAR;
    istr >> surf_area;
    istr >> movability;
    istr >> position;
    istr >> deltanucposn_sum;
    istr >> deltanucposn;
    istr >> torque;
    istr >> centerofmass;
    istr >> momentofinertia;

    inertia = 1/movability;

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

	if (NUCSHAPE==sphere)
	{

	// sphere

		for (double theta=-PI+2*PI/pointdensity; theta<PI-2*PI/pointdensity; theta+=2*PI/pointdensity)
			for (double phi=-PI+2*PI/pointdensity; phi<PI-2*PI/pointdensity; phi+=2*PI/pointdensity)
			{
				
				r = RADIUS * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RADIUS;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
			}


        cagepoints.push_back(vect(0,0,RADIUS));
        cagepoints.push_back(vect(0,0,-RADIUS));

	}
	else if (NUCSHAPE==ellipsoid)
	{

	// sphere

		for (double theta=-PI+2*PI/pointdensity; theta<PI-2*PI/pointdensity; theta+=2*PI/pointdensity)
			for (double phi=-PI+2*PI/pointdensity; phi<PI-2*PI/pointdensity; phi+=(2*PI/pointdensity) / ELLIPSOID_STRETCHFACTOR)
			{
				
				r = RADIUS * cos(phi);		// radius of circle
				
				xx = r * cos(theta);				// x and y of point
				yy = r * sin(theta);
				zz = sin(phi) * RADIUS * ELLIPSOID_STRETCHFACTOR;						// z just scaled by radius

				cagepoints.push_back(vect(xx,yy,zz));
			}


        cagepoints.push_back(vect(0,0,RADIUS));
        cagepoints.push_back(vect(0,0,-RADIUS));

	}
	else
	{


		for (double theta=-PI+2*PI/pointdensity; theta<PI-2*PI/pointdensity; theta+=2*PI/pointdensity)
			for (double phi=-PI+2*PI/pointdensity; phi<PI-2*PI/pointdensity; phi+=2*PI/pointdensity)
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

        cagepoints.push_back(vect(0,0,CAPSULE_HALF_LINEAR+RADIUS));
        cagepoints.push_back(vect(0,0,-CAPSULE_HALF_LINEAR-RADIUS));


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
