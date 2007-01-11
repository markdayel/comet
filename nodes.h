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

#ifndef nodes_H
#define nodes_H

//#include "stdafx.h"
#include "vect.h"
#include "Colour.h"
#include "links.h"

class links;
class Colour;

class nodes: public vect
{
public:

    NODEGRIDTYPE<nodes*>* nodegridptr;
	
	bool onseg;
    bool polymer;
    bool harbinger;
	bool stucktonucleator;
	bool move_harbinger_this_time;
    bool testnode;

    int nodenum;
    int creation_iter_num;
    

    short int testsurface;  // surface number that node is attached to

    short int threadnum;    

	short int gridx, gridy, gridz;

	//actin* ptheactin;
	//nodes* nextnode;
	//nodes* prevnode;


	vect unit_vec_posn;  // this is kept up-to-date in the updategrid() function
	//vect nearest_surface_point;

	vect link_force_vec;  
	
	vect rep_force_vec;
	vect nuc_repulsion_displacement_vec;
 	vect nucleator_stuck_position;
	vect delta;
	vect viscosity_velocity_sum;

    vect nucleator_link_force;	// force on nodes by the link to the nucleator 
								// again merely for display purposes
								// segmented in segments.cpp
								// nucleator attachements are dealt with in actin::nucleator_node_interactions()

    double pressure;

	double linkforce_transverse, linkforce_radial,  
			 repforce_transverse, repforce_radial,
			 links_broken;

    double linkenergy_transverse, linkenergy_radial,
			 repenergy_transverse, repenergy_radial,
			 linksenergy_broken;


	
	double nucleator_impacts;	// keeps track of ejections of node from nucleator
								// merely for display of ejection forces
								// segmented in segments.cpp
								// not used in moving node, nucleator etc.
								// that is done directly at the mo.

	//Colour colour;

	
	double dist_from_surface;

	double viscosity_velocity_unweight;

    vector <links> listoflinks;



    int savelinks(ofstream * outstream);
    
	nodes(void);
	~nodes(void);
	nodes(const double& set_x, const double& set_y,const double& set_z);

	bool depolymerize(void);
	bool polymerize(const double& set_x, const double& set_y, const double& set_z);
	int save_data(ofstream &ofstr);
	int load_data(ifstream &ifstr);

	void addlink(nodes& linkto, const double& dist);
	void removelink(const nodes* linkednode);

	void updategrid(void);
	void removefromgrid(void);
	void addtogrid();
	void setgridcoords(void);

    vect posnoflastgridupdate;

	inline void applyforces() 
	{	

	    if (VISCOSITY)
	    {
		    // simple average weighted by VISCOSITY_FACTOR
		    //
		    //delta = (delta + (viscosity_velocity_sum * VISCOSITY_FACTOR )) / 
		    //	                     ( 1 + VISCOSITY_FACTOR + mymax(VISCOSITY_EDGE_FACTOR, viscosity_velocity_unweight));

		    delta = ( (link_force_vec + rep_force_vec ) * NODE_FORCE_TO_DIST * NON_VISC_WEIGHTING + 
			        (viscosity_velocity_sum	* VISCOSITY_FACTOR )) / 
			                         ( NON_VISC_WEIGHTING + VISCOSITY_FACTOR * mymax(VISCOSITY_EDGE_FACTOR, viscosity_velocity_unweight));
        } else
        {
	        delta = (link_force_vec + rep_force_vec)				
                    * NODE_FORCE_TO_DIST;
        }

        *this += delta;

        setunitvec();  // we've moved the node, so reset the unit vector
        clearforces();   // and clear force sums we just used

        // note: remember to do this clearing too when switch clearing harbinger flag
        //       in crosslinknewnodes()
    }

	inline void clearforces()
	{
		
	    if (VISCOSITY)
	    {
		    viscosity_velocity_unweight = 0.0;
		    viscosity_velocity_sum.zero();
	    }

	    rep_force_vec.zero();
	    link_force_vec.zero();

	    pressure = 0;

	}



	inline void getdirectionalmags(const vect &displacement, double &dotmag, double &crossmag) const
	{
		dotmag = fabs(unit_vec_posn.dot(displacement));
		crossmag = displacement.length() - dotmag;
	}

	inline void adddirectionalmags(const vect &displacement, double &dotmag, double &crossmag) const
	{  
		double tmp_dotmag; 
		
		tmp_dotmag = fabs(unit_vec_posn.dot(displacement));

		dotmag   += tmp_dotmag;
		crossmag += displacement.length() - tmp_dotmag;
	}

	inline void clearstats()
	{
		// these stats are built up over many iterations, and cleared every readout frame

		linkforce_transverse = 
		linkforce_radial     = 
		repforce_transverse  = 
		repforce_radial      = 
		links_broken	     =
 		nucleator_impacts    = 0.0;

		nucleator_link_force.zero();

  //      linkenergy_transverse = 
  //      linkenergy_radial     =
		//repenergy_transverse  =
  //      repenergy_radial      =
		//linksenergy_broken    = 0.0;

	}

	void setunitvec(void)
	{	

		if (NUCSHAPE == nucleator::sphere)
		{
            dist_from_surface = this->length();	 // not really dist_from_surface yet, need to subtract radius
            unit_vec_posn = *this * (1/dist_from_surface);  // set unit vector position
			//nearest_surface_point = unit_vec_posn
			dist_from_surface -= RADIUS;

		}
		else
		{	// capsule
			if (fabs(z) < CAPSULE_HALF_LINEAR)
			{  // on cylinder, no z component

				dist_from_surface = calcdist(x,y);   // not really dist_from_surface yet, need to subtract radius
				unit_vec_posn = vect(x/dist_from_surface, y/dist_from_surface, 0);
				dist_from_surface -= RADIUS;
				//nearest_surface_point = unit_vec_posn + vect(0,0,z);
				onseg = true;

			}
			else
			{	// on ends

				onseg = false;

				if (z>0) // top
				{
					vect offsetvec = *this;
					offsetvec.z -= CAPSULE_HALF_LINEAR;

                    dist_from_surface = offsetvec.length();	// not really dist_from_surface yet, need to subtract radius

                    unit_vec_posn = offsetvec * (1/dist_from_surface);  // set unit vector position

					//nearest_surface_point = unit_vec_posn + vect(0,0,CAPSULE_HALF_LINEAR);

					dist_from_surface -= RADIUS;

				}
				else
				{
					vect offsetvec = *this;
					offsetvec.z += CAPSULE_HALF_LINEAR;

                    dist_from_surface = offsetvec.length();	 // not really dist_from_surface yet, need to subtract radius

                    unit_vec_posn = offsetvec * (1/dist_from_surface);  // set unit vector position

					//nearest_surface_point = unit_vec_posn + vect(0,0,CAPSULE_HALF_LINEAR);

					dist_from_surface -= RADIUS;

				}	
			}
		}
	}
};

#endif
