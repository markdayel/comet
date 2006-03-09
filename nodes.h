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

    bool polymer;
	nodes(void);
	~nodes(void);
	nodes(const double& set_x,const double& set_y,const double& set_z);
	nodes* nextnode;
	nodes* prevnode;
	bool depolymerize(void);
	bool polymerize(const double& set_x, const double& set_y, const double& set_z);
	int save_data(ofstream &ofstr);
	int load_data(ifstream &ifstr);

	vect unit_vec_posn;  // this is kept up-to-date in the updategrid() function

	vect link_force_vec;  
	
	vect rep_force_vec;
	vect nuc_repulsion_displacement_vec;
    //vect viscous_force_vec;
    //double viscous_force_recip_dist_sum;
	vector <links> listoflinks;

	double linkforce_transverse, linkforce_radial,  // index is the threadnum
			 repforce_transverse, repforce_radial,
			 links_broken;
	
	double nucleator_impacts;	// keeps track of ejections of node from nucleator
								// merely for display of ejection forces
								// segmented in segments.cpp
								// not used in moving node, nucleator etc.
								// that is done directly at the mo.

    vect nucleator_link_force;	// force on nodes by the link to the nucleator 
								// again merely for display purposes
								// segmented in segments.cpp
								// nucleator attachements are dealt with in actin::nucleator_node_interactions()

    //bool insidenucleator;

	int threadnum;    

	int gridx, gridy, gridz;
	//vect delta;
	void updategrid(void);
	void removefromgrid(void);
	void addtogrid(void);
	void setgridcoords(void);
	int nodenum;
	int	nodelinksbroken;
	int addlink(nodes* linkto, const double& dist);
	int removelink(nodes* link);
	Colour colour;
	actin *ptheactin;
	int creation_iter_num;
	bool harbinger;
	double theta;
	double phi;
	int savelinks(ofstream * outstream);
	bool onseg;
	double dist_from_surface;

    bool stucktonucleator;
    vect nucleator_stuck_position;

	bool move_harbinger_this_time;
    

    //bool dontupdate;

    //typedef void (nodes::*p_applyforces_fn)();

    //void (nodes::*p_applyforces_fn)(void);
    
	void applyforces();
    //void applyforces_visc();

	void getdirectionalmags(const vect &displacement, double &dotmag, double &crossmag)
	{
		dotmag = fabs(unit_vec_posn.dot(displacement));
		crossmag = displacement.length() - dotmag;
	}

	void adddirectionalmags(const vect &displacement, double &dotmag, double &crossmag)
	{  
		double tmp_dotmag; 
		
		tmp_dotmag = fabs(unit_vec_posn.dot(displacement));

		dotmag   += tmp_dotmag;
		crossmag += displacement.length() - tmp_dotmag;
	}

	inline void clearstats()
	{
		linkforce_transverse = 
		linkforce_radial     = 
		repforce_transverse  = 
		repforce_radial      = 
		links_broken	     =
 		nucleator_impacts    = 0;

		nucleator_link_force.zero();

	}

	void setunitvec(void)
	{	

		if (NUCSHAPE == nucleator::sphere)
		{
            dist_from_surface = this->length();	 // not really dist_from_surface yet, need to subtract radius
            unit_vec_posn = *this * (1/dist_from_surface);  // set unit vector position
			dist_from_surface -= RADIUS;

		}
		else
		{	// capsule
			if (fabs(z) < CAPSULE_HALF_LINEAR)
			{  // on cylinder, no z component

				dist_from_surface = calcdist(x,y);   // not really dist_from_surface yet, need to subtract radius
				unit_vec_posn = vect(x/dist_from_surface, y/dist_from_surface, 0);
				dist_from_surface -= RADIUS;

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


					dist_from_surface -= RADIUS;

				}
				else
				{
					vect offsetvec = *this;
					offsetvec.z += CAPSULE_HALF_LINEAR;

                    dist_from_surface = offsetvec.length();	 // not really dist_from_surface yet, need to subtract radius

                    unit_vec_posn = offsetvec * (1/dist_from_surface);  // set unit vector position

					dist_from_surface -= RADIUS;

				}	
			}
		}
	}
};

#endif
