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

#include "stdafx.h"
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
    vect viscous_force_vec;
    double viscous_force_recip_dist_sum;
	vector <links> listoflinks;

	double linkforce_transverse, linkforce_radial,  // index is the threadnum
			 repforce_transverse, repforce_radial,
			 links_broken;
	
	double nucleator_impacts;
    vect nucleator_link_force;

    bool insidenucleator;

	int gridx, gridy, gridz;
	vect delta;
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

    bool stucktonucleator;
    vect nucleator_stuck_position;
    

    //bool dontupdate;

    //typedef void (nodes::*p_applyforces_fn)();

    //void (nodes::*p_applyforces_fn)(void);
    
	void applyforces_novisc();
    void applyforces_visc();

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
	}

	void setunitvec(void)
	{	

		if (NUCSHAPE == nucleator::sphere)
		{
            double len = this->length();
            unit_vec_posn = *this * (1/len);  // set unit vector position
            if (len < RADIUS)
                insidenucleator = true;
            else
                insidenucleator = false;
		}
		else
		{	// capsule
			if (fabs(z) < CAPSULE_HALF_LINEAR)
			{  // on cylinder, no z component
				double len = calcdist(x,y);
				unit_vec_posn = vect(x/len, y/len, 0);

				onseg = true;

                if (len < RADIUS)
                    insidenucleator = true;
                else
                    insidenucleator = false;
			}
			else
			{	// on ends

				onseg = false;

				if (z>0) // top
				{
					vect offsetvec = *this;
					offsetvec.z -= CAPSULE_HALF_LINEAR;

                    double len = offsetvec.length();

                    unit_vec_posn = offsetvec * (1/len);  // set unit vector position

					if (len < RADIUS)
                        insidenucleator = true;
                    else
                        insidenucleator = false;
				}
				else
				{
					vect offsetvec = *this;
					offsetvec.z += CAPSULE_HALF_LINEAR;

                    double len = offsetvec.length();

                    unit_vec_posn = offsetvec * (1/len);  // set unit vector position

					if (len < RADIUS)
                        insidenucleator = true;
                    else
                        insidenucleator = false;
				}	
			}
		}
	}


    
};

#endif
