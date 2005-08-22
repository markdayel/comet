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

class links;
class Colour;

class nodes: public vect
{
public:
	//double x , y , z ;
	bool polymer;
	nodes(void);
	~nodes(void);
	nodes(const double& set_x,const double& set_y,const double& set_z);
	nodes* nextnode;
	nodes* prevnode;
	bool depolymerize(void);
	bool polymerize(const double& set_x, const double& set_y, const double& set_z);
	//int save(ofstream*);
	//int savedata(ofstream*);
	//int loaddata(ifstream *inputstream);
	int save_data(ofstream &ofstr);
	int load_data(ifstream &ifstr);
	//int xrank;
	//int yrank;
	//int zrank;
	//vect lastpos;
	vect unit_vec_posn;  // this is kept up-to-date in the updategrid() function
	vector <vect> link_force_vec;  // index is the threadnum
	//vector <vect> repulsion_displacement_vec;
	vector <vect> rep_force_vec;
	//vect momentum_vec;
	vector <links> listoflinks;

	vector <double> linkforce_transverse, linkforce_radial,  // index is the threadnum
			 repforce_transverse, repforce_radial,
			 links_broken;
			 //dispforce_transverse, dispforce_radial;  
	
	double nucleator_impacts;

	int gridx, gridy, gridz;
	vect delta;
	//double delta_x, delta_y, delta_z;
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
	bool dontupdate;

	//inline int applyforces(const int &threadnum);
	inline void applyforces(const int &threadnum)
	{	
		delta = (link_force_vec[threadnum] + rep_force_vec[threadnum]) * DELTA_T * FORCE_SCALE_FACT;
				//+ repulsion_displacement_vec[threadnum];

		*this+=delta;

		rep_force_vec[threadnum].zero();
		link_force_vec[threadnum].zero();
		//repulsion_displacement_vec[threadnum].zero();

		//return 0;
	}

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

	inline void clearstats(const int &threadnum)
	{
		linkforce_transverse[threadnum] = 
		linkforce_radial[threadnum]     = 
		repforce_transverse[threadnum]  = 
		repforce_radial[threadnum]      = 
		links_broken[threadnum]			= 0;

		nucleator_impacts  = 0;
//		dispforce_transverse[threadnum] = 
//		dispforce_radial[threadnum]     = 0;

		
	}

	void setunitvec(void)
	{	// TODO: fix this to bring in line with normal nuclator shape test

		if (NUCSHAPE == nucleator::sphere)
		{
			unit_vec_posn = this->unitvec();  // set unit vector position
		}
		else
		{	// capsule
			if (fabs(z) < CAPSULE_HALF_LINEAR)
			{  // on cylinder, no z component
				double len = calcdist(x,y);
				unit_vec_posn = vect(x/len, y/len, 0);
				onseg = true;
			}
			else
			{	// on ends

				onseg = false;

				if (z>0) // top
				{
					vect offsetvec = *this;
					offsetvec.z -= CAPSULE_HALF_LINEAR;
					unit_vec_posn = offsetvec.unitvec();
				}
				else
				{
					vect offsetvec = *this;
					offsetvec.z += CAPSULE_HALF_LINEAR;
					unit_vec_posn = offsetvec.unitvec();
				}	
			}
		}
	}
	bool onseg;
};

#endif
