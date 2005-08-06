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
	//MYDOUBLE x , y , z ;
	bool polymer;
	nodes(void);
	~nodes(void);
	nodes(const MYDOUBLE& set_x,const MYDOUBLE& set_y,const MYDOUBLE& set_z);
	nodes* nextnode;
	nodes* prevnode;
	bool depolymerize(void);
	bool polymerize(const MYDOUBLE& set_x, const MYDOUBLE& set_y, const MYDOUBLE& set_z);
	int save(ofstream*);
	int savedata(ofstream*);
	int loaddata(ifstream *inputstream);
	int save_data(ofstream &ofstr);
	int load_data(ifstream &ifstr);
	//int xrank;
	//int yrank;
	//int zrank;
	vect lastpos;
	vect unit_vec_posn;  // this is kept up-to-date in the updategrid() function
	vector <vect> link_force_vec;  // index is the threadnum
	vector <vect> repulsion_displacement_vec;
	vector <vect> rep_force_vec;
	//vect momentum_vec;
	vector <links> listoflinks;

	vector <MYDOUBLE> linkforce_transverse, linkforce_radial,
			 repforce_transverse, repforce_radial,
			 dispforce_transverse, dispforce_radial;  // index is the threadnum
	
	int gridx, gridy, gridz;
	vect delta;
	//MYDOUBLE delta_x, delta_y, delta_z;
	void updategrid(void);
	void removefromgrid(void);
	void addtogrid(void);
	void setgridcoords(void);
	int nodenum;
	int	nodelinksbroken;
	int addlink(nodes* linkto, const MYDOUBLE& dist);
	int removelink(nodes* link);
	Colour colour;
	actin *ptheactin;
	int creation_iter_num;
	bool harbinger;
	MYDOUBLE theta;
	MYDOUBLE phi;
	int savelinks(ofstream * outstream);
	bool dontupdate;

	//inline int applyforces(const int &threadnum);
	inline void applyforces(const int &threadnum)
	{	
		delta = (link_force_vec[threadnum] + rep_force_vec[threadnum]) * DELTA_T * FORCE_SCALE_FACT;

		*this+=delta;

		rep_force_vec[threadnum].zero();
		link_force_vec[threadnum].zero();
		repulsion_displacement_vec[threadnum].zero();

		//return 0;
	}

	inline void getdirectionalmags(const vect &displacement, MYDOUBLE &dotmag, MYDOUBLE &crossmag)
	{
		dotmag = fabs(unit_vec_posn.dot(displacement));
		crossmag = displacement.length() - dotmag;
	}

	inline void adddirectionalmags(const vect &displacement, MYDOUBLE &dotmag, MYDOUBLE &crossmag)
	{  
		MYDOUBLE tmp_dotmag = fabs(unit_vec_posn.dot(displacement));

		dotmag += tmp_dotmag;
		crossmag += displacement.length() - tmp_dotmag;
	}

	inline void clearstats(const int &threadnum)
	{
		linkforce_transverse[threadnum] = 
		linkforce_radial[threadnum]     = 
		repforce_transverse[threadnum]  = 
		repforce_radial[threadnum]      = 
		dispforce_transverse[threadnum] = 
		dispforce_radial[threadnum]     = 0;
	}
};

#endif
