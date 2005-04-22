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

class nodes
{
public:
	MYDOUBLE x , y , z ;
	bool polymer;
	nodes(void);
	~nodes(void);
	nodes(MYDOUBLE set_x, MYDOUBLE set_y, MYDOUBLE set_z);
	nodes* nextnode;
	nodes* prevnode;
	bool depolymerize(void);
	bool polymerize(MYDOUBLE set_x, MYDOUBLE set_y, MYDOUBLE set_z);
	int save(ofstream*);
	int savedata(ofstream*);
	int loaddata(ifstream *inputstream);
	//int xrank;
	//int yrank;
	//int zrank;
	vect lastpos;
	vector <vect> link_force_vec;
	vector <vect> repulsion_displacement_vec;
	vector <vect> rep_force_vec;
	//vect momentum_vec;
	vector <links> listoflinks;
	int applyforces(int threadnum);
	int gridx, gridy, gridz;
	int updategrid(void);
	int removefromgrid(void);
	int addtogrid(void);
	int setgridcoords(void);
	int nodenum;
	int	nodelinksbroken;
	int addlink(nodes* linkto, MYDOUBLE dist);
	int removelink(nodes* link);
	Colour colour;
	actin *ptheactin;
	int creation_iter_num;
	bool harbinger;
	MYDOUBLE theta;
	MYDOUBLE phi;
	int savelinks(ofstream * outstream);
};

#endif
