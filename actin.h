/*
Copyright (C) 2005 Mark J Dayel

You may not distribute, deploy, provide copies of, rent, loan, 
lease, transfer, or grant any rights in the software or derivative 
works thereof in any form to any person.  Reproduction, adaptation, 
or translation of this program without prior written permission 
from the author is prohibited.  Proofs of infringement of copyright 
violations include but not limited to similar code style and 
structures, similar layout and design, similar algorithm design, 
and containing parts of the original software source code.  
Copyright notice must remain intact and cannot be removed without 
prior written permission from the author.
*/

#ifndef actin_H 
#define actin_H

#include "nodes.h"
#include <vector>
#include "nucleator.h"
#include "Colour.h"
//#include "stdafx.h"

typedef vector<nodes*> Nodes1d;
typedef vector<Nodes1d> Nodes2d;
typedef vector<Nodes2d> Nodes3d;

//typedef vector<signed char> Bool1d;
//typedef vector<Bool1d> Bool2d;

typedef vector<MYDOUBLE> Dbl1d;
typedef vector<Dbl1d> Dbl2d;

class nodes;
class links;

struct int_vect 
{
	int x,y,z;
};

class linkform
{
public:
	linkform(){};
	linkform(int nn, MYDOUBLE dsqr)
	{
		nodenum = nn;
		distsqr = dsqr;
	}

	static int CompareDistance ( linkform elem1, linkform elem2 )
	{
		return (elem1.distsqr < elem2.distsqr);
	}

	int nodenum;
	MYDOUBLE distsqr;
};

class actin
{
public:
	enum projection 
	{
		xaxis,
		yaxis,
		zaxis
	};
	actin(void);
	~actin(void);
	int nucleate();
	int save(int);
	int saveinfo(void);
	int iterate(void);
	int addlinks(const int& linknode1,const int& linknode2);
	int savevrml(int filenum);
	ofstream opruninfo;
	ofstream opvelocityinfo;
	int highestnodecount;
	Dbl2d imageR, imageG, imageB;

	vector <Dbl2d> reportdat;

	static bool collisionthreaddone1;
	static bool collisionthreaddone2;
	static bool collisionthreaddone3;
	static bool collisionthreaddone4;

	int crosslinknewnodes(int numnewnodes);

	//links *link;
	//vector <nodes*> nodes;
	vector <int_vect> nucleatorgrid;
	vector <int> crosslinknodesdelay;
	vector <MYDOUBLE> link_radial_histogram;
	vector <MYDOUBLE> link_transverse_histogram;
	bool CompareDistance ( linkform* elem1, linkform* elem2 );
	
	int collisiondetection(void);
	int ejectfromnucleator(void);

	int applyforces(void);
	int savebmp(int filenum, projection proj);

	inline MYDOUBLE square(const MYDOUBLE& numb)
	{
		return numb*numb;
	}

	nucleator* nucleation_object;
	int linkforces(const bool& sumforces);
	Colour newnodescolour;
	static int iteration_num;
	int linksformed;
	int linksbroken;

	int doreportiteration;
	int doreportmaxnodes;
	
	int setnodecols(void);

	static vector <nodes> node;
	static vector <bool> donenode;	
	static vector <bool> repdonenode;
	static Nodes2d recti_near_nodes;
	static Nodes2d nodes_on_same_gridpoint;
	static Nodes1d nodes_within_nucleator;	
	static vector <int> nodesbygridpoint;
	inline static int findnearbynodes(const int& ournodenum,const int& adjgridpoints,const int& threadnum);
	inline static int dorepulsion(const int& node_i,const int& node_j, const MYDOUBLE& distsqr, const int& threadnum);
	static void *collisiondetectionthread(void*threadarg);
	inline static void *collisiondetectiondowork(thread_data* dat);
	static bool isinthread;
	//static Bool2d repulsedone;

	static void *applyforcesthread(void* threadarg);
	static vector <nodes*> linkremovefrom;
	static vector <nodes*> linkremoveto;
	//vector <int> linkformfrom;
	vector <linkform> linkformto;
	static void *linkforcesthread(void* threadarg);
	static void *repulsiveforcesthread(void* threadarg);
	static void *compressfilesthread(void* threadarg);

	int squash(MYDOUBLE thickness);
	//int repulsiveforces(void);
	int sortnodesbygridpoint(void);
	int nexttocrosslink;
	int find_centre(MYDOUBLE &centre_x, MYDOUBLE &centre_y, MYDOUBLE &centre_z);
	int save_linkstats(int filenum);
	void reportsnapshot(int filenum, int highestnode, int reportiteration);
	void savereport(int filenum, int highestnode);
	int savedata(int filenum);
	int loaddata(int filenum);
	void clear_nodegrid();
	int save_data(ofstream &ofstrm);
	int load_data(ifstream &ifstrm);
	void setdontupdates(void);
};

#endif
