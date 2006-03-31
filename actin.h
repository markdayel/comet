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

#include "rotationmatrix.h"
#include <vector>
#include <sstream>
#include "nodes.h"
#include "nucleator.h"
#include "Colour.h"
#include "threadedtaskqueue.h"
//#include <iterator>

//#include "stdafx.h"

typedef vector<nodes*> Nodes1d;
typedef vector<Nodes1d> Nodes2d;
typedef vector<Nodes2d> Nodes3d;

//typedef vector<signed char> Bool1d;
//typedef vector<Bool1d> Bool2d;

typedef vector<double> Dbl1d;
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
	linkform(int nn, double dsqr)
	{
		nodenum = nn;
		distsqr = dsqr;
	}

	static int CompareDistance ( linkform elem1, linkform elem2 )
	{
		return (elem1.distsqr < elem2.distsqr);
	}

	int nodenum;
	double distsqr;
};

class actin
{
public:
	enum projection 
	{
		xaxis = 0,
		yaxis = 1,
		zaxis = 2
	};
	enum processfgbg 
	{
		runfg = 0,
		runbg = 1
	};

	
	actin(void);
	~actin(void);
	int saveinfo(void);
	int iterate(void);
	int addlinks(const int& linknode1,const int& linknode2);
	int savevrml(int filenum);
	
	ofstream opruninfo;
	ofstream opvelocityinfo;
	ofstream outbmpfile_x,outbmpfile_y,outbmpfile_z;

    

    int lowestnodetoupdate;
	int highestnodecount;
	int lastsorthighestnode;

	Dbl2d imageR, imageG, imageB;

    Dbl1d imageRmax, imageGmax, imageBmax;

    bool BMP_intensity_scaling;

	rotationmatrix actin_rotation, camera_rotation, camera_rotation2,
			reverse_camera_rotation;

	//int num_rotate, num_displace;

	void crosslinknewnodes(const int &numnewnodes); 

	//vector <int_vect> nucleatorgrid;
	vector <int> crosslinknodesdelay;
	bool CompareDistance ( linkform* elem1, linkform* elem2 );
	
	int collisiondetection(void);
	void nucleator_node_interactions(void);
	void move_and_rotate(void);

	int applyforces(void);	

	ofstream::pos_type bitmap_start;

	void writebitmapfile(ofstream& outbmpfile, const Dbl2d& imageR, const Dbl2d& imageG, const Dbl2d& imageB);
	void writebitmapheader(ofstream& outbmpfile, const int & bitmapwidth, const int & bitmapheight);
	
	void savebmp(const int &filenum, projection proj, processfgbg fgbg, bool writefile);

	nucleator* p_nuc;
	int linkforces();
	Colour newnodescolour;
	static int iteration_num;
	int symbreakiter;
	int linksformed;
	int linksbroken;

	char temp_BMP_filename_x[255],
		 temp_BMP_filename_y[255],
		 temp_BMP_filename_z[255];

	int doreportiteration;
	int doreportmaxnodes;

	bool brokensymmetry;

	int attemptedpolrate, polrate;
	
	int setnodecols(void);

	vector <double> speckle_array;
	int speckle_array_size;

//#ifdef NODE_GRID_USE_ARRAYS
//	static nodes** __restrict nodegrid;
//#else
//	static Nodes3d nodegrid;
//#endif

	static vector <nodes> node;
	static vector <bool> donenode;	
    static Nodes2d nodes_by_thread;
	static Nodes2d recti_near_nodes;
	static Nodes2d nodes_on_same_gridpoint;
	static vector <int> nearby_collision_gridpoint_offsets;
	static vector <int>::iterator offset_begin;
	static vector <int>::iterator offset_end;

    //static vector <int> recti_near_nodes_size;
    //static vector <int> nodes_on_same_gridpoint_size;

	//static Nodes1d nodes_within_nucleator;	
	static int findnearbynodes(const nodes& ournode, const int& adjgridpoints, const int& threadnum);
    static int findnearbynodes_collision(const nodes& ournode, const int& threadnum);

	void findnearbynodes_collision_setup(const int& adjgridpoints);

    // -- Threading, comment these out
	// static void *collisiondetectionthread(void* threadarg);
	// static void *collisiondetectiondowork(thread_data* dat);
        // --
	static bool isinthread;
	// -- Threading add worker functions
	static void *collisiondetectiondowork(void* arg);//, pthread_mutex_t *mutex);
    //static void *collisiondetectiondoworkvisc(void* arg);//, pthread_mutex_t *mutex);
	static void *linkforcesdowork(void* arg);//, pthread_mutex_t *mutex);
	static void *applyforcesdowork(void* threadarg);//, pthread_mutex_t *mutex);
	// --
	//static Bool2d repulsedone;

	//static Node3d nodegrid;

	static void *applyforcesthread(void* threadarg);
	static Nodes2d linkremovefrom;
	static Nodes2d linkremoveto;
	vector <linkform> linkformto;
	static void *linkforcesthread(void* threadarg);
	//static void *compressfilesthread(void* threadarg);
	void compressfilesdowork(const int & filenum);   

	int squash(double thickness);
	void sortnodesbygridpoint(void);
	int nexttocrosslink;
	int find_center(vect &center);

	void clear_nodegrid();
	int save_data(ofstream &ofstrm);
	int load_data(ifstream &ifstrm);
	void setdontupdates(void);
	void set_sym_break_axes(void);
	bool load_sym_break_axes(void);
	void save_sym_break_axes(void);
	void clear_node_stats(void);

	inline int pixels(const double & coord) const
	{  // convert simulation distance into pixel distance
		return (int)((double) BMP_HEIGHT * ( (coord)/VIEW_HEIGHT) ); 
	}

	inline double dbl_pixels(const double & coord) const
	{  // convert simulation distance into pixel distance
		return ((double) BMP_HEIGHT * ( (coord)/VIEW_HEIGHT) ); 
	}
	
	inline double unpixels(const int & pix) const
	{  // convert pixel distance into simulation distance
		return ((double) pix / (double) BMP_HEIGHT) * (double) VIEW_HEIGHT;
	}

	//void removenodefromgrid(nodes* remnode)
	//{
	//	// are we on the grid?
	//	if (remnode->gridx==-1) return;  // return if not

	//	// are we the only grid node?
	//	if ((remnode->nextnode==remnode) &&
	//			(NODEGRID(remnode->gridx,remnode->gridy,remnode->gridz) == remnode))
	//		{	// if so, just delete the grid reference
	//			NODEGRID(remnode->gridx,remnode->gridy,remnode->gridz) = 0;
	//		}
	//		else
	//		{	// other nodes on grid
	//			if (NODEGRID(remnode->gridx,remnode->gridy,remnode->gridz) == remnode)
	//			{  // if we're the grid reference set ref to next node
	//				NODEGRID(remnode->gridx,remnode->gridy,remnode->gridz) = remnode->nextnode;
	//			}

	//			remnode->nextnode->prevnode = remnode->prevnode;  //  remove self from circular list
	//			remnode->prevnode->nextnode = remnode->nextnode;

	//		}

	//		remnode->gridx=remnode->gridy=remnode->gridz=-1;

	//	return;
	//}

	//void addnodetogrid(nodes* addnode)
	//{
	//	// are we already on the grid?
	//	//if (gridx!=-1) return 0;

	//	// is the new grid node empty?
	//	if ((NODEGRID(addnode->gridx,addnode->gridy,addnode->gridz) == 0))
	//	{	// if so, just add self to the grid reference
	//		NODEGRID(addnode->gridx,addnode->gridy,addnode->gridz) = addnode;
	//		addnode->nextnode = addnode->prevnode = addnode;  // and loop to self
	//	}
	//	else
	//	{	// otherwise sew into loop
	//		addnode->nextnode = NODEGRID(addnode->gridx,addnode->gridy,addnode->gridz);  // our next is the grid
	//		addnode->prevnode = addnode->nextnode->prevnode;  //our new previous is the new next's old previous

	//		addnode->nextnode->prevnode = addnode;  // and we are theirs
	//		addnode->prevnode->nextnode = addnode;
	//	}

	//return;
	//}


    void keep_mem_resident(void);
    void reservemorenodes(const int extranodes);
};

#endif
