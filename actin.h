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

class testnodeinfo
{
public:
    testnodeinfo(){};
    ~testnodeinfo(){};
    testnodeinfo(nodes *const & n, const vect& p):nodeptr(n),origunitvec(p){}  //nodes *

    nodes* nodeptr;
    vect origunitvec;
};

class tracknodeinfo
{
public:
    tracknodeinfo(){};
    ~tracknodeinfo(){};
    //tracknodeinfo(int const & setnodenum, const vect & setposn, const vect & setnucposn, const rotationmatrix & setrotation, int const & setframe ):
    tracknodeinfo(int const & setnodenum, const int & setx, const int & sety, int const & setframe ):
    nodenum(setnodenum),frame(setframe),x(setx),y(sety)
    {
        //posn = setposn;
        //nucposn = setnucposn;
        //rotation = setrotation;
    }  

    int nodenum;
    //vect posn, nucposn;
    //rotationmatrix rotation;
    int frame;
    int x,y;
};

class linkform
{
public:
	linkform(){};
    linkform(const int& nn, const double& dsqr):nodenum(nn),distsqr(dsqr){}

	static int CompareDistance (const linkform &elem1, const linkform &elem2)
	{
		return (elem1.distsqr < elem2.distsqr);
	}

	int nodenum;
	double distsqr;
};

/// One actin object is created per program run and it contains the nodes which make up the network, and associated functions.

/// The actin object is created in main().
/// etc etc

class actin
{
public:
	
	enum processfgbg 
	{
		runfg = 0,
		runbg = 1
	};

	
	actin(void);
	~actin(void);
	int saveinfo(void);
	void iterate(void);
	//int addlinks(const int& linknode1,const int& linknode2);
	bool addlinks(nodes& linknode1, nodes& linknode2) const;
	int savevrml(int filenum);
	
	ofstream opruninfo;
	ofstream opvelocityinfo, opinfo;
	ofstream outbmpfile_x,outbmpfile_y,outbmpfile_z;

    vector <int> nodes_to_track;
    vector <vector <tracknodeinfo> > node_tracks;

    static int lowestnodetoupdate;
	int highestnodecount;
	int lastsorthighestnode;

    int stationary_node_number;
    int stationary_node_xoffset, stationary_node_yoffset;;

	Dbl2d imageR, imageG, imageB;

    Dbl1d imageRmax, imageGmax, imageBmax;

    bool BMP_intensity_scaling;

	rotationmatrix world_to_nuc_rot, camera_rotation, camera_rotation2,
			reverse_camera_rotation;

	static rotationmatrix torque_rotate;
	static vect nuc_disp;

	//int num_rotate, num_displace;

	void crosslinknewnodes(const int &numnewnodes); 

	//vector <int_vect> nucleatorgrid;
	vector <int> crosslinknodesdelay;

    ofstream::pos_type bitmap_start;


	bool CompareDistance ( linkform* elem1, linkform* elem2 );
	
	void collisiondetection(void);
	void nucleator_node_interactions(void);
	//void move_and_rotate(void);

	void applyforces(void);
	void addapplyforcesthreads(const int threadnumoffset, const int &lowestnodenum, const int &highestnodenum);

    void set_nodes_to_track(const projection &  proj);
    void set_axisrotation(const projection &  proj, rotationmatrix & axisrotation);
	

	void writebitmapfile(ofstream& outbmpfile, const Dbl2d& imageR, const Dbl2d& imageG, const Dbl2d& imageB);
	void writebitmapheader(ofstream& outbmpfile, const int & bitmapwidth, const int & bitmapheight);
	
	void savebmp(const int &filenum, const projection & proj, const processfgbg& fgbg, bool writefile);



	nucleator* p_nuc;
	void linkforces();
	Colour newnodescolour;
	static int iteration_num;
	int symbreakiter;
	int linksformed;
	int linksbroken;

	bool currentlyusingthreads;

	char temp_BMP_filename_x[1024],
		 temp_BMP_filename_y[1024],
		 temp_BMP_filename_z[1024];

	//int doreportiteration;
	//int doreportmaxnodes;

	bool brokensymmetry;

	int attemptedpolrate, polrate;
	
	void setnodecols(void);

	vector <double> speckle_array;
	int speckle_array_size;



	static vector <nodes> node;
	static vector <bool> donenode;	
    static Nodes2d nodes_by_thread;
	static Nodes2d recti_near_nodes;
	static Nodes2d nodes_on_same_gridpoint;
	static vector <int> nearby_collision_gridpoint_offsets;
	static vector <int>::iterator nearby_collision_gridpoint_offset_begin;
	static vector <int>::iterator nearby_collision_gridpoint_offset_end;

    //static vector <vector <vector <nodes*>*>> gridpointsbythread;
    
    static vector<vector<NODEGRIDTYPE<nodes*>*> > gridpointsbythread;
    

    //static vector <int> recti_near_nodes_size;
    //static vector <int> nodes_on_same_gridpoint_size;

	//static Nodes1d nodes_within_nucleator;	
	static size_t findnearbynodes(const nodes& ournode, const int& adjgridpoints, const int& threadnum);
    //static size_t findnearbynodes_collision(const nodes& ournode, const int& threadnum);

	//void findnearbynodes_collision_setup(const int& adjgridpoints);

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

    vector <vector <testnodeinfo> > testnodes;
    double testsurfaceposn, lasttestsurfaceposn, lasttestsurfacesavedposn;
    double testsurfacerotation;
    double testforcemag;
    bool test_equilibrating;
    vect testdirection;     // defines the direction of the test surface
    double testangle;       // defines the arc

    void testforces_setup();
    void testforces_set_initial_surfaces();
    void testforces_cutlinks();   
    void testforces_remove_nontest_nodes();

    void testforces_addforces(const int &surface);
    void testforces_select_nodes(const double& testdist, const short int &setsurface);
    void testforces_saveiter();

    double sum_delta_movements();

	void squash(const double & thickness);
	void sortnodesbygridpoint(void);
    void sortgridpointsbythread(void);
	int nexttocrosslink;
	int find_center(vect &center);

    int currentsmallestgridthread;

	void clear_nodegrid();
	int save_data(ofstream &ofstrm);
	bool load_data(ifstream &ifstrm);
    void rebuildnodepointers();
	void setdontupdates(void);
	void set_sym_break_axes(void);
	bool load_sym_break_axes(void);
	void save_sym_break_axes(void);
	void clear_node_stats(void);

    void addbrownianforces();

    rotationmatrix nuc_to_world_rot;

	inline int pixels(const double & coord) const
	{  // convert simulation distance into pixel distance
        return int(coord > 0.0 ? dbl_pixels(coord) + 0.5 :
                                 dbl_pixels(coord) - 0.5);
	}

	inline double unpixels(const int & pix) const
	{  // convert pixel distance into simulation distance
		return ((double) pix / (double) BMP_HEIGHT) * (double) VIEW_HEIGHT;
	}

	inline double dbl_pixels(const double & coord) const
	{  // convert simulation distance into pixel distance
		return ((double) BMP_HEIGHT * ( (coord)/VIEW_HEIGHT) ); 
	}

    // convert positions from one frame of ref to the other

    inline void world_to_nuc_frame(vect & v)
    {
        world_to_nuc_rot.rotate(v); 
        v -= p_nuc->position;
    }

    inline void nuc_to_world_frame(vect & v)
    {
        v += p_nuc->position;
        nuc_to_world_rot.rotate(v); 
    }
	

    void reservemorenodes(const int extranodes);
};

#endif
