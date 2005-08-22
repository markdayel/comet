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

#ifndef stdafx_H
#define stdafx_H

#ifdef _WIN32

	#pragma warning(disable: 4511) // unable to generate copy constructor
	#pragma warning(disable: 4512)

	#pragma warning(disable: 4127)  // constant conditional expression

	#pragma inline_depth( 255 )
	#pragma inline_recursion( on )
	#pragma auto_inline( on )
	
	//#define inline __forceinline

#endif


#ifndef _WIN32

	//#include <unistd.h>
	#define VRMLDIR "vrml/"
	#define DATADIR "data/"
	#define REPORTDIR "reports/"
	#define BITMAPDIR "bitmaps/"
	#define TEMPDIR "temp/"


#else
	
	#define VRMLDIR "vrml\\"
	#define DATADIR "data\\"
	#define REPORTDIR "reports\\"
	#define BITMAPDIR "bitmaps\\"
	#define TEMPDIR "temp\\"

#endif

// defines

// #define FORCES_BOTH_WAYS 1
// #define NO_CALC_STATS 1
//#define NON_RANDOM 1 // keep nucleating from same places

#ifndef mymax
#define mymax(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef mymin
#define mymin(a,b)            (((a) < (b)) ? (a) : (b))
#endif
  


#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers
#include <stdio.h>

// includes

// standard headers

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <vector> 

#ifndef _WIN32
	#include <unistd.h>
#else
	//#include <Windows.h>
#endif
//#include <pstream>


// namespace:

using namespace std;

// pthreads stuff:

#include "pthread.h"
#include "semaphore.h"
//#include "sched.h"

//extern char temp_BMP_filename[];

extern bool USE_THREADS;

extern int NUM_THREADS;
extern pthread_attr_t thread_attr;

extern vector<pthread_t>  threads;

extern vector<struct thread_data>  collision_thread_data_array;
extern vector<sem_t> collision_thread_go;
extern vector<sem_t> collision_data_done;

extern vector<struct thread_data>  linkforces_thread_data_array;
extern vector<sem_t> linkforces_thread_go;
extern vector<sem_t> linkforces_data_done;

extern vector<struct thread_data>  applyforces_thread_data_array;
extern vector<sem_t> applyforces_thread_go;
extern vector<sem_t> applyforces_data_done;

extern vector<struct thread_data>  compressfiles_thread_data_array;
//extern vector<sem_t> compressfiles_thread_go;
//extern vector<sem_t> compressfiles_data_done;

extern pthread_mutex_t linkstoremove_mutex;

extern pthread_mutex_t filessavelock_mutex;
extern pthread_mutex_t filesdonelock_mutex;

extern pthread_mutex_t beadmovelock_mutex;

extern vector<pthread_mutex_t> collisiondetectiongolock_mutex;
extern vector<pthread_mutex_t> collisiondetectiondonelock_mutex;

struct thread_data
{
int startnode;
int endnode;
int threadnum;
};



// math related:

#if 0      // use fmath rather than math  _WIN32
	
	// don't use this part: 
	// for some reason it's *much* slower than doubles (RT typecasting?)	

	//#include <mathf.h>

    #define sin sinf
    #define cos cosf
    #define log logf
    #define floor floorf
    #define ceil ceilf
    #define sqrt sqrtf
	#define exp expf
	#define fabs fabsf
	#define atan2 atan2f

	//#define MYDOUBLE float	// set precision

	#define SQRT_ACCURACY_LOSS (double) 0.000000000001

#else

	//#define MYDOUBLE double	// set precision

	#define SQRT_ACCURACY_LOSS (double) 0.0000000000000000001

#endif

//typedef vector<nodes*> Node1d;
//typedef vector<Node1d> Node2d;
//typedef vector<Node2d> Node3d;

//extern Node3d nodegrid;

// global variables:

extern bool GRASS_IS_GREEN;  // this is a kludge to get rid of compiler warnings

extern double TOTAL_SIMULATION_TIME;  
extern double DELTA_T;
extern double MIN_TORQUE_TO_UPDATE;
extern double MIN_DISPLACEMENT_TO_UPDATE;
//extern double MAX_DISP_PERDT;
//extern double MAX_DISP_PERDT_DIVSQRTTWO;
extern int RECORDED_TIMESTEPS;		// number of recorded timesteps(data files)

extern bool NUCLEATOR_FORCES;

extern double FORCE_SCALE_FACT;  // convert forces (nom in pN) into node displacements (nom in uM)
										// this is related to effective viscosity and effective size of node
extern double FORCEBAR_SCALE;  // scale force for bars in output

extern double XLINK_NODE_RANGE;	// Limit crosslink to within this range
extern double NODE_INCOMPRESSIBLE_RADIUS;// repulsion is zero here
//extern double NODE_REPULSIVE_MAG;   // max repulsion (at dist=0)

extern double LINK_BREAKAGE_FORCE;  // breakage force per link
extern double LINK_FORCE;
extern double P_LINK_BREAK_IF_OVER;  // probablility that force will break link if over the link breakage force
extern double P_XLINK;
extern double P_NUC;

extern int BMP_WIDTH,BMP_HEIGHT;
extern bool NO_IMAGE_TEXT;
extern int BMP_COMPRESSION;
extern string BMP_OUTPUT_FILETYPE;

extern double GAUSSFWHM;
extern int SPECKLE_FACTOR;
extern double INIT_R_GAIN;
extern double INIT_G_GAIN;
extern double INIT_B_GAIN;

extern double RADIUS;   // radius and segment are the true radius and segment of nucleator
extern double CAPSULE_HALF_LINEAR;
extern double RAD_INCOMP; // RAD_INCOMP and SEG_INCOMP are the enlarged radius and segments
//extern double SEG_INCOMP; // to prevent putting nodes within the node NODE_INCOMPRESSIBLE_RADIUS of the nucleator
//extern double NODEMASS;
//extern double INERTIAL_DAMPING_HALFTIME;
//extern double DAMPING_FACTOR;
extern int CROSSLINKDELAY;
extern int NODES_TO_UPDATE;
extern double DISTANCE_TO_UPDATE;

extern double NODE_REPULSIVE_MAG;
extern double NODE_REPULSIVE_RANGE;
extern int InterRecordIterations;
extern unsigned int MAX_LINKS_PER_NODE;

extern int RADIAL_SEGMENTS;
extern int XLINK_NEAREST;

extern double VIEW_HEIGHT;

extern double LINK_TAUT_FORCE;
extern double LINK_TAUT_RATIO;

extern int ASYMMETRIC_NUCLEATION;

//extern int REPORT_AVERAGE_ITTERATIONS;

const int REPORT_NUM_VARIABLES = 8;

extern bool ROTATION;
extern double MofI;

extern bool FORCES_ON_SIDE;

extern double colormap[];


// compile-time options:

#define GRIDBOUNDS (double)40	  // size of grid in um
#define GRIDRES (double)0.2	  // low res grid range

const int MAXNODES = 50000;			// max nodes
//const double GLOBAL_DAMPING = 

#define HIST_MIN 1.0
#define HIST_MAX 3.0
#define HIST_BINS 20

extern int TOTAL_ITERATIONS;  // these variables are global and are calculated from others
extern int NODE_REPULSIVE_GRIDSEARCH;
extern int NODE_REPULSIVE_RANGE_GRIDSEARCH;
extern int NODE_XLINK_GRIDSEARCH;
//extern int GRIDSIZE;
const int GRIDSIZE =  (int) (GRIDBOUNDS/GRIDRES);

const double PI = (double) 3.141592653589793238462643383279502884197; // Pi
const double LN_TWO = (double) 0.69314718055995; // ln(2)

//inline double sqrt(const double &d);
inline double calcdist(const double & xdist, const double & ydist, const double & zdist);
inline double calcdist(const double & xdist, const double & ydist); 

inline void endian_swap(unsigned short& x);
inline void endian_swap(unsigned int& x);


// own headers
#include "comet.h"
#include "nucleator.h"

extern nucleator::shape NUCSHAPE;  //default to sphere

#include "nodes.h"
#include "links.h"
#include "actin.h"
#include "vect.h"



inline double calcdist(const vect & v1, const vect & v2);

/*inline float sqrt(float x) {
  float y;
  _asm {
    mov ecx, x
    and ecx, 0x7f800000
    fld x
    mov eax, ecx
    add ecx, 0x3f800000
    or  ecx, 0x00800000
    shr ecx, 1
    mov y, eax
    and ecx, 0x3f800000
    fld y
    fadd
    fstp y
    and [y], 0x007fffff
    or  [y], ecx
  }
  return y;
}  */


/*inline double sqrt(double d)
{
	double x,y, deltax;
	x = 1/d;
	do
	{
		y = x*x;
		deltax = x*(0.875 - ( (1.25*d) - (0.375*d*2) *y) *y);
		x+= deltax;
	} while (fabs(deltax) > (d/100));

	return (1/x);
}*/
/*
static inline double sqrt(double x)
{
	// "delta" is the acceptable error bound
	double delta = x/100;
	double y = (1+x)*.5,oldy;
	do {
		oldy = y;
		y = (x/y+y)*.5;
	}
	while(fabs(y-oldy) > delta);
	return y;
}
*/

//inline double sqrt(const double &d)
//{
//	return sqrt(d);
//}

//inline double calcdist(double xdist, double ydist, double zdist);
inline double calcdist(const double & xdist, const double & ydist, const double & zdist)
{
	double sqr = (xdist*xdist + ydist*ydist + zdist*zdist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 3D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return sqrt(sqr);
}

//inline double calcdist(double xdist, double ydist);
inline double calcdist(const double & xdist, const double & ydist) 
{
	double sqr = (xdist*xdist + ydist*ydist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 2D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return sqrt(sqr);
}

inline double calcdist(const vect & v1, const vect & v2) 
{
	return calcdist(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
}

inline void endian_swap(unsigned short& x)
{
    x = (unsigned short)( (x>>8) | (x<<8));
}

inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
	((x<<8) & 0x00FF0000) |
	((x>>8) & 0x0000FF00) |
	(x<<24);
}



/*
// __int64 for MSVC, "long long" for gcc
inline void endian_swap(unsigned long long& x)
{
    x = (x>>56) | 
	((x<<40) & 0x00FF000000000000) |
	((x<<24) & 0x0000FF0000000000) |
	((x<<8)  & 0x000000FF00000000) |
	((x>>8)  & 0x00000000FF000000) |
	((x>>24) & 0x0000000000FF0000) |
	((x>>40) & 0x000000000000FF00) |
        (x<<56);
}*/






// extern actin theactin;

//extern nodes* nodegrid[GRIDSIZE+1][GRIDSIZE+1][GRIDSIZE+1];

typedef vector<nodes*> Nodes1d;
typedef vector<Nodes1d> Nodes2d;
typedef vector<Nodes2d> Nodes3d;

typedef vector<signed char> Bool1d;
typedef vector<Bool1d> Bool2d;

typedef vector<double> Dbl1d;
typedef vector<Dbl1d> Dbl2d;

extern Nodes3d nodegrid;

//extern inline double fastsqrt(float n);
//extern float sse_sqrt(float n);
//extern void build_sqrt_table();


#endif

