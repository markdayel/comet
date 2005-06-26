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

// defines

// #define FORCES_BOTH_WAYS 1

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

#ifdef _WIN32

	#pragma inline_depth( 255 )
	#pragma inline_recursion( on )
	#pragma auto_inline( on )

#endif

// math related:

#if 1      // use fmath rather than math  _WIN32
	
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

	#define MYDOUBLE float	// set precision

	#define SQRT_ACCURACY_LOSS (MYDOUBLE) 0.000000000001

#else

	#define MYDOUBLE double	// set precision
	#define SQRT_ACCURACY_LOSS (MYDOUBLE) 0.0000000000000000001

#endif

//typedef vector<nodes*> Node1d;
//typedef vector<Node1d> Node2d;
//typedef vector<Node2d> Node3d;

//extern Node3d nodegrid;

// global variables:

extern MYDOUBLE TOTAL_SIMULATION_TIME;  
extern MYDOUBLE DELTA_T;
extern MYDOUBLE MAX_DISP_PERDT;
extern MYDOUBLE MAX_DISP_PERDT_DIVSQRTTWO;
extern int RECORDED_TIMESTEPS;		// number of recorded timesteps(data files)

extern MYDOUBLE FORCE_SCALE_FACT;  // convert forces (nom in pN) into node displacements (nom in uM)
										// this is related to effective viscosity and effective size of node

extern MYDOUBLE XLINK_NODE_RANGE;	// Limit crosslink to within this range
extern MYDOUBLE NODE_INCOMPRESSIBLE_RADIUS;// repulsion is zero here
//extern MYDOUBLE NODE_REPULSIVE_MAG;   // max repulsion (at dist=0)

extern MYDOUBLE LINK_BREAKAGE_FORCE;  // breakage force per link
extern MYDOUBLE LINK_FORCE;
extern MYDOUBLE P_LINK_BREAK_IF_OVER;  // probablility that force will break link if over the link breakage force
extern MYDOUBLE P_XLINK;
extern MYDOUBLE P_NUC;

extern MYDOUBLE RADIUS;   // radius and segment are the true radius and segment of nucleator
extern MYDOUBLE SEGMENT;
extern MYDOUBLE RAD_INCOMP; // RAD_INCOMP and SEG_INCOMP are the enlarged radius and segments
//extern MYDOUBLE SEG_INCOMP; // to prevent putting nodes within the node NODE_INCOMPRESSIBLE_RADIUS of the nucleator
//extern MYDOUBLE NODEMASS;
//extern MYDOUBLE INERTIAL_DAMPING_HALFTIME;
//extern MYDOUBLE DAMPING_FACTOR;
extern int CROSSLINKDELAY;
extern int NODES_TO_UPDATE;

extern MYDOUBLE NODE_REPULSIVE_MAG;
extern MYDOUBLE NODE_REPULSIVE_RANGE;
extern int InterRecordIterations;
extern unsigned int MAX_LINKS_PER_NODE;

extern int RADIAL_SEGMENTS;
extern int XLINK_NEAREST;

extern MYDOUBLE VIEW_HEIGHT;

extern MYDOUBLE LINK_TAUGHT_FORCE;
extern MYDOUBLE LINK_TAUGHT_RATIO;

extern int ASYMMETRIC_NUCLEATION;

extern int REPORT_AVERAGE_ITTERATIONS;

const int REPORT_NUM_VARIABLES = 8;


// compile-time options:

#define GRIDBOUNDS (MYDOUBLE)40	  // size of grid in um
#define GRIDRES (MYDOUBLE)0.2	  // low res grid range

const int MAXNODES = 50000;			// max nodes
//const MYDOUBLE GLOBAL_DAMPING = 

#define HIST_MIN 1.0
#define HIST_MAX 3.0
#define HIST_BINS 20

extern int TOTAL_ITERATIONS;  // these variables are global and are calculated from others
extern int NODE_REPULSIVE_GRIDSEARCH;
extern int NODE_REPULSIVE_RANGE_GRIDSEARCH;
extern int NODE_XLINK_GRIDSEARCH;
//extern int GRIDSIZE;
const int GRIDSIZE =  (int) (GRIDBOUNDS/GRIDRES);

const MYDOUBLE PI = (MYDOUBLE) 3.141592653589793238462643383279502884197; // Pi
const MYDOUBLE LN_TWO = (MYDOUBLE) 0.69314718055995; // ln(2)

// own headers
#include "comet.h"
#include "nucleator.h"
#include "nodes.h"
#include "links.h"
#include "actin.h"

// extern actin theactin;

//extern nodes* nodegrid[GRIDSIZE+1][GRIDSIZE+1][GRIDSIZE+1];

typedef vector<nodes*> Nodes1d;
typedef vector<Nodes1d> Nodes2d;
typedef vector<Nodes2d> Nodes3d;

typedef vector<signed char> Bool1d;
typedef vector<Bool1d> Bool2d;

typedef vector<MYDOUBLE> Dbl1d;
typedef vector<Dbl1d> Dbl2d;

extern Nodes3d nodegrid;

//extern inline MYDOUBLE fastsqrt(float n);
//extern float sse_sqrt(float n);
//extern void build_sqrt_table();

/*inline float mysqrt(float x) {
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


/*inline MYDOUBLE mysqrt(MYDOUBLE d)
{
	MYDOUBLE x,y, deltax;
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
static inline MYDOUBLE mysqrt(MYDOUBLE x)
{
	// "delta" is the acceptable error bound
	MYDOUBLE delta = x/100;
	MYDOUBLE y = (1+x)*.5,oldy;
	do {
		oldy = y;
		y = (x/y+y)*.5;
	}
	while(fabs(y-oldy) > delta);
	return y;
}
*/
inline MYDOUBLE mysqrt(MYDOUBLE d)
{
	return sqrt(d);
}

inline MYDOUBLE calcdist(MYDOUBLE xdist, MYDOUBLE ydist, MYDOUBLE zdist)
{
	MYDOUBLE sqr = (xdist*xdist + ydist*ydist + zdist*zdist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 3D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return mysqrt(sqr);
}

inline MYDOUBLE calcdist(MYDOUBLE xdist, MYDOUBLE ydist)
{
	MYDOUBLE sqr = (xdist*xdist + ydist*ydist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 2D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return mysqrt(sqr);
}

inline void endian_swap(unsigned short& x)
{
    x = (x>>8) | 
	(x<<8);
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

#endif

