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

// compile-time options:

// whether link or proximity viscosity
//#define LINK_VISCOSITY 1
//#define PROXIMITY_VISCOSITY 1
// note: also need to set VISCOSIY variable to true in cometparams.ini

// this doesn't seem to affect speed at all
#define NODE_GRID_USE_ARRAYS 1

//#define NOKBHIT 1
#define VTK_USE_ANSI_STDLIB

// NODEGRIDTYPE should be 'list' or 'vector'
// vector seems quite a bit faster
//#define NODEGRIDTYPELIST 1

//#define BMPS_USING_LINKS 1

#define NO_CALC_STATS 1 // calculating stats slows things considerably (prolly because of trig functions,
                        // but is required if we want force/energy graphs)

// #define NON_RANDOM 1   // keep nucleating from same places

#if defined (__SSE__) 
//#define USE_SSE_APPROX_SQRT 1
#endif

#define BMP_USE_FOCAL_DEPTH 1

//#define TIMERCPUTIME  // timers report cpu time instead of wall time

// mersennetwister.h is causing bus error for some reason...
//#define USE_MERSENNE 1


// test optimizations
#ifdef NODEGRIDTYPELIST
    #define NODEGRIDTYPE list
#else
    #define NODEGRIDTYPE vector
#endif

#ifdef NODE_GRID_USE_ARRAYS
	#define NODEGRID(i,j,k)  (*(nodegrid + ((((GRIDSIZE*GRIDSIZE*(i)) + (GRIDSIZE*(j)) + (k))))))
#else
	#define NODEGRID(i,j,k)	 nodegrid[(i)][(j)][(k)]
#endif

const unsigned int MAX_EXPECTED_LINKS = 20;   // reserves this no of links per node
                                              // OK if goes over, just slows things

#define SYM_BREAK_FILE "sym_break_axis.txt"
#define SEG_SCALE_FILE "segscalefactors.txt"
#define COMET_PARAMS_FILE "cometparams.ini"
#define NODESUPDATEFILE "nodesupdate.txt"
#define TESTNODESFILE "testnodes.txt"

#define USEBZIP2 1

#ifdef USEBZIP2
 #define COMPRESSCOMMAND "bzip2 -9 -f"
 #define DECOMPRESSCOMMAND "bzip2 -d"
 #define COMPRESSEDEXTENSION ".bz2"
#else
 #define COMPRESSCOMMAND "gzip -9 -f"                                     
 #define DECOMPRESSCOMMAND "gunzip"
 #define COMPRESSEDEXTENSION ".gz"
#endif


#ifdef _WIN32

    #define LINK_VTK 1    // link with vtk if on windows

    #define _CRT_SECURE_NO_DEPRECATE
    #define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

    // whether to use windows commands (/ or \ etc)
    // this is not *necessarily* a _WIN32 switch, since we may be using cygwin
	
    #define USEWINDOWSCOMMANDS

    #include <windows.h>
    #include <direct.h>      // windows name for dirent.h
	#include <process.h>
	#define getpid _getpid
	//#include <windows.h>

	//#pragma warning(disable: 4511) // unable to generate copy constructor
	//#pragma warning(disable: 4512)

	#pragma warning(disable: 4127)  // constant conditional expression

	#pragma warning(disable: 4996) // turn off depreciated warnings
	#pragma warning(disable: 4530) // turn off exception handling warning

	//#pragma auto_inline( on )

#else

    #include <dirent.h>
	#include <unistd.h>
    #include <sys/mman.h>

#endif

 
	extern char VRMLDIR[];
	extern char DATADIR[];
	extern char REPORTDIR[];
	extern char BITMAPDIR[];
	extern char TEMPDIR[];
	//extern char TEMPDIR2[];
	extern char VTKDIR[];
	extern char STATSDIR[];

    extern char IMAGEMAGICKCONVERT[];
    extern char IMAGEMAGICKMOGRIFY[];



#ifndef mymax
#define mymax(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef mymin
#define mymin(a,b)  (((a) < (b)) ? (a) : (b))
#endif

// includes

// standard headers


#include <stdio.h>
#include <assert.h>
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
#include <list>
#include <float.h>




// namespace:

using namespace std;

enum projection 
{
	xaxis = 0,
	yaxis = 1,
	zaxis = 2
};

// pthreads stuff:

#include "pthread.h"
#include "semaphore.h"
#include "threadedtaskqueue.h"

#ifdef USE_MERSENNE

#include "MersenneTwister.h"

extern vector <MTRand> mers_rand;

#endif

//class actin;


extern bool ALLOW_HARBINGERS_TO_MOVE;
extern bool CAGE_ON_SIDE;
extern bool SYM_BREAK_TO_RIGHT;

extern bool COMPRESSDATAFILES; 

extern double GRIDBOUNDS;
extern double GRIDRES;

extern bool VISCOSITY;
extern double VISCOSITY_EDGE_FACTOR;
extern double VISC_DIST;
extern double MAX_VISC_WEIGHTING;

extern bool USE_BROWNIAN_FORCES;
extern double BROWNIANFORCESCALE;

extern bool USE_THREADS;
extern int NUM_THREADS;
extern int NUM_THREAD_DATA_CHUNKS;

extern bool COL_NODE_BY_STRAIN;
extern bool COL_LINK_BY_DIRN;
extern bool COL_INDIVIDUAL_NODES;

//-- Threading
extern TaskQueue thread_queue;
extern bool USETHREAD_LINKFORCES;
extern bool USETHREAD_APPLYFORCES;
extern bool USETHREAD_COLLISION;

extern bool SEGMENT_BINS;
extern bool DRAW_CAGE;
//-- 

extern bool NOBGIMAGEMAGICK;

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

// global variables:

extern int MAXNODES;

extern bool REWRITESYMBREAK;
extern bool POST_PROCESS;

extern int POST_PROC_ORDER;


extern double TOTAL_SIMULATION_TIME;  
extern double DELTA_T;

extern bool STICK_TO_NUCLEATOR;
extern bool RESTICK_TO_NUCLEATOR;

extern bool NUCLEATOR_FORCES;


extern double NUC_LINK_FORCE;
extern double NUC_LINK_BREAKAGE_FORCE;

extern double FORCE_SCALE_FACT;  // convert forces (nom in pN) into node displacements (nom in uM)
										// this is related to effective viscosity and effective size of node
extern double FORCE_BAR_SCALE;  // scale force for bars in output

extern double XLINK_NODE_RANGE;	// Limit crosslink to within this range

extern double LINK_BREAKAGE_FORCE;  // breakage force per link
extern double LINK_FORCE;
extern bool USE_BREAKAGE_VISCOSITY;
extern double BREAKAGE_VISCOSITY_THRESHOLD;
extern double LINK_BREAKAGE_STRAIN;
extern double P_XLINK;
extern double P_NUC;
extern double MAX_POLYMERISATION_PRESSURE;

extern bool POST_VTK;
extern bool POST_STATS;

extern int BMP_WIDTH,BMP_HEIGHT;
extern int VTK_WIDTH;
extern int VTK_HEIGHT;
extern int VTK_AA_FACTOR;
extern int BMP_AA_FACTOR;

extern bool VTK_NUC_LINKS_ONLY;
extern double VTK_LINK_COLOUR_GAMMA;
extern bool VTK_MOVE_WITH_BEAD;
extern double VTK_VIEWANGLE;

extern bool NO_IMAGE_TEXT;
extern int BMP_COMPRESSION;
extern string BMP_OUTPUT_FILETYPE;

extern bool QUIET;

extern double LINK_POWER_SCALE;

extern bool X_BMP;										 
extern bool Y_BMP;
extern bool Z_BMP;

extern bool PLOTFORCES;

extern double FOCALDEPTH;

extern double BMP_INTENSITY_OFFSET;

extern double gridscanjitter;

extern double VISCOSITY_FACTOR;
extern double NON_VISC_WEIGHTING;

extern double GAUSSFWHM;
extern double SPECKLE_FACTOR;
extern bool   SPECKLE;

extern bool SPECKLEGRID;
extern double SPECKLEGRIDPERIOD;
extern double SPECKLEGRIDTIMEWIDTH;
extern double SPECKLEGRIDSTRIPEWIDTH;

extern bool POLY_FEEDBACK;
extern double POLY_FEEDBACK_DIST;
extern double POLY_FEEDBACK_MIN_PROB;
extern double POLY_FEEDBACK_FACTOR;

extern double BMP_INTENSITY_SCALE;
extern double INIT_R_GAIN;
extern double INIT_G_GAIN;
extern double INIT_B_GAIN;

extern double NODE_FORCE_TO_DIST;
extern double NODE_DIST_TO_FORCE;

extern double COVERSLIPGAP;

extern double IMPOSED_NUC_ROT_SPEED;
extern bool   IMPOSED_NUC_ROT;


extern bool   TEST_SQUASH;
extern double TEST_FORCE_INITIAL_MAG;
extern double TEST_FORCE_INCREMENT;
extern double TEST_DIST_EQUIL;

extern bool WRITE_BMPS_PRE_SYMBREAK;

extern double NUCLEATOR_INERTIA;
extern bool VARY_INERT_W_RAD;

extern double VTK_MIN_PLOT_LINK_FORCE_PCT;

extern double RADIUS;   // radius and segment are the true radius and segment of nucleator
extern double CAPSULE_HALF_LINEAR;
extern double RADIUS; // RADIUS and SEG_INCOMP are the enlarged radius and segments
extern int CROSSLINKDELAY;
extern int NODES_TO_UPDATE;
extern double DISTANCE_TO_UPDATE;

extern double NODE_REPULSIVE_POWER;
extern double NODE_REPULSIVE_MAG;
extern double NODE_REPULSIVE_RANGE;
extern double NODE_REPULSIVE_BUCKLE_RANGE;
extern double NODE_REPULSIVE_BUCKLE_MAG;
extern double NODE_REPULSIVE_BUCKLE_TO;
extern int InterRecordIterations;
extern unsigned int MAX_LINKS_PER_NEW_NODE;
extern unsigned int MAX_LINK_ATTEMPTS;

extern int RADIAL_SEGMENTS;
extern int XLINK_NEAREST;

extern double VIEW_HEIGHT;

extern int ASYMMETRIC_NUCLEATION;

const int REPORT_NUM_VARIABLES = 8;

extern bool ROTATION;
extern double MOFI;

extern bool FORCES_ON_SIDE;

extern int TOTAL_ITERATIONS;  // these variables are global and are calculated from others
extern int NODE_REPULSIVE_GRIDSEARCH;
extern int NODE_REPULSIVE_RANGE_GRIDSEARCH;
extern int NODE_XLINK_GRIDSEARCH;
extern short int GRIDSIZE;

extern bool BMP_TRACKS;

extern double TRACK_MIN_RANGE;
extern double TRACK_MAX_RANGE;

extern int TRACKFRAMESTEP;

extern bool BMP_FIX_BEAD_MOVEMENT;
extern bool BMP_FIX_BEAD_ROTATION;

const double RECIP_RAND_MAX =  (1/(double)RAND_MAX);
const double PI = (double) 3.141592653589793238462643383279502884197; // Pi
const double LN_TWO = (double) 0.69314718055995; // ln(2)

inline double calcdist(const double & xdist, const double & ydist, const double & zdist);
inline double calcdist(const double & xdist, const double & ydist); 

inline double recipcalcdist(const double & xdist, const double & ydist, const double & zdist);
inline double recipcalcdist(const double & xdist, const double & ydist);

inline void endian_swap(unsigned short& x);
inline void endian_swap(unsigned int& x);


// own headers
#include "comet.h"
#include "nucleator.h"



extern nucleator::shape NUCSHAPE;  //default to sphere

#include "nodes.h"

typedef vector<nodes*> Nodes1d;
typedef vector<Nodes1d> Nodes2d;
typedef vector<Nodes2d> Nodes3d;

typedef NODEGRIDTYPE<nodes*> NG1d;
typedef vector<NG1d> NG2d;
typedef vector<NG2d> NG3d;
typedef vector<NG3d> NG4d;

typedef vector<signed char> Bool1d;
typedef vector<Bool1d> Bool2d;

typedef vector<double> Dbl1d;
typedef vector<Dbl1d> Dbl2d;

struct thread_data
{
    vector <nodes>::iterator startnode;
    vector <nodes>::iterator endnode;
    int threadnum;
};

extern vector<struct thread_data>  collision_thread_data_array;
extern vector<struct thread_data>  linkforces_thread_data_array;
extern vector<struct thread_data>  applyforces_thread_data_array;


#ifdef NODE_GRID_USE_ARRAYS
	extern NG1d* __restrict nodegrid;
    // extern nodes** __restrict nodegrid;
#else
	//extern Nodes3d nodegrid;
    extern NG4d nodegrid;
#endif

extern bool VARY_P_XLINK;

#include "links.h"
#include "actin.h"
#include "vect.h"

extern actin *ptheactin;

//#include "consts.h"

//extern consts CONST;


#if defined(USE_SSE_APPROX_SQRT)

 inline float SSErsqrt(float x) {
 
         __asm {
            movss xmm0, x
            rsqrtss xmm0, xmm0
            movss x, xmm0
         }
         return x;
 }
  
 inline float SSEsqrt(float x) {

         __asm {
            movss xmm0, x
            sqrtss xmm0, xmm0
            movss x, xmm0
         }
         return x;
 }

 inline double SSErsqrt(double & x)
 {
     return (double) SSErsqrt((float)x);
 }

 inline double SSEsqrt(double & x)
 {
     return (double) SSEsqrt((float)x);
 }

#else

 inline float SSErsqrt(float x) {

         return 1 / sqrt(x);  
 }
  
 inline float SSEsqrt(float x) {
 
         return sqrt(x);  

 }

 inline double SSErsqrt(double & x)
 {
     return 1 / sqrt(x);
 }

 inline double SSEsqrt(double & x)
 {
     return sqrt(x);
 }


#endif


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
		return SSEsqrt(sqr);
}

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
		return SSEsqrt(sqr);
}

inline double calcdist(const vect & v1, const vect & v2) 
{
	return calcdist(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
}

inline double recipcalcdist(const double & xdist, const double & ydist, const double & zdist)
{
	double sqr = (xdist*xdist + ydist*ydist + zdist*zdist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 3D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return SSErsqrt(sqr);
}

inline double recipcalcdist(const double & xdist, const double & ydist) 
{
	double sqr = (xdist*xdist + ydist*ydist);
	if (sqr < SQRT_ACCURACY_LOSS)
	{
		//cout << "Accuracy loss: 2D dist close to zero. Increase math accuracy or reduce DELTA_T" << endl;
		//cout.flush();
		return SQRT_ACCURACY_LOSS;
	}
	else
		return SSErsqrt(sqr);
}

inline double recipcalcdist(const vect & v1, const vect & v2) 
{
	return recipcalcdist(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
}

inline void endian_swap(unsigned short& x)
{
    x = (unsigned short)((x>>8) | (x<<8));
}

inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
	    ((x<<8) & 0x00FF0000) |
	    ((x>>8) & 0x0000FF00) |
	    (x<<24);
}

#ifdef USE_MERSENNE

inline bool prob_to_bool(const double & prob)
{
    return (mers_rand[0].rand() < prob);
}

inline bool prob_to_bool(const double & prob, const unsigned int & thread)
{
    return (mers_rand[thread].rand() < prob);
}

inline double rand_0to1()
{
    return (mers_rand[0].rand());
}

#else

inline bool prob_to_bool(const double & prob)
{
    return (rand() < prob * RAND_MAX);
}

inline bool prob_to_bool(const double & prob, const unsigned int & )
{
    return (rand() < prob * RAND_MAX);
}

inline double rand_0to1()
{
    return ((double)rand() * (double) RECIP_RAND_MAX);
}

#endif

#endif

