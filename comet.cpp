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

#include "stdafx.h"
#include "comet.h"
#include "threadedtaskqueue.h"

// vtk
#ifdef LINK_VTK
  #include "comet_vtk.h"
#else
  #include "comet_novtk.h"
#endif

// #include "nucleator.h"
// #include "string.h"

double TOTAL_SIMULATION_TIME = 20000;  
double DELTA_T = 0.1;	
double MIN_TORQUE_TO_UPDATE = 0.2;
double MIN_DISPLACEMENT_TO_UPDATE = 0.001;
//double MAX_DISP_PERDT = 0.01;
//double MAX_DISP_PERDT_DIVSQRTTWO = 0.00707;
double GAUSSFWHM =  0.266;

bool NUCLEATOR_FORCES = true;

int BMP_WIDTH  = 800;
int BMP_HEIGHT = 600;

double INIT_R_GAIN = 20;
double INIT_G_GAIN = 20;
double INIT_B_GAIN = 20;

bool ALLOW_HARBINGERS_TO_MOVE = false;

bool   SPECKLE = false;
double SPECKLE_FACTOR = 1;

bool ROTATION = true;
bool GRASS_IS_GREEN = true;
bool DRAW_CAGE = false;

bool SEGMENT_BINS = true;

bool CAGE_ON_SIDE = true;

double MofI =  0.1;

bool VISCOSITY = false;

//bool FORCES_ON_SIDE = true;

bool NO_IMAGE_TEXT = false;
int BMP_COMPRESSION = 75;
string BMP_OUTPUT_FILETYPE = "png";

//int REPORT_AVERAGE_ITTERATIONS = 50;

int RECORDED_TIMESTEPS=200;			// number of recorded timesteps(data files)

int RESTORE_FROM_ITERATION = 0; // =0 don't load a checkpoint 
int RECORDING_INTERVAL = 0;
int NUMBER_RECORDINGS = 0;

double FORCE_SCALE_FACT = 0.001;	// convert forces (nom in pN) into node displacements (nom in uM)
					        // this is related to effective viscosity and effective size of node
double FORCE_BAR_SCALE   = 1;	// scale for relating force to image bar

double XLINK_NODE_RANGE =  1.0;		// Limit crosslink to within this range
//double NODE_INCOMPRESSIBLE_RADIUS = 0.2;	// repulsion is zero here
// double NODE_REPULSIVE_MAG = 1000;   // max repulsion (at dist=0)

double LINK_BREAKAGE_FORCE =  2;	 // breakage force per link
bool USE_BREAKAGE_STRAIN = false;
double LINK_BREAKAGE_STRAIN = 1.15;
double P_LINK_BREAK_IF_OVER =  0.25;  // probablility that force will break link if over the link breakage force
unsigned int MAX_LINKS_PER_NODE = 100;

double NUCLEATOR_INERTIA = 10;

double NUC_LINK_FORCE = 0.25;
double NUC_LINK_BREAKAGE_FORCE = 2;

double LINK_TAUT_FORCE =  5;
double LINK_TAUT_RATIO =  1.1;

double IMPOSED_NUC_ROT_SPEED = 1;
bool   IMPOSED_NUC_ROT = false;

bool WRITE_BMPS_PRE_SYMBREAK = false;

bool QUIET = false;

double VISCOSITY_FACTOR = 0.1;
//double VISCOSITY_UNWEIGHTING_FACTOR = 100;

//char temp_BMP_filename[255];

bool STICK_TO_NUCLEATOR = false;
bool RESTICK_TO_NUCLEATOR = true;

double LINK_FORCE = 0.1;
double P_XLINK =  0.5;
double P_NUC =  0.08;
double RADIUS =  1.0;
double CAPSULE_HALF_LINEAR =  6.0;

nucleator::shape NUCSHAPE = nucleator::sphere;  //default to sphere

//double SEG_INCOMP = 2*CAPSULE_HALF_LINEAR + NODE_INCOMPRESSIBLE_RADIUS/2;
//double RADIUS = RADIUS;// + NODE_INCOMPRESSIBLE_RADIUS/2;

//double NODEMASS = 1.0;
//double INERTIAL_DAMPING_HALFTIME = 50;

int TOTAL_ITERATIONS ;
int NODE_REPULSIVE_GRIDSEARCH ;
int NODE_XLINK_GRIDSEARCH ;
int NODE_REPULSIVE_RANGE_GRIDSEARCH;

int RADIAL_SEGMENTS = 12;
int NODES_TO_UPDATE = 5000;  //only update the NODES_TO_UPDATE newest nodes
double DISTANCE_TO_UPDATE = 0;

//double DAMPING_FACTOR = 10;
int CROSSLINKDELAY = 20;  // number of interations before crosslinking 
						 //  (to allow position to be equilibrated to something
						 //       reasonable before locking node in place)

double NODE_REPULSIVE_MAG =  0.00000001;
double NODE_REPULSIVE_RANGE = 1.0;

double NODE_FORCE_TO_DIST;
double NODE_DIST_TO_FORCE;

int ASYMMETRIC_NUCLEATION = 0;

int XLINK_NEAREST = 1;

double VIEW_HEIGHT = 12;

bool USE_THREADS;
int NUM_THREADS;
int NUM_THREAD_DATA_CHUNKS;
//-- ThreadTaskTeam
TaskQueue thread_queue;
//pthread_mutex_t removelinks_mutex;
//pthread_mutex_t nodedone_mutex;
bool USETHREAD_COLLISION   = true;
bool USETHREAD_APPLYFORCES = true;
bool USETHREAD_LINKFORCES  = true;
// --

#ifdef _NUMA
    //nsgid_t numa_group;
    radset_t radset;
    //cpuset_t cpuset;
#endif

pthread_attr_t thread_attr;
vector<pthread_t> threads;

vector<struct thread_data>  collision_thread_data_array;
vector<struct thread_data>  linkforces_thread_data_array;
vector<struct thread_data>  applyforces_thread_data_array;


// these variables need to be static/global for sharing across threads:

Nodes3d nodegrid;
vector <nodes>	actin::node;
vector <bool>   actin::donenode;	
Nodes2d actin::nodes_by_thread;
Nodes2d actin::recti_near_nodes;
Nodes2d actin::nodes_on_same_gridpoint;
//Nodes1d actin::nodes_within_nucleator;

//vector <int> actin::recti_near_nodes_size;
//vector <int> actin::nodes_on_same_gridpoint_size;

int actin::iteration_num;

bool actin::isinthread;
Nodes2d actin::linkremovefrom;
Nodes2d actin::linkremoveto;

bool REWRITESYMBREAK = false;
bool POST_PROCESS = false;

int InterRecordIterations = 0;

int load_data(actin &theactin, int iteration);
int save_data(actin &theactin, int iteration);
string get_datafilename(const int iteration);
void get_postprocess_iterations(const char *iterdesc, vector<int> &postprocess_iterations);
void postprocess(nucleator& nuc_object, actin &theactin, vector<int> &postprocess_iterations);
void rewrite_symbreak_bitmaps(nucleator& nuc_object, actin &theactin);

//#define NOKBHIT 1

#ifndef NOKBHIT
#ifndef _WIN32
	#include "kbhit.h"
	#define kbhit keyb.kbhit
    #define getch keyb.getch
#else
	#include <conio.h>
#endif
#endif

// main 

int main(int argc, char* argv[])
{
	cout << endl; 

    if(argc < 2 || argc > 4) 
	{
	    cerr << "Usage:" << endl << endl << argv[0] << " numThreads [R]" << endl << endl;
	    cerr << "where numThreads is the number of threads to use per calculation stage" << endl;
	    cerr << "Set numThreads to 1 to run in single threaded mode" << endl;
	    cerr << "and 'R' sets use of time-based random number seed" << endl << endl;
	    exit(EXIT_FAILURE);
	}

    if (argc == 2) 
	{
		if (strcmp(argv[1], "sym") == 0 ) 
			REWRITESYMBREAK = true;
	}
	else if (argc > 2)
	{
		if (strcmp(argv[1], "post") == 0 )
			POST_PROCESS = true;

	} else if (argc > 3) 
	{
		if (strcmp(argv[3], "q") == 0 || strcmp(argv[3], "Q") == 0)
			QUIET = true;
	}
 
  	NUM_THREADS = atoi(argv[1]);

	if (NUM_THREADS < 2)
	{
		cout << "Running in Single Threaded mode" << endl;
		NUM_THREADS = 1;
		USE_THREADS = false;
	}
	else
	{
		cout << "Running in multithreaded mode with " << NUM_THREADS << " threads (including parent)" << endl;
		USE_THREADS = true;
	}

    // create threads:
	// -- Threading TaskTeam, create and intialise the team
	if (USE_THREADS  && !POST_PROCESS && !REWRITESYMBREAK)
	{
#ifdef _NUMA
        radsetcreate(&radset);
        pid_t pid = getpid();
        rad_attach_pid(pid, radset, RAD_INSIST || RAD_SMALLMEM || RAD_MIGRATE || RAD_WAIT);
        //pid_t pid = getpid();
        //numa_group = nsg_init(pid, NSG_GETBYPID);

        //cpusetcreate(&cpuset);
        //cpuid_t cpu = cpu_get_current();
        //cpuaddset(cpuset, cpu);
#endif
	    thread_queue.create_threads(NUM_THREADS);
	}


	NUM_THREAD_DATA_CHUNKS = NUM_THREADS;  // was NUM_THREADS*4

	// data for threads (managed outside queue
	collision_thread_data_array.resize(NUM_THREAD_DATA_CHUNKS);
	linkforces_thread_data_array.resize(NUM_THREAD_DATA_CHUNKS);
	applyforces_thread_data_array.resize(NUM_THREAD_DATA_CHUNKS);


	// make directories

	char command1[255];

	cout << "Working Directory: ";
	cout.flush();

#ifndef _WIN32
    system("pwd");
#else
	system("chdir");
#endif

	
#ifndef NOKBHIT
#ifndef _WIN32
    keyboard keyb;
#endif
#endif

#ifndef _WIN32
	
	sprintf(command1, "mkdir %s 2>/dev/null", VRMLDIR  );
	system(command1);	
	sprintf(command1, "mkdir %s 2>/dev/null", DATADIR  );
	system(command1);
	sprintf(command1, "mkdir %s 2>/dev/null", REPORTDIR  );
	system(command1);
	sprintf(command1, "mkdir %s 2>/dev/null", BITMAPDIR  );
	system(command1);
	sprintf(command1, "mkdir %s 2>/dev/null", TEMPDIR  );
	system(command1);
#else
	
	sprintf(command1, "mkdir %s", VRMLDIR  );
	system(command1);
	sprintf(command1, "mkdir %s", DATADIR  );
	system(command1);
	sprintf(command1, "mkdir %s", REPORTDIR  );
	system(command1);
	sprintf(command1, "mkdir %s", BITMAPDIR  );
	system(command1);
	sprintf(command1, "mkdir %s", TEMPDIR  );
	system(command1);
#endif

	
if (!REWRITESYMBREAK)
{	// don't drop priority if re-writing bitmaps
	// because calling thread already has low priority
	// and this would drop it further
	// so process would halt on single cpu machine

#ifndef _WIN32
    char hostname[255];
    gethostname( hostname, 255);
    if (strcmp( hostname, "adenine.cgl.ucsf.edu") == 0 ||
        strcmp( hostname, "cytosine.cgl.ucsf.edu") == 0||
        strcmp( hostname, "thymine.cgl.ucsf.edu") == 0||
        strcmp( hostname, "uracil.cgl.ucsf.edu")== 0 )
    {
        nice(19);
    }
    else
    {
	    nice(15);
    }
#else
	//SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_LOWEST);
#endif
}


	// main parameters:

	// read the parameters file:

	cout << endl << "Parsing parameters file..." << endl << endl;

	ifstream param(COMET_PARAMS_FILE); 
	if(!param) 
	{
		cerr << "Cannot open " << COMET_PARAMS_FILE << endl << endl;
		exit(EXIT_FAILURE);
	}

	//double MAX_DISP = 1;

	string buffer;
	string unrecognisedlines;

	while (getline(param, buffer)) 
	{ 
		istringstream ss(buffer);

		string tag, buff2;

		ss >> tag >> std::ws;

		if (tag.size() == 0 || tag[0] == '#')
			// skip empty line or comment
            continue;

		if (tag == "TOTAL_SIMULATION_TIME") 
			{ss >> TOTAL_SIMULATION_TIME; } 

		else if (tag == "DELTA_T") 
			{ss >> DELTA_T;	}

		else if (tag == "CROSSLINKDELAY")	  
			{ss >> CROSSLINKDELAY;	}

		else if (tag == "ALLOW_HARBINGERS_TO_MOVE") 
			{ss >> buff2; if (buff2=="true") ALLOW_HARBINGERS_TO_MOVE = true; else ALLOW_HARBINGERS_TO_MOVE = false; }

		else if (tag == "CAGE_ON_SIDE") 
			{ss >> buff2; if (buff2=="true") CAGE_ON_SIDE = true; else CAGE_ON_SIDE = false; }

		else if (tag == "MIN_TORQUE_TO_UPDATE")	   
			{ss >> MIN_TORQUE_TO_UPDATE;	} 

		else if (tag == "MIN_DISPLACEMENT_TO_UPDATE") 
			{ss >> MIN_DISPLACEMENT_TO_UPDATE;	} 

		else if (tag == "RECORDING_INTERVAL") 
			{ss >> RECORDING_INTERVAL;	} 

		else if (tag == "RESTORE_FROM_ITERATION") 
			{ss >> RESTORE_FROM_ITERATION;	} 

		else if (tag == "NUCLEATOR_FORCES") 
			{ss >> NUCLEATOR_FORCES;	} 

		else if (tag == "FORCE_SCALE_FACT") 
			{ss >> FORCE_SCALE_FACT;	} 

		else if (tag == "FORCE_BAR_SCALE") 
			{ss >> FORCE_BAR_SCALE;	} 

		else if (tag == "NUCLEATOR_INERTIA") 
			{ss >> NUCLEATOR_INERTIA;	} 

		else if (tag == "XLINK_NODE_RANGE") 
			{ss >> XLINK_NODE_RANGE;	} 

		else if (tag == "SEGMENT_BINS") 
			{ss >> buff2; if(buff2=="true") SEGMENT_BINS = true; else SEGMENT_BINS = false; }

		else if (tag == "DRAW_CAGE")
			{ss >> buff2; if(buff2=="true") DRAW_CAGE = true; else DRAW_CAGE = false; }

		else if (tag == "P_XLINK") 
			{ss >> P_XLINK;	} 

		else if (tag == "NUC_LINK_FORCE") 
			{ss >> NUC_LINK_FORCE;	}

		else if (tag == "NUC_LINK_BREAKAGE_FORCE") 
			{ss >> NUC_LINK_BREAKAGE_FORCE;	}

		else if (tag == "LINK_BREAKAGE_FORCE") 
			{ss >> LINK_BREAKAGE_FORCE;	}

		else if (tag == "LINK_BREAKAGE_STRAIN") 
			{ss >> LINK_BREAKAGE_STRAIN;	} 

		else if (tag == "P_LINK_BREAK_IF_OVER") 
			{ss >> P_LINK_BREAK_IF_OVER;	} 
		
		else if (tag == "LINK_FORCE") 
			{ss >> LINK_FORCE;	}
        
		else if (tag == "STICK_TO_NUCLEATOR") 
			{ss >> buff2; if(buff2=="true") STICK_TO_NUCLEATOR = true; else STICK_TO_NUCLEATOR = false; }

		else if (tag == "RESTICK_TO_NUCLEATOR") 
			{ss >> buff2; if(buff2=="true") RESTICK_TO_NUCLEATOR = true; else RESTICK_TO_NUCLEATOR = false; }
		
		else if (tag == "USETHREAD_COLLISION") 
			{ss >> buff2; if(buff2=="true") USETHREAD_COLLISION = true; else USETHREAD_COLLISION = false; }
		
		else if (tag == "USETHREAD_LINKFORCES") 
			{ss >> buff2; if(buff2=="true") USETHREAD_LINKFORCES = true; else USETHREAD_LINKFORCES = false; }
		
		else if (tag == "USETHREAD_APPLYFORCES") 
			{ss >> buff2; if(buff2=="true") USETHREAD_APPLYFORCES = true; else USETHREAD_APPLYFORCES = false; } 
    	
		else if (tag == "USE_BREAKAGE_STRAIN") 
			{ss >> buff2; if(buff2=="true") USE_BREAKAGE_STRAIN = true; else USE_BREAKAGE_STRAIN = false; } 
		
		else if (tag == "P_NUC") 
			{ss >> P_NUC; }

		else if (tag == "VISCOSITY") 
			{ss >> buff2; if(buff2=="true") VISCOSITY = true; else VISCOSITY = false; } 
		
		else if (tag == "VISCOSITY_FACTOR") 
			{ss >> VISCOSITY_FACTOR; } 
        
		else if (tag == "IMPOSED_NUC_ROT_SPEED") 
			{ss >> IMPOSED_NUC_ROT_SPEED; } 

		else if (tag == "IMPOSED_NUC_ROT") 
			{ss >> buff2; if(buff2=="true") IMPOSED_NUC_ROT = true; else IMPOSED_NUC_ROT = false; }
		
		else if (tag == "WRITE_BMPS_PRE_SYMBREAK") 
			{ss >> buff2; if(buff2=="true") WRITE_BMPS_PRE_SYMBREAK = true; else WRITE_BMPS_PRE_SYMBREAK = false; }

		else if (tag == "RADIUS") 
			{ss >> RADIUS; } 
		
		else if (tag == "CAPSULE_HALF_LINEAR") 
			{ss >> CAPSULE_HALF_LINEAR; } 
		
		else if (tag == "MAX_LINKS_PER_NODE") 
			{ss >> MAX_LINKS_PER_NODE; } 
		
		else if (tag == "NODE_REPULSIVE_MAG") 
			{ss >> NODE_REPULSIVE_MAG; } 
		
		else if (tag == "NODE_REPULSIVE_RANGE") 
			{ss >> NODE_REPULSIVE_RANGE; } 
		
		else if (tag == "LINK_TAUT_FORCE") 
			{ss >> LINK_TAUT_FORCE; } 
		
		else if (tag == "LINK_TAUT_RATIO") 
			{ss >> LINK_TAUT_RATIO; }
		
		else if (tag == "ASYMMETRIC_NUCLEATION") 
			{ss >> ASYMMETRIC_NUCLEATION; } 
		
		else if (tag == "RADIAL_SEGMENTS") 
			{ss >> RADIAL_SEGMENTS; } 
		
		else if (tag == "XLINK_NEAREST") 
			{ss >> XLINK_NEAREST; } 
		
		else if (tag == "VIEW_HEIGHT") 
			{ss >> VIEW_HEIGHT;} 
		
		else if (tag == "NODES_TO_UPDATE") 
			{ss >> NODES_TO_UPDATE;} 
		
		else if (tag == "DISTANCE_TO_UPDATE") 
			{ss >> DISTANCE_TO_UPDATE;} 
		
		else if (tag == "GAUSSFWHM") 
			{ss >> GAUSSFWHM;} 
		
		else if (tag == "SPECKLE") 
			{ss >> buff2;if(buff2=="true") SPECKLE = true;else SPECKLE = false;} 
		
		else if (tag == "SPECKLE_FACTOR") 
			{ss >> SPECKLE_FACTOR;} 
		
		else if (tag == "INIT_R_GAIN") 
			{ss >> INIT_R_GAIN;} 
		
		else if (tag == "INIT_G_GAIN") 
			{ss >> INIT_G_GAIN;} 
		
		else if (tag == "INIT_B_GAIN") 
			{ss >> INIT_B_GAIN;} 
		
		else if (tag == "ROTATION") 
			{ss >> buff2;if(buff2=="true")ROTATION = true;else ROTATION = false; } 
		
		else if (tag == "MofI") 
			{ss >> MofI;} 
		
		else if (tag == "NO_IMAGE_TEXT") 
			{ss >> buff2;if(buff2=="true") NO_IMAGE_TEXT = true; else NO_IMAGE_TEXT = false; } 
		
		else if (tag == "BMP_COMPRESSION") 
			{ss >> BMP_COMPRESSION;	if (BMP_COMPRESSION > 100) BMP_COMPRESSION = 100; else if (BMP_COMPRESSION < 0)	BMP_COMPRESSION = 0;  } 
		
		else if (tag == "BMP_OUTPUT_FILETYPE") 
			{ss >> BMP_OUTPUT_FILETYPE;	} 
		
		else if (tag == "SHAPE") 
			{ss >> buff2; if (buff2 == "CAPSULE") NUCSHAPE = nucleator::capsule; else NUCSHAPE = nucleator::sphere; }

		else 
		{
			unrecognisedlines += buffer;
			unrecognisedlines += "\n";
			continue;
		}

		// got to here so line was recognised

		istringstream ss2(buffer);

		string tag2, buff3;

		ss2 >> tag2 >> buff3 >> std::ws;

		cout << setw(30) << tag2.c_str() << setw(7) << buff3.c_str() << endl;

	 }

	param.close();

	if (unrecognisedlines.length()!=0)
	{
		cerr << endl;
		cerr << "Warning---The following lines were unrecognised:" << endl << endl;
		cerr << unrecognisedlines << endl;
		cerr << endl << endl;
	}

	// calculate commonly used constants from parameters:

	TOTAL_ITERATIONS = (int) (((double)TOTAL_SIMULATION_TIME / (double)DELTA_T)+0.5);

	//NODE_REPULSIVE_GRIDSEARCH = (int) ceil(((double) NODE_INCOMPRESSIBLE_RADIUS )/GRIDRES) + 1;
	NODE_XLINK_GRIDSEARCH = (int) ceil(((double) XLINK_NODE_RANGE )/GRIDRES) + 1;
	NODE_REPULSIVE_RANGE_GRIDSEARCH = (int) ceil(((double) NODE_REPULSIVE_RANGE )/GRIDRES) + 1;

	// loop iterations per recorded timestep
	InterRecordIterations = RECORDING_INTERVAL;
	NUMBER_RECORDINGS = int(TOTAL_ITERATIONS / RECORDING_INTERVAL);

	NODE_FORCE_TO_DIST = DELTA_T * FORCE_SCALE_FACT;
	NODE_DIST_TO_FORCE = 1.0/(DELTA_T * FORCE_SCALE_FACT);


	//RADIUS = RADIUS;//+ NODE_INCOMPRESSIBLE_RADIUS/2;

	if (SPECKLE_FACTOR<0)
	{
        SPECKLE_FACTOR = 1;
		SPECKLE = false;
		cout << "Negative SPECKLE_FACTOR reset" << endl;
	}

	if (SPECKLE)
	{	// is speckle is on, scale the gain
		INIT_R_GAIN *= (1/SPECKLE_FACTOR);
		INIT_G_GAIN *= (1/SPECKLE_FACTOR);
		INIT_B_GAIN *= (1/SPECKLE_FACTOR);
	}


	
	vector<int> postprocess_iterations;
	postprocess_iterations.resize(0);

    if (POST_PROCESS)
    {
        cout << "Postprocessing iterations: ";
	    get_postprocess_iterations(argv[2], postprocess_iterations);
    }
        
    if (!REWRITESYMBREAK  && !POST_PROCESS)
	{	// not re-writing symmetry breaking bitmaps or post-processing
		// so this is a new calculation and no other process is using temp bmp files
		// and can clear them
		sprintf(command1, "rm -f %s*.bmp 2>/dev/null", TEMPDIR);
		system(command1);
	}




#ifdef NO_IMAGEMAGICK
	cerr << "Warning: compiled with ImageMagick support turned off" << endl << endl;
#else 
	if (NO_IMAGE_TEXT)
		cerr << "Warning: image text turned off" << endl << endl;
#endif





	// create main objects
	// create as static otherwise exit() doesn't call their destructors (!)
	static actin theactin;
	static nucleator nuc_object(NUCSHAPE, &theactin);

	if (REWRITESYMBREAK)
	{
		rewrite_symbreak_bitmaps(nuc_object, theactin);
		exit(EXIT_SUCCESS);
	}

	// write out parameters to screen

	cout << "Total simulation time:      " << TOTAL_SIMULATION_TIME << endl;
	cout << "Delta_t:                    " << DELTA_T << endl;

	// and file

	unsigned int rand_num_seed;



	theactin.opruninfo.flush();


	cout << "Total iterations: " << TOTAL_ITERATIONS << endl;
	cout << "Saving snapshot every " << InterRecordIterations  
		<< " iterations (" << NUMBER_RECORDINGS << " total)" << endl;

	unsigned int starttime, endtime, lasttime ,nowtime, lastitertime;

	lasttime=0;

	int lastlinksformed = 0;
	int lastlinksbroken = 0;

	// formatting

	cout.fill(' '); 
	cout.setf(ios::fixed);
	theactin.opruninfo.fill(' ');
	theactin.opruninfo.setf(ios::fixed);

	theactin.newnodescolour.setcol(0);

	double distfromorigin = 0;

	double x_angle, y_angle, z_angle, tot_rot;

	bool DISTANCE_TO_UPDATE_reached = false;

	vect last_center, center, delta_center;

	// Breakout if we are post processing
	if( !postprocess_iterations.empty() )
	{
		postprocess(nuc_object, theactin, postprocess_iterations );
	    exit(EXIT_SUCCESS); // FIXME: early finish, move to two fcns (ML)
	}

	// - - - - - - - - - - 
	// main iteration loop
	// - - - - - - - - - -
	// initialse from a checkpoint if requested
	int starting_iter = 1;

	if (RESTORE_FROM_ITERATION != 0)
	{
		cout << "Loading data...";
		cout.flush();

	    load_data(theactin, RESTORE_FROM_ITERATION);

	    cout << "restored from iteration "
			 << RESTORE_FROM_ITERATION << endl;

        // srand( (unsigned) 200 );
	    // cout << "reseeded: " << rand() << endl;

	    starting_iter = RESTORE_FROM_ITERATION; 
		//	starting_iter = RESTORE_FROM_ITERATION + 1; // don't overwrite

        if (theactin.highestnodecount > ((int)theactin.node.size() - 1000))
            theactin.reservemorenodes(10000);
	}
	else
	{
#ifndef _WIN32

		cout << "Deleting old save files...";
		cout.flush();
		// only if starting a new calculation from scratch, clear the directories
		sprintf(command1, "rm -f %s*_0*.%s 2>/dev/null", BITMAPDIR, BMP_OUTPUT_FILETYPE.c_str() );
		system(command1);
		sprintf(command1, "rm -f %s*.wrz 2>/dev/null", VRMLDIR );
		system(command1);
		sprintf(command1, "rm -f %s*.gz 2>/dev/null", REPORTDIR );
		system(command1);
		sprintf(command1, "rm -f %s*.gz 2>/dev/null", DATADIR );
		system(command1);
		sprintf(command1, "rm -f %s*.wrl %s*.txt 2>/dev/null", TEMPDIR, TEMPDIR );
		system(command1);
        sprintf(command1, "rm -f %s 2>/dev/null", SYM_BREAK_FILE );
		system(command1);        

		cout << "done." << endl;
#else

		cout << "Deleting old save files...";
		cout.flush();
		// only if starting a new calculation from scratch, clear the directories

		sprintf(command1, "del /q %s*_0*.%s", BITMAPDIR, BMP_OUTPUT_FILETYPE.c_str() );
		system(command1);
		sprintf(command1, "del /q %s*.wrz", VRMLDIR );
		system(command1);
		sprintf(command1, "del /q %s*.gz", REPORTDIR );
		system(command1);
		sprintf(command1, "del /q %s*.gz", DATADIR );
		system(command1);
		sprintf(command1, "del /q %s*.wrl %s*.txt", TEMPDIR, TEMPDIR );
		system(command1);
        sprintf(command1, "del /q %s", SYM_BREAK_FILE );
		system(command1);        

		cout << "done." << endl;

#endif

	}

	if (argc < 3) 
	{
		cerr << "Warning: Static random number seed used" <<  endl;
		theactin.opruninfo << "Warning: Static random number seed used" <<  endl;

		rand_num_seed = 200;
	} 
	else
	{
		rand_num_seed = (unsigned)time(NULL);
		theactin.opruninfo << "Time-based random number seed (" << rand_num_seed << ") used" << endl;
		cerr << "Time-based random number seed (" << rand_num_seed << ") used" << endl;
    }

	srand(rand_num_seed);

	int filenum = 0;

    char last_symbreak_bmp_filename[255] = "";

	cout << endl;
	cout << "Starting iterations..." << endl << endl; 

	cout << "(Press 'q' to abort run, or 'm' to make quicktime movie of frames so far)" << endl << endl;

	
	if (QUIET)
	{
		cout << "Running in quiet mode (no continual progress display)" << endl;
	}
	else
	{
		cout << "Itternum|TotNode|Links+|Links-|dist   |Rotn   |T   |SaveNum" << endl << endl;

	}

	lastitertime = starttime = (unsigned) time(NULL);

	for(int i=starting_iter;i<=TOTAL_ITERATIONS;i++)
	{

		filenum = (int)(i/InterRecordIterations);

		/*if ((i % 10) == 1)
			srand( rand_num_seed );*/

		nowtime = (unsigned) time(NULL);

		if ((nowtime > lasttime) || ((i % InterRecordIterations) == 1))
		{

            if (theactin.highestnodecount > ((int)theactin.node.size() - 1000))
                theactin.reservemorenodes(10000);

#ifndef NOKBHIT
            if (kbhit())  // taken out of main loop
		    {
			    int ch = getch();
			    if (( ch=='q' ) || ( ch=='Q') || ( ch=='m' ) || ( ch=='M'))
			    {

					if (( ch=='q' ) || ( ch=='Q'))
				    {
						cout << endl << "Abort run(y/n)?";
						cout.flush();
						int ch;
						do 
						{
							ch = getch();
						} 
						while (( ch=='q' ) || ( ch=='Q'));

						if (( ch=='y' ) || ( ch=='Y'))
						{
							cout << "y - Run aborted" << endl;
							cout.flush();
							break;
						}
						else
						{
							cout << "\r";
						}
					}
	 				else
					{
						cout << endl << endl << "Directory:";
						cout.flush();

#ifndef _WIN32
						system("pwd");
#else
						system("chdir");
#endif
						cout << endl << "Making movie of frames so far..." << endl;
						cout.flush();
						cout << endl;
						system("ffmpeg -v 1 -hq -me full -qscale 1 -y -i bitmaps//x_proj_%05d.png -vcodec mpeg4 -strict 100 x_proj.mov");
						system("ffmpeg -v 1 -hq -me full -qscale 1 -y -i bitmaps//y_proj_%05d.png -vcodec mpeg4 -strict 100 y_proj.mov");
						system("ffmpeg -v 1 -hq -me full -qscale 1 -y -i bitmaps//z_proj_%05d.png -vcodec mpeg4 -strict 100 z_proj.mov");
						cout << endl;
					}
				}
		    }

#endif
            //theactin.keep_mem_resident();            

			lasttime = nowtime;
			theactin.find_center(center);
			distfromorigin = center.length();
            nuc_object.nucleator_rotation.getangles(x_angle,y_angle,z_angle);
            tot_rot = fabs(x_angle) + fabs(y_angle) + fabs(z_angle);
		
			if ((!DISTANCE_TO_UPDATE_reached) && (DISTANCE_TO_UPDATE > 0.01) 
				&& (distfromorigin > (DISTANCE_TO_UPDATE*RADIUS)))
			{
				DISTANCE_TO_UPDATE_reached = true;
				NODES_TO_UPDATE = theactin.highestnodecount;
				cout << endl << "DISTANCE_TO_UPDATE distance reached at " << distfromorigin
					<< " updating only newest " << NODES_TO_UPDATE << " nodes" << endl;
			}

			if (!QUIET)
			{
				cout << "I" << setw(7) << i 
				<< "|N"<< setw(6)<< theactin.highestnodecount
				<< "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
				<< "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
				<< "|d" << setw(6) << setprecision(3) << distfromorigin
				<< "|R" << setw(6) << setprecision(1) << (180/PI) * tot_rot
				//<< "|NR" << setw(5) <<  theactin.num_rotate	
				<< "|T" <<  setw(3) <<((unsigned) time(NULL) - lastitertime) << "\r";
			}
			cout.flush();
		}

		theactin.iterate();
		//theactin.newnodescolour.setcol((double)i/(double)TOTAL_ITERATIONS);
		
		if (((i % InterRecordIterations) == 0) && (i>starting_iter))
		{

#ifdef NON_RANDOM

srand( rand_num_seed );

#endif

			theactin.find_center(center);
			distfromorigin = center.length();
            nuc_object.nucleator_rotation.getangles(x_angle,y_angle,z_angle);
            tot_rot = fabs(x_angle) + fabs(y_angle) + fabs(z_angle);
			
			delta_center = center - last_center;  
			last_center = center;

			nuc_object.nucleator_rotation.getangles(x_angle,y_angle,z_angle);

			theactin.opvelocityinfo 
				<< (i*DELTA_T) << "," 
				<< center.x << "," 
				<< center.y << "," 
				<< center.z << "," 
				<< delta_center.length() << endl;

			if (!QUIET)
			{
				cout << "I" << setw(7) << i 
				<< "|N"<< setw(6)<< theactin.highestnodecount
				<< "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
				<< "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
				<< "|d" << setw(6) << setprecision(3) << distfromorigin
				<< "|R" << setw(6) << setprecision(1) << (180/PI) * tot_rot
				//<< "|NR" << setw(5) <<  theactin.num_rotate	
				<< "|T" <<  setw(3) <<((unsigned) time(NULL) - lastitertime);

				cout << "|S " << setw(3) <<  (int)filenum  
					<< "/" << NUMBER_RECORDINGS;

				if ((!WRITE_BMPS_PRE_SYMBREAK) && 
					((strlen(last_symbreak_bmp_filename)!=0) || (!theactin.brokensymmetry)))
					cout << "*";
				
				cout << endl;

				cout.flush();

			}

			theactin.opruninfo << "I" << setw(7) << i 
			<< "|N"<< setw(6)<< theactin.highestnodecount
			<< "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
			<< "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
            << "|d" << setw(6) << setprecision(3) << distfromorigin
            << "|R" << setw(6) << setprecision(1) << (180/PI) * tot_rot
			//<< "|NR" << setw(5) <<  theactin.num_rotate	
			<< "|T" <<  setw(3) <<((unsigned) time(NULL) - lastitertime);
			//<< "|H" <<  setw(6) << theactin.harbinger;
			theactin.opruninfo << "|S" << setw(3) <<  (int)filenum  
				<< "/" << NUMBER_RECORDINGS << endl ;

            // we don't use vrml anymore, so don't bother with it for now
			// theactin.savevrml(filenum);

			if (WRITE_BMPS_PRE_SYMBREAK)
			{
				 if (!theactin.brokensymmetry)
				 {
					nuc_object.segs.addallnodes();  // put node data into segment bins
					nuc_object.segs.set_scale_factors();

					theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);
					theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);
					theactin.savebmp(filenum, actin::zaxis, actin::runfg, true);

					cout << "\r";
					cout.flush();
				}
			}
			else
			{  // calc every 5 for scaling

				 if (((i % (5 * InterRecordIterations)) == 0) &&
					  (!theactin.brokensymmetry) && 
					  (theactin.BMP_intensity_scaling == true))
				 {
					nuc_object.segs.addallnodes();  // put node data into segment bins
					nuc_object.segs.set_scale_factors();

					// calculate but don't write bitmaps to get scaling factors
					theactin.savebmp(filenum, actin::xaxis, actin::runfg, false);
					theactin.savebmp(filenum, actin::yaxis, actin::runfg, false);
					theactin.savebmp(filenum, actin::zaxis, actin::runfg, false);

					cout << "\r";
					cout.flush();
		        
				 }

			}

			if ((!theactin.brokensymmetry) && (distfromorigin > RADIUS))
			{	// symmetry broke: set directions, scale factors etc.
				
				theactin.brokensymmetry = true;

				// set camera rotation for save  
				theactin.set_sym_break_axes();
				theactin.save_sym_break_axes();

                // calculate but don't write bitmaps to get scaling factors
			    theactin.savebmp(filenum, actin::xaxis, actin::runfg, false);
			    theactin.savebmp(filenum, actin::yaxis, actin::runfg, false);
			    theactin.savebmp(filenum, actin::zaxis, actin::runfg, false);

				if (!QUIET)
				{
					cout << "\r";
					cout.flush();
				}

                // lock the bitmap scaling factors
                theactin.BMP_intensity_scaling = false;

                // use this point to set scale factors

                nuc_object.segs.addallnodes();  // put node data into segment bins
                nuc_object.segs.set_scale_factors();                
                nuc_object.segs.save_scalefactors();
                
				// reprocess bitmaps etc.

				// call another instance to write bitmaps
				// fix this to do in same process when threading working
				cout << "Spawning new background process to write bitmaps 1--" << filenum-1 
					<< " and continuing..." << endl;

				// look for existance of this file to tell when the process has finished
				sprintf(last_symbreak_bmp_filename, "%sz_proj_%05i.%s",BITMAPDIR, 
							1 ,BMP_OUTPUT_FILETYPE.c_str());

				if (WRITE_BMPS_PRE_SYMBREAK)
				{
					// if we have already written draft bitmaps, need to delete this one
					// so we know when it's re-written
					sprintf(command1, "rm %s 1>/dev/null 2>/dev/null", last_symbreak_bmp_filename);
					system(command1);
				}

				// call self with 'sym' argument to re-write the bitmaps
				sprintf(command1, "%s sym 1>/dev/null 2>/dev/null &", argv[0]);
				system(command1);

				// kludge to move the last save
				// so we can re-load to the point we left off
				//
				//	sprintf(command1 , "gzip -q -f -9 data*.txt 2>/dev/null" );
				//	system(command1);

				//	sprintf(command1 , "mv *data*.gz %s", DATADIR);
				//	system(command1);

				//	rewrite_symbreak_bitmaps(nuc_object, theactin);
				//	load_data(theactin, i);
			}  

			save_data(theactin, i);
			
			if (theactin.brokensymmetry)
			{	// only save bitmaps if symmetry broken, 
				// 'cause we'll write the others later

				// write the bitmaps, calling imagemagick in background
				// note: this *could* cause a problem if it gets around to the next one before this
				// has finished (because we're sharing x,y,z temp bitmap files)
				// bit this is probably impossible, since we're only calling this
				// once symmetry broken, and by then things are very slow


                nuc_object.segs.addallnodes();  // put node data into segment bins

                nuc_object.segs.savereport(filenum);
                nuc_object.segs.saveSDreport(filenum);                
		        nuc_object.segs.saveradialreport(filenum);
                nuc_object.segs.saveradialaxisreport(filenum, 0);
                nuc_object.segs.saveradialaxisreport(filenum, 1);
                nuc_object.segs.saveradialaxisreport(filenum, 2);

				theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);
				theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);
				theactin.savebmp(filenum, actin::zaxis, actin::runbg, true);

				cout << "\r";
				cout.flush();

				if (strlen(last_symbreak_bmp_filename)!=0)
				{
					ifstream ip_last_symbreak_bmp_file( last_symbreak_bmp_filename );

					if (ip_last_symbreak_bmp_file)
					{
						cout << "Finished background processing of symmetry breaking bitmaps               " << endl;

						*last_symbreak_bmp_filename = 0 ;

						ip_last_symbreak_bmp_file.close();
					}
				}

			}

			theactin.clear_node_stats();  // clear the stats data in the nodes

			theactin.compressfilesdowork(filenum);

			lastlinksformed = theactin.linksformed;
			lastlinksbroken = theactin.linksbroken;

			//theactin.num_rotate = 0;
			//theactin.num_displace = 0;

			theactin.opruninfo.flush();

			lastitertime = (unsigned) time(NULL);
		}
	}

	cout << endl << "Done " << endl << endl;

	theactin.saveinfo();

	endtime = (unsigned) time(NULL);

	cout << endl << "Time : " 
		<< ((endtime-starttime) / 3600) << "h " 
		<< ((endtime-starttime) / 60) % 60 << "m " 
		<< ((endtime-starttime) % 60) << "s " << endl;

	theactin.opruninfo << endl << "Time : " 
		<< ((endtime-starttime) / 3600) << "h " 
		<< ((endtime-starttime) / 60) % 60 << "m " 
		<< ((endtime-starttime) % 60) << "s " << endl;;	
	
	exit(EXIT_SUCCESS);
}


string get_datafilename(const int iteration)
{
	//stringstream filenamestr;
    char filename[255];
	sprintf(filename , "data_%07i.txt", iteration);
    return filename;
    
}

int load_data(actin &theactin, int iteration)
{

    char command1[255];
    string filename = get_datafilename(iteration);

	filename = DATADIR + filename;
    
    sprintf(command1, "gunzip %s",filename.c_str());
    system(command1);
    
    ifstream ifstrm( filename.c_str() );
    if(!ifstrm) {
	cout << "Unable to open file " << filename << " for input";
	return 1;
    }
    
    string str;    
    // load header
    ifstrm >> str;
    // ensure the identifier for the start of the actin
    if(str.compare("comet:") !=0 ){
	cout << "error in checkpoint file, 'comet:' expected" << endl;
	return 1;
    }
    
    int saved_iteration;
    ifstrm >> saved_iteration;
    theactin.load_data(ifstrm);
    
    ifstrm.close();
    
    sprintf(command1, "gzip %s",filename.c_str());
    system(command1);
    
    // check the iteration is correct
    if( saved_iteration != iteration ){
	cout << "error in saved file, saved iteration." << endl;
	return 1;
    }
    return iteration;
}

int save_data(actin &theactin, int iteration)
{
    string filename = get_datafilename(iteration);

	filename = TEMPDIR + filename;

    ofstream ofstrm( filename.c_str() );    
    if(!ofstrm) {
	cout << "Unable to open file " << filename;
	return 1;
    } 
    
    // write out a header  and save iteration
    ofstrm << "comet:" << endl
	   << iteration << endl;
    
    // actin does all the real work
    theactin.save_data(ofstrm);
    
    ofstrm.close();
    
    return 0;
}

void get_postprocess_iterations(const char *iterdesc, vector<int> &postprocess_iterations)
{
    // allow either:
    //  a single value        '12'
    //  a matlab style vector '1,2,3,4,5' | '0:5:50'
    char numbers[]    = "1234567890";
    char alldelim[]  = ",:"; // make big enough for the strcat (sorry)
    char rangedelim[] = ":";
    
    // Good old fashioned string handling...
    // why didn't I use std::string? (ML)
    
    // tokenise the description
    char* tokch = (char*) malloc(sizeof(char) * ( strlen(iterdesc)+1 ));
    strcpy(tokch, iterdesc);
    char* chptr = strtok(tokch, alldelim);
    while(chptr != NULL) {
	// ensure ints, no rubbish
	if(strspn(chptr, numbers) != strlen(chptr))
	{
	    cout << "Post processing iteration description not recognised (must be +ve int): "
		 << chptr << endl;
	    exit(EXIT_FAILURE);
	}
	postprocess_iterations.push_back( atoi(chptr) );
	chptr = strtok(NULL, alldelim);
    }
    free(tokch);
    
    // build the range if needed
    if(strcspn(iterdesc, rangedelim) != strlen(iterdesc) ){
	// : range
	int start = 0, step = 0, end = 0;
	if(postprocess_iterations.size() == 1) {
	    cerr << "ERROR: too few values in range request (2 min):" << iterdesc << endl;
	    exit(EXIT_FAILURE);
	} else if(postprocess_iterations.size() == 2) {
	    start = postprocess_iterations[0] * InterRecordIterations;
	    step = InterRecordIterations;
	    end = postprocess_iterations[1] * InterRecordIterations;
	} else if(postprocess_iterations.size() == 3) {
	    start = postprocess_iterations[0] * InterRecordIterations;
	    step = postprocess_iterations[1] * InterRecordIterations;
	    end = postprocess_iterations[2] * InterRecordIterations;	    
	} else if(postprocess_iterations.size() > 3) {
	    cerr << "ERROR: too many values in range request (3 max):" << iterdesc << endl;
	    exit(EXIT_FAILURE);
	}
	postprocess_iterations.resize(0);
	//postprocess_iterations.clear();
	for(int i= start; i<=end; i+=step){
	    postprocess_iterations.push_back(i);
	}
    } 
}

// We may want to have a vtk postprocess as a seperate option?
// otherwise I've lumped it with the bitmap processing here
// the frame counting is seperated, but consider unifying 
// this once we are happy with how it all works.
void postprocess(nucleator& nuc_object, actin &theactin, vector<int> &postprocess_iterations)
{
    int filenum;
    theactin.load_sym_break_axes();
    
    if(nuc_object.segs.load_scalefactors() ){
	// if we're able to load the scale factors
	// turn off the auto-scaling
	theactin.BMP_intensity_scaling = false;
    }
    
    // vtk
    CometVtkVis vtkvis(&theactin);
    
    cout << "Post processing " 
	 << postprocess_iterations.size()
	 <<  " data sets." << endl;
    
    int frame = (int) postprocess_iterations.size() - 1;
    for(vector<int>::reverse_iterator iteration = postprocess_iterations.rbegin(); 
	iteration != postprocess_iterations.rend(); ++iteration) {
	
	filenum = (int)(*iteration/InterRecordIterations);

	cout << "Post processing iteration: " << *iteration << " file " << filenum 
	     << " (" << (postprocess_iterations.size() - frame) << "/" 
	     << postprocess_iterations.size() << ")"
	     << endl; 
	
	load_data(theactin, *iteration);
	
	theactin.load_sym_break_axes();   // overwrite rotation matrixes	
	nuc_object.segs.addallnodes();  // put node data into segment bins
	
	cout << "- reports:" << endl;
        nuc_object.segs.savereport(filenum);
        nuc_object.segs.saveSDreport(filenum);
	nuc_object.segs.saveradialreport(filenum);       
        nuc_object.segs.saveradialaxisreport(filenum, 0);
        nuc_object.segs.saveradialaxisreport(filenum, 1);
        nuc_object.segs.saveradialaxisreport(filenum, 2);
	
	cout << "- bitmaps:" << endl;
	theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);  // was bg
	theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);  // was bg
	theactin.savebmp(filenum, actin::zaxis, actin::runfg, true);
	
	cout << endl << "- visualisation: " << filenum << endl;
	vtkvis.buildVTK(filenum);

	frame--;
    }
}

// void postprocess(nucleator& nuc_object, actin &theactin, 
// 		 vector<int> &postprocess_iterations)
// {
//     int filenum;
//     theactin.load_sym_break_axes();
//    
//     cout << endl
// 	 << "-- Post processing" << endl;
//    
//     if (nuc_object.segs.load_scalefactors())
//     {
// 	// if we're able to load the scale factors
// 	// turn off the auto-scaling
// 	theactin.BMP_intensity_scaling = false;
//     }
//     int i=1;
//     for(vector<int>::reverse_iterator iteration = postprocess_iterations.rbegin();
// 	iteration != postprocess_iterations.rend(); ++iteration, ++i)
//     {
// 	filenum = (int)(*iteration/InterRecordIterations);
//	
// 	cout << " iter " << i << "/" << postprocess_iterations.size() // Mark, change back if you prefer 
// 	     << " (frame " << filenum 
// 	     << "): " ;
//	
// 	load_data(theactin, *iteration);
//	
// 	theactin.load_sym_break_axes();   // overwrite rotation matcp rixes	
// 	nuc_object.segs.addallnodes();  // put node data into segment bins
//	
//         nuc_object.segs.savereport(filenum);
//         nuc_object.segs.saveSDreport(filenum);
// 	nuc_object.segs.saveradialreport(filenum);
//         nuc_object.segs.saveradialaxisreport(filenum, 0);
//         nuc_object.segs.saveradialaxisreport(filenum, 1);
//         nuc_object.segs.saveradialaxisreport(filenum, 2);
//	
// 	// run them in foreground to slow things down
// 	// so don't overload system with too many bg processes
//	
// 	theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);  // was bg
// 	theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);	// was bg
// 	theactin.savebmp(filenum, actin::zaxis, actin::runfg, true);
//        
// 	cout << "\r";
// 	cout.flush();
//     }
//     cout << endl << "-- done." << endl; 
// }

// r90 postprocess, Mark D remove if you are happy with minor changes above 
//
// void postprocess(nucleator& nuc_object, actin &theactin, vector<int> &postprocess_iterations)
// {
//
// 	int filenum;
//     theactin.load_sym_break_axes();
//
// 	if (nuc_object.segs.load_scalefactors())
//     {
//         // if we're able to load the scale factors
//         // turn off the auto-scaling
//         theactin.BMP_intensity_scaling = false;
//     }
//
//     for(vector<int>::iterator iteration = postprocess_iterations.end()-1; 
// 		iteration >= postprocess_iterations.begin(); ++iteration)
// 	{
// 		filenum = (int)(*iteration/InterRecordIterations);
//
// 		//cout << "Post processing iteration " << *iteration << ": ";
// 		cout << "Post processing frame " << filenum << "/" 
// 			 << ((*postprocess_iterations.end())/InterRecordIterations) << ": ";
//
// 		load_data(theactin, *iteration);
//
// 		theactin.load_sym_break_axes();   // overwrite rotation matcp rixes
//
// 		// cout << theactin.camera_rotation;
//				
// 		nuc_object.segs.addallnodes();  // put node data into segment bins
//
//         nuc_object.segs.savereport(filenum);
//         nuc_object.segs.saveSDreport(filenum);
// 		nuc_object.segs.saveradialreport(filenum);
//         nuc_object.segs.saveradialaxisreport(filenum, 0);
//         nuc_object.segs.saveradialaxisreport(filenum, 1);
//         nuc_object.segs.saveradialaxisreport(filenum, 2);
//
// 		// run them in foreground to slow things down
// 		// so don't overload system with too many bg processes
//
// 		theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);  // was bg
// 		theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);	// was bg
// 		theactin.savebmp(filenum, actin::zaxis, actin::runfg, true);
//      	
// 		cout << "\r";
// 		cout.flush();
//     }
// }

void rewrite_symbreak_bitmaps(nucleator& nuc_object, actin &theactin)
{

	int filenum;
    theactin.load_sym_break_axes();

	if (nuc_object.segs.load_scalefactors())
    {
        // if we're able to load the scale factors
        // turn off the auto-scaling
        theactin.BMP_intensity_scaling = false;
    }

    //cout << " InterRecordIterations " << InterRecordIterations
	//	 << " theactin.symbreakiter " << theactin.symbreakiter << endl;

    // go in reverse order, in case we're autoscaling

    for(int i = theactin.symbreakiter - InterRecordIterations; i > 0; i-=InterRecordIterations)
	{
		filenum = (int)(i/InterRecordIterations);

		//cout << "Post processing iteration " << *iteration << ": ";
		cout << "Writing bitmaps for frame " << filenum << "/" 
			 << (theactin.symbreakiter/InterRecordIterations)-1 << ": ";
		cout.flush();

		load_data(theactin, i);

		theactin.load_sym_break_axes();   // overwrite rotation matrixes

		// cout << theactin.camera_rotation;
					
		nuc_object.segs.addallnodes();  // put node data into segment bins

        nuc_object.segs.savereport(filenum);
        nuc_object.segs.saveSDreport(filenum);
		nuc_object.segs.saveradialreport(filenum);
        nuc_object.segs.saveradialaxisreport(filenum, 0);
        nuc_object.segs.saveradialaxisreport(filenum, 1);
        nuc_object.segs.saveradialaxisreport(filenum, 2);

		// run them in foreground to slow things down
		// so don't overload system with too many bg processes

		theactin.savebmp(filenum, actin::xaxis, actin::runbg, true);  // was bg
		theactin.savebmp(filenum, actin::yaxis, actin::runbg, true);	// was bg
		theactin.savebmp(filenum, actin::zaxis, actin::runfg, true);
		
		cout << "\r";
		cout.flush();
    }
	cout << endl;
}



/*


// Newtonian fluid drag
//
// F=6*Pi*a*eta*v
// where a = radius, eta = viscosity, v = velocity
// let it act on each node for now (acting on each
// inter-node link may be more appropriate)

// Elasticity of actin network
//
// ?? who knows?  Likely nonlinear and asymmetric.
// Probably very steep for tension, but shallow for compression
// up to a critical compressibility, then steep again
// This function will be critical.
//
// Start with very simple approximation:  
// forceonnode = E * ( actualdist - linkdist ) / linkdist


// linking:
//
// Let the linking be isotropic for now.
// argument can be made that entanglement will make effective
// linking largely isotropic even if filaments are aligned.
// However, alignment strongly implies that 
// the elasticity of the network will not
// be isotropic.  
// Deal with this later, after the first pass.

// Geometry
//
// let the bead be spherical for now

*/
