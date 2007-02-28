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

#include "mytimer.h"

//#include "consts.h"

class tracknodeinfo;

// vtk
#ifdef LINK_VTK
  #include "comet_vtk.h"
#else
  #include "comet_novtk.h"
#endif

const int default_nice_level = 15;   // v. low priority, but above screesaver :)

// these are fixed compile time limits, so make sure they are well above what will be used (else will core dump)

double GRIDBOUNDS =  50.0;	  // size of grid in um (i.e. bead can move half of this from origin)
double GRIDRES    =   0.8;	  // low res grid range
int MAXNODES = 50000;

// calc GRIDSIZE from above

short int GRIDSIZE = (int) (2*GRIDBOUNDS/GRIDRES);

// these are the default values for parameters that can be overridden by the cometparams.ini file

double TOTAL_SIMULATION_TIME = 20000;  
double DELTA_T = 0.1;	

bool NUCLEATOR_FORCES = true;

double GAUSSFWHM =  0.266;

int BMP_WIDTH  = 800;
int BMP_HEIGHT = 592; // divisible by 16 for movie compression
int BMP_AA_FACTOR = 1;

int VTK_WIDTH  = 1024;
int VTK_HEIGHT = 768;
int VTK_AA_FACTOR = 1;

double VTK_LINK_COLOUR_GAMMA = 1.8;
bool VTK_MOVE_WITH_BEAD = true;
double VTK_MIN_PLOT_LINK_FORCE_PCT = 0.0;

double BMP_INTENSITY_SCALE = 1.4;

double INIT_R_GAIN = 20;
double INIT_G_GAIN = 20;
double INIT_B_GAIN = 20;

double FOCALDEPTH = 2.0;

double BMP_INTENSITY_OFFSET = 0.24;

double VTK_VIEWANGLE = 40;

bool POST_BMP = true;
bool POST_VTK = false;
bool POST_STATS = false;
bool POST_REPORTS = false;
bool BMP_TRACKS = false;

bool PLOTFORCES = true;

bool ALLOW_HARBINGERS_TO_MOVE = false;

bool   SPECKLE = false;
double SPECKLE_FACTOR = 1;
bool SPECKLEGRID = false;
double SPECKLEGRIDPERIOD = 1.0;
double SPECKLEGRIDTIMEWIDTH = 0.1;
double SPECKLEGRIDSTRIPEWIDTH = 10.0;
bool POLY_FEEDBACK = false;
double POLY_FEEDBACK_DIST = 1.0;
double POLY_FEEDBACK_MIN_PROB = 0.05;
double POLY_FEEDBACK_FACTOR = 4;


bool ROTATION = true;

bool DRAW_CAGE = false;

bool SEGMENT_BINS = true;

bool CAGE_ON_SIDE = true;

bool BMP_FIX_BEAD_MOVEMENT = false;
bool BMP_FIX_BEAD_ROTATION = false;

bool COL_NODE_BY_STRAIN = false;
bool COL_LINK_BY_DIRN = false;
bool COL_INDIVIDUAL_NODES = false;

double MOFI =  0.1;

bool X_BMP = true;
bool Y_BMP = true;
bool Z_BMP = true;                 

double TRACK_MIN_RANGE = 0.0;  // frame numbers, can be floating point
double TRACK_MAX_RANGE = 1.0;

int TRACKFRAMESTEP = 5;

double gridscanjitter = 0.01;

int SAVE_DATA_PRECISION	= 4;


//bool FORCES_ON_SIDE = true;

bool NO_IMAGE_TEXT = false;
int BMP_COMPRESSION = 95;
string BMP_OUTPUT_FILETYPE = "png";

//int REPORT_AVERAGE_ITTERATIONS = 50;

int RESTORE_FROM_FRAME = 0; // =0 don't load a checkpoint 
int RECORDING_INTERVAL = 0;
int TOT_FRAMES = 400;

double FORCE_SCALE_FACT = 0.001;	// convert forces (nom in pN) into node displacements (nom in uM)
					        // this is related to effective viscosity and effective size of node
double FORCE_BAR_SCALE   = 1;	// scale for relating force to image bar

bool VARY_P_XLINK = false;

double XLINK_NODE_RANGE =  1.0;		// Limit crosslink to within this range

double LINK_BREAKAGE_FORCE =  2;	 // breakage force per link
bool USE_BREAKAGE_VISCOSITY = false;
double BREAKAGE_VISCOSITY_THRESHOLD = 1000;
double LINK_BREAKAGE_STRAIN = 1.15;
//double P_LINK_BREAK_IF_OVER =  0.25;  // probablility that force will break link if over the link breakage force
unsigned int MAX_LINKS_PER_NEW_NODE = 100;
unsigned int MAX_LINK_ATTEMPTS = 20;

bool USE_BROWNIAN_FORCES = false;
double BROWNIANFORCESCALE = 0.01;

double NUCLEATOR_INERTIA = 10;
bool VARY_INERT_W_RAD = false;

double NUC_LINK_FORCE = 0.25;
double NUC_LINK_BREAKAGE_FORCE = 2;

//double LINK_TAUT_FORCE =  5;
//double LINK_TAUT_RATIO =  1.1;

double IMPOSED_NUC_ROT_SPEED = 1;
bool   IMPOSED_NUC_ROT = false;

bool   TEST_SQUASH = false;
double TEST_FORCE_INITIAL_MAG = 0;
double TEST_FORCE_INCREMENT = 10;
double TEST_DIST_EQUIL = 0.0001;

bool WRITE_BMPS_PRE_SYMBREAK = false;
bool SYM_BREAK_TO_RIGHT = false;

bool QUIET = false;

bool VISCOSITY = false;
double NON_VISC_WEIGHTING = 1.0;
double MAX_VISC_WEIGHTING = 10.0;
double VISCOSITY_FACTOR = 0.1;
double VISCOSITY_EDGE_FACTOR = 4.0;
double VISC_DIST = 0.7;
//double VISCOSITY_UNWEIGHTING_FACTOR = 100;

//char temp_BMP_filename[1024];

bool STICK_TO_NUCLEATOR = false;
bool RESTICK_TO_NUCLEATOR = true;



double LINK_FORCE = 0.1;
double P_XLINK =  0.5;
double P_NUC =  0.08;
double MAX_POLYMERISATION_PRESSURE = 1000;

double RADIUS =  1.0;
double CAPSULE_HALF_LINEAR =  6.0;
double COVERSLIPGAP = 10000;

nucleator::shape NUCSHAPE = nucleator::sphere;  //default to sphere

int TOTAL_ITERATIONS ;
int NODE_REPULSIVE_GRIDSEARCH ;
int NODE_XLINK_GRIDSEARCH ;
int NODE_REPULSIVE_RANGE_GRIDSEARCH;

char VRMLDIR[2048];
char DATADIR[2048];
char REPORTDIR[2048];
char BITMAPDIR[2048];
char TEMPDIR[2048];
char VTKDIR[2048];
char STATSDIR[2048];

char IMAGEMAGICKCONVERT[1024];
char IMAGEMAGICKMOGRIFY[1024];

double LINK_POWER_SCALE = 0;

bool COMPRESSDATAFILES = true;

int RADIAL_SEGMENTS = 12;
int NODES_TO_UPDATE = 300000;  //only update the NODES_TO_UPDATE newest nodes
double DISTANCE_TO_UPDATE = 0;

int CROSSLINKDELAY = 200;  // number of interations before crosslinking 
						  //  (to allow position to be equilibrated to something
						  //       reasonable before locking node in place)

double NODE_REPULSIVE_POWER = 2.7;
double NODE_REPULSIVE_MAG =  1.5;
double NODE_REPULSIVE_RANGE = 1.0;
double NODE_REPULSIVE_BUCKLE_RANGE = NODE_REPULSIVE_RANGE * 0.3;
double NODE_REPULSIVE_BUCKLE_MAG = NODE_REPULSIVE_MAG * 0.7;
double NODE_REPULSIVE_BUCKLE_TO	= 0.2;

double NODE_FORCE_TO_DIST;
double NODE_DIST_TO_FORCE;

bool NOBGIMAGEMAGICK = false;

int ASYMMETRIC_NUCLEATION = 0;

int XLINK_NEAREST = 1;

double VIEW_HEIGHT = 12;

bool USE_THREADS;
int NUM_THREADS;
int NUM_THREAD_DATA_CHUNKS;

//-- ThreadTaskTeam
TaskQueue thread_queue;
bool USETHREAD_COLLISION   = true;
bool USETHREAD_APPLYFORCES = true;
bool USETHREAD_LINKFORCES  = true;
// --

bool CLUSTER = false;

pthread_attr_t thread_attr;
vector<pthread_t> threads;

vector<struct thread_data>  collision_thread_data_array;
vector<struct thread_data>  linkforces_thread_data_array;
vector<struct thread_data>  applyforces_thread_data_array;


// these variables need to be static/global for sharing across threads:


#ifdef NODE_GRID_USE_ARRAYS
	NG1d* __restrict nodegrid;
#else
	//Nodes3d nodegrid;
    NG4d nodegrid;
#endif

//unsigned int GRIDSIZE;

actin *ptheactin;

#ifdef USE_MERSENNE

#include "MersenneTwister.h"

vector <MTRand> mers_rand;

#endif

vector <nodes>	actin::node;
vector <bool>   actin::donenode;	
Nodes2d actin::nodes_by_thread;
Nodes2d actin::recti_near_nodes;
Nodes2d actin::nodes_on_same_gridpoint;
vector<vector<NODEGRIDTYPE<nodes*>*> > actin::gridpointsbythread;
vector <int> actin::nearby_collision_gridpoint_offsets;


rotationmatrix actin::torque_rotate;
//vect actin::nuc_disp;
int actin::lowestnodetoupdate;

vector<int>::iterator actin::nearby_collision_gridpoint_offset_begin;
vector<int>::iterator actin::nearby_collision_gridpoint_offset_end;

//Nodes1d actin::nodes_within_nucleator;

//vector <int> actin::recti_near_nodes_size;
//vector <int> actin::nodes_on_same_gridpoint_size;

int actin::iteration_num;

bool actin::isinthread;
Nodes2d actin::linkremovefrom;
Nodes2d actin::linkremoveto;

// set per instance, for when forks

bool REWRITESYMBREAK = false;
bool POST_PROCESS = false;
bool POST_PROCESS4CPU = false;
bool VIEW_VTK = false;

int POST_PROC_ORDER = +1;  // +1 = forward, -1 = reverse;

int InterRecordIterations = 0;

bool DISTANCE_TO_UPDATE_reached = false;
bool finished_writing_sym_bitmaps = false;

// functions

bool load_data(actin &theactin, int iteration, const bool &loadscale);
int save_data(actin &theactin, int iteration);
string get_datafilename(const int iteration);
void get_postprocess_iterations(const char *iterdesc, vector<int> &postprocess_iterations, const int& lastframedone);
void postprocess(nucleator& nuc_object, actin &theactin, vector<int> &postprocess_iterations, char* argv[]);
void rewrite_symbreak_bitmaps(nucleator& nuc_object, actin &theactin);

// keyboard monitoring function

#ifndef NOKBHIT
    #ifndef _WIN32
	    #include "kbhit.h"
	    #define kbhit keyb.kbhit
        #define getch keyb.getch
    #else
	    #include <conio.h>
    #endif
#endif



string strtoupper(string str)
{
   for(unsigned int i=0; i != str.length(); i++)
   {
      str[i] = toupper(str[i]);
   }
   return str;//return the converted string
}
 
string strtolower(string str)
{
   for(unsigned int i=0; i != str.length(); i++)
   {
      str[i] = tolower(str[i]);
   }
   return str;//return the converted string
}

// main 

int main(int argc, char* argv[])
{

    cout << endl; 

    if(argc < 2 || argc > 4) 
	{
	    cerr << "Usage:" << endl << endl << argv[0] << " numThreads" << endl << endl;
	    cerr << "where numThreads is the number of threads to use" << endl;
	    cerr << "Set numThreads to 1 to run in single threaded mode" << endl;
	    cerr << "or to postprocess, specify files in terms of filenumber" 
		 << endl << "./comet post 1:2:10" << endl
         << "./comet post 0:0  <- use 0:0 to process all frames" << endl
         << "./comet p4 0:0  <- use p4 to run 4 post-process threads simultaneously" << endl;
	    cerr << endl;

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

		if (strcmp(argv[1], "p4") == 0 )
		{
			POST_PROCESS = true;
			POST_PROCESS4CPU = true;
		}

		if (strcmp(argv[2], "q") == 0 || strcmp(argv[2], "Q") == 0)
			QUIET = true;

	} 
		
	if (argc > 3) 
	{
		if (strcmp(argv[3], "q") == 0 || strcmp(argv[3], "Q") == 0)
			QUIET = true;
	}

    bool NOBITMAPS = false;  

    sprintf(IMAGEMAGICKCONVERT,"convert");
    sprintf(IMAGEMAGICKMOGRIFY,"mogrify");



    

	char hostname[1024]="";

#ifndef _WIN32

	gethostname( hostname, 255);

#endif

    // is the machine part of a rocks cluster?  (if so, we'll set output to quiet etc. later on)
    if (    (strcmp( hostname, "ec3.ucsf.edu")       == 0) ||
            (strcmp( hostname, "sc1.ucsf.edu")       == 0) ||
            (strcmp( hostname, "compute-0-0.local")  == 0) ||
            (strcmp( hostname, "compute-0-1.local")  == 0) ||
            (strcmp( hostname, "compute-0-2.local")  == 0) ||
            (strcmp( hostname, "compute-0-3.local")  == 0) ||
            (strcmp( hostname, "compute-0-4.local")  == 0) ||
            (strcmp( hostname, "compute-0-5.local")  == 0) ||
            (strcmp( hostname, "compute-0-6.local")  == 0) ||
            (strcmp( hostname, "compute-0-7.local")  == 0) ||                  
            (strcmp( hostname, "compute-0-8.local")  == 0) ||
            (strcmp( hostname, "compute-0-9.local")  == 0) ||
            (strcmp( hostname, "compute-0-10.local") == 0) ||
            (strcmp( hostname, "compute-0-11.local") == 0) ||
            (strcmp( hostname, "compute-0-12.local") == 0) ||
            (strcmp( hostname, "compute-0-13.local") == 0) ||
            (strcmp( hostname, "compute-0-14.local") == 0) ||
            (strcmp( hostname, "compute-0-15.local") == 0) )
    {
        CLUSTER=true;
    }


    int nicelevel = 0;

	if (!POST_PROCESS && !REWRITESYMBREAK)
	{	// don't drop priority if re-writing bitmaps
		// because calling thread already has low priority
		// and this would drop it further
		// so process would halt on a single cpu machine
        // while the main process was running

		if (strcmp( hostname, "guanine.ucsg.edu") == 0)
		{
			nicelevel = 0;
		} else if (CLUSTER)
        {
            nicelevel = default_nice_level;
        }
		else
		{
			nicelevel = default_nice_level;
        }
    }

#ifndef _WIN32

	nice(nicelevel);

#endif


  	NUM_THREADS = atoi(argv[1]);

	if (NUM_THREADS < 2)
	{
		USE_THREADS = false;
		NUM_THREADS = 1;		
	}
	else
	{
		USE_THREADS = true;
	}

	NUM_THREAD_DATA_CHUNKS = NUM_THREADS;

    // create threads:
	// -- Threading TaskTeam, create and intialise the team
	if (USE_THREADS  && !POST_PROCESS && !REWRITESYMBREAK)
	{
	    thread_queue.create_threads(NUM_THREADS);
	}

	if (USE_THREADS)
		cout << "Running in multithreaded mode with " << NUM_THREADS << " threads (including parent)" << endl;
	else
		cout << "Running in Single Threaded mode" << endl;


	  // was NUM_THREADS*4

	// data for threads (managed outside queue
	collision_thread_data_array.resize(NUM_THREAD_DATA_CHUNKS);
	linkforces_thread_data_array.resize(NUM_THREAD_DATA_CHUNKS);
	applyforces_thread_data_array.resize(4 * NUM_THREAD_DATA_CHUNKS);



	// make directories

	char command1[1024];

	cout << "Working Directory: ";
	cout.flush();

#ifndef USEWINDOWSCOMMANDS
    system("pwd");
#else
	system("chdir");
#endif

	
#ifndef NOKBHIT
#ifndef _WIN32
    keyboard keyb;
#endif
#endif

    
    char path[2048];
    getcwd(path, 2048);
    
#ifndef USEWINDOWSCOMMANDS

	sprintf(VRMLDIR,"%s/vrml/", path);
    sprintf(DATADIR,"%s/data/", path);
    sprintf(REPORTDIR,"%s/reports/", path);
    sprintf(BITMAPDIR,"%s/bitmaps/", path);
    //sprintf(TEMPDIR,"%s/temp/", path);
    sprintf(VTKDIR,"%s/vtk/", path);
    sprintf(STATSDIR,"%s/statistics/", path);
    
    sprintf(TEMPDIR,"/tmp/comet.XXXXXX");
    mkdtemp(TEMPDIR);
    sprintf(TEMPDIR,"%s/",TEMPDIR);

    cout << "Temp dir: " << TEMPDIR << endl;

#else
	
    sprintf(VRMLDIR,"%s\\vrml\\", path);
    sprintf(DATADIR,"%s\\data\\", path);
    sprintf(REPORTDIR,"%s\\reports\\", path);
    sprintf(BITMAPDIR,"%s\\bitmaps\\", path);
    sprintf(TEMPDIR,"%s\\temp\\", path);
    sprintf(VTKDIR,"%s\\vtk\\", path);
    sprintf(STATSDIR,"%s\\statistics\\", path);

    

#endif


#ifndef USEWINDOWSCOMMANDS
	
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
	sprintf(command1, "mkdir %s 2>/dev/null", VTKDIR  );
	system(command1);
	sprintf(command1, "mkdir %s 2>/dev/null", STATSDIR  );
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
	sprintf(command1, "mkdir %s", VTKDIR  );
	system(command1);
	sprintf(command1, "mkdir %s", STATSDIR  );
	system(command1);
#endif

    bool ABORT = false;
			
    int RAND_SEED = -1; 

	// main parameters:

	// read the parameters file:

	cout << endl << "Parsing parameters file..." << endl << endl;

	ifstream param(COMET_PARAMS_FILE); 
	if(!param) 
	{
		cerr << "Cannot open " << COMET_PARAMS_FILE << endl << endl;
		exit(EXIT_FAILURE);
	}

	string buffer;
	string unrecognisedlines;

    double tempdbl;

	while (getline(param, buffer)) 
	{ 
        istringstream ss(strtoupper(buffer));

		string tag, buff2, buff3;

		ss >> tag >> std::ws;

		if (tag.size() == 0 || tag[0] == '#')
			// skip empty line or comment
            continue;

		if (tag == "TOTAL_SIMULATION_TIME") 
			{ss >> TOTAL_SIMULATION_TIME;}

        else if (tag == "TOT_FRAMES") 
			{ss >> TOT_FRAMES;}

		else if (tag == "RAND_SEED") 
			{ss >> RAND_SEED;}

		else if (tag == "SAVE_DATA_PRECISION") 
			{ss >> SAVE_DATA_PRECISION;}

		else if (tag == "DELTA_T") 
			{ss >> DELTA_T;}

		else if (tag == "GRIDBOUNDS") 
			{ss >> GRIDBOUNDS;}  

		else if (tag == "MAXNODES") 
			{ss >> MAXNODES;}

        else if (tag == "COMPRESSDATAFILES") 
			{ss >> buff2; if (buff2=="TRUE") COMPRESSDATAFILES = true; else COMPRESSDATAFILES = false;}

		else if (tag == "BMP_WIDTH")	  
			{ss >> BMP_WIDTH;}

		else if (tag == "BMP_HEIGHT")	  
			{ss >> BMP_HEIGHT;}

        else if (tag == "BMP_INTENSITY_SCALE")	  
			{ss >> BMP_INTENSITY_SCALE;}

		else if (tag == "VTK_WIDTH")	  
			{ss >> VTK_WIDTH;}

		else if (tag == "VTK_HEIGHT")	  
			{ss >> VTK_HEIGHT;}         

        else if (tag == "SYM_BREAK_TO_RIGHT") 
			{ss >> buff2; if (buff2=="TRUE") SYM_BREAK_TO_RIGHT = true; else SYM_BREAK_TO_RIGHT = false;}

		else if (tag == "VTK_MOVE_WITH_BEAD") 
			{ss >> buff2; if (buff2=="TRUE") VTK_MOVE_WITH_BEAD = true; else VTK_MOVE_WITH_BEAD = false;}

        else if (tag == "VTK_MIN_PLOT_LINK_FORCE_PCT") 
			{ss >> VTK_MIN_PLOT_LINK_FORCE_PCT;}

        else if (tag == "VIEW_VTK")  
			{ss >> buff2; if (buff2=="TRUE") VIEW_VTK = true; else VTK_MOVE_WITH_BEAD = false;}

		else if (tag == "VTK_AA_FACTOR")	  
			{ss >> VTK_AA_FACTOR;}

        else if (tag == "BMP_AA_FACTOR")	  
			{ss >> BMP_AA_FACTOR;}

		else if (tag == "VTK_LINK_COLOUR_GAMMA")	  
			{ss >> VTK_LINK_COLOUR_GAMMA;}

		else if (tag == "VTK_VIEWANGLE")	  
			{ss >> VTK_VIEWANGLE;}

		else if (tag == "CROSSLINKDELAY")	  
			{ss >> CROSSLINKDELAY;}
		
		else if (tag == "COVERSLIPGAP")	  
			{ss >> COVERSLIPGAP;}

		else if (tag == "ALLOW_HARBINGERS_TO_MOVE") 
			{ss >> buff2; if (buff2=="TRUE") ALLOW_HARBINGERS_TO_MOVE = true; else ALLOW_HARBINGERS_TO_MOVE = false;}

        else if (tag == "POST_PROC_ORDER") 
			{ss >> POST_PROC_ORDER;}

		else if (tag == "POST_BMP") 
			{ss >> buff2; if (buff2=="TRUE") POST_BMP = true; else POST_BMP = false;}   

		else if (tag == "POST_VTK") 
			{ss >> buff2; if (buff2=="TRUE") POST_VTK = true; else POST_VTK = false;}

		else if (tag == "POST_STATS") 
			{ss >> buff2; if (buff2=="TRUE") POST_STATS = true; else POST_STATS = false;}

		else if (tag == "POST_REPORTS") 
			{ss >> buff2; if (buff2=="TRUE") POST_REPORTS = true; else POST_REPORTS = false;}

        else if (tag == "BMP_TRACKS") 
			{ss >> buff2; if (buff2=="TRUE") BMP_TRACKS = true; else BMP_TRACKS = false;}

        else if (tag == "TRACK_MIN_RANGE") 
			{ss >> TRACK_MIN_RANGE;} 

        else if (tag == "TRACK_MAX_RANGE") 
			{ss >> TRACK_MAX_RANGE;}

        else if (tag == "TRACKFRAMESTEP") 
			{ss >> TRACKFRAMESTEP;}

		else if (tag == "X_BMP")     
			{ss >> buff2; if (buff2=="TRUE") X_BMP = true; else X_BMP = false;}

		else if (tag == "Y_BMP") 
			{ss >> buff2; if (buff2=="TRUE") Y_BMP = true; else Y_BMP = false;}

		else if (tag == "Z_BMP") 
			{ss >> buff2; if (buff2=="TRUE") Z_BMP = true; else Z_BMP = false;}

		else if (tag == "CAGE_ON_SIDE") 
			{ss >> buff2; if (buff2=="TRUE") CAGE_ON_SIDE = true; else CAGE_ON_SIDE = false;}

        else if (tag == "BMP_FIX_BEAD_MOVEMENT") 
			{ss >> buff2; if (buff2=="TRUE") BMP_FIX_BEAD_MOVEMENT = true; else BMP_FIX_BEAD_MOVEMENT = false;}

        else if (tag == "BMP_FIX_BEAD_ROTATION") 
			{ss >> buff2; if (buff2=="TRUE") BMP_FIX_BEAD_ROTATION = true; else BMP_FIX_BEAD_ROTATION = false;} 

		else if (tag == "RESTORE_FROM_FRAME") 
			{ss >> RESTORE_FROM_FRAME;} 

		else if (tag == "NUCLEATOR_FORCES") 
			{ss >> NUCLEATOR_FORCES;} 

		else if (tag == "FORCE_SCALE_FACT") 
			{ss >> FORCE_SCALE_FACT;} 

		else if (tag == "FORCE_BAR_SCALE") 
			{ss >> FORCE_BAR_SCALE;} 

		else if (tag == "NUCLEATOR_INERTIA") 
			{ss >> NUCLEATOR_INERTIA;} 

        else if (tag == "VARY_INERT_W_RAD")
			{ss >> buff2; if(buff2=="TRUE") VARY_INERT_W_RAD = true; else VARY_INERT_W_RAD = false;}

		else if (tag == "XLINK_NODE_RANGE") 
			{ss >> XLINK_NODE_RANGE;} 

		else if (tag == "SEGMENT_BINS") 
			{ss >> buff2; if(buff2=="TRUE") SEGMENT_BINS = true; else SEGMENT_BINS = false;}

		else if (tag == "DRAW_CAGE")
			{ss >> buff2; if(buff2=="TRUE") DRAW_CAGE = true; else DRAW_CAGE = false;}

		else if (tag == "P_XLINK") 
			{ss >> P_XLINK;}

        else if (tag == "VARY_P_XLINK")
			{ss >> buff2; if(buff2=="TRUE") VARY_P_XLINK = true; else VARY_P_XLINK = false;}

		else if (tag == "NUC_LINK_FORCE")    
			{ss >> NUC_LINK_FORCE;}

		else if (tag == "NUC_LINK_BREAKAGE_FORCE") 
			{ss >> NUC_LINK_BREAKAGE_FORCE;}

        else if (tag == "NUC_LINK_BREAKAGE_DIST")   // note must be called after NUC_LINK_FORCE
			{ss >> tempdbl; NUC_LINK_BREAKAGE_FORCE = NUC_LINK_FORCE * tempdbl; }

		else if (tag == "LINK_BREAKAGE_FORCE") 
			{ss >> LINK_BREAKAGE_FORCE;}

		else if (tag == "LINK_BREAKAGE_STRAIN") 
			{ss >> LINK_BREAKAGE_STRAIN;}

		else if (tag == "LINK_POWER_SCALE") 
			{ss >> LINK_POWER_SCALE;}

		else if (tag == "LINK_FORCE") 
			{ss >> LINK_FORCE;}	

		else if (tag == "MAX_POLYMERISATION_PRESSURE") 
			{ss >> MAX_POLYMERISATION_PRESSURE;}
        
		else if (tag == "STICK_TO_NUCLEATOR") 
			{ss >> buff2; if(buff2=="TRUE") STICK_TO_NUCLEATOR = true; else STICK_TO_NUCLEATOR = false;}

		else if (tag == "RESTICK_TO_NUCLEATOR") 
			{ss >> buff2; if(buff2=="TRUE") RESTICK_TO_NUCLEATOR = true; else RESTICK_TO_NUCLEATOR = false;}
		
		else if (tag == "USETHREAD_COLLISION") 
			{ss >> buff2; if(buff2=="TRUE") USETHREAD_COLLISION = true; else USETHREAD_COLLISION = false;}
		
		else if (tag == "USETHREAD_LINKFORCES") 
			{ss >> buff2; if(buff2=="TRUE") USETHREAD_LINKFORCES = true; else USETHREAD_LINKFORCES = false;}
		
		else if (tag == "USETHREAD_APPLYFORCES") 
			{ss >> buff2; if(buff2=="TRUE") USETHREAD_APPLYFORCES = true; else USETHREAD_APPLYFORCES = false;} 
    	
		else if (tag == "USE_BREAKAGE_VISCOSITY") 
			{ss >> buff2; if(buff2=="TRUE") USE_BREAKAGE_VISCOSITY = true; else USE_BREAKAGE_VISCOSITY = false;} 

        else if (tag == "BREAKAGE_VISCOSITY_THRESHOLD") 
			{ss >> BREAKAGE_VISCOSITY_THRESHOLD;}
		
		else if (tag == "P_NUC")      
			{ss >> P_NUC;}

		else if (tag == "POLY_FEEDBACK") 
			{ss >> buff2; if(buff2=="TRUE") POLY_FEEDBACK = true; else POLY_FEEDBACK = false;} 
		
		else if (tag == "POLY_FEEDBACK_DIST") 
			{ss >> POLY_FEEDBACK_DIST;}

        else if (tag == "POLY_FEEDBACK_MIN_PROB") 
			{ss >> POLY_FEEDBACK_MIN_PROB;}

        else if (tag == "POLY_FEEDBACK_FACTOR") 
			{ss >> POLY_FEEDBACK_FACTOR;}

		else if (tag == "VISCOSITY")    
			{ss >> buff2; if(buff2=="TRUE") VISCOSITY = true; else VISCOSITY = false;} 
		
		else if (tag == "VISCOSITY_FACTOR") 
			{ss >> VISCOSITY_FACTOR;} 

		else if (tag == "VISCOSITY_EDGE_FACTOR") 
			{ss >> VISCOSITY_EDGE_FACTOR;}

		else if (tag == "VISC_DIST") 
			{ss >> VISC_DIST;}    

        else if (tag == "BROWNIANFORCESCALE") 
			{ss >> BROWNIANFORCESCALE;}

        else if (tag == "USE_BROWNIAN_FORCES") 
			{ss >> buff2; if(buff2=="TRUE") USE_BROWNIAN_FORCES = true; else USE_BROWNIAN_FORCES = false;}

		else if (tag == "NON_VISC_WEIGHTING") 
			{ss >> NON_VISC_WEIGHTING;}		 

		else if (tag == "MAX_VISC_WEIGHTING") 
			{ss >> MAX_VISC_WEIGHTING;}
        
		else if (tag == "IMPOSED_NUC_ROT_SPEED")  
			{ss >> IMPOSED_NUC_ROT_SPEED;} 

		else if (tag == "IMPOSED_NUC_ROT") 
			{ss >> buff2; if(buff2=="TRUE") IMPOSED_NUC_ROT = true; else IMPOSED_NUC_ROT = false;}

        else if (tag == "TEST_SQUASH") 
			{ss >> buff2; if(buff2=="TRUE") TEST_SQUASH = true; else TEST_SQUASH = false;}

        else if (tag == "TEST_FORCE_INITIAL_MAG")  
			{ss >> TEST_FORCE_INITIAL_MAG;}

        else if (tag == "TEST_FORCE_INCREMENT")  
			{ss >> TEST_FORCE_INCREMENT;}

        else if (tag == "TEST_DIST_EQUIL")  
			{ss >> TEST_DIST_EQUIL;}

		else if (tag == "WRITE_BMPS_PRE_SYMBREAK") 
			{ss >> buff2; if(buff2=="TRUE") WRITE_BMPS_PRE_SYMBREAK = true; else WRITE_BMPS_PRE_SYMBREAK = false;}

		else if (tag == "RADIUS") 
			{ss >> RADIUS;} 
		
		else if (tag == "CAPSULE_HALF_LINEAR") 
			{ss >> CAPSULE_HALF_LINEAR;} 
		
		else if (tag == "MAX_LINKS_PER_NEW_NODE") 
			{ss >> MAX_LINKS_PER_NEW_NODE;}

        else if (tag == "MAX_LINK_ATTEMPTS") 
			{ss >> MAX_LINK_ATTEMPTS;}
		
		else if (tag == "NODE_REPULSIVE_MAG") 
			{ss >> NODE_REPULSIVE_MAG;}
		
		else if (tag == "NODE_REPULSIVE_POWER") 
			{ss >> NODE_REPULSIVE_POWER;}
		
		else if (tag == "NODE_REPULSIVE_RANGE") 
			{ss >> NODE_REPULSIVE_RANGE;}
		
		else if (tag == "ASYMMETRIC_NUCLEATION") 
			{ss >> ASYMMETRIC_NUCLEATION;} 
		
		else if (tag == "RADIAL_SEGMENTS") 
			{ss >> RADIAL_SEGMENTS;} 
		
		else if (tag == "XLINK_NEAREST") 
			{ss >> XLINK_NEAREST;} 
		
		else if (tag == "VIEW_HEIGHT") 
			{ss >> VIEW_HEIGHT;} 
		
		//else if (tag == "NODES_TO_UPDATE") 
		//	{ss >> NODES_TO_UPDATE;} 
		
		else if (tag == "DISTANCE_TO_UPDATE") 
			{ss >> DISTANCE_TO_UPDATE;} 
		
		else if (tag == "GAUSSFWHM") 
			{ss >> GAUSSFWHM;}         
#ifdef BMP_USE_FOCAL_DEPTH
        else if (tag == "FOCALDEPTH") 
			{ss >> FOCALDEPTH;}         
#endif
        else if (tag == "BMP_INTENSITY_OFFSET") 
			{ss >> BMP_INTENSITY_OFFSET;}
		
		else if (tag == "SPECKLE") 
			{ss >> buff2;if(buff2=="TRUE") SPECKLE = true;else SPECKLE = false;}

        else if (tag == "COL_NODE_BY_STRAIN") 
			{ss >> buff2;if(buff2=="TRUE") COL_NODE_BY_STRAIN = true;else COL_NODE_BY_STRAIN = false;}

        else if (tag == "COL_LINK_BY_DIRN") 
			{ss >> buff2;if(buff2=="TRUE") COL_LINK_BY_DIRN = true;else COL_LINK_BY_DIRN = false;}

        else if (tag == "COL_INDIVIDUAL_NODES") 
			{ss >> buff2;if(buff2=="TRUE") COL_INDIVIDUAL_NODES = true;else COL_INDIVIDUAL_NODES = false;}

        else if (tag == "PLOTFORCES")         
			{ss >> buff2;if(buff2=="TRUE") PLOTFORCES = true;else PLOTFORCES = false;}
		
		else if (tag == "SPECKLE_FACTOR") 
			{ss >> SPECKLE_FACTOR;}

        else if (tag == "SPECKLEGRID") 
			{ss >> buff2;if(buff2=="TRUE") SPECKLEGRID = true;else SPECKLEGRID = false;}

    	else if (tag == "SPECKLEGRIDPERIOD") 
		    {ss >> SPECKLEGRIDPERIOD;}

    	else if (tag == "SPECKLEGRIDTIMEWIDTH") 
		    {ss >> SPECKLEGRIDTIMEWIDTH;}

        else if (tag == "SPECKLEGRIDSTRIPEWIDTH") 
		    {ss >> SPECKLEGRIDSTRIPEWIDTH;}

		else if (tag == "INIT_R_GAIN") 
			{ss >> INIT_R_GAIN;} 
		
		else if (tag == "INIT_G_GAIN") 
			{ss >> INIT_G_GAIN;} 
		
		else if (tag == "INIT_B_GAIN") 
			{ss >> INIT_B_GAIN;} 
		
		else if (tag == "ROTATION") 
			{ss >> buff2;if(buff2=="TRUE") ROTATION = true; else ROTATION = false;} 
		
		else if (tag == "MOFI") 
			{ss >> MOFI;} 
		
		else if (tag == "NO_IMAGE_TEXT") 
			{ss >> buff2;if(buff2=="TRUE") NO_IMAGE_TEXT = true; else NO_IMAGE_TEXT = false;} 
		
		else if (tag == "BMP_COMPRESSION") 
			{ss >> BMP_COMPRESSION;	if (BMP_COMPRESSION > 100) BMP_COMPRESSION = 100; else if (BMP_COMPRESSION < 0)	BMP_COMPRESSION = 0;} 
		
		else if (tag == "BMP_OUTPUT_FILETYPE") 
			{ss >> buff2; istringstream ss3(strtolower(buff2)); ss3 >> BMP_OUTPUT_FILETYPE;} 
		
		else if (tag == "SHAPE") 
			{ss >> buff2; if (buff2 == "CAPSULE") NUCSHAPE = nucleator::capsule; else NUCSHAPE = nucleator::sphere;}

        else if (tag == "CLUSTER") 
			{ss >> buff2;if(buff2=="TRUE") CLUSTER = true; else CLUSTER = false;}

		else if (tag.find("VIS") == 0)
		{
			// VTK stuff, ignore for now---perhaps put VTK stuff in here?
		}
		else
		{
			unrecognisedlines += buffer;
			unrecognisedlines += "\n";
			continue;
		}

		// got to here so line was recognised, add it to the unrecognisedlines string

		istringstream ss2(buffer);

		string tag2;

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
		
 #ifndef NOKBHIT
		if (!REWRITESYMBREAK && !POST_PROCESS && !QUIET && !CLUSTER)
		{        
			cout << "Ignore errors and start run anyway(y/n)?";
		    
			cout.flush();
			int ch;
			do 
			{
				ch = getch();
			} 
			while (( ch!='y' ) && ( ch!='Y') && ( ch!='n') && ( ch!='N'));

			if (( ch=='n' ) || ( ch=='N'))
			{
				cout << "n - Run aborted" << endl;
				cout.flush();
				system("stty sane");
				exit(1);
			}
			else
			{
				cout << "\r";
			}	

		}
#endif

	}

    if (CLUSTER)
    {   // override some parameters if running on a cluster machine

        QUIET = true;   // quiet still produces output, but only once per written file, rather than per second
        NO_IMAGE_TEXT = true;  //remove image text for matrix display
        SYM_BREAK_TO_RIGHT = true;  // offset and fix break direction to right for matrix display
    }


    // calculate commonly used constants from parameters:

    if (NUCSHAPE == nucleator::capsule)
    {
        // capsule goes all over the place, so make sure we're doing all 3 bitmaps
        X_BMP = Y_BMP = Z_BMP = true;
    }

    if (NOBITMAPS)
       X_BMP = Y_BMP = Z_BMP = false;

	double grid_scan_range;

#ifdef PROXIMITY_VISCOSITY 
	
	// if using proximity-based viscosity, we might need to increase the collision-detection scan dist

	if (VISCOSITY) 
		grid_scan_range =  mymax(NODE_REPULSIVE_RANGE,VISC_DIST);
	else
	    grid_scan_range =  NODE_REPULSIVE_RANGE;

#else

	    grid_scan_range =  NODE_REPULSIVE_RANGE;

#endif

	GRIDRES = 1.05 * grid_scan_range;	// make the grid res a bit bigger so that only scan one GP in either dirn
	
    gridscanjitter = GRIDRES - grid_scan_range; // how far a node can move before we update its grid posn

    GRIDSIZE = (int) ( 2 * GRIDBOUNDS / GRIDRES );

	NODE_XLINK_GRIDSEARCH = (int) floor( XLINK_NODE_RANGE / GRIDRES ) + 1; 

#ifdef PROXIMITY_VISCOSITY 
	
	// if using proximity-based viscosity, we might need to increase the collision-detection scan dist
	NODE_REPULSIVE_RANGE_GRIDSEARCH = (int) floor( grid_scan_range  / GRIDRES ) + 1;

#else

	NODE_REPULSIVE_RANGE_GRIDSEARCH = (int) floor(((double) NODE_REPULSIVE_RANGE )/GRIDRES) + 1;

#endif

	TOTAL_ITERATIONS = (int) (((double)TOTAL_SIMULATION_TIME / (double)DELTA_T)+0.5);
    RECORDING_INTERVAL = int(TOTAL_ITERATIONS / TOT_FRAMES);
    TOTAL_ITERATIONS = RECORDING_INTERVAL * TOT_FRAMES; // just to make sure it's an exact muliple
	
	InterRecordIterations = RECORDING_INTERVAL;  // redundant---combine these two variables at some point
    
	//TOT_FRAMES = int(TOTAL_ITERATIONS / RECORDING_INTERVAL);

	NODE_FORCE_TO_DIST = DELTA_T * FORCE_SCALE_FACT;
	NODE_DIST_TO_FORCE = 1.0 / NODE_FORCE_TO_DIST;

    BMP_WIDTH  *= BMP_AA_FACTOR;
    BMP_HEIGHT *= BMP_AA_FACTOR;

	if (SPECKLE_FACTOR<0)
	{
        SPECKLE_FACTOR = 1;
		SPECKLE = false;
		cout << "Negative SPECKLE_FACTOR reset" << endl;
	}

	if (SPECKLE)
	{	// is speckle is on, scale the gain
		INIT_R_GAIN *= 1.0 / SPECKLE_FACTOR;
		INIT_G_GAIN *= 1.0 / SPECKLE_FACTOR;
		INIT_B_GAIN *= 1.0 / SPECKLE_FACTOR;
	}
	

    vector <double> distmoved(TOT_FRAMES+1,0.0);  // keep track of the distance moved and num nodes
    vector <int> numnodes(TOT_FRAMES+1,0);      // for the NODES_TO_UPDATE
    vector <vect> velocities(TOT_FRAMES+1);
    vector <vect> posn(TOT_FRAMES+1);

    int lastframedone = 0;

    if ((RESTORE_FROM_FRAME != 0) || (POST_PROCESS))
    {   // if restoring or post-processing, 
        // restore the NODESUPDATE list
        // post-processing uses this to set lastframedone
        // in case we're processing all the frames done

        ifstream ipnodesupdate(NODESUPDATEFILE);

        if (!ipnodesupdate) 
            { cout << "Unable to open file " << NODESUPDATEFILE << " for input"; }
        else
        {
            for (int fn = 1; fn != TOT_FRAMES+1; ++fn)
            {
                ipnodesupdate 
	                >> distmoved[fn]
	                >> numnodes[fn]
                    >> posn[fn]
                    >> velocities[fn] ;

                if (numnodes[fn]>0)
                    lastframedone = fn;
            }
            ipnodesupdate.close();
        }
  
    }

    vector<int> postprocess_iterations;
	postprocess_iterations.resize(0);

    if (POST_PROCESS)
    {
        cout << "Postprocessing iterations: ";
	    get_postprocess_iterations(argv[2], postprocess_iterations, lastframedone);
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

    if (RESTORE_FROM_FRAME == 0 && !REWRITESYMBREAK && !POST_PROCESS)
    {  // new run, clear directories

#ifndef NOKBHIT
		if (!QUIET)
		{        
			cout << "Clear directories and start run(y/n)?";
		    
			cout.flush();
			int ch;
			do 
			{
				ch = getch();
			} 
			while (( ch!='y' ) && ( ch!='Y') && ( ch!='n') && ( ch!='N'));

			if (( ch=='n' ) || ( ch=='N'))
			{
				cout << "n - Run aborted" << endl;
				cout.flush();
				system("stty sane");
				exit(1);
			}
			else
			{
				cout << endl;
			}	

		}
#endif  
        // only if starting a new calculation from scratch, clear the directories
        cout << "Deleting old save files...";
		cout.flush();

#ifndef USEWINDOWSCOMMANDS

		sprintf(command1, "rm -f %s*_0*.%s 2>/dev/null", BITMAPDIR, BMP_OUTPUT_FILETYPE.c_str() );
		system(command1);
		sprintf(command1, "rm -f %s*.wrz 2>/dev/null", VRMLDIR );
		system(command1);
		sprintf(command1, "rm -f %s*%s 2>/dev/null", REPORTDIR, COMPRESSEDEXTENSION);
		system(command1);
		sprintf(command1, "rm -f %s*%s 2>/dev/null", DATADIR, COMPRESSEDEXTENSION );
		system(command1);
		sprintf(command1, "rm -f %s*.wrl %s*.txt 2>/dev/null", TEMPDIR, TEMPDIR );
		system(command1);
        sprintf(command1, "rm -f %s 2>/dev/null", SYM_BREAK_FILE );
		system(command1);
		sprintf(command1, "rm -f %s*.* 2>/dev/null", VTKDIR );
		system(command1);        
		sprintf(command1, "rm -f %s*.* 2>/dev/null", STATSDIR );
		system(command1);        
#else
		sprintf(command1, "del /q %s*_0*.%s", BITMAPDIR, BMP_OUTPUT_FILETYPE.c_str() );
		system(command1);
		sprintf(command1, "del /q %s*.wrz", VRMLDIR );
		system(command1);
		sprintf(command1, "del /q %s*%s", REPORTDIR, COMPRESSEDEXTENSION );
		system(command1);
		sprintf(command1, "del /q %s*%s", DATADIR, COMPRESSEDEXTENSION );
		system(command1);
		sprintf(command1, "del /q %s*.wrl %s*.txt", TEMPDIR, TEMPDIR );
		system(command1);
		sprintf(command1, "del /q %s", SYM_BREAK_FILE );
		system(command1);
		sprintf(command1, "del /q %s*.*", VTKDIR );
		system(command1);
		sprintf(command1, "del /q %s*.*", STATSDIR );
		system(command1);
#endif

        cout << "done." << endl;
    }


	// create main objects
	// create as static otherwise exit() doesn't call their destructors (?!)

	//static consts CONST;

	static actin theactin;

    ptheactin = &theactin;  // ugly global pointer for access from nodes,nucleator and segments

    static nucleator nuc_object(NUCSHAPE);//, &theactin);

	if (REWRITESYMBREAK)
	{
		rewrite_symbreak_bitmaps(nuc_object, theactin);   
		exit(EXIT_SUCCESS);
	}

	// write out parameters to screen and file

	cout << "Total simulation time:      " << TOTAL_SIMULATION_TIME << endl;
	cout << "Delta_t:                    " << DELTA_T << endl;

	theactin.opruninfo.flush();

	cout << "Total iterations: " << TOTAL_ITERATIONS << endl;
	cout << "Frames every " << InterRecordIterations  
		<< " iterations (" << TOT_FRAMES << " total)" << endl;

	

	unsigned int starttime, endtime, lasttime ,nowtime;//, lastitertime;

    starttime = (unsigned) time(NULL);

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

	vect last_center, center, delta_center;
        
    



	// Breakout if we are post processing
	if( !postprocess_iterations.empty() )
	{
		postprocess(nuc_object, theactin, postprocess_iterations,  argv );
	    exit(EXIT_SUCCESS); // FIXME: early finish, move to two fcns (ML)
	}

	// - - - - - - - - - - 
	// main iteration loop
	// - - - - - - - - - -
	// initialse from a checkpoint if requested
	int starting_iter = 1;



	if (RESTORE_FROM_FRAME != 0)
	{
		cout << "Loading data...";
		cout.flush();

        starting_iter = RESTORE_FROM_FRAME * InterRecordIterations;

	    load_data(theactin, starting_iter, true);

        cout << "Restored from frame "
			 << RESTORE_FROM_FRAME << " (iteration " << starting_iter << ")" << endl;

        // srand( (unsigned) 200 );
 
		//	starting_iter = RESTORE_FROM_FRAME + 1; // don't overwrite

        if (theactin.highestnodecount > ((int)theactin.node.size() - 1000))
            theactin.reservemorenodes(10000);



        if (TEST_SQUASH)
        {
            theactin.testforces_setup();
        }
	}


	unsigned int rand_num_seed;

	if (RAND_SEED != -1) 
	{
		cerr << "Warning: Static random number seed used" <<  endl;
		theactin.opruninfo << "Warning: Static random number seed used" <<  endl;

		rand_num_seed = RAND_SEED;
	} 
	else
	{			  
		rand_num_seed = (unsigned int)( time(NULL) * getpid());
		theactin.opruninfo << "Time and PID based random number seed (" << rand_num_seed << ") used" << endl;
		cerr << "Time and PID based random number seed (" << rand_num_seed << ") used" << endl;
    }

	srand(rand_num_seed);

#ifdef USE_MERSENNE

    mers_rand.resize(NUM_THREAD_DATA_CHUNKS);

    for (int i=0; i !=  NUM_THREAD_DATA_CHUNKS; ++i)
    {
        mers_rand[NUM_THREAD_DATA_CHUNKS].seed( rand_num_seed + NUM_THREAD_DATA_CHUNKS );
    }

#endif

	int filenum = 0;
	double polrate = 0;

    char last_symbreak_bmp_filename[1024] = "";

    CometVtkVis vtkvis(VIEW_VTK);   // create vtk object

	cout << endl;
	cout << "Starting iterations..." << endl << endl; 

	if (!QUIET)
	{
		cout << "(Press 'q' to abort run, or 'm' to make quicktime movies of frames so far)" << endl << endl;
	}
	
	if (QUIET)
	{
		cout << "Running in quiet mode (no continual progress display)" << endl
			<< "Use kill " << getpid() << " to terminate"<< endl;
	}
	else
	{
		cout << "Itternum|TotNode|Pol  |Links+|Links-|dist   |Rotn   |T     |SaveNum" << endl << endl;
	}

    mytimer frametimer;
    frametimer.setprecision(1);
    frametimer.start();


    // this info file is opened here, because there's one per frame
	char infofilename[1024];
	sprintf ( infofilename , "%sinfo%05i.txt",REPORTDIR, 0 );
	theactin.opinfo.open(infofilename, ios::out | ios::trunc);
	if (!theactin.opinfo) 
	{ cout << "Unable to open file " << infofilename << " for output";}

    


	for(int i=starting_iter; i<=TOTAL_ITERATIONS; i++)
	{

		filenum = (int)(i/InterRecordIterations);

		/*if ((i % 10) == 1)
			srand( rand_num_seed );*/

		nowtime = (unsigned) time(NULL);

		if ((nowtime > lasttime) || ((i % InterRecordIterations) == 1))
		{

            if (theactin.highestnodecount > ((int)theactin.node.size() - 1000))
                theactin.reservemorenodes(10000);
			
			if (theactin.attemptedpolrate >0)
				polrate = (double) theactin.polrate/ (double) theactin.attemptedpolrate;// / (double) (i % InterRecordIterations);
            else
                polrate = 0;

			if (!QUIET)
			{
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
                                ABORT = true;
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
                            endtime = (unsigned) time(NULL);

                            cout << endl << "Time so far: " 
		                        << ((endtime-starttime) / 3600) << "h " 
		                        << ((endtime-starttime) / 60) % 60 << "m " 
		                        << ((endtime-starttime) % 60) << "s " << endl;

							cout << endl << endl << "Directory:";
							cout.flush();

	#ifndef USEWINDOWSCOMMANDS
							system("pwd");
	#else
							system("chdir");
	#endif
							cout << endl << "Making movie of frames so far..." << endl;
							cout.flush();
							cout << endl;
					        

							if (X_BMP)
                            {
                                sprintf(command1, "~/bin/ffmpeg -v 1 -b 1800000  -y -i bitmaps//x_proj_%%05d.%s -vcodec mpeg4 x_proj.mov", BMP_OUTPUT_FILETYPE.c_str());
                                system(command1);
                            }
							//if (Y_BMP)
							//	system("~/bin/ffmpeg -v 1 -b 1800000  -y -i bitmaps//y_proj_%05d.jpeg -vcodec mpeg4 y_proj.mov");
							//if (Z_BMP)
							//	system("~/bin/ffmpeg -v 1 -b 1800000  -y -i bitmaps//z_proj_%05d.jpeg -vcodec mpeg4 z_proj.mov");
							
							cout << endl;
						}
					}
				}
#endif
			}           

			lasttime = nowtime;
			theactin.find_center(center);

            if (!TEST_SQUASH)
			    distfromorigin = center.length();
            else
                distfromorigin = theactin.testsurfaceposn;

            theactin.nuc_to_world_rot.getangles(x_angle,y_angle,z_angle);
            tot_rot = fabs(x_angle) + fabs(y_angle) + fabs(z_angle);
		


			if (!QUIET)
			{   // once per second screen output
				cout<< "I"   << setw(7) << i 
				    << "|N"  << setw(6) << theactin.highestnodecount
				    << "|P"  << setw(4) << setprecision(2) << polrate
				    << "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
				    << "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
				    << "|d"  << setw(6) << setprecision(3) << distfromorigin
				    << "|R"  << setw(6) << setprecision(1) << (180/PI) * tot_rot
                    << "|T"  << setw(5) << frametimer << "\r"; //((unsigned) time(NULL) - lastitertime) << "\r";
				cout.flush();
			}			
		}

		theactin.iterate();
		


		if (((i % InterRecordIterations) == 0) && (i>starting_iter))
		{   // we've reached a frame save point
            
            if (theactin.attemptedpolrate >0)
				    polrate = (double) theactin.polrate/ (double) theactin.attemptedpolrate;// / (double) (i % InterRecordIterations);
                else
                    polrate = 0;

#ifdef NON_RANDOM

srand( rand_num_seed );

#endif

			theactin.find_center(center);

            posn[filenum] = center;
            distmoved[filenum] = center.length();
            numnodes[filenum] = theactin.highestnodecount;

			if (!TEST_SQUASH)
			    distfromorigin = distmoved[filenum];
            else
                distfromorigin = theactin.testsurfaceposn;

            theactin.nuc_to_world_rot.getangles(x_angle,y_angle,z_angle);
            tot_rot = fabs(x_angle) + fabs(y_angle) + fabs(z_angle);
			
			delta_center = center - last_center;  
			last_center = center;

            // velocity is how far moved divided by the frame time (i.e. DELTA_T * InterRecordIterations) 
            velocities[filenum] = delta_center * ( 1 / (DELTA_T * InterRecordIterations));

            //cout << "test nodes: "<< theactin.testnodes.size() << endl;

			theactin.opinfo.close();
			
            sprintf ( infofilename , "%sinfo%05i.txt",REPORTDIR, filenum );
			theactin.opinfo.open(infofilename, ios::out | ios::trunc);
			
            if (!theactin.opinfo) 
			{ cout << "Unable to open file " << infofilename << " for output"; }

			theactin.opvelocityinfo 
				<< (i*DELTA_T) << "," 
				<< center.x << "," 
				<< center.y << "," 
				<< center.z << "," 
				<< velocities[filenum].length() << endl;

            theactin.opvelocityinfo.flush();

            // once per save screen and log output

			if (!QUIET)
			{
				cout << "I" << setw(7) << i;
			}
			else
			{
				cout << "PID:" << setw(6) << getpid();    // write out pid if quiet so we know which proc to kill
			}

            theactin.opruninfo 
                    << "|N"  << setw(6) << theactin.highestnodecount
				    << "|P"  << setw(4) << setprecision(2) << polrate
				    << "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
				    << "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
				    << "|d"  << setw(6) << setprecision(3) << distfromorigin
				    << "|R"  << setw(6) << setprecision(1) << (180/PI) * tot_rot	
				    << "|T"  << setw(5) << frametimer 
				    << "|S " << setw(3) << (int)filenum  
				    << "/"   << TOT_FRAMES;

            cout    << "|N"  << setw(6) << theactin.highestnodecount
				    << "|P"  << setw(4) << setprecision(2) << polrate
				    << "|L+" << setw(4) << (theactin.linksformed-lastlinksformed)/2 
				    << "|L-" << setw(4) << (theactin.linksbroken-lastlinksbroken)/2
				    << "|d"  << setw(6) << setprecision(3) << distfromorigin
				    << "|R"  << setw(6) << setprecision(1) << (180/PI) * tot_rot	
				    << "|T"  << setw(5) << frametimer 
				    << "|S " << setw(3) << (int)filenum  
				    << "/"   << TOT_FRAMES;

			if ( !WRITE_BMPS_PRE_SYMBREAK && 
				 !finished_writing_sym_bitmaps)
				cout << "*";
            
            theactin.opruninfo << endl;
            theactin.opruninfo.flush();
            cout << endl;
			cout.flush();

			// we don't use vrml anymore, so don't bother with it for now
			// theactin.savevrml(filenum);

			if (WRITE_BMPS_PRE_SYMBREAK)
			{
				 if (!theactin.brokensymmetry)
				 {
					nuc_object.segs.addallnodes();  // put node data into segment bins
					nuc_object.segs.set_scale_factors();

					theactin.savebmp(filenum, xaxis, actin::runbg, true);
					theactin.savebmp(filenum, yaxis, actin::runbg, true);
					theactin.savebmp(filenum, zaxis, actin::runfg, true);

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
					theactin.savebmp(filenum, xaxis, actin::runfg, false);
					theactin.savebmp(filenum, yaxis, actin::runfg, false);
					theactin.savebmp(filenum, zaxis, actin::runfg, false);

					cout << "\r";
					cout.flush();
		        
				 }

			}

			if ((!TEST_SQUASH) && (!theactin.brokensymmetry) && (distfromorigin > RADIUS))
			{	// symmetry broke: set directions, scale factors etc.
				
				theactin.brokensymmetry = true;

				// set camera rotation for save  
				theactin.set_sym_break_axes();
				theactin.save_sym_break_axes();

                // calculate but don't write bitmaps to get scaling factors
			    theactin.savebmp(filenum, xaxis, actin::runfg, false);
			    theactin.savebmp(filenum, yaxis, actin::runfg, false);
			    theactin.savebmp(filenum, zaxis, actin::runfg, false);

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
				if (Z_BMP)
				{
					sprintf(last_symbreak_bmp_filename, "%sz_proj_%05i.%s",BITMAPDIR, 
							1 ,BMP_OUTPUT_FILETYPE.c_str());
				}
				else 
				{
					if (Y_BMP)
					{
						sprintf(last_symbreak_bmp_filename, "%sy_proj_%05i.%s",BITMAPDIR, 
								1 ,BMP_OUTPUT_FILETYPE.c_str());
					} else
					{
						sprintf(last_symbreak_bmp_filename, "%sx_proj_%05i.%s",BITMAPDIR, 
								1 ,BMP_OUTPUT_FILETYPE.c_str());
					}
				}
				

				if (WRITE_BMPS_PRE_SYMBREAK)
				{
					// if we have already written draft bitmaps, need to delete the last one
					// so we know when it's re-written
					sprintf(command1, "rm %s 1>/dev/null 2>/dev/null", last_symbreak_bmp_filename);
					system(command1);
				}

				// call self with 'sym' argument to re-write the bitmaps
				sprintf(command1, "%s sym 1>rewritesymbreak.txt 2>rewritesymbreak.err &", argv[0]);
				system(command1);

			}  

            save_data(theactin, i);

            

            // save the nodes per dist file for the nodesupdate
			ofstream opnodesupdate(NODESUPDATEFILE, ios::out | ios::trunc);
			
            if (!opnodesupdate) 
			    { cout << "Unable to open file " << NODESUPDATEFILE << " for output"; }
            else
            {
                // save the nodes update file
                for (int frame = 1; frame != TOT_FRAMES+1; ++frame)
                {
                    opnodesupdate 
				        << distmoved[frame]  << " " 
				        << numnodes[frame]   << " "
                        << posn[frame]       << " "
                        << velocities[frame] << endl;
                }

                opnodesupdate.close();
            }

            // find the symmetry break time



            // find max velocity

            double maxvel = 0;

            for (int frame = 1; frame != TOT_FRAMES+1; ++frame)
            {
                if (velocities[frame].length() > maxvel)
                    maxvel = velocities[frame].length();
            }

            
            // find where crosses half max

            int framemaxvelmoved = 0;                
            vect maxvelmoved;

            for (int frame = 1; frame != TOT_FRAMES+1; ++frame)
            {
                if (velocities[frame].length() > maxvel / 2.0) 
                {
                    maxvelmoved = velocities[frame];
                    framemaxvelmoved = frame;
                    break;
                }
            }

            //int SYMBREAKSEARCHWINDOWHALFWIDTH = 3; // look over double this width to find max speed
            //for (int frame = 1; frame != TOT_FRAMES+1; ++frame)
            //{
            //    if ((frame > SYMBREAKSEARCHWINDOWHALFWIDTH * 2) && (frame <= filenum))
            //    {   // we need length separate from the vectors because the mag is not just sum of delta lengths of the vectors,
            //        // but the sum of the lengths i.e. it's negative if it's slowing down
            //        vect deltav;
            //        deltav.zero();
            //        double deltavmag = 0.0;

            //        for (int fr = frame - SYMBREAKSEARCHWINDOWHALFWIDTH * 2; fr <= frame; ++fr)
            //        {  // find delta v over range
            //            deltavmag += velocities[fr].length() - velocities[fr-1].length();
            //            deltav    += velocities[fr]          - velocities[fr-1];
            //        }

            //        deltavmag *= 1.0/(SYMBREAKSEARCHWINDOWHALFWIDTH * 2); // rescale by window width
            //        deltav    *= 1.0/(SYMBREAKSEARCHWINDOWHALFWIDTH * 2);

            //        //vect vel = (posn[frame] - posn[frame - SYMBREAKSEARCHWINDOWHALFWIDTH * 2]) * ( (double) SYMBREAKSEARCHWINDOWHALFWIDTH / (DELTA_T * InterRecordIterations)) ;
            //        
            //        //cout << "Frame " << frame << setprecision(4) << " deltav: " << deltavmag << " " << endl ;
            //        if (deltavmag > maxvelmovedmag)
            //        {     
            //            maxvelmovedmag = deltavmag;
            //            maxvelmoved = deltav;
            //            framemaxvelmoved = frame - SYMBREAKSEARCHWINDOWHALFWIDTH; 
            //            //cout << setprecision(4) << "Max vel set to: " << maxvelmoved ;
            //            //cout << " At frame: "<< frame << endl;
            //        }
            //        
            //    }                 
            //}

            // write symmetry break time (and direction) to file

            ofstream opsymbreaktime("symbreaktime.txt", ios::out | ios::trunc);
            if (!opsymbreaktime) 
		    { 
                cout << "Unable to open file symbreaktime.txt for output"; 
            }
            else
            {
                opsymbreaktime << "Frame Time Speed Velocity" << endl;
                opsymbreaktime << framemaxvelmoved << " " << framemaxvelmoved * DELTA_T * InterRecordIterations << " " << " " << maxvelmoved << endl;
                opsymbreaktime.close();
            }

   
			
            if ( (DISTANCE_TO_UPDATE > 0.01) &&
				 (distmoved[filenum] > DISTANCE_TO_UPDATE*RADIUS) )
			{

                NODES_TO_UPDATE = theactin.highestnodecount;

                for (int fn = filenum; fn != 0; --fn)
                {
                    if ((distmoved[fn] < (distmoved[filenum] - DISTANCE_TO_UPDATE*RADIUS)) &&
                        (numnodes[fn] != 0))
                    {
                        NODES_TO_UPDATE = theactin.highestnodecount - numnodes[fn];
                        break;
                    }
                }

                if (!DISTANCE_TO_UPDATE_reached)
                {
				    DISTANCE_TO_UPDATE_reached = true;
				    cout << endl << "NODES_TO_UPDATE set to " << NODES_TO_UPDATE << endl;
                }
			}

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

				theactin.savebmp(filenum, xaxis, actin::runbg, true);
				theactin.savebmp(filenum, yaxis, actin::runbg, true);
				theactin.savebmp(filenum, zaxis, actin::runbg, true);

				cout << "\r";
				cout.flush();

				if (*last_symbreak_bmp_filename != 0)
				{
                    //check if last file exists yet
					ifstream ip_last_symbreak_bmp_file( last_symbreak_bmp_filename );

					if (ip_last_symbreak_bmp_file)
					{
						cout << "Finished background processing of symmetry breaking bitmaps               " << endl;

						*last_symbreak_bmp_filename = 0 ;

						ip_last_symbreak_bmp_file.close();

                        finished_writing_sym_bitmaps = true;
					}
				}

			}

            if (TEST_SQUASH && theactin.testsurfaceposn < 0)
            {
                cout << endl << endl << "Test position reached nucleator surface---finished." << endl;
                ABORT = true;
                break;
            }

            if (VIEW_VTK)
                vtkvis.buildVTK(filenum);

            theactin.setdontupdates();

			theactin.clear_node_stats();  // clear the cumulative stats data in the nodes *after setdontupdates*
                                          // so we preserve the numbers in the ones we're no longer calculating

			theactin.compressfilesdowork(filenum);

			lastlinksformed = theactin.linksformed;
			lastlinksbroken = theactin.linksbroken;

			theactin.opruninfo.flush();

            frametimer.start();

		}
	}

	cout << endl << "Done " << endl << endl;

	theactin.saveinfo();

	endtime = (unsigned) time(NULL);

    theactin.opruninfo 
        << endl << "Time : " 
		<< ((endtime-starttime) / 3600) << "h " 
		<< ((endtime-starttime) / 60) % 60 << "m " 
		<< ((endtime-starttime) % 60) << "s " << endl;
        
    cout<< endl << "Time : " 
		<< setprecision(1) << ((endtime-starttime) / 3600) << "h " 
		<< ((endtime-starttime) / 60) % 60 << "m " 
		<< ((endtime-starttime) % 60) << "s " << endl;
 
    // allow imagemagick to catch up before calling actin destructor and 
    // deleting the temp bitmap files
    if (!ABORT)
    {
        cout << "Pausing for background jobs to complete" << endl;
#ifdef _WIN32
        Sleep(10000);
#else
        usleep(10000);
#endif
    }
          
	
	exit(EXIT_SUCCESS);
}


string get_datafilename(const int iteration)
{
	//stringstream filenamestr;
    char filename[1024];
	sprintf(filename , "data_%07i.txt", iteration);
    return filename;
    
}

bool load_data(actin &theactin, int iteration, const bool &loadscale)
{

    string filename = get_datafilename(iteration);

	char tmpdatafile[1024];

	sprintf(tmpdatafile, "%stempdata_%u_%u.txt", TEMPDIR,
		(unsigned int) getpid(), (unsigned int) time(NULL) );

	filename = DATADIR + filename;

    // try to open uncompressed file first, then compressed file
    // (when working repeatedly on same dataset, prolly best to unzip the whole data directory)

    ifstream ifstrm;
    ifstrm.open( filename.c_str() );

    bool usingtmpdatafile = false;
    
    if (!ifstrm)
    {
        // uncompressed file open failed, so try decompressing the compressed file to a temp file and opening that
        ifstrm.close();

#ifndef _WIN32
	    char command1[1024];
        sprintf(command1, "%s -c %s%s > %s", DECOMPRESSCOMMAND, filename.c_str(), COMPRESSEDEXTENSION, tmpdatafile);
        system(command1);

        ifstrm.open( tmpdatafile );
        usingtmpdatafile = true;

#endif
    }
                                                    
    if(!ifstrm) {
	cout << "Unable to open file " << filename.c_str() << " or " << tmpdatafile << " for input" << endl;
    #ifdef _WIN32
    cout << "Note that on windows, you must manually unzip any gzipped data files" << endl;
    #endif
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

    ifstrm  >> saved_iteration
            >> theactin.brokensymmetry
            >> finished_writing_sym_bitmaps
            >> DISTANCE_TO_UPDATE_reached
            >> NODES_TO_UPDATE;

    if (!theactin.load_data(ifstrm))
    {
        cout << "Load data failed" << endl;
        cout.flush();
        abort();
    }
    theactin.setdontupdates();
    
    ifstrm.close();

 
    if (usingtmpdatafile)
    {
#ifndef _WIN32
        char command1[1024];
        sprintf(command1, "rm  %s",tmpdatafile);
        system(command1);
#endif    
    }


    // check the iteration is correct
    if( saved_iteration != iteration ){
	cout << "error in saved file, saved iteration." << endl;
	return 1;
    }

    // try to load sym break stuff if present

    if (theactin.load_sym_break_axes())
        theactin.brokensymmetry = true;  // set this true if able to load sym break file
                                         // since may be restoring to a pre sym break time
                                         // and we would still want same direction etc.

    if (loadscale)
    {
        // if we're able to load the scale factors
        // turn off the auto-scaling

        if (theactin.p_nuc->segs.load_scalefactors())
            theactin.BMP_intensity_scaling = false;
    }



    return true;
}

int save_data(actin &theactin, int iteration)
{
    string filename = get_datafilename(iteration);

    if (COMPRESSDATAFILES) // if compressing, put in temp directory, else save directly
	    filename = TEMPDIR + filename;
    else
        filename = DATADIR + filename;

    ofstream ofstrm( filename.c_str() );    
    if(!ofstrm) {
	cout << "Unable to open file " << filename << " for output" << endl;
	return 1;
    } 
    
    // write out a header  and save iteration
    ofstrm << "comet:" << endl
	   << iteration << " " 
       << theactin.brokensymmetry << " "
       << finished_writing_sym_bitmaps << " "
       << DISTANCE_TO_UPDATE_reached << " "
       << NODES_TO_UPDATE << endl;
    

    // what's the best way to set the data format (so keeps enough accuracy, but compresses well)
	ofstrm << setprecision(SAVE_DATA_PRECISION);
    //ofstrm.setf(ios_base::fixed);

    // actin does all the real work
    theactin.save_data(ofstrm);
    
    ofstrm.close();
    
    return 0;
}

void get_postprocess_iterations(const char *iterdesc, vector<int> &postprocess_iterations, const int& lastframedone)
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
	  cerr << "ERROR: too few values in range request (2 min):" 
	       << iterdesc << endl;
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
	  cerr << "ERROR: too many values in range request (3 max):" 
	       << iterdesc << endl;
	  exit(EXIT_FAILURE);
	}

    if (start == 0)
        start = InterRecordIterations;

    if (step == 0)
        step = InterRecordIterations;

    if (end == 0) // if user set last frame to 0, use the last one found from the nodesupdate file
    {
        end = lastframedone * InterRecordIterations;
    }

	postprocess_iterations.resize(0);
	//postprocess_iterations.clear();
	for(int i= start; i<=end; i+=step){
	    postprocess_iterations.push_back(i);
	}
    } 
}

// STUB function, doesn't do anything, but fill implementation out here
// 
void post_stats(actin &, const int filenum)
{
  cout << "  saving statistics for " << filenum << endl;  

  char nstats_filename[1024];
  sprintf(nstats_filename , "%snodestats_%05i.txt", STATSDIR, filenum);
  //ofstream nout;
  //nout.open(nstats_filename);

  char lstats_filename[1024];
  sprintf(lstats_filename , "%slinkstats_%05i.txt", STATSDIR, filenum);
  //ofstream lout;
  //lout.open(lstats_filename);

  // STUB for the function to get out statistics of current interest
  // fill in calculations here
}

// We may want to have a vtk postprocess as a seperate option?
// otherwise I've lumped it with the bitmap processing here
// the frame counting is seperated, but consider unifying 
// this once we are happy with how it all works.
void postprocess(nucleator& nuc_object, actin &theactin, 
		 vector<int> &postprocess_iterations, char* argv[])
{
	if (POST_PROCESS4CPU)
	{
		int threads = 4;
		int firstframe =  *postprocess_iterations.begin()  / InterRecordIterations;
		int lastframe  = *(postprocess_iterations.end()-1) / InterRecordIterations;

		//lastframe = lastframe - (lastframe % 4);  // make # frames divisible by 4 for quicktime compat

		char command1[1024];

		for (int i = 1; i != threads; ++i)
		{
			sprintf(command1, "%s post %i:%i:%i &", argv[0], firstframe + i, threads, lastframe);
			system(command1);
			//cout << command1 << endl;
		}

        // first one not background:

        sprintf(command1, "%s post %i:%i:%i", argv[0], firstframe , threads, lastframe);
		system(command1);

	} else
    {

		int filenum;

		// vtk
		
	    
		cout << "Post processing " 
		 << postprocess_iterations.size()
		 <<  " data sets." << endl;
	    
        const int maxframes = (int) postprocess_iterations.size();
		int frame = 0;

        theactin.node_tracks[xaxis].resize(0);
        theactin.node_tracks[yaxis].resize(0);
        theactin.node_tracks[zaxis].resize(0);

        // load the *last* image to set bmp calib
        int framemaxvelmoved;
        vect maxvelmoved;
        
        int iter = *(postprocess_iterations.end()-1);  // default to last frame if we can't find sym break time

        ifstream ipsymbreaktime("symbreaktime.txt");
        if (!ipsymbreaktime) 
	    { 
            cout << "Unable to open file symbreaktime.txt for input";
        }
        else
        {
            string buff;
            
            ipsymbreaktime >> buff >> buff >> buff >> buff;
            ipsymbreaktime >> framemaxvelmoved >> buff >> maxvelmoved;
            ipsymbreaktime.close();

            iter =  framemaxvelmoved * InterRecordIterations;
        }

        


        if (POST_BMP) // if we're doing the bitmaps, need to scale the intensities and segments
        {

            cout << "Loading frame " << (iter/InterRecordIterations) << " for scaling" << endl;

            load_data(theactin, iter, false);

            theactin.BMP_intensity_scaling = true;

            theactin.savebmp(filenum, xaxis, actin::runfg, false); 
		    theactin.savebmp(filenum, yaxis, actin::runfg, false); 
		    theactin.savebmp(filenum, zaxis, actin::runfg, false);

            cout << endl;

            nuc_object.segs.addallnodes();                                  
            nuc_object.segs.set_scale_factors();

        }

        if (BMP_TRACKS)
        {
            // load the first tracking frame
            
            iter = InterRecordIterations * (int)ceil(TRACK_MAX_RANGE);

            cout << "Loading frame " << (int)ceil(TRACK_MAX_RANGE) << " for track node selection" << endl;

            load_data(theactin, iter, false);

            theactin.set_nodes_to_track(xaxis); // only track x-axis for now

            theactin.savebmp(filenum, xaxis, actin::runfg, false);  // sets the initial offset position

            POST_PROC_ORDER = 1;  // must go forwards for tracks
        }

        // we can't create the vkt object here because we can't reuse the renderer without
        // causing segfaults for some reason
        // workaround is to create the vtk object within the loop, so that
        // the destructor clears out and we recreate the object every time
        //CometVtkVis vtkvis(true); // should this be false? (i.e. no render window)


        vector<int>::iterator start = postprocess_iterations.begin();
        vector<int>::iterator end = postprocess_iterations.end();

        if (POST_PROC_ORDER == -1)
        {
            start = postprocess_iterations.end();
            end = postprocess_iterations.begin();
        }

        for(vector<int>::iterator iteration  = start; 
		                          iteration != end;
                                  iteration += POST_PROC_ORDER)
		//for(vector<int>::reverse_iterator iteration  = postprocess_iterations.rbegin(); 
		//                                  iteration != postprocess_iterations.rend();
        //                                  iteration += order)
		//for(vector<int>::iterator iteration  = postprocess_iterations.begin(); 
		//						  iteration != postprocess_iterations.end();
		//						++iteration)
		{

		    filenum = (int)(*iteration/InterRecordIterations);
            frame++;

		    cout << "Post processing iteration: " << *iteration << " frame " << filenum 
			     << " (" << frame << "/" 
			     << maxframes << ")"
			     << endl;
            cout.flush();
    		
            //mytimer loadtimer;
            //loadtimer.start();

            //unsigned int nowtime = (unsigned) clock();

		    load_data(theactin, *iteration, false);

            //cout << " Loadtime: " << loadtimer << endl;

            //cout << " Loadtime: " << setprecision(2) <<  ((double)( (unsigned)clock() - nowtime) ) / (double) CLOCKS_PER_SEC << " seconds " << endl;
    		
            cout.flush();

		    if (POST_BMP || POST_REPORTS)
		    {		  
			    nuc_object.segs.addallnodes();  // put node data into segment bins
		    }

		    if (POST_REPORTS)
		    {
			    nuc_object.segs.savereport(filenum);
			    nuc_object.segs.saveSDreport(filenum);
			    nuc_object.segs.saveradialreport(filenum);       
			    nuc_object.segs.saveradialaxisreport(filenum, 0);
			    nuc_object.segs.saveradialaxisreport(filenum, 1);
			    nuc_object.segs.saveradialaxisreport(filenum, 2);
		    }

		    if (POST_BMP)
		    {	
			    // one has to be run in foreground, else will catch up 
                // with Imagemagick and overwrite temp files in use

			    actin::processfgbg xfg, yfg, zfg;

			    xfg = yfg = zfg = actin::runbg;
    			
			    if (Z_BMP)
			    {
				    zfg = actin::runfg;
			    }
			    else 
			    {
				    if (Y_BMP)
					    yfg = actin::runfg;
                    else
   					    xfg = actin::runfg;
			    }

			    theactin.savebmp(filenum, xaxis, xfg, true); 
			    theactin.savebmp(filenum, yaxis, yfg, true); 
			    theactin.savebmp(filenum, zaxis, zfg, true);

			    cout << endl;

		    }

		    if (POST_VTK)
		    {
                CometVtkVis vtkvis(true);  // true?

                mytimer rendertimer;
                rendertimer.start();
			    //cout << "- visualisation: " << filenum << endl;
                //unsigned int nowtime = (unsigned) clock();

			    vtkvis.buildVTK(filenum);
                
                //vtkvis.RestartRenderWindow(); // see if this gets rid of memory leak, slowdown etc.  this segfaults

                cout << " Rendertime: " << rendertimer  << " seconds " << endl;

		    }

		    if (POST_STATS)
		    {
		      // cout << "- statistics: " << filenum << endl;
		      post_stats(theactin, filenum);
		    }
    	}

	}

	if (POST_REPORTS)
	{
		char command1[1024];
		sprintf(command1, "(%s %s*report*.txt 2>/dev/null; mv %s*report*%s %s 2>/dev/null) &"
			,COMPRESSCOMMAND, TEMPDIR,TEMPDIR, COMPRESSEDEXTENSION, REPORTDIR);
		system(command1);
	}

}


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

		load_data(theactin, i, true);

		theactin.load_sym_break_axes();   // overwrite rotation matrixes

		// cout << theactin.sym_break_rotation_to_xy_plane;
					
		nuc_object.segs.addallnodes();  // put node data into segment bins

        nuc_object.segs.savereport(filenum);
        nuc_object.segs.saveSDreport(filenum);
		nuc_object.segs.saveradialreport(filenum);
        nuc_object.segs.saveradialaxisreport(filenum, 0);
        nuc_object.segs.saveradialaxisreport(filenum, 1);
        nuc_object.segs.saveradialaxisreport(filenum, 2);

		// run them in foreground to slow things down
		// so don't overload system with too many bg processes
        // (do we want to change to conditional fg/bg, if not need y, z etc. as above)

		theactin.savebmp(filenum, xaxis, actin::runfg, true);  // was bg
		theactin.savebmp(filenum, yaxis, actin::runfg, true);	// was bg
		theactin.savebmp(filenum, zaxis, actin::runfg, true);
		
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
