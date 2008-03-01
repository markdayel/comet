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
double MAX_NUC_MOVE = 100;

bool NO_SYMBREAK_ROTATION = false;

double GAUSSFWHM =  0.266;

int BMP_WIDTH  = 800;
int BMP_HEIGHT = 592; // divisible by 16 for movie compression
int BMP_AA_FACTOR = 1;

int VTK_WIDTH  = 1024;
int VTK_HEIGHT = 768;
int VTK_AA_FACTOR = 1;

bool VTK_RAND_NODE_COL = false;

int REFERENCEFRAME = 0;

unsigned int TEXT_POINTSIZE = 20;

bool VTK_NUC_LINKS_ONLY = false;
double VTK_NUC_LINKS_ONLY_AMPLIFY = 20;
double VTK_MAX_NODE_DIST = 10000000;
bool VTK_PLOT_NUC_IMPACTS = false;

bool POST_VTK_VIEW = false;

double COLOUR_GAMMA = 1.6;
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

unsigned int MAX_NODES_TO_TRACK = 30;

bool TRACKS_NO_STATIONARY_NODE = true;
bool BMP_LINKS_BROKEN = false;
bool BMP_TRANSVERSELINKSONLY = false;
bool BMP_RADIALLINKSONLY = false;

bool BLANK_BMP = false;

bool PLOTFORCES = true;

bool ALLOW_HARBINGERS_TO_MOVE = false;

bool NO_X_MOTION = false;

bool   SPECKLE = false;
double SPECKLE_FACTOR = 1;
bool SPECKLEGRID = false;
double SPECKLEGRIDPERIOD = 1.0;
double SPECKLEGRIDTIMEWIDTH = 0.1;
double SPECKLEGRIDSTRIPEWIDTH = 10.0;
bool SPECKLE_NO_ROTATE = false;

bool POLY_FEEDBACK = false;
double POLY_FEEDBACK_DIST = 1.0;
double POLY_FEEDBACK_MIN_PROB = 0.05;
double POLY_FEEDBACK_FACTOR = 4;

double ELLIPSOID_STRETCHFACTOR = 2.0;


bool ROTATION = true;

bool DRAW_CAGE = false;

bool SEGMENT_BINS = true;

bool CAGE_ON_SIDE = true;

bool BMP_FIX_BEAD_MOVEMENT = false;
bool BMP_FIX_BEAD_ROTATION = false;

bool COL_NODE_BY_STRAIN = false;
bool COL_LINK_BY_DIRN = false;
bool COL_INDIVIDUAL_NODES = false;
double NODE_SCALE_GAMMA = 1.0;
bool COL_GREY_BGND = false;
bool NO_BGND = false;
double COL_INDIVIDUAL_SCALE = 8000;

double MOFI =  0.1;

bool X_BMP = true;
bool Y_BMP = true;
bool Z_BMP = true;                 

double TRACK_MIN_RANGE = 0.0;  // frame numbers, can be floating point
double TRACK_MAX_RANGE = 1.0;

int TRACKFRAMESTEP = 5;

double gridscanjitter = 0.01;

int SAVE_DATA_PRECISION	= 4;


bool FORCES_ON_SIDE = false;

bool NO_IMAGE_TEXT = false;
bool NO_COLBAR = false;

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
//bool USE_BREAKAGE_VISCOSITY = false;
//double BREAKAGE_VISCOSITY_THRESHOLD = 1000;
double LINK_BREAKAGE_STRAIN = 1.15;
//double P_LINK_BREAK_IF_OVER =  0.25;  // probablility that force will break link if over the link breakage force
unsigned int MAX_LINKS_PER_NEW_NODE = 100;
unsigned int MAX_LINK_ATTEMPTS = 20;

bool USE_BROWNIAN_FORCES = false;
double BROWNIANFORCESCALE = 0.01;

double NUCLEATOR_INERTIA = 10;
bool VARY_INERT_W_RAD = false;

double DASHPOT_IMPEDANCE = 0.01;

double NUC_LINK_FORCE = 0.25;
double NUC_LINK_BREAKAGE_FORCE = 2;

//double LINK_TAUT_FORCE =  5;
//double LINK_TAUT_RATIO =  1.1;

double IMPOSED_NUC_ROT_SPEED = 1;
bool   IMPOSED_NUC_ROT = false;

double IMPOSED_NUC_DISP_SPEED = 1;
bool   IMPOSED_NUC_DISP = false;


bool   TEST_SQUASH = false;
double TEST_FORCE_INITIAL_MAG = 0;
double TEST_FORCE_INCREMENT = 10;
double TEST_DIST_EQUIL = 0.0001;

bool WRITE_BMPS_PRE_SYMBREAK = false;
bool SYM_BREAK_TO_RIGHT = false;
bool BMP_CENTER_ON_LEFT = false;
bool SYM_IN_COVERSLIP_PLANE = false;

double FINPITCH =  0.1; // 0.1  ; // turns per length
double FINWIDTHANGLE = 30.0; // the angular width of the 
double FINRATIO = 0.1; // 0 to 1, ratio of fin to non-fin

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

bool VTK_SCALECOLOPACITY = false;

double VTK_RIP_Z_ANGLE = 80;

double LINK_FORCE = 0.1;
double P_XLINK =  0.5;
double P_NUC =  0.08;
double MAX_POLYMERISATION_PRESSURE = 1000;

double RADIUS =  1.0;
double CAPSULE_HALF_LINEAR =  6.0;
double COVERSLIPGAP = 10000;    

double NUC_FRICTION_COEFF = 0.1;

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

bool DRAW_COMPRESSIVE_FORCES = true;

double LINK_POWER_SCALE = 0;

bool COMPRESSDATAFILES = true;

int RADIAL_SEGMENTS = 12;
int NODES_TO_UPDATE = 300000;  //only update the NODES_TO_UPDATE newest nodes
double DISTANCE_TO_UPDATE = 0;

int CROSSLINKDELAY = 200;  // number of interations before crosslinking 
						  //  (to allow position to be equilibrated to something
						  //       reasonable before locking node in place)

//double NODE_REPULSIVE_POWER = 2.7;
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

#ifdef USE_GSL_RANDOM

#include <gsl/gsl_rng.h>

//const gsl_rng_type * randomnumT;
gsl_rng * randomnum;

#endif


vector <nodes>	actin::node;
vector <bool>   actin::donenode;	
Nodes2d actin::nodes_by_thread;
Nodes2d actin::recti_near_nodes;
Nodes2d actin::nodes_on_same_gridpoint;
vector<vector<NODEGRIDTYPE<nodes*>*> > actin::gridpointsbythread;
vector <int> actin::nearby_collision_gridpoint_offsets;

vector <nodes> referencenodes;   // used for bitmaps/vtk to normalise energies etc.


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
bool READ_RAW = false;
bool POST_PROCESSSINGLECPU = false;
unsigned int POST_PROCESS_CPUS = 1;
bool VIEW_VTK = false;

bool VTK_NUC_WIREFRAME = true;

int POST_PROC_ORDER = +1;  // +1 = forward, -1 = reverse;

int InterRecordIterations = 0;

bool DISTANCE_TO_UPDATE_reached = false;
bool finished_writing_sym_bitmaps = false;

// functions

bool load_data(actin &theactin, int iteration, const bool &loadscale);
int save_data(actin &theactin, int iteration);
string get_datafilename(const int iteration);
void get_postprocess_iterations(const char *iterdesc, vector<int> &postprocess_iterations, const int& lastframedone);
void postprocess(nucleator& nuc_object, actin &theactin, 
		 vector<int> &postprocess_iterations, char* argv[], const vector<vect> & nodeposns);
void render_raw(actin &theactin);
// void rewrite_symbreak_bitmaps(nucleator& nuc_object, actin &theactin);

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
       str[i] = toupper(str[i]);

   return str;
}
 
string strtolower(string str)
{
   for(unsigned int i=0; i != str.length(); i++)
      str[i] = tolower(str[i]);

   return str;
}

// main 

int main(int argc, char* argv[])
{

    cout << endl; 

    cout << "               comet --- a bead motility simulator" << endl << endl;

    //char BOLD[255];
    //sprintf(BOLD,"%c%c%c",27,'[', 'A'

    if(argc < 2 || argc > 4) 
	{
        cerr << "For a new simulation setup the parameter file 'cometparams.ini'" << endl 
             << "in current directory and type:" << endl 
             //<< endl 
             << "       " << argv[0] << " <numThreads>" << endl 
             //<< endl
	         << "where <numThreads> is the number of CPUs to use" << endl
             << endl
             << "To process an existing dataset type:" << endl 
             <<  "       " << argv[0] << " <command> <frame range>" << endl 
             //<< endl 
		     << "where <command> is 'post' to write bitmap images" << endl 
             << "                   'vtk' to write 3D images " << endl
             << "                   'view' to enter 3D interactive mode" << endl
             << "                   'read_raw' to read 3D raw data (experimental)" << endl
             //<< "e.g.   " << endl 
             << "e.g    " <<argv[0] << " post 1:300      <- write bitmaps for frame 1--300" << endl
             << "       " <<argv[0] << " view 300:300    <- 3D interactive view for frame 300"   << endl
             << "[the range '0:0' can be used to process all frames]" << endl
	         << endl;

	    exit(EXIT_FAILURE);
	}

    if (argc > 1) 
	{
		if (strcmp(argv[1], "sym") == 0 ) 
			REWRITESYMBREAK = true;
        else
		if (strcmp(argv[1], "post") == 0 )
        {
			POST_PROCESS = true;
            POST_BMP = true;
            POST_VTK = false;
        }                                         
        else
        if (strcmp(argv[1], "view") == 0 )
        {
			POST_PROCESS = true;
            POST_BMP = false;
            POST_VTK = true;
            POST_VTK_VIEW = true;
        }
        else
        if (strcmp(argv[1], "vtk") == 0 )
        {
			POST_PROCESS = true;
            POST_BMP = false;
            POST_VTK = true;
            POST_VTK_VIEW = false;
        }
        else
		if (strcmp(argv[1], "psingleb") == 0 )
		{
			POST_PROCESS = true;
			POST_PROCESSSINGLECPU = true;  // this is used for when the multicpu post process calls the worker threads
            POST_BMP = true;
            POST_VTK = false;
		}
        else
        if (strcmp(argv[1], "psinglev") == 0 )
		{
			POST_PROCESS = true;
			POST_PROCESSSINGLECPU = true;  // this is used for when the multicpu post process calls the worker threads
            POST_BMP = false;
            POST_VTK = true;
		}
        else
        if (strcmp(argv[1], "read_raw") == 0 )   
		{
			POST_PROCESS = true;
			POST_PROCESSSINGLECPU = true;  // this is used for when the multicpu post process calls the worker threads
            POST_BMP = false;
            POST_VTK = true;
            READ_RAW = true;
		}
        else
        {

        }


        if (argc > 2)
        {
		if (strcmp(argv[2], "q") == 0 || strcmp(argv[2], "Q") == 0)
			QUIET = true;
        }

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
        POST_PROCESS_CPUS = 2;  // the cluster machines have 2 cpus
    }

    if (    (strcmp( hostname, "medusa.local")       == 0) ||
            (strcmp( hostname, "medusa.kinglab")     == 0) )
    {
        if (POST_VTK)  
            POST_PROCESS_CPUS = 4;
        else
            POST_PROCESS_CPUS = 4;
    }

    if (    (strcmp( hostname, "machaonbook.local")       == 0) )
    {
        if (POST_VTK)  
            POST_PROCESS_CPUS = 1;
        else
            POST_PROCESS_CPUS = 2;
    }


    int nicelevel = 0;

	if (!POST_PROCESS && !REWRITESYMBREAK)
	{	// don't drop priority if re-writing bitmaps
		// because calling thread already has low priority
		// and this would drop it further
		// so process would halt on a single cpu machine
        // while the main process was running

		if (strcmp( hostname, "guanine.ucsf.edu") == 0)
		{
			nicelevel = 0;
		} else
		{
			nicelevel = default_nice_level;
        }
    }


    if (CLUSTER)
    {   // override some defaults if running on a cluster machine

        QUIET = true;   // quiet still produces output, but only once per written file, rather than per second
        NO_IMAGE_TEXT = true;  //remove image text for matrix display
        SYM_BREAK_TO_RIGHT = true;  // offset and fix break direction to right for matrix display
        BMP_CENTER_ON_LEFT = true;
    }

#ifndef _WIN32

	nice(nicelevel);

#endif

    #ifdef USE_GSL_RANDOM

    gsl_rng_env_setup();

    //randomnumT = gsl_rng_default;  // gsl_rng_mt19937 is the default anyway
    //randomnumT = gsl_rng_mt19937;  // Mersenne Twister
    //randomnum = gsl_rng_alloc (randomnumT);

    randomnum = gsl_rng_alloc (gsl_rng_mt19937);

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

    // how many chunks to divide data into (these then get passed to threads)
    // probably best to keep equal to thread #
    // multiples of thread # would be better, except for collisiondeteciton code
    // which restricts nearby nodes to same thread data chunk
	NUM_THREAD_DATA_CHUNKS = NUM_THREADS; 

    // create threads:
	// -- Threading TaskTeam, create and intialise the team
	if (USE_THREADS  && !POST_PROCESS && !REWRITESYMBREAK)
	{
	    thread_queue.create_threads(NUM_THREADS);
	}

	if (USE_THREADS)
		cout << "Running in multithreaded mode with " << NUM_THREADS << " threads" << endl;
	else
		cout << "Running in Single Threaded mode" << endl;


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

	//cout << endl << "Parsing parameters file..." << endl;

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

        bool commandrecognized = true; // workaround for visual c++ compiler bug

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

        else if (tag == "TEXT_POINTSIZE")	  
			{ss >> TEXT_POINTSIZE;}

        else if (tag == "VTK_NUC_LINKS_ONLY") 
			{ss >> buff2; if (buff2=="TRUE") VTK_NUC_LINKS_ONLY = true; else VTK_NUC_LINKS_ONLY = false;} 

        else if (tag == "VTK_NUC_WIREFRAME") 
			{ss >> buff2; if (buff2=="TRUE") VTK_NUC_WIREFRAME = true; else VTK_NUC_WIREFRAME = false;} 

        else if (tag == "VTK_NUC_LINKS_ONLY_AMPLIFY")	  
			{ss >> VTK_NUC_LINKS_ONLY_AMPLIFY;} 

        else if (tag == "VTK_MAX_NODE_DIST")	  
			{ss >> VTK_MAX_NODE_DIST;} 
                                                  
        else if (tag == "VTK_PLOT_NUC_IMPACTS") 
			{ss >> buff2; if (buff2=="TRUE") VTK_PLOT_NUC_IMPACTS = true; else VTK_PLOT_NUC_IMPACTS = false;}

        else if (tag == "SYM_BREAK_TO_RIGHT") 
			{ss >> buff2; if (buff2=="TRUE") SYM_BREAK_TO_RIGHT = true; else SYM_BREAK_TO_RIGHT = false;}

        else if (tag == "BMP_CENTER_ON_LEFT") 
			{ss >> buff2; if (buff2=="TRUE") BMP_CENTER_ON_LEFT = true; else BMP_CENTER_ON_LEFT = false;}

        else if (tag == "BLANK_BMP")  
			{ss >> buff2; if (buff2=="TRUE") BLANK_BMP = true; else BLANK_BMP = false;}

        else if (tag == "SYM_IN_COVERSLIP_PLANE") 
			{ss >> buff2; if (buff2=="TRUE") SYM_IN_COVERSLIP_PLANE = true; else SYM_IN_COVERSLIP_PLANE = false;}   

		else if (tag == "VTK_MOVE_WITH_BEAD")   
			{ss >> buff2; if (buff2=="TRUE") VTK_MOVE_WITH_BEAD = true; else VTK_MOVE_WITH_BEAD = false;}

        else if (tag == "VTK_MIN_PLOT_LINK_FORCE_PCT") 
			{ss >> VTK_MIN_PLOT_LINK_FORCE_PCT;}    

        else if (tag == "VIEW_VTK")  
			{ss >> buff2; if (buff2=="TRUE") VIEW_VTK = true; else VIEW_VTK = false;}

		else if (tag == "VTK_AA_FACTOR")	  
			{ss >> VTK_AA_FACTOR;}

        else if (tag == "VTK_RIP_Z_ANGLE")	  
			{ss >> VTK_RIP_Z_ANGLE;}

        else if (tag == "BMP_AA_FACTOR")	  
			{ss >> BMP_AA_FACTOR;}

		else if (tag == "COLOUR_GAMMA")	  
			{ss >> COLOUR_GAMMA;}     

		else if (tag == "VTK_VIEWANGLE")	  
			{ss >> VTK_VIEWANGLE;}    

        else if (tag == "VTK_RAND_NODE_COL")  
			{ss >> buff2; if (buff2=="TRUE") VTK_RAND_NODE_COL = true; else VTK_RAND_NODE_COL = false;}

		else if (tag == "CROSSLINKDELAY")	  
			{ss >> CROSSLINKDELAY;}
		
		else if (tag == "COVERSLIPGAP")	  
			{ss >> COVERSLIPGAP;}

		else if (tag == "ALLOW_HARBINGERS_TO_MOVE") 
			{ss >> buff2; if (buff2=="TRUE") ALLOW_HARBINGERS_TO_MOVE = true; else ALLOW_HARBINGERS_TO_MOVE = false;}

        else if (tag == "POST_PROC_ORDER") 
			{ss >> POST_PROC_ORDER;}

        else if (tag == "POST_PROCESS_CPUS") 
			{ss >> POST_PROCESS_CPUS;}

		//else if (tag == "POST_BMP") 
		//	{ss >> buff2; if (buff2=="TRUE") POST_BMP = true; else POST_BMP = false;}   

		//else if (tag == "POST_VTK") 
		//	{ss >> buff2; if (buff2=="TRUE") POST_VTK = true; else POST_VTK = false;}

		else if (tag == "POST_STATS") 
			{ss >> buff2; if (buff2=="TRUE") POST_STATS = true; else POST_STATS = false;}

		else if (tag == "POST_REPORTS") 
			{ss >> buff2; if (buff2=="TRUE") POST_REPORTS = true; else POST_REPORTS = false;}

        else if (tag == "BMP_TRACKS") 
			{ss >> buff2; if (buff2=="TRUE") BMP_TRACKS = true; else BMP_TRACKS = false;}

        else if (tag == "TRACKS_NO_STATIONARY_NODE") 
			{ss >> buff2; if (buff2=="TRUE") TRACKS_NO_STATIONARY_NODE = true; else TRACKS_NO_STATIONARY_NODE = false;}

        else if (tag == "TRACK_MIN_RANGE")       
			{ss >> TRACK_MIN_RANGE;} 

        else if (tag == "TRACK_MAX_RANGE")   
			{ss >> TRACK_MAX_RANGE;}

        else if (tag == "TRACKFRAMESTEP")     
			{ss >> TRACKFRAMESTEP;}

        else if (tag == "MAX_NODES_TO_TRACK")       
			{ss >> MAX_NODES_TO_TRACK;} 

        else if (tag == "REFERENCEFRAME")     
			{ss >> REFERENCEFRAME;} 

        else if (tag == "NO_SYMBREAK_ROTATION") 
			{ss >> buff2; if (buff2=="TRUE") NO_SYMBREAK_ROTATION = true; else NO_SYMBREAK_ROTATION = false;}

		else if (tag == "X_BMP")     
			{ss >> buff2; if (buff2=="TRUE") X_BMP = true; else X_BMP = false;}

		else if (tag == "Y_BMP") 
			{ss >> buff2; if (buff2=="TRUE") Y_BMP = true; else Y_BMP = false;}

		else if (tag == "Z_BMP") 
			{ss >> buff2; if (buff2=="TRUE") Z_BMP = true; else Z_BMP = false;}

		else if (tag == "CAGE_ON_SIDE") 
			{ss >> buff2; if (buff2=="TRUE") CAGE_ON_SIDE = true; else CAGE_ON_SIDE = false;}

        else if (tag == "FORCES_ON_SIDE") 
			{ss >> buff2; if (buff2=="TRUE") FORCES_ON_SIDE = true; else FORCES_ON_SIDE = false;}

        else if (tag == "BMP_FIX_BEAD_MOVEMENT")         
			{ss >> buff2; if (buff2=="TRUE") BMP_FIX_BEAD_MOVEMENT = true; else BMP_FIX_BEAD_MOVEMENT = false;}

        else if (tag == "BMP_FIX_BEAD_ROTATION") 
			{ss >> buff2; if (buff2=="TRUE") BMP_FIX_BEAD_ROTATION = true; else BMP_FIX_BEAD_ROTATION = false;} 

		else if (tag == "RESTORE_FROM_FRAME") 
			{ss >> RESTORE_FROM_FRAME;} 

		else if (tag == "NUCLEATOR_FORCES") 
			{ss >> NUCLEATOR_FORCES;} 

        else if (tag == "MAX_NUC_MOVE") 
			{ss >> MAX_NUC_MOVE;}

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

        else
            commandrecognized = false;

        if (!commandrecognized)
        {   
            commandrecognized = true;
		    if (tag == "SEGMENT_BINS") 
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

            else if (tag == "DASHPOT_IMPEDANCE") 
			    {ss >> DASHPOT_IMPEDANCE;}

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
        	
		    //else if (tag == "USE_BREAKAGE_VISCOSITY") 
			//    {ss >> buff2; if(buff2=="TRUE") USE_BREAKAGE_VISCOSITY = true; else USE_BREAKAGE_VISCOSITY = false;} 

//            else if (tag == "BREAKAGE_VISCOSITY_THRESHOLD") 
//			    {ss >> BREAKAGE_VISCOSITY_THRESHOLD;}
    		
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

            else if (tag == "IMPOSED_NUC_DISP_SPEED")  
			    {ss >> IMPOSED_NUC_DISP_SPEED;} 

		    else if (tag == "IMPOSED_NUC_DISP") 
			    {ss >> buff2; if(buff2=="TRUE") IMPOSED_NUC_DISP = true; else IMPOSED_NUC_DISP = false;}

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
    		
		    //else if (tag == "NODE_REPULSIVE_POWER") 
			//    {ss >> NODE_REPULSIVE_POWER;}
    		
		    else if (tag == "NODE_REPULSIVE_RANGE") 
			    {ss >> NODE_REPULSIVE_RANGE;}
    		
		    else if (tag == "ASYMMETRIC_NUCLEATION") 
			    {ss >> ASYMMETRIC_NUCLEATION;} 

            else if (tag == "FINPITCH") 
			    {ss >> FINPITCH;} 

            else if (tag == "FINWIDTHANGLE") 
			    {ss >> FINWIDTHANGLE;} 

            else if (tag == "FINRATIO") 
			    {ss >> FINRATIO;} 
 		
		    else if (tag == "RADIAL_SEGMENTS") 
			    {ss >> RADIAL_SEGMENTS;} 
    		
		    else if (tag == "XLINK_NEAREST") 
			    {ss >> XLINK_NEAREST;} 
    		
		    else if (tag == "VIEW_HEIGHT") 
			    {ss >> VIEW_HEIGHT;} 
    		
		    else if (tag == "NODES_TO_UPDATE") 
		    	{ss >> NODES_TO_UPDATE;} 
    		
		    else if (tag == "DISTANCE_TO_UPDATE") 
			    {ss >> DISTANCE_TO_UPDATE;} 
    		
		    else if (tag == "GAUSSFWHM") 
			    {ss >> GAUSSFWHM;}         

            else if (tag == "FOCALDEPTH") 
			    {ss >> FOCALDEPTH;}         

            else if (tag == "BMP_INTENSITY_OFFSET") 
			    {ss >> BMP_INTENSITY_OFFSET;}
    		
		    else if (tag == "SPECKLE") 
			    {ss >> buff2;if(buff2=="TRUE") SPECKLE = true;else SPECKLE = false;} 

            else if (tag == "SPECKLE_NO_ROTATE") 
			    {ss >> buff2;if(buff2=="TRUE") SPECKLE_NO_ROTATE = true;else SPECKLE_NO_ROTATE = false;} 

            else if (tag == "COL_NODE_BY_STRAIN") 
			    {ss >> buff2;if(buff2=="TRUE") COL_NODE_BY_STRAIN = true;else COL_NODE_BY_STRAIN = false;}

            else if (tag == "COL_LINK_BY_DIRN") 
			    {ss >> buff2;if(buff2=="TRUE") COL_LINK_BY_DIRN = true;else COL_LINK_BY_DIRN = false;}

            else if (tag == "COL_INDIVIDUAL_NODES") 
			    {ss >> buff2;if(buff2=="TRUE") COL_INDIVIDUAL_NODES = true;else COL_INDIVIDUAL_NODES = false;}

            else if (tag == "BMP_LINKS_BROKEN") 
			    {ss >> buff2;if(buff2=="TRUE") BMP_LINKS_BROKEN = true;else BMP_LINKS_BROKEN = false;}

            else if (tag == "BMP_TRANSVERSELINKSONLY") 
			    {ss >> buff2;if(buff2=="TRUE") BMP_TRANSVERSELINKSONLY = true;else BMP_TRANSVERSELINKSONLY = false;}

            else if (tag == "BMP_RADIALLINKSONLY") 
			    {ss >> buff2;if(buff2=="TRUE") BMP_RADIALLINKSONLY = true;else BMP_RADIALLINKSONLY = false;}

            else if (tag == "NODE_SCALE_GAMMA")      
			    {ss >> NODE_SCALE_GAMMA;}

            else if (tag == "VTK_SCALECOLOPACITY") 
			    {ss >> buff2;if(buff2=="TRUE") VTK_SCALECOLOPACITY = true;else VTK_SCALECOLOPACITY = false;}

            else if (tag == "NO_X_MOTION") 
			    {ss >> buff2;if(buff2=="TRUE") NO_X_MOTION = true;else NO_X_MOTION = false;} 

            else if (tag == "COL_GREY_BGND") 
			    {ss >> buff2;if(buff2=="TRUE") COL_GREY_BGND = true;else COL_GREY_BGND = false;}  

            else if (tag == "NO_BGND") 
			    {ss >> buff2;if(buff2=="TRUE") NO_BGND = true;else NO_BGND = false;}  
            
            else if (tag == "COL_INDIVIDUAL_SCALE") 
			    {ss >> COL_INDIVIDUAL_SCALE;}

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
    		
            else if (tag == "DRAW_COMPRESSIVE_FORCES") 
			    {ss >> buff2;if(buff2=="TRUE") DRAW_COMPRESSIVE_FORCES = true; else DRAW_COMPRESSIVE_FORCES = false;}

		    else if (tag == "MOFI") 
			    {ss >> MOFI;}

            else if (tag == "NUC_FRICTION_COEFF") 
			    {ss >> NUC_FRICTION_COEFF;}
    		
		    else if (tag == "NO_IMAGE_TEXT") 
			    {ss >> buff2;if(buff2=="TRUE") NO_IMAGE_TEXT = true; else NO_IMAGE_TEXT = false;} 

            else if (tag == "NO_COLBAR") 
			    {ss >> buff2;if(buff2=="TRUE") NO_COLBAR = true; else NO_COLBAR = false;}
    		
		    else if (tag == "BMP_COMPRESSION") 
			    {ss >> BMP_COMPRESSION;	if (BMP_COMPRESSION > 100) BMP_COMPRESSION = 100; else if (BMP_COMPRESSION < 0)	BMP_COMPRESSION = 0;} 
    		
		    else if (tag == "BMP_OUTPUT_FILETYPE") 
			    {ss >> buff2; istringstream ss3(strtolower(buff2)); ss3 >> BMP_OUTPUT_FILETYPE;} 
    		
		    else if (tag == "SHAPE") 
			    {ss >> buff2; 
                if (buff2 == "CAPSULE") 
                    NUCSHAPE = nucleator::capsule; 
                else if (buff2 == "ELLIPSOID")
                    NUCSHAPE = nucleator::ellipsoid;
                else
                    NUCSHAPE = nucleator::sphere;}

            else if (tag == "ELLIPSOID_STRETCHFACTOR") 
			    {ss >> ELLIPSOID_STRETCHFACTOR;}

            else if (tag == "CLUSTER") 
			    {ss >> buff2;if(buff2=="TRUE") CLUSTER = true; else CLUSTER = false;}

		    else if (tag.find("VIS") == 0)
		    {
			    // VTK stuff, ignore for now---perhaps put VTK stuff in here?
		    }
        else
            commandrecognized = false;
        }
        
        if (!commandrecognized)
		{
			unrecognisedlines += buffer;
			unrecognisedlines += "\n";
			continue;
		}

		// got to here so line was recognised, add it to the unrecognisedlines string

		istringstream ss2(buffer);

		string tag2;

		ss2 >> tag2 >> buff3 >> std::ws;

		//if (!POST_PROCESSSINGLECPU)
        //    cout << setw(30) << tag2.c_str() << setw(7) << buff3.c_str() << endl;

	 }

	param.close();

	if (unrecognisedlines.length()!=0)
	{
		cerr << endl << "Warning---The following lines were unrecognised:" << endl << endl;
		cerr << unrecognisedlines << endl;
		
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


    if (BMP_TRACKS && POST_PROCESS) 
    {
        //POST_BMP = true;
        //POST_VTK = false;
        
        VTK_MOVE_WITH_BEAD = true;
    
        // check if nodetrack file exists, if so use it (and so we don't need to go single thread)
        // else switch to single thread for post processing

        ifstream nodetrackfile("nodetracks.txt", ios::in);
        
        if (!nodetrackfile)
        {
            cout << "Creating Bitmap Tracks.  Post processing will use one thread only." << endl;
		    POST_PROCESSSINGLECPU = true;  // this is used for when the multicpu post process calls the worker threads
            POST_BMP = true;  
        }
        else
        {
            nodetrackfile.close();
  
        }
    }
    
    // change bitmap directory if outputting links
    if (BMP_LINKS_BROKEN)
    {
        sprintf(BITMAPDIR,"%s/bitmaps_links/", path);
    	sprintf(command1, "mkdir %s 2>/dev/null", BITMAPDIR  );
	    system(command1);
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

    if (REFERENCEFRAME && !POST_PROCESS)
    {
        cerr << "REFERENCEFRAME set, but not post-processing.  Ignored." << endl;
        REFERENCEFRAME = 0; // force ignore of reference frame if calculating
    }

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

    if ((RESTORE_FROM_FRAME != 0) || POST_PROCESS || REWRITESYMBREAK)
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

    if (POST_PROCESS || REWRITESYMBREAK)
    {
        cout << "Postprocessing iterations: ";
        if (argc > 2)
	        get_postprocess_iterations(argv[2], postprocess_iterations, lastframedone);
        else
            get_postprocess_iterations("0:0", postprocess_iterations, lastframedone);
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
			cout << endl << "Clear directories and start run(y/n)?";
		    
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

    static nucleator nuc_object;//, &theactin);

    // formatting

	cout.fill(' '); 
	cout.setf(ios::fixed);
	theactin.opruninfo.fill(' ');
	theactin.opruninfo.setf(ios::fixed);

    if (READ_RAW)
    {
        render_raw(theactin);
        exit(EXIT_SUCCESS);
    }


	if (REWRITESYMBREAK)
	{
		//rewrite_symbreak_bitmaps(nuc_object, theactin);
        POST_BMP = true;
        POST_VTK = false;
        POST_PROCESSSINGLECPU = true;
        POST_REPORTS = true;
        POST_PROC_ORDER = 1; // must go forward because we test presence of last bitmap to know we're done
        

		postprocess(nuc_object, theactin, postprocess_iterations,  argv , posn);

        exit(EXIT_SUCCESS);
	}


    // Breakout if we are post processing
	if (POST_PROCESS) // ( !postprocess_iterations.empty() )
	{

		postprocess(nuc_object, theactin, postprocess_iterations,  argv , posn);
        
	    exit(EXIT_SUCCESS); // FIXME: early finish, move to two fcns (ML)
	}
	

	unsigned int starttime, endtime, lasttime ,nowtime;//, lastitertime;

    starttime = (unsigned) time(NULL);

	lasttime=0;

	int lastlinksformed = 0;
	int lastlinksbroken = 0;

    int framemaxvelmoved = 1;                
    vect maxvelmoved;

	theactin.newnodescolour.setcol(0);

	double distfromorigin = 0;

	double x_angle, y_angle, z_angle, tot_rot;

	vect last_center, center, delta_center;
       
    
    // write out parameters to screen and file

	cout << "Total simulation time:      " << TOTAL_SIMULATION_TIME << endl;
	cout << "Delta_t:                    " << DELTA_T << endl;

	theactin.opruninfo.flush();

	cout << "Total iterations: " << TOTAL_ITERATIONS << endl;
	cout << "Frames every " << InterRecordIterations  
		<< " iterations (" << TOT_FRAMES << " total)" << endl;

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

        lastlinksformed = ptheactin->linksformed;
        lastlinksbroken = ptheactin->linksbroken;

        theactin.clear_node_stats();
        theactin.sortnodesbygridpoint();

        // srand( (unsigned) 200 );
 
		//	starting_iter = RESTORE_FROM_FRAME + 1; // don't overwrite

        if (theactin.highestnodecount > ((int)theactin.node.size() - 1000))
            theactin.reservemorenodes(10000);

        theactin.find_center(last_center);

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

#ifdef USE_GSL_RANDOM

    cerr << "Using Mersenne Twister (Gnu Scientific Library) for random number generation" << endl;

gsl_rng_set(randomnum, rand_num_seed);

#else

    cerr << "Warning: Gnu Scientific Library not enabled.  Using rand() for random number generation." << endl;

#endif



//#ifdef USE_MERSENNE_THREADS
//
//    mers_rand.resize(NUM_THREAD_DATA_CHUNKS);
//
//    for (int i=0; i !=  NUM_THREAD_DATA_CHUNKS; ++i)
//    {
//        mers_rand[NUM_THREAD_DATA_CHUNKS].seed( rand_num_seed + NUM_THREAD_DATA_CHUNKS );
//    }
//
//#endif

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
					        
                            sprintf(command1, "movhere"); // call movie writing script
                            system(command1);
							
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
            else
                cout << " ";

            if ( ptheactin->nuc_struck_coverslip)
				cout << "H";
            else
                cout << " ";
            
            
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

                //vect sym_break_direction = posn[filenum] - posn[filenum - 10];

				// set camera rotation for save
    //            if ((SYM_IN_COVERSLIP_PLANE) && (COVERSLIPGAP > 0))
				//    theactin.set_sym_break_axes(true,sym_break_direction);  // constrain to zy plane
    //            else
    //                theactin.set_sym_break_axes(false,sym_break_direction);

				//theactin.save_sym_break_axes();

                // post-process one bitmap to set the symmetry breaking direction                                                   
                sprintf(command1, "%s post %d:%d 1>postsetsymbreak.txt 2>postsetsymbreak.err", argv[0], framemaxvelmoved , framemaxvelmoved  );
				system(command1);

                ptheactin->load_sym_break_axes();

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

				cout << "Spawning new background process to write bitmaps 1--" << filenum-1 
					<< " and continuing..." << endl;

				// look for existance of this file to tell when the process has finished
				if (Z_BMP)
				{
					sprintf(last_symbreak_bmp_filename, "%sz_proj_%05i.%s",BITMAPDIR, 
							filenum-1 ,BMP_OUTPUT_FILETYPE.c_str());
				}
				else 
				{
					if (Y_BMP)
					{
						sprintf(last_symbreak_bmp_filename, "%sy_proj_%05i.%s",BITMAPDIR, 
								filenum-1 ,BMP_OUTPUT_FILETYPE.c_str());
					} else
					{
						sprintf(last_symbreak_bmp_filename, "%sx_proj_%05i.%s",BITMAPDIR, 
								filenum-1 ,BMP_OUTPUT_FILETYPE.c_str());
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
                sprintf(command1, "%s sym 1:%d 1>rewritesymbreak.txt 2>rewritesymbreak.err &", argv[0], filenum-1);
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

            

            for (int frame = 1; frame != TOT_FRAMES+1; ++frame)
            {
                if (velocities[frame].length() > maxvel / 2.0) 
                {
                    maxvelmoved = velocities[frame];
                    framemaxvelmoved = frame;
                    break;
                }
            }


            // write symmetry break time (and direction) to file

            ofstream opsymbreaktime("symbreaktime.txt", ios::out | ios::trunc);
            if (!opsymbreaktime) 
		    { 
                cout << "Unable to open file symbreaktime.txt for output"; 
            }
            else
            {
                opsymbreaktime << "Frame Time Speed Velocity" << endl;
                opsymbreaktime << framemaxvelmoved << " " << framemaxvelmoved * DELTA_T * InterRecordIterations << " " << maxvelmoved << endl;
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

            nuc_object.segs.addallnodes();  // put node data into segment bins

            nuc_object.segs.savereport(filenum);
            nuc_object.segs.saveSDreport(filenum);                
	        nuc_object.segs.saveradialreport(filenum);
            nuc_object.segs.saveradialaxisreport(filenum, 0);
            nuc_object.segs.saveradialaxisreport(filenum, 1);
            nuc_object.segs.saveradialaxisreport(filenum, 2);


			if (theactin.brokensymmetry)
			{	// only save bitmaps if symmetry broken, 
				// 'cause we'll write the others later

				// write the bitmaps, calling imagemagick in background
				// note: this *could* cause a problem if it gets around to the next one before this
				// has finished (because we're sharing x,y,z temp bitmap files)
				// bit this is probably impossible, since we're only calling this
				// once symmetry broken, and by then things are very slow

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

            vect origin(0,0,0);
            if (VIEW_VTK)
                vtkvis.buildVTK(filenum, origin, origin);

            theactin.setdontupdates(); // note this is called *before* the clear_nodes_stats()
                                       // so the stats won't be cleared as they go out of the update range

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

	#ifdef USE_GSL_RANDOM

	gsl_rng_env_setup();

	gsl_rng_free (randomnum);

	#endif
			  
	
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
        ifstrm.clear();
#ifndef _WIN32
	    char command1[1024];
        sprintf(command1, "%s -c %s%s > %s", DECOMPRESSCOMMAND, filename.c_str(), COMPRESSEDEXTENSION, tmpdatafile);
        //cout << command1 << endl;
        system(command1);

        ifstrm.open( tmpdatafile );
        usingtmpdatafile = true;

#endif
    }
                                                    
    if(!ifstrm) 
    {
	    cout << "Unable to open file " << filename.c_str() << " or " << tmpdatafile << " for input" << endl;
#ifdef _WIN32
        cout << "Note that on windows, you must manually unzip any gzipped data files" << endl;
#endif
	    return false;
    }

    string str;    
    // load header
    ifstrm >> str;
    // ensure the identifier for the start of the actin
    if(str.compare("comet:") !=0 )
    {
	    cout << "error in checkpoint file, 'comet:' expected" << endl;
	    return false;
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

    //theactin.sortnodesbygridpoint();   // don't need to do this here
    
    ifstrm.close();

 
    if (usingtmpdatafile)
    {   // delete the temp data file
#ifndef _WIN32
        char command1[1024];
        sprintf(command1, "rm  %s",tmpdatafile);           
        system(command1);
#endif    
    }


    // check the iteration is correct
    if( saved_iteration != iteration )
    {
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
    if(!ofstrm) 
    {
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
    while(chptr != NULL) 
    {
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
    if(strcspn(iterdesc, rangedelim) != strlen(iterdesc) )
    {
        // : range
        int start = 0, step = 0, end = 0;
        if(postprocess_iterations.size() == 1) 
        {
          cerr << "ERROR: too few values in range request (2 min):" 
               << iterdesc << endl;
          exit(EXIT_FAILURE);
        } else if(postprocess_iterations.size() == 2) 
        {
          start = postprocess_iterations[0] * InterRecordIterations;
          step = InterRecordIterations;
          end = postprocess_iterations[1] * InterRecordIterations;
        } else if(postprocess_iterations.size() == 3) 
        {
          start = postprocess_iterations[0] * InterRecordIterations;
          step = postprocess_iterations[1] * InterRecordIterations;
          end = postprocess_iterations[2] * InterRecordIterations;	    
        } else if(postprocess_iterations.size() > 3) 
        {
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
        for(int i= start; i<=end; i+=step)
        {
            postprocess_iterations.push_back(i);                          
        }
    } 
}

// STUB function, doesn't do anything, but fill implementation out here
// 
//void post_stats(actin &, const int filenum)
//{
//  cout << "  saving statistics for " << filenum << endl;  
//
//  char nstats_filename[1024];
//  sprintf(nstats_filename , "%snodestats_%05i.txt", STATSDIR, filenum);
//  //ofstream nout;
//  //nout.open(nstats_filename);
//
//  char lstats_filename[1024];
//  sprintf(lstats_filename , "%slinkstats_%05i.txt", STATSDIR, filenum);
//  //ofstream lout;
//  //lout.open(lstats_filename);
//
//  // STUB for the function to get out statistics of current interest
//  // fill in calculations here
//}

// We may want to have a vtk postprocess as a seperate option?
// otherwise I've lumped it with the bitmap processing here
// the frame counting is seperated, but consider unifying 
// this once we are happy with how it all works.

void postprocess(nucleator& nuc_object, actin &theactin, 
		 vector<int> &postprocess_iterations, char* argv[], const vector<vect> & nodeposns)
{
    if (postprocess_iterations.empty())
        return; // skip if nothing to do

    // load the sym break info

    int framemaxvelmoved = 1;
    vect maxvelmoved;
    int symframe = 1;

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

        symframe = framemaxvelmoved;

        iter = symframe * InterRecordIterations;
    
    }                            
    cout << "Sym break at frame " << symframe  << endl;

    const int maxframe = (int) nodeposns.size() - 1;

    // calculate camera movement

    vector<vect> camerapos(maxframe + 1);
    vector<vect> cameratarget(maxframe + 1);
    vect cameravel;
    cameravel.zero();

    int accelwidth = 100; // how many frames to accelerate over
    int maxaccelframe = (int) (framemaxvelmoved * 1.25);

    cout << "Max camera acceleration will be at frame " << maxaccelframe << endl;

    vect meanvel = nodeposns[maxframe] / (double)(maxframe - framemaxvelmoved) ;

    for (int i = 1; i <= maxframe; ++i)
    {
        camerapos[i]    = camerapos[i-1] + cameravel;
        cameratarget[i] = camerapos[i];
        //cout << cameravel.length() << endl;

        cameravel = meanvel * 1.0 / ( 1.0 + exp( - 6.0 * (  2 * (double)(i - maxaccelframe) / (double)accelwidth) ));
    }

    if (REFERENCEFRAME)
    {
        cout << "Loading frame " << REFERENCEFRAME << " for reference" << endl;

        if (!load_data(theactin, REFERENCEFRAME * InterRecordIterations, false))
        {
            cerr << "Unable to open reference frame, ignoring" << endl;
        }

        REFERENCEFRAME = 0;

        // copy the nodes data
        referencenodes.assign(theactin.node.begin(), theactin.node.end());
    }



    //iter =  180 * InterRecordIterations;

    cout << "Loading frame " << (iter/InterRecordIterations) << " for scaling and to set sym break axis" << endl;

    load_data(theactin, iter, false);

    theactin.BMP_intensity_scaling = true;

    theactin.savebmp(framemaxvelmoved, xaxis, actin::runfg, false); 
    theactin.savebmp(framemaxvelmoved, yaxis, actin::runfg, false); 
    theactin.savebmp(framemaxvelmoved, zaxis, actin::runfg, false);

    cout << endl;

    nuc_object.segs.addallnodes();                                  
    nuc_object.segs.set_scale_factors();

    vect sym_break_direction = nodeposns[framemaxvelmoved] - nodeposns[framemaxvelmoved - 20];

    if ((SYM_IN_COVERSLIP_PLANE) && (COVERSLIPGAP > 0))
	    ptheactin->set_sym_break_axes(true, sym_break_direction);  // constrain to zy plane
    else
        ptheactin->set_sym_break_axes(false,sym_break_direction);

    if (!POST_PROCESSSINGLECPU)
    {   // prevent race condition when doing the mulit cpu post process---the caller thread does the symmetry breaking save
        ptheactin->save_sym_break_axes();
    }

    char command1[1024];

    if (!POST_VTK_VIEW && !POST_PROCESSSINGLECPU && (postprocess_iterations.size() > POST_PROCESS_CPUS) )
	{   // divide up tasks and spawn worker processes

		int firstframe =  *postprocess_iterations.begin()  / InterRecordIterations;
		int lastframe  = *(postprocess_iterations.end()-1) / InterRecordIterations;

	    for (unsigned int i = 1; i != POST_PROCESS_CPUS; ++i)
	    {
            if (POST_BMP)
		        sprintf(command1, "%s psingleb %i:%i:%i &", argv[0], firstframe + i, POST_PROCESS_CPUS, lastframe);
            else
                sprintf(command1, "%s psinglev %i:%i:%i &", argv[0], firstframe + i, POST_PROCESS_CPUS, lastframe);
		    system(command1);
		    //cout << command1 << endl;
	    }

        // do first one in foreground:

        if (POST_BMP)
	        sprintf(command1, "%s psingleb %i:%i:%i", argv[0], firstframe , POST_PROCESS_CPUS, lastframe);
        else
            sprintf(command1, "%s psinglev %i:%i:%i", argv[0], firstframe , POST_PROCESS_CPUS, lastframe);

        
	    system(command1);
 
	} else
    {

		int filenum;

		cout << "Post processing " 
		 << postprocess_iterations.size()
		 <<  " data sets." << endl;
	    
        const int maxframes = (int) postprocess_iterations.size();
		int frame = 0;

        theactin.node_tracks[xaxis].resize(0);
        theactin.node_tracks[yaxis].resize(0);
        theactin.node_tracks[zaxis].resize(0);



        theactin.savenodetracks = true;  // save new nodetracks if not exist


        if (BMP_TRACKS)
        {
            // try to load existing tracks file

            theactin.savenodetracks = !theactin.load_nodetracks();  // if load fails, we will make new tracks

            //cout << "Tot Nodes to track... " << ptheactin->node_tracks[xaxis].size() << endl;

            if (theactin.savenodetracks)    // make new tracks
            {
                // load the first tracking frame to choose new nodes
                
                iter = InterRecordIterations * (int)ceil(TRACK_MAX_RANGE);

                cout << "Loading frame " << (int)ceil(TRACK_MAX_RANGE) << " for track node selection" << endl;

                load_data(theactin, iter, false);

                theactin.set_nodes_to_track(xaxis); // only track x-axis for now

                theactin.savebmp(filenum, xaxis, actin::runfg, false);  // sets the initial offset position

                theactin.node_tracks[xaxis].resize(0); // clear the node tracks (after the last bmp call)

                POST_PROC_ORDER = 1;  // must go forwards for tracks
            }
            else
            {
                cout << "'nodetracks.txt' file found,  old tracks loaded" << endl;
            }                                          
    
        }

        //cout << "Tot Nodes to track... " << ptheactin->node_tracks[xaxis].size() << endl;

        // we can't create the vkt object here because we can't reuse the renderer without
        // causing segfaults for some reason
        // workaround is to create the vtk object within the loop, so that
        // the destructor clears out and we recreate the object every time

        CometVtkVis vtkvis(false); // should this be false? (i.e. no render window)


        vector<int>::iterator start = postprocess_iterations.begin();
        vector<int>::iterator end = postprocess_iterations.end();

        if (POST_PROC_ORDER == -1)
        {   // if going backwards, reverse the start and end
            start = postprocess_iterations.end();
            end = postprocess_iterations.begin();
        }

        for(vector<int>::iterator iteration  = start; 
		                          iteration != end;
                                  iteration += POST_PROC_ORDER)
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

		    load_data(theactin, *iteration, false);

            //cout << " Loadtime: " << loadtimer << endl;
    		
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

			    //xfg = yfg = zfg = actin::runbg;
    			//
			    //if (Z_BMP)
			    //{
				   // zfg = actin::runfg;
			    //}
			    //else 
			    //{
				   // if (Y_BMP)
					  //  yfg = actin::runfg;
       //             else
   				//	    xfg = actin::runfg;
			    //}

                xfg = yfg = zfg = actin::runfg;

			    theactin.savebmp(filenum, xaxis, xfg, true); 
			    theactin.savebmp(filenum, yaxis, yfg, true); 
			    theactin.savebmp(filenum, zaxis, zfg, true);

			    cout << endl;

                if (BMP_TRACKS && theactin.savenodetracks)
                {
                    theactin.save_nodetracks();  // save the node tracks
                }

		    }

		    if (POST_VTK)
		    {
                //CometVtkVis vtkvis(true);  // true?

                mytimer rendertimer;
                rendertimer.start();
			    //cout << "- visualisation: " << filenum << endl;
                //unsigned int nowtime = (unsigned) clock();

			    vtkvis.buildVTK(filenum, camerapos[filenum], cameratarget[filenum]);
                
                //vtkvis.RestartRenderWindow(); // see if this gets rid of memory leak, slowdown etc.  this segfaults

                cout << " Rendertime: " << rendertimer  << " seconds " << endl;

		    }

		    //if (POST_STATS)
		    //{
		      // cout << "- statistics: " << filenum << endl;
		    //  post_stats(theactin, filenum);
		    //}
    	}

	}

	if (POST_REPORTS)
	{
		char command5[1024];
		sprintf(command5, "(%s %s*report*.txt 2>/dev/null; mv %s*report*%s %s 2>/dev/null) &"
			,COMPRESSCOMMAND, TEMPDIR,TEMPDIR, COMPRESSEDEXTENSION, REPORTDIR);
		system(command1);
	}

}



void render_raw(actin &theactin)
{   
    // bit ugly right now---just read in images from 'raw' dir

#define RAWTYPE float

    

    vector<vector<vector<RAWTYPE> > > raw_data;

    int raw_width, 
	    raw_height;

	double	nm_per_pixel,
            nm_per_slice;

    int raw_firstframe,
        raw_lastframe;

    int totframes;

    char raw_filepattern[2048], filename[2048], buff[2048];


    int filenum;

    // read data settings in from raw files instead of previously calculated data

 	ifstream raw_info("raw/frameinfo.txt", ios::in);

    if (!raw_info) 
    { 
	    cout << "Unable to open file 'frameinfo.txt' for input" << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        string line;
        getline(raw_info,line);

        //cout << line;

	    raw_info    >> raw_width 
				    >> raw_height
				    >> nm_per_pixel
                    >> nm_per_slice
                    >> raw_firstframe
                    >> raw_lastframe
                    >> raw_filepattern;

	    raw_info.close();  

        cout << setprecision(2);
        cout << "File info: "
            << "Width:       " << raw_width << endl
			<< "Height:      " << raw_height << endl
			<< "nm/pixel:    " << nm_per_pixel << endl
            << "nm/slice:    " << nm_per_slice << endl
            << "First frame: " << raw_firstframe << endl
            << "Last Frame:  " << raw_lastframe << endl
            << "Filepattern: " << raw_filepattern << endl;

        
        
        if ((( raw_width < 1 ) || (raw_width > 2000)) ||
            (( raw_height < 1 ) || (raw_height > 2000)) || 
            (( raw_lastframe < 1 ) || (raw_lastframe > 2000)))
        {
            cout << " Parameters out of range.  Aborting read. " << endl;
            exit(EXIT_FAILURE);
        }

        cout << endl;


        totframes = (raw_lastframe - raw_firstframe) + 1;

        cout << "Raw type is set to " << sizeof(RAWTYPE) * 8 << " bits per pixel" << endl;

        double memreq=sizeof(RAWTYPE) * raw_width * raw_height * raw_lastframe;

        cout << "Allocating memory required for image : " <<  ((memreq / 1024) / 1024) << " MB" << endl;

        
        raw_data.resize(totframes);

        for (int z=0; z!=totframes; z++)
	    {
            raw_data[z].resize(raw_height);

            for (int x=0; x!=raw_height; x++)  // allocate data grid
            {
	            raw_data[z][x].resize(raw_width);
    	        
                //for (int y=0; y!=raw_height; y++)
	            //{
		        //    raw_data[z][x][y].resize(raw_lastframe);
	            //}
            }
        }

        cout << "Memory allocated" << endl;

    }

        CometVtkVis vtkvis(false); // should this be false? (i.e. no render window)


        for(filenum = raw_firstframe; filenum <= raw_lastframe ; filenum++ )
		{

            sprintf(buff, raw_filepattern, filenum);
            sprintf(filename, "raw/%s", buff);


            cout << "Reading file " << filename << "...";
            cout.flush();

            ifstream raw_in(filename, ios::in);

            // find filesize

            raw_in.seekg(0,ios::end);
            int file_size=raw_in.tellg();
            raw_in.seekg(0);

            int z = filenum - raw_firstframe;  // index from zero

            if ((z < 0) || (z > raw_lastframe))
            {
                cout << "Z out of range" << endl;
                exit(EXIT_FAILURE);
            }


            for (int y=0; y!=raw_height; y++)  // allocate data grid
            {
                //for (int y=0; y!=raw_height; y++)
	            //{
                    //for (int z=0; z!=raw_height; z++)
	                //{
                        raw_in.read((char*)(&raw_data[z][y][0]),raw_width*sizeof(RAWTYPE));
	                //}
                //}
            }

            

            int posn = (int)raw_in.tellg();
            if (posn == -1)
            {
                cout << endl << "Warning: End of file reached before data read in!" << endl;
            }
            else
            if  (posn != file_size)
            {
                cout << endl << "Warning: Finished reading " <<  file_size - posn << " bytes short of end of file " << endl;
            }
            else
            {
                cout << "read OK." << endl;
            }

            cout.flush();

            raw_in.close();
        }


    //CometVtkVis vtkvis(false); // should this be false? (i.e. no render window)


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
