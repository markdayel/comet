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

MYDOUBLE TOTAL_SIMULATION_TIME = 20000;  
MYDOUBLE DELTA_T = (MYDOUBLE)0.1;				
MYDOUBLE MAX_DISP_PERDT = (MYDOUBLE)0.01;
MYDOUBLE MAX_DISP_PERDT_DIVSQRTTWO = (MYDOUBLE)0.00707;

int RECORDED_TIMESTEPS=200;			// number of recorded timesteps(data files)

int RESTORE_FROM_ITERATION = 0; // =0 don't load a checkpoint 
int RECORDING_INTERVAL = 0;
int NUMBER_RECORDINGS = 0;

MYDOUBLE FORCE_SCALE_FACT = (MYDOUBLE)0.001;	// convert forces (nom in pN) into node displacements (nom in uM)
									// this is related to effective viscosity and effective size of node

MYDOUBLE XLINK_NODE_RANGE = (MYDOUBLE) 1.0;		// Limit crosslink to within this range
MYDOUBLE NODE_INCOMPRESSIBLE_RADIUS = (MYDOUBLE)0.2;	// repulsion is zero here
// MYDOUBLE NODE_REPULSIVE_MAG = 1000;   // max repulsion (at dist=0)

MYDOUBLE LINK_BREAKAGE_FORCE = (MYDOUBLE) 100;	 // breakage force per link
MYDOUBLE P_LINK_BREAK_IF_OVER = (MYDOUBLE) 0.25;  // probablility that force will break link if over the link breakage force
unsigned int MAX_LINKS_PER_NODE = 100;

MYDOUBLE LINK_TAUGHT_FORCE = (MYDOUBLE) 5;
MYDOUBLE LINK_TAUGHT_RATIO = (MYDOUBLE) 1.1;


MYDOUBLE LINK_FORCE = (MYDOUBLE)0.1;
MYDOUBLE P_XLINK = (MYDOUBLE) 0.5;
MYDOUBLE P_NUC = (MYDOUBLE) 0.08;
MYDOUBLE RADIUS = (MYDOUBLE) 1.0;
MYDOUBLE SEGMENT = (MYDOUBLE) 3.0;

//MYDOUBLE SEG_INCOMP = SEGMENT + NODE_INCOMPRESSIBLE_RADIUS/2;
MYDOUBLE RAD_INCOMP = RADIUS;// + NODE_INCOMPRESSIBLE_RADIUS/2;

//MYDOUBLE NODEMASS = 1.0;
//MYDOUBLE INERTIAL_DAMPING_HALFTIME = 50;

int TOTAL_ITERATIONS ;
int NODE_REPULSIVE_GRIDSEARCH ;
int NODE_XLINK_GRIDSEARCH ;
int NODE_REPULSIVE_RANGE_GRIDSEARCH;

int RADIAL_SEGMENTS = 12;
int NODES_TO_UPDATE = 5000;  //only update the NODES_TO_UPDATE newest nodes

//MYDOUBLE DAMPING_FACTOR = 10;
int CROSSLINKDELAY = 20;  // number of interations before crosslinking 
						 //  (to allow position to be equilibrated to something
						 //       reasonable before locking node in place)

MYDOUBLE NODE_REPULSIVE_MAG = (MYDOUBLE) 0.00000001;
MYDOUBLE NODE_REPULSIVE_RANGE = NODE_INCOMPRESSIBLE_RADIUS*2;

int ASYMMETRIC_NUCLEATION = 0;

int XLINK_NEAREST = 1;

MYDOUBLE VIEW_HEIGHT = 12;

bool USE_THREADS;

int NUM_THREADS;
pthread_attr_t thread_attr;

vector<pthread_t>  threads;

vector<struct thread_data>  collision_thread_data_array;
vector<sem_t> collision_thread_go;
vector<sem_t> collision_data_done;

vector<struct thread_data>  linkforces_thread_data_array;
vector<sem_t> linkforces_thread_go;
vector<sem_t> linkforces_data_done;

vector<struct thread_data>  applyforces_thread_data_array;
vector<sem_t> applyforces_thread_go;
vector<sem_t> applyforces_data_done;

vector<struct thread_data>  compressfiles_thread_data_array;
//vector<sem_t> compressfiles_thread_go;
//vector<sem_t> compressfiles_data_done;

pthread_mutex_t linkstoremove_mutex;

pthread_mutex_t filessavelock_mutex;
pthread_mutex_t filesdonelock_mutex;

// these variables need to be static/global for sharing across threads:

Nodes3d nodegrid;

vector <nodes>	actin::node;
vector <bool>   actin::donenode;	
vector <bool>   actin::repdonenode;	
//Bool2d actin::repulsedone;
Nodes2d actin::recti_near_nodes;
Nodes2d actin::nodes_on_same_gridpoint;
Nodes1d actin::nodes_within_nucleator;
vector <int> actin::nodesbygridpoint;
int actin::iteration_num;

bool actin::isinthread;
vector <nodes*> actin::linkremovefrom;
vector <nodes*> actin::linkremoveto;

int InterRecordIterations = 0;


int load_data(actin &theactin, int iteration);
int save_data(actin &theactin, int iteration);

// main 

int main(int argc, char* argv[])
{

cout << endl;

	if (argc < 2) 
	{
		cerr << "Usage:" << endl << endl << argv[0] << " numThreads [R]" << endl << endl;
		cerr << "where numThreads is the number of threads to use per calculation stage" << endl;
		cerr << "Set numThreads to 0 to run in single threaded mode" << endl;
		cerr << "and 'R' sets use of time-based random number seed" << endl << endl;
	    exit(EXIT_FAILURE);
	}

	ifstream param("cometparams.ini"); 
	if (!param) 
	{
		cerr << "Cannot open cometparams.ini" << endl << endl;
		exit(EXIT_FAILURE);
	}

	 srand( (unsigned) 200 );

	if (argc < 3) 
	{
		cerr << "Static random number seed used" <<  endl;
		srand( (unsigned) 200 );
	} else
	{
		cerr << "Time-based random number seed used" << endl;
		srand( (unsigned)time( NULL ) );
	}

#ifndef _WIN32
	nice(10);
#else
	//SetThreadPriority(GetCurrentThread(),THREAD_PRIORITY_LOWEST);
#endif

  	NUM_THREADS = atoi(argv[1]);

	if (NUM_THREADS < 1)
	{
		cout << "Running in Single Threaded mode" << endl;
		NUM_THREADS = 1;
		USE_THREADS = false;
	}
	else
	{
		if (NUM_THREADS<2)
			cout << "Running in multithreaded mode with 1 thread per stage" << endl;
		else
			cout << "Running in multithreaded mode with " << NUM_THREADS << " threads per stage" << endl;
		USE_THREADS = true;
	}

	threads.resize(NUM_THREADS*5);

	collision_thread_data_array.resize(NUM_THREADS);
	collision_thread_go.resize(NUM_THREADS);
	collision_data_done.resize(NUM_THREADS);

	linkforces_thread_data_array.resize(NUM_THREADS);
	linkforces_thread_go.resize(NUM_THREADS);
	linkforces_data_done.resize(NUM_THREADS);

	applyforces_thread_data_array.resize(NUM_THREADS);
	applyforces_thread_go.resize(NUM_THREADS);
	applyforces_data_done.resize(NUM_THREADS);

	compressfiles_thread_data_array.resize(NUM_THREADS);
//	compressfiles_thread_go.resize(NUM_THREADS);
//	compressfiles_data_done.resize(NUM_THREADS);

	pthread_attr_init(&thread_attr);
	pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_DETACHED);

	pthread_mutex_init(&linkstoremove_mutex,NULL);

	pthread_mutex_init(&filessavelock_mutex,NULL);
	pthread_mutex_init(&filesdonelock_mutex,NULL);
	
	pthread_mutex_lock(&filessavelock_mutex);

	int truethreadnum = 0;

	for (int i = 0; i < NUM_THREADS; i++)
		{
			collision_thread_data_array[i].threadnum = i;
			sem_init(&collision_thread_go[i],0,0);
			sem_init(&collision_data_done[i],0,0);
			
			linkforces_thread_data_array[i].threadnum = i;
			sem_init(&linkforces_thread_go[i],0,0);
			sem_init(&linkforces_data_done[i],0,0);

			applyforces_thread_data_array[i].threadnum = i;
			sem_init(&applyforces_thread_go[i],0,0);
			sem_init(&applyforces_data_done[i],0,0);

			if (USE_THREADS)  // only create threads if more than one cpu
			{
				// start the threads
				pthread_create(&threads[truethreadnum++], &thread_attr, actin::collisiondetectionthread, 
									(void *) &collision_thread_data_array[i]);

				pthread_create(&threads[truethreadnum++], &thread_attr, actin::linkforcesthread, 
									(void *) &linkforces_thread_data_array[i]);

				pthread_create(&threads[truethreadnum++], &thread_attr, actin::applyforcesthread, 
									(void *) &applyforces_thread_data_array[i]);
			}
		}

		// always create compress files thread

		compressfiles_thread_data_array[0].threadnum = 0;
		//sem_init(&compressfiles_thread_go[0],0,0);
		//sem_init(&compressfiles_data_done[0],0,0);

		// make sure semaphores are stopped
		//sem_trywait(&compressfiles_thread_go[0]);
        //sem_trywait(&compressfiles_data_done[0]);

		pthread_create(&threads[truethreadnum++], &thread_attr, actin::compressfilesthread, 
					(void *) &compressfiles_thread_data_array[0]);

	// main parameters:

	nucleator::shape nucshape=nucleator::sphere;  //default to sphere

	// read the parameters file:

	MYDOUBLE MAX_DISP = 1;

	string buffer;
       while (getline(param, buffer)) { 
               istringstream ss(buffer);
               string tag, buff2;
               ss >> tag >> std::ws;
               if (tag.size() == 0 || tag[0] == '#')
                       // skip empty line or comment
                       continue;

               if (tag == "TOTAL_SIMULATION_TIME") {
                       ss >> TOTAL_SIMULATION_TIME;
                       continue;
               } else if (tag == "DELTA_T") {
                       ss >> DELTA_T;
                       continue;
	       } else if (tag == "RECORDING_INTERVAL") {
		   ss >> RECORDING_INTERVAL;
		   continue;
	       } else if (tag == "RESTORE_FROM_ITERATION") {
		   ss >> RESTORE_FROM_ITERATION;
		   continue;
                           } else if (tag == "FORCE_SCALE_FACT") {
                       ss >> FORCE_SCALE_FACT;
                       continue;
			   } else if (tag == "XLINK_NODE_RANGE") {
                       ss >> XLINK_NODE_RANGE;
                       continue;
               } else if (tag == "NODE_INCOMPRESSIBLE_RADIUS") {
                       ss >> NODE_INCOMPRESSIBLE_RADIUS;
                       continue;
               } else if (tag == "P_XLINK") {
                       ss >> P_XLINK;
                       continue;
               } else if (tag == "LINK_BREAKAGE_FORCE") {
                       ss >> LINK_BREAKAGE_FORCE;
                       continue;
               } else if (tag == "P_LINK_BREAK_IF_OVER") {
                       ss >> P_LINK_BREAK_IF_OVER;
                       continue;
			   } else if (tag == "LINK_FORCE") {
                       ss >> LINK_FORCE;
                       continue;
				} else if (tag == "P_NUC") {
                       ss >> P_NUC;
                       continue;
				} else if (tag == "RADIUS") {
                       ss >> RADIUS;
                       continue;
				} else if (tag == "SEGMENT") {
                       ss >> SEGMENT;
                       continue;
				} else if (tag == "MAX_LINKS_PER_NODE") {
                       ss >> MAX_LINKS_PER_NODE;
                       continue;
				} else if (tag == "NODE_REPULSIVE_MAG") {
                       ss >> NODE_REPULSIVE_MAG;
                       continue;
				} else if (tag == "NODE_REPULSIVE_RANGE") {
                       ss >> NODE_REPULSIVE_RANGE;
                       continue;
			   } else if (tag == "LINK_TAUGHT_FORCE") {
                       ss >> LINK_TAUGHT_FORCE;
                      continue;
			   } else if (tag == "ASYMMETRIC_NUCLEATION") {
                       ss >> ASYMMETRIC_NUCLEATION;
                      continue;
			   } else if (tag == "RADIAL_SEGMENTS") {
                       ss >> RADIAL_SEGMENTS;
                      continue;
			   } else if (tag == "XLINK_NEAREST") {
                       ss >> XLINK_NEAREST;
                      continue;
			   } else if (tag == "VIEW_HEIGHT") {
                       ss >> VIEW_HEIGHT;
                      continue;
			   } else if (tag == "NODES_TO_UPDATE") {
                       ss >> NODES_TO_UPDATE;
                      continue;
				} else if (tag == "MAX_DISP") {
                       ss >> MAX_DISP;  // not true global, used to calc MAX_DISP_PERDT
                      continue;
			   } else if (tag == "SHAPE") 
					{
				       ss >> buff2;
					   if (buff2 == "CAPSULE") 
						   nucshape = nucleator::capsule;
					   else
						   nucshape = nucleator::sphere;
                       continue;
					}
       }
 
	param.close();

	// calculate commonly used constant's from parameters:

	TOTAL_ITERATIONS = (int) (((MYDOUBLE)TOTAL_SIMULATION_TIME / (MYDOUBLE)DELTA_T)+0.5);

	NODE_REPULSIVE_GRIDSEARCH = (int) ceil(((MYDOUBLE) NODE_INCOMPRESSIBLE_RADIUS )/GRIDRES) + 1;
	NODE_XLINK_GRIDSEARCH = (int) ceil(((MYDOUBLE) XLINK_NODE_RANGE )/GRIDRES) + 1;
	NODE_REPULSIVE_RANGE_GRIDSEARCH = (int) ceil(((MYDOUBLE) NODE_REPULSIVE_RANGE )/GRIDRES) + 1;

	// loop iterations per recorded timestep
	InterRecordIterations = RECORDING_INTERVAL;
	NUMBER_RECORDINGS = int(TOTAL_ITERATIONS / RECORDING_INTERVAL);
	// InterRecordIterations = (int)
	// (((MYDOUBLE)TOTAL_ITERATIONS / (MYDOUBLE) RECORDED_TIMESTEPS)+0.5 );	
	//SEG_INCOMP = SEGMENT + NODE_INCOMPRESSIBLE_RADIUS/2;
	RAD_INCOMP = RADIUS;//+ NODE_INCOMPRESSIBLE_RADIUS/2;

	MAX_DISP_PERDT = MAX_DISP * DELTA_T;
	MAX_DISP_PERDT_DIVSQRTTWO = MAX_DISP_PERDT / sqrt(2);

	//DAMPING_FACTOR = (DELTA_T * INERTIAL_DAMPING_HALFTIME) / (LN_TWO * NUM_THREADS);



	// create main objects

	actin theactin;
	nucleator nuc_object(nucshape, &theactin);


	cout << "Total simulation time:      " << TOTAL_SIMULATION_TIME << endl;
	cout << "Delta_t:                    " << DELTA_T << endl;
	cout << "MAX_DISP:                   " << MAX_DISP << endl;
	cout << "Nucleator radius:           " << RADIUS << endl;
if (nucshape == nucleator::capsule)
    cout << "Nucleator Segment:          " << SEGMENT << endl;
	cout << "P(nuc):                     " << P_NUC << endl;
	cout << "Force scale factor:         " << FORCE_SCALE_FACT << endl;
	cout << "Crosslink node range:       " << XLINK_NODE_RANGE << endl;
	cout << "Node Incompressible Radius: " << NODE_INCOMPRESSIBLE_RADIUS << endl;
	cout << "Node Repulsion Range:       " << NODE_REPULSIVE_RANGE << endl;
	cout << "Node Repulsion Magnitude:   " << NODE_REPULSIVE_MAG << endl;
	cout << "Max links per node:         " << MAX_LINKS_PER_NODE << endl;
	cout << "Link Taught Force:          " << LINK_TAUGHT_FORCE << endl;
	cout << "Link Taught Ratio:          " << LINK_TAUGHT_RATIO << endl;
	cout << "Link breakage force:        " << LINK_BREAKAGE_FORCE << endl;
	cout << "Link force:                 " << LINK_FORCE << endl;
	cout << "Max P(link):                " << P_XLINK << endl;
	cout << "P(link break):              " << P_LINK_BREAK_IF_OVER << endl << endl;

	theactin.opruninfo << "Total simulation time:      " << TOTAL_SIMULATION_TIME << endl;
	theactin.opruninfo << "Delta_t:                    " << DELTA_T << endl;
	theactin.opruninfo << "MAX_DISP:                   " << MAX_DISP << endl;
	theactin.opruninfo << "Nucleator radius:           " << RADIUS << endl;
if (nucshape == nucleator::capsule)
    theactin.opruninfo << "Nucleator Segment:          " << SEGMENT << endl;
	theactin.opruninfo << "P(nuc):                     " << P_NUC << endl;
	theactin.opruninfo << "Force scale factor:         " << FORCE_SCALE_FACT << endl;
	theactin.opruninfo << "Crosslink node range:       " << XLINK_NODE_RANGE << endl;
	theactin.opruninfo << "Node Incompressible Radius: " << NODE_INCOMPRESSIBLE_RADIUS << endl;
	theactin.opruninfo << "Node Repulsion Range:       " << NODE_REPULSIVE_RANGE << endl;
	theactin.opruninfo << "Node Repulsion Magnitude:   " << NODE_REPULSIVE_MAG << endl;	
	theactin.opruninfo << "Max links per node:         " << MAX_LINKS_PER_NODE << endl;
	theactin.opruninfo << "Link Taught Force:          " << LINK_TAUGHT_FORCE << endl;
	theactin.opruninfo << "Link Taught Ratio:          " << LINK_TAUGHT_RATIO << endl;
	theactin.opruninfo << "Max P(link):                " << P_XLINK << endl;
	theactin.opruninfo << "Link breakage force:        " << LINK_BREAKAGE_FORCE << endl;
	theactin.opruninfo << "Link force:                 " << LINK_FORCE << endl;
	theactin.opruninfo << "P(link break):              " << P_LINK_BREAK_IF_OVER << endl << endl;

	theactin.opruninfo.flush();

	cout << "Total iterations: " << TOTAL_ITERATIONS << endl;
	cout << "Saving snapshot every " << InterRecordIterations  
		<< " iterations (" << NUMBER_RECORDINGS << " total)" << endl;

	cout << "Starting iterations..." << endl << endl;

	unsigned int starttime, endtime, lasttime ,nowtime, lastitertime;

	lastitertime = starttime = (unsigned) time( NULL );

	lasttime=0;

	int lastlinksformed = 0;
	int lastlinksbroken = 0;

	cout.fill(' '); 
	cout.setf(ios::fixed);
	theactin.opruninfo.fill(' ');
	theactin.opruninfo.setf(ios::fixed);

	theactin.newnodescolour.setcol(0);

	MYDOUBLE centre_x,centre_y,centre_z;
	MYDOUBLE last_centre_x, last_centre_y , last_centre_z;
	MYDOUBLE delta_centre_x , delta_centre_y, delta_centre_z;

    last_centre_x = last_centre_y = last_centre_z = 0;
	centre_x = centre_y = centre_z = 0;

	// - - - - - - - - - - 
	// main iteration loop
	// - - - - - - - - - -
	// initialse from a checkpoint if requested
	int starting_iter = 1;
	if(RESTORE_FROM_ITERATION != 0){
	    load_data(theactin, RESTORE_FROM_ITERATION);
	    cout << "* Restored from iteration   "
		 << RESTORE_FROM_ITERATION << endl;
	    srand( (unsigned) 200 );
	    //cout << "reseeded: " << rand() << endl;

	    starting_iter = RESTORE_FROM_ITERATION; 
//	starting_iter = RESTORE_FROM_ITERATION + 1; // don't overwrite
	}
	
	cout << "Itternum|TotNode|NewLinks|RemLinks|Center_X|Center_Y|Center_Z|SnpshTm|SaveNum" << endl << endl;
	
	for(int i=starting_iter;i<=TOTAL_ITERATIONS;i++){
	        nowtime = (unsigned) time( NULL );
		if (nowtime > lasttime)
		{
			lasttime = nowtime;
			theactin.find_centre(centre_x,centre_y,centre_z);
			//delta_centre_x = centre_x - last_centre_x;
			//delta_centre_y = centre_y - last_centre_y;
			//delta_centre_z = centre_z - last_centre_z;

			cout << "I" << setw(7) << i 
			<< "|N"<< setw(6)<< theactin.highestnodecount
			<< "|L+" << setw(6) << (theactin.linksformed-lastlinksformed)/2 << "|L-"
			<< setw(6) << (theactin.linksbroken-lastlinksbroken)/2
			<< "|x " << setw(6) << setprecision(3) << centre_x
			<< "|y " << setw(6) << setprecision(3) << centre_y
			<< "|z " << setw(6) << setprecision(3) << centre_z
			<< "|T" <<  setw(6) <<((unsigned) time( NULL ) - lastitertime) << "\r";
			cout.flush();
		}
		theactin.iterate();
		//theactin.newnodescolour.setcol((MYDOUBLE)i/(MYDOUBLE)TOTAL_ITERATIONS);
		
		if ((i % InterRecordIterations) == 0)
		{
			theactin.setdontupdates();
			theactin.find_centre(centre_x, centre_y, centre_z);
			delta_centre_x = centre_x - last_centre_x;
			delta_centre_y = centre_y - last_centre_y;
			delta_centre_z = centre_z - last_centre_z;
			last_centre_x = centre_x;
			last_centre_y = centre_y;
			last_centre_z = centre_z;

			theactin.opvelocityinfo 
				<< (i*DELTA_T) << "," 
				<< centre_x << "," 
				<< centre_y << "," 
				<< centre_z << "," 
				<< calcdist(delta_centre_x, delta_centre_y, delta_centre_z) << endl;

			cout << "I" << setw(7) << i 
			<< "|N"<< setw(6)<< theactin.highestnodecount
			<< "|L+" << setw(6) << (theactin.linksformed-lastlinksformed)/2 << "|L-"
			<< setw(6) << (theactin.linksbroken-lastlinksbroken)/2
			<< "|x " << setw(6) << setprecision(3) << centre_x
			<< "|y " << setw(6) << setprecision(3) << centre_y
			<< "|z " << setw(6) << setprecision(3) << centre_z
			<< "|T" <<  setw(6) <<((unsigned) time( NULL ) - lastitertime);
			cout.flush();

			cout << "|S " << setw(3) <<  (int)(i/InterRecordIterations)  
				<< "/" << NUMBER_RECORDINGS << endl ;

			theactin.opruninfo << "I" << setw(7) << i 
			<< "|N"<< setw(6)<< theactin.highestnodecount
			<< "|L+" << setw(6) << (theactin.linksformed-lastlinksformed)/2 << "|L-"
			<< setw(6) << (theactin.linksbroken-lastlinksbroken)/2
			<< "|x " << setw(6) << setprecision(3) << centre_x
			<< "|y " << setw(6) << setprecision(3) << centre_y
			<< "|z " << setw(6) << setprecision(3) << centre_z
			<< "|T" <<  setw(6) <<((unsigned) time( NULL ) - lastitertime);
			//<< "|H" <<  setw(6) << theactin.harbinger;
			theactin.opruninfo << "|S" << setw(3) <<  (int)(i/InterRecordIterations)  
				<< "/" << NUMBER_RECORDINGS << endl ;

			//theactin.setnodecols();
			theactin.save((i/InterRecordIterations));
			theactin.savevrml((i/InterRecordIterations));
			theactin.savedata((i/InterRecordIterations));
			
			// test the load/save of data:
			//theactin.loaddata((i/InterRecordIterations));

			save_data(theactin, i);
			//load_data(theactin, i);
			//save_data(theactin, i+1);			
			srand( (unsigned) 200 );
			//cout << "reseeded: " << rand() << endl;

			theactin.savebmp((i/InterRecordIterations), actin::xaxis);
			theactin.savebmp((i/InterRecordIterations), actin::yaxis);
			theactin.savebmp((i/InterRecordIterations), actin::zaxis);

			nuc_object.clearradialsegments();

			if ((i/InterRecordIterations)>1)
			{  // wait for last save to complete...
				pthread_mutex_lock(&filesdonelock_mutex);
				//sem_wait(&compressfiles_data_done[0]);  // wait for the last data write
			}

			compressfiles_thread_data_array[0].startnode = (i/InterRecordIterations);
			//sem_post(&compressfiles_thread_go[0]);
			pthread_mutex_unlock(&filesdonelock_mutex);  // allow thread to grab done lock
			pthread_mutex_unlock(&filessavelock_mutex);  // start the thread

			lastlinksformed = theactin.linksformed;
			lastlinksbroken = theactin.linksbroken;

			lastitertime = (unsigned) time( NULL );

			theactin.opruninfo.flush();
		}
	}

	//sem_wait(&compressfiles_data_done[0]); // wait for the last data write
	pthread_mutex_lock(&filesdonelock_mutex);

	cout << endl << "Done " << endl << endl;

	theactin.saveinfo();

	endtime = (unsigned) time( NULL );

	cout << endl << "Time : " << (endtime-starttime) << " seconds" << endl;	
	theactin.opruninfo << endl << "Time : " << (endtime-starttime) << " seconds" << endl;	

	// Clean up threads
	pthread_attr_destroy(&thread_attr);

	pthread_mutex_destroy(&filesdonelock_mutex);  // allow thread to grab done lock
	pthread_mutex_destroy(&filessavelock_mutex);
	pthread_mutex_destroy(&linkstoremove_mutex);

	for (int i = 0; i < NUM_THREADS; ++i)
	{
		sem_destroy(&collision_thread_go[i]);
		sem_destroy(&collision_data_done[i]);

		sem_destroy(&linkforces_thread_go[i]);
		sem_destroy(&linkforces_data_done[i]);

		sem_destroy(&applyforces_thread_go[i]);
		sem_destroy(&applyforces_data_done[i]);

		//sem_destroy(&compressfiles_thread_go[i]);
		//sem_destroy(&compressfiles_data_done[i]);
	}

	//pthread_exit (NULL);

	exit(EXIT_SUCCESS);
}


string get_datafilename(const int iteration)
{
    stringstream filename;
    filename << "data-" << iteration << ".txt";
    return filename.str();
    
}

int load_data(actin &theactin, int iteration)
{
    string filename = get_datafilename(iteration);
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

/*
#define FP_BITS(fp) (*(DWORD *)&(fp))

static unsigned int fast_sqrt_table[0x10000];  // declare table of square roots 



typedef union FastSqrtUnion
{
  float f;
  unsigned int i;
} FastSqrtUnion;

void build_sqrt_table()
{
  unsigned int i;
  FastSqrtUnion s;
  
  for (i = 0; i <= 0x7FFF; i++)
  {
    
    // Build a float with the bit pattern i as mantissa
    //  and an exponent of 0, stored as 127
    
    s.i = (i << 8) | (0x7F << 23);
    s.f = (float)sqrt(s.f);
    
    // Take the square root then strip the first 7 bits of
    //  the mantissa into the table
    
    fast_sqrt_table[i + 0x8000] = (s.i & 0x7FFFFF);
    
    // Repeat the process, this time with an exponent of 1, 
    //  stored as 128
    
    s.i = (i << 8) | (0x80 << 23);
    s.f = (float)sqrt(s.f);
    
    fast_sqrt_table[i] = (s.i & 0x7FFFFF);
  }
}


inline MYDOUBLE fastsqrt(float n)
{
  if (FP_BITS(n) == 0)
    return 0.0;                 // check for square root of 0
  
  FP_BITS(n) = fast_sqrt_table[(FP_BITS(n) >> 8) & 0xFFFF] | ((((FP_BITS(n) - 0x3F800000) >> 1) + 0x3F800000) & 0x7F800000);
  
  return n;
}

float sse_sqrt(float n)
{
	__asm {
		movd        mm0,[n]
        pfrsqrt     mm1,mm0
        movq        mm2,mm1         // save first 1/sqrt approximation for later
        pfmul       mm1,mm1         // compute 1/(sqrt(x)^2) for pfrsqit1
        pfrsqit1    mm1,mm0         // iterate for accuracy
        pfrcpit2    mm1,mm2         // mm1 = 1 / sqrt(x)
        pfmul       mm0,mm1         // sqrt(x) = x / sqrt(x)
		movd        [n],mm0
	}
		return n;
}

*/
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
