/*

    comet - an actin-based motility simulator
    Copyright (C) 2005 Mark J Dayel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>

*/

#include "stdafx.h"
#include "actin.h"
#include "rotationmatrix.h"
#include <assert.h>


actin::actin(void)
{

    if (!REWRITESYMBREAK && !POST_PROCESS)
    {
	    opruninfo.open("comet_run_info.txt", ios::out | ios::trunc);
	    if (!opruninfo) 
	    { cout << "Unable to open file " << "comet_run_info.txt" << " for output"; return;}
    }
    else
    {
	    opruninfo.open("/dev/null", ios::out | ios::trunc);
    }

	opruninfo << "comet " << endl << "Mark J Dayel" << endl << "(" << __DATE__ << " " << __TIME__ << ") " << endl << endl;

	opruninfo << "PID: " << getpid() << endl;

	cout << endl;
	cout << "comet " << endl << "Mark J Dayel" << endl << "(" << __DATE__ << " " << __TIME__ << ") " << endl << endl;
	cout.flush();

    if (!REWRITESYMBREAK && !POST_PROCESS) // don't clobber if just doing post-process
    {
	    opvelocityinfo.open( VELOCITIESFILE , ios::out | ios::app);  // (doesn't work) don't truncate---append (in case continuing a run)
	    if (!opvelocityinfo) 
	    {
	        cout << "Unable to open file " << "velocities.txt" << " for output";
	        return;
	    }

		//opinfo.open("info.txt", ios::out | ios::trunc);
	 //   if (!opinfo) 
	 //   {
	 //       cout << "Unable to open file " << "velocities.txt" << " for output";
	 //       return;
	 //   }

    }
    else
    {
	    opvelocityinfo.open("/dev/null", ios::out | ios::trunc);
		//opinfo.open("/dev/null", ios::out | ios::trunc);
    }

//    testaverageposn.zero();
//    lasttestaverageposn.zero();
    
    if (!REWRITESYMBREAK && !POST_PROCESS)
    {
	    opvelocityinfo << "time,x,y,z,vel" << endl;
    }

	crosslinknodesdelay.resize(CROSSLINKDELAY);

	//doreportiteration = 9999999;

	donenode.resize(MAXNODES);

	node.resize(MAXNODES);// ,nodes(0,0,0, this));

	for (int i=0; i<MAXNODES; i++)
	{
		node[i].nodenum=i;
		//node[i].ptheactin=this;
	}

	cout << "     GridExtent : +/-" << GRIDBOUNDS << " uM" << endl;
	//cout << " GridResolution : " << GRIDRES << " uM" << endl;
	//cout << "TotalGridpoints : " << (GRIDSIZE*GRIDSIZE*GRIDSIZE) << endl;

	//cout << "Zeroing Grid...";
	//cout.flush();
    
    	     
    #ifdef NODE_GRID_USE_ARRAYS

	    nodegrid = new NG1d[(GRIDSIZE+1)*(GRIDSIZE+1)*(GRIDSIZE+1)];

	    clear_nodegrid();

    #else

        

        //#define NODEGRID(i,j,k)	 nodegrid[(i)][(j)][(k)]
	    nodegrid.resize(GRIDSIZE+1);
	    for (int i=0; i!=(GRIDSIZE+1); i++)  // allocate nodegrid and set nodegrid pointers to null
	    {
		    nodegrid[i].resize(GRIDSIZE+1);
		    for (int j=0; j!=(GRIDSIZE+1); j++)
		    {
			    nodegrid[i][j].resize(GRIDSIZE+1);
		    }
	    }



    #endif

        //reserve vectors for area close to nucleator (only for vector type)

    #ifndef NODEGRIDTYPELIST


        for (int i=0; i!=(GRIDSIZE+1); i++)  
	    {
            double x = (i - (GRIDSIZE/2)) * GRIDRES;
		    for (int j=0; j!=(GRIDSIZE+1); j++)
		    {
                double y = (j - (GRIDSIZE/2)) * GRIDRES;
                for (int k=0; k!=(GRIDSIZE+1); k++)
                {
                    double z = (k - (GRIDSIZE/2)) * GRIDRES;

                    if ((!REWRITESYMBREAK && !POST_PROCESS) &&
                        ((fabs(x) < RADIUS * 8) ||
                         (fabs(y) < RADIUS * 8) ||
                         (fabs(z) < RADIUS * 8)))
                    {
                        NODEGRID(i,j,k).reserve(64);
                    }
                }
		    }
	    }

    #endif

	//cout << "Done" << endl << endl;

	speckle_array_size = 10000;
	speckle_array.resize(speckle_array_size);


	double x1, x2, w, y1; //, y2;

	for (int i=0; i !=speckle_array_size; ++i)
	{   // gaussian random numbers for intensities
		do
		{
			do
			{
				x1 = 2.0 * rand_0to1() - 1.0;
				x2 = 2.0 * rand_0to1() - 1.0;
				w = x1 * x1 + x2 * x2;
			} while ( w >= 1.0 );

			w = sqrt( (-2.0 * log( w ) ) / w );

			y1 = x1 * w;

			y1 *= (1.0/7.0);
			y1 += 0.5;

		} while ( (y1<0) || (y1>1) );

		if (prob_to_bool(SPECKLE_FACTOR))
			speckle_array[i] = y1;
		else
			speckle_array[i] = 0;
	}


	//nucleatorgrid.reserve(5000);
	
	recti_near_nodes.resize(NUM_THREAD_DATA_CHUNKS);
    nodes_on_same_gridpoint.resize(NUM_THREAD_DATA_CHUNKS);
    nodes_by_thread.resize(NUM_THREAD_DATA_CHUNKS);
	

//    recti_near_nodes_size.resize(NUM_THREAD_DATA_CHUNKS);
//    nodes_on_same_gridpoint_size.resize(NUM_THREAD_DATA_CHUNKS);

    linkremovefrom.resize(NUM_THREAD_DATA_CHUNKS);
    linkremoveto.resize(NUM_THREAD_DATA_CHUNKS);

	for (int i = 0; i < NUM_THREAD_DATA_CHUNKS; ++i)
	{
        recti_near_nodes[i].resize(0);
        nodes_on_same_gridpoint[i].resize(0);
        nodes_by_thread[i].resize(0);
        linkremovefrom[i].resize(0);
        linkremoveto[i].resize(0);
        
        recti_near_nodes[i].reserve(1024);
		nodes_on_same_gridpoint[i].reserve(1024);
        nodes_by_thread[i].reserve(2 * MAXNODES / NUM_THREAD_DATA_CHUNKS);  // poss should / by NUM_THREAD_DATA_CHUNKS
        linkremovefrom[i].reserve(128);
        linkremoveto[i].reserve(128);
	}

	//nodes_within_nucleator.reserve(10000);

    linkformto.resize(0);
    linkformto.reserve(128);


	//cout << "Memory for Grid : " << (sizeof(nodes*)*GRIDSIZE*GRIDSIZE*GRIDSIZE/(1024*1024)) << " MB" << endl;

	//cout << "Memory for nodes: " << (node.size() * sizeof(node[0])/(1024*1024)) << " MB" << endl;


	//findnearbynodes_collision_setup(NODE_REPULSIVE_RANGE_GRIDSEARCH);  


	highestnodecount = 0;
	iteration_num = 0;
	linksformed = 0;
	linksbroken = 0;
	nexttocrosslink = 0;
	symbreakiter = 0;
    lowestnodetoupdate = 0;

	lastsorthighestnode = 0;
	attemptedpolrate = polrate = 0;

	newnodescolour.setcol(0);

	currentlyusingthreads = false;

    nuc_struck_coverslip = false;

    sym_break_x_angle = 0.0;
    sym_break_y_angle = 0.0;
    sym_break_z_angle = 0.0;

	//debug:

	//num_rotate = 0;
	//num_displace = 0;

    imageR.resize(BMP_WIDTH);
	imageG.resize(BMP_WIDTH);
	imageB.resize(BMP_WIDTH);

	for (int x = 0; x<BMP_WIDTH; x++)
	{
		imageR[x].resize(BMP_HEIGHT);
		imageG[x].resize(BMP_HEIGHT);
		imageB[x].resize(BMP_HEIGHT);

		for (int y = 0; y<BMP_HEIGHT; y++)
		{
			imageR[x][y]=0;
			imageG[x][y]=0;
			imageB[x][y]=0;
		}
	}

    //GetProcessId(GetCurrentThread());

	// setup the temp bitmap files with unique names:

	unsigned int filepidnum = (unsigned int)getpid();
	unsigned int filetidnum = (unsigned int)time(NULL);

	sprintf(temp_BMP_filename_x, "%sx_temp_%u_%u.bmp", TEMPDIR, filepidnum, filetidnum);
	sprintf(temp_BMP_filename_y, "%sy_temp_%u_%u.bmp", TEMPDIR, filepidnum, filetidnum);
	sprintf(temp_BMP_filename_z, "%sz_temp_%u_%u.bmp", TEMPDIR, filepidnum, filetidnum);

	outbmpfile_x.open(temp_BMP_filename_x, ios::out | ios::binary | ios::trunc);
	outbmpfile_y.open(temp_BMP_filename_x, ios::out | ios::binary | ios::trunc);
	outbmpfile_z.open(temp_BMP_filename_z, ios::out | ios::binary | ios::trunc);

	if ((!outbmpfile_x) || (!outbmpfile_y) || (!outbmpfile_z)) 
	{ 
        cout << "Unable to open temp bitmap file for output" << endl;
        cout << temp_BMP_filename_x << " " << temp_BMP_filename_y << " "<< temp_BMP_filename_z << endl;
        return;
    }

	// we write the header now, and just keep changing the pixel part as we write the frames

	writebitmapheader(outbmpfile_x, BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_y, BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_z, BMP_WIDTH, BMP_HEIGHT);
    
    imageRmax.resize(3);
    imageGmax.resize(3);
    imageBmax.resize(3);

    for (int proj = 0; proj != 3; ++proj)
    {
	    imageRmax[proj] = 1000/INIT_R_GAIN;
	    imageGmax[proj] = 1000/INIT_G_GAIN;
	    imageBmax[proj] = 1000/INIT_B_GAIN;  // prevent over sensitivity
    }

    brokensymmetry = false;

    BMP_intensity_scaling = true;

    world_to_nuc_rot.settoidentity();
    nuc_to_world_rot.settoidentity();

    gridpointsbythread.resize(NUM_THREAD_DATA_CHUNKS);

    for(int i = 0; i != NUM_THREAD_DATA_CHUNKS; ++i)
    {
        gridpointsbythread[i].resize(0);
        gridpointsbythread[i].reserve(2048);
        
    }
    currentsmallestgridthread = 0;

    for (int i=0; i != NUM_THREAD_DATA_CHUNKS; ++i)
    {
        collision_thread_data_array[i].endnode = node.begin();
        collision_thread_data_array[i].startnode = node.begin();
        collision_thread_data_array[i].threadnum = i;

        linkforces_thread_data_array[i].endnode = node.begin();
        linkforces_thread_data_array[i].startnode = node.begin();
        linkforces_thread_data_array[i].threadnum = i;

        applyforces_thread_data_array[i].endnode = node.begin();
        applyforces_thread_data_array[i].startnode = node.begin();
        applyforces_thread_data_array[i].threadnum = i;

    }


    lasttestsurfacesavedposn = testsurfaceposn = DBL_MAX;
    testsurfacerotation = 0;
    testangle = 0;
    testforcemag = 0;
     
    testdirection = vect(0,0,1);
    testangle = cos(30 * PI/360);

    test_equilibrating = false;

    node_tracks.resize(3);

    node_tracks[xaxis].resize(0);
    node_tracks[yaxis].resize(0);
    node_tracks[zaxis].resize(0);

    node_tracks[xaxis].reserve(10240);
    node_tracks[yaxis].reserve(10240);
    node_tracks[zaxis].reserve(10240);

    savenodetracks = true;


    collisiondetection_setrepforcelookup();

} 

actin::~actin(void)
{
	opruninfo.close();
	opvelocityinfo.close();
	//opinfo.close();
	outbmpfile_x.close();
	outbmpfile_y.close();
	outbmpfile_z.close();

#ifdef NODE_GRID_USE_ARRAYS
	delete [] nodegrid;
#endif

	// delete the temp bitmap files

	char command1[1024];

    // start a background process to clear the temp directory after a delay
    // (delay required because background imagemagick processes may still be accessing the bitmaps)
	sprintf(command1, "(sleep 120 ; rm %s*.bmp %s*.png %s*.txt 2>/dev/null ; rmdir %s ) &", TEMPDIR,TEMPDIR,TEMPDIR,TEMPDIR);  
	system(command1);

    system("stty sane 2>/dev/null");   // fix for something that messes the terminal up (kbhit?)  pipe error to /dev/null in case running in background
}


void actin::crosslinknewnodes(const int &numnewnodes)
{	
	if (numnewnodes == 0)
		return;  // :)
	
	vect nodeposvec;
	double distsqr;

	vect disp;

	// and link the new nodes...
	for (int i = nexttocrosslink; i != nexttocrosslink + numnewnodes; ++i)
	{
		attemptedpolrate++;

		if (!node[i].harbinger ||
            !node[i].polymer)
			continue;   // look for surviving harbingers to crosslink
                        // those that have been killed are no longer harbingers (check this?)
                        // so will abort here (so attemptedpolrate will
                        // be greater than polrate)

		//cout << "linking node " << i << endl;
		//cout << "created node at " << node[i].x << "," << node[i].y << "," << node[i].z << endl;

		node[i].harbinger = false;  // now our node exists, it's no longer a harbinger
		polrate++;     // node is created, so increase the actual pol rate monitoring variable

		if (findnearbynodes(node[i],NODE_XLINK_GRIDSEARCH,0)==0) // find nodes close by
			continue;	// skip if zero											

		// collect the nodes within link range (euclidian dist)
        // ready for sorting

        nodeposvec = node[i];

		linkformto.resize(0);

		for (vector <nodes*>::iterator nearnode  = recti_near_nodes[0].begin();
                                       nearnode != recti_near_nodes[0].end(); 
                                     ++nearnode)
		{

			if ( ((*nearnode)->harbinger) ||
                (!(*nearnode)->polymer) )
				continue;  // only crosslink to real nodes

			disp = (*(*nearnode)) - nodeposvec;

			distsqr = disp.sqrlength();

			if (distsqr < SQRT_ACCURACY_LOSS)
				continue;

			if ( distsqr < XLINK_NODE_RANGE * XLINK_NODE_RANGE)
			{
				linkformto.push_back(linkform( (*nearnode)->nodenum , distsqr ));
			}
		}


		// put them in order:
		if (XLINK_NEAREST)  // crosslink in order of distance?
			sort(linkformto.begin(), linkformto.end(), linkform::CompareDistance);
		else				// or in random order
			random_shuffle(linkformto.begin(), linkformto.end());

        unsigned int linkattempts = 0, successfullinks = 0;

		// and crosslink:
		for (vector <linkform>::iterator linkto  = linkformto.begin();
                                         linkto != linkformto.end();
                                       ++linkto)
		{
            if (( linkattempts++ == MAX_LINK_ATTEMPTS) || 
                (successfullinks == MAX_LINKS_PER_NEW_NODE))
                break; 

            if ( addlinks(node[i], node[linkto->nodenum]) )
            {   
                successfullinks++;
                opinfo << ( (node[i]-node[linkto->nodenum]).length() ) << " ";
            }
        }

        if (successfullinks > 0)
			opinfo << endl;

		if (!ALLOW_HARBINGERS_TO_MOVE)
		{   // this was a harbinger---if they were not having applyforces called
            // then force vectors will have built up so zero them now
			node[i].setunitvec();
			node[i].clearforces();  
		}

        node[i].nucleator_stuck_position = node[i];

        if (STICK_TO_NUCLEATOR)
        {
            node[i].stucktonucleator = true;    
        }

	}

	nexttocrosslink += numnewnodes;

	return;
}


void actin::getnodeextents(double & minx,double & miny,double & minz,
	                       double & maxx,double & maxy,double & maxz, int & existantnodes) const
{

	minx = miny = minz =   GRIDBOUNDS * 2;
	maxx = maxy = maxz = - GRIDBOUNDS * 2;
	existantnodes = 0;

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].polymer)
		{
			if (minx > node[i].x) minx = node[i].x;
			if (miny > node[i].y) miny = node[i].y;
			if (minz > node[i].z) minz = node[i].z;

			if (maxx < node[i].x) maxx = node[i].x;
			if (maxy < node[i].y) maxy = node[i].y;
			if (maxz < node[i].z) maxz = node[i].z;

			existantnodes++;
		}
	}

}

int actin::saveinfo()
{
	//char filename[1024];
 
	//char time[1024], date[1024];

    //_strtime( time );
    //_strdate( date );

	//sprintf ( filename , "nodes%05i.txt", filenum );

	// write header
	//opactindata << "Data Output " << date << " " << time << endl;
	//opruninfo << "Nodes Data Itteration Number: " << filenum << endl;
	opruninfo << "Highest Node Count: " << highestnodecount << endl;
	cout << "Highest Node Count: " << highestnodecount << endl;
	//opruninfo << "x,y,z,polymer " << endl;
	
	// find grid range:


    double minx, miny, minz;
	double maxx, maxy, maxz;
	int existantnodes;

    getnodeextents(minx, miny, minz,maxx, maxy, maxz, existantnodes);

    // write info

	opruninfo << "Final existant nodes: " << existantnodes << endl;
	cout << "Final existant nodes: " << existantnodes << endl;

    opruninfo << setprecision(1);
    opruninfo << "x bounds: " << setw(5) << minx << " to " << setw(5) << maxx << endl;
    opruninfo << "y bounds: " << setw(5) << miny << " to " << setw(5) << maxy << endl;
    opruninfo << "z bounds: " << setw(5) << minz << " to " << setw(5) << maxz << endl;
	
    cout << "Gridrange: +/- " << GRIDBOUNDS << " uM" << endl;
    
    cout << setprecision(1);
    cout << "x bounds: " << setw(5) << minx << " to " << setw(5) << maxx << endl;
    cout << "y bounds: " << setw(5) << miny << " to " << setw(5) << maxy << endl;
    cout << "z bounds: " << setw(5) << minz << " to " << setw(5) << maxz << endl;

	return 0;
}

int actin::savevrml(int filenum)
{
	char filename[1024];
	//char time[1024], date[1024];

    //_strtime( time );
    //_strdate( date );

	sprintf ( filename , "%snodes%05i.wrl", TEMPDIR,filenum );

	ofstream opvrml(filename, ios::out | ios::trunc);
	if (!opvrml) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	// write header
	opvrml << "#VRML V2.0 utf8" << endl;
	opvrml << "#comet (compiled " << __DATE__ << " " << __TIME__ << ") Mark J Dayel" << endl;
	//opactindata << "Data Output " << date << " " << time << endl;
	opvrml << "#Nodes Data Itteration Number: " << filenum << endl;
	opvrml << "#Highest Node Count: " << highestnodecount << endl;

	opvrml << "Shape {" << endl;
	opvrml << "    geometry PointSet {" << endl;
	opvrml << "        coord DEF nucleator Coordinate { " << endl;
	opvrml << "            point [ ";

	p_nuc->savevrml(&opvrml);

	opvrml << "] }" << endl;
	opvrml << "    }" << endl;
	opvrml << "}" << endl;

	opvrml << "Shape {" << endl;
	opvrml << "    geometry PointSet {" << endl;
	opvrml << "        coord DEF actin Coordinate { " << endl;
	opvrml << "            point [ ";

	// write data
	for (int i=0; i<highestnodecount; i++)
	{
		if ((node[i].polymer) &&			// don't plot depolymerised and
			(!node[i].listoflinks.empty()))	// only plot nodes with links
		{
			if (i==0)
				opvrml << node[i].x << "," << node[i].y << "," << node[i].z;
			else
				opvrml << "," << node[i].x << "," << node[i].y << "," << node[i].z;
		}
	}

	opvrml << "] }" << endl;
	opvrml << "        color Color { color [ ";

	for (int i=0; i<highestnodecount; i++)
	{
			if (i==0)
				opvrml	<< 1 << " " << 1 << " " << 1;   // RGB
			else
				opvrml	<< "," << 1 << " " << 1 << " " << 1;   // RGB
	}



	opvrml << "] }" << endl;
	opvrml << "    }" << endl;
	opvrml << "}" << endl;

	opvrml.close();



	return 0;
}


void actin::iterate()  // this is the main iteration loop call
{
	if (USE_THREADS)   // only use threads when we have enough nodes
	{
		if (highestnodecount > 0)
			currentlyusingthreads = true;
		else
			currentlyusingthreads = false;
	}
	
    if (!TEST_SQUASH)
    {   // only add new nodes if not testing forces
    	
	    // crosslink, but only after time to equilibrate position
	    // in mean time it is a 'harbinger'
	    // affected only by node collision repulsion and nucleator repulsion

        // add new nodes and store number created in rolling array
	    crosslinknodesdelay[iteration_num % CROSSLINKDELAY] = p_nuc->addnodes();

        // crosslink nodes after they are CROSSLINKDELAY iterations old:
	    crosslinknewnodes(crosslinknodesdelay[(iteration_num + 1) % CROSSLINKDELAY]);
    }
    else
    {
        if (!test_equilibrating)
        {   // this is for pushing on the network to test how strong it is

            testforces_select_nodes(testsurfaceposn,1);
            // if testing forces apply the force here

            testforces_addforces(1);

            if ((((lasttestsurfaceposn - testsurfaceposn) < TEST_DIST_EQUIL) &&
                  (lasttestsurfaceposn - testsurfaceposn) > 0 ))
            {   // has come to equilibrium
                // so save the point

                if (fabs(testsurfaceposn - lasttestsurfacesavedposn) > TEST_DIST_EQUIL)
                    testforces_saveiter();  // only save point if we've actually moved
                
                testforcemag += TEST_FORCE_INCREMENT;
            }

            lasttestsurfaceposn = testsurfaceposn;  // store the old position
        }
        else
        {
            if (sum_delta_movements() < 0.0005)
            {
                cout << endl << "Equilibrium reached, starting test" << endl;

                test_equilibrating = false;
                testforcemag = TEST_FORCE_INITIAL_MAG;

                testforces_set_initial_surfaces();
            }

        }
    }

    if (highestnodecount == lowestnodetoupdate)
    {
        iteration_num++;
        return;  // nothing to do
    }


	if ((lastsorthighestnode != highestnodecount)) // && ((iteration_num % 10) == 0))
	{    // only do full sort periodically, if we have new nodes
		sortnodesbygridpoint(); // sort the nodes so they can be divided sanely between threads
		lastsorthighestnode = highestnodecount;
	} 


    // by far the most time (>80%) is spent in collisiondetection() and a some (~10%) in linkforces()
    // so these functions are multithreaded

	// in multithreaded mode, these two functions just start the threads
    // otherwise they do the actual work
    // these functions depend on node positions 
    // so we don't update them here, but put store forces for later

	collisiondetection();	    // calc node-to-node repulsion
	linkforces();			    // and link forces

    if  (currentlyusingthreads && (USETHREAD_COLLISION || USETHREAD_LINKFORCES))
    {   
        // wait for threads to finish before updating node positions etc.
        thread_queue.complete_queued_tasks();
    }
    

    if (USE_BROWNIAN_FORCES)
        addbrownianforces();


    // applyforces() moves the nodes and nucleator and updates the grid
    // also calls setunitvec() for the nodes, which notes the last node position (for friction)
    applyforces();              
	
    if (!TEST_SQUASH)
        nucleator_node_interactions();	    // do forcable node ejection	
	
    if (COVERSLIPGAP > 0)   // skip if less than diameter (so set to 0 to disable) 
	    squash(COVERSLIPGAP);	 

	iteration_num++;

    //cout << "Iternum " << iteration_num << "   nodecount " << highestnodecount << endl;
    //cout.flush();

	return;
}

void actin::testforces_setup()
{   
    testdirection = vect(0,0,-1);
    double testangledeg = 60;
    testangle = cos(testangledeg * PI/360);  // note not PI/180 since both directions

    const int maxtestsurfaces = 3;

    testnodes.resize(maxtestsurfaces+1);

    // surface 0 is all the nodes in the test zone

    for (int i = 1; i != maxtestsurfaces+1; ++i)
    {
        testnodes[i].resize(0);
        testnodes[i].reserve(1024);
        
    }

    testforces_select_nodes(0,0);  // select *all* nodes within test section, and put into 'surface 0'
    testforces_cutlinks(); // cut them out from main network
    testforces_remove_nontest_nodes(); // remove the non-test nodes

    test_equilibrating = true;

}

void actin::testforces_set_initial_surfaces()
{

    // find the furthest node out
    testsurfaceposn = 0;
    lasttestsurfaceposn = DBL_MAX;
    nodes* p_furthestnode = NULL;

    for(vector <testnodeinfo>::iterator	i_test  = testnodes[0].begin(); 
							            i_test != testnodes[0].end();
					                  ++i_test)
    {   
        if (i_test->nodeptr->dist_from_surface > testsurfaceposn)
        {
            p_furthestnode = i_test->nodeptr;
            testsurfaceposn = p_furthestnode->dist_from_surface;
        }
    }

    // start with testsurfaceposn equal to the furthest out, 
    // and only the furthest node in testnodes

    
    if (p_furthestnode) // if there are any nodes in the test region
    {    
        testnodes[1].push_back(testnodeinfo(p_furthestnode,(*p_furthestnode).unitvec()));
        p_furthestnode->testsurface = 1;

        cout << endl << "Starting surface at " << testsurfaceposn << endl;
    }
    else
    {
        cout << endl << "Error: no nodes in test region" << endl;
        exit(1);
    }

    ofstream optest(TESTNODESFILE, ios::out | ios::trunc); // wipe the nodes test file
    optest << "Itteration Testnodes Force Position" << endl;
    optest.close();

}


void actin::testforces_select_nodes(const double& testdist, const short int &setsurface)
{
    const size_t oldnodes = testnodes[setsurface].size(); 

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        if ((i_node->testnode) && (i_node->testsurface == setsurface))
            continue;   // skip those that are already in this set
        
        if (i_node->dist_from_surface < testdist)
            continue;   // skip those with dist greater than testdist

        // make 2d versions of unit vector

        vect tempunitvector = i_node->unit_vec_posn;
        sym_break_rotation_to_xy_plane.rotate(tempunitvector);

        vect unitvecno_x = tempunitvector;
        unitvecno_x.x = 0;
        unitvecno_x = unitvecno_x.unitvec();

        vect unitvecno_y = tempunitvector;
        unitvecno_y.y = 0;
        unitvecno_y = unitvecno_y.unitvec();

        //if (i_node->unit_vec_posn.dot(testdirection) > testangle)
        if ((unitvecno_x.dot(testdirection) > testangle) &&
            (unitvecno_y.dot(testdirection) > testangle))
        {   
            // store the actual position of the node 
            // so we can tell where on the surface to keep the node
            testnodes[setsurface].push_back(testnodeinfo(&*i_node,(*i_node).unitvec()));
            
            // flag the node so don't add twice and so bitmap plots different colour
            i_node->testnode = true;  
            i_node->testsurface = setsurface;
        }

    }

    if (testnodes[setsurface].size() > oldnodes)  // do a save point if we have more nodes now
        testforces_saveiter();

}

void actin::testforces_cutlinks()
{   /// cut the test region out of the main network

    // this is separate from testforces_remove_nontest_nodes() in case
    // we want to cut the chunk out and see what happens to the rest of the network

    for(vector <nodes>::iterator    i_node  = node.begin(); 
								    i_node != node.begin() + highestnodecount;
						          ++i_node)
    {
        for (vector <links>::reverse_iterator i_link  = i_node->listoflinks.rbegin();
                                      i_link != i_node->listoflinks.rend();
                                    ++i_link )
		{
            if ( ( ( i_link->linkednodeptr->testnode) && (!i_node->testnode) ) ||
                 ( (!i_link->linkednodeptr->testnode) && ( i_node->testnode) ))
            {
                i_node->removelink(i_link->linkednodeptr);
            }
        }

    }
    
}

void actin::testforces_remove_nontest_nodes()
{
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        if(!i_node->testnode)
            i_node->depolymerize();
    }

}


void actin::testforces_addforces(const int &surface)
{

    vect tomove;

    vect lever_arm;

    vect forceontestsurface; //, torqueonsurface;

    forceontestsurface.zero();
    //torqueonsurface.zero();

    for(vector <testnodeinfo>::iterator	i_test  = testnodes[surface].begin(); 
							            i_test != testnodes[surface].end();
					                  ++i_test)
    {   
        // find the displacement to move the node to bring it to 
        // the testsurfaceposn distance

        tomove = (i_test->origunitvec * (RADIUS + testsurfaceposn)) -  // original position of node projected onto test surface
                *(i_test->nodeptr);  // the current node position

        forceontestsurface -= tomove;  // add movement to forceontestsurface (convert to force later)

        //lever_arm = *(i_test->nodeptr) - (testdirection * (RADIUS + testsurfaceposn));

        //torqueonsurface += lever_arm.cross(tomove);

        // move the node to the surface
         *(i_test->nodeptr) += tomove;

        i_test->nodeptr->setunitvec();
        i_test->nodeptr->updategrid();

    }

    forceontestsurface *= NODE_DIST_TO_FORCE;  // convert distance to force
                                                      

    //move the test surface

    testsurfaceposn += (-forceontestsurface.dot(testdirection) - testforcemag) 
                            * NODE_FORCE_TO_DIST / (double) testnodes[1].size();
                                                      // scale by # nodes on the surface,
                                                      // since moving the surface moves this many nodes
                                                      // so the inertia is proportional to node number

    //testsurfacerotation += torqueonsurface / testsurfaceinertia;


}

void actin::testforces_saveiter()
{
    if (testnodes[1].size() == 0)
        return;


    ofstream optest(TESTNODESFILE, ofstream::out | ofstream::app); // add to end

    if (!optest) 
        { cout << "Unable to open file " << TESTNODESFILE << " for output"; }
    else
    {
        cout << endl << setw(6) << setprecision(5)
            << iteration_num << " "
            << testnodes[1].size() << " "
            << testforcemag << " "
            << testsurfaceposn << endl;

        optest 
            << iteration_num << " "
            << testnodes[1].size() << " "
            << testforcemag << " "
            << testsurfaceposn << endl;

        optest.close();
    }

    lasttestsurfacesavedposn = testsurfaceposn;
    
}

double actin::sum_delta_movements()
{
    double delta_sum = 0.0;

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        if (i_node->polymer)
            delta_sum += i_node->delta.length();
    }

    return delta_sum;
}


bool actin::addlinks(nodes& linknode1, nodes& linknode2) const
{
	// crosslink a new node
	// returns 1 if link added, 0 if not

	if (&linknode1 == &linknode2) return false; // can't link to self

	if ((!linknode1.polymer) ||  // make sure both are polymers
		(!linknode2.polymer))  
		return false;

	// make sure the link doesn't go *through* the nucleator
	// i.e. that dot product of normal and link is not < - 0.1
    // actually not 0 else gets rid of *all* the surface links

	if (linknode1.unit_vec_posn.dot( linknode2 - linknode1 ) < - 0.1)
	   return false;

    // find distance

	double dist = calcdist(linknode1,linknode2);

	// calculate probability whether to link or not
	
    double pxlink;
    
    if (VARY_P_XLINK)
        pxlink = P_XLINK * (1 - dist / XLINK_NODE_RANGE);
    else
        pxlink = P_XLINK;


	if ( prob_to_bool(pxlink) )
    {
		linknode1.addlink(linknode2,dist); 
		linknode2.addlink(linknode1,dist);
        
		//cout << "linked " << node[linknode1].nodenum << " to " << node[linknode2].nodenum  << endl;
	
		return true;
	}

	return false;
}

void actin::nucleator_node_interactions()
{

	vect disp,forcevec;
	double dist, force;

#ifndef NO_CALC_STATS
    vect energyvec;
#endif

    for(vector <nodes>::iterator	i_node  = node.begin()+lowestnodetoupdate; 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
	{
        if (!i_node->polymer)
            continue;

        if (i_node->dist_from_surface < 0) // if node is inside nucleator           
        {   
            if (p_nuc->collision(*i_node)) // (*i)->x,(*i)->y,(*i)->z)==0)  
			    i_node->updategrid();  // ejected OK
		    //else
			//    i_node->depolymerize();  // not ejected OK, depolymerize?  (no, this is worse than not doing it.---just report it)
        }

		// do node-nucleator links:

        if (i_node->stucktonucleator)
		{
       
            vect nodepos = *i_node;      // get the node position in world coords  
            world_to_nuc_frame(nodepos); // change to nucleator frame

            disp = nodepos - i_node->nucleator_stuck_position; // this vector is in the nucleator frame
	        dist = disp.length();

            force = NUC_LINK_FORCE * dist;
    	    	        
            if ((force > NUC_LINK_BREAKAGE_FORCE) || (i_node->nodenum < lowestnodetoupdate))
            {
                //node[n].nucleator_link_force.zero();
                i_node->stucktonucleator = false;   // no longer stuck
            }
            else
            {
                // this is in the nucleator frame

                forcevec = disp * (force/dist);  // '(disp/dist)' is just to get the unit vector
                                                 // note this vector in nucleator frame!
                
                vect tomove = forcevec * -NODE_FORCE_TO_DIST  ;  // convert force to distance

				p_nuc->move_nuc(nodepos,tomove);		// add to nucleator movement vector in the nucleator frame



                nuc_to_world_rot.rotate(forcevec);      // rotate to world co-ords before adding to the node vector
                nuc_to_world_rot.rotate(tomove);

                // these are in the world frame

#ifndef NO_CALC_STATS

                // add to node segment stats (tension only since default dist is zero)
		        i_node->adddirectionalmags(forcevec, i_node->linkforce_radial, 
                                                     i_node->linkforce_transverse);

                //energyvec = disp * ( dist * NUC_LINK_FORCE / 2.0);

                //i_node->adddirectionalmags(energyvec, i_node->linkenergy_radial, 
                //                                      i_node->linkenergy_transverse);

#endif

				i_node->nucleator_link_force += tomove;	// add to nuc link force stats

                i_node->link_force_vec -= forcevec;			// add to move node
            }
            
        }

	}

}

void actin::collisiondetection(void)
{  

    if (highestnodecount == lowestnodetoupdate)
        return; // nothing to do

    fill(donenode.begin(), donenode.begin()+highestnodecount, false);
    
    // do collision detection

    if(currentlyusingthreads && USETHREAD_COLLISION)
    {
	    for (int i = 0; i < NUM_THREAD_DATA_CHUNKS; i++)
	    {
            // threads use the nodes_by_thread array
            // so don't need to pass work here, just the threadnum
	        collision_thread_data_array[i].threadnum = i;

            thread_queue.queue_task(&collisiondetectiondowork, &collision_thread_data_array[i]);
    
	    }

	//thread_queue.complete_current_tasks();
    } 
    else 
    {
        // if not using threads, do in one go:
	    collision_thread_data_array[0].threadnum = 0;
	    
        collisiondetectiondowork(&collision_thread_data_array[0]);//, NULL);

    }
    
    return;
}



size_t actin::findnearbynodes(const nodes& ournode, const int& adjgridpoints, const int& threadnum)
{
    // create list of nodes on the same gridpoint, and on adjacent gridpoints (including same GP)
    // 	
    // N.B. a good fraction of the total CPU time is spent in this function
    // may be worth linearizing the loop or poss amalgamating with the calling
    // functions so that creating and reading the list of gridpoints is not necessary
        
   
// save repeatedly dereferencing pointer by threadnum in inner loops:
// (maybe compiler does this anyway?)  prolly does, but loop is slow, so no harm making sure

    Nodes1d * const p_recti_near_nodeslist = &recti_near_nodes[threadnum];
    Nodes1d * const p_nodes_on_same_gridpoint = &nodes_on_same_gridpoint[threadnum];

    
    p_recti_near_nodeslist->resize(0);
    p_nodes_on_same_gridpoint->resize(0);
    
	//const int gridx = ournode.gridx;
    //const int gridy = ournode.gridy;
    //const int gridz = ournode.gridz;
    
    //int minx,miny,minz,maxx,maxy,maxz;
    
    const int minx = ournode.gridx - adjgridpoints;
    const int miny = ournode.gridy - adjgridpoints;
    const int minz = ournode.gridz - adjgridpoints;
    
    const int maxx = ournode.gridx + adjgridpoints + 1;
    const int maxy = ournode.gridy + adjgridpoints + 1;
    const int maxz = ournode.gridz + adjgridpoints + 1;



    p_nodes_on_same_gridpoint->insert(p_nodes_on_same_gridpoint->end(),
                                      NODEGRID(ournode.gridx,ournode.gridy,ournode.gridz).begin(),
                                      NODEGRID(ournode.gridx,ournode.gridy,ournode.gridz).end());

   
	int x,y,z;

	for (x = minx; x != maxx; ++x) 
    {
	    //ny = nodegrid[x];

	    for (y = miny; y != maxy; ++y) 
            {
				//nz = ny[y];	  

	            for (z = minz; z != maxz; ++z) 
                {
                    p_recti_near_nodeslist->insert(p_recti_near_nodeslist->end(),
                                      NODEGRID(x,y,z).begin(),
                                      NODEGRID(x,y,z).end());
            	    //nodeptr = nz[z]; //nodegrid[x][y][z];
	            }
	        }
    }

    return p_recti_near_nodeslist->size();
}

inline double actin::collisiondetection_getrepforce(const double & dist)
{
    return 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE / dist, NODE_REPULSIVE_POWER ) - 1 );
}

void actin::collisiondetection_setrepforcelookup()
{   	
    
    repulsiveforcelookup.resize(REPULSIVE_LOOKUP_DIVISIONS+1);

    for (int i = 0; i < REPULSIVE_LOOKUP_DIVISIONS+1; ++i)
    {
        double dist = MIN_REPULSION_DIST + i * ( (NODE_REPULSIVE_RANGE - MIN_REPULSION_DIST) / REPULSIVE_LOOKUP_DIVISIONS );
        repulsiveforcelookup[i]=collisiondetection_getrepforce( dist );

        //cout << i << "    " << dist << "   " << repulsiveforcelookup[i] << "   " << endl; 

    }

    maxrepulsivemag = collisiondetection_getrepforce( MIN_REPULSION_DIST );

    // // check
    //for (int i = 0; i < REPULSIVE_LOOKUP_DIVISIONS+1; ++i)
    //{
    //    double dist = MIN_REPULSION_DIST + i * ( (NODE_REPULSIVE_RANGE - MIN_REPULSION_DIST) / REPULSIVE_LOOKUP_DIVISIONS );
    //    cout << i << "    " << dist << "   " << repulsiveforcelookup[i] << "   " << collisiondetection_repforcelookup(dist) << endl; 
    //}

}




void * actin::collisiondetectiondowork(void* arg)//, pthread_mutex_t *mutex)
{	
    // cast arg
    const thread_data* const dat = (thread_data*) arg;

    if (nodes_by_thread[dat->threadnum].size() == 0)
        return NULL; // nothing to do

	const double local_NODE_REPULSIVE_RANGEsqared = NODE_REPULSIVE_RANGE * NODE_REPULSIVE_RANGE;
//    const double local_NODE_REPULSIVE_MAGscaled = 0.06 * NODE_REPULSIVE_MAG;

    vect nodeposvec;
    double distsqr;
    //double dist;
    double recipdist;
    //double recip_dist;
	double rep_force_mag;
    vect rep_force_vect;
    vect disp;

#ifndef NO_CALC_STATS
    vect repenergy_vect;

    // this is so that energy = 0 for dist = NODE_REPULSIVE_RANGE (greater than that is 0 since we don't do the calc'n)
//    const double energyoffset = - 0.06 * NODE_REPULSIVE_MAG * NODE_REPULSIVE_RANGE * ( (1.0 / ( 1 - NODE_REPULSIVE_POWER )) -1  );

#endif
    int x,y,z;
    //int minx,miny,minz;
    //int maxx,maxy,maxz;

#ifdef PROXIMITY_VISCOSITY
	double viscfactor;
#endif

    nodes *p_i_node, *p_nearnode, *p_sameGPnode;


    int sameGPnodenum;

    double dist;

    // loop over all the nodes given to this thread
    for(vector <nodes*>::iterator	i_node  = nodes_by_thread[dat->threadnum].begin(); 
									i_node != nodes_by_thread[dat->threadnum].end();
							      ++i_node) 
    {	
        p_i_node = *i_node;

        if(donenode[p_i_node->nodenum])     
	        continue;  // skip nodes already done

        //assert (p_i_node->nodegridptr);

        if(!p_i_node->nodegridptr)
            continue;  // skip if not on grid

        const int minx = p_i_node->gridx - NODE_REPULSIVE_RANGE_GRIDSEARCH;
        const int miny = p_i_node->gridy - NODE_REPULSIVE_RANGE_GRIDSEARCH;
        const int minz = p_i_node->gridz - NODE_REPULSIVE_RANGE_GRIDSEARCH;
        
        const int maxx = p_i_node->gridx + NODE_REPULSIVE_RANGE_GRIDSEARCH + 1;
        const int maxy = p_i_node->gridy + NODE_REPULSIVE_RANGE_GRIDSEARCH + 1;
        const int maxz = p_i_node->gridz + NODE_REPULSIVE_RANGE_GRIDSEARCH + 1;

	    // loop over nodes on same gridpoint
        for(NODEGRIDTYPE <nodes*>::iterator sameGPnode  = p_i_node->nodegridptr->begin(); 
									        sameGPnode != p_i_node->nodegridptr->end();
								          ++sameGPnode) 
        {
            p_sameGPnode = *sameGPnode;  

            sameGPnodenum = p_sameGPnode->nodenum;
    	    
            // do coverslip repulsion if necessary
            if (COVERSLIPGAP > 0)
            {

#ifdef FORCE_REPULSIVE_POWER_TO_TWO

                if (p_sameGPnode->x >   COVERSLIPGAP - NODE_REPULSIVE_RANGE)
                {
                    recipdist = 1.0 / ( p_sameGPnode->x - COVERSLIPGAP );

                    if ( recipdist < 1.0/0.05 )  // 1.0/0.05 (i.e. if less than dist of 0.05)
                       rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( local_NODE_REPULSIVE_RANGEsqared * recipdist * recipdist - 1 );
                    else
                       rep_force_mag = local_NODE_REPULSIVE_MAGscaled * (  NODE_REPULSIVE_RANGE * NODE_REPULSIVE_RANGE * 20 * 20  - 1 ) ;

                    p_sameGPnode->rep_force_vec.x -= 2*rep_force_mag ;
                }

                if (p_sameGPnode->x < - COVERSLIPGAP + NODE_REPULSIVE_RANGE)
                {
                    recipdist = 1.0 / ( - p_sameGPnode->x - COVERSLIPGAP );

                    if ( recipdist < 1.0/0.05 )   // 1.0/0.05 (i.e. if less than dist of 0.05)
                       rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( NODE_REPULSIVE_RANGE * recipdist * NODE_REPULSIVE_RANGE * recipdist - 1 );
                    else
                       rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( NODE_REPULSIVE_RANGE * NODE_REPULSIVE_RANGE * 20 * 20 - 1 ) ;

                    p_sameGPnode->rep_force_vec.x += 2*rep_force_mag ;
                }


#else

                if (p_sameGPnode->x >   COVERSLIPGAP - NODE_REPULSIVE_RANGE)
                {
                    recipdist = 1.0 / ( p_sameGPnode->x - COVERSLIPGAP );

                    if ( recipdist < 1.0/0.05 )  // 1.0/0.05 (i.e. if less than dist of 0.05)
                       rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) - 1 );
                    else
                       rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * 20 , NODE_REPULSIVE_POWER ) - 1 ) ;

                    p_sameGPnode->rep_force_vec.x -= 2*rep_force_mag ;
                }

                if (p_sameGPnode->x < - COVERSLIPGAP + NODE_REPULSIVE_RANGE)
                {
                    recipdist = 1.0 / ( - p_sameGPnode->x - COVERSLIPGAP );

                    if ( recipdist < 1.0/0.05 )   // 1.0/0.05 (i.e. if less than dist of 0.05)
                       rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) - 1 );
                    else
                       rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * 20 , NODE_REPULSIVE_POWER ) - 1 ) ;

                    p_sameGPnode->rep_force_vec.x += 2*rep_force_mag ;
                }


                if (p_sameGPnode->x >  COVERSLIPGAP - NODE_REPULSIVE_RANGE)
                {
                    dist = p_sameGPnode->x - COVERSLIPGAP ;

                    if (dist < MIN_REPULSION_DIST )
                        rep_force_mag = maxrepulsivemag;
                    else
                    {
#ifdef USE_REPULSIVE_LOOKUP
                                   rep_force_mag = collisiondetection_repforcelookup(dist); 
#else
                                   rep_force_mag = collisiondetection_getrepforce(dist);
#endif
                    }

                    p_sameGPnode->rep_force_vec.x -= 2*rep_force_mag ;
                }

                if (p_sameGPnode->x < - COVERSLIPGAP + NODE_REPULSIVE_RANGE)
                {
                    dist = - p_sameGPnode->x - COVERSLIPGAP ;

                    if (dist < MIN_REPULSION_DIST) 
                        rep_force_mag = maxrepulsivemag;
                    else
                    {
#ifdef USE_REPULSIVE_LOOKUP
                                   rep_force_mag = collisiondetection_repforcelookup(dist); 
#else
                                   rep_force_mag = collisiondetection_getrepforce(dist);
#endif
                    }

                    p_sameGPnode->rep_force_vec.x += 2*rep_force_mag ;
                }


#endif


            }

	        if (donenode[sameGPnodenum])
		      continue;  // skip if done
    	    
			donenode[sameGPnodenum] = true;  // mark node as done
    		
			// of these nodes, calculate euclidian dist
			nodeposvec = *p_sameGPnode;	// get xyz of our node

		

            // loop over adjacent gridpoints
	        for (x = minx; x != maxx; ++x)
            {
	            for (y = miny; y != maxy; ++y)
                {

	                for (z = minz; z != maxz; ++z) 
                    {
                        // this loop is where all the computational time is taken...


						const NODEGRIDTYPE <nodes*>::iterator nearnode_begin  = NODEGRID(x,y,z).begin();
						const NODEGRIDTYPE <nodes*>::iterator nearnode_end    = NODEGRID(x,y,z).end();
						
			            for(NODEGRIDTYPE <nodes*>::iterator nearnode  = nearnode_begin; 
									                        nearnode != nearnode_end;
		                                                  ++nearnode) 
			            {
                            p_nearnode = *nearnode;

			                if ( p_sameGPnode == p_nearnode)
				                continue;  // skip if self

                            if ((p_sameGPnode->nodenum < lowestnodetoupdate)  &&
                                (p_nearnode->nodenum < lowestnodetoupdate))
                                continue; // if we're not updating either node, then skip

			                disp = *p_nearnode - nodeposvec;   // find relative position of node
			                
                            
                            distsqr = disp.sqrlength();    // square of distance
                			                               // avoid sqrt at this point, in case it's out of range and we don't need to do the repulsion
               			
			                if (distsqr < local_NODE_REPULSIVE_RANGEsqared)
			                {

				                //dist = SSEsqrt(distsqr);  
                                recipdist = SSErsqrt(distsqr);
                                dist = 1.0/recipdist;
                                // calculate repulsive force
  
                                // this was the old function
					            //rep_force_mag = 13 * NODE_REPULSIVE_MAG * (exp ((-7.0*dist/NODE_REPULSIVE_RANGE)) - exp (-7.0));

                              
#ifdef FORCE_REPULSIVE_POWER_TO_TWO

                                // n.b. if you change this function, also change the energy function below to be the integral
                               // F_R = M_R (( \frac{d_R}{d} )^{P_R} -1)
                               if (recipdist < 1.0/0.05 )   // 1.0/0.05 (i.e. if less than dist of 0.05)
                                   rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( local_NODE_REPULSIVE_RANGEsqared * recipdist * recipdist - 1 );
                                   //rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) - 1 );
                               else
                                   rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( local_NODE_REPULSIVE_RANGEsqared / (0.05 * 0.05) - 1 ) ;
                                   //rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE / 0.05 , NODE_REPULSIVE_POWER ) - 1 ) ;
 
#else                          
                               if (dist < MIN_REPULSION_DIST) 
                                   rep_force_mag = maxrepulsivemag;
                               else
                               {
#ifdef USE_REPULSIVE_LOOKUP
                                   rep_force_mag = collisiondetection_repforcelookup(dist); 
#else
                                   rep_force_mag = collisiondetection_getrepforce(dist);
#endif
                               }

  
                               // if (recipdist < 1.0/0.05 )   // 1.0/0.05 (i.e. if less than dist of 0.05)
                               //    //rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( local_NODE_REPULSIVE_RANGEsqared * recipdist * recipdist - 1 );
                               //    rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) - 1 );
                               //else
                               //    //rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( local_NODE_REPULSIVE_RANGEsqared / (0.05 * 0.05) - 1 ) ;
                               //    rep_force_mag = local_NODE_REPULSIVE_MAGscaled * ( pow( NODE_REPULSIVE_RANGE / 0.05 , NODE_REPULSIVE_POWER ) - 1 ) ;

#endif


						   // if harbinger, ramp up the repulsive forces gradually (linearly with itter no):

					            if (p_sameGPnode->harbinger)
					            {	// harbinger being repelled
						            rep_force_mag *= (double) ( (iteration_num - p_sameGPnode->creation_iter_num) % CROSSLINKDELAY) / (double) CROSSLINKDELAY;
					            }

					            if (p_nearnode->harbinger)
					            {	// being repelled from harbinger
						            rep_force_mag *= (double) ( (iteration_num - p_nearnode->creation_iter_num) % CROSSLINKDELAY) / (double) CROSSLINKDELAY;
            					
						            if (p_sameGPnode->harbinger)
						            {	// two harbingers, so move them
							            p_sameGPnode->move_harbinger_this_time = true;
						            }
					            }

					            // convert force magnitude to vector

       
					            rep_force_vect = disp * rep_force_mag * recipdist;

                                // add to the actual repulsion

				                p_sameGPnode->rep_force_vec -= rep_force_vect ;


#ifndef NO_CALC_STATS

                                // add to the statistics
				                p_sameGPnode->adddirectionalmags(rep_force_vect, p_sameGPnode->repforce_radial,
								                                                 p_sameGPnode->repforce_transverse);

                                // this is the integral of the force function 
//                                repenergy_vect = disp * recipdist * ( 0.06 * SSEsqrt(distsqr) * NODE_REPULSIVE_MAG 
//                                    *  ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) / ( 1- NODE_REPULSIVE_POWER) - 1 ) - energyoffset);

//                                p_sameGPnode->adddirectionalmags(repenergy_vect, p_sameGPnode->repenergy_radial,
//								                                                 p_sameGPnode->repenergy_transverse);



#endif
                                // pressure can be used for polymerizationstalling---not part of the stats 
                                p_sameGPnode->pressure += rep_force_mag;

            #ifdef PROXIMITY_VISCOSITY
                                // do proximity-based viscosity

					            if ((VISCOSITY) && (distsqr < VISC_DIST * VISC_DIST))
					            {
						            viscfactor = mymin(MAX_VISC_WEIGHTING,recipdist);

						            p_sameGPnode->viscosity_velocity_sum += p_nearnode->delta * viscfactor;
						            p_sameGPnode->viscosity_velocity_unweight += viscfactor;
					            }
            #endif 
                            }
                        }
                    }
                }
            }
	    }
    }

    return (void*) NULL;
}



// ApplyForces
void actin::applyforces(void) 
{   /// this applys the previously calculated forces

    /// sets the rotation etc. up then calls threads to do the work

	if (highestnodecount == 0)
		return;

	// set the rotation matrix

	double x_angle = 0.0;
	double y_angle = 0.0;			
	double z_angle = 0.0;

	torque_rotate.settoidentity();

    if (IMPOSED_NUC_DISP)
    {
        double speed = IMPOSED_NUC_DISP_SPEED * DELTA_T;
        p_nuc->deltanucposn = vect( 0, speed, 0);
    }

	if (IMPOSED_NUC_ROT)
	{
		// rotate const speed

		double angle = IMPOSED_NUC_ROT_SPEED * 2 * PI * DELTA_T;

		torque_rotate.rotatematrix( angle, 0, 0);


	}
	else if (ROTATION)
	{
        // torque vector is in nucleator frame, so transform to world frame to make the torque rotation matrix

        vect world_torque_vec = p_nuc->torque;

        nuc_to_world_rot.rotate(world_torque_vec);

		x_angle = world_torque_vec.x / p_nuc->momentofinertia.x;
		y_angle = world_torque_vec.y / p_nuc->momentofinertia.y;
		z_angle = world_torque_vec.z / p_nuc->momentofinertia.z;

		torque_rotate.rotatematrix( x_angle, y_angle, z_angle);    
	}


    if (!TEST_SQUASH)
    {   // only move nucleator if not testing forces

	    if (IMPOSED_NUC_ROT || ROTATION)
	    {
            // rotate the nucleator

            world_to_nuc_rot.rotatematrix(torque_rotate); // rotate nucleator by (torque in world co-ords)

            nuc_to_world_rot = world_to_nuc_rot.inverse();  // set the inverse rotation

            p_nuc->last_torque_rotate = torque_rotate;    // store the torque matrix (used in nucleator friction)
        }

        // deltanucposn is in nucleator frame, so change to world frame and move the nucleator

        nuc_to_world_rot.rotate(p_nuc->deltanucposn);

        // if the nucleator is moving too fast, slow it down
        // this is to simulate the drag of the medium 
        // when the bead pops out in reality it's not in a vacuum
        if (p_nuc->deltanucposn.length() > MAX_NUC_MOVE * DELTA_T)
            p_nuc->deltanucposn *= MAX_NUC_MOVE * DELTA_T / p_nuc->deltanucposn.length();

        // update the nucleator position
        //p_nuc->last_position = p_nuc->position;
        p_nuc->position += p_nuc->deltanucposn;

        p_nuc->deltanucposn_sum +=  p_nuc->deltanucposn;

        p_nuc->last_delta_position = p_nuc->deltanucposn;

        // rotate the nucleator displacement vector
        //torque_rotate.rotate(p_nuc->position);
	    

	    //nuc_disp = p_nuc->deltanucposn; // store nucleator movement in static for threads

    }
    //else
    //{   // TEST_SQUASH is true, so lock the nucleator position
    //    nuc_disp.zero();
    //}

    

	// and zero
	p_nuc->deltanucposn.zero();

	// clear the torque vector now that the rotation matrix is set
	p_nuc->torque.zero();
	

#ifndef SEED_INSIDE

    if (currentlyusingthreads && USETHREAD_APPLYFORCES)	
    {

		// add the ones up to lowestnodetoupdate
        // do this separately since they should go faster---less to do 
        // and we want the threads to finish at the same time
        // note: applyforcesdowork() identifies these by these
        // dat->threadnum being >= NUM_THREAD_DATA_CHUNKS

		if (lowestnodetoupdate > 0)
			addapplyforcesthreads(NUM_THREAD_DATA_CHUNKS, 0, lowestnodetoupdate);

		addapplyforcesthreads(0, lowestnodetoupdate, highestnodecount);
        
        //addapplyforcesthreads(0, highestnodecount);

	    thread_queue.complete_queued_tasks();
    } 
    else 
    {
        if (lowestnodetoupdate > 0)
        {
			addapplyforcesthreads(NUM_THREAD_DATA_CHUNKS, 0, lowestnodetoupdate);

            applyforces_thread_data_array[0].startnode = node.begin();
	        applyforces_thread_data_array[0].endnode = node.begin() + lowestnodetoupdate;
	        applyforces_thread_data_array[0].threadnum = NUM_THREAD_DATA_CHUNKS;

            applyforcesdowork(&applyforces_thread_data_array[0]);
        }

        applyforces_thread_data_array[0].startnode = node.begin() + lowestnodetoupdate;
        applyforces_thread_data_array[0].endnode = node.begin() + highestnodecount;
        applyforces_thread_data_array[0].threadnum = 0;

        applyforcesdowork(&applyforces_thread_data_array[0]);

    }
    
    // note: we can't do it in threads because
    // of the linked list in the nodegrid.
	// not much time is spent here

	for(vector <nodes>::iterator	i_node  = node.begin() + lowestnodetoupdate; 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
	    i_node->updategrid(); // move the point on the grid if need to
    }

#endif
    return;
}

void actin::addapplyforcesthreads(const int threadnumoffset, const int &lowestnodenum, const int &highestnodenum)
{   ///  Splits node range between NUM_THREAD_DATA_CHUNKS and sets applyforces_thread_data_array for each

	if ((highestnodenum - lowestnodenum) <= 0)
		return;

	int numthreadnodes, start, end;

	numthreadnodes = (highestnodenum-lowestnodenum) / NUM_THREAD_DATA_CHUNKS;

	// if there is a remainder, distribute it between the threads
	if (((highestnodenum-lowestnodenum) % NUM_THREAD_DATA_CHUNKS) != 0)
		numthreadnodes += 1;

    for (int i = 0; i < NUM_THREAD_DATA_CHUNKS; i++)
    {
        start = i * numthreadnodes + lowestnodenum;
        end = (i+1) * numthreadnodes + lowestnodenum;
	    
        // if there was a remainder, then last thread will be short, and
        // allocation will overrun end of nodes, so truncate
        if (end > highestnodenum)
            end = highestnodenum;

        if (start >= end)
            continue;

        applyforces_thread_data_array[threadnumoffset + i].startnode = node.begin() + start;
        applyforces_thread_data_array[threadnumoffset + i].endnode = node.begin() + end;
        applyforces_thread_data_array[threadnumoffset + i].threadnum = threadnumoffset + i;
    
        thread_queue.queue_task(&applyforcesdowork, &applyforces_thread_data_array[threadnumoffset + i]);
    }

}


void* actin::applyforcesdowork(void* arg)//, pthread_mutex_t *mutex)
{   /// run concurrently by threads, applies the forces and rotations for the nodes

    // cast arg
    const thread_data* const dat = (thread_data*) arg;

    // if we're in the range above lowestnodetoupdate then
    // dat->threadnum will be < NUM_THREAD_DATA_CHUNKS
    // this is a bit of a cludge---may fix later

    if (dat->threadnum < NUM_THREAD_DATA_CHUNKS)
    {
        for(vector <nodes>::iterator i_node  = dat->startnode;
                                     i_node != dat->endnode;
                                   ++i_node)
	    {
		    if (!i_node->polymer)
			    continue;                       

		    if (i_node->harbinger)
		    {   // is a harbinger, check pressure---do we depolymerise it?
			    if (i_node->pressure > MAX_POLYMERISATION_PRESSURE)
			    {
				    i_node->depolymerize();
				    continue;
			    }

			    // special case: move harbinger if it is repelled by another harbinger, or if flag set
			    if (ALLOW_HARBINGERS_TO_MOVE || i_node->move_harbinger_this_time)
			    {	
				    i_node->applyforces();
				    i_node->move_harbinger_this_time = false;
			    }

                i_node->rep_force_vec.zero();  // clear the repulsive forces every time
                                                 // the possibility of move_harbinger_this_time
                                                 // means that repulsive forces could otherwise build up
		    }
		    else
		    {	// move if not harbinger
			    //torque_rotate.rotate(*i_node);	 // rotate
			    //*i_node -= nuc_disp;	 // move wrt nucleator frame of ref
                i_node->applyforces();
			    //if (i>=lowestnodetoupdate)
				   // i_node->applyforces();	         // move according to forces
		    }
        }
    }
    //else
    //{   // just move/rotate, don't apply forces
    //    for(vector <nodes>::iterator i_node  = dat->startnode;
    //                                 i_node != dat->endnode;
    //                               ++i_node)
	//    {
	//	    torque_rotate.rotate(*i_node);	 // rotate
	//	    *i_node -= nuc_disp;	             // move wrt nucleator frame of ref
    //    }
    //}

    return NULL;
}

void actin::linkforces()
{
   

    // remove the links that were broken last time
    for(int threadnum=0; threadnum<NUM_THREAD_DATA_CHUNKS; ++threadnum)
    {
        for (unsigned int i=0; i!=linkremovefrom[threadnum].size() ; ++i )
        {
            linkremovefrom[threadnum][i]->removelink(linkremoveto[threadnum][i]);  // remove the back link
	        linkremoveto[threadnum][i]->removelink(linkremovefrom[threadnum][i]);  // and remove from the list
        }
        // reset the lists of broken links:
        linkremovefrom[threadnum].resize(0);
        linkremoveto[threadnum].resize(0);
    }
   
    if (highestnodecount == lowestnodetoupdate)
        return;  // nothing to do


    // parcel out the work for the threads

    int numthreadnodes, start, end;


    numthreadnodes = (highestnodecount-lowestnodetoupdate) / NUM_THREAD_DATA_CHUNKS;

    // if there is a remainder, distribute it between the threads
    if (((highestnodecount-lowestnodetoupdate) % NUM_THREAD_DATA_CHUNKS) != 0)
        numthreadnodes += 1;



    if (currentlyusingthreads && USETHREAD_LINKFORCES)
    {
	    for (int i = 0; i < NUM_THREAD_DATA_CHUNKS; i++)
	    {
	        start = i * numthreadnodes + lowestnodetoupdate;
	        end = (i+1) * numthreadnodes + lowestnodetoupdate;
    	    
            // if there was a remainder, then last thread will be short, and
            // allocation will overrun end of nodes, so truncate
            if (end > highestnodecount)
                end = highestnodecount;

            if (start >= end)
                continue; // nothing to do

            //if (i==NUM_THREAD_DATA_CHUNKS-1)
	        //{	// put remainder in last thread (cludge for now)
		       // end += highestnodecount % NUM_THREAD_DATA_CHUNKS;
	        //}

            //cout << "Adding Links job " << start << " to " << end << " in thread " << i << endl;

	        linkforces_thread_data_array[i].startnode = node.begin() + start;
	        linkforces_thread_data_array[i].endnode = node.begin() + end;
	        linkforces_thread_data_array[i].threadnum = i;
    	    
		thread_queue.queue_task(&linkforcesdowork, &linkforces_thread_data_array[i]);
	    }

	    //thread_queue.complete_current_tasks();
    }
    else
    {
        // if not using threads, do in one go:
	    linkforces_thread_data_array[0].startnode = node.begin() + lowestnodetoupdate;
	    linkforces_thread_data_array[0].endnode = node.begin() + highestnodecount;
	    linkforces_thread_data_array[0].threadnum = 0;

	    linkforcesdowork(&linkforces_thread_data_array[0]);//, NULL);
    }

    return;
}

// LinkForces
void * actin::linkforcesdowork(void* arg)//, pthread_mutex_t *mutex)
{
    // cast arg
    const thread_data* const dat = (thread_data*) arg;
 

    if (dat->startnode == dat->endnode)
        return NULL;    // nothing to do

    double dist, force;
    vect nodeposvec, disp;
    vect forcevec;

#ifdef LINK_VISCOSITY
	double viscfactor;
#endif

#ifndef NO_CALC_STATS
//    vect energyvec;
#endif

    assert( &(*(dat->startnode)) != 0);

    // go through all nodes
    for(vector <nodes>::iterator i_node  = dat->startnode;
                                 i_node != dat->endnode;
                               ++i_node)
    {   
        assert( &(*(i_node)) != 0);

        // make sure they have links and are polymers
        if ((i_node->listoflinks.empty() ) || (!i_node->polymer))
	        continue;

        //cout << "Thread " << dat->threadnum << " List size " << (int)i_node->listoflinks.size() << endl;
        //cout.flush();

    	// store the node position for later calculation of link lengths
	    nodeposvec = *i_node;
   	
  		// go through links for each node
		for (vector <links>::iterator i_link  = i_node->listoflinks.begin();
                                      i_link != i_node->listoflinks.end();
                                    ++i_link )
		{	
            //cout << "Link Thread " << dat->threadnum << " " << endl;
            //cout.flush();

            assert( &(*(i_link)) != 0);
            assert( i_link->linkednodeptr != 0);
            assert( !i_link->broken ); // if link not broken  (shouldn't be here if broken anyway)

			disp = nodeposvec - *(i_link->linkednodeptr);
			dist = disp.length();
    	    
            if (i_link->getlinkforces(dist,force)) // returns false if broken
            {
				forcevec = disp * (force/dist); // convert force to vector
        		
				i_node->link_force_vec += forcevec;

        		i_link->forcesum += force;  // this is one of the stats, but calc'n is quick since no trig, so leave in

#ifndef NO_CALC_STATS


//                energyvec = disp * ( LINK_FORCE * i_link->linkforcescalefactor * (2 - (dist * i_link->orig_dist_recip)  )) / 2.0;


                if (force < 0.0) // put tension into link forces
				{
					i_node->adddirectionalmags(forcevec, i_node->linkforce_radial, i_node->linkforce_transverse);

//                    i_node->adddirectionalmags(energyvec, i_node->linkenergy_radial, 
//                                                          i_node->linkenergy_transverse);
				} 
				else 
				{	   // but put compression into repulsive forces
					i_node->adddirectionalmags(forcevec, i_node->repforce_radial , i_node->repforce_transverse );

//                    i_node->adddirectionalmags(energyvec, i_node->repenergy_radial, 
//                                                          i_node->repenergy_transverse);		
				}


                

#endif

                if (force > 0) // put repulsion into pressure
				{
                    i_node->pressure += force;
                }

#ifdef LINK_VISCOSITY
				if (VISCOSITY)
				{
					//vel_sum += node[n].delta;

					viscfactor = mymin(MAX_VISC_WEIGHTING,1/dist);

					node[n].viscosity_velocity_sum += i->linkednodeptr->delta * viscfactor;
					node[n].viscosity_velocity_unweight += viscfactor;
				}
#endif 

            }
            else 
			{
				// broken link: store which ones to break:
				linkremovefrom[dat->threadnum].push_back(&(*i_node));
				linkremoveto[dat->threadnum].push_back(i_link->linkednodeptr);
			} 

        }
		

    }

    return NULL;
}

//bool comparenodenum(tracknodeinfo & value, int & nodenum)
//{
//    return (value.nodenum == nodenum);
//}

void actin::set_axisrotation(const projection &  proj, rotationmatrix & axisrotation)
{
    axisrotation.settoidentity();

    if (proj == xaxis)
	{        
        return;
    } else
    if (proj == yaxis)
	{	

        axisrotation.rotatematrix(PI/2, zaxis);  // x to y is by rotation around z
	}
	else
	{	
        axisrotation.rotatematrix(-PI/2, xaxis);
        axisrotation.rotatematrix(PI/2, yaxis); // x to z is by rotation around y and z
        axisrotation.rotatematrix(PI/2, zaxis);
	}
}



void actin::set_nodes_to_track(const projection & proj)
{

    nodes_to_track.resize(0);

    vector <int> temp_nodes_to_track;
    temp_nodes_to_track.resize(0);

    rotationmatrix axisrotation;
    set_axisrotation(proj, axisrotation);

    rotationmatrix projection_rotation;
    projection_rotation.settoidentity();

    // projection_rotation gets us from the bead 
    // frame-of-ref to the bitmap projection frame-of-ref
    projection_rotation.rotatematrix(axisrotation);           // rotates for the projection axis 
    projection_rotation.rotatematrix(sym_break_rotation_to_zaxis);       // rotates for the symmetry breaking direction
                                                              // using sym_break_rotation_to_zaxis so we can tell the furthest node back
                                                              // from the sym break direction
    
    //if (!BMP_FIX_BEAD_ROTATION)
    //    projection_rotation.rotatematrix(nuc_to_world_rot); // compensates for bead rotation

    vect nucposn = ptheactin->p_nuc->position;
    //projection_rotation.rotate(nucposn);

    cout << "Selecting nodes from frame " << TRACK_MIN_RANGE << " to " << TRACK_MAX_RANGE << endl;

    vect tempposn;

    stationary_node_number = 0;    // node to lock the bitmap to (opposite the sym break direction)
    double furthest_node_posn = 0;
    double restrictrange = DBL_MAX;

    

    if (!POST_VTK) // only restrict if not plotting vtk tracks
    {
        restrictrange = RADIUS / 2;
        //MAX_NODES_TO_TRACK = 30;
    }
                      

    // first we find all the nodes within the frame range and put them into temp_nodes_to_track


    for(int i=0; i != highestnodecount; ++i)
    {       
        //if (temp_nodes_to_track.size() == MAX_NODES_TO_TRACK)   // this should be in later loop (if we do the redistribute loop)
        //        break;         // limit to 20 nodes
        
        if (node[i].dist_from_surface < 0.1)  // skip ones in contact with surface
            continue;

        if ((node[i].creation_iter_num < TRACK_MIN_RANGE * InterRecordIterations) || 
            (node[i].creation_iter_num > TRACK_MAX_RANGE * InterRecordIterations))  // is within range?
            continue;

        tempposn = node[i] - nucposn;
        projection_rotation.rotate(tempposn);

        if (fabs(tempposn.x) < restrictrange)  // if within range
        {

            temp_nodes_to_track.push_back(i);
            if (furthest_node_posn > tempposn.z)
            {
                furthest_node_posn = tempposn.z;
                stationary_node_number = i;
             }

        }

    }

    cout << "There are " << (int) temp_nodes_to_track.size() << " nodes with the frame range given" << endl;

    // start with the middle node furthest back:

    nodes_to_track.push_back(stationary_node_number);

    //vect lastpos = node[stationary_node_number];
    //projection_rotation.rotate(lastpos);



    // now we select a subset of these nodes distributed as evenly as we can over the sphere

    for(unsigned int i = 0; i != temp_nodes_to_track.size(); ++i)
    {   // this number of nodes to track
        
        if (i == MAX_NODES_TO_TRACK)
            break;         // limit to certain number of nodes


//        int max_dist_posn = 0;
       
        vector <int> nodeslist;
        vector <double> closestdist;

        //nodeslist.resize(temp_nodes_to_track.size());
        closestdist.resize(temp_nodes_to_track.size());

        for(unsigned int j = 0; j != temp_nodes_to_track.size(); ++j)
        {  
            // make sure node not already in list
            if ( find(nodes_to_track.begin(), nodes_to_track.end(), temp_nodes_to_track[j]) != nodes_to_track.end() )
                continue;

            vect curpos = node[temp_nodes_to_track[j]];
            projection_rotation.rotate(curpos);

            // find how far (dist squared) from all the current nodes

//            double tot_dist_sq = 0.0;

            double closest_dist = DBL_MAX;
            //vector<int>::iterator closest = nodes_to_track.begin(); 

            for(vector<int>::iterator i_nodenum2  = nodes_to_track.begin(); 
			                          i_nodenum2 != nodes_to_track.end();
			                        ++i_nodenum2)
            {

                if ( *i_nodenum2 == temp_nodes_to_track[j])
                    continue;   // make sure not self

                vect topos = node[*i_nodenum2];
                projection_rotation.rotate(topos);

                double dist = (curpos-topos).length();

                if (dist < closest_dist)
                {
                    closest_dist = dist;
                    //closest = i_nodenum2;
                }

            }

            closestdist[j] = closest_dist;   // distance
            //nodeslist[j] = *closest;      // nodenum
            
        }

        // find node with the highest closest distance (i.e. the furthest from the others)

        double high_closest_dist = 0;
        int   high_closest_dist_inx = 0;

        for(unsigned int j = 0; j != temp_nodes_to_track.size(); ++j)
        {
            if (closestdist[j] > high_closest_dist)
            {
                high_closest_dist = closestdist[j];
                high_closest_dist_inx = j;
            }

        }
         
        // add the node to the list

        nodes_to_track.push_back(temp_nodes_to_track[high_closest_dist_inx]);

        //cout << "Adding node " << high_closest_dist_inx << " ";

    }




    // Now we have the nodes to track within the main selection shell
    // optionally we also want to track corresponding nodes within the original shell (frames 10--30)
    // so we go through these nodes and find the nearest node (by radial angle) to the ones already selected

    // fixme: this is a clunky copy-paste-edit from above

    

    
//    double TRACK_MIN_RANGE2 = 10;
//    double TRACK_MAX_RANGE2 = 30;

//    double TRACK_TARGET_DIST = RADIUS * 1.5;

    if  (TRACKS_LENGTHS)
    {

     // first we find all the nodes within the frame range and put them into temp_nodes_to_track

        if (SECOND_SHELL)
        {   // overwrite the temp nodes array with the second shell if we're using it

            temp_nodes_to_track.resize(0);

            for(int i=0; i != highestnodecount; ++i)
            {       
                //if (temp_nodes_to_track.size() == MAX_NODES_TO_TRACK)   // this should be in later loop (if we do the redistribute loop)
                //        break;         // limit to 20 nodes
                
                if (node[i].dist_from_surface < 0.1)  // skip ones in contact with surface
                    continue;

                if ((node[i].creation_iter_num < TRACK_MIN_RANGE2 * InterRecordIterations) || 
                    (node[i].creation_iter_num > TRACK_MAX_RANGE2 * InterRecordIterations))  // is within range?
                    continue;

                tempposn = node[i] - nucposn;
                projection_rotation.rotate(tempposn);

                if (fabs(tempposn.x) < restrictrange)  // if within range
                {

                    temp_nodes_to_track.push_back(i);                                         
                }

            }

            cout << "There are " << (int) temp_nodes_to_track.size() << " nodes with the second frame range given" << endl;


        }

        size_t num_nodes = nodes_to_track.size(); // how many nodes have we got to find pairs for?



        // now we select the closest one in this shell to each of the first selected nodes

        for(unsigned int i = 0; i != num_nodes; ++i)
        {   // this number of nodes to track
            
            vect curpos = node[nodes_to_track[i]];
            projection_rotation.rotate(curpos);


            double closest_dist = DBL_MAX;
            double best_score = DBL_MAX;
            double score = 0;
            unsigned int best_nodenum = 0;

            // go through temp list, looking for best match

            for(unsigned int j = 0; j != temp_nodes_to_track.size(); ++j)
            {  
                // make sure node not already in list
                if ( find(nodes_to_track.begin(), nodes_to_track.end(), temp_nodes_to_track[j]) != nodes_to_track.end() )
                    continue;

                vect topos = node[temp_nodes_to_track[j]];
                projection_rotation.rotate(topos);

                if (SECOND_SHELL)
                {   // the experiment only measures in the symmetry breaking plane, 
                    // so we weight our node choices to favor in plane distance measures

                    score = fabs(curpos.x - topos.x) / 10 + // weak weighing towards the plane
                        fabs( acos( curpos.unitvec().dot( topos.unitvec() ) ) );  // find min angle

                    if ( score < best_score)
                    {
                        best_score = score;
                        best_nodenum = temp_nodes_to_track[j];
                        closest_dist = fabs( acos( curpos.unitvec().dot( topos.unitvec() ) ) );
                    }
                }
                else
                {
                    score = (fabs(curpos.x - topos.x) / RADIUS ) * 2 * TRACK_TARGET_DIST   // weighing towards the plane
                        + fabs( node[nodes_to_track[i]].dist_from_surface - node[temp_nodes_to_track[j]].dist_from_surface ) * 2 * TRACK_TARGET_DIST // weighting towards const radius
                        + fabs( (curpos-topos).length() - TRACK_TARGET_DIST);  // find min distance from targetdist range

                    if ( score < best_score)
                    {
                        best_score = score;
                        best_nodenum = temp_nodes_to_track[j];
                        closest_dist = (curpos-topos).length();
                    }

                }


            }

            

            if (best_nodenum)
            {
                // found closest, so add it to the list

                nodes_to_track.push_back(best_nodenum);

                cout << "Initial Length: " << closest_dist << endl;

            }

        }

    }

    //nodes_to_track.assign(temp_nodes_to_track.begin(), temp_nodes_to_track.end());

    cout << "Tracking " << (int) nodes_to_track.size() << " nodes" << endl;

    if (TRACKS_NO_STATIONARY_NODE)
        stationary_node_number = 0;
    else
        cout << "Node #" << stationary_node_number << " chosen to be stationary" << endl;

    cout << "Node track distances from surface: " ;
        
    for (unsigned int i=0; i != nodes_to_track.size(); ++i)
    {
        cout << node[nodes_to_track[i]].dist_from_surface << " ";
    }

    cout << endl;





}

bool actin::load_nodetracks()
{
    /// load the node tracks already chosen

    ifstream nodetrackfile("nodetracks.txt", ios::in);
    
    if (!nodetrackfile)
        return false;


    unsigned int numnodestotrack, numnodetracks;

    nodetrackfile >> numnodestotrack;

    nodes_to_track.resize(numnodestotrack);

    cout << "Loading " << numnodestotrack << " nodes to track" << endl;
    cout.flush();

    for (unsigned int i=0; i != numnodestotrack; ++i)
    {
        nodetrackfile >> nodes_to_track[i];
    }

    //cout << "Loading node_tracks" << endl;
    //cout.flush();

    for (int axis = xaxis; axis <= zaxis; ++axis)
    {
        nodetrackfile >> numnodetracks;

        node_tracks[axis].resize(numnodetracks);

        for (unsigned int i=0; i != numnodetracks; ++i)
        {
            nodetrackfile >> node_tracks[axis][i];
        }
    }

    nodetrackfile.close();

    return true;
}

bool actin::save_nodetracks()
{
     /// save the node tracks already chosen

    ofstream nodetrackfile("nodetracks.txt", ios::out | ios::trunc);
    
    if (!nodetrackfile)
        return false;


    nodetrackfile << nodes_to_track.size() << endl;

    for (unsigned int i=0; i != nodes_to_track.size(); ++i)
    {
        nodetrackfile << nodes_to_track[i] << " ";
    }

    nodetrackfile << endl;

    for (int axis = xaxis; axis <= zaxis; ++axis)
    {
        nodetrackfile << node_tracks[axis].size() << endl;

        for (unsigned int i=0; i != node_tracks[axis].size(); ++i)
        {
            nodetrackfile << node_tracks[axis][i] << endl;
        }
    }

    nodetrackfile.close();

    return true;

}



void actin::savebmp(const int &filenum, const projection & proj, const processfgbg& fgbg, bool writefile)
{ 
	// choose projection letter for filename etc.

	char projletter[] = "z";
	ofstream *p_outbmpfile;
	char *temp_BMP_filename;

    rotationmatrix axisrotation;
    
    set_axisrotation(proj, axisrotation);

    if (proj == xaxis)
	{
		if (!X_BMP)
			return;
        
		*projletter = 'x';
		p_outbmpfile = &outbmpfile_x;    
		temp_BMP_filename = temp_BMP_filename_x;
	}
	else if (proj == yaxis)
	{	
		if (!Y_BMP)
			return;

		*projletter = 'y';
		p_outbmpfile = &outbmpfile_y;
		temp_BMP_filename = temp_BMP_filename_y;     
	}
	else
	{	
		if (!Z_BMP)
			return;
		p_outbmpfile = &outbmpfile_z;
		temp_BMP_filename = temp_BMP_filename_z;
	}

    if (!QUIET)
	{
		cout << projletter << "-axis BMP";
		cout.flush();
	}

    rotationmatrix projection_rotation;
    projection_rotation.settoidentity();

    // projection_rotation gets us from the bead 
    // frame-of-ref to the bitmap projection frame-of-ref
    projection_rotation.rotatematrix(axisrotation);           // rotates for the projection axis 
    projection_rotation.rotatematrix(sym_break_rotation_to_xy_plane);        // rotates for the symmetry breaking direction
    
    //if (!BMP_FIX_BEAD_ROTATION)
    //    projection_rotation.rotatematrix(nuc_to_world_rot); // compensates for bead rotation
    

    //cout << "Projection Rotation :" << endl << projection_rotation << endl;


    // create copy of the node positions

    vector <vect> rotatednodepositions;
    rotatednodepositions.resize(highestnodecount);
    rotatednodepositions.assign(node.begin(), node.begin()+highestnodecount);

    // rotate the copy into the projection
	for(vector<vect>::iterator i_node  = rotatednodepositions.begin(); 
						       i_node != rotatednodepositions.end();
						     ++i_node)
    {
        projection_rotation.rotate(*i_node);
    }

	
	int x,y;

	// precalculate gaussian for psf

	double gaussmax;
    
    //if (BMP_LINKS_BROKEN)  
    //    gaussmax = (double) GAUSSFWHM * 5 / 2;  // wider for the links broken thing, so we don't notice the circles
    //else
        gaussmax = (double) GAUSSFWHM * 3 / 2; // full extent of gaussian radius -  fwhm is 2/3 this
	
	Dbl2d GaussMat, GaussMat2;

	const int xgmax = pixels(gaussmax);  // was BMP_WIDTH
	const int ygmax = pixels(gaussmax);

	int xg,yg;

	GaussMat.resize(2*xgmax+1);
	GaussMat2.resize(2*xgmax+1);

	for(xg = -xgmax; xg != xgmax+1; xg++)
	{
		GaussMat[xg+xgmax].resize(2*ygmax+1,0.0);
		GaussMat2[xg+xgmax].resize(2*ygmax+1,0.0);
                                                                
		for(yg = -ygmax; yg != ygmax+1; yg++)
		{
			if ( (xg*xg + yg*yg) > (xgmax*ygmax) )
				continue;  // don't do corners

			GaussMat[xg+xgmax][yg+ygmax] 
					= exp(-3*((double)(xg*xg + yg*yg)) / (double)(xgmax*ygmax));

			GaussMat2[xg+xgmax][yg+ygmax] 
					= exp(-12*((double)(xg*xg + yg*yg)) / (double)(xgmax*ygmax));
		}
	}

    int bmpcenterx = BMP_WIDTH  / 2;
    int bmpcentery = BMP_HEIGHT / 2;

    if (BMP_CENTER_ON_LEFT)
        bmpcenterx = BMP_HEIGHT/3; // note, this is meant to be height

	// clear the image array

	for (x = 0; x != BMP_WIDTH; ++x)
	{
        fill(imageR[x].begin(),imageR[x].end(),0.0);
        fill(imageG[x].begin(),imageG[x].end(),0.0);
        fill(imageB[x].begin(),imageB[x].end(),0.0);
	}

	// find extents of network

    double minx, miny, minz;
	double maxx, maxy, maxz; 

    minx = miny = minz = 0.0;
	maxx = maxy = maxz = 0.0;

    for(int i=0; i != highestnodecount; ++i)
    {
		if ((node[i].polymer) && (!node[i].listoflinks.empty()))
		{
			if (minx > rotatednodepositions[i].x) minx = rotatednodepositions[i].x;
			if (miny > rotatednodepositions[i].y) miny = rotatednodepositions[i].y;
			if (minz > rotatednodepositions[i].z) minz = rotatednodepositions[i].z;

			if (maxx < rotatednodepositions[i].x) maxx = rotatednodepositions[i].x;
			if (maxy < rotatednodepositions[i].y) maxy = rotatednodepositions[i].y;
			if (maxz < rotatednodepositions[i].z) maxz = rotatednodepositions[i].z;
		}
	}

    // find the nucleator position, rotated into this view

    vect nucposn = ptheactin->p_nuc->position;
    projection_rotation.rotate(nucposn); 

    vect origin;    // origin is where the bitmap should put the centre
                    // normally 0,0,0 but by putting it at the nucleator position
                    // we make the nuceator stationary

    if (BMP_FIX_BEAD_MOVEMENT)
        origin = nucposn;
    else
        origin.zero();

    // determine offset (movex, movey) needed to keep bead in picture

    int movex = 0;
	int movey = 0;
    
    int beadmaxx, beadmaxy;
	int beadminx, beadminy;

	double keep_within_border;

	if (NUCSHAPE == nucleator::sphere)
		 keep_within_border = 2 * RADIUS;
	else
		 keep_within_border = CAPSULE_HALF_LINEAR + ( 2 * RADIUS );


	beadmaxx = pixels(   keep_within_border - origin.y + nucposn.y) + bmpcenterx;
	beadmaxy = pixels(   keep_within_border - origin.z + nucposn.z) + bmpcentery;

	beadminx = pixels( - keep_within_border - origin.y + nucposn.y) + bmpcenterx;
	beadminy = pixels( - keep_within_border - origin.z + nucposn.z) + bmpcentery;	

	if (beadmaxx > BMP_WIDTH)
		movex = BMP_WIDTH - beadmaxx;

	if (beadmaxy > BMP_HEIGHT)
		movey = BMP_HEIGHT - beadmaxy;

	if (beadminx < 0)
		movex = -beadminx;

	if (beadminy < 0)
		movey = -beadminy;

   

    if ((stationary_node_number != 0) && (stationary_node_number < highestnodecount))
    {   
        // note this is set up initially with the first bitmap call after selecting nodes in comet.cpp
        stationary_node_xoffset =- pixels(rotatednodepositions[stationary_node_number].y);                      
		stationary_node_yoffset =- pixels(rotatednodepositions[stationary_node_number].z);
    }

    if (stationary_node_number != 0)
    {   // if we're locking the bitmap position to this node, add the displacement into movex, movey
        movex += stationary_node_xoffset;
        movey += stationary_node_yoffset;
    }

    // add the node gaussians to the double picture

	vect rot;
	double mult = 0;
    vect originalpos;

    int    SPECKLEGRIDPERIODiter = (int) (SPECKLEGRIDPERIOD    / DELTA_T);
    int SPECKLEGRIDTIMEWIDTHiter = (int) (SPECKLEGRIDTIMEWIDTH / DELTA_T);

    vect tmp_nodepos;

    Colour nodecol;
    nodecol.setwhite();

    // write out the track lengths

    if (BMP_TRACKS && TRACKS_LENGTHS)
    {
    
        ofstream nodetrackdist("nodedistances.txt", ios::out | ios::app);

        for (unsigned int i=0; i != (nodes_to_track.size() / 2)  ; ++i)
        { 

            int firstnode = nodes_to_track[i];
            int secondnode = nodes_to_track[i + (nodes_to_track.size() / 2) ];

            vect firstpos = node[firstnode];
            vect secondpos = node[secondnode];

            if ( ( firstnode > highestnodecount ) ||
                 ( secondnode > highestnodecount ))
            {   // if nodes don't exist yet, then set to origin
                firstpos = secondpos = vect(0,0,0);
            }

            nodetrackdist << (iteration_num / InterRecordIterations) << " " 
                //<< nodes_to_track[i] << " "
                << i+1 << " "
                << (firstpos-secondpos).length() << " "
                << firstpos << " "
                << secondpos << " ";

            nodetrackdist << endl;
        
        }

        nodetrackdist.close();
    }

   


    //  main loop over nodes:

    if (!VECTOR_NOT_BITMAP)
    {

        for(int i=0; i != highestnodecount; ++i)
        {
           if (!node[i].polymer)
			    continue;

            // calculate position in pixels 
            // with displacement to bring bead back in bounds

            x = pixels(rotatednodepositions[i].y - origin.y) + bmpcenterx + movex;                      
		    y = pixels(rotatednodepositions[i].z - origin.z) + bmpcentery + movey;


            // skip if out of focal depth
            if (fabs(rotatednodepositions[i].x - nucposn.x) > FOCALDEPTH)
                continue;

            //mult *= exp(-3.0 * fabs(rotatednodepositions[i].x) / FOCALDEPTH);
           
            // skip if off of the edge of the picture

            if ((x+xgmax < 0) || 
                (x-xgmax >= BMP_WIDTH) ||
                (y+ygmax < 0) || 
                (y-ygmax >= BMP_HEIGHT))
                continue;

            // set speckle magnitude for this point

           mult = 1.0;
           double specmult = 0.0;

		    if (SPECKLE)
            {
                if (SPECKLEGRID)
                {
                    // find position of node on original surface
                    // and rotate into observation frame, (without the bead rotation of course)
                    // to define thickness for speckle slice
                    originalpos = node[i].nucleator_stuck_position;
                    
                    
                    if (!SPECKLE_NO_ROTATE)
                        nuc_to_world_rot.rotate(originalpos);

                    vect posincameraframe = originalpos;

                      
                    sym_break_rotation_to_xy_plane.rotate(originalpos);
                    sym_break_rotation_to_xy_plane.rotate(posincameraframe);

                    axisrotation.rotate(posincameraframe);

                    double segnum;

                    if (fabs(posincameraframe.x) > RADIUS /2 ) // FOCALDEPTH )  // if outside focal depth then black
                    {
                        specmult = 0;
                    }
                    else
                    {
                        if ( ((node[i].creation_iter_num + SPECKLEGRIDPERIODiter/2) % SPECKLEGRIDPERIODiter) < SPECKLEGRIDTIMEWIDTHiter)
                        {
                            specmult = 1.0;  // if within time stripe, then white
                        }
                        else
                        {
                            //if ( ((int)(10 * 360 * (atan2(originalpos.y,originalpos.z)+PI)) % ((int)(10 * 360 * 2 * PI / RADIAL_SEGMENTS))) 
                            //      < (10 * 360 * 2 * PI * (SPECKLEGRIDSTRIPEWIDTH/360)))


                            //if (SPECKLE_NO_ROTATE)
                            //    segnum = ptheactin->p_nuc->segs.getsegmentnum(node[i].nucleator_stuck_position, proj);
                            //else
                            //    segnum = ptheactin->p_nuc->segs.getsegmentnum(originalpos, proj);


                            // note should we pass projection::xaxis, because we've already done the rotation above
                            if (SPECKLE_NO_ROTATE)
                                segnum = ptheactin->p_nuc->segs.getsegmentnum(node[i].nucleator_stuck_position, proj);
                            else
                                segnum = ptheactin->p_nuc->segs.getsegmentnum(originalpos, proj);


                            if ( fabs(segnum + 0.5 - (double)((int)segnum + 0.5)) < SPECKLEGRIDSTRIPEWIDTH )
                                specmult = 1.0;  // if on spoke then white
                            else
                                specmult = 0.0;  // else black
                        }
                    }

                }
                else                                                           
                {
                    specmult = speckle_array[node[i].creation_iter_num % speckle_array_size];
                }
            }
            else
            {
                specmult = 1.0;
            }

            specmult *= 10.0;




    #ifdef BMPS_USING_LINKS   

            for(vector<links>::iterator link_i  = node[i].listoflinks.begin(); 
	                                    link_i != node[i].listoflinks.end();
	                                  ++link_i)
            {

            // for the links, take the position as the midpoint of the link

            double posy = (rotatednodepositions[i].y + rotatednodepositions[link_i->linkednodenumber].y) / 2.0;
            double posz = (rotatednodepositions[i].z + rotatednodepositions[link_i->linkednodenumber].z) / 2.0;

            x = pixels(posy - origin.y) + bmpcenterx;                      
		    y = pixels(posz - origin.z) + bmpcentery;


    #endif 
      
            double value;

    #ifdef BMPS_USING_LINKS

            if (REFERENCEFRAME)
                value = getvaluetoplot(node[i], *link_i) - getvaluetoplot(referencenodes[i], *link_i);
            else
                value = getvaluetoplot(node[i], *link_i);

    #else
            if (REFERENCEFRAME)
                value = getvaluetoplot(node[i]) - getvaluetoplot(referencenodes[i]);
            else
                value = getvaluetoplot(node[i]);

    #endif

            //if (value < 0.01)   // do we want this? messes up the reference stuff
            //    continue;

            if (COL_NODE_BY_STRAIN)
            {
                //if (prob_to_bool(0.01))
                //    cout << " " << value;

                if (COL_INDIVIDUAL_NODES)
                {
                    nodecol.setcol(value);
                }
                else
                {
                    mult = value;
                       
                }
            
            }

            

      
            // add the gaussian
		    for(xg = -xgmax; xg != xgmax+1; ++xg)
            {
			    for(yg = -ygmax; yg != ygmax+1; ++yg)
			    {
				    if ((xg*xg+yg*yg)>(xgmax*ygmax))
					    continue;  // don't do corners

                    if ((x+xg < 0) || 
                        (x+xg >= BMP_WIDTH) ||
                        (y+yg < 0) || 
                        (y+yg >= BMP_HEIGHT))
                        continue;   // skip if outside image bounds
    				
                    

                    if (COL_NODE_BY_STRAIN)
                    {   // just store magnitudes in R, color at the end
                        if (!COL_INDIVIDUAL_NODES)
                        {
                        imageR[x+xg][y+yg] += 10 * mult * GaussMat[xg+xgmax][yg+ygmax];
                        }

                        if (COL_GREY_BGND)
                        {
                            imageG[x+xg][y+yg] += GaussMat[xg+xgmax][yg+ygmax];
                        }

                        if (COL_INDIVIDUAL_NODES)
                        {
                        imageR[x+xg][y+yg] += mult * 10 * nodecol.r * GaussMat[xg+xgmax][yg+ygmax];
                        imageG[x+xg][y+yg] += mult * 10 * nodecol.g * GaussMat[xg+xgmax][yg+ygmax];
                        imageB[x+xg][y+yg] += mult * 10 * nodecol.b * GaussMat[xg+xgmax][yg+ygmax];
                        }
                    }
                    else
                    {
				        imageG[x+xg][y+yg] += mult *
						            GaussMat[xg+xgmax][yg+ygmax];  // amount of actin
        				
				        if (SPECKLE)
                        {
					        imageR[x+xg][y+yg] += specmult *
						            GaussMat2[xg+xgmax][yg+ygmax];  // GaussMat2 is narrower than GaussMat
                        }

                        if (node[i].testnode)
                        {
                            imageB[x+xg][y+yg] += mult *
							        GaussMat2[xg+xgmax][yg+ygmax];

                            imageR[x+xg][y+yg] += mult * node[i].testsurface *
						            GaussMat2[xg+xgmax][yg+ygmax];
                        }
                    }
			    }
            }

           

    #ifdef BMPS_USING_LINKS

        }

    #endif

            // add the tracks data:                                           

        if ((node[i].dist_from_surface > 0.01) && savenodetracks && POST_PROCESS && BMP_TRACKS && ( stationary_node_number != i) &&
                ( find(nodes_to_track.begin(), nodes_to_track.end(), i) != nodes_to_track.end() ))  // is node in the track list?
            {
               
                vect temppos=node[i];// - p_nuc->position;

                world_to_nuc_frame(temppos); // we add the node positions relative to the nucleator

                //if (filenum > TRACK_MIN_RANGE) // start plotting after this point
                    node_tracks[proj].push_back(tracknodeinfo(i, temppos.x, temppos.y, temppos.z,
                                                            x, y, filenum));
            }

	    }

    }

    if (BMP_intensity_scaling)
    {

	    // normalize image

	    for (x = 0; x != BMP_WIDTH; ++x)
	    {
		    for (y = 0; y != BMP_HEIGHT; ++y)
		    {
			    if (imageR[x][y] > imageRmax[proj]) imageRmax[proj]=imageR[x][y];
			    if (imageG[x][y] > imageGmax[proj]) imageGmax[proj]=imageG[x][y];
			    if (imageB[x][y] > imageBmax[proj]) imageBmax[proj]=imageB[x][y];
		    }
	    }
    }

    double Rscale = BMP_INTENSITY_SCALE / imageRmax[proj];
    double Gscale = BMP_INTENSITY_SCALE / imageGmax[proj];
    double Bscale = BMP_INTENSITY_SCALE / imageBmax[proj];

    // scale
    if (( ( COL_NODE_BY_STRAIN && !COL_INDIVIDUAL_NODES ) || SPECKLE ) && fabs(1.0 - NODE_SCALE_GAMMA) > 0.001 ) 
    {
        for (x = 0; x != BMP_WIDTH; ++x)
	        {
		        for (y = 0; y != BMP_HEIGHT; ++y)
		        {
                    imageR[x][y] = pow( imageR[x][y] / imageRmax[proj] , 1.0 / NODE_SCALE_GAMMA);
		        }
	        }

        Rscale = 1.0;
    }

    
    if (SPECKLEGRID)
        Rscale *= 2.0;
     

    double minscale = mymin(mymin(Rscale,Gscale),Bscale);

    if (COL_NODE_BY_STRAIN)
    {   
        if (COL_INDIVIDUAL_NODES)
         { // scale RGB by same scale to keep colors
            
            for (x = 0; x != BMP_WIDTH; ++x)
            {
	            for (y = 0; y != BMP_HEIGHT; ++y)
	            {
		            imageR[x][y] = mymin( (imageR[x][y] * minscale) , 1.0); 
		            imageG[x][y] = mymin( (imageG[x][y] * minscale) , 1.0);
		            imageB[x][y] = mymin( (imageB[x][y] * minscale) , 1.0);
	            }
            }
        
        }
        else
        {
            for (x = 0; x != BMP_WIDTH; ++x)
            {
	            for (y = 0; y != BMP_HEIGHT; ++y)
	            {   
                    double colscale;
                    
                    if (imageR[x][y] > 0.00001)
                    {
                        nodecol.setcol( imageR[x][y] * Rscale );
                        colscale = mymin ( (imageR[x][y]* Rscale) / 0.6, 1.0);   // fade to black for bottom part
                    }
                    else 
                    {
                        colscale = 0.0;
                    }
	                

                    double greyscale = imageG[x][y] * Gscale ;

                    imageR[x][y] = mymin( mymax( greyscale , colscale * nodecol.r), 1.0);
                    imageG[x][y] = mymin( mymax( greyscale , colscale * nodecol.g), 1.0); 
                    imageB[x][y] = mymin( mymax( greyscale , colscale * nodecol.b), 1.0);
                    
	            }
            }
        }
    }
    else
    {   
        if (COL_GREY_BGND)
        {    
            if (NO_BGND)
            {

                // scale RGB separately
                for (x = 0; x != BMP_WIDTH; ++x)
                {
	                for (y = 0; y != BMP_HEIGHT; ++y)
	                {
                        
                        imageR[x][y] = mymin( imageR[x][y] * Rscale, 1.0);
                        imageG[x][y] = mymin( 0                    , 1.0); // green is the greyscale value
                        imageB[x][y] = mymin( imageB[x][y] * Bscale, 1.0);
	                }
                }
            }
            else
            {   
                double greyscale;
                // scale RGB separately
                for (x = 0; x != BMP_WIDTH; ++x)
                {
	                for (y = 0; y != BMP_HEIGHT; ++y)
	                {
                        greyscale = 0.7 * imageG[x][y] * Gscale ;

                        imageR[x][y] = mymin( mymax( greyscale , imageR[x][y] * Rscale), 1.0);
                        imageG[x][y] = mymin( mymax( greyscale , 0                    ), 1.0); // green is the greyscale value
                        imageB[x][y] = mymin( mymax( greyscale , imageB[x][y] * Bscale), 1.0);
	                }
            }
            }

        }
        else
        {
            // scale RGB separately
            for (x = 0; x != BMP_WIDTH; ++x)
            {
	            for (y = 0; y != BMP_HEIGHT; ++y)
	            {
		            imageR[x][y] = mymin( ((imageR[x][y] * Rscale) + BMP_INTENSITY_OFFSET) / (1+BMP_INTENSITY_OFFSET), 1.0);
		            imageG[x][y] = mymin( ((imageG[x][y] * Gscale) + BMP_INTENSITY_OFFSET) / (1+BMP_INTENSITY_OFFSET), 1.0);
		            imageB[x][y] = mymin( ((imageB[x][y] * Bscale) + BMP_INTENSITY_OFFSET) / (1+BMP_INTENSITY_OFFSET), 1.0);
	            }
            }
        }
 
    }



    // bail out if only doing scaling:

    if (!writefile)
        return;

    


	double cagedispy = nucposn.y;
	double cagedispz = nucposn.z;
	int cagemovex = movex;
	int cagemovey = movey;

    if ((CAGE_ON_SIDE) && (NUCSHAPE == nucleator::sphere))
	{   // move cage to side of image
		cagedispy = cagedispz = 0.0;
		cagemovex = p_nuc->segs.centerx - bmpcenterx - xgmax;
		cagemovey = p_nuc->segs.centery - bmpcentery - ygmax + BMP_HEIGHT/4;
	}

    if (BMP_FIX_BEAD_MOVEMENT)
    {
        cagedispy = cagedispz = 0.0;
    }

    const int CAGE_POINT_EXTENT = 4;
                                                       
	if (DRAW_CAGE && !TEST_SQUASH)
	{ // draw the nucleator points cage only if rotation is turned on

	    for (vector <vect>::iterator point  = p_nuc->cagepoints.begin(); 
	                                 point != p_nuc->cagepoints.end(); 
                                   ++point)
	    {
		    // rotate point	

		    rot = *point;

            if (!BMP_FIX_BEAD_ROTATION)     
                nuc_to_world_rot.rotate(rot);

            //sym_break_rotation_to_xy_plane.rotate(rot);        // rotates for the symmetry breaking direction
            //axisrotation.rotate(rot);           // and the projection

            projection_rotation.rotate(rot);

		    double dblx = dbl_pixels(rot.y + cagedispy) + (double) bmpcenterx;
		    double dbly = dbl_pixels(rot.z + cagedispz) + (double) bmpcentery;

            x = (int) (dblx + 0.5);
            y = (int) (dbly + 0.5);

            double xfractpxl = (double) x - dblx;
            double yfractpxl = (double) y - dbly;

		    x += cagemovex;  // displace to bring bead back in bounds
		    y += cagemovey;

            const double point_intensity = 0.1 + 0.7 * ( fabs( rot.x - RADIUS ) / ( 2.0 * RADIUS ) );

            //const double intensitysum = 10.0;

            for (int i  = - CAGE_POINT_EXTENT * BMP_AA_FACTOR; 
                     i !=   CAGE_POINT_EXTENT * BMP_AA_FACTOR + 1; ++i)
            {
                for (int j  = - CAGE_POINT_EXTENT * BMP_AA_FACTOR; 
                         j !=   CAGE_POINT_EXTENT * BMP_AA_FACTOR + 1; ++j)
                {
                    if ((i*i+j*j) > (CAGE_POINT_EXTENT * BMP_AA_FACTOR)*(CAGE_POINT_EXTENT * BMP_AA_FACTOR))  // corners
                        continue;

                    if ((x+i<0) || (x+i >= BMP_WIDTH) ||
			            (y+j<0) || (y+j >= BMP_HEIGHT))  // only plot if point in bounds
			            continue;

                    double dist = calcdist(xfractpxl + (double) i, yfractpxl + (double) j);

                    double intensity = point_intensity * exp( - 3.0 * dist * dist / 
                                 (double) (BMP_AA_FACTOR * CAGE_POINT_EXTENT));
                    
	                imageR[x+i][y+j] = mymin(imageR[x+i][y+j] + intensity , 1.0);
	                imageG[x+i][y+j] = mymin(imageG[x+i][y+j] + intensity , 1.0);
	                imageB[x+i][y+j] = mymin(imageB[x+i][y+j] + intensity , 1.0);
                    
                }
            }

	    }

	}


	
	if (!QUIET)
	{
		cout << ".";
		cout.flush();
	}

	// write the bins graphics

	if (SEGMENT_BINS)
	{
		p_nuc->segs.write_bins_bitmap(imageR, imageG, imageB,
				    p_nuc->segs.link_transverse, p_nuc->segs.link_transverse_scalefactor, proj);
	}
	
    if (COL_NODE_BY_STRAIN && !NO_COLBAR)
        p_nuc->segs.write_colourmap_bitmap(imageR, imageG, imageB); // BMP_AA_FACTOR);
    
    // write the bitmap file
		
	writebitmapfile(*p_outbmpfile, imageR, imageG, imageB);

	
	if (!QUIET)
	{
		cout << ".";
		cout.flush();
	}
	// add the imagemagick overlays
	
	char command1[1024000] = "", command2[1024000] = "" , command3[2048000] = "";
    //char command2[10240];

    // todo figure out how to clear a stringstream so don't have
	// to use as disposable objects
	// .clear() or .str("") don't seem to work.

    double strokewidth = 4;

    stringstream  drawcmd;

    if ( PLOTFORCES || ((POST_PROCESS || POST_PROCESSSINGLECPU) && BMP_TRACKS) )
    {
        // if we're drawing, begin the imagemagick draw command
	    drawcmd << "-fill none -stroke grey -strokewidth " << strokewidth * BMP_AA_FACTOR << " -draw \"";
    }

    if (PLOTFORCES)
    {
        stringstream tmp_drawcmd1,tmp_drawcmd2,tmp_drawcmd3; 
	    
	    p_nuc->segs.drawoutline(drawcmd,proj);				// draw outline of nucleator
    	
        // for these scales, larger number = longer lines

	    // check have lines to draw before adding them
	    // else ImageMagick intermittently crashes

        // scaled impact forces
        // note the values are already scaled to have a max of 1 from the symmetry breaking
        // so we're scaling the 3 colors mainly to show the range of small values
    	


        if (!VECTOR_NOT_BITMAP)
        {
            if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd1,proj,axisrotation, 4 * FORCE_BAR_SCALE) > 0)	
		        drawcmd << "\" -stroke blue   -strokewidth " << BMP_AA_FACTOR * strokewidth / 4 << " -draw \"" << tmp_drawcmd1.str();

            if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd2,proj,axisrotation, 0.8 * FORCE_BAR_SCALE) > 0)	
		        drawcmd << "\" -stroke red    -strokewidth " << BMP_AA_FACTOR * strokewidth / 2 << " -draw \"" << tmp_drawcmd2.str();
        }

        if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd3,proj,axisrotation, 0.16 * FORCE_BAR_SCALE) > 0)	
		        drawcmd << "\" -stroke yellow -strokewidth " << BMP_AA_FACTOR * strokewidth     << " -draw \"" << tmp_drawcmd3.str();
    }

    vect temp_nuc_posn;

    // create the tracks

    if (POST_PROCESS && BMP_TRACKS)
    {

        int lastx = 0, lasty = 0;

        //const int TRACKFRAMESTEP = 5;

        // not the most efficient algorythm, but good enough...

        // find the track node numbers, make sure more than two of each

        vector <int> temptracknodenumbers, tracknodenumbers;
        temptracknodenumbers.resize(0);
        tracknodenumbers.resize(0); 
                
        //cout << "Tot Nodes to track: " << node_tracks[proj].size() << endl;

        for (vector <tracknodeinfo>::iterator i_trackpoint  = node_tracks[proj].begin(); 
                                              i_trackpoint != node_tracks[proj].end(); 
                                            ++i_trackpoint)
        {
            if ((i_trackpoint->frame % TRACKFRAMESTEP == 0) &&   // only plot points *added* every trackframestep frames
                (filenum > i_trackpoint->frame))
                temptracknodenumbers.push_back(i_trackpoint->nodenum);
        }

        //cout << "Selected Nodes to track: " << temptracknodenumbers.size() << endl;

        for (vector <int>::iterator i_pointnodenum  = temptracknodenumbers.begin(); 
                                    i_pointnodenum != temptracknodenumbers.end(); 
                                  ++i_pointnodenum)
        {
            if (count(temptracknodenumbers.begin(), temptracknodenumbers.end(), *i_pointnodenum) > 2)
            { // make sure more than 3 points for the line
                if (find(tracknodenumbers.begin(), tracknodenumbers.end(), *i_pointnodenum) == tracknodenumbers.end())
                {   // if node not already there, add it
                    tracknodenumbers.push_back(*i_pointnodenum);  
                }
            }
        }

        //cout << "ReSelected Nodes to track: " << tracknodenumbers.size() << endl;

        // go through node by node

        for (vector <int>::iterator i_pointnodenum  = tracknodenumbers.begin(); 
                                    i_pointnodenum != tracknodenumbers.end(); 
                                  ++i_pointnodenum)
        {
            // and draw a line for each

            drawcmd << "\" -stroke gray -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \" polyline ";

            int linepoints = 0;

            for (vector <tracknodeinfo>::iterator i_trackpoint  = node_tracks[proj].begin(); 
                                                  i_trackpoint != node_tracks[proj].end(); 
                                                ++i_trackpoint)
            {   
                if (i_trackpoint->nodenum == *i_pointnodenum)
                {

                    // get the old node position in our current frame of ref
                    //tmp_nodepos = i_trackpoint->posn;
                    //temp_nuc_posn = i_trackpoint->nucposn;
                    //i_trackpoint->rotation.inverse().rotate(temp_nuc_posn);
                    //projection_rotation.rotate(temp_nuc_posn);
                    //projection_rotation.rotate(tmp_nodepos);

                    // convert to pixels

                    x = i_trackpoint->BMPx;
                    y = i_trackpoint->BMPy;

                    //x = pixels(i_trackpoint->y - origin.y) + bmpcenterx + movex;                      
		            //y = pixels(i_trackpoint->z - origin.z) + bmpcentery + movey;

                    //x = xgmax + movex + pixels(tmp_nodepos.y - temp_nuc_posn.y) +  bmpcenterx;                      
		            //y = ygmax + movey + pixels(tmp_nodepos.z - temp_nuc_posn.z) + bmpcentery;

                    if ((lastx != x) ||
                        (lasty != y))  // only plot if point different
                    {
                        if (((linepoints % 40) == 0) &&
                             (linepoints > 0))
                        {   // initial draw command at beginning, and every 40 point pairs
                            drawcmd << "\" -stroke magenta -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \" polyline ";
                            drawcmd << lastx << ","        // plot
			                        << lasty << " ";
                        }

                        drawcmd << x << ","        // plot
			                    << y << " ";

                        lastx = x;
                        lasty = y;

                        linepoints++;
                    }
                }
            }    
        }

    }

    if ( PLOTFORCES || (POST_PROCESS && BMP_TRACKS))
    {
        drawcmd << "\"";
    }

    // cout << drawcmd.str() << endl;

	int scalebarmicrons = (int)ceil(RADIUS*2);

    int scalebarlength = pixels(scalebarmicrons);

    // drawing image text takes forever on alpha and OSX
    // so have option to bypass
    // (also changed to get rid of quotes in imagemagick call)

    int text_height_pixels = int ( 25 * (double) TEXT_POINTSIZE / 20.0);

    char sourceimage[1024];

    if (VECTOR_NOT_BITMAP)
    {
        BMP_AA_FACTOR = 1 ; // set to 1 since no bitmap, no point in resizing
        sprintf(sourceimage, " -size %ix%i xc:white ", BMP_WIDTH , BMP_HEIGHT); // or 'transparent' or 'black'
        BMP_OUTPUT_FILETYPE = "svg";
    }
    else
    {
        sprintf(sourceimage, "%s", temp_BMP_filename );
    }

    if (NO_IMAGE_TEXT)
    {	
		//sprintf(command1,
		//"%s -quality %i -fill white -draw \"rectangle 5 576 %i 573\" %s %s %s%s_proj_%05i.%s", 
        //IMAGEMAGICKCONVERT, BMP_COMPRESSION, scalebarlength+5,drawcmd.str().c_str(), temp_BMP_filename, BITMAPDIR, 
		//	 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());

        sprintf(command1,
		"%s -quality %i -density 300x300 \"%s\" %s \"%s%s_proj_%05i.%s\"", 
        IMAGEMAGICKCONVERT, BMP_COMPRESSION, sourceimage, drawcmd.str().c_str(), BITMAPDIR, 
			 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());
    }
    else 
    {

		sprintf(command1,
		"%s \"%s\" -quality %i -font helvetica -fill white -pointsize %i -draw \"text %i %i '%ium' rectangle %i %i %i %i text +%i+%i '%s-projection' text +%i+%i 'Frame % 6i' text +%i+%i 'Time % 6i'\" -density 300x300 %s \"%s%s_proj_%05i.%s\"", 
            IMAGEMAGICKCONVERT, sourceimage, BMP_COMPRESSION,   // command and bitmap quality
            TEXT_POINTSIZE * BMP_AA_FACTOR,                           // font size
            5 * BMP_AA_FACTOR, BMP_HEIGHT - 5 * BMP_AA_FACTOR, scalebarmicrons,   // scale bar text
            5 * BMP_AA_FACTOR, BMP_HEIGHT - (text_height_pixels + 4 ) * BMP_AA_FACTOR, scalebarlength + 5 * BMP_AA_FACTOR,  BMP_HEIGHT - (text_height_pixels + 7) * BMP_AA_FACTOR,
            5 * BMP_AA_FACTOR,       text_height_pixels   * BMP_AA_FACTOR,  // first line of text
            projletter, 
            5 * BMP_AA_FACTOR, ( 2 * text_height_pixels ) * BMP_AA_FACTOR , 
            filenum,
			5 * BMP_AA_FACTOR, ( 3 * text_height_pixels ) * BMP_AA_FACTOR , 
            int(filenum * InterRecordIterations * DELTA_T), drawcmd.str().c_str(), BITMAPDIR,
			 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());
    }


    if (BMP_AA_FACTOR == 1)
    {   
        if ((fgbg == runbg) && (!NOBGIMAGEMAGICK))
            sprintf(command3, "%s &", command1);
        else
            sprintf(command3, "%s", command1);
    }
    else
    {
        sprintf(command2,"%s -quality %i -resize %f%% -density 300x300 \"%s%s_proj_%05i.%s\" ",
             IMAGEMAGICKMOGRIFY, BMP_COMPRESSION, 100/(double)BMP_AA_FACTOR, BITMAPDIR,
		     projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());

        // only command1 refers to the temp bmp, so we can always safely call the antialias resize command in the background

        

        if (!NOBGIMAGEMAGICK)
        {
            if (fgbg == runbg)
	            sprintf(command3, "(%s ; %s ) &", command1, command2);
            else
	            sprintf(command3, "%s ; %s &", command1, command2);
        }
        else
        {
            sprintf(command3, "%s ; %s", command1, command2);
        }

       
    }

    //opruninfo << command3 << endl;

	p_outbmpfile->flush();  // must flush the bitmap file buffer before calling imagemagick

    //cout << command3 << endl;

	system(command3);

	if (!QUIET)
	{
	cout << ". ";
	cout.flush();
	}

	return;
}

double actin::getvaluetoplot(nodes & mynode)
{

    if (BMP_LINKS_BROKEN)
        return mynode.links_broken / COL_INDIVIDUAL_SCALE;
    else if (BMP_TRANSVERSELINKSONLY)
        return (- mynode.linkforce_transverse + mynode.repforce_transverse)/ COL_INDIVIDUAL_SCALE;
    else if (BMP_RADIALLINKSONLY)
        return (- mynode.linkforce_radial + mynode.repforce_radial) / COL_INDIVIDUAL_SCALE;
    else  // energy
        return  (mynode.linkforce_radial * mynode.linkforce_radial +
                 mynode.linkforce_transverse * mynode.linkforce_transverse) / COL_INDIVIDUAL_SCALE;

}

double actin::getvaluetoplot(nodes & mynode, links & mylink)
{
    double value; 

    if (BMP_LINKS_BROKEN)
        value = mynode.links_broken / COL_INDIVIDUAL_SCALE;
    else
        value = mylink.forcesum / (InterRecordIterations * LINK_BREAKAGE_FORCE) ;   // scale by 1/LINK_BREAKAGE_FORCE

    vect linkvec = mynode - *(mylink.linkednodeptr);

    double dirfactor = fabs ( linkvec.unitvec().dot(mynode.unit_vec_posn) );  // we don't care about sign, just angle

    if (COL_LINK_BY_DIRN)
        value = dirfactor ;  // put the actual direction in the links
    else
        value = value * ( 1 - dirfactor );

    return value;

}

void actin::writebitmapheader(ofstream& outbmpfile, const int & bitmapwidth, const int & bitmapheight)
{

	// bitmap headers etc (see microsoft documentation---large chunks pasted from given code):

	// define data structures for bitmap header:

	#define DWORD unsigned int
	#define LONG unsigned int
	#define WORD unsigned short int
	#define BYTE unsigned char
	//#define FOURCC unsigned int

	#pragma pack(push,1)  // align the structs to byte boundaries

	typedef struct tagRGBQUAD {
		BYTE    rgbBlue; 
		BYTE    rgbGreen; 
		BYTE    rgbRed; 
		BYTE    rgbReserved; 
	} RGBQUAD; 

	typedef struct tagRGB {
		BYTE    B; 
		BYTE    G; 
		BYTE    R; 
	} RGB; 


	typedef struct tagBITMAPFILEHEADER { 
		WORD    bfType; 
		DWORD   bfSize; 
		WORD    bfReserved1; 
		WORD    bfReserved2; 
		DWORD   bfOffBits; 
	} BITMAPFILEHEADER, *PBITMAPFILEHEADER; 

	typedef struct tagBITMAPINFOHEADER{
		DWORD  biSize; 
		LONG   biWidth; 
		LONG   biHeight; 
		WORD   biPlanes; 
		WORD   biBitCount; 
		DWORD  biCompression; 
		DWORD  biSizeImage; 
		LONG   biXPelsPerMeter; 
		LONG   biYPelsPerMeter; 
		DWORD  biClrUsed; 
		DWORD  biClrImportant; 
	} BITMAPINFOHEADER, *PBITMAPINFOHEADER; 

	typedef struct tagBITMAPINFO { 
	BITMAPINFOHEADER bmiHeader; 
	RGBQUAD          bmiColors[1]; 
	} BITMAPINFO, *PBITMAPINFO; 

	#pragma pack(pop)

	
	// the header for saving the bitmaps...

	BITMAPFILEHEADER *fileHeader;
	BITMAPINFO       *fileInfo;

	fileHeader = (BITMAPFILEHEADER*)calloc(1, sizeof( BITMAPFILEHEADER ));
	fileInfo = (BITMAPINFO*)calloc(1, sizeof( BITMAPINFO ) + 256 * sizeof(RGBQUAD));

	fileHeader->bfType = 0x4d42;
	fileHeader->bfSize = (DWORD)(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFO) + (bitmapwidth * bitmapheight* 3));
	fileHeader->bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFO) + (256*sizeof(RGBQUAD));

	fileInfo->bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
	fileInfo->bmiHeader.biWidth         = bitmapwidth;
	fileInfo->bmiHeader.biHeight        = bitmapheight;
	fileInfo->bmiHeader.biPlanes        = 1;
	//fileInfo->bmiHeader.biBitCount      = 8;
	fileInfo->bmiHeader.biBitCount      = 24;
	fileInfo->bmiHeader.biCompression   = 0;  // BI_RGB = 0
	fileInfo->bmiHeader.biXPelsPerMeter = 1000;
	fileInfo->bmiHeader.biYPelsPerMeter = 1000;

#ifdef __BIG_ENDIAN__

    // if we're running on a big endian machine (e.g. ppc mac), we need to swap the bytes:

#pragma pack(push,1)  // align the structs to byte boundaries
	endian_swap(fileHeader->bfType);
	endian_swap(fileHeader->bfSize);
	endian_swap(fileHeader->bfOffBits);
	endian_swap(fileInfo->bmiHeader.biSize);
	endian_swap(fileInfo->bmiHeader.biWidth);
	endian_swap(fileInfo->bmiHeader.biHeight);
	endian_swap(fileInfo->bmiHeader.biPlanes);
	endian_swap(fileInfo->bmiHeader.biBitCount);
	endian_swap(fileInfo->bmiHeader.biCompression);
	endian_swap(fileInfo->bmiHeader.biXPelsPerMeter);
	endian_swap(fileInfo->bmiHeader.biYPelsPerMeter);
#pragma pack(pop)

#endif



	for (int i=0;i!=255;++i)
	{ // the bitmap palette  (not used for 24 bit, of course)
		fileInfo->bmiColors[i].rgbBlue = (BYTE) i;
		fileInfo->bmiColors[i].rgbGreen = (BYTE) i;
		fileInfo->bmiColors[i].rgbRed = (BYTE) i;
		fileInfo->bmiColors[i].rgbReserved = (BYTE) 0;
	}

	outbmpfile.seekp(0); // move to start

    // save headers

	outbmpfile.write((char*)fileHeader,sizeof(BITMAPFILEHEADER));
	outbmpfile.write((char*)fileInfo,sizeof(BITMAPINFO) + 256*sizeof(RGBQUAD));

	bitmap_start = outbmpfile.tellp();  // keep a note of where the picture data starts

	free(fileHeader);
	free(fileInfo);

}

void actin::writebitmapfile(ofstream& outbmpfile, const Dbl2d& myimageR, const Dbl2d& myimageG, const Dbl2d& myimageB)
{  // re-scale for byte output and save image data

//#define BYTE unsigned char

#pragma pack(push,1)  // align the structs to byte boundaries

	typedef struct tagRGB 
	{
		unsigned char    B; 
		unsigned char    G; 
		unsigned char    R; 
	} RGB; 

#pragma pack(pop)

    const int width  = (int) myimageR.size();
    const int height = (int) myimageR[0].size();

    //cout << "Width :" << width << "  Height: " << height << endl;

	RGB *line;
    line = new RGB[width];

	outbmpfile.seekp(bitmap_start); // move to start                          

	// write out the data, line by line (note: y is backwards)
	for (int y = (height-1); y>=0; y--)
		{ //cout << "Line : " << y << endl;
		//outbmpfile.write(picbuff + (bitmapwidth*y), bitmapwidth);
		for (int x = 0; x < width; x++)
			{
			line[x].B=(unsigned char)(255 * myimageB[x][y]);
			line[x].G=(unsigned char)(255 * myimageG[x][y]);
			line[x].R=(unsigned char)(255 * myimageR[x][y]);
			}
			outbmpfile.write((char*)line,width*3);
		}

	outbmpfile.flush();

	delete [] line;

}


void actin::squash(const double & thickness)
{  
	// squash with 'coverslip'

    const double halfthickness = thickness/2.0;

    nuc_struck_coverslip = false;

    // squash the network:

	for(vector <nodes>::iterator	i_node  = node.begin() + lowestnodetoupdate; 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        if (i_node->x >  halfthickness)  // above
            i_node->x =  halfthickness; 
        else 
        if (i_node->x < -halfthickness)  // below
			i_node->x = -halfthickness;
	}

    if (NO_X_MOTION)
    {
        p_nuc->position.x = 0.0;
        return;
    }

    // make sure the nucleator is in bounds, too
    // todo: do the capsule shape, too

    if (p_nuc->position.x >  halfthickness - RADIUS)
    {
        p_nuc->position.x =  halfthickness - RADIUS;
        //cout << "Warning - Nucleator has hit the coverslip" << endl;
        nuc_struck_coverslip = true;
    }

    if (p_nuc->position.x < - halfthickness + RADIUS)
    {
        p_nuc->position.x = - halfthickness + RADIUS;
        //cout << "Warning - Nucleator has hit the coverslip" << endl;
        nuc_struck_coverslip = true;
    }

    


	return;
}




void actin::sortnodesbygridpoint(void)
{
    // todo: change this to gridpoints-by-thread

//	nodes* nodeptr, *startnodeptr;


    fill(donenode.begin()+lowestnodetoupdate, donenode.begin()+highestnodecount, false);


    int threadnum;

    for (threadnum = 0; threadnum < NUM_THREAD_DATA_CHUNKS; ++threadnum)
    {
        nodes_by_thread[threadnum].resize(0);
    }

    threadnum = 0;


	for (int i=lowestnodetoupdate; i != highestnodecount; ++i)
	{	// collect the nodes in gridpoint order...

		if ((donenode[i]) || (!node[i].polymer))
			continue;

        nodes_by_thread[threadnum].insert(nodes_by_thread[threadnum].end(),
                                      NODEGRID(node[i].gridx, node[i].gridy, node[i].gridz).begin(),
                                      NODEGRID(node[i].gridx, node[i].gridy, node[i].gridz).end());

        for(NODEGRIDTYPE <nodes*>::iterator 
                    i_node  = NODEGRID(node[i].gridx, node[i].gridy, node[i].gridz).begin(); 
                    i_node != NODEGRID(node[i].gridx, node[i].gridy, node[i].gridz).end();
                  ++i_node)
        {
            donenode[(*i_node)->nodenum] = true;
        }

        // just scatter randomly between thread queues

        //threadnum = (threadnum + 1) % NUM_THREAD_DATA_CHUNKS;

		// find smallest thread queue and put next one in there

	    size_t minsize,threadsize;
        minsize = MAXNODES;

        for (int tn = 0; tn != NUM_THREAD_DATA_CHUNKS; ++tn)
        {
            threadsize = nodes_by_thread[tn].size();
            if (threadsize < minsize)
            {
                minsize = threadsize;
                threadnum = tn;
            }
        }

	}

	return;
}



void actin::sortgridpointsbythread(void)
{
    // todo: change this to gridpoints-by-thread

    // find the total number of nodes in each thread 

    vector <int> nodesinthread(NUM_THREAD_DATA_CHUNKS,0);

    int averagelength = 0;

    for (int i = 0; i != NUM_THREAD_DATA_CHUNKS; ++i)
    {
        for (vector<NODEGRIDTYPE<nodes*>*>::iterator 
                                                i_gp  = gridpointsbythread[i].begin();
		                                        i_gp != gridpointsbythread[i].end() ;
							                  ++i_gp )
        {
            nodesinthread[i] += (int) (*i_gp)->size();
        }

        averagelength += nodesinthread[i];
    }

    averagelength /= NUM_THREAD_DATA_CHUNKS;

}


void actin::compressfilesdowork(const int & filenum)
{

#ifndef _WIN32

	char command1[1024];
	// report files
		
	//sprintf(command1, "(%s \"%s*report*.txt\" 2>/dev/null; mv \"%s*report*%s\" \"%s\" 2>/dev/null) &"
	//	,COMPRESSCOMMAND,TEMPDIR,TEMPDIR, COMPRESSEDEXTENSION, REPORTDIR);

    // the 2>/dev/null on the move command is to
    // supress an odd  "set owner/group (was: 501/0): Operation not permitted"
    // warning
	sprintf(command1, "(%s \"%s\"*report*.txt ; mv \"%s\"*report*%s \"%s\" 2>/dev/null) &"
		,COMPRESSCOMMAND,TEMPDIR,TEMPDIR, COMPRESSEDEXTENSION, REPORTDIR);
    //cout << command1 << endl;
    system(command1);

	// save data file

    if (COMPRESSDATAFILES)
    {
	sprintf(command1 , "(%s \"%s\"*data*.txt ; mv \"%s\"*data*%s \"%s\" 2>/dev/null) &",
	     COMPRESSCOMMAND, TEMPDIR, TEMPDIR, COMPRESSEDEXTENSION, DATADIR);
	system(command1);
    }

                                       
	// wrl file

	//sprintf(command1 , "(gzip -9 -f -c \"%snodes%05i.wrl\" > \"%snodes%05i.wrz\" 2>/dev/null) &",
	//					TEMPDIR, filenum,  VRMLDIR,filenum);
	//system(command1);

#else

	char command2[1024];

	//sprintf(command1, "%s \"%s*report*.txt\"",COMPRESSCOMMAND, TEMPDIR);
	sprintf(command2, "move \"%s*report%05i.txt\" \"%s\"",TEMPDIR,filenum,REPORTDIR);
	//system(command1);
	system(command2);

	// save data file

	//sprintf(command1, "%s \"%s*data*.txt\"",COMPRESSCOMMAND, TEMPDIR2);
	//sprintf(command2, "move \"%s*data*%s\" \"%s\" >NUL",TEMPDIR,COMPRESSEDEXTENSION, DATADIR);
	sprintf(command2, "move \"%s*data*.txt\" \"%s\" >NUL",TEMPDIR,DATADIR);
	//system(command1);
	system(command2);
	//cout << command1 << endl << command2 << endl;

	// wrl file

	//sprintf(command1, "%s %snodes%05i.wrl",COMPRESSCOMMAND, TEMPDIR);
	//sprintf(command2, "move %snodes%05i.wrl%s %snodes%05i.wrz",TEMPDIR, COMPRESSEDEXTENSION, VRMLDIR);
	//system(command1);
	//system(command2);


#endif


}



int actin::find_center(vect &center)
{
	center = p_nuc->position;

	return 0;	  
}
		

void actin::clear_nodegrid()
{
    // clear the nodegrid
	
	for(int i=0; i!=GRIDSIZE; ++i)
		for(int j=0; j!=GRIDSIZE; ++j)
			for(int k = 0; k !=GRIDSIZE; ++k)
				NODEGRID(i,j,k).clear();

}

int actin::save_data(ofstream &ofstrm)
{
    // write out all stateful information for this object    
    // save actin
    ofstrm << "actin: "  << endl 
	   << highestnodecount << " "
	   << nexttocrosslink << " " 
	   << iteration_num << " " 
	   << linksbroken  << " " 
	   << linksformed << " "
       << world_to_nuc_rot << " "
	   << sym_break_rotation_to_xy_plane << " "
	   << sym_break_rotation_to_zaxis << endl;

    // save nucleator (must be before nodes, else unit vectors are wrong when read in)
    ofstrm << "nucleator:" << endl;
    p_nuc->save_data(ofstrm);
    
    // save nodes
    ofstrm << "nodes-links:" << endl;
    for(int i=0; i < highestnodecount; i++) {
	node[i].save_data(ofstrm);
    }
    
    // save linkdelays
    ofstrm << "xlinkdelays:" << endl 
	   << (unsigned int)crosslinknodesdelay.size() << ")";
    
    for(vector<int>::iterator i = crosslinknodesdelay.begin(); 
	i < crosslinknodesdelay.end();i++) {
	ofstrm << (*i);
	if(i < crosslinknodesdelay.end()-1)
	    ofstrm << " ";
	else
	    ofstrm << endl;
    }


    
    return 0;
}

bool actin::load_data(ifstream &ifstr)
{
	
    // clear the nodegrid
    clear_nodegrid();
    
    string str;
    unsigned char ch;
    
    // load actin
    ifstr >> str;
    // ensure the identifier for the start of the actin
    if(str.compare("actin:") !=0 )
    {
	    cout << "error in checkpoint file, 'actin:' expected '" << str << "' found" << endl;
	    return false;
    }

	ifstr >> highestnodecount 
	  >> nexttocrosslink  
	  >> iteration_num    
	  >> linksbroken      
	  >> linksformed      
	  >> world_to_nuc_rot   
	  >> sym_break_rotation_to_xy_plane  
	  >> sym_break_rotation_to_zaxis;

    nuc_to_world_rot = world_to_nuc_rot.inverse();

 
  
        // load nucleator
    ifstr >> str;
    if(str.compare("nucleator:") !=0 )
    {
	    cout << "error in checkpoint file, 'nucleator:' expected '" << str << "' found" << endl;
	    return false;
    }

    p_nuc->load_data(ifstr);

     // load nodes

    ifstr >> str;
    
    if(str.compare("nodes-links:") !=0 )
    {
	    cout << "error in data file, 'nodes-links:' expected '" << str << "' found" << endl;
	    cout << "'" << str <<"' last read." << endl;
	    return true;
    }



    // ** Remember the node vector is preallocated to MAXNODES
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
    {
		if (!i_node->load_data(ifstr))  // false if fails
        {
            return false;
        }

    }


    //if (POST_STATS || POST_VTK || (!REWRITESYMBREAK && !POST_PROCESS))
    rebuildnodepointers();

    // check node list
    /*
    for(int i=0; i < highestnodecount; i++) {
	cout << "node [" << i << ":" << &node[i] << "] " << node[i].listoflinks.size() << " ";
	for(vector<links>::iterator l=node[i].listoflinks.begin(); 
	    l<node[i].listoflinks.end(); ++l) {

	    cout << l->linkednodenumber << ":" << l->linkednodeptr << ", ";

	}
	cout << endl;
    }
    */
    
    // load xlinkdelays
    ifstr >> str;
    if(str.compare("xlinkdelays:") !=0 )
    {
	    cout << "error in checkpoint file, 'xlinkdelays:' expected" 
	         << endl;
	    return false;
    }
    
    int numcrosslinkdelay;
    ifstr >> numcrosslinkdelay >> ch;
    if(ch!=')' )
    {
	    cout << "error in checkpoint file, xlinkdelays 'NN)' expected" 
	     << endl;
	    return false;
    }
    //crosslinknodesdelay.clear();
    crosslinknodesdelay.resize(numcrosslinkdelay);
    // load each delay
    for(vector<int>::iterator i  = crosslinknodesdelay.begin(); 
	                          i != crosslinknodesdelay.end(); ++i) 
    {
	    ifstr >> (*i) ;
    }



    
    return true;
}

void actin::rebuildnodepointers()
{
    // This is not neat, and we have a nasty reliance on the data being
    // public all the way down to the nodes linklist link objects.
    // When loading the link has been stored as an index, and we now
    // convert that to a pointer.
    // Method is 
    //  Run through and rebuild the node ptrs for the linkednode indices
    //  important this is done after the full node list has been built
    //    - iterate over nodes
    //    - then iterate over links in the node linklist
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
	{
	    for(vector<links>::iterator i_link  = i_node->listoflinks.begin(); 
	                                i_link != i_node->listoflinks.end();
                                  ++i_link) 
        {
	        if(i_link->linkednodenumber>=0) // ensure the indx was explicitly set
				i_link->linkednodeptr = &node[i_link->linkednodenumber];
			else
				cout << "Link node < 0" << endl;
	    }
    }
}

void actin::setdontupdates(void)
{ 
    if (highestnodecount > NODES_TO_UPDATE)
	{
        lowestnodetoupdate = highestnodecount - NODES_TO_UPDATE;
	} else
	{
		lowestnodetoupdate = 0;
	}

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+lowestnodetoupdate;
							      ++i_node)
    {
        if (i_node->stucktonucleator) // must not allow the non-updated nodes to stay attached!
            i_node->stucktonucleator = false;

    }
}

void actin::set_sym_break_axes(bool constrain_to_zy_plane, vect sym_break_direction)
{
    vect tmp_nodepos;
	//vect CofM;
	//vect sym_break_direction;

	rotationmatrix tmp_rotation, tmp_rotation2;	
	rotationmatrix final_rotation;

    tmp_rotation.settoidentity();
 
    sym_break_x_angle = 0.0;
    sym_break_y_angle = 0.0;
    sym_break_z_angle = 0.0; 

    if (NUCSHAPE!=nucleator::capsule)
    {   // only rotate by x and y if not capsule, otherwise just by 90 degrees in z axis (after z rot calc)

        //find_center(sym_break_direction);  // which way did the bead go?

         
        //if (constrain_to_zy_plane) // if we're using coverslip, flatten direction to zy plane
        //   sym_break_direction.x=0.0;

        tmp_rotation.settoidentity();
        sym_break_x_angle = atan2(sym_break_direction.y,sym_break_direction.z);

        tmp_rotation.rotatematrix(sym_break_x_angle, 0 , 0);
        tmp_rotation.rotate(sym_break_direction);


        sym_break_y_angle = -atan2(sym_break_direction.x,sym_break_direction.z);

        tmp_rotation.rotatematrix(0, sym_break_y_angle , 0); // now pointing upwards
        tmp_rotation.rotate(sym_break_direction);


        // construct x & y rotation matrix

        tmp_rotation.settoidentity();

        tmp_rotation.rotatematrix(sym_break_y_angle , yaxis);
        tmp_rotation.rotatematrix(sym_break_x_angle , xaxis);

    }


	// we now have tmp_rotation which will transform the bead direction to be along the z axis
    // to align the break, we 
	// need to determine the z rotation by the principlal axis

    const double threshold = 1.005; // angle is average between edges.  Edge determined by only resetting if new value greater than this times old

	double chi, maxchi;

	maxchi = 0.0;
	sym_break_z_angle = 0.0;

	for(double theta = -PI; theta < PI; theta+=PI/360)    // PI/360 i.e. 0.5 degree intervals
    {
		chi = 0.0;

		tmp_rotation2.settoidentity();
		tmp_rotation2.rotatematrix(0,0,theta);

		for (int i=0; i != highestnodecount; ++i)
		{
			if ((!node[i].polymer) || (node[i].harbinger))
				continue;

			tmp_nodepos = node[i] - p_nuc->position; // node position relative to bead

			tmp_rotation.rotate(tmp_nodepos);   // rotate points to bring in line with sym break dir'n
			tmp_rotation2.rotate(tmp_nodepos);	// rotate z

            //tmp_nodepos = tmp_nodepos.unitvec();  // normalize so we don't weight further ones more


			
			//if ( tmp_nodepos.z <  0 ) // only look at front of bead (i.e. new z<0, in direction of movement)
			//{
				chi += fabs(tmp_nodepos.y);// - fabs(tmp_nodepos.x); // sum the distances from the axis
			//}
			
		}

		if (chi > (maxchi * threshold))
		{
			maxchi = chi;
			sym_break_z_angle = theta;
		}

	}

    double angle1=sym_break_z_angle;
    // other way
    maxchi = 0.0;
	sym_break_z_angle = 0.0;

    for(double theta = PI; theta > -PI; theta-=PI/360)    // PI/360 i.e. 0.5 degree intervals
    {
		chi = 0.0;

		tmp_rotation2.settoidentity();
		tmp_rotation2.rotatematrix(0,0,theta);

		for (int i=0; i != highestnodecount; ++i)
		{
			if ((!node[i].polymer) || (node[i].harbinger))
				continue;

			tmp_nodepos = node[i] - p_nuc->position; // node position relative to bead

			tmp_rotation.rotate(tmp_nodepos);   // rotate points to bring in line with sym break dir'n
			tmp_rotation2.rotate(tmp_nodepos);	// rotate z

            //tmp_nodepos = tmp_nodepos.unitvec();  // normalize so we don't weight further ones more


			
			//if ( tmp_nodepos.z <  0 ) // only look at front of bead (i.e. new z<0, in direction of movement)
			//{
				chi += fabs(tmp_nodepos.y);// - fabs(tmp_nodepos.x); // sum the distances from the axis
			//}
			
		}

		if (chi > (maxchi * threshold))
		{
			maxchi = chi;
			sym_break_z_angle = theta;
		}

	}

    if (sym_break_z_angle < 0 && angle1 > 0)
        sym_break_z_angle += 2*PI;  // make sure both angles are on same side

    sym_break_z_angle = (sym_break_z_angle + angle1) / 2;

    sym_break_z_angle += PI/2;


    sym_break_rotation_to_xy_plane.settoidentity();

    // now have all the angles
	// need to assemble them in the right (reverse) order z,y,x
	// to get the rotation matrix

	// sym_break_rotation_to_xy_plane includes a pre-rotation with the old x angle
	// to prevent it moving always onto the z-axis
	// sym_break_rotation_to_zaxis does not have this
    
    if (constrain_to_zy_plane) // if we're using coverslip, flatten direction to zy plane
    {
        cout << "Symbreak constrained to zy plane" << endl;
        sym_break_z_angle=0.0;
        sym_break_y_angle=0.0;
    }

    if (NUCSHAPE==nucleator::capsule)
    {
        sym_break_z_angle += PI/2;   // rotate into x/y plane for capsule
    }
    else
    {
        if (SYM_BREAK_TO_RIGHT)   // only do this for sphere
        {
            sym_break_rotation_to_xy_plane.rotatematrix(-PI/2,xaxis); // rotate so that always breaks to the right

        } else // restore the component of the sym break in the yz plane
        {
            sym_break_rotation_to_xy_plane.rotatematrix(-sym_break_x_angle, xaxis);
        }
    }

    

	sym_break_rotation_to_xy_plane.rotatematrix(sym_break_z_angle, zaxis);
	sym_break_rotation_to_xy_plane.rotatematrix(sym_break_y_angle, yaxis);
	sym_break_rotation_to_xy_plane.rotatematrix(sym_break_x_angle, xaxis);

    

	sym_break_rotation_to_zaxis.settoidentity();

	sym_break_rotation_to_zaxis.rotatematrix(sym_break_z_angle, zaxis);
	sym_break_rotation_to_zaxis.rotatematrix(sym_break_y_angle, yaxis);
	sym_break_rotation_to_zaxis.rotatematrix(sym_break_x_angle, xaxis);

    
    if (NO_SYMBREAK_ROTATION)
    {
        sym_break_rotation_to_xy_plane.settoidentity();
        sym_break_rotation_to_zaxis.settoidentity();
        sym_break_x_angle = sym_break_y_angle = sym_break_z_angle = 0.0;

    }


    //reverse_sym_break_rotation_to_xy_plane = sym_break_rotation_to_xy_plane.inverse();

	cout << setprecision(1) << "Camera rotation angles: " << sym_break_x_angle*180/PI << ", " << sym_break_y_angle*180/PI << ", " << sym_break_z_angle*180/PI << endl;

	symbreakiter = iteration_num;

}


void actin::save_sym_break_axes(void)
{
	ofstream opsymbreak(SYM_BREAK_FILE, ios::out | ios::trunc);
	if (!opsymbreak) 
	{ cout << "Unable to open file 'sym_break_axis.txt' for output"; return;}

	opsymbreak  << symbreakiter << endl
				<< sym_break_rotation_to_xy_plane << endl
				<< sym_break_rotation_to_zaxis << endl
                << sym_break_x_angle << " " << sym_break_y_angle << " " << sym_break_z_angle << endl;

	opsymbreak.close();

	cout << "'sym_break_axis.txt' file written" << endl;
}

bool actin::load_sym_break_axes(void)
{
	ifstream ipsymbreak("sym_break_axis.txt", ios::in);

	if (!ipsymbreak) 
	{ 
		cout << "Unable to open file 'sym_break_axis.txt' for input," << endl
			<< "symmetry breaking direction will be incorrect" << endl;
        return false;
	}
	else
	{
		ipsymbreak  >> symbreakiter 
					>> sym_break_rotation_to_xy_plane
					>> sym_break_rotation_to_zaxis
                    >> sym_break_x_angle >> sym_break_y_angle >> sym_break_z_angle ;

		ipsymbreak.close();  
        
        //reverse_sym_break_rotation_to_xy_plane = sym_break_rotation_to_xy_plane.inverse();
        
        return true;
	}
}


void actin::clear_node_stats(void)
{
	attemptedpolrate = polrate	= 0;

	for (int i=lowestnodetoupdate; i != highestnodecount; ++i)
	{
		node[i].clearstats();
	}

    p_nuc->deltanucposn_sum.zero();

}


void actin::reservemorenodes(const int extranodes)
{
    // have to do all this in case reallocating vectors causes it to move
    // the pointers

    const int oldsize = (int) node.size();

    // remove nodes from nodegrid
    // todo: make this more efficient (maybe, but it's not called very often)

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        i_node->removefromgrid();
    }


    // convert the link pointers to link node numbers
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
	{
	    for(vector<links>::iterator i_link  = i_node->listoflinks.begin(); 
	                                i_link != i_node->listoflinks.end();
                                  ++i_link) 
        {
            i_link->linkednodenumber = i_link->linkednodeptr->nodenum;
	    }
    }


    // resize the vectors
	node.resize(oldsize + extranodes);
	donenode.resize(oldsize + extranodes);

    // set the pointers within nodes
    for(int i=0; i != highestnodecount; ++i)
    {
        node[i].nodenum = i;
    }

    // set the link pointers from the link numbers
    rebuildnodepointers();

    // re-add to the grid
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
    {
        i_node->addtogrid();
    }

    MAXNODES = (int) node.size();

}

void actin::addbrownianforces()
{
    

    const double brownianscale = BROWNIANFORCESCALE * DELTA_T;

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
    {
        if (i_node->harbinger)
            continue;

        // qn: is it worth making this spherically symmetric?

        i_node->rep_force_vec.x += brownianscale * ( rand_0to1() - 0.5 );
        i_node->rep_force_vec.y += brownianscale * ( rand_0to1() - 0.5 );
        i_node->rep_force_vec.z += brownianscale * ( rand_0to1() - 0.5 );
    }
}
