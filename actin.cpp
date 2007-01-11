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

    if (!REWRITESYMBREAK)
    {
	    opvelocityinfo.open("velocities.txt", ios::out | ios::trunc);
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
    

	opvelocityinfo << "time,x,y,z,vel" << endl;

	crosslinknodesdelay.resize(CROSSLINKDELAY);

	//doreportiteration = 9999999;

	donenode.resize(MAXNODES);

	node.resize(MAXNODES);// ,nodes(0,0,0, this));

	for (int i=0; i<MAXNODES; i++)
	{
		node[i].nodenum=i;
		//node[i].ptheactin=this;
	}

	cout << "     GridExtent : " << GRIDBOUNDS << " uM" << endl;
	cout << " GridResolution : " << GRIDRES << " uM" << endl;
	cout << "TotalGridpoints : " << (GRIDSIZE*GRIDSIZE*GRIDSIZE) << endl;

	cout << "Zeroing Grid...";
	cout.flush();
    
    	     
    #ifdef NODE_GRID_USE_ARRAYS

	    nodegrid = new NG1d[(Node_Grid_Dim+1)*(Node_Grid_Dim+1)*(Node_Grid_Dim+1)];

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

        int i,j,k;
        double x,y,z;

        for (i=0; i!=(GRIDSIZE+1); i++)  
	    {
            x = (i - (GRIDSIZE/2)) * GRIDRES;
		    for (j=0; j!=(GRIDSIZE+1); j++)
		    {
                y = (j - (GRIDSIZE/2)) * GRIDRES;
                for (k=0; k!=(GRIDSIZE+1); k++)
                {
                    z = (k - (GRIDSIZE/2)) * GRIDRES;
                    if ((!REWRITESYMBREAK && !POST_PROCESS) &&
                        ((fabs(x) < RADIUS * 8) ||
                         (fabs(y) < RADIUS * 8) ||
                         (fabs(z) < RADIUS * 8)))
                        NODEGRID(i,j,k).reserve(64);
                }
		    }
	    }

    #endif

	cout << "Done" << endl << endl;

	speckle_array_size = 10000;
	speckle_array.resize(speckle_array_size);


	double x1, x2, w, y1; //, y2;

	for (int i=0; i<speckle_array_size; i++)
	{
		// gaussian random numbers (see http://www.taygeta.com/random/gaussian.html)

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
			//y2 = x2 * w;

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
		recti_near_nodes[i].reserve(1024);
		nodes_on_same_gridpoint[i].reserve(1024);
        nodes_by_thread[i].reserve(2 * MAXNODES / NUM_THREAD_DATA_CHUNKS);  // poss should / by NUM_THREAD_DATA_CHUNKS
        linkremovefrom[i].reserve(128);
        linkremoveto[i].reserve(128);

        recti_near_nodes[i].resize(0);
        nodes_on_same_gridpoint[i].resize(0);
        nodes_by_thread[i].resize(0);
        linkremovefrom[i].resize(0);
        linkremoveto[i].resize(0);
	}

	//nodes_within_nucleator.reserve(10000);
	linkformto.reserve(128);

	cout << "Memory for Grid : " << (sizeof(nodes*)*GRIDSIZE*GRIDSIZE*GRIDSIZE/(1024*1024)) << " MB" << endl;

	cout << "Memory for nodes: " << (node.size() * sizeof(node[0])/(1024*1024)) << " MB" << endl;


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
	outbmpfile_y.open(temp_BMP_filename_y, ios::out | ios::binary | ios::trunc);
	outbmpfile_z.open(temp_BMP_filename_z, ios::out | ios::binary | ios::trunc);

	if ((!outbmpfile_x) || (!outbmpfile_y) || (!outbmpfile_z)) 
	{ cout << "Unable to open file temp bitmap file for output"; return;}

	// we write the header now, and just keep changing the pixel part as we write the frames

	writebitmapheader(outbmpfile_x,BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_y,BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_z,BMP_WIDTH, BMP_HEIGHT);
    
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

    actin_rotation.settoidentity();
    inverse_actin_rotation.settoidentity();

    gridpointsbythread.resize(NUM_THREAD_DATA_CHUNKS);

    for(int i = 0; i != NUM_THREAD_DATA_CHUNKS; ++i)
    {
        gridpointsbythread[i].reserve(2048);
    }
    currentsmallestgridthread = 0;


    lasttestsurfacesavedposn = testsurfaceposn = DBL_MAX;
    testsurfacerotation = 0;
    testangle = 0;
    testforcemag = 0;
     
    testdirection = vect(0,0,1);
    testangle = cos(30 * PI/360);

    test_equilibrating = false;

    node_tracks.resize(3);

    node_tracks[xaxis].reserve(10240);
    node_tracks[yaxis].reserve(10240);
    node_tracks[zaxis].reserve(10240);

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

	char command1[255];

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_x );
	system(command1);

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_y );
	system(command1);

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_z );
	system(command1);

    system("stty sane");
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


int actin::saveinfo()
{
	//char filename[255];
	double minx, miny, minz;
	double maxx, maxy, maxz; 
	//char time[255], date[255];

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

	minx = miny = minz =   GRIDBOUNDS * 2;
	maxx = maxy = maxz = - GRIDBOUNDS * 2;
	int existantnodes = 0;

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


	opruninfo << "Final existant nodes: " << existantnodes << endl;
	cout << "Final existant nodes: " << existantnodes << endl;
	// write info

	opruninfo << "min x, min y, min z, max x, max y, max z" << endl;

	opruninfo << minx << "," << miny << "," << minz << "," 
					  << maxx << "," << maxy << "," << maxz << endl;



	cout << "min x, min y, min z, max x, max y, max z" << endl;

	cout << minx << ", " << miny << ", " << minz << ", " 
					  << maxx << ", " << maxy << ", " << maxz << endl;

	return 0;
}

int actin::savevrml(int filenum)
{
	char filename[255];
	//char time[255], date[255];

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
		if (highestnodecount > 400)
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
        {

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


	if ((lastsorthighestnode != highestnodecount) && ((iteration_num % 10) == 0))
	{    // only do full sort periodically, and if we have new nodes
		sortnodesbygridpoint(); // sort the nodes so they can be divided sanely between threads
		lastsorthighestnode = highestnodecount;
	} else
	{	// else just throw the new nodes in whichever node list
		for (int i=lastsorthighestnode; i!=highestnodecount; ++i)
		{
			nodes_by_thread[iteration_num % NUM_THREAD_DATA_CHUNKS].push_back(&node[i]);
		}
		lastsorthighestnode	= highestnodecount;
	}

	// in multithreaded mode, these functions just start the threads
    // else they do the work

	collisiondetection();	    // calc node-to-node repulsion
	linkforces();			    // and link forces

    if  (currentlyusingthreads && (USETHREAD_COLLISION || USETHREAD_LINKFORCES))
    {   // must wait for threads to finish before updating node positions
        thread_queue.complete_queued_tasks();
    }

	if (!TEST_SQUASH)
        nucleator_node_interactions();	    // do forcable node ejection
	
	applyforces();              // move the nodes and update the grid

    if (COVERSLIPGAP > RADIUS)   // skip if less than RADIUS (it's set to 0 to disable) 
	    squash(COVERSLIPGAP);	 

	iteration_num++;

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
        camera_rotation.rotate(tempunitvector);

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

	// crosslink poisson distribution (max probability at XLINK_NODE_RANGE/5)

	//pxlink = ((double)2.7*dist/(XLINK_NODE_RANGE/(double)5))
	//						*exp(-dist/(XLINK_NODE_RANGE/(double)5));

	// divide by dist*dist to compensate for increase numbers of nodes at larger distances



	//pxlink = ((( (double) 2.7 / (XLINK_NODE_RANGE/(double)5) )
	//						*exp(-dist/(XLINK_NODE_RANGE/(double)5) ) ));

	// normal distrib with center around 1/2 of XLINK_NODE_RANGE
	// and scale magnitude by 1/(XLINK_NODE_RANGE^2) to compensate for increased
	// number of nodes at this distance shell

	//if ((node[linknode1].z <0) || (node[linknode2].z <0))
 	//	return 0;

	// was using this one before 28 march:
	// (poissonian with max at XLINK_NODE_RANGE/5):
	//pxlink = ((( (double) 2.7 / (XLINK_NODE_RANGE/(double)5) )
	//						*exp(-dist/(XLINK_NODE_RANGE/(double)5) ) )/dist);


	// was using this (gaussian with max at XLINK_NODE_RANGE/2) for a bit:
	//pxlink = exp( -40*(dist-(XLINK_NODE_RANGE/2))*(dist-(XLINK_NODE_RANGE/2)))/(dist*dist*XLINK_NODE_RANGE*XLINK_NODE_RANGE);
	
	//double pxlink = 1;

	double dist = calcdist(linknode1,linknode2);

	//double pxlink = P_XLINK * exp(-4*dist/XLINK_NODE_RANGE);
	
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

    // now using insidenucleator flag set in setunitvec() of node

    for(vector <nodes>::iterator	i_node  = node.begin()+lowestnodetoupdate; 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
	{
        if (i_node->dist_from_surface < 0) // if node is inside nucleator           
        {   
            if (p_nuc->collision(*i_node)) // (*i)->x,(*i)->y,(*i)->z)==0)  
			    i_node->updategrid();  // ejected OK
		    else
			    i_node->depolymerize();  // not ejected OK, depolymerize
        }

		// do node-nucleator links:

        if (i_node->stucktonucleator)
		{
            disp = *i_node - i_node->nucleator_stuck_position;
	        dist = disp.length();

            force = NUC_LINK_FORCE * dist;
    	    	        
            if (force > NUC_LINK_BREAKAGE_FORCE)
            {
                //node[n].nucleator_link_force.zero();
                i_node->stucktonucleator = false;   // no longer stuck
            }
            else
            {

                forcevec = disp * (force/dist);  // '(disp/dist)' is just to get the unit vector

                i_node->link_force_vec -= forcevec;			// add to move node
				
                // add to node segment stats (tension only since default dist is zero)
		        i_node->adddirectionalmags(forcevec, i_node->linkforce_radial, i_node->linkforce_transverse);

				vect tomove = -forcevec * DELTA_T * FORCE_SCALE_FACT;  // convert force to distance

				p_nuc->move_nuc(*i_node,tomove);		// add to nucleator movement vector

				i_node->nucleator_link_force += tomove;	// add to nuc link force stats

            }
            
        }

	}

}

void actin::collisiondetection(void)
{
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

//void actin::findnearbynodes_collision_setup(const int& adjgridpoints)
//{
//	nearby_collision_gridpoint_offsets.reserve(128);
//	nearby_collision_gridpoint_offsets.resize(0);
//
//	int x,y,z;
//	
//	for (x = -adjgridpoints; x != adjgridpoints + 1; ++x) 
//    {
//	    for (y = -adjgridpoints; y != adjgridpoints + 1; ++y) 
//        {
//            for (z = -adjgridpoints; z != adjgridpoints + 1; ++z)	 
//            {
//				// skip the 0,0,0 gridpoint since will be added separately
//				if ((x!=0) || (y!=0) || (z!=0))
//					nearby_collision_gridpoint_offsets.push_back((Node_Grid_Dim*Node_Grid_Dim*x) + (Node_Grid_Dim*y) + z);
//            }
//        }
//    }
//
//	nearby_collision_gridpoint_offset_begin = nearby_collision_gridpoint_offsets.begin();
//	nearby_collision_gridpoint_offset_end   = nearby_collision_gridpoint_offsets.end();
//
//}


void * actin::collisiondetectiondowork(void* arg)//, pthread_mutex_t *mutex)
{	
    // cast arg
    const thread_data* const dat = (thread_data*) arg;

	const double local_NODE_REPULSIVE_RANGEsqared = NODE_REPULSIVE_RANGE * NODE_REPULSIVE_RANGE;

    vect nodeposvec;
    double distsqr;
    //double dist;
    double recipdist;
    //double recip_dist;
	double rep_force_mag;
    vect rep_force_vect;
    vect disp;

    int x,y,z;
    //int minx,miny,minz;
    //int maxx,maxy,maxz;

#ifdef PROXIMITY_VISCOSITY
	double viscfactor;
#endif

    nodes *p_i_node, *p_nearnode, *p_sameGPnode;


    int sameGPnodenum;

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
			            for(NODEGRIDTYPE <nodes*>::iterator nearnode  = NODEGRID(x,y,z).begin(); 
									                        nearnode != NODEGRID(x,y,z).end();
		                                                  ++nearnode) 
			            {
                            p_nearnode = *nearnode;

			                if ( p_sameGPnode == p_nearnode)
				                continue;  // skip if self


			                disp = *p_nearnode - nodeposvec;   // relative position of node
			                
                            
                            distsqr = disp.sqrlength();    // square of distance
                			                               // avoid sqrt at this point, in case it's out of range and we don't need to do the repulsion
               			
			                if (distsqr < local_NODE_REPULSIVE_RANGEsqared)
			                {

				                //dist = SSEsqrt(distsqr);  
                                recipdist = SSErsqrt(distsqr);
                                // calculate repulsive force
  
                                // this was the old function
					            //rep_force_mag = 13 * NODE_REPULSIVE_MAG * (exp ((-7.0*dist/NODE_REPULSIVE_RANGE)) - exp (-7.0));

                                if (recipdist < (1.0/0.05) )
                                   rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE * recipdist, NODE_REPULSIVE_POWER ) - 1 );
                               else
                                   rep_force_mag = 0.06 * NODE_REPULSIVE_MAG * ( pow( NODE_REPULSIVE_RANGE / 0.05 , NODE_REPULSIVE_POWER ) - 1 ) ;
 

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

                                // add to the statistics
				                p_sameGPnode->adddirectionalmags(rep_force_vect, p_sameGPnode->repforce_radial,
								                                                 p_sameGPnode->repforce_transverse);

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
	if (highestnodecount == 0)
		return;

	// set the rotation matrix

	double x_angle = 0.0;
	double y_angle = 0.0;			
	double z_angle = 0.0;

	torque_rotate.settoidentity();

	if (IMPOSED_NUC_ROT)
	{
		// rotate const speed

		x_angle = IMPOSED_NUC_ROT_SPEED * 2 * PI * DELTA_T;

		torque_rotate.rotatematrix( x_angle, 0, 0);

	}
	else if (ROTATION)
	{
		x_angle = p_nuc->torque.x / p_nuc->momentofinertia.x;
		y_angle = p_nuc->torque.y / p_nuc->momentofinertia.y;
		z_angle = p_nuc->torque.z / p_nuc->momentofinertia.z;

		torque_rotate.rotatematrix( x_angle, y_angle, z_angle);
	}

    if (!TEST_SQUASH)
    {   // don't move nucleator if testing forces

	    if (IMPOSED_NUC_ROT || ROTATION)
	    {
            // rotate the actin reference frame:
            actin_rotation.rotatematrix(torque_rotate);
            inverse_actin_rotation = actin_rotation.inverse();

            // rotate the nucleator displacement vector
	        torque_rotate.rotate(p_nuc->position);

            // todo: is this right?

            // update the nucleator position
	        p_nuc->position += p_nuc->deltanucposn;

	    }

	    nuc_disp = p_nuc->deltanucposn; // store nucleator movement in static for threads

    }
    else
    {   // TEST_SQUASH is true, so lock the nucleator position
        nuc_disp.zero();
    }

	// and zero
	p_nuc->deltanucposn.zero();

	// clear the torque vector now that the torque rotation matrix is set
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
    
    // note: we have to updategrid() for *all* the nodes, because of
    // the rotation and translation, and we can't do it in threads because
    // of the linked list in the nodegrid.
	// this is relatively slow, but not too bad.  could possibly thread this based on 
	// xyz co-ords of nodes, if the removefromgrid function is OK
    
    // note: is this necessary?  only use grid for repulsion anyway
    // and we're not calculating that

	for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
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
			    torque_rotate.rotate(*i_node);	 // rotate
			    *i_node -= nuc_disp;	 // move wrt nucleator frame of ref
                i_node->applyforces();
			    //if (i>=lowestnodetoupdate)
				   // i_node->applyforces();	         // move according to forces
		    }
        }
    }
    else
    {   // just move/rotate
        for(vector <nodes>::iterator i_node  = dat->startnode;
                                     i_node != dat->endnode;
                                   ++i_node)
	    {
		    torque_rotate.rotate(*i_node);	 // rotate
		    *i_node -= nuc_disp;	             // move wrt nucleator frame of ref
        }
    }

    return NULL;
}

void actin::linkforces()
{

    // remove the links for ones that were broken last time
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

            //if (i==NUM_THREAD_DATA_CHUNKS-1)
	        //{	// put remainder in last thread (cludge for now)
		       // end += highestnodecount % NUM_THREAD_DATA_CHUNKS;
	        //}

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
 
    double dist, force;
    vect nodeposvec, disp;
    vect forcevec;

#ifdef LINK_VISCOSITY
	double viscfactor;
#endif


    // go through all nodes
    for(vector <nodes>::iterator i_node  = dat->startnode;
                                 i_node != dat->endnode;
                               ++i_node)
    {   
        // make sure they have links and are polymers
	    if ((i_node->listoflinks.size() == 0) || (!i_node->polymer))
	        continue;

    	// store the node position for later calculation of link lengths
	    nodeposvec = *i_node;
   	
  		// go through links for each node
		for (vector <links>::iterator i_link  = i_node->listoflinks.begin();
                                      i_link != i_node->listoflinks.end();
                                    ++i_link )
		{	 
			assert( !i_link->broken ); // if link not broken  (shouldn't be here if broken anyway)

			disp = nodeposvec - *(i_link->linkednodeptr);
			dist = disp.length();
    	    
            if (i_link->getlinkforces(dist,force)) // returns false if broken
            {
				forcevec = disp * (force/dist); // convert force to vector
        		
				i_node->link_force_vec += forcevec;
        		
				if (force < 0) // put tension into link forces
				{
					i_node->adddirectionalmags(forcevec, i_node->linkforce_radial, i_node->linkforce_transverse);
				} 
				else 
				{	   // but put compression into repulsive forces
					i_node->adddirectionalmags(forcevec, i_node->repforce_radial , i_node->repforce_transverse );
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
        axisrotation.rotatematrix(PI/2, yaxis); // x to z is by rotation around y and z
        axisrotation.rotatematrix(PI/2, zaxis);
	}
}



void actin::set_nodes_to_track(const projection & proj)
{
    nodes_to_track.resize(0);

    rotationmatrix axisrotation;
    set_axisrotation(proj, axisrotation);

    rotationmatrix projection_rotation;
    projection_rotation.settoidentity();

    // projection_rotation gets us from the bead 
    // frame-of-ref to the bitmap projection frame-of-ref
    projection_rotation.rotatematrix(axisrotation);           // rotates for the projection axis 
    projection_rotation.rotatematrix(camera_rotation2);       // rotates for the symmetry breaking direction
                                                              // using camera_rotation2 so we can tell the furthest node back
                                                              // from the sym break direction
    
    if (!BMP_FIX_BEAD_ROTATION)
        projection_rotation.rotatematrix(inverse_actin_rotation); // compensates for bead rotation

    cout << "Selecting nodes from frame " << TRACK_MIN_RANGE << " to " << TRACK_MAX_RANGE << endl;

    vect tempposn;

    stationary_node_number = 0;    // node to lock the bitmap to (opposite the sym break direction)
    double furthest_node_posn = 0;

    for(int i=0; i != highestnodecount; ++i)
    {
        tempposn = node[i];
        projection_rotation.rotate(tempposn);

        if (fabs(tempposn.x) * 4 < RADIUS)  // if within RADIUS/2
        {
            if ((node[i].creation_iter_num > TRACK_MIN_RANGE * InterRecordIterations) && 
                (node[i].creation_iter_num < TRACK_MAX_RANGE * InterRecordIterations))  // is within range?
            {
                nodes_to_track.push_back(i);
                if (furthest_node_posn < tempposn.z)
                {
                    furthest_node_posn = tempposn.z;
                    stationary_node_number = i;
                 }
            }
        }

    }

    cout << "Tracking " << (int) nodes_to_track.size() << " nodes" << endl;
    cout << "Node #" << stationary_node_number << " chosen to be stationary" << endl;

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
    projection_rotation.rotatematrix(camera_rotation);        // rotates for the symmetry breaking direction
    
    if (!BMP_FIX_BEAD_ROTATION)
        projection_rotation.rotatematrix(inverse_actin_rotation); // compensates for bead rotation
    

    // create vector with rotated node positions

    vector <vect> rotatednodepositions;
    rotatednodepositions.resize(highestnodecount);

    for(int i=0; i != highestnodecount; ++i)
    {
        rotatednodepositions[i] = node[i];
        projection_rotation.rotate(rotatednodepositions[i]);
	}


	double minx, miny, minz;
	double maxx, maxy, maxz; 

	int beadmaxx, beadmaxy;
	int beadminx, beadminy;
	
	int x,y;

    int bmpcenterx = BMP_WIDTH/2;
    int bmpcentery = BMP_HEIGHT/2;

    if (SYM_BREAK_TO_RIGHT)
        bmpcenterx = BMP_WIDTH/3; 

	// precalculate gaussian for psf

	double gaussmax = (double) GAUSSFWHM * 3 / 2;  // full extent of gaussian radius -  fwhm is 2/3 this
	
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
					= exp(-9*((double)(xg*xg + yg*yg)) / (double)(xgmax*ygmax));
		}
	}

	// clear the image array

	for (x = 0; x != BMP_WIDTH; ++x)
	{
        fill(imageR[x].begin(),imageR[x].end(),0.0);
        fill(imageG[x].begin(),imageG[x].end(),0.0);
        fill(imageB[x].begin(),imageB[x].end(),0.0);
	}


	// determine size

	minx = miny = minz = 0;
	maxx = maxy = maxz = 0;

    double meanx = 0.0, meany = 0.0, meanz = 0.0;

    int countnodes = 0;

	// find extents of network and mean position

    for(int i=0; i != highestnodecount; ++i)
    {
		if ((node[i].polymer) && (!node[i].listoflinks.empty()))
		{
            //meanx += rotatednodepositions[i].x;
            //meany += rotatednodepositions[i].y;
            //meanz += rotatednodepositions[i].z;

            countnodes++;

			if (minx > rotatednodepositions[i].x) minx = rotatednodepositions[i].x;
			if (miny > rotatednodepositions[i].y) miny = rotatednodepositions[i].y;
			if (minz > rotatednodepositions[i].z) minz = rotatednodepositions[i].z;

			if (maxx < rotatednodepositions[i].x) maxx = rotatednodepositions[i].x;
			if (maxy < rotatednodepositions[i].y) maxy = rotatednodepositions[i].y;
			if (maxz < rotatednodepositions[i].z) maxz = rotatednodepositions[i].z;
		}
	}

    //meanx /= (double) count + 100; // + 100 buffers the instability when # nodes is v. small
    //meany /= (double) count + 100;
    //meanz /= (double) count + 100;

    vect nucposn = - ptheactin->p_nuc->position;
    projection_rotation.rotate(nucposn); 

    if (!BMP_FIX_BEAD_MOVEMENT)
    {
        meanx = nucposn.x;
        meany = nucposn.y;
        meanz = nucposn.z;
    }

    // determine offset (if any) needed to keep bead in picture

	double keep_within_border;

	if (p_nuc->geometry == nucleator::sphere)
		 keep_within_border = 2* RADIUS;
	else
		 keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);


	beadmaxx = pixels(  keep_within_border - meany) +  bmpcenterx;
	beadmaxy = pixels(  keep_within_border - meanz) +  bmpcentery;
	beadminx = pixels(- keep_within_border - meany) +  bmpcenterx;
	beadminy = pixels(- keep_within_border - meanz) +  bmpcentery;

	int movex = 0;
	int movey = 0;

	if (beadmaxx > BMP_WIDTH)
		movex =- (beadmaxx - BMP_WIDTH);
	if (beadmaxy > BMP_HEIGHT)
		movey =- (beadmaxy - BMP_HEIGHT);
	if (beadminx < 0)
		movex =- beadminx;
	if (beadminy < 0)
		movey =- beadminy;


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
    vect origposn;

    int    SPECKLEGRIDPERIODiter = (int) (SPECKLEGRIDPERIOD    / DELTA_T);
    int SPECKLEGRIDTIMEWIDTHiter = (int) (SPECKLEGRIDTIMEWIDTH / DELTA_T);

    vect tmp_nodepos;

    Colour nodecol;
    nodecol.setwhite();


    //  main loop over nodes:

    for(int i=0; i != highestnodecount; ++i)
    {
		if (!node[i].polymer)
			continue;

        // set speckle magnitude for this point

		if (SPECKLE)
        {
            if (SPECKLEGRID)
            {
                // find position of node on original surface
                // and rotate into observation frame, (without the bead rotation of course)
                // to define thickness for speckle slice
                origposn = node[i].nucleator_stuck_position;
                camera_rotation.rotate(origposn);  
                axisrotation.rotate(origposn);  
                
                if (fabs(origposn.x) * 2 > RADIUS )  // if outside RADIUS/2 then black
                {
                    mult = 0;
                }
                else
                {
                    if ( (node[i].creation_iter_num % SPECKLEGRIDPERIODiter) < SPECKLEGRIDTIMEWIDTHiter)
                    {
                        mult = 1.0;  // if within time stripe, then white
                    }
                    else
                    {
                        //if ( ((int)(10 * 360 * (atan2(origposn.y,origposn.z)+PI)) % ((int)(10 * 360 * 2 * PI / RADIAL_SEGMENTS))) 
                        //      < (10 * 360 * 2 * PI * (SPECKLEGRIDSTRIPEWIDTH/360)))

                        double segnum = ptheactin->p_nuc->segs.getsegmentnum(node[i].nucleator_stuck_position,proj); 

                        if ( fabs(segnum + 0.5 - (double)((int)segnum + 0.5)) < SPECKLEGRIDSTRIPEWIDTH )
                            mult = 1.0;  // if on spoke then white
                        else
                            mult = 0.0;  // else black
                    }
                }

            }
            else                                                           
            {
                mult = speckle_array[node[i].creation_iter_num % speckle_array_size];
            }
        }
        else
        {
            mult = 1.0;
        }

        x = pixels(rotatednodepositions[i].y - meany) +  bmpcenterx;                      
		y = pixels(rotatednodepositions[i].z - meanz) + bmpcentery;

#ifdef BMP_USE_FOCAL_DEPTH
            mult *= exp(-3.0 * fabs(rotatednodepositions[i].x) / FOCALDEPTH);
#endif

        // displace to bring bead back in bounds
        // and add the offset for the gaussian

		x += movex + xgmax;  
		y += movey + ygmax;

        //if (node[i].testnode)
        //    cout << "Plotting test node: " << i << endl;

        // add the gaussian

        if (COL_NODE_BY_STRAIN)
        {
            //double force;
            //double sum_force = 0;
            //double distance=0;

            //for(vector<links>::iterator l=node[i].listoflinks.begin(); l!=node[i].listoflinks.end(); ++l)
            //{
	           // l->getlinkforces(distance, force);
            //    sum_force += fabs(force);
            //}

            double value = node[i].linkforce_transverse + node[i].linkforce_radial;

	        double y = value / (4500 * LINK_BREAKAGE_FORCE);    
	        //y = pow( y , 1/VTK_LINK_COLOUR_GAMMA);	    
	        //y = y*0.9+0.1;
            //if (prob_to_bool(0.01))
            //    cout << " " << y;
            nodecol.setcol(y);
        }
  

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
				
				//imageR[x+xg+xgmax][y+yg+ygmax]+=		// link forces
				//	node[i].linkforce_transverse[0] * GaussMat[xg+xgmax][yg+ygmax];
				
                if (COL_NODE_BY_STRAIN)
                {
                    // just store magnitudes in R, color at the end
                    imageR[x+xg][y+yg] += 10 * nodecol.r * GaussMat[xg+xgmax][yg+ygmax];
                    imageG[x+xg][y+yg] += 10 * nodecol.g * GaussMat[xg+xgmax][yg+ygmax];
                    imageB[x+xg][y+yg] += 10 * nodecol.b * GaussMat[xg+xgmax][yg+ygmax];
                }
                else
                {
				    imageG[x+xg][y+yg] += mult *
						        GaussMat[xg+xgmax][yg+ygmax];  // amount of actin
    				
				    if (SPECKLE)
                    {
					    imageR[x+xg][y+yg] += mult *
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

        // add the tracks data:

        if (POST_PROCESS && BMP_TRACKS && ( stationary_node_number != i) &&
            ( find(nodes_to_track.begin(), nodes_to_track.end(), i) != nodes_to_track.end() ))  // is node in the track list?
        {
            node_tracks[proj].push_back(tracknodeinfo(i, x, y, filenum));
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

    const double Rscale = BMP_INTENSITY_SCALE / imageRmax[proj];
    const double Gscale = BMP_INTENSITY_SCALE / imageGmax[proj];
    const double Bscale = BMP_INTENSITY_SCALE / imageBmax[proj];

    if (COL_NODE_BY_STRAIN)
    {   // scale RGB by same scale to keep colors
        double minscale = mymin(mymin(Rscale,Gscale),Bscale);
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
    {   // scale RGB separately
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



    // bail out if only doing scaling:

    if (!writefile)
        return;

    //double cagedispx = meanx;
	double cagedispy = meany;
	double cagedispz = meanz;
	int cagemovex = movex;
	int cagemovey = movey;

	if ((CAGE_ON_SIDE) && (p_nuc->is_sphere()))
	{   // move cage to side of image
		cagedispy = cagedispz = 0.0;
		cagemovex = p_nuc->segs.centerx -  bmpcenterx - xgmax;
		cagemovey = p_nuc->segs.centery - bmpcentery - ygmax + BMP_HEIGHT/4;
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

            //camera_rotation.rotate(rot);        // rotates for the symmetry breaking direction
            //axisrotation.rotate(rot);           // and the projection

            projection_rotation.rotate(rot);

		    double dblx = dbl_pixels(rot.y - cagedispy) + (double) bmpcenterx;
		    double dbly = dbl_pixels(rot.z - cagedispz) + (double) bmpcentery;

            x = (int) (dblx + 0.5);
            y = (int) (dbly + 0.5);

            double xfractpxl = (double) x - dblx;
            double yfractpxl = (double) y - dbly;

		    x += cagemovex + xgmax;  // displace to bring bead back in bounds
		    y += cagemovey + ygmax;

            //double intensitysum = 0.0;  // normalise intensity on point-by-point basis (why do we need to do this??)

            //for (int i  = - CAGE_POINT_EXTENT * BMP_AA_FACTOR; 
            //         i !=   CAGE_POINT_EXTENT * BMP_AA_FACTOR + 1; ++i)
            //{
            //    for (int j  = - CAGE_POINT_EXTENT * BMP_AA_FACTOR; 
            //             j !=   CAGE_POINT_EXTENT * BMP_AA_FACTOR + 1; ++j)
            //    {
            //        if ((i*i+j*j) > (CAGE_POINT_EXTENT * BMP_AA_FACTOR)*(CAGE_POINT_EXTENT * BMP_AA_FACTOR))  // corners
            //            continue;

            //        if ((x+i<0) || (x+i >= BMP_WIDTH) ||
			         //   (y+j<0) || (y+j >= BMP_HEIGHT))  // only plot if point in bounds
			         //   continue;

            //        double dist = calcdist(xfractpxl + (double) i, yfractpxl + (double) j);

            //        double intensity =  exp( - 2.0 * dist * dist / 
            //                     (double) (BMP_AA_FACTOR * CAGE_POINT_EXTENT));
            //        
	           //     intensitysum += intensity;
            //                                                       
            //    }
            //}

            const double intensitysum = 10.0;

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

                    double intensity = (10.0/intensitysum) * exp( - 3.0 * dist * dist / 
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
	
    if (COL_NODE_BY_STRAIN)
        p_nuc->segs.write_colourmap_bitmap(imageR, imageG, imageB);
    
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

    stringstream tmp_drawcmd1,tmp_drawcmd2,tmp_drawcmd3, drawcmd;

    if ( PLOTFORCES || (POST_PROCESS && BMP_TRACKS))
    {
        // if we're drawing, begin the imagemagick draw command
	    drawcmd << "-fill none -stroke grey -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \"";
    }

    if (PLOTFORCES)
    {

	    
	    p_nuc->segs.drawoutline(drawcmd,proj);				// draw outline of nucleator
    					
	    // check have lines to draw before adding them
	    // else ImageMagick intermittently crashes

        // scaled impact forces
        // note the values are already scaled to have a max of 1 from the symmetry breaking
        // so we're scaling the 3 colors mainly to show the range of small values
    	
        //if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd1,proj, 1 * FORCE_BAR_SCALE) > 0)		
	    //	drawcmd << "\" -stroke blue -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \"" << tmp_drawcmd1.str();

	    if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd2,proj,0.5 * FORCE_BAR_SCALE) > 0)	
		    drawcmd << "\" -stroke red -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \"" << tmp_drawcmd2.str();	

	    if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd3,proj, 0.1 * FORCE_BAR_SCALE) > 0)	
		    drawcmd << "\" -stroke yellow -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \"" << tmp_drawcmd3.str();
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

        for (vector <tracknodeinfo>::iterator i_trackpoint  = node_tracks[proj].begin(); 
                                              i_trackpoint != node_tracks[proj].end(); 
                                            ++i_trackpoint)
        {
            if (i_trackpoint->frame % TRACKFRAMESTEP == 0)   // only plot points *added* every trackframestep frames
                temptracknodenumbers.push_back(i_trackpoint->nodenum);
        }

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


        // go through node by node

        for (vector <int>::iterator i_pointnodenum  = tracknodenumbers.begin(); 
                                    i_pointnodenum != tracknodenumbers.end(); 
                                  ++i_pointnodenum)
        {
            // and draw a line for each

            drawcmd << "\" -stroke magenta -strokewidth " << BMP_AA_FACTOR + 1 << " -draw \" polyline ";

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

                    x = i_trackpoint->x;
                    y = i_trackpoint->y;

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


    if (NO_IMAGE_TEXT)
    {	
		//sprintf(command1,
		//"%s -quality %i -fill white -draw \"rectangle 5 576 %i 573\" %s %s %s%s_proj_%05i.%s", 
        //IMAGEMAGICKCONVERT, BMP_COMPRESSION, scalebarlength+5,drawcmd.str().c_str(), temp_BMP_filename, BITMAPDIR, 
		//	 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());

        sprintf(command1,
		"%s -quality %i %s %s%s_proj_%05i.%s", 
        IMAGEMAGICKCONVERT, BMP_COMPRESSION, temp_BMP_filename, BITMAPDIR, 
			 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());
    }
    else 
    {
		sprintf(command1,
		"%s %s -quality %i -font helvetica -fill white -pointsize %i -draw \"text %i %i '%iuM' rectangle %i %i %i %i text +%i+%i '%s-projection' text +%i+%i 'Frame % 6i' text +%i+%i 'Time % 6.1f'\" %s %s%s_proj_%05i.%s", 
        IMAGEMAGICKCONVERT, temp_BMP_filename, BMP_COMPRESSION, 
            20 * BMP_AA_FACTOR, 5 * BMP_AA_FACTOR, 595 * BMP_AA_FACTOR, scalebarmicrons, 
            5 * BMP_AA_FACTOR, 576 * BMP_AA_FACTOR, scalebarlength + 5 * BMP_AA_FACTOR,  573 * BMP_AA_FACTOR,
            5 * BMP_AA_FACTOR, 20 * BMP_AA_FACTOR,  // first line of text
            projletter, 
            5 * BMP_AA_FACTOR, ( 20 + 25 )* BMP_AA_FACTOR , 
            filenum,
			5 * BMP_AA_FACTOR, ( 20 + 50 )* BMP_AA_FACTOR , 
            filenum * InterRecordIterations * DELTA_T, drawcmd.str().c_str(), BITMAPDIR,
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
        sprintf(command2,"%s -quality %i -resize %f%% %s%s_proj_%05i.%s ",
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

void actin::writebitmapheader(ofstream& outbmpfile, const int & bitmapwidth, const int & bitmapheight)
{

	// bitmap headers etc (see microsoft website---large chunks pasted from given code):

	// define data structures for bitmap header:

	#ifdef _WIN32
		#define QUADWORD __int64
	#else 
		#define QUADWORD long long
	#endif

	#define DWORD unsigned int
	#define LONG unsigned int
	#define WORD unsigned short int
	#define BYTE unsigned char
	#define FOURCC unsigned int

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

void actin::writebitmapfile(ofstream& outbmpfile, const Dbl2d& imageR, const Dbl2d& imageG, const Dbl2d& imageB)
{  // re-scale for byte output and save image data

#define BYTE unsigned char

#pragma pack(push,1)  // align the structs to byte boundaries

	typedef struct tagRGB 
	{
		BYTE    B; 
		BYTE    G; 
		BYTE    R; 
	} RGB; 

#pragma pack(pop)


	RGB *line;
	line = new RGB[BMP_WIDTH];

	outbmpfile.seekp(bitmap_start); // move to start

	// write out the data, line by line (note: y is backwards)
	for (int y = (BMP_HEIGHT-1); y>=0; y--)
		{
		//outbmpfile.write(picbuff + (bitmapwidth*y), bitmapwidth);
		for (int x = 0; x<BMP_WIDTH; x++)
			{
			line[x].B=(unsigned char)(255 * imageB[x][y]);
			line[x].G=(unsigned char)(255 * imageG[x][y]);
			line[x].R=(unsigned char)(255 * imageR[x][y]);
			}
			outbmpfile.write((char*)line,BMP_WIDTH*3);
		}

	outbmpfile.flush();

	delete [] line;

}


void actin::squash(const double & thickness)
{   // todo: check and fix this

	// squash with 'coverslip'

    vect rotpos;
    const double halfthickness = thickness/2.0;

	for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin() + highestnodecount;
							      ++i_node)
    {
        rotpos = *i_node;
        rotpos += p_nuc->position;

        inverse_actin_rotation.rotate(rotpos); // rotate

        if (rotpos.x >  halfthickness)  // above
			rotpos.x =  halfthickness * 0.999; 
        else 
        if (rotpos.x < -halfthickness)  // below
			rotpos.x = -halfthickness * 0.999;
        else 
            continue;                   // within coverslip, skip
        
        actin_rotation.rotate(rotpos);  // rotate back
        rotpos -= p_nuc->position;

        i_node->x = rotpos.x;
        i_node->y = rotpos.y;
        i_node->z = rotpos.z;

	}

	return;
}


void actin::sortnodesbygridpoint(void)
{
    // todo: change this to gridpoints-by-thread

//	nodes* nodeptr, *startnodeptr;

    int i;
	int tn;
	size_t minsize,threadsize;

	for (i=lowestnodetoupdate; i != highestnodecount; ++i)
	{
		donenode[i] = false;
	}

    int threadnum;

    for (threadnum = 0; threadnum < NUM_THREAD_DATA_CHUNKS; ++threadnum)
    {
        nodes_by_thread[threadnum].resize(0);
    }

    threadnum = 0;


	for (i=lowestnodetoupdate; i != highestnodecount; ++i)
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

        minsize = MAXNODES;

        for (tn = 0; tn != NUM_THREAD_DATA_CHUNKS; ++tn)
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

	char command1[255];
	// report files
		
	sprintf(command1, "(gzip -q -f -9 %s*report*.txt 2>/dev/null; mv %s*report*.gz %s 2>/dev/null) &"
		,TEMPDIR,TEMPDIR, REPORTDIR);
	system(command1);

	// save data file

	sprintf(command1 , "(gzip -q -f -9 %s*data*.txt 2>/dev/null; mv %s*data*.gz %s 2>/dev/null) &",
		TEMPDIR, TEMPDIR, DATADIR);
	system(command1);

	// wrl file

	sprintf(command1 , "(gzip -f -9 %snodes%05i.wrl 2>/dev/null; mv %snodes%05i.wrl.gz %snodes%05i.wrz 2>/dev/null) &",
						TEMPDIR, filenum, TEMPDIR,filenum, VRMLDIR,filenum );
	system(command1);

#else

	char command2[255];

	//sprintf(command1, "gzip -q -f -9 %s*report*.txt",TEMPDIR);
	sprintf(command2, "move %s*report%05i.txt %s",TEMPDIR,filenum,REPORTDIR);
	//system(command1);
	system(command2);

	// save data file

	//sprintf(command1, "gzip -q -f -9 %s*data*.txt",TEMPDIR2);
	//sprintf(command2, "move %s*data*.gz %s >NUL",TEMPDIR,DATADIR);
	sprintf(command2, "move %s*data*.txt %s >NUL",TEMPDIR,DATADIR);
	//system(command1);
	system(command2);
	//cout << command1 << endl << command2 << endl;

	// wrl file

	//sprintf(command1, "gzip -f -9 %snodes%05i.wrl",TEMPDIR);
	//sprintf(command2, "move %snodes%05i.wrl.gz %snodes%05i.wrz",TEMPDIR,VRMLDIR);
	//system(command1);
	//system(command2);


#endif


}



int actin::find_center(vect &center)
{
	center = - p_nuc->position;

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
       << actin_rotation << " "
	   << camera_rotation << " "
	   << camera_rotation2 << endl;

    
    
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

    // save nucleator
    ofstrm << "nucleator:" << endl;
    p_nuc->save_data(ofstrm);
    
    return 0;
}

int actin::load_data(ifstream &ifstr)
{
	

    // clear the nodegrid
    clear_nodegrid();
    
    string str;
    unsigned char ch;
    
    // load actin
    ifstr >> str;
    // ensure the identifier for the start of the actin
    if(str.compare("actin:") !=0 ){
	cout << "error in checkpoint file, 'actin:' expected" << endl;
	return 1;
    }

	ifstr >> highestnodecount 
	  >> nexttocrosslink  
	  >> iteration_num    
	  >> linksbroken      
	  >> linksformed      
	  >> actin_rotation   
	  >> camera_rotation  
	  >> camera_rotation2;

    inverse_actin_rotation = actin_rotation.inverse();

    ifstr >> str;

    // load nodes
    
    if(str.compare("nodes-links:") !=0 ){
	cout << "error in data file, 'nodes-links:' expected" << endl;
	cout << "'" << str <<"' last read." << endl;
	return 1;
    }
    // ** Remember the node vector is preallocated to MAXNODES
    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
    {
		i_node->load_data(ifstr);
    }

    // only add to grid if doing normal run or a post-process that needs links:
    // REVISIT this is really unexpected behaviour if you add a new post process operation
    // is it a usefule optimisation?
    if (POST_STATS || POST_VTK || (!REWRITESYMBREAK && !POST_PROCESS))
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
    if(str.compare("xlinkdelays:") !=0 ){
	cout << "error in checkpoint file, 'xlinkdelays:' expected" 
	     << endl;
	return 1;
    }
    
    int numcrosslinkdelay;
    ifstr >> numcrosslinkdelay >> ch;
    if(ch!=')' ){
	cout << "error in checkpoint file, xlinkdelays 'NN)' expected" 
	     << endl;
	return 1;
    }
    //crosslinknodesdelay.clear();
    crosslinknodesdelay.resize(numcrosslinkdelay);
    // load each delay
    for(vector<int>::iterator i = crosslinknodesdelay.begin(); 
	i != crosslinknodesdelay.end(); ++i) {
	ifstr >> (*i) ;
    }

    // load nucleator
    ifstr >> str;
    if(str.compare("nucleator:") !=0 ){
	cout << "error in checkpoint file, 'nucleator:' expected" 
	     << endl;
	return 1;
    }
    p_nuc->load_data(ifstr);
    
    return 0;
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
}

void actin::set_sym_break_axes(void)
{
	vect tmp_nodepos;
	//vect CofM;
	vect sym_break_direction;

	rotationmatrix tmp_rotation, tmp_rotation2;	
	rotationmatrix final_rotation;

    tmp_rotation.settoidentity();

	double x_angle = 0, y_angle = 0; 

    if (p_nuc->geometry==nucleator::sphere)
    {   // only rotate by x and y if sphere

	    find_center(sym_break_direction);  // which way did the bead go?

	    x_angle = atan2(sym_break_direction.y,sym_break_direction.z);

	    tmp_rotation.rotatematrix(x_angle, 0 , 0);
	    tmp_rotation.rotate(sym_break_direction);

	    tmp_rotation.settoidentity();

	    y_angle = -atan2(sym_break_direction.x,sym_break_direction.z);

	    tmp_rotation.rotatematrix(0, y_angle , 0); // now pointing upwards
	    tmp_rotation.rotate(sym_break_direction);

        // construct x & y rotation matrix

	    tmp_rotation.settoidentity();

	    tmp_rotation.rotatematrix(y_angle , yaxis);
	    tmp_rotation.rotatematrix(x_angle , xaxis);

    }


	// we now have tmp_rotation which will transform the bead direction to be along the z axis
	// need to determine the z rotation by the principlal axis

	double theta, chi, maxchi, maxchiangle;

	maxchi = 0;
	maxchiangle = 0;

	for(theta = -PI; theta < PI; theta+=PI/180)
	{
		chi = 0;

		tmp_rotation2.settoidentity();
		tmp_rotation2.rotatematrix(0,0,theta);

		for (int i=0; i != highestnodecount; ++i)
		{
			if ((!node[i].polymer) || (node[i].harbinger))
				continue;

			tmp_nodepos = node[i];
			tmp_rotation.rotate(tmp_nodepos);  // rotate points to bring in line with sym break dir'n
			tmp_rotation2.rotate(tmp_nodepos);	// rotate z

			// only look at front of bead (i.e. new z<0, in direction of movement)
			//if (tmp_nodepos.z<0)
			{
				chi += tmp_nodepos.y * tmp_nodepos.y; // sum the squares of distance from the axis
			}
			
		}

		if (chi > maxchi)
		{
			maxchi = chi;
			maxchiangle = theta;
		}

	}

    if (p_nuc->geometry==nucleator::capsule)
    {
        maxchiangle += PI/4;
    }

	// now have all the angles
	// need to assemble them in the right (reverse) order z,y,x
	// to get the rotation matrix

	// camera_rotation includes a pre-rotation with the old x angle
	// to prevent it moving always onto the z-axis
	// camera_rotation2 does not have this

	camera_rotation.settoidentity();


    if (SYM_BREAK_TO_RIGHT)
    {
        camera_rotation.rotatematrix(PI/2,xaxis); // rotate so that always breaks to the right

    } else // restore the component of the sym break in the yz plane
    {
        camera_rotation.rotatematrix(-x_angle, xaxis);
    }

	camera_rotation.rotatematrix(maxchiangle, zaxis);
	camera_rotation.rotatematrix(y_angle, yaxis);
	camera_rotation.rotatematrix(x_angle, xaxis);



	camera_rotation2.settoidentity();

	camera_rotation2.rotatematrix(maxchiangle, zaxis);
	camera_rotation2.rotatematrix(y_angle, yaxis);
	camera_rotation2.rotatematrix(x_angle, xaxis);


    reverse_camera_rotation = camera_rotation.inverse();

	cout << "Symmetry broken.  Camera rotation angles: " << x_angle*180/PI << "," << y_angle*180/PI << "," << maxchiangle*180/PI << endl;

	symbreakiter = iteration_num;

}


void actin::save_sym_break_axes(void)
{
	ofstream opsymbreak(SYM_BREAK_FILE, ios::out | ios::trunc);
	if (!opsymbreak) 
	{ cout << "Unable to open file 'sym_break_axis.txt' for output"; return;}

	opsymbreak  << symbreakiter << endl
				<< camera_rotation << endl
				<< camera_rotation2 << endl;

	opsymbreak.close();

	//cout << "'sym_break_axis.txt' file written" << endl;
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
					>> camera_rotation
					>> camera_rotation2;

		ipsymbreak.close();  
        
        reverse_camera_rotation = camera_rotation.inverse();
        
        return true;
	}
}


void actin::clear_node_stats(void)
{
	attemptedpolrate = polrate	= 0;

	for (int i=lowestnodetoupdate; i<highestnodecount; ++i)
	{
		node[i].clearstats();
	}

}

void actin::keep_mem_resident(void)
{

#ifdef _NUMA

 nmadvise(&nodegrid, sizeof(nodegrid), MADV_WILLNEED, NULL);
    //madvise((caddr_t)&nodegrid, sizeof(nodegrid), MADV_WILLNEED);
	for (int i=0; i!=(GRIDSIZE+1); i++)  
	{
        nmadvise(&nodegrid[i], sizeof(nodegrid[i]), MADV_WILLNEED, NULL);
        //madvise((caddr_t)&nodegrid[i], sizeof(nodegrid[i]), MADV_WILLNEED);
        for (int j=0; j!=(GRIDSIZE+1); j++)
		{
            nmadvise(&nodegrid[i][j], sizeof(nodegrid[i][j]), MADV_WILLNEED, NULL);
            //madvise((caddr_t)&nodegrid[i][j], sizeof(nodegrid[i][j]), MADV_WILLNEED);
        }
	}

#endif
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
    const double BROWNIANFORCESCALE = 0.01;

    const double brownianscale = BROWNIANFORCESCALE * DELTA_T;

    for(vector <nodes>::iterator	i_node  = node.begin(); 
									i_node != node.begin()+highestnodecount;
							      ++i_node)
    {
        i_node->rep_force_vec.x += brownianscale * rand_0to1();
        i_node->rep_force_vec.y += brownianscale * rand_0to1();
        i_node->rep_force_vec.z += brownianscale * rand_0to1();
    }
}
