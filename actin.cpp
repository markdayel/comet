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

actin::actin(void)
{

	opruninfo.open("comet_run_info.txt", ios::out | ios::trunc);
	if (!opruninfo) 
	{ cout << "Unable to open file " << "comet_run_info.txt" << " for output"; return;}

	opruninfo << "comet " << endl << "Mark J Dayel" << endl << "(" << __DATE__ << " " << __TIME__ << ") " << endl << endl;

	cout << endl;
	cout << "comet " << endl << "Mark J Dayel" << endl << "(" << __DATE__ << " " << __TIME__ << ") " << endl << endl;
	cout.flush();

	opvelocityinfo.open("velocities.txt", ios::out | ios::trunc);
	if (!opvelocityinfo) 
	{
	    cout << "Unable to open file " << "velocities.txt" << " for output";
	    return;
	}

	opvelocityinfo << "time,x,y,z,vel" << endl;

	// node = new nodes[MAXNODES];

	node.resize(MAXNODES);
	donenode.resize(MAXNODES);
	repdonenode.resize(MAXNODES);
	nodesbygridpoint.resize(MAXNODES);
	nodesbygridpoint_temp.resize(MAXNODES);
	crosslinknodesdelay.resize(CROSSLINKDELAY);

	doreportiteration = 9999999;

	//reportdat.resize(REPORT_NUM_VARIABLES);
	//for (int i=0; i!=REPORT_NUM_VARIABLES; i++)  
	//{
	//	reportdat[i].resize(REPORT_AVERAGE_ITTERATIONS+1);
	//}


	for (int i=0; i<MAXNODES; i++)
	{
		node[i].nodenum=i;
		node[i].ptheactin=this;
		nodesbygridpoint[i] = i;
		nodesbygridpoint_temp[i] = i;
		if (i<10000)
			node[i].listoflinks.reserve(64);
	}

	cout << "GridExtent : " << GRIDBOUNDS << " uM" << endl;
	cout << "GridResolution : " << GRIDRES << " uM" << endl;
	cout << "TotalGridpoints : " << (GRIDSIZE*GRIDSIZE*GRIDSIZE) << endl;

	cout << "Zeroing Grid...";
	cout.flush();

    nodegrid.resize(GRIDSIZE+1);
	for (int i=0; i!=(GRIDSIZE+1); i++)  // allocate nodegrid and set nodegrid pointers to null
	{
		nodegrid[i].resize(GRIDSIZE+1);
		for (int j=0; j!=(GRIDSIZE+1); j++)
		{
			nodegrid[i][j].resize(GRIDSIZE+1,0);
		}
	}

	cout << "Done" << endl << endl;

	nucleatorgrid.reserve(5000);
	
	// these are reserved 4x the NUM_THREADS because they may run concurrently 
	// between stages and between threads inside stages

	recti_near_nodes.resize(NUM_THREADS*4);
    nodes_on_same_gridpoint.resize(NUM_THREADS*4);

	for (int i = 0; i < NUM_THREADS; ++i)
	{
		recti_near_nodes[i].reserve(5000);
		nodes_on_same_gridpoint[i].reserve(MAXNODES);
	}

	nodes_within_nucleator.reserve(10000);
	linkformto.reserve(100*MAX_LINKS_PER_NODE);

	/*repulsedone.resize(MAXNODES);
	for (int i=0; i!=(MAXNODES); i++)   // this crashes:  why??
	{
		repulsedone[i].resize(MAXNODES);
	}*/

	cout << "Memory for Grid : " << (sizeof(nodes*)*GRIDSIZE*GRIDSIZE*GRIDSIZE/(1024*1024)) << " MB" << endl;

	highestnodecount = 0;
	iteration_num = 0;
	linksformed = 0;
	linksbroken = 0;
	nexttocrosslink = 0;
	symbreakiter = 0;

	newnodescolour.setcol(0);

	brokensymmetry = false;

	//debug:

	debug_num_rotate = 0;
	debug_num_displace = 0;

	int x,y;

	// bitmap size:
	//BMP_WIDTH = 800;
	//BMP_HEIGHT = 600;

	imageR.resize(BMP_WIDTH);
	imageG.resize(BMP_WIDTH);
	imageB.resize(BMP_WIDTH);

	for (x = 0; x<BMP_WIDTH; x++)
		{
			imageR[x].resize(BMP_HEIGHT);
			imageG[x].resize(BMP_HEIGHT);
			imageB[x].resize(BMP_HEIGHT);

			for (y = 0; y<BMP_HEIGHT; y++)
			{
				imageR[x][y]=0;
				imageG[x][y]=0;
				imageB[x][y]=0;
			}
		}

//GetProcessId(GetCurrentThread());

	// setup the bitmap file:

	sprintf(temp_BMP_filename_x, "%sx_temp_%i.bmp", TEMPDIR, (unsigned)time( NULL ));
	sprintf(temp_BMP_filename_y, "%sy_temp_%i.bmp", TEMPDIR, (unsigned)time( NULL ));
	sprintf(temp_BMP_filename_z, "%sz_temp_%i.bmp", TEMPDIR, (unsigned)time( NULL ));

	outbmpfile_x.open(temp_BMP_filename_x, ios::out | ios::binary | ios::trunc);
	outbmpfile_y.open(temp_BMP_filename_y, ios::out | ios::binary | ios::trunc);
	outbmpfile_z.open(temp_BMP_filename_z, ios::out | ios::binary | ios::trunc);

	if ((!outbmpfile_x) || (!outbmpfile_y) || (!outbmpfile_z)) 
	{ cout << "Unable to open file temp bitmap file for output"; return;}

	writebitmapheader(outbmpfile_x,BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_y,BMP_WIDTH, BMP_HEIGHT);
	writebitmapheader(outbmpfile_z,BMP_WIDTH, BMP_HEIGHT);

} 

actin::~actin(void)
{
	opruninfo.close();
	opvelocityinfo.close();
	outbmpfile_x.close();
	outbmpfile_y.close();
	outbmpfile_z.close();

	// delete the temp bitmap file

	char command1[255];

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_x );
	system(command1);

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_y );
	system(command1);

	sprintf(command1, "rm -f %s 2>/dev/null", temp_BMP_filename_z );
	system(command1);
}

int actin::nucleate()
{
	int numnewnodes;

	// add new nodes
	numnewnodes = p_nuc->addnodes();
	
return numnewnodes;
}

int actin::crosslinknewnodes(int numnewnodes)
{	
	int threadnum = 0;
	
	vect nodeposvec;
	double distsqr;

	vect disp;

	if (numnewnodes == 0)
		return 0;

	// and link the new nodes...
	for (int i=nexttocrosslink; i < (nexttocrosslink+numnewnodes); i++)
	{
		if (!node[i].harbinger)
			continue;  // look for harbingers to crosslink

		//cout << "linking node " << i << endl;
		//cout << "created node at " << node[i].x << "," << node[i].y << "," << node[i].z << endl;

		if (findnearbynodes(node[i],NODE_XLINK_GRIDSEARCH,threadnum)==0) continue;	// find nodes within grid point
												// skip if zero
							 
		// of these nodes, calculate euclidian dist

		nodeposvec = node[i];

		linkformto.resize(0);

		for (vector <nodes*>::iterator nearnode=recti_near_nodes[threadnum].begin(); nearnode<recti_near_nodes[threadnum].end() ; nearnode++ )
		{

			if ( ((*nearnode)->harbinger) || (!(*nearnode)->polymer) )
				continue;  // only crosslink to real nodes

			disp = (*(*nearnode)) - nodeposvec;

			distsqr = disp.sqrlength();

			if (distsqr < SQRT_ACCURACY_LOSS)
				continue;

			if ( distsqr < XLINK_NODE_RANGE*XLINK_NODE_RANGE)
			{
				linkformto.push_back(linkform( (*nearnode)->nodenum , distsqr ));
				//(*(linkformto.end()-1)).nodenum = (*nearnode)->nodenum;
				//(*(linkformto.end()-1)).distsqr = distsqr;

				//addlinks(i,(*nearnode)->nodenum,distsqr);  // add link if within the node link range
			}
		}

		if (XLINK_NEAREST)  // crosslink in order of distance?
			sort(linkformto.begin(), linkformto.end(), linkform::CompareDistance);
		else				// or in random order
			random_shuffle(linkformto.begin(),linkformto.end());

		// and crosslink:
		for (vector <linkform>::iterator linkto=linkformto.begin(); linkto<linkformto.end() ; linkto++ )
			{
				//if (i%100==0)
				//{
				//	cout << (*linkto).distsqr << " " << node[i].listoflinks.size() 
				//		<< " " << node[(*linkto).nodenum].listoflinks.size() << endl;
				//}
				if (node[i].listoflinks.size()<MAX_LINKS_PER_NODE)
				{
					if (node[(*linkto).nodenum].listoflinks.size()<MAX_LINKS_PER_NODE)
					{
						addlinks(i,(*linkto).nodenum);
					}
				}
				else
					break;
			}

		node[i].harbinger = false;  // crosslinked: now node exists

	}

	nexttocrosslink += numnewnodes;

	return numnewnodes;
}
//
//int actin::save(int filenum)
//{
//	char filename[255];
//	//char time[255], date[255];
//
//    //_strtime( time );
//    //_strdate( date );
//
//	sprintf ( filename , "nodes%05i.txt", filenum );
//
//	ofstream opactindata(filename, ios::out | ios::trunc);
//	if (!opactindata) 
//	{ cout << "Unable to open file " << filename << " for output"; return 1;}
//
//	// write header
//	opactindata << "comet (compiled " << __DATE__ << " " << __TIME__ << ") Mark J Dayel" << endl;
//	//opactindata << "Data Output " << date << " " << time << endl;
//	opactindata << "Nodes Data Itteration Number: " << filenum << endl;
//	opactindata << "Highest Node Count: " << highestnodecount << endl;
//		opactindata << "nodenum,x,y,z,harbinger,polymer,r,g,b,creationiter,numlinks " << endl;
//	
//	p_nuc->save(&opactindata);
//
//
//
//	// write data
//	for (int i=0; i<highestnodecount; i++)
//	{
//		if (node[i].polymer)
//		{
//			node[i].save(&opactindata);
//		}
//	}
//
//	opactindata.close();
//
//
//
//	// now save links
//
//	sprintf ( filename , "links%05i.txt", filenum );
//
//	ofstream oplinkdata(filename, ios::out | ios::trunc);
//
//	if (!oplinkdata) 
//	{ cout << "Unable to open file " << filename << " for output"; return 1;}
//
//	// write header
//	oplinkdata << "comet (compiled " << __DATE__ << " " << __TIME__ << ") Mark J Dayel" << endl;
//	//oplinkdata << "Data Output " << date << " " << time << endl;
//	oplinkdata << "Links Data Itteration Number: " << filenum << endl;
//	oplinkdata << "Highest Node Count: " << highestnodecount << endl;
//	oplinkdata << "nodenum,numlinks,linknum1,linknum2,..." << endl;
//
//	for (int i=0; i<highestnodecount; i++)
//	{
//		if (node[i].polymer)
//		{
//			node[i].savelinks(&oplinkdata);
//		}
//	}
//	oplinkdata.close();
//
//	return 0;
//}

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

	minx = miny = minz =   GRIDBOUNDS;
	maxx = maxy = maxz = - GRIDBOUNDS;
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
				opvrml	<< node[i].colour.r << " " << node[i].colour.g << " " << node[i].colour.b;
			else
				opvrml	<< "," << node[i].colour.r << " " << node[i].colour.g << " " << node[i].colour.b;
	}



	opvrml << "] }" << endl;
	opvrml << "    }" << endl;
	opvrml << "}" << endl;

	opvrml.close();



	return 0;
}


int actin::iterate()  // this is the main iteration loop call
{
	int numnewnodes;

	numnewnodes = nucleate();
	
	// crosslink, but only after time to equilibrate position
	// in mean time it is a 'harbinger'
	// affected only by node collision repulsion and nucleator repulsion

	crosslinknodesdelay[iteration_num%CROSSLINKDELAY] = numnewnodes;
	crosslinknewnodes(crosslinknodesdelay[(iteration_num+1)%CROSSLINKDELAY]);
	
	sortnodesbygridpoint();	// sort the nodes so they can be divided sanely between threads

	collisiondetection();	// calc node-to-node repulsion
	linkforces();			// and link forces

	ejectfromnucleator();	// do forcable node ejection
	
	move_and_rotate();		// move and rotation of the actin based on 
							// nucleator ejection data
							// Note: whole grid update (in applyforces()) 
							// must come immediately after move/rotate

	applyforces(); // move the nodes by the forces and update the grid

	iteration_num++;

	return 0;
}

int actin::addlinks(const int& linknode1,const int& linknode2)
{
	// crosslink a new node
	// returns 1 if link added, 0 if not

	double pxlink;

	if (linknode1==linknode2) return 0;

	if ((!node[linknode1].polymer)  || (!node[linknode2].polymer) ||
		(!node[linknode1].harbinger))
		return 0;

	//if ((linknode1==178) || (linknode1==180) ||
	//	(linknode2==178) || (linknode2==180) )
	//	cout << "stop here";

	//double dist = calcdist(node[linknode1].x-node[linknode2].x,
	//						 node[linknode1].y-node[linknode2].y,
	//						 node[linknode1].z-node[linknode2].z);

	double dist = calcdist(node[linknode1],node[linknode2]);

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
	
	pxlink = 1;

	if ( (P_XLINK*pxlink) > ( ((double) rand()) / (double)RAND_MAX ) )
	{
		node[linknode1].addlink(&node[linknode2],dist); 
		node[linknode2].addlink(&node[linknode1],dist);
        
		//cout << "linked " << node[linknode1].nodenum << " to " << node[linknode2].nodenum  << endl;
	
		return 1;
	}

	return 0;
}

void actin::ejectfromnucleator()
{
	vect nodeposvec, oldnucposn;
	nodes *nodeptr, *startnodeptr;

//  check node-nucleator repulsion
// collect the nodes:

	nodes_within_nucleator.resize(0);

	for (vector <int_vect>::iterator i=nucleatorgrid.begin(); i<nucleatorgrid.end() ; i++ )
	{	 
		nodeptr=nodegrid[i->x][i->y][i->z];
		startnodeptr=nodeptr;
			if (nodeptr!=0) 
			{
				do
				{
					if ((p_nuc->iswithinnucleator(nodeptr->x,nodeptr->y,nodeptr->z))
						&& (nodeptr->polymer))
					{	// inside nucleator
						nodes_within_nucleator.push_back(nodeptr);						
					}

					nodeptr = nodeptr->nextnode;
				}
				while (nodeptr!=startnodeptr);  //until back to start
				
			}
	}

	oldnucposn = p_nuc->position;

	// do the ejection:

	for (vector <nodes*>::iterator i=nodes_within_nucleator.begin(); i<nodes_within_nucleator.end() ; i++ )
	{  //eject them
		if (p_nuc->collision(*(*i))) // (*i)->x,(*i)->y,(*i)->z)==0)  
			(*i)->updategrid();  // ejected OK
		else
			(*i)->depolymerize();  // not ejected OK, depolymerize
	}

}

void actin::move_and_rotate()
{

	// do rotation

//bool toupdategrid = false;


	if (ROTATION && (p_nuc->torque.length() > MIN_TORQUE_TO_UPDATE))
	{

		
		debug_num_rotate++;

		double x_angle = p_nuc->torque.x / 
				p_nuc->momentofinertia.x;

		double y_angle = p_nuc->torque.y / 
				p_nuc->momentofinertia.y;

		double z_angle = p_nuc->torque.z / 
				p_nuc->momentofinertia.z;

		//double x_angle = 0;  // for debugging
		//double y_angle = 0;
		//double z_angle = 0.00002;

		rotationmatrix torque_rotate; // ( x_angle, rotationmatrix::xaxis);

		torque_rotate.rotatematrix( x_angle, y_angle, z_angle);


		// rotate the actin:

		for (int i=0; i<highestnodecount; ++i)
		{
			torque_rotate.rotate(node[i]);
		}

        // rotate the nucleator:

		p_nuc->nucleator_rotation.rotatematrix( -x_angle, -y_angle, -z_angle);


		// rotate the actin reference frame:

		actin_rotation.rotatematrix( x_angle, y_angle, z_angle);


		// clear the torque vectors:

		p_nuc->torque.zero();

		//toupdategrid = true;

	}


// do displacment

	if ((p_nuc->deltanucposn.length() > MIN_DISPLACEMENT_TO_UPDATE))
	{

		debug_num_displace++;

		// move the nucleator by moving nodes in opposite direction
		// note this is *before* rotating the deltanucposn vector, since the
		// nodes are in the nucleator reference frame
		
		for (int i = 0; i<highestnodecount; i++)
		{
			node[i]-=p_nuc->deltanucposn;
		}


		// rotate the nucleator displacement vector

		p_nuc->nucleator_rotation.rotate(p_nuc->deltanucposn);

		// update the nucleator position with the rotated vector

		p_nuc->position+=p_nuc->deltanucposn;

		// and zero

		p_nuc->deltanucposn.zero();

		//toupdategrid = true;
	}

	return;
}

int actin::collisiondetection(void)
{

	for (int i=0; i<highestnodecount; i++)
	{
		donenode[i] = false;
	}

	
	int numthreadnodes, start, end, higestorderednode;

	higestorderednode = highestnodecount; // (int) nodesbygridpoint.size();

	numthreadnodes = higestorderednode / NUM_THREADS;

// do collision detection

//int numtodo = 0;

	if (USE_THREADS)
	{

		for (int i = 0; i < NUM_THREADS; i++)
		{
			start = i * numthreadnodes;
			end = (i+1) * numthreadnodes;

			collision_thread_data_array[i].startnode = start;
			collision_thread_data_array[i].endnode = end;
			collision_thread_data_array[i].threadnum = i;

			if (i<NUM_THREADS-1)
			{
				collision_thread_data_array[i].endnode = end;
				//numtodo+=end-start;
			}
			else
			{	// put remainder in last thread (cludge for now)
				collision_thread_data_array[i].endnode = end + higestorderednode % NUM_THREADS;
				//numtodo+=end-start + highestnodecount % NUM_THREADS;
			};			
		}
// check
		/*
		if (rand() < RAND_MAX * 0.01)
		{
			for (int i = 0; i < NUM_THREADS; i++)
			{
				cout << collision_thread_data_array[i].startnode << " to " <<
					collision_thread_data_array[i].endnode << "...";
			}
			cout << "higestorderednode " << higestorderednode << endl;
		}
*/


	// check
		/*if (numtodo != highestnodecount)
		{
			cout << "numtodo not match :" << numtodo << " should be" << highestnodecount << endl;
			for (int i = 0; i < NUM_THREADS; i++)
			{
				cout << "Thread " << i << " start: " << collision_thread_data_array[i].startnode
					<< " end " << collision_thread_data_array[i].endnode << endl;
			}
			cout.flush();
		}*/

// thread test: call sequentially within this thread

		for (int i = 0; i < NUM_THREADS; i++)
		{
			collisiondetectiondowork(&collision_thread_data_array[i]);
		}


/*
		// release thread blocks to start threads:

	for (int i = 0; i < NUM_THREADS; i++)
		{
			//sem_post(&collision_thread_go[i]);
			// start the thread:
			pthread_mutex_unlock(&collisiondetectiongolock_mutex[i]);
		}


	// re-sync when done:
	for (int i = 0; i < NUM_THREADS; i++)
		{
			// grab done locks for all threads:
			pthread_mutex_lock(&collisiondetectiondonelock_mutex[i]);
		}

	// relock 'go'

	for (int i = 0; i < NUM_THREADS; i++)
		{
			// grab done locks for all threads:
			pthread_mutex_lock(&collisiondetectiongolock_mutex[i]);
		}

	// and let children take the 'done' locks again

	for (int i = 0; i < NUM_THREADS; i++)
		{
			// grab done locks for all threads:
			pthread_mutex_unlock(&collisiondetectiondonelock_mutex[i]);
		}
*/
/*
if (!collisionthreaddone1)
	cout << "Thread 0 still going!" << endl;
if (!collisionthreaddone2)
	cout << "Thread 1 still going!" << endl;
if (!collisionthreaddone3)
	cout << "Thread 2 still going!" << endl;
if (!collisionthreaddone4)
	cout << "Thread 3 still going!" << endl;
*/
	}
	else
	{  // if not using threads, do in one go:
		collision_thread_data_array[0].startnode = 0;
		collision_thread_data_array[0].endnode = higestorderednode;
		collision_thread_data_array[0].threadnum = 0;
		collisiondetectiondowork(&collision_thread_data_array[0]);
	}





	return 0;
}
/*
void * actin::collisiondetectionthread(void* threadarg)
{

	struct thread_data *dat;
	dat = (struct thread_data *) threadarg;

	cout << "Starting collision detection thread " << dat->threadnum << endl;

	while (true)
	{	
		pthread_mutex_lock(&collisiondetectiondonelock_mutex[dat->threadnum]);  // go only if released by main thread
		pthread_mutex_lock(&collisiondetectiongolock_mutex[dat->threadnum]);  // lock main thread
		// release the go lock so can be taken by main thread again:
		pthread_mutex_unlock(&collisiondetectiongolock_mutex[dat->threadnum]);

		if (dat->threadnum == 0)
			collisionthreaddone1 = false;
		if (dat->threadnum == 1)
			collisionthreaddone2 = false;
		if (dat->threadnum == 2)
			collisionthreaddone3 = false;
		if (dat->threadnum == 3)
			collisionthreaddone4 = false;

//	cout << "Thread running " << dat->threadnum << endl;

		collisiondetectiondowork(dat);

		if (dat->threadnum == 0)
			collisionthreaddone1 = true;
		if (dat->threadnum == 1)
			collisionthreaddone2 = true;
		if (dat->threadnum == 2)
			collisionthreaddone3 = true;
		if (dat->threadnum == 3)
			collisionthreaddone4 = true;


		pthread_mutex_unlock(&collisiondetectiondonelock_mutex[dat->threadnum]);  // release main thread block waiting for data
	}  

	return (void*) NULL;
}
*/

void * actin::collisiondetectiondowork(thread_data* dat)
{
	vect nodeposvec;
	double distsqr,dist;
	double force;
	vect tomove;

	vect disp;

	for (int i=dat->startnode; i<dat->endnode; i++)
	{	
		if ((node[nodesbygridpoint[i]].dontupdate) ||
			(donenode[nodesbygridpoint[i]]) ||
			(!node[nodesbygridpoint[i]].polymer))
		{
			continue;  // skip nodes already done
		}
		
		// find nodes on same gridpoint, and nodes within repulsive range

		if (findnearbynodes(node[nodesbygridpoint[i]],NODE_REPULSIVE_RANGE_GRIDSEARCH,dat->threadnum)==0) continue;	// find nodes within 1 grid point
												
		// skip if zero

		// loop over nodes on same gridpoint:
		for (vector <nodes*>::iterator sameGPnode=nodes_on_same_gridpoint[dat->threadnum].begin(); sameGPnode<nodes_on_same_gridpoint[dat->threadnum].end() ; sameGPnode++ )
		{

			if (donenode[(*sameGPnode)->nodenum])
				continue;  // skip if done
			
			if (( (*sameGPnode)->nodenum < dat->startnode) ||
				( (*sameGPnode)->nodenum >= dat->endnode) )
				 continue;    // skip if out of range of this thread

			donenode[(*sameGPnode)->nodenum] = true;  // mark node as done

			// of these nodes, calculate euclidian dist

			nodeposvec = (*(*sameGPnode));	// get xyz of our node

			for (vector <nodes*>::iterator nearnode=recti_near_nodes[dat->threadnum].begin(); nearnode<recti_near_nodes[dat->threadnum].end() ; nearnode++ )
			{
				if ((*sameGPnode)==(*nearnode)) continue;  // skip if self
				//if (repulsedone[(*nearnode)->nodenum][(*sameGPnode)->nodenum])
				//	continue;  // skip if done opposite pair
				
				// this check should not be necessary?:

				if (( (*sameGPnode)->nodenum < dat->startnode) ||
					( (*sameGPnode)->nodenum >= dat->endnode) )
					continue;    // skip if out of range of this thread  

				
				disp = (*(*nearnode)) - nodeposvec; 

				distsqr = disp.sqrlength();

				if (distsqr < 2 * SQRT_ACCURACY_LOSS)
				{	
					continue;
				}

				//dorepulsion((*sameGPnode)->nodenum,(*nearnode)->nodenum,dist,dat->threadnum);
				//dorepulsion(*(*sameGPnode),*(*nearnode),dist,dat->threadnum);
				
				if ( distsqr < NODE_REPULSIVE_RANGE*NODE_REPULSIVE_RANGE)
					{
						// calc dist between nodes

						dist = sqrt(distsqr); 

						force =  NODE_REPULSIVE_MAG - (NODE_REPULSIVE_MAG / NODE_REPULSIVE_RANGE) * dist;

					//if (dist > NODE_INCOMPRESSIBLE_RADIUS)
					//{
					//	force =  NODE_REPULSIVE_MAG - (NODE_REPULSIVE_MAG / 
					//		(NODE_REPULSIVE_RANGE - NODE_INCOMPRESSIBLE_RADIUS) * (dist-NODE_INCOMPRESSIBLE_RADIUS )) ;
					//}
					//else
					//{
					//	force =  NODE_REPULSIVE_MAG ;
					//}  // fix force if below node incompressible distance

					tomove = disp * ( 2 * force / dist);
					
					(*sameGPnode)->rep_force_vec[0] -= tomove ;
//
//#ifndef NO_CALC_STATS
//
					(*sameGPnode)->adddirectionalmags(tomove, (*sameGPnode)->repforce_radial[0], (*sameGPnode)->repforce_transverse[0]);
//

//#endif

//					if ( dist < NODE_INCOMPRESSIBLE_RADIUS)
//					{
//					                            
//						// how far to repulse (quarter for each node) (half each, but we will do i to j and j to i)
//						scale = (double)NODE_INCOMPRESSIBLE_RADIUS / (double)4.0*dist;  
//
//						tomove = disp * (2*scale);	// these are the vector half components (in direction from i to j)
//
//
//						(*sameGPnode)->repulsion_displacement_vec[0] = -tomove;
//
////#ifndef NO_CALC_STATS
////
////						(*sameGPnode)->adddirectionalmags(tomove, (*sameGPnode)->dispforce_radial[0], (*sameGPnode)->dispforce_transverse[0]);
////
////#endif
//					}
				}

			}
		}
	}
return (void*) NULL;
}

//int actin::findnearbynodes(const int& ournodenum, const int& adjgridpoints, const int& threadnum)
int actin::findnearbynodes(const nodes& ournode, const int& adjgridpoints, const int& threadnum)
{	// create list of nodes on the same gridpoint, and on adjacent gridpoints (including same GP)
	// 	
	// N.B. a good fraction of the total CPU time is spent in this function
	// may be worth linearizing the loop or poss amalgamating with the calling
	// functions so that creating and reading the list of gridpoints is not necessary

nodes *nodeptr, *startnodeptr;

// save repeatedly dereferencing pointer by threadnum in inner loops:
// (maybe compiler does this anyway?)
Nodes1d *p_recti_near_nodeslist = &recti_near_nodes[threadnum];
Nodes1d *p_nodes_on_same_gridpoint = &nodes_on_same_gridpoint[threadnum];

p_recti_near_nodeslist->resize(0);
p_nodes_on_same_gridpoint->resize(0); 

if (!ournode.polymer)
		return 0;

int gridx = ournode.gridx;
int gridy = ournode.gridy;
int gridz = ournode.gridz;

if ((gridx < 0) || (gridx > GRIDSIZE) ||
	(gridy < 0) || (gridy > GRIDSIZE) ||
	(gridz < 0) || (gridz > GRIDSIZE))
		return 0;

int minx,miny,minz,maxx,maxy,maxz;

minx = gridx-adjgridpoints;
miny = gridy-adjgridpoints;
minz = gridz-adjgridpoints;

maxx = gridx+adjgridpoints;
maxy = gridy+adjgridpoints;
maxz = gridz+adjgridpoints;

// truncate if out of grid bounds:

if (minx < 0) minx = 0;
if (miny < 0) miny = 0;
if (minz < 0) minz = 0;

if (maxx > GRIDSIZE) maxx = GRIDSIZE;
if (maxy > GRIDSIZE) maxy = GRIDSIZE;
if (maxz > GRIDSIZE) maxz = GRIDSIZE;

// do adjgridpoints by adjgridpoints scan on grid

for (int x = minx; x != maxx; ++x)  
	for (int y = miny; y != maxy; ++y)
		for (int z = minz; z != maxz; ++z)
		{
		nodeptr=nodegrid[x][y][z];
		startnodeptr=nodeptr;
			if (nodeptr!=0) 
			{
				do
				{
					p_recti_near_nodeslist->push_back(nodeptr);
					nodeptr = nodeptr->nextnode;					
				}
				while (nodeptr!=startnodeptr);  //until back to start
				
			}
		}

	// nodes on same gridpoint:

	nodeptr=nodegrid[gridx][gridy][gridz];
	startnodeptr=nodeptr;
	if (nodeptr!=0) 
	{
		do
		{
			p_nodes_on_same_gridpoint->push_back(nodeptr);
			nodeptr = nodeptr->nextnode;					
		}
		while (nodeptr!=startnodeptr);  //until back to start
		
	}

	return (int) p_recti_near_nodeslist->size();
}

//inline int actin::dorepulsion(nodes& node_i,nodes& node_j,
//							  const double& dist,const int& threadnum)
//{						// defunct, now in collisiondetectiondowork
//	if (&node_i==&node_j)
//		return 0;
//
//	double scale;
//
//	vect disp;
//
//	disp = node_j - node_i;
//
//	scale = (NODE_INCOMPRESSIBLE_RADIUS / 4)*dist;  
//
//	disp *= scale;
//
//
//#ifdef FORCES_BOTH_WAYS
//
//	node_i.repulsion_displacement_vec[threadnum] -= disp;
//	node_j.repulsion_displacement_vec[threadnum] += disp;
//
//#else
//
//	node_i.repulsion_displacement_vec[threadnum] -= disp*2;
//
//#endif
//
//
//	return 0;
//}

int actin::applyforces(void)  // this just applys previously calculated forces (in dorepulsion())
{
#ifndef SEED_INSIDE
	int numthreadnodes, start, end;

	numthreadnodes = highestnodecount / NUM_THREADS;

if (!GRASS_IS_GREEN) // (USE_THREADS)  // switch this off for now...
	{
		for (int i = 0; i < NUM_THREADS; i++)
		{
			start = i * numthreadnodes;
			end = (i+1) * numthreadnodes;

			applyforces_thread_data_array[i].startnode = start;
			applyforces_thread_data_array[i].endnode = end;
			applyforces_thread_data_array[i].threadnum = i;

			if (i<NUM_THREADS-1)
			{
				applyforces_thread_data_array[i].endnode = end;
			}
			else
			{	// put remainder in last thread (cludge for now)
				applyforces_thread_data_array[i].endnode = end + highestnodecount % NUM_THREADS;
			};
		}

		// release thread blocks to start threads:
	for (int i = 0; i < NUM_THREADS; i++)
		{
			sem_post(&applyforces_thread_go[i]);
		}

		// and wait 'till complete:
	for (int i = 0; i < NUM_THREADS; i++)
		{
			sem_wait(&applyforces_data_done[i]);
		}

	}
	else
	{
		for (int i=0; i<highestnodecount; i++)
			{
				for (int threadnum = 0; threadnum < NUM_THREADS; threadnum++)
				{
					if ((!node[i].dontupdate) && node[i].polymer)
						node[i].applyforces(threadnum);
				}
			}
	}

	for (int i=0; i<highestnodecount; i++)
	{
		if ((!node[i].dontupdate) && node[i].polymer)
			node[i].updategrid(); // move the point on the grid if need to
								  // and update the unit vector position
	}
#endif
	return 0;
}

void * actin::applyforcesthread(void* threadarg)
{
	struct thread_data *dat;
	dat = (struct thread_data *) threadarg;

	//cout << "Thread " << dat->threadnum << " started" << endl;
	//cout.flush();

	while (GRASS_IS_GREEN)
	{	

		sem_wait(&applyforces_thread_go[dat->threadnum]);  // go only if released by main thread

		//applyforcesdetectiondowork(dat);

		// sum forces etc. over all threads
		for (int i=dat->startnode; i<dat->endnode; i++)
			{
			for (int threadnum = 0; threadnum < NUM_THREADS; threadnum++)
			{
				node[i].applyforces(threadnum);
			}
			//node[i].updategrid(); // move the point on the grid if need to
		}

		sem_post(&applyforces_data_done[dat->threadnum]);  // release main thread block waiting for data

	}  

	return (void*) NULL;
}

int actin::linkforces()
{
	double dist;
	double force;
	vect nodeposvec, disp;
	vect tomove;
	
	// remove the links for ones that were broken last time
	for (unsigned int i=0; i<linkremovefrom.size() ; i++ )
	{
		linkremovefrom[i]->removelink(linkremoveto[i]);  // remove the back link
		linkremoveto[i]->removelink(linkremovefrom[i]);  // and remove from the list
	}

	// reset the lists of broken links:
	linkremovefrom.resize(0);
	linkremoveto.resize(0);
	
	int numthreadnodes, start, end;

	numthreadnodes = highestnodecount / NUM_THREADS;

if (!GRASS_IS_GREEN) // USE_THREADS)
	{
		for (int i = 0; i < NUM_THREADS; i++)
		{
			start = i * numthreadnodes;
			end = (i+1) * numthreadnodes;

			linkforces_thread_data_array[i].startnode = start;
			linkforces_thread_data_array[i].endnode = end;
			linkforces_thread_data_array[i].threadnum = i;

			if (i<NUM_THREADS-1)
			{
				linkforces_thread_data_array[i].endnode = end;
			}
			else
			{	// put remainder in last thread (cludge for now)
				linkforces_thread_data_array[i].endnode = end + highestnodecount % NUM_THREADS;
			};
		}

	// release thread blocks to start threads:
	for (int i = 0; i < NUM_THREADS; i++)
		{
			sem_post(&linkforces_thread_go[i]);
		}

	}
	else
	{
		// go through all nodes
		for (int n=0; n<highestnodecount; n++)
		{  	
			if ((node[n].dontupdate) || (!node[n].polymer))
				continue;

			nodeposvec = node[n];

			// go through links for each node
			for (vector <links>::iterator i=node[n].listoflinks.begin(); i<node[n].listoflinks.end() ; i++ )
			{	 
				if (!(i->broken))  // if link not broken  (shouldn't be here if broken anyway)
				{			
					
					disp = nodeposvec - *(i->linkednodeptr);

					dist = disp.length();;

					//if (dist < 0)
					//	continue;

					force = i->getlinkforces(dist);

					if (i->broken)
					{	// broken link: store which ones to break:
						linkremovefrom.push_back(&node[n]);
						linkremoveto.push_back(i->linkednodeptr);
						node[n].links_broken[0]++;
					}
					else
					{
						tomove = disp * (2* force/dist);  
												  // we're scaling the force by vector disp
												  // to get the direction, so we need to
												  // divide by the length of disp (i.e. dist)
												  // to prevent the length amplifying the force

//#ifdef FORCES_BOTH_WAYS
//						node[n].link_force_vec[0] += tomove/2;
//						i->linkednodeptr->link_force_vec[0] -= tomove/2;
//#else
						node[n].link_force_vec[0] += tomove;

//#ifndef NO_CALC_STATS
//
						if (force < 0) // put tension into link forces
						{
							node[n].adddirectionalmags(tomove, node[n].linkforce_radial[0], node[n].linkforce_transverse[0]);
						}
						else	// but put compression into link forces
						{
							node[n].adddirectionalmags(tomove, node[n].repforce_radial[0], node[n].repforce_transverse[0]);
						}



//						//node[n].linkforce_radial[0]     += fabs(  node[n].unit_vec_posn.dot(disp)  * (2 * force * scale) ) ;
//
//						//node[n].linkforce_transverse[0] += fabs( (node[n].unit_vec_posn.cross(disp)).length()  * (2 * force * scale));
//
//						//cout << node[n].linkforce_radial[0] << " " <<  node[n].linkforce_transverse[0] << endl;
//#endif
						

					}
				}
			}
		}
	}



	return 0;
}

/*
void * actin::linkforcesthread(void* threadarg)
{
	struct thread_data *dat;
	dat = (struct thread_data *) threadarg;

	double xdist, ydist, zdist, dist;
	double force;
	vect nodeposvec;

	//cout << "Thread " << dat->threadnum << " started" << endl;
	//cout.flush();

	while (true)
	{	

		sem_wait(&linkforces_thread_go[dat->threadnum]);  // go only if released by main thread

		for (int n=dat->startnode; n<dat->endnode; n++)
			{  	
				if (!node[n].polymer)
					continue;
				
				nodeposvec.x = node[n].x;	// get xyz of our node
				nodeposvec.y = node[n].y;
				nodeposvec.z = node[n].z;

				// go through links for each node
				for (vector <links>::iterator i=node[n].listoflinks.begin(); i<node[n].listoflinks.end() ; i++ )
				{	 
					if (!(i->broken))  // if link not broken  (shouldn't be here if broken anyway)
					{			
						xdist = nodeposvec.x - i->linkednodeptr->x;
						ydist = nodeposvec.y - i->linkednodeptr->y;
						zdist = nodeposvec.z - i->linkednodeptr->z;

						dist = calcdist(xdist,ydist,zdist);

						if (dist < 0)
							continue;

						force = i->getlinkforces(dist);

						if (i->broken)
						{	// broken link
							pthread_mutex_lock(&linkstoremove_mutex);  // lock the mutex
							linkremovefrom.push_back(&node[n]);
							linkremoveto.push_back(i->linkednodeptr);
							pthread_mutex_unlock(&linkstoremove_mutex);  //unlock
						}
						else
						{
							node[n].link_force_vec[dat->threadnum].x += force * (xdist/dist);
							node[n].link_force_vec[dat->threadnum].y += force * (ydist/dist);
							node[n].link_force_vec[dat->threadnum].z += force * (zdist/dist);
						}
					}
				}
			}

		sem_post(&linkforces_data_done[dat->threadnum]);  // release main thread block waiting for data

	}  

	return (void*) NULL;
}
*/

int actin::setnodecols(void)
{

	double maxcol,mincol;
	maxcol = 0;
	mincol = 0;	

	double val;

	maxcol = (double) 0.001;

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].link_force_vec[0].x>0.0001)
			val = sqrt(node[i].linkforce_transverse[0]);
		else
			val = 0;
	
		if ((node[i].polymer) && (!node[i].listoflinks.empty()) && (val > maxcol))
		{
			//maxcol = node[i].nodelinksbroken;
			maxcol = val;
		}
	}

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].link_force_vec[0].x>0.0001)
			val = sqrt(node[i].linkforce_transverse[0]);
		else
			val = 0;
		//node[i].colour.setcol((double)node[i].creation_iter_num/(double)TOTAL_ITERATIONS);
		//node[i].colour.setcol(((double)node[i].nodelinksbroken-mincol)/(maxcol-mincol));
		if ((node[i].polymer)&& (!node[i].listoflinks.empty()) && val > 0.001)
			{
				node[i].colour.setcol((val-mincol)/(maxcol-mincol));
			}
		else
		{
			node[i].colour.setcol(0);
		}
	}


	//for (int i=0; i<highestnodecount; i++)
	//{
	//node[i].colour.setcol((double)node[i].creation_iter_num/(double)iteration_num);
	//}
/*
	double maxcol,mincol;
	maxcol = 1;
	mincol = 0;	
	double val;

	for (int i=0; i<highestnodecount; i++)
	{
		val = (double)node[i].nodelinksbroken;
		if ((node[i].polymer) && (val > maxcol))
		{
			maxcol = val;
		}
	} 

	for (int i=0; i<highestnodecount; i++)
	{
	node[i].colour.setcol((double)node[i].nodelinksbroken/maxcol);
	}

	for (int i=0; i<highestnodecount;i++)
	{
		//i = (int) ((double) highestnodecount * ( (double) rand()) / (double)RAND_MAX);
		if (i%300==0)
		{
		node[i].colour.setcol((double)i/(double)highestnodecount);
		for (vector <links>::iterator l=node[i].listoflinks.begin(); l<node[i].listoflinks.end() ; l++ )
			{	 
				l->linkednodeptr->colour.setcol((double)i/(double)highestnodecount);
			}
		}
		else
		{
		node[i].colour.r = node[i].colour.g = node[i].colour.b = (double)0.4;
		}
	}
*/
	return 0;
}

int actin::savebmp(int filenum, projection proj, processfgbg fgbg)
{ 
	// choose projection letter for filename etc.

	char projletter[] = "z";
	ofstream *p_outbmpfile;
	char *temp_BMP_filename;

    if (proj == xaxis)
	{
		sprintf ( projletter , "x");
		p_outbmpfile = &outbmpfile_x;
		temp_BMP_filename = temp_BMP_filename_x;
	}
	else if (proj == yaxis)
	{
		sprintf ( projletter , "y");
		p_outbmpfile = &outbmpfile_y;
		temp_BMP_filename = temp_BMP_filename_y;
	}
	else
	{
		p_outbmpfile = &outbmpfile_z;
		temp_BMP_filename = temp_BMP_filename_z;
	}

	cout << "Processing "<< projletter << "-axis BMP";
	cout.flush();


	double minx, miny, minz;
	double maxx, maxy, maxz; 
	//double VIEW_HEIGHT, VIEW_HEIGHT, VIEW_HEIGHT;

	int beadmaxx, beadmaxy;
	int beadminx, beadminy;
	//double centerx, centery, centerz;
	
	int x,y;

	double gaussmax = (double) GAUSSFWHM * 3 / 2;  // full extent of gaussian radius -  fwhm is 2/3 this

	double imageRmax, imageGmax, imageBmax;

	for (x = 0; x<BMP_WIDTH; x++)
		{
			for (y = 0; y<BMP_HEIGHT; y++)
			{
				imageR[x][y]=0;
				imageG[x][y]=0;
				imageB[x][y]=0;
			}
		}

	minx = miny = minz = 0;
	maxx = maxy = maxz = 0;

	// determine size

	for (int i=0; i<highestnodecount; i++)
	{
		if ((node[i].polymer) && (!node[i].listoflinks.empty()))
		{
			if (minx > node[i].x) minx = node[i].x;
			if (miny > node[i].y) miny = node[i].y;
			if (minz > node[i].z) minz = node[i].z;

			if (maxx < node[i].x) maxx = node[i].x;
			if (maxy < node[i].y) maxy = node[i].y;
			if (maxz < node[i].z) maxz = node[i].z;

		}
	} 

	double meanx, meany, meanz;

	// temp: reset center:
	meanx = -p_nuc->position.x; 
	meany = -p_nuc->position.y; 
	meanz = -p_nuc->position.z; 

	p_nuc->nucleator_rotation.rotate(meanx,meany,meanz);
	camera_rotation.rotate(meanx,meany,meanz); 
	//camera_rotation2.rotate(meanx,meany,meanz); 

	
	// precalculate gaussian

	Dbl2d GaussMat;

	const int xgmax = pixels(gaussmax);  // was BMP_WIDTH
	const int ygmax = pixels(gaussmax);

	int xg,yg;

	GaussMat.resize(2*xgmax);

	for(xg = -xgmax; xg<xgmax; xg++)
	{
		GaussMat[xg+xgmax].resize(2*ygmax);
		for(yg = -ygmax; yg<ygmax; yg++)
		{
			if ((xg*xg+yg*yg)>(xgmax*ygmax))
				continue;  // don't do corners

			GaussMat[xg+xgmax][yg+ygmax] 
					= exp(-3*((double)(xg*xg+yg*yg))/(double)(xgmax*ygmax));
		}
	}

	double keep_within_border;

	if (p_nuc->geometry == nucleator::sphere)
		 keep_within_border = 2* RADIUS;
	else
		 keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);

	
	if (proj == xaxis)  // choose projection
		{
			beadmaxx = pixels(   keep_within_border - meany) +  BMP_WIDTH/2;
			beadmaxy = pixels(   keep_within_border - meanz) +  BMP_HEIGHT/2;
			beadminx = pixels(-  keep_within_border - meany) +  BMP_WIDTH/2;
			beadminy = pixels(-  keep_within_border - meanz) +  BMP_HEIGHT/2;
		}
		else if (proj == yaxis)
		{
			beadmaxx = pixels(   keep_within_border - meanx) +  BMP_WIDTH/2;
			beadmaxy = pixels(   keep_within_border - meanz) +  BMP_HEIGHT/2;
			beadminx = pixels(-  keep_within_border - meanx) +  BMP_WIDTH/2;
			beadminy = pixels(-  keep_within_border - meanz) +  BMP_HEIGHT/2;
		}
		else 
		{
			beadmaxx = pixels(   keep_within_border - meanx) +  BMP_WIDTH/2;
			beadmaxy = pixels(   keep_within_border - meany) +  BMP_HEIGHT/2;
			beadminx = pixels(-  keep_within_border - meanz) +  BMP_WIDTH/2;
			beadminy = pixels(-  keep_within_border - meany) +  BMP_HEIGHT/2;
		}

	int movex=0;
	int movey=0;

	if (beadmaxx > BMP_WIDTH)
		movex =- (beadmaxx - BMP_WIDTH);
	if (beadmaxy > BMP_HEIGHT)
		movey =- (beadmaxy - BMP_HEIGHT);
	if (beadminx < 0)
		movex =- beadminx;
	if (beadminy < 0)
		movey =- beadminy;


	vect rot;

	for (int i=0; i<highestnodecount; i++)
	{
		if (!node[i].polymer)
			continue;

		rot = node[i];

		p_nuc->nucleator_rotation.rotate(rot);
		camera_rotation.rotate(rot);
		//camera_rotation2.rotate(rot);

		if (proj == xaxis)  // choose projection
		{
			x = pixels(rot.y - meany) +  BMP_WIDTH/2; 
			y = pixels(rot.z - meanz) + BMP_HEIGHT/2;
		}
		else if (proj == yaxis)
		{
			x = pixels(rot.x - meanx) +  BMP_WIDTH/2; 
			y = pixels(rot.z - meanz) + BMP_HEIGHT/2;
		}
		else 
		{
			x = pixels(rot.x - meanx) +  BMP_WIDTH/2; 
			y = pixels(rot.y - meany) + BMP_HEIGHT/2;
		}

		x += movex;  // displace to bring bead back in bounds
		y += movey;

		if ((x<0) || (x >= BMP_WIDTH  - (2*xgmax+1)) ||
			(y<0) || (y >= BMP_HEIGHT - (2*ygmax+1)) )  // only plot if point in bounds
		{
			//cout << "point out of bounds " << x << "," << y << endl;
		}
		else
		{
			if ((i%SPECKLE_FACTOR)==0)
				for(xg = -xgmax; xg<xgmax; ++xg)
					for(yg = -ygmax; yg<ygmax; ++yg)
					{
						if ((xg*xg+yg*yg)>(xgmax*ygmax))
							continue;  // don't do corners
						
						//imageR[x+xg+xgmax][y+yg+ygmax]+=		// link forces
						//	node[i].linkforce_transverse[0] * GaussMat[xg+xgmax][yg+ygmax];
						
						imageG[x+xg+xgmax][y+yg+ygmax]+=
								1 * GaussMat[xg+xgmax][yg+ygmax];  // amount of actin

					}
		}
	}

	// normalize image

	imageRmax = 1000/INIT_R_GAIN;
	imageGmax = 1000/INIT_G_GAIN;
	imageBmax = 1000/INIT_B_GAIN;  // prevent over sensitivity

	for (x = 0; x<BMP_WIDTH; x++)
	{
		for (y = 0; y<BMP_HEIGHT; y++)
		{
			if (imageR[x][y]>imageRmax)
					imageRmax=imageR[x][y];
			if (imageG[x][y]>imageGmax)
					imageGmax=imageG[x][y];
			if (imageB[x][y]>imageBmax)
					imageBmax=imageB[x][y];
		}
	}

	// re-scale from 0 to 1

	for (x = 0; x<BMP_WIDTH; x++)
	{
		for (y = 0; y<BMP_HEIGHT; y++)
		{
			imageR[x][y] /= imageRmax;
			imageG[x][y] /= imageGmax;
			imageB[x][y] /= imageBmax;
		}
	}


// draw the nucleator points cage:

	for (vector <vect>::iterator point=p_nuc->cagepoints.begin(); 
	      point<p_nuc->cagepoints.end() ; ++point )
	{
		// rotate point	

		rot = *point;
		p_nuc->nucleator_rotation.rotate(rot);
		camera_rotation.rotate(rot); 
		//camera_rotation2.rotate(rot); 

		if (proj == xaxis)  // choose projection
		{
			x = pixels(rot.y - meany) +  BMP_WIDTH/2; // was BMP_WIDTH
			y = pixels(rot.z - meanz) + BMP_HEIGHT/2;
		}
		else if (proj == yaxis)
		{
			x = pixels(rot.x - meanx)+  BMP_WIDTH/2; // was BMP_WIDTH
			y = pixels(rot.z - meanz)+ BMP_HEIGHT/2;
		}
		else 
		{
			x = pixels(rot.x - meanx) +  BMP_WIDTH/2; // was BMP_WIDTH
			y = pixels(rot.y - meany) + BMP_HEIGHT/2;
		}

		x+=movex;  // displace to bring bead back in bounds
		y+=movey;

		if (((x+xgmax)<0) || ((x+xgmax)>=BMP_WIDTH) ||
			(((y+ygmax)<0) || ((y+ygmax)>=BMP_HEIGHT)))  // only plot if point in bounds
			continue;
		
		imageR[x+xgmax][y+ygmax] = p_nuc->colour.r;
		imageG[x+xgmax][y+ygmax] = p_nuc->colour.g;
		imageB[x+xgmax][y+ygmax] = p_nuc->colour.b;
		//imageG[x+xgmax][y+ygmax] = (imageG[x+xgmax][y+ygmax] + imageGmax / 2.0) / 2.0;
		//imageB[x+xgmax][y+ygmax] = (imageB[x+xgmax][y+ygmax] + imageBmax / 2.0) / 2.0;

	}

	cout << ".";
	cout.flush();

	// write the bins graphics

	p_nuc->segs.write_bins_bitmap(imageR, imageG, imageB,
					   p_nuc->segs.link_transverse, proj);

	//char filename[255];
	
	//sprintf ( filename , "temp.bmp");

	// sprintf ( filename , "%s_proj_%05i.bmp", projletter , filenum );


	// write the bitmap file
		
	writebitmapfile(*p_outbmpfile, imageR, imageG, imageB);

	cout << ".";
	cout.flush();

	// add the imagemagick overlays
	
	char command1[10240], command2[10240], command3[20480];

    stringstream tmp_drawcmd1,tmp_drawcmd2,tmp_drawcmd3, drawcmd;

	// todo figure out how to clear a stringstream so don't have
	// to use as disposable objects
	// .clear() or .str("") don't seem to work.

	int scalebarlength = pixels(1.0);

	drawcmd << "-fill none -stroke grey -draw \"";
	p_nuc->segs.drawoutline(drawcmd,proj);				// draw outline of nucleator
					
	// check have lines to draw before adding them
	// else ImageMagick intermittently crashes

	if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd1,proj,1) > 0)		
		drawcmd << "\" -stroke blue -draw \"" << tmp_drawcmd1.str();

	if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd2,proj,0.2) > 0)	// scaled impact forces
		drawcmd << "\" -stroke red -draw \"" << tmp_drawcmd2.str();	

	if (p_nuc->segs.drawsurfaceimpacts(tmp_drawcmd3,proj,0.05) > 0)	// scaled impact forces
		drawcmd << "\" -stroke yellow -draw \"" << tmp_drawcmd3.str();

	drawcmd << "\"";

// drawing image text takes forever on alpha and OSX
// so have option to bypass

if (NO_IMAGE_TEXT)
{	
		sprintf(command1,
		"convert -quality %i -fill white -draw \"rectangle 5 576 %i 573\" %s %s%s_proj_%05i.%s",
		     BMP_COMPRESSION, scalebarlength+5, temp_BMP_filename, BITMAPDIR, 
			 projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());
}
else 
{
		sprintf(command1,
		"convert -quality %i -font helvetica -fill white -pointsize 20 -draw \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 '%s-projection  \\nFrame % 4i\\nG-gain % 4i'\" %s %s%s_proj_%05i.%s",
			BMP_COMPRESSION, scalebarlength+5, projletter, filenum,
			(int)((1000/(double)imageGmax)+0.5),  temp_BMP_filename, BITMAPDIR,
			projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());
}

	sprintf(command2,
		"convert -quality %i %s %s%s_proj_%05i.%s %s%s_forces_%05i.%s", 
		BMP_COMPRESSION, drawcmd.str().c_str(), BITMAPDIR, 
		projletter, filenum, BMP_OUTPUT_FILETYPE.c_str(), BITMAPDIR, 
		projletter, filenum, BMP_OUTPUT_FILETYPE.c_str());

if (fgbg == runbg)
	sprintf(command3, "(%s ; %s )&", command1, command2);
else
	sprintf(command3, "(%s ; %s )", command1, command2);

//cout << command3 << endl;;

	system(command3);

    //system(command1);
    //system(command2);

	// opruninfo << endl << command1 << endl << command2 << endl;

//#ifdef _WIN32  // use 'del' on windows, 'rm' on unix:
//	sprintf(command2 , "del %s 2>/dev/null", filename);
//#else
//	sprintf(command2 , "rm -f %s 2>/dev/null", filename);
//#endif
//
//	system(command2);	

	cout << ". ";
	cout.flush();

	return filenum;
}

void actin::writebitmapheader(ofstream& outbmpfile, const int & bitmapwidth, const int & bitmapheight)
{

	// bitmap headers etc (see microsoft website):

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

	for (int i=0;i<=255;i++)
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


int actin::squash(double thickness)
{
	// squash with 'coverslip'

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].x > thickness/2)
			node[i].x = thickness/2;

		if (node[i].x < -thickness/2)
			node[i].x = -thickness/2;
	}

	return 0;
}


int actin::sortnodesbygridpoint(void)
{
	nodes* nodeptr, *startnodeptr;

	for (int i=0; i<highestnodecount; i++)
	{
		donenode[i] = false;
	}

	//nodesbygridpoint.swap(nodesbygridpoint_temp);  // store the old order
	//nodesbygridpoint.resize(0);	// blank the new

	// mark: todo: alter so that not restarting from scratch (use lowestnodeupdate)

	// int i;
	int nodenumber = 0;

	for (int i=0; i<highestnodecount; ++i)
	{	// collect the nodes in gridpoint order...

		//i = nodesbygridpoint_temp[j]; // get reference for old set
		if (donenode[i]) //|| (!node[i].polymer))
			continue;

		nodeptr=nodegrid[node[i].gridx][node[i].gridy][node[i].gridz];
		startnodeptr=nodeptr;
			if (nodeptr!=0) 
			{
				do
				{
					nodesbygridpoint[nodenumber++]=nodeptr->nodenum;
					donenode[nodeptr->nodenum] = true;  // mark node as done
					nodeptr = nodeptr->nextnode;					
				}
				while (nodeptr!=startnodeptr);  //until back to start
			}
		}


	return 0;
}

void * actin::compressfilesthread(void* threadarg)
{
	struct thread_data *dat;
	dat = (struct thread_data *) threadarg;

	while (GRASS_IS_GREEN)
	{	
		//sem_wait(&compressfiles_thread_go[0]);  // go only if released by main thread

		pthread_mutex_lock(&filessavelock_mutex);
		pthread_mutex_lock(&filesdonelock_mutex);

		char command1[255];//,command2[255];

		// gzip the text data:

		//sprintf(command1, "gzip -f -9 nodes%05i.txt",dat->startnode);
		//system(command1);

		//sprintf(command1, "gzip -f -9 links%05i.txt",dat->startnode);
		//system(command1);

		sprintf(command1, "gzip -q -f -9 %sreport*.txt 2>/dev/null" ,TEMPDIR);
		system(command1);

		sprintf(command1 , "gzip -q -f -9 %sdata*.txt 2>/dev/null",TEMPDIR);
		system(command1);

		// gzip the wrl file:

		sprintf(command1 , "gzip -f -9 %snodes%05i.wrl",TEMPDIR, dat->startnode );
		system(command1);

		sprintf(command1 , "mv %snodes%05i.wrl.gz %snodes%05i.wrz", TEMPDIR,dat->startnode, VRMLDIR, dat->startnode );
		system(command1);

		sprintf(command1 , "mv %s*report*.gz %s",TEMPDIR, REPORTDIR);
		system(command1);

		sprintf(command1 , "mv %s*data*.gz %s", TEMPDIR, DATADIR);
		system(command1);


//		// convert bmps to jpgs:
//#ifndef NO_IMAGEMAGICK
//			
//		sprintf(command1 , "convert -quality 75 x_proj_%05i.bmp x_proj_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//		sprintf(command1 , "convert -quality 75 y_proj_%05i.bmp y_proj_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//		sprintf(command1 , "convert -quality 75 z_proj_%05i.bmp z_proj_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//
//		sprintf(command1 , "convert -quality 75 x_forces_%05i.bmp x_forces_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//		sprintf(command1 , "convert -quality 75 y_forces_%05i.bmp y_forces_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//		sprintf(command1 , "convert -quality 75 z_forces_%05i.bmp z_forces_%05i.png", dat->startnode, dat->startnode );
//		system(command1);
//
//#ifdef _WIN32  // use 'del' on windows, 'rm' on unix:
//		sprintf(command2 , "del x_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		sprintf(command2 , "del y_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		sprintf(command2 , "del z_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		//sprintf(command2 , "del x_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		//system(command2);
//		//sprintf(command2 , "del y_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		//system(command2);
//		//sprintf(command2 , "del z_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//#else
//		sprintf(command2 , "rm -f x_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		sprintf(command2 , "rm -f y_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		sprintf(command2 , "rm -f z_proj_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//		//sprintf(command2 , "rm -f x_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		//system(command2);
//		//sprintf(command2 , "rm -f y_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		//system(command2);
//		//sprintf(command2 , "rm -f z_forces_%05i.bmp 2>/dev/null", dat->startnode);
//		system(command2);
//#endif
	
//#endif

		pthread_mutex_unlock(&filesdonelock_mutex);

		//sem_post(&compressfiles_data_done[0]);  // release main thread block waiting for data
	}	  
	return (void*) NULL;
}

void actin::compressfilesdowork(const int & filenum)
{
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

	sprintf(command1 , "( gzip -f -9 %snodes%05i.wrl ; mv %snodes%05i.wrl.gz %snodes%05i.wrz) &",
						TEMPDIR, filenum, TEMPDIR,filenum, VRMLDIR,filenum );
	system(command1);

}



int actin::find_center(vect &center)
{
	center = - p_nuc->position;

	return 0;
}

//
//int actin::savedata(int filenum)
//{
//	char filename[255];
//
//	sprintf ( filename , "data%05i.txt", filenum );
//
//	ofstream opdata(filename, ios::out | ios::trunc);
//	if (!opdata) 
//	{ cout << "Unable to open file " << filename << " for output"; return 1;}
//
//	opdata << highestnodecount << "," << iteration_num 
//		<< "," << linksbroken << "," << linksformed << endl;
//
//	for (int i=0; i<highestnodecount; i++)
//	{
//		node[i].savedata(&opdata);
//	}
//
//	opdata << endl;
//	opdata << (unsigned int) crosslinknodesdelay.size();
//
//	for (vector <int>::iterator i=crosslinknodesdelay.begin(); i<crosslinknodesdelay.end() ; i++ )
//	{	 
//		opdata << "," << (*i);
//	}
//
//	opdata.close();
//
//	return 0;
//}
////
////int actin::loaddata(int filenum)
////{
////
////	for (int i=0; i!=(GRIDSIZE+1); i++) // clear the nodegrid
////		for (int j=0; j!=(GRIDSIZE+1); j++)
////			for (int k=0; k!=(GRIDSIZE+1); k++)
////				nodegrid[i][j][k]=0;;
////
////	char filename[255];
////	char delim;
////	unsigned int numcrosslinkdelay;
////
////	sprintf ( filename , "data%05i.txt", filenum );
////
////	ifstream ipdata(filename);
////	if (!ipdata) 
////	{ cout << "Unable to open file " << filename << " for output"; return 1;}
////
////	ipdata >> highestnodecount >> delim >> iteration_num 
////		>> delim >> linksbroken >> delim >> linksformed;
////
////	for (int i=0; i<highestnodecount; i++)
////	{
////		node[i].loaddata(&ipdata);
////	}
////
////	ipdata >> numcrosslinkdelay;
////	crosslinknodesdelay.resize(numcrosslinkdelay);
////
////	for (vector <int>::iterator i=crosslinknodesdelay.begin(); i<crosslinknodesdelay.end() ; i++ )
////	{	 
////		ipdata >> delim >> (*i);
////	}
////
////	return 0;
////}

void actin::clear_nodegrid()
{
    // clear the nodegrid
    for(int i=0; i<GRIDSIZE; i++)
		for(int j=0; j<GRIDSIZE; j++)
			for(int k = 0; k <GRIDSIZE; k++)
				nodegrid[i][j][k] = NULL;
}

int actin::save_data(ofstream &ofstrm)
{
    // write out all stateful information for this object    
    // save actin
    ofstrm << "actin: "  << endl 
	   << highestnodecount << ","
	   << nexttocrosslink << "," 
	   << iteration_num << "," 
	   << linksbroken  << "," 
	   << linksformed << endl
    	   << actin_rotation << ","
	   << camera_rotation << endl;
    
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
	    ofstrm << ",";
	else
	    ofstrm << "." << endl;
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
    char ch;
    
    // load actin
    ifstr >> str;
    // ensure the identifier for the start of the actin
    if(str.compare("actin:") !=0 ){
	cout << "error in checkpoint file, 'actin:' expected" << endl;
	return 1;
    }
    
    ifstr >> highestnodecount >> ch
	  >> nexttocrosslink  >> ch
	  >> iteration_num    >> ch 
	  >> linksbroken      >> ch 
	  >> linksformed      
	  >> actin_rotation   >> ch
	  >> camera_rotation;

    ifstr >> str;
    // load nodes
    if(str.compare("nodes-links:") !=0 ){
	cout << "error in data file, 'nodes-links:' expected" << endl;
	cout << "'" << str <<"' last read." << endl;
	return 1;
    }
    // ** Remember the node vector is preallocated to MAXNODES
    for(int i=0; i < highestnodecount; i++) 
	{
		node[i].load_data(ifstr);
    }

    // This is not neat, and we have a nasty reliance on the data being
    // public all the way down to the nodes linklist link objects.
    // When loading the link has been stored as an index, and we now
    // convert that to a pointer.
    // Method is 
    //  Run through and rebuild the node ptrs for the linkednode indices
    //  important this is done after the full node list has been built
    //    - iterate over nodes
    //    - then iterate over links in the node linklist
    for(int i=0; i < highestnodecount; i++) {
	for(vector<links>::iterator l=node[i].listoflinks.begin(); 
	    l<node[i].listoflinks.end(); ++l) {
	    // if(l->linkednodenumber>=0) // ensure the indx was explicitly set
		l->linkednodeptr = &node[l->linkednodenumber];
	}
    }

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
	i < crosslinknodesdelay.end(); i++) {
	ifstr >> (*i) >> ch;
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

void actin::setdontupdates(void)
{
	if (highestnodecount > NODES_TO_UPDATE)
	{
	for (int i=0; i<(highestnodecount-NODES_TO_UPDATE); i++)
		{
			node[i].dontupdate = true;
		}
	}
}

void actin::set_sym_break_axes(void)
{

	vect tmp_nodepos;
	vect CofM;
	vect sym_break_direction;

	rotationmatrix tmp_rotation, tmp_rotation2;	
	rotationmatrix final_rotation;

	double x_angle, y_angle; //, z_angle;
	//double x_inertia = 0, y_inertia = 0;
	//double z_in_angle;

	//int numnodes;

	find_center(sym_break_direction);  // which way did the bead go?

	//cout << endl << "old sym_break_direction: " 
	//<< sym_break_direction.x << "," 
	//<< sym_break_direction.y << "," 
	//<< sym_break_direction.z << endl;

	x_angle = atan2(sym_break_direction.y,sym_break_direction.z);

	tmp_rotation.rotatematrix(x_angle, 0 , 0);
	tmp_rotation.rotate(sym_break_direction);

	tmp_rotation.settoidentity();

	y_angle = -atan2(sym_break_direction.x,sym_break_direction.z);

	tmp_rotation.rotatematrix(0, y_angle , 0); // now pointing upwards
	tmp_rotation.rotate(sym_break_direction);

		//cout << "new sym_break_direction: " 
		//<< sym_break_direction.x << "," 
		//<< sym_break_direction.y << "," 
		//<< sym_break_direction.z << endl;


	tmp_rotation.settoidentity();

	tmp_rotation.rotatematrix(y_angle , rotationmatrix::yaxis);
	tmp_rotation.rotatematrix(x_angle , rotationmatrix::xaxis);



	//tmp_rotation.rotatematrix(x_angle,y_angle,0);

	//find_center(sym_break_direction);  // which way did the bead go?
	//tmp_rotation.rotate(sym_break_direction);

	//cout << "combined xy new sym_break_direction: " 
	//	<< sym_break_direction.x << "," 
	//	<< sym_break_direction.y << "," 
	//	<< sym_break_direction.z << endl;

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

		for (int i=0; i<highestnodecount; ++i)
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

	// now have all the angles
	// need to assemble them in the right (reverse) order z,y,x
	// to get the rotation matrix

	camera_rotation.settoidentity();

	camera_rotation.rotatematrix(maxchiangle, rotationmatrix::zaxis);
	camera_rotation.rotatematrix(y_angle, rotationmatrix::yaxis);
	camera_rotation.rotatematrix(x_angle, rotationmatrix::xaxis);

	reverse_camera_rotation.settoidentity();
	//reverse_camera_rotation2.settoidentity();

	reverse_camera_rotation.rotatematrix(-maxchiangle, rotationmatrix::zaxis);
	reverse_camera_rotation.rotatematrix(-y_angle, rotationmatrix::yaxis);
	reverse_camera_rotation.rotatematrix(-x_angle, rotationmatrix::xaxis);

	cout << "Symmetry broken.  Camera rotation angles: " << x_angle*180/PI << "," << y_angle*180/PI << "," << maxchiangle*180/PI << endl;

	symbreakiter = iteration_num;

}
//
//void actin::set_sym_break_axes(void)
//{
//
//	vect tmp_nodepos;
//
//	//vect  CofM;
//
//	vect sym_break_direction;
//
//	rotationmatrix tmp_rotation, tmp_rotation2;	
//
//	double x_angle, y_angle;//, z_angle;
//
//	//double x_inertia = 0, y_inertia = 0;
//	
//	//double z_in_angle;
//
//	//int numnodes;
//
//	find_center(sym_break_direction);  // which way did the bead go?
//
//	x_angle = atan2(sym_break_direction.y,sym_break_direction.z);
//	y_angle = atan2(sym_break_direction.x,sym_break_direction.z);
//	//z_angle = atan2(sym_break_direction.x,sym_break_direction.y);
//
//	// make a temp rotation matrix in that direction:
//
//	//tmp_rotation.rotatematrix(0, y_angle, 0); // just rotate y
//	//tmp_rotation.rotatematrix(x_angle, y_angle, 0); // now pointing upwards
//
//	//cout << "old sym_break_direction: " 
//	//	<< sym_break_direction.x << "," 
//	//	<< sym_break_direction.y << "," 
//	//	<< sym_break_direction.z << endl;
//
//	//tmp_rotation.rotate(sym_break_direction);
//
//	//cout << "new sym_break_direction: " 
//	//	<< sym_break_direction.x << "," 
//	//	<< sym_break_direction.y << "," 
//	//	<< sym_break_direction.z << endl;
//
//
//	// find CofM
//
//	//numnodes=0;
//
//	//for (int i=0; i<highestnodecount; ++i)
//	//{
//	//	if (!node[i].polymer)
//	//		continue;
//
//	//	CofM += node[i];
//
//	//	++numnodes;
//
//	//}
//
//	//CofM *= (1/ (double) numnodes);
//
//	//tmp_rotation.rotate(CofM);
//
//	double theta, chi, maxchi, maxchiangle;
//
//	maxchi = 0;
//
//	for(theta = -PI; theta < PI; theta+=PI/180)
//	{
//		chi = 0;
//
//		tmp_rotation.settoidentity();
//		tmp_rotation.rotatematrix(x_angle,y_angle,0);
//		tmp_rotation2.rotatematrix(0,0,theta);
//
//		for (int i=0; i<highestnodecount; ++i)
//		{
//			if ((!node[i].polymer) || (node[i].harbinger))
//				continue;
//
//			tmp_nodepos = node[i];
//			tmp_rotation.rotate(tmp_nodepos);  // rotate points to bring in line with sym break dir'n
//			tmp_rotation2.rotate(tmp_nodepos);	// rotate z
//
//			// only look at front of bead (i.e. new z<0, in direction of movement)
//			if (tmp_nodepos.z<0)
//			{
//				chi += tmp_nodepos.x * tmp_nodepos.x; // sum the squares pf distance from the axis
//			}
//			
//		}
//
//		if (chi > maxchi)
//		{
//			maxchi = chi;
//			maxchiangle = theta;
//		}
//
//	}
//
//	
//
//	camera_rotation.settoidentity();
//	camera_rotation2.settoidentity();
//
//	camera_rotation.rotatematrix(x_angle, y_angle, 0);
//	camera_rotation2.rotatematrix(0, 0, maxchiangle);
//
//	reverse_camera_rotation.settoidentity();
//	reverse_camera_rotation2.settoidentity();
//
//	reverse_camera_rotation.rotatematrix(-x_angle, -y_angle, 0);
//	reverse_camera_rotation2.rotatematrix(0, 0, -maxchiangle);
//
//
//	cout << "Symmetry broken.  Camera rotation angles: " << x_angle*180/PI << "," << y_angle*180/PI << "," << maxchiangle*180/PI << endl;
//
//	symbreakiter = iteration_num;
//
//}

void actin::save_sym_break_axes(void)
{
	ofstream opsymbreak("sym_break_axis.txt", ios::out | ios::trunc);
	if (!opsymbreak) 
	{ cout << "Unable to open file 'sym_break_axis.txt' for output"; return;}

	opsymbreak  << symbreakiter << " "
				<< camera_rotation << " "
				<< reverse_camera_rotation << endl;

	opsymbreak.close();

	//cout << "'sym_break_axis.txt' file written" << endl;
}

void actin::load_sym_break_axes(void)
{
	ifstream ipsymbreak("sym_break_axis.txt", ios::in);

	if (!ipsymbreak) 
	{ 
		cout << "Unable to open file 'sym_break_axis.txt' for input, skipping." << endl;
	}
	else
	{
		ipsymbreak  >> symbreakiter 
					>> camera_rotation 
					>> reverse_camera_rotation;

		ipsymbreak.close();
	}
}


void actin::clear_node_stats(void)
{
	for (int threadnum = 0; threadnum < NUM_THREADS; ++threadnum)
	{
		for (int i=0; i<highestnodecount; ++i)
		{
			if (!node[i].dontupdate)
				node[i].clearstats(threadnum);
		}
	}

}



