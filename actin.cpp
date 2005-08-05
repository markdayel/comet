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
	{ cout << "Unable to open file " << "velocities.txt" << " for output"; return;}

	opvelocityinfo << "time,x,y,z,vel" << endl;

	// node = new nodes[MAXNODES];

	node.resize(MAXNODES);
	donenode.resize(MAXNODES);
	repdonenode.resize(MAXNODES);
	nodesbygridpoint.resize(MAXNODES);
	nodesbygridpoint_temp.resize(MAXNODES);
	crosslinknodesdelay.resize(CROSSLINKDELAY);

	doreportiteration = 9999999;

	reportdat.resize(REPORT_NUM_VARIABLES);
	for (int i=0; i!=REPORT_NUM_VARIABLES; i++)  
	{
		reportdat[i].resize(REPORT_AVERAGE_ITTERATIONS+1);
	}


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

	newnodescolour.setcol(0);

	brokensymmetry = false;

	//debug:

	debug_num_rotate = 0;
	debug_num_displace = 0;


} 

actin::~actin(void)
{
	opruninfo.close();
	opvelocityinfo.close();

}

int actin::nucleate()
{
	int numnewnodes;

	// add new nodes
	numnewnodes = nucleation_object->addnodes();
	
return numnewnodes;
}

int actin::crosslinknewnodes(int numnewnodes)
{	
	int threadnum = 0;
	
	vect nodeposvec;
	MYDOUBLE distsqr;

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

int actin::save(int filenum)
{
	char filename[255];
	//char time[255], date[255];

    //_strtime( time );
    //_strdate( date );

	sprintf ( filename , "nodes%05i.txt", filenum );

	ofstream opactindata(filename, ios::out | ios::trunc);
	if (!opactindata) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	// write header
	opactindata << "comet (compiled " << __DATE__ << " " << __TIME__ << ") Mark J Dayel" << endl;
	//opactindata << "Data Output " << date << " " << time << endl;
	opactindata << "Nodes Data Itteration Number: " << filenum << endl;
	opactindata << "Highest Node Count: " << highestnodecount << endl;
		opactindata << "nodenum,x,y,z,harbinger,polymer,r,g,b,creationiter,numlinks " << endl;
	
	nucleation_object->save(&opactindata);



	// write data
	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].polymer)
		{
			node[i].save(&opactindata);
		}
	}

	opactindata.close();



	// now save links

	sprintf ( filename , "links%05i.txt", filenum );

	ofstream oplinkdata(filename, ios::out | ios::trunc);

	if (!oplinkdata) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	// write header
	oplinkdata << "comet (compiled " << __DATE__ << " " << __TIME__ << ") Mark J Dayel" << endl;
	//oplinkdata << "Data Output " << date << " " << time << endl;
	oplinkdata << "Links Data Itteration Number: " << filenum << endl;
	oplinkdata << "Highest Node Count: " << highestnodecount << endl;
	oplinkdata << "nodenum,numlinks,linknum1,linknum2,..." << endl;

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].polymer)
		{
			node[i].savelinks(&oplinkdata);
		}
	}
	oplinkdata.close();

	return 0;
}

int actin::saveinfo()
{
	//char filename[255];
	MYDOUBLE minx, miny, minz;
	MYDOUBLE maxx, maxy, maxz; 
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

	sprintf ( filename , "nodes%05i.wrl", filenum );

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

	nucleation_object->savevrml(&opvrml);

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

	ejectfromnucleator();

	numnewnodes = nucleate();
	
	// crosslink, but only after time to equilibrate position
	// in mean time it is in limbo
	// affected only by node collision repulsion and nucleator repulsion

	crosslinknodesdelay[iteration_num%CROSSLINKDELAY] = numnewnodes;
	crosslinknewnodes(crosslinknodesdelay[(iteration_num+1)%CROSSLINKDELAY]);
	
	sortnodesbygridpoint();
	collisiondetection();	// start collision detection threads
	linkforces(false);			// start link forces threads
	//repulsiveforces();		// start the repulsive forces thread

	/*
	//  wait 'till complete before applying forces
	if (USE_THREADS)
	{
		for (int i = 0; i < NUM_THREADS; i++)
		{
			sem_wait(&collision_data_done[i]);
			//sem_wait(&linkforces_data_done[i]);
			//sem_wait(&compressfiles_data_done[i]);
		}
	}
	*/


	//  debugging:
	if (((iteration_num+1)%InterRecordIterations)==0)
	{
		doreportiteration = 0;
		doreportmaxnodes = highestnodecount;

		// clear the array
		for (int i=0; i!=REPORT_NUM_VARIABLES; ++i)  
		{
			for (int j=0; j!=REPORT_AVERAGE_ITTERATIONS; ++j)
			{
				reportdat[i][j].resize(doreportmaxnodes+1,0);
			}
		}
	}

	if (doreportiteration < REPORT_AVERAGE_ITTERATIONS)
	{
		reportsnapshot(iteration_num,doreportmaxnodes,doreportiteration);
		++doreportiteration;

		if (doreportiteration==REPORT_AVERAGE_ITTERATIONS)
			savereport(iteration_num,doreportmaxnodes);
	}

	applyforces(); // move the nodes


	//squash(2.4f); // squash with 'coverslip' this dist is slide-coverslip dist
			
	iteration_num++;

	/*  debugging:
	if ((iteration_num%5000)==0)
	{
		for (int i = 0; i<highestnodecount;i++)
		{
			node[i].depolymerize();
		}
	}
  */  
	return 0;
}

int actin::addlinks(const int& linknode1,const int& linknode2)
{
	// crosslink a new node

	MYDOUBLE pxlink;

	if (linknode1==linknode2) return 0;

	if ((!node[linknode1].polymer)  || (!node[linknode2].polymer) ||
		(!node[linknode1].harbinger))
		return 0;

	//if ((linknode1==178) || (linknode1==180) ||
	//	(linknode2==178) || (linknode2==180) )
	//	cout << "stop here";

	MYDOUBLE dist = calcdist(node[linknode1].x-node[linknode2].x,
							 node[linknode1].y-node[linknode2].y,
							 node[linknode1].z-node[linknode2].z);
	
	// crosslink poisson distribution (max probability at XLINK_NODE_RANGE/5)

	//pxlink = ((MYDOUBLE)2.7*dist/(XLINK_NODE_RANGE/(MYDOUBLE)5))
	//						*exp(-dist/(XLINK_NODE_RANGE/(MYDOUBLE)5));

	// divide by dist*dist to compensate for increase numbers of nodes at larger distances



	//pxlink = ((( (MYDOUBLE) 2.7 / (XLINK_NODE_RANGE/(MYDOUBLE)5) )
	//						*exp(-dist/(XLINK_NODE_RANGE/(MYDOUBLE)5) ) ));

	// normal distrib with center around 1/2 of XLINK_NODE_RANGE
	// and scale magnitude by 1/(XLINK_NODE_RANGE^2) to compensate for increased
	// number of nodes at this distance shell

	//if ((node[linknode1].z <0) || (node[linknode2].z <0))
 	//	return 0;

	// was using this one before 28 march:
	// (poissonian with max at XLINK_NODE_RANGE/5):
	//pxlink = ((( (MYDOUBLE) 2.7 / (XLINK_NODE_RANGE/(MYDOUBLE)5) )
	//						*exp(-dist/(XLINK_NODE_RANGE/(MYDOUBLE)5) ) )/dist);


	// was using this (gaussian with max at XLINK_NODE_RANGE/2) for a bit:
	//pxlink = exp( -40*(dist-(XLINK_NODE_RANGE/2))*(dist-(XLINK_NODE_RANGE/2)))/(dist*dist*XLINK_NODE_RANGE*XLINK_NODE_RANGE);
	
	pxlink = 1;

	if ( (P_XLINK*pxlink) > ( ((MYDOUBLE) rand()) / (MYDOUBLE)RAND_MAX ) )
	{
		node[linknode1].addlink(&node[linknode2],dist); 
		node[linknode2].addlink(&node[linknode1],dist);
        
		//cout << "linked " << node[linknode1].nodenum << " to " << node[linknode2].nodenum  << endl;
	}

	return 0;
}

int actin::ejectfromnucleator()
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
					if ((nucleation_object->iswithinnucleator(nodeptr->x,nodeptr->y,nodeptr->z))
						&& (nodeptr->polymer))
					{	// inside nucleator
						nodes_within_nucleator.push_back(nodeptr);						
					}

					nodeptr = nodeptr->nextnode;
				}
				while (nodeptr!=startnodeptr);  //until back to start
				
			}
	}

	oldnucposn = nucleation_object->position;

	// do the ejection:

	for (vector <nodes*>::iterator i=nodes_within_nucleator.begin(); i<nodes_within_nucleator.end() ; i++ )
	{  //eject them
		if (nucleation_object->collision(*(*i))) // (*i)->x,(*i)->y,(*i)->z)==0)  
			(*i)->updategrid();  // ejected OK
		else
			(*i)->depolymerize();  // not ejected OK, depolymerize
	}


	// do rotation

bool toupdategrid = false;


	if (ROTATION && (nucleation_object->torque.length() > MIN_TORQUE_TO_UPDATE))
	{

		
		debug_num_rotate++;


		MYDOUBLE x_angle = nucleation_object->torque.x / 
				nucleation_object->momentofinertia.x;

		MYDOUBLE y_angle = nucleation_object->torque.y / 
				nucleation_object->momentofinertia.y;

		MYDOUBLE z_angle = nucleation_object->torque.z / 
				nucleation_object->momentofinertia.z;

		//MYDOUBLE x_angle = 0;  // for debugging
		//MYDOUBLE y_angle = 0;
		//MYDOUBLE z_angle = 0.00002;

		rotationmatrix torque_rotate; // ( x_angle, rotationmatrix::xaxis);

		torque_rotate.rotatematrix( x_angle, y_angle, z_angle);


		// rotate the actin:

		for (int i=0; i<highestnodecount; ++i)
		{
			torque_rotate.rotate(node[i]);
		}

        // rotate the nucleator:

		nucleation_object->nucleator_rotation.rotatematrix( -x_angle, -y_angle, -z_angle);


		// rotate the actin reference frame:

		actin_rotation.rotatematrix( x_angle, y_angle, z_angle);


		// clear the torque vectors:

		nucleation_object->torque.zero();

		toupdategrid = true;

	}


// do displacment

if ((nucleation_object->deltanucposn.length() > MIN_DISPLACEMENT_TO_UPDATE))
{

	debug_num_displace++;

	// move the nucleator by moving nodes in opposite direction
	// note this is *before* rotating the deltanucposn vector, since the
	// nodes are in the nucleator reference frame
	
	for (int i = 0; i<highestnodecount; i++)
	{
		node[i]-=nucleation_object->deltanucposn;
	}


	// rotate the nucleator displacement vector

	nucleation_object->nucleator_rotation.rotate(nucleation_object->deltanucposn);

	// update the nucleator position with the rotated vector

	nucleation_object->position+=nucleation_object->deltanucposn;

	// and zero

	nucleation_object->deltanucposn.zero();

	toupdategrid = true;
}

if (toupdategrid)
{
	for (int i = 0; i<highestnodecount; i++)
	{
		node[i].updategrid();
	}
}


	return 0;
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

inline void * actin::collisiondetectiondowork(thread_data* dat)
{
	vect nodeposvec;
	MYDOUBLE distsqr,dist;
	MYDOUBLE scale,force, forcescale;

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

				if (distsqr < SQRT_ACCURACY_LOSS)
				{	
					continue;
				}

				//dorepulsion((*sameGPnode)->nodenum,(*nearnode)->nodenum,dist,dat->threadnum);
				//dorepulsion(*(*sameGPnode),*(*nearnode),dist,dat->threadnum);
				
				if ( distsqr < NODE_REPULSIVE_RANGE*NODE_REPULSIVE_RANGE)
					{
						// calc dist between nodes

						dist = mysqrt(distsqr); 

					if (dist > NODE_INCOMPRESSIBLE_RADIUS)
					{
						force =  NODE_REPULSIVE_MAG - (NODE_REPULSIVE_MAG / 
							(NODE_REPULSIVE_RANGE - NODE_INCOMPRESSIBLE_RADIUS) * (dist-NODE_INCOMPRESSIBLE_RADIUS )) ;
					}
					else
					{
						force =  NODE_REPULSIVE_MAG ;
					}  // fix force if below node incompressible distance

					forcescale = force / dist;
					

#ifdef FORCES_BOTH_WAYS
					
					(*sameGPnode)->rep_force_vec[0] -= disp * forcescale;
					(*nearnode)->rep_force_vec[0] += disp * forcescale;

#else
					(*sameGPnode)->rep_force_vec[0] -= disp * (2.0 * forcescale) ;

					(*sameGPnode)->repforce_radial[0]     += fabs(  (*sameGPnode)->unitvec().dot(disp)  * (2 * forcescale) ) ;

					(*sameGPnode)->repforce_transverse[0] += fabs( ((*sameGPnode)->unitvec().cross(disp)).length()  * (2 * forcescale));


#endif

					if ( dist < NODE_INCOMPRESSIBLE_RADIUS)
					{
					                            
						// how far to repulse (quarter for each node) (half each, but we will do i to j and j to i)
						scale = (MYDOUBLE)NODE_INCOMPRESSIBLE_RADIUS / (MYDOUBLE)4.0*dist;  

						disp*=scale;	// these are the vector half components (in direction from i to j)

#ifdef FORCES_BOTH_WAYS

						(*sameGPnode)->repulsion_displacement_vec[0] -= disp;
						(*nearnode)->repulsion_displacement_vec[0]   += disp;

#else
						(*sameGPnode)->repulsion_displacement_vec[0] -= disp * 2;

						(*sameGPnode)->dispforce_radial[0]     += fabs(  (*sameGPnode)->unitvec().dot(disp)  * (2 * forcescale) ) ;

						(*sameGPnode)->dispforce_transverse[0] += fabs( ((*sameGPnode)->unitvec().cross(disp)).length()  * (2 * forcescale));

#endif
					}
				}

			}
		}
	}
return (void*) NULL;
}

//int actin::findnearbynodes(const int& ournodenum, const int& adjgridpoints, const int& threadnum)
inline int actin::findnearbynodes(const nodes& ournode, const int& adjgridpoints, const int& threadnum)
{

Nodes1d *p_recti_near_nodeslist = &recti_near_nodes[threadnum];
Nodes1d *p_nodes_on_same_gridpoint = &nodes_on_same_gridpoint[threadnum];

nodes *nodeptr, *startnodeptr;

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
//							  const MYDOUBLE& dist,const int& threadnum)
//{						// defunct, now in collisiondetectiondowork
//	if (&node_i==&node_j)
//		return 0;
//
//	MYDOUBLE scale;
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

if (false) // (USE_THREADS)  // switch this off for now...
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

	while (true)
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

int actin::linkforces(const bool& sumforces)
{
	MYDOUBLE scale,dist;
	MYDOUBLE force;
	vect nodeposvec, disp;
	
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

if (false) // USE_THREADS)
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
			if ((node[n].dontupdate && !sumforces) || (!node[n].polymer))
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
					}
					else
					{
						scale = 1/dist;

						if (sumforces)
						{ //used to calculate stored link energy
							node[n].link_force_vec[0].x += fabs(force);
							i->linkednodeptr->link_force_vec[0].x += fabs(force);
						}
						else
						{
#ifdef FORCES_BOTH_WAYS
						node[n].link_force_vec[0] += disp * (force * scale);
						i->linkednodeptr->link_force_vec[0] -= disp * (force * scale);
#else
						node[n].link_force_vec[0] += disp * (2 * force * scale);

						node[n].linkforce_radial[0]     += fabs(  node[n].unitvec().dot(disp)  * (2 * force * scale) ) ;

						node[n].linkforce_transverse[0] += fabs( (node[n].unitvec().cross(disp)).length()  * (2 * force * scale));

						//cout << node[n].linkforce_radial[0] << " " <<  node[n].linkforce_transverse[0] << endl;
#endif
						}

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

	MYDOUBLE xdist, ydist, zdist, dist;
	MYDOUBLE force;
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

	MYDOUBLE maxcol,mincol;
	maxcol = 0;
	mincol = 0;	

	linkforces(true);

	MYDOUBLE val;

	maxcol = (MYDOUBLE) 0.001;

	for (int i=0; i<highestnodecount; i++)
	{
		if (node[i].link_force_vec[0].x>0.0001)
			val = mysqrt(node[i].link_force_vec[0].x);
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
			val = mysqrt(node[i].link_force_vec[0].x);
		else
			val = 0;
		//node[i].colour.setcol((MYDOUBLE)node[i].creation_iter_num/(MYDOUBLE)TOTAL_ITERATIONS);
		//node[i].colour.setcol(((MYDOUBLE)node[i].nodelinksbroken-mincol)/(maxcol-mincol));
		if ((node[i].polymer)&& (!node[i].listoflinks.empty()) && val > 0.001)
			{
				node[i].colour.setcol((val-mincol)/(maxcol-mincol));
			}
		else
		{
			node[i].colour.setcol(0);
		}
		for (int threadnum=0; threadnum < NUM_THREADS; threadnum++)
				{  // clear the dummy values
					node[i].link_force_vec[threadnum].x = 0;
					node[i].link_force_vec[threadnum].y = 0;
					node[i].link_force_vec[threadnum].z = 0;
				}
	}


	//for (int i=0; i<highestnodecount; i++)
	//{
	//node[i].colour.setcol((MYDOUBLE)node[i].creation_iter_num/(MYDOUBLE)iteration_num);
	//}
/*
	MYDOUBLE maxcol,mincol;
	maxcol = 1;
	mincol = 0;	
	MYDOUBLE val;

	for (int i=0; i<highestnodecount; i++)
	{
		val = (MYDOUBLE)node[i].nodelinksbroken;
		if ((node[i].polymer) && (val > maxcol))
		{
			maxcol = val;
		}
	} 

	for (int i=0; i<highestnodecount; i++)
	{
	node[i].colour.setcol((MYDOUBLE)node[i].nodelinksbroken/maxcol);
	}

	for (int i=0; i<highestnodecount;i++)
	{
		//i = (int) ((MYDOUBLE) highestnodecount * ( (MYDOUBLE) rand()) / (MYDOUBLE)RAND_MAX);
		if (i%300==0)
		{
		node[i].colour.setcol((MYDOUBLE)i/(MYDOUBLE)highestnodecount);
		for (vector <links>::iterator l=node[i].listoflinks.begin(); l<node[i].listoflinks.end() ; l++ )
			{	 
				l->linkednodeptr->colour.setcol((MYDOUBLE)i/(MYDOUBLE)highestnodecount);
			}
		}
		else
		{
		node[i].colour.r = node[i].colour.g = node[i].colour.b = (MYDOUBLE)0.4;
		}
	}
*/
	return 0;
}

int actin::savebmp(int filenum, projection proj)
{ // todo: change to do 3 axes in one go to save calculation time

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

	char filename[255];

	if (proj == xaxis)  // choose projection
	{
		sprintf ( filename , "x_proj_%05i.bmp", filenum );
	}
	else if (proj == yaxis)
	{
		sprintf ( filename , "y_proj_%05i.bmp", filenum );
	}
	else 
	{
		sprintf ( filename , "z_proj_%05i.bmp", filenum );
	}

	ofstream outbmpfile(filename, ios::out | ios::binary | ios::trunc);
	if (!outbmpfile) 
	{ cout << "Unable to open file '" << filename << "' for output"; return 1;}

	MYDOUBLE minx, miny, minz;
	MYDOUBLE maxx, maxy, maxz; 
	MYDOUBLE dimx, dimy, dimz;

	int beadmaxx, beadmaxy;
	int beadminx, beadminy;
	//MYDOUBLE centerx, centery, centerz;
	
	int width = 800;
	int height = 600;

	int x,y;

	MYDOUBLE gaussmax = (MYDOUBLE) GAUSSFWHM * 3 / 2;  // full extent of gaussian radius -  fwhm is 2/3 this

	MYDOUBLE imageRmax, imageGmax, imageBmax;

	imageR.resize(width);
	imageG.resize(width);
	imageB.resize(width);

	for (x = 0; x<width; x++)
		{
			imageR[x].resize(height);
			imageG[x].resize(height);
			imageB[x].resize(height);

			for (y = 0; y<height; y++)
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

	//dimx = mymax(((maxx - minx) + 4 * gaussmax),12);
	//dimy = mymax(((maxy - miny) + 4 * gaussmax),12);
	//dimz = mymax(((maxz - minz) + 4 * gaussmax),12);

 	dimx = dimy = dimz = VIEW_HEIGHT;  //these are the x,y and z scales (should be equal)

	MYDOUBLE meanx, meany, meanz;

/*

	MYDOUBLE val;

	Dbl1d linkforces;

	linkforces.reserve(MAXNODES);

	meanx = meany = meanz = 0.0;
	
	for (int i=0; i<highestnodecount; i++)
	{
		val =0;
		for (int threadnum=0; threadnum < NUM_THREADS; threadnum++)
		{
			val+= sqrt(node[i].link_force_vec[threadnum].x*node[i].link_force_vec[threadnum].x +
					node[i].link_force_vec[threadnum].y*node[i].link_force_vec[threadnum].y +
					node[i].link_force_vec[threadnum].z*node[i].link_force_vec[threadnum].z);
		}
		if ((node[i].polymer) && (!node[i].listoflinks.empty()))
		{
			//maxcol = node[i].nodelinksbroken;
			linkforces[i] = val;
		}
		meanx+= node[i].x;
		meany+= node[i].y;
		meanz+= node[i].z;
	}



	meanx/=highestnodecount;
	meany/=highestnodecount;
	meanz/=highestnodecount;

*/

	// temp: reset center:
	meanx = -nucleation_object->position.x; 
	meany = -nucleation_object->position.y; 
	meanz = -nucleation_object->position.z; 

	nucleation_object->nucleator_rotation.rotate(meanx,meany,meanz);
	camera_rotation.rotate(meanx,meany,meanz); 
	
	//meanx = (maxx+minx)/2;
	//meany = (maxy+miny)/2;
	//meanz = (maxz+minz)/2;

	// precalculate gaussian

	Dbl2d GaussMat;

	const int xgmax = (int)((height *( ((gaussmax)/dimy) ))-0.5);  // was width
	const int ygmax = (int)((height *( ((gaussmax)/dimz) ))-0.5);

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
					= exp(-3*((MYDOUBLE)(xg*xg+yg*yg))/(MYDOUBLE)(xgmax*ygmax));
		}
	}

	MYDOUBLE keep_witin_border;

	if (nucleation_object->geometry == nucleator::sphere)
		keep_witin_border = 2* RADIUS;
	else
		keep_witin_border = SEGMENT/2+(2*RADIUS);

	
	if (proj == xaxis)  // choose projection
		{
			beadmaxx =(int)(((height *(((  keep_witin_border - meany)/dimy) ))-0.5) +  width/2);
			beadmaxy =(int)(((height *(((  keep_witin_border - meanz)/dimz) ))-0.5) +  height/2);
			beadminx =(int)(((height *(((- keep_witin_border - meany)/dimy) ))-0.5) +  width/2);
			beadminy =(int)(((height *(((- keep_witin_border - meanz)/dimz) ))-0.5) +  height/2);
		}
		else if (proj == yaxis)
		{
			beadmaxx =(int)(((height *(((  keep_witin_border - meanx)/dimx) ))-0.5) +  width/2);
			beadmaxy =(int)(((height *(((  keep_witin_border - meanz)/dimz) ))-0.5) +  height/2);
			beadminx =(int)(((height *(((- keep_witin_border - meanx)/dimx) ))-0.5) +  width/2);
			beadminy =(int)(((height *(((- keep_witin_border - meanz)/dimz) ))-0.5) +  height/2);
		}
		else 
		{
			beadmaxx =(int)(((height *(((  keep_witin_border - meanx)/dimx) ))-0.5) +  width/2);
			beadmaxy =(int)(((height *(((  keep_witin_border - meany)/dimy) ))-0.5) +  height/2);
			beadminx =(int)(((height *(((- keep_witin_border - meanz)/dimx) ))-0.5) +  width/2);
			beadminy =(int)(((height *(((- keep_witin_border - meany)/dimy) ))-0.5) +  height/2);
		}

	int movex=0;
	int movey=0;

	if (beadmaxx>width)
		movex=-(beadmaxx-width);
	if (beadmaxy>height)
		movey=-(beadmaxy-height);
	if (beadminx<0)
		movex=-(beadminx);
	if (beadminy<0)
		movey=-(beadminy);


//int harbingers = 0;

//MYDOUBLE xx,yy,zz;
//MYDOUBLE rotx, roty, rotz;

vect rot;

	for (int i=0; i<highestnodecount; i++)
	{
		if (!node[i].polymer)
			continue;

		rot = node[i];

		nucleation_object->nucleator_rotation.rotate(rot);
		camera_rotation.rotate(rot);

		if (proj == xaxis)  // choose projection
		{
			x = (int)(((height *(((rot.y - meany)/dimy) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.z - meanz)/dimz) ))-0.5) + height/2);
		}
		else if (proj == yaxis)
		{
			x = (int)(((height *(((rot.x - meanx)/dimx) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.z - meanz)/dimz) ))-0.5) + height/2);
		}
		else 
		{
			x = (int)(((height *(((rot.x - meanx)/dimx) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.y - meany)/dimy) ))-0.5) + height/2);
		}

		x+=movex;  // displace to bring bead back in bounds
		y+=movey;

		if ((x<0) || (x>=width- (2*xgmax+1)) ||
			(y<0) || (y>=height-(2*ygmax+1)))  // only plot if point in bounds
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
						//	node[i].nodelinksbroken * GaussMat[xg+xgmax][yg+ygmax];
						//imageR[x+xg+xgmax][y+yg+ygmax]+=		// link forces
						//		linkforces[i] * GaussMat[xg+xgmax][yg+ygmax];
						
						imageR[x+xg+xgmax][y+yg+ygmax]+=		// link forces
							node[i].linkforce_transverse[0] * GaussMat[xg+xgmax][yg+ygmax];
						
						imageG[x+xg+xgmax][y+yg+ygmax]+=
								1 * GaussMat[xg+xgmax][yg+ygmax];  // amount of actin
						//imageB[x+xg+xgmax][y+yg+ygmax]+=          // Blue: number of links 
						//		node[i].listoflinks.size() * GaussMat[xg+xgmax][yg+ygmax];
						//imageB[x+xg+xgmax][y+yg+ygmax]+=          // Blue: number of links 
						//		node[i].listoflinks.size() * GaussMat[xg+xgmax][yg+ygmax];
					}



		}
//	if (node[i].harbinger)
//		harbingers++;
	}
//cout << endl << "Harbingers:" << harbingers << endl;


		// normalize image

	

	imageRmax = 1000/INIT_R_GAIN;
	imageGmax = 1000/INIT_G_GAIN;
	imageBmax = 1000/INIT_B_GAIN;  // prevent over sensitivity

	// imageRmax = 1;  // max sensitivity for red channel

	for (y = 0; y<height; y++)
		{
			for (x = 0; x<width; x++)
			{
				if (imageR[x][y]>imageRmax)
						imageRmax=imageR[x][y];
				if (imageG[x][y]>imageGmax)
						imageGmax=imageG[x][y];
				if (imageB[x][y]>imageBmax)
						imageBmax=imageB[x][y];
			}
		}


// draw the nucleator points cage:

for (vector <vect>::iterator point=nucleation_object->cagepoints.begin(); 
	      point<nucleation_object->cagepoints.end() ; ++point )
	{
		// rotate point	

		rot = *point;
		nucleation_object->nucleator_rotation.rotate(rot);
		camera_rotation.rotate(rot); 

		if (proj == xaxis)  // choose projection
		{
			x = (int)(((height *(((rot.y - meany)/dimy) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.z - meanz)/dimz) ))-0.5) + height/2);
		}
		else if (proj == yaxis)
		{
			x = (int)(((height *(((rot.x - meanx)/dimx) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.z - meanz)/dimz) ))-0.5) + height/2);
		}
		else 
		{
			x = (int)(((height *(((rot.x - meanx)/dimx) ))-0.5) +  width/2); // was width
			y = (int)(((height *(((rot.y - meany)/dimy) ))-0.5) + height/2);
		}

		x+=movex;  // displace to bring bead back in bounds
		y+=movey;

		if (((x+xgmax)<0) || ((x+xgmax)>=width) ||
			(((y+ygmax)<0) || ((y+ygmax)>=height)))  // only plot if point in bounds
			continue;
		
		imageR[x+xgmax][y+ygmax] = nucleation_object->colour.r * imageRmax;
		imageG[x+xgmax][y+ygmax] = nucleation_object->colour.g * imageGmax;
		imageB[x+xgmax][y+ygmax] = nucleation_object->colour.b * imageBmax;
		//imageG[x+xgmax][y+ygmax] = (imageG[x+xgmax][y+ygmax] + imageGmax / 2.0) / 2.0;
		//imageB[x+xgmax][y+ygmax] = (imageB[x+xgmax][y+ygmax] + imageBmax / 2.0) / 2.0;

	}

	
	// the header for saving the bitmaps...

	BITMAPFILEHEADER *fileHeader;
	BITMAPINFO       *fileInfo;

	fileHeader = (BITMAPFILEHEADER*)calloc(1, sizeof( BITMAPFILEHEADER ));
	fileInfo = (BITMAPINFO*)calloc(1, sizeof( BITMAPINFO ) + 256 * sizeof(RGBQUAD));

	fileHeader->bfType = 0x4d42;
	fileHeader->bfSize = (DWORD)(sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFO) + (width * height* 3));
	fileHeader->bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFO) + (256*sizeof(RGBQUAD));

	fileInfo->bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
	fileInfo->bmiHeader.biWidth         = width;
	fileInfo->bmiHeader.biHeight        = height;
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

	// save headers

	outbmpfile.write((char*)fileHeader,sizeof(BITMAPFILEHEADER));
	outbmpfile.write((char*)fileInfo,sizeof(BITMAPINFO) + 256*sizeof(RGBQUAD));

	// re-scale for byte output and save image data

	RGB *line;
	line = new RGB[width];

	// write out the data, line by line (note: y is backwards)
	for (y = (height-1); y>=0; y--)
		{
		//outbmpfile.write(picbuff + (width*y), width);
		for (x = 0; x<width; x++)
			{
			line[x].B=(unsigned char)(255 * (MYDOUBLE)(imageB[x][y])/(MYDOUBLE)imageBmax);
			line[x].G=(unsigned char)(255 * (MYDOUBLE)(imageG[x][y])/(MYDOUBLE)imageGmax);
			line[x].R=(unsigned char)(255 * (MYDOUBLE)(imageR[x][y])/(MYDOUBLE)imageRmax);
			}
			outbmpfile.write((char*)line,width*3);
		}

	outbmpfile.close();

	delete [] line;
	free(fileHeader);
	free(fileInfo);

	if (FORCES_ON_SIDE)
	{
		movex = 50-(width/2);
		meanx=meany=meanz=0.0;
	}

	if( nucleation_object->is_sphere() ) 
	{
	    draw_bead_forces(filenum, proj,
			     width, height,
			     xgmax, ygmax,
			     movex, movey,
			     meanx, meany, meanz,
			     dimx, dimy, dimz, imageGmax,
				 false, 1, 1, "white");

		// bead overlay not implemented yet
		//if (NUCLEATOR_FORCES)
		//{
		//	draw_bead_forces(filenum, proj,
		//			width, height,
		//			xgmax, ygmax,
		//			movex, movey,
		//			meanx, meany, meanz,
		//			dimx, dimy, dimz, imageGmax,
		//			true, 0.2, 3, "red");
		//}

	} 
	else 
	{
	    draw_capsule_forces(filenum, proj,
			     width, height,
			     xgmax, ygmax,
			     movex, movey,
			     meanx, meany, meanz,
			     dimx, dimy, dimz, imageGmax,
				 false, 1, 1, "white");

		if (NUCLEATOR_FORCES)
		{
			draw_capsule_forces(filenum, proj,
					width, height,
					xgmax, ygmax,
					movex, movey,
					meanx, meany, meanz,
					dimx, dimy, dimz, imageGmax,
					true, 0.2, 3, "red");
		}
	}

char command2[255];

#ifdef _WIN32  // use 'del' on windows, 'rm' on unix:

	if (proj == xaxis)  // choose projection
	{
		sprintf(command2 , "del x_proj_%05i.bmp 2>/dev/null", filenum);
	}
	else if (proj == yaxis)
	{
		sprintf(command2 , "del y_proj_%05i.bmp 2>/dev/null", filenum);
	}
	else 
	{
		sprintf(command2 , "del z_proj_%05i.bmp 2>/dev/null", filenum);
	}

#else

	if (proj == xaxis)  // choose projection
	{
		sprintf(command2 , "rm -f x_proj_%05i.bmp 2>/dev/null", filenum);
	}
	else if (proj == yaxis)
	{
		sprintf(command2 , "rm -f y_proj_%05i.bmp 2>/dev/null", filenum);
	}
	else 
	{
		sprintf(command2 , "rm -f z_proj_%05i.bmp 2>/dev/null", filenum);
	}

system(command2);		


#endif

	return filenum;
}

void actin::draw_capsule_forces(int filenum, projection proj,
				int width, int height,
				int xgmax, int ygmax,
				int movex, int movey,
				MYDOUBLE meanx, MYDOUBLE meany, MYDOUBLE meanz,
				MYDOUBLE dimx, MYDOUBLE dimy, MYDOUBLE dimz,
				MYDOUBLE imageGmax, 
				bool nuc_forces_overlay, MYDOUBLE scale, MYDOUBLE thickness,  char*  color)
{
    // FIXME to stringstream
    char command1[2048], command2[2048];

    stringstream drawcmd;
    
    // bead center:
    int c_centerx, c_centery;
    int c_radius = 0;
    int c_segment = 0;
    
    // draw the outline of the capsule
    if (proj == xaxis) {  // choose projection
	
	c_centerx = xgmax+movex + (int)(((height *(((0 - meany)/dimy) ))-0.5) +  width/2); // was width
	c_centery = xgmax+movey + (int)(((height *(((0 - meanz)/dimz) ))-0.5) + height/2);
	
	c_radius  = (int)(height * ( (RADIUS)/dimz) ); 
	c_segment = (int)(height * ( (SEGMENT/2.0)/dimz) ); 
	
    } else if (proj == yaxis) {
	
	c_centerx = xgmax+movex + (int)(((height *(((0 - meanx)/dimx) ))-0.5) +  width/2); // was width
	c_centery = xgmax+movey + (int)(((height *(((0 - meanz)/dimz) ))-0.5) + height/2);
	
	c_radius  = (int)(height * ( (RADIUS)/dimx) ); 
	c_segment = (int)(height * ( (SEGMENT/2.0)/dimz) );
	
    } else {
	
	c_centerx = xgmax+movex + (int)(((height *(((0 - meanx)/dimx) ))-0.5) +  width/2); // was width
	c_centery = xgmax+movey + (int)(((height *(((0 - meany)/dimy) ))-0.5) + height/2);

	c_radius  = (int)(height * ( (RADIUS)/dimx) );
	
    }


    
    if( (proj == xaxis) || (proj == yaxis) ) 
	{
		drawcmd << "ellipse "
			<< c_centerx << "," << c_centery + c_segment << " "
			<< c_radius << "," << c_radius << " "
			<< "0,180";
		
		drawcmd << " line "
			<< c_centerx - c_radius << "," << c_centery + c_segment << " "
			<< c_centerx - c_radius << "," << c_centery - c_segment;
		
		drawcmd << " line "
			<< c_centerx + c_radius << "," << c_centery + c_segment << " "
			<< c_centerx + c_radius << "," << c_centery - c_segment;
		
		drawcmd << " ellipse "
			<< c_centerx << "," << c_centery - c_segment << " "
			<< c_radius << "," << c_radius << " "
			<< "180,360";
    } 
	else 
	{
		// z projection
		drawcmd << "circle "
			<< c_centerx << "," << c_centery
			<< " "
			<< c_centerx << "," << c_centery+c_radius;
    }
    
    
    // draw force indicators
    int lineoriginx,lineoriginy, lineendx, lineendy;
    
    /*
    MYDOUBLE radial_rep_tot_x=(MYDOUBLE) 0.01;
    MYDOUBLE radial_rep_tot_y=(MYDOUBLE) 0.01;
    MYDOUBLE radial_rep_tot_z=(MYDOUBLE) 0.01;
    
    MYDOUBLE radial_rep_mean_x, radial_rep_mean_y, radial_rep_mean_z;
    */
    
    // scale by the number of iterations between outputs
    for(int i=0; i<nucleation_object->n_force_segments(); i++)
    {
		nucleation_object->radial_rep_distrib_x[i]/=InterRecordIterations;
		nucleation_object->radial_rep_distrib_y[i]/=InterRecordIterations;
		nucleation_object->radial_rep_distrib_z[i]/=InterRecordIterations;
    }

    /*
    // increment the totals
    for(int i=0; i<RADIAL_SEGMENTS; i++)
    {
	radial_rep_tot_x+=nucleation_object->radial_rep_distrib_x[i];
	radial_rep_tot_y+=nucleation_object->radial_rep_distrib_y[i];
	radial_rep_tot_z+=nucleation_object->radial_rep_distrib_z[i];
    }

    radial_rep_mean_x = radial_rep_tot_x / (MYDOUBLE) nucleation_object->n_force_segments();
    radial_rep_mean_y = radial_rep_tot_y / (MYDOUBLE) nucleation_object->n_force_segments();
    radial_rep_mean_z = radial_rep_tot_z / (MYDOUBLE) RADIAL_SEGMENTS;
    */
    
    // Find and plot the force from the radial_rep_totals
    MYDOUBLE linescale;
    const MYDOUBLE ticklength = 0.01;

    // body
    for(int i=0; i<nucleation_object->nbdy_segs; i++)
    {
	/*
	cout << (nucleation_object->radial_rep_distrib_x[i]) << ", " 
	     << (nucleation_object->radial_rep_distrib_y[i]) << ", "
	     << (nucleation_object->radial_rep_distrib_z[i]) << endl;
	*/
	if (proj == xaxis)  // choose projection
	{
	    lineoriginx = xgmax+movex + (int)(((height *(((nucleation_object->fbar_bdy_x[i]-meany)/dimy) ))-0.5) + width/2);
	    lineoriginy = xgmax+movey + (int)(((height *(((nucleation_object->fbar_bdy_y[i]-meanz)/dimz) ))-0.5) + height/2);
	    
	    // length 0-1
	    linescale = ((nucleation_object->radial_rep_distrib_x[i])* scale * FORCEBAR_SCALE) + ticklength;
	    
		if (linescale > 0.9*RADIUS)
			linescale = RADIUS;

	    lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
	    lineendy = int(lineoriginy);	    
	
	    drawcmd << " line "
		    << lineoriginx << "," << lineoriginy << " "
		    << lineendx    << "," << lineendy;
	}
	else if (proj == yaxis)  // choose projection
	{
	    lineoriginx = xgmax+movex + (int)(((height*(((nucleation_object->fbar_bdy_x[i]-meanx)/dimx) ))-0.5) + width/2);
	    lineoriginy = xgmax+movey + (int)(((height*(((nucleation_object->fbar_bdy_y[i]-meanz)/dimz) ))-0.5) + height/2);
	    
	    // length 0-1
	    linescale = ((nucleation_object->radial_rep_distrib_y[i])* scale * FORCEBAR_SCALE) + ticklength;
	    
		if (linescale > 0.9*RADIUS)
			linescale = RADIUS;

	    lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
	    lineendy = int(lineoriginy);	    
	
	    drawcmd << " line "
		    << lineoriginx << "," << lineoriginy << " "
		    << lineendx    << "," << lineendy;
	}
	// do nothing additional for z
    }


    //  end caps
    for(int i=0; i< nucleation_object->ncap_segs; i++)
    {
	// choose projection
	if (proj == xaxis) {
	    
	    lineoriginx = xgmax+movex + (int)(((height*(((nucleation_object->fbar_cap_x[i]-meany)/dimy) ))-0.5) + width/2);
	    lineoriginy = xgmax+movey + (int)(((height*(((nucleation_object->fbar_cap_y[i]-meanz)/dimz) ))-0.5) + height/2);

	    // length 0-1
	    linescale = ((nucleation_object->radial_rep_distrib_x[nucleation_object->nbdy_segs + i])* scale * FORCEBAR_SCALE) + ticklength;
	    
		if (linescale > 0.9*RADIUS)
			linescale = RADIUS;
		
		// cout << linescale << endl;
	    if(nucleation_object->fbar_cap_y[i]>0) {
		// upper cap
		lineoriginy += c_segment;
		lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
		lineendy = int(lineoriginy - linescale*(lineoriginy - c_centery - c_segment));
	    } else {
		// lower cap
		lineoriginy -= c_segment;
		lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
		lineendy = int(lineoriginy - linescale*(lineoriginy - c_centery + c_segment));
	    }
	    drawcmd << " line "
		    << lineoriginx << "," << lineoriginy << " "
		    << lineendx    << "," << lineendy;
	    
	} else if (proj == yaxis) {
	    
	    lineoriginx = xgmax+movex + (int)(((height*(((nucleation_object->fbar_cap_x[i]-meanx)/dimx) ))-0.5) + width/2);
	    lineoriginy = xgmax+movey + (int)(((height*(((nucleation_object->fbar_cap_y[i]-meanz)/dimz) ))-0.5) + height/2);
	    
	    // length 0-1
	    linescale = ((nucleation_object->radial_rep_distrib_y[nucleation_object->nbdy_segs + i])* scale * FORCEBAR_SCALE) + ticklength;
	    if (linescale > 0.9*RADIUS)
			linescale = RADIUS;
		
		//cout << "y endcap: " << linescale << endl;
	    if(nucleation_object->fbar_cap_y[i]>0) {
		// upper cap
		lineoriginy += c_segment;
		lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
		lineendy = int(lineoriginy - linescale*(lineoriginy - c_centery - c_segment));
	    } else {
		// lower cap
		lineoriginy -= c_segment;
		lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
		lineendy = int(lineoriginy - linescale*(lineoriginy - c_centery + c_segment));
	    }
	    drawcmd << " line "
		    << lineoriginx << "," << lineoriginy << " "
		    << lineendx    << "," << lineendy;
	} else {
	    
	    // z projection
	    lineoriginx = xgmax+movex + (int)(((height *(((nucleation_object->fbar_cap_x[i]-meanx)/dimx) ))-0.5) + width/2);
	    lineoriginy = xgmax+movey + (int)(((height *(((nucleation_object->fbar_cap_y[i]-meany)/dimy) ))-0.5) + height/2);
	    
	    // length 0-1
	    linescale = ((nucleation_object->radial_rep_distrib_z[i])* scale * FORCEBAR_SCALE) + ticklength;
	    // cout << linescale << endl;
		if (linescale > 0.9*RADIUS)
			linescale = RADIUS;

	    lineendx = int(lineoriginx - linescale*(lineoriginx - c_centerx));
	    lineendy = int(lineoriginy - linescale*(lineoriginy - c_centery));
	    
	    drawcmd << " line "
		    << lineoriginx << "," << lineoriginy << " "
		    << lineendx    << "," << lineendy;
	}
	
    }

    // scalebar
    int scalebarlength;
    scalebarlength = (int)(height * 1.0 /dimx);
    
    // call imagemagick to write text on image
    if (proj == xaxis)
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" x_forces_%05i.png x_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'X-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						x_proj_%05i.bmp x_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" x_proj_%05i.png x_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}

    }
    else if (proj == yaxis)
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" y_forces_%05i.png y_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Y-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						y_proj_%05i.bmp y_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" y_proj_%05i.png y_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}
    }
    else 
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" z_forces_%05i.png z_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Z-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						z_proj_%05i.bmp z_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" z_proj_%05i.png z_forces_%05i.png",
				color, drawcmd.str().c_str(), filenum, filenum);
		}
    }
    //cout << endl << command1 << endl << endl << command2 << endl;
    system(command1);
    system(command2);

    for(int i=0; i<nucleation_object->n_force_segments(); i++)
    {
	nucleation_object->radial_rep_distrib_x[i]*=InterRecordIterations;
	nucleation_object->radial_rep_distrib_y[i]*=InterRecordIterations;
	nucleation_object->radial_rep_distrib_z[i]*=InterRecordIterations;
    }

}

void actin::draw_bead_forces(int filenum, projection proj,
				int width, int height,
				int xgmax, int ygmax,
				int movex, int movey,
				MYDOUBLE meanx, MYDOUBLE meany, MYDOUBLE meanz,
				MYDOUBLE dimx, MYDOUBLE dimy, MYDOUBLE dimz,
				MYDOUBLE imageGmax, 
				bool nuc_forces_overlay, MYDOUBLE scale, MYDOUBLE thickness,  char*  color)
{


    char command1[2048], command2[2048];
    
    char drawstring[2048] = "";
    char drawtemp[256];
    
    // bead center:
    int x, y;
    int beadcenterx, beadcentery;
    if (proj == xaxis) { // choose projection

	beadcenterx = xgmax+movex + (int)(((height *(((- meany)/dimy) ))-0.5) +  width/2); // was width
	beadcentery = xgmax+movey + (int)(((height *(((- meanz)/dimz) ))-0.5) + height/2);
	
	x = xgmax+movex + (int)(((height *(((RADIUS - meany)/dimy) ))-0.5) +  width/2); // was width
	y = ygmax+movey + (int)(((height *(((- meanz)/dimz) ))-0.5) + height/2);

    } else if (proj == yaxis) {
	
	beadcenterx = xgmax+movex + (int)(((height *(((- meanx)/dimx) ))-0.5) +  width/2); // was width
	beadcentery = xgmax+movey + (int)(((height *(((- meanz)/dimz) ))-0.5) + height/2);
	
	x = xgmax+movex + (int)(((height *(((RADIUS - meanx)/dimx) ))-0.5) +  width/2); // was width
	y = ygmax+movey + (int)(((height *(((       - meanz)/dimz) ))-0.5) + height/2);

    } else {
	
	beadcenterx = xgmax+movex + (int)(((height *(((- meanx)/dimx) ))-0.5) +  width/2); // was width
	beadcentery = xgmax+movey + (int)(((height *(((- meany)/dimy) ))-0.5) + height/2);
	
	x = xgmax+movex + (int)(((height *(((RADIUS - meanx)/dimx) ))-0.5) +  width/2); // was width
	y = ygmax+movey + (int)(((height *(((- meany)/dimy) ))-0.5) + height/2);
    }

    sprintf(drawtemp, "circle %i,%i %i,%i ",beadcenterx, beadcentery, x, beadcentery);
    strcat(drawstring,drawtemp);
    
    MYDOUBLE xscale,yscale;
    int lineoriginx,lineoriginy;
    
    MYDOUBLE radial_rep_tot_x=(MYDOUBLE) 0.01; 
    MYDOUBLE radial_rep_tot_y=(MYDOUBLE) 0.01;
    MYDOUBLE radial_rep_tot_z=(MYDOUBLE) 0.01;
    

    
    for (int i=0; i<RADIAL_SEGMENTS; i++)
    {
	nucleation_object->radial_rep_distrib_x[i]/=InterRecordIterations;
	nucleation_object->radial_rep_distrib_y[i]/=InterRecordIterations;
	nucleation_object->radial_rep_distrib_z[i]/=InterRecordIterations;
    }
    
    for (int i=0; i<RADIAL_SEGMENTS; i++)
    {
	radial_rep_tot_x+=nucleation_object->radial_rep_distrib_x[i];
	radial_rep_tot_y+=nucleation_object->radial_rep_distrib_y[i];
	radial_rep_tot_z+=nucleation_object->radial_rep_distrib_z[i];
    }

 //   MYDOUBLE radial_rep_mean_x,radial_rep_mean_y,radial_rep_mean_z;      
 //   radial_rep_mean_x = radial_rep_tot_x / (MYDOUBLE) RADIAL_SEGMENTS;
 //   radial_rep_mean_y = radial_rep_tot_y / (MYDOUBLE) RADIAL_SEGMENTS;
 //   radial_rep_mean_z = radial_rep_tot_z / (MYDOUBLE) RADIAL_SEGMENTS;
    
    for (int i=0; i<RADIAL_SEGMENTS; i++)
    {
	//cout << nucleation_object->radial_rep_distrib_x[i] << " " ;
	xscale = RADIUS * -sin((2*PI*i/(MYDOUBLE)RADIAL_SEGMENTS));
	yscale = RADIUS *  cos((2*PI*i/(MYDOUBLE)RADIAL_SEGMENTS));
	
	if (proj == xaxis){  // choose projection

	    lineoriginx = xgmax+movex + (int)(((height *(((xscale - meany)/dimy) ))-0.5) +  width/2); // was width
	    lineoriginy = ygmax+movey + (int)(((height *(((yscale - meanz)/dimz) ))-0.5) + height/2);
	    
	    x = xgmax+movex
		+(int)(((height *(((((-nucleation_object->radial_rep_distrib_x[i])*100* scale+1)* xscale - meany)/dimy) ))-0.5)
		       + width/2); // was width
	    y = ygmax+movey
		+ (int)(((height *(((((-nucleation_object->radial_rep_distrib_x[i])*100* scale+1)* yscale - meanz)/dimz) ))-0.5)
			+ height/2);
	    
	} else if (proj == yaxis) {

	    lineoriginx = xgmax+movex + (int)(((height *(((xscale - meanx)/dimx) ))-0.5) +  width/2); // was width
	    lineoriginy = ygmax+movey + (int)(((height *(((yscale - meanz)/dimz) ))-0.5) + height/2);
	    
		x = xgmax+movex
		+ (int)(((height *(((((-nucleation_object->radial_rep_distrib_y[i])*100* scale+1)* xscale - meanx)/dimx) ))-0.5)
			+  width/2); // was width
	    y = ygmax+movey
		+ (int)(((height *(((((-nucleation_object->radial_rep_distrib_y[i])*100* scale+1)* yscale - meanz)/dimz) ))-0.5)
			+ height/2);

	} else {
	    
	    lineoriginx = xgmax+movex + (int)(((height *(((xscale - meanx)/dimx) ))-0.5) +  width/2); // was width
	    lineoriginy = ygmax+movey + (int)(((height *(((yscale - meany)/dimy) ))-0.5) + height/2);
	    
	    x = xgmax+movex
		+ (int)(((height *(((((-nucleation_object->radial_rep_distrib_z[i])*100* scale+1)* xscale - meanx)/dimx) ))-0.5)
			+  width/2); // was width
	    y = ygmax+movey
		+ (int)(((height *(((((-nucleation_object->radial_rep_distrib_z[i])*100* scale+1)* yscale - meany)/dimy) ))-0.5)
			+ height/2);
	}

	sprintf(drawtemp, "line %i,%i %i,%i ",lineoriginx,lineoriginy,x,y);
	strcat(drawstring,drawtemp);
    }
    
    int scalebarlength;
    scalebarlength = (int)(height * 1.0 /dimx);

	//if (!nuc_forces_overlay)  // overlay not coded yet
	//{
	//// #ifdef NO_IMAGE_TEXT
	//	// call imagemagick to write text on image
	//	if (proj == xaxis)  // call imagemagick to write text on image
	//	{
	//	sprintf ( command1 , "mogrify -font helvetica -fill white -pointsize 20 -draw \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'X-Projection  \\nFrame % 4i\\nG-gain % 4i'\" x_proj_%05i.bmp",scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5), filenum );
	//	sprintf ( command2 , "convert -fill none -stroke white -draw \"%s\" x_proj_%05i.bmp x_forces_%05i.bmp", drawstring, filenum, filenum);

	//	} else if (proj == yaxis) {
	//	sprintf ( command1 , "mogrify -font helvetica -fill white -pointsize 20 -draw \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Y-Projection  \\nFrame % 4i\\nG-gain % 4i'\" y_proj_%05i.bmp",scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5), filenum );
	//	sprintf ( command2 , "convert -fill none -stroke white -draw \"%s\" y_proj_%05i.bmp y_forces_%05i.bmp", drawstring, filenum, filenum);

	//	} else {
	//	
	//	sprintf ( command1 , "mogrify -font helvetica -fill white -pointsize 20 -draw \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Z-Projection  \\nFrame % 4i\\nG-gain % 4i'\" z_proj_%05i.bmp",scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5), filenum );
	//	sprintf ( command2 , "convert -fill none -stroke white -draw \"%s\" z_proj_%05i.bmp z_forces_%05i.bmp", drawstring, filenum, filenum);
	//	
	//	}

	//}

// #endif

    //cout << endl << command1 << endl << endl << command2 << endl << endl;
    
    /*if (proj == xaxis)
    {
	sprintf(command1,
		"convert -font helvetica -fill %s -pointsize 20 -draw \
                   \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'X-Projection \\nFrame % 4i\\nG-gain % 4i' \"\
                   -fill none -stroke %s -draw \"%s\" x_proj_%05i.bmp x_forces_%05i.bmp", color,
		scalebarlength, filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5), color, drawstring, filenum, filenum);
	
	sprintf(command2,
		"mogrify -font helvetica -fill %s -pointsize 20 -draw \
                 \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'X-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
                 x_proj_%05i.bmp", color,
		scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum );

    } else if (proj == yaxis) {
	
	sprintf(command1,
		"convert -font helvetica -fill %s -pointsize 20 -draw \
                 \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Y-Projection \\nFrame % 4i\\nG-gain % 4i' \" \
                 -fill none -stroke %s -draw \"%s\" y_proj_%05i.bmp y_forces_%05i.bmp", color,
		scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),color, drawstring, filenum, filenum);
	sprintf(command2,
		"mogrify -font helvetica -fill %s -pointsize 20 -draw \" \
                 text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Y-Projection \\nFrame % 4i\\nG-gain % 4i'\" \
                 y_proj_%05i.bmp", color,
		scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum );

    } else {
	sprintf (command1,
		 "convert -font helvetica -fill %s -pointsize 20 -draw \
                  \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Z-Projection \\nFrame % 4i\\nG-gain % 4i' \" \
                  -fill none -stroke %s -draw \"%s\" z_proj_%05i.bmp z_forces_%05i.bmp", color,
		 scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),color, drawstring, filenum, filenum);
	
	sprintf(command2,
		"mogrify -font helvetica -fill %s -pointsize 20 -draw \
                 \"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Z-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
                 z_proj_%05i.bmp", color,
		scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5), filenum );
    }*/
    
    // call imagemagick to write text on image
    if (proj == xaxis)
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" x_forces_%05i.png x_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'X-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						x_proj_%05i.bmp x_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" x_proj_%05i.png x_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}

    }
    else if (proj == yaxis)
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" y_forces_%05i.png y_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Y-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						y_proj_%05i.bmp y_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" y_proj_%05i.png y_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}
    }
    else 
    {

		if (nuc_forces_overlay)
		{
			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" z_forces_%05i.png z_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}
		else
		{
			sprintf(command1,
				"convert -font helvetica -fill white -pointsize 20 -draw \
						\"text 5 595 '1uM' rectangle 5 576 %i 573 text +5+20 'Z-Projection  \\nFrame % 4i\\nG-gain % 4i'\" \
						z_proj_%05i.bmp z_proj_%05i.png",  
				scalebarlength,filenum,(int)((1000/(MYDOUBLE)imageGmax)+0.5),  filenum, filenum );

			sprintf(command2,
				"convert -fill none -stroke %s -draw \"%s\" z_proj_%05i.png z_forces_%05i.png",
				color, drawstring, filenum, filenum);
		}
    }
    //cout << endl << command1 << endl << endl << command2 << endl;

#ifndef NO_IMAGEMAGICK
    system(command1);
    system(command2);
#endif
    

	// put the radalsegments data back (for next time if doing multiple projections...)

	for (int i=0; i<RADIAL_SEGMENTS; i++)
	{
		nucleation_object->radial_rep_distrib_x[i]*=InterRecordIterations;
		nucleation_object->radial_rep_distrib_y[i]*=InterRecordIterations;
		nucleation_object->radial_rep_distrib_z[i]*=InterRecordIterations;
	}

//	return filenum;

}

int actin::squash(MYDOUBLE thickness)
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

/*
int actin::repulsiveforces(void)  // this is defunct
{
	MYDOUBLE xdist, ydist, zdist, dist;
	MYDOUBLE force, distsqr;
	vect nodeposvec;
	
	int numthreadnodes, start, end;

	numthreadnodes = highestnodecount / NUM_THREADS;

	for (int i=0; i<highestnodecount; i++)
	{
		repdonenode[i] = false;
	}

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
			}
			else
			{	// put remainder in last thread (cludge for now)
				collision_thread_data_array[i].endnode = end + highestnodecount % NUM_THREADS;
			}
		}

		// release thread blocks to start threads:
	for (int i = 0; i < NUM_THREADS; i++)
		{
//			sem_post(&compressfiles_thread_go[i]);
		}

	}
	else
	{
		// go through all nodes
		for (int i=0; i<highestnodecount; i++)
		{  	
			if ((!node[i].polymer) || (node[i].harbinger))
				continue;  // no forces act on harbingers (or depolymerized nodes)

			if (findnearbynodes(nodesbygridpoint[i],NODE_REPULSIVE_RANGE_GRIDSEARCH,0)==0) continue;	// find nodes within 1 grid point
												// skip if zero
			//tmpnodeptr = &node[i];  // for debugging
			//sameGPnode= &tmpnodeptr;

			// loop over nodes on same gridpoint:
			for (vector <nodes*>::iterator sameGPnode=nodes_on_same_gridpoint[0].begin(); sameGPnode<nodes_on_same_gridpoint[0].end() ; sameGPnode++ )
			{
				if ((!(*sameGPnode)->polymer) || ((*sameGPnode)->harbinger))
				continue;  // no forces act on harbingers (or depolymerized nodes)

				for (vector <nodes*>::iterator nearnode=recti_near_nodes[0].begin(); nearnode<recti_near_nodes[0].end() ; nearnode++ )
				{
					if (((*sameGPnode)==(*nearnode)) || (repdonenode[(*sameGPnode)->nodenum]==true))
						 continue;  // skip if self or if done
					
					// this check should not be necessary?:

					xdist = (*nearnode)->x - (*sameGPnode)->x;
					ydist = (*nearnode)->y - (*sameGPnode)->y;
					zdist = (*nearnode)->z - (*sameGPnode)->z;

					distsqr = xdist*xdist + ydist*ydist + zdist*zdist;

					if (distsqr < SQRT_ACCURACY_LOSS)
					{	
						repdonenode[(*sameGPnode)->nodenum] = true;  // mark node as done
						continue;
					}

					if (( distsqr < NODE_REPULSIVE_RANGE*NODE_REPULSIVE_RANGE) &&
						 (distsqr > NODE_INCOMPRESSIBLE_RADIUS*NODE_INCOMPRESSIBLE_RADIUS))
								// must be between NODE_INCOMPRESSIBLE_RANGE and NODE_INCOMPRESSIBLE_RADIUS
					{
						dist = calcdist(xdist,ydist,zdist);  // calc dist between nodes

						force =  NODE_REPULSIVE_MAG - (NODE_REPULSIVE_MAG / 
							(NODE_REPULSIVE_RANGE - NODE_INCOMPRESSIBLE_RADIUS) * dist) ;

						(*sameGPnode)->rep_force_vec[0].x -= force * (xdist/dist);
						(*sameGPnode)->rep_force_vec[0].y -= force * (ydist/dist);
						(*sameGPnode)->rep_force_vec[0].z -= force * (zdist/dist);

					}

					repdonenode[(*sameGPnode)->nodenum] = true;  // mark node as done
				}
			}
			}
	
	}
	return 0;
}
void * actin::repulsiveforcesthread(void* threadarg)
{
	struct thread_data *dat;
	dat = (struct thread_data *) threadarg;

	MYDOUBLE xdist, ydist, zdist;
	MYDOUBLE distsqr, dist, force;
	vect nodeposvec;

	//cout << "Thread " << dat->threadnum << " started" << endl;
	//cout.flush();

	while (true)
	{	

//		sem_wait(&compressfiles_thread_go[dat->threadnum]);  // go only if released by main thread

		//applyforcesdetectiondowork(dat);

		for (int i=dat->startnode; i<dat->endnode; i++)
			{  	
			if ((!node[i].polymer) || (node[i].harbinger))
				continue;  // no forces act on harbingers (or depolymerized nodes)

			if (findnearbynodes(nodesbygridpoint[i],NODE_REPULSIVE_RANGE_GRIDSEARCH,dat->threadnum+NUM_THREADS)==0) continue;	// find nodes within 1 grid point
												// skip if zero
			//tmpnodeptr = &node[i];  // for debugging
			//sameGPnode= &tmpnodeptr;

			// loop over nodes on same gridpoint:
			for (vector <nodes*>::iterator sameGPnode=nodes_on_same_gridpoint[dat->threadnum+NUM_THREADS].begin(); sameGPnode<nodes_on_same_gridpoint[dat->threadnum+NUM_THREADS].end() ; sameGPnode++ )
			{
				if (( (*sameGPnode)->nodenum < dat->startnode) ||
					( (*sameGPnode)->nodenum >= dat->endnode) )
					continue;    // skip if out of range of this thread

				// of these nodes, calculate euclidian dist

				nodeposvec.x = (*sameGPnode)->x;	// get xyz of our node
				nodeposvec.y = (*sameGPnode)->y;
				nodeposvec.z = (*sameGPnode)->z;

				for (vector <nodes*>::iterator nearnode=recti_near_nodes[dat->threadnum].begin(); nearnode<recti_near_nodes[dat->threadnum].end() ; nearnode++ )
				{
					if ((*sameGPnode)==(*nearnode)) continue;  // skip if self
					//if (repulsedone[(*nearnode)->nodenum][(*sameGPnode)->nodenum])
					//	continue;  // skip if done opposite pair
					
					// this check should not be necessary?:

					if (( (*sameGPnode)->nodenum < dat->startnode) ||
						( (*sameGPnode)->nodenum >= dat->endnode) )
						continue;    // skip if out of range of this thread  

					xdist = (*nearnode)->x - nodeposvec.x;
					ydist = (*nearnode)->y - nodeposvec.y;
					zdist = (*nearnode)->z - nodeposvec.z;

					distsqr = xdist*xdist + ydist*ydist + zdist*zdist;

					if (distsqr < SQRT_ACCURACY_LOSS)
					{	
						repdonenode[(*sameGPnode)->nodenum] = true;  // mark node as done
						continue;
					}

					if ( distsqr < NODE_REPULSIVE_RANGE*NODE_REPULSIVE_RANGE)
					{
						dist = calcdist(xdist,ydist,zdist);  // calc dist between nodes

						force =  NODE_REPULSIVE_MAG 
						- (NODE_REPULSIVE_MAG / (NODE_REPULSIVE_RANGE - NODE_INCOMPRESSIBLE_RADIUS) ) * dist;

						(*sameGPnode)->rep_force_vec[dat->threadnum].x += force * (xdist/dist);
						(*sameGPnode)->rep_force_vec[dat->threadnum].y += force * (ydist/dist);
						(*sameGPnode)->rep_force_vec[dat->threadnum].z += force * (zdist/dist);

					}

					repdonenode[(*sameGPnode)->nodenum] = true;  // mark node as done
				}
			}
			}
		}

//		sem_post(&compressfiles_data_done[dat->threadnum]);  // release main thread block waiting for data

	  
	return (void*) NULL;
}
*/

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

	while (true)
	{	
		//sem_wait(&compressfiles_thread_go[0]);  // go only if released by main thread

		pthread_mutex_lock(&filessavelock_mutex);
		pthread_mutex_lock(&filesdonelock_mutex);

		char command1[255];//,command2[255];

		// gzip the text data:

		sprintf(command1, "gzip -f -9 nodes%05i.txt",dat->startnode);
		system(command1);

		sprintf(command1, "gzip -f -9 links%05i.txt",dat->startnode);
		system(command1);

		sprintf(command1, "gzip -q -f -9 report*.txt 2>/dev/null" );
		system(command1);

		sprintf(command1 , "gzip -q -f -9 data*.txt 2>/dev/null" );
		system(command1);

		// gzip the wrl file:

		sprintf(command1 , "gzip -f -9 nodes%05i.wrl",dat->startnode );
		system(command1);

		sprintf(command1 , "mv nodes%05i.wrl.gz nodes%05i.wrz",dat->startnode, dat->startnode );
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

int actin::find_center(vect &center)
{
	center = - nucleation_object->position;

	return 0;
}

int actin::save_linkstats(int filenum)
{
	char filename[255];
//	MYDOUBLE radial_force, transverse_force;

	sprintf ( filename , "linkstats%05i.txt", filenum );

	ofstream oplinkstats(filename, ios::out | ios::trunc);
	if (!oplinkstats) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	// write header
	
	oplinkstats << "Link Radius,Radial Force,Transverse Force" << endl;
	for (int i=0; i<highestnodecount; i++)
	{
		if ((node[i].polymer) && (!node[i].harbinger))
		{
            //link_force_vec[0].x
		}
	}


	return 0;
}
void actin::reportsnapshot(int filenum, int highestnode, int reportiteration)
{

 	MYDOUBLE dist,compx,compy,compz;

	//int varnum;

    for (int i=0; i<highestnode; i++)
	{
		if ((node[i].polymer) && (!node[i].harbinger))  // is point valid?
		{

			reportdat[0][reportiteration][i]=1;  // valid point
			
			//varnum = 0;

			// calculate distance from center of bead
			dist = calcdist(node[i].x,node[i].y,node[i].z);

			// calculate components of normalized radial vector
			compx = node[i].x / dist;
			compy = node[i].y / dist;
			compz = node[i].z / dist;

			// distance from center of bead
			reportdat[1][reportiteration][i]=dist;

			// dot product of repulsive force vector and normalised radial vector
			reportdat[2][reportiteration][i]=(node[i].rep_force_vec[0].x * compx +
						node[i].rep_force_vec[0].y * compy +
						node[i].rep_force_vec[0].z * compz );

			// cross product of repulsive force vector and normalised radial vector
			reportdat[3][reportiteration][i]=calcdist(node[i].rep_force_vec[0].y * compz - node[i].rep_force_vec[0].z * compy,
						node[i].rep_force_vec[0].z * compx - node[i].rep_force_vec[0].x * compz,
						node[i].rep_force_vec[0].x * compy - node[i].rep_force_vec[0].y * compx);

			// dot product of repulsive displacement vector and normalised radial vector
			reportdat[4][reportiteration][i]=(node[i].repulsion_displacement_vec[0].x * compx +
						node[i].repulsion_displacement_vec[0].y * compy +
						node[i].repulsion_displacement_vec[0].z * compz );

			// cross product of repulsive displacement vector and normalised radial vector
			reportdat[5][reportiteration][i]=calcdist(node[i].repulsion_displacement_vec[0].y * compz - node[i].repulsion_displacement_vec[0].z * compy,
						node[i].repulsion_displacement_vec[0].z * compx - node[i].repulsion_displacement_vec[0].x * compz,
						node[i].repulsion_displacement_vec[0].x * compy - node[i].repulsion_displacement_vec[0].y * compx);

			// dot product of link force vector and normalised radial vector
			reportdat[6][reportiteration][i]=(node[i].link_force_vec[0].x * compx +
						node[i].link_force_vec[0].y * compy +
						node[i].link_force_vec[0].z * compz );

			// cross product of link force vector and normalised radial vector
			reportdat[7][reportiteration][i]=calcdist(node[i].link_force_vec[0].y * compz - node[i].link_force_vec[0].z * compy,
						node[i].link_force_vec[0].z * compx - node[i].link_force_vec[0].x * compz,
						node[i].link_force_vec[0].x * compy - node[i].link_force_vec[0].y * compx);
		}
		else
		{
			reportdat[0][reportiteration][i]=-1;  // not valid point
		}
	}

	return;

}

void actin::savereport(int filenum, int highestnode)
{

	char filename[255];

	sprintf ( filename , "report%05i.txt", filenum );

	ofstream opreport(filename, ios::out | ios::trunc);
	if (!opreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

	// write header
	
	opreport << "Nodenum,NodeRadius,NodeRadiusSD,RepForceRadial,RepForceRadialSD,RepForceTrans,RepForceTransSD,RepDisplRadial,RepDisplRadialSD,RepDisplTrans,RepDisplTransSD,LinkForceRadial,LinkForceRadialSD,LinkForceTrans,LinkForceTransSD" << endl;
	
	double sum, mean, sd, varsum, count;

	for(int i=1; i<highestnode; ++i)
	{
		opreport << i << ",";

		for(int j=1; j<REPORT_NUM_VARIABLES; ++j)  // skip i=0 since this is the 'valid' switch
		{
			sum = 0;
			varsum = 0;
			count = 0;

			for (int k=0; k<REPORT_AVERAGE_ITTERATIONS; k++)
			{

				if (reportdat[0][k][i] > 0)
				{ // valid point
					++count;
					sum+=reportdat[j][k][i];
				}
			}

			mean = sum/count;

			for (int k=0; k<REPORT_AVERAGE_ITTERATIONS; k++)
			{

				if (reportdat[0][k][i] > 0)
				{ // valid point
					varsum+=pow((mean - reportdat[j][k][i]),2);
				}
			}

			sd = varsum/count;

			opreport << mean << "," << sd << "," ;
		}

		opreport << endl;
	}

	opreport.close();

}

int actin::savedata(int filenum)
{
	char filename[255];

	sprintf ( filename , "data%05i.txt", filenum );

	ofstream opdata(filename, ios::out | ios::trunc);
	if (!opdata) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	opdata << highestnodecount << "," << iteration_num 
		<< "," << linksbroken << "," << linksformed << endl;

	for (int i=0; i<highestnodecount; i++)
	{
		node[i].savedata(&opdata);
	}

	opdata << endl;
	opdata << (unsigned int) crosslinknodesdelay.size();

	for (vector <int>::iterator i=crosslinknodesdelay.begin(); i<crosslinknodesdelay.end() ; i++ )
	{	 
		opdata << "," << (*i);
	}

	opdata.close();

	return 0;
}

int actin::loaddata(int filenum)
{

	for (int i=0; i!=(GRIDSIZE+1); i++) // clear the nodegrid
		for (int j=0; j!=(GRIDSIZE+1); j++)
			for (int k=0; k!=(GRIDSIZE+1); k++)
				nodegrid[i][j][k]=0;;

	char filename[255];
	char delim;
	unsigned int numcrosslinkdelay;

	sprintf ( filename , "data%05i.txt", filenum );

	ifstream ipdata(filename);
	if (!ipdata) 
	{ cout << "Unable to open file " << filename << " for output"; return 1;}

	ipdata >> highestnodecount >> delim >> iteration_num 
		>> delim >> linksbroken >> delim >> linksformed;

	for (int i=0; i<highestnodecount; i++)
	{
		node[i].loaddata(&ipdata);
	}

	ipdata >> numcrosslinkdelay;
	crosslinknodesdelay.resize(numcrosslinkdelay);

	for (vector <int>::iterator i=crosslinknodesdelay.begin(); i<crosslinknodesdelay.end() ; i++ )
	{	 
		ipdata >> delim >> (*i);
	}

	return 0;
}

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
	   << linksformed << endl;
    
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
    nucleation_object->save_data(ofstrm);
    
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
	  >> nexttocrosslink >> ch
	  >> iteration_num >> ch 
	  >> linksbroken  >> ch 
	  >> linksformed;
    
    // load nodes
    ifstr >> str;
    if(str.compare("nodes-links:") !=0 ){
	cout << "error in data file, 'nodes-links:' expected" << endl;
	return 1;
    }
    // ** Remember the node vector is preallocated to MAXNODES
    for(int i=0; i < highestnodecount; i++) {
	node[i].load_data(ifstr);
	node[i].setgridcoords(); // now we are okay to do pointer
	node[i].addtogrid();
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
    crosslinknodesdelay.clear();
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
    nucleation_object->load_data(ifstr);
    
    
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

rotationmatrix actin::get_sym_break_axes(void)
{

	vect tmp_nodepos;

	vect  CofM;

	vect sym_break_direction;

	rotationmatrix tmp_rotation, final_rotation;	

	MYDOUBLE x_angle, y_angle, z_angle;

	MYDOUBLE x_inertia = 0, y_inertia = 0;
	
	MYDOUBLE z_in_angle;

	int numnodes;

	find_center(sym_break_direction);  // which way did the bead go?

	x_angle = atan2(sym_break_direction.y,sym_break_direction.z);
	y_angle = atan2(sym_break_direction.x,sym_break_direction.z);
	z_angle = atan2(sym_break_direction.x,sym_break_direction.y);

	// make a temp rotation matrix in that direction:

	//tmp_rotation.rotatematrix(0, y_angle, 0); // just rotate y
	tmp_rotation.rotatematrix(x_angle, y_angle, 0); // now pointing upwards

	//cout << "old sym_break_direction: " 
	//	<< sym_break_direction.x << "," 
	//	<< sym_break_direction.y << "," 
	//	<< sym_break_direction.z << endl;

	//tmp_rotation.rotate(sym_break_direction);

	//cout << "new sym_break_direction: " 
	//	<< sym_break_direction.x << "," 
	//	<< sym_break_direction.y << "," 
	//	<< sym_break_direction.z << endl;


	// find CofM

	numnodes=0;

	for (int i=0; i<highestnodecount; ++i)
	{
		if (!node[i].polymer)
			continue;

		CofM += node[i];

		++numnodes;

	}

	CofM *= (1/ (MYDOUBLE) numnodes);

	tmp_rotation.rotate(CofM);


	for (int i=0; i<highestnodecount; ++i)
	{
		if (!node[i].polymer)
			continue;

		tmp_nodepos = node[i];
		tmp_rotation.rotate(tmp_nodepos);  // rotate points to bring in line with sym break dir'n


		// determine principal axis

		x_inertia+=(tmp_nodepos.x-CofM.x)*(tmp_nodepos.x-CofM.x);
		y_inertia+=(tmp_nodepos.y-CofM.y)*(tmp_nodepos.y-CofM.y);
		
	}

	z_in_angle = atan2(sqrt(x_inertia),sqrt(y_inertia));

	cout << "sqrt(x_inertia): " << sqrt(x_inertia) 
		 << " sqrt(y_inertia): " << sqrt(y_inertia) 
		 << " z_in_angle: " << 180*z_in_angle << endl; 

	//final_rotation.rotatematrix(x_angle, y_angle, -z_in_angle);
	final_rotation.rotatematrix(0, 0, z_in_angle);
	final_rotation.rotatematrix(x_angle, y_angle, 0);

	
	return final_rotation;
}

void actin::clearstats(void)
{
	for (int threadnum = 0; threadnum < NUM_THREADS; ++threadnum)
	{
		for (int i=0; i<highestnodecount; ++i)
		{
			node[i].clearstats(threadnum);
		}
	}

}
