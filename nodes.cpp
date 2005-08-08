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
#include "nodes.h"

nodes::nodes(void)
{
	nextnode = this;  // initialise to point to self
	prevnode = this;
	
	x = y = z = 0.0;
	lastpos.zero();
	gridx = gridy = gridz = -1;  //note this is meaningless unless polymer==true
	//rep_force_vec.x  = rep_force_vec.y  = rep_force_vec.z  = 0.0;
	
	//repulsion_displacement_vec = new vect[NUM_THREADS+1];
	repulsion_displacement_vec.resize(NUM_THREADS+1);
	link_force_vec.resize(NUM_THREADS+1);
	rep_force_vec.resize(NUM_THREADS+1);

	for (int i = 0; i < NUM_THREADS; ++i)
	{
		repulsion_displacement_vec[i].zero();
		link_force_vec[i].zero();
		rep_force_vec[i].zero();
	}

	unit_vec_posn.zero();

	
//	momentum_vec.x = momentum_vec.y = momentum_vec.z = 0.0;	
	polymer = false;
	colour.setcol(0);
	creation_iter_num = 0;
	nodelinksbroken = 0;
	harbinger = false;
	dontupdate = false;

	linkforce_transverse.resize(NUM_THREADS);
	linkforce_radial.resize(NUM_THREADS);
	repforce_transverse.resize(NUM_THREADS);
	repforce_radial.resize(NUM_THREADS);
	dispforce_transverse.resize(NUM_THREADS);
	dispforce_radial.resize(NUM_THREADS);

	for (int threadnum = 0; threadnum < NUM_THREADS; ++threadnum)
	{
		clearstats(threadnum);
	}
}

nodes::nodes(const double& set_x, const double& set_y,const double& set_z)
{
	gridx = gridy = gridz = -1;	
    nextnode = this;  // initialise to point to self
	prevnode = this;
	//rep_force_vec.x  = rep_force_vec.y  = rep_force_vec.z  = 0.0;
	lastpos.zero();
	//repulsion_displacement_vec =  new vect[NUM_THREADS+1];

	repulsion_displacement_vec.resize(NUM_THREADS+1);
	link_force_vec.resize(NUM_THREADS+1);
	rep_force_vec.resize(NUM_THREADS+1);

	for (int i = 0; i < NUM_THREADS; ++i)
	{
		repulsion_displacement_vec[i].zero();
		link_force_vec[i].zero();
		rep_force_vec[i].zero();
	}

	for (int threadnum = 0; threadnum < NUM_THREADS; ++threadnum)
	{
		clearstats(threadnum);
	}

//	momentum_vec.x = momentum_vec.y = momentum_vec.z = 0.0;	
	colour.setcol(0);
	creation_iter_num = ptheactin->iteration_num;
	nodelinksbroken =0;
	harbinger = false;
	dontupdate = false;

	polymerize(set_x,  set_y,  set_z);
	
}

nodes::~nodes(void)
{
//	delete[] repulsion_displacement_vec;
}

bool nodes::depolymerize(void) 
{
	removefromgrid();
	x = y = z = 0.0;
	polymer = false;

	// remove all back links
	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		i->linkednodeptr->removelink(this);
	}

	// and remove our own links to other nodes
	ptheactin->linksbroken+=(int) listoflinks.size();
	listoflinks.resize(0);

	return false;
}

bool nodes::polymerize(const double& set_x, const double& set_y, const double& set_z)
{
	x = set_x;
	y = set_y;
	z = set_z;
	
	lastpos=*this; 

	polymer = true ;

	nextnode = this;  // initialise to point to self
	prevnode = this;

	setgridcoords(); // set grid by x,y,z
	addtogrid();     // add node to the grid
	unit_vec_posn=this->unitvec();

	colour=ptheactin->newnodescolour;

	creation_iter_num = ptheactin->iteration_num;

	harbinger = true;

	//colour.setcol(x);

	return true;
}

int nodes::save(ofstream *outputstream) 
{
	*outputstream	<< nodenum << "," << x << "," << y << "," << z 
					<< "," << harbinger << "," << polymer << "," 
					<< colour.r << "," << colour.g << "," << colour.b << ","
					<< creation_iter_num << (unsigned int) listoflinks.size() 
				 << "," << endl;

	return 0;
}

int nodes::savedata(ofstream *outputstream) 
{
	*outputstream	<< nodenum << "," << x << "," << y << "," << z 
					<< "," << harbinger << "," << polymer << "," 
					<< colour.r << "," << colour.g << "," << colour.b << ","
					<< creation_iter_num << "," << (unsigned int) listoflinks.size();
					
	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		*outputstream << ",";
		i->savedata(outputstream);
	}
	
	*outputstream << endl;

	return 0;
}

int nodes::loaddata(ifstream *inputstream) 
{
	char delim;
	unsigned int linklistsize;

	*inputstream	>> nodenum >> delim >> x >> delim >> y >> delim >> z 
					>> delim >> harbinger >> delim >> polymer >> delim 
					>> colour.r >> delim >> colour.g >> delim >> colour.b >> delim
					>> creation_iter_num >> delim >> linklistsize ;
					
	listoflinks.resize(linklistsize);  // is resize right? or should it be reserve?

	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		*inputstream >> delim;
		i->loaddata(inputstream);
	}

	setgridcoords();
	addtogrid(); 

	return 0;
}


int nodes::save_data(ofstream &ostr) 
{
    // save the nodes
    ostr << nodenum << "," 
	 << x << "," << y << "," << z << "," 
	 << harbinger << "," << polymer << "," 
	 << colour.r << "," << colour.g << "," << colour.b << ","
	 << creation_iter_num << ":";
    
    // now the links
    ostr << (unsigned int) listoflinks.size() << ")";
    for(vector<links>::iterator l=listoflinks.begin(); l<listoflinks.end(); ++l)
    {
	l->save_data(ostr);
	ostr << ".";
    }
    
    // done
    ostr << endl;
    
    return 0;
}

int nodes::load_data(ifstream &istrm) 
{
    // read in from the stream to our private data
    char ch;    
    istrm >> nodenum >> ch 
	  >> x >> ch >> y >> ch >> z >> ch
	  >> harbinger >> ch >> polymer >> ch
	  >> colour.r >> ch >> colour.g >> ch >> colour.b >> ch
	  >> creation_iter_num >> ch;
    
    // check we are ready to read links
    if(ch!=':' ){
	cout << "error in checkpoint file, end of node ':' expected" 
	     << endl;
	return 1;
    }
    
    int linklistsize;
    istrm >> linklistsize >> ch;
    if( ch!=')' ){
	cout << "error in checkpoint file, xlinkdelays 'NN)' expected" 
	     << endl;
	return 1;
    }
  
    listoflinks.clear(); // note doesn't free memory
    listoflinks.resize(linklistsize);

    for(int i=0; i<linklistsize; ++i)
    {	 
	// construct the links and add to vector
	links link;
	link.load_data(istrm);
	listoflinks[i] = link;
	istrm >> ch;
    }
    // note we don't set pointer or build the grid here
    // because we need to be sure this is the node
    // stored by the actin.  The actin knows that.
    return 0;
}


void nodes::updategrid(void)
{
	int gridtmpx, gridtmpy, gridtmpz;
	int oldgridx = gridx;  // store old grid pos'n
	int oldgridy = gridy;
	int oldgridz = gridz;

	setgridcoords();  // set gridx,y,z by x,y,z position
	unit_vec_posn=this->unitvec();  //set unit vector position

	if	(((gridx!=oldgridx) ||     // has the node moved gridpoints?
		  (gridy!=oldgridy) || 
		  (gridz!=oldgridz)))
	{
		// node moved, check not out of grid bounds
		if ((gridx>GRIDSIZE) ||
			(gridy>GRIDSIZE) ||
			(gridz>GRIDSIZE) ||
			(gridx<0) ||
			(gridy<0) ||
			(gridz<0)) 
		{  
			cout << "Node out of grid bounds, deleted:" << x << " " << y << " " << z << endl;
			gridx=oldgridx;
			gridy=oldgridy;
			gridz=oldgridz;
			depolymerize();
			return;
		}
	}

	else
	{  // not moved - return
		return;
	}


	// move the node in the grid array

	// is the node on the grid:
	if (oldgridx!=-1)
	{
		// if so, remove node from old grid pos'n

		gridtmpx = gridx; gridx = oldgridx;  // need to remove from old pos'n
		gridtmpy = gridy; gridy = oldgridy;
		gridtmpz = gridz; gridz = oldgridz;

		removefromgrid();

		gridx = gridtmpx;  // restore new grid co-ords
		gridy = gridtmpy;
		gridz = gridtmpz;
		
	}

	// add self to new grid position if not depolymerized

	if (polymer) 
		addtogrid();
			
	return;
}

void nodes::removefromgrid(void)
{
	// are we on the grid?
	if (gridx==-1) return;  // return if now

	// are we the only grid node?
		if ((nextnode==this) &&
			(nodegrid[gridx][gridy][gridz] == this))
		{	// if so, just delete the grid reference
			nodegrid[gridx][gridy][gridz] = 0;
		}
		else
		{	// other nodes on grid
			if (nodegrid[gridx][gridy][gridz] == this)
			{  // if we're the grid reference set ref to next node
				nodegrid[gridx][gridy][gridz] = nextnode;
			}

            nextnode->prevnode = prevnode;  //  remove self from circular list
			prevnode->nextnode = nextnode;

		}

		gridx=gridy=gridz=-1;

	return;
}

void nodes::addtogrid(void)
{
	// are we already on the grid?
	//if (gridx!=-1) return 0;

	// is the new grid node empty?
	if ((nodegrid[gridx][gridy][gridz] == 0))
	{	// if so, just add self to the grid reference
		nodegrid[gridx][gridy][gridz] = this;
		nextnode = prevnode = this;  // and loop to self
	}
	else
	{	// otherwise sew into loop
		nextnode = nodegrid[gridx][gridy][gridz];  // our next is the grid
		prevnode = nextnode->prevnode;  //our new previous is the new next's old previous

		nextnode->prevnode = this;  // and we are theirs
		prevnode->nextnode = this;
	}

	return;
}

void nodes::setgridcoords(void)
{  

	gridx = (((int)(x / GRIDRES)) + (GRIDSIZE/2) ); // find new grid pos'n
	gridy = (((int)(y / GRIDRES)) + (GRIDSIZE/2) );
	gridz = (((int)(z / GRIDRES)) + (GRIDSIZE/2) );

	return;
} 

int nodes::addlink(nodes* linkto, const double& dist)
{

	//if (listoflinks.size()<MAX_LINKS_PER_NODE)
	//{
		listoflinks.push_back(links(linkto,dist));
		ptheactin->linksformed++;
		return true;
	//}

	//return false;
}

int nodes::removelink(nodes* link)
{
	//  check node-nucleator repulsion

	//int templinksbroken = ptheactin->linksbroken;

	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		if (i->linkednodeptr==link)
		{
			listoflinks.erase(i);
			ptheactin->linksbroken++;
			continue;
		}
	}

nodelinksbroken++;

//if (templinksbroken == ptheactin->linksbroken)
//{
//	// didn't remove the link for some reason
//	cout << "Tried to remove but link not in list " << link->nodenum << " " <<endl; 
//}
	
return 0;
}

int nodes::savelinks(ofstream * outputstream)
{
	if (listoflinks.empty())
		return 0;

	*outputstream << nodenum << "," << (unsigned int) listoflinks.size();

	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		*outputstream << "," << i->linkednodeptr->nodenum;
	}

	*outputstream << endl;

	return 0;
}
