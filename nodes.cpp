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

nodes::nodes(void) :
onseg(false),
polymer(false),
harbinger(true),
stucktonucleator(false),
move_harbinger_this_time(false),
threadnum(0),
gridx(-1),gridy(-1),gridz(-1),
nodelinksbroken(0),
creation_iter_num(0)
//, ptheactin(NULL)
//: local_NODE_FORCE_TO_DIST(DELTA_T * FORCE_SCALE_FACT)
{
	nextnode = this;  // initialise to point to self
	prevnode = this;
	
	x = y = z = 0.0;
	//gridx = gridy = gridz = -1;  //note this is meaningless unless polymer==true
	
    link_force_vec.zero();
	rep_force_vec.zero();
    //viscous_force_vec.zero();
	nuc_repulsion_displacement_vec.zero();
    //viscous_force_recip_dist_sum = 0;

	unit_vec_posn.zero();

    nucleator_link_force.zero();

	//threadnum = 0;

	listoflinks.reserve(MAX_LINKS_PER_NEW_NODE*2);

    //insidenucleator = false;
	//polymer = false;
	colour.setcol(0);
	//creation_iter_num = 0;
	//nodelinksbroken = 0;
	//harbinger = true;
//	dontupdate = false;

	//move_harbinger_this_time = false;

    //stucktonucleator = false;
    
	nucleator_stuck_position.zero();

	clearforces();
	clearstats();

}

nodes::nodes(const double& set_x, const double& set_y,const double& set_z)
//: ptheactin(pactin)
//: local_NODE_FORCE_TO_DIST(DELTA_T * FORCE_SCALE_FACT)
{
	gridx = gridy = gridz = -1;	
    nextnode = this;  // initialise to point to self
	prevnode = this;

	link_force_vec.zero();
	rep_force_vec.zero();
    //viscous_force_vec.zero();
	nuc_repulsion_displacement_vec.zero();
    //viscous_force_recip_dist_sum = 0;

    nucleator_link_force.zero();

	clearstats();

	threadnum = 0;

	listoflinks.reserve(MAX_LINKS_PER_NEW_NODE*2);

	colour.setcol(0);
	creation_iter_num = ptheactin->iteration_num;
	nodelinksbroken =0;
	harbinger = true;
//	dontupdate = false;
//    insidenucleator = false;

	polymerize(set_x,  set_y,  set_z);

    stucktonucleator = false;
    nucleator_stuck_position.zero();

	clearforces();

	move_harbinger_this_time = false;

}											    

nodes::~nodes(void)
{
}





bool nodes::depolymerize(void) 
{	
	removefromgrid();
	//ptheactin->removenodefromgrid(this);
	
	x = y = z = 0.0;
	polymer = false;

	// remove all back links
	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		i->linkednodeptr->removelink(this);
	}

	// and remove our own links to other nodes
	ptheactin->linksbroken += (int) listoflinks.size();
    links_broken += (int) listoflinks.size();
	listoflinks.resize(0);

	return false;
}

bool nodes::polymerize(const double& set_x, const double& set_y, const double& set_z)
{
	x = set_x;
	y = set_y;
	z = set_z;
	
	//lastpos=*this; 

	polymer = true ;

	nextnode = this;  // initialise to point to self
	prevnode = this;

	setgridcoords(); // set grid by x,y,z
	addtogrid();     // add node to the grid
	//ptheactin->addnodetogrid(this);     // add node to the grid
	unit_vec_posn=this->unitvec();

	colour=ptheactin->newnodescolour;

	creation_iter_num = ptheactin->iteration_num;

	harbinger = true;

    setunitvec();

	//colour.setcol(x);

	return true;
}

int nodes::save_data(ofstream &ostr) 
{
    // save the nodes
    ostr << nodenum << " " 
	 << x << " " << y << " " << z << " " 
	 << harbinger << " " 
	 << polymer << " " 
	 << colour.r << " " << colour.g << " " << colour.b << " "
	 << delta.x << " "
	 << delta.y << " "
	 << delta.z << " "
	 << linkforce_transverse << " " 
	 << linkforce_radial << " " 
	 << repforce_transverse << " " 
	 << repforce_radial << " " 
	 << links_broken << " " 
	 << nucleator_impacts << " "
	 << stucktonucleator << " "
	 << nucleator_stuck_position.x << " "
     << nucleator_stuck_position.y << " "
     << nucleator_stuck_position.z << " "
     << nucleator_link_force.x << " "
     << nucleator_link_force.y << " "
     << nucleator_link_force.z << " "
	 << creation_iter_num << ":";
    
    // now the links
    ostr << (unsigned int) listoflinks.size() << ")";
    for(vector<links>::iterator l=listoflinks.begin(); l<listoflinks.end(); ++l)
    {
	l->save_data(ostr);
	ostr << " ";
    }
    
    // done
    ostr << endl;
    
    return 0;
}

int nodes::load_data(ifstream &istrm) 
{
    // read in from the stream to our private data
    char ch;    
    istrm >> nodenum 
	  >> x  >> y  >> z 
	  >> harbinger  
	  >> polymer 
	  >> colour.r  >> colour.g  >> colour.b 
	  >> delta.x 
      >> delta.y 
      >> delta.z 
	  >> linkforce_transverse  
	  >> linkforce_radial  
	  >> repforce_transverse  
	  >> repforce_radial  
	  >> links_broken  
	  >> nucleator_impacts 
	  >> stucktonucleator 
	  >> nucleator_stuck_position.x 
      >> nucleator_stuck_position.y 
      >> nucleator_stuck_position.z 
      >> nucleator_link_force.x 
      >> nucleator_link_force.y 
      >> nucleator_link_force.z 
	  >> creation_iter_num >> ch;

	//if (nucleator_impacts>0.00001)
	//	cout << "Loaded nucleator impact " <<nucleator_impacts << endl;
    
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
	//istrm >> ch;
    }
    // note we don't set pointer or build the grid here
    // because we need to be sure this is the node
    // stored by the actin.  The actin knows that.

    setunitvec();
	updategrid();

    return 0;
}


void nodes::updategrid(void)
{
	int gridtmpx, gridtmpy, gridtmpz;
	int oldgridx = gridx;			// store old grid pos'n
	int oldgridy = gridy;
	int oldgridz = gridz;

	setgridcoords();				// set gridx,y,z by x,y,z position
	//setunitvec();	

	if	((gridx != oldgridx) ||		// has the node moved gridpoints?
		 (gridy != oldgridy) || 
		 (gridz != oldgridz))
	{
		// node moved, check not out of grid bounds
		if ((gridx > GRIDSIZE) ||
			(gridy > GRIDSIZE) ||
			(gridz > GRIDSIZE) ||
			(gridx < 0) ||
			(gridy < 0) ||
			(gridz < 0)) 
		{  
			cout << "Node out of grid bounds, deleted:" << x << " " << y << " " << z << endl;

			gridx = oldgridx;
			gridy = oldgridy;
			gridz = oldgridz;

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
		//ptheactin->removenodefromgrid(this);

		gridx = gridtmpx;  // restore new grid co-ords
		gridy = gridtmpy;
		gridz = gridtmpz;
		
	}

	// add self to new grid position if not depolymerized

	if (polymer)
	{
		addtogrid();
		//ptheactin->addnodetogrid(this);
	}
			
	return;
}

void nodes::removefromgrid(void)
{
	// are we on the grid?
	if (gridx==-1) return;  // return if now

	// are we the only grid node?
		if ((nextnode==this) &&
			(NODEGRID(gridx,gridy,gridz) == this))
		{	// if so, just delete the grid reference
			NODEGRID(gridx,gridy,gridz) = 0;
		}
		else
		{	// other nodes on grid
			if (NODEGRID(gridx,gridy,gridz) == this)
			{  // if we're the grid reference set ref to next node
				NODEGRID(gridx,gridy,gridz) = nextnode;
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
	if ((NODEGRID(gridx,gridy,gridz) == 0))
	{	// if so, just add self to the grid reference
		NODEGRID(gridx,gridy,gridz) = this;
		nextnode = prevnode = this;  // and loop to self
	}
	else
	{	// otherwise sew into loop
		nextnode = NODEGRID(gridx,gridy,gridz);  // our next is the grid
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

	//if (listoflinks.size()<MAX_LINKS_PER_NEW_NODE)
	//{
		listoflinks.push_back(links(linkto,dist));
		(ptheactin->linksformed)++;
		return true;
	//}

	//return false;
}

int nodes::removelink(nodes* link)
{
	//  check node-nucleator repulsion

	//int templinksbroken = ptheactin->linksbroken;

	if (listoflinks.size()==0)
		return 0;

	for (vector<links>::iterator i  = listoflinks.begin();
		                         i != listoflinks.end() ;
								  ++i )
	{	 
		if (i->linkednodeptr==link)
		{
			listoflinks.erase(i);
			ptheactin->linksbroken++;
			break;
			//continue;
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

	*outputstream << nodenum << " " << (unsigned int) listoflinks.size();

	for (vector <links>::iterator i=listoflinks.begin(); i<listoflinks.end() ; i++ )
	{	 
		*outputstream << " " << i->linkednodeptr->nodenum;
	}

	*outputstream << endl;

	return 0;
}
