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
#include "actin.h"

class Colour;

nodes::nodes(void) :
nodegridptr(NULL),
onseg(false),
polymer(false),
harbinger(true),
stucktonucleator(false),
move_harbinger_this_time(false),
threadnum(0),
gridx(-1),gridy(-1),gridz(-1),
links_broken(0)
{   /// Create a blank node

    /// This is the default constructor, called in the main node allocation.  
    /// Before the node joins the network, it must be initialised properly with call to polymerize.
    /// alternatively call the other nodes(x,y,z) constructor.

    threadnum = 0;
    testnode = false;
    testsurface = 0;
    creation_iter_num = 0;

//    initial_repulsive_energy = -1.0; // start -ve before set

  	x = y = z = 0.0;
    link_force_vec.zero();
	rep_force_vec.zero();
    //unit_vec_correct = false;

	nuc_repulsion_displacement_vec.zero();

	unit_vec_posn.zero();

    nucleator_link_force.zero();

    posnoflastgridupdate.zero();

    listoflinks.resize(0);
    listoflinks.reserve(MAX_EXPECTED_LINKS);
    

	//colour.setcol(0);
 
	nucleator_stuck_position.zero();

	clearforces();
	clearstats();

}

nodes::nodes(const double& set_x, const double& set_y,const double& set_z):
onseg(false),
harbinger(true),
stucktonucleator(false),
move_harbinger_this_time(false),
links_broken(0)
{
    threadnum = 0;
    testnode = false;
    testsurface = 0;
    creation_iter_num = 0;
    //unit_vec_correct = false;
//    initial_repulsive_energy = -1.0; // start -ve before set

	link_force_vec.zero();
	rep_force_vec.zero();

    nuc_repulsion_displacement_vec.zero();

    nucleator_link_force.zero();

	clearstats();

    listoflinks.resize(0);
	listoflinks.reserve(MAX_EXPECTED_LINKS);
    

	//colour.setcol(0);
    
    nucleator_stuck_position.zero();
    clearforces();

    nodegridptr = NULL;
    gridx = gridy = gridz = -1;
    creation_iter_num = 0;
    polymer = false;

	polymerize(set_x,  set_y,  set_z);
}											    

nodes::~nodes(void)
{
}



bool nodes::depolymerize(void) 
{	
	removefromgrid();
	
	x = y = z = 0.0;
	polymer = false;
    stucktonucleator = false;

	// remove all back links
	for (vector <links>::iterator i=listoflinks.begin(); i!=listoflinks.end() ; ++i )
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
{   /// Initialises the node at given position (world co-ords) before adding to the network

    /// Effectively this creates a harbinger.  The node isn't crosslinked here, just created so it
    /// can interact with the network (added to the grid, position set etc.)
    /// it is not stuck to the nucleator here, since it's passed world co-ords, and the nucleator stuck
    /// position is relative to the nucleator.

	x = set_x;
	y = set_y;
	z = set_z;

	polymer = true ;

	setgridcoords(); // set grid by x,y,z
	addtogrid();     // add node to the grid

    //unit_vec_posn = this->unitvec();

	//colour=ptheactin->newnodescolour;

	creation_iter_num = ptheactin->iteration_num;

    //nucleator_stuck_position = *this;

	harbinger = true;
    //unit_vec_correct = false;
    setunitvec();

    //previous_pos_in_nuc_frame = pos_in_nuc_frame;

	return true;       
}

bool nodes::polymerize(const vect& v)
{
    return polymerize(v.x, v.y, v.z);
}

int nodes::save_data(ofstream &ostr) 
{
    Colour dummycol;
    // save the nodes
    ostr << nodenum << " " 
	 << harbinger << " " 
	 << polymer << " " 
     << stucktonucleator << " "
     << testnode << " "
     << testsurface << " "
	 << links_broken << " "     
     << creation_iter_num << " "    
     << dist_from_surface << " "     
	 << x << " " << y << " " << z << " " 
     << delta << " "
     << delta_sum << " "
//     << unit_vec_posn << " "
//     << pos_in_nuc_frame << " "
//     << previous_pos_in_nuc_frame << " "
	 << linkforce_transverse << " " 
	 << linkforce_radial << " " 
	 << repforce_transverse << " " 
	 << repforce_radial << " "  
	 << nucleator_impacts << " " 
	 << nucleator_stuck_position << " " 
     << nucleator_link_force
	 << ":";
    
    // now the links
    ostr << (unsigned int) listoflinks.size() << ")";
    for(vector<links>::iterator l=listoflinks.begin(); l!=listoflinks.end(); ++l)
    {
	    l->save_data(ostr);
	    ostr << " ";
    }
    
    // done
    ostr << endl;
    
    return 0;
}

bool nodes::load_data(ifstream &istrm) 
{
    //unit_vec_correct = false;
    Colour dummycol;
    // read in from the stream to our private data
    char ch;

    ifstream::pos_type filepos = istrm.tellg();

    istrm >> nodenum 
	 >> harbinger 
	 >> polymer 
     >> stucktonucleator
     >> testnode
     >> testsurface
	 >> links_broken     
     >> creation_iter_num    
     >> dist_from_surface
     >> x >> y >> z 
     >> delta
     >> delta_sum
//     >> unit_vec_posn
//     >> pos_in_nuc_frame
//     >> previous_pos_in_nuc_frame
	 >> linkforce_transverse 
	 >> linkforce_radial 
	 >> repforce_transverse 
	 >> repforce_radial  
	 >> nucleator_impacts 
	 >> nucleator_stuck_position 
     >> nucleator_link_force
	 >> ch;

   // istrm >> nodenum 
	  //>> x >> y >> z 
	  //>> harbinger  
	  //>> polymer
   //   >> testnode
   //   >> testsurface
	  //>> delta 
	  //>> linkforce_transverse  
	  //>> linkforce_radial  
	  //>> repforce_transverse  
	  //>> repforce_radial                           
	  //>> links_broken
	  //>> nucleator_impacts 
	  //>> stucktonucleator 
	  //>> nucleator_stuck_position 
   //   >> nucleator_link_force 
	  //>> creation_iter_num >> ch;

    if (stucktonucleator && !STICK_TO_NUCLEATOR)
        stucktonucleator = false; // if we're continuing run and turning off nucleator attachments, we delete them here

	//if (nucleator_impacts>0.00001)
	//	cout << "Loaded nucleator impact " <<nucleator_impacts << endl;
    
    //cout << "Loaded node " << nodenum << endl;

    if (nodenum < 0)
    {
 

        cout << "Nodenum < 0 : " << nodenum << endl;
        cout << "Read: " << nodenum << " " 
	         << x << " " << y << " " << z << " " 
	         << harbinger << " " 
	         << polymer << " " 
             << testnode << " "
             << testsurface << " "
	         << delta << " " 
	         << linkforce_transverse << " " 
	         << linkforce_radial << " " 
	         << repforce_transverse << " " 
	         << repforce_radial << " " 
	         << links_broken << " " 
	         << nucleator_impacts << " "
	         << stucktonucleator << " "
	         << nucleator_stuck_position << " " 
             << nucleator_link_force << " "
	         << creation_iter_num << endl;
        return false;
    }

    // check we are ready to read links
    if(ch!=':' ){
	cout << endl << "error in checkpoint file at pos " << filepos << ", end of node ':' expected" 
	     << endl;


            if (!istrm.seekg(filepos, ios::beg))
                cout << "Unable to report text of line (seek failed)" << endl;
            else
            {

                char line1[16480];//,line2[16480],line3[16480];

                istrm.getline(line1, 16480, '\n');
                ifstream::pos_type pos1= istrm.tellg();
                //istrm.getline(line2, 16480, '\n');
                //ifstream::pos_type pos2= istrm.tellg();
                //istrm.getline(line3, 16480, '\n');
                //ifstream::pos_type pos3= istrm.tellg();

                cout << "Read the following line:" << endl
                     << pos1 << " " << line1 << endl;
                     //<< pos2 << " " << line2 << endl
                     //<< pos3 << " " << line3 << endl;
            
            }

            cout << "Interpreted as follows: " << nodenum << " " 
	         << x << " " << y << " " << z << " " 
	         << harbinger << " " 
	         << polymer << " " 
             << testnode << " "
             << testsurface << " "
	         << delta << " " 
	         << linkforce_transverse << " " 
	         << linkforce_radial << " " 
	         << repforce_transverse << " " 
	         << repforce_radial << " " 
	         << links_broken << " " 
	         << nucleator_impacts << " "
	         << stucktonucleator << " "
	         << nucleator_stuck_position << " " 
             << nucleator_link_force << " "
	         << creation_iter_num << endl;
	return false;
    }
    
    int linklistsize;
    istrm >> linklistsize >> ch;
    if( ch!=')' ){
	cout << "error in checkpoint file, xlinkdelays 'NN)' expected" 
	     << endl;

    setunitvec();

	return false;
    }

    // use resize(0) if we can, to save reallocating memory
#ifdef NODEGRIDTYPELIST
    listoflinks.clear();
#else
    listoflinks.resize(0);
#endif

    for(int i=0; i != linklistsize; ++i)
    {	 
        listoflinks.push_back(links(istrm));
    }

    // note we don't set pointer or build the grid here
    // because we need to be sure this is the node
    // stored by the actin.  The actin knows that.

    setunitvec();

//    previous_pos_in_nuc_frame=pos_in_nuc_frame; // not worth saving this?

    // only add to grid if doing normal run or a post-process that needs links:
    if (POST_VTK || (!REWRITESYMBREAK && !POST_PROCESS))
	    updategrid();

    return true;
}


void nodes::updategrid(void)
{
    // check if node has moved far enough since last update to warrent another update
    if ((fabs(posnoflastgridupdate.x - x) < gridscanjitter) &&
        (fabs(posnoflastgridupdate.y - y) < gridscanjitter) &&
        (fabs(posnoflastgridupdate.z - z) < gridscanjitter))
        return;  // no, then just skip the update

	short int oldgridx = gridx;			// store old grid pos'n
	short int oldgridy = gridy;
	short int oldgridz = gridz;

	setgridcoords();				// set gridx,y,z by x,y,z position

	if	((gridx == oldgridx) &&		// has the node moved gridpoints?
		 (gridy == oldgridy) && 
		 (gridz == oldgridz))
         return;  // no, then skip the update

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

 	short int gridtmpx, gridtmpy, gridtmpz;

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
	{
		addtogrid();
	}
			
	return;
}

void nodes::removefromgrid(void)
{
	// are we on the grid?
	if (nodegridptr==NULL) return;  // return if not on grid

    // find the node
    for (NODEGRIDTYPE <nodes*>::iterator i_node  = nodegridptr->begin();
		                                 i_node != nodegridptr->end() ;
							           ++i_node )
	{	 
		if (this == *i_node)
		{
            // remove node
			nodegridptr->erase(i_node);
			break;
		}
	}

	gridx = gridy = gridz = -1; // undo grid co-ords
    
    posnoflastgridupdate.zero();

    nodegridptr = NULL;

	return;
}

void nodes::addtogrid()
{
    nodegridptr = &NODEGRID(gridx,gridy,gridz);

    nodegridptr->push_back(this); // add to nodegrid

    posnoflastgridupdate = *this;

	return;
}

void nodes::setgridcoords(void)
{  

	gridx = (((short int)(x / GRIDRES)) + (GRIDSIZE/2) ); // find new grid pos'n
	gridy = (((short int)(y / GRIDRES)) + (GRIDSIZE/2) );
	gridz = (((short int)(z / GRIDRES)) + (GRIDSIZE/2) );

	return;
} 

void nodes::addlink(nodes& linkto, const double& dist)
{
    assert( &linkto != 0);

	listoflinks.push_back(links(linkto,dist));
	ptheactin->linksformed++;

	return;
}


void nodes::removelink(const nodes* linkednode)
{
	if (listoflinks.size()==0)
		return;

    // find the link to erase

	for (vector <links>::iterator i  = listoflinks.begin();
		                          i != listoflinks.end() ;
							    ++i )
	{	 
		if ( (i->linkednodeptr == linkednode) || ( i->linkednodeptr == 0) )
		{
			listoflinks.erase(i);
            
            links_broken++;
			ptheactin->linksbroken++;

			break;
		}
	}
}

int nodes::savelinks(ofstream * outputstream)
{
	if (listoflinks.empty())
		return 0;

	*outputstream << nodenum << " " << (unsigned int) listoflinks.size();

	for (vector <links>::iterator i=listoflinks.begin(); i!=listoflinks.end() ; ++i )
	{	 
		*outputstream << " " << i->linkednodeptr->nodenum;
	}

	*outputstream << endl;

	return 0;
}



    void nodes::setunitvec(void)
	{	
        //if (unit_vec_correct)
        //    return;
        //unit_vec_correct = true;  // prevent unnecessary recalculation if node not moved, etc.

//        previous_pos_in_nuc_frame = pos_in_nuc_frame; // old position
        //ptheactin->torque_rotate.rotate(previous_pos_in_nuc_frame);  // (don't need this!) need to rotate in case nucleator rotated

        pos_in_nuc_frame=*this;
        ptheactin->world_to_nuc_frame(pos_in_nuc_frame);

		if (NUCSHAPE == nucleator::sphere)
		{
            dist_from_surface = pos_in_nuc_frame.length();	 // not really dist_from_surface yet, need to subtract radius
            unit_vec_posn = pos_in_nuc_frame / dist_from_surface;  // set unit vector position
			//nearest_surface_point = unit_vec_posn
			dist_from_surface -= RADIUS;

		}
		else if (NUCSHAPE == nucleator::ellipsoid)
        {   // TODO: note the vector isn't right yet!                                                                    
            vect posonsphere(pos_in_nuc_frame.x, pos_in_nuc_frame.y, pos_in_nuc_frame.z / ELLIPSOID_STRETCHFACTOR);
            
            dist_from_surface = posonsphere.length();	 // not really dist_from_surface yet, need to subtract radius
            unit_vec_posn = posonsphere / dist_from_surface;  // normalise to set unit vector position
			
            
            //nearest_surface_point = unit_vec_posn
			dist_from_surface -= RADIUS;

		}
        else
		{	// capsule
			if (fabs(pos_in_nuc_frame.z) < CAPSULE_HALF_LINEAR)
			{  // on cylinder, no z component

				dist_from_surface = calcdist(pos_in_nuc_frame.x,pos_in_nuc_frame.y);   // not really dist_from_surface yet, need to subtract radius
				unit_vec_posn.x = pos_in_nuc_frame.x / dist_from_surface;
                unit_vec_posn.y = pos_in_nuc_frame.y / dist_from_surface;
                unit_vec_posn.z = 0.0;
				dist_from_surface -= RADIUS;
				//nearest_surface_point = unit_vec_posn + vect(0,0,z);
				onseg = true;

			}
			else
			{	// on ends

				onseg = false;

				if (pos_in_nuc_frame.z>0) // top
				{
					vect offsetvec = pos_in_nuc_frame;
					offsetvec.z -= CAPSULE_HALF_LINEAR;

                    dist_from_surface = offsetvec.length();	// not really dist_from_surface yet, need to subtract radius

                    unit_vec_posn = offsetvec / dist_from_surface;  // set unit vector pos

					//nearest_surface_point = unit_vec_posn + vect(0,0,CAPSULE_HALF_LINEAR);

					dist_from_surface -= RADIUS;

				}
				else
				{
					vect offsetvec = pos_in_nuc_frame;
					offsetvec.z += CAPSULE_HALF_LINEAR;

                    dist_from_surface = offsetvec.length();	 // not really dist_from_surface yet, need to subtract radius

                    unit_vec_posn = offsetvec / dist_from_surface;  // set unit vector pos

					//nearest_surface_point = unit_vec_posn + vect(0,0,CAPSULE_HALF_LINEAR);

					dist_from_surface -= RADIUS;

				}	
			}
        }
    }
