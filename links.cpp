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
#include "links.h"

//#define LINK_POWER_SCALE 1.5

links::links(void)
{
	orig_distsqr = 0.0;
	orig_dist = 0.0;
	broken = true;
	linkednodeptr = 0;
	linkednodenumber = -1;
}

links::~links(void)
{
}

links::links(nodes& linknode, const double& dist)
{
	orig_distsqr = dist*dist;
	orig_dist = dist;
	orig_dist_recip = 1/orig_dist;
	strengthscalefactor = pow(orig_dist_recip,LINK_POWER_SCALE);
	broken = false;
	linkednodeptr = &linknode;
}

double links::getlinkforces(const double & dist)
{  // return force (nominally in pN)

    // calculate forces:

    double force = - LINK_FORCE * (dist - orig_dist) * orig_dist_recip * strengthscalefactor;

    // decide if link broken:

    if (USE_BREAKAGE_STRAIN)
    {	
		// since strainlimit = distlimit / orig_dist,
		// then distlimit = strainlimit * orig_dist
	    if (dist > LINK_BREAKAGE_STRAIN * orig_dist)
		{
			broken = true;
			force = 0;  
		}
    }
    else
    { // just using force
		if (-force > LINK_BREAKAGE_FORCE * strengthscalefactor)
		{   
			broken = true;
            force = 0;
		}

    }

	return force;
}


int links::save_data(ofstream &ostr) 
{
    ostr << broken << " " 
		 << orig_dist << " "
		 << linkednodeptr->nodenum;
	
    return 0;
}

int links::load_data(ifstream &istr) 
{
    istr >> broken 
		 >> orig_dist  
		 >> linkednodenumber;

	orig_dist_recip = 1/orig_dist;
	strengthscalefactor = pow(orig_dist_recip,LINK_POWER_SCALE);

    return 0;
}

void links::reset_link(void)
{
    broken = false;
}
