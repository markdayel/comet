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

links::links(void)
{
	orig_distsqr = 0.0;
	orig_dist = 0.0;
	breakcount = 0;
	breaklastiter = 0;
	broken = true;
	linkednodeptr = 0;
	linkednodenumber = -1;
}

links::~links(void)
{
}

links::links(nodes* linknodep, double dist)
{
	orig_distsqr = dist*dist;
	orig_dist = dist;
	orig_dist_recip = 1/dist;
	broken = false;
	breakcount = 0;
	breaklastiter = 0;
	linkednodeptr = linknodep;
}

double links::getlinkforces(const double
                            
                            & dist)
{  // return force (nominally in pN)
	double force;//=0.0;
    double strain = dist / orig_dist;
	double stress_over_breakage;
	// is link loose or taut?

    static const bool local_USE_BREAKAGE_STRAIN = USE_BREAKAGE_STRAIN;
    // calculate forces:

    force = - LINK_FORCE * (dist - orig_dist) * orig_dist_recip;

	if (dist > (orig_dist * LINK_TAUT_RATIO))
	{  // filaments taut:  go to high strain regime

		force += - ( LINK_TAUT_FORCE  * 
            (dist - (orig_dist * LINK_TAUT_RATIO)) * orig_dist_recip );
    }

    // decide if link broken:

    if (local_USE_BREAKAGE_STRAIN)
    {
	    if (strain > LINK_BREAKAGE_STRAIN)
		{
			breakcount++;

			if ( breakcount * P_LINK_BREAK_IF_OVER * DELTA_T * RAND_MAX > rand() )
			{
				broken = true;
				force = 0;  
			}

		}
		else
		{
			breakcount = 0;
		}
    }
    else
    { // just using force
		if (-force > LINK_BREAKAGE_FORCE)
		{
			stress_over_breakage = -force / LINK_BREAKAGE_FORCE;
			breakcount++;

			if ( (breakcount*P_LINK_BREAK_IF_OVER*DELTA_T*stress_over_breakage) * RAND_MAX > 
					rand() )
			//if ((++breakcount>MAX_LINK_BREAKCOUNT) && breaklastiter)
			{
				broken = true;
                force = 0;
			}

		}
		else
		{
			breakcount = 0;
		}
    }

	return force;
}

int links::save_data(ofstream &ostr) 
{
    ostr << broken << "," 
	 << breakcount << "," 
	 << breaklastiter << "," 
	 << orig_dist << ","
	 << linkednodeptr->nodenum;

    //if(linkednodeptr != NULL)
    //	ostr << linkednodeptr->nodenum;
    //  else
    //ostr << -1;
    return 0;
}

int links::load_data(ifstream &istr) 
{
    char ch;
    istr >> broken >> ch
	 >> breakcount >> ch 
	 >> breaklastiter >> ch 
	 >> orig_dist >> ch 
	 >> linkednodenumber;
    return 0;
}

void links::reset_link(void)
{
    breakcount = 0;
    breaklastiter = false;
    broken = false;
}
