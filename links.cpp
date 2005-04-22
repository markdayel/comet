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
}

links::~links(void)
{
}

links::links(nodes* linknodep, MYDOUBLE dist)
{
	orig_distsqr = dist*dist;
	orig_dist = dist;
	broken = false;
	breakcount = 0;
	breaklastiter = 0;
	linkednodeptr = linknodep;
}


MYDOUBLE links::getlinkforces(MYDOUBLE dist)
{  // return force (nominally in pN)
	MYDOUBLE force=0.0;
	MYDOUBLE stress_over_breakage;
	// is link loose or taught?

	if (dist > (orig_dist*LINK_TAUGHT_RATIO))
	{  // filaments taught:  go to high strain regime

		force =		- ( LINK_FORCE * (dist - orig_dist) +
				 LINK_TAUGHT_FORCE * (dist - (orig_dist*LINK_TAUGHT_RATIO)))
					       / orig_dist;

		if ((-force) > LINK_BREAKAGE_FORCE)
		{
			stress_over_breakage = (-force)/LINK_BREAKAGE_FORCE;
			breakcount++;

			if ( (breakcount*P_LINK_BREAK_IF_OVER*DELTA_T*stress_over_breakage) > 
				      ( ((MYDOUBLE) rand()) / (MYDOUBLE)(RAND_MAX) ) )
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
	else
	{  // loose: entropic spring

		force = - LINK_FORCE * (dist - orig_dist) / orig_dist;

	}

	return force;
}

int links::savedata(ofstream *outputstream) 
{
	*outputstream << broken << "," << breakcount << "," << breaklastiter
		<< "," << orig_dist << "," << linkednodeptr->nodenum;

	return 0;
}

int links::loaddata(ifstream *inputstream) 
{
	char delim;

	*inputstream >> broken >> delim >> breakcount >> delim >> breaklastiter 
		>> delim >> orig_dist >> delim >>  linkednodeptr->nodenum;

	orig_distsqr = orig_dist*orig_dist;

	return 0;
}
