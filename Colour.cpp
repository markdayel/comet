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
#include "Colour.h"

Colour::Colour(void)
{
	setcol(0.0);
}

Colour::~Colour(void)
{
}

int Colour::setcol(MYDOUBLE magnitude)
{  // colour scale 0--1
	
	mag=magnitude;

	r =  (MYDOUBLE) mymin(  mymax( (4*    (mag-(MYDOUBLE)0.25) ) , 0 ), 1 );
	b =  (MYDOUBLE) mymin(  mymax( (4*    ((MYDOUBLE)0.75-mag) ) , 0 ), 1 );
	g =  (MYDOUBLE) mymin(  mymax( (4*fabs(mag-(MYDOUBLE)0.5)-1) , 0 ), 1 );

	return 0;
}
