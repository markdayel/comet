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

Colour::Colour(double val)
{
	setcol(val);
}

Colour::~Colour(void)
{
}

int Colour::setcol(double magnitude)
{  // colour scale 0--1
	
	mag=magnitude;

	r =  (double) mymin(  mymax( ((double)4*    (mag-(double)0.25) ) , (double)0 ), (double)1 );
	b =  (double) mymin(  mymax( ((double)4*    ((double)0.75-mag) ) , (double)0 ), (double)1 );
	g =  (double) mymin(  mymax( ((double)4*fabs(mag-(double)0.5)-1) , (double)0 ), (double)1 );

	R = (unsigned char) (r*255);
	G = (unsigned char) (G*255);
	B = (unsigned char) (B*255);

	return 0;
}
