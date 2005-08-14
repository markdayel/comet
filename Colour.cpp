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

int Colour::setcol(const double & magnitude)
{  // colour scale 0--1
	
	mag = magnitude;

	// hot colormap
	//r =  mymin(  mymax( 8*mag / 3    , 0 ), 1 );
	//g =  mymin(  mymax( 8*mag / 3 - 1, 0 ), 1 );
	//b =  mymin(  mymax( 4*mag     - 3, 0 ), 1 );

	// rainbow colormap, 0 = black
	r = mymin( mymax( -2 * fabs( mag - 1   ) + 1.3  ,0) ,1);
	g = mymin( mymax( -2 * fabs( mag - 0.55) + 1.05 ,0) ,1);
	b = mymin( mymax( -5 * fabs( mag - 0.25) + 1.2  ,0) ,1);

	// hls, sort of
	//r =   mymin(  mymax( (4*    (mag-0.25) ) , 0 ), 1 );
	//b =   mymin(  mymax( (4*    (0.75-mag) ) , 0 ), 1 );
	//g =   mymin(  mymax( (4*fabs(mag-0.5)-1) , 0 ), 1 );

	R = (unsigned char) (r*255);
	G = (unsigned char) (g*255);
	B = (unsigned char) (b*255);

	return 0;
}
