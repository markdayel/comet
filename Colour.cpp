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

Colour::Colour(const double val)
{
	setcol(val);
}

Colour::Colour(const double rr, const double gg, const double bb)
{
	r=rr;
    g=gg;
    b=bb;
}

Colour::~Colour(void)
{
}

int Colour::setcol(const double & magnitude)
{  // colour scale 0--1
	
	mag = mymax(0, mymin( 1, magnitude )); // make sure it's between 0 and 1


	// hot colormap
	//r =  mymin(  mymax( 8*mag / 3    , 0 ), 1 );
	//g =  mymin(  mymax( 8*mag / 3 - 1, 0 ), 1 );
	//b =  mymin(  mymax( 4*mag     - 3, 0 ), 1 );

	// rainbow colormap, 0 = black
//	r = mymin( mymax( -2 * fabs( mag - 1   ) + 1.3  ,0) ,1);
//	g = mymin( mymax( -2 * fabs( mag - 0.55) + 1.05 ,0) ,1);
	////b = mymin( mymax( -5 * fabs( mag - 0.25) + 1.2  ,0) ,1);
//    b = mymin( mymax( -2 * fabs( mag - 0.25) + 1.2  ,0) ,1);

	// hls, sort of
	//r =   mymin(  mymax( (4*    (mag-0.25) ) , 0 ), 1 );
	//b =   mymin(  mymax( (4*    (0.75-mag) ) , 0 ), 1 );
	//g =   mymin(  mymax( (4*fabs(mag-0.5)-1) , 0 ), 1 );

    //r = mymin( mymax( -2 * abs( mag - 1   ) + 1.3  ,0) ,1);
    //g = mymin( mymax( -2 * abs( mag - 0.4) + 1.15 ,0) ,1);
    //b = mymin( mymax( -2 * abs( mag - 0.05) + 1.2  ,0) ,1);

    

//r = mymin( mymax( -1.8* fabs( mag - 1   ) + 1.3  ,0) ,1);
//g = mymin( mymax( -2.6 * fabs( mag - 0.58) + 1.05 ,0) ,1);

//b = mymin( mymax( -3 * fabs( mag - 0.25) + 1.2  ,0) ,1);


    r = mymin( mymax( -1.8* fabs( mag - 1   ) + 1.3  ,0) ,1);
if (mag < 0.58) 
    g = mymin( mymax( -2.0 * fabs( mag - 0.58) + 1.05 ,0) ,1);
else
    g = mymin( mymax( -2.6 * fabs( mag - 0.58) + 1.05 ,0) ,1);

if (mag < 0.25) 
    b = mymin( mymax( -2.5 * fabs( mag - 0.25) + 1.2  ,0) ,1);
else
    b = mymin( mymax( -3 * fabs( mag - 0.25) + 1.2  ,0) ,1);


// gamma adjustment --- note this is applied both to the data and the scale bar, so we don't misrepresent anything
    r = pow( r , 1 / COLOUR_GAMMA);
    g = pow( g , 1 / COLOUR_GAMMA);
    b = pow( b , 1 / COLOUR_GAMMA);

	R = (unsigned char) (r*255);
	G = (unsigned char) (g*255);
	B = (unsigned char) (b*255);

	return 0;
}
