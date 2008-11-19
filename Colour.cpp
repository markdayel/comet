/*
    comet - a actin-based motility simulator
    Copyright (C) 2005 Mark J Dayel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
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
