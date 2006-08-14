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

#ifndef Colour_H
#define Colour_H

class Colour
{
public:
	Colour(void);
	Colour(double val);
	~Colour(void);
	double r,g,b;
	unsigned char R,G,B;
	double mag;
	int setcol(const double & mag);
    void setwhite()
    {
        r=g=b=1.0;
        R=G=B=255;
    };
};

#endif

