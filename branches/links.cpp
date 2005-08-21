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


//
//
//int links::savedata(ofstream *outputstream) 
//{
//	*outputstream << broken << "," << breakcount << "," << breaklastiter
//		<< "," << orig_dist << "," << linkednodeptr->nodenum;
//
//	return 0;
//}
//
//int links::loaddata(ifstream *inputstream) 
//{
//	char delim;
//
//	*inputstream >> broken >> delim >> breakcount >> delim >> breaklastiter 
//		>> delim >> orig_dist >> delim >>  linkednodeptr->nodenum;
//
//	orig_distsqr = orig_dist*orig_dist;
//	orig_dist_recip = 1/orig_dist;
//
//	return 0;
//}

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
