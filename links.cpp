/*

    comet - an actin-based motility simulator
    Copyright (C) 2009 Mark J Dayel

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

    When using the code in academic work please cite:

      Dayel MJ, Akin O, Landeryou M, Risca V, Mogilner A, et al. (2009) 
      In Silico Reconstitution of Actin-Based Symmetry Breaking and Motility. 
      PLoS Biol 7(9):e1000201. doi:10.1371/journal.pbio.1000201

    and include any modifications to the source code with the published work.

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
    last_link_length = 0.0;

    forcesum=0.0;
}

links::~links(void)
{
}

links::links(nodes& linknode, const double& linkdist)
{
	linkednodeptr = &linknode;
    linkednodenumber = linknode.nodenum;
    orig_dist = linkdist;
    last_link_length = linkdist;
    //lastbutone_link_length = linkdist;
    broken = false;
    //last_force_set = false;

    orig_distsqr = orig_dist*orig_dist;
	orig_dist_recip = 1/orig_dist;
	linkforcescalefactor = pow(orig_dist_recip,LINK_POWER_SCALE);
    forcesum=0.0;
}

links::links(ifstream &istr)
{
    load_data(istr);
}


int links::save_data(ofstream &ostr) 
{
    ostr << broken << " " 
		 << linkednodeptr->nodenum << " "
		 << orig_dist << " "
         //<< last_force << " "
         //<< last_force_set << " "
         << last_link_length << " "
         //<< lastbutone_link_length << " "
         << forcesum;
	
    return 0;
}

int links::load_data(ifstream &istr) 
{
    //double dummy;

    istr >> broken 
		 >> linkednodenumber
         >> orig_dist  
         //>> last_force 
         //>> last_force_set
         >> last_link_length
         //>> lastbutone_link_length
         >> forcesum;                  

   

	orig_dist_recip = 1/orig_dist;
	linkforcescalefactor = pow(orig_dist_recip,LINK_POWER_SCALE);
    orig_distsqr = orig_dist*orig_dist;
    broken = false;   // todo: are we catching broken links for removal before/after load?
    //last_force_set = true;
    // note cannot set the node pointer yet
    // linkednodeptr = &linknode;

    return 0;
}

