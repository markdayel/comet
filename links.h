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

#ifndef links_H
#define links_H

class links
{
public:
	links(void);
	links(nodes& linknodep, const double &dist);
    links(ifstream &istr);
	 ~links(void);
	//int savedata(ofstream *outputstream);
	//int loaddata(ifstream *inputstream);
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	nodes* linkednodeptr;
	int linkednodenumber;
	double orig_distsqr;
	double orig_dist;
	double orig_dist_recip;
	double linkforcescalefactor;
    double last_link_length;//, lastbutone_link_length;

    double forcesum;

    //double last_force;

	//double last_dist;
	//double last_but_one_dist;

	//int breakcount;
	bool broken;//, last_force_set;	
	//bool breaklastiter;
	//double theta;
	//double phi;

    inline void clearstats()
	{
        forcesum=0.0;
    }

    inline bool getlinkforces(const double & dist, double &force)
    {   /// calculate link forces
        /// break link if above LINK_BREAKAGE_FORCE
        /// returns false if link broken

        // force is on the node in the direction of the link,
        // i.e. as the link is stretched longer, the force becomes more negative

        // calculate forces:

        force = -LINK_FORCE * (dist - orig_dist) * orig_dist_recip * linkforcescalefactor; 

        if (-force > LINK_BREAKAGE_FORCE)
        {   
            broken = true;
            return false;
        }

        force -= NODE_DIST_TO_FORCE * DASHPOT_IMPEDANCE * (dist - last_link_length);

        //force -= NODE_DIST_TO_FORCE * DASHPOT_IMPEDANCE * ( (dist - last_link_length) - 
        //                     ((dist - last_link_length) - (last_link_length - lastbutone_link_length)) / 2 );

        //lastbutone_link_length = last_link_length;
        last_link_length = dist;

        return true;
    }
	
};

#endif
