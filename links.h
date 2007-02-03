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

    double forcesum;

    double last_force;

	//double last_dist;
	//double last_but_one_dist;

	//int breakcount;
	bool broken, last_force_set;	
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


        // calculate forces:

        force = - LINK_FORCE * (dist - orig_dist) * orig_dist_recip * linkforcescalefactor;

        if (force < 0)
        {
            // decide if link broken:

            if (USE_BREAKAGE_VISCOSITY)
            {	
		        // since strainlimit = distlimit / orig_dist,
		        // then distlimit = strainlimit * orig_dist
    	        
                //double p_break = DELTA_T * 1 * exp( (-LINK_BREAKAGE_FORCE * 12) / (force * force));

                if (-force > LINK_BREAKAGE_FORCE)
		        {
			        broken = true;
			        return false;  
		        }

                if (last_force_set)
                {
                    if ( fabs(last_force - force) > BREAKAGE_VISCOSITY_THRESHOLD * DELTA_T )
                    {
                        broken = true;
                        return false;
                    }
                }

            }
            else
            { // just using force
		        if (-force > LINK_BREAKAGE_FORCE)
		        {   
			        broken = true;
                    return false;
		        }

            }

        }

        last_force = force;
        last_force_set = true;

	    return true;
    }
	
};

#endif
