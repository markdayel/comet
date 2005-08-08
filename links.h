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
	links(nodes* linknodep, double dist);
	 ~links(void);
	int savedata(ofstream *outputstream);
	int loaddata(ifstream *inputstream);
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	nodes* linkednodeptr;
	int linkednodenumber;
	double orig_distsqr;
	double orig_dist;
	double orig_dist_recip;

	bool broken;
	int breakcount;
	bool breaklastiter;
	//inline double getlinkforces(const double& distsq);
	double theta;
	double phi;

	inline const double getlinkforces(const double& dist) 
	{  // return force (nominally in pN)
		double force;//=0.0;
		double stress_over_breakage;
		// is link loose or taut?

		if (dist > (orig_dist*LINK_TAUT_RATIO))
		{  // filaments taut:  go to high strain regime

			force =		- ( LINK_FORCE * (dist - orig_dist) +
					LINK_TAUT_FORCE * (dist - (orig_dist*LINK_TAUT_RATIO)))
							* orig_dist_recip;

			if ((-force) > LINK_BREAKAGE_FORCE)
			{
				stress_over_breakage = (-force)/LINK_BREAKAGE_FORCE;
				breakcount++;

				if ( (breakcount*P_LINK_BREAK_IF_OVER*DELTA_T*stress_over_breakage) > 
						( ((double) rand()) / (double)(RAND_MAX) ) )
				//if ((++breakcount>MAX_LINK_BREAKCOUNT) && breaklastiter)
				{
					broken = true;
					force = 0;  
				}

			}
			else
			{
				breakcount = 0;
			}

		}
		else
		{  // loose: entropic spring

			force = - LINK_FORCE * (dist - orig_dist) * orig_dist_recip;

		}

		return force;
	}
};

#endif
