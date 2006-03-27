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
	//int savedata(ofstream *outputstream);
	//int loaddata(ifstream *inputstream);
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	nodes* linkednodeptr;
	int linkednodenumber;
	double orig_distsqr;
	double orig_dist;
	double orig_dist_recip;

	//double last_dist;
	//double last_but_one_dist;

	int breakcount;
	bool broken;	
	bool breaklastiter;
	//double theta;
	//double phi;

	double getlinkforces(const double& dist);
	
    void reset_link(void);
};

#endif
