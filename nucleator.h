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

#ifndef nucleator_H
#define nucleator_H

#include "vect.h"
#include "rotationmatrix.h"
#include "Colour.h"

class nodes;
class actin;

class nucleator 
{
public:
	enum shape {
		sphere,
		capsule };
	nucleator(void);
	~nucleator(void);
	nucleator(shape set_geometry, actin *ptheactin);
	shape geometry;
	MYDOUBLE radius;
	MYDOUBLE segment;
	MYDOUBLE surf_area;
	MYDOUBLE movability;
	actin *ptheactin;

	vect position, deltanucposn;
	vect torque, centerofmass;
	vect momentofinertia;

	rotationmatrix nucleator_rotation;

	Colour colour;



	vector <vect> cagepoints;

	vector <MYDOUBLE> radial_rep_distrib_x;
	vector <MYDOUBLE> radial_rep_distrib_y;
	vector <MYDOUBLE> radial_rep_distrib_z;
	
	int nbdy_segs;
	int ncap_segs;

	vector <MYDOUBLE> fbar_cap_x;
	vector <MYDOUBLE> fbar_cap_y;
	vector <MYDOUBLE> fbar_cap_ang;
	vector <MYDOUBLE> fbar_bdy_x;
	vector <MYDOUBLE> fbar_bdy_y;

	int addnodes(void);
	int addnodessphere(void);
	int addnodescapsule(void);
	int definenucleatorgrid(void);
	int save(ofstream *outputstream) ;
	int savevrml(ofstream *outputstream) ;
	int saveradialsegments(ofstream *outputstream);
	int clearradialsegments();
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	void set_rep_bins();
//	int  get_rep_bin(MYDOUBLE angle);
	int  get_zbin(const MYDOUBLE x, const MYDOUBLE y);
	int  get_angbin(const MYDOUBLE x, const MYDOUBLE y);
//	MYDOUBLE get_rep_angle(MYDOUBLE x, MYDOUBLE y);
	bool is_sphere();
	bool is_capsule();

	bool iswithinnucleator(const MYDOUBLE& x, const MYDOUBLE& y, const MYDOUBLE& z);
	bool collision(nodes &node);//(MYDOUBLE &x, MYDOUBLE &y, MYDOUBLE &z);
	int n_force_segments();
	void definecagepoints(void);
};

#endif

