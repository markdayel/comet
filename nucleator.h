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
#include "segments.h"

class nodes;
class actin;

class nucleator 
{
public:
	enum shape {
		sphere = 0,
		capsule = 1};
	nucleator(void);
	~nucleator(void);
	nucleator(shape set_geometry);//, actin *ptheactin);
	shape geometry;

	double surf_area;
	double movability;
	//actin *ptheactin;

	segments segs;

	vect position, deltanucposn;
	vect torque, centerofmass;
	vect momentofinertia;

	//rotationmatrix nucleator_rotation;

	Colour colour;

	vector <vect> cagepoints;


	inline bool is_sphere()
	{
		return (geometry == sphere);
	}

    void move_nuc(vect& origin_of_movement, vect& tomove);

	int addnodes(void);
	int addnodessphere(void);
	int addnodescapsule(void);
	//int definenucleatorgrid(void);
	//int save(ofstream *outputstream) ;
	int savevrml(ofstream *outputstream) ;
	//int saveradialsegments(ofstream *outputstream);
	//int clearradialsegments();
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	//void set_rep_bins();
	//int  get_rep_bin(double angle);
	//int  get_zbin(const double x, const double y);
	//int  get_angbin(const double x, const double y);
	//double get_rep_angle(double x, double y);
	//bool is_sphere();
	//bool is_capsule();

	//bool iswithinnucleator(const double& x, const double& y, const double& z);
	bool collision(nodes &node);//(double &x, double &y, double &z);
	//int n_force_segments();
	void definecagepoints(void);
};

#endif

