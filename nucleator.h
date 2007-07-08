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
		capsule = 1,
        ellipsoid = 2};
	nucleator(void);
	~nucleator(void);
	//shape geometry;

	double surf_area;
	double inertia;

	segments segs;

	vect position, deltanucposn;
	vect torque, centerofmass;
	vect momentofinertia;


	Colour colour;

	vector <vect> cagepoints;


	//inline bool is_sphere()
	//{
	//	return (geometry == sphere);
	//}

    void move_nuc(const vect& origin_of_movement, const vect& tomove);

	int addnodes(void) const;
	int addnodessphere(void) const;
	int addnodescapsule(void) const;
    int addnodesellipsoid(void) const;
    double polyfeedbackprob(const double& x, const double& y, const double& z) const;

	int savevrml(ofstream *outputstream);
	int save_data(ofstream &ostr);
	int load_data(ifstream &istr);
	//bool is_sphere();
	//bool is_capsule();

	bool collision(nodes &node);

    vect eject_sphere(const vect& inputpoint) const;
    vect eject_ellipsoid(const vect& inputpoint) const;   
    vect eject_capsule(const vect& inputpoint) const;
    vect eject_point(const vect& inputpoint) const;
	
    void definecagepoints(void);
};

#endif

