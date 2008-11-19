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

	vect position, last_delta_position, deltanucposn_sum, deltanucposn;
	vect torque, centerofmass;
	vect momentofinertia;

    rotationmatrix last_torque_rotate;


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

