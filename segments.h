/*
    comet - an actin-based motility simulator
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

    When using the code in academic work please cite:

      Dayel MJ, Akin O, Landeryou M, Risca V, Mogilner A, et al. (2009) 
      In Silico Reconstitution of Actin-Based Symmetry Breaking and Motility. 
      PLoS Biol 7(9):e1000201. doi:10.1371/journal.pbio.1000201

    and include any modifications to the source code with the published work.
*/

#ifndef segments_H 
#define segments_H



//#include "stdafx.h"
//#include "nucleator.h"

typedef vector<double> Dbl1d;
typedef vector<Dbl1d> Dbl2d;
typedef vector<Dbl2d> Dbl3d;

typedef vector<int> Int1d;
typedef vector<Int1d> Int2d;
typedef vector<Int2d> Int3d;
typedef vector<Int3d> Int4d;
typedef vector<Int4d> Int5d;

class nodes;
class nucleator;
class actin;

#include "stdafx.h"

class segments
{
public:
	//segments(nucleator::shape set_geometry);
	segments(void);
	~segments(void);

	void setupsegments(nucleator *pn, actin * pactin);

	void drawoutline(ostream& out, const projection &axis) const;
	int drawsurfaceimpacts(ostream& out, const projection &axis,const rotationmatrix& axisrotation, const double scale) const;

    int drawline(ostream& drawcmd,
                    const double& startx, const double& starty, 
                    double forcevecx, double forcevecy,
                    const double& maxlen, const double& scalefactor) const;

	void addnode(const nodes& node);
	void addallnodes();
    void calcSD(const Dbl3d& data, Dbl2d& SD);
	//void addsurfaceimpact(nodes& node, const double& mag);

	void getsegmentnum(const vect& node, int& xseg, int& yseg, int& zseg) const;
    void getsegmentnum(const vect& node, double& xseg, double& yseg, double& zseg) const;
    double getsegmentnum(const vect& node, const projection &axis) const;

	void getsegmentdist(const nodes& node,int& xdist, int& ydist, int& zdist) const;

	double getcapsuleseg(const double & x, const double & y) const;
	void getsegmentposition(double& x, double& y, double& z, const int & seg,
								  const int & dist, const int & axis) const;


	int  dist_to_seg(const double & dist) const;

	void savereport(const int& filenum) const;
	void saveradialreport(const int& filenum) const;
    void saveSDreport(const int& filenum) const;
    void saveradialaxisreport(const int& filenum, const int axis) const;

	void write_bins_bitmap(Dbl2d &imageR, Dbl2d &imageG, Dbl2d &imageB,
					   const Dbl3d & var, const double& scale, const projection &axis);

    void write_colourmap_bitmap(Dbl2d &imageR, Dbl2d &imageG, Dbl2d &imageB); //, double scalefactor);

	void setbitmapcoords();

    bool load_scalefactors(void);
    void save_scalefactors(void);

	//void clearsurfaceimpacts(void);
	void clearbins(void);

	nucleator *  p_nuc;
	actin *  ptheactin;

	int num_segs;					// number of segments all the way round
	int num_cap_segs;				// num segs on *one* cap
	int num_straight_segs;			// num segs on *one* side
	int num_dist_segs;				// number of segs going outwards from nuc
	
	double curved_seg_area;
	double straight_seg_area;

	double curved_length;			// length of *one* cap
	double straight_length;			// length of *one* side

	double straight_seg_len;		// distance between transverse segs on straight part of nuc
	double cap_seg_len;				// distance between transverse segs on caps of nuc
	double dist_step;				// distance between outward segs

	int centerx, centery;			// co-ords of center of nucleator drawing


	//vector <double> linestartx;	// x position of force line origin
	//vector <double> linestarty;	// y position of force line origin
	//vector <double> lineunitvecx;	// unit vector of force line x 
	//vector <double> lineunitvecy;	// unit vector of force line y
	
	// for position vectors, first index is axis num

	Dbl2d linestartx;	// x position of force line origin
	Dbl2d linestarty;	// y position of force line origin
	Dbl2d lineunitvecx;	// unit vector of force line x 
	Dbl2d lineunitvecy;	// unit vector of force line y

									// for axes, 0=x, 1=y, 2=z

	Dbl3d 							// 3d vectors for the rest [axis][segnum][dist]
		numnodes,
		rep_radial,
		rep_transverse,
		link_radial,
		link_transverse, 
		links_broken;

    Dbl2d							// 2d vector for surface [axis][segnum]
		surfaceimpacts;

	Dbl3d
		surfacestuckforce;  // [axis][segment][direction] where direction: 0=radial, 1=transverse

	Dbl2d 							// 2d vectors for SD [axis][dist]
		numnodes_SD,
		rep_radial_SD,
		rep_transverse_SD,
		link_radial_SD,
		link_transverse_SD, 
		links_broken_SD;

    double 	numnodes_scalefactor,
		    rep_radial_scalefactor,
		    rep_transverse_scalefactor,
		    link_radial_scalefactor,
		    link_transverse_scalefactor, 
		    links_broken_scalefactor,
            surfaceimpacts_scalefactor;


	double radialdist;
	int num_radial_bins;

	int bins_bitmap_width, bins_bitmap_height;

	Dbl1d						// radial average
		radial_numnodes,				
		radial_rep_radial,
		radial_rep_transverse,
		radial_link_radial,
		radial_link_transverse,
		radial_links_broken;

    Dbl1d 
        radial_repenergy_radial,
        radial_repenergy_transverse, 
        radial_linkenergy_radial, 
        radial_linkenergy_transverse, 
        radial_linksenergy_broken;

	//Int4d segment_pixel;

    void set_scale_factors(void);
};

#endif
