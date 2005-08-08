#ifndef segments_H 
#define segments_H

//#include "stdafx.h"
//#include "nucleator.h"

typedef vector<MYDOUBLE> Dbl1d;
typedef vector<Dbl1d> Dbl2d;
typedef vector<Dbl2d> Dbl3d;

class nodes;
class nucleator;
class actin;

class segments
{
public:
	//segments(nucleator::shape set_geometry);
	segments(void);
	~segments(void);

	void setupsegments(nucleator *pn, actin * pactin);

	void drawoutline(ostream& out, const int& axis) const;
	void drawsurfaceimpacts(ostream& out, const int& axis, const MYDOUBLE scale) const;

	void addnode(const nodes& node);
	void addsurfaceimpact(const nodes& node, const MYDOUBLE& mag);

	void getsegmentnum(const nodes& node, int& xseg, int& yseg, int& zseg) const;
	int getcapsuleseg(const MYDOUBLE & x, const MYDOUBLE & y) const;

	void clearsurfaceimpacts(void);
	void clearnodes(void);

	nucleator * p_nuc;
	actin * p_actin;

	int num_segs;					// number of segments all the way round
	int num_cap_segs;				// num segs on *one* cap
	int num_straight_segs;			// num segs on *one* side
	int num_dist_segs;				// number of segs going outwards from nuc
	

	MYDOUBLE curved_length;			// length of *one* cap
	MYDOUBLE straight_length;		// length of *one* side

	MYDOUBLE straight_seg_len;		// distance between transverse segs on straight part of nuc
	MYDOUBLE cap_seg_len;			// distance between transverse segs on caps of nuc
	MYDOUBLE dist_step;				// distance between outward segs


	vector <MYDOUBLE> linestartx;	// x position of force line origin
	vector <MYDOUBLE> linestarty;	// y position of force line origin
	vector <MYDOUBLE> lineunitvecx;	// unit vector of force line x 
	vector <MYDOUBLE> lineunitvecy;	// unit vector of force line y

									// for axes, 0=x, 1=y, 2=z

	Dbl2d							// 2d vector for surface [axis][segnum]
		surfaceimpacts;

	Dbl3d 							// 3d vectors for the rest [axis][segnum][dist]
		numnodes,
		rep_radial,
		rep_transverse,
		link_radial,
		link_transverse,
		disp_radial,
		disp_transverse;

};

#endif
