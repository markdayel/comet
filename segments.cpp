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

#include "stdafx.h"
#include "nucleator.h"
#include "segments.h"


//segments::segments(nucleator::shape set_geometry) : geometry (set_geometry)
//{	// the constructor is really setupsegments
//	// it can't be done here because we need the nucleator pointer
//  // and I can't get it to work in the constructor (not possible?)
//}

segments::segments(void)
{
}

segments::~segments(void)
{
}

void segments::setupsegments(nucleator *pnuc, actin * pactin)
{

	p_nuc = pnuc;
	ptheactin = pactin;
	
	num_dist_segs = 20;		// normal bins radial distance
	dist_step = RADIUS * 0.1;

	radialdist = RADIUS * 0.1;		// for radial-only averaging
	num_radial_bins = 20;

	centerx = BMP_WIDTH  / 7;  // specifies center of nucleator on forces and seg maps
	centery = BMP_HEIGHT / 2;

	bins_bitmap_width  = centerx * 2;			// width and height of pixels to plot
	bins_bitmap_height = (int) ((double) BMP_HEIGHT * 2/3);

	// calc segment lengths etc.

	curved_length = PI * RADIUS;	// one cap
	straight_length = 2 * CAPSULE_HALF_LINEAR;		// one side

	num_cap_segs = (int) (RADIAL_SEGMENTS/2);
	cap_seg_len = curved_length / (double) num_cap_segs;

	// => will need to scale straight numbers by (cap_seg_len/straight_seg_len) before readout

	// add up segments

	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		num_straight_segs = 0;
		straight_length = 0;

		num_segs = num_cap_segs * 2;
	}
	else
	{	// capsule

		num_straight_segs = (int)(straight_length / cap_seg_len);
		straight_seg_len = straight_length / num_straight_segs;

		num_segs = num_cap_segs * 2 + num_straight_segs * 2;
	}

	// allocate and zero the 3d vectors

	       numnodes.resize(3);			
	     rep_radial.resize(3);
	 rep_transverse.resize(3);
	    link_radial.resize(3);
	link_transverse.resize(3);
	   links_broken.resize(3);
  surfacestuckforce.resize(3);

	for (int axis = 0; axis != 3; ++axis)
	{
		       numnodes[axis].resize(num_segs);			
		     rep_radial[axis].resize(num_segs);
		 rep_transverse[axis].resize(num_segs);
		    link_radial[axis].resize(num_segs);
		link_transverse[axis].resize(num_segs);
		   links_broken[axis].resize(num_segs);
	  surfacestuckforce[axis].resize(num_segs);

        for(int i = 0; i != num_segs; ++i)
		{
			       numnodes[axis][i].resize(num_dist_segs);
			     rep_radial[axis][i].resize(num_dist_segs);
			 rep_transverse[axis][i].resize(num_dist_segs);
			    link_radial[axis][i].resize(num_dist_segs);
			link_transverse[axis][i].resize(num_dist_segs);
			   links_broken[axis][i].resize(num_dist_segs);
		  surfacestuckforce[axis][i].resize(2);  // this is the direction, 0=radial, 1=transverse

		}
	}

	radial_numnodes.resize(num_radial_bins);				
	radial_rep_radial.resize(num_radial_bins);	
	radial_rep_transverse.resize(num_radial_bins);	
	radial_link_radial.resize(num_radial_bins);	
	radial_link_transverse.resize(num_radial_bins);	
	radial_links_broken.resize(num_radial_bins);


    radial_repenergy_radial.resize(num_radial_bins);
    radial_repenergy_transverse.resize(num_radial_bins); 
    radial_linkenergy_radial.resize(num_radial_bins); 
    radial_linkenergy_transverse.resize(num_radial_bins); 
    radial_linksenergy_broken.resize(num_radial_bins);


	// allocate surfaceimpacts vector

	surfaceimpacts.resize(3);
	for (int axis = 0; axis < 3; ++axis)
	{
		surfaceimpacts[axis].resize(num_segs);
	}


    // allocate stddev vectors

    numnodes_SD.resize(3);
	rep_radial_SD.resize(3);
	rep_transverse_SD.resize(3);
	link_radial_SD.resize(3);
	link_transverse_SD.resize(3);
	links_broken_SD.resize(3);

    for (int axis = 0; axis < 3; ++axis)
	{
        numnodes_SD[axis].resize(num_dist_segs);
		rep_radial_SD[axis].resize(num_dist_segs);
		rep_transverse_SD[axis].resize(num_dist_segs);
		link_radial_SD[axis].resize(num_dist_segs);
		link_transverse_SD[axis].resize(num_dist_segs);
		links_broken_SD[axis].resize(num_dist_segs);
    }

    // clear these vectors

	clearbins();


	// define force lines co-ords

	linestartx.resize(3);
	linestarty.resize(3);

	lineunitvecx.resize(3);
	lineunitvecy.resize(3);

	for (int axis = 0; axis < 3; ++axis)
	{
		linestartx[axis].resize(num_segs);
		linestarty[axis].resize(num_segs);

		lineunitvecx[axis].resize(num_segs);
		lineunitvecy[axis].resize(num_segs);
	}


	double startx,starty, theta;
	double dist;

	for (int axis = 0; axis != 3; ++axis)
	{
		// top cap

		for(int i = 0; i != num_cap_segs; ++i)
		{
			//       offset     segnum        angle of 1 seg
			theta = (0.5 + (double) i) * (PI / (double) num_cap_segs);

			startx = - cos (theta) * RADIUS;  // unit components
			starty =   sin (theta) * RADIUS;

			linestartx[axis][i] = startx;	// line start
	        
			if ((p_nuc->geometry == nucleator::sphere) || (axis == 2))
			{	// sphere or capsule z
				linestarty[axis][i] = starty;
			}
			else
			{
				linestarty[axis][i] = starty + CAPSULE_HALF_LINEAR;
			}

			lineunitvecx[axis][i] = - startx;
			lineunitvecy[axis][i] = - starty;
		}

		// bottom cap

		for(int i = num_cap_segs; i != (2 *num_cap_segs); ++i)
		{
			//       offset     segnum                      angle of 1 seg
			theta = (0.5 + (double) i) * (PI / (double) num_cap_segs);

			startx = - cos (theta) * RADIUS;  // unit components
			starty =   sin (theta) * RADIUS;

			if ((p_nuc->geometry == nucleator::sphere) || (axis == zaxis))
			{	// sphere or capsule z
				linestartx[axis][i] = startx;  // line start
				linestarty[axis][i] = starty;
		        
				lineunitvecx[axis][i] = - startx;
				lineunitvecy[axis][i] = - starty;
			}
			else
			{	
				// capsule x or y
				linestartx[axis][i+num_straight_segs] = startx ;  // line start
				linestarty[axis][i+num_straight_segs] = starty - CAPSULE_HALF_LINEAR;
		        
				lineunitvecx[axis][i+num_straight_segs] = - startx;
				lineunitvecy[axis][i+num_straight_segs] = - starty;
			}
		}

		// sides

		if  ((p_nuc->geometry == nucleator::capsule) && (axis != zaxis))
		{
			
			for(int i = 0; i != num_straight_segs; ++i)
			{

				dist = (0.5 + (double) i) * straight_seg_len;

				startx = RADIUS;  // unit components
				starty = CAPSULE_HALF_LINEAR - dist ;

				linestartx[axis][i+num_cap_segs] = startx;
				linestarty[axis][i+num_cap_segs] = starty;

				// N.B. this is scaled by (cap_seg_len/straight_seg_len) to compenate for 
				// difference in segments lengths
		       
				lineunitvecx[axis][i+num_cap_segs] = - startx;
				lineunitvecy[axis][i+num_cap_segs] = 0;
			}

			for(int i = 0; i != num_straight_segs; ++i)
			{

				dist = (0.5 + (double) i) * straight_seg_len;

				startx = - RADIUS;  // unit components
				starty = dist - CAPSULE_HALF_LINEAR ;

				linestartx[axis][i+2*num_cap_segs+num_straight_segs] = startx;
				linestarty[axis][i+2*num_cap_segs+num_straight_segs] = starty;

				// N.B. this is scaled by (cap_seg_len/straight_seg_len) to compenate for 
				// difference in segments lengths
		        
				lineunitvecx[axis][i+2*num_cap_segs+num_straight_segs] = - startx;
				lineunitvecy[axis][i+2*num_cap_segs+num_straight_segs] = 0;
			}
			
		}
	}
// N.B. TODO (important): need to derive scaled surface area
// (integral of surface multiplied by orthoganal scaling factor used
// during the impacts)
// but use lengths for now:

	straight_seg_area = straight_seg_len;// CAPSULE_HALF_LINEAR * 2 * 2 * PI * RADIUS /  (2 * num_straight_segs) ;

	curved_seg_area   = cap_seg_len;// (3 * PI * RADIUS * RADIUS * RADIUS / 4) /	(2 * num_cap_segs);

    numnodes_scalefactor = 0.1;
	rep_radial_scalefactor = 0.1;
	rep_transverse_scalefactor = 0.1;
	link_radial_scalefactor = 0.1;
	link_transverse_scalefactor = 0.1;
	links_broken_scalefactor = 0.1;
    surfaceimpacts_scalefactor = 0.1;
}

void segments::getsegmentnum(const vect& node, int& xseg, int& yseg, int& zseg) const
{  // set the segment number for each axis, dependent on node position

    double dblxseg, dblyseg, dblzseg;

    getsegmentnum(node, dblxseg, dblyseg, dblzseg);

	xseg = (int) dblxseg;
	yseg = (int) dblyseg;
	zseg = (int) dblzseg;

}

double segments::getsegmentnum(const vect& node, const projection &axis)
{
    double dblxseg, dblyseg, dblzseg;

    getsegmentnum(node, dblxseg, dblyseg, dblzseg);

    if (axis == xaxis)
    {
        return dblxseg;
    }
    else if (axis == yaxis)
    {
        return dblyseg;
    }

    return dblzseg;

}

void segments::getsegmentnum(const vect& node, double& xseg, double& yseg, double& zseg) const
{  // set the segment number for each axis, dependent on node position

	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		xseg = (double) num_cap_segs * (1+ (atan2(-node.z,node.y) / PI));
		yseg = (double) num_cap_segs * (1+ (atan2(-node.z,node.x) / PI));
		zseg = (double) num_cap_segs * (1+ (atan2(-node.y,node.x) / PI));
		
		return;
	}
	else
	{	// capsule
	
		xseg = getcapsuleseg(node.y,node.z);
		yseg = getcapsuleseg(node.x,node.z);
		zseg = (double)num_cap_segs * (1+ (atan2(-node.y,node.x) / PI));
	}

}

double segments::getcapsuleseg(const double & x, const double & y) const
{	// get the segment num for the capsule (fractional)
	// segments are numbered clockwise
	// starting at upper cap on left most edge

	
	if (fabs(y) < CAPSULE_HALF_LINEAR)
	{	// on cylinder

		if (x > 0)	// RHS
		{
			return (double) num_cap_segs + ((CAPSULE_HALF_LINEAR - y) / straight_seg_len);
		}
		else		// LHS
		{
			return 2 * (double) num_cap_segs + (double)num_straight_segs 
						        + ((y + CAPSULE_HALF_LINEAR) / straight_seg_len);
		}
	} 
	else
	{	
		if (y > 0)
		{
			// top cap
			return ((double)num_cap_segs * (1+ (atan2( -(y - CAPSULE_HALF_LINEAR),  x) / PI)) );
 		}
        else
		{	// bottom cap
			return num_straight_segs + num_cap_segs
				+  ((double)num_cap_segs * (1+ (atan2(  (y + CAPSULE_HALF_LINEAR), -x) / PI)) );
		}
	}

}

void segments::getsegmentdist(const nodes& node,int& xdist, int& ydist, int& zdist) const
{

	vect rot_pos, rot_unit;

	rot_pos = node;

	//ptheactin->camera_rotation.rotate(rot_pos); 
	
	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		xdist = dist_to_seg(calcdist(rot_pos.y, rot_pos.z) - RADIUS); 
		ydist = dist_to_seg(calcdist(rot_pos.x, rot_pos.z) - RADIUS); 
		zdist = dist_to_seg(calcdist(rot_pos.x, rot_pos.y) - RADIUS); 
		return;
	}
	else
	{	// capsule

		zdist = dist_to_seg(calcdist(rot_pos.x, rot_pos.y) - RADIUS); 

		if (node.onseg)
		{	// on cylinder

			xdist = dist_to_seg(fabs(rot_pos.y) - RADIUS); 
			ydist = dist_to_seg(fabs(rot_pos.x) - RADIUS); 
		
		} 
		else if (node.z > 0)
		{	// top cap
	
			xdist = dist_to_seg(calcdist(rot_pos.y, rot_pos.z - CAPSULE_HALF_LINEAR) - RADIUS); 
			ydist = dist_to_seg(calcdist(rot_pos.x, rot_pos.z - CAPSULE_HALF_LINEAR) - RADIUS); 

		}
		else	//  (y < CAPSULE_HALF_LINEAR)
		{	// bottom cap

			xdist = dist_to_seg(calcdist(rot_pos.y, rot_pos.z + CAPSULE_HALF_LINEAR) - RADIUS); 
			ydist = dist_to_seg(calcdist(rot_pos.x, rot_pos.z + CAPSULE_HALF_LINEAR) - RADIUS); 
	
		}

	}

}


int segments::dist_to_seg(const double & dist) const
{	// convert distance to seg num, and return -1 if out of range

	if (dist < 0) return -1;
	
	int dist_seg = (int) (dist / dist_step);

	return (dist_seg>=num_dist_segs)?(-1):(dist_seg);
}



void segments::addnode(const nodes& node)
{	// add the node's data to the segments

	// TODO:  fix this for capsule rotation
    
	int xseg, yseg, zseg;
	int xdist, ydist, zdist;
	int dist;
	double xfactor,yfactor,zfactor;
	//double radius;

	nodes           rot_pos = node;
	vect           rot_unit = node.unit_vec_posn;
	vect rot_nuc_link_force = node.nucleator_link_force;
	
	vect xfacvec, yfacvec, zfacvec;

	ptheactin->camera_rotation.rotate(rot_pos); 
	ptheactin->camera_rotation.rotate(rot_unit); 
	ptheactin->camera_rotation.rotate(rot_nuc_link_force);

    // there is an annoying rounding error problem for 
    // capsule when fabs(rot_pos.z) == CAPSULE_HALF_LINEAR
    // [caused by the atan2() functions in getsegmentnum()]
    // cludge:

    if (fabs(rot_pos.z) - CAPSULE_HALF_LINEAR < 0.000001)
    {
        rot_pos.z += 0.000001;
    }

	getsegmentnum(rot_pos, xseg, yseg, zseg);
	getsegmentdist(rot_pos, xdist, ydist, zdist);

	// find scaling factor if off-axis, using rotated vector
	// and radius using normal pos'n

	xfacvec = vect(0,rot_unit.y,rot_unit.z); 
	yfacvec = vect(rot_unit.x,0,rot_unit.z); 

	//xfactor = calcdist(rot_unit.y,rot_unit.z);
	//yfactor = calcdist(rot_unit.x,rot_unit.z);

	if (p_nuc->is_sphere())
	{
		zfacvec = vect(rot_unit.x,rot_unit.y,0);
		//zfactor = calcdist(rot_unit.x,rot_unit.y);	
		//radius = node.length() - RADIUS;
	}
	else
	{
		if (node.onseg)
		{  // on cylinder, don't scale by z component
			zfacvec = vect(rot_unit.x,rot_unit.y,0);
			zfactor = 1;
			//radius = calcdist(node.x,node.y) - RADIUS;
		}
		else
		{	// on ends
			zfacvec = vect(rot_unit.x,rot_unit.y,(fabs(node.z) - CAPSULE_HALF_LINEAR) - RADIUS); 

			zfactor = calcdist(rot_unit.x,rot_unit.y);
			//radius = calcdist(node.x,node.y,fabs(node.z) - CAPSULE_HALF_LINEAR) - RADIUS;
		}

	}

	xfactor = xfacvec.length();
	yfactor = yfacvec.length();
	zfactor = zfacvec.length();

	//radius = node.dist_from_surface;

	surfaceimpacts[0][xseg] += xfactor * node.nucleator_impacts; 
	surfaceimpacts[1][yseg] += yfactor * node.nucleator_impacts; 
	surfaceimpacts[2][zseg] += zfactor * node.nucleator_impacts;

	//if ((xseg==11))
	//	cout << "node.nucleator_impacts: "  << setw(30) << setprecision(20) << node.nucleator_impacts << endl;

	// check this, may not be right

	surfacestuckforce[0][xseg][0] += xfactor * rot_nuc_link_force.y; 
	surfacestuckforce[0][xseg][1] += xfactor * rot_nuc_link_force.z;

	surfacestuckforce[1][yseg][0] += yfactor * rot_nuc_link_force.x;
	surfacestuckforce[1][yseg][1] += yfactor * rot_nuc_link_force.z; 

	surfacestuckforce[2][zseg][0] += zfactor * rot_nuc_link_force.x;
	surfacestuckforce[2][zseg][1] += zfactor * rot_nuc_link_force.y;


	if (xdist!=-1)
	{
              numnodes[0][xseg][xdist] += xfactor;
			rep_radial[0][xseg][xdist] += xfactor * node.repforce_radial;
		rep_transverse[0][xseg][xdist] += xfactor * node.repforce_transverse;
	   	   link_radial[0][xseg][xdist] += xfactor * node.linkforce_radial;
	   link_transverse[0][xseg][xdist] += xfactor * node.linkforce_transverse;
		  links_broken[0][xseg][xdist] += xfactor * node.links_broken;
	}


	if (ydist!=-1)
	{
              numnodes[1][yseg][ydist] += yfactor;
			rep_radial[1][yseg][ydist] += yfactor * node.repforce_radial;
		rep_transverse[1][yseg][ydist] += yfactor * node.repforce_transverse;
		   link_radial[1][yseg][ydist] += yfactor * node.linkforce_radial;
	   link_transverse[1][yseg][ydist] += yfactor * node.linkforce_transverse;
		  links_broken[1][yseg][ydist] += yfactor * node.links_broken;
	}


	if (zdist!=-1)
	{
              numnodes[2][zseg][zdist] += zfactor;
			rep_radial[2][zseg][zdist] += zfactor * node.repforce_radial;
		rep_transverse[2][zseg][zdist] += zfactor * node.repforce_transverse;
		   link_radial[2][zseg][zdist] += zfactor * node.linkforce_radial;
	   link_transverse[2][zseg][zdist] += zfactor * node.linkforce_transverse;
		  links_broken[2][zseg][zdist] += zfactor * node.links_broken;

    }

	dist = (int) (node.dist_from_surface / radialdist);

	if ((dist < num_radial_bins) && (dist >= 0))
    {
	    radial_numnodes[dist]++;
	    radial_rep_radial[dist] += node.repforce_radial; 
	    radial_rep_transverse[dist] += node.repforce_transverse; 
	    radial_link_radial[dist] += node.linkforce_radial;
	    radial_link_transverse[dist] += node.linkforce_transverse;
	    radial_links_broken[dist] += node.links_broken;

        radial_linkenergy_transverse[dist] += node.linkenergy_transverse;
        radial_linkenergy_radial[dist] += node.linkenergy_radial;
		radial_repenergy_transverse[dist] += node.repenergy_transverse;
        radial_repenergy_radial[dist] += node.repenergy_radial;
		radial_linksenergy_broken[dist] += node.linksenergy_broken;
    }

}

void segments::clearbins(void)
{
	for (int axis = 0; axis != 3; ++axis)
	{
		for(int seg = 0; seg != num_segs; ++seg)
		{
			surfaceimpacts[axis][seg]=0;

			for(int dist = 0; dist<num_dist_segs; ++dist)
			{
					   numnodes[axis][seg][dist]=0;
					 rep_radial[axis][seg][dist]=0;
				 rep_transverse[axis][seg][dist]=0;
					link_radial[axis][seg][dist]=0;
				link_transverse[axis][seg][dist]=0;
   				   links_broken[axis][seg][dist]=0;
			  
			}

			for(int dirn = 0; dirn<2; ++dirn)
			{
				surfacestuckforce[axis][seg][dirn]=0;
			}

		}
	}
	
	for(int dist = 0; dist<num_dist_segs; ++dist)
	{
		radial_numnodes[dist] =			
		radial_rep_radial[dist] =
		radial_rep_transverse[dist] =
		radial_link_radial[dist] =
		radial_link_transverse[dist] =
		radial_links_broken[dist] = 0.0;

        radial_linkenergy_transverse[dist] =
        radial_linkenergy_radial[dist] =
		radial_repenergy_transverse[dist] =
        radial_repenergy_radial[dist] =
		radial_linksenergy_broken[dist] =0.0;
	}

    for (int axis = 0; axis != 3; ++axis)
	{
    	for(int dist = 0; dist != num_dist_segs; ++dist)
		{
            numnodes_SD[axis][dist] =
		    rep_radial_SD[axis][dist] =
		    rep_transverse_SD[axis][dist] =
		    link_radial_SD[axis][dist] =
		    link_transverse_SD[axis][dist] =
		    links_broken_SD[axis][dist] = 0.0;
        }
    }
}

void segments::drawoutline(ostream& drawcmd, const projection & axis) const
{

    int radius;
    int segment;

	radius  = ptheactin->pixels(RADIUS); 
	segment = ptheactin->pixels(CAPSULE_HALF_LINEAR);
	
	// draw outline

    if	((p_nuc->geometry == nucleator::sphere) || (axis == zaxis))
	{	
		// sphere

			drawcmd << " circle "
			<< centerx << "," << centery
			<< " "
			<< centerx << "," << centery+radius;
    } 
	else 
	{	
		// capsule

		drawcmd << " ellipse "
			<< centerx << "," << centery + segment << " "
			<< radius << "," << radius << " "
			<< "0,180";
		
		drawcmd << " line "
			<< centerx - radius << "," << centery + segment << " "
			<< centerx - radius << "," << centery - segment;
		
		drawcmd << " line "
			<< centerx + radius << "," << centery + segment << " "
			<< centerx + radius << "," << centery - segment;
		
		drawcmd << " ellipse "
			<< centerx << "," << centery - segment << " "
			<< radius << "," << radius << " "
			<< "180,360";
    }
 
}

int segments::drawsurfaceimpacts(ostream& drawcmd, const projection & axis, const double scale) const
{

	double linelen, linex, liney;
	double startx, starty;

	double seg_area, unscaledlen;

	int numlinesplotted = 0;

	int segstodraw = num_segs;

	if ((p_nuc->geometry == nucleator::capsule) && (axis == zaxis))
	{	// capsule z axis, no linear section to plot
		segstodraw = 2 * num_cap_segs;
	}
	
	for (int i=0; i<segstodraw; ++i)
	{

		startx = linestartx[axis][i];
		starty = linestarty[axis][i];

		if ((p_nuc->geometry == nucleator::capsule) && 
			(axis != 2) &&
			(fabs(starty) < CAPSULE_HALF_LINEAR))
		{	// on capsule side
			seg_area = straight_seg_area;
		}
		else
		{
			seg_area = curved_seg_area;
		}
	
		// draw the surface impacts:

		unscaledlen = surfaceimpacts[axis][i] * scale
                            / (seg_area * surfaceimpacts_scalefactor);

		//cout << setw(10) << setprecision(5) << unscaledlen << endl;

		//cout << "surfaceimpacts[" << axis <<"][" <<i<<"]" << surfaceimpacts[axis][i] << endl;

		linex = unscaledlen * lineunitvecx[axis][i];
		liney = unscaledlen * lineunitvecy[axis][i];

		linelen = calcdist(linex,liney);

		// truncate if too long
		if (linelen > 0.9 * RADIUS)
		{
			linex *= 0.9 * RADIUS / linelen;
			liney *= 0.9 * RADIUS / linelen;			
		}

		// don't plot zero length lines
		if ( (ptheactin->pixels(startx) == ptheactin->pixels(startx + linex)) &&
			 (ptheactin->pixels(starty) == ptheactin->pixels(starty + liney)) )
		{   // do nothing
		}
		else
		{

		drawcmd << " line "
				<< centerx + ptheactin->pixels(startx) << "," 
				<< centery + ptheactin->pixels(starty) << " "

				<< centerx + ptheactin->pixels(startx + linex) << ","
				<< centery + ptheactin->pixels(starty + liney);
				
				numlinesplotted++;

		}

		if (STICK_TO_NUCLEATOR)
		{

			// draw the nucleator link forces:

			unscaledlen = calcdist(surfacestuckforce[axis][i][0],surfacestuckforce[axis][i][1]) * scale
								/ (seg_area * surfaceimpacts_scalefactor);  // same scale factor as above

			linex = - surfacestuckforce[axis][i][0] * scale / (seg_area * surfaceimpacts_scalefactor);
			liney = - surfacestuckforce[axis][i][1] * scale / (seg_area * surfaceimpacts_scalefactor);

			linelen = calcdist(linex,liney);

			// truncate if too long (this can be longer than 0.9 RADIUS, if want)
			if (linelen > 0.9 * RADIUS)
			{
				linex *= 0.9 * RADIUS / linelen;
				liney *= 0.9 * RADIUS / linelen;			
			}

			// don't plot zero length lines
			if ( (ptheactin->pixels(startx) == ptheactin->pixels(startx + linex)) &&
				 (ptheactin->pixels(starty) == ptheactin->pixels(starty + liney)) )
			{	// do nothing
			}
			else
			{

			drawcmd << " line "
					<< centerx + ptheactin->pixels(startx) << "," 
					<< centery + ptheactin->pixels(starty) << " "

					<< centerx + ptheactin->pixels(startx + linex) << ","
					<< centery + ptheactin->pixels(starty + liney);
					
					numlinesplotted++;

			}
		}

	}

return numlinesplotted;
}

void segments::addallnodes()
{
    clearbins();

	for (int i=0; i<ptheactin->highestnodecount; i++)
	{
		if ((ptheactin->node[i].polymer))  // is point valid?
		{
			addnode(ptheactin->node[i]);
		}
	}

    // calculate stddev's from bins

    calcSD(numnodes,numnodes_SD);
    calcSD(rep_radial,rep_radial_SD);
    calcSD(rep_transverse,rep_transverse_SD);
    calcSD(link_radial,link_radial_SD);
    calcSD(link_transverse,link_transverse_SD);
    calcSD(links_broken,links_broken_SD);

}

void segments::calcSD(const Dbl3d& data, Dbl2d& SD)
{
    Dbl2d mean;
    mean.resize(3);
    for (int axis = 0; axis != 3; ++axis)
    {
	    mean[axis].resize(num_dist_segs);
        fill(mean[axis].begin(), mean[axis].end(), 0);
    }

    // calc sums for means

    for (int axis = 0; axis != 3; ++axis)
	{
	    for(int seg = 0; seg != num_segs; ++seg)
		{
		    for(int dist = 0; dist<num_dist_segs; ++dist)
    		{
                mean[axis][dist] += data[axis][seg][dist];
			}
		}
	}

    // convert sums to means
    for (int axis = 0; axis != 3; ++axis)
	{
	    for(int dist = 0; dist != num_dist_segs; ++dist)
        {
            mean[axis][dist] /= num_segs;
	    }
    }
    
    // sum the squares

    for (int axis = 0; axis != 3; ++axis)
	{
	    for(int seg = 0; seg != num_segs; ++seg)
		{
		    for(int dist = 0; dist<num_dist_segs; ++dist)
    		{
                SD[axis][dist] += data[axis][seg][dist] * data[axis][seg][dist];
			}
		}
	}

    // convert to means of squares
    
    for (int axis = 0; axis != 3; ++axis)
	{
        for(int dist = 0; dist != num_dist_segs; ++dist)
        {
            SD[axis][dist] /= num_segs;
	    }
    }

    // subtract square of mean

    for (int axis = 0; axis != 3; ++axis)
	{
        for(int dist = 0; dist != num_dist_segs; ++dist)
        {
            SD[axis][dist] -= mean[axis][dist] * mean[axis][dist];
	    }
    }

    // and take squareroot
    for (int axis = 0; axis != 3; ++axis)
	{
        for(int dist = 0; dist != num_dist_segs; ++dist)
        {
            if (SD[axis][dist]>SQRT_ACCURACY_LOSS)
                SD[axis][dist] = sqrt(SD[axis][dist]);
	    }
    }

}


void segments::savereport(const int& filenum) const
{

	char filename[255];

	double x,y,z;
	double radius;

	bool capsuleside;

	sprintf ( filename , "%sreport%05i.txt",TEMPDIR, filenum );

	ofstream opreport(filename, ios::out | ios::trunc);
	if (!opreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

	// write header
	
	opreport << setprecision(4);

	opreport << "Axis,segment,distseg,x,y,z,Radius,area,capsuleside,numnodes,RepForceRadial,RepForceTrans,LinkForceRadial,LinkForceTrans,LinksBroken" << endl;


	for(int dist = 0; dist<num_dist_segs; ++dist)
	{
		for (int axis = 0; axis != 3; ++axis)
		{
			for(int seg = 0; seg != num_segs; ++seg)
			{
				if ((axis == zaxis) && (seg >= 2 * num_cap_segs))
					continue;  // no linear segment on z axis

				getsegmentposition(x, y, z, seg, dist, axis);

				capsuleside = false;

                if (axis == xaxis)
                {
                    if (p_nuc->geometry == nucleator::sphere)
				    {	// sphere or capsule z
					    radius = calcdist(y,z) - RADIUS;
				    }
				    else
				    {
					    if (fabs(z) < CAPSULE_HALF_LINEAR)
					    {
						    capsuleside = true;
						    radius = fabs(x) - RADIUS;
					    }
					    else
					    {
						    radius = calcdist(x , fabs(z) - CAPSULE_HALF_LINEAR) - RADIUS;
					    }
				    }
                }
                else if (axis == yaxis)
                {
                    if ((p_nuc->geometry == nucleator::sphere) || (axis == zaxis))
				    {	// sphere or capsule z
					    radius = calcdist(x,z) - RADIUS;
				    }
				    else
				    {
					    if (fabs(z) < CAPSULE_HALF_LINEAR)
					    {
						    capsuleside = true;
						    radius = fabs(x) - RADIUS;
					    }
					    else
					    {
						    radius = calcdist(x , fabs(z) - CAPSULE_HALF_LINEAR) - RADIUS;
					    }
				    }
                }
                else
                {
				    radius = calcdist(x,y) - RADIUS;
		        }


                // rotate *after* determining radius

                //ptheactin->camera_rotation.rotate(x,y,z);

// Axis,segment,distseg,x,y,z,Radius,numnodes,area,capsuleside,RepForceRadial,RepForceTrans,RepDisplRadial,RepDisplTrans,LinkForceRadial,LinkForceTrans" << endl;

				opreport 
					<< axis << ","
					<< seg << ","
					<< dist << ","
					<< x << ","
					<< y << ","
					<< z << ","
					<< radius << ","
					<< 0 << ","
					<< capsuleside << ","
					<< numnodes[axis][seg][dist] << ","
					<< rep_radial[axis][seg][dist] << ","
					<< rep_transverse[axis][seg][dist] << ","
					<< link_radial[axis][seg][dist] << ","
					<< link_transverse[axis][seg][dist] << ","
					<< links_broken[axis][seg][dist] << endl;

			}
		}
//	opradialreport << "RadialSegDist,numnodes,RepForceRadial,RepForceTrans,RepDisplRadial,RepDisplTrans,LinkForceRadial,LinkForceTrans" << endl;



	}

	opreport << endl;

	opreport.close();

}

void segments::saveradialreport(const int& filenum) const
{

	char filename [255];

	sprintf ( filename , "%sreport_radial%05i.txt", TEMPDIR, filenum );

	ofstream opradialreport(filename, ios::out | ios::trunc);
	if (!opradialreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

    opradialreport << "Distance,Nodes,F_Rep(radial),F_Rep(trans),F_Link(radial),F_Link(trans),LinksBroken"
                   << ",E_Rep(radial),E_Rep(trans),E_Link(radial),E_Link(trans),E_LinksBroken" << endl;

	for (int dist = 0; dist < num_radial_bins; ++dist)
	{

		opradialreport 
			<< ((double)dist + 0.5) * radialdist << ","
			<< radial_numnodes[dist] << ","
			<< radial_rep_radial[dist] << ","
			<< radial_rep_transverse[dist] << ","
			<< radial_link_radial[dist] << ","
			<< radial_link_transverse[dist] << ","
			<< radial_links_broken[dist] << ","
            << radial_repenergy_radial[dist] << ","
            << radial_repenergy_transverse[dist] << ","
            << radial_linkenergy_radial[dist] << ","
            << radial_linkenergy_transverse[dist] << ","
		    << radial_linksenergy_broken[dist]
            << endl;
	}

	opradialreport << endl;
	opradialreport.close();
}

void segments::saveradialaxisreport(const int& filenum, const int axis) const
{

	char filename [255];

    char projletter[] = "z";

    if (axis == xaxis)
	{
		sprintf ( projletter , "x");
	}
	else if (axis == yaxis)
	{
		sprintf ( projletter , "y");
	}

	sprintf ( filename , "%s%s_report_radial_%05i.txt", TEMPDIR,projletter, filenum );

	ofstream opradialreport(filename, ios::out | ios::trunc);
	if (!opradialreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

	opradialreport << "RadialSegDist,numnodes,RepForceRadial,RepForceTrans,LinkForceRadial,LinkForceTrans,LinksBroken" << endl;

    double temp_radial_numnodes,
        temp_radial_rep_radial,
        temp_radial_rep_transverse,
        temp_radial_link_radial,
        temp_radial_link_transverse,
        temp_radial_links_broken;


	for (int dist = 0; dist < num_radial_bins; ++dist)
	{
        temp_radial_numnodes = 0;
        temp_radial_rep_radial = 0;
        temp_radial_rep_transverse = 0;
        temp_radial_link_radial = 0;
        temp_radial_link_transverse = 0;
        temp_radial_links_broken = 0;

        for(int seg = 0; seg<num_segs; ++seg)
	    {
            temp_radial_numnodes += numnodes[axis][seg][dist];
            temp_radial_rep_radial += rep_radial[axis][seg][dist];
            temp_radial_rep_transverse += rep_transverse[axis][seg][dist];
            temp_radial_link_radial += link_radial[axis][seg][dist];
            temp_radial_link_transverse += link_transverse[axis][seg][dist];
            temp_radial_links_broken += links_broken[axis][seg][dist];
        }

		opradialreport 
			<< ((double)dist + 0.5) * radialdist << ","
			<< temp_radial_numnodes << ","
			<< temp_radial_rep_radial << ","
			<< temp_radial_rep_transverse << ","
			<< temp_radial_link_radial << ","
			<< temp_radial_link_transverse << ","
			<< temp_radial_links_broken << endl;
	}

	opradialreport << endl;
	opradialreport.close();
}

void segments::saveSDreport(const int& filenum) const
{

	char filename [255];

	sprintf ( filename , "%sreport_SD%05i.txt", TEMPDIR, filenum );

	ofstream opSDreport(filename, ios::out | ios::trunc);
	if (!opSDreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

	opSDreport << "Axis,RadialSegDist,numnodesSD,RepForceRadialSD,RepForceTransSD,LinkForceRadialSD,LinkForceTransSD,LinksBrokenSD" << endl;

    for (int axis = 0; axis != 3; ++axis)
	{
	    for (int dist = 0; dist != num_radial_bins; ++dist)
	    {

		    opSDreport << axis << ","
			    << ((double)dist + 0.5) * radialdist << ","
			    << numnodes_SD[axis][dist] << ","
			    << rep_radial_SD[axis][dist] << ","
			    << rep_transverse_SD[axis][dist] << ","
			    << link_radial_SD[axis][dist] << ","
			    << link_transverse_SD[axis][dist] << ","
			    << links_broken_SD[axis][dist] << endl;
	    }
    }

	opSDreport << endl;
	opSDreport.close();
}

void segments::getsegmentposition(double& x, double& y, double& z, const int & seg,
								  const int & dist, const int & axis) const
{
	if (axis == xaxis)
	{
		x = 0;
		y = linestartx[axis][seg] - lineunitvecx[axis][seg] * dist * dist_step;
		z = linestarty[axis][seg] - lineunitvecy[axis][seg] * dist * dist_step;
	}
	else if (axis == yaxis)
	{
		x = linestartx[axis][seg] - lineunitvecx[axis][seg] * dist * dist_step;
		y = 0;
		z = linestarty[axis][seg] - lineunitvecy[axis][seg] * dist * dist_step;
	}
	else
	{
		x = linestartx[axis][seg] - lineunitvecx[axis][seg] * dist * dist_step;
		y = linestarty[axis][seg] - lineunitvecy[axis][seg] * dist * dist_step;
		z = 0;
	}

}

//void segments::setbitmapcoords()
//{
//
//	for (int axis = 0; axis < 3; ++axis)
//		for(int i = 0; i<num_segs; ++i)
//			for(int j = 0; j<num_dist_segs; ++j)
//			{
//				segment_pixel[axis][i][j].resize(0);
//			}
//
//	const int offsetx = bins_bitmap_width / 2;			// center these pixels
//	const int offsety = bins_bitmap_height / 2;
//
//	int pix_x, pix_y;	// these are bitmap co-ords
//
//	int xseg, yseg, zseg;
//	int xdist, ydist, zdist;
//	int seg,dist;
//
//	nodes dummynode;
//
//	// go through bitmap
//
//	for (pix_x = 0; pix_x < bins_bitmap_width; ++pix_x)
//		for (pix_y = 0; pix_y < bins_bitmap_height; ++pix_y)
//			for(int axis = 0; axis < 3; ++axis)
//			{
//				if (axis == 0)
//				{
//					dummynode.y = ptheactin->unpixels(pix_x - offsetx);
//					dummynode.z = ptheactin->unpixels(pix_y - offsety);
//				}
//				else if (axis == 1)
//				{
//					dummynode.x = ptheactin->unpixels(pix_x - offsetx);
//					dummynode.z = ptheactin->unpixels(pix_y - offsety);
//				}
//				else
//				{
//					dummynode.x = ptheactin->unpixels(pix_x - offsetx);
//					dummynode.y = ptheactin->unpixels(pix_y - offsety);
//				}
//				
//				//ptheactin->reverse_camera_rotation.rotate(dummynode); 
//		
//				dummynode.setunitvec();
//
//				// get bin numbers:
//
//				getsegmentnum( dummynode, xseg,  yseg,  zseg);
//				getsegmentdist(dummynode, xdist, ydist, zdist);
//
//				if (axis == 0)
//				{
//					seg = xseg;
//					dist = xdist;
//				}
//				else if (axis == 1)
//				{
//					seg = yseg;
//					dist = ydist;
//				}
//				else
//				{
//					seg = zseg;
//					dist = zdist;
//				}
//
//				if (dist == -1)
//					continue;
//
//				segment_pixel[axis][seg][dist].push_back(pix_x + centerx - offsetx);
//				segment_pixel[axis][seg][dist].push_back(pix_y + centery - offsety);
//
//			}
//
//}

void segments::write_bins_bitmap(Dbl2d &imageR, Dbl2d &imageG, Dbl2d &imageB,
					   const Dbl3d & var, const double& scale, const projection &axis)
{

	const int offsetx = bins_bitmap_width / 2;			// center these pixels
	const int offsety = bins_bitmap_height / 2;

	//int x, y;	// these are bitmap co-ords

	int xseg, yseg, zseg;
	int xdist, ydist, zdist;
	int seg,dist;

	int pix_x,pix_y;

	nodes dummynode;
    Colour colour;

	double value;

	// go through bitmap

	for (pix_x = 0; pix_x < bins_bitmap_width; ++pix_x)
		for (pix_y = 0; pix_y < bins_bitmap_height; ++pix_y)
		{
			if (axis == xaxis)
			{
				dummynode.y = ptheactin->unpixels(pix_x - offsetx);
				dummynode.z = ptheactin->unpixels(pix_y - offsety);
			}
			else if (axis == yaxis)
			{
				dummynode.x = ptheactin->unpixels(pix_x - offsetx);
				dummynode.z = ptheactin->unpixels(pix_y - offsety);
			}
			else
			{
				dummynode.x = ptheactin->unpixels(pix_x - offsetx);
				dummynode.y = ptheactin->unpixels(pix_y - offsety);
			}
			
			//ptheactin->reverse_camera_rotation.rotate(dummynode); 
	
			dummynode.setunitvec();

			// get bin numbers:

			getsegmentnum( dummynode, xseg,  yseg,  zseg);
			getsegmentdist(dummynode, xdist, ydist, zdist);

			if (axis == xaxis)
			{
				seg = xseg;
				dist = xdist;
			}
			else if (axis == yaxis)
			{
				seg = yseg;
				dist = ydist;
			}
			else
			{
				seg = zseg;
				dist = zdist;
			}

			if (dist == -1)
				continue;

			value = var[axis][seg][dist] / scale;  // scaled between 0 and 1

			if (value < 0.00001)	// skip if zero
				continue;

            

			colour.setcol(value);

			imageR[pix_x + centerx - offsetx][pix_y + centery - offsety] = colour.r;
			imageG[pix_x + centerx - offsetx][pix_y + centery - offsety] = colour.g;
			imageB[pix_x + centerx - offsetx][pix_y + centery - offsety] = colour.b;

		}

        write_colourmap_bitmap(imageR, imageG, imageB);

}

void segments::write_colourmap_bitmap(Dbl2d &imageR, Dbl2d &imageG, Dbl2d &imageB)
{

    // write out colormap key

    const int keyheight = bins_bitmap_height;
    const int keywidth  = 10 * BMP_AA_FACTOR;

    const int keyxorig = 1 * BMP_AA_FACTOR; //centerx - offsetx;
    const int keyyorig = centery - (keyheight / 2);

    Colour colour;
    int pix_x,pix_y;

	    for (pix_y = 0; pix_y < keyheight; ++pix_y)
	    {
		    colour.setcol(1 - (double)pix_y / (double)keyheight);

		    for (pix_x = 0; pix_x < keywidth; ++pix_x)
		    {
			    imageR[pix_x + keyxorig][pix_y + keyyorig] = colour.r;
			    imageG[pix_x + keyxorig][pix_y + keyyorig] = colour.g;
			    imageB[pix_x + keyxorig][pix_y + keyyorig] = colour.b;
		    }
	    }

}

void segments::set_scale_factors(void)
{
    for (int ax = 0; ax < 3; ++ax)
	    for(int seg = 0; seg<num_segs; ++seg)
        {

            if (surfaceimpacts[ax][seg] > surfaceimpacts_scalefactor)
				surfaceimpacts_scalefactor = surfaceimpacts[ax][seg];

		    for(int dist = 0; dist<num_dist_segs; ++dist)
		    {
			    if (numnodes[ax][seg][dist] > numnodes_scalefactor)
				    numnodes_scalefactor = numnodes[ax][seg][dist];

                if (rep_radial[ax][seg][dist] > rep_radial_scalefactor)
				    rep_radial_scalefactor = rep_radial[ax][seg][dist];

                if (rep_transverse[ax][seg][dist] > rep_transverse_scalefactor)
				    rep_transverse_scalefactor = rep_transverse[ax][seg][dist];

                if (link_radial[ax][seg][dist] > link_radial_scalefactor)
				    link_radial_scalefactor = link_radial[ax][seg][dist];

                if (link_transverse[ax][seg][dist] > link_transverse_scalefactor)
				    link_transverse_scalefactor = link_transverse[ax][seg][dist];

                if (links_broken[ax][seg][dist] > links_broken_scalefactor)
				    links_broken_scalefactor = links_broken[ax][seg][dist];

			}
        }

}


void segments::save_scalefactors(void)
{
	ofstream opscalefact("segscalefactors.txt", ios::out | ios::trunc);
	if (!opscalefact) 
	{ cout << "Unable to open file 'sym_break_axis.txt' for output"; return;}

	opscalefact  << surfaceimpacts_scalefactor << " "
				<< numnodes_scalefactor << " "
				<< rep_radial_scalefactor << " "
                << rep_transverse_scalefactor << " "
                << link_radial_scalefactor << " "
                << link_transverse_scalefactor << " "
                << links_broken_scalefactor << endl
                << ptheactin->imageRmax[0] << " "
                << ptheactin->imageGmax[0] << " "
                << ptheactin->imageBmax[0] << endl
                << ptheactin->imageRmax[1] << " "
                << ptheactin->imageGmax[1] << " "
                << ptheactin->imageBmax[1] << endl
                << ptheactin->imageRmax[2] << " "
                << ptheactin->imageGmax[2] << " "
                << ptheactin->imageBmax[2] << endl;

	opscalefact.close();

	//cout << "'sym_break_axis.txt' file written" << endl;
}

bool segments::load_scalefactors(void)
{
	ifstream ipscalefact(SEG_SCALE_FILE, ios::in);

	if (!ipscalefact) 
	{ 
		cout << "Unable to open file 'segscalefactors.txt' for input, skipping." << endl;
        return false;
	}
	else
	{
	    ipscalefact  >> surfaceimpacts_scalefactor 
				    >> numnodes_scalefactor 
				    >> rep_radial_scalefactor 
                    >> rep_transverse_scalefactor 
                    >> link_radial_scalefactor 
                    >> link_transverse_scalefactor 
                    >> links_broken_scalefactor
                    >> ptheactin->imageRmax[0]
                    >> ptheactin->imageGmax[0]
                    >> ptheactin->imageBmax[0]
                    >> ptheactin->imageRmax[1]
                    >> ptheactin->imageGmax[1]
                    >> ptheactin->imageBmax[1]
                    >> ptheactin->imageRmax[2]
                    >> ptheactin->imageGmax[2]
                    >> ptheactin->imageBmax[2];

		ipscalefact.close();
        return true;
	}
}

