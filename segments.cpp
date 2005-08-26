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
	p_actin = pactin;
	
	num_dist_segs = 20;		// normal bins radial distance
	dist_step = 0.2;

	radialdist = 0.1;		// for radial-only averaging
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

	// allocate and zero 1d surfaceimpacts vector

	surfaceimpacts.resize(3);
	for (int axis = 0; axis < 3; ++axis)
	{
		surfaceimpacts[axis].resize(num_segs);
	}

	//clearsurfaceimpacts();


	// allocate and zero the 2d vectors

	       numnodes.resize(3);			
	     rep_radial.resize(3);
	 rep_transverse.resize(3);
	    link_radial.resize(3);
	link_transverse.resize(3);
	   links_broken.resize(3);
//	    disp_radial.resize(3);
//	disp_transverse.resize(3);

//	segment_pixel.resize(3);

	for (int axis = 0; axis < 3; ++axis)
	{
		       numnodes[axis].resize(num_segs);			
		     rep_radial[axis].resize(num_segs);
		 rep_transverse[axis].resize(num_segs);
		    link_radial[axis].resize(num_segs);
		link_transverse[axis].resize(num_segs);
		   links_broken[axis].resize(num_segs);

		//segment_pixel[axis].resize(num_segs);

		for(int i = 0; i<num_segs; ++i)
		{
			       numnodes[axis][i].resize(num_dist_segs);
			     rep_radial[axis][i].resize(num_dist_segs);
			 rep_transverse[axis][i].resize(num_dist_segs);
			    link_radial[axis][i].resize(num_dist_segs);
			link_transverse[axis][i].resize(num_dist_segs);
			   links_broken[axis][i].resize(num_dist_segs);

			//segment_pixel[axis][i].resize(num_dist_segs);

			//for(int j = 0; j<num_dist_segs; ++j)
			//{
			//	segment_pixel[axis][i][j].resize(0);

			//}

		}
	}


	radial_numnodes.resize(num_radial_bins);				
	radial_rep_radial.resize(num_radial_bins);	
	radial_rep_transverse.resize(num_radial_bins);	
	radial_link_radial.resize(num_radial_bins);	
	radial_link_transverse.resize(num_radial_bins);	
	radial_links_broken.resize(num_radial_bins);	

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

	for (int axis = 0; axis < 3; ++axis)
	{
		// top cap

		for(int i = 0; i< num_cap_segs; ++i)
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

		for(int i = num_cap_segs; i< (2 *num_cap_segs); ++i)
		{
			//       offset     segnum                      angle of 1 seg
			theta = (0.5 + (double) i) * (PI / (double) num_cap_segs);

			startx = - cos (theta) * RADIUS;  // unit components
			starty =   sin (theta) * RADIUS;

			if ((p_nuc->geometry == nucleator::sphere) || (axis == 2))
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

		if  ((p_nuc->geometry == nucleator::capsule) && (axis != 2))
		{
			
			for(int i = 0; i< num_straight_segs; ++i)
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

			for(int i = 0; i< num_straight_segs; ++i)
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
	curved_seg_area   = cap_seg_len;
	straight_seg_area =  straight_seg_len;
	// capsule_z_seg_area = 

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

	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		xseg = (int)((double)num_cap_segs * (1+ (atan2(-node.z,node.y) / PI)) );
		yseg = (int)((double)num_cap_segs * (1+ (atan2(-node.z,node.x) / PI)) );
		zseg = (int)((double)num_cap_segs * (1+ (atan2(-node.y,node.x) / PI)) );
		
		return;
	}
	else
	{	// capsule
	
		xseg = getcapsuleseg(node.y,node.z);
		yseg = getcapsuleseg(node.x,node.z);
		zseg = (int)((double)num_cap_segs * (1+ (atan2(-node.y,node.x) / PI)) );
	}

}

int segments::getcapsuleseg(const double & x, const double & y) const
{	// get the segment num for the capsule
	// segments are numbered clockwise
	// starting at upper cap on left most edge
	
	if (fabs(y) < CAPSULE_HALF_LINEAR)
	{	// on cylinder

		if (x>0)	// RHS
		{
			return num_cap_segs + (int) (((CAPSULE_HALF_LINEAR) - y)/straight_seg_len);
		}
		else		// LHS
		{
			return 2 * num_cap_segs + num_straight_segs 
						        + (int) ((y + (CAPSULE_HALF_LINEAR))/straight_seg_len);
		}
	} 
	else if (y > 0)
	{	// top cap

		return (int)((double)num_cap_segs * (1+ (atan2(-(y - (CAPSULE_HALF_LINEAR)),x) / PI)) );

	}
	else	//  (y < CAPSULE_HALF_LINEAR)
	{	// bottom cap

		return num_straight_segs
			+ (int)((double)num_cap_segs * (1+ (atan2((y + (CAPSULE_HALF_LINEAR)),-x) / PI)) );
		
	}

}

void segments::getsegmentdist(const nodes& node,int& xdist, int& ydist, int& zdist) const
{

	vect rot_pos, rot_unit;

	rot_pos = node;

	//p_actin->camera_rotation.rotate(rot_pos); 
	
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
	double radius;

	nodes rot_pos;
	vect rot_unit;

	rot_pos = node;
	rot_unit = node.unit_vec_posn;

	p_actin->camera_rotation.rotate(rot_pos); 
	p_actin->camera_rotation.rotate(rot_unit); 

	getsegmentnum(rot_pos, xseg, yseg, zseg);
	getsegmentdist(rot_pos, xdist, ydist, zdist);

	if (xdist!=-1)
		numnodes[0][xseg][xdist] ++;

	if (ydist!=-1)	
		numnodes[1][yseg][ydist] ++;

	if (zdist!=-1)
		numnodes[2][zseg][zdist] ++;

	// find scaling factor if off-axis, using rotated vector
	// and radius using normal pos'n

	xfactor = calcdist(rot_unit.y,rot_unit.z);
	yfactor = calcdist(rot_unit.x,rot_unit.z);

	if (p_nuc->geometry == nucleator::sphere)
	{
		zfactor = calcdist(rot_unit.x,rot_unit.y);	
		radius = node.length() - RADIUS;
	}
	else
	{
		if (node.onseg)
		{  // on cylinder, don't scale by z component
			zfactor = 1;
			radius = calcdist(node.x,node.y) - RADIUS;
		}
		else
		{	// on ends

			zfactor = calcdist(rot_unit.x,rot_unit.y);	
			radius = calcdist(node.x,node.y,fabs(node.z) - CAPSULE_HALF_LINEAR) - RADIUS;
		}

	}

	surfaceimpacts[0][xseg] += xfactor * node.nucleator_impacts;
	surfaceimpacts[1][yseg] += yfactor * node.nucleator_impacts;
	surfaceimpacts[2][zseg] += zfactor * node.nucleator_impacts;


	if (xdist!=-1)
	{

			rep_radial[0][xseg][xdist] += xfactor * node.repforce_radial;
		rep_transverse[0][xseg][xdist] += xfactor * node.repforce_transverse;
		link_radial[0][xseg][xdist] += xfactor * node.linkforce_radial;
	link_transverse[0][xseg][xdist] += xfactor * node.linkforce_transverse;
		links_broken[0][xseg][xdist] += xfactor * node.links_broken;

            
//			disp_radial[0][xseg][xdist] += node.dispforce_radial;
//		disp_transverse[0][xseg][xdist] += node.dispforce_transverse;
	}


	if (ydist!=-1)
	{
			rep_radial[1][yseg][ydist] += yfactor * node.repforce_radial;
		rep_transverse[1][yseg][ydist] += yfactor * node.repforce_transverse;
		link_radial[1][yseg][ydist] += yfactor * node.linkforce_radial;
	link_transverse[1][yseg][ydist] += yfactor * node.linkforce_transverse;
		links_broken[1][yseg][ydist] += yfactor * node.links_broken;


//			disp_radial[1][yseg][ydist] += node.dispforce_radial;
//		disp_transverse[1][yseg][ydist] += node.dispforce_transverse;
	}


	if (zdist!=-1)
	{
			rep_radial[2][zseg][zdist] += zfactor * node.repforce_radial;
		rep_transverse[2][zseg][zdist] += zfactor * node.repforce_transverse;
		link_radial[2][zseg][zdist] += zfactor * node.linkforce_radial;
	link_transverse[2][zseg][zdist] += zfactor * node.linkforce_transverse;
		links_broken[2][zseg][zdist] += zfactor * node.links_broken;

    }

//			disp_radial[2][zseg][zdist] += node.dispforce_radial[threadnum];
//		disp_transverse[2][zseg][zdist] += node.dispforce_transverse[threadnum];
	
	


	dist = (int) (radius / radialdist);

	if ((dist < num_radial_bins) && (dist > 0))
    {
	    radial_numnodes[dist]++;
	    radial_rep_radial[dist] += node.repforce_radial; 
	    radial_rep_transverse[dist] += node.repforce_transverse; 
	    radial_link_radial[dist] += node.linkforce_radial;
	    radial_link_transverse[dist] += node.linkforce_transverse;
	    radial_links_broken[dist] += node.links_broken;
    }
	


}

//void segments::addsurfaceimpact(nodes& node, const double& mag)
//{	// add a surface impact
//	
//	
//	node.nucleator_impacts += mag;

 //   int xseg,yseg,zseg;
	//vect rot_pos, rot_unit;

	//rot_pos = node;
	//rot_unit = node.unit_vec_posn;

	////p_actin->camera_rotation.rotate(rot_pos); 
	////p_actin->camera_rotation.rotate(rot_unit); 

	//getsegmentnum(rot_pos, xseg, yseg, zseg);

 //   surfaceimpacts[0][xseg]+= mag * calcdist(rot_unit.x,rot_unit.y);
	//surfaceimpacts[1][yseg]+= mag * calcdist(rot_unit.x,rot_unit.z);

	//if (p_nuc->geometry == nucleator::sphere)
	//{
	//	surfaceimpacts[2][zseg]+= mag * calcdist(rot_unit.x,rot_unit.y);	
	//}
	//else
	//{
	//	if (node.onseg)
	//	{  // on cylinder, don't scale by z component
	//		surfaceimpacts[2][zseg]+= mag;
	//	}
	//	else
	//	{	// on ends

	//		surfaceimpacts[2][zseg]+= mag * calcdist(rot_unit.x,rot_unit.y);	
	//	}

	//}

//}

//void segments::clearsurfaceimpacts(void)
//{
//	for (int axis = 0; axis < 3; ++axis)
//	{
//		for(int i = 0; i<num_segs; ++i)
//		{
//			surfaceimpacts[axis][i]=0;
//		}
//	}
//}

void segments::clearbins(void)
{
	for (int axis = 0; axis < 3; ++axis)
	{
		for(int i = 0; i<num_segs; ++i)
		{
		surfaceimpacts[axis][i]=0;

			for(int j = 0; j<num_dist_segs; ++j)
			{
				       numnodes[axis][i][j]=0;
				     rep_radial[axis][i][j]=0;
				 rep_transverse[axis][i][j]=0;
				    link_radial[axis][i][j]=0;
				link_transverse[axis][i][j]=0;
	   			   links_broken[axis][i][j]=0;
//				    disp_radial[axis][i][j]=0;
//			    disp_transverse[axis][i][j]=0;
			}
		}
	}

	
	for(int i = 0; i<num_radial_bins; ++i)
	{
		radial_numnodes[i] =			
		radial_rep_radial[i] =
		radial_rep_transverse[i] =
		radial_link_radial[i] =
		radial_link_transverse[i] =
		radial_links_broken[i] = 0;
	}
}

void segments::drawoutline(ostream& drawcmd, const int& axis) const
{

    int radius;
    int segment;

	radius  = p_actin->pixels(RADIUS); 
	segment = p_actin->pixels(CAPSULE_HALF_LINEAR);
	
	// draw outline

    if	((p_nuc->geometry == nucleator::sphere) || (axis == 2))
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
 

	// draw segment divisions?



}

int segments::drawsurfaceimpacts(ostream& drawcmd, const int& axis, const double scale) const
{

	double linelen, linex, liney;
	double startx, starty;

	double seg_area, unscaledlen;

	int numlinesplotted = 0;

	int segstodraw = num_segs;

	if ((p_nuc->geometry == nucleator::capsule) && (axis == 2))
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
	
		unscaledlen = surfaceimpacts[axis][i] * scale
                            / (seg_area * surfaceimpacts_scalefactor);

		//cout << "surfaceimpacts[" << axis <<"][" <<i<<"]" << surfaceimpacts[axis][i] << endl;

		linex = unscaledlen * lineunitvecx[axis][i];
		liney = unscaledlen * lineunitvecy[axis][i];

		linelen = calcdist(linex,liney);

		if (linelen > 0.9 * RADIUS)
		{
			linex *= 0.9 * RADIUS / linelen;
			liney *= 0.9 * RADIUS / linelen;			
		}

		// don't plot zero length lines
		if ( (p_actin->pixels(startx) == p_actin->pixels(startx + linex)) &&
			 (p_actin->pixels(starty) == p_actin->pixels(starty + liney)) )
			continue;

		numlinesplotted++;

		drawcmd << " line "
				<< centerx + p_actin->pixels(startx) << "," 
				<< centery + p_actin->pixels(starty) << " "

				<< centerx + p_actin->pixels(startx + linex) << ","
				<< centery + p_actin->pixels(starty + liney);

	}

return numlinesplotted;
}

void segments::addallnodes()
{
	for (int i=0; i<p_actin->highestnodecount; i++)
	{
		if ((p_actin->node[i].polymer))  // is point valid?
		{
			addnode(p_actin->node[i]);
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
	
	opreport << "Axis,segment,distseg,x,y,z,Radius,area,capsuleside,numnodes,RepForceRadial,RepForceTrans,LinkForceRadial,LinkForceTrans,LinksBroken" << endl;


	for(int dist = 0; dist<num_dist_segs; ++dist)
	{
		for (int axis = 0; axis < 3; ++axis)
		{
			for(int seg = 0; seg<num_segs; ++seg)
			{
				if ((axis == 2) && (seg >= 2 * num_cap_segs))
					continue;  // no linear segment on z axis

				getsegmentposition(x, y, z, seg, dist, axis);

				capsuleside = false;

				if ((p_nuc->geometry == nucleator::sphere) || (axis == 2))
				{	// sphere or capsule z
					radius = calcdist(x,y) - RADIUS;
				}
				else
				{
					if (fabs(y) < CAPSULE_HALF_LINEAR)
					{
						capsuleside = true;
						radius = fabs(x) - RADIUS;
					}
					else
					{
						radius = calcdist(x , fabs(y) - CAPSULE_HALF_LINEAR) - RADIUS;
					}
				}

                // rotate *after* determining radius

                p_actin->camera_rotation.rotate(x,y,z);


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
//					<< disp_radial[axis][seg][dist] << ","
//					<< disp_transverse[axis][seg][dist] << ",";



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

	opradialreport << "RadialSegDist,numnodes,RepForceRadial,RepForceTrans,LinkForceRadial,LinkForceTrans,LinksBroken" << endl;

	for (int dist = 0; dist < num_radial_bins; ++dist)
	{

		opradialreport 
			<< ((double)dist + 0.5) * radialdist << ","
			<< radial_numnodes[dist] << ","
			<< radial_rep_radial[dist] << ","
			<< radial_rep_transverse[dist] << ","
			<< radial_link_radial[dist] << ","
			<< radial_link_transverse[dist] << ","
			<< radial_links_broken[dist] << endl;
	}

	opradialreport << endl;
	opradialreport.close();
}

void segments::getsegmentposition(double& x, double& y, double& z, const int & seg,
								  const int & dist, const int & axis) const
{
	if (axis == 0)
	{
		x = 0;
		y = linestartx[axis][seg] - lineunitvecx[axis][seg] * dist * dist_step;
		z = linestarty[axis][seg] - lineunitvecy[axis][seg] * dist * dist_step;
	}
	else if (axis == 1)
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
//					dummynode.y = p_actin->unpixels(pix_x - offsetx);
//					dummynode.z = p_actin->unpixels(pix_y - offsety);
//				}
//				else if (axis == 1)
//				{
//					dummynode.x = p_actin->unpixels(pix_x - offsetx);
//					dummynode.z = p_actin->unpixels(pix_y - offsety);
//				}
//				else
//				{
//					dummynode.x = p_actin->unpixels(pix_x - offsetx);
//					dummynode.y = p_actin->unpixels(pix_y - offsety);
//				}
//				
//				//p_actin->reverse_camera_rotation.rotate(dummynode); 
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
					   const Dbl3d & var, const double& scale, const int& axis)
{

	const int offsetx = bins_bitmap_width / 2;			// center these pixels
	const int offsety = bins_bitmap_height / 2;

	//int x, y;	// these are bitmap co-ords

	int xseg, yseg, zseg;
	int xdist, ydist, zdist;
	int seg,dist;

	int pix_x,pix_y;


	nodes dummynode;

	double value;

//	double maxval = SQRT_ACCURACY_LOSS;

	// set max value:

	//for (int ax = 0; ax < 3; ++ax)
	//	for(int seg = 0; seg<num_segs; ++seg)
	//		for(int dist = 0; dist<num_dist_segs; ++dist)
	//		{
	//			if (var[ax][seg][dist] > maxval)
	//				maxval = var[ax][seg][dist];
	//		}

	//for(int seg = 0; seg<num_segs; ++seg)
	//	for(int dist = 0; dist<num_dist_segs; ++dist)
	//	{
	//		value = var[axis][seg][dist] / maxval;  // scaled between 0 and 1

	//		if (value < 0.00001)	// skip if zero
	//			continue;

	//		dummynode.colour.setcol(value);

	//		for (vector <int>::iterator i=segment_pixel[axis][seg][dist].begin(); i<segment_pixel[axis][seg][dist].end() ; i++ )
	//		{	 

	//		x = *i;

	//		i++;
	//		y = *i;

	//		imageR[x][y] = dummynode.colour.r;
	//		imageG[x][y] = dummynode.colour.g;
	//		imageB[x][y] = dummynode.colour.b;

	//		}
		//}

	// go through bitmap

	for (pix_x = 0; pix_x < bins_bitmap_width; ++pix_x)
		for (pix_y = 0; pix_y < bins_bitmap_height; ++pix_y)
		{
			if (axis == 0)
			{
				dummynode.y = p_actin->unpixels(pix_x - offsetx);
				dummynode.z = p_actin->unpixels(pix_y - offsety);
			}
			else if (axis == 1)
			{
				dummynode.x = p_actin->unpixels(pix_x - offsetx);
				dummynode.z = p_actin->unpixels(pix_y - offsety);
			}
			else
			{
				dummynode.x = p_actin->unpixels(pix_x - offsetx);
				dummynode.y = p_actin->unpixels(pix_y - offsety);
			}
			
			//p_actin->reverse_camera_rotation.rotate(dummynode); 
	
			dummynode.setunitvec();

			// get bin numbers:

			getsegmentnum( dummynode, xseg,  yseg,  zseg);
			getsegmentdist(dummynode, xdist, ydist, zdist);

			if (axis == 0)
			{
				seg = xseg;
				dist = xdist;
			}
			else if (axis == 1)
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

			dummynode.colour.setcol(value);

			imageR[pix_x + centerx - offsetx][pix_y + centery - offsety] = dummynode.colour.r;
			imageG[pix_x + centerx - offsetx][pix_y + centery - offsety] = dummynode.colour.g;
			imageB[pix_x + centerx - offsetx][pix_y + centery - offsety] = dummynode.colour.b;

		}

// write out colormap key

const int keyheight = bins_bitmap_height;
const int keywidth  = 10;

const int keyxorig = 0; //centerx - offsetx;
const int keyyorig = centery - (keyheight / 2);


	for (pix_y = 0; pix_y < keyheight; ++pix_y)
	{
		dummynode.colour.setcol(1 - (double)pix_y / (double)keyheight);

		for (pix_x = 0; pix_x < keywidth; ++pix_x)
		{
			imageR[pix_x + keyxorig][pix_y + keyyorig] = dummynode.colour.r;
			imageG[pix_x + keyxorig][pix_y + keyyorig] = dummynode.colour.g;
			imageB[pix_x + keyxorig][pix_y + keyyorig] = dummynode.colour.b;
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
                << links_broken_scalefactor << endl;



	opscalefact.close();

	//cout << "'sym_break_axis.txt' file written" << endl;
}

void segments::load_scalefactors(void)
{
	ifstream ipscalefact("segscalefactors.txt", ios::in);

	if (!ipscalefact) 
	{ 
		cout << "Unable to open file 'segscalefactors.txt' for input, skipping." << endl;
	}
	else
	{
	    ipscalefact  >> surfaceimpacts_scalefactor 
				    >> numnodes_scalefactor 
				    >> rep_radial_scalefactor 
                    >> rep_transverse_scalefactor 
                    >> link_radial_scalefactor 
                    >> link_transverse_scalefactor 
                    >> links_broken_scalefactor;

		ipscalefact.close();
	}
}

