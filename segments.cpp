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
	
	num_dist_segs = 20;
	dist_step = 0.2;

	// calc segment lengths etc.

	curved_length = PI * RADIUS;	// one cap
	straight_length = 2*CAPSULE_HALF_LINEAR;		// one side

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

	clearsurfaceimpacts();


	// allocate and zero the 2d vectors

	       numnodes.resize(3);			
	     rep_radial.resize(3);
	 rep_transverse.resize(3);
	    link_radial.resize(3);
	link_transverse.resize(3);
	    disp_radial.resize(3);
	disp_transverse.resize(3);

	for (int axis = 0; axis < 3; ++axis)
	{
		       numnodes[axis].resize(num_segs);			
		     rep_radial[axis].resize(num_segs);
		 rep_transverse[axis].resize(num_segs);
		    link_radial[axis].resize(num_segs);
		link_transverse[axis].resize(num_segs);
		    disp_radial[axis].resize(num_segs);
		disp_transverse[axis].resize(num_segs);

		for(int i = 0; i<num_segs; ++i)
		{
			       numnodes[axis][i].resize(num_dist_segs);
			     rep_radial[axis][i].resize(num_dist_segs);
			 rep_transverse[axis][i].resize(num_dist_segs);
			    link_radial[axis][i].resize(num_dist_segs);
			link_transverse[axis][i].resize(num_dist_segs);
			    disp_radial[axis][i].resize(num_dist_segs);
			disp_transverse[axis][i].resize(num_dist_segs);
		}
	}

	clearnodes();

	// define force lines co-ords

	linestartx.resize(3);
	linestarty.resize(3);

	lineunitvecx.resize(3);
	lineunitvecy.resize(3);

	for (int axis = 0; axis < 2; ++axis)
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

}

void segments::getsegmentnum(const nodes& node, int& xseg, int& yseg, int& zseg) const
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
	
	if ((y <= CAPSULE_HALF_LINEAR) && (y >= -CAPSULE_HALF_LINEAR))
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

		return num_cap_segs + num_straight_segs
			+ (int)((double)num_cap_segs * (1+ (atan2((y + (CAPSULE_HALF_LINEAR)),-x) / PI)) );
		
	}

}

void segments::getsegmentdist(const nodes& node,int& xdist, int& ydist, int& zdist) const
{
	
	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		xdist = dist_to_seg(calcdist(node.y, node.z) - RADIUS); 
		ydist = dist_to_seg(calcdist(node.x, node.z) - RADIUS); 
		zdist = dist_to_seg(calcdist(node.x, node.y) - RADIUS); 
		return;
	}
	else
	{	// capsule

		zdist = dist_to_seg(calcdist(node.x, node.y) - RADIUS); 

		if ((node.z <= CAPSULE_HALF_LINEAR) && (node.z >= -CAPSULE_HALF_LINEAR))
		{	// on cylinder

			xdist = dist_to_seg(fabs(node.y) - RADIUS); 
			ydist = dist_to_seg(fabs(node.x) - RADIUS); 
		
		} 
		else if (node.z > 0)
		{	// top cap
	
			xdist = dist_to_seg(calcdist(node.y, node.z - CAPSULE_HALF_LINEAR) - RADIUS); 
			ydist = dist_to_seg(calcdist(node.x, node.z - CAPSULE_HALF_LINEAR) - RADIUS); 

		}
		else	//  (y < CAPSULE_HALF_LINEAR)
		{	// bottom cap

			xdist = dist_to_seg(calcdist(node.y, node.z + CAPSULE_HALF_LINEAR) - RADIUS); 
			ydist = dist_to_seg(calcdist(node.x, node.z + CAPSULE_HALF_LINEAR) - RADIUS); 
			
		}

	}

}


int segments::dist_to_seg(const double dist) const
{	// convert distance to seg num, and return -1 if out of range
	int dist_seg = (int) (dist / dist_step);

	return ((dist_seg>num_dist_segs)||(dist_seg < 0))?(-1):(dist_seg);
}



void segments::addnode(const nodes& node)
{	// add the node's data to the segments
    
	int xseg, yseg, zseg;
	int xdist, ydist, zdist;

	getsegmentnum(node, xseg, yseg, zseg);
	getsegmentdist(node, xdist, ydist, zdist);

	if (xdist!=-1)
		numnodes[2][zseg][zdist] ++;

	if (ydist!=-1)	
		numnodes[1][yseg][ydist] ++;

	if (zdist!=-1)
		numnodes[0][xseg][xdist] ++;
	
	for (int threadnum = 0; threadnum < NUM_THREADS; ++threadnum)
	{
		if (xdist!=-1)
		{

			 rep_radial[0][xseg][xdist] += node.repforce_radial[threadnum];
		 rep_transverse[0][xseg][xdist] += node.repforce_transverse[threadnum];
			link_radial[0][xseg][xdist] += node.linkforce_radial[threadnum];
		link_transverse[0][xseg][xdist] += node.linkforce_transverse[threadnum];
			disp_radial[0][xseg][xdist] += node.dispforce_radial[threadnum];
		disp_transverse[0][xseg][xdist] += node.dispforce_transverse[threadnum];
		}


		if (ydist!=-1)
		{
			 rep_radial[1][yseg][ydist] += node.repforce_radial[threadnum];
		 rep_transverse[1][yseg][ydist] += node.repforce_transverse[threadnum];
			link_radial[1][yseg][ydist] += node.linkforce_radial[threadnum];
		link_transverse[1][yseg][ydist] += node.linkforce_transverse[threadnum];
			disp_radial[1][yseg][ydist] += node.dispforce_radial[threadnum];
		disp_transverse[1][yseg][ydist] += node.dispforce_transverse[threadnum];
		}


		if (zdist!=-1)
		{
			 rep_radial[2][zseg][zdist] += node.repforce_radial[threadnum];
		 rep_transverse[2][zseg][zdist] += node.repforce_transverse[threadnum];
			link_radial[2][zseg][zdist] += node.linkforce_radial[threadnum];
		link_transverse[2][zseg][zdist] += node.linkforce_transverse[threadnum];
			disp_radial[2][zseg][zdist] += node.dispforce_radial[threadnum];
		disp_transverse[2][zseg][zdist] += node.dispforce_transverse[threadnum];
		}


	}
}

void segments::addsurfaceimpact(const nodes& node, const double& mag)
{	// add a surface impact
	
    int xseg,yseg,zseg;

	getsegmentnum(node, xseg, yseg, zseg);

    surfaceimpacts[0][xseg]+= mag * fabs(RADIUS - node.x);
	surfaceimpacts[1][yseg]+= mag * fabs(RADIUS - node.y);

	if (p_nuc->geometry == nucleator::sphere)
	{
		surfaceimpacts[2][zseg]+= mag * fabs(1 - node.z);	
	}
	else
	{
		if ((node.z < CAPSULE_HALF_LINEAR) && (node.z > -CAPSULE_HALF_LINEAR))
		{  // on cylinder, don't scale by z component
			surfaceimpacts[2][zseg]+= mag;
		}
		else
		{	// on ends
			if (node.z > 0)
			{
				surfaceimpacts[2][zseg]+= mag * fabs(RADIUS - (node.z-CAPSULE_HALF_LINEAR));	
			}
			else
			{
				surfaceimpacts[2][zseg]+= mag * fabs(RADIUS - (node.z+CAPSULE_HALF_LINEAR));
			}

		}

	}

}

void segments::clearsurfaceimpacts(void)
{
	for (int axis = 0; axis < 3; ++axis)
	{
		for(int i = 0; i<num_segs; ++i)
		{
			surfaceimpacts[axis][i]=0;
		}
	}
}

void segments::clearnodes(void)
{
	for (int axis = 0; axis < 3; ++axis)
	{
		for(int i = 0; i<num_segs; ++i)
		{
			for(int j = 0; j<num_dist_segs; ++j)
			{
				       numnodes[axis][i][j]=0;
				     rep_radial[axis][i][j]=0;
				 rep_transverse[axis][i][j]=0;
				    link_radial[axis][i][j]=0;
				link_transverse[axis][i][j]=0;
				    disp_radial[axis][i][j]=0;
			    disp_transverse[axis][i][j]=0;
			}
		}
	}
}

void segments::drawoutline(ostream& drawcmd, const int& axis) const
{

	int centerx, centery;
    int radius;
    int segment;

	centerx = BMP_WIDTH / 8;
	centery = BMP_HEIGHT / 2;
	
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

void segments::drawsurfaceimpacts(ostream& drawcmd, const int& axis, const double scale) const
{

	int centerx, centery;
    //int radius;
    //int segment;

	centerx = BMP_WIDTH / 8;
	centery = BMP_HEIGHT / 2;
	
	double linelen, linex, liney;
	double startx, starty;

	double seg_area, unscaledlen;

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
	
		unscaledlen = surfaceimpacts[axis][i] * FORCEBAR_SCALE * scale / seg_area;

		linex = unscaledlen * lineunitvecx[axis][i];
		liney = unscaledlen * lineunitvecy[axis][i];

		linelen = calcdist(linex,liney);

		if (linelen > 0.9 * RADIUS)
		{
			linex*= 0.9 * RADIUS / linelen;
			liney*= 0.9 * RADIUS / linelen;			
		}

		drawcmd << " line "
				<< centerx + p_actin->pixels(startx) << "," 
				<< centery + p_actin->pixels(starty) << " "

				<< centerx + p_actin->pixels(startx + linex) << ","
				<< centery + p_actin->pixels(starty + liney);

	}

}

void segments::addallnodes()
{
	for (int i=0; i<p_actin->highestnodecount; i++)
	{
		if ((p_actin->node[i].polymer) && (!p_actin->node[i].harbinger))  // is point valid?
		{
			addnode(p_actin->node[i]);
		}
	}

}

void segments::drawnodestats(ostream& drawcmd, const int& axis) const
{

}

void segments::savereport(const int& filenum) const
{

	char filename[255];

	double x,y,z;
	double radius;

	sprintf ( filename , "report%05i.txt", filenum );

	ofstream opreport(filename, ios::out | ios::trunc);
	if (!opreport) 
	{ cout << "Unable to open file " << filename << " for output"; return;}

	// write header
	
	opreport << "Axis,segment,distseg,x,y,z,NodeRadius,RepForceRadial,RepForceTrans,RepDisplRadial,RepDisplTrans,LinkForceRadial,LinkForceTrans" << endl;
	
	for (int axis = 0; axis < 3; ++axis)
	{
		for(int seg = 0; seg<num_segs; ++seg)
		{
			if ((axis == 2) && (seg >= 2 * num_cap_segs))
				continue;  // no linear segment on z axis

			if (axis == 0)
			{
                x = 0;
				y = linestartx[axis][seg];
				z = linestarty[axis][seg];
			}
			else if (axis == 1)
			{
				x = linestartx[axis][seg];
				y = 0;
				z = linestarty[axis][seg];
			}
			else
			{
				x = linestartx[axis][seg];
				y = linestarty[axis][seg];
				z = 0;
			}


			if ((p_nuc->geometry == nucleator::sphere) || (axis == 2))
			{	// sphere or capsule z
				radius = calcdist(x,y) - RADIUS;
			}
			else
			{
				if (fabs(y) < CAPSULE_HALF_LINEAR)
					radius = fabs(x) - RADIUS;
				else
					radius = calcdist(x,fabs(y)-CAPSULE_HALF_LINEAR) - RADIUS;
			}

			for(int dist = 0; dist<num_dist_segs; ++dist)
			{
				opreport 
					<< axis << ","
					<< seg << ","
					<< dist << ","
					<< x << ","
					<< y << ","
					<< z << ","
					<< radius << ",";

					 

			}
		}
	}




	opreport << endl;

	opreport.close();
}
