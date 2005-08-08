

#include "stdafx.h"
#include "nucleator.h"
#include "segments.h"


//segments::segments(nucleator::shape set_geometry) : geometry (set_geometry)
//{	// the constructor is really setupsegments
//	// it can't be done here because we need the nucleator pointer
//  // and I can't get it to work in the constructor (not possible?)
//}

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
	cap_seg_len = curved_length / (MYDOUBLE) num_cap_segs;



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

	// define forces lines co-ords

	linestartx.resize(num_segs);
	linestarty.resize(num_segs);

	lineunitvecx.resize(num_segs);
	lineunitvecy.resize(num_segs);

	MYDOUBLE startx,starty, theta;
	MYDOUBLE dist;

	// top cap
	
	for(int i = 0; i< num_cap_segs; ++i)
	{
		//       offset     segnum                      angle of 1 seg
		theta = (0.5 + (MYDOUBLE) i) * (PI / (MYDOUBLE) num_cap_segs);

		startx = - cos (theta);  // unit components
		starty =   sin (theta);

		linestartx[i] = startx * RADIUS;	// line start
		linestarty[i] = starty * RADIUS;
        
		lineunitvecx[i] = - startx * FORCEBAR_SCALE;
		lineunitvecy[i] = - starty * FORCEBAR_SCALE;

		if (p_nuc->geometry == nucleator::capsule)
		{
			linestarty[i] += CAPSULE_HALF_LINEAR;
		}
	}

	// bottom cap

	for(int i = num_cap_segs; i< (2 *num_cap_segs); ++i)
	{
		//       offset     segnum                      angle of 1 seg
		theta = (0.5 + (MYDOUBLE) i) * (PI / (MYDOUBLE) num_cap_segs);

		startx = - cos (theta);  // unit components
		starty =   sin (theta);

		linestartx[i+num_straight_segs] = startx * RADIUS;  // line start
		linestarty[i+num_straight_segs] = starty * RADIUS;
        
		lineunitvecx[i+num_straight_segs] = - startx * FORCEBAR_SCALE;
		lineunitvecy[i+num_straight_segs] = - starty * FORCEBAR_SCALE;

		if (p_nuc->geometry == nucleator::capsule)
		{
			linestarty[i+num_straight_segs] -= CAPSULE_HALF_LINEAR;
		}
	}

	// sides

	if (p_nuc->geometry == nucleator::capsule)
	{
		
		for(int i = 0; i< num_straight_segs; ++i)
		{

			dist = (0.5 + (MYDOUBLE) i) * straight_seg_len;

			startx = 1;  // unit components
			starty = CAPSULE_HALF_LINEAR - dist ;

			linestartx[i+num_cap_segs] = startx * RADIUS;
			linestarty[i+num_cap_segs] = starty;

			// N.B. this is scaled by (cap_seg_len/straight_seg_len) to compenate for 
			// difference in segments lengths
	       
			lineunitvecx[i+num_cap_segs] = - startx * FORCEBAR_SCALE * (cap_seg_len/straight_seg_len);
			lineunitvecy[i+num_cap_segs] = 0;
		}

		for(int i = 0; i< num_straight_segs; ++i)
		{

			dist = (0.5 + (MYDOUBLE) i) * straight_seg_len;

			startx = - 1;  // unit components
			starty = dist - CAPSULE_HALF_LINEAR ;

			linestartx[i+2*num_cap_segs+num_straight_segs] = startx * RADIUS;
			linestarty[i+2*num_cap_segs+num_straight_segs] = starty;

			// N.B. this is scaled by (cap_seg_len/straight_seg_len) to compenate for 
			// difference in segments lengths
	        
			lineunitvecx[i+2*num_cap_segs+num_straight_segs] = - startx * FORCEBAR_SCALE * (cap_seg_len/straight_seg_len);
			lineunitvecy[i+2*num_cap_segs+num_straight_segs] = 0;
		}
		
	}

}

void segments::getsegmentnum(const nodes& node, int& xseg, int& yseg, int& zseg) const
{  // set the segment number for each axis, dependent on node position

	if (p_nuc->geometry == nucleator::sphere)
	{	// sphere

		xseg = (int)((MYDOUBLE)num_cap_segs * (1+ (atan2(-node.z,node.y) / PI)) );
		yseg = (int)((MYDOUBLE)num_cap_segs * (1+ (atan2(-node.z,node.x) / PI)) );
		zseg = (int)((MYDOUBLE)num_cap_segs * (1+ (atan2(-node.y,node.x) / PI)) );
		
		return;
	}
	else
	{	// capsule
	
		xseg = getcapsuleseg(node.y,node.z);
		yseg = getcapsuleseg(node.x,node.z);
		zseg = (int)((MYDOUBLE)num_cap_segs * (1+ (atan2(-node.y,node.x) / PI)) );
	}

}

int segments::getcapsuleseg(const MYDOUBLE & x, const MYDOUBLE & y) const
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

		return (int)((MYDOUBLE)num_cap_segs * (1+ (atan2(-(y - (CAPSULE_HALF_LINEAR)),x) / PI)) );

	}
	else	//  (y < CAPSULE_HALF_LINEAR)
	{	// bottom cap

		return num_cap_segs + num_straight_segs
			+ (int)((MYDOUBLE)num_cap_segs * (1+ (atan2((y + (CAPSULE_HALF_LINEAR)),-x) / PI)) );
		
	}

}

segments::segments(void)
{
}

segments::~segments(void)
{
}

void segments::addnode(const nodes& node)
{	// add the node's data to the segments

}

void segments::addsurfaceimpact(const nodes& node, const MYDOUBLE& mag)
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

void segments::drawsurfaceimpacts(ostream& drawcmd, const int& axis, const MYDOUBLE scale) const
{

	int centerx, centery;
    int radius;
    int segment;

	centerx = BMP_WIDTH / 8;
	centery = BMP_HEIGHT / 2;
	
	radius  = p_actin->pixels(RADIUS); 
	segment = p_actin->pixels(CAPSULE_HALF_LINEAR);

	MYDOUBLE linelen, linex, liney;
	MYDOUBLE startx, starty;

	int segstodraw;

	if ((p_nuc->geometry == nucleator::sphere) || (axis == 2))
	{
		segstodraw = 2 * num_cap_segs;
	}
	else
	{
		segstodraw = 2 * num_cap_segs + 2 * num_straight_segs;
	}

	for (int i=0; i<segstodraw; ++i)
	{

		startx = linestartx[i];
		starty = linestarty[i];
	
		linex = scale*surfaceimpacts[axis][i]*lineunitvecx[i];
		liney = scale*surfaceimpacts[axis][i]*lineunitvecy[i];

		if ((p_nuc->geometry == nucleator::capsule) && (axis == 2))
		{	// capsule z axis, do ends only and bring them to center
			if (i>=num_cap_segs)
			{
				starty = linestarty[i+num_straight_segs] + CAPSULE_HALF_LINEAR;
				startx = linestartx[i+num_straight_segs];

				linex = scale*surfaceimpacts[axis][i]*lineunitvecx[i+num_straight_segs];
				liney = scale*surfaceimpacts[axis][i]*lineunitvecy[i+num_straight_segs];
			}
			else
			{
				starty = linestarty[i] - CAPSULE_HALF_LINEAR;
			}
		}

		linelen = calcdist(linex,liney);

		if (linelen > 0.9)
		{
			linex*= 0.9/linelen;
			liney*= 0.9/linelen;			
		}

		linex*= RADIUS;
		liney*= RADIUS;

		drawcmd << " line "
				<< centerx + p_actin->pixels(startx) << "," 
				<< centery + p_actin->pixels(starty) << " "

				<< centerx + p_actin->pixels(startx + linex) << ","
				<< centery + p_actin->pixels(starty + liney);

	}

}
