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

#ifdef LINK_VTK

/// win32 vtk 4.2 libaries need floats not doubles
#ifdef _WIN32
  #define VTK_FLOAT_PRECISION float
#else
  #define VTK_FLOAT_PRECISION double
#endif

#ifndef NOKBHIT
#ifndef _WIN32
	#include "kbhit.h"
	#define kbhit keyb.kbhit
    #define getch keyb.getch
#else
	#include <conio.h>
#endif
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "comet_vtk.h"

// vtk
#include "vtkCubeSource.h"
//
#include "vtkMath.h"
#include "vtkProperty.h"
#include "vtkSphereSource.h"
#include "vtkCylinderSource.h"
#include "vtkLineSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkGaussianSplatter.h"
#include "vtkImageData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkTubeFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkAxes.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkContourFilter.h"
#include "vtkLight.h"
#include "vtkCamera.h"
#include "vtkLookupTable.h"
#include "vtkInteractorObserver.h"
#include "vtkImageData.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkIdType.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"


#include "vtkLightCollection.h"
// vol rend
#include "vtkImageImport.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMIPFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkImageCast.h"
#include "vtkShepardMethod.h"

#include "vtkSmoothPolyDataFilter.h"

// output
#include "vtkDataWriter.h"
#include "vtkStructuredPointsWriter.h"

// ImageWriter
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"

using namespace std; // REVISIT only iostream probably temporary -- remove later

// Simple visualisation uing the current 'actin' simulation object to get at the results
// Makes extensive use of the VTK visulisation toolkit (public.kitware.com/VTK/) 
// tested with v4.4 (current stable Dec 2005), see examples for some of the uses below.
// This code is called by 'post' option, and reads 'cometparams.ini' to set options for the
// visualisation.  Some of this is a bit clumsy, as inconsistent options can be set, by design
// whatever is requested will be added to the visualisation, it's up to the user to set 
// sensible and compatible options.

CometVtkVis::CometVtkVis(actin * theactin)
{
    p_actin = theactin;

    voxel_scale = 4.0;
    p_scale = 100;
    nuc_opacity = 0.4;
    vx_intensity_scale = 1.0;
    file_prefix = "vtk";
        
    // find out data size
    ni = int(BMP_HEIGHT / voxel_scale);
    nj = int(BMP_HEIGHT / voxel_scale);
    nk = int(BMP_HEIGHT / voxel_scale);
    
    OptsInteractive         = true;
    OptsRenderNucleator     = true;
    OptsRenderNodes         = false;
    OptsRenderLinks         = false;
    OptsShadeLinks          = false;
    OptsVolumeRenderNodes   = false;
    OptsIsoRenderNodes      = true;
    OptsRenderAxes          = false;
    OptsRenderText          = false;
    OptsNormaliseFrames     = false;
    OptsRenderProjection    = RIP;
    renderwin_npx = 0;
    renderwin_npy = 0;

    setOptions();

	radius_pixels = p_actin->dbl_pixels(0.99 * RADIUS)/voxel_scale;
    
    if(renderwin_npx == 0 || renderwin_npy == 0) {
	if(OptsInteractive) {
	    renderwin_npx = VTK_WIDTH;
	    renderwin_npy = VTK_HEIGHT;
	} else {
	    renderwin_npx = VTK_WIDTH * VTK_AA_FACTOR;
	    renderwin_npy = VTK_HEIGHT * VTK_AA_FACTOR;
	}
    }
    cout << "  render win (nx, ny) = " << renderwin_npx << " " << renderwin_npy << endl;
    cout << "  voxel    (ni nj nk) = " << ni<<" " << nj << " " << nk << endl;

    // REVISIT: Make a renderer each time? or clear out and rebuild?
    // a renderer and render window
    //renderer = vtkRenderer::New();
    //render_win = vtkRenderWindow::New();
    //render_win->AddRenderer(renderer);	
}

CometVtkVis::~CometVtkVis()
{
    // finished
    // REVISIT: deleting renderer kills the pipeline?
    // iren->Delete();
    //render_win->Delete();
    //renderer->Delete();
}

void CometVtkVis::buildVTK(int framenumber)
{
    char filename[255];
    sprintf(filename , "%s%s_%05i.png", VTKDIR, file_prefix.c_str(), framenumber);
  
    renderer = vtkRenderer::New();
    render_win = vtkRenderWindow::New();
    renderer->SetBackground(0, 0, 0);		

	render_win->LineSmoothingOn();
	render_win->PointSmoothingOn();
	//render_win->PolygonSmoothingOn();

	if(!OptsInteractive)	 // increase quality for non-interactive
	{
		render_win->SetAAFrames(10);
	}

	//render_win->SetStereoRender(2);


    render_win->SetSize(renderwin_npx, renderwin_npy);  
    if(!OptsInteractive) {
	render_win->OffScreenRenderingOn();
    }												
    render_win->AddRenderer(renderer);    
    
    // set projection early
    setProjection(); 
    
    // add objects to renderer
    if(OptsRenderNucleator)
	addNucleator();
  
    if(OptsRenderNodes)
	addNodes();
  
    if(OptsRenderLinks)        
	addLinks();
  
    if(OptsVolumeRenderNodes)
	addGuassianNodeVolume(false);        
    
    if(OptsIsoRenderNodes)
	addGuassianNodeVolume(true);
    
    if(OptsRenderAxes)
	addAxes();
    
    // if(OptsRenderText)        
    // addVoxelBound(renderer);    
    addLight();
    
    // -- rendering
    if(OptsInteractive) {

	render_win-> Render();

	// allow interaction
	iren = vtkRenderWindowInteractor::New();
	iren->SetRenderWindow(render_win);
	iren->Initialize();
	iren->Start();
	
	iren->Delete();
    } else 
	{
		render_win->Render();
		saveImage(filename);
		if (VTK_AA_FACTOR!=1)
		{
			char command1[255];
			sprintf(command1, 
 "%s -resize %f%%  -font helvetica -fill white -pointsize 20 -draw \"text +5+20 'Frame % 6i\\nTime % 6.1f'\" %s &",
				IMAGEMAGICKMOGRIFY, 100/(double)VTK_AA_FACTOR,  
			    framenumber, framenumber * InterRecordIterations * DELTA_T, filename);
			//cout << command1 << endl;
			system(command1);
		}
    }
    
    // clear all actors from the renderer
    // renderer->RemoveAllProps();
    // this reuse of renderer seems to cause segfaults in vtk.
    
    // CHECK: order matters in deletion of these objects
    // otherwise too many X Servers error.
    renderer->Delete();
    render_win->Delete();
}

// -- Helper functions
void CometVtkVis::addAxes()
{
    // see finance.cxx in vtk examples
    vtkAxes *axes = vtkAxes::New();
    VTK_FLOAT_PRECISION origin[3] = {0.0, 0.0, 0.0};
    axes->SetOrigin( origin );
    axes->SetScaleFactor(1.3*p_actin->dbl_pixels(RADIUS)/voxel_scale);
    
    vtkTubeFilter *axes_tubes = vtkTubeFilter::New();
    axes_tubes->SetInput( axes->GetOutput() );
    axes_tubes->SetRadius(axes->GetScaleFactor()/50.0);
    axes_tubes->SetNumberOfSides(10);
    
    vtkPolyDataMapper *axes_mapper = vtkPolyDataMapper::New();
    axes_mapper->SetInput(axes_tubes->GetOutput());
    vtkActor *axes_actor = vtkActor::New();
    axes_actor->SetMapper(axes_mapper);
    axes_mapper->Delete();
    
    renderer->AddActor(axes_actor);
    axes_actor->Delete();
}

void CometVtkVis::addVoxelBound()
{
    // Add a crude voxel bounding box, sometimes useful
    // create lines
    VTK_FLOAT_PRECISION pt[3];
    
    vtkLineSource *line1 = vtkLineSource::New();
    vtkPolyDataMapper *map1 = vtkPolyDataMapper::New();
    vtkActor *line_actor1 = vtkActor::New();
    pt[0] = -ni/2;
    pt[1] = -nj/2;
    pt[2] = -nk/2;
    line1->SetPoint1( pt );
    pt[0] = -ni/2;
    pt[1] = +nj/2;
    pt[2] = -nk/2;
    line1->SetPoint2( pt );

    // add the actor to the scene
    map1->SetInput(line1->GetOutput());  
    line_actor1->SetMapper(map1);
    line_actor1->GetProperty ()->SetColor(1, 1, 0);
    renderer->AddActor(line_actor1);
    line1->Delete();
    map1->Delete();
    line_actor1->Delete();

    vtkLineSource *line2 = vtkLineSource::New();
    vtkPolyDataMapper *map2 = vtkPolyDataMapper::New();
    vtkActor *line_actor2 = vtkActor::New();
    pt[0] = -ni/2;
    pt[1] = -nj/2;
    pt[2] = -nk/2;
    line2->SetPoint1( pt );
    pt[0] = -ni/2;
    pt[1] = -nj/2;
    pt[2] = +nk/2;
    line2->SetPoint2( pt );
    // add the actor to the scene
    map2->SetInput(line2->GetOutput());  
    line_actor2->SetMapper(map2);
    line_actor2->GetProperty ()->SetColor(0, 1, 0);
    renderer->AddActor(line_actor2);
    line2->Delete();
    map2->Delete();
    line_actor2->Delete();

    vtkLineSource *line3 = vtkLineSource::New();
    vtkPolyDataMapper *map3 = vtkPolyDataMapper::New();
    vtkActor *line_actor3 = vtkActor::New();
    pt[0] = -ni/2;
    pt[1] = -nj/2;
    pt[2] = -nk/2;
    line3->SetPoint1( pt );
    pt[0] = +ni/2;
    pt[1] = -nj/2;
    pt[2] = -nk/2;
    line3->SetPoint2( pt );

    // add the actor to the scene
    map3->SetInput(line3->GetOutput());  
    line_actor3->SetMapper(map3);
    line_actor3->GetProperty ()->SetColor(1, 0, 0);
    renderer->AddActor(line_actor3);
    line3->Delete();
    map3->Delete();
    line_actor3->Delete();
}

void CometVtkVis::saveImage(char* image_filename)
{
    cout << "  saving image to file: " << image_filename << endl;
  
    vtkWindowToImageFilter *rwin_to_image = vtkWindowToImageFilter::New();
    rwin_to_image->SetInput(render_win);
  
    vtkPNGWriter *png_writer = vtkPNGWriter::New();
    png_writer->SetInput( rwin_to_image->GetOutput() );
    png_writer->SetFileName( image_filename );
    png_writer->Write();    
    png_writer->Delete();
  
    rwin_to_image->Delete();
}

void CometVtkVis::addGuassianNodeVolume(bool do_iso)
{	
    vtkStructuredPoints *spoints = vtkStructuredPoints::New();
    createGaussianVoxelRepresentation(spoints);
  
    if(do_iso){
	// render the volume using an iso surface
	addStructuredPointIsoRender(spoints);    
    } else {
	// volume render the points
	addStructuredPointVolumeRender(spoints);
    }
    
    // Delete spoints
    spoints->Delete();
}

void CometVtkVis::addNucleator()
{
    if(p_actin->p_nuc->geometry == nucleator::sphere){
	addSphericalNucleator();
    } else {
	addCapsuleNucleator();
    }
}

void CometVtkVis::addCapsuleNucleator()
{
    // -- Sources
    // -- endcap source
    // create top sphere geometry
    vtkSphereSource *endcap = vtkSphereSource::New();
    endcap->SetRadius(0.95*RADIUS); // can't we get radius through nucleator:: ??
    endcap->SetThetaResolution(18);
    endcap->SetPhiResolution(18);
    endcap->SetStartPhi(0);
    endcap->SetEndPhi(90);
    // endcap mapper
    vtkPolyDataMapper *endcap_mapper = vtkPolyDataMapper::New();
    endcap_mapper->SetInput( endcap->GetOutput() );
    endcap->Delete();
    
    // - body source
    vtkCylinderSource *body = vtkCylinderSource::New();
    body->SetRadius(0.95*RADIUS); // why not through nucleator:: ??
    body->SetHeight(2*CAPSULE_HALF_LINEAR);
    body->SetResolution(18);
    body->CappingOff();
    // body mapper
    vtkPolyDataMapper *body_mapper = vtkPolyDataMapper::New();
    body_mapper->SetInput( body->GetOutput() );
    body->Delete();

    // -- Actors
    // rotate the nucleator sections
    double nrotation[3];
    p_actin->p_nuc->nucleator_rotation.getangles(nrotation[0], 
						 nrotation[1], 
						 nrotation[2]);
    nrotation[0] *= vtkMath::DoubleRadiansToDegrees();
    nrotation[1] *= vtkMath::DoubleRadiansToDegrees();
    nrotation[2] *= vtkMath::DoubleRadiansToDegrees();
    
    // - upper endcap
    // actor coordinates geometry, properties, transformation
    vtkActor *endcap1_actor = vtkActor::New();
    endcap1_actor->SetMapper(endcap_mapper);
    VTK_FLOAT_PRECISION centre[3];
    centre[0] = 0;
    centre[1] = 0;
    centre[2] = +CAPSULE_HALF_LINEAR;
    endcap1_actor->SetOrientation( nrotation[0],
				   nrotation[1],
				   nrotation[2]);
    endcap1_actor->SetPosition( centre );
    endcap1_actor->GetProperty ()->SetColor(0, 0, 1);	// sphere color blue
    endcap1_actor->GetProperty()->SetOpacity(0.4);
    
    // - lower endcap
    vtkActor *endcap2_actor = vtkActor::New();
    endcap2_actor->SetMapper(endcap_mapper);
    centre[0] = 0;
    centre[1] = 0;
    centre[2] = -CAPSULE_HALF_LINEAR;
    // endcap2_actor->RotateX( 180 );
    endcap2_actor->SetOrientation( nrotation[0] + 180,
				   nrotation[1],
				   nrotation[2]);
    endcap2_actor->SetPosition( centre );
    endcap2_actor->GetProperty ()->SetColor(0, 0, 1);	// sphere color blue
    endcap2_actor->GetProperty()->SetOpacity(nuc_opacity);
    
    endcap_mapper->Delete();
    
    // - body
    vtkActor *body_actor = vtkActor::New();
    body_actor->SetMapper(body_mapper);
    centre[0] = 0;
    centre[1] = 0;
    centre[2] = 0;
    body_actor->SetPosition( centre );
    body_actor->SetOrientation( nrotation[0] + 90,
				nrotation[1],
				nrotation[2]);
    body_actor->GetProperty()->SetColor(0, 0, 1);	// sphere color blue
    body_actor->GetProperty()->SetOpacity(0.4);
    body_mapper->Delete();

    // add the actors to the scene
    renderer->AddActor(endcap1_actor);
    endcap1_actor->Delete();

    renderer->AddActor(endcap2_actor);
    endcap2_actor->Delete();

    renderer->AddActor(body_actor);
    body_actor->Delete();
}

void CometVtkVis::addSphericalNucleator()
{   
    // create sphere geometry
    vtkSphereSource *sphere = vtkSphereSource::New();
    
    //cout << "  voxel_scale: " << voxel_scale << endl;
    //cout << "  nucleator radius: " << radius_pixels << endl;
    sphere->SetRadius(radius_pixels);
    sphere->SetThetaResolution(64);
    sphere->SetPhiResolution(64);

    double nx, ny, nz;
    
    // temp: reset center:
	if (!VTK_MOVE_WITH_BEAD)
	{
		nx = -p_actin->p_nuc->position.x; 
		ny = -p_actin->p_nuc->position.y; 
		nz = -p_actin->p_nuc->position.z; 
	}
	else
	{
		nx = ny = nz = 0.0;

	}
    p_actin->p_nuc->nucleator_rotation.rotate(nx, ny, nz);
    vtk_cam_rot.rotate(nx, ny, nz); 

    // stops bead 
    double keep_within_border;
    if( p_actin->p_nuc->geometry == nucleator::sphere )
	keep_within_border = 2*RADIUS;
    else
	keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);

    int beadminx = int(p_actin->pixels(-keep_within_border - nx)/voxel_scale) + ni/2 + 1;
    int beadminy = int(p_actin->pixels(-keep_within_border - ny)/voxel_scale) + nj/2 + 1; 
    int beadminz = int(p_actin->pixels(-keep_within_border - nz)/voxel_scale) + nk/2 + 1; 
    int beadmaxx = int(p_actin->pixels( keep_within_border - nx)/voxel_scale) + ni/2 + 1; 
    int beadmaxy = int(p_actin->pixels( keep_within_border - ny)/voxel_scale) + nj/2 + 1; 
    int beadmaxz = int(p_actin->pixels( keep_within_border - nz)/voxel_scale) + nk/2 + 1; 
    
    int movex = 0;
    int movey = 0;
    int movez = 0;

    if(beadminx < 0)
	movex = -(beadminx);
    if(beadminy < 0)
	movey = -(beadminy);
    if(beadminz < 0)
	movez = -(beadminz);    

    if(beadmaxx > ni)
	movex = -(beadmaxx - ni);
    if(beadmaxy > nj)
	movey = -(beadmaxy - nj);
    if(beadmaxz > nk)
	movez = -(beadmaxz - nk);
    
    nx = p_actin->dbl_pixels(-nx)/voxel_scale; 
    ny = p_actin->dbl_pixels(-ny)/voxel_scale; 
    nz = p_actin->dbl_pixels(-nz)/voxel_scale;
    
    // displace to bring nucleator back in bounds
    nx += movex;
    ny += movey;
    nz += movez;

    // map
    vtkPolyDataMapper *map = vtkPolyDataMapper::New();
    map->SetInput(sphere->GetOutput());
    sphere->Delete();
    
    // actor coordinates geometry, properties, transformation
    vtkActor *nuc_actor = vtkActor::New();
    nuc_actor->SetMapper(map);
    nuc_actor->SetPosition(nx, ny, nz);
    map->Delete();

    nuc_actor->GetProperty()->SetColor(0.7, 0.7, 0.7); // sphere color 
    nuc_actor->GetProperty()->SetOpacity(nuc_opacity);
	nuc_actor->GetProperty()->SetDiffuse(0.25);
	nuc_actor->GetProperty()->SetSpecular(1.0);
	nuc_actor->GetProperty()->SetSpecularPower(5.0);

	
    // add the actor to the scene
    renderer->AddActor(nuc_actor);
    nuc_actor->Delete();
}

void CometVtkVis::fillVoxelSetFromActinNodes(vector< vector< vector<double > > >  &vx)
					    // double vd, double *min)
{
    // zero the voxel set
    for(int i=0; i<ni; ++i) {
	for(int j=0; j<nj; ++j) {
	    for(int k=0; k<nk; ++k) {		
		vx[i][j][k] = 0;
	    }
	}
    }
  
    double meanx, meany, meanz;
  
    // temp: reset center:
    meanx = -p_actin->p_nuc->position.x; 
    meany = -p_actin->p_nuc->position.y; 
    meanz = -p_actin->p_nuc->position.z; 
  
    p_actin->p_nuc->nucleator_rotation.rotate(meanx, meany, meanz);
    vtk_cam_rot.rotate(meanx, meany, meanz); 
  
    // move to static
    // Gaussian splat
    double gaussmax = (double) GAUSSFWHM * 3/2.0;  
    // full extent of gaussian radius -  fwhm is 2/3 this
  
    const int splat_sz = (int)(p_actin->dbl_pixels(gaussmax) / voxel_scale); // as in actin 
    cout << "  gaussian splat extent: "<< splat_sz << endl;
  
    vector< vector< vector<double > > >  splat;
    splat.resize(2*splat_sz+1);
    for(int i=0; i<2*splat_sz+1 ;i++) {
	splat[i].resize(2*splat_sz+1);
	for(int j=0; j<2*splat_sz+1; j++) {
	    splat[i][j].resize(2*splat_sz+1); // extent x2
	}
    }
  

    // create splat to reuse
    for(int i=-splat_sz; i<=splat_sz; i++) {
	for(int j=-splat_sz; j<=splat_sz; j++) {
	    for(int k=-splat_sz; k<=splat_sz; k++) {
	
		// cut the corners
		if(i*i + j*j + k*k > splat_sz*splat_sz )
		    continue;
	
		splat[i+splat_sz][j+splat_sz][k+splat_sz]
		    = exp(-3*((double)(i*i+j*j+k*k)) / 
			  (double)(splat_sz*splat_sz*splat_sz) );
	    }
	}
    }
  
    // stops bead moving out of visualisation
    double keep_within_border;
    if( p_actin->p_nuc->geometry == nucleator::sphere )
	keep_within_border = 2*RADIUS;
    else
	keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);
    
    int beadminx = int(p_actin->pixels(-keep_within_border - meanx)/voxel_scale) + ni/2 + 1;
    int beadminy = int(p_actin->pixels(-keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
    int beadminz = int(p_actin->pixels(-keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
    int beadmaxx = int(p_actin->pixels( keep_within_border - meanx)/voxel_scale) + ni/2 + 1; 
    int beadmaxy = int(p_actin->pixels( keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
    int beadmaxz = int(p_actin->pixels( keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
  
    int movex = 0;
    int movey = 0;
    int movez = 0;
  
    if(beadminx < 0)
	movex = -(beadminx);
    if(beadminy < 0)
	movey = -(beadminy);
    if(beadminz < 0)
	movez = -(beadminz);    
  
    if(beadmaxx > ni)
	movex = -(beadmaxx - ni);
    if(beadmaxy > nj)
	movey = -(beadmaxy - nj);
    if(beadmaxz > nk)
	movez = -(beadmaxz - nk);
  
    // loop over the nodes, and splat them into the volume
    vect node_pos;
    int x,y,z;
    int vi,vj,vk;
    double vx_imax = 0;
    double vx_isum = 0;
    int n_isums = 0;
    for(int n=0; n<p_actin->highestnodecount; n++) {
	if (!p_actin->node[n].polymer)
	    continue;
      
	node_pos = p_actin->node[n]; // copy x,y,z coords from node[n]
    
	p_actin->p_nuc->nucleator_rotation.rotate(node_pos); 
	vtk_cam_rot.rotate(node_pos); // bring rip to y-axis
    
	// node centre in local voxel coords (nuc at centre) 
	// REVISIT: check why we need to add one here
	x = int(p_actin->pixels(node_pos.x - meanx)/voxel_scale) + ni/2 + 1; 
	y = int(p_actin->pixels(node_pos.y - meany)/voxel_scale) + nj/2 + 1;
	z = int(p_actin->pixels(node_pos.z - meanz)/voxel_scale) + nk/2 + 1;
    
	// displace nucleator to bring back in bounds
	x += movex;
	y += movey;
	z += movez;
    
	if((x<0) || (x>=ni) ||
	   (y<0) || (y>=nj) ||
	   (z<0) || (z>=nk) )  // only plot if point in bounds
	{
	    //  cout << "point out of bounds " << x << "," << y << endl;
	} else {
	    // if((n%SPECKLE_FACTOR)==0)
      
	    for(int i=-splat_sz; i<=splat_sz; ++i) {
		for(int j=-splat_sz; j<=splat_sz; ++j) {
		    for(int k=-splat_sz; k<=splat_sz; ++k) {
	    
			if((i*i+j*j+k*k) > (splat_sz*splat_sz))
			    continue;  // don't do corners
	    
			vi = x+i;
			vj = y+j;
			vk = z+k;
	    
			if((vi<0) || (vi>=ni) ||
			   (vj<0) || (vj>=nj) ||
			   (vk<0) || (vk>=nk) )  // only plot if point in bounds		  
			    continue;
	    
			// amount of actin		  
			vx[vi][vj][vk] += 
			    splat[i+splat_sz][j+splat_sz][k+splat_sz];
	    
			// keep tabs on the maximum
			if(vx[vi][vj][vk] > vx_imax)
			    vx_imax = vx[vi][vj][vk];
			// calc sum
			vx_isum += vx[vi][vj][vk];
			++n_isums;
		    }
		}
	    } // loop over splat

	} // else (add intensity)
    
    } // loop over nodes
    cout << "  max gausian intensity: "  << vx_imax << endl;
    cout << "  mean gausian intensity: " << vx_isum/n_isums << endl;
  
    // if(p_actin->BMP_intensity_scaling){
    if(OptsNormaliseFrames){
    
	cout << "  normalising Frame, vx_max: " << vx_imax << endl;
	// scale the data
	for(int i=0; i<ni; ++i) {
	    for(int j=0; j<nj; ++j) {
		for(int k=0; k<nk; ++k) {		
		    vx[i][j][k] = vx[i][j][k]/vx_imax;
		}
	    }
	}
    
    } // NormaliseFrames 
  
}


void CometVtkVis::createGaussianVoxelRepresentation(vtkStructuredPoints *sp)
{
    // Accumulate a gaussian splat
  
    // find out extent of the data
    // initialise extents
  
    // create the voxel set as vtk structured point data
    // create a structured grid for the nodes, and sample them
    // simple binning of the node data into a voxel set, with
    // gaussian contributions.
  
    // create the voxel 3d vector
    vector< vector< vector<double > > >  vx;
    vx.resize(ni);  
    for(int i=0; i< ni; i++) {
	vx[i].resize(nj);
	for (int j=0; j<nj; j++) {
	    vx[i][j].resize(nk);
	}
    }
  
    // fill voxel set from the actin node data
    cout << "  filling vtk voxel data " << endl;
    const double vd = 1;
    VTK_FLOAT_PRECISION vx_orig[3] = {-ni/2.0-1, -nj/2.0-1, -nk/2.0-1};
    fillVoxelSetFromActinNodes(vx); //, vd, vx_orig);
    

    // set up the spoints for volume rendering
    int dim[3] = {ni, nj, nk};
    sp->SetDimensions(dim);
    sp->SetSpacing(vd, vd, vd);
    sp->SetOrigin( vx_orig ); 
  
    // create the underlying voxel data (only UnsignedChar or Short)
    // see vtk examples for voxel data
    vtkUnsignedCharArray *vxdata = vtkUnsignedCharArray::New();    
    double iscale = 1.0;
    if(OptsNormaliseFrames) {
	cout << "  normalised frame" << endl;
    } else {
	cout << "  using fixed intensity scale:" << vx_intensity_scale << endl;
	iscale = vx_intensity_scale;
    }
    
    // loop over voxel set converting to structured points suitable for 
    // vtk volume rendering
    int offset = 0;
    for(int k=0 ; k<nk ; ++k){
	for(int j=0 ; j<nj ; ++j){
	    for(int i=0 ; i<ni; ++i){
		
		double value = vx[i][j][k]/iscale; 
		if(value > 0) {
		    if(value > 1) {
			vxdata->InsertValue(offset, 255 );
		    } else {
			vxdata->InsertValue(offset, (int)(value*255) );
		    }
		} else {
		    vxdata->InsertValue(offset, 0);
		}
		// show voxel volume border, to aid alignment 
		//if(i==0 || j==0 || k==0) {
		//  data->InsertValue(offset, 255);	  
		//}		
		offset++;		
	    }
	}
    }
    sp->GetPointData()->SetScalars(vxdata);
    vxdata->Delete();  
}

void CometVtkVis::addStructuredPointVolumeRender(vtkStructuredPoints *vx)
{
    // Standard transfer functions
    vtkPiecewiseFunction *opacity_tf = vtkPiecewiseFunction::New();
    opacity_tf->AddPoint(0,  0.0);
    opacity_tf->AddPoint(255, 1.0);
    
    // green for volume
    vtkColorTransferFunction *colour_tf = vtkColorTransferFunction::New();
    colour_tf->AddRGBPoint(  0.0, 0.0, 0.0, 0.0); // black
    colour_tf->AddRGBPoint(255.0, 0.0, 1.0, 0.0); // green
    
    
    vtkVolumeProperty *volume_property = vtkVolumeProperty::New();
    volume_property->SetScalarOpacity(opacity_tf);
    volume_property->SetColor(colour_tf);
    volume_property->ShadeOn();
    volume_property->SetInterpolationTypeToLinear();
    
    vtkVolumeRayCastMIPFunction *mip = vtkVolumeRayCastMIPFunction::New();
    vtkVolumeRayCastMapper *volume_map = vtkVolumeRayCastMapper::New();
    volume_map->SetVolumeRayCastFunction(mip);
    // really important for small spacing
    // adjust to be comparable to vd
    volume_map->SetSampleDistance(0.1); // important
    
    // map the actual data
    // using Gaussian cast
    volume_map->SetInput( vx );
    
    vtkVolume *volume = vtkVolume::New();
    volume->SetMapper(volume_map);
    volume->SetProperty(volume_property);
    
    renderer->AddVolume(volume);
    
    volume->Delete();
    volume_map->Delete();
    mip->Delete();
    volume_property->Delete();
    colour_tf->Delete();
    opacity_tf->Delete();
}

void CometVtkVis::addStructuredPointIsoRender(vtkStructuredPoints *sp)
{
    // render using an iso-surface
    vtkContourFilter *iso_surf = vtkContourFilter::New();
    iso_surf->SetInput( sp );
    iso_surf->SetValue(0, 50);

	vtkSmoothPolyDataFilter *smooth = vtkSmoothPolyDataFilter::New();
	smooth->SetInput(iso_surf->GetOutput());
	smooth->SetNumberOfIterations( 100 );
	smooth->SetRelaxationFactor( 0.05 );

    vtkPolyDataMapper *iso_mapper = vtkPolyDataMapper::New();
    //iso_mapper->SetInput( iso_surf->GetOutput() );
	iso_mapper->SetInput( smooth->GetOutput() );
    iso_mapper->ScalarVisibilityOff();
    iso_surf->Delete();
    
    vtkActor *iso_actor = vtkActor::New();
    iso_actor->SetMapper(iso_mapper);
    iso_mapper->Delete();
    
    iso_actor->GetProperty()->SetOpacity(1.0);
    iso_actor->GetProperty()->SetColor(0.0, 1.0, 0.0);
    
    renderer->AddActor(iso_actor);
    iso_actor->Delete();
}

// ML:REVISIT I don't think this works
//
//
void CometVtkVis::addGaussianSplatterToNodes()
{
    // -- set up data (see vtk finance.cxx example)
    // data is on an unstructured grid
    
    vtkPoints *pts = vtkPoints::New();
    vtkFloatArray *scalars = vtkFloatArray::New();
    
    // Gaussian splatter accumulate mode...
    for(int i=0; i<p_actin->highestnodecount; i++)
    {
	pts->InsertPoint(i,
			 p_actin->node[i].x - p_actin->p_nuc->position.x,
			 p_actin->node[i].y - p_actin->p_nuc->position.y,
			 p_actin->node[i].z - p_actin->p_nuc->position.z);
	scalars->InsertValue(i, 100);
    }
    // put point and scalar data onto the grid
    vtkUnstructuredGrid *grid_data = vtkUnstructuredGrid::New();
    grid_data->SetPoints(pts);
    grid_data->GetPointData()->SetScalars(scalars);    
    pts->Delete();
    scalars->Delete();
    
    // gaussian splatter for nodes
    vtkGaussianSplatter *nodes_splatter = vtkGaussianSplatter::New();
    nodes_splatter->SetInput( grid_data );
    nodes_splatter->ScalarWarpingOff();
    nodes_splatter->SetSampleDimensions(100, 100, 100);
    nodes_splatter->SetRadius(0.05);
    
    // render using an iso-surface 
    vtkContourFilter *iso_surf = vtkContourFilter::New();
    iso_surf->SetInput( nodes_splatter->GetOutput() );
    iso_surf->SetValue(0, 64);
    
    vtkPolyDataMapper *iso_mapper = vtkPolyDataMapper::New();
    iso_mapper->SetInput( iso_surf->GetOutput() );
    iso_mapper->ScalarVisibilityOff();
    
    vtkActor *iso_actor = vtkActor::New();
    iso_actor->SetMapper( iso_mapper );
    iso_mapper->Delete();
    iso_actor->GetProperty()->SetOpacity(0.8);
    iso_actor->GetProperty()->SetColor(0.0, 1.0, 0.0);
    
    renderer->AddActor(iso_actor);
    // delete actor etc...
}

void CometVtkVis::addNodes()
{
    // represent nodes by spheres
    
    // create sphere geometry
    vtkSphereSource *sphere = vtkSphereSource::New();
    sphere->SetRadius( 0.05*p_actin->dbl_pixels(RADIUS)/voxel_scale ); // scale with the nucleator
    sphere->SetThetaResolution(10); // low res is fine
    sphere->SetPhiResolution(10);
    
    // map
    vtkPolyDataMapper *map = vtkPolyDataMapper::New();
    map->SetInput(sphere->GetOutput());
    sphere->Delete();

    // temp: reset center:
    double meanx, meany, meanz;
    meanx = -p_actin->p_nuc->position.x; 
    meany = -p_actin->p_nuc->position.y; 
    meanz = -p_actin->p_nuc->position.z; 
    
    p_actin->p_nuc->nucleator_rotation.rotate(meanx, meany, meanz);
    vtk_cam_rot.rotate(meanx, meany, meanz); 

    // stops bead 
    double keep_within_border;
    if( p_actin->p_nuc->geometry == nucleator::sphere )
	keep_within_border = 2*RADIUS;
    else
	keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);

    int beadminx = int(p_actin->pixels(-keep_within_border - meanx)/voxel_scale) + ni/2 + 1;
    int beadminy = int(p_actin->pixels(-keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
    int beadminz = int(p_actin->pixels(-keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
    int beadmaxx = int(p_actin->pixels( keep_within_border - meanx)/voxel_scale) + ni/2 + 1; 
    int beadmaxy = int(p_actin->pixels( keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
    int beadmaxz = int(p_actin->pixels( keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
    
    int movex = 0;
    int movey = 0;
    int movez = 0;

    if(beadminx < 0)
	movex = -(beadminx);
    if(beadminy < 0)
	movey = -(beadminy);
    if(beadminz < 0)
	movez = -(beadminz);    

    if(beadmaxx > ni)
	movex = -(beadmaxx - ni);
    if(beadmaxy > nj)
	movey = -(beadmaxy - nj);
    if(beadmaxz > nk)
	movez = -(beadmaxz - nk);

    // loop over the nodes
    double ncx, ncy, ncz;
    for(int i=0; i<p_actin->highestnodecount; i++)
    {
	ncx = p_actin->node[i].x;
	ncy = p_actin->node[i].y;
	ncz = p_actin->node[i].z;

	p_actin->p_nuc->nucleator_rotation.rotate(ncx, ncy, ncz); 
	vtk_cam_rot.rotate(ncx, ncy, ncz); // bring rip to y-axis
	
	ncx = p_actin->dbl_pixels(ncx - meanx)/voxel_scale; 
	ncy = p_actin->dbl_pixels(ncy - meany)/voxel_scale;
	ncz = p_actin->dbl_pixels(ncz - meanz)/voxel_scale;
	
	// displace to bring node back in bounds
	ncx += movex;
	ncy += movey;
	ncz += movez;
	
	// actor coordinates geometry, properties, transformation
	vtkActor *node_actor = vtkActor::New();
	node_actor->SetMapper(map);
	node_actor->SetPosition(ncx, ncy, ncz);
	
	if( p_actin->node[i].polymer ) { 
	    node_actor->GetProperty()->SetColor(1.0, 0, 0);     // plot polymerised red
	} else if( p_actin->node[i].harbinger ){
	    node_actor->GetProperty()->SetColor(0.5, 0.5, 0.0); // plot harbinger purple
	} else if( !p_actin->node[i].listoflinks.empty() ){
	    node_actor->GetProperty()->SetColor(0.5, 0.5, 0.5); // plot empty grey
	} else {
	    node_actor->GetProperty()->SetColor(0, 1, 0);       // otherwise green
	}	
	// add the actor to the scene	
	renderer->AddActor(node_actor);
	node_actor->Delete();
    }
    
    map->Delete();
}

void CometVtkVis::addLinks()
{	
  // Move to fcn
  // temp: reset center:
  double meanx = 0.0, meany=0.0, meanz=0.0;
  
  if(!VTK_MOVE_WITH_BEAD) {
    meanx = -p_actin->p_nuc->position.x; 
    meany = -p_actin->p_nuc->position.y; 
    meanz = -p_actin->p_nuc->position.z;
  }
  
  p_actin->p_nuc->nucleator_rotation.rotate(meanx, meany, meanz);
  vtk_cam_rot.rotate(meanx, meany, meanz); 
  
  // stops bead 
  double keep_within_border;
  if( p_actin->p_nuc->geometry == nucleator::sphere )
    keep_within_border = 2*RADIUS;
  else
    keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);
  
  int beadminx = int(p_actin->pixels(-keep_within_border - meanx)/voxel_scale) + ni/2 + 1;
  int beadminy = int(p_actin->pixels(-keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
  int beadminz = int(p_actin->pixels(-keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
  int beadmaxx = int(p_actin->pixels( keep_within_border - meanx)/voxel_scale) + ni/2 + 1; 
  int beadmaxy = int(p_actin->pixels( keep_within_border - meany)/voxel_scale) + nj/2 + 1; 
  int beadmaxz = int(p_actin->pixels( keep_within_border - meanz)/voxel_scale) + nk/2 + 1; 
  
  int movex = 0;
  int movey = 0;
  int movez = 0;
  
  if(beadminx < 0)
    movex = -(beadminx);
  if(beadminy < 0)
    movey = -(beadminy);
  if(beadminz < 0)
    movez = -(beadminz);    
  
  if(beadmaxx > ni)
    movex = -(beadmaxx - ni);
  if(beadmaxy > nj)
    movey = -(beadmaxy - nj);
  if(beadmaxz > nk)
    movez = -(beadmaxz - nk);
  // -- ^^ Move to fcn
  
  // only used for strain coloring
  vtkLookupTable *lut = vtkLookupTable::New();
  lut->SetRange(0.7, 1.0);
  lut->SetNumberOfColors(100);
  lut->Build();
  // ML REVISIT, crude linear ramp, 
  // do this better way to do this via vtkLut calls
  
  Colour col;
  
  VTK_FLOAT_PRECISION rgba[4];
  for(int i=0; i<100; i++) {
    col.setcol((double)i/100.0);
    
    rgba[0] = col.r;
    rgba[1] = col.g;
    rgba[2] = col.b;
    rgba[3] = 1.0;	
    lut->SetTableValue(i, rgba);
  }
  
  // loop over the nodes
  VTK_FLOAT_PRECISION n_pt[3];
  VTK_FLOAT_PRECISION l_pt[3];
  vect nodeposvec;
  
  // use a cell array to store line data 
  vtkCellArray *cellarray = vtkCellArray::New();
  // with scalars
  vtkFloatArray *cellscalars = vtkFloatArray::New();  
  // made up of lines
  vtkPolyData *linedata = vtkPolyData::New();
  vtkPoints *linepts = vtkPoints::New();
  
  for(int i=0; i<p_actin->highestnodecount; i++) {      
    nodeposvec = p_actin->node[i];
    
    n_pt[0] = p_actin->node[i].x;
    n_pt[1] = p_actin->node[i].y;
    n_pt[2] = p_actin->node[i].z;		
    
    // -- Move to fcn
    p_actin->p_nuc->nucleator_rotation.rotate(n_pt[0], n_pt[1], n_pt[2]); 
    vtk_cam_rot.rotate(n_pt[0], n_pt[1], n_pt[2]); // bring rip to y-axis
    
    n_pt[0] = p_actin->dbl_pixels(n_pt[0] - meanx)/voxel_scale; 
    n_pt[1] = p_actin->dbl_pixels(n_pt[1] - meany)/voxel_scale;
    n_pt[2] = p_actin->dbl_pixels(n_pt[2] - meanz)/voxel_scale;
    
    // displace to bring node back in bounds
    n_pt[0] += movex;
    n_pt[1] += movey;
    n_pt[2] += movez;
    // -- ^^ Move to fcn
    
    if(!p_actin->node[i].listoflinks.empty()) {
      // nodes thisnode = p_actin->node[i];
      
      // loop over linked nodes
      // check this out if we use this function
      for(vector<links>::iterator link_i=p_actin->node[i].listoflinks.begin(); 
	  link_i != p_actin->node[i].listoflinks.end();
	  ++link_i) {
	
	// assume links to nodes 'less' than us have already been added
	// cout << i << "==" << p_actin->node[i].nodenum << endl;
	if( (*link_i).linkednodenumber < i)
	  continue;
	
	// create line for the link
	l_pt[0] = (*link_i).linkednodeptr->x;
	l_pt[1] = (*link_i).linkednodeptr->y;
	l_pt[2] = (*link_i).linkednodeptr->z;
	
	// -- Move to fcn
	p_actin->p_nuc->nucleator_rotation.rotate(l_pt[0], l_pt[1], l_pt[2]); 
	vtk_cam_rot.rotate(l_pt[0], l_pt[1], l_pt[2]); // bring rip to y-axis
	
	l_pt[0] = p_actin->dbl_pixels(l_pt[0] - meanx)/voxel_scale; 
	l_pt[1] = p_actin->dbl_pixels(l_pt[1] - meany)/voxel_scale;
	l_pt[2] = p_actin->dbl_pixels(l_pt[2] - meanz)/voxel_scale;
	
	// displace to bring node back in bounds
	l_pt[0] += movex;
	l_pt[1] += movey;
	l_pt[2] += movez;
	// -- ^^ Move to fcn

	vtkIdType linept_ids[2];	
	linept_ids[0] = linepts->InsertNextPoint( n_pt );
	linept_ids[1] = linepts->InsertNextPoint( l_pt );
	vtkIdType cell_id = cellarray->InsertNextCell(2, linept_ids);	
	
	Colour col;	  
	// Set scalar for the line
	if(OptsShadeLinks) {
	  vect displacement = nodeposvec - *(link_i->linkednodeptr);
	  double distance = displacement.length();      
	  double strain = fabs(distance-link_i->orig_dist) / link_i->orig_dist;    
	  double y = strain / LINK_BREAKAGE_STRAIN;
	  //double y = fabs( link_i->getlinkforces( distance) ) / LINK_BREAKAGE_FORCE;
	  //y = (5+log(y))/5;	 // log transform	    
	  //double VTK_LINK_COLOUR_GAMMA = 1.8;	    
	  y = pow( y , 1/VTK_LINK_COLOUR_GAMMA);	    
	  y = y*0.9+0.1; // prevent zeros  because colorscheme makes them black
	  col.setcol(y);
	  
	  cellscalars->InsertValue(cell_id, y);
	} else {
	  cellscalars->InsertValue(cell_id, 1.0);
	}	
      } // links loop
    }
  } // node loop
  linedata->SetPoints(linepts);
  linepts->Delete();
  linedata->SetLines(cellarray);
  cellarray->Delete();
  linedata->GetCellData()->SetScalars(cellscalars);
  
  // mapper
  vtkPolyDataMapper *map = vtkPolyDataMapper::New();
  map->SetLookupTable( lut );
  map->SetInput( linedata );
  linedata->Delete();
  // actor
  vtkActor *lines_actor = vtkActor::New();
  lines_actor->SetMapper(map);
  //render
  renderer->AddActor(lines_actor);
  
  lut->Delete();
  lines_actor->Delete();
}
// ML FIXME:  This is wrong, sort out a better method for mean strain
double CometVtkVis::getMeanNodeLinkForce(const int id) 
{
    // loop over the links  
    double strain = 0;
    int n_links = 0;
    double distance;
    vect displacement;
    vect nodeposvec = p_actin->node[id];
    for(vector<links>::iterator i=p_actin->node[id].listoflinks.begin(); 
	i<p_actin->node[id].listoflinks.end(); i++ ) {	 
	displacement = nodeposvec - *(i->linkednodeptr);
	distance = displacement.length();      
	strain = distance / i->orig_dist;
	++n_links;
    }
    if(n_links != 0)
	strain = strain/n_links;
  
    return strain; 
}

void CometVtkVis::addLight()
{
    //vtkLight *light = vtkLight::New();
    ////light->SetLightTypeToCameraLight();
    ////light->SetFocalPoint(0, 0, 0);
    //light->SetPosition(0, 0, 0);
    //renderer->AddLight(light);
    //light->Delete();

	// remove random lights that may have appeared from nowhere

	renderer->GetLights()->RemoveAllItems();

	//renderer->SetAmbient(0); // this causes segfault for some reason

 // vtkLightCollection *lights = renderer->GetLights();
 // int numlights = lights->  //GetNumberOfItems();
 //   vtkLight *light_to_remove;
 //   lights->InitTraversal();
 //   for(int i=0; i<numlights; ++i) 
 //   {
 //       light_to_remove = lights->GetNextItem();
 //	      renderer->RemoveLight(light_to_remove);
 //   }

	// add our own lights

	vtkLight *light = vtkLight::New();
	
	light->SetIntensity(0.7);
	light->SetColor(1,0.8,0.8);
    light->SetPosition(radius_pixels*5, -radius_pixels*5, -radius_pixels*5);
    renderer->AddLight(light);
    light->Delete();

	vtkLight *light2 = vtkLight::New();
	light2->SetIntensity(0.7);
	light2->SetColor(0.9,0.9,1.0);
    light2->SetPosition(radius_pixels*5, +radius_pixels*5, -radius_pixels*5);
    renderer->AddLight(light2);
    light2->Delete();

}

void CometVtkVis::setProjection()
{
    if(OptsRenderProjection==X){
		// x
		vtk_cam_rot = p_actin->camera_rotation;
		renderer->GetActiveCamera()->SetPosition(-radius_pixels*7, 0, 0);
		renderer->GetActiveCamera()->SetViewUp(0, 0, -1);
    } else if(OptsRenderProjection==Y){
	// y
		vtk_cam_rot = p_actin->camera_rotation;
		renderer->GetActiveCamera()->SetPosition(0, radius_pixels*7, 0);
		renderer->GetActiveCamera()->SetViewUp(0, 0, -1);
    } else if(OptsRenderProjection==Z){
	// z
		vtk_cam_rot = p_actin->camera_rotation;
		renderer->GetActiveCamera()->SetPosition(0, 0, -radius_pixels*7);
		renderer->GetActiveCamera()->SetViewUp(0, -1, 0);
    } else if(OptsRenderProjection==RIP){

		vtk_cam_rot = p_actin->camera_rotation2;	  // note using different rotation matrix!
	// rip
	//cout << "  projection: rip" << endl;
	renderer->GetActiveCamera()->SetPosition(radius_pixels*4, 0, -radius_pixels*7);
	renderer->GetActiveCamera()->SetViewUp(1, 0, 0);    
    } else {
	cout << "!ERROR: unknown projection:" << OptsRenderProjection << endl;
    }
    
    renderer->GetActiveCamera()->SetFocalPoint(0, 0, 0);
	renderer->GetActiveCamera()->ParallelProjectionOff(); // ParallelProjectionOn();
	renderer->GetActiveCamera()->SetViewAngle(VTK_VIEWANGLE);
    // FIXME: ML
    // should scale properly here to a value linked to the render setup
    renderer->GetActiveCamera()->SetParallelScale(p_scale);
}

void CometVtkVis::setOptions()
{
    ifstream param("cometparams.ini"); 
    if(!param) 
    {
	cerr << "Cannot open cometparams.ini" << endl << endl;
	exit(EXIT_FAILURE);
    }
    
    string buffer;
    while (getline(param, buffer)) 
    { 
	istringstream ss(buffer);
	string tag, value;
	ss >> tag >> std::ws;
	if(tag.size() == 0 || tag[0] == '#')
	    // skip empty line or comment
	    continue;
      
	if(tag == "VIS_INTERACTIVE") 
	{ 
	    ss >> value;
	    OptsInteractive = getBoolOpt(value);
	    continue;
	}
      
	if(tag == "VIS_NUCLEATOR") 
	{
	    ss >> value;
	    OptsRenderNucleator = getBoolOpt(value);
	    continue;
	}
      
	if(tag == "VIS_NODES") 
	{
	    ss >> value;
	    OptsRenderNodes = getBoolOpt(value);
	    continue;	  
	}
	if(tag == "VIS_LINKS") 
	{
	    ss >> value;
	    OptsRenderLinks = getBoolOpt(value);
	    continue;
	}
	if(tag == "VIS_SHADELINKS") 
	{
	    ss >> value;
	    OptsShadeLinks = getBoolOpt(value);
	    continue;
	}      
	if(tag == "VIS_VOLUMENODES") 
	{
	    ss >> value;
	    OptsVolumeRenderNodes = getBoolOpt(value);
	    continue;
	}
	if(tag == "VIS_ISONODES") 
	{
	    ss >> value;
	    OptsIsoRenderNodes = getBoolOpt(value);
	    continue;
	}    
	if(tag == "VIS_AXES") 
	{
	    ss >> value;
	    OptsRenderAxes = getBoolOpt(value);
	    continue;
	}
	if(tag == "VIS_TEXT") 
	{
	    ss >> value;
	    OptsRenderText = getBoolOpt(value);
	    continue;
	}
	if(tag == "VIS_PROJECTION") 
	{
	    ss >> value;
	    OptsRenderProjection = getProjectionOpt(value);
	    continue;
	}
	if(tag == "VIS_VOXELSCALE") 
	{
	    ss >> voxel_scale;
	    continue;
	}
	if(tag == "VIS_NI-NJ-NK") 
	{
	    ss >> ni >> nj >> nk;
	    continue;
	}
	if(tag == "VIS_NUCOPACITY") 
	{
	    ss >> nuc_opacity;
	    continue;
	}
	if(tag == "VIS_VXISCALE") 
	{
	    ss >> vx_intensity_scale;
	    continue;
	}
	if(tag == "VIS_FILEPREFIX") 
	{
	    ss >> file_prefix;
	    continue;
	}
	if(tag == "VIS_NORMALISEFRAMES") 
	{
	    ss >> value;
	    OptsNormaliseFrames = getBoolOpt(value);
	    continue;
	}
	if(tag == "VIS_PSCALE") 
	{
	    ss >> p_scale;
	    continue;
	}

    }

	
  
}

bool CometVtkVis::getBoolOpt(const string &value) {
    // if(value=="true" || value == "TRUE" || value == "T" || value = "1")
    if(value.compare("true")==0)
	return true;
    else
	return false;
}

CometVtkVis::OptProjectionType CometVtkVis::getProjectionOpt(const string &value)
{
    if(value.compare("x")==0)
	return X;
    else if(value.compare("y")==0)
	return Y;
    else if(value.compare("z")==0)
	return Z;
    else if(value.compare("rip")==0)
	return RIP;
    else {
	cerr << "unknown projection type:" << value << endl;
	return NONE;
    }
}

void CometVtkVis::reportOptions()
{   
    cout << "OptsInteractive       = " << OptsInteractive       << endl;
    cout << "OptsRenderNucleator   = " << OptsRenderNucleator   << endl;
    cout << "OptsRenderNodes       = " << OptsRenderNodes       << endl;
    cout << "OptsRenderLinks       = " << OptsRenderLinks       << endl;
    cout << "OptsShadeLinks        = " << OptsShadeLinks        << endl;
    cout << "OptsVolumeRenderNodes = " << OptsVolumeRenderNodes << endl;
    cout << "OptsIsoRenderNodes    = " << OptsIsoRenderNodes    << endl;
    cout << "OptsRenderAxes        = " << OptsRenderAxes        << endl;
    cout << "OptsRenderText        = " << OptsRenderText        << endl;
    cout << "OptsRenderProjection  = " << OptsRenderProjection  << endl;
    cout << "OptsNormaliseFrames   = " << OptsNormaliseFrames   << endl;
    cout << "voxel_scale           = " << voxel_scale           << endl;
    cout << "ni nj nk              = " << ni<<" "<<nj<<" "<< nk << endl;
    cout << "vx intensity scale    = " << vx_intensity_scale    << endl;
    cout << "file_prefix           = " << file_prefix           << endl;
    cout << "projection_scale      = " << p_scale               << endl;
}


#endif

