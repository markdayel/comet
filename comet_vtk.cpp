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

#include "Colour.h"

#ifdef LINK_VTK


#define NUCTEXTUREJPEG1 "nuctex.jpg"
#define NUCTEXTUREJPEG2 "~/comet/nuctex.jpg"

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

#include "tracknodeinfo.h"

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
//#include "vtkIdType.h"
#include "vtkCellData.h"
#include "vtkDataSet.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

//#include "vtkSmartPointer.h"

#include "vtkCardinalSpline.h"

#include "vtkStructuredGrid.h"
#include "vtkHedgeHog.h"
#include "vtkStreamLine.h"
#include "vtkStreamTracer.h"
#include "vtkTubeFilter.h"
#include "vtkRungeKutta4.h"

#include "vtkRenderLargeImage.h"
#include "vtkAssembly.h"

#include "vtkExporter.h"

//#include "vtkPlanes.h"
#include "vtkPlane.h"
//#include "vtkPoints.h"
//#include "vtkNormals.h"

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
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkPolyDataNormals.h"

// output
#include "vtkDataWriter.h"
#include "vtkStructuredPointsWriter.h"

// ImageWriter
#include "vtkWindowToImageFilter.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkBMPWriter.h"
#include "vtkVRMLExporter.h"

// Texture map
#include "vtkTexture.h"
#include "vtkTextureMapToSphere.h"
#include "vtkTransformTextureCoords.h"
#include "vtkJPEGReader.h"

using namespace std; // REVISIT only iostream probably temporary -- remove later

// Simple visualisation uing the current 'actin' simulation object to get at the results
// Makes extensive use of the VTK visulisation toolkit (public.kitware.com/VTK/) 
// tested with v4.4 (current stable Dec 2005), see examples for some of the uses below.
// This code is called by 'post' option, and reads 'cometparams.ini' to set options for the
// visualisation.  Some of this is a bit clumsy, as inconsistent options can be set, by design
// whatever is requested will be added to the visualisation, it's up to the user to set 
// sensible and compatible options.

CometVtkVis::CometVtkVis(bool VIEW_VTK) // this parameter *should* be whether to create a viewing window
                                        // does it need to be true for off-screen as well?
{

    voxel_scale = 4.0;
    p_scale = 100;
    nuc_opacity = 0.4;
    vx_intensity_scale = 1.0;
    file_prefix = "vtk";
        
    // find out data size
    ni = int(BMP_HEIGHT/voxel_scale) ;
    nj = int(BMP_HEIGHT/voxel_scale) ;
    nk = int(BMP_HEIGHT/voxel_scale) ;
    
    //OptsInteractive         = true;
    OptsRenderNucleator     = true;
    OptsRenderNodes         = false;
    OptsRenderLinks         = false;
    OptsShadeLinks          = false;
    OptsVolumeRenderNodes   = false;
    OptsIsoRenderNodes      = true;
    OptsRenderAxes          = false;
    OptsRenderText          = false;
    OptsNormaliseFrames     = false;
    OptsSkipOutOfFocusPoints = true;
    OptsRenderProjection    = RIP;
    OptsUseNucTextureMap    = true;
    renderwin_npx = 0;
    renderwin_npy = 0;
    VIS_PARALLELPROJECTION = true;

    VTK_HIGHQUAL = false;

    OptsCameraDistMult = 7;
    VIS_LINETHICKNESS = 1.0;

    //voxelscalefactor = ptheactin->dbl_pixels(1.0)/voxel_scale;

    voxelscalefactor = (800 / 30) / voxel_scale;   // fix this---the voxel scale stuff is a mess (and unnecessary?)

    setOptions();

    // kludge for os x rendering
//#ifdef __APPLE__ 
//    VTK_AA_FACTOR = 1;
//#endif

    // this is *.99 so that links on the surface don't look like they go inside the sphere
    radius_pixels = voxelscalefactor * 0.99 * RADIUS;
    capsule_length_pixels = voxelscalefactor * CAPSULE_HALF_LINEAR;
    
    //if(renderwin_npx == 0 || renderwin_npy == 0) {
    //  if(POST_VTK_VIEW) {
	renderwin_npx = VTK_WIDTH;
	renderwin_npy = VTK_HEIGHT;
    //  } else {
	//renderwin_npx = VTK_WIDTH * VTK_AA_FACTOR;
	//renderwin_npy = VTK_HEIGHT * VTK_AA_FACTOR;
    //  }
    //}
    viewupvect=vect(0,0,-1);
    
    // REVISIT: Make a renderer each time? or clear out and rebuild?
    // a renderer and render window
    //renderer = vtkRenderer::New();
    //render_win = vtkRenderWindow::New();
    //render_win->AddRenderer(renderer);

    renderer = vtkRenderer::New();
	renderer->SetBackground(0, 0, 0);

    OptVIEW_VTK = VIEW_VTK;

    //cout << "  render win (nx, ny) = " << renderwin_npx << " " << renderwin_npy << endl;
    //cout << "  voxel    (ni nj nk) = " << ni<<" " << nj << " " << nk << endl;

    render_win = vtkRenderWindow::New();

    //if(!(POST_VTK_VIEW || VIEW_VTK))
    //    render_win->OffScreenRenderingOn();    // this seems to do nothing...

    render_win->SetSize(renderwin_npx, renderwin_npy);
    
	//render_win->LineSmoothingOn();
    //render_win->PointSmoothingOn();

    // smoothing seems to cause lines on sphere for some reason
    //if (!OptsRenderNodes)   // smoothing is too slow if rendering the nodes
    //   render_win->PolygonSmoothingOn();   // seem to need this else rendering artifact on sphere
    
//    const bool VTK_STEREO = true;
//
//    if (VTK_STEREO)
//    {
//        render_win->StereoCapableWindowOn();
//#if VTK_MAJOR_VERSION > 4
//        render_win->SetStereoTypeToAnaglyph();
//#else
//        render_win->SetStereoTypeToRedBlue();
//#endif
//    }

    render_win->AddRenderer(renderer);
    
    if(POST_VTK_VIEW || VIEW_VTK)	 // only create interactive window if we're actually using it
    {
        iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(render_win);
        iren->Initialize();   
    } 
    else
    {
        if(VTK_HIGHQUAL)	 // increase quality for non-interactive
        {
            render_win->SetAAFrames(7);
        }

    }

    // scale linewidths by aa factor unless rendering to screen
    if (VIEW_VTK)
        linewidths = (double)VTK_AA_FACTOR * VIS_LINETHICKNESS;
    else
        linewidths = 1.0;


    // write the colourmap file
    // this is kind of slow, so only write if doesn't already exist
   
    sprintf(VTK_colmap_filename, "%stempVTKcolmap.png", TEMPDIR);

    ifstream isthere(VTK_colmap_filename);
    
    if (isthere)
    {
         isthere.close();
    }
    else
    {   // colourmap file doesn't exist, so create
        
        char tmpfilename[1024];

        sprintf(tmpfilename,         "%stempVTKcolmap.bmp", TEMPDIR);
        ofstream outbmpfile;
        // write out file with colourmap

        Dbl2d imageR, imageG, imageB;

        const int width  = VTK_WIDTH * VTK_AA_FACTOR;
        const int height = VTK_HEIGHT * VTK_AA_FACTOR;

        imageR.resize(width);
	    imageG.resize(width);
	    imageB.resize(width);

	    for (int x = 0; x<width; x++)
	    {
		    imageR[x].resize(height);
		    imageG[x].resize(height);
		    imageB[x].resize(height);

            fill(imageR[x].begin(),imageR[x].end(),0.0);
            fill(imageG[x].begin(),imageG[x].end(),0.0);
            fill(imageB[x].begin(),imageB[x].end(),0.0);
	    }

        ptheactin->p_nuc->segs.write_colourmap_bitmap(imageR, imageG, imageB); //, VTK_AA_FACTOR); 
    
	    outbmpfile.open(tmpfilename, ios::out | ios::binary | ios::trunc);

        ptheactin->writebitmapheader(outbmpfile, width, height);
	    ptheactin->writebitmapfile(outbmpfile, imageR, imageG, imageB);

        outbmpfile.close();

        // change the black background to transparent, and convert to png:

        char command1[1024];
        sprintf(command1, 
            "%s -transparent black %s %s",
            IMAGEMAGICKCONVERT, tmpfilename, VTK_colmap_filename);
        system(command1);

        // remove the temp bitmap
        sprintf(command1, "rm %s", tmpfilename);
        system(command1);
    }

    texturereadOK = false;

#if VTK_MAJOR_VERSION > 4	
    
    if (OptsUseNucTextureMap)
    {
        ifstream texturefile1(NUCTEXTUREJPEG1, ios::in);
        ifstream texturefile2(NUCTEXTUREJPEG2, ios::in);   
        
        if (texturefile1)
        {
            texturefile1.close();

	        // read texture image
	        tx_reader = vtkJPEGReader::New();
	        cout << "Reading texture file '"<< NUCTEXTUREJPEG1 << "'...";
            cout.flush();
            tx_reader->SetFileName(NUCTEXTUREJPEG1);
            cout << " done." << endl;
            cout.flush();
            texturereadOK = true; // todo: how to test if vtk read the jpeg properly?
        }
        else if (texturefile2)
        {
            texturefile2.close();

            // read texture image
	        tx_reader = vtkJPEGReader::New();
	        cout << "Reading texture file '"<< NUCTEXTUREJPEG2 << "'...";
            cout.flush();
            tx_reader->SetFileName(NUCTEXTUREJPEG2);
            cout << " done." << endl;
            cout.flush();
            texturereadOK = true;
        }
        else
        {
            texturereadOK = false;
            cerr << "Can't read texture file '"<< NUCTEXTUREJPEG1 << "' or '" << NUCTEXTUREJPEG2 << "'" << endl;
        }
    }


#endif
 
}

CometVtkVis::~CometVtkVis()
{
    // finished
    // REVISIT: deleting renderer kills the pipeline?
    // iren->Delete();
    //render_win->Delete();
    
    render_win->RemoveRenderer(renderer);

    render_win->Delete();

    renderer->Delete();


    if(POST_VTK_VIEW || OptVIEW_VTK)	 // increase quality for non-interactive
    {
        iren->Delete();
    }

#if VTK_MAJOR_VERSION > 4
    if (texturereadOK)
        tx_reader->Delete();
#endif
    // remove the scalebar png file

    //char command1[1024];
    //sprintf(command1, "rm %s", VTK_colmap_filename);
    //system(command1);

}

void CometVtkVis::RestartRenderWindow()
{
    render_win->Delete();

    if(POST_VTK_VIEW || OptVIEW_VTK)	 // increase quality for non-interactive
    {
        iren->Delete();

        iren = vtkRenderWindowInteractor::New();
        iren->SetRenderWindow(render_win);
        iren->Initialize();
    }

}

void CometVtkVis::buildVTK(const int &framenumber, vect & cameraposition, vect & cameratarget)
{
    
    setProjection(cameraposition,cameratarget);  // also sets the rotation of the objects, so comes before actors:
    
    // add objects to renderer
                 
    if(OptsRenderNodes)
	    addNodes();
  
    if(OptsRenderLinks)        
	    addLinks();

    if(OptsRenderNucleator)
    {   
        addNucleator(false); // render a solid nucleator
        
        if (VTK_NUC_WIREFRAME)
	        addNucleator(true);  // if wireframe, add another nucleator, in wireframe and slightly bigger
    }

    if(OptsVolumeRenderNodes || OptsIsoRenderNodes)
    {
	    //addGuassianNodeVolume(true);

        vtkStructuredPoints *spoints = vtkStructuredPoints::New();
        createGaussianVoxelRepresentation(spoints);

        if(OptsVolumeRenderNodes)
	        addStructuredPointVolumeRender(spoints);        
                            
        if(OptsIsoRenderNodes)
        {
	        addStructuredPointIsoRender(spoints, 250.0, Colour( 0.3, 0.6, 0.3), 1.0 );
            //addStructuredPointIsoRender(spoints, 150.0, Colour( 0.6, 0.3, 0.2), 0.5 );
            //addStructuredPointIsoRender(spoints, 100.0, Colour( 0.2, 0.3, 0.6), 0.5 );
        }
        
        spoints->Delete();

    }

    if(OptsRenderAxes)
	    addTracks(framenumber);
        //addAxes();
       
    // if(OptsRenderText)        
    // addVoxelBound(renderer);

    addLight();
    
    //renderer->GetActiveCamera()->ViewingRaysModified();

	
	//cout << "Rendering..." << endl;
    //cout.flush();

   
    
    // -- rendering
    if(POST_VTK_VIEW) 
    {
        cout << "Starting interactive mode" << endl;
        
        iren->Start();
    } 
    else 
    {
        if (OptVIEW_VTK)
        {
            // render to screen
            iren->Render();
        }

        // render to file
        saveImage(framenumber);
        //saveVRML(framenumber); // save vrml last (have to change lighting before vrml save for some reason)

    }


    
    // clear all actors from the renderer

#if VTK_MAJOR_VERSION < 5
    renderer->RemoveAllProps();
#else
    renderer->RemoveAllViewProps();
#endif
        
    // this reuse of renderer seems to cause segfaults in vtk.
    // CHECK: order matters in deletion of these objects
    // otherwise too many X Servers error.

    //render_win->RemoveRenderer(renderer);
    
    //renderer->Delete();
}

// -- Helper functions
void CometVtkVis::addAxes()
{
    // see finance.cxx in vtk examples
    vtkAxes *axes = vtkAxes::New();
    VTK_FLOAT_PRECISION origin[3] = {0.0, 0.0, 0.0};
    axes->SetOrigin( origin );
    axes->SetScaleFactor(1.3*voxelscalefactor * RADIUS);
    
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

void CometVtkVis::saveImage(const int &framenumber)
{
    char tmpfilename[1024],filename[1024];
    sprintf(tmpfilename , "%s%s_%05i.bmp", TEMPDIR, file_prefix.c_str(), framenumber );
    sprintf(filename    , "%s%s_%05i.%s", VTKDIR, file_prefix.c_str(), framenumber, BMP_OUTPUT_FILETYPE.c_str());

    //cout << "Saving " << filename << endl;
  

    //vtkWindowToImageFilter *rwin_to_image = vtkWindowToImageFilter::New();
    //rwin_to_image->SetInput(render_win);

    vtkRenderLargeImage *renderLarge = vtkRenderLargeImage::New();
    renderLarge->SetInput(renderer);
    renderLarge->SetMagnification(VTK_AA_FACTOR);

    // vtkBMPWriter seems a bit faster than vtkPNGWriter, and we're compressing with IM to png anyway later
    // so use BMP here
    vtkBMPWriter *imagewriter = vtkBMPWriter::New();
    
    imagewriter->SetInput( renderLarge->GetOutput() ); //rwin_to_image->GetOutput() );
    imagewriter->SetFileName( tmpfilename );
    imagewriter->Write();    
    
    imagewriter->Delete();

    //rwin_to_image->Delete();
    renderLarge->Delete();
    
    char textstring[1024];
    
    int text_height_pixels = int ( 25 * (double) TEXT_POINTSIZE / 20.0);

    if (NO_IMAGE_TEXT)
    {
        *textstring=0;
    }
    else
    {
        sprintf(textstring, "-font helvetica -fill white -pointsize %i -draw \"text +%i+%i 'Frame % 6i' text +%i+%i 'Time % 6i'\"",
            TEXT_POINTSIZE, 5 , text_height_pixels,
            framenumber,
            5 , text_height_pixels * 2 ,
            int(framenumber * InterRecordIterations * DELTA_T) );
    }

    if (NO_COLBAR)
    {
        if (VTK_AA_FACTOR!=1)
        {
            char command1[1024];
            sprintf(command1, 
                "(%s -quality %i -resize %f%% %s %s %s ; rm %s ) &",
                IMAGEMAGICKCONVERT, BMP_COMPRESSION, 100/(double)VTK_AA_FACTOR, textstring,  
                  tmpfilename, filename, tmpfilename);
            //cout << command1 << endl;
            system(command1);
        }
        else
        {
            char command1[1024];
            
            sprintf(command1, 
                "(%s -quality %i %s %s %s ; rm %s ) &",
                IMAGEMAGICKCONVERT, BMP_COMPRESSION, textstring,
                tmpfilename, filename, tmpfilename);

            //cout << command1 << endl;
            system(command1);
        }
    }
    else
    {
        if (VTK_AA_FACTOR!=1)
        {
            char command1[1024];
            sprintf(command1, 
                "(%s -compose Dst_Over -composite -quality %i -resize %f%% %s %s %s %s ; rm %s ) &",
                IMAGEMAGICKCONVERT, BMP_COMPRESSION, 100/(double)VTK_AA_FACTOR, textstring,  
                 VTK_colmap_filename, tmpfilename, filename, tmpfilename);
            //cout << command1 << endl;
            system(command1);
        }
        else
        {
            char command1[1024];
            sprintf(command1, 
                "(%s -compose Dst_Over -composite -quality %i %s %s %s %s ; rm %s ) &",
                IMAGEMAGICKCONVERT, BMP_COMPRESSION, textstring,
                VTK_colmap_filename, tmpfilename, filename, tmpfilename);
            //cout << command1 << endl;
            system(command1);
        }
    }

}

void CometVtkVis::saveVRML(const int &framenumber)
{   // for some reason this produces *huge* files (>16MB) when nucleator texture is turned on
    // it also doesn't seem to work at all...(lighting problem?)

    char vrmltmpfilename[1024],vrmlfilename[1024];
    sprintf(vrmltmpfilename , "%s%s_%05i.wrl", TEMPDIR, file_prefix.c_str(), framenumber );
    sprintf(vrmlfilename , "%s%s_%05i.wrz", VRMLDIR, file_prefix.c_str(), framenumber );


    // for some reason directional lights don't seem to work for vrml
    // so remove them and switch on ambient lighting:

    renderer->GetLights()->RemoveAllItems();

    renderer->SetAmbient(0.8,0.8,0.8);


    cout << "Saving " << vrmlfilename << endl;

	vtkVRMLExporter *VRMLExporter = vtkVRMLExporter::New();

	VRMLExporter->SetInput( render_win );
    VRMLExporter->SetFileName( vrmltmpfilename );
    VRMLExporter->Write();    

    VRMLExporter->Delete();
  
    char command1[1024];
    sprintf(command1 , "gzip -f -c %s > %s ", vrmltmpfilename, vrmlfilename );
    system(command1);
}



void CometVtkVis::addNucleator(bool wireframe)
{

    // creat assembly and add nucleator
    vtkAssembly *nucleator_actor=vtkAssembly::New();

    if (NUCSHAPE == nucleator::capsule) 
        add_capsule_to_assembly(nucleator_actor, wireframe);   
    else 
	    add_sphere_to_assembly(nucleator_actor, wireframe);  // sphere or ellipsoid


    // set nucleator rotation and translation matrix

    vtkMatrix4x4 *transformationmatrix = vtkMatrix4x4::New();
    set_transform_matrix(transformationmatrix, ptheactin->nuc_to_world_rot);
    nucleator_actor->SetUserMatrix(transformationmatrix);

    
    // make wireframe slightly bigger, so lines appear on outside
    if (wireframe)
        nucleator_actor->SetScale(1.015);


    // add nucleator
    renderer->AddActor(nucleator_actor);

    nucleator_actor->Delete();
    transformationmatrix->Delete();
}

void CometVtkVis::add_capsule_to_assembly(vtkAssembly* nucleator_actor, bool wireframe)
{
    // -- Sources
    // -- endcap source
    // create top sphere geometry
    vtkSphereSource *endcap = vtkSphereSource::New();
    endcap->SetRadius(radius_pixels); // can't we get radius through nucleator:: ??
    endcap->SetThetaResolution(32);
    endcap->SetPhiResolution(16);    // half the phi resolution because it's a hemisphere
    endcap->SetStartPhi(0);
    endcap->SetEndPhi(90);

    // endcap mapper
    vtkPolyDataMapper *endcap_mapper = vtkPolyDataMapper::New();
    SetFocalDepthPlanes(endcap_mapper);
    endcap_mapper->SetInput( endcap->GetOutput() );
    endcap->Delete();
    
    // - body source
    vtkCylinderSource *body = vtkCylinderSource::New();
    body->SetRadius(radius_pixels); // why not through nucleator:: ??
    body->SetHeight(2*capsule_length_pixels);
    body->SetResolution(32);
    body->CappingOff();
   

    // body mapper
    vtkPolyDataMapper *body_mapper = vtkPolyDataMapper::New();
    SetFocalDepthPlanes(body_mapper);
    body_mapper->SetInput( body->GetOutput() );
    body->Delete();
  
    // -- Actors
    
    // - upper endcap
    // actor coordinates geometry, properties, transformation
    vtkActor *endcap1_actor = vtkActor::New();
    endcap1_actor->SetMapper(endcap_mapper);
 
    endcap1_actor->SetOrientation( 0, 0, 0);

    endcap1_actor->SetPosition( 0,0,  CAPSULE_HALF_LINEAR*voxelscalefactor);
    endcap1_actor->GetProperty()->SetColor(0.7, 0.7, 0.7);	
    endcap1_actor->GetProperty()->SetOpacity(nuc_opacity);
    
    // - lower endcap
    vtkActor *endcap2_actor = vtkActor::New();
    endcap2_actor->SetMapper(endcap_mapper);

    endcap2_actor->SetOrientation( 180, 0, 0);

    endcap2_actor->SetPosition( 0, 0, -CAPSULE_HALF_LINEAR*voxelscalefactor);
    endcap2_actor->GetProperty()->SetColor(0.7, 0.7, 0.7);	
    endcap2_actor->GetProperty()->SetOpacity(nuc_opacity);
    
    endcap_mapper->Delete();
    
    // - body
    vtkActor *body_actor = vtkActor::New();
    body_actor->SetMapper(body_mapper);
    
    body_actor->SetOrientation(  90, 0, 0 );
    body_actor->SetPosition( 0,0,0 );
    body_actor->GetProperty()->SetColor(0.7, 0.7, 0.7);	
    body_actor->GetProperty()->SetOpacity(nuc_opacity);
    body_mapper->Delete();


    if (wireframe)
    {
        endcap1_actor->GetProperty()->SetRepresentationToWireframe();
        endcap2_actor->GetProperty()->SetRepresentationToWireframe();
        body_actor->GetProperty()->SetRepresentationToWireframe();

        endcap1_actor->GetProperty()->SetColor(0, 0, 0);
        endcap2_actor->GetProperty()->SetColor(0, 0, 0);
        body_actor->GetProperty()->SetColor(0, 0, 0);

        endcap1_actor->GetProperty()->SetLineWidth(linewidths / 2);
        endcap2_actor->GetProperty()->SetLineWidth(linewidths / 2);
        body_actor->GetProperty()->SetLineWidth(linewidths / 2);
    }
    
    //assemble nucleator from parts

    nucleator_actor->AddPart(endcap1_actor);
    nucleator_actor->AddPart(endcap2_actor);
    nucleator_actor->AddPart(body_actor);
    endcap1_actor->Delete();
    endcap2_actor->Delete();
    body_actor->Delete();
}

void CometVtkVis::add_sphere_to_assembly(vtkAssembly* nucleator_actor, bool wireframe)
{

    // create sphere geometry
    vtkSphereSource *sphere = vtkSphereSource::New();
    
    //cout << "  voxel_scale: " << voxel_scale << endl;
    //cout << "  nucleator radius: " << radius_pixels << endl;

    sphere->SetRadius(radius_pixels);
    sphere->SetThetaResolution(32);
    sphere->SetPhiResolution(32);
    //sphere->LatLongTessellationOn();

    vtkActor *nuc_actor = vtkActor::New();
    
    nuc_actor->GetProperty()->SetOpacity(nuc_opacity);

    nuc_actor->GetProperty()->SetDiffuse(0.25);
    nuc_actor->GetProperty()->SetAmbient(0.3);
    nuc_actor->GetProperty()->SetSpecular(1.0);
    nuc_actor->GetProperty()->SetSpecularPower(5.0);
    
    vtkTransform *ellipseTransform = vtkTransform::New();
    vtkTransformPolyDataFilter *ellipseFilter = vtkTransformPolyDataFilter::New();
    ellipseFilter->SetInput( sphere->GetOutput() );

    if (NUCSHAPE == nucleator::ellipsoid)
        ellipseTransform->Scale(1,1,ELLIPSOID_STRETCHFACTOR);     

    ellipseFilter->SetTransform( ellipseTransform );

     // mapper
    vtkPolyDataMapper *mapper = vtkPolyDataMapper::New(); 

    SetFocalDepthPlanes(mapper);


    if((OptsUseNucTextureMap) && (VTK_MAJOR_VERSION > 4) && texturereadOK) 
    {
      // texture commands not in VTK 4 or below, so skip:
 #if VTK_MAJOR_VERSION > 4	
        
        // add texture map to the nucleator,
	    // see vtk example: 'GenerateTextureCoords.tcl'

	    // create texturemap to sphere
	    vtkTextureMapToSphere *tx_mapper = vtkTextureMapToSphere::New();
	    tx_mapper->PreventSeamOff();
        tx_mapper->SetInput( ellipseFilter->GetOutput() );
	    vtkTransformTextureCoords *tx_xfm =  vtkTransformTextureCoords::New();
	    tx_xfm->SetInput( tx_mapper->GetOutput() );
    	
        // add the texture to the nucleator
	    vtkTexture *tx = vtkTexture::New();

        tx->SetInputConnection( tx_reader->GetOutputPort() );
        nuc_actor->SetTexture( tx );

	    // set the mapper
	    mapper->SetInputConnection( tx_xfm->GetOutputPort() );
 	
	    // clear up

	    tx_mapper->Delete();
	    tx_xfm->Delete();
	    tx->Delete();

 #endif


    } else 
    {
	    mapper->SetInput(ellipseFilter->GetOutput());
	    nuc_actor->GetProperty()->SetColor(0.7, 0.7, 0.7); // sphere color 
    }

    if (wireframe)
    {
      nuc_actor->GetProperty()->SetRepresentationToWireframe();
      nuc_actor->GetProperty()->SetColor(0, 0, 0);
    }

    nuc_actor->SetMapper(mapper);

    sphere->Delete();
    mapper->Delete();              

    nucleator_actor->AddPart(nuc_actor);
    nuc_actor->Delete();
}

void CometVtkVis::set_transform_matrix(vtkMatrix4x4 * vtkmat, const rotationmatrix & rotmat) const
{   // convert 3x3 rotationmatrix into 4x4 vtk transform matrix

    // set rotation

    vtkmat->SetElement(0,0,rotmat.xx);
    vtkmat->SetElement(1,0,rotmat.xy);
    vtkmat->SetElement(2,0,rotmat.xz);
    vtkmat->SetElement(3,0,0);

    vtkmat->SetElement(0,1,rotmat.yx);
    vtkmat->SetElement(1,1,rotmat.yy);
    vtkmat->SetElement(2,1,rotmat.yz);
    vtkmat->SetElement(3,1,0);
    
    vtkmat->SetElement(0,2,rotmat.zx);
    vtkmat->SetElement(1,2,rotmat.zy);
    vtkmat->SetElement(2,2,rotmat.zz);
    vtkmat->SetElement(3,2,0);

    // default zero translation
    //vtkmat->SetElement(0,3,0);
    //vtkmat->SetElement(1,3,0);
    //vtkmat->SetElement(2,3,0);
    //vtkmat->SetElement(3,3,1);

    // set translation

    vtkmat->SetElement(0,3,ptheactin->p_nuc->position.x * voxelscalefactor);
    vtkmat->SetElement(1,3,ptheactin->p_nuc->position.y * voxelscalefactor);
    vtkmat->SetElement(2,3,ptheactin->p_nuc->position.z * voxelscalefactor);
    vtkmat->SetElement(3,3,1);

}


void CometVtkVis::fillVoxelSetFromActinNodes(vector< vector< vector<double > > >  &vx)
					    // double vd, double *min)
{
    // zero the voxel set
    for(int i=0; i<ni; ++i) 
    {
	    for(int j=0; j<nj; ++j) 
        {
            fill(vx[i][j].begin(),vx[i][j].end(),0);
	        //for(int k=0; k<nk; ++k) 
            //{		            
		    //    vx[i][j][k] = 0;
	        //}
	    }
    }
   
  
    // move to static
    // Gaussian splat
    double gaussmax = (double) GAUSSFWHM * 3/2.0;  
    // full extent of gaussian radius -  fwhm is 2/3 this
  
    const int splat_sz = (int)(voxelscalefactor * gaussmax) ; // as in actin 
    cout << "  gaussian splat extent: "<< splat_sz << endl;
  
    vector< vector< vector<double > > >  splat;
    splat.resize(2*splat_sz+1);
    for(int i=0; i<2*splat_sz+1 ;i++) 
    {
	    splat[i].resize(2*splat_sz+1);
	    for(int j=0; j<2*splat_sz+1; j++) 
        {
	        splat[i][j].resize(2*splat_sz+1); // extent x2
	    }
    }
  

    // create splat to reuse
    for(int i=-splat_sz; i<=splat_sz; i++) 
    {
	    for(int j=-splat_sz; j<=splat_sz; j++) 
        {
	        for(int k=-splat_sz; k<=splat_sz; k++) 
            {
	
		    // cut the corners
		    if(i*i + j*j + k*k > splat_sz*splat_sz )
		        continue;
    	
		    splat[i+splat_sz][j+splat_sz][k+splat_sz]
		        = exp(-3*((double)(i*i+j*j+k*k)) / 
			      (double)(splat_sz*splat_sz*splat_sz) );
	        }
	    }
    }
  
  
    // loop over the nodes, and splat them into the volume
    vect node_pos;
    int x,y,z;
    int vi,vj,vk;
    double vx_imax = 0;
    double vx_isum = 0;
    int n_isums = 0;
    for(int n=0; n<ptheactin->highestnodecount; n++) 
    {
	    if (!ptheactin->node[n].polymer)
	        continue;
          
	    node_pos = ptheactin->node[n]; // copy x,y,z coords from node[n]
        
	    //ptheactin->nuc_to_world_rot.rotate(node_pos); 
	    //vtk_cam_rot.rotate(node_pos); // bring rip to y-axis
        
	    // node centre in local voxel coords (nuc at centre) 
	    // REVISIT: check why we need to add one here
	    x = int(voxelscalefactor * node_pos.x + ni/2 + 1); 
	    y = int(voxelscalefactor * node_pos.y + nj/2 + 1);
	    z = int(voxelscalefactor * node_pos.z + nk/2 + 1);
        
        
	    if((x<0) || (x>=ni) ||
	       (y<0) || (y>=nj) ||
	       (z<0) || (z>=nk) )  // only plot if point in bounds
	    {
	        //  cout << "point out of bounds " << x << "," << y << endl;
	    } 
        else
        {
	        // if((n%SPECKLE_FACTOR)==0)
          
	        for(int i=-splat_sz; i<=splat_sz; ++i) 
            {
		        for(int j=-splat_sz; j<=splat_sz; ++j) 
                {
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
  
    if(OptsNormaliseFrames)
    {
	    cout << "  normalising Frame, vx_max: " << vx_imax << endl;
	    // scale the data
	    for(int i=0; i<ni; ++i) 
        {
	        for(int j=0; j<nj; ++j) 
            {
		        for(int k=0; k<nk; ++k) 
                {		
		            vx[i][j][k] = vx[i][j][k]/vx_imax;
		        }
	        }
	    }
        
    } 
  
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
    for(int i=0; i< ni; i++) 
    {
	    vx[i].resize(nj);
	    for (int j=0; j<nj; j++) 
        {
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
    
    if(OptsNormaliseFrames) 
    {
	    cout << "  normalised frame" << endl;
    } 
    else 
    {
	    cout << "  using fixed intensity scale:" << vx_intensity_scale << endl;
	    iscale = vx_intensity_scale;
    }
    
    // loop over voxel set converting to structured points suitable for 
    // vtk volume rendering
    int offset = 0;
    for(int i=0; i<ni; ++i) 
    {
        for(int j=0; j<nj; ++j) 
        {
	        for(int k=0; k<nk; ++k) 
            {
		
		    double value = vx[i][j][k]/iscale; 
		    if(value > 0) 
            {
		        if(value > 1) 
                {
			        vxdata->InsertValue(offset, 255 );
		        } 
                else 
                {
			        vxdata->InsertValue(offset, (int)(value*255) );
		        }
		    } 
            else 
            {
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


void CometVtkVis::addflow()
{
    // Accumulate a gaussian splat
  
    // find out extent of the data
    // initialise extents
  
    // create the voxel set as vtk structured point data
    // create a structured grid for the nodes, and sample them
    // simple binning of the node data into a voxel set, with
    // gaussian contributions.
  
    // create the voxel 3d vector
    vector< vector< vector<vect > > >  flow;
    flow.resize(ni);  
    for(int i=0; i< ni; i++) 
    {
	    flow[i].resize(nj);
	    for (int j=0; j<nj; j++) 
        {
	        flow[i][j].resize(nk);
	    
            fill(flow[i][j].begin(),flow[i][j].end(),vect(0,0,0));
        }
    }
  
    // fill voxel set from the actin node data
    cout << "  filling vtk voxel flow data " << endl;
//    const double vd = 1;
//    VTK_FLOAT_PRECISION vx_orig[3] = {-ni/2.0-1, -nj/2.0-1, -nk/2.0-1};
//    fillVoxelSetFromActinNodes(vx); //, vd, vx_orig);
    

    
    // loop over the nodes, and splat them into the volume
    //vect node_pos;
    int x,y,z;
    //int vi,vj,vk;
    double flow_imax = 0;
//    double flow_isum = 0;
//   int n_isums = 0;

    for(vector <nodes>::iterator	i_node  = ptheactin->node.begin(); 
									i_node != ptheactin->node.begin() + ptheactin->highestnodecount;
							      ++i_node)
    {

	    if (!i_node->polymer)
	        continue;
   
        
	    //ptheactin->nuc_to_world_rot.rotate(node_pos); 
	    //vtk_cam_rot.rotate(node_pos); // bring rip to y-axis
        
	    // node centre in local voxel coords (nuc at centre) 
	    // REVISIT: check why we need to add one here
	    x = int(voxelscalefactor * i_node->x + ni/2 + 1); 
	    y = int(voxelscalefactor * i_node->y + nj/2 + 1);
	    z = int(voxelscalefactor * i_node->z + nk/2 + 1);
        
        
	    if((x<0) || (x>=ni) ||
	       (y<0) || (y>=nj) ||
	       (z<0) || (z>=nk) )  // only plot if point in bounds
	    {
	        //  cout << "point out of bounds " << x << "," << y << endl;
	    } 
        else
        {
      
    	    
            // amount of actin		  
            flow[x][y][z] += i_node->delta_sum;
	        
            // keep tabs on the maximum
	        if(flow[x][y][z].length() > flow_imax)
	            flow_imax = flow[x][y][z].length();


	    } // else (add intensity)
        
    } // loop over nodes

    //cout << "  max gausian intensity: "  << flow_imax << endl;
    //cout << "  mean gausian intensity: " << flow_isum/n_isums << endl;
  
    if(OptsNormaliseFrames)
    {
	    cout << "  normalising Frame, flow_imax: " << flow_imax << endl;
	    // scale the data
	    for(int i=0; i<ni; ++i) 
        {
	        for(int j=0; j<nj; ++j) 
            {
		        for(int k=0; k<nk; ++k) 
                {		
		            flow[i][j][k] = flow[i][j][k] * flow_imax;
		        }
	        }
	    }
        
    } 




    // set up the spoints for volume rendering
    int dim[3] = {ni, nj, nk};

 
  
    // create the underlying voxel data (only UnsignedChar or Short)
    // see vtk examples for voxel data

    vtkFloatArray *vectors = vtkFloatArray::New();
    vectors->SetNumberOfComponents(3);
    vectors->SetNumberOfTuples(dim[0]*dim[1]*dim[2]);

    vtkPoints *points = vtkPoints::New();
    points->Allocate(dim[0]*dim[1]*dim[2]);

    //points->SetDimensions(dim);
    //points->SetSpacing(vd, vd, vd);
    //points->SetOrigin( vx_orig );

    //vtkUnsignedCharArray *vxdata = vtkUnsignedCharArray::New();    

    
    VTK_FLOAT_PRECISION vel[3], pos[3];

    double vd = 1.0;

    // loop over voxel set converting to structured points suitable for 
    // vtk volume rendering
    int offset = 0;
    for(int i=0 ; i<ni ; ++i)
    {   
        pos[0] = (i - ni/2 - 1) * vd ; // ((VTK_FLOAT_PRECISION)(i - ni/2 - 1)) / voxelscalefactor ;

	    for(int j=0 ; j<nj ; ++j)
        {
            pos[1] = (j - nj/2 - 1) * vd ; //((VTK_FLOAT_PRECISION)(j - nj/2 - 1)) / voxelscalefactor;

	        for(int k=0 ; k<nk; ++k)
            {
		        pos[2] = (k - nk/2 - 1) * vd ; //((VTK_FLOAT_PRECISION)(k - nk/2 - 1)) / voxelscalefactor;

                vel[0] = flow[i][j][k].x;
                vel[1] = flow[i][j][k].y;
                vel[2] = flow[i][j][k].z;

                points->InsertPoint(offset, pos);
                vectors->InsertTuple(offset, vel);

		        offset++;		
            }
        }
    }



    vtkStructuredGrid *sgrid = vtkStructuredGrid::New();
    sgrid->SetDimensions(dim);

    sgrid->SetPoints(points);

    sgrid->GetPointData()->SetVectors(vectors);
    vectors->Delete();

    vtkHedgeHog *hedgehog = vtkHedgeHog::New();
    hedgehog->SetInput(sgrid);
    hedgehog->SetScaleFactor(10);   



    //vtkRungeKutta4 *integ = vtkRungeKutta4::New();

    //vtkStreamTracer *flowStreamer = vtkStreamTracer::New();

    //flowStreamer->SetInput(sgrid->GetO);
    //flowStreamer->SetSource(hedgehog->GetOutput());
    //flowStreamer->SetStartPosition(RADIUS, RADIUS, RADIUS);
    //flowStreamer->SetMaximumPropagationTime(500);
    //flowStreamer->SetStepLength(0.25);
    //flowStreamer->SetIntegrationStepLength(0.025);
    //flowStreamer->SetIntegrationDirectionToIntegrateBothDirections();
    //flowStreamer->SetIntegrator(integ);

    //vtkTubeFilter *flowtube = vtkTubeFilter::New();

    //flowtube->SetInput(flowStreamer->GetOutput());
    //flowtube->SetRadius(0.01);
    //flowtube->SetNumberOfSides(12);
    //flowtube->SetVaryRadiusToVaryRadiusByVector();

    //vtkPolyDataMapper *flowtubeMapper = vtkPolyDataMapper::New();
    //flowtubeMapper->SetInput(flowtube->GetOutput());

    //vtkActor *flowtubeactor = vtkActor::New();
    //flowtubeactor->SetMapper(flowtubeMapper);

    //renderer->AddActor(flowtubeactor);

    //flowtubeactor->Delete();
    
    vtkPolyDataMapper *sgridMapper = vtkPolyDataMapper::New();
    sgridMapper->SetInput(hedgehog->GetOutput());

    vtkActor *sgridActor = vtkActor::New();
    sgridActor->SetMapper(sgridMapper);
    sgridActor->GetProperty()->SetColor(1,1,1);

    renderer->AddActor(sgridActor);

    sgridActor->Delete();
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

void CometVtkVis::addStructuredPointIsoRender(vtkStructuredPoints *sp, const double threshold, const Colour col, const double opacity)
{
    // render using an iso-surface
    vtkContourFilter *iso_surf = vtkContourFilter::New();
    iso_surf->SetInput( sp );
    iso_surf->SetValue(0, threshold);
    iso_surf->ComputeNormalsOn();

	vtkWindowedSincPolyDataFilter *smooth = vtkWindowedSincPolyDataFilter::New();
    smooth->SetInput(iso_surf->GetOutput());
    smooth->SetNumberOfIterations( 200 );
//	smooth->SetRelaxationFactor( 0.02 );
    smooth->SetFeatureAngle(45);

    vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
    normals->SetInput( smooth->GetOutput());
    normals->SetFeatureAngle(45);
    normals->FlipNormalsOn();

    //smooth->BoundarySmoothingOn();

    vtkPolyDataMapper *iso_mapper = vtkPolyDataMapper::New();
	iso_mapper->SetInput( normals->GetOutput() );
    iso_mapper->ScalarVisibilityOff();
    iso_surf->Delete();
    
    vtkActor *iso_actor = vtkActor::New();
    iso_actor->SetMapper(iso_mapper);
    iso_mapper->Delete();
    
    iso_actor->GetProperty()->SetOpacity(opacity);
    iso_actor->GetProperty()->SetColor(col.r, col.g, col.b);
    
    renderer->AddActor(iso_actor);
    iso_actor->Delete();
    smooth->Delete();
    normals->Delete();
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
    for(int i=0; i<ptheactin->highestnodecount; i++)
    {
	pts->InsertPoint(i,
			 ptheactin->node[i].x - ptheactin->p_nuc->position.x,
			 ptheactin->node[i].y - ptheactin->p_nuc->position.y,
			 ptheactin->node[i].z - ptheactin->p_nuc->position.z);
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
    sphere->SetRadius( 0.01*voxelscalefactor * (RADIUS) ); // scale with the nucleator
    sphere->SetThetaResolution(10); // low res is fine
    sphere->SetPhiResolution(10);
    
    // map
    vtkPolyDataMapper *map = vtkPolyDataMapper::New();

    SetFocalDepthPlanes(map);
    map->SetInput(sphere->GetOutput());
    sphere->Delete();

    Colour col;

    //const double VTK_MAX_NODE_DIST = 10000000; // 0.1;

    // loop over the nodes
    VTK_FLOAT_PRECISION ncx, ncy, ncz;
    for(int i=0; i < ptheactin->highestnodecount; i++)
    {
        if (ptheactin->node[i].dist_from_surface > VTK_MAX_NODE_DIST)
            continue;

	    ncx = ptheactin->node[i].x;
	    ncy = ptheactin->node[i].y;
	    ncz = ptheactin->node[i].z;

        convert_to_vtkcoord(ncx, ncy, ncz);
    	
	    // actor coordinates geometry, properties, transformation
	    vtkActor *node_actor = vtkActor::New();
	    node_actor->SetMapper(map);
	    node_actor->SetPosition(ncx, ncy, ncz);
    	
	    
        if (VTK_RAND_NODE_COL)
        {
            col.setcol( double(ptheactin->node[i].creation_iter_num % 200) / 200);
            
            if(ptheactin->node[i].harbinger)
            {
                col.r = 1.0;
                col.g = 1.0;
                col.b = 1.0;
            }
            
            node_actor->GetProperty()->SetColor(col.r, col.g, col.b); 
        }
        else 
        {              
            if ( i < ptheactin->lowestnodetoupdate )
            {
	            node_actor->GetProperty()->SetColor(1.0,1.0,0.0); // plot harbinger red
	        } 
            else if ( ptheactin->node[i].harbinger )
            {
	            node_actor->GetProperty()->SetColor(1.0,0.0,0.0); // plot harbinger red
	        } 
            else if ( ptheactin->node[i].listoflinks.empty() )
            {
	            node_actor->GetProperty()->SetColor(0.5, 0.5, 0.5); // plot empty grey
	        } 
            else if ( ptheactin->node[i].polymer ) 
            { 
                node_actor->GetProperty()->SetColor(0.2, 1.0, 0.2);     // plot polymerised green
	        } 
            else
                node_actor->GetProperty()->SetColor(0.2, 0.2, 1.0);
        
            //node_actor->GetProperty()->SetColor(1.0, 0.2, 0.2);
        
        }

        

	    // add the actor to the scene	
	    renderer->AddActor(node_actor);
	    node_actor->Delete();
    }
    
    map->Delete();
}

void CometVtkVis::addTracks(const int & framenumber)
{
    


    //double lastx = 0, lasty = 0, lastz = 0;
    VTK_FLOAT_PRECISION x,y,z;

    //const int TRACKFRAMESTEP = 5;

    // not the most efficient algorythm, but good enough...

    // find the track node numbers, make sure more than two of each

    vector <int> temptracknodenumbers, tracknodenumbers;
    temptracknodenumbers.resize(0);
    tracknodenumbers.resize(0); 
            
    //cout << "Tot Nodes to track: " << ptheactin->node_tracks[xaxis].size() << endl;

    for (vector <tracknodeinfo>::iterator i_trackpoint  = ptheactin->node_tracks[xaxis].begin(); 
                                          i_trackpoint != ptheactin->node_tracks[xaxis].end(); 
                                        ++i_trackpoint)
    {
        if ((i_trackpoint->frame % TRACKFRAMESTEP == 0) &&           // only plot points *added* every trackframestep frames
            (framenumber > i_trackpoint->frame))  // and only if we're in a frame in which they exist!  
        {        
            temptracknodenumbers.push_back(i_trackpoint->nodenum);
            //cout << "Adding track point: frame " << i_trackpoint->frame << " itternum " << ptheactin->iteration_num << endl;
        }
    }

    //cout << "Selected Nodes to track: " << temptracknodenumbers.size() << endl;

    for (vector <int>::iterator i_pointnodenum  = temptracknodenumbers.begin(); 
                                i_pointnodenum != temptracknodenumbers.end(); 
                              ++i_pointnodenum)
    {
        if (count(temptracknodenumbers.begin(), temptracknodenumbers.end(), *i_pointnodenum) > 2)
        { // make sure more than 3 points for the line
            if (find(tracknodenumbers.begin(), tracknodenumbers.end(), *i_pointnodenum) == tracknodenumbers.end())
            {   // if node not already there, add it
                tracknodenumbers.push_back(*i_pointnodenum);  
            }
        }
    }

    //cout << "ReSelected Nodes to track: " << tracknodenumbers.size() << endl;

    // go through node by node

    for (vector <int>::iterator i_pointnodenum  = tracknodenumbers.begin(); 
                                i_pointnodenum != tracknodenumbers.end(); 
                              ++i_pointnodenum)
    {

        // we do this for each line...
     
        // add the track points

        vtkCardinalSpline *aSplineX = vtkCardinalSpline::New();
        vtkCardinalSpline *aSplineY = vtkCardinalSpline::New();
        vtkCardinalSpline *aSplineZ = vtkCardinalSpline::New();

        vtkPoints *inputPoints = vtkPoints::New();

        vtkPoints *outputPoints = vtkPoints::New();


        int time = 0;

        for (vector <tracknodeinfo>::iterator i_trackpoint  = ptheactin->node_tracks[xaxis].begin(); 
                                              i_trackpoint != ptheactin->node_tracks[xaxis].end(); 
                                            ++i_trackpoint)
        {   
            if (i_trackpoint->nodenum == *i_pointnodenum)
            {
                x=i_trackpoint->x + ptheactin->p_nuc->position.x;
                y=i_trackpoint->y + ptheactin->p_nuc->position.y;
                z=i_trackpoint->z + ptheactin->p_nuc->position.z;

                convert_to_vtkcoord(x, y, z);

                aSplineX->AddPoint(time, x);
                aSplineY->AddPoint(time, y);
                aSplineZ->AddPoint(time, z);

                inputPoints->InsertPoint(time, x, y, z); 

                ++time;

            }
        }

        int numberOfInputPoints = time;

        // now have the track points for the line, interpolate...

        int numberOfOutputPoints = 5 * numberOfInputPoints; // spline will have 20 times as many points

        vtkPolyData *profileData = vtkPolyData::New();

        vtkCellArray *lines = vtkCellArray::New();
        lines->InsertNextCell(numberOfOutputPoints);
        
        //cout << "Line" << endl;

        for (int i = 0; i < numberOfOutputPoints; ++i)
        {
            double time = (double)(numberOfInputPoints - 1) / (double)(numberOfOutputPoints - 1) * (double)i;
            outputPoints->InsertPoint(i, aSplineX->Evaluate(time), aSplineY->Evaluate(time), aSplineZ->Evaluate(time));
            
            //cout << "Point:" << aSplineX->Evaluate(time) << "," << aSplineY->Evaluate(time) << "," << aSplineZ->Evaluate(time) << endl;

            lines->InsertCellPoint(i);
        }

        profileData->SetPoints(outputPoints);
        profileData->SetLines(lines);
            
        // Add thickness to the resulting line.
        vtkTubeFilter *profileTubes = vtkTubeFilter::New();
        profileTubes->SetNumberOfSides(6);
        profileTubes->SetInput(profileData);
        profileTubes->SetRadius(0.2);
  

        vtkPolyDataMapper *profileMapper = vtkPolyDataMapper::New();
        profileMapper->SetInput(profileTubes->GetOutput());

        vtkActor *profile = vtkActor::New();

        profile->SetMapper(profileMapper);
        profile->GetProperty()->SetDiffuseColor(1.0,0.0,0.0);
        profile->GetProperty()->SetSpecular(.3);
        profile->GetProperty()->SetSpecularPower(30);

        renderer->AddActor(profile);

        profile->Delete();
        profileMapper->Delete();
        profileTubes->Delete();
        profileData->Delete();
        lines->Delete();

        aSplineX->Delete();
        aSplineY->Delete();
        aSplineZ->Delete();

        inputPoints->Delete();
        outputPoints->Delete();
    }




}





bool CometVtkVis::convert_to_vtkcoord(VTK_FLOAT_PRECISION &x, VTK_FLOAT_PRECISION &y, VTK_FLOAT_PRECISION &z)
{

    //ptheactin->nuc_to_world_rot.rotate(x, y, z); 
    //vtk_cam_rot.rotate(x, y, z); // bring rip to y-axis

    bool infocus = !(OptsSkipOutOfFocusPoints && fabs(x) > FOCALDEPTH);

    x *= voxelscalefactor; 
    y *= voxelscalefactor;
    z *= voxelscalefactor;

    return infocus;
}

void CometVtkVis::set_mean_posns()
{
  
  meanx = meany = meanz =0.0;

  if(!VTK_MOVE_WITH_BEAD) {
    meanx = -ptheactin->p_nuc->position.x; 
    meany = -ptheactin->p_nuc->position.y; 
    meanz = -ptheactin->p_nuc->position.z;
  }
  
  ptheactin->nuc_to_world_rot.rotate(meanx, meany, meanz);
  vtk_cam_rot.rotate(meanx, meany, meanz); 
  
  // stops bead 
  double keep_within_border;
  if( NUCSHAPE == nucleator::sphere )
    keep_within_border = 2*RADIUS;
  else
    keep_within_border = CAPSULE_HALF_LINEAR+(2*RADIUS);
  
  int beadminx = int(voxelscalefactor * (-keep_within_border - meanx) + ni/2 + 1);
  int beadminy = int(voxelscalefactor * (-keep_within_border - meany) + nj/2 + 1); 
  int beadminz = int(voxelscalefactor * (-keep_within_border - meanz) + nk/2 + 1); 
  int beadmaxx = int(voxelscalefactor * ( keep_within_border - meanx) + ni/2 + 1); 
  int beadmaxy = int(voxelscalefactor * ( keep_within_border - meany) + nj/2 + 1); 
  int beadmaxz = int(voxelscalefactor * ( keep_within_border - meanz) + nk/2 + 1); 
  
  movex=movey=movez=0;
  
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
}


void CometVtkVis::addLinks()
{	
  // Move to fcn
  // temp: reset center:
  
  set_mean_posns();
  
  const int numcol = 256;

  // only used for strain coloring
  vtkLookupTable *lut = vtkLookupTable::New();
  lut->SetNumberOfTableValues(numcol);
  lut->SetRange(0.7, 1.0);
  lut->Build();
  
  Colour col;
  
  VTK_FLOAT_PRECISION rgba[4];
  
  

  if(OptsShadeLinks)
  {
      for(int i=0; i<numcol; i++) 
      {
        double value = (double)i/(double)(numcol-1);
        col.setcol(value);
        
        rgba[0] = col.r;
        rgba[1] = col.g;
        rgba[2] = col.b;

        if (VTK_SCALECOLOPACITY)
            rgba[3] = value;
        else
            rgba[3] = 1.0;	

        lut->SetTableValue(i, rgba);
      }
  }
  else
  {
    rgba[0] = 0.75;
    rgba[1] = 0.75;
    rgba[2] = 0.75;
    rgba[3] = 1.0;

    for(int i=0; i<numcol; i++) 
    {	
        lut->SetTableValue(i, rgba);
    }
  }
  
  // loop over the nodes
  VTK_FLOAT_PRECISION ncx,ncy,ncz; //n_pt[3];
  VTK_FLOAT_PRECISION lptx, lpty, lptz; //l_pt[3];
  vect testpos;
  vect nodeposvec;
  
  // use a cell array to store line data 
  vtkCellArray *cellarray = vtkCellArray::New();
  // with scalars
  vtkFloatArray *cellscalars = vtkFloatArray::New();  
  // made up of lines
  vtkPolyData *linedata = vtkPolyData::New();
  vtkPoints *linepts = vtkPoints::New();

  double force;

  //bool firstpointinfocus, secondpointinfocus; // used to reject links if both points are out of focal planes
  
  for(int i=0; i != ptheactin->highestnodecount; i++) 
  { 
    if (!ptheactin->node[i].polymer)
          continue;
    
    nodeposvec = ptheactin->node[i];

    ncx = ptheactin->node[i].x;
    ncy = ptheactin->node[i].y;
    ncz = ptheactin->node[i].z;		
    
    //firstpointinfocus = 
    convert_to_vtkcoord(ncx, ncy, ncz);
 
    if(!ptheactin->node[i].listoflinks.empty() && !VTK_NUC_LINKS_ONLY) 
    {
      // nodes thisnode = ptheactin->node[i];
      
      // loop over linked nodes
      // check this out if we use this function                          
      for(vector<links>::iterator link_i = ptheactin->node[i].listoflinks.begin(); 
	                             link_i != ptheactin->node[i].listoflinks.end();
	                           ++link_i) 
      {
        // assume links to nodes 'less' than us have already been added
        if( (*link_i).linkednodeptr->nodenum < i)
            continue;

        lptx = link_i->linkednodeptr->x;
        lpty = link_i->linkednodeptr->y;
        lptz = link_i->linkednodeptr->z;

        //secondpointinfocus = 
        convert_to_vtkcoord(lptx, lpty, lptz);

        //if (!firstpointinfocus && !secondpointinfocus)
        //    continue;  // both points out of focus, so skip link

        vect linkvec = nodeposvec - ptheactin->node[link_i->linkednodenumber];//*(link_i->linkednodeptr) ;

        force = link_i->forcesum / InterRecordIterations;
        //double distance = linkvec.length();     
        //link_i->getlinkforces(distance, force);

        double forcescaled = -force / LINK_BREAKAGE_FORCE;;

        //if (prob_to_bool(0.01)) cout << " " << forcescaled;

        if ( forcescaled * 100.0 < VTK_MIN_PLOT_LINK_FORCE_PCT)  // don't plot links with less than this mag
                continue;


        double magcol;

        if (COL_LINK_BY_DIRN)
        {   
            // dot product, and we don't care about sign, just angle                          
            double dirfactor = fabs ( linkvec.unitvec().dot(ptheactin->node[i].unit_vec_posn) );  
            magcol = dirfactor;

        }
        else
        {
            magcol = forcescaled;

            //double strain = fabs(distance-link_i->orig_dist) / link_i->orig_dist;    
            //double y = strain / LINK_BREAKAGE_STRAIN;
        }

        if ( magcol < 0.0 ) magcol = 0.0;
        if ( magcol > 1.0 ) magcol = 1.0;

        //magcol *= magcol; // square to get energy    
        magcol = magcol * 0.9 + 0.1; // prevent zeros because colorscheme makes them black
        //col.setcol(magcol);

        // create line for the link

        vtkIdType linept_ids[2];	
        linept_ids[0] = linepts->InsertNextPoint( ncx, ncy, ncz );
        linept_ids[1] = linepts->InsertNextPoint( lptx, lpty, lptz );
        vtkIdType cell_id = cellarray->InsertNextCell(2, linept_ids);


        // Set scalar for the line
        if(OptsShadeLinks && !ptheactin->node[i].testnode) // colour testnodes white
        {    
            cellscalars->InsertValue(cell_id, magcol);
        } 
        else 
        {
            cellscalars->InsertValue(cell_id, 0.7);   // also see the colour mapping above
        }

      } // links loop
    
    }

     if (ptheactin->node[i].stucktonucleator)
     {
        vect stuck_pos_world_frame = ptheactin->node[i].nucleator_stuck_position;;
        ptheactin->nuc_to_world_frame(stuck_pos_world_frame);

        vect displacement = ptheactin->node[i] - stuck_pos_world_frame;
        double distance = displacement.length();

        force = NUC_LINK_FORCE * distance;

        // note this is scaled by LINK_BREAKAGE_FORCE not NUC_LINK_BREAKAGE_FORCE, to keep colours consistant

        double magcol = force / LINK_BREAKAGE_FORCE;    // note: colour relative to normal link scale
        if ( (force / NUC_LINK_BREAKAGE_FORCE )* 100.0 < VTK_MIN_PLOT_LINK_FORCE_PCT) // only plot ones with this force or greater
            continue;

        //magcol = pow( magcol , 1 / COLOUR_GAMMA);	    
        magcol = magcol * 0.9 + 0.1; // prevent zeros because colorscheme makes them black
        //col.setcol(y);

        // create line for the link
        lptx = stuck_pos_world_frame.x;
        lpty = stuck_pos_world_frame.y;
        lptz = stuck_pos_world_frame.z;

        convert_to_vtkcoord(lptx, lpty, lptz);

        if (VTK_NUC_LINKS_ONLY)
        {   // if monitoring attachment only, amplify the length
            vect linkendpos = stuck_pos_world_frame + displacement * VTK_NUC_LINKS_ONLY_AMPLIFY;

            // and change colour to relative to nuc link breakage force

            magcol = force / NUC_LINK_BREAKAGE_FORCE;
            magcol = magcol * 0.9 + 0.1;

            ncx = linkendpos.x;
            ncy = linkendpos.y;
            ncz = linkendpos.z;

            convert_to_vtkcoord(ncx, ncy, ncz);

        }

        vtkIdType linept_ids[2];		
        linept_ids[0] = linepts->InsertNextPoint( ncx, ncy, ncz );
        linept_ids[1] = linepts->InsertNextPoint( lptx, lpty, lptz );
        vtkIdType cell_id = cellarray->InsertNextCell(2, linept_ids);

        if(OptsShadeLinks && !ptheactin->node[i].testnode) // colour testnodes white
        {
            cellscalars->InsertValue(cell_id, magcol);
        } 
        else 
        {
            cellscalars->InsertValue(cell_id, 0.7);   // also see the colour mapping above
        }	
            
     }

     //const bool VTK_PLOT_NUC_IMPACTS = false; // true;

     if (VTK_PLOT_NUC_IMPACTS)
     {


        vect displacement = ptheactin->node[i].nucleator_impacts * NODE_DIST_TO_FORCE / (double) InterRecordIterations ;
        double distance = displacement.length();

        force = distance;

        // note this is scaled by LINK_BREAKAGE_FORCE not NUC_LINK_BREAKAGE_FORCE, to keep colours consistant

        double magcol = force / LINK_BREAKAGE_FORCE;    // note: colour relative to normal link scale
        if ( (force / NUC_LINK_BREAKAGE_FORCE )* 100.0 < VTK_MIN_PLOT_LINK_FORCE_PCT) // only plot ones with this force or greater
            continue;

        //magcol = pow( magcol , 1 / COLOUR_GAMMA);	    
        magcol = magcol * 0.9 + 0.1; // prevent zeros because colorscheme makes them black
        //col.setcol(y);

        // create line for the link
        lptx = ptheactin->node[i].x;
        lpty = ptheactin->node[i].y;
        lptz = ptheactin->node[i].z;

        convert_to_vtkcoord(lptx, lpty, lptz);

        if (VTK_NUC_LINKS_ONLY)
        {   // if monitoring attachment only, amplify the length
            vect linkendpos = ptheactin->node[i] + displacement * VTK_NUC_LINKS_ONLY_AMPLIFY;

            // and change colour to relative to nuc link breakage force

            magcol = force / NUC_LINK_BREAKAGE_FORCE;
            magcol = magcol * 0.9 + 0.1;

            ncx = linkendpos.x;
            ncy = linkendpos.y;
            ncz = linkendpos.z;

            convert_to_vtkcoord(ncx, ncy, ncz);

        }

        vtkIdType linept_ids[2];		
        linept_ids[0] = linepts->InsertNextPoint( ncx, ncy, ncz );
        linept_ids[1] = linepts->InsertNextPoint( lptx, lpty, lptz );
        vtkIdType cell_id = cellarray->InsertNextCell(2, linept_ids);

        if(OptsShadeLinks && !ptheactin->node[i].testnode) // colour testnodes white
        {
            cellscalars->InsertValue(cell_id, magcol);
        } 
        else 
        {
            cellscalars->InsertValue(cell_id, 0.7);   // also see the colour mapping above
        }	
            
     }

    } // node loop

  linedata->SetPoints(linepts);
  linepts->Delete();
  linedata->SetLines(cellarray);
  
  cellarray->Delete();
  linedata->GetCellData()->SetScalars(cellscalars);
  // mapper
  vtkPolyDataMapper *map = vtkPolyDataMapper::New();
  
  SetFocalDepthPlanes(map);
 
  map->SetLookupTable( lut );
  map->SetInput( linedata );

  linedata->Delete();
  // actor
  vtkActor *lines_actor = vtkActor::New();
  //lines_actor->GetProperty()->SetAmbient(2.0);
  //lines_actor->GetProperty()->SetDiffuse(2.0);
  //lines_actor->GetProperty()->SetSpecular(2.0);
  lines_actor->GetProperty()->SetLineWidth(linewidths);
  lines_actor->SetMapper(map);
  map->Delete();
  //render
  renderer->AddActor(lines_actor);
  
  lut->Delete();
  lines_actor->Delete();

}

void CometVtkVis::SetFocalDepthPlanes(vtkPolyDataMapper *map)
{    // sets the mapper clipping planes to reject data outside of distance of FOCALDEPTH
     // above and below the nucleator (in x direction)
  

  vect focdepposvec = vect(-FOCALDEPTH*voxelscalefactor,0,0);

  ptheactin->sym_break_rotation_to_xy_plane.inverse().rotate(focdepposvec);

  //// the planes are focdepposx above and below the nucleator
  vect nucpos=ptheactin->p_nuc->position * voxelscalefactor; 

  vect planeorigin1 = nucpos + focdepposvec;
  vect planeorigin2 = nucpos - focdepposvec;


  map->RemoveAllClippingPlanes();

  vtkPlane* newplane = vtkPlane::New();
  newplane->SetOrigin(  planeorigin1.x, planeorigin1.y, planeorigin1.z );
  newplane->SetNormal( -focdepposvec.x,-focdepposvec.y,-focdepposvec.z );
  map->AddClippingPlane( newplane );
  newplane->Delete();

  vtkPlane* newplane2 = vtkPlane::New();
  newplane2->SetOrigin(  planeorigin2.x, planeorigin2.y, planeorigin2.z );
  newplane2->SetNormal(  focdepposvec.x, focdepposvec.y, focdepposvec.z );
  map->AddClippingPlane( newplane2 );
  newplane2->Delete();
}

// ML FIXME:  This is wrong, sort out a better method for mean strain
double CometVtkVis::getMeanNodeLinkForce(const int id) 
{
    // loop over the links  
    double strain = 0;
    int n_links = 0;
    double distance;
    vect displacement;
    vect nodeposvec = ptheactin->node[id];
    for(vector<links>::iterator i=ptheactin->node[id].listoflinks.begin(); 
	i!=ptheactin->node[id].listoflinks.end(); i++ ) {	 
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

    renderer->SetAmbient(0,0,0);

	renderer->GetLights()->RemoveAllItems();

    vtkLightCollection *lights = renderer->GetLights();
    int numlights = lights->GetNumberOfItems();
    vtkLight *light_to_remove;
    lights->InitTraversal();
    for(int i=0; i<numlights; ++i) 
    {
        light_to_remove = lights->GetNextItem();
 	    renderer->RemoveLight(light_to_remove);
    }
    // and still there's a ******* headlight somewhere that I can't delete!

    // N.B. also there is an ambient setting for the object---set in the add nucleator function

	// add our own lights

    vect campos;

    renderer->GetActiveCamera()->GetPosition(campos.x,campos.y,campos.z);

    campos-= ptheactin->p_nuc->position;

    vect l1pos=campos * 1.3;
    vect l2pos=campos * 1.3;

    rotationmatrix rot1, rot2;
    
    rot1.rotatematrix(-40 * PI / 180 , zaxis);
    rot1.rotatematrix(50 * PI / 180 , yaxis);
    rot1.rotate(l1pos);
    
    rot2.rotatematrix(40 * PI / 180 , zaxis);
    rot2.rotatematrix(-10 * PI / 180 , yaxis);
    rot2.rotate(l2pos);

    l1pos += ptheactin->p_nuc->position;
    l2pos += ptheactin->p_nuc->position;

	vtkLight *light1 = vtkLight::New();
	light1->SetIntensity(0.8);
	light1->SetColor(1,0.9,0.8);
    light1->SetPosition(l1pos.x, l1pos.y, l1pos.z);
    renderer->AddLight(light1);
    light1->Delete();

    vtkLight *light2 = vtkLight::New();                                                       
	light2->SetIntensity(0.8);
	light2->SetColor(0.8,0.8,1.0);
    light2->SetPosition(l2pos.x, l2pos.y, l2pos.z);
    renderer->AddLight(light2);
    light2->Delete();

}

void CometVtkVis::setProjection()
{
    vect origin(0,0,0);
    setProjection(origin,origin);
}

void CometVtkVis::setProjection(vect & cameraposition,vect & cameratarget)
{
    
    vect camerapos;
    vect focalpoint; 

    if (VTK_MOVE_WITH_BEAD)
    {
		camerapos  = ptheactin->p_nuc->position * voxelscalefactor;
        focalpoint = ptheactin->p_nuc->position * voxelscalefactor; 
	} else
    {
		camerapos  = cameraposition * voxelscalefactor;
        focalpoint = cameratarget   * voxelscalefactor;
	}

    

    // renderer->ResetCamera();

    renderer->GetActiveCamera()->SetFocalPoint(focalpoint.x, focalpoint.y, focalpoint.z);
  
    if (VIS_PARALLELPROJECTION)
        renderer->GetActiveCamera()->ParallelProjectionOn(); // ParallelProjectionOn();
    else
        renderer->GetActiveCamera()->ParallelProjectionOff();

    renderer->GetActiveCamera()->SetViewAngle(VTK_VIEWANGLE);
    // FIXME: ML
    // should scale properly here to a value linked to the render setup
    renderer->GetActiveCamera()->SetParallelScale(p_scale);
    //renderer->ResetCamera();
  
    
    // this sets where camera is relative to the bead

    viewupvect = vect(0,0,-1);

    if (OptsRenderProjection == X )
    {
        vtk_cam_rot = ptheactin->sym_break_rotation_to_xy_plane.inverse();
    }   
    else if (OptsRenderProjection == Y )
    {
        vtk_cam_rot = ptheactin->sym_break_rotation_to_xy_plane.inverse();
        vtk_cam_rot.rotatematrix(0,PI/2,0);
    }
    else if ( OptsRenderProjection == Z )
    {
        vtk_cam_rot = ptheactin->sym_break_rotation_to_xy_plane.inverse();
        vtk_cam_rot.rotatematrix(0,0,-PI/2);

    }
    else if ( OptsRenderProjection == RIP )
    {
        vtk_cam_rot = ptheactin->sym_break_rotation_to_zaxis.inverse();	  // note using different rotation matrix!

        vtk_cam_rot.rotatematrix(0,VTK_RIP_Z_ANGLE*(PI/180),0); // angle of view

    }
    else 
    {
      cout << "!ERROR: unknown projection:" << OptsRenderProjection << endl;
    }

    // set position of the camera relative to the origin
    vect camera_posn_vect(-radius_pixels*OptsCameraDistMult,0,0);

    vtk_cam_rot.rotate(camera_posn_vect);
    camera_posn_vect += camerapos;

    vtk_cam_rot.rotate(viewupvect);

    renderer->GetActiveCamera()->SetPosition(camera_posn_vect.x,camera_posn_vect.y,camera_posn_vect.z);
    renderer->GetActiveCamera()->SetViewUp(viewupvect.x, viewupvect.y, viewupvect.z);
    renderer->GetActiveCamera()->ComputeViewPlaneNormal();
    renderer->ResetCameraClippingRange();

    //renderer->ResetCamera();

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
     
    if(tag == "VIS_VTK_HIGHQUAL") 
	{ 
	    ss >> value;
	    VTK_HIGHQUAL = getBoolOpt(value);
	    continue;
	}

	//if(tag == "VIS_INTERACTIVE") 
	//{ 
	//    ss >> value;
	//    OptsInteractive = getBoolOpt(value);
	//    continue;
	//}
      
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

    if(tag == "VIS_CAMERADISTMULT") 
	{
	    ss >> OptsCameraDistMult;
	    continue;
	}
    if(tag == "VIS_LINETHICKNESS") 
	{
	    ss >> VIS_LINETHICKNESS;
	    continue;
	}
    if(tag == "VIS_SKIPOUTOFFOCUS") 
        {
	    ss >> value;
	    OptsSkipOutOfFocusPoints = getBoolOpt(value);
	    continue;
	} 
    if(tag == "VIS_PARALLELPROJECTION") 
        {
	    ss >> value;
	    VIS_PARALLELPROJECTION = getBoolOpt(value);
	    continue;
	} 
    if(tag == "VIS_USENUCTEXMAP") 
        {
	    ss >> value;
	    OptsUseNucTextureMap = getBoolOpt(value);
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
    cout << "OptsInteractive          = " << POST_VTK_VIEW          << endl;
    cout << "OptsRenderNucleator      = " << OptsRenderNucleator      << endl;
    cout << "OptsRenderNodes          = " << OptsRenderNodes          << endl;
    cout << "OptsRenderLinks          = " << OptsRenderLinks          << endl;
    cout << "OptsShadeLinks           = " << OptsShadeLinks           << endl;
    cout << "OptsVolumeRenderNodes    = " << OptsVolumeRenderNodes    << endl;
    cout << "OptsIsoRenderNodes       = " << OptsIsoRenderNodes       << endl;
    cout << "OptsRenderAxes           = " << OptsRenderAxes           << endl;
    cout << "OptsRenderText           = " << OptsRenderText           << endl;
    cout << "OptsRenderProjection     = " << OptsRenderProjection     << endl;
    cout << "OptsNormaliseFrames      = " << OptsNormaliseFrames      << endl;
    cout << "voxel_scale              = " << voxel_scale              << endl;
    cout << "ni nj nk                 = " << ni<<" "<<nj<<" "<< nk    << endl;
    cout << "vx intensity scale       = " << vx_intensity_scale       << endl;
    cout << "file_prefix              = " << file_prefix              << endl;
    cout << "projection_scale         = " << p_scale                  << endl;
    cout << "CameraDistMult           = " << OptsCameraDistMult       << endl;
    cout << "OptsSkipOutOfFocusPoints = " << OptsSkipOutOfFocusPoints << endl;
    cout << "OptsUseNucTextureMap = " << OptsUseNucTextureMap << endl;
}

#endif

