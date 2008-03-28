#include "vtkPolyDataMapper.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkProperty.h"

#include "vtkPointData.h"
#include "vtkDataSet.h"
#include "vtkStructuredPoints.h"
#include "vtkLookupTable.h"
#include "vtkPiecewiseFunction.h"
#include "vtkColorTransferFunction.h"
#include "vtkVolumeProperty.h"
#include "vtkVolumeRayCastCompositeFunction.h"
#include "vtkVolumeRayCastMIPFunction.h"
#include "vtkVolumeRayCastMapper.h"
#include "vtkImageCast.h"
#include "vtkShepardMethod.h"

#include "vtkRenderLargeImage.h"
#include "vtkJPEGWriter.h"
#include "vtkPNGWriter.h"
#include "vtkBMPWriter.h"

// iso
#include "vtkContourFilter.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkPolyDataNormals.h"

// call back
#include "vtkCallbackCommand.h"

// file writer
#include "vtkJPEGWriter.h"
#include "vtkWindowToImageFilter.h"

#include <vector>
using namespace std;

#define PI 3.14159

#define RAWTYPE float

// REVISIT, move to header file
void createStructuredPointRepresentation(const vector<vector<vector<RAWTYPE> > > 
					 &raw_data, 
					 const double vspacing[3],
					 vtkStructuredPoints *sp);
void addVolumeRender(vtkRenderer *renderer,
		     vtkStructuredPoints *vx);
void addIsoRender(vtkRenderer *renderer, vtkStructuredPoints *sp, 
		  const double threshold, const double opacity);
int loadRaw(std::vector<std::vector<std::vector<RAWTYPE > > > &raw_data,
	    double &nm_per_pixel, double &nm_per_slice);
void getMinMax(const std::vector<std::vector<std::vector<RAWTYPE> > > 
	       &raw_data, 
	       double &min, double &max);
void renderBead(bool render_vol, bool render_iso, double isothreshold);
void filterRaw(std::vector<std::vector<std::vector<RAWTYPE > > > &raw_data,
	       int ksz);

// ---
// ---
void saveCurrentImage(vtkRenderWindow* ren_win)
{
   // might want to use large image methods?
   vtkWindowToImageFilter *to_image = vtkWindowToImageFilter::New();
   vtkJPEGWriter *writer = vtkJPEGWriter::New();
   
   to_image->SetInput(ren_win);
   writer->SetInput(to_image->GetOutput());
   ren_win->Render();
   to_image->Modified();
   writer->SetFileName("grab.jpeg");
   writer->Write();
   cout << "wrote current window to 'grab.jpeg'" << endl;
   
   writer->Delete();
   to_image->Delete();
}    

void saveCurrentImage(vtkRenderWindow* ren_win, int framenum)
{
   // might want to use large image methods?
   vtkWindowToImageFilter *to_image = vtkWindowToImageFilter::New();
   vtkJPEGWriter *writer = vtkJPEGWriter::New();
   
   char filename[1024];
   sprintf(filename  , "vtk_%05i.jpg",  framenum );

   to_image->SetInput(ren_win);
   writer->SetInput(to_image->GetOutput());
   ren_win->Render();
   to_image->Modified();
   writer->SetFileName(filename);
   writer->Write();
   
   writer->Delete();
   to_image->Delete();
}

void saveBigImage(vtkRenderer* ren, int framenum)
{

    int MAGFACTOR=2;

    char filename[1024];
	char tmpfilename[1024];
	
	

	#define IMAGEMAGICKCONVERT "convert"
    #define BMP_COMPRESSION 95

    sprintf(filename  , "vtk_%05i.jpg",  framenum );
	sprintf(tmpfilename  , "tmp_vtk_%05i.jpg",  framenum );

    //cout << "Saving " << filename << endl;
  

    //vtkWindowToImageFilter *rwin_to_image = vtkWindowToImageFilter::New();
    //rwin_to_image->SetInput(render_win);

    vtkRenderLargeImage *renderLarge = vtkRenderLargeImage::New();
    renderLarge->SetInput(ren);
    renderLarge->SetMagnification(MAGFACTOR);



    // vtkBMPWriter seems a bit faster than vtkPNGWriter, and we're compressing with IM to png anyway later
    // so use BMP here
    vtkBMPWriter *imagewriter = vtkBMPWriter::New();
    
    imagewriter->SetInput( renderLarge->GetOutput() ); //rwin_to_image->GetOutput() );
	if (MAGFACTOR==1)    
		imagewriter->SetFileName( filename );
	else
		imagewriter->SetFileName( tmpfilename );

	//renderLarge->Modified();
	//renderLarge->Update();

    imagewriter->Write();    
    
    imagewriter->Delete();

    //rwin_to_image->Delete();
    renderLarge->Delete();

    if (MAGFACTOR!=1)
    {
        char command1[1024];
        sprintf(command1, 
            "(%s -quality %i -resize %f%% %s %s ; rm %s ) &",
            IMAGEMAGICKCONVERT, BMP_COMPRESSION, 100/(double)MAGFACTOR,  
              tmpfilename, filename, tmpfilename);
        //cout << command1 << endl;
        system(command1);
    }


}

void renderCallback(vtkObject *caller, unsigned long eid, 
		    void * cldata, void* cdata)
{
    vtkRenderWindowInteractor *iren = 
    	reinterpret_cast<vtkRenderWindowInteractor *>(caller);

    vtkRenderWindow* ren_win = iren->GetRenderWindow();
    
    if(iren->GetKeyCode()=='s' || iren->GetKeyCode()=='S' )
	saveCurrentImage(ren_win);
}

int main(int argc, char *argv[])
{
    bool render_vol = false;
    bool render_iso = false;    
    double isothreshold = 0.5;
    
    if(argc<2){
	cout << "Usage: " << endl
	     << "   " << argv[0]
	     << " (vol) (iso iso_threshold) " << endl
	     << "    where iso_threshold is 0-1"
	     << " eg 'vol iso 0.5'." <<  endl;
	exit(EXIT_FAILURE);
    }

    if (strcmp(argv[1], "vol") == 0 )
    {
	render_vol = true;
    }
    for(int i=0; i<2; ++i)
    {
	if( (i+1)<argc && strcmp(argv[i+1], "iso") == 0 )
	{
	    render_iso = true;
	    if( (i+2)<argc )
		isothreshold = atof(argv[i+2]);
	}
    }
    renderBead(render_vol, render_iso, isothreshold);
}

void renderBead(bool render_vol, bool render_iso, double isothreshold)
{
    std::vector<std::vector<std::vector<RAWTYPE> > > raw_data;
    double nm_per_pixel, nm_per_slice;

    loadRaw(raw_data, nm_per_pixel, nm_per_slice);    

    filterRaw(raw_data, 3);

    // convert to structured points to render via vtk
    vtkStructuredPoints *spoints = vtkStructuredPoints::New();
    double vspacing[3];
    vspacing[0] = nm_per_slice;
    vspacing[1] = nm_per_pixel;
    vspacing[2] = nm_per_pixel;

	double x_center = (nm_per_slice * raw_data.size()) / 2.0;
	double y_center = (nm_per_pixel * raw_data[0].size()) / 2.0;
	double z_center = (nm_per_pixel * raw_data[0][0].size()) / 2.0;

    createStructuredPointRepresentation(raw_data, vspacing, spoints);
    
    // create render context
    vtkRenderer *ren = vtkRenderer::New();
    ren->SetBackground( 0.0, 0.0, 0.0 );

    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->SetSize( 800, 800 );
	renWin->AddRenderer( ren );

	renWin->Modified();

    if(render_iso) {
	cout << "Rendering iso surface, threshold:" 
	     << isothreshold << endl;
	addIsoRender(ren, spoints, isothreshold, 1.0 );
    }
    if(render_vol) {
	cout << "Rendering volume." << endl; 
	addVolumeRender(ren, spoints);  
    }
    
	//ren->SetAmbient(0.8,0.8,0.8);
    
    vtkRenderWindowInteractor *iren =  vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

	//renWin->Render();
    

    //vtkCallbackCommand *keypress = vtkCallbackCommand::New();
    //keypress->SetCallback(renderCallback);
    //iren->AddObserver(vtkCommand::KeyPressEvent,keypress);
    //cout << "Rendering... press 's' to grab image." << endl;     
    
	//iren->Initialize();
    //iren->Start();
	//iren->Enable();
	//ren->Modified();
	//iren->Modified();
	//iren->Render();
	
	
    ren->GetActiveCamera()->SetFocalPoint(x_center,y_center,z_center);
	
	ren->GetActiveCamera()->SetViewUp(1,0,0);
	
	ren->GetActiveCamera()->Modified();
	ren->GetActiveCamera()->ParallelProjectionOff();
	
	ren->GetActiveCamera()->ViewingRaysModified();
	
	//ren->GetActiveCamera()->SetViewAngle(60);
	
	double steps = 180;
	
	double camerax = x_center*2;
	double cameralen = y_center * 2.5;
	
	double step = 360 / steps;
	
	
	for (int frame = 1; frame <= steps ; ++frame)
	{
		double angle = 2 * PI * step * (double) frame / 360.0;
		
		ren->GetActiveCamera()->SetPosition(x_center + cameralen * cos(angle), y_center + cameralen * cos(angle), z_center + cameralen * sin(angle)  );
	
		ren->GetActiveCamera()->ComputeViewPlaneNormal();
		ren->ResetCameraClippingRange();

		cout << "Rendering Frame " << frame << " of " << steps << endl;
		cout.flush();

		//saveBigImage(ren, frame);
		saveCurrentImage(renWin, frame);
	
	}

    spoints->Delete();       
    ren->Delete();
    renWin->Delete();
    iren->Delete();
}

void createStructuredPointRepresentation(
    const vector<vector<vector<RAWTYPE> > > &raw_data,
    const double vspacing[3],
    vtkStructuredPoints *sp)
{
    // create the voxel set as vtk structured point data
    // set up the spoints for volume rendering
    int ni, nj, nk;
    ni = raw_data.size(); // z
    nj = raw_data[0].size(); // height
    nk = raw_data[0][0].size(); // width

    int dim[3] = {ni, nj, nk};
    sp->SetDimensions(dim);
    sp->SetSpacing(vspacing[0], 
		   vspacing[1],
		   vspacing[2]);
    
    // VTK_FLOAT_PRECISION 
    double vx_orig[3] = {-ni/2.0-1, -nj/2.0-1, -nk/2.0-1};  
    sp->SetOrigin( vx_orig ); 
    
    // create the underlying voxel data (only UnsignedChar or Short)
    // see vtk examples for voxel data
    vtkUnsignedCharArray *vxdata = vtkUnsignedCharArray::New();    
    double iscale = 1.0;
    
    double max_val, min_val;
    getMinMax(raw_data, min_val, max_val);
    
    // loop over voxel set converting to structured points suitable for 
    // vtk volume rendering
    // NB: min val is assumed 0
    int offset = 0;
    for(int k=0; k<nk; ++k) {	
	for(int j=0; j<nj; ++j) {	
	    for(int i=0; i<ni; ++i) {
		double val = fabs( raw_data[i][j][k] - 1.5 * min_val) ;
		int vval = int(255*val/max_val);
		vxdata->InsertValue(offset, int(255*val/max_val));	  
		offset++;		
	    }
	}
    }
    sp->GetPointData()->SetScalars(vxdata);
    vxdata->Delete();  
}

void addIsoRender(vtkRenderer *renderer,
		  vtkStructuredPoints *sp, 
		  const double threshold, 
		  const double opacity)
{
    // render using an iso-surface
    vtkContourFilter *iso_surf = vtkContourFilter::New();
    iso_surf->SetInput( sp );
    iso_surf->SetValue(0, threshold*255);
    iso_surf->ComputeNormalsOn();
    
    vtkWindowedSincPolyDataFilter *smooth 
	= vtkWindowedSincPolyDataFilter::New();
    smooth->SetInput(iso_surf->GetOutput());
    smooth->SetNumberOfIterations( 40 );
    smooth->SetFeatureAngle(45);
    smooth->BoundarySmoothingOn();  
    
    vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
    normals->SetInput( smooth->GetOutput());
    //normals->SetInput( iso_surf->GetOutput());
    //normals->SetFeatureAngle(45);
    //normals->FlipNormalsOn();
   
    
    vtkPolyDataMapper *iso_mapper = vtkPolyDataMapper::New();
    iso_mapper->SetInput( normals->GetOutput() );
    iso_mapper->ScalarVisibilityOff();
    iso_surf->Delete();
    
    vtkActor *iso_actor = vtkActor::New();
    iso_actor->SetMapper(iso_mapper);
    iso_mapper->Delete();
    
    iso_actor->GetProperty()->SetOpacity(opacity);
    iso_actor->GetProperty()->SetColor(0.3, 0.6, 0.3);
	iso_actor->SetPosition(0,0,0);

	//iso_actor->SetOrientation(PI/2,PI/2,PI/2);

    renderer->AddActor(iso_actor);
    iso_actor->Delete();
    smooth->Delete();
    normals->Delete();
}

void addVolumeRender(vtkRenderer *renderer, vtkStructuredPoints *vx)
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
    volume_map->SetSampleDistance(50); // important
    
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

void writeKernel(vector<vector<vector<RAWTYPE> > > &kernel)
{
  int ni = kernel.size(); // z
  int nj = kernel[0].size(); // height
  int nk = kernel[0][0].size(); // width
  int mj = (nj-1)/2;
  int mk = (nk-1)/2;
  std::ofstream os("k.dat" );
  for(int i=0; i<ni; ++i){
    os << i << " " << kernel[i][mj][mk] << endl;
  }
  os.close();
  
}

void gaussianKernel(vector<vector<vector<RAWTYPE> > > &kernel, int sz)
{
    // REVISIT: check and adjust values to get desired extent
    int ex = sz*2+1;
    double sigma  = ex/(2*2.3); // FWHM
    
    kernel.resize(ex);    
    for(int i=0; i<ex; i++) {
	kernel[i].resize(ex);	
	for(int j=0; j<ex; j++){
	    kernel[i][j].resize(ex);	    
	}
    }
    // loop calculating kernel value
    double x2;
    double ksum = 0;
    for(int i=0; i<ex; i++) {
	for(int j=0; j<ex; j++) {
	    for(int k=0; k<ex; k++) {
		x2 = ( (i-sz)*(i-sz) + (j-sz)*(j-sz) + (k-sz)*(k-sz) );
		kernel[i][j][k] = (1/(sigma*sqrt(2*M_PI)))
		    *exp( -(x2)/(2.0*sigma*sigma) );
		ksum += kernel[i][j][k];
	    }
	}
    }
    // normalise
    for(int i=0; i<ex; i++) {
	for(int j=0; j<ex; j++) {
	    for(int k=0; k<ex; k++) {
		kernel[i][j][k] = kernel[i][j][k]/ksum;
	    }
	}
    }
    writeKernel(kernel);
    
}

RAWTYPE medianFilter(vector<RAWTYPE> &kernel_data)
{
    sort(kernel_data.begin(), kernel_data.end());
    return kernel_data[(kernel_data.size()/2) + 1];
}

// ksz is kernel size, ie [-ksz..0..+ksz]
void filterRaw(std::vector<std::vector<std::vector<RAWTYPE > > > &raw_data,
	       int ksz)
{
    int nkn = ksz*2 + 1; // kernel is odd
    cout << "extent of the kernel, nkn: " << nkn << endl;

    vector<RAWTYPE> kernel_data;
    kernel_data.resize(nkn*nkn*nkn);
    
    int ni = raw_data.size(); // z
    int nj = raw_data[0].size(); // height
    int nk = raw_data[0][0].size(); // width
    
    // create temporary volume for filtered data
    std::vector<std::vector<std::vector<RAWTYPE> > > filtered_data;
    filtered_data.resize(ni);    
    for(int i=0; i<ni; i++) {
	filtered_data[i].resize(nj);	
	for(int j=0; j<nj; j++){
	    filtered_data[i][j].resize(nk);	    
	}
    }
    
    // create kernel
    std::vector<std::vector<std::vector<RAWTYPE> > > kernel;
    gaussianKernel(kernel, ksz);

    // apply kernel
    // loop over raw data filling filtered data
    for(int i=0; i<ni; i++) {
	for(int j=0; j<nj; j++){
	    for (int k=0; k<nk; k++){

	 
		double fd = 0;
		for(int ki=0; ki<kernel.size(); ki++) {
		    for(int kj=0; kj<kernel[0].size(); kj++){
			for (int kk=0; kk<kernel[0][0].size(); kk++){
			    int dki= (i-ksz)+ki;
			    int dkj= (j-ksz)+kj;			    
			    int dkk= (k-ksz)+kk;
			    
			    if(dki<0)
				dki=0;
			    else if(dki>(ni-1))
				dki = ni-1;
			    
			    if(dkj<0)
				dkj=0;
			    else if(dkj>(nj-1))
				dkj = nj-1;
			    
			    if(dkk<0)
				dkk=0;
			    else if(dkk>(nk-1))
				dkk = nk-1;

			    // kernel_data[kdi++] = raw_data[dki][dkj][dkk];
			    fd += raw_data[dki][dkj][dkk]*kernel[ki][kj][kk];
			}
		    }
		}
		// filter the kernel data
		// filtered_data[i][j][k] = medianFilter(kernel_data);
		filtered_data[i][j][k] = fd;
	    }
	}
    }

    // copy filtered data back to raw data
    for(int i=0; i<ni; i++) {
	for(int j=0; j<nj; j++){
	    for (int k=0; k<nk; k++){
		raw_data[i][j][k] = filtered_data[i][j][k];
	    }
	}
    }


    return;
}

int loadRaw(std::vector<std::vector<std::vector<RAWTYPE > > > &raw_data,
	    double &nm_per_pixel, double &nm_per_slice)
{   
    // bit ugly right now---just read in images from 'raw' dir  
    int raw_width, raw_height;
    int raw_firstframe, raw_lastframe;
    int totframes;
    char raw_filepattern[2048], filename[2048], buff[2048];
    int filenum;
    
    // read data settings in from raw files instead 
    // of previously calculated data
    
    ifstream raw_info("frameinfo.txt", ios::in);
    
    if (!raw_info) { 
	cout << "Unable to open file 'frameinfo.txt' for input" << endl;
	exit(EXIT_FAILURE);
    } else {
	std::string line;
	getline(raw_info,line);
	//cout << line;
	raw_info    >> raw_width 
		    >> raw_height
		    >> nm_per_pixel
		    >> nm_per_slice
		    >> raw_firstframe
		    >> raw_lastframe
		    >> raw_filepattern;
	
	raw_info.close();  
	
	cout << setprecision(2);
	cout << "File info: "
	     << "Width:       " << raw_width << endl
	     << "Height:      " << raw_height << endl
	     << "nm/pixel:    " << nm_per_pixel << endl
	     << "nm/slice:    " << nm_per_slice << endl
	     << "First frame: " << raw_firstframe << endl
	     << "Last Frame:  " << raw_lastframe << endl
	     << "Filepattern: " << raw_filepattern << endl;
	
	
	if ((( raw_width < 1 ) || (raw_width > 2000)) ||
	    (( raw_height < 1 ) || (raw_height > 2000)) || 
	    (( raw_lastframe < 1 ) || (raw_lastframe > 2000))) {
	    cout << " Parameters out of range.  Aborting read. " << endl;
	    exit(EXIT_FAILURE);
	}
	
	cout << endl;
	
	
	totframes = (raw_lastframe - raw_firstframe) + 1;
	
	cout << "Raw type is set to " << sizeof(RAWTYPE) * 8 
	     << " bits per pixel" << endl;
	
	double memreq=sizeof(RAWTYPE) * raw_width * raw_height * raw_lastframe;
	
	cout << "Allocating memory required for image : " 
	     <<  ((memreq / 1024) / 1024) << " MB" << endl;
	
	
	raw_data.resize(totframes);
	
	for (int z=0; z!=totframes; z++)
	{
	    raw_data[z].resize(raw_height);
	    
	    for (int x=0; x!=raw_height; x++)  // allocate data grid
	    {
		raw_data[z][x].resize(raw_width);
		
		//for (int y=0; y!=raw_height; y++)
		//{
		//    raw_data[z][x][y].resize(raw_lastframe);
		//}
	    }
	}
	
	cout << "Memory allocated" << endl;
	
    }
    
    for(filenum = raw_firstframe; filenum <= raw_lastframe ; filenum++ ) {
	
	sprintf(buff, raw_filepattern, filenum);
	sprintf(filename, "%s", buff);
	
	//cout << "Reading file " << filename << "...";
	//cout.flush();
	
	ifstream raw_in(filename, ios::in);
	
	// find filesize	
	raw_in.seekg(0,ios::end);
	int file_size = raw_in.tellg();
	raw_in.seekg(0);
	
	int z = filenum - raw_firstframe;  // index from zero	
	if ((z < 0) || (z > raw_lastframe))
	{
	    cout << "Z out of range" << endl;
	    exit(EXIT_FAILURE);
	}
	
	
	for (int y=0; y!=raw_height; y++)  // allocate data grid
	{
	    //for (int y=0; y!=raw_height; y++)
	    //{
	    //for (int z=0; z!=raw_height; z++)
	    //{
	    raw_in.read((char*)(&raw_data[z][y][0]),raw_width*sizeof(RAWTYPE));
	    //}
	    //}
	}	
	int posn = (int)raw_in.tellg();
	if(posn == -1)
	{
	    cout << endl 
		 << "Warning: End of file reached before data read in!" 
		 << endl;
	}
	else
	    if(posn != file_size)
	    {
	//	cout << endl 
	//	     << "Warning: Finished reading " 
	//	     <<  file_size - posn 
	//	     << " bytes short of end of file " << endl;
	    }
	    else
	    {
		cout << "read OK." << endl;
	    }
	
	cout.flush();
	
	raw_in.close();
    }   
}
    
void getMinMax(const vector<vector<vector<RAWTYPE> > > &raw_data, 
	       double &min_val, double &max_val)
{
    min_val = raw_data[0][0][0];
    max_val = raw_data[0][0][0];

    int ni = raw_data.size(); // z
    int nj = raw_data[0].size(); // height
    int nk = raw_data[0][0].size(); // width

    for(int f=0; f<ni; f++)
    {
	for(int i=0; i<nj; i++)
	{
	    for(int j=0; j<nk; j++)
	    {
		if(raw_data[f][i][j]>max_val)
		    max_val = raw_data[f][i][j];
		
		if(raw_data[f][i][j]<min_val)
		    min_val = raw_data[f][i][j];
	    }
	}
    }
    cout << "min val:"<< min_val << ", "
	 << "max val:"<< max_val << endl;
}
