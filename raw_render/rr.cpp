// First include the required header files for the VTK classes we are using.
#include "vtkConeSource.h"
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


#define RAWTYPE float

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
    
    vtkStructuredPoints *spoints = vtkStructuredPoints::New();
    double vspacing[3];
    vspacing[0] = nm_per_slice;
    vspacing[1] = nm_per_pixel;
    vspacing[2] = nm_per_pixel;
    createStructuredPointRepresentation(raw_data, vspacing, spoints);
    
    vtkRenderer *ren = vtkRenderer::New();
    ren->SetBackground( 1.0, 1.0, 1.0 );

    if(render_iso) {
	cout << "Rendering iso surface, threshold:" 
	     << isothreshold << endl;
	addIsoRender(ren, spoints, isothreshold, 1.0 );
    }
    if(render_vol) {
	cout << "Rendering volume." << endl; 
	addVolumeRender(ren, spoints);  
    }
    
    vtkRenderWindow *renWin = vtkRenderWindow::New();
    renWin->AddRenderer( ren );
    renWin->SetSize( 400, 400 );
    
    vtkRenderWindowInteractor *iren =  vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    vtkCallbackCommand *keypress = vtkCallbackCommand::New();
    keypress->SetCallback(renderCallback);
    iren->AddObserver(vtkCommand::KeyPressEvent,keypress);

    cout << "Rendering... press 's' to grab image." << endl;     
    iren->Initialize();
    iren->Start();

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
    int offset = 0;
    for(int k=0; k<nk; ++k) {	
	for(int j=0; j<nj; ++j) {	
	    for(int i=0; i<ni; ++i) {
		double val = raw_data[i][j][k];
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
    smooth->SetNumberOfIterations( 20 );
    //smooth->SetFeatureAngle(45);
    smooth->BoundarySmoothingOn();  
    
    vtkPolyDataNormals *normals = vtkPolyDataNormals::New();
    normals->SetInput( smooth->GetOutput());
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
    volume_map->SetSampleDistance(100); // important
    
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
    
    ifstream raw_info("raw/frameinfo.txt", ios::in);
    
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
	sprintf(filename, "raw/%s", buff);
	
	cout << "Reading file " << filename << "...";
	cout.flush();
	
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
		cout << endl 
		     << "Warning: Finished reading " 
		     <<  file_size - posn 
		     << " bytes short of end of file " << endl;
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
	    }
	}
    }
    cout << "min val:"<< min_val << ", "
	 << "max val:"<< max_val << endl;
}
