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

#include <string>
#include "stdafx.h"

#include "vtkVersion.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkStructuredPoints.h"

#include "vtkPolyDataMapper.h"

#include "vtkJPEGReader.h"

//#include "vtkMesaRenderer.h"
//#include "vtkMesaRenderWindow.h"

//#include "vtkOpenGLRenderWindow.h"
//#include "vtkOpenGLRenderer.h"

// vtk 4.2 libaries need floats not doubles
#if VTK_MAJOR_VERSION < 5
  #define VTK_FLOAT_PRECISION float
#else
  #define VTK_FLOAT_PRECISION double
#endif

//#define vtkRenderWindow vtkOpenGLOffscreenRenderWindow
//#define vtkRenderer vtkOpenGLRenderer

class CometVtkVis {
 private:
  //actin* ptheactin;
  enum OptProjectionType {X, Y, Z, RIP, NONE};
  double voxel_scale; // const
  double p_scale;
  int ni, nj, nk;
  double nuc_opacity;
  int renderwin_npx;
  int renderwin_npy;
  double vx_intensity_scale;
  std::string file_prefix;

  VTK_FLOAT_PRECISION meanx, meany, meanz;

  int movex,movey,movez;

  rotationmatrix vtk_cam_rot;

  vtkRenderer *renderer;
  vtkRenderWindow *render_win;
  vtkRenderWindowInteractor *iren;

  vtkJPEGReader *tx_reader;

  double radius_pixels, capsule_length_pixels;
  
  bool OptsInteractive;
  bool OptsRenderNucleator;
  bool OptsRenderNodes;
  bool OptsRenderLinks;
  bool OptsShadeLinks;
  bool OptsVolumeRenderNodes;
  bool OptsIsoRenderNodes;
  bool OptsRenderAxes;
  bool OptsRenderText;
  bool OptsNormaliseFrames;
  bool OptsUseNucTextureMap;
  double OptsCameraDistMult;
  OptProjectionType OptsRenderProjection;
  bool OptsSkipOutOfFocusPoints;
  bool VTK_HIGHQUAL;
  double VIS_LINETHICKNESS;
  bool VIS_PARALLELPROJECTION;

  bool OptVIEW_VTK;

  bool texturereadOK;
  
  bool getBoolOpt(const std::string &value);
  OptProjectionType getProjectionOpt(const std::string &value);
  
 public:
  CometVtkVis(bool VIEW_VTK);//actin * theactin);
  ~CometVtkVis();
  
  char VTK_colmap_filename[1024];

  bool convert_to_vtkcoord(VTK_FLOAT_PRECISION &x, VTK_FLOAT_PRECISION &y, VTK_FLOAT_PRECISION &z);

  void set_mean_posns();

  void buildVTK(int framenumber, vect & cameraposition, vect & cameratarget);

  void RestartRenderWindow();
  void addNucleator();
  void addCapsuleNucleator();
  void addSphericalNucleator();
  void addNodes();
  void addLinks();
  void addGaussianSplatterToNodes();
  void addGuassianNodeVolume(bool do_iso);
  void addLight();
  void addVolumeTest();
  void addStructuredPointVolumeRender(vtkStructuredPoints *vx);
  void addStructuredPointIsoRender(vtkStructuredPoints *sp);
  void addAxes();
  void addVoxelBound();

  void SetFocalDepthPlanes(vtkPolyDataMapper *map);
  void saveVRML(int framenumber);
  
  void setOptions();
  void reportOptions();
  void createGaussianVoxelRepresentation(vtkStructuredPoints *vx);
  void setProjection();
  void setProjection(vect & cameraposition,vect & cameratarget);
  void fillVoxelSetFromActinNodes(vector< vector< vector<double > > >  &vx);// ,double vd, double* min);
  void saveImage(int framenumber);
  void renderProjections(vtkRenderWindow *render_win, string base_filename);
  double getMeanNodeLinkForce(const int node_i);
};


