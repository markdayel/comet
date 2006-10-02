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

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkStructuredPoints.h"

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

  rotationmatrix vtk_cam_rot;

  vtkRenderer *renderer;
  vtkRenderWindow *render_win;
  vtkRenderWindowInteractor *iren;

  double radius_pixels;
  
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
  double OptsCameraDistMult;
  OptProjectionType OptsRenderProjection;
  bool OptsSkipOutOfFocusPoints;
  
  bool getBoolOpt(const std::string &value);
  OptProjectionType getProjectionOpt(const std::string &value);
  
 public:
  CometVtkVis();//actin * theactin);
  ~CometVtkVis();
  
  void buildVTK(int framenumber);
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
  
  void setOptions();
  void reportOptions();
  void createGaussianVoxelRepresentation(vtkStructuredPoints *vx);
  void setProjection();
  void fillVoxelSetFromActinNodes(vector< vector< vector<double > > >  &vx);// ,double vd, double* min);
  void saveImage(char *filename);
  void renderProjections(vtkRenderWindow *render_win, string base_filename);
  double getMeanNodeLinkForce(const int node_i);
};


