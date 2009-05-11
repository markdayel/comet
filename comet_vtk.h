/*
    comet - an actin-based motility simulator
    Copyright (C) 2005 Mark J Dayel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

//#include <string>
#include "stdafx.h"

#include "Colour.h"
                                                                   
#include "vtkVersion.h"

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkStructuredPoints.h"

#include "vtkPolyDataMapper.h"

#include "vtkJPEGReader.h"

#include "vtkAssembly.h"

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

  bool vtk_dummy_object;

  int movex,movey,movez;

  rotationmatrix vtk_cam_rot;

  vtkRenderer *renderer;
  vtkRenderWindow *render_win;
  vtkRenderWindowInteractor *iren;

  vtkJPEGReader *tx_reader;

  double radius_pixels, capsule_length_pixels;
  double linewidths;
  
  //bool OptsInteractive;    // changed to POST_VTK_VIEW
  bool OptsRenderNucleator;
  bool OptsRenderNodes;
  bool OptsRenderLinks;
  bool OptsShadeLinks;
  bool OptsVolumeRenderNodes;
  bool OptsIsoRenderNodes;
  bool OptsRenderTracks;
  bool OptsRenderText;
  bool OptsNormaliseFrames;
  bool OptsUseNucTextureMap;
  double OptsCameraDistMult;
  OptProjectionType OptsRenderProjection;
  bool OptsSkipOutOfFocusPoints;
  bool VTK_HIGHQUAL;
  double VIS_LINETHICKNESS;
  bool VIS_PARALLELPROJECTION;

  vect viewupvect;

  double voxelscalefactor;

  bool OptVIEW_VTK;

  bool texturereadOK;
  
  bool getBoolOpt(const std::string &value);
  OptProjectionType getProjectionOpt(const std::string &value);
  
 public:
  CometVtkVis(bool VIEW_VTK, bool dummy_vtk);//actin * theactin);
  ~CometVtkVis();
  
  char VTK_colmap_filename[1024];

  bool convert_to_vtkcoord(VTK_FLOAT_PRECISION &x, VTK_FLOAT_PRECISION &y, VTK_FLOAT_PRECISION &z);

  void set_mean_posns();

  void buildVTK(const int &framenumber, vect & cameraposition, vect & cameratarget);

  void RestartRenderWindow();
  void addNucleator(bool wireframe);
  //void addCapsuleNucleator(bool wireframe, double scale);
  //void addSphericalNucleator(bool wireframe, double scale);
  void add_capsule_to_assembly(vtkAssembly* nucleator_actor, bool wireframe);
  void add_sphere_to_assembly(vtkAssembly* nucleator_actor, bool wireframe);
  void addNodes();
  void addLinks();
  void addGaussianSplatterToNodes();
  //void addGuassianNodeVolume(bool do_iso);
  void addLight();
  void addVolumeTest();
  void addStructuredPointVolumeRender(vtkStructuredPoints *vx);
  void addStructuredPointIsoRender(vtkStructuredPoints *sp, const double threshold, const Colour col, const double Opacity);
  void addAxes();
  void addVoxelBound();

  void addTracks(const int &framenumber);
  void addTrackDistances();

  void addflow();

  void set_transform_matrix(vtkMatrix4x4 * vtkmat, const rotationmatrix & rotmat) const;

  void SetFocalDepthPlanes(vtkPolyDataMapper *map);
  void saveVRML(const int &framenumber);
  
  void setOptions();
  void reportOptions();
  void createGaussianVoxelRepresentation(vtkStructuredPoints *vx);
  void setProjection();
  void setProjection(vect & cameraposition,vect & cameratarget);
  void fillVoxelSetFromActinNodes(vector< vector< vector<double > > >  &vx);// ,double vd, double* min);
  void saveImage(const int &framenumber);
  void saveImageRotationSet(const int &framenumber);
  void saveImage(const int &framenumber, char* filename);
  void renderProjections(vtkRenderWindow *render_win, string base_filename);
  double getMeanNodeLinkForce(const int node_i);
};


