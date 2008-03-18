PROG = rr
VTKSRCS = rr.cpp

VTKOBJS = $(VTKSRCS:.cpp=.o)

LIBS = -lm `Magick++-config --cppflags --cxxflags --ldflags --libs`


# vtk 
VTKINCLUDES = -I/usr/include/vtk-5.0/
VTKLIBPATH  = -L/usr/lib/vtk-5.0/
VTKLIBS     = -lvtkRendering -lvtkImaging -lvtkCommon -lvtkGraphics -lvtkVolumeRendering

CXXFLAGS = -O3 -ffast-math

CXX = g++

$(PROG):$(VTKOBJS)
	$(CXX) $(CXXFLAGS) \
	-o $(PROG) $(VTKOBJS) $(VTKLIBPATH) $(VTKINCLUDES) $(LIBS) \
	$(VTKLIBS)

rr.o:rr.cpp
	$(CXX) $(CXXFLAGS) $(VTKINCLUDES) \
	-c -o $@ $<

