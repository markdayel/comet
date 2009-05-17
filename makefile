PROG = /cluster/bin/comet
SRCS = actin.cpp Colour.cpp links.cpp \
       nodes.cpp nucleator.cpp rotationmatrix.cpp \
       segments.cpp kbhit.cpp threadedtaskqueue.cpp \

VTKSRCS = comet_vtk.cpp comet.cpp 

OBJS = $(SRCS:.cpp=.o)
VTKOBJS = $(VTKSRCS:.cpp=.o)
ALLOBJS = $(OBJS) $(VTKOBJS) 

# note: '-framework AppKit -framework OpenGL' are required for vtk on os x

# LIBS = -lm -I/usr/local/include/stlport -L/usr/local/lib/stlport
#LIBS = -lm -lpthread -L"/System/Library/Frameworks/OpenGL.framework/Libraries" -framework AppKit -framework OpenGL -lGL -lGLU

# vtk 
#VTKINCLUDES = -I/opt/local/include/vtk/
#VTKLIBPATH  = -L/opt/local/lib/vtk/

#VTKINCLUDES = -I/Users/mark/VTKcvs/VTKBuild/include/vtk-5.1/
#VTKLIBPATH  = -L/Users/mark/VTKcvs/VTKBuild/lib/

# libraries to use for vtk link (also set -DLINK_VTK in CXXFLAGS)
# assumes vtk libraries etc built with cmake and in /Users/mark/VTKcvs/VTKBuild.compat/
#LIBS = -lpthread -lgsl -lgslcblas -framework AppKit -framework OpenGL
#VTKLIBPATH = /usr/lib
#VTKINCLUDES = -I/Users/mark/VTKcvs/VTKBuild.compat/include/vtk-5.1/
#VTKLIBPATH  = -L/Users/mark/VTKcvs/VTKBuild.compat/lib/

# vtk with ubuntu (this doesn't work---default vtk package install is missing libraries---need to use cmake)
#LIBS = -lpthread -lgsl -lgslcblas
#VTKINCLUDES = -I/usr/include/vtk-5.0/
#VTKLIBS     = -lvtkRendering -lvtkImaging -lvtkCommon -lvtkGraphics \
              -lvtkpng -lvtkIO -lvtkexpat -lvtkFiltering -lvtkftgl \
              -lvtkHybrid -lvtkVolumeRendering -lvtkjpeg -lvtktiff -lvtksys -lvtkzlib#

# libraries to use for non-vtk (and also unset -DLINK_VTK)

LIBS = -lpthread  -lgsl -lgslcblas 
VTKLIBS = 
VTKLIBPATH =
VTKINCLUDES =

CXXFLAGS =  -O3 -march=native -mtune=native\
            -ffast-math \
            -Wall -Wextra \
            -DNDEBUG -g

# N.B use -DLINK_VTK if comiling with VTK
#	    -DLINK_VTK \ 

# -arch x86_64 is probably the way to go with 64 bit once it's supported in os x
# -m64 makes things ~20% faster, but can't complie  vtk in 64 bit 
# on Tiger because Tiger libraries are 32 bit


CXX = c++

$(PROG):$(ALLOBJS)
	$(CXX) $(CXXFLAGS) $(VTKLIBPATH) -o $(PROG) $(ALLOBJS) $(LIBS) $(VTKLIBS)


comet.o:comet.cpp
	$(CXX) $(FRAMEWORK) $(CXXFLAGS) $(VTKINCLUDES) -c -o $@ $< 

comet_vtk.o:comet_vtk.cpp
	$(CXX) $(FRAMEWORK) $(CXXFLAGS) $(VTKINCLUDES) -c -o $@ $< 

clean:
	rm -f $(ALLOBJS)

distclean: clean
	rm -f $(PROG)

dataclean:
	rm -f ../bin/*.png
	rm -f ../bin/*.gz
	rm -f ../bin/*.txt
	rm -f ../bin/*.wrz
	rm -f ../bin/*.bmp

