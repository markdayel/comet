PROG = /cluster/bin/comet
SRCS = actin.cpp Colour.cpp links.cpp \
       nodes.cpp nucleator.cpp rotationmatrix.cpp \
       segments.cpp kbhit.cpp threadedtaskqueue.cpp \

VTKSRCS = comet_vtk.cpp comet.cpp 

OBJS = $(SRCS:.cpp=.o)
VTKOBJS = $(VTKSRCS:.cpp=.o)
ALLOBJS = $(OBJS) $(VTKOBJS) 

# LIBS = -lm -I/usr/local/include/stlport -L/usr/local/lib/stlport
#LIBS = -lm -lpthread -L"/System/Library/Frameworks/OpenGL.framework/Libraries" -framework AppKit -framework OpenGL -lGL -lGLU

# vtk 
#VTKINCLUDES = -I/opt/local/include/vtk/
#VTKLIBPATH  = -L/opt/local/lib/vtk/

#VTKINCLUDES = -I/Users/mark/VTKcvs/VTKBuild/include/vtk-5.1/
#VTKLIBPATH  = -L/Users/mark/VTKcvs/VTKBuild/lib/

# libraries to use for vtk link (and also set -DLINK_VTK)

#LIBS = -lpthread -lgsl -lgslcblas -framework AppKit -framework OpenGL
#VTKLIBPATH = /usr/lib
#VTKINCLUDES = -I/Users/mark/VTKcvs/VTKBuild.compat/include/vtk-5.1/
#VTKLIBPATH  = -L/Users/mark/VTKcvs/VTKBuild.compat/lib/

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

CXXFLAGS =  -O3 \
            -ffast-math -mfpmath=387 -mfpmath=sse \
            -Wall -Wno-deprecated -Wextra \
            -g \
            -DNDEBUG

#-march=amdfam10 -mtune=amdfam10 \
#CXXFLAGS =  -O3 \
#            -march=amdfam10 -mtune=amdfam10 \
#            -Wall -Wno-deprecated -Wextra \
#            -ftree-vectorize\
#            -g \
#            -DNDEBUG\
#            -fprefetch-loop-arrays 

#            -mfpmath=sse -msse2 -msse3 \
#\
#            -L/opt/local/lib/ -I/opt/local/include/ 

#	    -mdynamic-no-pic -fno-pic\
#            -fasm-blocks\


#CXXFLAGS =  -march=nocona -mtune=nocona \
#            -mfpmath=sse -msse2 -msse3 \
#            -Wall -Wno-deprecated -Wextra \
#            -g \
#	    -DLINK_VTK \ 
#            -arch x86_64 \


#            -m64

# -arch x86_64 is probably the way to go with 64 bit once it's supported in os x
# -m64 makes things ~20% faster, but can't complie  vtk in 64 bit 
# on Tiger because Tiger libraries are 32 bit

# -fomit-frame-pointer
# -g

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

