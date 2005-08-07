PROG = comet
SRCS = actin.cpp Colour.cpp comet.cpp links.cpp nodes.cpp nucleator.cpp stdafx.cpp rotationmatrix.cpp
OBJS = $(SRCS:.cpp=.o)
LIBS =  -lm -lpthread
#CXXFLAGS = -O5 -ffast-math -fsingle-precision-constant -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays
#CXXFLAGS = -O5 -ffast-math -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays -g -Wall
#CXXFLAGS = -O5 -ffast-math -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays -Wall
#CXXFLAGS = -O5 -march=athlon-xp -mfpmath=sse -fmove-all-movables \
#                -freduce-all-givs -Winline -Wall -g
#CXXFLAGS = -O5 -march=athlon-xp -ffast-math -mfpmath=sse -Wall -g 
#CXXFLAGS = -O3 -ffast-math -Wall -g -DSEED_INSIDE
CXXFLAGS = -O3 -ffast-math -Wall -g
#CXXFLAGS = -Wall
CXX = c++

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o ../bin/$(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(PROG)

dataclean:
	rm -f ../bin/*.png
	rm -f ../bin/*.gz
	rm -f ../bin/*.txt
	rm -f ../bin/*.wrz
	rm -f ../bin/*.bmp
	rm -f ../bin/*.wrl


