PROG = comet
SRCS = actin.cpp Colour.cpp comet.cpp links.cpp nodes.cpp nucleator.cpp stdafx.cpp 
OBJS = $(SRCS:.cpp=.o)
LIBS =  -lm -lpthread
#CXXFLAGS = -O5 -ffast-math -fsingle-precision-constant -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays
#CXXFLAGS = -O5 -ffast-math -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays -g -Wall
#CXXFLAGS = -O5 -ffast-math -march=athlon-xp -mfpmath=sse -fomit-frame-pointer -malign-double -fprefetch-loop-arrays -Wall
#CXXFLAGS = -O5 -march=athlon-xp -mfpmath=sse -fmove-all-movables \
#                -freduce-all-givs -Winline -Wall -g
CXXFLAGS = -O5 -march=athlon-xp -ffast-math -mfpmath=sse -Wall -g 
#CXXFLAGS = -O5 -march=athlon-xp -ffast-math -mfpmath=sse -Wall
#CXXFLAGS = -Wall
CXX = c++

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(PROG)

