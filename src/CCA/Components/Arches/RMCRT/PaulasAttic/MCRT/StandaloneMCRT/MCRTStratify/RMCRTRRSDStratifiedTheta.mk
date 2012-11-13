CXX = g++
# CFLAGS = -g -MD
CFLAGS = -pg -O3 -MD

SRCS = RMCRTRRSDStratifiedTheta.cc Surface.cc RealSurface.cc TopRealSurface.cc BottomRealSurface.cc \
	   FrontRealSurface.cc BackRealSurface.cc LeftRealSurface.cc RightRealSurface.cc \
	   VirtualSurface.cc ray.cc VolElement.cc MakeTableFunction.cc \
	

OBJS := $(patsubst %.cc,%.o,$(filter %.cc,$(SRCS)))

RMCRTRRSDStratifiedTheta : $(OBJS) 
			$(CXX) $(CFLAGS) $(OBJS) -o RMCRTRRSDStratifiedTheta

.cc.o: $<
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o *.d RMCRTRRSDStratifiedTheta *.out

-include $(SRCS:.cc=.d)
