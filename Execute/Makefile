CXX	= g++
RM	= rm -f
INCLUDES = -I/usr/local/include
LDFLAGS	=
LDLIBS	= -L/usr/local/lib -lgsl -lgslcblas -lm

SRCS	= ../Modules/main.cpp
OBJS	= $(SRCS:.cpp=.o)


DUST_SCAT: $(SRCS)
	$(CXX) $(INCLUDES) -o MAIN $(SRCS) $(LDLIBS)

clean:
	$(RM) $(OBJS)

clobber:
	$(RM) MAIN ../Results/*.dat ../Results/*/*.dat
