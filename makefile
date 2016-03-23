#  USAGE:
#     make          ... to build the program

CC       = g++-4.9
CLINKER  = $(CC)
OPTFLAGS = -Wall -fopenmp -DAPPLE -O3 -std=c++11
CFLAGS	 = $(OPTFLAGS)
OCL_LIBS = -framework OpenCL
BASEDIR  = base
INCDIR   = ./$(BASEDIR)
OBJ      = o
SRC      = cpp
EXE      =
RM       = rm -f
EXES     = vertexer$(EXE)
VER_OBJS = Event.$(OBJ) VertexerFast.$(OBJ) VertexCandidate.$(OBJ) 

all: $(EXES)
	$(RM) $(VER_OBJS)

vertexer$(EXE): $(VER_OBJS)
	$(CLINKER) $(CFLAGS) -o vertexer$(EXE) main.cc -I$(INCDIR) $(VER_OBJS) $(OCL_LIBS)

Event.$(OBJ):
	$(CC) $(CFLAGS) -c $(BASEDIR)/Event.cpp

VertexCandidate.$(OBJ):
	$(CC) $(CFLAGS) -c $(BASEDIR)/VertexCandidate.cpp

VertexerFast.$(OBJ):	VertexCandidate.$(OBJ)
	$(CC) $(CFLAGS) -c $(BASEDIR)/VertexerFast.cpp

debug: 	CFLAGS += -DDEBUG
debug:	all

clean:
	$(RM) $(EXES) *.$(OBJ)
