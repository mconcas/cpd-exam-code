#  USAGE:
#     make          ... to build the program

CC          = g++-4.9
CLINKER     = $(CC)
OPTFLAGS    = -fopenmp -DAPPLE -O3 -std=c++11
LIBS        = -lm
PRE         = ./
CFLAGS	  = $(OPTFLAGS)
OCL_LIBS 	= -framework OpenCL
OBJ=o
EXE=
RM=rm -f

EXES     =  vertexer$(EXE)
VER_OBJS	=  Event.$(OBJ) vertexer.$(OBJ) 

all: $(EXES)

vertexer$(EXE): $(VER_OBJS)
	$(CLINKER) $(CFLAGS) -o vertexer$(EXE) $(VER_OBJS) $(LIBS) $(OCL_LIBS)

clean:
	$(RM) $(EXES) *.$(OBJ)

.SUFFIXES:
.SUFFIXES: .c .cpp .$(OBJ)

.c.$(OBJ):
	$(CC) $(CFLAGS) -c $<

.cpp.$(OBJ):
	$(CC) $(CFLAGS) -c $<
