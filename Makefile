INCDIR= $(shell pwd)/inc
SRCDIR= $(shell pwd)/src
OBJDIR= $(shell pwd)/obj
BINDIR= $(shell pwd)/bin

VPATH = $(SRCDIR)

CXX=g++
#CFLAGS=-c -O3 -Xclang -fopenmp -g -Wall `root-config --cflags` -I${INCDIR}
CFLAGS=-c -O3 -g -Wall -DDEBUG_PRINT_LEVEL=4 -DROOT_OUTPUT_LEVEL=4 `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter -lMinuit2 -lGenVector

TARGET=run_ssd_reconstruction.cpp

EXECUTABLE=$(TARGET:%.cpp=$(BINDIR)/%)

FILES= $(wildcard $(SRCDIR)/*.cpp)
SOURCES=$(FILES)

OBJECTS = $(FILES:$(SRCDIR)/%.cpp=${OBJDIR}/%.o)
INCLUDES = $(wildcard $(INCDIR)/*.h)

OBJ=$(TARGET:%.cpp=${OBJDIR}/%.o) $(OBJECTS)

all: $(TARGET) $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)
	
$(OBJDIR)/%.o: %.cpp $(INCLUDES)
	$(CXX) $(CFLAGS) $< -o $@

print-%  : ; @echo $* = $($*)

clean:
	- $(RM) $(BINDIR)/* $(OBJDIR)/*
