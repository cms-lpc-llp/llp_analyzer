CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

AUX = ../../Aux/

INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include -I$(AUX)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR)
CPPFLAGS += -g -std=c++11


TARGET  = CompareTrees
TARGET1 = PrintNtp1 
TARGET2 = PrintRazorEvents
TARGET3 = CompareEvents
TARGET4 = ListSelectedSusyHggTree

SRC = app/CompareTrees.cc src/ntp1.cc src/aux.cc
SRC1 = app/PrintNtp1.cc src/ntp1.cc
SRC2 = app/PrintRazorEvents.cc src/RazorEvents.cc
SRC3 = app/CompareEvents.cc src/aux.cc
SRC4 = app/ListSelectedSusyHggTree.cc src/SusyHggTree.cc src/HighRes.cc $(AUX)/src/CommandLineInput.cc

OBJ  = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)

all : $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET3) : $(OBJ3)
	$(LD) $(CPPFLAGS) -o $(TARGET3) $(OBJ3) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET4) : $(OBJ4)
	$(LD) $(CPPFLAGS) -o $(TARGET4) $(OBJ4) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) *~
