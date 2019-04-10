CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lMathMore

CPPFLAGS += -g

TARGET = PlotOnOutputs
TARGET1 = GetTableOfSignalEfficiency
TARGET2 = 2016vs2017
TARGET3 = MakeRoot
TARGET4 = CalculateYield
TARGET5 = RazorVsMt2
TARGET6 = ReducedTree
TARGET7 = EcalPrefiring
TARGET8 = MetRecipe

SRC = app/PlotOnOutputs.cc
SRC1 = app/GetTableOfSignalEfficiency.cc
SRC2 = app/2016vs2017.cc
SRC3 = app/MakeRoot.cc
SRC4 = app/CalculateYield.cc  
SRC5 = app/RazorVsMt2.cc  
SRC6 = app/ReducedTree.cc  
SRC7 = app/EcalPrefiring.cc
SRC8 = app/MetRecipe.cc

OBJ = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ3 = $(SRC3:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)
OBJ5 = $(SRC5:.cc=.o)
OBJ6 = $(SRC6:.cc=.o)
OBJ7 = $(SRC7:.cc=.o)
OBJ8 = $(SRC8:.cc=.o)

all : $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8)


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

$(TARGET5) : $(OBJ5)
	$(LD) $(CPPFLAGS) -o $(TARGET5) $(OBJ5) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET6) : $(OBJ6)
	$(LD) $(CPPFLAGS) -o $(TARGET6) $(OBJ6) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET7) : $(OBJ7)
	$(LD) $(CPPFLAGS) -o $(TARGET7) $(OBJ7) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET8) : $(OBJ8)
	$(LD) $(CPPFLAGS) -o $(TARGET8) $(OBJ8) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc	
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@	
	@echo $<
clean :
	rm -f *.o src/*.o $(TARGET) $(TARGET1)  $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8) *~

