include ../Makefile.inc

INCLUDEDIR= src
SHARED = lib/libRazorRun2.so
INCLUDELIST= RooRazor2DTail_SYS.h RooRazor3DBinPdf.h RooRazor3DBinReparamPdf.h RooRazor3DBinErrorPdf.h RooRazor3DBinNumericRsqTurnOffPdf.h RooRazor3DBinNumericMRSlopePdf.h RooRazor3DBinMRSlopePdf.h RooRazor3DBinMRSlopeErrorPdf.h MLMixedSample.h RooBTagMult.h RooSigmoid.h Linkdef.h

.PHONY: clean all

all: $(SHARED) 

$(INCLUDEDIR)/rootdict.cc: 
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST) 

$(INCLUDEDIR)/RooRazor2DTail_SYS.o:
	$(CXX) $(INCLUDEDIR)/RooRazor2DTail_SYS.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c
	
$(INCLUDEDIR)/RooRazor3DBinReparamPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinReparamPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinErrorPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinErrorPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinNumericRsqTurnOffPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinNumericRsqTurnOffPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinNumericMRSlopePdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinNumericMRSlopePdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinMRSlopePdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinMRSlopePdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinMRSlopeErrorPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinMRSlopeErrorPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooBTagMult.o:
	$(CXX) $(INCLUDEDIR)/RooBTagMult.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooSigmoid.o:
	$(CXX) $(INCLUDEDIR)/RooSigmoid.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/MLMixedSample.o:
	$(CXX) $(INCLUDEDIR)/MLMixedSample.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(SHARED): $(INCLUDEDIR)/rootdict.o $(INCLUDEDIR)/RooRazor2DTail_SYS.o $(INCLUDEDIR)/RooRazor3DBinPdf.o $(INCLUDEDIR)/RooRazor3DBinReparamPdf.o $(INCLUDEDIR)/RooRazor3DBinErrorPdf.o $(INCLUDEDIR)/RooRazor3DBinNumericRsqTurnOffPdf.o $(INCLUDEDIR)/RooRazor3DBinNumericMRSlopePdf.o $(INCLUDEDIR)/RooRazor3DBinMRSlopePdf.o $(INCLUDEDIR)/RooRazor3DBinMRSlopeErrorPdf.o $(INCLUDEDIR)/RooBTagMult.o $(INCLUDEDIR)/RooSigmoid.o $(INCLUDEDIR)/MLMixedSample.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^

clean:
	@rm -f $(INCLUDEDIR)/*.o
	@rm -f $(INCLUDEDIR)/rootdict.*
	@rm -f $(SHARED) 
