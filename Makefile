include Makefile.inc

DIRS = python
SRCDIR = src
INCLUDEDIR = include
ANADIR = analyzers
BINDIR = bin
INCLUDELIST= SimpleTable.h Linkdef.h

ANALYZERS = $(wildcard $(ANADIR)/*.cc)
ANALYZERSH = $(ANALYZERS:cc=h)
ANALYZERSOBJ = $(ANALYZERS:cc=o)
RUNNERS = $(addprefix $(BINDIR)/Run,$(notdir $(basename $(ANALYZERS))))
RUNNERSCC = $(addsuffix .cc,$(addprefix $(ANADIR)/,$(notdir $(RUNNERS))))
UTILS = $(SRCDIR)/JetCorrectorParameters.cc $(SRCDIR)/SimpleJetCorrectionUncertainty.cc  $(SRCDIR)/JetCorrectionUncertainty.cc $(SRCDIR)/SimpleJetCorrector.cc $(SRCDIR)/FactorizedJetCorrector.cc $(SRCDIR)/SimpleJetResolution.cc $(SRCDIR)/BTagCalibrationStandalone.cc $(SRCDIR)/EnergyScaleCorrection_class.cc $(SRCDIR)/EnergyScaleCorrection_class_2017.cc $(SRCDIR)/Hemisphere.cc $(SRCDIR)/Davismt2.cc $(SRCDIR)/RazorHelper.cc $(SRCDIR)/DBSCAN.cc $(SRCDIR)/LiteTreeMuonSystem.cc $(SRCDIR)/LiteLiteTreeMuonSystem.cc $(SRCDIR)/HNLMuonSystemTree.cc
UTILSOBJ = $(UTILS:cc=o)
EXECUTABLES = NormalizeNtuple SkimNtuple $(RUNNERS)
HELPERSCRIPT = python/MakeAnalyzerCode.py


.PHONY: clean all lxplus copy_runners

all: copy_runners $(EXECUTABLES)

lxplus: all

clean:
	@-rm $(EXECUTABLES)
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o

copy_runners:
		@for d in $(subst Run,,$(notdir $(basename $(RUNNERSCC)))); do ( if [ ! -f "src/Run"$$d".cc" ]; then echo $$d" file does not exists, copying"; $(HELPERSCRIPT) $$d; fi ) ; done



$(INCLUDEDIR)/rootdict.cc:
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST)

$(SRCDIR)/SimpleTable.o: $(SRCDIR)/SimpleTable.cc
	$(CXX) -c $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/llp_event.o: $(SRCDIR)/llp_event.C $(INCLUDEDIR)/llp_event.h
	$(CXX) $(SRCDIR)/llp_event.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(RUNNERS): $(BINDIR)/Run%: $(SRCDIR)/llp_event.o $(SRCDIR)/RazorAnalyzer.o $(UTILSOBJ) $(ANADIR)/%.o $(SRCDIR)/Run%.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

SkimNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/SkimNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

MakePlots: $(SRCDIR)/SimpleTable.o ./macros/BackgroundStudies/OverlayKinematicPlots_Selected.C $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
