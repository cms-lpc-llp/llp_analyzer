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
RUNNERSCC = $(addsuffix .cc,$(addprefix $(ANADIR),$(notdir RUNNERS)))
UTILS = $(SRCDIR)/JetCorrectorParameters.cc $(SRCDIR)/SimpleJetCorrectionUncertainty.cc  $(SRCDIR)/JetCorrectionUncertainty.cc $(SRCDIR)/SimpleJetCorrector.cc $(SRCDIR)/FactorizedJetCorrector.cc $(SRCDIR)/SimpleJetResolution.cc $(SRCDIR)/BTagCalibrationStandalone.cc $(SRCDIR)/EnergyScaleCorrection_class.cc $(SRCDIR)/EnergyScaleCorrection_class_2017.cc $(SRCDIR)/Hemisphere.cc $(SRCDIR)/Davismt2.cc $(SRCDIR)/RazorHelper.cc
UTILSOBJ = $(UTILS:cc=o)
EXECUTABLES = NormalizeNtuple SkimNtuple $(RUNNERS)
HELPERSCRIPT = python/MakeAnalyzerCode.py

.PHONY: clean all lxplus 

all: $(EXECUTABLES)
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
lxplus: all

clean:
	@-rm $(EXECUTABLES)
	@rm -f $(SRCDIR)/*.o $(ANADIR)/*.o
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

$(INCLUDEDIR)/rootdict.cc:
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST)

$(SRCDIR)/SimpleTable.o: $(SRCDIR)/SimpleTable.cc 
	$(CXX) -c $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorEvents.o: $(SRCDIR)/RazorEvents.C $(INCLUDEDIR)/RazorEvents.h
	$(CXX) $(SRCDIR)/RazorEvents.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorEventsRun1.o: $(SRCDIR)/RazorEventsRun1.C $(INCLUDEDIR)/RazorEventsRun1.h
	$(CXX) $(SRCDIR)/RazorEventsRun1.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorEventsUpgradeTiming.o: $(SRCDIR)/RazorEventsUpgradeTiming.C $(INCLUDEDIR)/RazorEventsUpgradeTiming.h
	$(CXX) $(SRCDIR)/RazorEventsUpgradeTiming.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)


$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzerRun1.o: $(SRCDIR)/RazorEventsRun1.o $(SRCDIR)/RazorAnalyzerRun1.cc
	$(CXX) $(SRCDIR)/RazorAnalyzerRun1.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(SRCDIR)/RazorAnalyzerUpgradeTiming.o: $(SRCDIR)/RazorEventsUpgradeTiming.o $(SRCDIR)/RazorAnalyzerUpgradeTiming.cc
	$(CXX) $(SRCDIR)/RazorAnalyzerUpgradeTiming.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

$(UTILSOBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSOBJ): $(ANADIR)/%.o: $(ANADIR)/%.cc $(ANADIR)/%.h
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) $<

$(ANALYZERSH): 
	$(HELPERSCRIPT) $(notdir $(basename $@))

$(RUNNERSCC): 
	$(HELPERSCRIPT) $(notdir $(basename $($@:Run=)))

$(RUNNERS): $(BINDIR)/Run%: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorEventsRun1.o $(SRCDIR)/RazorEventsUpgradeTiming.o $(SRCDIR)/RazorAnalyzer.o $(SRCDIR)/RazorAnalyzerRun1.o $(SRCDIR)/RazorAnalyzerUpgradeTiming.o $(UTILSOBJ) $(ANADIR)/%.o $(SRCDIR)/Run%.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) -I$(ANADIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS) 

NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

SkimNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/SkimNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)

MakePlots: $(SRCDIR)/SimpleTable.o ./macros/BackgroundStudies/OverlayKinematicPlots_Selected.C $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX14FLAGS)
