include ../Makefile.inc

SHARED = libRazor.so
INCLUDELIST= RooCruijff.h Linkdef.h

all: $(SHARED) 

rootdict.cc: $(INCLUDELIST)
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) $^ 

$(SHARED): rootdict.o RooCruijff.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	mv $(SHARED) ../lib

clean:
	@rm -f *.o
	@rm -f rootdict.*
	@rm -f $(SHARED) 
	@rm -f ../lib/$(SHARED)
