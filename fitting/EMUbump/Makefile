include Makefile.inc

DIRS = src

all:
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
clean:
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

force_look :
	true