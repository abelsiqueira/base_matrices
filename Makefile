default: error

BASEMATRICES_DIR = ./
include make.inc

error:
	@echo '------------------------------------------'
	@echo '|             base_matrices              |'
	@echo '------------------------------------------'
	@echo '| make all       | make library and test |'
	@echo '| make library   | make library          |'
	@echo '| make test      | test                  |'
	@echo '| make clean     | clean                 |'
	@echo '| make purge     | purge                 |'
	@echo '------------------------------------------'

all: library test

library:
	(cd Lib; make )

test:
	(cd Tests; make )

clean:
	(cd Lib; make clean )
	(cd Tests; make clean )

purge:
	(cd Lib; make purge )
	(cd Tests; make purge )

##### External Libraries
chollibrary:
	( cd $(METISDIR)/Lib ; $(MAKE) )
	( cd $(GOTOBLASDIR) ; $(MAKE) libs netlib shared )
	( cd $(TIMDAVISDIR)/CHOLMOD ; $(MAKE) library )

choltests:
	( cd $(METISDIR) ; $(MAKE) )
	( cd $(GOTOBLASDIR) ; $(MAKE) tests )

cholclean:
	( cd $(METISDIR) ; $(MAKE) clean )
	( cd $(GOTOBLASDIR) ; $(MAKE) clean )
	( cd $(TIMDAVISDIR)/CHOLMOD ; $(MAKE) clean )

cholpurge:
	( cd $(METISDIR) ; $(MAKE) realclean )
	( cd $(GOTOBLASDIR) ; $(MAKE) clean )
	( cd $(TIMDAVISDIR) ; $(MAKE) purge )
