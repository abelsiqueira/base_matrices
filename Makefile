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

install:
	$(MAKE) -C Lib install
