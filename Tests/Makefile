#
# Tests
#

default: all

BASEMATRICES_DIR = ..
include $(BASEMATRICES_DIR)/make.inc

LIBSHORT = $(BASEMATRICES_DIR)/Lib/$(LIBNAME)
EXEC = test1.tst test2.tst test3.tst test4.tst test_dense.tst
#VAL = valgrind
#VALFLAGS = --leak-check=full -v
VAL =
VALFLAGS =

all: $(EXEC)
	( cd Matrix; make )
	$(VAL) $(VALFLAGS) ./test1.tst < Matrix/fin_elem_mm10.mtx
	$(VAL) $(VALFLAGS) ./test1.tst < Matrix/fin_elem_triinf_mm10.mtx
	$(VAL) $(VALFLAGS) ./test2.tst < Matrix/test.mtxnvec
	$(VAL) $(VALFLAGS) ./test3.tst < Matrix/ret_mm10.mtx
	$(VAL) $(VALFLAGS) ./test4.tst < Matrix/test.mtxnvec
	$(VAL) $(VALFLAGS) ./test_dense.tst

ldl.tst: ldl.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(LIBS) $(I) $(CXXFLAGS)

%.tst: %.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(LIBS) $(I) $(CXXFLAGS)

clean:
	( cd Matrix; make clean )
	$(RM) $(EXEC) ldl

purge: clean
