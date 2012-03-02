#
# Tests
#

default: all

BASEMATRICES_DIR = ..
include $(BASEMATRICES_DIR)/make.inc

LIBSHORT = $(BASEMATRICES_DIR)/Lib/$(BASEMATRICESLIBNAME)
EXEC = test1 test2 test3 test4 test_dense
#VAL = valgrind
#VALFLAGS = --leak-check=full -v
VAL = 
VALFLAGS =

all: $(EXEC)
	( cd Matrix; make )
	$(VAL) $(VALFLAGS) ./test1 < Matrix/fin_elem_mm10.mtx
	$(VAL) $(VALFLAGS) ./test1 < Matrix/fin_elem_triinf_mm10.mtx
	$(VAL) $(VALFLAGS) ./test2 < Matrix/test.mtxnvec
	$(VAL) $(VALFLAGS) ./test3 < Matrix/ret_mm10.mtx
	$(VAL) $(VALFLAGS) ./test4 < Matrix/test.mtxnvec
	$(VAL) $(VALFLAGS) ./test_dense

ldl: ldl.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS)

test1: test1.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS) 

test2: test2.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS) 
	
test3: test3.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS) 

test4: test4.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS) 

test_dense: test_dense.cpp $(LIBSHORT)
	$(CXX) -o $@ $< $(BASEMATRICESLIB) $(BASEMATRICESI) $(CXXFLAGS)

clean:
	( cd Matrix; make clean )
	$(RM) $(EXEC) ldl

purge: clean
