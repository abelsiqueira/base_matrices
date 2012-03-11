#
# Matrices
#

default: finelem retgen finelemtriinf

BASEMATRICES_DIR = ../..
include $(BASEMATRICES_DIR)/make.inc

finelem: fin_elem.tst
	./fin_elem.tst 10

retgen: ret_gen.tst
	./ret_gen.tst 10

finelemtriinf: fin_elem_triinf.tst
	./fin_elem_triinf.tst 10

fin_elem.tst: fin_elem.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

ret_gen.tst: ret_gen.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

fin_elem_triinf.tst: fin_elem_triinf.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

clean:
	$(RM) fin_elem.tst ret_gen.tst fin_elem_triinf.tst *.mtx

