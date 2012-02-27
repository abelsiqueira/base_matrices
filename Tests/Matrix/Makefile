#
# Matrices
#

default: finelem retgen finelemtriinf

BASEMATRICES_DIR = ../..
include $(BASEMATRICES_DIR)/make.inc

finelem: fin_elem
	./fin_elem 10

retgen: ret_gen
	./ret_gen 10

finelemtriinf: fin_elem_triinf
	./fin_elem_triinf 10

fin_elem: fin_elem.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

ret_gen: ret_gen.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

fin_elem_triinf: fin_elem_triinf.cpp
	$(CXX) -o $@ $< $(CXXFLAGS)

clean:
	$(RM) fin_elem ret_gen fin_elem_triinf *.mtx

