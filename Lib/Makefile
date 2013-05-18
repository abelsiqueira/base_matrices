#
# Library
#

default: all

BASEMATRICES_DIR = ..
include $(BASEMATRICES_DIR)/make.inc

OBJ = base_common.o base_sparse.o base_factor.o base_dense.o base_triplet.o

all: $(BASEMATRICESLIBNAME) $(OBJ)

$(BASEMATRICESLIBNAME): $(OBJ)
	$(AR) $(BASEMATRICESLIBNAME) $(OBJ)

clean:
	$(RM) $(OBJ)

purge: clean
	$(RM) $(BASEMATRICESLIBNAME)

%.o: $(BASEMATRICES_DIR)/Src/%.cpp $(BASEMATRICES_DIR)/Include/%.h
	$(CXX) -c -o $@ $< $(I) $(CXXFLAGS) 

