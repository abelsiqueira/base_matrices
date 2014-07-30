#
# Library
#

default: all

BASEMATRICES_DIR = ..
include $(BASEMATRICES_DIR)/make.inc

OBJ = base_common.o base_sparse.o base_factor.o base_dense.o base_triplet.o

all: $(LIBNAME) $(OBJ)

$(LIBNAME): $(OBJ)
	$(AR) $(LIBNAME) $(OBJ)

clean:
	$(RM) $(OBJ)

purge: clean
	$(RM) $(LIBNAME)

%.o: $(BASEMATRICES_DIR)/Src/%.cpp $(BASEMATRICES_DIR)/Include/%.h
	$(CXX) -c -o $@ $< $(I) $(CXXFLAGS)

install: $(LIBNAME)
	$(CP) $(LIBNAME) /usr/lib/
	$(CP) ../Include/*.h /usr/include/
