# Before this file in included, there must be a definition of
# BASEMATRICES_DIR
# Also, $HOME/.bashrc must have the definitions of the variables as shown in the INSTALL

CXX = g++

#BMLONG = -DLONG
#PROFILE = -pg
#DEBUG = -ggdb

CXXFLAGS = -O4 -fexceptions -Wall -Werror -Wextra $(DEBUG) $(BMLONG) $(PROFILE)

# Other commands
RM = rm -f
CP = cp -f
AR = ar rv

BASEMATRICESLIBNAME = libbasematrices.a

LIBDIR = $(HOME)/Libraries

BLAS = $(LIBDIR)/GotoBLAS2/libgoto2.a -lgfortran -lgfortranbegin -lpthread
METIS = $(LIBDIR)/metis-4.0.3/libmetis.a
CHOLMOD = $(LIBDIR)/CHOLMOD/Lib/libcholmod.a $(LIBDIR)/AMD/Lib/libamd.a \
  $(LIBDIR)/COLAMD/Lib/libcolamd.a $(LIBDIR)/CAMD/Lib/libcamd.a \
  $(LIBDIR)/CCOLAMD/Lib/libccolamd.a \
  $(LIBDIR)/SuiteSparse_config/libsuitesparseconfig.a

BASEMATRICESLIB = -L$(BASEMATRICES_DIR)/Lib -lbasematrices $(CHOLMOD) $(METIS) $(BLAS) -lrt

CHOLMODI = -I$(LIBDIR)/CHOLMOD/Include -I$(LIBDIR)/SuiteSparse_config
CUTERI = -I$(CUTER)/common/include
I = -I$(LIBDIR)/base_matrices/Include $(CHOLMODI) $(CUTERI)
