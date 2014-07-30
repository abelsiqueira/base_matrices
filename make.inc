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

LIBNAME = libbasematrices.a

BLAS = -lopenblas -lgfortran -lgfortranbegin -lpthread
METIS = -lmetis
CHOLMOD = -lcholmod -lamd -lcolamd -lcamd -lccolamd -lsuitesparseconfig

LIBS = -L$(BASEMATRICES_DIR)/Lib -lbasematrices $(CHOLMOD) $(METIS) $(BLAS) -lrt
I = -I$(BASEMATRICES_DIR)/Include/
