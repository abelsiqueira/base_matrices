# Before this file in included, there must be a definition of
# BASEMATRICES_DIR
# Also, $HOME/.bashrc must have the definitions of the variables as shown in the INSTALL

CXX = g++

#BMLONG = -DLONG
PROFILE = -pg
DEBUG = -ggdb

CXXFLAGS = -O4 -fexceptions -Wall -Werror -Wextra $(DEBUG) $(BMLONG) $(PROFILE)

# Other commands
RM = rm -f
CP = cp -f
AR = ar rv

BASEMATRICESLIBNAME = libbasematrices.a

BASEMATRICESLIB = $(BASEMATRICES_DIR)/Lib/$(BASEMATRICESLIBNAME) -lcholmod \
     -lamd -lcolamd -lccolamd -lcamd -lsuitesparseconfig -lmetis -lgoto2 \
     -lgfortran -lgfortranbegin -lpthread -lrt
