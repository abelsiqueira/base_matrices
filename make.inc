# Before this file in included, there must be a definition of
# BASEMATRICES_DIR
# Also, $HOME/.bashrc must have the definitions of the variables as shown in the INSTALL

include $(BASEMATRICES_DIR)/config.mk

# Compilation
CXX = g++
CXXFLAGS = -O4 -fexceptions -Wall -Werror -Wextra -ggdb $(BMLONG) -m32

# Other commands
RM = rm -f
CP = cp -f
AR = ar rv

BASEMATRICESLIBNAME = libbasematrices.a

BASEMATRICESLIB = $(BASEMATRICES_DIR)/Lib/$(BASEMATRICESLIBNAME) -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -lgoto2 -lgfortran -lgfortranbegin
