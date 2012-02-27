# Before this file in included, there must be a definition of
# BASEMATRICES_DIR

include $(BASEMATRICES_DIR)/config.mk

# Compilation
CXX = g++
CXXFLAGS = -O4 -fexceptions -Wall -Werror -Wextra -ggdb $(BMLONG) -m32

# Other commands
RM = rm -f
CP = cp -f
AR = ar rv

##############################################
###### Adjust this variables accordingly #####
##############################################

##### METIS #####
METIS_PATH = $(HOME)/Libraries/metis
METIS = $(METIS_PATH)/libmetis.a

##### GotoBLAS2 #####
GOTODIR = $(HOME)/Libraries/GotoBLAS2
GOTO = $(GOTODIR)/libgoto2.a $(GOTODIR)/libgoto2.so -lgfortran -lgfortranbegin

##### CHOLMOD, AMD, CAMD, COLAMD, CCOLAMD, UFconfig #####
CHOLMOD_ROOT = $(HOME)/Libraries
AMD = $(CHOLMOD_ROOT)/AMD
CAMD = $(CHOLMOD_ROOT)/CAMD
COLAMD = $(CHOLMOD_ROOT)/COLAMD
CHOLMOD = $(CHOLMOD_ROOT)/CHOLMOD
CCOLAMD = $(CHOLMOD_ROOT)/CCOLAMD
UFCONFIG = $(CHOLMOD_ROOT)/UFconfig
CHOLMODI = -I$(CHOLMOD)/Include -I$(UFCONFIG)
CHOLMODLIB = $(CHOLMOD)/Lib/libcholmod.a $(AMD)/Lib/libamd.a $(CAMD)/Lib/libcamd.a $(COLAMD)/Lib/libcolamd.a $(CCOLAMD)/Lib/libccolamd.a 

###########################################
##### STOP editing anything past here #####
###########################################
BASEMATRICESLIBNAME = libbasematrices.a
INCLUDES = -I$(BASEMATRICES_DIR)/Include -I$(METIS_PATH)/Lib $(CHOLMODI)
LIBS = $(BASEMATRICES_DIR)/Lib/$(BASEMATRICESLIBNAME) $(CHOLMODLIB) $(METIS) $(GOTO) 
BASEMATRICESLIB = $(LIBS)
BASEMATRICESI = $(INCLUDES)

#BASEMATRICESLIB = $(BASEMATRICES_DIR)/Lib/$(BASEMATRICESLIBNAME) -lcholmod -lamd -lcolamd -lccolamd -lcamd -lmetis -lgoto2 -lgfortran -lgfortranbegin
