# $Id: Makefile,v 1.5 2011-06-24 20:06:52 copi Exp $

# HEALPix.  Use the healpix-config I have written to make life easier.
HEALPIX_INC=`healpix-config --cppflags`
HEALPIX_LIBS=`healpix-config --cpplibs`

DOXYGEN = doxygen

USE_LIB_HEALPIX = create_twopt_table

override INCLUDES += -I.
# Set to the appropriate flag for openmp compilation, for
# g++ this is -fopenmp
OPENMP =

OPTIMIZE = -O3 -ffast-math -fomit-frame-pointer
CPPFLAGS = $(INCLUDES) $(OPTIMIZE) $(OPENMP)

# Sort also removes duplicates which is what we really want.
ALL_TARGETS=$(sort $(USE_LIB_HEALPIX))

# Rule for linking all the targets
$(ALL_TARGETS) :
	$(LINK.cc) -o $@.out $^ $(LIBS)

doc :
	$(DOXYGEN) Npoint_Functions.doxy

clean :
	$(RM) *.out *.o *~

clean-doc :
	$(RM) -r html

# Library dependencies.  Set the libraries and include paths.
$(USE_LIB_HEALPIX) : override LDFLAGS+=$(HEALPIX_LIBS)
$(USE_LIB_HEALPIX) : override CPPFLAGS+=$(HEALPIX_INC)

# Individual target dependencies
create_twopt_table : create_twopt_table.o

# Individual file dependencies
create_twopt_table.o : create_twopt_table.cpp \
	myRange.h buffered_pair_binary_file.h Twopt_Table.h
