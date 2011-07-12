# $Id: Makefile,v 1.7 2011-07-10 03:07:06 copi Exp $

# HEALPix.  Use the healpix-config I have written to make life easier.
HEALPIX_INC=`healpix-config --cppflags`
HEALPIX_LIBS=`healpix-config --cpplibs`

DOXYGEN = doxygen

DEFINES=
override INCLUDES += -I.
# Set to the appropriate flag for openmp compilation, for
# g++ this is -fopenmp.  See OPENMP_DEFAULT for targets that use this.
OPENMP=-fopenmp

OPTIMIZE = -O3 -ffast-math -fomit-frame-pointer -Wall -W

# Special handling of targets
USE_LIB_HEALPIX = create_twopt_table calculate_twopt_correlation_function \
	calculate_equilateral_threept_correlation_function
ifdef USE_LZMA_COMPRESSION
	override DEFINES+=-DUSE_LZMA_COMPRESSION
	USE_LIB_LZMA = create_twopt_table \
		calculate_twopt_correlation_function \
		calculate_equilateral_threept_correlation_function
	USE_LIB_Z=
else
	USE_LIB_Z = create_twopt_table calculate_twopt_correlation_function \
		calculate_equilateral_threept_correlation_function
	USE_LIB_LZMA=
endif
# Targets that are built with openmp by default.  To turn this off for a
# compilation invoke make as
# make target OPENMP=
OPENMP_DEFAULT = create_twopt_table calculate_twopt_correlation_function \
	calculate_equilateral_threept_correlation_function
# Targets that don't need anything special.
EXTRA_TARGETS =

# Sort also removes duplicates which is what we really want.
ALL_TARGETS=$(sort $(USE_LIB_HEALPIX) $(USE_LIB_LZMA) \
                   $(OPENMP_DEFAULT) $(EXTRA_TARGETS) )

CPPFLAGS = $(INCLUDES) $(OPTIMIZE) $(DEFINES)

all :
	@echo Available targets: $(ALL_TARGETS)
	@echo Use a command like: make target USE_LZMA_COMPRESSION=1
	@echo to use LZMA compression instead of libz.
	@echo "  [Note that libz is about 5 times faster in creating two"
	@echo "   point tables and slightly faster in calculating the two point"
	@echo "  correlation function.]"

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
$(USE_LIB_LZMA) : override LDFLAGS+=-llzma
$(USE_LIB_Z) : override LDFLAGS+=-lz
$(OPENMP_DEFAULT) : override CPPFLAGS+=$(OPENMP)

# Individual target dependencies
create_twopt_table : create_twopt_table.o
calculate_twopt_correlation_function : calculate_twopt_correlation_function.o
calculate_equilateral_threept_correlation_function : \
	calculate_equilateral_threept_correlation_function.o

# Individual file dependencies
create_twopt_table.o : create_twopt_table.cpp \
	myRange.h buffered_pair_binary_file.h Twopt_Table.h \
	ZLIB_Wrapper.h LZMA_Wrapper.h
calculate_twopt_correlation_function.o : \
	calculate_twopt_correlation_function.cpp \
	Twopt_Table.h \
	ZLIB_Wrapper.h LZMA_Wrapper.h
calculate_equilateral_threept_correlation_function.o : \
	calculate_equilateral_threept_correlation_function.cpp \
	Twopt_Table.h Pixel_Triangles.h \
	ZLIB_Wrapper.h LZMA_Wrapper.h
