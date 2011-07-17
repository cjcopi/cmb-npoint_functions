# $Id: Makefile,v 1.13 2011-07-17 03:58:25 copi Exp $

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
	calculate_equilateral_threept_correlation_function \
	calculate_isosceles_threept_correlation_function
# Targets that may use compression
USE_COMPRESSION = create_twopt_table \
		calculate_twopt_correlation_function \
		calculate_equilateral_threept_correlation_function \
		calculate_isosceles_threept_correlation_function
ifdef USE_NO_COMPRESSION
	override DEFINES+=-DUSE_NO_COMPRESSION
	COMPRESSION_WRAPPER=No_Compression_Wrapper.h
else ifdef USE_LZMA_COMPRESSION
	override DEFINES+=-DUSE_LZMA_COMPRESSION
	USE_LIB_LZMA = $(USE_COMPRESSION)
	USE_LIB_Z=
	COMPRESSION_WRAPPER=LZMA_Wrapper.h
else
	USE_LIB_Z = $(USE_COMPRESSION)
	USE_LIB_LZMA=
	COMPRESSION_WRAPPER=ZLIB_Wrapper.h
endif
# Targets that are built with openmp by default.  To turn this off for a
# compilation invoke make as
# make target OPENMP=
OPENMP_DEFAULT = create_twopt_table calculate_twopt_correlation_function \
	calculate_equilateral_threept_correlation_function \
	calculate_isosceles_threept_correlation_function
# Targets that don't need anything special.
EXTRA_TARGETS =

# Sort also removes duplicates which is what we really want.
ALL_TARGETS=$(sort $(USE_LIB_HEALPIX) $(USE_COMPRESSION) \
                   $(OPENMP_DEFAULT) $(EXTRA_TARGETS) )

CPPFLAGS = $(INCLUDES) $(OPTIMIZE) $(DEFINES)

all :
	@echo
	@echo Available targets:
# Print the targets one per line
	@echo $(ALL_TARGETS) | sed 's/\s/\n/g' | sed 's/^/  /g'
	@echo Use a command like: make target USE_LZMA_COMPRESSION=1
	@echo to use LZMA compression instead of libz.
	@echo Use a command like: make target USE_NO_COMPRESSION=1
	@echo to use no compression.
	@echo "  [Note that libz is about 5 times faster in creating two"
	@echo "   point tables and slightly faster in calculating the two point"
	@echo "  correlation function.]"
	@echo

# Rule for linking all the targets
$(ALL_TARGETS) :
	$(LINK.cc) -o $@.out $^ $(LIBS)

doc :
	$(DOXYGEN) Npoint_Functions.doxy

clean :
	$(RM) *.out *.o *~

clean-doc :
	$(RM) -r html

clean-all : clean clean-doc

# Targets to always be run when asked for.
.PHONY : all clean clean-doc clean-all doc

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
calculate_isosceles_threept_correlation_function : \
	calculate_isosceles_threept_correlation_function.o
# Individual file dependencies
create_twopt_table.o : create_twopt_table.cpp \
	buffered_pair_binary_file.h Twopt_Table.h \
	$(COMPRESSION_WRAPPER)
calculate_twopt_correlation_function.o : \
	calculate_twopt_correlation_function.cpp \
	Twopt_Table.h \
	$(COMPRESSION_WRAPPER)
calculate_equilateral_threept_correlation_function.o : \
	calculate_equilateral_threept_correlation_function.cpp \
	Twopt_Table.h Pixel_Triangles.h \
	$(COMPRESSION_WRAPPER) \
	Npoint_Functions_Utils.h
calculate_isosceles_threept_correlation_function.o : \
	calculate_isosceles_threept_correlation_function.cpp \
	Twopt_Table.h Pixel_Triangles.h \
	$(COMPRESSION_WRAPPER) \
	Npoint_Functions_Utils.h
