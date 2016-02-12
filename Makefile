# $Id: Makefile,v 1.28 2016/02/09 21:44:09 copi Exp $

# HEALPix.  Use the healpix-config I have written to make life easier.
HEALPIX_INC=`healpix-config --cppflags` -I$(HOME)/opt/myHealpix/include
HEALPIX_LIBS=`healpix-config --cpplibs`

DOXYGEN=doxygen

DEFINES=
override INCLUDES+=-I.
# Set to the appropriate flag for openmp compilation, for
# g++ this is -fopenmp.  See OPENMP_DEFAULT for targets that use
# this. Invoke make as
# make target OPENMP=
# to avoid using openmp for the default targets.  The -DOMP is included so
# that programs can use this for conditional inclusion of code.
OPENMP=-DOMP -fopenmp

OPTIMIZE=-O3 -ffast-math -fomit-frame-pointer -Wall -Wextra -Wno-unknown-pragmas

# Special handling of targets
USE_LIB_HEALPIX=create_twopt_table calculate_twopt_correlation_function \
	calculate_equilateral_threept_correlation_function \
	calculate_isosceles_threept_correlation_function \
	calculate_fourpt_correlation_function \
	calculate_LCDM_fourpt_correlation_function \
	calculate_constrained_fourpt_correlation_function \
	test_rhombic_quadrilaterals \
	create_rhombic_quadrilaterals_list \
	create_rhombic_quadrilaterals_list_parallel

# Targets that are built with openmp by default.  To turn this off for a
# compilation invoke make as
# make target OPENMP=
OPENMP_DEFAULT=create_twopt_table calculate_twopt_correlation_function \
	calculate_equilateral_threept_correlation_function \
	calculate_isosceles_threept_correlation_function \
	calculate_fourpt_correlation_function \
	calculate_LCDM_fourpt_correlation_function \
	calculate_constrained_fourpt_correlation_function \
	create_rhombic_quadrilaterals_list_parallel
# Targets that don't need anything special.
EXTRA_TARGETS=

# Sort also removes duplicates which is what we really want.
ALL_TARGETS=$(sort $(USE_LIB_HEALPIX) \
                   $(OPENMP_DEFAULT) $(EXTRA_TARGETS) )

CPPFLAGS=$(INCLUDES) $(OPTIMIZE) $(DEFINES)

all :
	@echo
	@echo Available targets:
# Print the targets one per line
	@echo $(ALL_TARGETS) | sed 's/\s/\n/g' | sed 's/^/  /g'
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
$(USE_LIB_HEALPIX) : override LIBS+=$(HEALPIX_LIBS)
$(USE_LIB_HEALPIX) : override CPPFLAGS+=$(HEALPIX_INC)
$(USE_LIB_LZMA) : override LIBS+=-llzma
$(USE_LIB_Z) : override LIBS+=-lz
$(OPENMP_DEFAULT) : override CPPFLAGS+=$(OPENMP) 

# Individual target dependencies
create_twopt_table : create_twopt_table.o
calculate_twopt_correlation_function : calculate_twopt_correlation_function.o
calculate_equilateral_threept_correlation_function : \
	calculate_equilateral_threept_correlation_function.o
calculate_isosceles_threept_correlation_function : \
	calculate_isosceles_threept_correlation_function.o
calculate_fourpt_correlation_function : \
	calculate_fourpt_correlation_function.o
calculate_LCDM_fourpt_correlation_function : \
	calculate_LCDM_fourpt_correlation_function.o
calculate_constrained_fourpt_correlation_function : \
	calculate_constrained_fourpt_correlation_function.o
test_rhombic_quadrilaterals : \
	test_rhombic_quadrilaterals.o
create_rhombic_quadrilaterals_list : \
	create_rhombic_quadrilaterals_list.o
create_rhombic_quadrilaterals_list_parallel : \
	create_rhombic_quadrilaterals_list_parallel.o

# Individual file dependencies
create_twopt_table.o : create_twopt_table.cpp \
	buffered_pair_binary_file.h Twopt_Table.h
calculate_twopt_correlation_function.o : \
	calculate_twopt_correlation_function.cpp \
	Twopt_Table.h
calculate_equilateral_threept_correlation_function.o : \
	calculate_equilateral_threept_correlation_function.cpp \
	Twopt_Table.h Pixel_Triangles.h \
	Npoint_Functions_Utils.h
calculate_isosceles_threept_correlation_function.o : \
	calculate_isosceles_threept_correlation_function.cpp \
	Twopt_Table.h Pixel_Triangles.h \
	Npoint_Functions_Utils.h
calculate_fourpt_correlation_function.o : \
	calculate_fourpt_correlation_function.cpp \
	Quadrilateral_List_File.h \
	Npoint_Functions_Utils.h
calculate_LCDM_fourpt_correlation_function.o : \
	calculate_LCDM_fourpt_correlation_function.cpp \
	Quadrilateral_List_File.h \
	Npoint_Functions_Utils.h
calculate_constrained_fourpt_correlation_function.o : \
	calculate_constrained_fourpt_correlation_function.cpp \
	Quadrilateral_List_File.h \
	Npoint_Functions_Utils.h
test_rhombic_quadrilaterals.o : \
	test_rhombic_quadrilaterals.cpp \
	Twopt_Table.h Pixel_Triangles.h Pixel_Quadrilaterals.h \
	Npoint_Functions_Utils.h
create_rhombic_quadrilaterals_list.o : \
	create_rhombic_quadrilaterals_list.cpp \
	Twopt_Table.h Pixel_Triangles.h Pixel_Quadrilaterals.h
create_rhombic_quadrilaterals_list_parallel.o : \
	create_rhombic_quadrilaterals_list_parallel.cpp \
	Twopt_Table.h Pixel_Triangles.h Pixel_Quadrilaterals.h
