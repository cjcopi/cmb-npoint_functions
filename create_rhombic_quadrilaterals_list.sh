#!/bin/sh

# $Id$

# This is the companion to the C++ code in
# create_rhombic_quadrilaterals_list.cpp
# that builds/runs that program and performs the necessary processing to
# get a compact quadrilateral table.

if [ $# -ne 2 ]; then
    echo "Usage: $0 <two point table filename> <output quad file>" 2>&1
    exit 1
fi

TPT_FILENAME=$1
OUTPUT_QUAD_FILE=$2
if [ ! -f ${TPT_FILENAME} ]; then
    echo "Two point table file does not exist : ${TPT_FILENAME}" 2>&1
    exit 1
fi

(make create_rhombic_quadrilaterals_list 2>&1 1> /dev/null) \
    || (echo Build failed; exit 1)

TMPFILE=data/tmp_rhombic.$$
# Create the sorted, unique table.
./create_rhombic_quadrilaterals_list.out ${TPT_FILENAME} \
    | perl -lane '@F=sort{$a<=>$b} @F; print "@F"' \
    | sort -n -k1 -k2 -k3 -k4 -u > ${TMPFILE}

# Now that we have the sorted table of unique quadrilaterals we want to
# compress it into a very terse, easy to use format. Use a python program
# for this.
./create_compressed_quadrilateral_table.py \
	${TPT_FILENAME} ${TMPFILE} ${OUTPUT_QUAD_FILE}
