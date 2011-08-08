#!/bin/sh

# Simple script to generate the rhombic quadrilaterals for a particular
# case and test them against a known good result.

# $Id$

# Known good results
DATAFILE=data/test_rhombic_quadrilaterals.dat

(make test_rhombic_quadrilaterals 2>&1 1> /dev/null) \
    || (echo Build failed; exit 1)
TMPFILE=test_tmp.$$
./test_rhombic_quadrilaterals.out | \
    perl -lane '@F=sort{$a<=>$b}@F; print "@F"' > ${TMPFILE}
RES=`sort -n ${TMPFILE} ${DATAFILE} | sort -n | uniq -c | egrep '^\s*1'`
if [ "X$RES" != "X" ]; then
    echo "Test failed:"
    echo "$TMPFILE contains the data with errors"
else
    /bin/rm $TMPFILE
fi
