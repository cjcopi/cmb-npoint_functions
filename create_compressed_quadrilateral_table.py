#!/usr/bin/env python

import sys, re
import struct
import numpy as np
import Twopt_Table

def usage (progname) :
    sys.stderr.write(
        """Usage: %s <arguments>
Arguments:
  <two point table filename>
  <input quadrilaterals list filename>
  <output quadrilaterals filename>
"""%(progname))
    sys.exit(0)


class quadrilateral_list_file (object) :
    def __init__ (self, filename, twopt_table) :
        self.maxbytes = -1
        self.fp = open (filename, "wb")
        self.tpt = twopt_table
        self.write_header()

    def __del__ (self) :
        self.fp.close()

    def write_maxbytes (self) :
        currpos = self.fp.tell()
        self.fp.seek(self.maxbytes_loc, 0)
        self.fp.write(struct.pack ("L", self.maxbytes))
        self.fp.seek(currpos, 0)

    def write_header (self) :
        version = '\x01'
        if self.tpt.isNest() :
            scheme = '\x00'
        else :
            scheme = '\x01'
        self.fp.write(version)
        self.fp.write(struct.pack ("L", self.tpt.Nside()))
        self.fp.write(scheme)
        self.fp.write(struct.pack ("d", self.tpt.bin_value()))
        # Padding for max bytes
        self.fp.write(struct.pack ("8x"))
        self.maxbytes_loc = 1+8+1+8

    def add_entry (self, res) :
        n = len(res)
        bytes = 4*n
        if bytes > self.maxbytes :
            self.maxbytes = bytes
            self.write_maxbytes()
        self.fp.write(struct.pack ("L", bytes))
        self.fp.write(struct.pack ("%di"%n, *res))

class compress_quad_list (object) :
    def __init__ (self) :
        self.prev = - np.ones(4)
        self.res = [ np.array([], dtype=np.int),
                     np.array([], dtype=np.int),
                     np.array([], dtype=np.int),
                     np.array([], dtype=np.int)
                     ]
        self.count = np.zeros(4, dtype=np.int)

    def _split_line (self, line) :
        s = re.split ('\s+', line.strip())
        if s is None or len(s) != 4 : return None
        return map (lambda x : int(x), s)

    def initialize (self, line) :
        pts = self._split_line (line)
        if pts is None : return False
        self.prev = np.asarray(pts, dtype=np.int)
        for j in range(len(self.res)) :
            self.res[j] = np.array([], dtype=np.int)
            self.count[j] = 1
        # Yes, this one is treated special
        self.res[3] = np.array([pts[3]], dtype=np.int)
	return True

    def next_line (self, line) :
        pts = self._split_line (line)
        if pts[2] == self.prev[2] \
                and pts[1] == self.prev[1] \
                and pts[0] == self.prev[0] :
            self.res[3] = np.concatenate((self.res[3], [pts[3]]))
            self.count[3] += 1
            return None

        # Otherwise we start concatenating lists.  If we get here we KNOW
        # that we need to process pts[2]
        self.res[2] = np.concatenate((self.res[2], [self.prev[2]],
                                      [self.count[3]], self.res[3]))
        self.res[3] = np.array([pts[3]], dtype=np.int)
        self.count[3] = 1
        self.count[2] += 1
        self.prev[2] = pts[2]

        if pts[1] != self.prev[1] or pts[0] != self.prev[0] :
            self.res[1] = np.concatenate((self.res[1], [self.prev[1]],
                                          [self.count[2]-1], self.res[2]))
            self.res[2] = np.array([], dtype=np.int)
            self.count[2] = 1
            self.count[1] += 1
            self.prev[1] = pts[1]

        if pts[0] != self.prev[0] :
            # Concatentate and return
            ret = np.concatenate((self.res[0], [self.prev[0]],
                                  [self.count[1]-1], self.res[1]))
            self.res[1] = np.array([], dtype=np.int)
            self.count[1] = 1
            self.res[0] = np.array([], dtype=np.int)
            self.prev[0] = pts[0]
            return ret
        return None

    def finalize (self) :
        return self.next_line ("-1 -1 -1 -1")


if len(sys.argv) != 4 :
    usage (sys.argv[0])

tpt = Twopt_Table.Twopt_Table();
tpt.read_file (sys.argv[1])

fp_quad_in = open (sys.argv[2], "r")
qlf = quadrilateral_list_file (sys.argv[3], tpt)

cql = compress_quad_list()
if cql.initialize (fp_quad_in.readline()) :
    for line in fp_quad_in :
        res = cql.next_line (line)
        if res is not None :
            qlf.add_entry (res)
    res = cql.finalize()
    qlf.add_entry (res)

fp_quad_in.close()



