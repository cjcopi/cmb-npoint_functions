import sys
import struct
import zlib
import numpy as np
import healpy

# $Id: Twopt_Table.py,v 1.2 2011-07-30 02:29:08 copi Exp $

class Twopt_Table :
    """Python class for the storage of a single bin of a two point table.  See
    the C++ version for details.  Only zlib compression and integer tables
    have been implemented."""

    def __init__ (self) :
        self.table = []
        self.pixlist = []
        self.cosbin = 0
        self.nside = 0
        self.nmax = 0
        self.scheme = '\x00'

    def _read_table_from_stream (self, fd) :
        """Read the table from the stream with compression.
        The Nmax() and Npix() MUST be set correctly before calling.
        Here \a fd is the file stream.  The read table is set after a
        successful call to this this method.
        
        This is used internally by read_file().
        It should be used with caution.
        """
        Nelem = self.Nmax()*self.Npix()
        buf = fd.read()
        buf = zlib.decompress(buf)
        self.table = struct.unpack ("%di"%Nelem, buf)
        
    def _read_header_from_stream (self, fd) :
        """Read the header from the stream.
        It is assumed the header starts at the current stream position.
        On success thee stream position is left immediately after the
        header. On failure the stream is left in an undefined state.

        This is used internally by read_file().
        It should be used with caution.
        """
        version = fd.read (1)
        if version != '\x03' :
            sys.stderr.write ("Twopt_Table only supports file format version 3\n")
            return False

        self.cosbin = struct.unpack ("d", fd.read(8))[0]
        self.nside = struct.unpack ("L", fd.read(8))[0]
        npix = struct.unpack ("L", fd.read(8))[0]
        self.pixlist = struct.unpack ("%di"%(npix),
                                      fd.read(4*npix)) 
        self.scheme = fd.read (1)
        self.nmax = struct.unpack ("L", fd.read(8))[0]


    def _write_table_to_stream (self, fd) :
        """Write the table to the stream with compression.
        Here \a fd is the file stream.
        
        This is used internally by write_file().
        It should be used with caution.
        """
        Nelem = len(self.table)
        buf = struct.pack ("%di"%Nelem, *(self.table))
        bufc = zlib.compress(buf, 6)
        print len(buf), len(bufc)
        fd.write (bufc)

    def _write_header_to_stream (self, fd) :
        """Write the header from the stream.
        It is assumed the header is to start at the current stream position.
        On success the stream position is left immediately after the
        header. On failure the stream is left in an undefined state.

        This is used internally by write_file().
        It should be used with caution.
        """
        version = '\x03' # Write version 3
        fd.write (version)
        fd.write (struct.pack ("dLL",
                               self.cosbin, self.nside, self.Npix()))
        fd.write (struct.pack ("%di"%(self.Npix()), *(self.pixlist)))
        fd.write (self.scheme)
        fd.write (struct.pack ("L", self.nmax))


    def read_file (self, filename) :
        """Read the table from a binary file.
        At present version 3 of the file format is supported. 
        """
        try :
            fd = open (filename, "rb")
        except :
            return False
        stat = True
        try :
            self._read_header_from_stream (fd)
            self._read_table_from_stream (fd)
        except Exception as e :
            print e
        #except :
            stat = False
        finally :
            fd.close()
        return stat

    def read_file_header (self, filename) :
        """Read the table header from a binary file.
        Only the header is read, not the table.  This is useful for getting
        information about the two point tables, such as the pixels in them,
        the bin value, without having to read and decompress the whole file.
        At present version 3 of the file format is supported.
        """
        try :
            fd = open (filename, "rb")
            self._read_header_from_stream (fd)
        except Exception as e :
            print e
            return False
        finally :
            fd.close()
        return True

    def write_file (self, filename) :
        """Write the table from a binary file.
        The pixel list must already be set.
        At present version 3 of the file format is supported. 
        """
        try :
            fd = open (filename, "wb")
        except :
            return False
        stat = True
        try :
            self._write_header_to_stream (fd)
            self._write_table_to_stream (fd)
        # except :
        except Exception as e :
            print "Write failed:", e
            stat = False
        finally :
            fd.close()
        return stat

    def swap_scheme (self) :
        """Swap the HEALPix ordering scheme.  The pixel list will be sorted
        after this operation regardless of its initial state."""
        newpixlist = np.asarray(self.pixlist).copy()
        if self.isNest() :
            pixconvert = healpy.nest2ring
        else :
            pixconvert = healpy.ring2nest
        newpixlist = pixconvert (self.Nside(), newpixlist)
        newpixlist.sort()
        # Make the new table the maximum size.  This makes it easier to
        # convert.  In the end we pull out the rows we want.
        bad_value = 12*self.Nside()**2 + 1 # Impossible value
        newtab = np.ones((12*self.Nside()**2, self.Nmax()), dtype=np.int) \
            * bad_value
        # Loop over existing pixels and do conversion
        for i in range(self.Npix()) :
            p = pixconvert (self.Nside(), self.pixel_list()[i])
            for j in range(self.Nmax()) :
                if self.get(i,j) == -1 : break
                newtab[p,j] = pixconvert (self.Nside(),
                                          self.pixel_list()[self.get(i,j)])
        # Make sure rows are sorted
        newtab.sort(axis=1)
        # Replace the bad_value with -1
        newtab[np.where(newtab == bad_value)] = -1
        # Now reset everything in the class
        self.table = tuple(newtab[newpixlist].flatten())
        self.pixlist = newpixlist.copy()
        if self.isNest() :
            self.scheme = '\x01'
        else :
            self.scheme = '\x00'

    def bin_value (self) :
        return self.cosbin

    def pixel_list (self) :
        return self.pixlist

    def Npix (self) :
        return len(self.pixlist)

    def Nside (self) :
        return self.nside

    def Nmax (self) :
        return self.nmax

    def isNest (self) :
        return self.scheme == '\x00'

    def isRing (self) :
        return not self.isNest()

    def get (self, i, j) :
        return self.table[i*self.Nmax()+j]
