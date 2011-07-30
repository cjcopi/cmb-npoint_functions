import sys
import struct
import zlib

# $Id$

class Twopt_Table :
    """Python class for the storage of a single bin of a two point table.  See
    the c++ version for details.  Only zlib compression integer tables, and
    reading has been implemented."""

    def __init__ (self) :
        self.table_read = []
        self.pixlist = []
        self.cosbin = 0
        self.nside = 0
        self.nmax = 0

    def _read_table_from_stream (self, fd) :
        """Read the table from the stream with compression.
        The Nmax() and Npix() MUST be set correctly before calling.
        Here \a fd is the file stream.  The read table is set after a
        successful call to this this method.
        
        This is used internally by read_file().  It should be used with caution.
        """
        Nelem = self.Nmax()*self.Npix();
        buf = fd.read()
        buf = zlib.decompress(buf)
        self.table_read = struct.unpack ("%di"%Nelem, buf)
        
    def _read_header_from_stream (self, fd) :
        """Read the header from the stream.
        It is assumed the header starts at the current stream position.
        On success thee stream position is left immediately after the
        header. On failure the stream is left in an undefined state.

        This is used internally by read_file().  It should be used with caution.
        """
        version = fd.read (1)
        if version != '\x02' :
            sys.stderr.write ("Twopt_Table only supports file format version 2\n")
            return False

        self.cosbin = struct.unpack ("d", fd.read(8))[0]
        self.nside = struct.unpack ("L", fd.read(8))[0]
        npix = struct.unpack ("L", fd.read(8))[0]
        self.pixlist = struct.unpack ("%di"%(npix),
                                      fd.read(4*npix)) 
        self.nmax = struct.unpack ("L", fd.read(8))[0]


    def read_file (self, filename) :
        """Read the table from a binary file.
        At present version 2 of the file format is supported. 
        """
        try :
            fd = open (filename, "rb")
        except :
            return False
        stat = True
        try :
            self._read_header_from_stream (fd)
            self._read_table_from_stream (fd)
        except :
            stat = False
        finally :
            fd.close()
        return stat

    def read_file_header (self, filename) :
        """Read the table header from a binary file.
        Only the header is read, not the table.  This is useful for getting
        information about the two point tables, such as the pixels in them,
        the bin value, without having to read and decompress the whole file.
        At present version 2 of the file format is supported.
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

    def bin_value (self) :
        return self.cosbin

    def pixel_list (self) :
        return self.pixlist

    def pixel_list (self, ind) :
        return self.pixlist[ind]

    def pixel_list (self) :
        return self.pixlist

    def Npix (self) :
        return len(self.pixlist)

    def Nside (self) :
        return self.nside

    def Nmax (self) :
        return self.nmax

    def get (self, i, j) :
        return self.table_read[i*self.Nmax()+j]

    def get (self) :
        return self.table_read
