#ifndef TWOPT_TABLE_H
#define TWOPT_TABLE_H

#include <vector>
#include <string>
#include <fstream>
#include <tr1/memory> // For std::tr1::shared_ptr

#include <lzma.h>
#include <inttypes.h>

namespace {
  /// @cond IDTAG
  const std::string TWOPT_TABLE_RCSID
  ("$Id: Twopt_Table.h,v 1.4.2.1 2011-07-09 05:04:35 copi Exp $");
  /// @endcond
}

/** Storage for a single bin of a two point table.
 *
 *  A two point table consists of a list of pixels in the NEST scheme, the
 *  value of the lower edge of the bin, and a rectangle table of pixel
 *  indices in the bin.  The size of the table is Npix() x Nmax() where
 *  Nmax() is the maximum number of entries in a row.  The table is "-1"
 *  padded to make it rectangular.
 *
 *  Note that the pixel \b index is stored in the table, not the pixel
 *  number itself.  For a full sky map with the pixels in order these two
 *  are the same, however, for a masked sky or for the pixels not in order
 *  (for some reason) then the pixel index is \b not the same as the pixel
 *  number.  To get the pixel number use the appropriate entry from pixel_list().
 *
 *  Reading and writing two point tables are different processes and are
 *  internally treated differently.  You should not mix reading and writing
 *  of tables.  The intention is to have one code create the tables and
 *  other codes use them.
 *
 *  Internally the table data is stored using liblzma for compression.  The
 *  user does not need to know this since the reading and writing routines
 *  handle it transparently.  However this does mean the files are much
 *  smaller than they would be otherwise.  Also this trades off significant
 *  file io latency for uncompressed files (by far the slowest part of
 *  calculating a two point correlation function) with the necessity for
 *  more CPU power/memory to decompress the data.
 */
template<typename T>
class Twopt_Table {
private :
  // The write table has to be allowed to grow.
  std::vector<std::vector<T> > table_write;
  // The read table is a known size.
  std::tr1::shared_ptr<T> table_read;
  std::vector<T> pixlist;
  double cosbin;
  size_t nmax;

  // LZMA compression options
  static const int compression_level = 6; // 0 to 9
  // Chunk size in bytes to use for compression and decompression
  //static const size_t compression_chunk = 1024*1024;
  /** Write the output table to the stream with compression.
   *  The table and nmax MUST be set correctly before calling.
   */
  bool write_table_to_stream (std::ofstream& out)
  {
    size_t Nelem = Nmax() * Npix();

    // Create the rectangular buffer and -1 fill.
    std::tr1::shared_ptr<T> buf_full (new T [Nelem]);
    std::fill (&buf_full.get()[0], &buf_full.get()[Nelem], -1);
    // Fill in values
    for (size_t p=0; p < Npix(); ++p) {
      for (size_t j=0; j < table_write[p].size(); ++j) {
	buf_full.get()[p*Nmax()+j] = table_write[p][j];
      }
    }

    /* Create the compression buffer.  We know the data CAN be compressed
     * so we set this to be the same size as the full buffer.  This is
     * overkill but should should be safe.
     */
    size_t Nbytes = Nelem * sizeof(T) / sizeof(uint8_t);
    std::tr1::shared_ptr<uint8_t> buf_comp (new uint8_t [Nbytes]);
    lzma_ret ret;
    lzma_stream strm = LZMA_STREAM_INIT;

    // Initialize the buffer
    ret = lzma_easy_encoder (&strm, compression_level, LZMA_CHECK_CRC64);
    if (ret != LZMA_OK) {
      std::cerr << "Error initializing compression buffer : "
                << ret << std::endl;
      return false;
    }
#if 0 // If we need to do the compression in chunks instead of all at once.
    size_t compression_chunk = 1024*1024;
    size_t out_offset;
    size_t Nbytes_compressed;
    out_offset = 0;
    Nbytes_compressed = 0;
    size_t Ncompress = std::min (compression_chunk,
				 Nbytes-Nbytes_compressed);
    strm.next_in  = reinterpret_cast<uint8_t*>(buf_full.get());
    strm.next_out = buf_comp.get();
    while (Ncompress > 0) {
      strm.avail_in = Ncompress;
      strm.avail_out = Nbytes - out_offset;

      ret = lzma_code (&strm, LZMA_RUN);
      if (ret != LZMA_OK) {
	std::cerr << "Error compressing buffer : " << ret << std::endl;
	return false;
      }
      Nbytes_compressed += Ncompress;
      out_offset += Nbytes - out_offset - strm.avail_out;
      Ncompress = std::min (compression_chunk,
			    Nbytes-Nbytes_compressed);
    };
#else
    strm.next_in  = reinterpret_cast<uint8_t*>(buf_full.get());
    strm.next_out = buf_comp.get();
    strm.avail_in = Nbytes;
    strm.avail_out = Nbytes;

    ret = lzma_code (&strm, LZMA_RUN);
    if (ret != LZMA_OK) {
      std::cerr << "Error compressing buffer : " << ret << std::endl;
      return false;
    }
#endif
    ret = lzma_code (&strm, LZMA_FINISH);
    if ((ret != LZMA_OK) && (ret != LZMA_STREAM_END)) {
      std::cerr << "Error cleaning up compression buffer : "
                << ret << std::endl;
      return false;
    }
    out.write (reinterpret_cast<char*>(buf_comp.get()), strm.total_out);
    lzma_end (&strm);
    return true;
  }

  /** Read the table from the stream with compression.
   *  The Nmax() and Npix() MUST be set correctly before calling.
   */
  bool read_table_from_stream (std::ifstream& in,
			       std::tr1::shared_ptr<T> *buf)
  {
    *buf = std::tr1::shared_ptr<T>(new T [Nmax()*Npix()]);
    // First read in compressed values from the file.
    std::streampos curpos = in.tellg();
    in.seekg (0, std::ios::end);
    size_t in_len = in.tellg() - curpos;
    in.seekg (curpos, std::ios::beg);
    std::tr1::shared_ptr<uint8_t> buf_comp (new uint8_t [ in_len ]);
    in.read (reinterpret_cast<char*>(buf_comp.get()), in_len);

    // Now decompress and fill in buf
    lzma_stream strm = LZMA_STREAM_INIT;
    lzma_ret ret;
    ret = lzma_stream_decoder (&strm, UINT64_MAX, 0);
    if (ret != LZMA_OK) {
      std::cerr << "Error initializing decompression : "
                << ret << std::endl;
      return false;
    }
    strm.next_in = buf_comp.get();
    strm.avail_in = in_len;
    strm.next_out = reinterpret_cast<uint8_t*>(buf->get());
    strm.avail_out = Nmax()*Npix() * sizeof(T)/sizeof(uint8_t);

    ret = lzma_code (&strm, LZMA_RUN);
    if ((ret != LZMA_OK) && (ret != LZMA_STREAM_END)) {
      std::cerr << "Error decompressing buffer : " << ret << std::endl;
      return false;
    }
    ret = lzma_code (&strm, LZMA_FINISH);
    if ((ret != LZMA_OK) && (ret != LZMA_STREAM_END)) {
      std::cerr << "Error cleaning up decompression buffer : "
                << ret << std::endl;
      return false;
    }
    lzma_end (&strm);
    return true;
  }
public :
  /** \name Constructors
   *  Construct a two point table.
   */
  //@{
  /// Generic constructor.
  Twopt_Table () : table_write(), table_read(), pixlist(), cosbin(0), nmax(0) {}
  /** Construct a table given the pixel list and the value of the left edge
   *  of the bin.
  */
  Twopt_Table (const std::vector<T>& pl, double binvalue)
    : table_write(pl.size()), table_read(), pixlist(pl), cosbin(binvalue), nmax(0) {}
  //@}

  /// Add an entry to the two point table.
  inline void add (const T& i, const T& j)
  { table_write[i].push_back(j); }

  /** Add a pair symmetrically to the two point table.
   *  This is equivalent to calling add() twice for the pairs \a i,\a j and
   *  \a j,\a i.
   */
  inline void add_pair (const T& i, const T& j)
  { add(i,j); add(j,i); }

  /** Write the table to a binary file.
   *  At present version 1 of the file format is written.  This format is
   * version number (char)
   * bin value (double)
   * Npix (size_t)
   * list of pixels (Npix of them of type T)
   * Nmax (size_t)
   * table values (Npix x Nmax of them of type T written in row major order)
  *
  * The table is -1 padded to make it rectangular.
  */
  void write_file (const std::string& filename)
  {
    char version = 1;
    size_t Npix = pixlist.size();
    std::ofstream out (filename.c_str(),
		       std::fstream::out | std::fstream::trunc
		       | std::fstream::binary);
    // First header
    out.write (&version, sizeof(version));
    out.write (reinterpret_cast<char*>(&cosbin), sizeof(cosbin));
    out.write (reinterpret_cast<char*>(&Npix), sizeof(Npix));
    for (size_t p=0; p < Npix; ++p) {
      out.write (reinterpret_cast<char*>(&pixlist[p]), sizeof(T));
    }
    // Now figure out what the maximum number of values in a pixel bin are
    nmax = 0;
    for (size_t p=0; p < Npix; ++p) {
      nmax = std::max (nmax, table_write[p].size());
    }
    out.write (reinterpret_cast<char*>(&nmax), sizeof(nmax));

    // Now write out the values.
    write_table_to_stream (out);
    out.close();
  }

  /** Read the table from a binary file.
   *  At present version 1 of the file format is supported.  See
   *  write_file() for details.
   */
  bool read_file (const std::string& filename)
  {
    char version;
    size_t Npix;
    std::ifstream in (filename.c_str(),
		      std::fstream::in | std::fstream::binary);
    if (! in) return false;

    int nmax_old = nmax;
    // First header
    in.read (&version, sizeof(version));
    if (version != 1) {
      std::cerr << "Twopt_Table only supports file format version 1\n";
      return false;
    }
    in.read (reinterpret_cast<char*>(&cosbin), sizeof(cosbin));
    in.read (reinterpret_cast<char*>(&Npix), sizeof(Npix));
    pixlist.resize(Npix);
    for (size_t p=0; p < Npix; ++p) {
      in.read (reinterpret_cast<char*>(&pixlist[p]), sizeof(T));
    }
    in.read (reinterpret_cast<char*>(&nmax), sizeof(nmax));
    read_table_from_stream (in, &table_read);
#if 0
    if (nmax != nmax_old)
      table_read = std::tr1::shared_ptr<T>(new T [Npix*nmax]);
    in.read (reinterpret_cast<char*>(table_read.get()), Npix*nmax*sizeof(T));
#endif

    in.close();
    return true;
  }

  /** Reset the two point table.
   *  This ONLY clears the table.  The pixel list and bin value are
   *  unchanged.  This is useful to call after writing a table with
   *  write_file() to prepare for filling in another bin in the two point
   *  table.
   */
  inline void reset()
  { for (size_t p=0; p < table_write.size(); ++p) table_write[p].clear(); }

  /** \name Accessors
   *  Access internal information.
   */
  //@{
  /// The value of the left edge of the bin.
  inline double bin_value () const { return cosbin; }
  /// The list of pixels.
  inline const std::vector<T>& pixel_list () const { return pixlist; }
  /// The number of pixels.
  inline size_t Npix() const { return pixlist.size(); }
  /// The maximum number of values in each row of the table.
  inline size_t Nmax() const { return nmax; }
  /// Value from the two point write table.
  inline T& get_write_value (T i, T j)
  { return table_write[i][j]; }
  /** Value from the two point read table.
   *  This value cannot be changed.
   */
  inline const T& get_read_value (T i, T j) const
  { return table_read.get()[i*Nmax()+j]; }
  //@}

  /// Assign the value of the left edge of the bin.
  inline void bin_value (double bv) { cosbin=bv; }
  /// Assign the list of pixels.
  inline void pixel_list (const std::vector<T>& pl) 
  { 
    pixlist.resize (pl.size());
    std::copy (pl.begin(), pl.end(), pixlist.begin());
  }
};

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
