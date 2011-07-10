#ifndef LZMA_WRAPPER_H
#define LZMA_WRAPPER_H

#include <fstream>
#include <tr1/memory> // For std::tr1::shared_ptr

#include <lzma.h>
#include <inttypes.h>

namespace {
  /// @cond IDTAG
  const std::string LZMA_WRAPPER_RCSID
  ("$Id: LZMA_Wrapper.h,v 1.1 2011-07-09 21:24:53 copi Exp $");
  /// @endcond
}

/** Wrapper for simple lzma compression.
 *  This is a generic interface that allows reading and writing chunks of
 *  data to a stream with compression.  This allows easy replacement of the
 *  compression routines by writting a wrapper with the same interface that
 *  uses a different compression library.
 *
 *  Lzma produces smaller files zlib but is much slower.  This wrapper can
 *  be used as a drop-in replacement for ZLIB_Wrapper if smaller files are
 *  very important.
 */
class LZMA_Wrapper {
private :
  // LZMA compression options
  static const int compression_level = 6; // 0 to 9
public :
  /// Generic constructor.
  LZMA_Wrapper() {}

  /** Write the buffer to the stream with compression.
   *   The provided buffer, \a buf_in, of size \a Nbytes is compressed and
   *  written to the output stream, \a out, at the current location in the
   *  file. 
   */
  bool write_buffer (std::ofstream& out,
		     const void *buf_in, size_t Nbytes)
  {
    /* Create the compression buffer.  We know the data CAN be compressed
     * so we set this to be the same size as the full buffer.  This is
     * overkill but should should be safe.
     */
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

    strm.next_in  = reinterpret_cast<const uint8_t*>(buf_in);
    strm.next_out = buf_comp.get();
    strm.avail_in = Nbytes;
    strm.avail_out = Nbytes;

    ret = lzma_code (&strm, LZMA_RUN);
    if (ret != LZMA_OK) {
      std::cerr << "Error compressing buffer : " << ret << std::endl;
      return false;
    }

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

  /** Read the buffer from the stream with compression.
   *  The provided buffer, \a buf_out, of size \a Nbytes is filled from
   *  the input stream, \a in.  The compressed bytes are read from the
   *  current location to the end of the file, decompressed, and returned
   *  in \a buf_out.
   */
  bool read_buffer (std::ifstream& in,
		    void *buf_out, size_t Nbytes)
  {
    // First read in compressed values from the file.
    std::streampos curpos = in.tellg();
    in.seekg (0, std::ios::end);
    size_t in_len = in.tellg() - curpos;
    in.seekg (curpos, std::ios::beg);
    std::tr1::shared_ptr<uint8_t> buf_comp (new uint8_t [ in_len ]);
    in.read (reinterpret_cast<char*>(buf_comp.get()), in_len);

    // Now decompress and fill in buf_out
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
    strm.next_out = reinterpret_cast<uint8_t*>(buf_out);
    strm.avail_out = Nbytes;

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
};

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
