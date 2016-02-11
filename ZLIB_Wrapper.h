#ifndef ZLIB_WRAPPER_H
#define ZLIB_WRAPPER_H

#include <fstream>
#include <tr1/memory> // For std::tr1::shared_ptr

#include <zlib.h>

namespace {
  /// @cond IDTAG
  const std::string ZLIB_WRAPPER_RCSID
  ("$Id: ZLIB_Wrapper.h,v 1.5 2016/02/09 20:31:44 copi Exp $");
  /// @endcond
}

namespace Npoint_Functions {
  /** Wrapper for simple zlib compression.
   *  This is a generic interface that allows reading and writing chunks of
   *  data to a stream with compression.  This allows easy replacement of the
   *  compression routines by writting a wrapper with the same interface that
   *  uses a different compression library.
   */
  class ZLIB_Wrapper {
  private :
    // ZLIB compression options
    static const int compression_level = 6; // 0 to 9
  public :
    /// Generic constructor.
    ZLIB_Wrapper() {}

    /** Write the buffer to the stream with compression.
     *   The provided buffer, \a buf_in, of size \a Nbytes is compressed and
     *  written to the output stream, \a out, at the current location in the
     *  file. 
     */
    bool write_buffer (std::ofstream& out,
                       void *buf_in, size_t Nbytes)
    {
      if (Nbytes == 0) return true;

      /* Create the compression buffer.  We know the data CAN be compressed
       * so we set this to be the same size as the full buffer.  This is
       * overkill but should should be safe.
       */
      std::tr1::shared_ptr<unsigned char> buf_comp (new unsigned char [Nbytes]);
      int ret;
      z_stream strm;

      // Initialize the buffer
      strm.zalloc = Z_NULL;
      strm.zfree = Z_NULL;
      strm.opaque = Z_NULL;
      ret = deflateInit (&strm, compression_level);
      if (ret != Z_OK) {
        std::cerr << "Error initializing compression buffer : "
                  << ret << std::endl;
        return false;
      }

      strm.next_in  = reinterpret_cast<unsigned char*>(buf_in);
      strm.next_out = buf_comp.get();
      strm.avail_in = Nbytes;
      strm.avail_out = Nbytes;

      ret = deflate (&strm, Z_NO_FLUSH);
      if ((ret != Z_OK) && (ret != Z_STREAM_END)) {
        std::cerr << "Error compressing buffer : " << ret << std::endl;
        deflateEnd (&strm);
        return false;
      }

      ret = deflate (&strm, Z_FINISH);
      if (ret != Z_STREAM_END) {
        std::cerr << "Error cleaning up compression buffer : "
                  << ret << std::endl;
        deflateEnd (&strm);
        return false;
      }
      out.write (reinterpret_cast<char*>(buf_comp.get()),
                 Nbytes - strm.avail_out);
      deflateEnd (&strm);
      return (! out.fail());
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
      if (in_len == 0) return true;
      std::tr1::shared_ptr<unsigned char> buf_comp (new unsigned char [ in_len ]);
      in.read (reinterpret_cast<char*>(buf_comp.get()), in_len);

      // Now decompress and fill in buf_out
      z_stream strm;
      int ret;
      strm.zalloc = Z_NULL;
      strm.zfree = Z_NULL;
      strm.opaque = Z_NULL;
      strm.avail_in = 0;
      strm.next_in = Z_NULL;
      ret = inflateInit(&strm);
      if (ret != Z_OK) {
        std::cerr << "Error initializing decompression : "
                  << ret << std::endl;
        return false;
      }
      strm.next_in = buf_comp.get();
      strm.avail_in = in_len;
      strm.next_out = reinterpret_cast<unsigned char*>(buf_out);
      strm.avail_out = Nbytes;

      ret = inflate (&strm, Z_NO_FLUSH);
      if (ret != Z_STREAM_END) {
        std::cerr << "Error decompressing buffer : " << ret << std::endl;
        inflateEnd (&strm);
        return false;
      }
      inflateEnd (&strm);
      return (! in.fail());
    }
  };
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
