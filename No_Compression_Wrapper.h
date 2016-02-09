#ifndef NO_COMPRESSION_WRAPPER_H
#define NO_COMPRESSION_WRAPPER_H

#include <fstream>

namespace {
  /// @cond IDTAG
  const std::string NO_COMPRESSION_WRAPPER_RCSID
  ("$Id$");
  /// @endcond
}

namespace Npoint_Functions {
  /** Wrapper for no compression.
   *  This implements the generic interface for reading and writing a chunk
   *  of data to a file, in this case with no compression.  See
   *  ZLIB_Wrapper for more details.
   */
  class No_Compression_Wrapper {
  public :
    /// Generic constructor.
    No_Compression_Wrapper() {}

    /** Write the buffer to the stream without compression.
     *  The provided buffer, \a buf_in, of size \a Nbytes is written to
     *  the output stream, \a out, at the current location in the file.
     */
    bool write_buffer (std::ofstream& out,
		       void *buf_in, size_t Nbytes)
    {
      out.write (reinterpret_cast<char*>(buf_in), Nbytes);
      return (! out.fail());
    }

    /** Read the buffer from the stream without compression.
     *  The provided buffer, \a buf_out, of size \a Nbytes is filled from
     *  the input stream, \a in.  The bytes are read from the current
     *  location to the end of the file and returned in \a buf_out.  It is
     *  \b assumed that the stream contains at least \a Nbytes of data.
     */
    bool read_buffer (std::ifstream& in,
		      void *buf_out, size_t Nbytes)
    {
      in.read (reinterpret_cast<char*>(buf_out), Nbytes);
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
