#ifndef BUFFERED_PAIR_BINARY_FILE_H
#define BUFFERED_PAIR_BINARY_FILE_H

#include <string>
#include <fstream>
#include <tr1/memory> // for std::tr1::shared_ptr

namespace {
  /// @cond IDTAG
  const std::string BUFFERED_PAIR_BINARY_FILE_RCSID
  ("$Id: buffered_pair_binary_file.h,v 1.5 2016/02/09 20:31:44 copi Exp $");
  /// @endcond
}

namespace Npoint_Functions {
  /** Buffered binary file for a pair of values.
   *
   *  A binary file is created that stores a sequence of pairs of values.
   *  The reads and writes are internally buffered to cut down on filesystem
   *  io.  The file is written in the byte order of the host machine, nothing
   *  special is don't to make the output portable.  The intent is to use
   *  these for temporary files. */
  template<typename T>
  class buffered_pair_binary_file {
  private :
    size_t buf_size; // Number of entries to buffer.
    std::tr1::shared_ptr<std::fstream> fd;
    std::string fname;
    /* Write buffer information. */
    size_t nbuf_write;
    /* Read buffer information.  These are curr num in the buf, total num in
     * buf, total number in the file. 
     */
    size_t nbuf_read, Nbuf, Ntotal;
    // Buffer, shared for read and write.
    std::tr1::shared_ptr<T> buf;

  public :
    /** Construct a  binary file with a buffer.
     *   The size of the buffer is specified by buf_pairs.  This is the
     *  number of pairs of values to store in the buffer.  This MUST be set
     *  during the intial construction as it cannot be changed.
     */
    buffered_pair_binary_file (const std::string& filename="",
                               size_t buf_pairs=1000000)
      : buf_size(2*buf_pairs), fd(new std::fstream), fname(filename), nbuf_write(0),
        nbuf_read(0), Nbuf(0), Ntotal(0), buf(new T[buf_size])
    {}
    /** Destruct the binary file.
     *  The write buffer is flushed, the file closed, and the buffer freed.
     */
    ~buffered_pair_binary_file ()
    {
      if (fd->is_open()) flush();
    }

    /** Create the buffered file.
     *  The buffered file is created, overwriting the file if it exists., and
     *  opened for writing.
     */
    void create ()
    {
      if (fd->is_open()) fd->close();
      fd->open (fname.c_str(),
                std::fstream::out
                | std::fstream::trunc | std::fstream::binary);
      nbuf_write = 0;
    }

    /** Open the buffered file for reading.
     *  The write buffer is flushed before opening for read.
     */
    void open_read ()
    {
      if (fd->is_open()) {
        flush();
        fd->close();
      }
      // Open the file
      fd->open (fname.c_str(), std::fstream::in | std::fstream::binary);
      // Get total number of entries in the file.
      fd->seekg (0, std::ios::end);
      Ntotal = fd->tellg() / sizeof(T);
      // Rewind
      fd->seekg (0, std::ios::beg);
      nbuf_read = 0;
      Nbuf = 0;
    }

    /** Append a pair of values to the binary file.
     *  The values are buffered internally and only written when the buffer
     *  fills.  To write the values to disk see flush() and close().
     */
    void append (T i, T j)
    {
      if (nbuf_write >= buf_size) flush();
      buf.get()[nbuf_write++] = i;
      buf.get()[nbuf_write++] = j;
    }

    /** Read the next pair of values from the binary file.
     */
    bool read_next_pair (T& i, T& j)
    {
      if (nbuf_read >= Nbuf) {
        // Read another buffer full
        Nbuf = Ntotal - fd->tellg()/sizeof(T);
        if (Nbuf > buf_size) Nbuf = buf_size;
        else if (Nbuf == 0) return false; // End of file
        fd->read (reinterpret_cast<char*>(buf.get()), Nbuf*sizeof(T));
        nbuf_read = 0;
      }
      i = buf.get()[nbuf_read++];
      j = buf.get()[nbuf_read++];
      return true;
    }
  
    /** Flush the write buffer to disk.
     *  The internal buffer is written to disk (which may also be buffered
     *  by the C++ iostream routines).  This routine is safe to call on
     *  files opened for reading (nothing will happen).
     */
    void flush()
    {
      if (nbuf_write > buf_size) {
        std::cerr << "Write buffer too big : " << nbuf_write << std::endl;
        nbuf_write = buf_size;
      }
      if (nbuf_write > 0) {
        fd->write (reinterpret_cast<char*>(buf.get()), nbuf_write*sizeof(T));
      }
      nbuf_write = 0;
    }

    /** Close the binary file.
     *  The write buffer is flushed and the file is closed.
     */
    void close()
    {
      if (fd->is_open()) {
        flush();
        fd->close();
      }
    }

    /** \name Filename
     *  Access and set the filename.
     */
    //@{
    /// Get the filename of the current binary file.
    inline std::string filename() const { return fname; }
    /** Set the name of the binary file.
     *  If an existing file was in use it is first closed.  The new filename
     *  is set but the file is NOT opened.  You must call create() or
     *  open_read() to use the new file.
     */
    inline void filename (std::string& newfile)
    { close(); fname = newfile; }
    //@}
  };
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
