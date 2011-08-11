#ifndef QUADRILATERAL_LIST_FILE_H
#define QUADRILATERAL_LIST_FILE_H

#include <vector>
#include <string>
#include <fstream>
#include <tr1/memory> // For std::tr1::shared_ptr

namespace {
  /// @cond IDTAG
  const std::string QUADRILATERAL_LIST_FILE_RCSID
  ("$Id: Pixel_Quadrilaterals.h,v 1.4 2011-08-11 16:23:48 copi Exp $");
  /// @endcond
}

namespace Npoint_Functions {
  /** List of pixels for quadrilaterals stored in a compressed format.  *
   *  This is a "raw" class providing a wrapper around the file format used
   *  to store lists of quadrilaterals.  It only provides read access to
   *  the file.
   */
  template<typename T>
  class Quadrilateral_List_File {
  private :
    size_t nside;
    Healpix_Ordering_Scheme scheme;
    double binval;
    std::tr1::shared_ptr<std::ifstream> fd;
    T *buf;
  public :  
    /** Constructor.
     *  If a filename is provided the class is initialized and ready for
     *  use. */
    Quadrilateral_List_File (const std::string& filename="")
      : nside(0), scheme(NEST), binval(0.0),
	fd(new std::ifstream), buf(0)
    { if (filename != "") initialize (filename); }

    /** Destructor.
     *  Close file and free memory. */
    ~Quadrilateral_List_File ()
    {
      if (fd->is_open()) fd->close();
      if (buf != 0) delete [] buf;
    }

    /** Initialize. 
     *  The file is opened, read, and prepared for use.  See next() for
     *  usage.  On error the file is left in an indeterminant state.
     */
    bool initialize (const std::string& filename)
    { 
      if (fd->is_open()) fd->close();
      fd->open (filename.c_str(), std::fstream::in | std::fstream::binary);
      if (! *fd) {
	std::cerr << "Failed to open file " << filename << std::endl;
	return false;
      }
      
      // Read header.
      char version, s;
      size_t maxbytes;
      fd->read (&version, sizeof(version));
      if (version != 1) {
	std::cerr << "Only version 1 supported\n";
	return false;
      }
      fd->read (reinterpret_cast<char*>(&nside), sizeof(nside));
      fd->read (&s, sizeof(s));
      if (s == 0) scheme = NEST;
      else scheme = RING;
      fd->read (reinterpret_cast<char*>(&binval), sizeof(binval));
      fd->read (reinterpret_cast<char*>(&maxbytes), sizeof(maxbytes));

      buf = new T [ maxbytes/sizeof(T) ];

      return true;
    }

    /** Get the next set of quadrilaterals to process.

     * A pointer to the memory is returned.  Do \b not free this memory or
     * problems will ensue (I warned you that this is a raw class).  The
     * quadrilaterals are stored as a list of values in a "recursive"
     * order.  The list of numbers is in the format:
     * p0 Np1 { p1 Np2 [ p2 Np3 (p3 p3 ...) p2 Np3 (p3 ...) ... ] p1 Np2
     * [p2 ...] ... }.
     *
     * When no more quadrilaterals are available "0" is returned.
     */
    T* next()
    {
      size_t bytes;

      fd->read (reinterpret_cast<char*>(&bytes), sizeof(bytes));
      if (! *fd) return 0;

      fd->read (reinterpret_cast<char*>(buf), bytes);
      return buf;
    }

    /// \name Accessors
    //@{
    /** Nside of the pixels in the quadrilateral list. */
    size_t Nside() const { return nside; }
    /** Scheme of the pixels in the quadrilateral list. */
    Healpix_Ordering_Scheme Scheme() const { return scheme; }
    /** Value at the center of the bin for this quadrilateral list.
     *  This is specific to rhombic quadrilaterals .... */
    double bin_value() const { return binval; }
    //@}
  };
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
