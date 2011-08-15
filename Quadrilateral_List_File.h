#ifndef QUADRILATERAL_LIST_FILE_H
#define QUADRILATERAL_LIST_FILE_H

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <tr1/memory> // For std::tr1::shared_ptr

#include <healpix_map.h>

namespace {
  /// @cond IDTAG
  const std::string QUADRILATERAL_LIST_FILE_RCSID
  ("$Id: Quadrilateral_List_File.h,v 1.2 2011-08-12 22:26:27 copi Exp $");
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

      if (buf != 0) delete [] buf;
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

  // Not sure this really belongs here, but...
  /** Calculate the four point function.
   *  Use a Quadrilateral_List_File to calculate the four point function
   *  for the provided HEALPix map.  It is \b assumed that the scheme of
   *  the map is the same as that of the quadrilateral list.
   *
   *  \relates Quadrilateral_List_File
   */
  template<typename TM, typename TL>
  TM calculate_fourpoint_function (const Healpix_Map<TM>& map,
				   Quadrilateral_List_File<TL>& qlf)
  {
    size_t ind;
    TL p[4];
    TL N[4];
    TL *arr;
    size_t Nquad = 0;
    TM C[4];
    C[0] = 0.0;

    while ((arr = qlf.next()) != 0) {
      ind = 0;
      p[0] = arr[ind++];
      N[1] = arr[ind++];
      C[1] = 0.0;
      for (int n1=0; n1 < N[1]; ++n1) {
	p[1] = arr[ind++];
	N[2] = arr[ind++];
	C[2] = 0.0;
	for (int n2=0; n2 < N[2]; ++n2) {
	  p[2] = arr[ind++];
	  N[3] = arr[ind++];
	  Nquad += N[3];
	  C[3] = 0.0;
	  for (int n3=0; n3 < N[3]; ++n3) {
	    C[3] += map[arr[ind++]];
	  }
	  C[2] += map[p[2]] * C[3];
	}
	C[1] += map[p[1]] * C[2];
      }
      C[0] += map[p[0]] * C[1];
    }

    if (Nquad > 0) C[0] /= Nquad;
    return C[0];
  }

  /** Calculate the four point function for a list of maps.
   *  Use a Quadrilateral_List_File to calculate the four point function
   *  for the provided list of HEALPix maps.  It is \b assumed that the
   *  scheme of the maps is the same as that of the quadrilateral list.
   *
   *  This is a specialized version of calculate_fourpoint_function()
   *  optimized for more than one map at a time.
   *
   *  \relates Quadrilateral_List_File
   */
  template<typename TM, typename TL>
  void calculate_fourpoint_function_list
  (const std::vector<Healpix_Map<TM> >& maps,
   Quadrilateral_List_File<TL>& qlf,
   std::vector<TM>& C4)
  {
    C4.resize(maps.size());
    std::fill (C4.begin(), C4.end(), 0.0);

    size_t ind;
    TL p[4];
    TL N[4];
    TL *arr;
    size_t Nquad = 0;
    // Temporary storage space
    std::vector<TM> C1(maps.size());
    std::vector<TM> C2(maps.size());
    std::vector<TM> C3(maps.size());

    while ((arr = qlf.next()) != 0) {
      ind = 0;
      p[0] = arr[ind++];
      N[1] = arr[ind++];
      std::fill (C1.begin(), C1.end(), 0.0);
      for (int n1=0; n1 < N[1]; ++n1) {
	p[1] = arr[ind++];
	N[2] = arr[ind++];
	std::fill (C2.begin(), C2.end(), 0.0);
	for (int n2=0; n2 < N[2]; ++n2) {
	  p[2] = arr[ind++];
	  N[3] = arr[ind++];
	  Nquad += N[3];
	  std::fill (C3.begin(), C3.end(), 0.0);
	  for (int n3=0; n3 < N[3]; ++n3) {
	    for (size_t j=0; j < C3.size(); ++j)
	      C3[j] += maps[j][arr[ind]];
	    ++ind;
	  }
	  for (size_t j=0; j < C2.size(); ++j)
	    C2[j] += maps[j][p[2]] * C3[j];
	}
	for (size_t j=0; j < C2.size(); ++j)
	  C1[j] += maps[j][p[1]] * C2[j];
      }
      for (size_t j=0; j < C4.size(); ++j)
	C4[j] += maps[j][p[0]] * C1[j];
    }
    if (Nquad > 0) {
      for (size_t j=0; j < C4.size(); ++j)
	C4[j] /= Nquad;
    }
  }
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
