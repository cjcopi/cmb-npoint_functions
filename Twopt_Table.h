#ifndef TWOPT_TABLE_H
#define TWOPT_TABLE_H

#include <vector>
#include <string>
#include <fstream>

namespace {
  /// @cond IDTAG
  const std::string TWOPT_TABLE_RCSID
  ("$Id: Twopt_Table.h,v 1.2 2011-07-07 04:06:19 copi Exp $");
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
 */
template<typename T>
class Twopt_Table {
private :
  std::vector<std::vector<T> > table;
  std::vector<T> pixlist;
  double cosbin;
  size_t nmax;
public :
  /** \name Constructors
   *  Construct a two point table.
   */
  //@{
  /// Generic constructor.
  Twopt_Table () : table(), pixlist(), cosbin(0), nmax(0) {}
  /** Construct a table given the pixel list and the value of the left edge
   *  of the bin.
  */
  Twopt_Table (const std::vector<T>& pl, double binvalue)
    : table(pl.size()), pixlist(pl), cosbin(binvalue), nmax(0) {}
  //@}

  /// Add an entry to the two point table.
  inline void add (const T& i, const T& j)
  { table[i].push_back(j); }

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
      nmax = std::max (nmax, table[p].size());
    }
    out.write (reinterpret_cast<char*>(&nmax), sizeof(nmax));

    // Now write out the values.  Pad with -1 to get all to be the same
    // length.
    T minusone = -1;
    for (size_t p=0; p < Npix; ++p) {
      size_t k;
      for (k=0; k < table[p].size(); ++k) {
	out.write (reinterpret_cast<char*>(&table[p][k]), sizeof(T));
      }
      for (; k < nmax; ++k) {
	out.write (reinterpret_cast<char*>(&minusone), sizeof(T));
      }
    }
    
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
    table.resize(Npix);
    for (size_t p=0; p < Npix; ++p) {
      table[p].resize(nmax);
      for (size_t k=0; k < nmax; ++k) {
	in.read (reinterpret_cast<char*>(&table[p][k]), sizeof(T));
      }
    }
    
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
  { for (size_t p=0; p < table.size(); ++p) table[p].clear(); }

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
  /// Value from the two point table.
  inline T operator() (T i, T j) const
  { return table[i][j]; }
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
