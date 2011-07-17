#ifndef NPOINT_FUNCTIONS_UTILS_H
#define NPOINT_FUNCTIONS_UTILS_H

#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace {
  /// @cond IDTAG
  const std::string NPOINT_FUNCTIONS_UTILS_RCSID
  ("$Id: Npoint_Functions_Utils.h,v 1.2 2011-07-16 23:52:38 copi Exp $");
  /// @endcond
}

/** Namespace for containing all Npoint function related objects.
 */
namespace Npoint_Functions {
  /** Make a numbered filename from a prefix and suffix.
   *  The number is zero padded to a fixed size.  No additional characters
   *  are added to the file name so the suffix MUST include the dot such as
   *  ".dat", if desired.  Similarly the prefix MUST include a separator
   *  character such as an underscore, "prefix_", if desired.
   */
  std::string make_filename (const std::string& prefix, int filenum,
			      int digits=5,
			      const std::string& suffix=".dat")
  {
    std::ostringstream sstr;
    sstr << prefix << std::setw(digits) << std::setfill('0') << filenum
	 << suffix;
    return sstr.str();
  }

  /** Find sequentially numbered files.
   *  A list of existing, sequentially numbered files is generated from the
   *  prefix, digits, and suffix provided.  The filenames are generated by
   *  make_filename() and numbered from \a start in increments of \a
   *  increment.  Existence here means that the file can be opened 
   *  for reading.
   */
  std::vector<std::string>
  get_sequential_file_list (const std::string& prefix,
			    int start, int increment,
			    int digits=5,
			    const std::string& suffix=".dat")
  {
    std::vector<std::string> files;
    std::ifstream in;
    int Nfile = start;
    while (true) {
      std::string fname = make_filename (prefix, Nfile, digits, suffix);
      in.open(fname.c_str());
      if (! in) break;
      in.close();
      files.push_back (fname);
      Nfile += increment;
    }
    return files;
  }

  /** Find sequentially numbered files.
   *  A short hand version for the common case when the file numbers start
   *  at 0 and are incremented by 1.
   */
  std::vector<std::string>
  get_sequential_file_list (const std::string& prefix,
			    int digits=5,
			    const std::string& suffix=".dat")
  { return get_sequential_file_list (prefix, 0, 1, digits, suffix); }


  /** Convert a string to any (valid) type.
   *  The conversion is done using a stringstream and the usual c++ io
   *  mechanism.  This is not the most robust way to do things and it
   *  doesn't allow for complete error checking, however, it is simple
   *  which is what we want here.
   */
  template<typename T>
  bool from_string (const std::string& instr, T& val)
  {
    std::istringstream iss(instr);
    return !(iss >> val).fail();
  }

  /** Generate a range of values.
   *  Sequentially generate a range of values from a starting value in
   *  increments of delta.  This can be used with STL algorithms such as
   *  std::generate().
   *
   *  A sequential list of real numbers from 1 to 2 (inclusive) in steps of
   *  0.1 may be generated as
   * \code
   *  std::vector<double> rseq(11);
   *  std::generate (rseq.begin(), rseq.end(), myRange<double>(1.0, 0.1));
   * \endcode
   */
  template<typename T>
  class myRange {
  private :
    T start, delta, curr, next;
  public :
    /** Construct a range.
     *  By default the range starts at zero and increases in steps of one.
     */
    myRange (T start_=0, T delta_=1)
      : start(start_), delta(delta_), curr(start_), next(start_) {}
    /// Get the next value in the range.
    inline T operator() () { curr=next; next += delta; return curr; }
    /** Reset the range.
     *  The range is reset to its initial value starting value.
     */
    inline void reset() { curr = next = start; }
  };
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
