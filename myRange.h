#ifndef MYRANGE_H
#define MYRANGE_H

#include <vector>
#include <string>
#include <fstream>

namespace {
  /// @cond IDTAG
  const std::string MYRANGE_RCSID
  ("$Id$");
  /// @endcond
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

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
