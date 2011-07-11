#ifndef PIXEL_TRIANGLES_H
#define PIXEL_TRIANGLES_H

#include <vector>
#include <string>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_TRIANGLES_RCSID
  ("$Id: Pixel_Triangles.h,v 1.1 2011-07-11 20:45:06 copi Exp $");
  /// @endcond
}

/** Storage for pixel triangles.
 *  All possible triangles are stored, including cyclic permutations of
 *  triangle with the same side lengths.  See Pixel_Triangles_Isosceles or
 *  Pixel_Triangles_Equilateral for specialized versions.
 */
template<typename T>
class Pixel_Triangles {
private :
  std::vector<std::vector<T> > triangles;
protected :
  /** Find matches in two lists and append them to a new list. */
  void append_matches (const T *L1, size_t NL1, const T *L2, size_t NL2,
		       std::vector<T>& res)
  {
    const T *it1, *it2; // "iterators"
    
    it1 = L1;
    it2 = L2;
    
    /* Now loop over the iterators storing matches
     * Since the lists are monotonically increasing and -1 padded at the end
     * a simple linear search is an efficient algorithm.
     */
    while ( (it1 != &L1[NL1]) && (it2 != &L2[NL2])
	    && (*it1 != -1) && (*it2 != -1) ) {
      if (*it1 == *it2) {
	res.push_back (*it1);
	++it1;
	++it2;
      } else if (*it1 < *it2) {
	++it1;
      } else {
	++it2;
      }
    }
  }
public :
  /// Generic constructor.
  Pixel_Triangles () : triangles() {}

  /// Add a triangle
  inline void add (const T& p1, const T& p2, const T& p3)
  {
    size_t n = triangles.size();
    triangles.push_back (std::vector<T>(3));
    triangles[n][0] = p1;
    triangles[n][1] = p2;
    triangles[n][2] = p3;
  }

  /** Reset the list of triangles.
   *  All triangles are erased.
   */
  inline void reset()
  { triangles.clear(); }

  /** Find all triangles.
   *  Find all the triangles that can be made up from the provided two
   *  point tables.  It is assumed that all the two point tables are
   *  different, thus if this is used to calculate isosceles or
   *  equilaterial triangles then all the cyclic permutations will be found
   *  as separate triangles.
   */
  void find_triangles (const Twopt_Table<T>& t1,
		       const Twopt_Table<T>& t2,
		       const Twopt_Table<T>& t3)
  {
    T p1, p2, p3;
    T i1, i2;
    std::vector<T> trip;
    size_t ntriplets=0; // Keep track of triplets, just for convenience

    for (size_t j1=0; j1 < t1.Npix(); ++j1) {
      i1 = j1; // to make the code look symmetric
      p1 = t1.pixel_list()[j1];
      for (size_t j2=0; (j2 < t1.Nmax()) && (t1(j1,j2) != -1); ++j2) {
	i2 = t1(j1,j2);
	p2 = t1.pixel_list()[i2];
	// Finally can search for and add appropriate pairs.
	trip.clear();
	append_matches (&t2(i1,0), t2.Nmax(), &t3(i2,0), t3.Nmax(),
			trip);
	// Now put all the triplets in the list
	for (size_t k=0; k < trip.size(); ++k) {
	  add (p1, p2, t1.pixel_list()[trip[k]]);
	}
      }
    }
  }

  /** \name Accessors
   *  Access internal information.
   */
  //@{
  /// Number of triangles in the list.
  inline size_t size() const
  { return triangles.size(); }
  /** The three pixels that are the corners of the requested triangle.
   *  This value cannot (should not) be changed.
   */
  inline const std::vector<T>& operator() (size_t j) const
  { return triangles[j]; }
  //@}
};


/** Storage for isosceles pixel triangles.
 *  Only the unique triangles are stored.  The angular distance between
 *  pixel pairs 1,2 and 1,3 are equal.  The angular distance between pixel
 *  pair 2,3 is different than the other two pairs.  Pixel 2 is always less
 *  than pixel 3 (these two pixels are interchangeable).
 *
 *  This is a specialized version of Pixel_Triangles.
 */
template<typename T>
class Pixel_Triangles_Isosceles : public Pixel_Triangles<T> {
protected :
  /** Find matches in two lists and append them to a new list.
   *  Here the minimum allowed value is provided.  All values appended to the
   *  list will be greater than or equal to this value. */
  void append_matches (T minval, const T *L1, size_t NL1,
		       const T *L2, size_t NL2,
		       std::vector<T>& res)
  {
    const T *it1, *it2; // "iterators"
    // For some reason this can fail?  Is it because of the -1 padding?
    //it1 = std::lower_bound (&L1[0], &L1[NL1], minval);
    //it2 = std::lower_bound (&L2[0], &L2[NL2], minval);
    it1 = L1;
    while ((it1 != &L1[NL1]) && (*it1 < minval)) ++it1;
    it2 = L2;
    while ((it2 != &L2[NL2]) && (*it2 < minval)) ++it2;

    /* Now loop over the iterators storing matches
     * Since the lists are monotonically increasing and -1 padded at the end
     * a simple linear search is an efficient algorithm.
     */
    while ( (it1 != &L1[NL1]) && (it2 != &L2[NL2])
	    && (*it1 != -1) && (*it2 != -1) ) {
      if (*it1 == *it2) {
	res.push_back (*it1);
	++it1;
	++it2;
      } else if (*it1 < *it2) {
	++it1;
      } else {
	++it2;
      }
    }
  }
public :
  /// Generic constructor.
  Pixel_Triangles_Isosceles () : Pixel_Triangles<T>() {}

  /** Find all isosceles triangles.
   *  Find all the triangles that can be made up from the provided two
   *  point tables.  The first table is the one for the two equal sides.
   *  The first two pixels are in monotonically increasing pixel index
   *  order.  The the angular distance between pixel pairs 1,2 and 1,3 are
   *  equal.  The angular distance between pixel pair 2,3 is different
   *  than the other two pairs.
   */
  void find_triangles (const Twopt_Table<T>& t1,
		       const Twopt_Table<T>& t2)
  {
    T p1, p2, p3;
    T i1, i2;
    std::vector<T> trip;
    size_t ntriplets=0; // Keep track of triplets, just for convenience

    for (size_t j1=0; j1 < t1.Npix(); ++j1) {
      i1 = j1; // to make the code look symmetric
      p1 = t1.pixel_list()[j1];
      for (size_t j2=0; (j2 < t1.Nmax()) && (t1(j1,j2) != -1); ++j2) {
	i2 = t1(j1,j2);
	p2 = t1.pixel_list()[i2];
	// Finally can search for and add appropriate pairs.
	trip.clear();
	append_matches (i2, &t1(i1,0), t1.Nmax(), &t2(i2,0), t2.Nmax(),
			trip);
	// Now put all the triplets in the list
	for (size_t k=0; k < trip.size(); ++k) {
	  add (p1, p2, t1.pixel_list()[trip[k]]);
	}
      }
    }
  }
};

/** Storage for equilateral pixel triangles.
 *  Only the unique triangles are stored.  The pixels are stored in
 *  monotonically increasing order.
 *
 *  This is a specialized version of Pixel_Triangles_Isosceles.
 */
template<typename T>
class Pixel_Triangles_Equilateral : public Pixel_Triangles_Isosceles<T> {
public :
  /// Generic constructor.
  Pixel_Triangles_Equilateral () : Pixel_Triangles_Isosceles<T>() {}

  /** Find all equilateral triangles.
   *  Find all the triangles that can be made up from the provided two
   *  point table.  The triangles are stored in monotonically increasing
   *  pixel index order.
   */
  void find_triangles (const Twopt_Table<T>& t)
  {
    T p1, p2, p3;
    T i1, i2;
    std::vector<T> trip;
    size_t ntriplets=0; // Keep track of triplets, just for convenience

    for (size_t j1=0; j1 < t.Npix(); ++j1) {
      i1 = j1; // to make the code look symmetric
      p1 = t.pixel_list()[j1];
      for (size_t j2=0; (j2 < t.Nmax()) && (t(j1,j2) != -1); ++j2) {
	i2 = t(j1,j2);
	p2 = t.pixel_list()[i2];
	if (p2 < p1) continue;
	// Finally can search for and add appropriate pairs.
	trip.clear();
	append_matches (i2, &t(i1,0), t.Nmax(), &t(i2,0), t.Nmax(),
			trip);
	// Now put all the triplets in the list
	for (size_t k=0; k < trip.size(); ++k) {
	  add (p1, p2, t.pixel_list()[trip[k]]);
	}
      }
    }
  }
};

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
