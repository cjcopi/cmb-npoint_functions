#ifndef PIXEL_QUADRILATERALS_H
#define PIXEL_QUADRILATERALS_H

#include <vector>
#include <string>

#include <Pixel_Triangles.h>
#include <healpix_tables.h>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_QUADRILATERALS_RCSID
  ("$Id: Pixel_Quadrilaterals.h,v 1.21 2011-07-27 19:50:00 copi Exp $");
  /// @endcond
}

namespace Npoint_Functions {
  /** Rhombic quadrilaterals.
   *
   *  Rhombic quadrilaterals are constructed from two equilateral triangles
   *  connected along one side, that is, all sides have the same length and
   *  one of the diagonals has the length of the sides.  We use the fact
   *  that the pixels in the triangle are stored in monotonically
   *  increasing order.
   *
   *  Even with this specialization the quad table can be huge.  For this
   *  reason we create a class that incrementally calculates sets of points.
   *  This costs more in overhead but requires significantly less memory.
   */
  template<typename T>
  class Pixel_Quadrilaterals_Rhombic {
  private :
    size_t ind_curr;
    T pixval_end;
    Pixel_Triangles_Equilateral<T> *t;
    std::vector<size_t> skiplist;
  public :  
    Pixel_Quadrilaterals_Rhombic () : ind_curr(0), pixval_end(0),
				      t(0), skiplist() {}
    /// \name Initialize Search
    //@{
    /** Initialize the rhombic quadrilateral search with a triangle.
     *  The given triangle will be used for subsequent searches.  See
     *  next().
     */
    void initialize (Pixel_Triangles_Equilateral<T>& triangle,
		     T pixel_value=-1)
    { 
      ind_curr = 0; t = &triangle;
      /* Create the skip list. Since the actual pixel numbers are stored in
       * triangle the skip list is indexed by pixel number at the relevant
       * Nside. */
      skiplist.assign (12*triangle.Nside()*triangle.Nside(), 0);
      T prev = triangle.get(0,0);
      for (size_t j=1; j < triangle.size(); ++j) {
	if (triangle.get(j,0) != prev) {
	  prev = triangle.get(j,0);
	  skiplist[prev] = j;
	}
      }
      /* There can be 0 in the skip list if pixels do not appear in
       * triangles, so we go through the list again, backwards and fill in
       * the zeros.  The beginning of the list can have zeros so find
       * the first non-zero value.  Also, if the end of the list has zeros
       * we want to fill them in with the max value. */
      size_t ind_start = 1; // [0] must be 0!
      while (skiplist[ind_start] == 0) ++ind_start;
      if (skiplist[skiplist.size()-1] == 0) 
	skiplist[skiplist.size()-1] = triangle.size();
      for (size_t ind=skiplist.size()-2; ind > ind_start; --ind) {
	if (skiplist[ind] == 0) skiplist[ind] = skiplist[ind+1];
      }

      initialize (pixel_value);
    }
    /** Initialize the search for a particular pixel value.
     *  A negative pixel value indicates we want next() to step through all
     *  quadrilaterals. */
    void initialize (T pixel_value=-1)
    {
      if (pixel_value < 0) {
	ind_curr = 0;
	pixval_end = skiplist.size();
      } else {
	ind_curr = skiplist[pixel_value];
	pixval_end = pixel_value;
      }
    }
    //@}

    /// \name Accessors
    //@{
    /** HEALPix Nside of the pixels in the quadrilaterals. */
    inline size_t Nside() const { return t->Nside(); }
    /** HEALPix ordering scheme of the pixels in the quadrilaterals. */
    inline Healpix_Ordering_Scheme Scheme() const { return t->Scheme(); }
    //@}

    /** Get the next set of rhombic quadrilaterals.
     *  The quadrilaterals are constructed for each triangle provided to
     *  initialize(). The quadrilaterals are then made up of the three
     *  points in the triangle, returned in \a pts, and a third
     *  point.  The list of all third points is returned in \a thirdpt.
     *  Note that the orientation of the quadrilateral is lost in this
     *  process.  We do \b not ensure that the quadrilaterals are right
     *  handed.
     *
     *  When the quadrilateral is initialized it is set to either find all
     *  quadrilaterals or only those with a particular pixel index, see
     *  initialize() for details.  Repeated calls to next() will return the
     *  \a true when there are more quadrilaterals to find.
     */
    bool next (std::vector<T>& pts, std::vector<T>& thirdpt)
    {
      if ((ind_curr >= t->size()-1)
	  || (t->get(ind_curr,0) > pixval_end)) return false;
      thirdpt.clear();
      pts.resize(3);
      // Points are not ordered in any special way.
      std::copy (t->get(ind_curr).begin(), t->get(ind_curr).end(),
		 pts.begin());
      // Shorthand
      Npoint_Functions::Orientation o = t->orientation(ind_curr);
      size_t j = ind_curr+1;
      // First loop over triangles with the first two points equal
      while ((j < t->size())
	     && (t->get(j,1) == pts[1])
	     && (t->get(j,0) == pts[0])) {
	if (o != t->orientation(j)) {
	  thirdpt.push_back (t->get(j,2));
	}
	++j;
      }
      /* Next loop over triangles looking for the case when the first and
       * third points are equal to the first and third points of our base
       * triangle. */
      while ((j < t->size())
	     && (t->get(j,1) < pts[2])
	     && (t->get(j,0) == pts[0])) {
	if ((o != t->orientation(j)) && (t->get(j,2) == pts[2])) {
	  thirdpt.push_back (t->get(j,1));
	}
	++j;
      }
      /* Next skip to triangles where the first and third points are equal
       * to the first and second points of our base triangle. */
      while ((j < t->size())
	     && (t->get(j,1) < pts[2])
	     && (t->get(j,0) == pts[0])) ++j;
      // Now loop over these triangles.
      while ((j < t->size())
	     && (t->get(j,1) == pts[2])
	     && (t->get(j,0) == pts[0])) {
	if (o == t->orientation(j)) {
	  thirdpt.push_back (t->get(j,2));
	}
	++j;
      }
      /* Next look for triangles with the second and third points equal to
       * the second and third points of our base triangle. Unfortunately I
       * don't know of a smarter way to do this. */
      T prev;
      while ((j < t->size()) && (t->get(j,0) < pts[1])) {
	prev = t->get(j,0);
	/* For the given value of the first point skip triangles until we
	 * get to one where the second point can possibly match the second
	 * point of the base triangle. */
	while ((j < t->size()) 
	       && (t->get(j,1) < pts[1])
	       && (t->get(j,0) == prev)) ++j;
	/* Now loop over the triangles where the second point matches. */
	while ((j < t->size()) && (t->get(j,1) == pts[1])
	       && (t->get(j,0) == prev)) {
	  if ((o != t->orientation(j)) && 
	      (t->get(j,2) == pts[2])) {
	    thirdpt.push_back(t->get(j,0));
	  }
	  ++j;
	}
	/* Now skip the rest of the triangles where the first point
	 * matches. */
	while ((j < t->size()) && (t->get(j,0) == prev)) ++j;
      }
      /* Now loop over the triangles while its second point is less
       * than the third point of our base triangle. Look for triangles in
       * which the third point equals the third point of our base
       * triangle. */
      while ((j < t->size()) 
	     && (t->get(j,1) < pts[2])
	     && (t->get(j,0) == pts[1])) {
	if ((o == t->orientation(j)) && (t->get(j,2) == pts[2])) {
	  thirdpt.push_back (t->get(j,1));
	}
	++j;
      }
      /* Finally loop over the triangles where the first and second points
       * are equal to the second and third points of our base triangle. */
      while ((j < t->size()) 
	     && (t->get(j,1) == pts[2])
	     && (t->get(j,0) == pts[1])) {
	if (o != t->orientation(j)) {
	  thirdpt.push_back (t->get(j,2));
	}
	++j;
      }
      ++ind_curr;
      return true;
    }
  };
}

#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
