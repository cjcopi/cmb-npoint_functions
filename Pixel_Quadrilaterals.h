#ifndef PIXEL_QUADRILATERALS_H
#define PIXEL_QUADRILATERALS_H

#include <vector>
#include <string>

#include <Pixel_Triangles.h>
#include <healpix_tables.h>

#include <pixel_ringinfo.h>
#include <pixel_tools.h>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_QUADRILATERALS_RCSID
  ("$Id: Pixel_Quadrilaterals.h,v 1.1 2011-08-09 21:46:47 copi Exp $");
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

  template<typename T>
  class Pixel_Quadrilaterals_Rhombic_Full 
    : public Pixel_Quadrilaterals_Rhombic<T> {
  private :
    // So we don't need to keep recreating it.
    myHealpix::pixel_ringinfo pri;
    // NEED Healpix_Base unless triangle is already in ring scheme
    Healpix_Base HBase;
    // Simple state engine information
    // The base pixel on which we are working.
    enum BasePix { BASE0, BASE4} basepix;
    // Operation to perform at the current pixel
    enum Operation { FINDQUADS, SHIFT, REFLECT1, REFLECT2, REFLECT3}
      operation;
    int optcount;
    // Stored information
    std::vector<T> pts_saved, thirdpt_saved, pts_latest, thirdpt_latest;
    std::vector<T> pixlist;
    size_t ind_curr;

    // Wrappers to handle the scheme for our pixels.
    inline void _pri_setpix (T p)
    {
      if (HBase.Scheme() == NEST) {
	pri.from_pixel (HBase.nest2ring(p));
      } else {
	pri.from_pixel (p);
      }
    }
    inline T _pri_frompix () const {
      if (HBase.Scheme() == NEST) {
	return HBase.nest2ring(pri.to_pixel());
      } else {
	return pri.to_pixel();
      }
    }
    inline T _shift_pix (T p)
    {
      _pri_setpix (p);
      pri.shift_by_base_pixel();
      return _pri_frompix();
    }      
    inline T _reflect_base0 (T p)
    {
      _pri_setpix (p);
      pri.reflect_through_base0();
      return _pri_frompix();
    }
    inline T _reflect_zaxis (T p)
    {
      _pri_setpix (p);
      pri.reflect_through_zaxis();
      return _pri_frompix();
    }
    inline T _reflect_z0 (T p)
    {
      _pri_setpix (p);
      pri.reflect_through_z0();
      return _pri_frompix();
    }

  public :
    Pixel_Quadrilaterals_Rhombic_Full () :
      Pixel_Quadrilaterals_Rhombic<T>(), pri(), HBase(),
      basepix(BASE0), operation(FINDQUADS), optcount(0),
      pts_saved(3), thirdpt_saved(), pts_latest(3), thirdpt_latest(),
      pixlist(), ind_curr(0) {}

    void initialize (Pixel_Triangles_Equilateral<T>& triangle)
    {
      myHealpix::base0_list (triangle.Nside(), pixlist);
      Pixel_Quadrilaterals_Rhombic<T>::initialize (triangle, pixlist[0]);
      basepix = BASE0;
      operation = FINDQUADS;
      optcount = 0;
      ind_curr = 0;
      pri.Nside = triangle.Nside();
      HBase.SetNside (triangle.Nside(), triangle.Scheme());
    }

    bool next (std::vector<T>& pts, std::vector<T>& thirdpt)
    {
      if (operation == SHIFT) {
	for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _shift_pix (pts_latest[j]);
	for (size_t j=0; j < thirdpt_latest.size(); ++j)
	  thirdpt_latest[j] = _shift_pix (thirdpt_latest[j]);
	++optcount;
	if (optcount == 4) operation = REFLECT1;
	else if (optcount == 8) operation = REFLECT2;
	else if (optcount == 12) operation = REFLECT3;
	else if (optcount == 16) {
	  operation = FINDQUADS;
	  optcount = 0;
	}
      } else if (operation == FINDQUADS) {
	if (! Pixel_Quadrilaterals_Rhombic<T>::next (pts, thirdpt)) {
	  // Find the next set of quadrilaterals and set up for all the transformations
	  // It is possible some pixels will not form any quadrilaterals
	  do {
	    ++ind_curr;
	    if (ind_curr >= pixlist.size()) {
	      if (basepix == BASE4) return false;
	      myHealpix::base4_list (this->Nside(), pixlist);
	      basepix = BASE4;
	      ind_curr = 0;
	    }
	    Pixel_Quadrilaterals_Rhombic<T>::initialize (pixlist[ind_curr]);
	  } while (! Pixel_Quadrilaterals_Rhombic<T>::next (pts, thirdpt));
	}
	// Save the pixel info
	std::copy (pts.begin(), pts.end(), pts_saved.begin());
	std::copy (pts.begin(), pts.end(), pts_latest.begin());
	thirdpt_saved.resize (thirdpt.size());
	thirdpt_latest.resize (thirdpt.size());
	std::copy (thirdpt.begin(), thirdpt.end(), thirdpt_saved.begin());
	std::copy (thirdpt.begin(), thirdpt.end(), thirdpt_latest.begin());
	operation = SHIFT;
	optcount = 0;
      } else if (operation == REFLECT1) {
	// Reflect through base0 or z-axis, as appropriate
	if (basepix == BASE0) {
	  for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_base0 (pts_saved[j]);
	  for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_base0 (thirdpt_saved[j]);
	  ++optcount;
	  operation = SHIFT;
	} else { // BASE4
	  for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_zaxis (pts_saved[j]);
	  for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_zaxis (thirdpt_saved[j]);
	  ++optcount;
	  operation = SHIFT;
	}
      } else if (operation == REFLECT2) {
	// Reflect through z=0 line
	for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_z0 (pts_saved[j]);
	for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_z0 (thirdpt_saved[j]);
	  ++optcount;
	  operation = SHIFT;
      } else if (operation == REFLECT3) {
	/* Reflect through z=0 line and through base0 or z-axis, as
	 * appropriate.   This could be further optimized. */
	if (basepix == BASE0) {
	  for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_base0 (pts_saved[j]);
	  for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_base0 (thirdpt_saved[j]);
	  ++optcount;
	  operation = SHIFT;
	} else { // BASE4
	  for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_zaxis (pts_saved[j]);
	  for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_zaxis (thirdpt_saved[j]);

	for (size_t j=0; j < pts_latest.size(); ++j) pts_latest[j] = _reflect_z0 (pts_saved[j]);
	for (size_t j=0; j < thirdpt_latest.size(); ++j)
	    thirdpt_latest[j] = _reflect_z0 (thirdpt_saved[j]);
	  ++optcount;
	  operation = SHIFT;
	}
      }

      std::copy (pts_latest.begin(), pts_latest.end(), pts.begin());
      thirdpt.resize (thirdpt_latest.size());
      std::copy (thirdpt_latest.begin(), thirdpt_latest.end(), thirdpt.begin());
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
