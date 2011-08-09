#ifndef PIXEL_TRIANGLES_H
#define PIXEL_TRIANGLES_H

#include <vector>
#include <string>

#include <healpix_base.h>
#include <healpix_tables.h>
#include <vec3.h>

#include <Npoint_Functions_Utils.h>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_TRIANGLES_RCSID
  ("$Id: Pixel_Triangles.h,v 1.21 2011-07-27 19:50:00 copi Exp $");
  /// @endcond
}

namespace {
    /** Find matches in two lists and append them to a new list. */
  template<class IteratorType, typename T>
  void append_matches (IteratorType it1, IteratorType it1end,
		       IteratorType it2, IteratorType it2end,
		       std::vector<T>& res)
  {
    /* Loop over the iterators storing matches.
     * Since the lists are monotonically increasing and -1 padded at the end
     * a simple linear search is an efficient algorithm.
     */
    while ( (it1 != it1end) && (it2 != it2end)
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
  
  /* Find matches in two lists and append them to a new list.
   *  Here the minimum allowed value is provided.  All values appended to the
   *  list will be greater than or equal to this value. */
  template<class IteratorType, typename T>
  void append_matches (T minval,
		       IteratorType it1, IteratorType it1end,
		       IteratorType it2, IteratorType it2end,
		       std::vector<T>& res)
  {
    while ((it1 != it1end) && (*it1 < minval)) ++it1;
    while ((it2 != it2end) && (*it2 < minval)) ++it2;
    append_matches (it1, it1end, it2, it2end, res);
  }
}

namespace Npoint_Functions {
    /** Allowed orientations of a triangle.
     *  See calculate_orientation() for conventions.
     *  \relates Pixel_Triangles
     */
    enum Orientation { RIGHTHANDED, LEFTHANDED };

    /** Calculate the Orientation from three vectors.
     *  The orientation is either righthanded or lefthanded.
     *  Righthanded is defined by 
     *  \f[ (\hat n_1\times\hat n_2)\cdot \hat n_3 > 0. \f]
     *  \relates Pixel_Triangles
     */
    Orientation calculate_orientation (const vec3& n1, const vec3& n2,
				       const vec3& n3)
    {
      double val = dotprod (crossprod (n1, n2), n3);
      return ((val > 0) ? RIGHTHANDED : LEFTHANDED);
    }

  /** Storage for pixel triangles.
   *  All possible triangles are stored, including cyclic permutations of
   *  triangle with the same side lengths.  See Pixel_Triangles_Isosceles or
   *  Pixel_Triangles_Equilateral for specialized versions.
   *
   *  The actual pixel values are stored, not the indices to the pixel list
   *  as is done in the two point table.
   */
  template<typename T>
  class Pixel_Triangles {
  private :
    std::vector<std::vector<T> > triangles; // List of pixels in triangle.
    std::vector<double> edge_length; // Length of triangle edges.
     /// Orientation of the triangles.
    std::vector<Orientation> orient;
    // Vectors to the center of HEALPix pixels.
    std::vector<vec3> v;
    // HEALPix Nside of the pixels in the triangles.
    size_t nside;
    // HEALPix ordering scheme for the pixels in the triangles.
    Healpix_Ordering_Scheme scheme;
  protected :
    /// Add a triangle to the list.
    inline void add (const T& p1, const T& p2, const T& p3)
    {
      size_t n = triangles.size();
      triangles.push_back (std::vector<T>(3));
      triangles[n][0] = p1;
      triangles[n][1] = p2;
      triangles[n][2] = p3;
      orient.push_back (calculate_orientation (v[p1], v[p2],
					       v[p3]));
    }

    /// Set the edge lengths of the triangle
    inline void set_edge_lengths (double l1, double l2, double l3)
    {
      this->edge_length[0] = l1;
      this->edge_length[1] = l2;
      this->edge_length[2] = l3;
    }


    /** Internal routine for initializing the state of the class for a set
     *  of two point tables.  The list of vectors to the HEALPix pixel
     *  centers is also calculated.
     */
    void initialize (const Twopt_Table<T>& t1,
		     const Twopt_Table<T>& t2,
		     const Twopt_Table<T>& t3)
    {
      this->reset();
      set_edge_lengths (t1.bin_value(), t2.bin_value(), t3.bin_value());
      this->nside = t1.Nside();
      this->scheme = t1.Scheme();
      fill_vector_list (t1.Nside(), t1.Scheme(), v);
    }

  public :
    /// Generic constructor.
    Pixel_Triangles () : triangles(), edge_length(3), orient(), v(),
			 nside(0), scheme(NEST) {}

    /** Reset the list of triangles.
     *  All triangles are erased.
     */
    inline void reset()
    { 
      triangles.clear();
      v.clear();
      orient.clear();
    }

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
      T p1, p2;
      T i2;
      std::vector<T> trip;

      this->initialize (t1, t2, t3);

      for (size_t i1=0; i1 < t1.Npix(); ++i1) {
	p1 = t1.pixel_list(i1);
	for (size_t j2=0; (j2 < t1.Nmax()) && (t1(i1,j2) != -1); ++j2) {
	  i2 = t1(i1,j2);
	  p2 = t1.pixel_list(i2);
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (&t2(i1,0), &t2(i1,t2.Nmax()),
			  &t3(i2,0), &t3(i2,t3.Nmax()), trip);
	  // Now put all the triplets in the list.
	  for (size_t k=0; k < trip.size(); ++k) {
	    add (p1, p2, t1.pixel_list(trip[k]));
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
    /** The three pixels that are the corners of the requested triangle.
     *  This value cannot (should not) be changed.
     */
    inline const std::vector<T>& get (size_t j) const
    { return triangles[j]; }
    /** Pixel index \a j of triangle \a i. */
    inline T get (size_t i, size_t j) const
    { return triangles[i][j]; }
    /** The Orientation of the triangle.
     *   See calculate_orientation() for more details. 
     */
    inline Orientation orientation (size_t j) const
    { return orient[j]; }
    /** The edge lengths of the triangles.
     *  The edge lengths are the dot products between the vectors to the
     *  points of the triangle in the order
     * \f[ \left\{ \hat n_1\cdot\hat n_2, \quad \hat n_2\cdot\hat n_3, \quad \hat
     *     n_3\cdot\hat n_1 \right\}. \f]
     */
    inline const std::vector<double>& lengths() const
    { return edge_length; }
    /** Nside of the pixels in the triangles.
     */
    inline int Nside() const { return nside; }
    /** HEALPix ordering scheme of the pixels in the triangles.
     */
    inline Healpix_Ordering_Scheme Scheme() const { return scheme; }
    //@}
  };


  /** Storage for isosceles pixel triangles.
   *  Only the unique triangles are stored.  The angular distance between
   *  pixel pairs 2,3 and 3,1 are equal.  The angular distance between pixel
   *  pair 1,2 is different than the other two pairs.
   *
   *  This is a specialized version of Pixel_Triangles.
   */
  template<typename T>
  class Pixel_Triangles_Isosceles : public Pixel_Triangles<T> {
  public :
    /// Generic constructor.
    Pixel_Triangles_Isosceles () : Pixel_Triangles<T>() {}

    /** Find all isosceles triangles.
     *  Find all the triangles that can be made up from the provided two
     *  point tables.  The first table, \a t1, is the one for the two equal
     *  sides.  The first two pixels in the table come from \a t2 and their
     *  order is set by ensuring that the triangle is righthanded.
     */
    void find_triangles (const Twopt_Table<T>& tequal,
			 const Twopt_Table<T>& tother)
    {
      T p1, p2;
      T i2;
      std::vector<T> trip;

      this->initialize (tother, tequal, tequal);

      for (size_t i1=0; i1 < tother.Npix(); ++i1) {
	p1 = tother.pixel_list(i1);
	for (size_t j2=0; (j2 < tother.Nmax()) && (tother(i1,j2) != -1);
	     ++j2) {
	  i2 = tother(i1,j2);
	  p2 = tother.pixel_list(i2);
	  if (p2 < p1) continue; // Don't double count triangles.
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (&tequal(i1,0), &tequal(i1,tequal.Nmax()),
			  &tequal(i2,0), &tequal(i2,tequal.Nmax()), trip);
	  // Now put all the triplets in the list.
	  for (size_t k=0; k < trip.size(); ++k) {
	    add (p1, p2, tequal.pixel_list(trip[k]));
	  }
	}
      }
    }
  };

  /** Storage for equilateral pixel triangles.
   *  Only the unique triangles are stored.  The pixels are stored as in
   *  monotonically increasing order.  This fact can be used to speed up
   *  searches through the triangles.
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
      T p1, p2;
      T i2;
      std::vector<T> trip;

      this->initialize(t, t, t);

      for (size_t i1=0; i1 < t.Npix(); ++i1) {
	p1 = t.pixel_list(i1);
	for (size_t j2=0; (j2 < t.Nmax()) && (t(i1,j2) != -1); ++j2) {
	  i2 = t(i1,j2);
	  p2 = t.pixel_list(i2);
	  if (p2 < p1) continue;
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (i2, &t(i1,0), &t(i1,t.Nmax()),
			  &t(i2,0), &t(i2,t.Nmax()), trip);
	  // Now put all the triplets in the list.
	  for (size_t k=0; k < trip.size(); ++k) {
	    add (p1, p2, t.pixel_list(trip[k]));
	  }
	}
      }
    }
  };

  // Temporary location?
  /** Calculate all quadrilaterals.  This is specialized to Equilateral
   * triangles and only calculates equilateral quadrilaterals.  We use the
   * fact that the pixels in the triangle are stored such that the triangles
   * are righthanded.  We further use the fact that min(pix1, pix2) will
   * be monotonically increasing, thus we can truncate the search.
   *
   * Even with this specialization the quad table can be huge.  For this
   * reason we create a class that incrementally calculates sets of points.
   * This costs more in overhead but requires significantly less memory.
   */
  template<typename T>
  class Quads {
  private :
    size_t ind_curr;
    Pixel_Triangles_Equilateral<T> *t;
    std::vector<T> pts; // So we don't have to keep recreating it.
  public :  
    Quads () : ind_curr(0), t(0), pts(4) {}
    void initialize (Pixel_Triangles_Equilateral<T>& triangle)
    { ind_curr = 0; t = &triangle; }
    bool next (std::vector<std::vector<T> >& quads)
    {
      if (ind_curr == t->size()) return false;
      quads.clear();
      // Order of pts doesn't matter.
      std::copy (t->get(ind_curr).begin(), t->get(ind_curr).end(),
		 pts.begin());
      T cmax = std::max (pts[1], pts[2]);
      for (size_t jj=ind_curr+1; jj < t->size(); ++jj) {
	if (t->get(jj,0) > cmax) break;
	if (t->get(jj,0) == pts[0]) {
	  if (pts[1] == t->get(jj,2)) {
	    pts[3] = t->get(jj,1);
	  } else if (pts[2] == t->get(jj,1)) {
	    pts[3] = t->get(jj,2);
	  } else {
	    continue;
	  }
	} else if (pts[1] == t->get(jj,0)) {
	  if (pts[2] == t->get(jj,2)) {
	    pts[3] = t->get(jj,1);
	  } else {
	    continue;
	  }
	} else if (pts[1] == t->get(jj,1)) {
	  if (pts[2] == t->get(jj,0)) {
	    pts[3] = t->get(jj,2);
	  } else {
	    continue;
	  }
	} else if (pts[1] == t->get(jj,2)) {
	  if (pts[2] == t->get(jj,1)) {
	    pts[3] = t->get(jj,0);
	  } else {
	    continue;
	  }
	} else {
	  continue;
	}
	quads.push_back (pts);
      }
      ++ind_curr;
      return true;
    }
  };


  /** Calculate all rhombic quadrilaterals.  Rhombic quadrilaterals are
   *  constructed from two equilateral triangles connected along one side,
   *  that is, all sides have the same length and one of the diagonals has
   *  the length of the sides.  We use the fact that the pixels in the
   *  triangle are stored in monotonically increasing order.
   *
   *  Even with this specialization the quad table can be huge.  For this
   *  reason we create a class that incrementally calculates sets of points.
   *  This costs more in overhead but requires significantly less memory.
   */
  template<typename T>
  class Quadrilaterals_Rhombic {
  private :
    size_t ind_curr;
    Pixel_Triangles_Equilateral<T> *t;
  public :  
    Quadrilaterals_Rhombic () : ind_curr(0), t(0) {}
    /** Initialize the rhombic quadrilateral search.
     *  The given triangle will be used for subsequent searches.  See
     *  next().
     */
    void initialize (Pixel_Triangles_Equilateral<T>& triangle)
    { 
      ind_curr = 0; t = &triangle;
    }
    /** Get the next set of rhombic quadrilaterals.
     *  The quadrilaterals are constructed for each triangle provided to
     *  initialize(). The quadrilaterals are then made up of the three
     *  points in the triangle, returned in \a pts, and a third
     *  point.  The list of all third points is returned in \a thirdpt.
     *  Note that the orientation of the quadrilateral is lost in this
     *  process.  We do \b not ensure that the quadrilaterals are right
     *  handed.
     */
    bool next (std::vector<T>& pts, std::vector<T>& thirdpt)
    {
      if (ind_curr >= t->size()-1) return false;
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
