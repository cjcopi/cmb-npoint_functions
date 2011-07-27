#ifndef PIXEL_TRIANGLES_H
#define PIXEL_TRIANGLES_H

#include <vector>
#include <string>

#include <healpix_base.h>
#include <vec3.h>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_TRIANGLES_RCSID
  ("$Id: Pixel_Triangles.h,v 1.16 2011-07-25 22:13:15 copi Exp $");
  /// @endcond
}

namespace {
  /* Helper function to create a list of vectors pointing to HEALPix pixel
   *  centers.  The vectors are labelled by the pixel INDEX in the
   *  two point table, not the actual pixel number.
   */
  template<typename T>
  void fill_vector_list (const Npoint_Functions::Twopt_Table<T>& t,
			 std::vector<vec3>& veclist)
  {
    Healpix_Base HBase (t.Nside(), NEST, SET_NSIDE);
    veclist.resize (t.Npix());
    for (size_t i=0; i < t.Npix(); ++i) {
      veclist[i] = HBase.pix2vec (t.pixel_list(i));
    }
  }

    /** Find matches in two lists and append them to a new list. */
  template<typename T>
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

  /* Find matches in two lists and append them to a new list.
   *  Here the minimum allowed value is provided.  All values appended to the
   *  list will be greater than or equal to this value. */
  template<typename T>
  void append_matches (T minval, const T *L1, size_t NL1,
		       const T *L2, size_t NL2,
		       std::vector<T>& res)
  {
    const T *it1, *it2; // "iterators"
    it1 = L1;
    while ((it1 != &L1[NL1]) && (*it1 < minval)) ++it1;
    it2 = L2;
    while ((it2 != &L2[NL2]) && (*it2 < minval)) ++it2;
    append_matches (it1, NL1, it2, NL2, res);
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
  protected :
     /// Orientation of the triangles.
    std::vector<Orientation> orient;

    /// Add a triangle to the list.
    inline void add (const T& p1, const T& p2, const T& p3)
    {
      size_t n = triangles.size();
      triangles.push_back (std::vector<T>(3));
      triangles[n][0] = p1;
      triangles[n][1] = p2;
      triangles[n][2] = p3;
    }

    /// Set the edge lengths of the triangle
    inline void set_edge_lengths (double l1, double l2, double l3)
    {
      this->edge_length[0] = l1;
      this->edge_length[1] = l2;
      this->edge_length[2] = l3;
    }

  public :
    /// Generic constructor.
    Pixel_Triangles () : triangles(), edge_length(3), orient() {}

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
      T p1, p2;
      T i1, i2;
      std::vector<T> trip;
      this->reset();
      set_edge_lengths (t1.bin_value(), t2.bin_value(), t3.bin_value());


      // Create list of vectors so we can get the orientation.
      std::vector<vec3> v;
      fill_vector_list (t1, v);
      orient.clear();

      for (size_t j1=0; j1 < t1.Npix(); ++j1) {
	i1 = j1; // to make the code look symmetric
	p1 = t1.pixel_list(i1);
	for (size_t j2=0; (j2 < t1.Nmax()) && (t1(j1,j2) != -1); ++j2) {
	  i2 = t1(j1,j2);
	  p2 = t1.pixel_list(i2);
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (&t2(i1,0), t2.Nmax(), &t3(i2,0), t3.Nmax(),
			  trip);
	  // Now put all the triplets in the list
	  for (size_t k=0; k < trip.size(); ++k) {
	    add (p1, p2, t1.pixel_list(trip[k]));
	    orient.push_back (calculate_orientation (v[i1], v[i2],
						     v[trip[k]])); 
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
    //@}
  };


  /** Storage for isosceles pixel triangles.
   *  Only the unique triangles are stored.  The angular distance between
   *  pixel pairs 2,3 and 3,1 are equal.  The angular distance between pixel
   *  pair 1,2 is different than the other two pairs.  The order of pixels
   *  1 and 2 is set to ensure that the orientation of the triangle is
   *  righthanded. 
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
      T i1, i2;
      std::vector<T> trip;
      this->reset();
      set_edge_lengths (tother.bin_value(), tequal.bin_value(),
			tequal.bin_value());

      // Create list of vectors so we can get the orientation.
      std::vector<vec3> v;
      fill_vector_list (tequal, v);

      for (size_t j1=0; j1 < tother.Npix(); ++j1) {
	i1 = j1; // to make the code look symmetric
	p1 = tother.pixel_list(i1);
	for (size_t j2=0; (j2 < tother.Nmax()) && (tother(j1,j2) != -1);
	     ++j2) {
	  i2 = tother(j1,j2);
	  p2 = tother.pixel_list(i2);
	  if (p2 < p1) continue; // Don't double count triangles.
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (&tequal(i1,0), tequal.Nmax(),
			  &tequal(i2,0), tequal.Nmax(), trip);
	  /* Now put all the triplets in the list making sure the triangles
	   * are righthanded.  The first two pixels are swapped if the
	   * given orientation is lefthanded.
	   */
	  for (size_t k=0; k < trip.size(); ++k) {
	    if (calculate_orientation (v[i1], v[i2], v[trip[k]]) 
		== RIGHTHANDED)
	      add (p1, p2, tequal.pixel_list(trip[k]));
	    else
	      add (p2, p1, tequal.pixel_list(trip[k]));
	  }
	}
      }
      // All triangles are righthanded.
      this->orient.assign (this->size(), RIGHTHANDED);
    }
  };

  /** Storage for equilateral pixel triangles.
   *  Only the unique triangles are stored.  The pixels are stored as to
   *  ensure that the triangles are righthanded and that the first pixel
   *  has the smallest index.
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
      T i1, i2;
      std::vector<T> trip;
      this->reset();
      set_edge_lengths (t.bin_value(), t.bin_value(), t.bin_value());

      // Create list of vectors so we can get the orientation.
      std::vector<vec3> v;
      fill_vector_list (t, v);

      for (size_t j1=0; j1 < t.Npix(); ++j1) {
	i1 = j1; // to make the code look symmetric
	p1 = t.pixel_list(i1);
	for (size_t j2=0; (j2 < t.Nmax()) && (t(j1,j2) != -1); ++j2) {
	  i2 = t(j1,j2);
	  p2 = t.pixel_list(i2);
	  if (p2 < p1) continue;
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  append_matches (i2, &t(i1,0), t.Nmax(), &t(i2,0), t.Nmax(),
			  trip);
	  /* Now put all the triplets in the list making sure the triangles
	   * are righthanded.  The last two pixels are swapped if the
	   * given orientation is lefthanded.
	   */
	  for (size_t k=0; k < trip.size(); ++k) {
	    if (calculate_orientation (v[i1], v[i2], v[trip[k]]) 
		== RIGHTHANDED)
	      add (p1, p2, t.pixel_list(trip[k]));
	    else
	      add (p1, t.pixel_list(trip[k]), p2);
	  }
	}
      }
      // All triangles are righthanded.
      this->orient.assign (this->size(), RIGHTHANDED);
    }
  };

  // Temporary location?
  /** Calculate all quadrilaterals.  This is specialized to Isosceles
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
    Pixel_Triangles_Isosceles<T> *t;
    std::vector<T> pts; // So we don't have to keep recreating it.
  public :  
    Quads () : ind_curr(0), t(0), pts(4) {}
    void initialize (Pixel_Triangles_Isosceles<T>& triangle)
    { ind_curr = 0; t = &triangle; }
    bool next (std::vector<std::vector<T> >& quads)
    {
      if (ind_curr == t->size()) return false;
      quads.clear();
      T min1, min2;
      // Order of pts doesn't matter.
      std::copy ((*t)(ind_curr).begin(), (*t)(ind_curr).end(), pts.begin());
      min1 = std::min ((*t)(ind_curr)[0], (*t)(ind_curr)[1]);
      for (size_t jj=ind_curr+1; jj < t->size(); ++jj) {
	min2 = std::min ((*t)(jj)[0], (*t)(jj)[1]);
	if (min2 > min1) break;
	if (((*t)(ind_curr)[0] == (*t)(jj)[1])
	    && ((*t)(ind_curr)[1] == (*t)(jj)[0])) {
	  pts[3] = (*t)(jj)[2];
	  quads.push_back (pts);
	}
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
