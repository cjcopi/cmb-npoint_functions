#ifndef PIXEL_TRIANGLES_H
#define PIXEL_TRIANGLES_H

#include <vector>
#include <string>

#include <healpix_base.h>
#include <vec3.h>

namespace {
  /// @cond IDTAG
  const std::string PIXEL_TRIANGLES_RCSID
  ("$Id: Pixel_Triangles.h,v 1.7 2011-07-23 01:00:31 copi Exp $");
  /// @endcond
}

namespace {
  /** Helper function to create a list of vectors pointing to HEALPix pixel
   * centers. 
   */
  template<typename T>
  void fill_vector_list (const Npoint_Functions::Twopt_Table<T>& t,
			 std::vector<vec3>& veclist)
  {
    Healpix_Base HBase (t.Nside(), NEST, SET_NSIDE);
    veclist.resize (t.Npix());
    for (size_t i=0; i < t.Npix(); ++i) {
      veclist[i] = HBase.pix2vec (t.pixel_list()[i]);
    }
  }
}

namespace Npoint_Functions {
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
    std::vector<int> orient; // Orientation of the triangle.

    /** Calculate the orientation from three vectors.
     *  The orientation is an integer:  +1 for righthanded and -1 for lefthanded.
     *  Righthanded is defined by 
     *  \f[ (\hat n_1\times\hat n_2)\cdot \hat n_3 > 0. \f]
     */
    int calculate_orientation (const vec3& n1, const vec3& n2,
			       const vec3& n3)
    {
      double val = dotprod (crossprod (n1, n2), n3);
      return ((val > 0) ? +1 : -1);
    }
  
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
	p1 = t1.pixel_list()[i1];
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
	    orient.push_back (calculate_orientation (v[i1], v[i2],
						     trip[k])); 
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
    /** The orientation of the triangle.
     *  A value of +1 represents righthanded triangles and -1 for
     *  lefthanded triangles.  See calculate_orientation() for more
     *  details. 
     */
    inline int orientation (size_t j) const
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
    void find_triangles (const Twopt_Table<T>& t1,
			 const Twopt_Table<T>& t2)
    {
      T p1, p2, p3;
      T i1, i2;
      std::vector<T> trip;
      this->reset();
      set_edge_lengths (t2.bin_value(), t1.bin_value(), t1.bin_value());

      // Create list of vectors so we can get the orientation.
      std::vector<vec3> v;
      fill_vector_list (t1, v);

      for (size_t j1=0; j1 < t2.Npix(); ++j1) {
	i1 = j1; // to make the code look symmetric
	p1 = t2.pixel_list()[i1];
	for (size_t j2=j1+1; (j2 < t2.Nmax()) && (t2(j1,j2) != -1); ++j2) {
	  i2 = t2(j1,j2);
	  p2 = t2.pixel_list()[i2];
	  // Finally can search for and add appropriate pairs.
	  trip.clear();
	  this->append_matches (&t1(i1,0), t1.Nmax(), &t1(i2,0), t1.Nmax(),
				trip);
	  // Now put all the triplets in the list
	  for (size_t k=0; k < trip.size(); ++k) {
	    p3 = t1.pixel_list()[trip[k]];
	    if (calculate_orientation (v[p1], v[p2], v[p3]) > 0)
	      add (p1, p2, p3);
	    else
	      add (p2, p1, p3);
	  }
	}
      }
      // All triangles are righthanded.
      this->orient.assign (this->size(), 1);
    }
  };

  /** Storage for equilateral pixel triangles.
   *  Only the unique triangles are stored.  The pixels are stored as to
   *  ensure that the triangles are righthanded.
   *
   *  This is a specialized version of Pixel_Triangles_Isosceles.
   */
  template<typename T>
  class Pixel_Triangles_Equilateral : public Pixel_Triangles_Isosceles<T> {
  protected :
    /** Find matches in two lists and append them to a new list.
     *  Here the minimum allowed value is provided.  All values appended to the
     *  list will be greater than or equal to this value. */
    void append_matches (T minval, const T *L1, size_t NL1,
                         const T *L2, size_t NL2,
                         std::vector<T>& res)
    {
      const T *it1, *it2; // "iterators"
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
	p1 = t.pixel_list()[i1];
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
	    if (calculate_orientation (v[i1], v[i2], v[trip[k]]) > 0)
	      add (p1, p2, t.pixel_list()[trip[k]]);
	    else
	      add (p2, p1, t.pixel_list()[trip[k]]);
	  }
	}
      }
      // All triangles are righthanded.
      this->orient.assign (this->size(), 1);
    }
  };
}
#endif

/* For emacs, this is a c++ header
 * Local Variables:
 * mode: c++
 * End:
 */
