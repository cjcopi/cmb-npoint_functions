#include <iostream>
#include <iomanip>
#include <string>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Pixel_Quadrilaterals.h>

#include <healpix_base.h>

#include <pixel_tools.h>
#include <pixel_ringinfo.h>

/* We pull a lot of the details out of Pixel_Quadrilaterals_Rhombic_Full so
 * that we can parallelize the code. */

namespace {
  const std::string CREATE_RHOMBIC_QUADRILATERALS_LIST_PARALLEL_RCSID
  ("$Id: create_rhombic_quadrilaterals_list.cpp,v 1.2 2011-08-10 18:56:28 copi Exp $");
}

// Simple wrapper
struct PixelInfo {
  enum BasePix { BASE0, BASE4} basepix;
  int pixnum;
};
  
// Simple class to handle transforming pixels.
class PixelTrans {
private :
  myHealpix::pixel_ringinfo pri;
  Healpix_Base HBase;

  // Wrappers to handle the scheme for our pixels.
  inline void pri_setpix (int p)
  {
    if (HBase.Scheme() == NEST) {
      pri.from_pixel (HBase.nest2ring(p));
    } else {
      pri.from_pixel (p);
    }
  }

  inline int pri_frompix ()
  {
    if (HBase.Scheme() == NEST) {
      return HBase.ring2nest(pri.to_pixel());
    } else {
      return pri.to_pixel();
    }
  }

  inline int shift_pix_by_base (int p)
  {
    pri_setpix (p);
    pri.shift_by_base_pixel();
    return pri_frompix();
  }      

  inline int reflect_pix_through_zaxis (int p)
  {
    pri_setpix (p);
    pri.reflect_through_zaxis();
    return pri_frompix();
  }

  inline int reflect_pix_through_z0 (int p)
  {
    pri_setpix (p);
    pri.reflect_through_z0();
    return pri_frompix();
  }

public :
  PixelTrans (size_t Nside, Healpix_Ordering_Scheme scheme)
    : pri(Nside), HBase (Nside, scheme, SET_NSIDE) {}

  inline void shift_by_base (std::vector<int>& pl)
  { for (size_t j=0; j < pl.size(); ++j) pl[j] = shift_pix_by_base (pl[j]); }

  inline void reflect_through_zaxis (std::vector<int>& pl)
  { for (size_t j=0; j < pl.size(); ++j) pl[j] = reflect_pix_through_zaxis (pl[j]); }

  inline void reflect_through_z0 (std::vector<int>& pl)
  { for (size_t j=0; j < pl.size(); ++j) pl[j] = reflect_pix_through_z0 (pl[j]); }
};


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <two point table name>\n";
  exit (0);
}

void write_quad_buffer (const std::vector<int>& quad_buf)
{
#pragma omp critical
  {
    for (size_t j=0; j < quad_buf.size(); j+=4) {
      for (int k=0; k < 4; ++k) {
	std::cout << quad_buf[j+k] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void add_quads (const std::vector<int>& tri,
		const std::vector<int>& thirdpt,
		std::vector<int>& quad_buf, 
		size_t Nbuf)
{
  for (size_t j=0; j < thirdpt.size(); ++j) {
    for (size_t i=0; i < tri.size(); ++i) quad_buf.push_back(tri[i]);
    quad_buf.push_back(thirdpt[j]);
    /* We really don't want the buffer to "overflow" since
     * it could end up needed a huge amount of memory.  Instead we
     * will just empty it if the need arises. */
    if (quad_buf.size() >= 4*Nbuf) {
      write_quad_buffer (quad_buf);
      quad_buf.clear();
    }
  }
}


int main (int argc, char *argv[])
{
  if (argc != 2) usage (argv[0]);

  std::string twopt_table_file = argv[1];

  Npoint_Functions::Pixel_Quadrilaterals_Rhombic<int> q;
  Npoint_Functions::Twopt_Table<int> twopt_table;
  twopt_table.read_file (twopt_table_file);
  Npoint_Functions::Pixel_Triangles_Equilateral<int> triangles;
  triangles.find_triangles (twopt_table);
  q.initialize (triangles);

  /* Build the list of pixels storing information about their base pixel as
   * this is needed for the transformations. */
  std::vector<PixelInfo> pixel_list;
  {
    std::vector<int> pl0, pl4;
    myHealpix::base0_list (q.Nside(), pl0);
    myHealpix::base4_list (q.Nside(), pl4);
    pixel_list.resize(pl0.size()+pl4.size());
    for (size_t j=0; j < pl0.size(); ++j) {
      pixel_list[j].pixnum = pl0[j];
      pixel_list[j].basepix = PixelInfo::BASE0;
    }      
    for (size_t j=0; j < pl4.size(); ++j) {
      pixel_list[j+pl0.size()].pixnum = pl4[j];
      pixel_list[j+pl0.size()].basepix = PixelInfo::BASE4;
    }      
  }

#pragma omp parallel shared (pixel_list) firstprivate (q)
  {
    std::vector<int> tri;
    std::vector<int> thirdpt;
    thirdpt.reserve(1000);
    int pix;
    
    PixelTrans pixtrans (q.Nside(), q.Scheme());

    /* Buffer space.  We save the quadrilaterals in a buffer and then
     * write them out all at once.  This is done so that we don't have
     * threads fighting each other to get write access.  We can't have them
     * all write at once.
     */
    // Number of quad space we want to reserve.  This is a bit under 500MB.
    const size_t Nbuf = 30000000;
    std::vector<int> quad_buf;
    quad_buf.reserve(4*Nbuf);
    
#pragma omp for schedule(dynamic,1)
    for (size_t j=0; j < pixel_list.size(); ++j) {
      pix = pixel_list[j].pixnum; // shorthand
      q.initialize(pix);
      while (q.next(tri, thirdpt)) {
	// First add the quads.
	add_quads (tri, thirdpt, quad_buf, Nbuf);
	// Next shift by base pixel 3 times.
	for (int n=0; n < 3; ++n) {
	  pixtrans.shift_by_base (tri);
	  pixtrans.shift_by_base (thirdpt);
	  add_quads (tri, thirdpt, quad_buf, Nbuf);
	}
	// Then reflect through z=0 line
	pixtrans.reflect_through_z0 (tri);
	pixtrans.reflect_through_z0 (thirdpt);
	add_quads (tri, thirdpt, quad_buf, Nbuf);
	// and shift by base pixel 3 times.
	for (int n=0; n < 3; ++n) {
	  pixtrans.shift_by_base (tri);
	  pixtrans.shift_by_base (thirdpt);
	  add_quads (tri, thirdpt, quad_buf, Nbuf);
	}

	// If we have a base0 pixel we are done
	if (pixel_list[j].basepix == PixelInfo::BASE0) continue;

	// Otherwise we have more transformations to do.
	// Reflect through z-axis
	pixtrans.reflect_through_zaxis (tri);
	pixtrans.reflect_through_zaxis (thirdpt);
	add_quads (tri, thirdpt, quad_buf, Nbuf);
	// and shift by base pixel 3 times.
	for (int n=0; n < 3; ++n) {
	  pixtrans.shift_by_base (tri);
	  pixtrans.shift_by_base (thirdpt);
	  add_quads (tri, thirdpt, quad_buf, Nbuf);
	}
	// Then reflect through z=0 line
	pixtrans.reflect_through_z0 (tri);
	pixtrans.reflect_through_z0 (thirdpt);
	add_quads (tri, thirdpt, quad_buf, Nbuf);
	// and shift by base pixel 3 times.
	for (int n=0; n < 3; ++n) {
	  pixtrans.shift_by_base (tri);
	  pixtrans.shift_by_base (thirdpt);
	  add_quads (tri, thirdpt, quad_buf, Nbuf);
	}
      }
      write_quad_buffer (quad_buf);
      quad_buf.clear();
    }
  }

  return 0;
}
