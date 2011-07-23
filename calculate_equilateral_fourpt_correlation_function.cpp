#include <iostream>
#include <iomanip>
#include <string>

#ifdef OMP
#include <omp.h>
#endif

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CALCULATE_QUADRILATERAL_FOURPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_isosceles_threept_correlation_function.cpp,v 1.3 2011-07-16 23:52:38 copi Exp $");
}

/* Calculate all quadrilaterals.  This is specialized to Isosceles
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
  Npoint_Functions::Pixel_Triangles_Isosceles<T> *t;
  std::vector<T> pts; // So we don't have to keep recreating it.
public :  
  Quads () : ind_curr(0), t(0), pts(4) {}
  void initialize (Npoint_Functions::Pixel_Triangles_Isosceles<T>&
		   triangle)
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


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <map fits file> "
            << "<twopt tables prefix>\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if (argc != 3) usage (argv[0]);
  std::string mapfile = argv[1];
  std::string twopt_prefix = argv[2];

  Healpix_Map<double> map;
  read_Healpix_map_from_fits (mapfile, map);
  if (map.Scheme() == RING) map.swap_scheme();

  // Figure out how many bins there are by trying to open files.
  std::vector<std::string> twopt_table_file
    = Npoint_Functions::get_sequential_file_list(twopt_prefix);

  std::vector<double> bin_list(twopt_table_file.size());
  std::vector<double> Corr(twopt_table_file.size());

#pragma omp parallel shared(Corr, bin_list, twopt_table_file)
  {
    double C4;
    size_t Nquads;
    Npoint_Functions::Twopt_Table<int> twopt_table_equal;
    Npoint_Functions::Twopt_Table<int> twopt_table;
    Npoint_Functions::Pixel_Triangles_Isosceles<int> triangles;
    std::vector<std::vector<int> > quads;
    Quads<int> q;

#pragma omp for schedule(dynamic,2)
    for (size_t k=0; k < twopt_table_file.size(); ++k) {
      /* So only one thread writes at a time.  Really only matters on
       * initial startup. */
#pragma omp critical
      {
	std::cerr 
#ifdef OMP
	  << omp_get_thread_num() << " "
#endif       
	  << k << std::endl;
      }

      twopt_table_equal.read_file (twopt_table_file[k]);
      C4 = 0;
      Nquads = 0;
      // Now loop over the bins again so we get the unequal bin length.
      for (size_t kk=0; kk < twopt_table_file.size(); ++kk) {
	twopt_table.read_file (twopt_table_file[kk]);
	triangles.find_triangles (twopt_table_equal, twopt_table);
	q.initialize (triangles);
	while (q.next(quads)) {
	  for (size_t j=0; j < quads.size(); ++j) {
	    C4 += map[quads[j][0]] * map[quads[j][1]]
	      * map[quads[j][2]] * map[quads[j][3]];
	  }
	  Nquads += quads.size();
	}
      }
      bin_list[k] = triangles.lengths()[1];
      if (Nquads != 0) C4 /= Nquads;
      Corr[k] = C4;
    }
  }

  for (size_t k=0; k < bin_list.size(); ++k) {
    // Same format as spice
    std::cout << std::acos(bin_list[k]) << " " << bin_list[k] << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
