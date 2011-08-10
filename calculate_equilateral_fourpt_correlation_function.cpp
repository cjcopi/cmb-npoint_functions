#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

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
  ("$Id: calculate_equilateral_fourpt_correlation_function.cpp,v 1.4 2011-07-27 04:16:53 copi Exp $");
}


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
    Npoint_Functions::Twopt_Table<int> twopt_table;
    Npoint_Functions::Pixel_Triangles_Equilateral<int> triangles;
    std::vector<std::vector<int> > quads;
    Npoint_Functions::Quads<int> q;

#pragma omp for schedule(dynamic,2)
    for (size_t k=0; k < twopt_table_file.size(); ++k) {
      twopt_table.read_file (twopt_table_file[k]);

      /* So only one thread writes at a time.  Really only matters on
       * initial startup. */
#pragma omp critical
      {
	std::cerr 
#ifdef OMP
	  << omp_get_thread_num() << " "
#endif       
	  << k << " " << twopt_table.Nmax() << std::endl;
      }

      bin_list[k] = twopt_table.bin_value();
      if (twopt_table.Nmax() == 0) {
	std::cerr << "Nmax = 0 for bin " << k << std::endl;
	continue;
      }
      C4 = 0;
      Nquads = 0;
      triangles.find_triangles (twopt_table);
      q.initialize (triangles);
      while (q.next(quads)) {
	for (size_t j=0; j < quads.size(); ++j) {
	  C4 += map[quads[j][0]] * map[quads[j][1]]
	    * map[quads[j][2]] * map[quads[j][3]];
	}
	Nquads += quads.size();
      }
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
