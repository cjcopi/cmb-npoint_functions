#include <iostream>
#include <string>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CALCULATE_EQUILATERAL_THREEPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_equilateral_threept_correlation_function.cpp,v 1.4 2011-07-14 17:35:45 copi Exp $");
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

  std::vector<std::string> twopt_table_list
    = Npoint_Functions::get_sequential_file_list (twopt_prefix);
  std::vector<double> bin_list(twopt_table_list.size());
  std::vector<double> Corr(twopt_table_list.size());

#pragma omp parallel shared(twopt_table_list, Corr, bin_list)
  {
    double C3;
    Npoint_Functions::Twopt_Table<int> twopt_table;
    Npoint_Functions::Pixel_Triangles_Equilateral<int> triangles;

#pragma omp for schedule(guided)
    for (size_t k=0; k < twopt_table_list.size(); ++k) {
      twopt_table.read_file (twopt_table_list[k]);
      // Cast to quiet the compiler about the signed/unsigned comparison.
      if ((size_t)map.Npix() < twopt_table.Npix()) {
      	std::cerr << "Map does not have enough pixels.\n";
	std::exit(1);
      }
      triangles.find_triangles (twopt_table);
      C3 = 0;
      for (size_t j=0; j < triangles.size(); ++j) {
	C3 += map[triangles(j)[0]] * map[triangles(j)[1]]
	  * map[triangles(j)[2]];
      }
      if (triangles.size() != 0) C3 /= triangles.size();
      bin_list[k] = triangles.lengths()[0];
      Corr[k] = C3;
    }
  }

  for (size_t k=0; k < bin_list.size(); ++k) {
    // Same format as spice
    std::cout << std::acos(bin_list[k]) << " " << bin_list[k] << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
