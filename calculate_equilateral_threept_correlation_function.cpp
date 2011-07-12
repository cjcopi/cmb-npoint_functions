#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>

namespace {
  const std::string CALCULATE_EQUILATERAL_THREEPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_twopt_correlation_function.cpp,v 1.6 2011-07-10 02:43:47 copi Exp $");
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
  int Nbin = 0;
  {
    std::ostringstream sstr;
    std::ifstream in;
    while (true) {
      sstr.str("");
      sstr << twopt_prefix << std::setw(5) << std::setfill('0') << Nbin
           << ".dat";
      ++Nbin;
      in.open (sstr.str().c_str());
      if (! in) break;
      in.close();
    }
  }

  std::vector<double> bin_list(Nbin);
  std::vector<double> Corr(Nbin);

  std::cerr << "Nbin = " << Nbin << std::endl;
#pragma omp parallel shared(Nbin, Corr, bin_list)
  {
    std::ostringstream sstr;
    double C3;
    Twopt_Table<int> twopt_table;
    Pixel_Triangles_Equilateral<int> triangles;

#pragma omp for schedule(guided)
    for (int k=0; k < Nbin; ++k) {
      sstr.str("");
      sstr << twopt_prefix << std::setw(5) << std::setfill('0') << k << ".dat";
      twopt_table.read_file (sstr.str());
      triangles.find_triangles (twopt_table);
      C3 = 0;
      for (size_t j=0; j < triangles.size(); ++j) {
	C3 += map[triangles(j)[0]] * map[triangles(j)[1]]
	  * map[triangles(j)[2]] ;
      }
      C3 /= triangles.size();
      bin_list[k] = twopt_table.bin_value();
      Corr[k] = C3;
    }
  }

  for (int k=0; k < Nbin; ++k) {
    // Same format as spice
    std::cout << std::acos(bin_list[k]) << " " << bin_list[k] << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
