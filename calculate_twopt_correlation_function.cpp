#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <Twopt_Table.h>

namespace {
  const std::string CALCULATE_TWOPT_CORRELATION_FUNCTION_RCSID
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
  std::vector<std::string> twopt_table_file;
  {
    int Nbin = 0;
    std::ostringstream sstr;
    Npoint_Functions::Twopt_Table<int> tp;
    while (true) {
      sstr.str("");
      sstr << twopt_prefix << std::setw(5) << std::setfill('0') << Nbin
           << ".dat";
      if (! tp.read_file_header (sstr.str())) break;
      twopt_table_file.push_back (sstr.str());
      ++Nbin;
    }
  }

  std::vector<double> bin_list(twopt_table_file.size());
  std::vector<double> Corr(twopt_table_file.size());

#pragma omp parallel shared(Corr, bin_list, twopt_table_file)
  {
    size_t Npair;
    double C2, Csum;
    int p1, p2;
    Npoint_Functions::Twopt_Table<int> twopt_table;
#pragma omp for schedule(guided)
    for (size_t k=0; k < twopt_table_file.size(); ++k) {
      twopt_table.read_file (twopt_table_file[k]);
      C2 = 0;
      Npair = 0;
      for (size_t i=0; i < twopt_table.Npix(); ++i) {
	Csum = 0;
	p1 = twopt_table.pixel_list()[i];
	for (size_t j=0;
	     ((j < twopt_table.Nmax())
	      && (twopt_table(i,j) != -1));
	     ++j) {
	  p2 = twopt_table.pixel_list()[twopt_table(i,j)];
	  if (p1 > p2) continue; // Avoid double counting.
	  ++Npair;
	  Csum += map[p2];
	}
	C2 += map[p1] * Csum;
      }
      C2 /= Npair;
      bin_list[k] = twopt_table.bin_value();
      Corr[k] = C2;
    }
  }

  for (size_t k=0; k < twopt_table_file.size(); ++k) {
    // Same format as spice
    std::cout << std::acos(bin_list[k]) << " " << bin_list[k] << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
