#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>
//#include <paramfile.h>

#include <Twopt_Table.h>

namespace {
  const std::string CALCULATE_TWOPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_twopt_correlation_function.cpp,v 1.2 2011-07-07 18:30:04 copi Exp $");
}


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <map fits file> <twopt tables prefix>\n";
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
      sstr << twopt_prefix << std::setw(5) << std::setfill('0') << Nbin << ".dat";
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
    size_t Npair;
    double C2, Csum;
    int p1, p2;
    Twopt_Table<int> twopt_table;
#pragma omp for schedule(guided,Nbin/10)
    for (int k=0; k < Nbin; ++k) {
      std::cerr << k << std::endl;
      sstr.str("");
      sstr << twopt_prefix << std::setw(5) << std::setfill('0') << k << ".dat";
      twopt_table.read_file (sstr.str());
      C2 = 0;
      Npair = 0;
      for (size_t i=0; i < twopt_table.Npix(); ++i) {
	Csum = 0;
	p1 = twopt_table.pixel_list()[i];
	for (size_t j=0;
	     ((j < twopt_table.Nmax())
	      && (twopt_table.get_read_value(i,j) != -1));
	     ++j) {
	  p2 = twopt_table.pixel_list()[twopt_table.get_read_value(i,j)];
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

  for (int k=0; k < Nbin; ++k) {
    // Same format as spice
    std::cout << std::acos(bin_list[k]) << " " << bin_list[k] << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
