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
  ("$Id: create_twopt_table.cpp,v 1.1 2011-07-07 02:41:27 copi Exp $");
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
  
  std::ostringstream sstr;
  int Nbin = 0;
  size_t Npair = 0;
  double C2, Csum;
  Twopt_Table<int> twopt_table;
  while (true) {
    sstr.str("");
    sstr << twopt_prefix << std::setw(5) << std::setfill('0') << Nbin << ".dat";
    if (! twopt_table.read_file (sstr.str())) break;
    C2 = 0;
    for (size_t i=0; i < twopt_table.Npix(); ++i) {
      Csum = 0;
      for (size_t j=0;
	   ((j < twopt_table.Nmax())
	    && (twopt_table.pixel_list()[j] != -1));
	   ++j) {
	if (twopt_table.pixel_list()[i] > twopt_table.pixel_list()[j])
	  continue;
	++Npair;
	Csum += map[twopt_table.pixel_list()[j]];
      }
      C2 += map[twopt_table.pixel_list()[i]] * Csum;
    }
    C2 /= Npair;
    std::cout << twopt_table.bin_value() << " " << C2 << std::endl;
    ++Nbin;
  }

  return 0;
}
