#include <iostream>
#include <string>
#include <cmath>

#ifdef OMP
#include <omp.h>
#endif

#include <healpix_map.h>
#include <healpix_map_fitsio.h>

#include <Quadrilateral_List_File.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CALCULATE_FOURPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_fourpt_correlation_function.cpp,v 1.2 2011-08-15 15:40:36 copi Exp $");
}


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <map fits file> "
            << "<quad list prefix>\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if (argc != 3) usage (argv[0]);
  std::string mapfile = argv[1];
  std::string quad_list_prefix = argv[2];

  Healpix_Map<double> map;
  read_Healpix_map_from_fits (mapfile, map);


  // Figure out how many bins there are by trying to open files.
  std::vector<std::string> quad_list_files
    = Npoint_Functions::get_range_file_list(quad_list_prefix, 0, 180);

  std::vector<double> bin_list(quad_list_files.size());
  std::vector<double> Corr(quad_list_files.size());

#pragma omp parallel shared(Corr, bin_list, quad_list_files)
  {
    Npoint_Functions::Quadrilateral_List_File<int> qlf;

#pragma omp for schedule(dynamic,2)
    for (size_t k=0; k < quad_list_files.size(); ++k) {

      if (! qlf.initialize (quad_list_files[k])) {
	std::cerr << "Error initializing quadrilateral list from "
		  << quad_list_files[k] << std::endl;
	std::exit(1);
      }
      if (static_cast<size_t>(map.Nside()) != qlf.Nside()) {
	std::cerr << "Map has Nside = " << map.Nside()
		  << " but quad list has Nside = " << qlf.Nside()
		  << "\nGiving up!\n";
	std::exit(1);
      }
      if (map.Scheme() != qlf.Scheme()) map.swap_scheme();

#pragma omp critical
      {
	std::cerr 
#ifdef OMP
	  << omp_get_thread_num() << " "
#endif       
	  << k << std::endl;
      }

      bin_list[k] = qlf.bin_value();
      Corr[k] = calculate_fourpoint_function (map, qlf);
    }
  }
  
  for (size_t k=0; k < bin_list.size(); ++k) {
    // Same format as spice
    std::cout << bin_list[k]*M_PI/180 << " " 
	      << cos(bin_list[k]*M_PI/180) << " "
	      << Corr[k] << std::endl;
  }

  return 0;
}
