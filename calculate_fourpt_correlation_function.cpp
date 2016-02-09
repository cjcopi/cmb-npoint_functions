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
  ("$Id$");
}


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <map fits file> "
            << "<quad list prefix> [<mask file>]\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if ((argc < 3) || (argc > 4)) usage (argv[0]);
  std::string mapfile = argv[1];
  std::string quad_list_prefix = argv[2];

  Healpix_Map<double> map;
  read_Healpix_map_from_fits (mapfile, map);
  bool have_mask = false;
  Healpix_Map<double> mask;
  if (argc == 4) {
    read_Healpix_map_from_fits (argv[3], mask);
    have_mask = true;
  }

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
      if (have_mask) {
	if (static_cast<size_t>(mask.Nside()) != qlf.Nside()) {
	  std::cerr << "Mask and quadrilateral lists do not have"
		    << " the same Nside: " << mask.Nside() 
		    << " != " << qlf.Nside() << std::endl;
	  std::exit(1);
	}
	if (mask.Scheme() != qlf.Scheme()) mask.swap_scheme();
      }

#pragma omp critical
      {
	std::cerr 
#ifdef OMP
	  << omp_get_thread_num() << " "
#endif       
	  << k << std::endl;
      }

      bin_list[k] = qlf.bin_value();
      if (have_mask) {
	Corr[k] = calculate_masked_fourpoint_function (map, mask, qlf);
      } else {
	Corr[k] = calculate_fourpoint_function (map, qlf);
      }
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
