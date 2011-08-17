#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#ifdef OMP
#include <omp.h>
#endif

#include <healpix_map.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <alm_powspec_tools.h>
#include <powspec.h>
#include <powspec_fitsio.h>
#include <planck_rng.h>

#include <Quadrilateral_List_File.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CALCULATE_LCDM_FOURPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_LCDM_fourpt_correlation_function.cpp,v 1.5 2011-08-16 02:32:32 copi Exp $");
}


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <cl fits file> "
            << "<quad list prefix> <num maps to generate>\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if (argc != 4) usage (argv[0]);
  std::string clfile = argv[1];
  std::string quad_list_prefix = argv[2];
  size_t Nmaps;
  if (! Npoint_Functions::from_string (argv[3], Nmaps)) {
    std::cerr << "Could not parse Nmaps\n";
    usage (argv[0]);
  }

  // Figure out how many bins there are by trying to open files.
  std::vector<std::string> quad_list_files
    = Npoint_Functions::get_range_file_list(quad_list_prefix, 0, 400);
  if (quad_list_files.size() == 0) {
    std::cerr << "No quad list files found!\n";
    usage (argv[0]);
  }

  int Lmax;
  Healpix_Ordering_Scheme qlf_scheme;
  std::vector<Healpix_Map<double> > maps (Nmaps);
  {
    // Figure out Lmax from the Nside and set up maps.
    Npoint_Functions::Quadrilateral_List_File<int> qlf;
    qlf.initialize (quad_list_files[0]);
    qlf_scheme = qlf.Scheme();
    Lmax = std::min (2000UL, 4*qlf.Nside()+1);
    for (size_t j=0; j < maps.size(); ++j) {
      // alm2map REQUIRES the map to be in RING order.
      maps[j].SetNside (qlf.Nside(), RING);
    }
  }  
  PowSpec cl;
  read_powspec_from_fits (clfile, cl, 1, Lmax);

  // Make the maps
#pragma omp parallel shared(cl, maps)
  {
    planck_rng rng;
    /* Seed with random values.  Make sure the threads don't stomp on each
     * other by making the seeding section critical. */
#pragma omp critical
    {
      unsigned int seed[4];
      std::ifstream inseed ("/dev/urandom",
			    std::fstream::in | std::fstream::binary);
      inseed.read (reinterpret_cast<char*>(seed), sizeof(seed));
      inseed.close();
      rng.seed (seed[0], seed[1], seed[2], seed[3]);
    }
    Alm<xcomplex<double> > alm (cl.Lmax(), cl.Lmax());
#pragma omp for schedule(static)
    for (size_t k=0; k < maps.size(); ++k) {
      create_alm (cl, alm, rng);
      alm2map (alm, maps[k]);
      if (maps[k].Scheme() != qlf_scheme) maps[k].swap_scheme();
    }
  }

  std::vector<double> bin_list(quad_list_files.size());
  /* We will generate this by bin for each map so make the bin number the
   * first index. */
  std::vector<std::vector<double> > Corr(quad_list_files.size());

#pragma omp parallel shared(Corr, bin_list, quad_list_files, maps)
  {
    Npoint_Functions::Quadrilateral_List_File<int> qlf;
    
#pragma omp for schedule(dynamic,2)
    for (size_t k=0; k < quad_list_files.size(); ++k) {
      if (! qlf.initialize (quad_list_files[k])) {
	std::cerr << "Error initializing quadrilateral list from "
		  << quad_list_files[k] << std::endl;
	std::exit(1);
      }

      bin_list[k] = qlf.bin_value();
      Npoint_Functions::calculate_fourpoint_function_list
	(maps, qlf, Corr[k]);
    }
  }
  
  std::cout << "# LCDM four point function from " << quad_list_prefix
	    << std::endl;
  std::cout << "# First line is bin values, rest are the four point function.\n";
  for (size_t k=0; k < bin_list.size(); ++k) {
    std::cout << bin_list[k] << " ";
  }
  std::cout << std::endl;

  for (size_t j=0; j < maps.size(); ++j) {
    for (size_t k=0; k < bin_list.size(); ++k) {
      std::cout << Corr[k][j] << " ";
    }
    std::cout << std::endl;
  }

  return 0;
}
