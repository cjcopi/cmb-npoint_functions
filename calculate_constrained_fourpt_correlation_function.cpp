#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <alm.h>
#include <alm_healpix_tools.h>
#include <alm_fitsio.h>

#include <dirtree.h>

#include <Quadrilateral_List_File.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CALCULATE_constrained_FOURPT_CORRELATION_FUNCTION_RCSID
  ("$Id: calculate_constrained_fourpt_correlation_function.cpp,v 1.3 2016/02/09 20:31:44 copi Exp $");
}


void usage (const char *progname)
{
  std::cerr << "Usage: " << progname
            << " <quad list prefix> <Alm dir>"
            << " <num alm start> <num alm end>"
            << " [<mask file>>]\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if ((argc < 5) || (argc > 6)) usage (argv[0]);
  std::string quad_list_prefix = argv[1];
  std::string alm_dir = argv[2];
  size_t Nstart, Nend;
  if (! Npoint_Functions::from_string (argv[3], Nstart)) {
    std::cerr << "Could not parse Nstart\n";
    usage (argv[0]);
  }
  if (! Npoint_Functions::from_string (argv[4], Nend)) {
    std::cerr << "Could not parse Nend\n";
    usage (argv[0]);
  }
  bool have_mask = false;
  Healpix_Map<double> mask;
  if (argc == 6) {
    read_Healpix_map_from_fits (argv[5], mask);
    have_mask = true;
  }

  // Figure out how many bins there are by trying to open files.
  std::vector<std::string> quad_list_files
    = Npoint_Functions::get_range_file_list(quad_list_prefix, 0, 400);
  if (quad_list_files.size() == 0) {
    std::cerr << "No quad list files found!\n";
    usage (argv[0]);
  }

  int Lmax;
  std::vector<Healpix_Map<double> > maps (Nend-Nstart);
  // Make maps
  {
    Npoint_Functions::Quadrilateral_List_File<int> qlf;
    qlf.initialize (quad_list_files[0]);
    if (have_mask) {
      if (static_cast<size_t>(mask.Nside()) != qlf.Nside()) {
        std::cerr << "Mask and quadrilateral lists do not have"
                  << " the same Nside: " << mask.Nside() 
                  << " != " << qlf.Nside() << std::endl;
        std::exit(1);
      }
      if (mask.Scheme() != qlf.Scheme()) mask.swap_scheme();
    }
    Lmax = std::min(200UL, 4*qlf.Nside()+1);
    //#pragma omp parallel shared(qlf, maps)
    {
      Alm<xcomplex<double> > alm (Lmax, Lmax);
      //#pragma omp for schedule(static)
      for (size_t k=0; k < maps.size(); ++k) {
        read_Alm_from_fits (dirtree::filename(alm_dir, "alm_T_", ".fits",
                                              k+Nstart),
                            alm, Lmax, Lmax);
        maps[k].SetNside (qlf.Nside(), RING);
        alm2map (alm, maps[k]);
        if (maps[k].Scheme() != qlf.Scheme()) maps[k].swap_scheme();
      }
    }
  }

  std::vector<double> bin_list(quad_list_files.size());
  /* We will generate this by bin for each map so make the bin number the
   * first index. */
  std::vector<std::vector<double> > Corr(quad_list_files.size());

#pragma omp parallel shared(Corr, bin_list, quad_list_files, maps, mask)
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
      if (have_mask) {
        Npoint_Functions::calculate_masked_fourpoint_function_list
          (maps, mask, qlf, Corr[k]);
      } else {
        Npoint_Functions::calculate_fourpoint_function_list
          (maps, qlf, Corr[k]);
      }
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
