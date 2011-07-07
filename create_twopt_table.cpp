#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <sstream>

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <paramfile.h>
#include <vec3.h>

#include <Twopt_Table.h>
#include <buffered_pair_binary_file.h>
#include <myRange.h>

namespace {
  const std::string CREATE_TWOPT_TABLE_RCSID
  ("$Id$");
}

void dtheta_to_cosbin (double dtheta, std::vector<double>& cosbin)
{
  dtheta *= M_PI / 180.0; // Convert from degrees to radians
  size_t Nbin = M_PI/dtheta;
  cosbin.resize(Nbin);
  double theta = M_PI;
  for (size_t j=0; j < cosbin.size(); ++j) {
    cosbin[j] = std::cos(theta);
    theta -= dtheta;
  }
  /* We MUST have the 1.0 bin for the bin finding algorithm below
   * even though nothing can be in it.
   */
  if (std::abs(cosbin[cosbin.size()-1] - 1) > 1.0e-12) cosbin.push_back(1.0);
}
void mask_to_pixlist (const Healpix_Map<double>& mask,
		      std::vector<int>& pixlist)
{
  pixlist.clear();
  for (size_t j=0; j < mask.Npix(); ++j) {
    if (mask[j] > 0.5) pixlist.push_back(j);
  }
}

void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <parameter file name>\n";
  exit (1);
}


int main (int argc, char *argv[])
{
  if (argc != 2) usage (argv[0]);
  std::vector<double> bin_list;
  std::vector<int> pixel_list;
  
  paramfile params (argv[1]);
  int Nside = params.find<int> ("Nside", -1);
  std::string maskfile = params.find<std::string> ("maskfile", "");
  double dtheta = params.find<double> ("dtheta");
  std::string tmpfile_prefix = params.find<std::string> ("tmpfile_prefix");
  std::string twoptfile_prefix = params.find<std::string> ("twoptfile_prefix");
  bool clean_tmpfiles = params.find<bool> ("clean_tmpfiles", false);

  if ((Nside == -1) && (maskfile == "")) {
    std::cerr << "Maskfile or Nside must be set in the parameter file.\n";
    return 1;
  }

  Healpix_Map<double> mask;
  if (maskfile != "") {
    read_Healpix_map_from_fits ("mask_r7_kq75y7.fits", mask);
    if (mask.Scheme() == RING) mask.swap_scheme();
    Nside = mask.Nside();
    mask_to_pixlist (mask, pixel_list);
  } else {
    pixel_list.resize (12*Nside*Nside);
    std::generate (pixel_list.begin(), pixel_list.end(), myRange<int>());
  }

  dtheta_to_cosbin (dtheta, bin_list);

  size_t Npix = pixel_list.size();
  size_t Nbin = bin_list.size() - 1;

  Healpix_Base HBase (Nside, NEST, SET_NSIDE);
  // Create list of vectors.
  std::vector<vec3> veclist (Npix);
  for (size_t i=0; i < Npix; ++i) {
    veclist[i] = HBase.pix2vec (pixel_list[i]);
  }

  std::vector<buffered_pair_binary_file<int> > binfiles;
  {
    std::ostringstream fstr;
    for (int k=0; k < Nbin; ++k) {
      fstr.str("");
      fstr << tmpfile_prefix << std::setw(5) << std::setfill('0') << k << ".dat";
      binfiles.push_back(buffered_pair_binary_file<int> (fstr.str()));
      binfiles[k].create();
    }
  }
  
  std::cout << "Creating temporary files.\n";
  // Now create values and write to temporary files the files.
  vec3_t<double> v0, v1;
  int ibin = 0;
  double dp;
  int dir;
  for (size_t i=0; i < Npix; ++i) {
    if (i%10000==0) std::cout << i << std::endl;
    for (size_t j=i+1; j < Npix; ++j) {
      dp = dotprod (veclist[i], veclist[j]);
      dp = std::max (dp, -1.0);
      dp = std::min (dp,  1.0);
      if (dp > bin_list[ibin]) dir = +1;
      else dir = -1;
      while ((dp < bin_list[ibin]) || (dp > bin_list[ibin+1])) ibin += dir;
      binfiles[ibin].append(i, j);
    }
  }
  /* Free memory.  This should flush buffers, close files, and release
   * allocated memory.
   */
  binfiles.clear();
  std::cout << "Temporary files created.\n";

  /* Now create the 2 point tables.  For each table these are 2 dim arrays
   * of size Npix x ?, indexed as (p,?) where p is the pixel at the center.
   * The entries in the table are the pixels paired with p that are in the
   * given bin. Note that this table may not be completely full if dtheta
   * is small compared to Nside. */
  std::cout << "Creating two point tables.\n";

#pragma omp parallel shared(binfiles, Nbin, Npix, pixel_list, bin_list, \
  tmpfile_prefix, twoptfile_prefix)
  {
    Twopt_Table<int> twopt_table (pixel_list, bin_list[0]);

    int i, j;
    std::ostringstream sstr;
#pragma omp for schedule(static)
    for (size_t k=0; k < Nbin; ++k) {
      twopt_table.reset();
      twopt_table.bin_value (bin_list[k]);
      sstr.str("");
      sstr << tmpfile_prefix << std::setw(5) << std::setfill('0') << k << ".dat";
      buffered_pair_binary_file<int> binfile(sstr.str());

      // Next open the file for reading
      binfile.open_read();
      // Now fill in the table by looping over all pairs of pixels.
      while (binfile.read_next_pair (i, j)) {
	twopt_table.add_pair (i, j);
      }

      if (clean_tmpfiles) unlink(sstr.str().c_str());

      sstr.str("");
      sstr << twoptfile_prefix << std::setw(5) << std::setfill('0') << k << ".dat";
      twopt_table.write_file (sstr.str());
    }
  }
  std::cout << "Two point tables created.\n";

  return 0;
}
