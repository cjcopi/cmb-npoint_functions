#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <string_utils.h> // For parsing a text file
#include <paramfile.h>
#include <vec3.h>

#include <Twopt_Table.h>
#include <buffered_pair_binary_file.h>
#include <Npoint_Functions_Utils.h>

// Documentation read by doxygen for the front page of the project.

/** \mainpage Npoint Functions
 *
 *  The most common npoint function calculated for use in the CMB and
 *  elsewhere is the two point function, also known as the correlation
 *  function, \f$C(\theta)\f$.  This function can be quickly calculated in
 *  a number of ways, in particular using PolSpice.  Higher point functions
 *  cannot efficiently be calculated using such clever techniques.
 *
 *  A pixel based method has been described by Eriksen, H.K., Lilje, P.B.,
 *  Banday, A.J, and Gorski, K.M., "Estimating N-point correlation
 *  functions from pixelized sky maps", ApJS, 151, 1 (2004).  The code
 *  presented here is an independent implementation of the basic idea with
 *  some of the optimizations presented in their paper.  Some other
 *  optimizations have also been included (such as compressing the two
 *  point tables and using temporary files when first calculating the two
 *  point tables).  The added optimizations allow for larger number of
 *  pixels to be employed.
 */
namespace {
  const std::string CREATE_TWOPT_TABLE_RCSID
  ("$Id: create_twopt_table.cpp,v 1.11 2011-07-20 18:34:01 copi Exp $");
}

void mask_to_pixlist (const Healpix_Map<double>& mask,
		      std::vector<int>& pixlist)
{
  pixlist.clear();
  for (int j=0; j < mask.Npix(); ++j) {
    if (mask[j] > 0.5) pixlist.push_back(j);
  }
}

bool read_text_file (const std::string& cosbinfile,
		     std::vector<double>& bin_list)
{
  /* Read the file line by line and extract the first column.
   *  Anything following a # is a comment.
   */
  std::string line;
  std::string::iterator it;
  std::vector<double> vals;
  std::ifstream in (cosbinfile.c_str());
  if (! in.is_open()) return false;

  bin_list.clear();
  while (in.good()) {
    std::getline (in, line);
    line = trim (line);
    it = std::find (line.begin(), line.end(), '#');
    if (it != line.end()) line.erase (it, line.end());
    if (line == "") continue;
    vals.clear();
    split (line, vals);
    bin_list.push_back(vals[0]);
  }
  in.close();
  return true;
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
  double dcosbin = params.find<double> ("dcosbin", -100);
  double dtheta = params.find<double> ("dtheta", -200);
  std::string cosbinfile = params.find<std::string> ("cosbinfile", "");
  std::string tmpfile_prefix = params.find<std::string> ("tmpfile_prefix");
  std::string twoptfile_prefix = params.find<std::string> ("twoptfile_prefix");
  bool clean_tmpfiles = params.find<bool> ("clean_tmpfiles", false);

  if ((Nside == -1) && (maskfile == "")) {
    std::cerr << "Maskfile or Nside must be set in the parameter file.\n";
    return 1;
  }

  if ((dcosbin == -100) && (cosbinfile == "") && (dtheta == -200)) {
    std::cerr << "cosbinfile or dcosbin or dtheta must be set in the parameter file.\n";
    return 1;
  }

  Healpix_Map<double> mask;
  if (maskfile != "") {
    read_Healpix_map_from_fits (maskfile, mask);
    if (mask.Scheme() == RING) mask.swap_scheme();
    Nside = mask.Nside();
    mask_to_pixlist (mask, pixel_list);
  } else {
    pixel_list.resize (12*Nside*Nside);
    std::generate (pixel_list.begin(), pixel_list.end(),
		   Npoint_Functions::myRange<int>());
  }

  std::vector<double> cosbin;
  if (cosbinfile != "") {
    if (! read_text_file (cosbinfile, bin_list)) {
      std::cerr << "Failed reading " << cosbinfile << std::endl;
      return 1;
    }
  } else if (dcosbin != -100) {
    int Nbin = 2/dcosbin;
    bin_list.resize(Nbin);
    std::generate (bin_list.begin(), bin_list.end(),
		   Npoint_Functions::myRange<double>(-1.0+dcosbin/2,
						     dcosbin));
  } else {
    int Nbin = 180/dtheta;
    bin_list.resize(Nbin);
    /* Run this "backward" since we will use bins in cos(theta) and
     * cos(180)=-1. */
    std::generate (bin_list.begin(), bin_list.end(),
		   Npoint_Functions::myRange<double>(180-dtheta/2,
						     -dtheta));
    /* We want equal spacing/width in theta (I guess) so
     * create cosbin here with this in mind. */
    cosbin.push_back(-1.1);
    for (size_t j=0; j < bin_list.size()-1; ++j) {
      cosbin.push_back(cos(0.5*(bin_list[j]+bin_list[j+1])*M_PI/180));
    }
    cosbin.push_back(1.1);
  }

  /* Convert bin_list to bin edges. Put the ends of the bins a little off
   * the min and max values so we don't have to deal with the dot product
   * numerically being a little too big or small.
   */  
  if (cosbin.size() == 0) {
    cosbin.push_back(-1.1);
    for (size_t j=0; j < bin_list.size()-1; ++j) {
      cosbin.push_back(0.5*(bin_list[j]+bin_list[j+1]));
    }
    cosbin.push_back(1.1);
  }

  size_t Npix = pixel_list.size();
  std::cout << "Generating for\n Nside = " << Nside
	    << "\n Npix = " << Npix
	    << "\n Nbin = " << bin_list.size()
	    << std::endl;

  Healpix_Base HBase (Nside, NEST, SET_NSIDE);
  // Create list of vectors.
  std::vector<vec3> veclist (Npix);
  for (size_t i=0; i < Npix; ++i) {
    veclist[i] = HBase.pix2vec (pixel_list[i]);
  }

  std::vector<Npoint_Functions::buffered_pair_binary_file<int> > binfiles;
  for (size_t k=0; k < bin_list.size(); ++k) {
    binfiles.push_back(Npoint_Functions::buffered_pair_binary_file<int>
		       (Npoint_Functions::make_filename
			(tmpfile_prefix, k))); 
    binfiles[k].create();
  }
  
  std::cout << "Creating temporary files.\n";
  /* Now create values and write to temporary files the files.
   * This could be parallelized but not easily.  We use the order we step
   * through the  pixels to ensure that the tables we create below are
   * sorted withouthaving to actually run a sorting algorithm on them.
   */
  vec3_t<double> v0, v1;
  size_t ibin = 0;
  double dp;
  int dir;
  for (size_t i=0; i < Npix; ++i) {
    for (size_t j=i+1; j < Npix; ++j) {
      dp = dotprod (veclist[i], veclist[j]);
      //dp = std::max (dp, -1.0);
      //dp = std::min (dp,  1.0);
      if (dp > cosbin[ibin]) dir = +1;
      else dir = -1;
      /* Find the bin.  Since we use the NEST scheme a simple linear search
       * is efficient; sequential pixels are near each other so it should
       * be a short walk between pixel pairs.  This is only true if the
       * pixel list is sorted.  If the pixel list is randomized then this
       * won't be efficient.  Even so, the number of bins is expected to be
       * small so a more sophisticated algorithm isn't warranted.
       *
       * We REQUIRE the bin list to be inclusive, ie start at -1 (or
       * smaller) and end at 1 (or larger).
       */
      while ((dp < cosbin[ibin]) || (dp > cosbin[ibin+1]))
	ibin += dir;
      binfiles[ibin].append(i, j);
    }
  }
  /* Free memory.  This should flush buffers, close files, and release
   * allocated memory.
   */
  binfiles.clear();
  std::cout << "Temporary files created.\n";

  std::cout << "Creating two point tables.\n";
  // Now create the 2 point tables.   This can trivially be parallelized.
#pragma omp parallel shared(binfiles, Npix, pixel_list, bin_list, \
  tmpfile_prefix, twoptfile_prefix)
  {
    Npoint_Functions::Twopt_Table<int>
      twopt_table (Nside, pixel_list, bin_list[0]);

    int i, j;
#pragma omp for schedule(guided)
    for (size_t k=0; k < bin_list.size(); ++k) {
      twopt_table.reset();
      twopt_table.bin_value (bin_list[k]);
      Npoint_Functions::buffered_pair_binary_file<int>
	binfile(Npoint_Functions::make_filename (tmpfile_prefix, k));

      // Next open the file for reading
      binfile.open_read();
      // Now fill in the table by looping over all pairs of pixels.
      while (binfile.read_next_pair (i, j)) {
	twopt_table.add_pair (i, j);
      }
      if (clean_tmpfiles) unlink(binfile.filename().c_str());

      twopt_table.write_file 
	(Npoint_Functions::make_filename (twoptfile_prefix, k));

    }
  }
  std::cout << "Two point tables created.\n";

  return 0;
}
