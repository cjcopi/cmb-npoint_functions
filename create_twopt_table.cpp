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
#include <string_utils.h> // For parsing a text file
#include <paramfile.h>
#include <vec3.h>

#include <Twopt_Table.h>
#include <buffered_pair_binary_file.h>
#include <myRange.h>

namespace {
  const std::string CREATE_TWOPT_TABLE_RCSID
  ("$Id: create_twopt_table.cpp,v 1.6 2011-07-09 05:07:01 copi Exp $");
}

void dtheta_to_cosbin (double dtheta, std::vector<double>& cosbin)
{
  dtheta *= M_PI / 180.0; // Convert from degrees to radians
  size_t Nbin = M_PI/dtheta;
  cosbin.clear();
  double theta = M_PI;
  while (theta > 0) {
    cosbin.push_back(std::cos(theta));
    theta -= dtheta;
  }
}
void dcostheta_to_cosbin (double dcostheta, std::vector<double>& cosbin)
{
  size_t Nbin = 2/dcostheta;
  cosbin.clear();
  double cb = -1;
  while (cb < 1) {
    cosbin.push_back(cb);
    cb += dcostheta;
  }
}
void mask_to_pixlist (const Healpix_Map<double>& mask,
		      std::vector<int>& pixlist)
{
  pixlist.clear();
  for (size_t j=0; j < mask.Npix(); ++j) {
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
  // Make sure the list ends with the 1.0 bin
  if (abs(bin_list[bin_list.size()-1] - 1) > 1.0e-12)
    bin_list.push_back (1.0);
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
  double dtheta = params.find<double> ("dtheta", -200);
  double dcostheta = params.find<double> ("dcostheta", -200);
  std::string cosbinfile = params.find<std::string> ("cosbinfile", "");
  std::string tmpfile_prefix = params.find<std::string> ("tmpfile_prefix");
  std::string twoptfile_prefix = params.find<std::string> ("twoptfile_prefix");
  bool clean_tmpfiles = params.find<bool> ("clean_tmpfiles", false);

  if ((Nside == -1) && (maskfile == "")) {
    std::cerr << "Maskfile or Nside must be set in the parameter file.\n";
    return 1;
  }

  if ((dtheta == -200) && (dcostheta == -200) && (cosbinfile == "")) {
    std::cerr << "A source of bins must be provided in the parameter file\n"
	      << "  Set cosbinfile, dtheta, or dcostheta.\n";
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
    std::generate (pixel_list.begin(), pixel_list.end(), myRange<int>());
  }

  if (cosbinfile != "") {
    if (! read_text_file (cosbinfile, bin_list)) {
      std::cerr << "Failed reading " << cosbinfile << std::endl;
      return 1;
    }
  } else if (dcostheta != -200) {
    dcostheta_to_cosbin (dcostheta, bin_list);
  } else {
    dtheta_to_cosbin (dtheta, bin_list);
  }

  size_t Npix = pixel_list.size();
  size_t Nbin = bin_list.size() - 1;
  std::cout << "Generating for\n Nside = " << Nside
	    << "\n Npix = " << Npix
	    << "\n Nbin = " << Nbin
	    << std::endl;

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
  /* Now create values and write to temporary files the files.
   * This could be parallelized but not easily.  We use the order we step
   * through the  pixels to ensure that the tables we create below are
   * sorted withouthaving to actually run a sorting algorithm on them.
   */
  vec3_t<double> v0, v1;
  int ibin = 0;
  double dp;
  int dir;
  for (size_t i=0; i < Npix; ++i) {
    for (size_t j=i+1; j < Npix; ++j) {
      dp = dotprod (veclist[i], veclist[j]);
      dp = std::max (dp, -1.0);
      dp = std::min (dp,  1.0);
      if (dp > bin_list[ibin]) dir = +1;
      else dir = -1;
      /* Find the bin.  Since we use the NEST scheme a simple linear search
       *  is efficient; sequential pixels are near each other so it should
       *  be a short walk between pixel pairs.  This is only true if the
       *  pixel list is sorted.  If the pixel list is randomized then this
       * won't be efficient.  Even so, the number of bins is expected to be
       * small so a more sophisticated algorithm isn't warranted.
       *
       *  We do NOT require the bin list to be inclusive, ie start at -1
       *  (or smaller) and end at 1 (or larger).  If a value is "off the
       *  end" of the list we stick the value in the first or last bin, as
       *  appropriate.
       */
      while ( (((ibin != 0) && (dir < 0)) || ((ibin != Nbin-1) && (dir > 0)))
	      && ((dp < bin_list[ibin]) || (dp > bin_list[ibin+1])))
	ibin += dir;
      binfiles[ibin].append(i, j);
    }
  }
  /* Free memory.  This should flush buffers, close files, and release
   * allocated memory.
   */
  binfiles.clear();
  std::cout << "Temporary files created.\n";

  // Now create the 2 point tables.   This can trivially be parallelized.
#pragma omp parallel shared(binfiles, Nbin, Npix, pixel_list, bin_list, \
  tmpfile_prefix, twoptfile_prefix)
  {
    Twopt_Table<int> twopt_table (pixel_list, bin_list[0]);

    int i, j;
    std::ostringstream sstr;
#pragma omp for schedule(guided)
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
