#include <iostream>
#include <iomanip>
#include <string>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Pixel_Quadrilaterals.h>
#include <Npoint_Functions_Utils.h>

namespace {
  const std::string CREATE_RHOMBIC_QUADRILATERALS_LIST_RCSID
  ("$Id: create_twopt_table.cpp,v 1.12 2011-08-09 21:08:04 copi Exp $");
}

void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <two point table prefix> "
	    << "<two point table file number>\n";
  exit (0);
}

int main (int argc, char *argv[])
{
  if (argc != 3) usage (argv[0]);

  std::string twopt_prefix = argv[1];
  int filenum;
  if (! Npoint_Functions::from_string (argv[2], filenum)) {
    std::cerr << "Could not read filenumber from '" << argv[2] << "'\n";
    usage (argv[0]);
  }

  Npoint_Functions::Twopt_Table<int> twopt_table;
  twopt_table.read_file (Npoint_Functions::make_filename (twopt_prefix,
							  filenum)); 
  Npoint_Functions::Pixel_Triangles_Equilateral<int> triangles;
  std::vector<int> tri;
  std::vector<int> thirdpt;
  thirdpt.reserve(100);
  Npoint_Functions::Pixel_Quadrilaterals_Rhombic_Full<int> q;
  triangles.find_triangles (twopt_table);
  q.initialize (triangles);
  while (q.next(tri, thirdpt)) {
    for (size_t j=0; j < thirdpt.size(); ++j) {
      for (size_t i=0; i < tri.size(); ++i)
	std::cout << tri[i] << " ";
      std::cout << thirdpt[j] << std::endl;
    }
  }
  return 0;
}
