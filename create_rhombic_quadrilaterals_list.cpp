#include <iostream>
#include <iomanip>
#include <string>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Pixel_Quadrilaterals.h>

namespace {
  const std::string CREATE_RHOMBIC_QUADRILATERALS_LIST_RCSID
  ("$Id: create_rhombic_quadrilaterals_list.cpp,v 1.3 2016/02/09 20:31:44 copi Exp $");
}

void usage (const char *progname)
{
  std::cerr << "Usage: " << progname << " <two point table name>\n";
  exit (0);
}

int main (int argc, char *argv[])
{
  if (argc != 2) usage (argv[0]);

  std::string twopt_table_file = argv[1];

  Npoint_Functions::Twopt_Table<int> twopt_table;
  twopt_table.read_file (twopt_table_file);                                        
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
