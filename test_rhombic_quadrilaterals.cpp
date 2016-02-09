#include <iostream>
#include <iomanip>
#include <string>

#include <Twopt_Table.h>
#include <Pixel_Triangles.h>
#include <Pixel_Quadrilaterals.h>
#include <Npoint_Functions_Utils.h>

// $Id$

int main ()
{
  std::string twopt_prefix ("data/twopt_Nside32_0.01_");

  Npoint_Functions::Twopt_Table<int> twopt_table;
  twopt_table.read_file (Npoint_Functions::make_filename (twopt_prefix,
							  150)); 
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
