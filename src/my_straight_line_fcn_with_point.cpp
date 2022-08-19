#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

#include "my_straight_line_fcn_with_point.h"

using namespace std;
using namespace ROOT::Math;

namespace ROOT
{
  namespace Minuit2
  {
    double my_straight_line_fcn_with_point::operator()(const vector<double>& parameters) const
    {
      double dummy_double;
      
      double px, py, pz; // unit vector of direction
      
      double chi2 = 0.;
      
      XYZVector dummy_vector;
      
      XYZVector shortest_distance_vector;
      
      XYZVector track_direction;

      assert(parameters.size() == 4);
      
      px = parameters[2];
      py = parameters[3];
      pz = 1.;
      
      track_direction.SetCoordinates(px, py, pz);
      track_direction = track_direction.Unit();
      
      dummy_vector.SetCoordinates(parameters[0], parameters[1], the_midpoint_z);
      
      for (int i = 0; i < cluster_lines.size(); i++)
      {
        shortest_distance_vector = cluster_lines[i].direction.Cross(track_direction);
        
        dummy_double = fabs(shortest_distance_vector.Dot(cluster_lines[i].point - dummy_vector)) / sqrt(shortest_distance_vector.Mag2());
        
        chi2 += dummy_double;
      }
      
      chi2 += the_guessed_vertex_weight * sqrt(track_direction.Cross(the_guessed_vertex.point - dummy_vector).Mag2());
 
      //cout << chi2 << endl;
      
      return chi2;
    }
    
    
  }  // namespace Minuit2
  
}  // namespace ROOT
