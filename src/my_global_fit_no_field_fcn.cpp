#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>

#include "my_global_fit_no_field_fcn.h"

using namespace std;
using namespace ROOT::Math;

namespace ROOT
{
  namespace Minuit2
  {
    double my_global_fit_no_field_fcn::operator()(const vector<double>& parameters) const
    {
      XYZVector vertex;
      XYZVector track_direction;

      double x0, y0, z0, x, y, z;
      double px, py, pz;
      
      double chi2 = 0.;

      double dummy_double;

      XYZVector shortest_distance_vector;

      x0 = the_upstream_track.track_line.point.X();
      y0 = the_upstream_track.track_line.point.Y();
      z0 = the_upstream_track.track_line.point.Z();

      z = parameters[0];
      
      x = x0 + (z - z0) * the_upstream_track.track_line.direction.X() / the_upstream_track.track_line.direction.Z();
      y = y0 + (z - z0) * the_upstream_track.track_line.direction.Y() / the_upstream_track.track_line.direction.Z();

      vertex.SetCoordinates(x, y, z);
      
      for (int i = 1; i < parameters.size() - 1; i += 2)
      {
        px = parameters[i];
        py = parameters[i + 1];
        pz = 1.;
      
        track_direction.SetCoordinates(px, py, pz);
        track_direction = track_direction.Unit();
      
        for (int j = 0; j < the_pattern.combined_tracks.at((i - 1) / 2).cluster_lines.size(); j++)
        {
          shortest_distance_vector = the_pattern.combined_tracks.at((i - 1) / 2).cluster_lines.at(j).direction.Cross(track_direction);
          
          dummy_double = fabs(shortest_distance_vector.Dot(the_pattern.combined_tracks.at((i - 1) / 2).cluster_lines.at(j).point - vertex)) / sqrt(shortest_distance_vector.Mag2());
          
          chi2 += dummy_double;
        }
        
      }
      
      return chi2;
    }
    
  }  // end namespace Minuit2
  
}  // end namespace ROOT
