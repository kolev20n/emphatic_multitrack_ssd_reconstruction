// This class follows Munuit 2 instructions to fit a
// straight line through a point and a set of skew lines.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef my_straight_line_fcn_with_point_h
#define my_straight_line_fcn_with_point_h

#include "Minuit2/FCNBase.h"
#include <TMath.h>
#include <Math/Vector3D.h>

#include <vector>

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

namespace ROOT
{
  namespace Minuit2
  {
    class my_straight_line_fcn_with_point : public FCNBase
    {
    public:
      my_straight_line_fcn_with_point(const std::vector<line_3d>& strips, const line_3d a_guessed_vertex, const double a_guessed_vertex_weight, const double midpoint_z)
      : cluster_lines(strips), the_guessed_vertex(a_guessed_vertex), the_guessed_vertex_weight(a_guessed_vertex_weight), the_midpoint_z(midpoint_z), fErrorDef(1.) {}
      
      ~my_straight_line_fcn_with_point() {}
      
      virtual double Up() const {return fErrorDef;}
      virtual double operator()(const std::vector<double>&) const;
      
      void SetErrorDef(double def) {fErrorDef = def;}
      
    private:
      std::vector<line_3d> cluster_lines;
      line_3d the_guessed_vertex;
      double the_guessed_vertex_weight;
      double the_midpoint_z;

      double fErrorDef;
    };
    
  }
  
}

#endif
