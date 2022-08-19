// This class follows Munuit 2 instructions to fit a
// straight line through a set of skew lines.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef my_straight_line_fcn_h
#define my_straight_line_fcn_h

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
    class my_straight_line_fcn : public FCNBase
    {
    public:
      my_straight_line_fcn(const std::vector<line_3d>& strips, const double midpoint_z)
      : cluster_lines(strips), the_midpoint_z(midpoint_z), fErrorDef(1.) {}
      
      ~my_straight_line_fcn() {}
      
      virtual double Up() const {return fErrorDef;}
      virtual double operator()(const std::vector<double>&) const;
      
      void SetErrorDef(double def) {fErrorDef = def;}
      
    private:
      std::vector<line_3d> cluster_lines;
      double the_midpoint_z;
      
      double fErrorDef;
    };
    
  }
  
}

#endif
